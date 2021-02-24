//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MPMRayStudyBase.h"

// MOOSE includes
#include "Assembly.h"
#include "AuxiliarySystem.h"

InputParameters
MPMRayStudyBase::validParams()
{
  auto params = RayTracingStudy::validParams();

  params.addRequiredCoupledVar("displacement",
                               "The displacement variables for moving the particles");
  params.addRequiredCoupledVar("mass", "The mass variable");

  // Ray tracing by default requires that we set names for each Ray - we don't need them named here
  params.set<bool>("_use_ray_registration") = false;

  // Whether or not the particles added during generateInitialParticles() are replicated
  params.addPrivateParam<bool>("_initial_particles_replicated", true);

  return params;
}

MPMRayStudyBase::MPMRayStudyBase(const InputParameters & parameters)
  : RayTracingStudy(parameters),
    CoupleableMooseVariableDependencyIntermediateInterface(this, false),
    _particle_mass_index(registerRayAuxData("particle_mass")),
    _initial_particles_replicated(getParam<bool>("_initial_particles_replicated")),
    _displacements(coupledValues("displacement")),
    _mass(coupledValue("mass")),
    _claim_rays(*this,
                _mesh,
                _initial_particles,
                _initial_local_particles,
                !getParam<bool>("_initial_particles_replicated")),
    _particles(
        libmesh_make_unique<std::unordered_map<const Elem *, std::vector<std::shared_ptr<Ray>>>>()),
    _old_particles(
        libmesh_make_unique<std::unordered_map<const Elem *, std::vector<std::shared_ptr<Ray>>>>()),
    _generated_initial_particles(false),
    _generating_initial_particles(false)
{
  if (coupledComponents("displacement") != _mesh.dimension())
    paramError("displacements", "Must provide as many displacements as the mesh dimension");

  // Get the displacement variables for reinit
  for (unsigned int i = 0; i < _mesh.dimension(); ++i)
  {
    auto var = getVar("displacement", i);
    auto cast_var = dynamic_cast<MooseVariableFEBase *>(var);
    mooseAssert(cast_var, "Cannot be cast to MooseVariableFEBase");
    _displacement_vars.insert(cast_var);
  }

  // Get the mass variable for reinit
  {
    auto var = getVar("mass", 0);
    _mass_var = dynamic_cast<MooseVariableFEBase *>(var);
    mooseAssert(_mass_var, "Cannot be cast to MooseVariableFEBase");
  }
}

void
MPMRayStudyBase::postExecuteStudy()
{
  // _particles contains the new binned rays, move them to _old_particles
  // so that we can store into _particles again on the next iteration
  std::swap(_particles, _old_particles);
  _particles->clear();

  // Zero the mass variable, as we will be filling into it
  std::vector<std::string> zero_names = {_mass_var->name()};
  _fe_problem.getAuxiliarySystem().zeroVariables(zero_names);

  // Set the mass variable as active so that we can reinit it
  _fe_problem.setActiveElementalMooseVariables({_mass_var}, 0);

  // Cast the mass variable to a field variable so we can access phi
  MooseVariableField<Real> * mass_var_fe = dynamic_cast<MooseVariableField<Real> *>(_mass_var);
  mooseAssert(mass_var_fe, "Not a field variable");

  // Loop through the elem -> rays (particles) map to set the new displacements
  // and move into the buffer to be traced
  for (auto & pair : *_old_particles)
  {
    const Elem * elem = pair.first;
    std::vector<std::shared_ptr<Ray>> & rays = pair.second;

    _fe_problem.setCurrentSubdomainID(elem, 0);
    _fe_problem.prepare(elem, 0);
    mass_var_fe->prepareAux();

    // The contribution to this element from the particles
    DenseVector<Real> contribution(_mass_var->phiSize(), 0);

    for (auto & ray : rays)
    {
      // For phi reinit at the particle point
      _fe_problem.reinitElemPhys(elem, {ray->currentPoint()}, 0);

      // Append to the contribution at each node for this particle
      for (std::size_t i = 0; i < _mass_var->phiSize(); ++i)
        contribution(i) += ray->auxData(_particle_mass_index) * mass_var_fe->phi()[i][0];
    }

    // Add the contribution to the mass variable
    mass_var_fe->setDofValues(contribution);
    mass_var_fe->add(_fe_problem.getAuxiliarySystem().solution());
  }

  // Done with the reinits and aux solution access
  _fe_problem.clearActiveElementalMooseVariables(0);
  _fe_problem.getAuxiliarySystem().solution().close();
}

void
MPMRayStudyBase::onCompleteRay(const std::shared_ptr<Ray> & ray)
{
  // When a Ray is complete, it has traversed through its displacement for this iteration
  // Now we bin it into the element it ended up in for postprocessing after all rays are done
  (*_particles)[ray->currentElem()].push_back(ray);
}

void
MPMRayStudyBase::generateRays()
{
  // On the first run around, we need to initialize the particles at the starting
  // positions determined by generateInitialParticles(), which fills the initial
  // rays into _old_particles
  if (!_generated_initial_particles)
  {
    _generated_initial_particles = true;
    _generating_initial_particles = true;

    // Overridden in derived classes to generate the initial particles into _initial_particles
    generateInitialParticles();
    _generating_initial_particles = false;

    // Claim the initial particles - at this point, they only have starting points set. We need
    // them on the right processors and with the right starting elems. This fills from
    // _initial_particles into _initial_local_particles (the particles this proc starts)
    _claim_rays.claim();

    // Take the initial rays (particles), and bin by element
    _old_particles->clear();
    for (auto & ray : _initial_local_particles)
      (*_old_particles)[ray->currentElem()].emplace_back(std::move(ray));

    // Done with these
    _initial_particles.clear();
    _initial_local_particles.clear();
  }

  // Helper for getting the displacement at a point within an element
  auto get_displacement = [&](const Elem * elem, const Point & point) {
    _fe_problem.setCurrentSubdomainID(elem, 0);
    if (_fe_problem.assembly(0).elem() != elem)
      _fe_problem.prepare(elem, 0);

    _fe_problem.reinitElemPhys(elem, {point}, 0);

    Point displacement(0, 0, 0);
    for (std::size_t d = 0; d < _mesh.dimension(); ++d)
      displacement(d) = (*_displacements[d])[0];

    return displacement;
  };

  // Set the active variables to the displacement variables so we can reinit them
  _fe_problem.setActiveElementalMooseVariables(_displacement_vars, 0);

  // Loop through the elem -> rays (particles) map to set the new displacements
  // and move into the buffer to be traced
  for (auto & pair : *_old_particles)
  {
    const Elem * elem = pair.first;
    std::vector<std::shared_ptr<Ray>> & rays = pair.second;

    for (auto & ray : rays)
    {
      const auto point = ray->currentPoint();
      const auto displacement = get_displacement(elem, point);

      // Must clear because we're reusing an old ray
      ray->resetCounters();
      ray->clearStartingInfo();

      // Set the start and displacement for the next iteration
      ray->setStart(point, elem);
      ray->setStartingEndPoint(point + displacement);

      // Move into the buffer to be traced
      moveRayToBuffer(ray);
    }
  }

  // Done reiniting the displacement variables
  _fe_problem.clearActiveElementalMooseVariables(0);
}

const std::shared_ptr<Ray> &
MPMRayStudyBase::generateInitialParticle(const Point & point, const Elem * elem /* = nullptr */)
{
  if (!_generating_initial_particles)
    mooseError("generateInitialParticle() can only be called during generateInitialParticles()");

  std::shared_ptr<Ray> ray = _initial_particles_replicated ? acquireReplicatedRay() : acquireRay();
  ray->setStart(point, elem);
  _initial_particles.emplace_back(std::move(ray));
  return _initial_particles.back();
}
