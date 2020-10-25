//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AQViewFactorRayStudy.h"

// MOOSE includes
#include "TimedPrint.h"
#include "ReflectRayBC.h"

// Local includes
#include "AQViewFactorRayBC.h"
#include "HCAngularQuadrature.h"

registerMooseObject("HeatConductionApp", AQViewFactorRayStudy);

InputParameters
AQViewFactorRayStudy::validParams()
{
  auto params = ViewFactorRayStudyBase::validParams();

  params.addParam<unsigned int>(
      "polar_quad_order",
      16,
      "Order of the polar quadrature [polar angle is between ray and normal]. Must be even.");
  params.addParam<unsigned int>(
      "azimuthal_quad_order",
      8,
      "Order of the azimuthal quadrature per quadrant [azimuthal angle is measured in "
      "a plan perpendicular to the normal].");

  params.set<bool>("_has_reciprocity") = false;

  return params;
}

AQViewFactorRayStudy::AQViewFactorRayStudy(const InputParameters & parameters)
  : ViewFactorRayStudyBase(parameters)
{
  // create angular quadrature
  if (_mesh.dimension() == 2)
    HCAngularQuadrature::getHalfRangeAQ2D(
        getParam<unsigned int>("polar_quad_order"), _aq_angles, _aq_weights);
  else if (_mesh.dimension() == 3)
    HCAngularQuadrature::getHalfRangeAQ3D(4 * getParam<unsigned int>("azimuthal_quad_order"),
                                          getParam<unsigned int>("polar_quad_order"),
                                          _aq_angles,
                                          _aq_weights);
}

void
AQViewFactorRayStudy::initialSetup()
{
  ViewFactorRayStudyBase::initialSetup();

  // RayBC coverage checks (at least one AQViewFactorRayBC and optionally a ReflectRayBC)
  std::vector<RayBoundaryConditionBase *> ray_bcs;
  RayTracingStudy::getRayBCs(ray_bcs, 0);
  unsigned int vf_bc_count = 0;
  for (RayBoundaryConditionBase * rbc : ray_bcs)
  {
    auto view_factor_bc = dynamic_cast<AQViewFactorRayBC *>(rbc);
    if (view_factor_bc)
    {
      ++vf_bc_count;

      if (!view_factor_bc->hasBoundary(_bnd_ids))
        mooseError(_error_prefix,
                   ": The boundary restriction of ",
                   rbc->type(),
                   " '",
                   rbc->name(),
                   "' does not match 'boundary'");
    }
    else
    {
      auto reflect_bc = dynamic_cast<ReflectRayBC *>(rbc);
      if (reflect_bc)
      {
        if (reflect_bc->hasBoundary(_bnd_ids))
          mooseError(_error_prefix,
                     ":\nThe boundaries applied in ReflectRayBC '",
                     rbc->name(),
                     "' cannot include any of the boundaries in ",
                     type());
      }
      else
        mooseError(
            _error_prefix,
            " does not support the ",
            rbc->type(),
            " ray boundary condition.\nSupported RayBCs: ReflectRayBC and AQViewFactorRayBC.");
    }
    if (vf_bc_count != 1)
      mooseError(_error_prefix, " requires one and only one AQViewFactorRayBC.");
  }
}

void
AQViewFactorRayStudy::generateRays()
{
  CONSOLE_TIMED_PRINT("AQViewFactorRayStudy generating rays");

  // Determine number of Rays and points to allocate space before generation and for output
  std::size_t num_local_rays = 0;
  std::size_t num_local_start_points = 0;
  for (const auto & start_elem : _start_elems)
  {
    num_local_start_points += start_elem._points.size();
    num_local_rays += start_elem._points.size() * _aq_angles.size();
  }

  // Print out totals while we're here
  std::size_t num_total_points = num_local_start_points;
  std::size_t num_total_rays = num_local_rays;
  _communicator.sum(num_total_points);
  _communicator.sum(num_total_rays);
  _console << "AQViewFactorRayStudy generated " << num_total_points
           << " points with an angular quadrature of " << _aq_angles.size()
           << " directions per point requiring " << num_total_rays << " rays" << std::endl;

  // Reserve space in the buffer ahead of time before we fill it
  reserveRayBuffer(num_local_rays);

  // loop through all starting points and spawn rays from each for each point and angle
  for (const auto & start_elem : _start_elems)
  {
    // Get normal for the element we're starting on
    auto inward_normal =
        getSideNormal(start_elem._start_elem, start_elem._incoming_side, /* tid = */ 0);
    // We actually want the normal of the original element (remember that we may swap starting
    // elements to the element on the other face per requirements of the ray tracer)
    if (start_elem._start_elem != start_elem._elem)
      inward_normal *= -1;
    // Lastly, if the boundary is external and the internal convention is positive, we must
    // switch the normal because our AQ uses the inward normal
    if (_internal_convention == 0 &&
        !start_elem._start_elem->neighbor_ptr(start_elem._incoming_side))
      inward_normal *= -1;

    // Create rotation matrix and rotate vector omega (the normal)
    const auto rotation_matrix =
        HCAngularQuadrature::aqRoationMatrix(inward_normal, _mesh.dimension());

    // Loop through all points and then all directions
    for (std::size_t start_i = 0; start_i < start_elem._points.size(); ++start_i)
      for (std::size_t l = 0; l < _aq_angles.size(); ++l)
      {
        const auto direction = HCAngularQuadrature::getAngularDirection(
            l, rotation_matrix, _aq_angles, _mesh.dimension());

        // angular weight function differs in 2D/3D
        // 2D: the quadrature abscissae are the angles between direction & normal.
        //     The weight is the cosine of that angle
        // 3D: the quadrature abscissae are the azimuthal angle phi and the cosine of the angle
        //     between normal and direction (== mu). The weight is mu in that case.
        const auto awf =
            _mesh.dimension() == 3 ? _aq_angles[l].second : std::cos(_aq_angles[l].second);
        const auto start_weight = start_elem._weights[start_i] * _aq_weights[l] * awf;

        // Acquire a Ray and fill with the starting information
        std::shared_ptr<Ray> ray = acquireRay(/* tid = */ 0, rayDataSize(), rayAuxDataSize());
        ray->setStartDirection(start_elem._points[start_i], direction);
        ray->setStartingElem(start_elem._start_elem);
        ray->setStartingIncomingSide(start_elem._incoming_side);
        ray->setID(generateUniqueRayID());
        ray->setAuxData(_ray_index_start_bnd_id, start_elem._bnd_id);
        ray->setAuxData(_ray_index_start_total_weight, start_weight);

        // Move the Ray into the buffer to be traced
        moveRayToBuffer(ray);
      }
  }
}
