//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ViewFactorRayBC.h"
#include "ViewFactorRayStudy.h"
#include "RayTracingStudy.h"

registerMooseObject("RayTracingApp", ViewFactorRayBC);

InputParameters
ViewFactorRayBC::validParams()
{
  InputParameters params = RayBC::validParams();
  return params;
}

ViewFactorRayBC::ViewFactorRayBC(const InputParameters & params)
  : RayBC(params),
    _ray_index_start_dot(_study.getRayAuxDataIndex("start_dot")),
    _ray_index_start_bnd_id(_study.getRayAuxDataIndex("start_bnd_id")),
    _ray_index_start_weight(_study.getRayAuxDataIndex("start_weight")),
    _ray_index_end_bnd_id(_study.getRayAuxDataIndex("end_bnd_id")),
    _ray_index_end_weight(_study.getRayAuxDataIndex("end_weight")),
    _vf_study(nullptr)
{
  // make sure we have a have a ViewFactorRayStudy
  _vf_study = dynamic_cast<ViewFactorRayStudy *>(&_study);
  if (!_vf_study)
    mooseError("ViewFactorRayBC must be paired with ViewFactorRayStudy.");
}

void
ViewFactorRayBC::apply(const Elem * /* elem */,
                       const unsigned short /* intersected_side */,
                       const BoundaryID bnd_id,
                       const Point & /* intersection_point */,
                       const std::shared_ptr<Ray> & ray,
                       const bool /* applying_at_corner */)
{
  const BoundaryID end_bnd_id = ray->auxData(_ray_index_end_bnd_id);

  // If we actually hit the boundary we were trying to get to: contribute
  if (bnd_id == end_bnd_id)
  {
    const Real dot = ray->auxData(_ray_index_start_dot);
    // NOTE: in reality this should use _normals[_qp] BUT
    // raytracing works only with first order elems so all
    // normals are the same; the absolute value is used to not
    // having to worry of the sign of the normal
    const Real dot_end = std::abs(ray->direction() * _normals[0]);
    const BoundaryID start_bnd_id = ray->auxData(_ray_index_start_bnd_id);
    const Real start_weight = ray->auxData(_ray_index_start_weight);
    const Real end_weight = ray->auxData(_ray_index_end_weight);

    Real denom = ray->distance();
    if (_mesh.dimension() == 3)
      denom *= ray->distance();
    Real value = start_weight * end_weight * dot * dot_end / denom;
    _vf_study->viewFactorInfo(start_bnd_id, bnd_id, _tid) += value;
  }

  // Kill the Ray if it's moved some
  // We guard this with a distance check because some Rays may start on an internal boundary, in
  // which case they would be immediately killed before moving
  if (ray->distance() > 0)
    ray->setShouldContinue(false);
}
