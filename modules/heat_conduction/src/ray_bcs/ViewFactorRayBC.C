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
    _ray_index_end_x(_study.getRayAuxDataIndex("end_x")),
    _ray_index_end_y(_study.getRayAuxDataIndex("end_y")),
    _ray_index_end_z(_study.getRayAuxDataIndex("end_z")),
    _vf_study(nullptr)
{
  // make sure we have a have a ViewFactorRayStudy
  _vf_study = dynamic_cast<ViewFactorRayStudy *>(&_study);
  if (!_vf_study)
    mooseError("ViewFactorRayBC must be paired with ViewFactorRayStudy.");
}

void
ViewFactorRayBC::apply(const Elem * elem,
                       const unsigned short intersected_side,
                       const BoundaryID bnd_id,
                       const Point & intersection_point,
                       const std::shared_ptr<Ray> & ray,
                       const bool /* applying_at_corner */)
{
  // Nothing to do if we haven't moved
  if (ray->distance() == 0)
    return;

  // The boundary ID this Ray wants to end on
  const BoundaryID end_bnd_id = ray->auxData(_ray_index_end_bnd_id);

  // If we actually hit the boundary we were trying to get to: possibly contribute
  if (bnd_id == end_bnd_id)
  {
    // Build the end point from the individual values
    const Point end_point = Point(ray->auxData(_ray_index_end_x),
                                  ray->auxData(_ray_index_end_y),
                                  ray->auxData(_ray_index_end_z));

    // Accumulate only if we're at the end point
    if (end_point.absolute_fuzzy_equals(intersection_point))
    {
      // Starting dot product
      const Real start_dot = ray->auxData(_ray_index_start_dot);

      // Ending dot product
      // NOTE: This should really be _normals[_qp] but because ray-tracing only works
      // with first order elements, we simply use the cached centroid normal from the
      // study
      const auto & normal = _study.getSideNormal(elem, intersected_side, _tid);
      const Real dot_end = std::abs(ray->direction() * normal);

      // The boundary ID this Ray started on
      const BoundaryID start_bnd_id = ray->auxData(_ray_index_start_bnd_id);

      // The start and end weights for this Ray
      const Real start_weight = ray->auxData(_ray_index_start_weight);
      const Real end_weight = ray->auxData(_ray_index_end_weight);

      // Accumulate into the view factor info
      Real denom = ray->distance();
      if (_mesh.dimension() == 3)
        denom *= ray->distance();
      Real value = start_weight * end_weight * start_dot * dot_end / denom;
      _vf_study->viewFactorInfo(start_bnd_id, bnd_id, _tid) += value;
    }
  }

  // Kill the Ray
  ray->setShouldContinue(false);
}
