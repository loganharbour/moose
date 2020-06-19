//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ViewFactorRayBC.h"

// Local includes
#include "ViewFactorRayStudy.h"

registerMooseObject("RayTracingApp", ViewFactorRayBC);

InputParameters
ViewFactorRayBC::validParams()
{
  InputParameters params = RayBC::validParams();
  return params;
}

ViewFactorRayBC::ViewFactorRayBC(const InputParameters & params)
  : RayBC(params),
    _vf_study(ViewFactorRayStudy::castFromStudy(params)),
    _ray_index_start_bnd_id(_study.getRayAuxDataIndex("start_bnd_id")),
    _ray_index_start_total_weight(_study.getRayAuxDataIndex("start_total_weight")),
    _ray_index_end_bnd_id(_study.getRayAuxDataIndex("end_bnd_id")),
    _ray_index_end_weight(_study.getRayAuxDataIndex("end_weight")),
    _ray_index_start_end_distance(_study.getRayAuxDataIndex("start_end_distance"))
{
}

void
ViewFactorRayBC::apply(const Elem * elem,
                       const unsigned short intersected_side,
                       const BoundaryID bnd_id,
                       const Point & /* intersection_point */,
                       const std::shared_ptr<Ray> & ray,
                       const bool applying_multiple)
{
  // If we hit the end boundary and are at the end point (determined by the correct distance), then
  // we can contribute to the view factor info
  if (ray->auxData(_ray_index_end_bnd_id) == bnd_id &&
      MooseUtils::absoluteFuzzyEqual(ray->distance(), ray->auxData(_ray_index_start_end_distance)))
  {
    if (applying_multiple)
      mooseError("Should not contribute while applying multiple ViewFactorRayBC");

    // The boundary ID this Ray started on
    const BoundaryID start_bnd_id = ray->auxData(_ray_index_start_bnd_id);
    // Starting total weight
    const Real start_total_weight = ray->auxData(_ray_index_start_total_weight);

    // Ending total weight
    // NOTE: This should really be _normals[_qp] but because ray-tracing only works
    // with first order elements, we simply use the centroid normal from the study
    const auto & normal = _study.getSideNormal(elem, intersected_side, _tid);
    const Real end_total_weight =
        std::abs(ray->direction() * normal) * ray->auxData(_ray_index_end_weight);

    // Accumulate into the view factor info
    Real denom = ray->distance();
    if (_mesh.dimension() == 3)
      denom *= ray->distance();
    const Real value = start_total_weight * end_total_weight / denom;
    _vf_study.viewFactorInfo(start_bnd_id, bnd_id, _tid) += value;
  }

  // Either hit an obstacle here or hit its end and contributed: done with this Ray
  ray->setShouldContinue(false);
}
