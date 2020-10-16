//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PointToPointViewFactorRayBC.h"

// Local includes
#include "PointToPointViewFactorRayStudy.h"

registerMooseObject("RayTracingApp", PointToPointViewFactorRayBC);

InputParameters
PointToPointViewFactorRayBC::validParams()
{
  InputParameters params = RayBC::validParams();
  return params;
}

PointToPointViewFactorRayBC::PointToPointViewFactorRayBC(const InputParameters & params)
  : RayBC(params),
    _vf_study(RayTracingStudy::castFromStudy<PointToPointViewFactorRayStudy>(params)),
    _ray_index_start_bnd_id(_vf_study.rayIndexStartBndID()),
    _ray_index_start_total_weight(_vf_study.rayIndexStartTotalWeight()),
    _ray_index_end_bnd_id(_vf_study.rayIndexEndBndID()),
    _ray_index_end_weight(_vf_study.rayIndexEndWeight()),
    _ray_index_end_elem_id(_vf_study.rayIndexEndElemID()),
    _ray_index_end_side(_vf_study.rayIndexEndSide())
{
}

void
PointToPointViewFactorRayBC::apply(const Elem * elem,
                                   const unsigned short intersected_side,
                                   const BoundaryID bnd_id,
                                   const Point & /* intersection_point */,
                                   const std::shared_ptr<Ray> & ray,
                                   const unsigned int /*num_applying*/)
{
  // Hit the end boundary and are on the correct elem and side -> contribute to view factor info
  if (ray->auxData(_ray_index_end_bnd_id) == bnd_id &&
      elem->id() == (dof_id_type)ray->auxData(_ray_index_end_elem_id) &&
      intersected_side == ray->auxData(_ray_index_end_side))
  {
    mooseAssert(num_applying == 1, "Should not contribute to multiple");

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
    mooseAssert(!std::isnan(value), "Encountered NaN");
    _vf_study.addToViewFactorInfo(value, start_bnd_id, bnd_id, _tid);
  }

  // Either hit an obstacle here or hit its end and contributed: done with this Ray
  ray->setShouldContinue(false);
}
