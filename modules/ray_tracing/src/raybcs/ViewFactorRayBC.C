//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ViewFactorRayBC.h"

#include "RayTracingStudy.h"

registerMooseObject("RayTracingApp", ViewFactorRayBC);

InputParameters
ViewFactorRayBC::validParams()
{
  return RayBC::validParams();
}

ViewFactorRayBC::ViewFactorRayBC(const InputParameters & params)
  : RayBC(params),
    _ray_index_start_dot(_study.getRayAuxDataIndex("start_dot")),
    _ray_index_start_bnd_id(_study.getRayAuxDataIndex("start_bnd_id")),
    _ray_index_start_weight(_study.getRayAuxDataIndex("start_weight")),
    _ray_index_end_weight(_study.getRayAuxDataIndex("end_weight"))
{
}

void
ViewFactorRayBC::apply(const Elem * /* elem */,
                       const unsigned short /* intersected_side */,
                       const BoundaryID bnd_id,
                       const Point & intersection_point,
                       const std::shared_ptr<Ray> & ray)
{
  const Real dot = ray->auxData(_ray_index_start_dot);
  const BoundaryID start_bnd_id = ray->auxData(_ray_index_start_bnd_id);
  const Real start_weight = ray->auxData(_ray_index_start_weight);
  const Real end_weight = ray->auxData(_ray_index_end_weight);

  // If we're at the right end point, we made it
  if (intersection_point.absolute_fuzzy_equals(ray->end()))
  {
    std::cerr << "at end from bnd " << start_bnd_id << " to bid " << bnd_id;
    std::cerr << ": dot = " << dot;
    std::cerr << ", distance = " << ray->distance() << std::endl;
  }

  // Die regardless
  ray->setShouldContinue(false);
}
