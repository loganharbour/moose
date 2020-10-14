//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AQViewFactorRayBC.h"

// Local includes
#include "AQViewFactorRayStudy.h"

registerMooseObject("RayTracingApp", AQViewFactorRayBC);

InputParameters
AQViewFactorRayBC::validParams()
{
  InputParameters params = RayBC::validParams();
  return params;
}

AQViewFactorRayBC::AQViewFactorRayBC(const InputParameters & params)
  : RayBC(params),
    _vf_study(RayTracingStudy::castFromStudy<AQViewFactorRayStudy>(params)),
    _ray_index_start_bnd_id(_vf_study.rayIndexStartBndID()),
    _ray_index_start_total_weight(_vf_study.rayIndexStartTotalWeight())
{
}

void
AQViewFactorRayBC::apply(const Elem * /* elem */,
                         const unsigned short /* intersected_side */,
                         const BoundaryID bnd_id,
                         const Point & /* intersection_point */,
                         const std::shared_ptr<Ray> & ray,
                         const unsigned int num_applying)
{
  // The boundary ID this Ray started on
  const BoundaryID start_bnd_id = ray->auxData(_ray_index_start_bnd_id);
  // Starting total weight
  const Real start_total_weight = ray->auxData(_ray_index_start_total_weight);
  // Value to append (divide by num_applying if we hit an edge or node)
  const Real value = start_total_weight / (Real)num_applying;
  mooseAssert(!std::isnan(value), "Encountered NaN");

  // Accumulate into the view factor info
  _vf_study.addToViewFactorInfo(value, start_bnd_id, bnd_id, _tid);

  // Either hit an obstacle here or hit its end and contributed: done with this Ray
  ray->setShouldContinue(false);
}
