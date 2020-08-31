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
    _vf_study(AQViewFactorRayStudy::castFromStudy(params)),
    _ray_index_start_bnd_id(_study.getRayAuxDataIndex("start_bnd_id")),
    _ray_index_start_total_weight(_study.getRayAuxDataIndex("start_total_weight"))
{
}

void
AQViewFactorRayBC::apply(const Elem * elem,
                       const unsigned short intersected_side,
                       const BoundaryID bnd_id,
                       const Point & /* intersection_point */,
                       const std::shared_ptr<Ray> & ray,
                       const bool applying_multiple)
{
  // Hit the end boundary and are on the correct elem and side -> contribute to view factor info
  //if (applying_multiple)
  //  mooseError("Should not contribute while applying multiple AQViewFactorRayBC\n\n",
  //             ray->getInfo(&_study));

  // The boundary ID this Ray started on
  const BoundaryID start_bnd_id = ray->auxData(_ray_index_start_bnd_id);
  // Starting total weight
  const Real start_total_weight = ray->auxData(_ray_index_start_total_weight);

  // Accumulate into the view factor info
  _vf_study.addToViewFactorInfo(start_total_weight, start_bnd_id, bnd_id, _tid);

  if (std::isnan(start_total_weight))
    mooseError("Encountered NaN in AQViewFactorRayBC\n", ray->getInfo(&_study));

  // Either hit an obstacle here or hit its end and contributed: done with this Ray
  ray->setShouldContinue(false);
}
