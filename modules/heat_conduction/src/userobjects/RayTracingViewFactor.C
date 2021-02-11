//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "RayTracingViewFactor.h"

#include "ViewFactorRayStudy.h"

registerMooseObject("HeatConductionApp", RayTracingViewFactor);

InputParameters
RayTracingViewFactor::validParams()
{
  InputParameters params = ViewFactorBase::validParams();
  params.addRequiredParam<UserObjectName>("ray_study_name",
                                          "Name of the view factor ray study UO.");
  params.addClassDescription("Computes view factors for arbitrary geometries using raytracing.");
  return params;
}

RayTracingViewFactor::RayTracingViewFactor(const InputParameters & parameters)
  : ViewFactorBase(parameters), _ray_study(getUserObject<ViewFactorRayStudy>("ray_study_name"))
{
  if (_mesh.dimension() == 1)
    mooseError("View factor calculations for 1D geometry makes no sense");
  else if (_mesh.dimension() == 2)
    _divisor = 2;
  else
    _divisor = libMesh::pi;
}

void
RayTracingViewFactor::execute()
{
  // compute areas
  auto current_boundary_name = _mesh.getBoundaryName(_current_boundary_id);
  auto it = _side_name_index.find(current_boundary_name);
  if (it == _side_name_index.end())
    mooseError("Current boundary name: ",
               current_boundary_name,
               " with id ",
               _current_boundary_id,
               " not in boundary parameter.");

  _areas[it->second] += _current_side_volume;
}

void
RayTracingViewFactor::initialize()
{
  // set view_factors to zero
  for (unsigned int j = 0; j < _n_sides; ++j)
    _areas[j] = 0;
}

void
RayTracingViewFactor::finalizeViewFactor()
{
  gatherSum(_areas);

  // get the _view_factors from ray study
  for (auto & from_name : boundaryNames())
    for (auto & to_name : boundaryNames())
    {
      unsigned int from_index = getSideNameIndex(from_name);
      unsigned int to_index = getSideNameIndex(to_name);
      BoundaryID from_id = _mesh.getBoundaryID(from_name);
      BoundaryID to_id = _mesh.getBoundaryID(to_name);
      _view_factors[from_index][to_index] = _ray_study.viewFactorInfo(from_id, to_id);
    }

  // divide view_factor Fij by Ai and pi
  for (unsigned int i = 0; i < _n_sides; ++i)
    for (auto & vf : _view_factors[i])
      vf /= (_areas[i] * _divisor);
}

void
RayTracingViewFactor::threadJoinViewFactor(const UserObject & y)
{
  const RayTracingViewFactor & pps = static_cast<const RayTracingViewFactor &>(y);
  for (unsigned int i = 0; i < _n_sides; ++i)
    _areas[i] += pps._areas[i];
}
