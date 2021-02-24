//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MPMBoundingBoxRayStudy.h"

// MOOSE includes
#include "Assembly.h"

registerMooseObject("RayTracingApp", MPMBoundingBoxRayStudy);

InputParameters
MPMBoundingBoxRayStudy::validParams()
{
  auto params = MPMRayStudyBase::validParams();

  params.addRequiredParam<Point>("bbox_min", "Minimum bounding box for the particle grid");
  params.addRequiredParam<Point>("bbox_max", "Maximum bounding box for the particle grid");
  params.addRequiredParam<std::vector<unsigned int>>("grid_intervals",
                                                     "Intervals for the particle grid");

  return params;
}

MPMBoundingBoxRayStudy::MPMBoundingBoxRayStudy(const InputParameters & parameters)
  : MPMRayStudyBase(parameters),
    _bbox_min(getParam<Point>("bbox_min")),
    _bbox_max(getParam<Point>("bbox_max")),
    _grid_intervals(getParam<std::vector<unsigned int>>("grid_intervals"))
{
  if (_grid_intervals.size() != _mesh.dimension())
    paramError("grid_intervals", "Must provide as many intervals as the mesh dimension");
  if (!boundingBox().contains_point(_bbox_min))
    paramError("bbox_min", "Not contained within the mesh");
  if (!boundingBox().contains_point(_bbox_max))
    paramError("bbox_max", "Not contained within the mesh");
}

void
MPMBoundingBoxRayStudy::generateInitialParticles()
{
  // Width of the grid in each dimension
  const auto grid_width = _bbox_max - _bbox_min;
  // Delta for each grid interval
  Point dgrid;
  for (unsigned int d = 0; d < _mesh.dimension(); ++d)
    dgrid(d) = grid_width(d) / (Real)_grid_intervals[d];

  // Create a particle (ray) for each point in the grid
  Point particle_point;
  for (unsigned int i = 0; i < (_grid_intervals[0] + 1); ++i)
  {
    particle_point(0) = _bbox_min(0) + dgrid(0) * (Real)i;
    for (unsigned int j = 0; j < (_mesh.dimension() > 1 ? _grid_intervals[1] + 1 : 1); ++j)
    {
      particle_point(1) = _bbox_min(1) + dgrid(1) * (Real)j;
      for (unsigned int k = 0; k < (_mesh.dimension() > 2 ? _grid_intervals[2] + 1 : 1); ++k)
      {
        particle_point(2) = _bbox_min(2) + dgrid(2) * (Real)k;

        auto & ray = generateInitialParticle(particle_point);
        ray->auxData(_particle_mass_index) = 1;
      }
    }
  }
}
