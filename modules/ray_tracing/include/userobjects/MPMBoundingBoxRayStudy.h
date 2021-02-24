//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MPMRayStudyBase.h"

class MPMBoundingBoxRayStudy : public MPMRayStudyBase
{
public:
  MPMBoundingBoxRayStudy(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  void generateInitialParticles() override;

  /// Minimum of the bounding box for spawning particles
  const Point _bbox_min;
  /// Maximum of the bounding box for spawning particles
  const Point _bbox_max;
  /// Intervals in the bounding box for spawning particles
  const std::vector<unsigned int> _grid_intervals;
};
