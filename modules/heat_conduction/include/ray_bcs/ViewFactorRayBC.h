//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "RayBC.h"

class ViewFactorRayStudy;

class ViewFactorRayBC : public RayBC
{
public:
  ViewFactorRayBC(const InputParameters & params);

  static InputParameters validParams();

  virtual void apply(const Elem * elem,
                     const unsigned short intersected_side,
                     const BoundaryID bnd_id,
                     const Point & intersection_point,
                     const std::shared_ptr<Ray> & ray,
                     const bool applying_at_corner) override;

protected:
  /// Index in the Ray aux data for the starting dot product
  const unsigned int _ray_index_start_dot;
  /// Index in the Ray aux data for the starting boundary ID
  const unsigned int _ray_index_start_bnd_id;
  /// Index in the Ray aux data for the starting weight
  const unsigned int _ray_index_start_weight;
  /// Index in the Ray aux data for the ending weight
  const unsigned int _ray_index_end_weight;

  ViewFactorRayStudy * _vf_study;
};
