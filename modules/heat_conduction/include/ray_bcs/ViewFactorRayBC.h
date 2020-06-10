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
  const RayDataIndex _ray_index_start_dot;
  /// Index in the Ray aux data for the starting boundary ID
  const RayDataIndex _ray_index_start_bnd_id;
  /// Index in the Ray aux data for the starting weight
  const RayDataIndex _ray_index_start_weight;
  /// Index in the Ray aux data for the ending boundary ID
  const RayDataIndex _ray_index_end_bnd_id;
  /// Index in the Ray aux data for the ending weight
  const RayDataIndex _ray_index_end_weight;
  /// The x-index for the point on the boundary where this Ray should end
  const RayDataIndex _ray_index_end_x;
  /// The y-index for the point on the boundary where this Ray should end
  const RayDataIndex _ray_index_end_y;
  /// The z-index for the point on the boundary where this Ray should end
  const RayDataIndex _ray_index_end_z;

  ViewFactorRayStudy * _vf_study;
};
