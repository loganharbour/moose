//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "RayTracingStudy.h"

class ViewFactorRayStudy : public RayTracingStudy
{
public:
  ViewFactorRayStudy(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual void generateRays() override;
  virtual void generatePoints();
  virtual void defineRays();

  void defineRay(const Elem * starting_elem,
                 const Point & start_point,
                 const Point & end_point,
                 const Point & normal,
                 const unsigned short side,
                 const BoundaryID bnd_id);

  Point sideNormal(const Elem * elem, const unsigned short side);

  /// The user supplied boundary IDs we need view factors on
  const std::vector<BoundaryID> _bnd_ids;

  /// Index in the Ray aux data for the starting dot product
  const unsigned int _ray_index_dot;
  /// Index in the Ray aux data for the starting boundary ID
  const unsigned int _ray_index_bnd_id;

  /// Face FE used for creating face normals
  std::unique_ptr<FEBase> _fe_face;
  /// Face quadrature used for creating face normals
  std::unique_ptr<QBase> _q_face;

  /// The objects that this proc needs to spawn Rays from (indexed by _bnd_ids)
  std::vector<std::vector<std::tuple<Point, const Elem *, unsigned short>>> _start_points;
  /// The objects that this proc needs to spawn Rays to (indexed by _bnd_ids)
  std::vector<std::vector<Point>> _end_points;

  /// The next available ID to assign to a Ray in defineRay()
  dof_id_type _next_id;
};
