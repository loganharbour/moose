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
  virtual void initialize() override;
  virtual void finalize() override;

  struct StartElem
  {
    StartElem(const Elem * elem,
              const unsigned short side,
              const Point & normal,
              const std::vector<Point> & points,
              const std::vector<Real> & weights)
      : _elem(elem), _side(side), _normal(normal), _points(points), _weights(weights)
    {
    }

    const Elem * _elem;
    const unsigned short int _side;
    const Point _normal;
    const std::vector<Point> _points;
    const std::vector<Real> _weights;
  };

  struct ViewFactorEntry
  {
    ViewFactorEntry(BoundaryID from, BoundaryID to)
    : from_bnd_id(from), to_bnd_id(to), view_factor(0)
    {
    }

    ViewFactorEntry(BoundaryID from, BoundaryID to, Real v)
    : from_bnd_id(from), to_bnd_id(to), view_factor(v)
    {
    }

    BoundaryID from_bnd_id;
    BoundaryID to_bnd_id;
    Real view_factor;
  };

  // returns a writeable reference to _vf_info pair from_bnd_id -> to_bnd_id
  Real & viewFactorInfo(BoundaryID from_id, BoundaryID to_id, THREAD_ID tid);

  // const accessor into view factor info
  Real viewFactorInfo(BoundaryID from_id, BoundaryID to_id) const;

protected:
  virtual void generateRays() override;
  virtual void generatePoints();
  virtual void defineRays();

  void defineRay(const Elem * starting_elem,
                 const Point & start_point,
                 const Point & end_point,
                 const Point & normal,
                 const unsigned short side,
                 const BoundaryID bnd_id,
                 const Real start_weight,
                 const Real end_weight);

  /// The user supplied boundary IDs we need view factors on
  const std::vector<BoundaryID> _bnd_ids;

  /// Index in the Ray aux data for the starting dot product
  const unsigned int _ray_index_start_dot;
  /// Index in the Ray aux data for the starting boundary ID
  const unsigned int _ray_index_start_bnd_id;
  /// Index in the Ray aux data for the starting weight
  const unsigned int _ray_index_start_weight;
  /// Index in the Ray aux data for the ending weight
  const unsigned int _ray_index_end_weight;

  /// Face FE used for creating face normals
  std::unique_ptr<FEBase> _fe_face;
  /// Face quadrature used for creating face normals
  std::unique_ptr<QBase> _q_face;

  /// The objects that this proc needs to spawn Rays from (indexed by _bnd_ids)
  std::vector<std::vector<StartElem>> _start_info;
  /// The objects (point and weight) that this proc needs to spawn Rays to (indexed by _bnd_ids)
  std::vector<std::vector<std::pair<Point, Real>>> _end_points;

  /// The next available ID to assign to a Ray in defineRay()
  dof_id_type _next_id;

  /// view factor information by tid and then from/to pair
  std::vector<std::vector<ViewFactorEntry>> _vf_info;
};
