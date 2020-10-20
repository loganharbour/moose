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

// libMesh includes
#include "libmesh/parallel.h"

/**
 * Base RayTracingStudy that view factor studies inherit from
 */
class ViewFactorRayStudyBase : public RayTracingStudy
{
public:
  ViewFactorRayStudyBase(const InputParameters & parameters);

  static InputParameters validParams();

  void initialSetup() override;

  /**
   * Data structure used for storing all of the information needed to spawn
   * Rays from a single element.
   */
  struct StartElem
  {
    StartElem() {}

    StartElem(const Elem * elem,
              const Elem * start_elem,
              const unsigned short int incoming_side,
              const BoundaryID bnd_id,
              const std::vector<Point> & points,
              const std::vector<Real> & weights)
      : _elem(elem),
        _start_elem(start_elem),
        _incoming_side(incoming_side),
        _bnd_id(bnd_id),
        _points(points),
        _weights(weights)
    {
      mooseAssert(_points.size() == _weights.size(), "Point and weight size not equal");
    }

    /// The element the points originate from
    const Elem * _elem;
    /// The element the trace will start from
    const Elem * _start_elem;
    /// The incoming side on start_elem that the trace will start from
    unsigned short int _incoming_side;
    /// The boundary ID associated with this start elem
    BoundaryID _bnd_id;
    /// The points on start_elem to spawn Rays from
    std::vector<Point> _points;
    /// The weights associated with each point
    std::vector<Real> _weights;
  };

  /**
   * Adds into the view factor info; to be used in ViewFactorRayBC
   * @param value The value to add
   * @param from_id The from boundary
   * @param to_id The to boundary
   * @param tid The thread
   */
  void addToViewFactorInfo(Real value,
                           const BoundaryID from_id,
                           const BoundaryID to_id,
                           const THREAD_ID tid);

  /**
   * Accessor for the finalized view factor info
   * @param from_id The from boundary
   * @param to_id The to boundary
   */
  Real viewFactorInfo(const BoundaryID from_id, const BoundaryID to_id) const;

  /**
   * Get the index in the Ray aux data for the starting boundary ID
   */
  RayDataIndex rayIndexStartBndID() const { return _ray_index_start_bnd_id; }
  /**
   * Get the index in the Ray aux data for the starting total weight (dot * qp weight)
   */
  RayDataIndex rayIndexStartTotalWeight() const { return _ray_index_start_total_weight; }

protected:
  void generateRays() override = 0;
  void preExecuteStudy() override;
  void postExecuteStudy() override;

  virtual void postAddStartElem(const Elem * elem,
                                const unsigned short side,
                                const BoundaryID bnd_id,
                                const std::vector<Point> & points,
                                const std::vector<Real> & weights);

  /// The vector of user supplied boundary IDs we need view factors on
  const std::vector<BoundaryID> _bnd_ids_vec;
  /// The user supplied boundary IDs we need view factors on
  const std::set<BoundaryID> _bnd_ids;

  /// The convention for spawning rays from internal sidesets
  const MooseEnum _internal_convention;

  /// Index in the Ray aux data for the starting boundary ID
  const RayDataIndex _ray_index_start_bnd_id;
  /// Index in the Ray aux data for the starting total weight (dot * qp weight)
  const RayDataIndex _ray_index_start_total_weight;

  /// Face FE used for creating face quadrature points and weights
  const std::unique_ptr<FEBase> _fe_face;
  /// Face quadrature used for _fe_face
  const std::unique_ptr<QBase> _q_face;

  /// The StartElem objects that this proc needs to spawn Rays from
  std::vector<StartElem> _start_elems;

  /// Whether or not this view factor method exhibits reciprocity
  const bool _has_reciprocity;

private:
  void generateStartElems();

  /// View factor information by tid and then from/to pair
  std::vector<std::unordered_map<BoundaryID, std::unordered_map<BoundaryID, Real>>>
      _threaded_vf_info;
  /// Cumulative view factor information
  std::map<BoundaryID, std::map<BoundaryID, Real>> _vf_info;
};

namespace libMesh
{
namespace Parallel
{
template <>
class Packing<ViewFactorRayStudyBase::StartElem>
{
public:
  typedef Real buffer_type;

  static unsigned int packing_size(const std::size_t num_points);

  static unsigned int packed_size(typename std::vector<Real>::const_iterator in);

  static unsigned int packable_size(const ViewFactorRayStudyBase::StartElem & start_elem,
                                    const void *);

  template <typename Iter, typename Context>
  static void
  pack(const ViewFactorRayStudyBase::StartElem & object, Iter data_out, const Context *);

  template <typename BufferIter, typename Context>
  static ViewFactorRayStudyBase::StartElem unpack(BufferIter in, Context *);
};

} // namespace Parallel

} // namespace libMesh
