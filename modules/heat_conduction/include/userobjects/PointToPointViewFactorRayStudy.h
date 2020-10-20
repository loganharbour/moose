//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ViewFactorRayStudyBase.h"

/**
 * RayTracingStudy used to generate Rays for view factor computation
 * using the point-to-point method.
 */
class PointToPointViewFactorRayStudy : public ViewFactorRayStudyBase
{
public:
  PointToPointViewFactorRayStudy(const InputParameters & parameters);

  static InputParameters validParams();

  void initialSetup() override;

  /**
   * This is called in RayStudyTraceRay to grab the RayBCs on a boundary.
   *
   * With view factor computation, we only have one RayBC. Therefore, cache it up front in
   * initialSetup() and return the cached object here.
   */
  void getRayBCs(std::vector<RayBoundaryConditionBase *> & result,
                 const std::vector<TraceRayBndElement> &,
                 THREAD_ID tid,
                 RayID) override
  {
    result = _threaded_cached_ray_bcs[tid];
  }

  /**
   * Data structure used for storing all of the information needed to spawn Rays to a single element
   */
  struct EndElem
  {
    EndElem() {}
    EndElem(const dof_id_type elem_id,
            const unsigned short int side,
            const BoundaryID bnd_id,
            const std::vector<Point> & points,
            const std::vector<Real> & weights)
      : _elem_id(elem_id), _side(side), _bnd_id(bnd_id), _points(points), _weights(weights)
    {
    }

    /// The end element ID
    dof_id_type _elem_id;
    /// The side on the end element
    unsigned short int _side;
    /// The end boundary ID
    BoundaryID _bnd_id;
    /// The points
    std::vector<Point> _points;
    /// The weights
    std::vector<Real> _weights;
  };

  /**
   * Gets the index in the Ray aux data for the ending boundary ID
   */
  RayDataIndex rayIndexEndBndID() const { return _ray_index_end_bnd_id; }
  /**
   * Gets the index in the Ray aux data for the ending weight
   */
  RayDataIndex rayIndexEndWeight() const { return _ray_index_end_weight; }
  /**
   * Gets the index in the Ray aux data for the ending elem ID
   */
  RayDataIndex rayIndexEndElemID() const { return _ray_index_end_elem_id; }
  /**
   * Gets the index in the Ray aux data for the ending side
   */
  RayDataIndex rayIndexEndSide() const { return _ray_index_end_side; }

protected:
  void preExecuteStudy() override;
  void generateRays() override;

  /// The ray direction * side normal tolerance for spawning rays
  const Real _dot_tol;

  /// Index in the Ray aux data for the ending boundary ID
  const RayDataIndex _ray_index_end_bnd_id;
  /// Index in the Ray aux data for the ending weight
  const RayDataIndex _ray_index_end_weight;
  /// Index in the Ray aux data for the ending elem ID
  const RayDataIndex _ray_index_end_elem_id;
  /// Index in the Ray aux data for the ending side
  const RayDataIndex _ray_index_end_side;

private:
  void postAddStartElem(const Elem * elem,
                        const unsigned short side,
                        const BoundaryID bnd_id,
                        const std::vector<Point> & points,
                        const std::vector<Real> & weights) override;

  /**
   * Fills _end_elems
   */
  void generateEndElems();

  /**
   * Helper for whether or not this processor should spawn a Ray with the given conditions
   *
   * This helps distribute work across processors
   */
  bool shouldGenerate(const StartElem & start_elem, const EndElem & end_elem) const;

  /// Used for caching the single RayBC per thread for use in getRayBCs()
  std::vector<std::vector<RayBoundaryConditionBase *>> _threaded_cached_ray_bcs;

  /// Stores the EndElem information for the elements we need to trace to
  std::vector<EndElem> _end_elems;
};

namespace libMesh
{
namespace Parallel
{

template <>
class Packing<PointToPointViewFactorRayStudy::EndElem>
{
public:
  typedef Real buffer_type;

  static unsigned int packing_size(const std::size_t num_points);

  static unsigned int packed_size(typename std::vector<Real>::const_iterator in);

  static unsigned int packable_size(const PointToPointViewFactorRayStudy::EndElem & end_elem,
                                    const void *);

  template <typename Iter>
  static void
  pack(const PointToPointViewFactorRayStudy::EndElem & object, Iter data_out, const void *);

  template <typename BufferIter>
  static PointToPointViewFactorRayStudy::EndElem unpack(BufferIter in, void *);
};

} // namespace Parallel

} // namespace libMesh
