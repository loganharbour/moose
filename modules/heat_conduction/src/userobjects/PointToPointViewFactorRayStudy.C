//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PointToPointViewFactorRayStudy.h"

// MOOSE includes
#include "TimedPrint.h"

// Local includes
#include "PointToPointViewFactorRayBC.h"

// libMesh includes
#include "libmesh/parallel_sync.h"

registerMooseObject("HeatConductionApp", PointToPointViewFactorRayStudy);

InputParameters
PointToPointViewFactorRayStudy::validParams()
{
  auto params = ViewFactorRayStudyBase::validParams();

  params.addRangeCheckedParam<Real>(
      "dot_tol",
      1.0e-3,
      "dot_tol > 0",
      "Tolerance for the dot product of a ray and a side normal, i.e., if the dot "
      "product is less than this, the ray will not be spawned");

  params.set<bool>("_has_reciprocity") = true;

  return params;
}

PointToPointViewFactorRayStudy::PointToPointViewFactorRayStudy(const InputParameters & parameters)
  : ViewFactorRayStudyBase(parameters),
    _dot_tol(getParam<Real>("dot_tol")),
    _ray_index_end_bnd_id(registerRayAuxData("end_bnd_id")),
    _ray_index_end_weight(registerRayAuxData("end_weight")),
    _ray_index_end_elem_id(registerRayAuxData("end_elem_id")),
    _ray_index_end_side(registerRayAuxData("end_side")),
    _threaded_cached_ray_bcs(libMesh::n_threads())
{
}

void
PointToPointViewFactorRayStudy::initialSetup()
{
  RayTracingStudy::initialSetup();

  std::vector<RayBC *> ray_bcs;

  // Make sure we have only one RayBC: a PointToPointViewFactorRayBC with the same boundaries
  RayTracingStudy::getRayBCs(ray_bcs, 0);
  if (ray_bcs.size() != 1 || !dynamic_cast<PointToPointViewFactorRayBC *>(ray_bcs[0]))
    mooseError(_error_prefix, ": Requires one and only one RayBC, a PointToPointViewFactorRayBC");
  if (!ray_bcs[0]->hasBoundary(_bnd_ids))
    mooseError("The PointToPointViewFactorRayBC '",
               ray_bcs[0]->name(),
               "' must be applied to the same boundaries as the ",
               _error_prefix);

  // Cache our one RayBC per thread so that we don't spend unnecessary time querying for RayBCs on
  // every boundary
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    RayTracingStudy::getRayBCs(ray_bcs, tid);
    _threaded_cached_ray_bcs[tid] = ray_bcs;
  }
}

void
PointToPointViewFactorRayStudy::preExecuteStudy()
{
  // Clear this before we start filling into it
  _end_elems.clear();

  ViewFactorRayStudyBase::preExecuteStudy();

  // Now that we have _start_elems filled, fill _end_elems
  generateEndElems();
}

void
PointToPointViewFactorRayStudy::postAddStartElem(const Elem * elem,
                                                 const unsigned short side,
                                                 const BoundaryID bnd_id,
                                                 const std::vector<Point> & points,
                                                 const std::vector<Real> & weights)
{
  // We also need to add starting information to an EndElem object so we can communicate
  // who everyone should spawn Rays to
  _end_elems.emplace_back(elem->id(), side, bnd_id, points, weights);
}

void
PointToPointViewFactorRayStudy::generateEndElems()
{
  // Decide on the processors that have StartingElems so we know who to send to
  std::vector<dof_id_type> have_points(_communicator.size());
  _communicator.allgather((dof_id_type)_start_elems.size(), have_points);

  // Fill the send data to each proc that needs it
  std::unordered_map<processor_id_type, std::vector<EndElem>> send_end_map;
  for (processor_id_type pid = 0; pid < _communicator.size(); ++pid)
    if (pid != _pid && have_points[pid])
      send_end_map.emplace(pid, _end_elems);

  // Communciate and fill into our _end_elems
  auto end_functor = [this](processor_id_type, const std::vector<EndElem> & end_elems) {
    _end_elems.reserve(_end_elems.size() + end_elems.size());
    for (const EndElem & end_elem : end_elems)
      _end_elems.emplace_back(end_elem);
  };
  Parallel::push_parallel_packed_range(_communicator, send_end_map, (void *)nullptr, end_functor);
}

bool
PointToPointViewFactorRayStudy::shouldGenerate(const StartElem & start_elem,
                                               const EndElem & end_elem) const
{
  const auto start_elem_id = start_elem._elem->id();

  if (start_elem_id == end_elem._elem_id)
    return start_elem._bnd_id < end_elem._bnd_id;

  const bool start_even = start_elem_id % 2 == 0;
  const bool end_even = end_elem._elem_id % 2 == 0;

  if (start_even && end_even)
    return start_elem_id > end_elem._elem_id;
  if (!start_even && !end_even)
    return end_elem._elem_id > start_elem_id;
  else
    return start_even;
}

void
PointToPointViewFactorRayStudy::generateRays()
{
  CONSOLE_TIMED_PRINT("PointToPointViewFactorRayStudy generating rays");

  // Determine number of Rays and points to allocate space before generation and for output
  std::size_t num_local_rays = 0;
  std::size_t num_local_start_points = 0;
  for (const auto & start_elem : _start_elems)
  {
    num_local_start_points += start_elem._points.size();

    for (const auto & end_elem : _end_elems)
      if (shouldGenerate(start_elem, end_elem))
        for (const auto & start_point : start_elem._points)
        {
          const auto & normal =
              getSideNormal(start_elem._start_elem, start_elem._incoming_side, /* tid = */ 0);

          for (const auto & end_point : end_elem._points)
            if (!start_point.absolute_fuzzy_equals(end_point))
            {
              const Point direction = (end_point - start_point).unit();
              if (normal * direction < -_dot_tol)
                ++num_local_rays;
            }
        }
  }

  // Print out totals while we're here
  std::size_t num_total_points = num_local_start_points;
  std::size_t num_total_rays = num_local_rays;
  _communicator.sum(num_total_points);
  _communicator.sum(num_total_rays);
  _console << "PointToPointViewFactorRayStudy generated " << num_total_points
           << " points requiring " << num_total_rays << " rays" << std::endl;

  // Reserve space in the buffer ahead of time before we fill it
  reserveRayBuffer(num_local_rays);

  // Loop over...
  for (const auto & start_elem : _start_elems)  // all start elements
    for (const auto & end_elem : _end_elems)    // all end elements
      if (shouldGenerate(start_elem, end_elem)) // we're the one to generate this combo
        for (std::size_t start_i = 0; start_i < start_elem._points.size(); ++start_i) // all points
        {
          const auto & normal =
              getSideNormal(start_elem._start_elem, start_elem._incoming_side, /* tid = */ 0);

          for (std::size_t end_i = 0; end_i < end_elem._points.size(); ++end_i)
          {
            // Don't spawn to myself
            if (start_elem._points[start_i].absolute_fuzzy_equals(end_elem._points[end_i]))
              continue;

            const auto direction = (end_elem._points[end_i] - start_elem._points[start_i]).unit();
            const Real dot = normal * direction;
            if (dot > -_dot_tol)
              continue;

            // Acquire a Ray and fill with the starting information
            std::shared_ptr<Ray> ray = acquireRay(/* tid = */ 0, rayDataSize(), rayAuxDataSize());
            ray->setStartingElem(start_elem._start_elem);
            ray->setIncomingSide(start_elem._incoming_side);
            ray->setStartDirection(start_elem._points[start_i], direction);
            ray->setID(generateUniqueRayID());
            ray->setAuxData(_ray_index_start_bnd_id, start_elem._bnd_id);
            ray->setAuxData(_ray_index_start_total_weight,
                            start_elem._weights[start_i] * std::abs(dot));
            ray->setAuxData(_ray_index_end_bnd_id, end_elem._bnd_id);
            ray->setAuxData(_ray_index_end_weight, end_elem._weights[end_i]);
            ray->setAuxData(_ray_index_end_elem_id, end_elem._elem_id);
            ray->setAuxData(_ray_index_end_side, end_elem._side);

            // Move the Ray to the buffer to be traced
            moveRayToBuffer(ray);
          }
        }
}

namespace libMesh
{
namespace Parallel
{

unsigned int
Packing<PointToPointViewFactorRayStudy::EndElem>::packing_size(const std::size_t num_points)
{
  // Number of points, elem_id, side, bnd_id
  unsigned int total_size = 4;
  // Points
  total_size += num_points * 3;
  // Weights
  total_size += num_points;

  return total_size;
}

unsigned int
Packing<PointToPointViewFactorRayStudy::EndElem>::packed_size(
    typename std::vector<buffer_type>::const_iterator in)
{
  unsigned int num_points = *in++;
  return packing_size(num_points);
}

unsigned int
Packing<PointToPointViewFactorRayStudy::EndElem>::packable_size(
    const PointToPointViewFactorRayStudy::EndElem & end_elem, const void *)
{
  mooseAssert(end_elem._points.size() == end_elem._weights.size(), "Size mismatch");
  return packing_size(end_elem._points.size());
}

template <>
PointToPointViewFactorRayStudy::EndElem
Packing<PointToPointViewFactorRayStudy::EndElem>::unpack(std::vector<Real>::const_iterator in,
                                                         void *)
{
  // EndElem to fill into
  PointToPointViewFactorRayStudy::EndElem end_elem;

  // Number of points
  const std::size_t num_points = *in++;

  // Elem id
  // Because the buffer type is a Real, we make sure to preserve this value by copying bytes
  std::memcpy(&end_elem._elem_id, &(*in++), sizeof(dof_id_type));

  // Side
  end_elem._side = static_cast<unsigned short>(*in++);

  // Boundary ID
  end_elem._bnd_id = static_cast<BoundaryID>(*in++);

  // Points
  end_elem._points.resize(num_points);
  for (std::size_t i = 0; i < num_points; ++i)
  {
    end_elem._points[i](0) = *in++;
    end_elem._points[i](1) = *in++;
    end_elem._points[i](2) = *in++;
  }

  // Weights
  end_elem._weights.resize(num_points);
  for (std::size_t i = 0; i < num_points; ++i)
    end_elem._weights[i] = *in++;

  return end_elem;
}

template <>
void
Packing<PointToPointViewFactorRayStudy::EndElem>::pack(
    const PointToPointViewFactorRayStudy::EndElem & end_elem,
    std::back_insert_iterator<std::vector<buffer_type>> data_out,
    const void *)
{
  // Number of points
  data_out = static_cast<buffer_type>(end_elem._points.size());

  // Elem id
  // Because the buffer type is a Real, we make sure to preserve this value by copying bytes
  static_assert(sizeof(dof_id_type) <= sizeof(buffer_type), "Elem id won't fit into buffer type");
  Real id_as_Real;
  std::memcpy(&id_as_Real, &end_elem._elem_id, sizeof(dof_id_type));
  data_out = id_as_Real;

  // Side
  data_out = static_cast<buffer_type>(end_elem._side);

  // Boundary id
  data_out = static_cast<buffer_type>(end_elem._bnd_id);

  // Points
  for (const auto & point : end_elem._points)
  {
    data_out = point(0);
    data_out = point(1);
    data_out = point(2);
  }

  // Weights
  std::copy(end_elem._weights.begin(), end_elem._weights.end(), data_out);
}

} // namespace Parallel

} // namespace libMesh
