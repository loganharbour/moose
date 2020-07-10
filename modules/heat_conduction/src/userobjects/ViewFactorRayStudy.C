//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ViewFactorRayStudy.h"

// MOOSE includes
#include "Conversion.h"
#include "TimedPrint.h"

// Local includes
#include "ViewFactorRayBC.h"

// libMesh includes
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel_sync.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature.h"

registerMooseObject("HeatConductionApp", ViewFactorRayStudy);

InputParameters
ViewFactorRayStudy::validParams()
{
  auto params = RayTracingStudy::validParams();

  params.addRequiredParam<std::vector<BoundaryName>>(
      "boundary", "The list of boundary IDs from the mesh where this boundary condition applies");

  MooseEnum qorders("CONSTANT FIRST SECOND THIRD FOURTH FIFTH SIXTH SEVENTH EIGHTH NINTH TENTH "
                    "ELEVENTH TWELFTH THIRTEENTH FOURTEENTH FIFTEENTH SIXTEENTH SEVENTEENTH "
                    "EIGHTTEENTH NINTEENTH TWENTIETH",
                    "CONSTANT");
  params.addParam<MooseEnum>("face_order", qorders, "The face quadrature rule order");

  MooseEnum qtypes("GAUSS GRID", "GRID");
  params.addParam<MooseEnum>("face_type", qtypes, "The face quadrature type");

  MooseEnum convention("positive=0 negative=1", "positive");
  params.addParam<MooseEnum>(
      "internal_convention",
      convention,
      "The convention for spawning rays from internal sidesets; denotes the sign of the dot "
      "product between a ray and the internal sideset side normal");

  // Shouldn't ever need RayKernels for view factors
  params.set<bool>("ray_kernel_coverage_check") = false;
  params.suppressParameter<bool>("ray_kernel_coverage_check");

  // So that the study executes before the RayTracingViewFactor
  params.set<bool>("force_preaux") = true;
  params.suppressParameter<bool>("force_preaux");

  // Need to use internal sidesets
  params.set<bool>("use_internal_sidesets") = true;
  params.suppressParameter<bool>("use_internal_sidesets");

  // May have duplicate external boundaries on one boundary so ignore
  params.set<bool>("external_ray_bc_coverage_check") = false;
  params.suppressParameter<bool>("external_ray_bc_coverage_check");

  // No need to use Ray registration
  params.set<bool>("_use_ray_registration") = false;
  // Do not need to bank Rays on completion
  params.set<bool>("_bank_rays_on_completion") = false;

  return params;
}

ViewFactorRayStudy::ViewFactorRayStudy(const InputParameters & parameters)
  : RayTracingStudy(parameters),
    _bnd_ids(_mesh.getBoundaryIDs(getParam<std::vector<BoundaryName>>("boundary"))),
    _internal_convention(getParam<MooseEnum>("internal_convention")),
    _ray_index_start_bnd_id(registerRayAuxData("start_bnd_id")),
    _ray_index_start_total_weight(registerRayAuxData("start_total_weight")),
    _ray_index_end_bnd_id(registerRayAuxData("end_bnd_id")),
    _ray_index_end_weight(registerRayAuxData("end_weight")),
    _ray_index_start_end_distance(registerRayAuxData("start_end_distance")),
    _fe_face(FEBase::build(_mesh.dimension(), FEType(CONSTANT, MONOMIAL))),
    _q_face(QBase::build(Moose::stringToEnum<QuadratureType>(getParam<MooseEnum>("face_type")),
                         _mesh.dimension() - 1,
                         Moose::stringToEnum<Order>(getParam<MooseEnum>("face_order")))),
    _threaded_vf_info(libMesh::n_threads()),
    _threaded_cached_ray_bcs(libMesh::n_threads())
{
  _fe_face->attach_quadrature_rule(_q_face.get());
  _fe_face->get_normals();
  _fe_face->get_xyz();
}

void
ViewFactorRayStudy::initialSetup()
{
  RayTracingStudy::initialSetup();

  // We optimized away RayKernels, so don't allow them
  std::vector<RayKernelBase *> ray_kernels;
  RayTracingStudy::getRayKernels(ray_kernels, 0);
  if (!ray_kernels.empty())
    mooseError("The ViewFactorRayStudy '", name(), "' is not compatible with RayKernels");

  // Make sure we have only one RayBC: a ViewFactorRayBC with the same boundaries
  std::vector<RayBC *> ray_bcs;
  RayTracingStudy::getRayBCs(ray_bcs, 0);
  if (ray_bcs.size() != 1)
    mooseError(
        "The ViewFactorRayStudy '", name(), "' requires one and only one RayBC, a ViewFactorRayBC");
  ViewFactorRayBC * vf_ray_bc = dynamic_cast<ViewFactorRayBC *>(ray_bcs[0]);
  if (!vf_ray_bc)
    mooseError(
        "The ViewFactorRayStudy '", name(), "' requires one and only one RayBC, a ViewFactorRayBC");
  if (!vf_ray_bc->hasBoundary(_bnd_ids))
    mooseError("The ViewFactorRayBC '",
               vf_ray_bc->name(),
               "' must be applied to the same boundaries as the ViewFactorRayStudy '",
               name(),
               "'");

  // Cache our one RayBC per thread so that we don't spend unnecessary time querying for RayBCs on
  // every boundary
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    RayTracingStudy::getRayBCs(ray_bcs, tid);
    _threaded_cached_ray_bcs[tid] = ray_bcs;
  }
}

void
ViewFactorRayStudy::preExecuteStudy()
{
  RayTracingStudy::preExecuteStudy();

  // Clear and zero the view factor maps we're about to accumulate into for each thread
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    _threaded_vf_info[tid].clear();
    for (const BoundaryID from_id : _bnd_ids)
      for (const BoundaryID to_id : _bnd_ids)
        _threaded_vf_info[tid][from_id][to_id] = 0;
  }

  generatePoints();

  generateRayIDs();
}

void
ViewFactorRayStudy::postExecuteStudy()
{
  RayTracingStudy::postExecuteStudy();

  // Finalize the cumulative _vf_info
  auto before = _threaded_vf_info[0];
  _vf_info.clear();
  for (const BoundaryID from_id : _bnd_ids)
    for (const BoundaryID to_id : _bnd_ids)
      if (from_id <= to_id)
      {
        Real & entry = _vf_info[from_id][to_id];

        // Zero before summing
        entry = 0;

        // Sum over threads
        for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
        {
          entry += _threaded_vf_info[tid][from_id][to_id];
          entry += _threaded_vf_info[tid][to_id][from_id];
        }

        // Sum over processors
        _communicator.sum(entry);

        // Apply reciprocity
        if (from_id != to_id)
          _vf_info[to_id][from_id] = entry;
      }
}

void
ViewFactorRayStudy::addToViewFactorInfo(Real value,
                                        const BoundaryID from_id,
                                        const BoundaryID to_id,
                                        const THREAD_ID tid)
{
  if (!currentlyPropagating() && !currentlyGenerating())
    mooseError("addToViewFactorInfo() can only be called during Ray generation and propagation");
  mooseAssert(_threaded_vf_info[tid].count(from_id),
              "Threaded view factor info does not have from boundary");
  mooseAssert(_threaded_vf_info[tid][from_id].count(to_id),
              "Threaded view factor info does not have from -> to boundary");

  _threaded_vf_info[tid][from_id][to_id] += value;
}

Real
ViewFactorRayStudy::viewFactorInfo(const BoundaryID from_id, const BoundaryID to_id) const
{
  auto it = _vf_info.find(from_id);
  if (it == _vf_info.end())
    mooseError("From boundary id ", from_id, " not in view factor map.");

  auto itt = it->second.find(to_id);
  if (itt == it->second.end())
    mooseError("From boundary id ", from_id, " to boundary_id ", to_id, " not in view factor map.");
  return itt->second;
}

ViewFactorRayStudy &
ViewFactorRayStudy::castFromStudy(const InputParameters & params)
{
  RayTracingStudy * study = params.getCheckedPointerParam<RayTracingStudy *>("_ray_tracing_study");
  ViewFactorRayStudy * vf_study = dynamic_cast<ViewFactorRayStudy *>(study);
  if (!vf_study)
    ::mooseError(params.get<std::string>("_type"),
                 " '",
                 params.get<std::string>("_object_name"),
                 "' must be paired with a ViewFactorRayStudy");
  return *vf_study;
}

void
ViewFactorRayStudy::generatePoints()
{
  _internal_bnd_ids.clear();

  const auto & normals = _fe_face->get_normals();
  const auto & points = _fe_face->get_xyz();
  const auto & weights = _fe_face->get_JxW();

  _start_info.clear();
  _end_info.clear();

  // Starting elements we have that are on the wrong side of an internal boundary
  std::vector<StartElem> send_start_info;
  // Our local end elements
  std::vector<EndElem> local_end_info;

  // Get all possible points on the user defined boundaries on this proc
  for (const BndElement * belem : *_mesh.getBoundaryElementRange())
  {
    const Elem * elem = belem->_elem;
    const auto side = belem->_side;
    const auto bnd_id = belem->_bnd_id;

    // Skip if we don't own you
    if (elem->processor_id() != processor_id())
      continue;

    // Skip if the boundary id isn't one we're looking for
    if (std::find(_bnd_ids.begin(), _bnd_ids.end(), bnd_id) == _bnd_ids.end())
      continue;

    // Sanity check on QGRID not working on some types
    if (_q_face->type() == QGRID && elem->type() == TET4)
      mooseError(
          "Cannot use GRID quadrature type with tetrahedral elements in ViewFactorRayStudy '",
          _name,
          "'");

    // Reinit this face
    _fe_face->reinit(elem, side);
    const auto & normal = normals[0];

    // See if this boundary is internal
    const Elem * neighbor = elem->neighbor_ptr(side);
    if (neighbor)
    {
      if (!neighbor->active())
        mooseError("ViewFactorRayStudy does not work with adaptivity");

      // Mark if this is an internal sideset for future use
      _internal_bnd_ids.insert(bnd_id);

      // With the positive convention, the Rays that we want to spawn from internal boundaries
      // have positive dot products with the outward normal on the side. This is not efficient
      // for the ray tracer, so set them up to be spawned from the other element on this
      // side instead
      if (_internal_convention == 0)
      {
        const unsigned short neighbor_side = neighbor->which_neighbor_am_i(elem);

        // If the neighbor is on another processor, ship this send elem to it instead
        if (neighbor->processor_id() != processor_id())
          send_start_info.emplace_back(neighbor, neighbor_side, bnd_id, -normal, points, weights);
        // If the neighbor is local, add it to our own send info
        else
          _start_info.emplace_back(neighbor, neighbor_side, bnd_id, -normal, points, weights);

        // Use the new elem id/side for the end elem as well
        local_end_info.emplace_back(neighbor->id(), bnd_id, points, weights);

        continue;
      }
    }

    // The boundary is external and we can add it like normal
    _start_info.emplace_back(elem, side, bnd_id, normal, points, weights);
    local_end_info.emplace_back(elem->id(), bnd_id, points, weights);
  }

  // Decide on the internal boundaries
  _communicator.set_union(_internal_bnd_ids);

  // Send the additional start info if any
  if (_internal_convention == 0)
  {
    // Fill the elements to send to each proc based on who owns the elem
    std::unordered_map<processor_id_type, std::vector<StartElem *>> send_start_map;
    for (StartElem & start_elem : send_start_info)
      send_start_map[start_elem._elem->processor_id()].push_back(&start_elem);

    // Communicate and fill into our _start_info
    auto start_functor = [this](processor_id_type, const std::vector<StartElem *> & start_elems) {
      _start_info.reserve(_start_info.size() + start_elems.size());
      for (StartElem * start_elem : start_elems)
        _start_info.emplace_back(std::move(*start_elem));
    };
    Parallel::push_parallel_packed_range(_communicator, send_start_map, this, start_functor);
  }

  // Decide on the processors that have boundary elements so that we know who to send to
  std::vector<dof_id_type> have_boundary(_communicator.size());
  _communicator.allgather((dof_id_type)_start_info.size(), have_boundary);

  // The data we're sending to each proc - here we have pointers because the pack/unpack
  // routines are optimized to copy as little data as possible
  std::vector<EndElem *> send_end_elems;
  send_end_elems.reserve(local_end_info.size());
  for (auto & end_elem : local_end_info)
    send_end_elems.push_back(&end_elem);

  // Fill the send data to each proc that needs it
  std::unordered_map<processor_id_type, std::vector<EndElem *>> send_end_map;
  for (processor_id_type pid = 0; pid < _communicator.size(); ++pid)
    if (pid != _pid && have_boundary[pid])
      send_end_map[pid] = send_end_elems;

  // Communciate and fill into our _end_info
  auto end_functor = [this](processor_id_type, const std::vector<EndElem *> & end_elems) {
    _end_info.reserve(_end_info.size() + end_elems.size());
    for (EndElem * end_elem : end_elems)
      _end_info.emplace_back(std::move(*end_elem));
  };
  Parallel::push_parallel_packed_range(_communicator, send_end_map, (void *)nullptr, end_functor);

  // _end_info at this point only contains remote info - add ours too. We add these last
  // so that they get generated last in order to prioritize Rays that need to leave our proc
  _end_info.reserve(_end_info.size() + local_end_info.size());
  for (EndElem & end_elem : local_end_info)
    _end_info.emplace_back(std::move(end_elem));
}

bool
ViewFactorRayStudy::shouldGenerate(const dof_id_type start_elem_id,
                                   const dof_id_type end_elem_id,
                                   const BoundaryID start_bnd_id,
                                   const BoundaryID end_bnd_id) const
{
  if (start_elem_id == end_elem_id)
    return start_bnd_id < end_bnd_id;

  const bool start_even = start_elem_id % 2 == 0;
  const bool end_even = end_elem_id % 2 == 0;

  if (start_even && end_even)
    return start_elem_id > end_elem_id;
  if (!start_even && !end_even)
    return end_elem_id > start_elem_id;
  else
    return start_even;
}

void
ViewFactorRayStudy::generateRayIDs()
{
  // Determine number of Rays and points
  std::size_t num_local_rays = 0;
  std::size_t num_local_start_points = 0;
  for (const auto & start_elem : _start_info)
  {
    const auto start_elem_id = start_elem._elem->id();
    const auto start_bnd_id = start_elem._bnd_id;
    const auto & start_normal = start_elem._normal;

    num_local_start_points += start_elem._points.size();

    for (const auto & end_elem : _end_info)
    {
      const auto end_elem_id = end_elem._elem_id;
      const auto end_bnd_id = end_elem._bnd_id;
      if (!shouldGenerate(start_elem_id, end_elem_id, start_bnd_id, end_bnd_id))
        continue;

      for (const auto & start_point : start_elem._points)
        for (const auto & end_point : end_elem._points)
        {
          if (start_point.absolute_fuzzy_equals(end_point))
            continue;

          const Point direction = (end_point - start_point).unit();
          const Real dot = start_normal * direction;
          if (dot < -1e-6)
            ++num_local_rays;
        }
    }
  }

  // Decide on the starting Ray IDs for each rank. Each rank will have a contiguous set of IDs for
  // its rays and there will be no holes in the IDs between ranks Number of rays this proc has
  // Vector of sizes across all processors (for rank 0 only)
  std::vector<std::size_t> proc_sizes;
  // Send this procesor's local size to rank 0
  _comm.gather(0, num_local_rays, proc_sizes);
  // Rank 0 has the proc sizes and will decide on a starting ID for every rank
  std::vector<RayID> proc_starting_id;
  if (processor_id() == 0)
  {
    proc_starting_id.resize(_comm.size());
    RayID current_starting_id = 0;
    for (processor_id_type pid = 0; pid < _comm.size(); ++pid)
    {
      proc_starting_id[pid] = current_starting_id;
      current_starting_id += proc_sizes[pid];
    }
  }

  // Send the starting ID to each rank
  _comm.scatter(proc_starting_id, _current_starting_id, 0);
  if (!num_local_rays)
    _current_starting_id = Ray::INVALID_RAY_ID;

  // Print out totals while we're here
  std::size_t num_total_points = num_local_start_points;
  std::size_t num_total_rays = num_local_rays;
  _communicator.sum(num_total_points);
  _communicator.sum(num_total_rays);
  _console << "View factor study generated " << num_total_points << " points requiring "
           << num_total_rays << " rays" << std::endl;
}

void
ViewFactorRayStudy::generateRays()
{
  CONSOLE_TIMED_PRINT("View factor study generating rays");

  for (const auto & start_elem : _start_info)
  {
    const Elem * elem = start_elem._elem;
    const auto start_elem_id = start_elem._elem->id();
    const auto side = start_elem._side;
    const auto start_bnd_id = start_elem._bnd_id;
    const auto & start_normal = start_elem._normal;

    for (const auto & end_elem : _end_info)
    {
      const auto end_elem_id = end_elem._elem_id;
      const auto end_bnd_id = end_elem._bnd_id;

      if (!shouldGenerate(start_elem_id, end_elem_id, start_bnd_id, end_bnd_id))
        continue;

      const bool end_is_internal = _internal_bnd_ids.count(end_bnd_id);

      for (std::size_t start_i = 0; start_i < start_elem._points.size(); ++start_i)
      {
        const auto & start_point = start_elem._points[start_i];
        const auto start_weight = start_elem._weights[start_i];

        for (std::size_t end_i = 0; end_i < end_elem._points.size(); ++end_i)
        {
          const auto & end_point = end_elem._points[end_i];
          const auto end_weight = end_elem._weights[end_i];

          if (start_point.absolute_fuzzy_equals(end_point))
            continue;

          const auto direction = (end_point - start_point).unit();
          const Real dot = start_normal * direction;
          if (dot > -1e-6)
            continue;

          Point mock_end_point = end_point;
          if (end_is_internal)
            mock_end_point += direction;

          // Create a Ray and add it to the buffer for future tracing
          std::shared_ptr<Ray> ray = _ray_pool.acquire(
              start_point, mock_end_point, rayDataSize(), rayAuxDataSize(), elem, side);
          addToWorkingBuffer(ray);

          ray->setID(_current_starting_id++);
          ray->setAuxData(_ray_index_start_bnd_id, start_bnd_id);
          ray->setAuxData(_ray_index_start_total_weight, start_weight * std::abs(dot));
          ray->setAuxData(_ray_index_end_bnd_id, end_bnd_id);
          ray->setAuxData(_ray_index_end_weight, end_weight);

          // Note that ray->end() is a fake end point. We must do this due to Rays that end on
          // internal sidesets. These Rays cannot have ray->end() = end_point, because in the
          // case where the Ray ends on the ending side but not on the ending internal side
          // (an internal side is associated with a single element, and the Ray may hit the
          // neighbor first), the Ray would never hit the internal side! It would die first.
          // Therefore, we will use the distance from start -> actual end to see if we've
          // hit the actual end point.
          ray->setAuxData(_ray_index_start_end_distance, (start_point - end_point).norm());

          ray.reset();

          chunkyTraceAndBuffer();
        }
      }
    }
  }
}

namespace libMesh
{
namespace Parallel
{

unsigned int
Packing<ViewFactorRayStudy::StartElem *>::packed_size(typename std::vector<Real>::const_iterator in)
{
  unsigned int total_size = 0;

  // Number of points
  unsigned int num_points = *in++;

  // Number of points, elem, side, bnd_id, normal
  total_size += 7;
  // Points and weights
  total_size += num_points * 4;

  return total_size;
}

unsigned int
Packing<ViewFactorRayStudy::StartElem *>::packable_size(
    const ViewFactorRayStudy::StartElem * const start_elem, const void *)
{
  unsigned int total_size = 0;

  // Number of points, elem, side, bnd_id, normal
  total_size += 7;
  // Points and weights
  total_size += start_elem->_points.size() * 4;

  return total_size;
}

template <>
ViewFactorRayStudy::StartElem *
Packing<ViewFactorRayStudy::StartElem *>::unpack(std::vector<Real>::const_iterator in,
                                                 ViewFactorRayStudy * study)
{
  // Number of points
  const std::size_t num_points = *in++;

  // Elem id
  const dof_id_type elem_id = *in++;
  const Elem * elem =
      elem_id != DofObject::invalid_id ? study->meshBase().query_elem_ptr(elem_id) : nullptr;

  // Side
  const unsigned short side = *in++;

  // Boundary ID
  const BoundaryID bnd_id = *in++;

  // Normal
  Point normal;
  normal(0) = *in++;
  normal(1) = *in++;
  normal(2) = *in++;

  // Points
  std::vector<Point> points(num_points);
  for (std::size_t i = 0; i < num_points; ++i)
  {
    points[i](0) = *in++;
    points[i](1) = *in++;
    points[i](2) = *in++;
  }

  // Weights
  std::vector<Real> weights(num_points);
  for (std::size_t i = 0; i < num_points; ++i)
    weights[i] = *in++;

  // We store this with std::move, so this is safe if we do it right!
  return new ViewFactorRayStudy::StartElem(
      elem, side, bnd_id, std::move(normal), std::move(points), std::move(weights));
}

template <>
void
Packing<ViewFactorRayStudy::StartElem *>::pack(
    const ViewFactorRayStudy::StartElem * const start_elem,
    std::back_insert_iterator<std::vector<Real>> data_out,
    const ViewFactorRayStudy *)
{
  // Number of points
  data_out = start_elem->_points.size();

  // Elem id
  data_out = start_elem->_elem->id();

  // Side
  data_out = start_elem->_side;

  // Boundary id
  data_out = start_elem->_bnd_id;

  // Normal
  data_out = start_elem->_normal(0);
  data_out = start_elem->_normal(1);
  data_out = start_elem->_normal(2);

  // Points
  for (const auto & point : start_elem->_points)
  {
    data_out = point(0);
    data_out = point(1);
    data_out = point(2);
  }

  // Weights
  std::copy(start_elem->_weights.begin(), start_elem->_weights.end(), data_out);
}

unsigned int
Packing<ViewFactorRayStudy::EndElem *>::packed_size(typename std::vector<Real>::const_iterator in)
{
  unsigned int total_size = 0;

  // Number of points
  unsigned int num_points = *in++;

  // Number of points, elem id, bnd_id
  total_size += 3;
  // Points and weights
  total_size += num_points * 4;

  return total_size;
}

unsigned int
Packing<ViewFactorRayStudy::EndElem *>::packable_size(
    const ViewFactorRayStudy::EndElem * const end_elem, const void *)
{
  unsigned int total_size = 0;

  // Number of points, elem_id, bnd_id
  total_size += 3;
  // Points and weights
  total_size += end_elem->_points.size() * 4;

  return total_size;
}

template <>
ViewFactorRayStudy::EndElem *
Packing<ViewFactorRayStudy::EndElem *>::unpack(std::vector<Real>::const_iterator in, void *)
{
  // Number of points
  const std::size_t num_points = *in++;

  // Elem id
  const dof_id_type elem_id = *in++;

  // Boundary ID
  const BoundaryID bnd_id = *in++;

  // Points
  std::vector<Point> points(num_points);
  for (std::size_t i = 0; i < num_points; ++i)
  {
    points[i](0) = *in++;
    points[i](1) = *in++;
    points[i](2) = *in++;
  }

  // Weights
  std::vector<Real> weights(num_points);
  for (std::size_t i = 0; i < num_points; ++i)
    weights[i] = *in++;

  // We store this with std::move, so this is safe if we do it right!
  return new ViewFactorRayStudy::EndElem(elem_id, bnd_id, std::move(points), std::move(weights));
}

template <>
void
Packing<ViewFactorRayStudy::EndElem *>::pack(const ViewFactorRayStudy::EndElem * const end_elem,
                                             std::back_insert_iterator<std::vector<Real>> data_out,
                                             const void *)
{
  // Number of points
  data_out = end_elem->_points.size();

  // Elem id
  data_out = end_elem->_elem_id;

  // Boundary id
  data_out = end_elem->_bnd_id;

  // Points
  for (const auto & point : end_elem->_points)
  {
    data_out = point(0);
    data_out = point(1);
    data_out = point(2);
  }

  // Weights
  std::copy(end_elem->_weights.begin(), end_elem->_weights.end(), data_out);
}

} // namespace Parallel

} // namespace libMesh
