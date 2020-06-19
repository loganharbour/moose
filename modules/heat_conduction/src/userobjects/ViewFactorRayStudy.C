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

  // Shouldn't ever need RayKernels for view factors
  params.set<bool>("ray_kernel_coverage_check") = false;
  params.suppressParameter<bool>("ray_kernel_coverage_check");

  // So that the study executes before the RayTracingViewFactor
  params.set<bool>("force_preaux") = true;
  params.suppressParameter<bool>("force_preaux");

  // Need to use internal sidesets
  params.set<bool>("use_internal_sidesets") = true;
  params.suppressParameter<bool>("use_internal_sidesets");

  // No need to use Ray registration
  params.set<bool>("_use_ray_registration") = false;
  // Do not need to bank Rays on completion
  params.set<bool>("_bank_rays_on_completion") = false;

  return params;
}

ViewFactorRayStudy::ViewFactorRayStudy(const InputParameters & parameters)
  : RayTracingStudy(parameters),
    _bnd_ids(_mesh.getBoundaryIDs(getParam<std::vector<BoundaryName>>("boundary"))),
    _ray_index_start_bnd_id(registerRayAuxData("start_bnd_id")),
    _ray_index_start_total_weight(registerRayAuxData("start_total_weight")),
    _ray_index_end_bnd_id(registerRayAuxData("end_bnd_id")),
    _ray_index_end_weight(registerRayAuxData("end_weight")),
    _ray_index_start_end_distance(registerRayAuxData("start_end_distance")),
    _fe_face(FEBase::build(_mesh.dimension(), FEType(CONSTANT, MONOMIAL))),
    _q_face(QBase::build(QGRID,
                         _mesh.dimension() - 1,
                         Moose::stringToEnum<Order>(getParam<MooseEnum>("face_order")))),
    _vf_info(libMesh::n_threads()),
    _threaded_cached_ray_bcs(libMesh::n_threads())
{
  _fe_face->attach_quadrature_rule(_q_face.get());
  _fe_face->get_normals();
  _fe_face->get_xyz();
}

void
ViewFactorRayStudy::initialize()
{
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    _vf_info[tid].clear();

    for (const BoundaryID from_id : _bnd_ids)
      for (const BoundaryID to_id : _bnd_ids)
        if (from_id <= to_id)
          _vf_info[tid][from_id][to_id] = 0;
  }
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

Real &
ViewFactorRayStudy::viewFactorInfo(BoundaryID from_id, BoundaryID to_id, THREAD_ID tid)
{
  mooseAssert(std::find(_bnd_ids.begin(), _bnd_ids.end(), from_id) != _bnd_ids.end(),
              "Invalid from_id");
  mooseAssert(std::find(_bnd_ids.begin(), _bnd_ids.end(), to_id) != _bnd_ids.end(),
              "Invalid to_id");

  return _vf_info[tid][from_id][to_id];
}

Real
ViewFactorRayStudy::viewFactorInfo(BoundaryID from_id, BoundaryID to_id) const
{
  // this function should be called after summation over all
  // threads
  auto it = _vf_info[0].find(from_id);
  if (it == _vf_info[0].end())
    mooseError("From id ", from_id, " not in view factor map.");

  auto itt = it->second.find(to_id);
  if (itt == it->second.end())
    mooseError("From boundary id ", from_id, " to boundary_id ", to_id, " not in view factor map.");
  return itt->second;
}

void
ViewFactorRayStudy::finalize()
{
  // accumulation of _vf_info for all threads into copy 0
  for (THREAD_ID tid = 1; tid < libMesh::n_threads(); ++tid)
    for (auto & p : _vf_info[tid])
      for (auto & pp : p.second)
        _vf_info[0][p.first][pp.first] += pp.second;

  // first all processors must have the same map_of_map entries
  std::set<std::pair<BoundaryID, BoundaryID>> bnd_id_pairs;
  for (auto & p : _vf_info[0])
    for (auto & pp : p.second)
      bnd_id_pairs.insert(std::pair<BoundaryID, BoundaryID>(p.first, pp.first));
  comm().set_union(bnd_id_pairs);

  for (auto & p : bnd_id_pairs)
  {
    BoundaryID from = p.first;
    BoundaryID to = p.second;
    auto it = _vf_info[0].find(from);
    if (it == _vf_info[0].end() || it->second.find(to) == it->second.end())
      _vf_info[0][from][to] = 0;
  }

  // now accumulate over all processors
  for (auto & p : _vf_info[0])
    for (auto & pp : p.second)
      gatherSum(_vf_info[0][p.first][pp.first]);

  // Apply reciprocity
  // NOTE we modify the map so we need to make a copy
  // here to avoid the problem of looping over entries
  // we just created
  auto vf = _vf_info[0];
  for (auto & p : vf)
  {
    BoundaryID from = p.first;
    for (auto & pp : p.second)
    {
      BoundaryID to = pp.first;
      Real v = pp.second;

      if (from > to)
        mooseError("This should never happen, from boundary id can never be larger than to "
                   "boundary id.");

      _vf_info[0][to][from] = v;
    }
  }
}

void
ViewFactorRayStudy::generatePoints()
{
  _start_info.resize(_bnd_ids.size());
  _end_points.resize(_bnd_ids.size());
  _internal_bnd_ids.clear();

  // Get all possible centroids on the user defined boundaries on this proc
  for (const auto & belem : *_mesh.getBoundaryElementRange())
  {
    const Elem * elem = belem->_elem;
    const auto side = belem->_side;
    const auto bnd_id = belem->_bnd_id;

    // Skip if we don't own you
    if (elem->processor_id() != processor_id())
      continue;

    // Skip if the boundary id isn't one we're looking for
    const auto find_it = std::find(_bnd_ids.begin(), _bnd_ids.end(), bnd_id);
    if (find_it == _bnd_ids.end())
      continue;

    // Mark if this is an internal sideset for future use
    if (elem->neighbor_ptr(side))
      _internal_bnd_ids.insert(bnd_id);

    // The location in the _bnd_ids vector that bnd_id is
    const auto bnd_id_index = std::distance(_bnd_ids.begin(), find_it);

    // Reinit this face
    _fe_face->reinit(elem, side);

    // Add to the starting info
    _start_info[bnd_id_index].emplace_back(
        elem, side, _fe_face->get_normals()[0], _fe_face->get_xyz(), _fe_face->get_JxW());
  }

  // Decide on the internal boundaries
  _communicator.set_union(_internal_bnd_ids);

  // The processors that have elements on each boundary (indexed by indexing in _bnd_ids)
  std::vector<std::vector<processor_id_type>> bnd_id_procs(_bnd_ids.size());
  /// To be filled by each proc for each boundary
  std::vector<processor_id_type> have;
  /// Find out which processors have elements on which boundary
  for (unsigned int bnd_id_index = 0; bnd_id_index < _bnd_ids.size(); ++bnd_id_index)
  {
    have.clear();

    // If we have this boundary, fill our entry with our pid
    if (!_start_info[bnd_id_index].empty())
      have.push_back(_pid);
    // And if we don't fill with an invalid pid
    else
      have.push_back(libMesh::DofObject::invalid_processor_id);
    _comm.allgather(have, true);

    // All entries that aren't invalid are someone that has elems on this bnd_id
    for (const auto & pid : have)
      if (pid != libMesh::DofObject::invalid_processor_id)
        bnd_id_procs[bnd_id_index].push_back(pid);
  }

  // Fill who we need to send points (boundary ID, point, and weight) to
  std::unordered_map<processor_id_type, std::vector<std::tuple<unsigned int, Point, Real>>>
      send_tuples;
  for (unsigned int my_bnd_id_index = 0; my_bnd_id_index < _bnd_ids.size(); ++my_bnd_id_index)
  {
    // No points to send
    if (_start_info[my_bnd_id_index].empty())
      continue;

    const auto my_bnd_id = _bnd_ids[my_bnd_id_index];
    const auto & start_elems = _start_info[my_bnd_id_index];

    // Find procs that we need to send these boundary points to
    std::set<processor_id_type> send_to;
    for (unsigned int to_bnd_id_index = 0; to_bnd_id_index < _bnd_ids.size(); ++to_bnd_id_index)
    {
      const auto to_bnd_id = _bnd_ids[to_bnd_id_index];
      if (to_bnd_id <= my_bnd_id)
        for (const auto to_pid : bnd_id_procs[to_bnd_id_index])
          send_to.insert(to_pid);
    }

    // Fill to send for each of these procs
    for (const auto & to_pid : send_to)
    {
      auto & send_entry = send_tuples[to_pid];
      send_entry.reserve(send_entry.size() + start_elems.size());
      for (const auto & start_elem : start_elems)
        for (unsigned int i = 0; i < start_elem._points.size(); ++i)
          send_entry.emplace_back(my_bnd_id_index, start_elem._points[i], start_elem._weights[i]);
    }
  }

  // And send those points
  _end_points.clear();
  auto receive_points_functor =
      [this](processor_id_type /* pid */,
             const std::vector<std::tuple<unsigned int, Point, Real>> & tuples) {
        for (const auto & tuple : tuples)
        {
          const auto bnd_id_index = std::get<0>(tuple);
          const auto & point = std::get<1>(tuple);
          const auto weight = std::get<2>(tuple);

          _end_points[bnd_id_index].emplace_back(point, weight);
        }
      };
  Parallel::push_parallel_vector_data(_comm, send_tuples, receive_points_functor);

  // Determine number of Rays and points
  std::size_t num_local_start_points = 0;
  std::size_t num_local_rays = 0;
  for (unsigned int start_bnd_id_i = 0; start_bnd_id_i < _bnd_ids.size(); ++start_bnd_id_i)
    for (const auto & start_elem : _start_info[start_bnd_id_i])
    {
      num_local_start_points += start_elem._points.size();
      for (unsigned int end_bnd_id_i = 0; end_bnd_id_i < _bnd_ids.size(); ++end_bnd_id_i)
        if (_bnd_ids[start_bnd_id_i] <= _bnd_ids[end_bnd_id_i])
          for (const auto & start_point : start_elem._points)
            for (const auto & end_point_pair : _end_points[end_bnd_id_i])
            {
              const auto & end_point = end_point_pair.first;
              if (start_point.absolute_fuzzy_equals(end_point))
                continue;
              const Point direction = (end_point - start_point).unit();
              if (start_elem._normal * direction < -TOLERANCE)
                ++num_local_rays;
            }
    }

  // While we're at it, set the work buffer capacity ahead of time
  setWorkingBufferCapacity(num_local_rays);

  // Print out totals while we're here
  std::size_t num_total_points = num_local_start_points;
  std::size_t num_total_rays = num_local_rays;
  _communicator.sum(num_total_points);
  _communicator.sum(num_total_rays);
  _console << "View factor study generated " << num_total_points << " points requiring "
           << num_total_rays << " rays" << std::endl;

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
}

void
ViewFactorRayStudy::generateRays()
{
  // Generate the start and end points
  generatePoints();

  CONSOLE_TIMED_PRINT("View factor study generating rays");

  for (unsigned int start_bnd_id_i = 0; start_bnd_id_i < _bnd_ids.size(); ++start_bnd_id_i)
  {
    // Don't have any start points on this bounary
    if (_start_info[start_bnd_id_i].empty())
      continue;

    const auto start_bnd_id = _bnd_ids[start_bnd_id_i];

    // Loop through all end boundaries to see which we need to send Rays to
    for (unsigned int end_bnd_id_i = 0; end_bnd_id_i < _bnd_ids.size(); ++end_bnd_id_i)
    {
      const auto end_bnd_id = _bnd_ids[end_bnd_id_i];
      const bool end_is_internal = _internal_bnd_ids.count(end_bnd_id);

      // Only send rays when start_bnd_id <= end_bnd_id
      if (start_bnd_id > end_bnd_id)
        continue;

      // Loop through all elems on said boundary
      for (const auto & start_elem : _start_info[start_bnd_id_i])
      {
        const auto elem = start_elem._elem;
        const auto side = start_elem._side;
        const auto normal = start_elem._normal;
        const auto & points = start_elem._points;
        const auto & weights = start_elem._weights;

        // Loop over all start points associated with this boundary
        for (unsigned int i = 0; i < points.size(); ++i)
        {
          const auto & start_point = points[i];
          const auto start_weight = weights[i];

          // And spawn from the start point to all of the end points
          for (const auto & end_point_pair : _end_points[end_bnd_id_i])
          {
            const auto & end_point = end_point_pair.first;
            const auto end_weight = end_point_pair.second;

            // Don't spawn when start == end
            if (start_point.absolute_fuzzy_equals(end_point))
              continue;

            // Direction from start -> end
            const Point direction = (end_point - start_point).unit();
            // Dot product with direction and normal: Only keep when dot < -TOLERANCE
            // dot = (0, 1] is the wrong direction and [-TOLERANCE, 0] implies in the same plane
            const Real dot = normal * direction;
            if (dot > -TOLERANCE)
              continue;

            // For Rays that could possibly end on an internal boundary, the end point we set on the
            // Ray must be past where it may end on the internal boundary. This is a requirement of
            // the ray tracer. Therefore, if the end boundary is internal, bump the end point a
            // little in the direction of the Ray.
            Point mock_end_point = end_point;
            if (end_is_internal)
              mock_end_point += direction;

            // Create a Ray and add it to the buffer for future tracing
            std::shared_ptr<Ray> ray = _ray_pool.acquire();
            addToWorkingBuffer(ray);

            ray->setID(_current_starting_id++);
            ray->setStartingElem(elem);
            ray->setIncomingSide(side);
            ray->setStart(start_point);
            ray->setEnd(mock_end_point);
            ray->setDirection(direction);

            ray->auxData().resize(rayAuxDataSize());
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
          }
        }
      }
    }

    // Done creating Rays from this start boundary
    _start_info[start_bnd_id_i].clear();
  }

  // Clear the point info now that we're done with it
  _start_info.clear();
  _end_points.clear();
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
