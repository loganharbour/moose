//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ViewFactorRayStudy.h"

// libMesh includes
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel_sync.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature.h"

#include "SerializerGuard.h"

registerMooseObject("RayTracingApp", ViewFactorRayStudy);

InputParameters
ViewFactorRayStudy::validParams()
{
  auto params = RayTracingStudy::validParams();

  params.addRequiredParam<std::vector<BoundaryName>>(
      "boundary", "The list of boundary IDs from the mesh where this boundary condition applies");

  MooseEnum qorders("CONSTANT FIRST SECOND", "CONSTANT");
  params.addParam<MooseEnum>("face_order", qorders, "The face quadrature rule order");

  return params;
}

ViewFactorRayStudy::ViewFactorRayStudy(const InputParameters & parameters)
  : RayTracingStudy(parameters),
    _bnd_ids(_mesh.getBoundaryIDs(getParam<std::vector<BoundaryName>>("boundary"))),
    _ray_index_start_dot(registerRayAuxData("start_dot")),
    _ray_index_start_bnd_id(registerRayAuxData("start_bnd_id")),
    _ray_index_start_weight(registerRayAuxData("start_weight")),
    _ray_index_end_weight(registerRayAuxData("end_weight")),
    _fe_face(FEBase::build(_mesh.dimension(), FEType(FIRST, LAGRANGE))),
    _q_face(QBase::build(QGAUSS,
                         _mesh.dimension() - 1,
                         Moose::stringToEnum<Order>(getParam<MooseEnum>("face_order")))),
    _vf_info(libMesh::n_threads())
{
  _use_ray_registration = false;

  _fe_face->attach_quadrature_rule(_q_face.get());
  _fe_face->get_normals();
  _fe_face->get_xyz();
}

void
ViewFactorRayStudy::initialize()
{
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
    _vf_info[tid].clear();
}

Real &
ViewFactorRayStudy::viewFactorInfo(BoundaryID from_id, BoundaryID to_id, THREAD_ID tid)
{
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
        mooseError(
            "This should never happen, from boundary id can never be larger than to boundary id.");

      _vf_info[0][to][from] = v;
    }
  }
}

void
ViewFactorRayStudy::generatePoints()
{
  _start_info.resize(_bnd_ids.size());
  _end_points.resize(_bnd_ids.size());

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

    // The location in the _bnd_ids vector that bnd_id is
    const auto bnd_id_index = find_it - _bnd_ids.begin();

    // Reinit this face
    _fe_face->reinit(elem, side);

    // Add to the starting info
    _start_info[bnd_id_index].emplace_back(
        elem, side, _fe_face->get_normals()[0], _fe_face->get_xyz(), _fe_face->get_JxW());
  }

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
}

void
ViewFactorRayStudy::defineRays()
{
  // Start with an ID for Rays that is unique to this processor
  _next_id = DofObject::invalid_id / (dof_id_type)_comm.size() * (dof_id_type)_pid;

  // Generate the start and end points
  generatePoints();

  for (unsigned int start_bnd_id_index = 0; start_bnd_id_index < _bnd_ids.size();
       ++start_bnd_id_index)
  {
    // Don't have any start points on this bounary
    if (_start_info[start_bnd_id_index].empty())
      continue;

    const auto start_bnd_id = _bnd_ids[start_bnd_id_index];

    // Loop through all end boundaries to see which we need to send Rays to
    for (unsigned int end_bnd_id_index = 0; end_bnd_id_index < _bnd_ids.size(); ++end_bnd_id_index)
    {
      const auto end_bnd_id = _bnd_ids[end_bnd_id_index];

      // Only send rays when start_bnd_id <= end_bnd_id
      if (start_bnd_id > end_bnd_id)
        continue;

      // Loop through all elems on said boundary
      for (const auto & start_elem : _start_info[start_bnd_id_index])
      {
        const auto elem = start_elem._elem;
        const auto side = start_elem._side;
        const auto normal = start_elem._normal;
        const auto & points = start_elem._points;
        const auto & weights = start_elem._weights;

        // Loop over all start points associated with this boundary
        for (unsigned int i = 0; i < points.size(); ++i)
        {
          const auto start_point = points[i];
          const auto start_weight = weights[i];

          // And spawn from the start point to all of the end points
          for (const auto & end_point_pair : _end_points[end_bnd_id_index])
          {
            const auto end_point = end_point_pair.first;
            const auto end_weight = end_point_pair.second;

            // Don't spawn when start == end
            if (start_bnd_id == end_bnd_id && start_point.absolute_fuzzy_equals(end_point))
              continue;

            defineRay(
                elem, start_point, end_point, normal, side, start_bnd_id, start_weight, end_weight);
          }
        }
      }
    }
  }
}

void
ViewFactorRayStudy::defineRay(const Elem * starting_elem,
                              const Point & start_point,
                              const Point & end_point,
                              const Point & normal,
                              const unsigned short side,
                              const BoundaryID bnd_id,
                              const Real start_weight,
                              const Real end_weight)
{
  std::shared_ptr<Ray> ray = _ray_pool.acquire();
  _working_buffer->push_back(ray);

  ray->setStartingElem(starting_elem);
  ray->setIncomingSide(side);
  ray->setID(_next_id++);

  const Point direction = (end_point - start_point).unit();
  ray->setStart(start_point);
  ray->setEnd(end_point);
  ray->setDirection(direction);

  ray->auxData().resize(rayAuxDataSize());
  // For computing the dot product, inward normals are assumed. We just
  // use the absolute value and don't have to worry about it
  ray->setAuxData(_ray_index_start_dot, std::abs(normal * direction));
  ray->setAuxData(_ray_index_start_bnd_id, bnd_id);
  ray->setAuxData(_ray_index_start_weight, start_weight);
  ray->setAuxData(_ray_index_end_weight, end_weight);
}

void
ViewFactorRayStudy::generateRays()
{
  defineRays();

  // Increment required variable when adding rays
  _local_rays_started += _working_buffer->size();

  // And spawn
  if (_method == SMART)
    chunkyTraceAndBuffer();
  else
  {
    traceAndBuffer(_working_buffer->begin(), _working_buffer->end());
    _working_buffer->clear();
  }
}
