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

  return params;
}

ViewFactorRayStudy::ViewFactorRayStudy(const InputParameters & parameters)
  : RayTracingStudy(parameters),
    _bnd_ids(_mesh.getBoundaryIDs(getParam<std::vector<BoundaryName>>("boundary"))),
    _ray_index_dot(registerRayAuxData("dot")),
    _ray_index_bnd_id(registerRayAuxData("bnd_id")),
    _fe_face(FEBase::build(_mesh.dimension(), FEType(FIRST, LAGRANGE))),
    _q_face(QBase::build(QGAUSS, _mesh.dimension() - 1, FIRST))
{
  _use_ray_registration = false;

  _fe_face->attach_quadrature_rule(_q_face.get());
  _fe_face->get_normals();
}

void
ViewFactorRayStudy::generatePoints()
{
  _start_points.resize(_bnd_ids.size());
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

    const auto centroid = elem->side_ptr(side)->centroid();
    _start_points[bnd_id_index].emplace_back(centroid, elem, side);
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
    if (!_start_points[bnd_id_index].empty())
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

  // Fill who we need to send points to
  std::unordered_map<processor_id_type, std::vector<std::pair<unsigned int, Point>>> send_pairs;
  for (unsigned int my_bnd_id_index = 0; my_bnd_id_index < _bnd_ids.size(); ++my_bnd_id_index)
  {
    // No points to send
    if (_start_points[my_bnd_id_index].empty())
      continue;

    const auto my_bnd_id = _bnd_ids[my_bnd_id_index];
    const auto & my_point_tups = _start_points[my_bnd_id_index];

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
      auto & send_entry = send_pairs[to_pid];
      send_entry.reserve(send_entry.size() + my_point_tups.size());
      for (const auto & my_point_tup : my_point_tups)
        send_entry.emplace_back(my_bnd_id_index, std::get<0>(my_point_tup));
    }
  }

  // And send those points
  _end_points.clear();
  auto receive_points_functor = [this](processor_id_type /* pid */,
                                       const std::vector<std::pair<unsigned int, Point>> & pairs) {
    for (const auto & pair : pairs)
    {
      const auto bnd_id_index = pair.first;
      const auto & point = pair.second;
      _end_points[bnd_id_index].push_back(point);
    }
  };
  Parallel::push_parallel_vector_data(_comm, send_pairs, receive_points_functor);
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
    if (_start_points[start_bnd_id_index].empty())
      continue;

    const auto start_bnd_id = _bnd_ids[start_bnd_id_index];

    // Loop through all end boundaries to see which we need to send Rays to
    for (unsigned int end_bnd_id_index = 0; end_bnd_id_index < _bnd_ids.size(); ++end_bnd_id_index)
    {
      const auto end_bnd_id = _bnd_ids[end_bnd_id_index];

      // Only send rays when start_bnd_id <= end_bnd_id
      if (start_bnd_id > end_bnd_id)
        continue;

      // Loop through all start points on said boundary
      for (const auto & start_point_tup : _start_points[start_bnd_id_index])
      {
        const auto start_point = std::get<0>(start_point_tup);
        const auto start_elem = std::get<1>(start_point_tup);
        const auto incoming_side = std::get<2>(start_point_tup);
        const auto normal = sideNormal(start_elem, incoming_side);

        // ...and all end points on said boundary
        for (const auto & end_point : _end_points[end_bnd_id_index])
        {
          // Don't spawn when start == end
          if (start_bnd_id == end_bnd_id && start_point.absolute_fuzzy_equals(end_point))
            continue;

          defineRay(start_elem, start_point, end_point, normal, incoming_side, start_bnd_id);
        }
      }
    }
  }
}

Point
ViewFactorRayStudy::sideNormal(const Elem * elem, const unsigned short side)
{
  _fe_face->reinit(elem, side);
  return _fe_face->get_normals()[0];
}

void
ViewFactorRayStudy::defineRay(const Elem * starting_elem,
                              const Point & start_point,
                              const Point & end_point,
                              const Point & normal,
                              const unsigned short side,
                              const BoundaryID bnd_id)
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
  ray->setAuxData(_ray_index_dot, -normal * direction);
  ray->setAuxData(_ray_index_bnd_id, bnd_id);
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
