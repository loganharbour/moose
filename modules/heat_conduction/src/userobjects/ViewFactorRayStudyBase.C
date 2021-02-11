// //* This file is part of the MOOSE framework
// //* https://www.mooseframework.org
// //*
// //* All rights reserved, see COPYRIGHT for full restrictions
// //* https://github.com/idaholab/moose/blob/master/COPYRIGHT
// //*
// //* Licensed under LGPL 2.1, please see LICENSE for details
// //* https://www.gnu.org/licenses/lgpl-2.1.html
//
// #include "ViewFactorRayStudyBase.h"
//
// // MOOSE includes
// #include "RayTracingPackingUtils.h"
//
// // libMesh includes
// #include "libmesh/parallel_algebra.h"
// #include "libmesh/parallel_sync.h"
// #include "libmesh/enum_quadrature_type.h"
// #include "libmesh/fe_base.h"
// #include "libmesh/quadrature.h"
//
// InputParameters
// ViewFactorRayStudyBase::validParams()
// {
//   auto params = RayTracingStudy::validParams();
//
//   params.addRequiredParam<std::vector<BoundaryName>>(
//       "boundary", "The list of boundaries where view factors are desired");
//
//   MooseEnum qorders("CONSTANT FIRST SECOND THIRD FOURTH FIFTH SIXTH SEVENTH EIGHTH NINTH TENTH "
//                     "ELEVENTH TWELFTH THIRTEENTH FOURTEENTH FIFTEENTH SIXTEENTH SEVENTEENTH "
//                     "EIGHTTEENTH NINTEENTH TWENTIETH",
//                     "CONSTANT");
//   params.addParam<MooseEnum>("face_order", qorders, "The face quadrature rule order");
//
//   MooseEnum qtypes("GAUSS GRID", "GRID");
//   params.addParam<MooseEnum>("face_type", qtypes, "The face quadrature type");
//
//   MooseEnum convention("positive=0 negative=1", "positive");
//   params.addParam<MooseEnum>(
//       "internal_convention",
//       convention,
//       "The convention for spawning rays from internal sidesets; denotes the sign of the dot "
//       "product between a ray and the internal sideset side normal");
//
//   // Shouldn't ever need RayKernels for view factors
//   params.set<bool>("ray_kernel_coverage_check") = false;
//   params.suppressParameter<bool>("ray_kernel_coverage_check");
//
//   // So that the study executes before the RayTracingViewFactor
//   params.set<bool>("force_preaux") = true;
//   params.suppressParameter<bool>("force_preaux");
//
//   // Need to use internal sidesets
//   params.set<bool>("use_internal_sidesets") = true;
//   params.suppressParameter<bool>("use_internal_sidesets");
//
//   // Don't verify Rays in opt mode by default - it's expensive
//   params.set<bool>("verify_rays") = false;
//
//   // No need to use Ray registration
//   params.set<bool>("_use_ray_registration") = false;
//   // Do not need to bank Rays on completion
//   params.set<bool>("_bank_rays_on_completion") = false;
//
//   // Whether or not the given method has reciprocity
//   params.addPrivateParam<bool>("_has_reciprocity", false);
//
//   return params;
// }
//
// ViewFactorRayStudyBase::ViewFactorRayStudyBase(const InputParameters & parameters)
//   : RayTracingStudy(parameters),
//     _bnd_ids_vec(_mesh.getBoundaryIDs(getParam<std::vector<BoundaryName>>("boundary"))),
//     _bnd_ids(_bnd_ids_vec.begin(), _bnd_ids_vec.end()),
//     _internal_convention(getParam<MooseEnum>("internal_convention")),
//     _ray_index_start_bnd_id(registerRayAuxData("start_bnd_id")),
//     _ray_index_start_total_weight(registerRayAuxData("start_total_weight")),
//     _fe_face(FEBase::build(_mesh.dimension(), FEType(CONSTANT, MONOMIAL))),
//     _q_face(QBase::build(Moose::stringToEnum<QuadratureType>(getParam<MooseEnum>("face_type")),
//                          _mesh.dimension() - 1,
//                          Moose::stringToEnum<Order>(getParam<MooseEnum>("face_order")))),
//     _has_reciprocity(getParam<bool>("_has_reciprocity")),
//     _threaded_vf_info(libMesh::n_threads())
// {
//   _fe_face->attach_quadrature_rule(_q_face.get());
//   _fe_face->get_xyz();
// }
//
// void
// ViewFactorRayStudyBase::initialSetup()
// {
//   RayTracingStudy::initialSetup();
//
//   // We optimized away RayKernels, so don't allow them
//   if (hasRayKernels(/* tid = */ 0))
//     mooseError(_error_prefix, " is not compatible with RayKernels.");
// }
//
// void
// ViewFactorRayStudyBase::preExecuteStudy()
// {
//   // Clear and zero the view factor maps we're about to accumulate into for each thread
//   for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
//   {
//     _threaded_vf_info[tid].clear();
//     for (const BoundaryID from_id : _bnd_ids)
//       for (const BoundaryID to_id : _bnd_ids)
//         _threaded_vf_info[tid][from_id][to_id] = 0;
//   }
//
//   generateStartElems();
// }
//
// void
// ViewFactorRayStudyBase::postExecuteStudy()
// {
//   RayTracingStudy::postExecuteStudy();
//
//   // Finalize the cumulative _vf_info;
//   _vf_info.clear();
//   for (const BoundaryID from_id : _bnd_ids)
//     for (const BoundaryID to_id : _bnd_ids)
//       if (!_has_reciprocity || from_id <= to_id)
//       {
//         Real & entry = _vf_info[from_id][to_id];
//
//         // Zero before summing
//         entry = 0;
//
//         // Sum over threads
//         for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
//         {
//           entry += _threaded_vf_info[tid][from_id][to_id];
//           if (_has_reciprocity)
//             entry += _threaded_vf_info[tid][to_id][from_id];
//         }
//
//         // Sum over processors
//         _communicator.sum(entry);
//
//         if (_has_reciprocity && from_id != to_id)
//           _vf_info[to_id][from_id] = entry;
//       }
// }
//
// void
// ViewFactorRayStudyBase::addToViewFactorInfo(Real value,
//                                             const BoundaryID from_id,
//                                             const BoundaryID to_id,
//                                             const THREAD_ID tid)
// {
//   mooseAssert(currentlyPropagating(), "Can only be called during Ray tracing");
//   mooseAssert(_threaded_vf_info[tid].count(from_id),
//               "Threaded view factor info does not have from boundary");
//   mooseAssert(_threaded_vf_info[tid][from_id].count(to_id),
//               "Threaded view factor info does not have from -> to boundary");
//
//   _threaded_vf_info[tid][from_id][to_id] += value;
// }
//
// Real
// ViewFactorRayStudyBase::viewFactorInfo(const BoundaryID from_id, const BoundaryID to_id) const
// {
//   auto it = _vf_info.find(from_id);
//   if (it == _vf_info.end())
//     mooseError(_error_prefix, ": From boundary id ", from_id, " not in view factor map.");
//
//   auto itt = it->second.find(to_id);
//   if (itt == it->second.end())
//     mooseError(_error_prefix,
//                ": From boundary id ",
//                from_id,
//                " to boundary_id ",
//                to_id,
//                " not in view factor map.");
//   return itt->second;
// }
//
// void
// ViewFactorRayStudyBase::generateStartElems()
// {
//   const auto & points = _fe_face->get_xyz();
//   const auto & weights = _fe_face->get_JxW();
//
//   // Clear before filling
//   _start_elems.clear();
//
//   // Starting elements we have that are on the wrong side of an internal boundary
//   std::unordered_map<processor_id_type, std::vector<StartElem>> send_start_map;
//
//   // Get all possible points on the user defined boundaries on this proc
//   for (const BndElement * belem : *_mesh.getBoundaryElementRange())
//   {
//     const Elem * elem = belem->_elem;
//     const auto side = belem->_side;
//     const auto bnd_id = belem->_bnd_id;
//
//     // Skip if we don't own you
//     if (elem->processor_id() != _pid)
//       continue;
//
//     // Skip if the boundary id isn't one we're looking for
//     if (!_bnd_ids.count(bnd_id))
//       continue;
//
//     // Sanity check on QGRID not working on some types
//     if (_q_face->type() == QGRID && elem->type() == TET4)
//       mooseError(
//           "Cannot use GRID quadrature type with tetrahedral elements in ViewFactorRayStudy '",
//           _name,
//           "'");
//
//     // The elem/side that we will actually start the trace from
//     // (this may change on internal sidesets)
//     const Elem * start_elem = elem;
//     auto start_side = side;
//
//     // Reinit this face for points
//     _fe_face->reinit(elem, side);
//
//     // See if this boundary is internal
//     const Elem * neighbor = elem->neighbor_ptr(side);
//     if (neighbor)
//     {
//       if (!neighbor->active())
//         mooseError(type(), " does not work with adaptivity");
//
//       // With the positive convention, the Rays that we want to spawn from internal boundaries
//       // have positive dot products with the outward normal on the side. The ray-tracer requires
//       // that we provide an element incoming side that is actually incoming (the dot product with
//       // the direction and the normal is negative). Therefore, switch the physical trace to start
//       // from the other element and the corresponding side
//       if (_internal_convention == 0)
//       {
//         start_elem = neighbor;
//         start_side = neighbor->which_neighbor_am_i(elem);
//       }
//     }
//
//     // If we own the true starting elem, add to our start info. Otherwise, package the
//     // start info to be sent to the processor that will actually start this trace
//     const auto start_pid = start_elem->processor_id();
//     auto & add_to = _pid ? _start_elems : send_start_map[start_pid];
//     add_to.emplace_back(elem, start_elem, start_side, bnd_id, points, weights);
//
//     // Entry point for other methods that need to do other things with this info
//     postAddStartElem(elem, side, bnd_id, points, weights);
//   }
//
//   // If the internal convention is positive, we may have points that we switched to another
//   // element for the actual trace, so communicate those to the processors that will
//   // actually be starting them
//   if (_internal_convention == 0)
//   {
//     // Functor that takes in StartElems and appends them to our local list
//     auto append_start_elems = [this](processor_id_type,
//                                      const std::vector<StartElem> & start_elems) {
//       _start_elems.reserve(_start_elems.size() + start_elems.size());
//       for (const StartElem & start_elem : start_elems)
//         _start_elems.emplace_back(start_elem);
//     };
//
//     // Communicate and act on data
//     Parallel::push_parallel_packed_range(_communicator, send_start_map, this, append_start_elems);
//   }
// }
//
// void
// ViewFactorRayStudyBase::postAddStartElem(const Elem * /* elem */,
//                                          const unsigned short /* side */,
//                                          const BoundaryID /* bnd_id */,
//                                          const std::vector<Point> & /* points */,
//                                          const std::vector<Real> & /* weights */)
// {
// }
//
// namespace libMesh
// {
// namespace Parallel
// {
//
// unsigned int
// Packing<ViewFactorRayStudyBase::StartElem>::packing_size(const std::size_t num_points)
// {
//   // Number of points, elem_id, start_elem_id, incoming_side, bnd_id
//   unsigned int total_size = 5;
//   // Points
//   total_size += num_points * 3;
//   // Weights
//   total_size += num_points;
//
//   return total_size;
// }
//
// unsigned int
// Packing<ViewFactorRayStudyBase::StartElem>::packed_size(
//     typename std::vector<Real>::const_iterator in)
// {
//   const std::size_t num_points = *in++;
//   return packing_size(num_points);
// }
//
// unsigned int
// Packing<ViewFactorRayStudyBase::StartElem>::packable_size(
//     const ViewFactorRayStudyBase::StartElem & start_elem, const void *)
// {
//   mooseAssert(start_elem._points.size() == start_elem._weights.size(), "Size mismatch");
//   return packing_size(start_elem._points.size());
// }
//
// template <>
// ViewFactorRayStudyBase::StartElem
// Packing<ViewFactorRayStudyBase::StartElem>::unpack(std::vector<Real>::const_iterator in,
//                                                    ViewFactorRayStudyBase * study)
// {
//   // StartElem to fill into
//   ViewFactorRayStudyBase::StartElem start_elem;
//
//   // Number of points
//   const std::size_t num_points = static_cast<std::size_t>(*in++);
//
//   // Elem id
//   RayTracingPackingUtils::unpack(start_elem._elem, *in++, &study->meshBase());
//
//   // Start elem id
//   RayTracingPackingUtils::unpack(start_elem._start_elem, *in++, &study->meshBase());
//
//   // Incoming side
//   start_elem._incoming_side = static_cast<unsigned short>(*in++);
//
//   // Boundary ID
//   start_elem._bnd_id = static_cast<BoundaryID>(*in++);
//
//   // Points
//   start_elem._points.resize(num_points);
//   for (std::size_t i = 0; i < num_points; ++i)
//   {
//     start_elem._points[i](0) = *in++;
//     start_elem._points[i](1) = *in++;
//     start_elem._points[i](2) = *in++;
//   }
//
//   // Weights
//   start_elem._weights.resize(num_points);
//   for (std::size_t i = 0; i < num_points; ++i)
//     start_elem._weights[i] = *in++;
//
//   return start_elem;
// }
//
// template <>
// void
// Packing<ViewFactorRayStudyBase::StartElem>::pack(
//     const ViewFactorRayStudyBase::StartElem & start_elem,
//     std::back_insert_iterator<std::vector<Real>> data_out,
//     const ViewFactorRayStudyBase * study)
// {
//   // Number of points
//   data_out = static_cast<buffer_type>(start_elem._points.size());
//
//   // Elem id
//   data_out = RayTracingPackingUtils::pack<buffer_type>(start_elem._elem, &study->meshBase());
//
//   // Start elem id
//   data_out = RayTracingPackingUtils::pack<buffer_type>(start_elem._start_elem, &study->meshBase());
//
//   // Incoming side
//   data_out = static_cast<buffer_type>(start_elem._incoming_side);
//
//   // Boundary id
//   data_out = static_cast<buffer_type>(start_elem._bnd_id);
//
//   // Points
//   for (const auto & point : start_elem._points)
//   {
//     data_out = point(0);
//     data_out = point(1);
//     data_out = point(2);
//   }
//
//   // Weights
//   std::copy(start_elem._weights.begin(), start_elem._weights.end(), data_out);
// }
//
// } // namespace Parallel
//
// } // namespace libMesh
