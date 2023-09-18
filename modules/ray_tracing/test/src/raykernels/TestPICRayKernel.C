// //* This file is part of the MOOSE framework
// //* https://www.mooseframework.org
// //*
// //* All rights reserved, see COPYRIGHT for full restrictions
// //* https://github.com/idaholab/moose/blob/master/COPYRIGHT
// //*
// //* Licensed under LGPL 2.1, please see LICENSE for details
// //* https://www.gnu.org/licenses/lgpl-2.1.html

// #include "TestPICRayKernel.h"

// #include "BasicPICRayStudy.h"

// registerMooseObject("RayTracingTestApp", TestPICRayKernel);

// InputParameters
// TestPICRayKernel::validParams()
// {
//   auto params = GeneralRayKernel::validParams();
//   return params;
// }

// TestPICRayKernel::TestPICRayKernel(const InputParameters & params)
//   : GeneralRayKernel(params), _ray_velocity_index(getStudy<BasicPICRayStudy>().rayVelocityIndex())
// {
// }

// void
// TestPICRayKernel::onSegment()
// {
//   // The distance that we should move during this trace
//   const auto max_distance =
//   // The current timestep
//   const auto dt = _fe_problem.dt();
//   // Writeable reference to the current ray's velocity
//   auto & velocity = currentRay()->data(_ray_velocity_index);
//   // The distance that we should move during the segment within this trace
//   const auto distance = velocity * dt;

// }
