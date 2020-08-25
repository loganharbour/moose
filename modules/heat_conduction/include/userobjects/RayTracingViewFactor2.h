//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ViewFactorBase.h"

// Forward Declarations
class RayTracingViewFactor2;
class ViewFactorRayStudy2;

template <>
InputParameters validParams<RayTracingViewFactor2>();

/**
 * Computes the view factors for planar faces in unobstructed radiative heat transfer
 */
class RayTracingViewFactor2 : public ViewFactorBase
{
public:
  static InputParameters validParams();

  RayTracingViewFactor2(const InputParameters & parameters);

  virtual void execute() override;
  virtual void initialize() override;

protected:
  virtual void threadJoinViewFactor(const UserObject & y) override;
  virtual void finalizeViewFactor() override;

  const ViewFactorRayStudy2 & _ray_study;
  Real _divisor;
};
