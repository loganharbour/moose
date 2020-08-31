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
class AQRayTracingViewFactor;
class AQViewFactorRayStudy;

template <>
InputParameters validParams<AQRayTracingViewFactor>();

/**
 * Computes the view factors for planar faces in unobstructed radiative heat transfer
 */
class AQRayTracingViewFactor : public ViewFactorBase
{
public:
  static InputParameters validParams();

  AQRayTracingViewFactor(const InputParameters & parameters);

  virtual void execute() override;
  virtual void initialize() override;

protected:
  virtual void threadJoinViewFactor(const UserObject & y) override;
  virtual void finalizeViewFactor() override;

  const AQViewFactorRayStudy & _ray_study;
  Real _divisor;
};
