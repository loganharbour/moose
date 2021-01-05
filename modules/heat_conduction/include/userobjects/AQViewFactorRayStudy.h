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

class RayTracingAngularQuadrature;

/**
 * RayTracingStudy used to generate Rays for view factor computation
 * using the angular quadrature method.
 */
class AQViewFactorRayStudy : public ViewFactorRayStudyBase
{
public:
  AQViewFactorRayStudy(const InputParameters & parameters);

  static InputParameters validParams();

  void initialSetup() override;

protected:
  void generateRays() override;

private:
  /// if problem has 3 dimensions
  bool _is3D;

  /// angular quadrature object
  std::unique_ptr<RayTracingAngularQuadrature> _aq;

  ///@{ angular quadrature info
  std::vector<Real> _aq_angles;
  std::vector<Real> _aq_weights;
  ///@}

  /// the number of angular directions
  unsigned int _ndir;
};
