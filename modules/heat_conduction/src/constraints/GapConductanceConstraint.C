//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GapConductanceConstraint.h"

registerADMooseObject("HeatConductionApp", GapConductanceConstraint);

defineADValidParams(
    GapConductanceConstraint,
    ADMortarConstraint,
    params.addClassDescription(
        "Computes the residual and Jacobian contributions for the 'Lagrange Multiplier' "
        "implementation of the thermal contact problem. For more information, see the "
        "detailed description here: http://tinyurl.com/gmmhbe9");

    params.addRequiredParam<Real>("k", "Gap conductance"););

template <ComputeStage compute_stage>
GapConductanceConstraint<compute_stage>::GapConductanceConstraint(
    const InputParameters & parameters)
  : ADMortarConstraint<compute_stage>(parameters), _k(getParam<Real>("k"))
{
}

template <ComputeStage compute_stage>
ADReal
GapConductanceConstraint<compute_stage>::computeQpResidual(Moose::MortarType mortar_type)
{
  switch (mortar_type)
  {
    case Moose::MortarType::Master:
      return _lambda[_qp] * _test_master[_i][_qp];
    case Moose::MortarType::Slave:
      return -_lambda[_qp] * _test_slave[_i][_qp];
    case Moose::MortarType::Lower:
    {
      auto l = (_phys_points_master[_qp] - _phys_points_slave[_qp]).norm();
      return (_k * (_u_master[_qp] - _u_slave[_qp]) / l - _lambda[_qp]) * _test[_i][_qp];
    }
    default:
      return 0;
  }
}
