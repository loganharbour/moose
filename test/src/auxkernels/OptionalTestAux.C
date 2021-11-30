//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "OptionalTestAux.h"

registerMooseObject("MooseTestApp", OptionalTestAux);

InputParameters
OptionalTestAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<MaterialPropertyName>("prop1", "first optional property");
  params.addRequiredParam<MaterialPropertyName>("prop2", "second optional property");
  params.addRequiredParam<bool>("expect1", "expect first property to exist");
  params.addRequiredParam<bool>("expect2", "expect second property to exist");
  return params;
}

OptionalTestAux::OptionalTestAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _prop1(getOptionalMaterialProperty<Real>("prop1")),
    _prop2(getOptionalMaterialProperty<Real>("prop2")),
    _expect1(getParam<bool>("expect1")),
    _expect2(getParam<bool>("expect2"))
{
}

Real
OptionalTestAux::computeValue()
{
  if (_expect1 == !_prop1)
    mooseError("Unexpected existence of prop1");
  if (_expect2 == !_prop2)
    mooseError("Unexpected existence of prop2");

  return (_prop1 ? (*_prop1)[_qp] : 0.0) + (_prop2 ? (*_prop2)[_qp] : 0.0);
}
