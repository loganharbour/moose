//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeshChangedInterface.h"

#include "MooseApp.h"

defineLegacyParams(MeshChangedInterface);

InputParameters
MeshChangedInterface::validParams()
{
  InputParameters params = emptyInputParameters();
  return params;
}

MeshChangedInterface::MeshChangedInterface(const InputParameters & params)
{
  params.getCheckedPointerParam<MooseApp *>("_moose_app")->registerInterfaceObject(*this);
}
