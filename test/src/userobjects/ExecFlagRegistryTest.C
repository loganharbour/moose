//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ExecFlagRegistryTest.h"

registerMooseObject("MooseTestApp", ExecFlagRegistryTest);

InputParameters
ExecFlagRegistryTest::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  params.addParam<bool>("add_all", false, "True to test adding the ALL execute flag");
  params.addParam<bool>("add_undefined", false, "True to test adding an undefined execute flag");

  return params;
}

ExecFlagRegistryTest::ExecFlagRegistryTest(const InputParameters & params)
  : GeneralUserObject(params)
{
  if (getParam<bool>("add_all"))
    _app.addExecFlag(EXEC_ALL);
  if (getParam<bool>("add_undefined"))
    _app.addExecFlag(ExecFlagType("FOO"));
}
