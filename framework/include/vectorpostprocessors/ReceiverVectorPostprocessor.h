//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralVectorPostprocessor.h"

class ReceiverVectorPostprocessor : public GeneralVectorPostprocessor
{
public:
  static InputParameters validParams();

  ReceiverVectorPostprocessor(const InputParameters & parameters);

  void initialize() override final {}
  void execute() override final {}
};
