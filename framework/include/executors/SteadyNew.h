//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Executor.h"

class SteadyExecutor;

class Steady : public Executor
{
public:
  static InputParameters validParams();

  Steady(const InputParameters & parameters);

  virtual void init() override;

protected:
  virtual Result run() override;

private:
  Executor & buildExecutors();

  Executor & _steady_executor;
};
