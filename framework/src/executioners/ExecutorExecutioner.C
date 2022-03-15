//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ExecutorExecutioner.h"
#include "Executor.h"

registerMooseObject("MooseApp", ExecutorExecutioner);

InputParameters
ExecutorExecutioner::validParams()
{
  return Executioner::validParams();
}

ExecutorExecutioner::ExecutorExecutioner(const Executor & executor)
  : Executioner(executor.parameters(), false), _executor(executor)
{
}

void
ExecutorExecutioner::execute()
{
}

bool
ExecutorExecutioner::lastSolveConverged() const
{
  return _executor.lastSolveConverged();
}
