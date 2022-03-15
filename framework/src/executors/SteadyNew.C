//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Steady.h"

#include "SteadyExecutor.h"
#include "FEProblemSolveExecutor.h"

registerMooseObject("MooseApp", Steady);

InputParameters
Steady::validParams()
{
  InputParameters params = Executor::validParams();

  params += SteadyExecutor::executorParams();
  params += FEProblemSolveExecutor::executorParams();

  return params;
}

Steady::Steady(const InputParameters & parameters)
  : Executor(parameters), _steady_executor(buildExecutors())
{
}

Executor &
Steady::buildExecutors()
{
  const auto solve_params = _app.getFactory().getValidParams("FEProblemSolveExecutor");
  const auto & solve_executor =
      addSubExecutor("FEProblemSolveExecutor", "solve", solve_params, /* apply_params = */ true);

  auto steady_params = _app.getFactory().getValidParams("SteadyExecutor");
  steady_params.set<std::vector<ExecutorName>>("executors") = {solve_executor.name()};
  return addSubExecutor("SteadyExecutor", "steady", steady_params, /* apply_params = */ true);
}

Executor::Result
Steady::run()
{
  Result & result = newResult();
  result.record(_steady_executor.name(), _steady_executor.execute());
  return result;
}

void
Steady::init()
{
  _steady_executor.init();
}
