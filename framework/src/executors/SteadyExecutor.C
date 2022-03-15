//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SteadyExecutor.h"

registerMooseObject("MooseApp", SteadyExecutor);

InputParameters
SteadyExecutor::executorParams()
{
  auto params = emptyInputParameters();
  params.addParam<Real>("time", 0, "System time");
  return params;
}

InputParameters
SteadyExecutor::validParams()
{
  InputParameters params = Executor::validParams();
  params += SteadyExecutor::executorParams();
  params.addRequiredParam<std::vector<ExecutorName>>("executors", "foo");
  return params;
}

SteadyExecutor::SteadyExecutor(const InputParameters & parameters)
  : Executor(parameters),
    _executors(getExecutors("executors")),
    _system_time(getParam<Real>("time")),
    _time_step(_fe_problem.timeStep()),
    _time(_fe_problem.time())
{
}

void
SteadyExecutor::init()
{
  if (_app.isRecovering())
  {
    _console << "\nCannot recover steady solves!\nExiting...\n" << std::endl;
    return;
  }

  if (_fe_problem.getNonlinearSystemBase().containsTimeKernel())
    mooseError("You have specified time kernels in your steady state simulation");

  _fe_problem.execute(EXEC_PRE_MULTIAPP_SETUP);
  _fe_problem.initialSetup();
}

Executor::Result
SteadyExecutor::run()
{
  Result & result = newResult();

  if (_app.isRecovering())
    return result;

  _time_step = 0;
  _time = _time_step;
  _fe_problem.outputStep(EXEC_INITIAL);
  _time = _system_time;

  _fe_problem.advanceState();

  // first step in any steady state solve is always 1 (preserving backwards compatibility)
  _time_step = 1;

  _fe_problem.timestepSetup();

  for (auto & executor : _executors)
  {
    Result executor_result = executor->execute();
    result.record(executor->name(), executor_result);
    if (!executor_result.convergedAll())
    {
      result.fail("Failed");
      break;
    }
  }

  // need to keep _time in sync with _time_step to get correct output
  _time = _time_step;
  _fe_problem.outputStep(EXEC_TIMESTEP_END);
  _time = _system_time;

  _fe_problem.execMultiApps(EXEC_FINAL);
  _fe_problem.finalizeMultiApps();
  _fe_problem.postExecute();
  _fe_problem.execute(EXEC_FINAL);
  _time = _time_step;
  _fe_problem.outputStep(EXEC_FINAL);
  _time = _system_time;

  return result;
}
