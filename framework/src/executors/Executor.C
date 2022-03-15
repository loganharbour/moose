//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Executor.h"
#include "MooseApp.h"
#include "FEProblem.h"

#include "ExecutorExecutioner.h"

InputParameters
Executor::validParams()
{
  InputParameters params = MooseObject::validParams();
  params += Reporter::validParams();
  params += ReporterInterface::validParams();

  params.addParam<ExecFlagType>(
      "begin_exec_flag", EXEC_NONE, "exec flag associated with the beginning of this executor");
  params.addParam<ExecFlagType>(
      "end_exec_flag", EXEC_NONE, "exec flag associated with the end of this executor");

  params.addParam<bool>("verbose", false, "Set to true to print additional information");

  params.registerBase("Executor");
  return params;
}

Executor::Executor(const InputParameters & parameters)
  : MooseObject(parameters),
    Reporter(this),
    ReporterInterface(this),
    ExecutorInterface(static_cast<MooseObject &>(*this)),
    UserObjectInterface(this),
    PostprocessorInterface(this),
    Restartable(this, "Executors"),
    PerfGraphInterface(this),
    _fe_problem(*getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")),
    _begin_flag(getParam<ExecFlagType>("begin_exec_flag")),
    _end_flag(getParam<ExecFlagType>("end_exec_flag")),
    _verbose(getParam<bool>("verbose")),
    _dummy_executioner(std::make_shared<ExecutorExecutioner>(*this))
{
  // if (!parameters.isParamSetByUser("begin_exec_flag"))
  // {
  //   _begin_flag = ExecFlagType("exec_" + _name + "_begin");
  //   _app.addExecFlag(_begin_flag);
  // }
  // if (!parameters.isParamSetByUser("end_exec_flag"))
  // {
  //   _end_flag = ExecFlagType("exec_" + _name + "_end");
  //   _app.addExecFlag(_end_flag);
  // }
}

Executor::Result
Executor::execute()
{
  TIME_SECTION("execute", 3);

  _fe_problem.executeAllObjects(_begin_flag);
  auto result = run();
  _fe_problem.executeAllObjects(_end_flag);

  _result = result;
  return result;
}

Executor &
Executor::addSubExecutor(const std::string & type,
                         const std::string & name,
                         InputParameters params,
                         const bool apply_params /* = false */)
{
  if (apply_params)
    params.applyParameters(_pars);
  return *_app.addExecutor(type, this->name() + "_" + name, params);
}
