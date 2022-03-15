//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ExecutorInterface.h"

#include "FEProblem.h"
#include "MooseObject.h"
#include "NullExecutor.h"
#include "Executor.h"
#include "MooseApp.h"

InputParameters
ExecutorInterface::validParams()
{
  return emptyInputParameters();
}

ExecutorInterface::ExecutorInterface(const MooseObject & moose_object)
  : _ei_moose_object(moose_object), _ei_app(moose_object.getMooseApp())
{
}

Executor &
ExecutorInterface::getExecutor(const std::string & param_name) const
{
  const auto & params = _ei_moose_object.parameters();
  if (!params.isParamValid(param_name))
    _ei_moose_object.mooseError("Failed to get a parameter with the name \"",
                                param_name,
                                "\" when getting an ExecutorName.",
                                "\n\nKnown parameters:\n",
                                _ei_moose_object.parameters());
  if (!params.isType<ExecutorName>(param_name))
    _ei_moose_object.paramError(param_name,
                                "Parameter of type \"",
                                params.type(param_name),
                                "\" is not an expected type for getting the name of an Executor.");

  const auto & executor_name = _ei_moose_object.getParam<ExecutorName>(param_name);
  return getExecutorByName(executor_name);
}

Executor &
ExecutorInterface::getExecutorByName(const ExecutorName & executor_name) const
{
  auto & executor = _ei_app.getExecutor(executor_name);
  _coupled_executors.insert(&executor);
  return executor;
}

std::vector<Executor *>
ExecutorInterface::getExecutors(const std::string & param_name) const
{
  const auto & params = _ei_moose_object.parameters();
  if (!params.isParamValid(param_name))
    _ei_moose_object.mooseError("Failed to get a parameter with the name \"",
                                param_name,
                                "\" when getting a vector of ExecutorName.",
                                "\n\nKnown parameters:\n",
                                _ei_moose_object.parameters());
  if (!params.isType<std::vector<ExecutorName>>(param_name))
    _ei_moose_object.paramError(param_name,
                                "Parameter of type \"",
                                params.type(param_name),
                                "\" is not an expected type for getting Executor name(s).");

  const auto & executor_names = _ei_moose_object.getParam<std::vector<ExecutorName>>(param_name);
  return getExecutorsByName(executor_names);
}

std::vector<Executor *>
ExecutorInterface::getExecutorsByName(const std::vector<ExecutorName> & executor_names) const
{
  std::vector<Executor *> executors;
  executors.reserve(executor_names.size());
  for (const auto & executor_name : executor_names)
    executors.push_back(&getExecutorByName(executor_name));
  return executors;
}
