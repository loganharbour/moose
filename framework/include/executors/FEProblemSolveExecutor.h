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

class FEProblemBase;

class FEProblemSolveExecutor : public Executor
{
public:
  FEProblemSolveExecutor(const InputParameters & parameters);

  static InputParameters validParams();
  static InputParameters executorParams();

  static const std::set<std::string> & mooseLineSearches() { return _moose_line_searches; }

protected:
  virtual Result run() override;

  /// Reference to FEProblem
  FEProblemBase & _fe_problem;

  /// The number of steps to perform in a grid sequencing algorithm. This is one
  /// less than the number of grids requested by the user in their input,
  /// e.g. if they requested num_grids = 1, then there won't be any steps
  /// because we only need to perform one solve per time-step. Storing this
  /// member in this way allows for easy boolean operations, e.g. if (_num_grid_steps)
  /// as opposed to if (_num_grids)
  const unsigned int _num_grid_steps;

  /// Moose provided line searches
  static std::set<std::string> const _moose_line_searches;
};
