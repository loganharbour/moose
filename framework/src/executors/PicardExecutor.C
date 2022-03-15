//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PicardExecutor.h"

#include "Executioner.h"
#include "FEProblemBase.h"
#include "NonlinearSystem.h"
#include "AllLocalDofIndicesThread.h"
#include "Console.h"

registerMooseObject("MooseApp", PicardExecutor);

InputParameters
PicardExecutor::executorParams()
{
  return FixedPointExecutor::executorParams();
}

InputParameters
PicardExecutor::validParams()
{
  return FixedPointExecutor::validParams();
}

PicardExecutor::PicardExecutor(const InputParameters & params) : FixedPointExecutor(params)
{
  allocateStorage(true);
}

void
PicardExecutor::allocateStorage(const bool primary)
{
  Real relaxation_factor;
  TagID old_tag_id;
  const std::vector<PostprocessorName> * transformed_pps;
  std::vector<std::vector<PostprocessorValue>> * transformed_pps_values;
  if (primary)
  {
    relaxation_factor = _relax_factor;
    old_tag_id = _problem.addVectorTag("xn_m1", Moose::VECTOR_TAG_SOLUTION);
    _old_tag_id = old_tag_id;
    transformed_pps = &_transformed_pps;
    transformed_pps_values = &_transformed_pps_values;
  }
  else
  {
    relaxation_factor = _secondary_relaxation_factor;
    old_tag_id = _problem.addVectorTag("secondary_xn_m1", Moose::VECTOR_TAG_SOLUTION);
    _secondary_old_tag_id = old_tag_id;
    transformed_pps = &_secondary_transformed_pps;
    transformed_pps_values = &_secondary_transformed_pps_values;
  }

  if (relaxation_factor != 1.)
  {
    // Store a copy of the previous solution
    _nl.addVector(old_tag_id, false, PARALLEL);

    // Allocate storage for the previous postprocessor values
    (*transformed_pps_values).resize((*transformed_pps).size());
    for (size_t i = 0; i < (*transformed_pps).size(); i++)
      (*transformed_pps_values)[i].resize(1);
  }
}

void
PicardExecutor::saveVariableValues(const bool primary)
{
  Real relaxation_factor;
  TagID old_tag_id;
  if (primary)
  {
    relaxation_factor = _relax_factor;
    old_tag_id = _old_tag_id;
  }
  else
  {
    relaxation_factor = _secondary_relaxation_factor;
    old_tag_id = _secondary_old_tag_id;
  }

  if (relaxation_factor != 1.)
  {
    // Save variable previous values
    NumericVector<Number> & solution = _nl.solution();
    NumericVector<Number> & transformed_old = _nl.getVector(old_tag_id);
    transformed_old = solution;
  }
}

void
PicardExecutor::savePostprocessorValues(const bool primary)
{
  Real relaxation_factor;
  const std::vector<PostprocessorName> * transformed_pps;
  std::vector<std::vector<PostprocessorValue>> * transformed_pps_values;
  if (primary)
  {
    relaxation_factor = _relax_factor;
    transformed_pps = &_transformed_pps;
    transformed_pps_values = &_transformed_pps_values;
  }
  else
  {
    relaxation_factor = _secondary_relaxation_factor;
    transformed_pps = &_secondary_transformed_pps;
    transformed_pps_values = &_secondary_transformed_pps_values;
  }

  if (relaxation_factor != 1.)
    // Save postprocessor previous values
    for (size_t i = 0; i < (*transformed_pps).size(); i++)
      (*transformed_pps_values)[i][0] = getPostprocessorValueByName((*transformed_pps)[i]);
}

bool
PicardExecutor::useFixedPointAlgorithmUpdateInsteadOfPicard(const bool primary)
{
  // unrelaxed Picard is the default update for fixed point iterations
  // old values are required for relaxation
  if (primary)
    return _relax_factor != 1. && _fixed_point_it > 0;
  else
    return _secondary_relaxation_factor != 1. && _main_fixed_point_it > 0;
}

void
PicardExecutor::transformPostprocessors(const bool primary)
{
  Real relaxation_factor;
  const std::vector<PostprocessorName> * transformed_pps;
  std::vector<std::vector<PostprocessorValue>> * transformed_pps_values;
  if (primary)
  {
    relaxation_factor = _relax_factor;
    transformed_pps = &_transformed_pps;
    transformed_pps_values = &_transformed_pps_values;
  }
  else
  {
    relaxation_factor = _secondary_relaxation_factor;
    transformed_pps = &_secondary_transformed_pps;
    transformed_pps_values = &_secondary_transformed_pps_values;
  }

  // Relax the postprocessors
  for (size_t i = 0; i < (*transformed_pps).size(); i++)
  {
    // Get new postprocessor value
    const Real current_value = getPostprocessorValueByName((*transformed_pps)[i]);
    const Real old_value = (*transformed_pps_values)[i][0];

    // Compute and set relaxed value
    Real new_value = current_value;
    new_value = relaxation_factor * current_value + (1 - relaxation_factor) * old_value;
    _problem.setPostprocessorValueByName((*transformed_pps)[i], new_value);
  }
}

void
PicardExecutor::transformVariables(const std::set<dof_id_type> & transformed_dofs,
                                   const bool primary)
{
  Real relaxation_factor;
  TagID old_tag_id;
  if (primary)
  {
    relaxation_factor = _relax_factor;
    old_tag_id = _old_tag_id;
  }
  else
  {
    relaxation_factor = _secondary_relaxation_factor;
    old_tag_id = _secondary_old_tag_id;
  }

  NumericVector<Number> & solution = _nl.solution();
  NumericVector<Number> & transformed_old = _nl.getVector(old_tag_id);

  for (const auto & dof : transformed_dofs)
    solution.set(dof,
                 (transformed_old(dof) * (1.0 - relaxation_factor)) +
                     (solution(dof) * relaxation_factor));

  solution.close();
  _nl.update();
}

void
PicardExecutor::printFixedPointConvergenceHistory()
{
  _console << "\n 0 Picard |R| = "
           << Console::outputNorm(std::numeric_limits<Real>::max(), _fixed_point_initial_norm)
           << '\n';

  for (unsigned int i = 0; i <= _fixed_point_it; ++i)
  {
    Real max_norm =
        std::max(_fixed_point_timestep_begin_norm[i], _fixed_point_timestep_end_norm[i]);
    _console << std::setw(2) << i + 1
             << " Picard |R| = " << Console::outputNorm(_fixed_point_initial_norm, max_norm)
             << '\n';
  }

  _console << std::flush;
}
