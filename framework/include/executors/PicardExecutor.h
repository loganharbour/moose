//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FixedPointExecutor.h"

class PicardExecutor : public FixedPointExecutor
{
public:
  PicardExecutor(const InputParameters & params);

  static InputParameters executorParams();
  static InputParameters validParams();

  /**
   * Allocate storage for the fixed point algorithm.
   * This creates the system vector of old (older, pre/post solve) variable values and the
   * array of old (older, pre/post solve) postprocessor values.
   *
   * @param primary Whether this routine is to allocate storage for the primary transformed
   *                quantities (as main app) or the secondary ones (as a subapp)
   */
  void allocateStorage(const bool primary);

private:
  /**
   * Saves the current values of the variables, and update the old(er) vectors.
   *
   * @param primary Whether this routine is to save the variables for the primary transformed
   *                quantities (as main app) or the secondary ones (as a subapp)
   */
  void saveVariableValues(const bool primary) override final;

  /**
   * Saves the current values of the postprocessors, and update the old(er) vectors.
   *
   * @param primary Whether this routine is to save the variables for the primary transformed
   *                quantities (as main app) or the secondary ones (as a subapp)
   */
  void savePostprocessorValues(const bool primary) override final;

  /**
   * Use the fixed point algorithm transform instead of simply using the Picard update
   *
   * @param primary Whether this routine is used for the primary transformed
   *                quantities (as main app) or the secondary ones (as a subapp)
   */
  bool useFixedPointAlgorithmUpdateInsteadOfPicard(const bool primary) override final;

  /**
   * Use the fixed point algorithm to transform the postprocessors.
   * If this routine is not called, the next value of the postprocessors will just be from
   * the unrelaxed Picard fixed point algorithm.
   *
   * @param primary Whether this routine is to save the variables for the primary transformed
   *                quantities (as main app) or the secondary ones (as a subapp)
   */
  void transformPostprocessors(const bool primary) override final;

  /**
   * Use the fixed point algorithm to transform the variables.
   * If this routine is not called, the next value of the variables will just be from
   * the unrelaxed Picard fixed point algorithm.
   *
   * @param transformed_dofs The dofs that will be affected by the algorithm
   * @param primary Whether this routine is to save the variables for the primary transformed
   *                quantities (as main app) or the secondary ones (as a subapp)
   */
  void transformVariables(const std::set<dof_id_type> & transformed_dofs,
                          const bool primary) override final;

  /// Print the convergence history of the coupling, at every coupling iteration
  void printFixedPointConvergenceHistory() override final;

  /// Vector tag id for the previous solution variable, as a main app
  TagID _old_tag_id;

  /// Vector tag id for the previous solution variable, as a sub app
  TagID _secondary_old_tag_id;
};
