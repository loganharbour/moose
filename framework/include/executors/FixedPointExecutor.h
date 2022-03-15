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
class MooseMesh;
class NonlinearSystemBase;
class AuxiliarySystem;

class FixedPointExecutor : public Executor
{
public:
  FixedPointExecutor(const InputParameters & parameters);

  static InputParameters executorParams();
  static InputParameters validParams();

  /// Enumeration for fixed point convergence reasons
  enum class MooseFixedPointConvergenceReason
  {
    UNSOLVED = 0,
    CONVERGED_NONLINEAR = 1,
    CONVERGED_ABS = 2,
    CONVERGED_RELATIVE = 3,
    CONVERGED_CUSTOM = 4,
    REACH_MAX_ITS = 5,
    DIVERGED_MAX_ITS = -1,
    DIVERGED_NONLINEAR = -2,
    DIVERGED_FAILED_MULTIAPP = -3
  };

protected:
  Result run() override final;

  /// Reference to FEProblem
  FEProblemBase & _problem;
  /// Displaced problem
  std::shared_ptr<DisplacedProblem> _displaced_problem;
  /// Mesh
  MooseMesh & _mesh;
  /// Displaced mesh
  MooseMesh * _displaced_mesh;
  /// Reference to nonlinear system base for faster access
  NonlinearSystemBase & _nl;
  /// Reference to auxiliary system for faster access
  AuxiliarySystem & _aux;

  /// Minimum fixed point iterations
  const unsigned int _min_fixed_point_its;
  /// Maximum fixed point iterations
  const unsigned int _max_fixed_point_its;
  /// Whether or not we activate fixed point iteration
  const bool _has_fixed_point_its;
  /// Whether or not to treat reaching maximum number of fixed point iteration as converged
  const bool _accept_max_it;
  /// Whether or not to use residual norm to check the fixed point convergence
  const bool _has_fixed_point_norm;
  /// Relative tolerance on residual norm
  const Real _fixed_point_rel_tol;
  /// Absolute tolerance on residual norm
  const Real _fixed_point_abs_tol;
  /// Whether or not we force evaluation of residual norms even without multiapps
  const bool _fixed_point_force_norms;
  /// Postprocessor value for user-defined fixed point convergence check
  const PostprocessorValue * const _fixed_point_custom_pp;
  /// Relaxation factor for fixed point Iteration
  const Real _relax_factor;
  /// Relative tolerance on postprocessor value
  const Real _custom_rel_tol;
  /// Absolute tolerance on postprocessor value
  const Real _custom_abs_tol;

  /// The variables (transferred or not) that are going to be relaxed
  const std::vector<std::string> & _transformed_vars;
  /// The postprocessors (transferred or not) that are going to be relaxed
  const std::vector<PostprocessorName> & _transformed_pps;
  /// Previous values of the relaxed postprocessors
  std::vector<std::vector<PostprocessorValue>> _transformed_pps_values;

  /// Relaxation factor outside of fixed point iteration (used as a subapp)
  Real _secondary_relaxation_factor;
  /// Variables to be relaxed outside of fixed point iteration (used as a subapp)
  std::vector<std::string> _secondary_transformed_variables;
  /// Postprocessors to be relaxed outside of fixed point iteration (used as a subapp)
  std::vector<PostprocessorName> _secondary_transformed_pps;
  /// Previous values of the postprocessors relaxed outside of the fixed point iteration (used as a subapp)
  std::vector<std::vector<PostprocessorValue>> _secondary_transformed_pps_values;

  ///@{ Variables used by the fixed point iteration
  /// fixed point iteration counter
  unsigned int _fixed_point_it;
  /// fixed point iteration counter for the main app
  unsigned int _main_fixed_point_it;
  /// Initial residual norm
  Real _fixed_point_initial_norm;
  /// Full history of residual norm after evaluation of timestep_begin
  std::vector<Real> _fixed_point_timestep_begin_norm;
  /// Full history of residual norm after evaluation of timestep_end
  std::vector<Real> _fixed_point_timestep_end_norm;
  /// Status of fixed point solve
  MooseFixedPointConvergenceReason _fixed_point_status;
  ///@}

  /// Save both the variable and postprocessor values
  virtual void saveAllValues(const bool primary);

  /**
   * Saves the current values of the variables, and update the old(er) vectors.
   *
   * @param primary Whether this routine is to save the variables for the primary transformed
   *                quantities (as main app) or the secondary ones (as a subapp)
   */
  virtual void saveVariableValues(const bool primary) = 0;

  /**
   * Saves the current values of the postprocessors, and update the old(er) vectors.
   *
   * @param primary Whether this routine is to save the variables for the primary transformed
   *                quantities (as main app) or the secondary ones (as a subapp)
   */
  virtual void savePostprocessorValues(const bool primary) = 0;

  /**
   * Use the fixed point algorithm transform instead of simply using the Picard update
   * This routine can be used to alternate Picard iterations and fixed point algorithm
   * updates based on the values of the variables before and after a solve / a Picard iteration.
   *
   * @param primary Whether this routine is used for the primary transformed
   *                quantities (as main app) or the secondary ones (as a subapp)
   */
  virtual bool useFixedPointAlgorithmUpdateInsteadOfPicard(const bool primary) = 0;

  /**
   * Use the fixed point algorithm to transform the postprocessors.
   * If this routine is not called, the next value of the postprocessors will just be from
   * the unrelaxed Picard fixed point algorithm.
   *
   * @param primary Whether this routine is to save the variables for the primary transformed
   *                quantities (as main app) or the secondary ones (as a subapp)
   */
  virtual void transformPostprocessors(const bool primary) = 0;

  /**
   * Use the fixed point algorithm to transform the variables.
   * If this routine is not called, the next value of the variables will just be from
   * the unrelaxed Picard fixed point algorithm.
   *
   * @param transformed_dofs The dofs that will be affected by the algorithm
   * @param primary Whether this routine is to save the variables for the primary transformed
   *                quantities (as main app) or the secondary ones (as a subapp)
   */
  virtual void transformVariables(const std::set<dof_id_type> & transformed_dofs,
                                  const bool primary) = 0;

  /// Print the convergence history of the coupling, at every fixed point iteration
  virtual void printFixedPointConvergenceHistory() = 0;

  virtual bool
  solveStep(Real & begin_norm, Real & end_norm, const std::set<dof_id_type> & transformed_dofs);

private:
  /// Computes and prints the user-specified postprocessor assessing convergence
  void computeCustomConvergencePostprocessor();

  /// Examine the various convergence metrics
  bool examineFixedPointConvergence(bool & converged);

  /// Print information about the fixed point convergence
  void printFixedPointConvergenceReason();

  /// Whether sub-applications are automatically advanced no matter what happens during their solves
  bool autoAdvance() const;

  Executor & _executor;

  /// Old value of the custom convergence check postprocessor
  Real _pp_old;
  /// Current value of the custom convergence check postprocessor
  Real _pp_new;
  /// Scaling of custom convergence check postprocessor (its initial value)
  Real _pp_scaling;
  /// Convergence history of the custom convergence check postprocessor
  std::ostringstream _pp_history;

  /// Maximum number of xfem updates per step
  const unsigned int _max_xfem_update;
  /// Controls whether xfem should update the mesh at the beginning of the time step
  const bool _update_xfem_at_timestep_begin;
  /// Counter for number of xfem updates that have been performed in the current step
  unsigned int _xfem_update_count;
  /// Whether step should be repeated due to xfem modifying the mesh
  bool _xfem_repeat_step;

  /// Time of previous fixed point solve as a subapp
  Real _old_entering_time;

  /// Console output for whether the solve is skipped or not
  const std::string _solve_message;

  /// force the current step to fail, triggering are repeat with a cut dt
  bool _fail_step;

  /// Whether the user has set the auto_advance parameter for handling advancement of
  /// sub-applications in multi-app contexts
  const bool _auto_advance_set_by_user;

  /// The value of auto_advance set by the user for handling advancement of sub-applications in
  /// multi-app contexts
  const bool _auto_advance_user_value;
};
