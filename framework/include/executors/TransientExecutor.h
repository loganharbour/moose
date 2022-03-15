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

class TimeStepper;

class TransientExecutor : public Executor
{
public:
  TransientExecutor(const InputParameters & parameters);

  static InputParameters validParams();
  static InputParameters executorParams();

  void init() override;

  /**
   * Get the maximum dt
   * @return The maximum dt
   */
  Real & dtMax() { return _dtmax; }

  /**
   * Get the minimal dt
   * @return The minimal dt
   */
  Real & dtMin() { return _dtmin; }

  /**
   * Return the start time
   * @return The start time
   */
  Real getStartTime() { return _start_time; }

  /**
   * Get the end time
   * @return The end time
   */
  Real & endTime() { return _end_time; }

  /**
   * Get the timestep tolerance
   * @return The timestep tolerance
   */
  Real & timestepTol() { return _timestep_tolerance; }

  /**
   * Is the current step at a sync point (sync times, time interval, target time, etc)?
   * @return Bool indicataing whether we are at a sync point
   */
  bool atSyncPoint() const { return _at_sync_point; }

  /**
   * Get the unconstrained dt
   * @return Value of dt before constraints were applied
   */
  Real unconstrainedDT() const { return _unconstrained_dt; }

  /**
   * Set the number of time steps
   * @param num_steps number of time steps
   */
  void forceNumSteps(const unsigned int num_steps) { _num_steps = num_steps; }

protected:
  virtual Result run() override;

  /// The nonlinear system
  NonlinearSystemBase & _nl;

  /// The aux system
  AuxiliarySystem & _aux;

  /// Whether to use the auxiliary system solution to determine steady-states
  const bool _check_aux;

  const Moose::TimeIntegratorType _time_scheme;

  std::shared_ptr<TimeStepper> _time_stepper;

  /// The current timestep
  int & _t_step;
  /// The current time
  Real & _time;
  /// The previous time
  Real & _time_old;
  /// The current timestep size
  Real & _dt;
  /// The old timestep size
  Real & _dt_old;

  Real & _unconstrained_dt;
  bool & _at_sync_point;

  /// Whether or not the last solve converged
  bool & _last_solve_converged;

  /// Whether step should be repeated due to xfem modifying the mesh
  bool _xfem_repeat_step;

  /// The end time
  Real _end_time;
  /// The minimum timestep size
  Real _dtmin;
  /// The maximum timstep size
  Real _dtmax;
  /// The number of timesteps
  unsigned int _num_steps;
  /// The number of timesteps during startup
  int _n_startup_steps;

  /// Whether or not to use steady state detection
  bool _steady_state_detection;
  /// The steady state detection tolerance
  Real _steady_state_tolerance;
  /// The start time for steady state detection
  Real _steady_state_start_time;

  /// The times to sync
  std::set<Real> & _sync_times;

  /// Whether or not to abort on failure
  bool _abort;

  /// This parameter controls how the system will deal with _dt <= _dtmin
  /// If true, the time stepper is expected to throw an error
  /// If false, the executioner will continue through EXEC_FINAL
  const bool _error_on_dtmin;

  bool & _time_interval;
  Real _next_interval_output_time;
  Real _time_interval_output_interval;

  Real _start_time;
  Real _timestep_tolerance;
  Real & _target_time;
  bool _use_multiapp_dt;

  Real & _solution_change_norm;

  NumericVector<Number> & _sln_diff;

  /// Whether to divide the solution difference norm by dt. If taking 'small' time steps this member
  /// should probably be true. If taking very 'large' timesteps in an attempt to reach a
  /// steady-state, this member should probably be be false.
  const bool _normalize_solution_diff_norm_by_dt;

private:
  void setupTimeIntegrator();

  const std::vector<Executor *> _executors;
};
