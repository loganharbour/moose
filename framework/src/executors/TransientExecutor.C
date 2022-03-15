//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TransientExecutor.h"

#include "AuxiliarySystem.h"
#include "NonlinearSystem.h"
#include "TimeStepper.h"

registerMooseObject("MooseApp", TransientExecutor);

InputParameters
TransientExecutor::executorParams()
{
  auto params = emptyInputParameters();

  /**
   * For backwards compatibility we'll allow users to set the TimeIntegration scheme inside of
   the
   * executioner block
   * as long as the TimeIntegrator does not have any additional parameters.
   */
  MooseEnum schemes("implicit-euler explicit-euler crank-nicolson bdf2 explicit-midpoint dirk "
                    "explicit-tvd-rk-2 newmark-beta",
                    "implicit-euler");

  params.addParam<Real>("start_time", 0.0, "The start time of the simulation");
  params.addParam<Real>("end_time", 1.0e30, "The end time of the simulation");
  params.addParam<Real>("dt", 1., "The timestep size between solves");
  params.addParam<Real>("dtmin", 1.0e-13, "The minimum timestep size in an adaptive run");
  params.addParam<Real>("dtmax", 1.0e30, "The maximum timestep size in an adaptive run");
  params.addParam<bool>(
      "reset_dt", false, "Use when restarting a calculation to force a change in dt.");
  params.addParam<unsigned int>("num_steps",
                                std::numeric_limits<unsigned int>::max(),
                                "The number of timesteps in a transient run");
  params.addParam<int>("n_startup_steps", 0, "The number of timesteps during startup");

  params.addDeprecatedParam<bool>("trans_ss_check",
                                  false,
                                  "Whether or not to check for steady state conditions",
                                  "Use steady_state_detection instead");
  params.addDeprecatedParam<Real>("ss_check_tol",
                                  1.0e-08,
                                  "Whenever the relative residual changes by less "
                                  "than this the solution will be considered to be "
                                  "at steady state.",
                                  "Use steady_state_tolerance instead");
  params.addDeprecatedParam<Real>(
      "ss_tmin",
      0.0,
      "Minimum amount of time to run before checking for steady state conditions.",
      "Use steady_state_start_time instead");

  params.addParam<bool>(
      "steady_state_detection", false, "Whether or not to check for steady state conditions");
  params.addParam<Real>("steady_state_tolerance",
                        1.0e-08,
                        "Whenever the relative residual changes by less "
                        "than this the solution will be considered to be "
                        "at steady state.");
  params.addParam<Real>(
      "steady_state_start_time",
      0.0,
      "Minimum amount of time to run before checking for steady state conditions.");
  params.addParam<bool>(
      "normalize_solution_diff_norm_by_dt",
      true,
      "Whether to divide the solution difference norm by dt. If taking 'small' "
      "time steps you probably want this to be true. If taking very 'large' timesteps in an "
      "attempt to *reach* a steady-state, you probably want this parameter to be false.");

  params.addParam<std::vector<std::string>>("time_periods", "The names of periods");
  params.addParam<std::vector<Real>>("time_period_starts", "The start times of time periods");
  params.addParam<std::vector<Real>>("time_period_ends", "The end times of time periods");
  params.addParam<bool>(
      "abort_on_solve_fail", false, "abort if solve not converged rather than cut timestep");
  params.addParam<bool>(
      "error_on_dtmin",
      true,
      "Throw error when timestep is less than dtmin instead of just aborting solve.");
  params.addParam<MooseEnum>("scheme", schemes, "Time integration scheme used.");
  params.addParam<Real>("timestep_tolerance",
                        1.0e-13,
                        "the tolerance setting for final timestep size and sync times");

  params.addParam<bool>("use_multiapp_dt",
                        false,
                        "If true then the dt for the simulation will be "
                        "chosen by the MultiApps.  If false (the "
                        "default) then the minimum over the master dt "
                        "and the MultiApps is used");

  params.addParam<bool>("check_aux",
                        false,
                        "Whether to check the auxiliary system for convergence to steady-state. "
                        "If false, then the nonlinear system is used.");

  params.addParamNamesToGroup(
      "steady_state_detection steady_state_tolerance steady_state_start_time check_aux",
      "Steady State Detection");

  params.addParamNamesToGroup("start_time dtmin dtmax n_startup_steps trans_ss_check ss_check_tol "
                              "ss_tmin abort_on_solve_fail timestep_tolerance use_multiapp_dt",
                              "Advanced");

  params.addParamNamesToGroup("time_periods time_period_starts time_period_ends", "Time Periods");

  return params;
}

InputParameters
TransientExecutor::validParams()
{
  InputParameters params = Executor::validParams();
  params += TransientExecutor::executorParams();
  params.addRequiredParam<std::vector<ExecutorName>>("executors", "foo");
  return params;
}

TransientExecutor::TransientExecutor(const InputParameters & parameters)
  : Executor(parameters),
    _nl(_fe_problem.getNonlinearSystemBase()),
    _aux(_fe_problem.getAuxiliarySystem()),
    _check_aux(getParam<bool>("check_aux")),
    _time_scheme(getParam<MooseEnum>("scheme").getEnum<Moose::TimeIntegratorType>()),
    _t_step(_fe_problem.timeStep()),
    _time(_fe_problem.time()),
    _time_old(_fe_problem.timeOld()),
    _dt(_fe_problem.dt()),
    _dt_old(_fe_problem.dtOld()),
    _unconstrained_dt(declareRecoverableData<Real>("unconstrained_dt", -1)),
    _at_sync_point(declareRecoverableData<bool>("at_sync_point", false)),
    _last_solve_converged(declareRecoverableData<bool>("last_solve_converged", true)),
    _xfem_repeat_step(false),
    _end_time(getParam<Real>("end_time")),
    _dtmin(getParam<Real>("dtmin")),
    _dtmax(getParam<Real>("dtmax")),
    _num_steps(getParam<unsigned int>("num_steps")),
    _n_startup_steps(getParam<int>("n_startup_steps")),
    _steady_state_detection(getParam<bool>("steady_state_detection")),
    _steady_state_tolerance(getParam<Real>("steady_state_tolerance")),
    _steady_state_start_time(getParam<Real>("steady_state_start_time")),
    _sync_times(_app.getOutputWarehouse().getSyncTimes()),
    _abort(getParam<bool>("abort_on_solve_fail")),
    _error_on_dtmin(getParam<bool>("error_on_dtmin")),
    _time_interval(declareRecoverableData<bool>("time_interval", false)),
    _start_time(getParam<Real>("start_time")),
    _timestep_tolerance(getParam<Real>("timestep_tolerance")),
    _target_time(declareRecoverableData<Real>("target_time", -std::numeric_limits<Real>::max())),
    _use_multiapp_dt(getParam<bool>("use_multiapp_dt")),
    _solution_change_norm(declareRecoverableData<Real>("solution_change_norm", 0.0)),
    _sln_diff(_check_aux ? _aux.addVector("sln_diff", false, PARALLEL)
                         : _nl.addVector("sln_diff", false, PARALLEL)),
    _normalize_solution_diff_norm_by_dt(getParam<bool>("normalize_solution_diff_norm_by_dt")),
    _executors(getExecutors("executors"))
{
  // Handle deprecated parameters
  if (!parameters.isParamSetByAddParam("trans_ss_check"))
    _steady_state_detection = getParam<bool>("trans_ss_check");
  if (!parameters.isParamSetByAddParam("ss_check_tol"))
    _steady_state_tolerance = getParam<Real>("ss_check_tol");
  if (!parameters.isParamSetByAddParam("ss_tmin"))
    _steady_state_start_time = getParam<Real>("ss_tmin");

  _t_step = 0;
  _dt = 0;
  _next_interval_output_time = 0.0;

  // Either a start_time has been forced on us, or we want to tell the App about what our start time
  // is (in case anyone else is interested.
  if (_app.hasStartTime())
    _start_time = _app.getStartTime();
  else if (parameters.isParamSetByUser("start_time"))
    _app.setStartTime(_start_time);

  _time = _time_old = _start_time;
  _fe_problem.transient(true);

  setupTimeIntegrator();

  if (_app.halfTransient()) // Cut timesteps and end_time in half...
  {
    _end_time /= 2.0;
    _num_steps /= 2.0;

    if (_num_steps == 0) // Always do one step in the first half
      _num_steps = 1;
  }
}

void
TransientExecutor::init()
{
  if (!_time_stepper.get())
  {
    InputParameters pars = _app.getFactory().getValidParams("ConstantDT");
    pars.set<SubProblem *>("_subproblem") = &_fe_problem;
    // pars.set<Transient *>("_executioner") = this;

    /**
     * We have a default "dt" set in the Transient parameters but it's possible for users to set
     * other parameters explicitly that could provide a better calculated "dt". Rather than provide
     * difficult to understand behavior using the default "dt" in this case, we'll calculate "dt"
     * properly.
     */
    if (!_pars.isParamSetByAddParam("end_time") && !_pars.isParamSetByAddParam("num_steps") &&
        _pars.isParamSetByAddParam("dt"))
      pars.set<Real>("dt") = (getParam<Real>("end_time") - getParam<Real>("start_time")) /
                             static_cast<Real>(getParam<unsigned int>("num_steps"));
    else
      pars.set<Real>("dt") = getParam<Real>("dt");

    pars.set<bool>("reset_dt") = getParam<bool>("reset_dt");
    _time_stepper = _app.getFactory().create<TimeStepper>("ConstantDT", "TimeStepper", pars);
  }

  _fe_problem.execute(EXEC_PRE_MULTIAPP_SETUP);
  _fe_problem.initialSetup();

  /**
   * If this is a restart run, the user may want to override the start time, which we already set in
   * the constructor. "_time" however will have been "restored" from the restart file. We need to
   * honor the original request of the developer now that the restore has been completed. This must
   * occur before we init the time stepper (since that prints out the start time). The multiapp case
   * is also bit complicated. If we didn't set a start time, the app won't have it yet, so we just
   * restart the old time from the current time.
   */
  if (_app.isRestarting())
  {
    if (_app.hasStartTime())
      _time = _time_old = _app.getStartTime();
    else
      _time_old = _time;
  }

  _time_stepper->init();

  if (_app.isRecovering()) // Recover case
  {
    if (_t_step == 0)
      mooseError("Internal error in Transient executioner: _t_step is equal to 0 while recovering "
                 "in init().");

    _dt_old = _dt;
  }
}

Executor::Result
TransientExecutor::run()
{
  Result & result = newResult();

  // if (_app.isRecovering())
  //   return result;
  //
  // _time_step = 0;
  // _time = _time_step;
  // _fe_problem.outputStep(EXEC_INITIAL);
  // _time = _system_time;
  //
  // _fe_problem.advanceState();
  //
  // // first step in any steady state solve is always 1 (preserving backwards compatibility)
  // _time_step = 1;
  //
  // _fe_problem.timestepSetup();
  //
  // for (auto & executor : _executors)
  // {
  //   Result executor_result = executor->execute();
  //   result.record(executor->name(), executor_result);
  //   if (!executor_result.convergedAll())
  //   {
  //     result.fail("Failed");
  //     break;
  //   }
  // }
  //
  // // need to keep _time in sync with _time_step to get correct output
  // _time = _time_step;
  // _fe_problem.outputStep(EXEC_TIMESTEP_END);
  // _time = _system_time;
  //
  // _fe_problem.execMultiApps(EXEC_FINAL);
  // _fe_problem.finalizeMultiApps();
  // _fe_problem.postExecute();
  // _fe_problem.execute(EXEC_FINAL);
  // _time = _time_step;
  // _fe_problem.outputStep(EXEC_FINAL);
  // _time = _system_time;

  return result;
}

void
TransientExecutor::setupTimeIntegrator()
{
  if (_pars.isParamSetByUser("scheme") && _fe_problem.hasTimeIntegrator())
    mooseError("You cannot specify time_scheme in the Executioner and independently add a "
               "TimeIntegrator to the system at the same time");

  if (!_fe_problem.hasTimeIntegrator())
  {
    // backwards compatibility
    std::string ti_str;
    using namespace Moose;

    switch (_time_scheme)
    {
      case TI_IMPLICIT_EULER:
        ti_str = "ImplicitEuler";
        break;
      case TI_EXPLICIT_EULER:
        ti_str = "ExplicitEuler";
        break;
      case TI_CRANK_NICOLSON:
        ti_str = "CrankNicolson";
        break;
      case TI_BDF2:
        ti_str = "BDF2";
        break;
      case TI_EXPLICIT_MIDPOINT:
        ti_str = "ExplicitMidpoint";
        break;
      case TI_LSTABLE_DIRK2:
        ti_str = "LStableDirk2";
        break;
      case TI_EXPLICIT_TVD_RK_2:
        ti_str = "ExplicitTVDRK2";
        break;
      case TI_NEWMARK_BETA:
        ti_str = "NewmarkBeta";
        break;
      default:
        mooseError("Unknown scheme: ", _time_scheme);
        break;
    }

    InputParameters params = _app.getFactory().getValidParams(ti_str);
    _fe_problem.addTimeIntegrator(ti_str, ti_str, params);
  }
}
