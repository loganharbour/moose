//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ReporterTransferInterface.h"

#include "UserObject.h"
#include "Reporter.h"
#include "Transfer.h"

InputParameters
ReporterTransferInterface::validParams()
{
  InputParameters params = emptyInputParameters();
  return params;
}

ReporterTransferInterface::ReporterTransferInterface(const Transfer * transfer)
  : _rti_transfer(*transfer)
{
}

void
ReporterTransferInterface::addReporterTransferMode(const ReporterName & name,
                                                   const ReporterMode & mode,
                                                   FEProblemBase & problem)
{
  check(name, problem, "addReporterTransferMode()");

  problem.getReporterData(ReporterData::WriteKey())
      .getReporterStateBase(name)
      .addConsumer(mode, _rti_transfer);
}

void
ReporterTransferInterface::transferReporter(const ReporterName & from_reporter,
                                            const ReporterName & to_reporter,
                                            const FEProblemBase & from_problem,
                                            FEProblemBase & to_problem,
                                            unsigned int time_index)
{
  check(from_reporter, to_reporter, from_problem, to_problem, "transferReporter()");

  const ReporterData & from_data = from_problem.getReporterData();
  ReporterData & to_data = to_problem.getReporterData(ReporterData::WriteKey());
  from_data.getReporterContextBase(from_reporter).transfer(to_data, to_reporter, time_index);
}

void
ReporterTransferInterface::transferToVectorReporter(const ReporterName & from_reporter,
                                                    const ReporterName & to_reporter,
                                                    const FEProblemBase & from_problem,
                                                    FEProblemBase & to_problem,
                                                    dof_id_type index,
                                                    unsigned int time_index)
{
  check(from_reporter, to_reporter, from_problem, to_problem, "transferToVectorReporter()");

  const ReporterData & from_data = from_problem.getReporterData();
  ReporterData & to_data = to_problem.getReporterData(ReporterData::WriteKey());
  from_data.getReporterContextBase(from_reporter)
      .transferToVector(to_data, to_reporter, index, time_index);
}

void
ReporterTransferInterface::declareClone(const ReporterName & from_reporter,
                                        const ReporterName & to_reporter,
                                        const FEProblemBase & from_problem,
                                        FEProblemBase & to_problem,
                                        const ReporterMode & mode)
{
  check(from_reporter, from_problem, "declareClone()");

  const ReporterData & from_data = from_problem.getReporterData();
  ReporterData & to_data = to_problem.getReporterData(ReporterData::WriteKey());
  from_data.getReporterContextBase(from_reporter)
      .declareClone(to_data, to_reporter, mode, _rti_transfer);

  // Hide variables (if requested in parameters) if name is associated with a reporter object
  if (to_problem.hasUserObject(to_reporter.getObjectName()))
  {
    UserObject & uo = to_problem.getUserObject<UserObject>(to_reporter.getObjectName());
    Reporter * rep = dynamic_cast<Reporter *>(&uo);
    if (rep)
      rep->buildOutputHideVariableList({to_reporter.getCombinedName()});
  }
}

void
ReporterTransferInterface::declareVectorClone(const ReporterName & from_reporter,
                                              const ReporterName & to_reporter,
                                              const FEProblemBase & from_problem,
                                              FEProblemBase & to_problem,
                                              const ReporterMode & mode)
{
  check(from_reporter, from_problem, "declareVectorClone()");

  const ReporterData & from_data = from_problem.getReporterData();
  ReporterData & to_data = to_problem.getReporterData(ReporterData::WriteKey());

  from_data.getReporterContextBase(from_reporter)
      .declareVectorClone(to_data, to_reporter, mode, _rti_transfer);

  // Hide variables (if requested in parameters) if name is associated with a reporter object
  if (to_problem.hasUserObject(to_reporter.getObjectName()))
  {
    UserObject & uo = to_problem.getUserObject<UserObject>(to_reporter.getObjectName());
    Reporter * rep = dynamic_cast<Reporter *>(&uo);
    if (rep)
      rep->buildOutputHideVariableList({to_reporter.getCombinedName()});
  }
}

void
ReporterTransferInterface::resizeReporter(const ReporterName & name,
                                          FEProblemBase & problem,
                                          dof_id_type n)
{
  check(name, problem, "resizeReporter()");

  problem.getReporterData(ReporterData::WriteKey()).getReporterContextBase(name).resize(n);
}

std::vector<ReporterName>
ReporterTransferInterface::getReporterNamesHelper(std::string prefix,
                                                  const std::string & obj_name,
                                                  const std::vector<ReporterName> & rep_names)
{
  if (!prefix.empty())
    prefix += ":";
  std::vector<ReporterName> rnames;
  rnames.reserve(rep_names.size());
  for (const auto & rn : rep_names)
    rnames.emplace_back(obj_name, prefix + rn.getObjectName() + ":" + rn.getValueName());
  return rnames;
}

void
ReporterTransferInterface::check(const ReporterName & reporter,
                                 const FEProblemBase & problem,
                                 const std::string & method) const
{
  if (!problem.getReporterData().hasReporterValue(reporter))
    _rti_transfer.mooseError("In ",
                             method,
                             ": Reporter with the name \"",
                             reporter,
                             "\" within app \"",
                             problem.getMooseApp().name(),
                             "\" was not found.");
}

void
ReporterTransferInterface::check(const ReporterName & from_reporter,
                                 const ReporterName & to_reporter,
                                 const FEProblemBase & from_problem,
                                 const FEProblemBase & to_problem,
                                 const std::string & method) const
{
  if (!from_problem.getReporterData().hasReporterValue(from_reporter))
    _rti_transfer.mooseError("In ",
                             method,
                             ": Reporter with the name \"",
                             from_reporter,
                             "\"\n was not found within app \"",
                             from_problem.getMooseApp().name(),
                             "\".");

  if (!to_problem.getReporterData().hasReporterValue(to_reporter))
    _rti_transfer.mooseError("In ",
                             method,
                             ": Reporter with the name \"",
                             to_reporter,
                             "\"\n was not found within app \"",
                             to_problem.getMooseApp().name(),
                             "\".");
}
