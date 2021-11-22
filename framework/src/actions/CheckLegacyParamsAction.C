//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CheckLegacyParamsAction.h"

#include "AppFactory.h"
#include "Registry.h"

registerMooseAction("MooseApp", CheckLegacyParamsAction, "check_legacy_params");

InputParameters
CheckLegacyParamsAction::validParams()
{
  auto params = Action::validParams();
  params.addClassDescription("Checks whether or not objects exist that are constructed with the "
                             "legacy input parameter construction");
  return Action::validParams();
}

CheckLegacyParamsAction::CheckLegacyParamsAction(InputParameters params) : Action(params) {}

void
CheckLegacyParamsAction::act()
{
  // Not a big fan of testing objects within the object itself... but this is a temp object
  // and this is the easiest way to do it. MooseTestApp uses this parameter for the test
  const bool for_test = _app.parameters().have_parameter<bool>("test_check_legacy_params") &&
                        _app.parameters().get<bool>("test_check_legacy_params");

  std::set<std::pair<std::string, std::string>> objects;

  // Get the MooseObjects and Actions whose input parameters are constructed
  // using the legacy method, skipping those registered to MooseApp
  // (which is required so that apps that use legacy params from framework
  // objects still pass while we deprecate)
  for (const auto & object_label_pair :
       moose::internal::getRegistry().getLegacyConstructedObjects())
    if (object_label_pair.second != "MooseApp" || for_test)
      objects.insert(object_label_pair);

  // Get the applications whose input parameters are constructed using the legacy
  // method, skipping only MooseApp (which we need for now for deprecation)
  std::stringstream apps_out;
  for (const auto & app_name : AppFactory::instance().getLegacyConstructedApps())
    if (app_name != "MooseApp")
      objects.insert(std::make_pair(app_name, ""));

  if (objects.size())
  {
    std::stringstream warning;
    warning << "The following object(s) are constructed using the legacy input parameter "
               "construction:\n\n";
    for (const auto & object_label_pair : objects)
    {
      warning << "  " << object_label_pair.first;
      if (object_label_pair.second.size())
        warning << " (" << object_label_pair.second << ")";
      warning << "\n";
    }

    warning
        << "\nConvert InputParameters validParams<T>() for each object into a static"
        << "\nmember function InputParameters T::validParams() and remove the old function."
        << "\n\nSee mooseframework.org/newsletter/2021_11.html#legacy-input-parameter-deprecation"
        << "\nfor more information.\n";
    mooseError(warning.str());
  }
}
