//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GenericRegistry.h"
#include "MooseEnumItem.h"

#include <shared_mutex>

namespace moose
{
namespace internal
{

class ExecFlagRegistry;

/**
 * Get the global ExecFlagRegistry singleton.
 */
ExecFlagRegistry & getExecFlagRegistry();

/**
 * Registry for statically defining execute flags with consistent numbering.
 */
class ExecFlagRegistry : private GenericRegistry<std::string, MooseEnumItem>
{
public:
  /**
   * Registers an execute flag with the given name.
   */
  const MooseEnumItem & registerFlag(const std::string & name);

  /**
   * \returns True if the execute flag \p flag is registered, false otherwise.
   */
  bool isFlagRegistered(const MooseEnumItem & flag) const;

private:
  ExecFlagRegistry();

  /// So it can be constructed
  friend ExecFlagRegistry & getExecFlagRegistry();
};

}
}

#define defineExecFlag(flag) moose::internal::getExecFlagRegistry().registerFlag(flag)
