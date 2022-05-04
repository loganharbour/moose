#include "ExecFlagRegistry.h"

#include "MooseUtils.h"

namespace moose
{
namespace internal
{

ExecFlagRegistry &
getExecFlagRegistry()
{
  static ExecFlagRegistry exec_flag_registry;
  return exec_flag_registry;
}

ExecFlagRegistry::ExecFlagRegistry()
  : GenericRegistry<std::string, MooseEnumItem>("ExecFlagRegistry")
{
}

const MooseEnumItem &
ExecFlagRegistry::registerFlag(const std::string & name)
{
  const auto create = [&name](const std::size_t id) { return MooseEnumItem(name, id); };
  const auto id = registerItem(MooseUtils::toUpper(name), create);
  return item(id);
}

bool
ExecFlagRegistry::isFlagRegistered(const MooseEnumItem & flag) const
{
  return keyExists(flag.name()) && item(id(flag.name())).id() == flag.id();
}

}
}
