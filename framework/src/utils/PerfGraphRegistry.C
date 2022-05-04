#include "PerfGraphRegistry.h"

#include "DataIO.h"

namespace moose
{
namespace internal
{

PerfGraphRegistry &
getPerfGraphRegistry()
{
  // In C++11 this is even thread safe! (Lookup "Static Initializers")
  static PerfGraphRegistry perf_graph_registry_singleton;

  return perf_graph_registry_singleton;
}

PerfGraphRegistry::PerfGraphRegistry()
  : GenericRegistry<std::string, PerfGraphSectionInfo>("PerfGraphRegistry")
{
  // Reserve space so that re-allocation doesn't need to happen much
  // This does not take much memory and, for most cases, will keep a single
  // reallocation from happening
  _id_to_section_info.reserve(5000);
}

unsigned int
PerfGraphRegistry::registerSection(const std::string & section_name, unsigned int level)
{
  return actuallyRegisterSection(section_name, level, "", false);
}

PerfID
PerfGraphRegistry::registerSection(const std::string & section_name,
                                   unsigned int level,
                                   const std::string & live_message,
                                   const bool print_dots)
{
  if (section_name.empty())
    mooseError("Section name not provided when registering timed section!");
  if (live_message.empty())
    mooseError("Live message not provided when registering timed section!");

  return actuallyRegisterSection(section_name, level, live_message, print_dots);
}

PerfID
PerfGraphRegistry::actuallyRegisterSection(const std::string & section_name,
                                           unsigned int level,
                                           const std::string & live_message,
                                           const bool print_dots)
{
  const auto create =
      [this, &section_name, &level, &live_message, &print_dots](const std::size_t id)
  {
    const PerfGraphSectionInfo info(id, section_name, level, live_message, print_dots);

    // Also register in _id_to_section_info
    {
      std::unique_lock write_lock(_id_to_section_info_mutex);
      mooseAssert(_id_to_section_info.size() == id, "Section is already inserted");
      _id_to_section_info.push_back(info);
    }

    return info;
  };

  return registerItem(section_name, create);
  ;
}

const PerfGraphSectionInfo &
PerfGraphRegistry::readSectionInfo(PerfID section_id)
{
  return _id_to_section_info[section_id];
}

}
}

void
dataStore(std::ostream & stream, moose::internal::PerfGraphSectionInfo & info, void * context)
{
  dataStore(stream, info._id, context);
  dataStore(stream, info._name, context);
  dataStore(stream, info._level, context);
  dataStore(stream, info._live_message, context);
  dataStore(stream, info._print_dots, context);
}

void
dataLoad(std::istream & stream, moose::internal::PerfGraphSectionInfo & info, void * context)
{
  dataLoad(stream, info._id, context);
  dataLoad(stream, info._name, context);
  dataLoad(stream, info._level, context);
  dataLoad(stream, info._live_message, context);
  dataLoad(stream, info._print_dots, context);
}
