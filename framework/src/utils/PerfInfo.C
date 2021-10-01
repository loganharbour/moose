//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PerfInfo.h"

#include "PerfNode.h"

PerfInfoBase::PerfInfoBase(const std::string & name,
                           const unsigned long int num_calls,
                           const Real self_time,
                           const Real total_time,
                           const Real root_time,
                           const long int self_memory,
                           const long int total_memory)
  : _name(name),
    _num_calls(num_calls),
    _self_time(self_time),
    _total_time(total_time),
    _root_time(root_time),
    _self_memory(self_memory),
    _total_memory(total_memory)
{
}

PerfTreeInfo::PerfTreeInfo(const PerfNode & node,
                           const std::string & name,
                           const unsigned int depth,
                           const Real root_time)
  : PerfInfoBase(name,
                 node.numCalls(),
                 node.selfTimeSec(),
                 node.totalTimeSec(),
                 root_time,
                 node.selfMemory(),
                 node.totalMemory()),
    _depth(depth)
{
}

PerfSectionInfo::PerfSectionInfo(const std::string & name,
                                 const unsigned long int num_calls,
                                 const Real self_time,
                                 const Real total_time,
                                 const Real root_time,
                                 const long int self_memory,
                                 const long int total_memory)
  : PerfInfoBase(name, num_calls, self_time, total_time, root_time, self_memory, total_memory)
{
}
