//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"

class PerfGraph;
class PerfNode;

/**
 * Base class for storing static information about the PerfGraph for use in output.
 */
class PerfInfoBase
{
public:
  PerfInfoBase(const std::string & name,
               const unsigned long int num_calls,
               const Real self_time,
               const Real total_time,
               const Real root_time,
               const long int self_memory,
               const long int total_memory);

  /**
   * @returns The name
   */
  const std::string & name() const { return _name; }

  /**
   * @returns The number of calls
   */
  unsigned long int numCalls() const { return _num_calls; }

  /**
   * @returns The time spent without children
   */
  Real selfTime() const { return _self_time; }
  /**
   * @returns The average time spent without children
   */
  Real selfTimeAvg() const { return selfTime() / static_cast<Real>(numCalls()); }
  /**
   * @returns The percentage of time spent (relative to total) without children
   */
  Real selfTimePercent() const { return 100. * selfTime() / _root_time; }
  /**
   * @returns The memory gained without children
   */
  long int selfMemory() const { return _self_memory; }

  /**
   * @returns The time spent with children
   */
  Real totalTime() const { return _total_time; }
  /**
   * @returns The average time spent with children
   */
  Real totalTimeAvg() const { return totalTime() / static_cast<Real>(totalTime()); }
  /**
   * @returns The percentage of time spent (relative to total) with children
   */
  Real totalTimePercent() const { return 100. * totalTime() / _root_time; }
  /**
   * @returns The memory gained with children
   */
  long int totalMemory() const { return _total_memory; }

private:
  /// The name
  const std::string _name;
  /// Number of calls
  const unsigned long int _num_calls;
  /// Time spent without children
  const Real _self_time;
  /// Time spent with children
  const Real _total_time;
  /// The total time of the root (used in percentages)
  const Real _root_time;
  /// Memory gained without children
  const long int _self_memory;
  /// Memory gained with children
  const long int _total_memory;
};

/**
 * Represents a section in the PerfGraph tree, used for output
 */
class PerfTreeInfo : public PerfInfoBase
{
public:
  PerfTreeInfo(const PerfNode & node,
               const std::string & name,
               const unsigned int depth,
               const Real root_time);

  /**
   * @returns The depth in the tree
   */
  unsigned int depth() const { return _depth; }

  /**
   * @returns The children
   */
  const std::vector<std::shared_ptr<const PerfTreeInfo>> & children() const { return _children; }

  /**
   * Adds a child
   */
  void addChild(const std::shared_ptr<const PerfTreeInfo> & child_info)
  {
    _children.push_back(child_info);
  }

private:
  /// Depth in the tree
  const unsigned int _depth;
  /// Children in the tree
  std::vector<std::shared_ptr<const PerfTreeInfo>> _children;
};

/**
 * Represents a section in the PerfGraph, used for output
 */
class PerfSectionInfo : public PerfInfoBase
{
public:
  PerfSectionInfo(const std::string & name,
                  const unsigned long int num_calls,
                  const Real self_time,
                  const Real total_time,
                  const Real root_time,
                  const long int self_memory,
                  const long int total_memory);
};
