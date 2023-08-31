//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <string>

/**
 * Represents a mesh generator in a tree for use in traversing said tree
 * for accessing mesh meta data.
 *
 * Needed because the mesh generators are not constructed in restart/recover,
 * so we need a way to access this data.
 */
class MeshGeneratorNode
{
public:
  MeshGeneratorNode(const std::string & name, const std::string & type);

  /**
   * Comparator for nodes that sorts by name
   */
  struct Comparator
  {
    bool operator()(const MeshGeneratorNode * const & a, const MeshGeneratorNode * const & b) const
    {
      return a->name() < b->name();
    }
  };

  /**
   * @return The name of the mesh generator
   */
  const std::string & name() const { return _name; }
  /**
   * @return The type of the mesh generator
   */
  const std::string & type() const { return _type; }

  /**
   * Adds the node \p node as a child to this mesh generator
   */
  void addChild(const MeshGeneratorNode & node);

  void addParent(const MeshGeneratorNode & node);

  bool isChild(const MeshGeneratorNode & node) const { return _children.count(&node); }
  bool isParent(const MeshGeneratorNode & node) const { return _parents.count(&node); }

private:
  /// The name of the mesh generator
  const std::string _name;
  /// The type of the mesh generator
  const std::string _type;
  /// The children of this mesh generator
  std::set<const MeshGeneratorNode *, Comparator> _children;
  /// The parents of this mesh generator
  std::set<const MeshGeneratorNode *, Comparator> _parents;
};
