//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeshGeneratorNode.h"

MeshGeneratorNode::MeshGeneratorNode(const std::string & name, const std::string & type)
  : _name(name), _type(type)
{
}

void
MeshGeneratorNode::addChild(const MeshGeneratorNode & node)
{
  _children.insert(&node);
}
void
MeshGeneratorNode::addParent(const MeshGeneratorNode & node)
{
  _parents.insert(&node);
}
