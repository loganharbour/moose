//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include <map>
#include <vector>

class MooseApp;

namespace moose
{
namespace internal
{

template <typename T>
class StaticWarehouse
{
public:
  const std::vector<T *> & get(const MooseApp & app) const;

private:
  void add(const MooseApp & app, T & object);

  friend T;

  std::map<const MooseApp *, std::vector<T *>> _objects;
};

template <typename T>
void
StaticWarehouse<T>::add(const MooseApp & app, T & object)
{
  _objects[&app].push_back(&object);
}

template <typename T>
const std::vector<T *> &
StaticWarehouse<T>::get(const MooseApp & app) const
{
  const auto it = _objects.find(&app);
  if (it != _objects.end())
    return it->second;
  static std::vector<T *> empty;
  return empty;
}

template <typename T>
StaticWarehouse<T> &
getStaticWarehouse()
{
  static StaticWarehouse<T> warehouse;
  return warehouse;
}

}
}
