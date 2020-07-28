//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ViewFactorBase.h"
#include "libmesh/quadrature.h"

#include <limits>

InputParameters
ViewFactorBase::validParams()
{
  InputParameters params = SideUserObject::validParams();
  params.addParam<Real>("view_factor_tol",
                        std::numeric_limits<Real>::max(),
                        "Tolerance for checking view factors. Default is to allow everything.");
  params.addParam<bool>(
      "print_view_factor_info", true, "Flag to print information about computed view factors.");
  params.addParam<bool>("normalize_view_factor",
                        true,
                        "Determines if view factors are normalized to sum to one (consistent with "
                        "their definition).");
  params.addClassDescription(
      "A base class for automatic computation of view factors between sidesets.");
  return params;
}

ViewFactorBase::ViewFactorBase(const InputParameters & parameters)
  : SideUserObject(parameters),
    _n_sides(boundaryIDs().size()),
    _areas(_n_sides),
    _view_factor_tol(getParam<Real>("view_factor_tol")),
    _normalize_view_factor(getParam<bool>("normalize_view_factor")),
    _print_view_factor_info(getParam<bool>("print_view_factor_info"))
{
  // sizing the view factor array
  _view_factors.resize(_n_sides);
  for (auto & v : _view_factors)
    v.resize(_n_sides);

  // set up the map from the side id to the local index & side name to local index
  std::vector<BoundaryName> boundary_names = getParam<std::vector<BoundaryName>>("boundary");
  for (unsigned int j = 0; j < boundary_names.size(); ++j)
    _side_name_index[boundary_names[j]] = j;
}

unsigned int
ViewFactorBase::getSideNameIndex(std::string name) const
{
  auto it = _side_name_index.find(name);
  if (it == _side_name_index.end())
    mooseError("Boundary ", name, " does not exist.");
  return it->second;
}

Real
ViewFactorBase::getViewFactor(BoundaryID from_id, BoundaryID to_id) const
{
  auto from_name = _mesh.getBoundaryName(from_id);
  auto to_name = _mesh.getBoundaryName(to_id);

  return getViewFactor(from_name, to_name);
}

Real
ViewFactorBase::getViewFactor(BoundaryName from_name, BoundaryName to_name) const
{
  auto from = _side_name_index.find(from_name);
  auto to = _side_name_index.find(to_name);
  if (from == _side_name_index.end())
    mooseError("Boundary id ",
               _mesh.getBoundaryID(from_name),
               " with name ",
               from_name,
               " not listed in boundary parameter.");

  if (to == _side_name_index.end())
    mooseError("Boundary id ",
               _mesh.getBoundaryID(to_name),
               " with name ",
               to_name,
               " not listed in boundary parameter.");

  return _view_factors[from->second][to->second];
}

void
ViewFactorBase::finalize()
{
  // do some communication before finalizing view_factors
  for (unsigned int i = 0; i < _n_sides; ++i)
    gatherSum(_view_factors[i]);

  finalizeViewFactor();
  checkAndNormalizeViewFactor();
}

void
ViewFactorBase::threadJoin(const UserObject & y)
{
  const ViewFactorBase & pps = static_cast<const ViewFactorBase &>(y);
  for (unsigned int i = 0; i < _n_sides; ++i)
  {
    for (unsigned int j = 0; j < _n_sides; ++j)
      _view_factors[i][j] += pps._view_factors[i][j];
  }
  threadJoinViewFactor(y);
}

void
ViewFactorBase::checkAndNormalizeViewFactor()
{
  // check view factors
  for (unsigned int from = 0; from < _n_sides; ++from)
  {
    Real s = 0;
    for (unsigned int to = 0; to < _n_sides; ++to)
      s += _view_factors[from][to];

    if (_print_view_factor_info)
      _console << "View factors from sideset " << boundaryNames()[from] << " sum to " << s
               << std::endl;

    if (std::abs(1 - s) > _view_factor_tol)
      mooseError("View factor from boundary ", boundaryNames()[from], " add to ", s);
  }

  // normalize view factors
  if (_normalize_view_factor)
  {
    // compute number of entries
    unsigned int ne = _n_sides * (_n_sides + 3) / 2;

    // allocate space
    DenseVector<Real> rhs(ne);
    DenseVector<Real> corrections(ne);
    DenseMatrix<Real> matrix(ne, ne);

    // equations for the Lagrange multiplier
    unsigned int row = 0;
    for (unsigned int i = 0; i < _n_sides; ++i)
    {
      rhs(row) = 1;
      for (unsigned int j = 0; j < _n_sides; ++j)
        rhs(row) -= _view_factors[i][j];

      matrix(row, indexHelper(i, i)) = 1;
      for (unsigned int j = 0; j < i; ++j)
        matrix(row, indexHelper(j, i)) = _areas[j] / _areas[i];

      for (unsigned int j = i + 1; j < _n_sides; ++j)
        matrix(row, indexHelper(i, j)) = 1;

      ++row;
    }

    // equations for the delta_ii, i.e. corrections for diagonal elements
    for (unsigned int i = 0; i < _n_sides; ++i)
    {
      matrix(row, i) = 1;
      matrix(row, indexHelper(i, i)) = 2;
      ++row;
    }

    // equations for the delta_ij, j > i, i.e. the corrections for the off-diagonal elements
    for (unsigned int i = 0; i < _n_sides; ++i)
      for (unsigned int j = i + 1; j < _n_sides; ++j)
      {
        Real ar = _areas[i] / _areas[j];
        matrix(row, i) = 1 + ar;
        matrix(row, indexHelper(i, j)) = 2 * (1 + ar * ar);
        ++row;
      }

    // solve the linear system
    matrix.lu_solve(rhs, corrections);

    // apply the corrections
    for (unsigned int from = 0; from < _n_sides; ++from)
      for (unsigned int to = 0; to < _n_sides; ++to)
      {
        if (from <= to)
          _view_factors[from][to] += corrections(indexHelper(from, to));
        else
          _view_factors[from][to] += corrections(indexHelper(to, from)) * _areas[to] / _areas[from];
      }
  }

  for (unsigned int from = 0; from < _n_sides; ++from)
  {
    std::string from_name;
    for (auto pair : _side_name_index)
      if (pair.second == from)
        from_name = pair.first;
    auto from_id = _mesh.getBoundaryID(from_name);

    for (unsigned int to = 0; to < _n_sides; ++to)
    {
      std::string to_name;
      for (auto pair : _side_name_index)
        if (pair.second == to)
          to_name = pair.first;
      auto to_id = _mesh.getBoundaryID(to_name);
      _console << from_name << " (" << from_id << ") -> " << to_name << " (" << to_id
               << ") = " << _view_factors[from][to] << std::endl;
    }
  }
}

unsigned int
ViewFactorBase::indexHelper(unsigned int i, unsigned int j) const
{
  mooseAssert(i <= j, "indexHelper requires i <= j");
  if (i == j)
    return _n_sides + i;
  unsigned int pos = 2 * _n_sides;
  for (unsigned int l = 0; l < _n_sides; ++l)
    for (unsigned int m = l + 1; m < _n_sides; ++m)
    {
      if (l == i && m == j)
        return pos;
      else
        ++pos;
    }
  mooseError("Should never get here");
  return 0;
}
