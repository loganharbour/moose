#pragma once

#include "GeneralUserObject.h"

template <typename T>
class MooseVariableFE;

class GavinRenameMe : public GeneralUserObject
{
public:
  GavinRenameMe(const InputParameters & params);

  static InputParameters validParams();

  void initialize() override final {}
  void finalize() override final {}
  void execute() override final;

private:
  AuxiliarySystem & _aux;
  const unsigned int _N;
  MooseVariableFE<Real> & _from_var;
  MooseVariableFE<Real> & _to_var;
};
