#include "GavinRenameMe.h"

#include "MooseVariableFE.h"

registerMooseObject("MooseApp", GavinRenameMe);

InputParameters
GavinRenameMe::validParams()
{
  auto params = GeneralUserObject::validParams();

  params.addRequiredParam<unsigned int>("N", "Foo");
  params.addRequiredParam<AuxVariableName>("from_var", "Foo");
  params.addRequiredParam<AuxVariableName>("to_var", "Foo");

  return params;
}

GavinRenameMe::GavinRenameMe(const InputParameters & params)
  : GeneralUserObject(params),
    _aux(_fe_problem.getAuxiliarySystem()),
    _N(getParam<unsigned int>("N")),
    _from_var(dynamic_cast<MooseVariableFE<Real> &>(
        _aux.getVariable(0, params.get<AuxVariableName>("from_var")))),
    _to_var(dynamic_cast<MooseVariableFE<Real> &>(
        _aux.getVariable(0, params.get<AuxVariableName>("to_var"))))
{
}

void
GavinRenameMe::execute()
{
  // get the current element values from from_var
  std::set<MooseVariableFEBase *> needed_moose_vars({&_from_var});
  _fe_problem.setActiveElementalMooseVariables(needed_moose_vars, 0);
  _from_var.prepareAux();
  for (const auto & elem : *_fe_problem.mesh().getActiveLocalElementRange())
  {
    _fe_problem.prepare(elem, 0);
    _fe_problem.setCurrentSubdomainID(elem, 0);
    _fe_problem.reinitElem(elem, 0);

    const auto from_var_val = _from_var.sln()[0];
    std::cerr << elem->id() << " from val = " << from_var_val << std::endl;
  }

  // collect the elements that need to be marked here

  // add to the aux; here I'm adding a value to all of them (the element ID)
  // but obviously you'll only add to a subset of them and will set 1 instead.
  // you may need to zero the other values... I don't remember
  needed_moose_vars = {&_to_var};
  _to_var.prepareAux();
  _fe_problem.setActiveElementalMooseVariables(needed_moose_vars, 0);
  for (const auto & elem : *_fe_problem.mesh().getActiveLocalElementRange())
  {
    _fe_problem.prepare(elem, 0);
    _fe_problem.setCurrentSubdomainID(elem, 0);
    _fe_problem.reinitElem(elem, 0);

    _to_var.setNodalValue(elem->id());
    _to_var.insert(_aux.solution());
  }

  _aux.solution().close();
  _fe_problem.clearActiveElementalMooseVariables(0);
}
