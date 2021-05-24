//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ReceiverVectorPostprocessor.h"

registerMooseObject("MooseApp", ReceiverVectorPostprocessor);

InputParameters
ReceiverVectorPostprocessor::validParams()
{
  InputParameters params = GeneralVectorPostprocessor::validParams();
  params.addClassDescription(
      "Reports the value stored in this vector postprocessor, which is usually filled "
      "in by another object. The receiver does not compute its own value.");

  params.addRequiredParam<std::vector<std::string>>("vector_names",
                                                    "Names of the column vectors in this object");

  params.addParam<std::vector<std::vector<Real>>>(
      "default_values",
      "Default values to fill. May be empty. Leading dimension must be equal to leading dimension "
      "of vector_names parameter.");

  return params;
}

ReceiverVectorPostprocessor::ReceiverVectorPostprocessor(const InputParameters & parameters)
  : GeneralVectorPostprocessor(parameters)
{
  const auto vectors = declareVectors(getParam<std::vector<std::string>>("vector_names"));

  if (isParamValid("default_values"))
  {
    const auto & default_values = getParam<std::vector<std::vector<Real>>>("default_values");
    if (default_values.size() != vectors.size())
      paramError("default_values",
                 "Leading dimension must be the equal to the number of declared vectors (",
                 vectors.size(),
                 ")");

    for (const auto i : index_range(vectors))
      (*vectors[i]) = default_values[i];
  }
}
