//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TestExtrudedConcentricCircleMeshGenerator.h"

#include "MooseApp.h"
#include "Factory.h"

registerMooseObject("MooseTestApp", TestExtrudedConcentricCircleMeshGenerator);

InputParameters
TestExtrudedConcentricCircleMeshGenerator::validParams()
{
  auto params = MeshGenerator::validParams();

  params.addRequiredParam<unsigned int>("num_sectors",
                                        "Number of azimuthal sectors in each quadrant");
  params.addRequiredParam<Real>("pitch", "The pitch for the concentric circle");

  params.addRequiredParam<std::vector<unsigned int>>("num_layers",
                                                     "The number of layers for the extrusion");
  params.addRequiredParam<std::vector<Real>>("heights", "The heights for the extrusion");

  return params;
}

TestExtrudedConcentricCircleMeshGenerator::TestExtrudedConcentricCircleMeshGenerator(
    const InputParameters & parameters)
  : MeshGenerator(parameters)
{
  // For the below, we first get the default parameters for each subgenerator that we
  // want to create. This allows us to then modify them as we please before creating
  // a subgenerator with said parameters
  {
    auto params = _app.getFactory().getValidParams("ConcentricCircleMeshGenerator");
    params.set<bool>("has_outer_square") = true;
    params.set<unsigned int>("num_sectors") = getParam<unsigned int>("num_sectors");
    params.set<bool>("preserve_volumes") = true;
    params.set<std::vector<unsigned int>>("rings") = {1, 2, 3, 4};
    params.set<std::vector<Real>>("radii") = {0.2, 0.3, 0.4};
    params.set<Real>("pitch") = getParam<Real>("pitch");

    // Here, we give it a name that is <name of this object>_circle, note that
    // this is used as the "input" for the extruder generator below
    addMeshSubgenerator("ConcentricCircleMeshGenerator", name() + "_circle", params);
  }

  {
    auto params = _app.getFactory().getValidParams("FancyExtruderGenerator");
    params.set<MeshGeneratorName>("input") = name() + "_circle";
    params.set<Point>("direction") = Point(0, 0, 1);
    params.set<std::vector<unsigned int>>("num_layers") =
        getParam<std::vector<unsigned int>>("num_layers");
    params.set<std::vector<Real>>("heights") = getParam<std::vector<Real>>("heights");

    // We store this mesh so that we can return it in generate()
    _extruded_mesh = &addMeshSubgenerator("FancyExtruderGenerator", name() + "_extruded", params);
  }
}

std::unique_ptr<MeshBase>
TestExtrudedConcentricCircleMeshGenerator::generate()
{
  // This generate() method will be called once the subgenerators that we depend on are
  // called, and we want the final product to be the extruded mesh so we will move
  // that now completed mesh and make it our final product
  return std::move(*_extruded_mesh);
}
