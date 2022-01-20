[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 10
  []
[]

[AuxVariables]
  [from]
    order = CONSTANT
    family = MONOMIAL
  []
  [to]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[UserObjects/gavin]
  type = GavinRenameMe
  from_var = from
  to_var = to
  N = 5
  execute_on = INITIAL
[]

[Problem]
  solve = False
[]

[Executioner]
  type = Steady
[]
