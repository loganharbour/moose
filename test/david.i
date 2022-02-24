[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 5
  []
[]

[Variables/u]
  order = CONSTANT
  family = MONOMIAL
[]

[DGKernels/dg]
  type = DGDiffusion
  sigma = 1
  epsilon = 1
  variable = u
[]

[UserObjects/ghosting]
  type = GhostingUserObject
[]

[Executioner]
  type = Steady
[]

[Problem]
  solve = false
[]
