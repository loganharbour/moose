[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 1
  []
[]

[UserObjects/test]
  type = ExecFlagRegistryTest
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]
