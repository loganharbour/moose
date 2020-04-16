[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 4
    ny = 4
    xmax = 8
    ymax = 8
  []
[]

# needed for a higher order quadrature for the UO
[Variables/u]
[]

[UserObjects]
  [./view_factor]
    type = UnobstructedPlanarViewFactor
    boundary = 'top left right bottom'
    execute_on = INITIAL
  [../]
[]

[Postprocessors]
  [./left_right]
    type = ViewFactorPP
    from_boundary = 'left'
    to_boundary = 'right'
    view_factor_object_name = view_factor
  [../]

  [./left_top]
    type = ViewFactorPP
    from_boundary = 'left'
    to_boundary = 'top'
    view_factor_object_name = view_factor
  [../]
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  csv = true
[]
