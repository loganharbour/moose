[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 100
    ny = 100
  []
[]

[UserObjects/study]
  type = MPMBoundingBoxRayStudy
  bbox_min = '0.05 0.8 0'
  bbox_max = '0.3 0.95 0'
  grid_intervals = '5 5'

  execute_on = TIMESTEP_END

  displacement = 'disp_x disp_y'
  mass = mass

  ray_kernel_coverage_check = false
[]


[AuxVariables]
  [disp_x]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.05
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE
    initial_condition = -0.1
  []
  [mass]
    order = FIRST
    family = LAGRANGE
  []
[]

[Problem]
  solve = False
[]

[Executioner]
  type = Transient
  num_steps = 5
[]

[Outputs]
  exodus = true
[]
