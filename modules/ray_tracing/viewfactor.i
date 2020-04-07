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

[AuxVariables/pid]
  order = CONSTANT
  family = MONOMIAL
[]

[AuxKernels/pid]
  type = ProcessorIDAux
  variable = pid
[]

[UserObjects/viewfactor]
  type = ViewFactorRayStudy
  execute_on = initial
  boundary = 'top left right bottom'
  ray_kernel_coverage_check = false
[]

[RayBCs/viewfactor]
  type = ViewFactorRayBC
  boundary = 'top left right bottom'
[]

[Executioner]
  type = Steady
[]

[Problem]
  solve = false
[]

[Outputs]
  exodus = true
[]
