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

[UserObjects]
  [./view_factor_study]
    type = ViewFactorRayStudy
    execute_on = INITIAL
    boundary = 'top left right bottom'
    face_order = TENTH
  [../]

  [./view_factor]
    type = RayTracingViewFactor
    boundary = 'top left right bottom'
    execute_on = INITIAL
    ray_study_name = view_factor_study
  [../]
[]

[RayBCs/viewfactor]
  type = ViewFactorRayBC
  boundary = 'top left right bottom'
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
