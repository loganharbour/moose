[Mesh]
  type = GeneratedMesh
  dim = 3
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  zmin = 0
  zmax = 1
  nx = 2
  ny = 2
  nz = 2
[]

[UserObjects]
  [./view_factor_study]
    type = PointToPointViewFactorRayStudy
    execute_on = initial
    boundary = 'left right front back bottom top'
    face_order = FIFTH
  [../]

  [./view_factor]
    type = RayTracingViewFactor
    boundary = 'left right front back bottom top'
    execute_on = INITIAL
    ray_study_name = view_factor_study
  [../]
[]

[RayBCs/viewfactor]
  type = PointToPointViewFactorRayBC
  boundary = 'left right front back bottom top'
[]

[Postprocessors]
  [./left_right]
    type = ViewFactorPP
    from_boundary = left
    to_boundary = right
    view_factor_object_name = view_factor
  [../]

  [./left_top]
    type = ViewFactorPP
    from_boundary = left
    to_boundary = top
    view_factor_object_name = view_factor
  [../]

  [./left_back]
    type = ViewFactorPP
    from_boundary = left
    to_boundary = back
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
