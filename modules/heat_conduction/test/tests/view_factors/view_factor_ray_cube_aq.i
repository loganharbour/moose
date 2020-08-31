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
    type = AQViewFactorRayStudy
    execute_on = initial
    boundary = 'left right front back bottom top'
    face_order = FIFTH
    polar_quad_order = 12
    azimuthal_quad_order = 4
    always_cache_traces = true
    aux_data_on_cache_traces = true
  [../]

  [./view_factor]
    type = AQRayTracingViewFactor
    boundary = 'left right front back bottom top'
    execute_on = INITIAL
    ray_study_name = view_factor_study
    normalize_view_factor = false
    print_view_factor_info = true
  [../]
[]

[RayBCs/viewfactor]
  type = AQViewFactorRayBC
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
