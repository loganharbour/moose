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
    type = ViewFactorRayStudy2
    execute_on = initial
    boundary = 'left right front back bottom top'
    face_order = SECOND
    polar_quad_order = 12
    azimuthal_quad_order = 24
    #face_type = Gauss
    always_cache_traces = true
    aux_data_on_cache_traces = true
  [../]

  [./view_factor]
    type = RayTracingViewFactor2
    boundary = 'left right front back bottom top'
    execute_on = INITIAL
    ray_study_name = view_factor_study
    normalize_view_factor = false
    print_view_factor_info = true
  [../]
[]

[RayBCs/viewfactor]
  type = ViewFactorRayBC2
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
  [rays]
    type = RayTracingExodus
    study = view_factor_study
    output_aux_data = true
    execute_on = final
  []
[]
