[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 2
  ymin = 0
  ymax = 1
  nx = 2
  ny = 2
  nz = 2
[]

[UserObjects]
  [./view_factor_study]
    type = AQViewFactorRayStudy
    execute_on = initial
    boundary = 'left right bottom top'
    face_order = TENTH
    polar_quad_order = 80
    always_cache_traces = true
    aux_data_on_cache_traces = true
  [../]

  [./view_factor]
    type = AQRayTracingViewFactor
    boundary = 'left right bottom top'
    execute_on = INITIAL
    normalize_view_factor = false
    ray_study_name = view_factor_study
  [../]
[]

[RayBCs/viewfactor]
  type = AQViewFactorRayBC
  boundary = 'left right bottom top'
[]

##
## Reference: bottom -> left/right = 0.19098
##            bottom -> top = 0.61803
## Result at spatial order 20, angular order 200 & -r2
##            bottom -> left/right = 0.1911
##            bottom -> top = 0.6177
##
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

  [./left_bottom]
    type = ViewFactorPP
    from_boundary = left
    to_boundary = bottom
    view_factor_object_name = view_factor
  [../]

  [./bottom_left]
    type = ViewFactorPP
    from_boundary = bottom
    to_boundary = left
    view_factor_object_name = view_factor
  [../]

  [./bottom_right]
    type = ViewFactorPP
    from_boundary = bottom
    to_boundary = right
    view_factor_object_name = view_factor
  [../]

  [./bottom_top]
    type = ViewFactorPP
    from_boundary = bottom
    to_boundary = top
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
