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

[UserObjects]
  [./viewfactor_study]
    type = ViewFactorRayStudy
    execute_on = initial
    boundary = 'top left right bottom'
    ray_kernel_coverage_check = false
    force_preaux = true
    face_order = SECOND
  [../]

  [./unobstructed_viewfactor_uo]
    type = UnobstructedPlanarViewFactor
    boundary = 'top left right bottom'
    execute_on = INITIAL
    normalize_view_factor = false
  [../]

  [./viewfactor_uo]
    type = RayTracingViewFactor
    boundary = 'top left right bottom'
    execute_on = INITIAL
    ray_study_name = viewfactor_study
    normalize_view_factor = false
  [../]
[]

[RayBCs/viewfactor]
  type = ViewFactorRayBC
  boundary = 'top left right bottom'
[]

[Executioner]
  type = Steady
[]

[Variables/u]
[]

[Postprocessors]
  [./vf_unob_left_right]
    type = ViewFactorPP
    from_boundary = 'left'
    to_boundary = 'right'
    view_factor_object_name = unobstructed_viewfactor_uo
  [../]

  [./vf_unob_left_top]
    type = ViewFactorPP
    from_boundary = 'left'
    to_boundary = 'top'
    view_factor_object_name = unobstructed_viewfactor_uo
  [../]

  [./vf_left_right]
    type = ViewFactorPP
    from_boundary = 'left'
    to_boundary = 'right'
    view_factor_object_name = viewfactor_uo
  [../]

  [./vf_left_top]
    type = ViewFactorPP
    from_boundary = 'left'
    to_boundary = 'top'
    view_factor_object_name = viewfactor_uo
  [../]
[]

[Problem]
  kernel_coverage_check = false
[]

[Outputs]
  exodus = true
[]
