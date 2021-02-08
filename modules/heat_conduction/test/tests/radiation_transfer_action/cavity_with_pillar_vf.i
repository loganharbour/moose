[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Mesh]
  [./cartesian]
    type = CartesianMeshGenerator
    dim = 3
    dx = '0.1 0.3 0.4 0.3 0.1'
    #ix = '  5  15  20  15   5'
    ix = '  1  3  4  3   1'
    dy = '0.1 0.3 0.4 0.3 0.1'
    #iy = '  5  15  20  15   5'
    iy = '  1  3  4  3   1'
    dz = '0.1 0.8 0.2 0.1'
    #iz = '  5  40  10   5'
    iz = ' 1   8   2   1'
    subdomain_id = '
                    1 1 1 1 1
                    1 15 15 15 1
                    1 15 1 15 1
                    1 15 15 15 1
                    1 1 1 1 1

                    1 12 12 12 1
                    11 0 103 0  14
                    11 104 2 102  14
                    11 0 101 0  14
                    1 13 13 13 1

                    1 12 12 12 1
                    11 0 0 0  14
                    11 0 105 0  14
                    11 0 0 0  14
                    1 13 13 13 1

                    1 1 1 1 1
                    1 16 16 16 1
                    1 16 16 16 1
                    1 16 16 16 1
                    1 1 1 1 1
                   '
  [../]

  [./left_interior]
    type = SideSetsBetweenSubdomainsGenerator
    master_block = 11
    paired_block = '0 101 102 103 104 105'
    new_boundary = left_interior_wall
    input = cartesian
  [../]

  [./right_interior]
    type = SideSetsBetweenSubdomainsGenerator
    master_block = 14
    paired_block = '0 101 102 103 104 105'
    new_boundary = right_interior_wall
    input = left_interior
  [../]

  [./bottom_interior]
    type = SideSetsBetweenSubdomainsGenerator
    master_block = 12
    paired_block = '0 101 102 103 104 105'
    new_boundary = bottom_interior_wall
    input = right_interior
  [../]

  [./top_interior]
    type = SideSetsBetweenSubdomainsGenerator
    master_block = 13
    paired_block = '0 101 102 103 104 105'
    new_boundary = top_interior_wall
    input = bottom_interior
  [../]

  [./front_interior]
    type = SideSetsBetweenSubdomainsGenerator
    master_block = 15
    paired_block = '0 101 102 103 104 105'
    new_boundary = front_interior_wall
    input = top_interior
  [../]

  [./back_interior]
    type = SideSetsBetweenSubdomainsGenerator
    master_block = 16
    paired_block = '0 101 102 103 104 105'
    new_boundary = back_interior_wall
    input = front_interior
  [../]

  [./pillar_left]
    type = SideSetsBetweenSubdomainsGenerator
    master_block = 2
    paired_block = 104
    new_boundary = pillar_left
    input = 'back_interior'
  [../]

  [./pillar_right]
    type = SideSetsBetweenSubdomainsGenerator
    master_block = 2
    paired_block = 102
    new_boundary = pillar_right
    input = 'pillar_left'
  [../]

  [./pillar_bottom]
    type = SideSetsBetweenSubdomainsGenerator
    master_block = 2
    paired_block = 103
    new_boundary = pillar_bottom
    input = 'pillar_right'
  [../]

  [./pillar_top]
    type = SideSetsBetweenSubdomainsGenerator
    master_block = 2
    paired_block = 101
    new_boundary = pillar_top
    input = 'pillar_bottom'
  [../]

  [./pillar_back]
    type = SideSetsBetweenSubdomainsGenerator
    master_block = 2
    paired_block = 105
    new_boundary = pillar_back
    input = 'pillar_top'
  [../]

  [./rename_block]
    type = RenameBlockGenerator
    old_block_id = '2 11 12 13 14 15 16 101 102 103 104 105'
    new_block_id = '2 1  1  1  1  1  1  0   0   0   0   0'
    input = 'pillar_back'
  [../]

  #[./delete_0]
  #  type = BlockDeletionGenerator
  #  block_id = 0
  #  input = rename_block
  #[../]
[]

[GrayDiffuseRadiation]
  [./cavity]
    sidesets = '6 7 8 9 10 11 12 13 14 15 16'
    emissivity = '0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8'
    n_patches = '5 5 5 5 5 5 5 5 5 5 5'
    partitioners = 'metis metis metis metis metis metis metis metis metis metis metis'
    #n_patches = '10 10 10 10 10 10 10 10 10 10 10'
    #partitioners = 'centroid centroid centroid centroid centroid centroid centroid centroid centroid centroid centroid'
    #centroid_partitioner_directions = 'x x x x x x x x x x x'
    temperature = temperature
    ray_tracing_face_order = SECOND
    #unobstructed_cavity = true
  [../]
[]

[Variables]
  [./temperature]
    initial_condition = 300
    block = '1 2'
  [../]
[]

[Kernels]
  #[./hctd]
  #  type = HeatConductionTimeDerivative
  #  variable = temperature
  #  block = '1 2'
  #[../]

  [./hc]
    type = HeatConduction
    variable = temperature
    block = '1 2'
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = temperature
    boundary = left
    value = 500
  [../]

  [./front]
    type = DirichletBC
    variable = temperature
    boundary = front
    value = 300
  [../]
[]

[Materials]
  [./hcmat]
    type = HeatConductionMaterial
    thermal_conductivity = 25.0
    specific_heat = 490.0
    block = '1 2'
  [../]

  [./density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '80'
    block = '1 2'
  [../]
[]

[Executioner]
  type = Steady #Transient
  #end_time = 5000
  #dt = 10

  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
[]

[Outputs]
  exodus = true
[]
