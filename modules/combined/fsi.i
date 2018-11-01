[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 1
  xmax = 2
  displacements = 'disp_x disp_y'
[]

[MeshModifiers]
  [./subdomain1]
    type = SubdomainBoundingBox
    bottom_left = '1.0 0 0'
    block_id = 1
    top_right = '2.0 1.0 0'
  [../]
  [./interface]
    type = SideSetsBetweenSubdomains
    depends_on = subdomain1
    master_block = '0'
    paired_block = '1'
    new_boundary = 'master0_interface'
  [../]
[]

[Variables]
  [./u]
    order = FIRST
    family = LAGRANGE
    block = '0'
  [../]


  [./disp_x]
    order = FIRST
    family = LAGRANGE
    # block = 1
  [../]
  [./vx]
    block = 1
  [../]
[]

[AuxVariables]
  [./disp_y]
    initial_condition = 0
    block = 1
  [../]
[]

[Kernels]
  [./diff_u]
    type = MatDiffusion
    variable = u
    D_name = 4
    block = 0
    use_displaced_mesh = true
  [../]
  [./diff_disp_x_left]
    type = Diffusion
    variable = disp_x
    block = 0
  [../]
  [./diff_disp_x]
    type = MatDiffusion
    variable = disp_x
    D_name = 10
    block = 1
  [../]
  [./accel]
    type = CoupledTimeDerivative
    variable = disp_x
    v = vx
    block = 1
  [../]
  [./vx_time_derivative_term]
    type = CoupledTimeDerivative
    variable = vx
    v = disp_x
    block = 1
  [../]
  [./source_vx]
    type = MatReaction
    variable = vx
    block = 1
    mob_name = 1
  [../]
[]

[InterfaceKernels]
  [./penalty_interface]
    type = CoupledPenaltyInterfaceDiffusion
    variable = u
    neighbor_var = disp_x
    slave_coupled_var = vx
    boundary = master0_interface
    penalty = 1e6
  [../]
[]

[BCs]
  [./left_velocity]
    type = FunctionDirichletBC
    variable = u
    boundary = 'left'
    function = '1 + cos(t)'
  [../]
  [./left_disp_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  [../]
  [./right_disp_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'right'
    value = 0
  [../]
  [./right_vx]
    type = DirichletBC
    variable = vx
    boundary = 'right'
    value = 0
  [../]
[]

[Materials]
  [./block0]
    type = GenericConstantMaterial
    block = '0'
    prop_names = 'D'
    prop_values = '4'
  [../]
  [./block1]
    type = GenericConstantMaterial
    block = '1'
    prop_names = 'D'
    prop_values = '2'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 100
  # num_steps = 3
  dt = 0.05
  dtmin = 0.05
  solve_type = PJFNK
  petsc_options = '-options_left -pc_svd_monitor'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = 'none'
[]

[Outputs]
  exodus = true
  print_linear_residuals = true
  [./dofmap]
    type = DOFMap
    execute_on = 'initial'
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]
