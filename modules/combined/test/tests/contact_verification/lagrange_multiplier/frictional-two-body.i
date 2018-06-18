[GlobalParams]
  displacements = 'disp_x disp_y'
  D_name = 1e6
[]

[Mesh]
  file = two-body-no-sep.e
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./disp_x]
    block = '1 2'
  [../]
  [./disp_y]
    block = '1 2'
  [../]
  [./lm]
    block = 3
  [../]
  [./tangent_lm]
    block = 3
  [../]
  [./vel_x]
    block = 2
  [../]
  [./vel_y]
    block = 2
  [../]
[]

# [ICs]
#   [./block2y]
#     block = 2
#     variable = disp_y
#     type = ConstantIC
#     value = -1
#   [../]
#   [./block2x]
#     block = 2
#     variable = disp_x
#     type = ConstantIC
#     value = 0
#   [../]
# []

[Kernels]
  [./disp_x]
    type = MatDiffusion
    variable = disp_x
  [../]
  [./disp_y]
    type = MatDiffusion
    variable = disp_y
  [../]
  [./accel_x]
    type = CoupledTimeDerivative
    variable = disp_x
    v = vel_x
    block = 2
  [../]
  [./accel_y]
    type = CoupledTimeDerivative
    variable = disp_y
    v = vel_y
    block = 2
  [../]
  [./coupled_time_velx]
    type = CoupledTimeDerivative
    variable = vel_x
    v = disp_x
    block = 2
  [../]
  [./coupled_time_vely]
    type = CoupledTimeDerivative
    variable = vel_y
    v = disp_y
    block = 2
  [../]
  [./source_velx]
    type = MatReaction
    variable = vel_x
    mob_name = 1
    block = 2
  [../]
  [./source_vely]
    type = MatReaction
    variable = vel_y
    mob_name = 1
    block = 2
  [../]
[]

[Constraints]
  [./lm]
    type = LMConstraint
    slave = 10
    master = 20
    variable = lm
    master_variable = disp_x
    disp_y = disp_y
  [../]
  [./tan_lm]
    type = TangentialLMConstraint
    slave = 10
    master = 20
    variable = tangent_lm
    contact_pressure = lm
    master_variable = vel_x
    vel_y = vel_y
    mu = 0.3
  [../]
[]

[BCs]
  [./botx]
    type = DirichletBC
    variable = disp_x
    boundary = 40
    value = 0.0
  [../]
  [./boty]
    type = DirichletBC
    variable = disp_y
    boundary = 40
    value = 0.0
  [../]
  [./topy]
    type = NeumannBC
    variable = disp_y
    boundary = 30
    value = -10
  [../]
  [./leftx]
    type = NeumannBC
    variable = disp_x
    boundary = 50
    value = 1
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 2
  dt = 1
  dtmin = 1
  solve_type = 'NEWTON'
  line_search = 'basic'
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_test_jacobian -snes_fd -pc_svd_monitor'# -snes_test_jacobian_view'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'svd'

  l_max_its = 100
  nl_max_its = 20
[]

[Outputs]
  exodus = true
[]

[Contact]
  [./leftright]
    master = 20
    slave = 10
    model = frictionless
    formulation = lagrange
    # penalty = 1e6
    system = constraint
    lm = lm
    tangent_lm = tangent_lm
    vel_x = vel_x
    vel_y = vel_y
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]
