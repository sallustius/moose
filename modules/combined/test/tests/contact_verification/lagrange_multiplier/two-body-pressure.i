[GlobalParams]
  displacements = 'disp_x disp_y'
  D_name = 1e0
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
  [./left_top_block_x]
    type = DirichletBC
    variable = disp_x
    boundary = 50
    value = 0
  [../]
  [./left_top_block_y]
    type = DirichletBC
    variable = disp_y
    boundary = 50
    value = 0
  [../]

  [./topx]
    type = DirichletBC
    variable = disp_x
    boundary = 30
    value = 0
  [../]
  # [./topy]
  #   type = DirichletBC
  #   variable = disp_y
  #   boundary = 30
  #   value = -0.004
  # [../]
  [./topy]
    type = Pressure
    variable = disp_y
    boundary = 30
    factor = 1e0
    component = 1
    use_displaced_mesh = false
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1
  dtmin = 1
  solve_type = 'NEWTON'
  line_search = 'bt'
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_test_jacobian'# -snes_test_jacobian_view'
  petsc_options_iname = '-pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'NONZERO	       1e-15'

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
