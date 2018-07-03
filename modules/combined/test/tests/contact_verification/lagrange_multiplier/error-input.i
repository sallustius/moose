[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./lm]
  [../]
[]

[ICs]
  [./block2y]
    variable = disp_y
    type = ConstantIC
    value = -1
  [../]
  [./block2x]
    variable = disp_x
    type = ConstantIC
    value = 0
  [../]
[]

[Kernels]
  [./disp_x]
    type = Diffusion
    variable = disp_x
  [../]
  [./disp_y]
    type = Diffusion
    variable = disp_y
  [../]
[]

[Constraints]
  [./lm]
    type = LMConstraint
    slave = 0
    master = 0
    variable = lm
    master_variable = disp_x
    disp_y = disp_y
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1
  dtmin = 1
  solve_type = 'NEWTON'
  line_search = 'bt'
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_test_jacobian -snes_test_jacobian_view'
  petsc_options_iname = '-pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'NONZERO	       1e-15'

  l_max_its = 100
  nl_max_its = 100
  l_tol = 1e-3
[]

[Outputs]
  exodus = true
[]

[Contact]
  [./leftright]
    master = 0
    slave = 0
    model = frictionless
    formulation = lagrange
    # penalty = 1e6
    system = constraint
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
