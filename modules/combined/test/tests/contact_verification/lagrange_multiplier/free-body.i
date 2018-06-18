[GlobalParams]
  displacements = 'disp_x disp_y'
  D_name = 1e6
  use_displaced_mesh = true
[]

[Mesh]
  file = one-body.e
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./vel_x]
  [../]
  [./vel_y]
  [../]
[]

[Kernels]
  [./disp_x_diffusion]
    type = MatDiffusion
    variable = disp_x
  [../]
  [./disp_y_diffusion]
    type = MatDiffusion
    variable = disp_y
  [../]
  [./accel_x]
    type = CoupledTimeDerivative
    variable = disp_x
    v = vel_x
  [../]
  [./accel_y]
    type = CoupledTimeDerivative
    variable = disp_y
    v = vel_y
  [../]
  [./coupled_time_velx]
    type = CoupledTimeDerivative
    variable = vel_x
    v = disp_x
  [../]
  [./coupled_time_vely]
    type = CoupledTimeDerivative
    variable = vel_y
    v = disp_y
  [../]
  [./source_velx]
    type = MatReaction
    variable = vel_x
    mob_name = 1
  [../]
  [./source_vely]
    type = MatReaction
    variable = vel_y
    mob_name = 1
  [../]
[]


[BCs]
  # [./topy]
  #   type = Pressure
  #   variable = disp_y
  #   boundary = 30
  #   factor = 1e0
  #   component = 1
  #   use_displaced_mesh = true
  # [../]
  [./topy]
    type = NeumannBC
    variable = disp_y
    boundary = 30
    value = -1
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 4
  dt = 1
  dtmin = 1
  solve_type = 'NEWTON'
  line_search = 'bt'
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_test_jacobian'# -snes_test_jacobian_view'
  petsc_options_iname = '-pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'NONZERO	       1e-15'
  # scheme = 'explicit-euler'

  l_max_its = 100
  nl_max_its = 20
[]

[Outputs]
  exodus = true
  checkpoint = true
  dofmap = true
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
