[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  # scaling = 1e-6
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = SMALL
    add_variables = false
    use_displaced_mesh = false
  [../]
[]

[BCs]
  [./symmy]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./symmx]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./tdisp]
    type = DirichletBC
    variable = disp_y
    boundary = top
    value = 0.1
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1
    poissons_ratio = 0
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  petsc_options = '-pc_svd_monitor -snes_test_jacobian -snes_test_jacobian_view -snes_converged_reason -ksp_converged_reason'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'svd'
  type = Transient
  dt = 0.05

  solve_type = 'NEWTON'

  dtmin = 0.05
  num_steps = 1
[]

[Outputs]
  exodus = true
[]
