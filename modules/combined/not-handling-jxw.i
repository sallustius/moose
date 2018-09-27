[GlobalParams]
  displacements = 'disp_x'
[]

[Mesh]
  type = GeneratedMesh
  nx = 2
  dim = 1
[]

[Variables]
  [./disp_x]
  [../]
  [./u]
  [../]
[]

[Kernels]
  [./disp_x]
    type = Diffusion
    variable = disp_x
  [../]
  [./u]
    type = Diffusion
    variable = u
    use_displaced_mesh = true
    disp_x = disp_x
  [../]
[]

[BCs]
  [./left_u]
    boundary = left
    variable = u
    type = DirichletBC
    value = 0
  [../]
  [./right_u]
    boundary = right
    variable = u
    type = DirichletBC
    value = 1
  [../]
  [./disp_x_left]
    boundary = left
    variable = disp_x
    type = DirichletBC
    value = 0
  [../]
  [./disp_x_right]
    boundary = right
    variable = disp_x
    type = FunctionDirichletBC
    function = 't'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  solve_type = NEWTON
  type = Transient
  num_steps = 5
  petsc_options = '-snes_test_jacobian -snes_test_jacobian_view'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  exodus = true
[]
