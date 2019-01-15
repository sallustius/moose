[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  type = GeneratedMesh
  nx = 10
  ny = 5
  dim = 2
[]

[Variables]
  [./disp_x]
  [../]
  [./u]
  [../]
[]

[AuxVariables]
  [./disp_y]
  [../]
[]

[Kernels]
  [./disp_x]
    type = Diffusion
    variable = disp_x
  [../]
[]

[ADKernels]
  [./u]
    type = WeirdTestKernel
    variable = u
    use_displaced_mesh = true
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
  num_steps = 2
  dtmin = 1
[]

[ICs]
  [./disp_x]
    type = FunctionIC
    variable = disp_x
    function = '1 + x'
  [../]
  [./u]
    type = FunctionIC
    variable = u
    function = '2 + 2 * x'
  [../]
[]

[Outputs]
  exodus = true
[]
