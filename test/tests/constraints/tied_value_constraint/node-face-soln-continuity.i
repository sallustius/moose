[Mesh]
  file = node-face-6.msh
[]

[Variables]
  [./T]
  [../]
[]

[BCs]
  [./neumann]
    type = FunctionDirichletBC
    function = exact
    variable = T
    boundary = 'left right top_left top_right bottom_left bottom_right'
  [../]
[]

[Functions]
  [exact]
    type = ParsedFunction
    value = 'sin(x) * cos(y)'
  []
[]

[Kernels]
  [./conduction]
    type = Diffusion
    variable = T
  [../]
  [./forcing_function]
    type = BodyForce
    variable = T
    function = '2*sin(x)*cos(y)'
  [../]
[]

[Constraints]
  [./mortar]
    type = TiedValueConstraint
    slave = 'left_interface'
    master = 'right_interface'
    variable = T
    master_variable = T
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  [l2]
    type = ElementL2Error
    variable = T
    function = exact
  []
[]

[Executioner]
  solve_type = NEWTON
  type = Steady
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre    boomeramg'
[]

[Outputs]
  exodus = true
[]
