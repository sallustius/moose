[Mesh]
  [file]
    file = node-face-6.msh
    type = FileMeshGenerator
  []
  [master]
    input = file
    type = LowerDBlockFromSidesetGenerator
    sidesets = 'right_interface'
    new_block_id = '200'
  []
  [slave]
    input = master
    type = LowerDBlockFromSidesetGenerator
    sidesets = 'left_interface'
    new_block_id = '100'
  []
[]

[Variables]
  [./T]
    block = domain
  [../]
  [lambda]
    block = 100
  []
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
  [flux]
    type = ParsedFunction
    value = 'cos(x) * cos(y)'
  []
[]

[Kernels]
  [./conduction]
    type = Diffusion
    variable = T
    block = 'domain'
  [../]
  [./forcing_function]
    type = BodyForce
    variable = T
    function = '2*sin(x)*cos(y)'
    block = 'domain'
  [../]
[]

[Constraints]
  [./mortar]
    type = EqualValueConstraint
    slave_boundary = 'left_interface'
    master_boundary = 'right_interface'
    slave_subdomain = 100
    master_subdomain = 200
    variable = lambda
    slave_variable = T
    compute_lm_residuals = false
  [../]
  [lm-nf]
    type = LMTiedValueConstraint
    variable = lambda
    master_variable = T
    primal_var = T
    slave = 'left_interface'
    master = 'right_interface'
  []
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
    block = domain
  []
  [h1]
    type = ElementH1Error
    variable = T
    function = exact
    block = domain
  []
  [lambda_l2]
    type = ElementL2Error
    variable = lambda
    function = flux
    block = 100
  []
[]

[Executioner]
  solve_type = PJFNK
  type = Steady
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
[]

[Outputs]
  exodus = true
  [dof]
    type = DOFMap
    execute_on = 'initial'
  []
[]
