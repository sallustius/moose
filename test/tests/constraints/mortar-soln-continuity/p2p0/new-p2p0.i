[Mesh]
  file = mortar-1.msh
[]

[MeshModifiers]
  [./master]
    type = LowerDBlockFromSideset
    sidesets = '106'
    new_block_id = '1000'
  [../]
  [./slave]
    type = LowerDBlockFromSideset
    sidesets = '107'
    new_block_id = '1001'
  [../]
[]

[Problem]
  kernel_coverage_check = false
[]

[AuxVariables]
  [exact]
    order = SECOND
    block = 'domain'
  []
[]

[AuxKernels]
  [function]
    type = FunctionAux
    variable = exact
    function = exact_soln_primal
    block = 'domain'
  []
[]

[Variables]
  [./T]
    block = 'domain'
    order = SECOND
  [../]
  [./lambda]
    block = '1001'
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

[BCs]
  [./neumann]
    type = FunctionGradientNeumannBC
    exact_solution = exact_soln_primal
    variable = T
    boundary = 'bottom_left bottom_right right top_right top_left left'
  [../]
[]

[Kernels]
  [./conduction]
    type = Diffusion
    variable = T
    block = 'domain'
  [../]
  [./sink]
    type = Reaction
    variable = T
    block = 'domain'
  [../]
  [./forcing_function]
    type = BodyForce
    variable = T
    function = forcing_function
    block = 'domain'
  [../]
[]

[Functions]
  [./forcing_function]
    type = ParsedFunction
    value = '-12*x^2 -24*y^2 - 2*x + x^4 + 2*y^4 + x*y^2'
  [../]
  [./exact_soln_primal]
    type = ParsedFunction
    value = 'x^4 + 2*y^4 + x*y^2'
  [../]
  [exact_soln_lambda]
    type = ParsedFunction
    value = '4*x^3 + y^2'
  []
[]

[Debug]
  show_var_residual_norms = 1
[]

[Constraints]
  [./mortar]
    type = SolutionContinuityTest
    master_boundary_id = 106
    slave_boundary_id = 107
    master_subdomain_id = 1000
    slave_subdomain_id = 1001
    variable = T
    lm_variable = lambda
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
  type = Steady
  petsc_options_iname = '-pc_type -snes_linesearch_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu       basic                 NONZERO               1e-15'
[]

[Outputs]
  exodus = true
  csv = true
  [dofmap]
    type = DOFMap
    execute_on = 'initial'
  []
[]

[Postprocessors]
  [L2lambda]
    type = ElementL2Error
    variable = lambda
    function = exact_soln_lambda
    outputs = 'console csv'
    execute_on = 'timestep_end'
    block = '1001'
  []
  [L2u]
    type = ElementL2Error
    variable = T
    function = exact_soln_primal
    outputs = 'console csv'
    execute_on = 'timestep_end'
    block = 'domain'
  []
  [H1u]
    type = ElementH1Error
    variable = T
    function = exact_soln_primal
    outputs = 'console csv'
    execute_on = 'timestep_end'
    block = 'domain'
  []
[]
