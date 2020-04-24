mu=1.5
rho=2.5

[GlobalParams]
  integrate_p_by_parts = true
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1.0
    ymin = 0
    ymax = 1.0
    elem_type = QUAD9
    nx = 2
    ny = 2
  []
[]

[Variables]
  [./velocity]
    family = LAGRANGE_VEC
    order = SECOND
  [../]
  [./p]
  [../]
[]

[Kernels]
  [./mass]
    type = INSADMass
    variable = p
  [../]
  [mass_forcing]
    type = BodyForce
    variable = p
    function =
  []

  [./x_momentum_space]
    type = INSADMomentumAdvection
    variable = velocity
  [../]
  [./momentum_viscous]
    type = INSADMomentumViscous
    variable = velocity
  [../]
  [./momentum_pressure]
    type = INSADMomentumPressure
    variable = velocity
    p = p
    integrate_p_by_parts = true
  [../]
  [./momentum_forcing]
    type = VectorBodyForce
    variable = velocity
    function_x =
    function_y =
  [../]
[]

[BCs]
  [./vel_x]
    type = VectorFunctionDirichletBC
    boundary = 'left right top bottom'
    function_x =
    function_y =
    variable = velocity
  [../]
  [./p]
    type = FunctionDirichletBC
    boundary = 'left right top bottom'
    function = p_func
    variable = p
  [../]
[]

[Functions]
  [./vel_x_source_func]
    type = ParsedFunction
    value =
  [../]
  [./vel_y_source_func]
    type = ParsedFunction
    value =
  [../]
  [./p_source_func]
    type = ParsedFunction
    value = '-0.06*x*pi*cos(0.3*x*y*pi) - 0.14*y*pi*cos(0.2*x*y*pi) - 0.2*pi*cos((1/2)*x*pi) - 0.09*pi*cos(0.3*y*pi)'
  [../]
  [./vel_x_func]
    type = ParsedFunction
    value = '0.4*sin(0.5*pi*x) + 0.4*sin(pi*y) + 0.7*sin(0.2*pi*x*y) + 0.5'
  [../]
  [./vel_y_func]
    type = ParsedFunction
    value = '0.6*sin(0.8*pi*x) + 0.3*sin(0.3*pi*y) + 0.2*sin(0.3*pi*x*y) + 0.3'
  [../]
  [./p_func]
    type = ParsedFunction
    value = '0.5*sin(0.5*pi*x) + 1.0*sin(0.3*pi*y) + 0.5*sin(0.2*pi*x*y) + 0.5'
  [../]
[]

[Materials]
  [./const]
    type = ADGenericConstantMaterial
    prop_names = 'rho mu'
    prop_values = '${rho}  ${mu}'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    solve_type = 'NEWTON'
  [../]
[]

[Executioner]
  type = Steady
  petsc_options = '-snes_converged_reason -ksp_converged_reason'
  line_search = 'none'
  nl_rel_tol = 1e-12
  nl_max_its = 6
[]

[Outputs]
  [./exodus]
    type = Exodus
  [../]
  [./csv]
    type = CSV
  [../]
[]

[Postprocessors]
  [./L2vel_x]
    type = ElementL2Error
    variable = vel_x
    function = vel_x_func
    execute_on = 'timestep_end'
  [../]
  [./L2vel_y]
    variable = vel_y
    function = vel_y_func
    type = ElementL2Error
    execute_on = 'timestep_end'
  [../]
  [./L2p]
    variable = p
    function = p_func
    type = ElementL2Error
    execute_on = 'timestep_end'
  [../]
  [./L2vxx]
    variable = vxx
    function = vxx_func
    type = ElementL2Error
    execute_on = 'timestep_end'
  [../]
[]
