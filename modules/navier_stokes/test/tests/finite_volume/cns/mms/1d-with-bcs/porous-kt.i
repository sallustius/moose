[GlobalParams]
  fp = fp
  limiter = 'upwind'
  two_term_boundary_expansion = false
[]

[Mesh]
  [cartesian]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = .1
    xmax = .6
    nx = 2
  []
[]

[Modules]
  [FluidProperties]
    [fp]
      type = IdealGasFluidProperties
    []
  []
[]

[Problem]
  kernel_coverage_check = false
  fv_bcs_integrity_check = false
[]

[Variables]
  [pressure]
    type = MooseVariableFVReal
  []
  [ud]
    type = MooseVariableFVReal
  []
  [temperature]
    type = MooseVariableFVReal
  []
[]

[ICs]
  [pressure]
    type = FunctionIC
    variable = pressure
    function = 'exact_p'
  []
  [ud]
    type = FunctionIC
    variable = ud
    function = 'exact_ud'
  []
  [temperature]
    type = FunctionIC
    variable = temperature
    function = 'exact_T'
  []
[]

[FVKernels]
  [mass_advection]
    type = PCNSFVLaxFriedrichs
    variable = pressure
    eqn = "mass"
  []
  [mass_fn]
    type = FVBodyForce
    variable = pressure
    function = 'forcing_p'
  []

  [momentum_x_advection]
    type = PCNSFVLaxFriedrichs
    variable = ud
    momentum_component = x
    eqn = "momentum"
  []
  [eps_grad]
    type = PNSFVPGradEpsilon
    variable = ud
    momentum_component = 'x'
    epsilon_function = 'eps'
  []
  [momentum_fn]
    type = FVBodyForce
    variable = ud
    function = 'forcing_ud'
  []

  [fluid_energy_advection]
    type = PCNSFVLaxFriedrichs
    variable = temperature
    eqn = "energy"
  []
  [energy_fn]
    type = FVBodyForce
    variable = temperature
    function = 'forcing_T'
  []
[]

[FVBCs]
  [mass_left]
    variable = pressure
    type = PCNSFVLaxFriedrichsBC
    boundary = left
    T_fluid = 'exact_T'
    superficial_velocity = 'vector_ud'
    eqn = 'mass'
  []
  [momentum_left]
    variable = ud
    type = PCNSFVLaxFriedrichsBC
    boundary = left
    T_fluid = 'exact_T'
    superficial_velocity = 'vector_ud'
    eqn = 'momentum'
    momentum_component = 'x'
  []
  [energy_left]
    variable = ud
    type = PCNSFVLaxFriedrichsBC
    boundary = left
    T_fluid = 'exact_T'
    superficial_velocity = 'vector_ud'
    eqn = 'energy'
  []
  [mass_right]
    variable = pressure
    type = PCNSFVLaxFriedrichsBC
    boundary = right
    pressure = 'exact_p'
    eqn = 'mass'
  []
  [momentum_right]
    variable = ud
    type = PCNSFVLaxFriedrichsBC
    boundary = right
    pressure = 'exact_p'
    eqn = 'momentum'
    momentum_component = 'x'
  []
  [energy_right]
    variable = temperature
    type = PCNSFVLaxFriedrichsBC
    boundary = right
    pressure = 'exact_p'
    eqn = 'energy'
  []

  # Use these to help create more accurate cell centered gradients for cells adjacent to boundaries
  [T_left]
    type = FVFunctionDirichletBC
    variable = temperature
    function = 'exact_T'
    boundary = 'left'
  []
  [sup_vel_left]
    type = FVFunctionDirichletBC
    variable = ud
    function = 'exact_ud'
    boundary = 'left'
  []
  [p_right]
    type = FVFunctionDirichletBC
    variable = pressure
    function = 'exact_p'
    boundary = 'right'
  []
[]

[Materials]
  [var_mat]
    type = PorousPrimitiveVarMaterial
    pressure = pressure
    superficial_vel_x = ud
    T_fluid = temperature
    porosity = porosity
  []
  [porosity]
    type = GenericFunctionMaterial
    prop_names = 'porosity'
    prop_values = 'eps'
  []
[]

[Functions]
[exact_p]
  type = ParsedFunction
  value = '101000.0*cos(x)'
[]
[forcing_p]
  type = ParsedFunction
  value = '-1.4186482097536*sin(x)*cos(1.2*x)/cos(1.1*x) + 1.56051303072896*sin(1.1*x)*cos(x)*cos(1.2*x)/cos(1.1*x)^2 - 1.70237785170432*sin(1.2*x)*cos(x)/cos(1.1*x)'
[]
[exact_ud]
  type = ParsedFunction
  value = '1.1*cos(1.2*x)'
[]
[forcing_ud]
  type = ParsedFunction
  value = '-101000.0*sin(x)*cos(1.3*x) - 1.56051303072896*sin(x)*cos(1.2*x)^2/(cos(1.1*x)*cos(1.3*x)) + 1.71656433380186*sin(1.1*x)*cos(x)*cos(1.2*x)^2/(cos(1.1*x)^2*cos(1.3*x)) - 3.74523127374952*sin(1.2*x)*cos(x)*cos(1.2*x)/(cos(1.1*x)*cos(1.3*x)) + 2.02866693994765*sin(1.3*x)*cos(x)*cos(1.2*x)^2/(cos(1.1*x)*cos(1.3*x)^2)'
[]
[exact_T]
  type = ParsedFunction
  value = '273.15*cos(1.1*x)'
[]
[forcing_T]
  type = ParsedFunction
  value = '-1.4186482097536*(274098.960775862*cos(1.1*x) + 0.605*cos(1.2*x)^2/cos(1.3*x)^2)*sin(x)*cos(1.2*x)/cos(1.1*x) + 1.56051303072896*(274098.960775862*cos(1.1*x) + 0.605*cos(1.2*x)^2/cos(1.3*x)^2)*sin(1.1*x)*cos(x)*cos(1.2*x)/cos(1.1*x)^2 - 1.70237785170432*(274098.960775862*cos(1.1*x) + 0.605*cos(1.2*x)^2/cos(1.3*x)^2)*sin(1.2*x)*cos(x)/cos(1.1*x) + 1.4186482097536*(-301508.856853448*sin(1.1*x) - 1.452*sin(1.2*x)*cos(1.2*x)/cos(1.3*x)^2 + 1.573*sin(1.3*x)*cos(1.2*x)^2/cos(1.3*x)^3)*cos(x)*cos(1.2*x)/cos(1.1*x)'
[]
[vector_ud]
  type = ParsedVectorFunction
  value_x = '1.1*cos(1.2*x)'
[]
[eps]
  type = ParsedFunction
  value = 'cos(1.3*x)'
[]
[]

[Executioner]
  solve_type = NEWTON
  type = Transient
  num_steps = 1
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_max_its = 50
  line_search = bt
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
  dtmin = 1
[]

[Outputs]
  exodus = true
  csv = true
[]

[Debug]
  show_var_residual_norms = true
[]

[Postprocessors]
  [h]
    type = AverageElementSize
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
  [L2pressure]
    type = ElementL2Error
    variable = pressure
    function = exact_p
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
  [L2ud]
    variable = ud
    function = exact_ud
    type = ElementL2Error
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
  [L2temperature]
    variable = temperature
    function = exact_T
    type = ElementL2Error
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
[]
