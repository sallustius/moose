mu=1.1
rho=1.1

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = -0.6
    xmax = 0.6
    # ymin = -0.6
    # ymax = 0.6
    nx = 2
    # ny = 2
  []
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [u]
    order = CONSTANT
    family = MONOMIAL
    fv = true
    initial_condition = 1
  []
  # [v]
  #   order = CONSTANT
  #   family = MONOMIAL
  #   fv = true
  #   initial_condition = 1
  # []
  [pressure]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
[]

# [AuxVariables]
#   [U]
#     order = CONSTANT
#     family = MONOMIAL
#     fv = true
#   []
# []

# [AuxKernels]
#   [mag]
#     type = VectorMagnitudeAux
#     variable = U
#     x = u
#     y = v
#   []
# []

[FVKernels]
  [mass]
    # type = NSFVKernel
    type = FVMatAdvection
    variable = pressure
    advected_quantity = 1
    advected_interp_method = 'average'
    vel = 'velocity'
    # pressure = pressure
    # u = v
    # v = v
    # mu = ${mu}
    # rho = ${rho}
  []
  [mass_forcing]
    type = FVBodyForce
    variable = pressure
    function = forcing_p
  []

  # [u_advection]
  #   type = NSFVKernel
  #   variable = u
  #   advected_quantity = 'rhou'
  #   pressure = pressure
  #   u = u
  #   v = v
  #   mu = ${mu}
  #   rho = ${rho}
  # []
  [u_advection]
    type = FVMatAdvection
    variable = u
    advected_quantity = 'rhou'
    vel = 'velocity'
    advected_interp_method = 'average'
  []
  [u_viscosity]
    type = FVDiffusion
    variable = u
    coeff = ${mu}
  []
  [u_pressure]
    # FVMomPressure inherits from FVMatAdvection and in FVMomPressure::validParams we set
    # 'advected_quantity = NS::prssure'
    type = FVMomPressure
    variable = u
    momentum_component = 'x'

    # these parameters shouldn't be used for anything but are still required
    vel = 'velocity'
    advected_interp_method = 'average'
  []
  [u_forcing]
    type = FVBodyForce
    variable = u
    function = forcing_u
  []

  # [v_advection]
  #   type = NSFVKernel
  #   variable = v
  #   advected_quantity = 'rhov'
  #   pressure = pressure
  #   u = u
  #   v = v
  #   mu = ${mu}
  #   rho = ${rho}
  # []
  # [v_viscosity]
  #   type = FVDiffusion
  #   variable = v
  #   coeff = ${mu}
  # []
  # [v_pressure]
  #   type = FVMomPressure
  #   variable = v
  #   momentum_component = 'y'
  # []
  # [v_forcing]
  #   type = FVBodyForce
  #   variable = v
  #   function = forcing_v
  # []
[]

[FVBCs]
  [u_advection]
    type = FVMatAdvectionFunctionBC
    # boundary = 'left right top bottom'
    boundary = 'left right'
    variable = u
    vel = 'velocity'
    flux_variable_exact_solution = 'exact_rhou'
    advected_quantity = 'rhou'
    advected_interp_method = 'average'
    vel_x_exact_solution = 'exact_u'
    # vel_y_exact_solution = 'exact_v'
  []
  [u_diffusion]
    type = FVDiffusionFunctionBC
    # boundary = 'left right top bottom'
    boundary = 'left right'
    variable = u
    exact_solution = 'exact_u'
    coeff = '${mu}'
    coeff_function = '${mu}'
  []
  [u_pressure]
    type = FVMomPressureFunctionBC
    # boundary = 'left right top bottom'
    boundary = 'left right'
    variable = u
    momentum_component = 'x'
    p = pressure
    pressure_exact_solution = 'exact_p'
  []

  # [v_advection]
  #   type = FVMatAdvectionFunctionBC
  #   boundary = 'left right top bottom'
  #   variable = v
  #   vel = 'velocity'
  #   flux_variable_exact_solution = 'exact_rhov'
  #   vel_x_exact_solution = 'exact_u'
  #   vel_y_exact_solution = 'exact_v'
  # []
  # [v_diffusion]
  #   type = FVDiffusionFunctionBC
  #   boundary = 'left right top bottom'
  #   variable = v
  #   exact_solution = 'exact_v'
  #   coeff = '${mu}'
  #   coeff_function = '${mu}'
  # []
  # [v_pressure]
  #   type = FVMomPressureFunctionBC
  #   boundary = 'left right top bottom'
  #   variable = v
  #   momentum_component = 'y'
  #   p = pressure
  #   pressure_exact_solution = 'exact_p'
  # []

  [mass_continuity_flux]
    type = FVMatAdvectionFunctionBC
    variable = pressure
    # boundary = 'top bottom left right'
    boundary = 'left right'
    advected_quantity = 1
    vel = 'velocity'
    flux_variable_exact_solution = 1
    vel_x_exact_solution = 'exact_u'
    # vel_y_exact_solution = 'exact_v'
  []
[]

[Materials]
  [rho]
    type = ADGenericConstantMaterial
    prop_names = 'rho'
    prop_values = ${rho}
  []
  [ins_fv]
    type = INSFVMaterial
    u = 'u'
    v = '1'
    # we need to compute this here for advection in FVMomPressure
    pressure = 'pressure'
  []
[]

[Functions]
[exact_u]
  type = ParsedFunction
  value = '1.1*sin(1.1*x)'
[]
[exact_rhou]
  type = ParsedFunction
  value = '1.1*rho*sin(1.1*x)'
  vars = 'rho'
  vals = '${rho}'
[]
[forcing_u]
  type = ParsedFunction
  value = '1.331*mu*sin(1.1*x) + 2.662*rho*sin(1.1*x)*cos(1.1*x) - 1.21*sin(1.1*x)'
  vars = 'mu rho'
  vals = '${mu} ${rho}'
[]
[exact_p]
  type = ParsedFunction
  value = '1.1*cos(1.1*x)'
[]
[forcing_p]
  type = ParsedFunction
  value = '1.21*cos(1.1*x)'
[]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      100                lu           NONZERO'
[]

[Outputs]
  exodus = true
  csv = true
  [dof]
    type = DOFMap
    execute_on = 'initial'
  []
[]

[Postprocessors]
  [h]
    type = AverageElementSize
    outputs = 'console csv'
    execute_on = 'timestep_end'
  []
  [./L2u]
    type = ElementL2Error
    variable = u
    function = exact_u
    outputs = 'console csv'
    execute_on = 'timestep_end'
  [../]
  # [./L2v]
  #   variable = v
  #   function = exact_v
  #   type = ElementL2Error
  #   outputs = 'console csv'
  #   execute_on = 'timestep_end'
  # [../]
  [./L2p]
    variable = pressure
    function = exact_p
    type = ElementL2Error
    outputs = 'console csv'
    execute_on = 'timestep_end'
  [../]
[]
