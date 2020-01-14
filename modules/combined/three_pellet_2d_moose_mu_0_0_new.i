[GlobalParams]
  displacements = 'disp_x disp_y'
  order = FIRST
  family = LAGRANGE
  temperature = temp
[]

[Problem]
  coord_type = RZ
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
  group_variables = 'disp_x disp_y'
[]

[Mesh]
  file = three_pellet_2D_smeared.e
  patch_size = 5
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./temp]
    initial_condition = 298
  [../]
[]

[AuxVariables]
  [./gap_conductance]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./coolant_press_ramp]
    type = PiecewiseLinear
    x = '0 1.0e5 5.1e5'
    y = '0  1.0  1.0'
  [../]
  [./coolant_temp_ramp]
    type = PiecewiseLinear
    x = '0 1.0e5 5.0e5 5.1e5'
    y = '300.0  680.0  680.0 300.0'
  [../]
  [./volumetric_heat]
    type = PiecewiseLinear
    x = '0 1.0e5 2.0e5 3.0e5 4.0e5 5.0e5 5.1e5'
    y = '0  5.0e8 1.0e9 2.0e9 5.0e9 0  0'
  [../]
[]

[Modules]
  [./TensorMechanics]
    [./Master]
      [./fuel]
        strain = FINITE
        block = '3'
        add_variables = true
        incremental = true
        eigenstrain_names = fuel_eigenstrain
        generate_output = 'stress_xx stress_yy vonmises_stress elastic_strain_yy'
        extra_vector_tags = 'ref'
      [../]
      [./clad]
        strain = FINITE
        block = '1'
        add_variables = true
        incremental = true
        generate_output = 'stress_xx stress_yy vonmises_stress elastic_strain_yy creep_strain_yy'
        extra_vector_tags = 'ref'
      [../]
    [../]
  [../]
[]

[Kernels]
  [./heat]
    type = HeatConduction
    variable = temp
    extra_vector_tags = 'ref'
  [../]
  [./heat_ie]
    type = HeatConductionTimeDerivative
    variable = temp
    extra_vector_tags = 'ref'
  [../]
  [./heat_source_fuel]
    type = HeatSource
    variable = temp
    block = '3'
    function = volumetric_heat
    extra_vector_tags = 'ref'
  [../]
[]

[AuxKernels]
  [./gap_cond]
    type = MaterialRealAux
    property = gap_conductance
    variable = gap_conductance
    boundary = 10
    execute_on = linear
  [../]
[]

[Contact]
  [./pellet_clad_mechanical]
    master = 5
    slave = 10
    system = Constraint
    penalty = 1e7
    model = frictionless
    formulation = kinematic
    normal_smoothing_distance = 0.1
  [../]
[]

[ThermalContact]
  [./pellet_clad_thermal]
    type = GapHeatTransfer
    variable = temp
    master = 5
    slave = 10
    gap_conductivity = 1.0
    quadrature = true
  [../]
[]

[BCs]
  [./no_x_all]
    type = PresetBC
    variable = disp_x
    boundary = 12
    value = 0.0
  [../]
# pin entire clad bottom in y
  [./no_y_clad_bottom]
    type = PresetBC
    variable = disp_y
    boundary = 1
    value = 0.0
  [../]
# pin fuel bottom in y
  [./no_y_fuel_bottom]
    type = PresetBC
    variable = disp_y
    boundary = '1020'
    value = 0.0
  [../]
# wedge bcs
  [./clad_temp]
    type = FunctionPresetBC
    variable = temp
    boundary = '1 2 3'
    function = coolant_temp_ramp
  [../]
  [./Pressure]
    # apply coolant pressure on clad outer walls
    [./coolantPressure]
      boundary = '1 2 3'
      factor = 15.5e6
      function = coolant_press_ramp
    [../]
  [../]
[]

[Materials]
  [./fuel_thermal]
    type = HeatConductionMaterial
    block = '3'
    thermal_conductivity = 20
    specific_heat = 500.0
  [../]
  [./fuel_elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    block = '3'
    youngs_modulus = 2e11
    poissons_ratio = 0.345
  [../]
  [./fuel_thermal_strain]
    type = ComputeThermalExpansionEigenstrain
    block = '3'
    stress_free_temperature = 298
    thermal_expansion_coeff = 50.0e-6
    eigenstrain_name = fuel_eigenstrain
  [../]
  [./fuel_stress]
    type = ComputeFiniteStrainElasticStress
    block = '3'
  [../]
  [./radial_return_stress]
    type = ComputeMultipleInelasticStress
    block = 1
    inelastic_models = 'power_law_creep'
    tangent_operator = elastic
  [../]
  [./power_law_creep]
    type = PowerLawCreepStressUpdate
    block = 1
    coefficient = 1.0e-15
    n_exponent = 1
    activation_energy = 0
  [../]
  [./clad_elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    block = 1
    youngs_modulus = 7.5e10
    poissons_ratio = 0.3
  [../]
  [./clad_thermal]
    type = HeatConductionMaterial
    block = 1
    thermal_conductivity = 16.0
    specific_heat = 330.0
  [../]
  [./clad_density]
    type = Density
    block = 1
    density = 6551.0
  [../]
  [./fuel_density]
    type = Density
    block = '3'
    density = 10000.0
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options = '-ksp_monitor_true_residual -pc_svd_monitor'
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  # petsc_options_value = 'lu       superlu_dist'
  petsc_options_iname = '-pc_type -ksp_gmres_restart'
  petsc_options_value = 'svd      200'

  line_search = 'none'

  l_max_its = 200
  l_tol = 1e-3
  nl_max_its = 10
  nl_rel_tol = 1e-4

  dtmax = 1.0e5
  dtmin = 1e3
  end_time = 2.0e5

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e3
    optimal_iterations = 7
    iteration_window = 2
    time_t ='0   5.0e4 1.0e5 5.1e5'
    time_dt ='2e3  1.2e4  7.5e3   7.5e3'
  [../]

  [./Quadrature]
     order = fifth
     side_order = seventh
  [../]

  verbose = true
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  [./ave_temp_interior]
    type = SideAverageValue
    boundary = 9
    variable = temp
    execute_on = 'initial linear'
  [../]
  [./clad_inner_vol]
   type = InternalVolume
    boundary = 7
    execute_on = 'initial timestep_end'
  [../]
  [./pellet_volume]
    type = InternalVolume
    boundary = 8
    execute_on = 'initial timestep_end'
  [../]
  [./gas_volume]
    type = InternalVolume
    boundary = 9
    execute_on = 'initial linear'
  [../]
  [./avg_clad_temp]
    type = SideAverageValue
    boundary = 7 #bws was 9
    variable = temp
    execute_on = 'initial timestep_end'
  [../]
  [./flux_from_clad]           # area integrated heat flux from the cladding
    type = SideFluxIntegral
    variable = temp
    boundary = 5
    diffusivity = thermal_conductivity
    execute_on = linear
  [../]
  [./flux_from_fuel]          # area integrated heat flux from the fuel
    type = SideFluxIntegral
    variable = temp
    boundary = 10
    diffusivity = thermal_conductivity
    execute_on = linear
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./residual]
    type = Residual
  [../]
  [./nl_its]
    type = NumNonlinearIterations
  [../]
  [./lin_its]
    type = NumLinearIterations
  [../]
  [cum_nl]
    type = CumulativeValuePostprocessor
    postprocessor = nl_its
  []
  [cum_lin]
    type = CumulativeValuePostprocessor
    postprocessor = lin_its
  []
  [./clad_elongation]
    type = NodalMaxValue
    variable = disp_y
    boundary = 2
  [../]
[]

[Outputs]
  csv = true
  interval = 1
  exodus = true
  color = false
  print_linear_residuals = false
  [./console]
    type = Console
    perf_log = true
    solve_log = true
  [../]
[]
