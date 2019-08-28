[Mesh]
  type = GeneratedMesh
  zmax = 3.6
  xmin = -0.004
  ymin = -0.004
  ymax = 0.004
  nx = 20
  ny = 20
  nz = 150
  dim = 3
  xmax = 0.004
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Materials]
  [bar_elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 219e9
    poissons_ratio = 0.345
  []
  [bar_thermal_expansion_strain]
    type = ComputeThermalExpansionEigenstrain
    stress_free_temperature = '298'
    temperature = 'temp'
    thermal_expansion_coeff = 10e-6
    eigenstrain_name = eigenstrain
  []
  [bar_stress]
    type = ComputeFiniteStrainElasticStress
  []
[]

[Modules]
  [TensorMechanics]
    [Master]
      [all]
        generate_output = 'vonmises_stress'
        add_variables = true
        temperature = temp
        eigenstrain_names = 'eigenstrain'
        strain = FINITE
        decomposition_method = EigenSolution
        use_displaced_mesh = true
      []
    []
  []
[]

[BCs]
  [Pressure]
    [boundary_pressure]
      function = 5e6
      boundary = 'left right bottom top front'
    []
  []
  [bottom_disp_x]
    type = PresetBC
    variable = disp_x
    boundary = 'back'
    value = 0.0
  []
  [bottom_disp_y]
    type = PresetBC
    variable = disp_y
    boundary = 'back'
    value = 0.0
  []
  [bottom_disp_z]
    type = PresetBC
    variable = disp_z
    boundary = 'back'
    value = 0.0
  []
[]

[Functions]
  [temp_space_func]
    type = ParsedFunction
    value = '1000+50*(x/0.004)'
  []
  [temp_time_factor]
    type = PiecewiseLinear
    y = '0.53 1 1'
    x = '0 10800 201830400'
  []
  [temp_func]
    type = CompositeFunction
    functions = 'temp_space_func temp_time_factor'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  petsc_options = '-ksp_monitor_true_residual'
  petsc_options_iname = '-snes_type -snes_ls -ksp_gmres_restart -pc_type -pc_hypre_type'
  petsc_options_value = 'newtonls   basic    201                hypre    boomeramg'
  line_search = none
  dt = 4036608
  solve_type = PJFNK
  end_time = 100915200
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-4
  l_tol = 8e-3
  l_max_its = 100
  [TimeStepper]
    type = TimeSequenceStepper
    time_sequence = '540 1080 1620 2160 2700 3240 3780 4320 4860 5400 5940 6480 7020 7560 8100 8640 9180 9720 10260 10800 11340 11880 12420 12960'
  []
[]

[Outputs]
  exodus = true
  csv = true
  print_linear_residuals = false
  checkpoint = true
[]

[AuxVariables]
  [temp]
  []
[]

[AuxKernels]
  [temp_auxkernel]
    type = FunctionAux
    function = temp_func
    variable = temp
    execute_on = 'LINEAR INITIAL'
  []
[]
