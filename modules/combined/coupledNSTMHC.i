# Attempt to couple Navier Stokes with thermal conduction and tensor mechanics.
# Using the bump.i test as a basis for this

[Mesh]
  [annular]
  type = AnnularMeshGenerator
  nt = 10
  rmax = 0.11
  rmin = 0.025
  nr = 10
  []
#  [./assignBase]
#  type = ParsedGenerateSideset
#  combinatorial_geometry = 'y < -0.105'
#  input = annular
#  new_sideset_name = base
#  replace = true
#  [../]
  [splitDomain]
  type = ParsedSubdomainMeshGenerator
  block_id = 1
  combinatorial_geometry = '(x*x + y*y) > 0.1^2'
  input = annular
  []
  [setOuter]
  type = ParsedGenerateSideset
  input = splitDomain
  combinatorial_geometry = '(x*x * y*y) > 0.109'
  new_sideset_name = 'outer'
  []
  [interfaceSideset]
  type = SideSetsBetweenSubdomainsGenerator
  input = setOuter
  new_boundary = midline
  paired_block = 1
  master_block = 0
  []
[]

[Variables]
  [./T]
    initial_condition = 300
  [../]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[]

[Modules]
  [./FluidProperties]
    [./ideal_gas]
      type = IdealGasFluidProperties
      gamma = 1.4
      #R = 287
    [../]
  [../]

  [./CompressibleNavierStokes]
    total_energy_scaling = '9.869232667160121e-6'
    family = LAGRANGE
    order = FIRST
    block = 0
    initial_pressure = 101325.
    initial_temperature = 300.
    initial_velocity = '17.3594354746921 0 0' # Mach 0.5: = 0.5*sqrt(gamma*R*T)
    fluid_properties = ideal_gas

    stagnation_boundary = 'rmin'
    stagnation_pressure = 120192.995549849 # Pa, Mach=0.5 at 1 atm
    stagnation_temperature = 315 # K, Mach=0.5 at 1 atm
    stagnation_flow_direction = '1 0'

    no_penetration_boundary = 'midline'
  [../]
[]



[Materials]
  [./fluid]
    type = Air
    block = 0
    rho = rho
    rhou = rhou
    rhov = rhov
    rhoE = rhoE
    vel_x = vel_x
    vel_y = vel_y
    temperature = temperature
    enthalpy = enthalpy
    # This value is not used in the Euler equations, but it *is* used
    # by the stabilization parameter computation, which it decreases
    # the amount of artificial viscosity added, so it's best to use a
    # realistic value.
    dynamic_viscosity = 0.0
    fluid_properties = ideal_gas
  [../]
  [./blockMaterial]
    type = HeatConductionMaterial
    specific_heat = 2250
    thermal_conductivity = 0.47
    block=1
 [../]
 [./density_]
   type = GenericConstantMaterial
   prop_names = density
   prop_values = 2250
   block=1
 [../]
 [./Elasticity_tensor]
   type = ComputeIsotropicElasticityTensor
   youngs_modulus = 500E9
   poissons_ratio = 0.42
   block=1
 [../]
 [./strain]
   type = ComputeSmallStrain
   displacements = 'disp_x disp_y'
   eigenstrain_names = eigenstrain
   block=1
 [../]
 [./stress]
   type = ComputeLinearElasticStress
   block = 1
 [../]
 [./thermal_strain]
   type= ComputeThermalExpansionEigenstrain
   thermal_expansion_coeff = 12E-5
   temperature = T
   stress_free_temperature = 293
   eigenstrain_name = eigenstrain
   block = 1
 [../]
[]

[Kernels]
  [./HeatConduction]
    type = HeatConduction
    variable = T
    block = 1
  [../]
  [./HeatConductionTimeDerivative]
    type = HeatConductionTimeDerivative
    variable=T
    block = 1
  [../]
  [./TensorMechanics]
    displacements = 'disp_x disp_y'
    block = 1
  [../]
#  [./gravity_y]
#    #Gravity is applied to the balloon
#    type = Gravity
#    variable = disp_y
#    value = -9.81
#  [../]
[]


#[Preconditioning]
#  [./SMP_PJFNK]
#    type = SMP
#    full = true
#    solve_type = 'PJFNK'
#  [../]
#[]



[Executioner]
  type = Transient
  dt = 5.e-5
  dtmin = 1.e-5
  start_time = 0.0
  num_steps = 10
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-9
  # nl_abs_step_tol = 1e-15
  nl_max_its = 5
  l_tol = 1e-4 # Relative linear tolerance for each Krylov solve
  l_max_its = 100 # Number of linear iterations for each Krylov solve

  # Specify the order as FIRST, otherwise you will get warnings in DEBUG mode...
  [./Quadrature]
    type = TRAP
    order = FIRST
  [../]
[]



[Outputs]
  file_base = step_out
  interval = 1
  exodus = true
  [dof]
    type = DOFMap
    execute_on = 'initial'
  []
[]
