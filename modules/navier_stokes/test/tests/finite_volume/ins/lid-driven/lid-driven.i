mu=.01
rho=1

[GlobalParams]
  vel = 'velocity'
  velocity_interp_method = 'rc'
  advected_interp_method = 'average'
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = .1
    ymin = 0
    ymax = .1
    nx = 20
    ny = 20
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
  []
  [v]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [pressure]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
[]

[AuxVariables]
  [U]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
[]

[AuxKernels]
  [mag]
    type = VectorMagnitudeAux
    variable = U
    x = u
    y = v
  []
[]

[FVKernels]
  [mass]
    type = INSFVMass
    variable = pressure
    constrain_pressure = true
    pressure = pressure
    u = u
    v = v
    mu = ${mu}
    rho = ${rho}
  []

  [u_advection]
    type = NSFVKernel
    variable = u
    advected_quantity = 'rhou'
    pressure = pressure
    u = u
    v = v
    mu = ${mu}
    rho = ${rho}
  []

  [u_viscosity]
    type = FVDiffusion
    variable = u
    coeff = ${mu}
  []

  [u_pressure]
    type = FVMomPressure
    variable = u
    momentum_component = 'x'
  []

  [v_advection]
    type = NSFVKernel
    variable = v
    advected_quantity = 'rhov'
    pressure = pressure
    u = u
    v = v
    mu = ${mu}
    rho = ${rho}
  []

  [v_viscosity]
    type = FVDiffusion
    variable = v
    coeff = ${mu}
  []

  [v_pressure]
    type = FVMomPressure
    variable = v
    momentum_component = 'y'
  []
[]

[FVBCs]
  [top_x]
    type = FVDirichletBC
    variable = u
    value = 1
    boundary = 'top'
  []

  [no_slip_x]
    type = FVDirichletBC
    variable = u
    value = 0
    boundary = 'left right bottom'
  []

  [no_slip_y]
    type = FVDirichletBC
    variable = v
    value = 0
    boundary = 'left right top bottom'
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
    v = 'v'
    pressure = 'pressure'
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
  [dof]
    type = DOFMap
    execute_on = 'initial'
  []
[]
