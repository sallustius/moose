mu=1
rho=1
advected_interp_method='average'
velocity_interp_method='rc'

[Mesh]
  [inlet]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = -1
    xmax = 1
    ymin = -1
    ymax = 1
    zmin = -9
    zmax = -1
    nx = 4
    ny = 4
    nz = 16
    boundary_name_prefix = inlet
  []

  [elbow]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = -1
    xmax = 1
    ymin = -1
    ymax = 1
    zmin = -1
    zmax = 1
    nx = 4
    ny = 4
    nz = 4
    boundary_name_prefix = elbow
    boundary_id_offset = 6
  []

  [outlet]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 1
    xmax = 9
    ymin = -1
    ymax = 1
    zmin = -1
    zmax = 1
    nx = 16
    ny = 4
    nz = 4
    boundary_name_prefix = outlet
    boundary_id_offset = 12
  []

  [combiner_1]
    type = StitchedMeshGenerator
    inputs = 'inlet elbow'
    clear_stitched_boundary_ids = true
    stitch_boundaries_pairs = 'inlet_front elbow_back'
  []
  [combiner_2]
    type = StitchedMeshGenerator
    inputs = 'combiner_1 outlet'
    clear_stitched_boundary_ids = true
    stitch_boundaries_pairs = 'elbow_right outlet_left'
  []

  uniform_refine = 2
[]

[Problem]
  kernel_coverage_check = false
  fv_bcs_integrity_check = true
[]

[Variables]
  [u]
    type = INSFVVelocityVariable
    initial_condition = 1e-15
  []
  [v]
    type = INSFVVelocityVariable
    initial_condition = 1e-15
  []
  [w]
    type = INSFVVelocityVariable
    initial_condition = 1e-15
  []
  [pressure]
    type = INSFVPressureVariable
  []
[]

[AuxVariables]
  [viz_x]
    family = MONOMIAL
    order = FIRST
  []
  [viz_y]
    family = MONOMIAL
    order = FIRST
  []
  [viz_z]
    family = MONOMIAL
    order = FIRST
  []
[]

[AuxKernels]
  [viz_x]
    type = ParsedAux
    function = 'u'
    args = 'u'
    execute_on = 'timestep_end'
    variable = viz_x
  []
  [viz_y]
    type = ParsedAux
    function = 'v'
    args = 'v'
    execute_on = 'timestep_end'
    variable = viz_y
  []
  [viz_z]
    type = ParsedAux
    function = 'w'
    args = 'w'
    execute_on = 'timestep_end'
    variable = viz_z
  []
[]


[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    vel = 'velocity'
    pressure = pressure
    u = u
    v = v
    w = w
    mu = ${mu}
    rho = ${rho}
  []

  [u_advection]
    type = INSFVMomentumAdvection
    variable = u
    advected_quantity = 'rhou'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    pressure = pressure
    u = u
    v = v
    w = w
    mu = ${mu}
    rho = ${rho}
  []
  [u_viscosity]
    type = FVDiffusion
    variable = u
    coeff = ${mu}
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = u
    momentum_component = 'x'
    p = pressure
  []

  [v_advection]
    type = INSFVMomentumAdvection
    variable = v
    advected_quantity = 'rhov'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    pressure = pressure
    u = u
    v = v
    w = w
    mu = ${mu}
    rho = ${rho}
  []
  [v_viscosity]
    type = FVDiffusion
    variable = v
    coeff = ${mu}
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = v
    momentum_component = 'y'
    p = pressure
  []

  [w_advection]
    type = INSFVMomentumAdvection
    variable = w
    advected_quantity = 'rhow'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    pressure = pressure
    u = u
    v = v
    w = w
    mu = ${mu}
    rho = ${rho}
  []
  [w_viscosity]
    type = FVDiffusion
    variable = w
    coeff = ${mu}
  []
  [w_pressure]
    type = INSFVMomentumPressure
    variable = w
    momentum_component = 'z'
    p = pressure
  []
[]

[FVBCs]
  [inlet-u]
    type = INSFVInletVelocityBC
    boundary = 'inlet_back'
    variable = u
    function = '0'
  []
  [inlet-v]
    type = INSFVInletVelocityBC
    boundary = 'inlet_back'
    variable = v
    function = '0'
  []
  [inlet-w]
    type = INSFVInletVelocityBC
    boundary = 'inlet_back'
    variable = w
    function = '1'
  []
  [walls-u]
    type = INSFVNoSlipWallBC
    boundary = 'elbow_bottom elbow_front elbow_left elbow_top inlet_bottom inlet_left inlet_right inlet_top outlet_back outlet_bottom outlet_front outlet_top'
    variable = u
    function = 0
  []
  [walls-v]
    type = INSFVNoSlipWallBC
    boundary = 'elbow_bottom elbow_front elbow_left elbow_top inlet_bottom inlet_left inlet_right inlet_top outlet_back outlet_bottom outlet_front outlet_top'
    variable = v
    function = 0
  []
  [walls-w]
    type = INSFVNoSlipWallBC
    boundary = 'elbow_bottom elbow_front elbow_left elbow_top inlet_bottom inlet_left inlet_right inlet_top outlet_back outlet_bottom outlet_front outlet_top'
    variable = w
    function = 0
  []
  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = 'outlet_right'
    variable = pressure
    function = '0'
  []
[]

[Materials]
  [ins_fv]
    type = INSFVMaterial
    u = 'u'
    v = 'v'
    w = 'w'
    pressure = 'pressure'
    rho = ${rho}
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      200                lu           NONZERO'
  line_search = 'none'
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
  [in]
    type = SideIntegralVariablePostprocessor
    variable = w
    boundary = 'inlet_back'
  []
  [out]
    type = SideIntegralVariablePostprocessor
    variable = u
    boundary = 'outlet_right'
  []
[]
