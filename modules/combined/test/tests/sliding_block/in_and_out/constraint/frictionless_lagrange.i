#  This is a benchmark test that checks constraint based frictionless
#  contact using the penalty method.  In this test a sinusoidal
#  displacement is applied in the horizontal direction to simulate
#  a small block come in and out of contact as it slides down a larger block.
#
#  The sinusoid is of the form 0.4sin(4t)+0.2. The gold file is run
#  on one processor and the benchmark
#  case is run on a minimum of 4 processors to ensure no parallel variability
#  in the contact pressure and penetration results.  Further documentation can
#  found in moose/modules/contact/doc/sliding_block/
#

[Mesh]
  file = sliding_elastic_blocks_2d_with_lagrange.e
  patch_size = 80
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./disp_x]
    block = '1 2'
  [../]
  [./disp_y]
    block = '1 2'
  [../]
  [./lm]
    block = 30
  [../]
[]

[AuxVariables]
  [./penetration]
  [../]
  [./inc_slip_x]
  [../]
  [./inc_slip_y]
  [../]
  [./accum_slip_x]
  [../]
  [./accum_slip_y]
  [../]
[]

[Functions]
  [./vertical_movement]
    type = ParsedFunction
    value = -t
  [../]
  [./horizontal_movement]
    type = ParsedFunction
    value = -0.04*sin(4*t)+0.02
  [../]
[]

[SolidMechanics]
  [./solid]
    disp_x = disp_x
    disp_y = disp_y
    block = '1 2'
  [../]
[]

[NodalKernels]
  [./lm]
    type = LM
    slave = 3
    master = 2
    block = 30
    variable = lm
  [../]
[]

[AuxKernels]
  [./zeroslip_x]
    type = ConstantAux
    variable = inc_slip_x
    boundary = 3
    execute_on = timestep_begin
    value = 0.0
  [../]
  [./zeroslip_y]
    type = ConstantAux
    variable = inc_slip_y
    boundary = 3
    execute_on = timestep_begin
    value = 0.0
  [../]
  [./accum_slip_x]
    type = AccumulateAux
    variable = accum_slip_x
    accumulate_from_variable = inc_slip_x
    execute_on = timestep_end
  [../]
  [./accum_slip_y]
    type = AccumulateAux
    variable = accum_slip_y
    accumulate_from_variable = inc_slip_y
    execute_on = timestep_end
  [../]
  [./penetration]
    type = PenetrationAux
    variable = penetration
    boundary = 3
    paired_boundary = 2
  [../]
[]

[Postprocessors]
  [./nonlinear_its]
    type = NumNonlinearIterations
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  [../]
  [./left_y]
    type = DirichletBC
    variable = disp_y
    boundary = 1
    value = 0.0
  [../]
  [./right_x]
    type = FunctionPresetBC
    variable = disp_x
    boundary = 4
    function = horizontal_movement
  [../]
  [./right_y]
    type = FunctionPresetBC
    variable = disp_y
    boundary = 4
    function = vertical_movement
  [../]
[]

[Materials]
  [./left]
    type = LinearIsotropicMaterial
    block = 1
    disp_y = disp_y
    disp_x = disp_x
    poissons_ratio = 0.3
    youngs_modulus = 1e3
  [../]
  [./right]
    type = LinearIsotropicMaterial
    block = 2
    disp_y = disp_y
    disp_x = disp_x
    poissons_ratio = 0.3
    youngs_modulus = 1e3
  [../]
  [./dummy]
    type = GenericConstantMaterial
    block = 30
    prop_names = 'dumb'
    prop_values = '1'
  [../]
[]

[Preconditioning]
  [./smp]
     type = SMP
     full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'

  petsc_options = '-options_left'

  l_max_its = 30
  dt = 0.1
  end_time = 3
  dtmin = 0.01
  nl_max_its = 10

  # [./Predictor]
  #   type = SimplePredictor
  #   scale = 1.0
  # [../]
[]

[Outputs]
  file_base = frictionless_lagrange_out
  checkpoint = true
  print_linear_residuals = false
  [./exodus]
    type = Exodus
    elemental_as_nodal = true
  [../]
  [./dof_map]
    execute_on = 'initial'
    type = DOFMap
  [../]
[]

[Contact]
  [./leftright]
    slave = 3
    master = 2
    model = frictionless
    formulation = lagrange
    system = constraint
    lm = lm
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]
