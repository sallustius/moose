[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 1
[]

[Variables]
  [./u]
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  [../]
  [./right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  [../]
[]

[Preconditioning]
  [./fdp]
    type = FDP
    full = true
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options = '-pc_svd_monitor -ksp_view_pmat'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'svd'
[]

[Outputs]
  exodus = true
[]

[AuxVariables]
  [./random]
  [../]
[]

[ICs]
  [./random]
    type = RandomIC
    variable = random
  [../]
[]

[Adaptivity]
  marker = box
  max_h_level = 1
  initial_steps = 1

  [./Indicators]
    [./random_indicator]
      type = ElementIntegralIndicator
      variable = random
    [../]
  [../]

  [./Markers]
    [./errorfrac_random]
      type = ErrorFractionMarker
      refine = 0.1
      coarsen = 0
      indicator = random_indicator
    [../]
    [./uniform]
      type = UniformMarker
      mark = 'refine'
    [../]
    [./box]
      type = BoxMarker
      bottom_left = '0.5 0 0'
      top_right = '1 1 0'
      inside = 'refine'
      outside = 'do_nothing'
    [../]
  [../]
[]
