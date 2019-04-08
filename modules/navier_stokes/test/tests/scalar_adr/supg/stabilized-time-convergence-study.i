[Mesh]
  type = GeneratedMesh
  dim  = 2
  xmin = -1
  xmax = 1
  ymin = -1
  ymax = 1
  nx   = 4
  ny   = 4
  elem_type = QUAD9
[]


[Variables]
  [./u]
    order = SECOND
  [../]
[]

[Kernels]
  [time]
    type = TimeDerivative
    variable = u
  []
  [adv]
    type = Advection
    variable = u
    u = 1
    p = 0
    tau_type = mod
    alpha = 0
  []
  # [diff]
  #   type = Diffusion
  #   variable = u
  # []
  [mms_forcing_function]
    type = BodyForce
    variable = u
    function = forcing_function
  []
  # [supg]
  #   type = SUPG
  #   variable = u
  #   forcing_func = forcing_function
  #   velocity = '1 0 0'
  #   include_transient_term = true
  # []
  # [time_supg]
  #   type = TimeDerivativeSUPG
  #   variable = u
  #   velocity = '1 0 0'
  # []
[]

[BCs]
  [./left]
    type = FunctionDirichletBC
    variable = u
    # boundary = 'left right top bottom'
    boundary = 'left'
    function = exact_soln
  [../]
[]

[Materials]
  [./mat]
    # These Materials are required by the INSBase class; we don't use them for anything.
    type = GenericConstantMaterial
    prop_names = 'mu rho'
    prop_values = '0 1'
  [../]
[]

[Functions]
  [./forcing_function]
    type = ParsedFunction
    value = ''
  [../]
  [exact_soln]
    type = ParsedFunction
    value = ''
  []
[]

[Executioner]
  type = Transient
  end_time = 1
  solve_type = PJFNK
[]

[Outputs]
  exodus = true
  csv = true
[]

[Postprocessors]
  [l2_err]
    type = ElementL2Error
    variable = u
    function = exact_soln
    outputs = 'csv console'
  []
[]

[ICs]
  [ini]
    type = FunctionIC
    variable = u
    function = exact_soln
  []
[]
