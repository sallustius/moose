# velocity=1

# [GlobalParams]
#   u = ${velocity}
#   p = 0
#   tau_type = opt
#   alpha = 0
# []

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

# [AuxVariables]
#   [exact]
#     order = SECOND
#   []
# []

# [AuxKernels]
#   [exact]
#     type = FunctionAux
#     function = exact_soln
#     variable = exact
#   []
# []

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
  # [./adv]
  #   type = Advection
  #   variable = u
  # [../]
  [diff]
    type = Diffusion
    variable = u
  []
  [mms_forcing_function]
    type = BodyForce
    variable = u
    function = forcing_function
  []
[]

[BCs]
  [./left]
    type = FunctionDirichletBC
    variable = u
    boundary = 'left right top bottom'
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
