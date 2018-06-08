[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  file = two-body.e
[]

[Problem]
  kernel_coverage_check = false
[]

[Variables]
  [./disp_x]
    block = '1 2'
  [../]
  [./disp_y]
    block = '1 2'
  [../]
  [./lm]
    block = 3
  [../]
[]

[ICs]
  [./block2]
    block = 2
    variable = disp_y
    type = ConstantIC
    value = -1
  [../]
[]

[Kernels]
  [./disp_x]
    type = Diffusion
    variable = disp_x
  [../]
  [./disp_y]
    type = Diffusion
    variable = disp_y
  [../]
[]

[Constraints]
  [./lm]
    type = LMConstraint
    slave = 10
    master = 20
    variable = lm
    master_variable = disp_x
    disp_y = disp_y
  [../]
[]

[BCs]
  [./botx]
    type = DirichletBC
    variable = disp_x
    boundary = 40
    value = 0.0
  [../]
  [./boty]
    type = DirichletBC
    variable = disp_y
    boundary = 40
    value = 0.0
  [../]
  [./topx]
    type = DirichletBC
    variable = disp_x
    boundary = 30
    value = 0
  [../]
  [./topy]
    type = DirichletBC
    variable = disp_y
    boundary = 30
    value = -1.5
  [../]
  # [./topy]
  #   type = NeumannBC
  #   variable = disp_y
  #   boundary = 30
  #   value = 1
  # [../]
[]

# [SolidMechanics]
#   [./solid]
#     disp_x = disp_x
#     disp_y = disp_y
#     block = '1 2'
#   [../]
# []

# [Materials]
#   [./bodies]
#     type = LinearIsotropicMaterial
#     block = '1 2'
#     disp_x = disp_x
#     disp_y = disp_y
#     poissons_ratio = 0.3
#     youngs_modulus = 1
#   [../]
#   [./dummy]
#     type = GenericConstantMaterial
#     block = '3'
#     prop_names = 'dumb'
#     prop_values = '1'
#   [../]
# []

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1
  dtmin = 1
  solve_type = 'NEWTON'
  line_search = 'bt'

  l_max_its = 100
  nl_max_its = 100
  l_tol = 1e-3
[]

[Outputs]
  exodus = true
  checkpoint = true
  [./dof_map]
    type = DOFMap
  [../]
[]

[Contact]
  [./leftright]
    master = 20
    slave = 10
    model = frictionless
    formulation = lagrange
    # penalty = 1e6
    system = constraint
    lm = lm
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]
