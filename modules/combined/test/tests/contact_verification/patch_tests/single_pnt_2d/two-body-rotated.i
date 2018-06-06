[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Mesh]
  file = two-body-rotated.e
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

[Constraints]
  [./lm]
    type = LMConstraint
    slave = 10
    master = 20
    variable = lm
    master_variable = disp_x
  [../]
[]

# [Kernels]
#   [./TensorMechanics]
#     use_displaced_mesh = true
#     block = '1 2'
#   [../]
# []

# [BCs]
#   [./botx]
#     type = DirichletBC
#     variable = disp_x
#     boundary = 40
#     value = 0.0
#   [../]
#   [./boty]
#     type = DirichletBC
#     variable = disp_y
#     boundary = 40
#     value = 0.0
#   [../]
#   [./topx]
#     type = DirichletBC
#     variable = disp_x
#     boundary = 30
#     value = -1.0606601717798212
#   [../]
#   [./topy]
#     type = DirichletBC
#     variable = disp_y
#     boundary = 30
#     value = -1.0606601717798212
#   [../]
# []

[Materials]
  [./bot_elas_tens]
    type = ComputeIsotropicElasticityTensor
    block = '1'
    youngs_modulus = 1
    poissons_ratio = 0.3
  [../]
  [./bot_strain]
    type = ComputePlaneFiniteStrain
    block = '1'
  [../]
  [./bot_stress]
    type = ComputeFiniteStrainElasticStress
    block = '1'
  [../]
  [./top_elas_tens]
    type = ComputeIsotropicElasticityTensor
    block = '2'
    youngs_modulus = 1
    poissons_ratio = 0.3
  [../]
  [./top_strain]
    type = ComputePlaneFiniteStrain
    block = '2'
  [../]
  [./top_stress]
    type = ComputeFiniteStrainElasticStress
    block = '2'
  [../]
  [./dummy]
    type = GenericConstantMaterial
    block = '3'
    prop_names = 'dumb'
    prop_values = '1'
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 2
  dt = 1
  dtmin = 1
  solve_type = 'NEWTON'
  line_search = 'bt'

  l_max_its = 100
  nl_max_its = 200
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
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
