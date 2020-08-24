[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 20
  []
  [subdomain]
    type = SubdomainBoundingBoxGenerator
    block_id = 1
    bottom_left = '0.5 0 0'
    top_right = '1 1 0'
    input = gmg
  []
[]

[Problem]
  type = FEProblem
[]

[Variables]
  [u]
    initial_condition = 1
    block = 0
  []
  [v]
    initial_condition = 10
    block = 1
  []
  [w]
    initial_condition = 1
  []
[]

[Kernels]
  [diff_u]
    type = Diffusion
    variable = u
    block = 0
  []
  [force_u]
    type = BodyForce
    variable = u
    block = 0
    function = 1
  []

  [diff_v]
    type = Diffusion
    variable = v
    block = 1
  []
  [force_v]
    type = BodyForce
    variable = v
    block = 1
    function = 1
  []

  [diff_w]
    type = Diffusion
    variable = w
  []
  [force_w]
    type = BodyForce
    variable = w
    function = 1
  []
[]

[AuxVariables]
  [u_aux]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
    block = 0
  []
  [v_aux]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
    block = 1
  []

  # w_aux should exist in both blocks
  [w_aux]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
    block = '0 1'
  []
  [b0_var1]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
    block = 0
  []
  [b0_var2]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
    block = 0
  []
  [b0_var3]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
    block = 0
  []
  [b0_var4]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
    block = 0
  []
  [b1_var1]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
    block = 1
  []
  [b1_var2]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
    block = 1
  []
  [b1_var3]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
    block = 1
  []
  [b1_var4]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
    block = 1
  []

  # This should also exist in both blocks
  [b12_var1]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
  []
[]

[BCs]
  [u_left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [v_right]
    type = DirichletBC
    variable = v
    boundary = right
    value = 0
  []
  [w_left]
    type = DirichletBC
    variable = w
    boundary = left
    value = 0
  []
  [w_right]
    type = DirichletBC
    variable = w
    boundary = right
    value = 0
  []
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Steady
  solve_type = newton
[]

[Outputs]
  [out_test]
    type = Exodus
  []
[]
