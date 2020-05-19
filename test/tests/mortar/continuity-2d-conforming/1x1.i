[Mesh]
  second_order = true
  [./left_block]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = -1
    xmax = 0
    ymin = 0
    ymax = 1
    nx = 1
    ny = 2
    elem_type = QUAD4
  [../]
  [./left_block_sidesets]
    type = RenameBoundaryGenerator
    input = left_block
    old_boundary_id = '0 1 2 3'
    new_boundary_name = 'lb_bottom lb_right lb_top lb_left'
  [../]
  [./left_block_id]
    type = SubdomainIDGenerator
    input = left_block_sidesets
    subdomain_id = 1
  [../]
  [./right_block]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    nx = 1
    ny = 2
    elem_type = QUAD4
  [../]
  [./right_block_id]
    type = SubdomainIDGenerator
    input = right_block
    subdomain_id = 2
  [../]
  [./combined]
    type = MeshCollectionGenerator
    inputs = 'left_block_id right_block_id'
  [../]
  [./block_rename]
    type = RenameBlockGenerator
    input = combined
    old_block_id = '1 2'
    new_block_name = 'left_block right_block'
  [../]
  [./right_block_sidesets]
    type = SideSetsFromPointsGenerator
    input = block_rename
    points = '0.5 0 0
              1 0.5 0
              0.5 1 0
              0 0.5 0'
    new_boundary = 'rb_bottom rb_right rb_top rb_left'
  [../]
  [slave]
    input = right_block_sidesets
    type = LowerDBlockFromSidesetGenerator
    sidesets = 'lb_right'
    new_block_id = '10001'
    new_block_name = 'slave_lower'
  []
  [master]
    input = slave
    type = LowerDBlockFromSidesetGenerator
    sidesets = 'rb_left'
    new_block_id = '10000'
    new_block_name = 'master_lower'
  []
[]

[Functions]
  [./exact_sln]
    type = ParsedFunction
    value = y
  [../]
  [./ffn]
    type = ParsedFunction
    value = 0
  [../]
[]

[Variables]
  [./u]
    order = SECOND
    family = LAGRANGE
    block = 'left_block right_block'
  [../]

  [./lm]
    order = FIRST
    family = LAGRANGE
    block = 'slave_lower'
    # use_dual = true
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
    block = 'left_block right_block'
  [../]
  [./ffn]
    type = BodyForce
    variable = u
    function = ffn
    block = 'left_block right_block'
  [../]
[]

[Constraints]
  [./ced]
    type = EqualValueConstraint
    variable = lm
    slave_variable = u
    master_boundary = rb_left
    master_subdomain = master_lower
    slave_boundary = lb_right
    slave_subdomain = slave_lower
  [../]
[]

[BCs]
  [./all]
    type = FunctionDirichletBC
    variable = u
    boundary = 'lb_left lb_top lb_bottom rb_top rb_bottom rb_right'
    # boundary = 'lb_left rb_right'
    function = exact_sln
  [../]
[]

[Postprocessors]
  [./l2_error]
    type = ElementL2Error
    variable = u
    function = exact_sln
    block = '1 2'
    execute_on = 'initial timestep_end'
  [../]
[]

[Preconditioning]
  [./fmp]
    type = SMP
    full = true
    solve_type = 'NEWTON'
  [../]
[]

[Executioner]
  type = Steady

  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_ksp_ew'
  nl_rel_tol = 1e-12
[]

[Outputs]
  exodus = true
[]
