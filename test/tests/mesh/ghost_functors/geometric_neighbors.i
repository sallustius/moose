[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 4
  ny = 4

  # We are testing geometric ghosted functors
  # so we have to use distributed mesh
  parallel_type = distributed
[]

[Variables]
  [./u]
  [../]

  [./v]
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]

  [./diff_v]
    type = Diffusion
    variable = v
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
  [./left_v]
    type = DirichletBC
    variable = v
    boundary = left
    value = 0
  [../]
  [./right_v]
    type = DirichletBC
    variable = v
    boundary = right
    value = 10
  [../]
[]


[AuxVariables]
  [./ghosted_elements]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./proc]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./random_elemental]
    type = GhostAux
    variable = ghosted_elements
    ghost_user_object = ghost_uo
    execute_on = timestep_end
  [../]
  [./proc]
    type = ProcessorIDAux
    variable = proc
    execute_on = initial
  [../]
[]

[UserObjects]
  [./ghost_uo]
    type = GhostUserObject
    element_side_neighbor_layers = 2
    some_variable = v
    execute_on = timestep_end
  [../]
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]
