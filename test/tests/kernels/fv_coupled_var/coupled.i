[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
[]

[Variables]
  [u]
  []
  [v]
    family = MONOMIAL
    order = CONSTANT
    fv = true
    initial_condition = 1
  []
  [w]
    family = MONOMIAL
    order = CONSTANT
    fv = true
    initial_condition = 1
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
  [rxn]
    type = Reaction
    variable = u
  []
[]

[FVKernels]
  [diff]
    type = FVDiffusion
    variable = v
    coeff = coeff
  []
  [rxn]
    type = FVReaction
    variable = v
  []
  [diffw]
    type = FVDiffusion
    variable = w
    coeff = coeff
  []
  # [rxnw]
  #   type = FVReaction
  #   variable = w
  # []
[]

[FVBCs]
  [left]
    type = FVDirichletBC
    variable = v
    boundary = left
    value = 0
  []
  [right]
    type = FVDirichletBC
    variable = v
    boundary = right
    value = 1
  []
  [leftw]
    type = FVDirichletBC
    variable = w
    boundary = left
    value = 0
  []
  [rightw]
    type = FVDirichletBC
    variable = w
    boundary = right
    value = 1
  []
[]

[Materials]
  [diff]
    type = ADGenericConstantMaterial
    prop_names = 'coeff'
    prop_values = '1'
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  []
[]

[Problem]
  kernel_coverage_check = off
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
[]

[Outputs]
  exodus = true
  [dof]
    type = DOFMap
    execute_on = 'initial'
  []
[]
