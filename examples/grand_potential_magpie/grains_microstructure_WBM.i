# 2 grains in UO2 under irradiation
# missing: nucleation of voids/bubbles, sinks (GB, dislocations)

[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  type = GeneratedMesh
  dim = 2 # Problem dimension
  nx = 10 # Number of elements in the x-direction
  ny = 10 # Number of elements in the y-direction
  xmin = 0 # minimum x-coordinate of the mesh
  xmax = 10 # maximum x-coordinate of the mesh
  ymax = 10 # maximum y-coordinate of the mesh
  elem_type = QUAD4 # Type of elements used in the mesh
  uniform_refine = 3 # Initial uniform refinement of the mesh
[]

[GlobalParams]
  op_num = 10
  var_name_base = gr
  time_scale = 1
  grain_num = 10
[]

[Variables]
  [./PolycrystalVariables]
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
[]
#
#[Modules]
#  [./PhaseField]
#    [./Conserved]
#      [./cv]
#        kappa = kappa
#        free_energy = F1
#        mobility = M
#        solve_type = FORWARD_SPLIT
#        args = 'ci cg'
#      [../]
#      [./ci]
#        kappa = kappa
#        free_energy = F1
#        mobility = M
#        solve_type = FORWARD_SPLIT
#        args = 'cv cg'
#      [../]
#      [./cg]
#        kappa = kappa
#        free_energy = F1
#        mobility = M
#        solve_type = FORWARD_SPLIT
#        args = 'ci cv'
#      [../]
#    [../]
#  [../]
#[]

[ICs]
  #[./cg]
  #  type = ConstantIC
  #  variable = cg
  #  value = 1e-5
  #[../]
  #[./cv]
  #  type = ConstantIC
  #  variable = cv
  #  value = 1e-5
  #[../]
  #[./ci]
  #  type = ConstantIC
  #  variable = ci
  #  value = 1e-5
  #[../]
[./PolycrystalICs]
  [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
   [../]
[../]
[]

[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = bt
  [../]
[]

[AuxVariables]
  [./local_energy]
    family = MONOMIAL
    order =  CONSTANT
  [../]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  # AuxKernel block, defining the equations used to calculate the auxvars
  #[./local_energy]
  #  type = TotalFreeEnergy
  #  variable = local_energy
  #  f_name = F1
  #  interfacial_vars = 'cg   ci     cv'
  #  kappa_names = 'kappa   kappa   kappa'
  #  execute_on = timestep_end
  #[../]
  [./bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
[]

[BCs]
  [./Periodic]
    [./all]
      auto_direction = 'x y'
    [../]
  [../]

[]

[Materials]
  #[./floc] # Ideal solution model
  #  type = DerivativeParsedMaterial
  #  f_name = F1
  #  args = 'ci cv cg'
  #  function = '1/Vm * (Na*Ef*(cv + ci + cg) + R*T*(cg*plog((cg),tol) + (1-cg)*plog((1-cg),tol) + cv*plog((cv),tol) + (1-cv)*plog((1-cv),tol) + ci*plog((ci),tol) + (1-ci)*plog((1-ci),tol) ))'
  #  constant_names =       'Ef         Vm          R           Na     tol    T'
  #  constant_expressions = '3          2.4389e22  5.1896e19  6.022e23   1e-7 1200'
  #  derivative_order = 2
  #  #outputs = exodus
  #[../]
  #[./constants]
  #  type = GenericConstantMaterial
  #  prop_names  = 'kappa    M'
  #  prop_values = '1        1'
  #[../]
  [./UO2GrGr]
    # Material properties
    type = UO2GrGr # Quantitative material properties for copper grain growth.  Dimensions are nm and ns
    block = 0 # Block ID (only one block in this problem)
    T = 1200 # Constant temperature of the simulation (for mobility calculation)
    wGB = 0.5 # Width of the diffuse GB
    length_scale = 1e-06
  [../]
[]


[Postprocessors]
  [./dt]
    # Outputs the current time step
    type = TimestepSize
  [../]
  [./total_free_energy]
    type = ElementIntegralVariablePostprocessor
    variable = local_energy
  [../]
  #[./cgaverage]
  #  type = Receiver
  #[../]

[]

[Preconditioning]
    active = 'SMP'
    [./SMP]
      type = SMP
      full = true
    [../]
[]

#[MultiApps]
#  [./radiation_damage_app]
#    type = TransientMultiApp
#    input_files = bubbles_nanostructure.i
#    execute_on = timestep_begin
#
#  []
#[]
#
#[Transfers]
#  [./radiation_damage_transfer]
#    type = MultiAppPostprocessorTransfer
#    multi_app = radiation_damage_app
#    direction = from_multiapp
#    to_postprocessor = cgaverage
#    from_postprocessor = cg_average
#    reduction_type =  average
#  [../]
#[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = 'NEWTON'
 l_max_its = 15
  nl_max_its = 15
  l_tol = 1.0e-6
  nl_rel_tol = 1.0e-10
  nl_abs_tol = 1.0e-10
  start_time = 0.0
  num_steps = 120
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  petsc_options_value = 'hypre boomeramg 101 ds'
  dtmax = 1e+6
  [./TimeStepper]
    type = IterationAdaptiveDT
        dt = 1000
        optimal_iterations = 8
        iteration_window = 2
  [../]
  [./Adaptivity]
  #  Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
    initial_adaptivity = 2 # Number of times mesh is adapted to initial condition
    refine_fraction = 0.7 # Fraction of high error that will be refined
    coarsen_fraction = 0.1 # Fraction of low error that will coarsened
    max_h_level = 4 # Max number of refinements used, starting from initial mesh (before uniform refinement)
  [../]
[]
#
[Debug]
  show_var_residual_norms = true
[]

[Outputs]
  exodus = true
  #csv = true
[]
