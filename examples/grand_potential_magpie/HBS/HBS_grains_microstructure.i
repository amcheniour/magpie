# This is based on Larry's grandpotential.i file in marmot-problems/Magpie/examples
# The domain is 2D 400 nm x 400 nm with a xenon bubble in the middle of UO2 matrix, radius = 44 nm
# The variables are chemical potentials of U vacancies and interstitials, Xe atoms
# This is solving the rate theory equation for wi, wv, wg
# It is using a source term from MAGPIE based on BCMC

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 100
  ny = 100
  xmin = 0
  xmax = 3000
  ymin = 0
  ymax = 3000
  # uniform_refine = 2
[]

[GlobalParams]
  op_num = 2
  grain_num = 2
  var_name_base = gr
[]

[Variables]
  # chemical potentials for vacancies, interstitials and gas atoms
  [./wv]
  [../]
  [./wi]
  [../]
  [./wg]
  [../]

  [./PolycrystalVariables]
  [../]
[]

[Functions]
  # burnup function
  [./fburnup] # in GWD/tHM
    type = ParsedFunction
    value = 'q*sigmaf_u235*flux*t*950'
    vars = 'q sigmaf_u235 flux'
    vals = '0.045 5.906e-22 1.5e13'
  [../]
[]

[AuxVariables]
  # grain boundaries
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./burnup_var]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./local_energy]
    family = MONOMIAL
    order =  CONSTANT
  [../]
  # number density variables created for the rasterizer in Magpie
  [./rho_U235]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rho_U238]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rho_O]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rhog_var]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # concentration of Xe
  [./c_Xe]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  [../]
  # vacancies
  [./cv]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # interstitials
  [./ci]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # number density of interstitials for the rasterizer
  [./rhoi_var]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # temperature
  [./T]
    initial_condition = 600
  [../]


[]

[ICs]
  [./PolycrystalICs]
    [./BicrystalCircleGrainIC]
      radius = 1300
       x = 1500
       y = 1500
    [../]
  [../]

  [./IC_wv]
    type = ConstantIC
    value = 0.0
    variable = wv
  [../]
  [./IC_wi]
    type = ConstantIC
    value = 0.0
    variable = wi
  [../]
  [./IC_wg]
    type = ConstantIC
    value = 0.0
    variable = wg
  [../]
[]


[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = bt
  [../]
[]

[Modules]
  [./PhaseField]
    [./GrandPotential]
      switching_function_names = 1

      chemical_potentials = 'wv wi wg'
      mobilities = 'Dchiv Dchii Dchig'
      susceptibilities = 'chiv chii chig'
      free_energies_w = 'rhov rhoi rhog'

      gamma_gr = 'gmm'
      free_energies_gr = 'omegam'
      kappa_gr = kappa
      mobility_name_gr = L
      var_name_base = gr

    [../]
  [../]
[]

[Kernels]
  [./vac_source_U238] # Magpie source term for U238
    type = MyTRIMElementSource
    variable = wv
    runner = runner
    ivar = 1
    defect = VAC
    prefactor = 0.04092
  [../]
  [./vac_source_U235] # Magpie source term for U235
    type = MyTRIMElementSource
    variable = wv
    runner = runner
    ivar = 0
    defect = VAC
    prefactor = 0.04092
  [../]
  [./recombination_v] # Recombination term = Recombination rate * rhoi * rhov
    type = GrandPotentialRecombination
    variable = wv
    rho = rhov
    rho_r = rhoi
    value = 1
    D = Di
    omega = 0.04092
    a0 = 0.25
    Z = 50
    hm = 1
    args = wi
  [../]
  [./dislocation_sink_v] # Dislocation sink = sink strength * Diffusion coefficient * rho
     type = GrandPotentialSink
     variable = wv
     rho = rhov
     D = D
     mask = 1
     value = 1
     sink_strength = dislocation_density
 [../]

  [./source_i_U238]
    type = MyTRIMElementSource
    variable = wi
    runner = runner
    ivar = 1
    defect = INT
    prefactor = 0.04092
  [../]
  [./source_i_U235]
    type = MyTRIMElementSource
    variable = wi
    runner = runner
    ivar = 0
    defect = INT
    prefactor = 0.04092
  [../]
  [./source_i_i]
    type = MyTRIMElementSource
    variable = wi
    runner = runner
    ivar = 4
    defect = INT
    prefactor = 0.04092
  [../]
  [./source_i_v]
    type = MyTRIMElementSource
    variable = wi
    runner = runner
    ivar = 4
    defect = VAC
    prefactor = -0.04092
  [../]
  [./recombination_i]
    type = GrandPotentialRecombination
    variable = wi
    rho = rhoi
    rho_r = rhov
    value = 1
    D = Di
    omega = 0.04092
    a0 = 0.25
    Z = 50
    hm = 1
    args = wv
  [../]
  [./dislocation_sink_i]
     type = GrandPotentialSink
     variable = wi
     rho = rhoi
     mask = 1
     D = Di
     value = 1.02
     sink_strength = dislocation_density
   [../]
  [./source_g]
    type = MyTRIMElementSource
    variable = wg
    runner = runner
    ivar = 3
    defect = INT
    prefactor = 0.04092
  [../]
[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x y' # 2D
    [../]
  [../]
[]

[AuxKernels]

  # Grain boundaries
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]

  # Number densities for different concentrations (U,O,Xe,Uint)
  [./rho_U235_aux] # U235
    type = MaterialRealAux
    property = rho_U235_matl
    variable = rho_U235
    execute_on = 'initial timestep_begin'
  [../]
  [./rho_U238_aux] # U238
    type = MaterialRealAux
    property = rho_U238_matl
    variable = rho_U238
    execute_on = 'initial timestep_begin'
  [../]
  [./rho_O] # Oxygen
    type = MaterialRealAux
    property = rho_O_matl
    variable = rho_O
    execute_on = 'initial timestep_begin'
  [../]
  [./rhog_aux] # Xenon
    type = MaterialRealAux
    property = rhog
    variable = rhog_var
    execute_on = 'initial timestep_begin'
  [../]
  [./rhoi_aux] # U interstitial
    type = MaterialRealAux
    property = rhoi
    variable = rhoi_var
    execute_on = 'initial timestep_begin'
  [../]
  [./local_energy]
    type = TotalFreeEnergy
    variable = local_energy
    f_name = free_energy
    interfacial_vars = 'cv   ci     c_Xe'
    kappa_names = 'kappa_c   kappa_c   kappa_c'
    execute_on = timestep_end
  [../]
  # Number fractions of point defects
  [./cv_aux]
    type = MaterialRealAux
    property = cv_mat
    variable = cv
    execute_on = 'initial timestep_begin'
  [../]
  [./ci_aux]
    type = MaterialRealAux
    property = ci_mat
    variable = ci
    execute_on = 'initial timestep_begin'
  [../]
  [./c_Xe]
    type = MaterialRealAux
    property = c_Xe_matl
    variable = c_Xe
    execute_on = 'initial timestep_begin'
  [../]
[]

[Materials]
  [./omegam]
    type = DerivativeParsedMaterial
    args = 'wv wi wg'
    f_name = omegam
    material_property_names = 'Va kvmatrix cvmatrixeq kimatrix cimatrixeq kgmatrix cgmatrixeq'
    function = '-0.5*wv^2/Va^2/kvmatrix-wv/Va*cvmatrixeq-0.5*wi^2/Va^2/kimatrix-wi/Va*cimatrixeq-0.5*wg^2/Va^2/kgmatrix-wg/Va*cgmatrixeq'
    derivative_order = 2
    outputs = exodus
  [../]

  [./free_energy]
    type = ParsedMaterial
    args = 'wv wi wg'
    material_property_names = 'rhov rhoi rhog omegam'
    f_name = free_energy
    function = 'omegam + rhov*wv + rhoi*wi + rhog*wg'
    outputs = exodus
  [../]

  [./rhov]
    type = DerivativeParsedMaterial
    args = 'wv'
    f_name = rhov
    material_property_names = 'Va kvmatrix cvmatrixeq'
    function = 'wv/Va^2/kvmatrix + cvmatrixeq/Va'
    derivative_order = 2
    # outputs = exodus
  [../]

  [./rhoi]
    type = DerivativeParsedMaterial
    args = 'wi'
    f_name = rhoi
    material_property_names = 'Va kimatrix cimatrixeq'
    function = 'wi/Va^2/kimatrix + cimatrixeq/Va'
    derivative_order = 2
    #outputs = exodus
  [../]

  [./rhog]
    type = DerivativeParsedMaterial
    args = 'wg'
    f_name = rhog
    material_property_names = 'Va kgmatrix cgmatrixeq'
    function = 'wg/Va^2/kgmatrix + cgmatrixeq/Va'
    derivative_order = 2
    #outputs = exodus
  [../]


  # U and O number densities
  [./rho_U235_matl]
    type = ParsedMaterial
    f_name = rho_U235_matl
    constant_names = Va
    constant_expressions = 0.04092
    function = '0.01485/Va'
    # outputs = exodus
  [../]
  [./rho_U238_matl]
    type = ParsedMaterial
    f_name = rho_U238_matl
    constant_names = Va
    constant_expressions = 0.04092
    function = '0.31515/Va'
    # outputs = exodus
  [../]
  [./rho_O_matl]
    type = ParsedMaterial
    f_name = rho_O_matl
    constant_names = Va
    constant_expressions = 0.04092
    function = '0.66/Va/2'
    # outputs = exodus
  [../]

  # Constants
  [./const]
    type = GenericConstantMaterial
    prop_names =  ' D    Va      cvbubeq cgbubeq cibubeq  kgbub  kvbub kibub gmb     gmm T    Efvbar    Efgbar    kTbar     f0     tgrad_corr_mult  kappa_c kappa_op gamma_asymm Di'
    prop_values = ' 0.01 0.04092 0.5459  0.4541  0.0      1.41   1.41  1.41  0.9218 1.5 1200 7.505e-3  7.505e-3  2.588e-4  0.0    0.0              1.0     0.5273   1.5         1 '
  [../]

  [./kappa]
    type = ParsedMaterial
    f_name = kappa
    material_property_names = 'T'
    constant_names = 'l_int sigma_GB'
    constant_expressions = '100 9.36226'
    function = '3/4 * sigma_GB * l_int'
  [../]

  [./mu]
    type = ParsedMaterial
    f_name = mu
    material_property_names = 'T'
    constant_names = 'l_int sigma_GB'
    constant_expressions = '100 9.36226'
    function = '6 * sigma_GB / l_int'
  [../]

  [./L]
    type = ParsedMaterial
    f_name = L
    material_property_names = 'T'
    constant_names = 'M0 l_int Q kB'
    constant_expressions = '1.475989e9 100 2.77 8.617343e-5'
    function = '4/3 * M0 * exp(-Q/(kB*T)) / l_int'
  [../]
  # Equilibrium concentrations of defects
  [./cvmatrixeq]    #For values, see Li et al., Nuc. Inst. Methods in Phys. Res. B, 303, 62-27 (2013).
    type = ParsedMaterial
    f_name = cvmatrixeq
    material_property_names = 'T'
    constant_names        = 'kB           Efv'
    constant_expressions  = '8.6173324e-5 3.0'
    function = 'exp(-Efv/(kB*T))'
  [../]
  [./cimatrixeq]
    type = ParsedMaterial
    f_name = cimatrixeq
    material_property_names = 'T'
    constant_names        = 'kB           Efi'
    constant_expressions  = '8.6173324e-5 3.0'
    function = 'exp(-Efi/(kB*T))'
  [../]
  [./cgmatrixeq]
    type = ParsedMaterial
    f_name = cgmatrixeq
    material_property_names = 'T'
    constant_names        = 'kB           Efg'
    constant_expressions  = '8.6173324e-5 3.0'
    function = 'exp(-Efg/(kB*T))'
  [../]

  # See free energy expression
  [./kvmatrix_parabola]
    type = ParsedMaterial
    f_name = kvmatrix
    material_property_names = 'T  cvmatrixeq'
    constant_names        = 'c0v  c0g  a1                                               a2'
    constant_expressions  = '0.01 0.01 0.178605-0.0030782*log(1-c0v)+0.0030782*log(c0v) 0.178605-0.00923461*log(1-c0v)+0.00923461*log(c0v)'
    function = '((-a2+3*a1)/(4*(c0v-cvmatrixeq))+(a2-a1)/(2400*(c0v-cvmatrixeq))*T)'
    # outputs = exodus
  [../]
  [./kimatrix_parabola]
    type = ParsedMaterial
    f_name = kimatrix
    material_property_names = 'kvmatrix'
    function = 'kvmatrix'
  [../]
  [./kgmatrix_parabola]
    type = ParsedMaterial
    f_name = kgmatrix
    material_property_names = 'kvmatrix'
    function = 'kvmatrix'
  [../]
  [./c_Xe_matl]
    type = ParsedMaterial
    f_name = c_Xe_matl
    material_property_names = 'Va rhog'
    function = 'Va * rhog'
  [../]
  [./cv_mat]
    type = ParsedMaterial
    f_name = cv_mat
    material_property_names = 'Va rhov'
    function = 'Va * rhov'
    outputs = exodus
  [../]
  [./ci_mat]
    type = ParsedMaterial
    f_name = ci_mat
    material_property_names = 'Va rhoi'
    function = 'Va * rhoi'
    outputs = exodus
  [../]

  # Mobilities for Grand Potential Model
  [./Mobility_v]
    type = DerivativeParsedMaterial
    f_name = Dchiv
    material_property_names = 'D chiv'
    function = 'D*chiv'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./Mobility_i]
    type = DerivativeParsedMaterial
    f_name = Dchii
    material_property_names = 'Di chii'
    function = 'Di*chii'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./Mobility_g]
    type = DerivativeParsedMaterial
    f_name = Dchig
    material_property_names = 'Dtot chig'
    function = 'Dtot*chig'
    derivative_order = 2
    #outputs = exodus
  [../]

  # Susceptibilities
  [./chiv]
    type = DerivativeParsedMaterial
    f_name = chiv
    material_property_names = 'Va kvmatrix '
    function = '(kvmatrix) / Va^2'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./chii]
    type = DerivativeParsedMaterial
    f_name = chii
    material_property_names = 'Va kimatrix '
    function = '(kimatrix) / Va^2'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./chig]
    type = DerivativeParsedMaterial
    f_name = chig
    material_property_names = 'Va kgmatrix '
    function = '(kgmatrix ) / Va^2'
    derivative_order = 2
    # outputs = exodus
  [../]

  # Dislocation density per nm2
  [./dislocation_density]
    type = ParsedMaterial
    f_name = dislocation_density
    args = 'burnup_var'
    function = '(10^(2.2e-2*burnup_var+13.8))*1e-18' #This is an empirical expression  https://www.sciencedirect.com/science/article/pii/S002231151200637X
    # outputs = exodus
  [../]

[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
 # time step
  [./dt]
    type = TimestepSize
  [../]
  [./total_free_energy]
    type = ElementIntegralVariablePostprocessor
    variable = local_energy
  [../]
[]

[Executioner]
  type = Transient
  nl_max_its = 15
  scheme = bdf2
  #solve_type = NEWTON
   solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type'
  petsc_options_value = 'asm       1               ilu'
  l_max_its = 15
  l_tol = 1.0e-4
  nl_rel_tol = 1.0e-8
  start_time = 0.0
  #num_steps = 1500
  end_time = 63072000
  dtmax = 1e6
  nl_abs_tol = 1e-10
  [./TimeStepper]
    type = IterationAdaptiveDT
        dt = 1000
        optimal_iterations = 8
        iteration_window = 2
  [../]
[]

[UserObjects]
  active = 'rasterizer runner neutronics_fission_generator'
  [./neutronics_fission_generator]
    type = PKAFissionFragmentEmpirical
    relative_density = '1'
    fission_rate = 1.5e-08
  [../]
  [./xenon]
    type = PKAConstant
    pka_rate = 3e-8
    m = 131
    Z = 54
    E = 70e6
  [../]
  [./rasterizer]
    type = MyTRIMRasterizer
    var = 'rho_U235  rho_U238   rho_O  rhog_var rhoi_var'
    M   = '235       238        16      135     238'
    Z   = '92        92         8       54      92'
    site_volume = 1 # nm^3 per UO2 unit
    periodic_var = wv
    pka_generator = neutronics_fission_generator
    length_unit = NANOMETER
    max_pka_count = 1000
    recoil_rate_scaling = 1
    r_rec = 5.45
  [../]
  [./runner]
    type = MyTRIMElementRun
    rasterizer = rasterizer
  [../]
[]

[Adaptivity]
 marker = errorfrac
 max_h_level = 2
 [./Indicators]
   [./error]
     type = GradientJumpIndicator
     variable = bnds
   [../]
 [../]
 [./Markers]
   [./errorfrac]
     type = ErrorFractionMarker
     coarsen = 0.1
     indicator = error
     refine = 0.7
   [../]
 [../]
[]

[Outputs]
  exodus = true
  csv = true
[]
