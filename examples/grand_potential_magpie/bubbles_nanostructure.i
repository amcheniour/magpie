# This is based on Larry's grandpotential.i file in marmot-problems/Magpie/examples
# The domain is 2D 400 nm x 400 nm with a xenon bubble in the middle of UO2 matrix, radius = 44 nm
# The variables are chemical potentials of U vacancies and interstitials, Xe atoms
# This is solving the rate theory equation for wi, wv, wg
# It is using a source term from MAGPIE based on BCMC

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 30
  ny = 30
  xmin = 0
  xmax = 400
  ymin = 0
  ymax = 400
  # uniform_refine = 2
[]

[GlobalParams]
  op_num = 2
  var_name_base = etam
  bubspac = 70
  numbub = 8
  rand_seed = 3147315
[]

[Variables]
  # chemical potentials for vacancies, interstitials and gas atoms
  [./wv]
  [../]
  [./wi]
  [../]
  [./wg]
  [../]

  # order parameters: etab0 is the bubble, etam0 is the grain, etam1 is currently not used
  [./etab0]
  [../]
  [./etam0]
  [../]
  [./etam1]
  [../]
[]

[AuxVariables]
  # grain boundaries
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]

  # concentration of Xe
  [./c_Xe]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
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


  # concentrations of (other) defects
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

  # number of defects generated each time step
  [./interstitial_rate_Xe]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vacancy_rate_U]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./interstitial_rate_U]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./interstitial_rate_i]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # burnup in GWd/tHM
  [./burnup_var]
    order = CONSTANT
    family = MONOMIAL
  [../]

  # temperature
  [./T]
    initial_condition = 1200
  [../]


[]

[ICs]
  [./IC_etab0]
    type = MultiSmoothSuperellipsoidIC
    variable = etab0
    invalue = 1
    outvalue = 0
    semiaxis_a = 30
    semiaxis_b = 30
    semiaxis_c = 1
    exponent = 2
    semiaxis_a_variation = 0
    semiaxis_b_variation = 0
    semiaxis_c_variation = 0
    prevent_overlap = true
    vary_axes_independently = false
    int_width = 40

  [../]
  [./IC_etam0]
    type = MultiSmoothSuperellipsoidIC
    variable = etam0
    #bubspac = 100
    outvalue = 1
    invalue = 0
    semiaxis_a = 30
    semiaxis_b = 30
    semiaxis_c = 1
    exponent = 2
    semiaxis_a_variation = 0
    semiaxis_b_variation = 0
    semiaxis_c_variation = 0
    prevent_overlap = true
    vary_axes_independently = false
    int_width = 40
  [../]
  [./IC_etam1]
    type = MultiSmoothSuperellipsoidIC
    variable = etam1
    #bubspac = 100
    invalue = 0
    outvalue = 0
    semiaxis_a = 30
    semiaxis_b = 30
    semiaxis_c = 1
    exponent = 2
    semiaxis_a_variation = 0
    semiaxis_b_variation = 0
    semiaxis_c_variation = 0
    prevent_overlap = true
    vary_axes_independently = false
    int_width = 40
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

[Functions]
  # bubble IC
  [./ic_func_etab0]
    type = ParsedFunction
    vars = 'kappa   mu'
    vals = '0.5273  0.004688'
    value = 'r:=sqrt((x-200)^2+(y-200)^2);0.5*(1.0-tanh((r-44)*sqrt(mu/2.0/kappa)))'
  [../]
  # grain IC
  [./ic_func_etam0]
    type = ParsedFunction
    vars = 'kappa   mu'
    vals = '0.5273  0.004688'
    value = 'r:=sqrt((x-200)^2+(y-200)^2);0.5*(1.0+tanh((r-44)*sqrt(mu/2.0/kappa)))'
  [../]
  [./ic_func_etam1]
    type = ParsedFunction
    vars = 'kappa   mu'
    vals = '0.5273  0.004688'
    value = '0'
  [../]

  # burnup function
  [./fburnup] # in GWD/tHM
    type = ParsedFunction
    value = 'q*sigmaf_u235*flux*t*950'
    vars = 'q sigmaf_u235 flux'
    vals = '0.045 5.906e-22 1.5e13'
  [../]

[]

[Modules]
  [./PhaseField]
    [./GrandPotential]
      switching_function_names = 'hb hm'

      chemical_potentials = 'wv wi wg'
      mobilities = 'Dchiv Dchii Dchig'
      susceptibilities = 'chiv chii chig'
      free_energies_w = ' rhovbub rhovmatrix rhoibub rhoimatrix rhogbub rhogmatrix'

      gamma_gr = 'gmm'
      op_num = 2
      var_name_base = etam
      free_energies_gr = 'omegab omegam'
      kappa_gr = kappa
      mobility_name_gr = L

      additional_ops = etab0
      free_energies_op = 'omegab omegam'
      mobility_name_op = L
      kappa_op = kappa_op
      gamma_op = gmb
      gamma_grxop = gmb
    [../]
  [../]
[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x y' # 2D
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
 # [./recombination_v] # Recombination term = Recombination rate * rhoi * rhov
 #   type = GrandPotentialRecombination
 #   variable = wv
 #   rho = rhov
 #   rho_r = rhoi
 #   value = 1
 #   D = Di
 #   omega = 0.04092
 #   a0 = 0.25
 #   Z = 50
 #   hm = hm
 #   args = 'wi etab0 etam0'
 # [../]
 # [./dislocation_sink_v] # Dislocation sink = sink strength * Diffusion coefficient * rho
 #    type = GrandPotentialSink
 #    variable = wv
 #    rho = rhov
 #    D = D
 #    mask = hm
 #    value = 1
 #    sink_strength = dislocation_density
 #    args = 'etab0 etam0 etam1'
 #[../]

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
  #[./recombination_i]
  #  type = GrandPotentialRecombination
  #  variable = wi
  #  rho = rhoi
  #  rho_r = rhov
  #  value = 1
  #  D = Di
  #  omega = 0.04092
  #  a0 = 0.25
  #  Z = 50
  #  hm = hm
  #  args = 'wv etab0 etam0'
  #[../]
  #[./dislocation_sink_i]
  #   type = GrandPotentialSink
  #   variable = wi
  #   rho = rhoi
  #   mask = hm
  #   D = Di
  #   value = 1.02
  #   sink_strength = dislocation_density
  #   args = 'etab0 etam0 etam1'
  # [../]
  [./source_g]
    type = MyTRIMElementSource
    variable = wg
    runner = runner
    ivar = 3
    defect = INT
    prefactor = 0.04092
  [../]
[]

[AuxKernels]
  # Numbers of defects created each time step
  [./interstitial_rate_Xe]
    type = MyTRIMElementResultAux
    runner = runner
    defect = INT
    variable = interstitial_rate_Xe
    ivar = 3
    execute_on = timestep_end
  [../]
  [./vacancy_rate_U]
    type = MyTRIMElementResultAux
    runner = runner
    defect = VAC
    variable = vacancy_rate_U
    ivar = 1
    execute_on = timestep_end
  [../]
  [./interstitial_rate_U]
    type = MyTRIMElementResultAux
    runner = runner
    defect = INT
    variable = interstitial_rate_U
    ivar = 1
    execute_on = timestep_end
  [../]
  [./interstitial_rate_i]
    type = MyTRIMElementResultAux
    runner = runner
    defect = INT
    variable = interstitial_rate_i
    ivar = 4
    execute_on = timestep_end
  [../]

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
  # Auxkernels for sink rates
  # The value of absorbed_int/vac is integrated over the domain to get the total number of atoms absorbed (see postprocessor)
  # It is also multiplied by the time step because the absorption rate is per second
  #[./disloc_absorption_int]
  #  type = SinkAbsorptionRate
  #  variable = absorbed_int
  #  rho = rhoi
  #  D = Di
  #  mask = hm
  #  sink_strength = disloc_density_dt
  #[../]
  #[./disloc_absorption_vac]
  #  type = SinkAbsorptionRate
  #  variable = absorbed_vac
  #  rho = rhov
  #  D = D
  #  mask = hm
  #  sink_strength = disloc_density_dt
  #[../]

  # Burnup auxkernel defined by the function fburnup
  [./burnup_aux]
    type = FunctionAux
    function = fburnup
    variable = burnup_var
  [../]
[]

[Materials]

  # Switching functions for bubble and matrix
  [./hb]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = hb
    all_etas = 'etab0 etam0 etam1'
    phase_etas = 'etab0'
    # outputs = exodus
  [../]
  [./hm]
    type = SwitchingFunctionMultiPhaseMaterial
    h_name = hm
    all_etas = 'etab0 etam0 etam1'
    phase_etas = 'etam0 etam1'
    outputs = exodus
  [../]

  # Grand potential densities for bubble and matrix
  [./omegab]
    type = DerivativeParsedMaterial
    args = 'wv wi wg'
    f_name = omegab
    material_property_names = 'Va kvbub cvbubeq kibub cibubeq kgbub cgbubeq f0'
    function = '-0.5*wv^2/Va^2/kvbub-wv/Va*cvbubeq-0.5*wi^2/Va^2/kibub-wi/Va*cibubeq-0.5*wg^2/Va^2/kgbub-wg/Va*cgbubeq+f0'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./omegam]
    type = DerivativeParsedMaterial
    args = 'wv wi wg'
    f_name = omegam
    material_property_names = 'Va kvmatrix cvmatrixeq kimatrix cimatrixeq kgmatrix cgmatrixeq'
    function = '-0.5*wv^2/Va^2/kvmatrix-wv/Va*cvmatrixeq-0.5*wi^2/Va^2/kimatrix-wi/Va*cimatrixeq-0.5*wg^2/Va^2/kgmatrix-wg/Va*cgmatrixeq'
    derivative_order = 2
    #outputs = exodus
  [../]

  # Number densities defind in each phase for each defect (vacancy, interstitial, gas atoms)
  [./rhovbub]
    type = DerivativeParsedMaterial
    args = 'wv'
    f_name = rhovbub
    material_property_names = 'Va kvbub cvbubeq'
    function = 'wv/Va^2/kvbub + cvbubeq/Va'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./rhovmatrix]
    type = DerivativeParsedMaterial
    args = 'wv'
    f_name = rhovmatrix
    material_property_names = 'Va kvmatrix cvmatrixeq'
    function = 'wv/Va^2/kvmatrix + cvmatrixeq/Va'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./rhoibub]
    type = DerivativeParsedMaterial
    args = 'wi'
    f_name = rhoibub
    material_property_names = 'Va kibub cibubeq'
    function = 'wi/Va^2/kibub + cibubeq/Va'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./rhoimatrix]
    type = DerivativeParsedMaterial
    args = 'wi'
    f_name = rhoimatrix
    material_property_names = 'Va kimatrix cimatrixeq'
    function = 'wi/Va^2/kimatrix + cimatrixeq/Va'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./rhogbub]
    type = DerivativeParsedMaterial
    args = 'wg'
    f_name = rhogbub
    material_property_names = 'Va kgbub cgbubeq'
    function = 'wg/Va^2/kgbub + cgbubeq/Va'
    derivative_order = 2
    #outputs = exodus
  [../]
  [./rhogmatrix]
    type = DerivativeParsedMaterial
    args = 'wg'
    f_name = rhogmatrix
    material_property_names = 'Va kgmatrix cgmatrixeq'
    function = 'wg/Va^2/kgmatrix + cgmatrixeq/Va'
    derivative_order = 2
    #outputs = exodus
  [../]

  # Total number densities in the whole domain for each defect
  [./rhov]
    type = DerivativeParsedMaterial
    args = 'wv'
    f_name = rhov
    material_property_names = 'hm hb rhovmatrix(wv) rhovbub(wv)'
    function = 'hm * rhovmatrix + hb * rhovbub'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./rhoi]
    type = DerivativeParsedMaterial
    args = 'wi'
    f_name = rhoi
    material_property_names = 'hm hb rhoimatrix(wi) rhoibub(wi)'
    function = 'hm * rhoimatrix + hb * rhoibub'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./rhog]
    type = DerivativeParsedMaterial
    args = 'wg'
    f_name = rhog
    material_property_names = 'hm hb rhogmatrix(wg) rhogbub(wg)'
    function = 'hm * rhogmatrix + hb * rhogbub'
    derivative_order = 2
    # outputs = exodus
  [../]

  # U and O number densities
  [./rho_U235_matl]
    type = ParsedMaterial
    material_property_names = c_U235_matl
    f_name = rho_U235_matl
    constant_names = Va
    constant_expressions = 0.04092
    function = 'c_U235_matl/Va'
    # outputs = exodus
  [../]
  [./rho_U238_matl]
    type = ParsedMaterial
    material_property_names = c_U238_matl
    f_name = rho_U238_matl
    constant_names = Va
    constant_expressions = 0.04092
    function = 'c_U238_matl/Va'
    # outputs = exodus
  [../]
  [./rho_O_matl]
    type = ParsedMaterial
    material_property_names = 'c_O_matl'
    f_name = rho_O_matl
    constant_names = Va
    constant_expressions = 0.04092
    function = 'c_O_matl/Va/2'
    # outputs = exodus
  [../]

  # Constants
  [./const]
    type = GenericConstantMaterial
    prop_names =  'kappa   mu       L   D    Va      cvbubeq cgbubeq cibubeq  kgbub  kvbub kibub gmb     gmm T    Efvbar    Efgbar    kTbar     f0     tgrad_corr_mult  kappa_c kappa_op gamma_asymm Di'
    prop_values = '0.5273  0.004688 0.1 0.01 0.04092 0.5459  0.4541  0.0      1.41   1.41  1.41  0.9218 1.5 1200 7.505e-3  7.505e-3  2.588e-4  0.0    0.0              1.0     0.5273   1.5         1 '
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

  # Material properties to define the number fractions of U, O, Xe, U vacancy and interstitial
  [./c_U235_matl]
    type = ParsedMaterial
    f_name = c_U235_matl
    material_property_names = 'hm'
    function = 'hm*0.01485'
    outputs = exodus
  [../]
  [./c_U238_matl]
    type = ParsedMaterial
    f_name = c_U238_matl
    material_property_names = 'hm'
    function = 'hm*0.31515'
    outputs = exodus
  [../]
  [./c_O_matl]
    type = ParsedMaterial
    f_name = c_O_matl
    material_property_names = 'hm'
    function = '0.66*hm'
    outputs = exodus
  [../]
  [./c_Xe_matl]
    type = ParsedMaterial
    f_name = c_Xe_matl
    material_property_names = 'hm hb Va rhogmatrix rhogbub'
    function = 'Va * (hm * rhogmatrix + hb * rhogbub)'
  [../]
  [./cv_mat]
    type = ParsedMaterial
    f_name = cv_mat
    material_property_names = 'hm hb Va rhovmatrix rhovbub'
    function = 'Va * (hm * rhovmatrix + hb * rhovbub)'
    outputs = exodus
  [../]
  [./ci_mat]
    type = ParsedMaterial
    f_name = ci_mat
    material_property_names = 'hm hb Va rhoimatrix rhoibub'
    function = 'Va * (hm * rhoimatrix + hb * rhoibub)'
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
    material_property_names = 'Va hb kvbub hm kvmatrix '
    function = '(hm/kvmatrix + hb/kvbub) / Va^2'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./chii]
    type = DerivativeParsedMaterial
    f_name = chii
    material_property_names = 'Va hb kibub hm kimatrix '
    function = '(hm/kimatrix + hb/kibub) / Va^2'
    derivative_order = 2
    # outputs = exodus
  [../]
  [./chig]
    type = DerivativeParsedMaterial
    f_name = chig
    material_property_names = 'Va hb kgbub hm kgmatrix '
    function = '(hm/kgmatrix + hb/kgbub) / Va^2'
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
  # Dislocation density per nm2 * time step in sec
  [./disloc_dens_dt]
    type = ParsedMaterial
    f_name = disloc_density_dt
    material_property_names = 'dislocation_density dt'
    function = 'dislocation_density * dt'
    # outputs = exodus
  [../]



  # Time step material
  [./dt]
    type = TimeStepMaterial
  [../]

[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Postprocessors]
  # total number of defects in each time step
  [./total_Xe_interstitial_production]
    type = ElementIntegralVariablePostprocessor
    variable = interstitial_rate_Xe
  [../]

  [./total_U_vacancy_production]
    type = ElementIntegralVariablePostprocessor
    variable = vacancy_rate_U
  [../]

  [./total_U_interstitial_production]
    type = ElementIntegralVariablePostprocessor
    variable = interstitial_rate_U
  [../]

  [./total_i_interstitial_production]
    type = ElementIntegralVariablePostprocessor
    variable = interstitial_rate_i
  [../]
  #
  # [./total_rep_in_235_production]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = rep_in_235
  # [../]
  #
  # [./total_rep_out_235_production]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = rep_out_235
  # [../]
  #
  # [./total_rep_in_238_production]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = rep_in_238
  # [../]
  #
  # [./total_rep_out_238_production]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = rep_out_238
  # [../]

 # average number fraction of defect
  [./cv_average]
    type = ElementAverageValue
    variable = cv
  [../]
  [./ci_average]
    type = ElementAverageValue
    variable = ci
  [../]
  [./cg_average]
    type = ElementAverageValue
    variable = c_Xe
  [../]

 # initial matrix volume
  [./UO2_volume]
    type = ElementIntegralMaterialProperty
    mat_prop = hm
    execute_on = initial
  [../]

 # time step
  [./dt]
    type = TimestepSize
  [../]

[]

[Executioner]
  type = Transient
  nl_max_its = 15
  scheme = bdf2
  solve_type = NEWTON
  # solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_asm_overlap -sub_pc_type'
  petsc_options_value = 'asm       1               ilu'
  l_max_its = 15
  l_tol = 1.0e-4
  nl_rel_tol = 1.0e-8
  start_time = 0.0
  num_steps = 1500
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
  file_base = density
[]
