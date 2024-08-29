!------------------------------------------------------------------------------!
!                              PincFlow Namelists                              !
!------------------------------------------------------------------------------!

&domain
  sizeX                      =                      32 ! Cells in x
  sizeY                      =                       1 ! Cells in y
  sizeZ                      =                     100 ! Cells in z
  nbx                        =                       3 ! Halo/ghost cells in x
  nby                        =                       3 ! Halo/ghost cells in y
  nbz                        =                       3 ! Ghost cells in z
  lx_dim(0)                  =                  0.0E+0 ! Minimum of x
  lx_dim(1)                  =                  9.0E+6 ! Maximum of x
  ly_dim(0)                  =                  0.0E+0 ! Minimum of y
  ly_dim(1)                  =                  3.0E+5 ! Maximum of y
  lz_dim(0)                  =                  0.0E+0 ! Minimum of z
  lz_dim(1)                  =                  1.0E+5 ! Maximum of z
  nprocx                     =                {nprocx} ! Processors in x
  nprocy                     =                {nprocy} ! Processors in y
&end

&variables
  include_ice                =                       F ! Ice microphysics
                                                       ! parameterization
  include_ice2               =                       F ! ...
  include_tracer             =                       T ! Tracer equation
  include_testoutput         =                       F ! ...
&end

&outputList
  dimOut(1)                  =                       T ! ...
  dimOut(2)                  =                       T ! ...
  dimOut(3)                  =                       T ! ...
  varOut(1)                  =                   'rho' ! ...
  varOut(2)                  =                     'u' ! ...
  varOut(3)                  =                     'v' ! ...
  varOut(4)                  =                     'w' ! ...
  varOut(5)                  =                    'pi' ! ...
  varOut(6)                  =                   'chi' ! ...
  varIn(1)                   =                   'rho' ! ...
  varIn(2)                   =                     'u' ! ...
  varIn(3)                   =                     'v' ! ...
  varIn(4)                   =                     'w' ! ...
  varIn(5)                   =                    'pi' ! ...
  varIn(6)                   =                   'chi' ! ...
  offset(1)                  =                  0.0E+0 ! ...
  offset(2)                  =                  0.0E+0 ! ...
  offset(3)                  =                  0.0E+0 ! ...
  offset(4)                  =                  0.0E+0 ! ...
  offset(5)                  =                  0.0E+0 ! ...
  offset(6)                  =                  0.0E+0 ! ...
  iIn                        =                      23 ! ...
  optVarOut(1)               =                       0 ! ...
  optVarOut(2)               =                       0 ! ...
  optVarOut(3)               =                       0 ! ...
  optVarOut(4)               =                       0 ! ...
  optVarOut(5)               =                       0 ! ...
  optVarOut(6)               =                       0 ! ...
  wkbVarOut(1)               =                       0 ! ...
  wkbVarOut(2)               =                       0 ! ...
  wkbVarOut(3)               =                       0 ! ...
  wkbVarOut(4)               =                       0 ! ...
  wkbVarOut(5)               =                       0 ! ...
  wkbVarOut(6)               =                       0 ! ...
  wkbVarOut(7)               =                       0 ! ...
  wkbVarOut(8)               =                       0 ! ...
  wkbVarOut(9)               =                       0 ! ...
  wkbVarOut(10)              =                       0 ! ...
  wkbVarOut(11)              =                       0 ! ...
  wkbVarOut(12)              =                       0 ! ...
  wkbVarOut(13)              =                       0 ! ...
  outputType                 =                  'time' ! 'timeStep' or 'time'
  nOutput                    =                       1 ! Output every nOutput
                                                       ! time steps for
                                                       ! outputType =
                                                       ! 'timeStep'
  maxIter                    =                       1 ! Stop after maxIter
                                                       ! time steps for
                                                       ! outputType =
                                                       ! 'timeStep'
  outputTimeDiff             =                  9.0E+3 ! Output every
                                                       ! outputTimeDiff seconds
                                                       ! for outputType =
                                                       ! 'time'
  maxTime                    =                  9.0E+4 ! Stop after maxTime
                                                       ! seconds for outputType
                                                       ! = 'time'
  restart                    =                       F ! Restart the model from
                                                       ! state in previous
                                                       ! simulation
  thetaOffset                =                       F ! Subtract background
                                                       ! potential temperature
  rhoOffset                  =                       T ! Subtract background
                                                       ! density
  detailedinfo               =                       F ! Provide info on the
                                                       ! final state of Poisson
                                                       ! solver
  RHS_diagnostics            =                       T ! Provide info about the
                                                       ! RHS of Poisson
                                                       ! equation
  fancy_namelists            =                       T ! Write all namelists
                                                       ! with comments
&end

&debuggingList
  verbose                    =                       F ! ...
  dtMin_dim                  =                  1.0E-5 ! Stop if dt < dtMin
&end

&testCaseList
  testCase                   =             'raytracer' ! Predefined test case
&end

&monochromeWave
  lambda_dim                 =                  6.0E+3 ! Vertical wavelength
&end

&wavePacket
  wavePacketType             =                       1 ! 1 = Gaussian, 2 =
                                                       ! Cosine
  wavePacketDim              =                       1 ! 1 = 1D, 2 = 2D, 3 = 3D
  lambdaX_dim                =                  1.0E+3 ! Wavelength in x (0.0
                                                       ! for infinite
                                                       ! wavelength)
  lambdaY_dim                =                  0.0E+0 ! Wavelength in y (0.0
                                                       ! for infinite
                                                       ! wavelength)
  lambdaZ_dim                =                  1.0E+3 ! Wavelength in z
  amplitudeFactor            =                  9.0E-1 ! Normalized buoyancy
                                                       ! amplitude
  x0_dim                     =                  5.0E+2 ! Wave-packet center in x
  y0_dim                     =                  1.5E+4 ! Wave-packet center in y
  z0_dim                     =                  1.0E+4 ! Wave-packet center in z
  sigma_dim                  =                  2.0E+3 ! Vertical width of
                                                       ! Gaussian wave packet
  sigma_hor_dim              =                  0.0E+0 ! Cosine distribution
                                                       ! width in x (0.0 for
                                                       ! infinite width)
  amp_mod_x                  =                  1.0E+0 ! Fractional amplitude
                                                       ! modulation in x
  sigma_hor_yyy_dim          =                  0.0E+0 ! Cosine distribution
                                                       ! width in y (0.0 for
                                                       ! infinite width)
  amp_mod_y                  =                  1.0E+0 ! Fractional amplitude
                                                       ! modulation in y
  L_cos_dim                  =                  2.0E+4 ! Half width of vertical
                                                       ! cosine profile of wave
                                                       ! packet
  omiSign                    =                      -1 ! Frequency branch
  u0_jet_dim                 =                  0.0E+0 ! Maximum amplitude of
                                                       ! jet
  z0_jet_dim                 =                  5.0E+4 ! Center of jet
  L_jet_dim                  =                  5.0E+3 ! Half width of vertical
                                                       ! cosine profile of jet
  inducedwind                =                       F ! ...
&end

&LagrangeRayTracing
  xrmin_dim                  =                  0.0E+0 ! Left bound of initial
                                                       ! rays
  xrmax_dim                  =                  9.0E+6 ! Right bound of initial
                                                       ! rays
  yrmin_dim                  =                  0.0E+0 ! Backward bound of
                                                       ! initial rays
  yrmax_dim                  =                  3.0E+5 ! Forward bound of
                                                       ! initial rays
  zrmin_dim                  =                  0.0E+0 ! Bottom bound of
                                                       ! initial rays
  zrmax_dim                  =                  1.0E+5 ! Top bound of initial
                                                       ! rays
  nrxl                       =                      26 ! Initial ray volumes
                                                       ! within dx
  nryl                       =                       1 ! Initial ray volumes
                                                       ! within dy
  nrzl                       =                      26 ! Initial ray volumes
                                                       ! within dz
  fac_dk_init                =                  1.0E-1 ! Initial fraction dk/k
  fac_dl_init                =                  1.0E-1 ! Initial fraction dl/l
  fac_dm_init                =                  1.0E-4 ! Initial fraction dm/m
  nrk_init                   =                       1 ! Initial ray volumes
                                                       ! within dk
  nrl_init                   =                       1 ! Initial ray volumes
                                                       ! within dl
  nrm_init                   =                       1 ! Initial ray volumes
                                                       ! within dm
  nsmth_wkb                  =                       2 ! Width of smoothing
                                                       ! operator for mean flow
                                                       ! tendencies
  lsmth_wkb                  =                       T ! Smoothing operator for
                                                       ! mean flow tendencies
  sm_filter                  =                       2 ! 1 = Box filter, 2 =
                                                       ! Shapiro filter
  lsaturation                =                       F ! Switch for saturation
                                                       ! scheme
  alpha_sat                  =                  1.4E+0 ! Saturation threshold
  steady_state               =                       F ! Steady-state mode
  case_wkb                   =                       1 ! 1 = Gaussian wave
                                                       ! packet, 2 = cosine
                                                       ! wave packet, 3 =
                                                       ! mountain wave
  amp_wkb                    =                  5.0E-1 ! Relative amplitude of
                                                       ! wave packet
  wlrx_init                  =                  0.0E+0 ! Wavelength in x of
                                                       ! wave packet
  wlry_init                  =                  3.0E+5 ! Wavelength in y of
                                                       ! wave packet
  wlrz_init                  =                  1.0E+3 ! Wavelength in z of
                                                       ! wave packet
  xr0_dim                    =                  4.5E+6 ! Wave packet center in x
  yr0_dim                    =                  1.5E+4 ! Wave packet center in y
  zr0_dim                    =                  3.0E+4 ! Wave packet center in z
  sigwpx_dim                 =                  1.5E+6 ! Wave packet width in x
                                                       ! (0.0 for infinite
                                                       ! width)
  sigwpy_dim                 =                  0.0E+0 ! Wave packet width in y
                                                       ! (0.0 for infinite
                                                       ! width)
  sigwpz_dim                 =                  5.0E+3 ! Wave packet width in z
                                                       ! (0.0 for infinite
                                                       ! width)
  branchr                    =                       1 ! Frequency branch
  lindUinit                  =                       F ! Induced wind at
                                                       ! initial time
  topographyTime_wkb         =                  0.0E+0 ! WKB topography growth
                                                       ! time
  mountainHeight_wkb_dim     =                  0.0E+0 ! WKB mountain height
  mountainWidth_wkb_dim      =                  1.0E+0 ! WKB mountain half width
  mountain_case_wkb          =                       6 ! WKB topography shape
  range_factor_wkb           =                  1.0E+1 ! Ratio between large
                                                       ! and small WKB scales
  blocking                   =                       F ! Simple blocked-layer
                                                       ! scheme
  nwm                        =                       1 ! Number of initial wave
                                                       ! modes
  launch_level               =                       0 ! Ray-volume launch level
  launch_algorithm           =                  'clip' ! Ray-volume launch
                                                       ! algorithm
  zmin_wkb_dim               =                  0.0E+0 ! Minimum altitude for
                                                       ! wave-mean-flow
                                                       ! interaction
  nray_fac                   =                       1 ! Maximum multiplication
                                                       ! factor per spectral
                                                       ! dimension
  cons_merge                 =                    'wa' ! Conserved quantity in
                                                       ! ray-volume merging
                                                       ! ('wa' = wave action,
                                                       ! 'en' = wave energy)
  nRayOutput                 =                      10 ! Number of dominant ray
                                                       ! volumes in output
&end

&bubble
  dTheta0_dim                =                  6.0E+0 ! ...
  xRadius_dim                =                  2.5E+3 ! ...
  zRadius_dim                =                  2.5E+3 ! ...
  xCenter_dim                =                  0.0E+0 ! ...
  yCenter_dim                =                  1.5E+5 ! ...
  zCenter_dim                =                  0.0E+0 ! ...
  zExcentricity              =                  1.0E+0 ! ...
  rhoCenter_dim              =                  1.0E+3 ! ...
&end

&robert_bubble
  dTheta1_dim                =                  5.0E-1 ! Potential temperature
                                                       ! offset
  a1_dim                     =                  1.5E+2 ! Radius of plateau
  sigma1_dim                 =                  5.0E+1 ! Gaussian edge profile
  xCenter1_dim               =                  0.0E+0 ! ...
  zCenter1_dim               =                  2.5E+4 ! ...
  dTheta2_dim                =                 -1.5E-1 ! ...
  a2_dim                     =                  0.0E+0 ! ...
  sigma2_dim                 =                  5.0E+1 ! ...
  xCenter2_dim               =                  6.0E+1 ! ...
  zCenter2_dim               =                  6.4E+2 ! ...
&end

&mountainwavelist
  u_relax                    =                  7.5E+1 ! Zonal relaxation wind
  v_relax                    =                  0.0E+0 ! Meridional relaxation
                                                       ! wind
  w_relax                    =                  0.0E+0 ! Vertical relaxation
                                                       ! wind
  t_relax                    =                1.728E+5 ! Relaxation time
  t_ramp                     =                  3.6E+3 ! Not used at the moment
  xextent_relax              =                  3.0E+6 ! Zonal extent of
                                                       ! relaxation region
  yextent_relax              =                  0.0E+0 ! Meridional extent of
                                                       ! relaxation region
  wind_relaxation            =                       F ! Relaxation switch
  surface_layer_depth        =                  0.0E+0 ! Surface-layer depth
&end

&baroclinic_LC
  zero_initial_state         =                       F ! ...
  z_trpp0_dim                =                  8.0E+3 ! Mean tropopause height
  z_baro_dim                 =                  2.4E+4 ! Altitude above which
                                                       ! the atmosphere is
                                                       ! barotropic
  thet0_dim                  =                 2.73E+2 ! Characteristic
                                                       ! potential temperature
  ntrp_dim                   =                  1.0E-2 ! Buoyancy frequency of
                                                       ! the troposphere
  nstr_dim                   =                 2.45E-2 ! Buoyancy frequency of
                                                       ! the stratosphere
  jwdth_dim                  =                 8.25E+5 ! Jet half width
  kaptpp                     =                  7.0E+1 ! Jet slope
  add_ptptb                  =                       T ! Add local potential-
                                                       ! temperature
                                                       ! perturbation
  ptptb_x_dim                =                  5.0E+6 ! Location in x of local
                                                       ! potential-temperature
                                                       ! perturbation
  ptptb_y_dim                =                  4.0E+6 ! Location in y of local
                                                       ! potential-temperature
                                                       ! perturbation
  ptptb_z_dim                =                  9.0E+3 ! Location in z of local
                                                       ! potential-temperature
                                                       ! perturbation
  ptptb_dh_dim               =                  5.0E+5 ! Horizontal width of
                                                       ! local potential-
                                                       ! temperature
                                                       ! perturbation
  ptptb_dz_dim               =                  2.0E+3 ! Vertical width of
                                                       ! local potential-
                                                       ! temperature
                                                       ! perturbation
  ptptb_amp_dim              =                  3.0E+0 ! Amplitude of local
                                                       ! potential-temperature
                                                       ! perturbation
  add_noise                  =                       F ! Add noise to the
                                                       ! initial potential
                                                       ! temperature
  proc_noise                 =                  1.0E-2 ! Relative amplitude of
                                                       ! potential-temperature
                                                       ! noise
  tau_relax                  =                1.296E+6 ! ...
  tau_relax_low              =                  0.0E+0 ! ...
  sigma_tau                  =                  0.0E+0 ! Relative thickness of
                                                       ! relaxation profile
  tau_jet                    =                  7.2E+3 ! Time for jet formation
  Sponge_Rel_Bal_Type        =                   'env' ! Relaxation state:
                                                       ! 'hyd' = hydrostatic
                                                       ! balance, 'env' =
                                                       ! geostrophic balance
  ta_hs_dim                  =                  0.0E+0 ! Thermal-relaxation
                                                       ! time scale outside
                                                       ! tropical boundary
                                                       ! layer (zero means
                                                       ! infinity)
  ts_hs_dim                  =                  0.0E+0 ! Thermal-relaxation
                                                       ! time scale in tropical
                                                       ! boundary layer (zero
                                                       ! means infinity)
  tf_hs_dim                  =                  0.0E+0 ! Boundary-layer
                                                       ! Rayleigh-damping  time
                                                       ! scale
  sigb_hs                    =                  0.0E+0 ! Sigma of boundary-
                                                       ! layer top
&end

&baroclinic_ID
  bar_sigma_y                =                  1.0E-1 ! ...
  u_strength                 =                  3.0E+1 ! Jet strength
  dTh_atm                    =                  2.0E+1 ! Potential-temperature
                                                       ! difference between
                                                       ! poles and tropics
  init_2Dto3D                =                       F ! Initialize 3D state
                                                       ! from 2D balanced state
  init_bal                   =             'geostr_id' ! ...
  lastrecordnum              =                      88 ! Last record in 2D
                                                       ! output
  fileinitstate2D            =    'pf_all_2D3D_in.dat' ! File name with 2D
                                                       ! balanced state
  output_theta_bgr           =                       F ! Output environmental
                                                       ! potential temperature
  output_br_vais_sq          =                       T ! Output environmental
                                                       ! squared buoyancy
                                                       ! frequency
  output_heat                =                       T ! ...
  balance_eq                 =                 'PI_an' ! ...
  output_rho_bgr             =                       T ! Output environmental
                                                       ! density
&end

&modelList
  model                      = 'pseudo_incompressible' ! Dynamic equations
  vert_theta                 =                  9.0E+1 ! Rotation about x
  vert_alpha                 =                  0.0E+0 ! Rotation about y
&end

&solverList
  cfl                        =                  5.0E-1 ! CFL number
  cfl_wave                   =                  2.5E-1 ! WKB CFL number
  dtMax_dim                  =                  1.0E+2 ! Maximum time step
  tStepChoice                =                   'cfl' ! 'fix' or 'cfl'
  timeScheme                 =          'semiimplicit' ! 'LS_Will_RK3' or
                                                       ! 'semiimplicit'
  auxil_equ                  =                       F ! Buoyancy equation
  fluxType                   =                'upwind' ! 'ILES', 'central' or
                                                       ! 'upwind'
  reconstType                =                 'MUSCL' ! 'MUSCL', 'constant',
                                                       ! 'SALD' or 'ALDM'
  musclType                  =                'muscl1' ! 'muscl1' or 'muscl2'
  limiterType1               =             'MCVariant' ! 'minmod', 'MCVariant'
                                                       ! or 'Cada'
  fluctuationMode            =                       T ! Subtract background
                                                       ! density
  TurbScheme                 =                       F ! Turbulence scheme
  turb_dts                   =                  5.0E+3 ! Turbulent damping time
  DySmaScheme                =                       T ! Dynamic Smagorinsky
                                                       ! scheme
  dtWave_on                  =                       T ! Limit time step by
                                                       ! inverse buoyancy
                                                       ! frequency
  heatingONK14               =                       F ! Heating as in O'Neill
                                                       ! and Klein (2014)
  dens_relax                 =                       F ! Heating by density
                                                       ! relaxation
  shap_dts_fac               =                 -1.0E+0 ! Shapiro-filter damping
                                                       ! time
  n_shap                     =                       1 ! Order of Shapiro filter
&end

&poissonSolverList
  tolPoisson                 =                  1.0E-4 ! Abort criterion
  abs_tol                    =                  0.0E+0 ! Lower bound for
                                                       ! tolerance
  tolCond                    =                 1.0E-23 ! Preconditioner
                                                       ! tolerance
  maxIterPoisson             =                   10000 ! Maximum iterations
  poissonSolverType          =              'bicgstab' ! 'bicgstab', 'gcr',
                                                       ! 'adi' or 'hypre'
  preconditioner             =                   'yes' ! 'no' or 'yes'
  dtau                       =                  4.0E-4 ! Time parameter for
                                                       ! preconditioner
  maxIterADI                 =                       2 ! Preconditioner
                                                       ! iterations
  initialCleaning            =                       T ! Enforce initial non-
                                                       ! divergence
  pressureScaling            =                       F ! Scale by P
  correctMomentum            =                       T ! Correct momentum so
                                                       ! that divergence
                                                       ! constraint is
                                                       ! fulfilled
  correctDivError            =                       F ! Subtract divergence
  tolcrit                    =                   'abs' ! 'abs' or 'rel'
&end

&atmosphereList
  referenceQuantities        =                 'Klein' ! 'Klein', 'WKB', 'SI'
                                                       ! or 'general'
  specifyReynolds            =                       F ! Use inverse Reynolds
                                                       ! number
  ReInv                      =                  0.0E+0 ! Inverse Reynolds number
  mu_viscous_dim             =                  0.0E+0 ! Kinematic viscosity
  mu_conduct_dim             =                  0.0E+0 ! Heat conductivity
  background                 =            'isothermal' ! 'realistic',
                                                       ! 'isothermal',
                                                       ! 'isentropic', 'const-
                                                       ! N', 'diflapse' or
                                                       ! 'HeldSuarez'
  N_BruntVaisala_dim         =                  1.8E-2 ! Buoyancy frequency for
                                                       ! 'const-N'
  theta0_dim                 =                  3.0E+2 ! Background potential
                                                       ! temperature for
                                                       ! 'isentropic'
  Temp0_dim                  =                  3.0E+2 ! Background temperature
                                                       ! for 'isothermal'
  press0_dim                 =                6.888E+4 ! Ground pressure
  backgroundFlow_dim(1)      =                  0.0E+0 ! Initial wind
  backgroundFlow_dim(2)      =                  0.0E+0 ! Initial wind
  backgroundFlow_dim(3)      =                  0.0E+0 ! Initial wind
  f_Coriolis_dim             =                  1.0E-4 ! Coriolis frequency
  corset                     =              'constant' ! 'constant' or
                                                       ! 'periodic'
  z_tr_dim                   =                  1.2E+4 ! Tropopause height
  theta_tr_dim               =                  3.0E+2 ! Potential temperature
                                                       ! in troposphere
  gamma_t                    =                  0.0E+0 ! Lapse rate in
                                                       ! troposphere
  gamma_s                    =                  0.0E+0 ! Lapse rate in
                                                       ! stratosphere
  tp_strato_dim              =                  2.0E+2 ! Temperature in
                                                       ! stratosphere for
                                                       ! 'HeldSuarez'
  tp_srf_trp_dim             =                  3.0E+2 ! Tropical surface
                                                       ! temperature for
                                                       ! 'HeldSuarez'
  tpdiffhor_tropo_dim        =                  5.0E+1 ! Temperature difference
                                                       ! between poles and
                                                       ! tropics for
                                                       ! 'HeldSuarez'
  ptdiffvert_tropo_dim       =                  1.0E+1 ! Vertical potential
                                                       ! temperature difference
                                                       ! in troposphere for
                                                       ! 'HeldSuarez'
&end

&topographyList
  topography                 =                       F ! Terrain-following
                                                       ! coordinates
  ipolTFC                    =                       2 ! Interpolation in the
                                                       ! transformation of w
  freeSlipTFC                =                       F ! Transformed free-slip
                                                       ! condition
  testTFC                    =                       F ! Various TFC tests
  topographyTime             =                  0.0E+0 ! Topography growth time
  mountainHeight_dim         =                  4.0E+2 ! Maximum height
  mountainWidth_dim          =                  1.0E+3 ! Half width
  mountain_case              =                       3 ! Predefined topography
  range_factor               =                  1.0E+1 ! Ratio between large
                                                       ! and small scales
&end

&boundaryList
  rhoFluxCorr                =                       F ! ...
  iceFluxCorr                =                       F ! ...
  uFluxCorr                  =                       F ! ...
  vFluxCorr                  =                       F ! ...
  wFluxCorr                  =                       F ! ...
  thetaFluxCorr              =                       F ! ...
  nbCellCorr                 =                       1 ! ...
  spongeLayer                =                       F ! General sponge layer
                                                       ! switch
  sponge_uv                  =                       F ! Sponge layer for
                                                       ! horizontal wind if
                                                       ! unifiedSponge =
                                                       ! .false.
  spongeHeight               =                  3.3E-1 ! Relative height of
                                                       ! lower sponge layer
                                                       ! edge (scale height for
                                                       ! unifiedSponge = .true.
                                                       ! and spongeType =
                                                       ! 'exponential')
  spongeAlphaZ_dim           =                  2.0E-4 ! Maximum relaxation
                                                       ! rate for unifiedSponge
                                                       ! = .true.
  spongeAlphaZ_fac           =                  1.0E-2 ! Sponge layer factor
                                                       ! for unifiedSponge =
                                                       ! .false.
  unifiedSponge              =                       F ! Unified sponge for
                                                       ! both time schemes,
                                                       ! applied to wind and
                                                       ! density
  lateralSponge              =                       F ! Lateral sponge for
                                                       ! unifiedSponge = .true.
  spongeType                 =            'polynomial' ! Sponge layer profile
                                                       ! for unifiedSponge =
                                                       ! .true.
  spongeOrder                =                       1 ! Order of polynomial
                                                       ! sponge
  cosmoSteps                 =                       1 ! Relative strength of
                                                       ! COSMO sponge
&end

&boundaryList2
  xBoundary                  =              'periodic' ! Boundary conditions in
                                                       ! x ('periodic' only)
  yBoundary                  =              'periodic' ! Boundary conditions in
                                                       ! y ('periodic' only)
  zBoundary                  =            'solid_wall' ! Boundary conditions in
                                                       ! z ('periodic' or
                                                       ! 'solid_wall')
&end

&wkbList
  rayTracer                  =                       T ! Ray-tracer switch
&end

&tracerList
  tracerSetup                =  'increase_in_z_tracer' ! ...
  include_gw_tracer_forcing  =                       T ! ...
  include_tracer_mixing      =                       T ! ...
  tracerdifference           =                       T ! ...
  include_prime              =                       T ! ...
  include_env_tracer_forcing =                       T ! ...
&end

&iceList
  iceTestcase                =             'ice_cloud' ! ...
  init_SIce                  =                  1.0E+5 ! ...
  init_nAer                  =                  1.0E+8 ! ...
  init_qv                    =                  1.0E-4 ! ...
  init_m_ice                 =                 1.0E-16 ! ...
  radius_solution            =                  1.0E-8 ! ...
  sigma_r                    =                  1.0E+0 ! ...
  T_nuc                      =                  2.0E+2 ! ...
  p_nuc                      =                  2.0E+4 ! ...
  dt_ice                     =                  1.0E-3 ! ...
  NUC_approx_type            =                'linFit' ! ...
  ISSR_top                   =                  0.0E+0 ! ...
  super_simplified           =                       F ! ...
  kT_linFit                  =                       F ! ...
  dv_exp2                    =                       F ! ...
  cm_dryAir                  =                       F ! ...
  mu_linFit                  =                       F ! ...
  sedimentation_on           =                       F ! ...
  nucleation_on              =                       F ! ...
  evaporation_on             =                       F ! ...
  awi_type                   =                 'const' ! ...
  SIce_threshold_type        =                'linFit' ! ...
&end

&iceList2
  inN                        =                       0 ! ...
  inQ                        =                       0 ! ...
  inQv                       =                       0 ! ...
  nVarIce                    =                       0 ! ...
  dt_ice2                    =                  0.0E+0 ! ...
  no_ice_source              =                       F ! ...
&end
