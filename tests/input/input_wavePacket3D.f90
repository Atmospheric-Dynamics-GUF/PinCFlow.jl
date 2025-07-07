!------------------------------------------------------------------------------!
!                              PincFlow Namelists                              !
!------------------------------------------------------------------------------!

&domain
  sizeX                    =                     512 ! Cells in x
  sizeY                    =                      16 ! Cells in y
  sizeZ                    =                    1000 ! Cells in z
  nbx                      =                       3 ! Halo/ghost cells in x
  nby                      =                       3 ! Halo/ghost cells in y
  nbz                      =                       3 ! Ghost cells in z
  lx_dim(0)                =                  0.0E+0 ! Minimum of x
  lx_dim(1)                =                  9.0E+6 ! Maximum of x
  ly_dim(0)                =                  0.0E+0 ! Minimum of y
  ly_dim(1)                =                  3.0E+5 ! Maximum of y
  lz_dim(0)                =                  0.0E+0 ! Minimum of z
  lz_dim(1)                =                  1.0E+5 ! Maximum of z
  nprocx                   =                {nprocx} ! Processors in x
  nprocy                   =                {nprocy} ! Processors in y
&end

&variables
  include_ice              =                       F ! Ice microphysics
                                                     ! parameterization
  include_tracer           =                       F ! Tracer equation
  include_testoutput       =                       F ! ...
&end

&outputList
  atmvarOut(1)             =                  'rhop' ! Atmospheric output
                                                     ! variable
  rayvarOut(1)             =                      '' ! Raytracer output variable
  icevarOut(1)             =                      '' ! Ice output variable
  saverayvols              =                       F ! Save ray volumes
  prepare_restart          =                       F ! Save everything needed
                                                     ! for a restart
  restart                  =                       F ! Restart the model from
                                                     ! state in previous
                                                     ! simulation
  iIn                      =                      -1 ! Restart at time step iIn
  runName                  =          'wavePacket3D' ! Run name for netCDF file
  outputType               =                  'time' ! 'timeStep' or 'time'
  nOutput                  =                       1 ! Output every nOutput
                                                     ! time steps for
                                                     ! outputType = 'timeStep'
  maxIter                  =                       1 ! Stop after maxIter time
                                                     ! steps for outputType =
                                                     ! 'timeStep'
  outputTimeDiff           =                 1.08E+4 ! Output every
                                                     ! outputTimeDiff seconds
                                                     ! for outputType = 'time'
  maxTime                  =                 1.08E+4 ! Stop after maxTime
                                                     ! seconds for outputType =
                                                     ! 'time'
  detailedinfo             =                       F ! Provide info on the
                                                     ! final state of Poisson
                                                     ! solver
  RHS_diagnostics          =                       T ! Provide info about the
                                                     ! RHS of Poisson equation
  fancy_namelists          =                       T ! Write all namelists with
                                                     ! comments
&end

&debuggingList
  verbose                  =                       F ! ...
  dtMin_dim                =                  1.0E-5 ! Stop if dt < dtMin
&end

&testCaseList
  testCase                 =            'wavePacket' ! Predefined test case
&end

&monochromeWave
  lambda_dim               =                  6.0E+3 ! Vertical wavelength
&end

&wavePacket
  wavePacketType           =                       1 ! 1 = Gaussian, 2 = Cosine
  wavePacketDim            =                       2 ! 1 = 1D, 2 = 2D, 3 = 3D
  lambdaX_dim              =                  3.0E+5 ! Wavelength in x (0.0 for
                                                     ! infinite wavelength)
  lambdaY_dim              =                  3.0E+5 ! Wavelength in y (0.0 for
                                                     ! infinite wavelength)
  lambdaZ_dim              =                 -1.0E+3 ! Wavelength in z
  amplitudeFactor          =                  5.0E-1 ! Normalized buoyancy
                                                     ! amplitude
  x0_dim                   =                  4.5E+6 ! Wave-packet center in x
  y0_dim                   =                  1.5E+5 ! Wave-packet center in y
  z0_dim                   =                  3.0E+4 ! Wave-packet center in z
  sigma_dim                =                  5.0E+3 ! Vertical width of
                                                     ! Gaussian wave packet
  sigma_hor_dim            =                  1.5E+6 ! Cosine distribution
                                                     ! width in x (0.0 for
                                                     ! infinite width)
  amp_mod_x                =                  1.0E+0 ! Fractional amplitude
                                                     ! modulation in x
  sigma_hor_yyy_dim        =                  1.5E+6 ! Cosine distribution
                                                     ! width in y (0.0 for
                                                     ! infinite width)
  amp_mod_y                =                  1.0E+0 ! Fractional amplitude
                                                     ! modulation in y
  L_cos_dim                =                  1.0E+4 ! Half width of vertical
                                                     ! cosine profile of wave
                                                     ! packet
  omiSign                  =                       1 ! Frequency branch
  u0_jet_dim               =                  0.0E+0 ! Maximum amplitude of jet
  z0_jet_dim               =                  5.0E+4 ! Center of jet
  L_jet_dim                =                  5.0E+3 ! Half width of vertical
                                                     ! cosine profile of jet
  inducedwind              =                       F ! ...
&end

&LagrangeRayTracing
  xrmin_dim                =                  0.0E+0 ! Left bound of initial
                                                     ! rays
  xrmax_dim                =                  4.0E+6 ! Right bound of initial
                                                     ! rays
  yrmin_dim                =                  0.0E+0 ! Backward bound of
                                                     ! initial rays
  yrmax_dim                =                  4.0E+4 ! Forward bound of initial
                                                     ! rays
  zrmin_dim                =                  3.0E+3 ! Bottom bound of initial
                                                     ! rays
  zrmax_dim                =                  7.0E+4 ! Top bound of initial rays
  nrxl                     =                       1 ! Initial ray volumes
                                                     ! within dx
  nryl                     =                       1 ! Initial ray volumes
                                                     ! within dy
  nrzl                     =                       1 ! Initial ray volumes
                                                     ! within dz
  fac_dk_init              =                  1.0E-1 ! Initial fraction dk/kh
  fac_dl_init              =                  1.0E-1 ! Initial fraction dl/kh
  fac_dm_init              =                  1.0E-1 ! Initial fraction dm/m
  nrk_init                 =                       2 ! Initial ray volumes
                                                     ! within dk
  nrl_init                 =                       1 ! Initial ray volumes
                                                     ! within dl
  nrm_init                 =                       2 ! Initial ray volumes
                                                     ! within dm
  nsmth_wkb                =                       2 ! Width of smoothing
                                                     ! operator for mean flow
                                                     ! tendencies
  lsmth_wkb                =                       T ! Smoothing operator for
                                                     ! mean flow tendencies
  sm_filter                =                       2 ! 1 = Box filter, 2 =
                                                     ! Shapiro filter
  lsaturation              =                       T ! Switch for saturation
                                                     ! scheme
  alpha_sat                =                  1.0E+0 ! Saturation threshold
  single_column            =                       F ! Single-column mode
  steady_state             =                       F ! Steady-state mode
  case_wkb                 =                       3 ! 1 = Gaussian wave
                                                     ! packet, 2 = cosine wave
                                                     ! packet, 3 = mountain
                                                     ! wave
  amp_wkb                  =                  5.0E-1 ! Relative amplitude of
                                                     ! wave packet
  wlrx_init                =                  1.0E+5 ! Wavelength in x of wave
                                                     ! packet
  wlry_init                =                  0.0E+0 ! Wavelength in y of wave
                                                     ! packet
  wlrz_init                =                  1.0E+3 ! Wavelength in z of wave
                                                     ! packet
  xr0_dim                  =                  2.0E+6 ! Wave packet center in x
  yr0_dim                  =                  2.0E+4 ! Wave packet center in y
  zr0_dim                  =                  3.0E+4 ! Wave packet center in z
  sigwpx_dim               =                  1.0E+6 ! Wave packet width in x
                                                     ! (0.0 for infinite width)
  sigwpy_dim               =                  0.0E+0 ! Wave packet width in y
                                                     ! (0.0 for infinite width)
  sigwpz_dim               =                  5.0E+3 ! Wave packet width in z
                                                     ! (0.0 for infinite width)
  branchr                  =                      -1 ! Frequency branch
  lindUinit                =                       F ! Induced wind at initial
                                                     ! time
  blocking                 =                       F ! Blocked-layer scheme
  long_threshold           =                  2.5E-1 ! Long-number threshold of
                                                     ! the blocked-layer scheme
  drag_coefficient         =                  1.0E+0 ! Drag coefficient of the
                                                     ! blocked-layer scheme
  nwm                      =                       1 ! Number of initial wave
                                                     ! modes
  zmin_wkb_dim             =                  0.0E+0 ! Minimum altitude for
                                                     ! wave-mean-flow
                                                     ! interaction
  nray_fac                 =                      20 ! Maximum multiplication
                                                     ! factor per spectral
                                                     ! dimension
  cons_merge               =                    'en' ! Conserved quantity in
                                                     ! ray-volume merging ('wa'
                                                     ! = wave action, 'en' =
                                                     ! wave energy)
&end

&bubble
  dTheta0_dim              =                  6.0E+0 ! ...
  xRadius_dim              =                  2.5E+3 ! ...
  zRadius_dim              =                  2.5E+3 ! ...
  xCenter_dim              =                  0.0E+0 ! ...
  yCenter_dim              =                  1.5E+5 ! ...
  zCenter_dim              =                  0.0E+0 ! ...
  zExcentricity            =                  1.0E+0 ! ...
  rhoCenter_dim            =                  1.0E+3 ! ...
&end

&robert_bubble
  dTheta1_dim              =                  5.0E-1 ! Potential temperature
                                                     ! offset
  a1_dim                   =                  1.5E+2 ! Radius of plateau
  sigma1_dim               =                  5.0E+1 ! Gaussian edge profile
  xCenter1_dim             =                  0.0E+0 ! ...
  zCenter1_dim             =                  2.5E+4 ! ...
  dTheta2_dim              =                 -1.5E-1 ! ...
  a2_dim                   =                  0.0E+0 ! ...
  sigma2_dim               =                  5.0E+1 ! ...
  xCenter2_dim             =                  6.0E+1 ! ...
  zCenter2_dim             =                  6.4E+2 ! ...
&end

&mountainwavelist
  u_relax                  =                  7.5E+1 ! Zonal relaxation wind
  v_relax                  =                  0.0E+0 ! Meridional relaxation
                                                     ! wind
  w_relax                  =                  0.0E+0 ! Vertical relaxation wind
  t_relax                  =                1.728E+5 ! Relaxation time
  t_ramp                   =                  3.6E+3 ! Not used at the moment
  xextent_relax            =                  3.0E+6 ! Zonal extent of
                                                     ! relaxation region
  yextent_relax            =                  0.0E+0 ! Meridional extent of
                                                     ! relaxation region
  wind_relaxation          =                       F ! Relaxation switch
  surface_layer_depth      =                  0.0E+0 ! Surface-layer depth
&end

&baroclinic_LC
  zero_initial_state       =                       F ! ...
  z_trpp0_dim              =                  8.0E+3 ! Mean tropopause height
  z_baro_dim               =                  2.4E+4 ! Altitude above which the
                                                     ! atmosphere is barotropic
  thet0_dim                =                 2.73E+2 ! Characteristic potential
                                                     ! temperature
  ntrp_dim                 =                  1.0E-2 ! Buoyancy frequency of
                                                     ! the troposphere
  nstr_dim                 =                 2.45E-2 ! Buoyancy frequency of
                                                     ! the stratosphere
  jwdth_dim                =                 8.25E+5 ! Jet half width
  kaptpp                   =                  7.0E+1 ! Jet slope
  add_ptptb                =                       T ! Add local potential-
                                                     ! temperature perturbation
  ptptb_x_dim              =                  5.0E+6 ! Location in x of local
                                                     ! potential-temperature
                                                     ! perturbation
  ptptb_y_dim              =                  4.0E+6 ! Location in y of local
                                                     ! potential-temperature
                                                     ! perturbation
  ptptb_z_dim              =                  9.0E+3 ! Location in z of local
                                                     ! potential-temperature
                                                     ! perturbation
  ptptb_dh_dim             =                  5.0E+5 ! Horizontal width of
                                                     ! local potential-
                                                     ! temperature perturbation
  ptptb_dz_dim             =                  2.0E+3 ! Vertical width of local
                                                     ! potential-temperature
                                                     ! perturbation
  ptptb_amp_dim            =                  3.0E+0 ! Amplitude of local
                                                     ! potential-temperature
                                                     ! perturbation
  add_noise                =                       F ! Add noise to the initial
                                                     ! potential temperature
  proc_noise               =                  1.0E-2 ! Relative amplitude of
                                                     ! potential-temperature
                                                     ! noise
  tau_relax                =                1.296E+6 ! ...
  tau_relax_low            =                  0.0E+0 ! ...
  sigma_tau                =                  0.0E+0 ! Relative thickness of
                                                     ! relaxation profile
  tau_jet                  =                  7.2E+3 ! Time for jet formation
  Sponge_Rel_Bal_Type      =                   'env' ! Relaxation state: 'hyd'
                                                     ! = hydrostatic balance,
                                                     ! 'env' = geostrophic
                                                     ! balance
  ta_hs_dim                =                  0.0E+0 ! Thermal-relaxation time
                                                     ! scale outside tropical
                                                     ! boundary layer (zero
                                                     ! means infinity)
  ts_hs_dim                =                  0.0E+0 ! Thermal-relaxation time
                                                     ! scale in tropical
                                                     ! boundary layer (zero
                                                     ! means infinity)
  tf_hs_dim                =                  0.0E+0 ! Boundary-layer Rayleigh-
                                                     ! damping  time scale
  sigb_hs                  =                  0.0E+0 ! Sigma of boundary-layer
                                                     ! top
&end

&baroclinic_ID
  bar_sigma_y              =                  1.0E-1 ! ...
  u_strength               =                  3.0E+1 ! Jet strength
  dTh_atm                  =                  2.0E+1 ! Potential-temperature
                                                     ! difference between poles
                                                     ! and tropics
  init_2Dto3D              =                       F ! Initialize 3D state from
                                                     ! 2D balanced state
  init_bal                 =             'geostr_id' ! ...
  lastrecordnum            =                      88 ! Last record in 2D output
  fileinitstate2D          =    'pf_all_2D3D_in.dat' ! File name with 2D
                                                     ! balanced state
  output_theta_bgr         =                       F ! Output environmental
                                                     ! potential temperature
  output_br_vais_sq        =                       T ! Output environmental
                                                     ! squared buoyancy
                                                     ! frequency
  output_heat              =                       T ! ...
  balance_eq               =                 'PI_an' ! ...
  output_rho_bgr           =                       T ! Output environmental
                                                     ! density
&end

&modelList
  model                    = 'pseudo_incompressible' ! Dynamic equations
  vert_theta               =                  9.0E+1 ! Rotation about x
  vert_alpha               =                  0.0E+0 ! Rotation about y
&end

&solverList
  cfl                      =                  5.0E-1 ! CFL number
  cfl_wave                 =                  5.0E-1 ! WKB CFL number
  dtMax_dim                =                  1.2E+3 ! Maximum time step
  tStepChoice              =                   'cfl' ! 'fix' or 'cfl'
  timeScheme               =          'semiimplicit' ! 'LS_Will_RK3' or
                                                     ! 'semiimplicit'
  auxil_equ                =                       F ! Buoyancy equation
  limiterType1             =             'MCVariant' ! 'minmod', 'MCVariant' or
                                                     ! 'Cada'
  TurbScheme               =                       F ! Turbulence scheme
  turb_dts                 =                  5.0E+3 ! Turbulent damping time
  DySmaScheme              =                       T ! Dynamic Smagorinsky
                                                     ! scheme
  dtWave_on                =                       T ! Limit time step by
                                                     ! inverse buoyancy
                                                     ! frequency
  heatingONK14             =                       F ! Heating as in O'Neill
                                                     ! and Klein (2014)
  dens_relax               =                       F ! Heating by density
                                                     ! relaxation
  shap_dts_fac             =                 -1.0E+0 ! Shapiro-filter damping
                                                     ! time
  n_shap                   =                       1 ! Order of Shapiro filter
&end

&poissonSolverList
  tolPoisson               =                  1.0E-8 ! Abort criterion
  abs_tol                  =                  0.0E+0 ! Lower bound for tolerance
  tolCond                  =                 1.0E-23 ! Preconditioner tolerance
  maxIterPoisson           =                    1000 ! Maximum iterations
  preconditioner           =                   'yes' ! 'no' or 'yes'
  dtau                     =                  4.0E+0 ! Time parameter for
                                                     ! preconditioner
  maxIterADI               =                       2 ! Preconditioner iterations
  initialCleaning          =                       T ! Enforce initial non-
                                                     ! divergence
  correctMomentum          =                       T ! Correct momentum so that
                                                     ! divergence constraint is
                                                     ! fulfilled
  correctDivError          =                       F ! Subtract divergence
  tolcrit                  =                   'abs' ! 'abs' or 'rel'
&end

&atmosphereList
  referenceQuantities      =                 'Klein' ! 'Klein', 'WKB', 'SI' or
                                                     ! 'general'
  specifyReynolds          =                       F ! Use inverse Reynolds
                                                     ! number
  ReInv                    =                  0.0E+0 ! Inverse Reynolds number
  mu_viscous_dim           =                  0.0E+0 ! Kinematic viscosity
  mu_conduct_dim           =                  0.0E+0 ! Heat conductivity
  background               =            'isothermal' ! 'realistic',
                                                     ! 'isothermal',
                                                     ! 'isentropic', 'const-N',
                                                     ! 'diflapse' or
                                                     ! 'HeldSuarez'
  N_BruntVaisala_dim       =                  1.8E-2 ! Buoyancy frequency for
                                                     ! 'const-N'
  theta0_dim               =                  3.0E+2 ! Background potential
                                                     ! temperature for
                                                     ! 'isentropic'
  Temp0_dim                =                  3.0E+2 ! Background temperature
                                                     ! for 'isothermal'
  press0_dim               =              1.01325E+5 ! Ground pressure
  backgroundFlow_dim(1)    =                  0.0E+0 ! Initial wind
  backgroundFlow_dim(2)    =                  0.0E+0 ! Initial wind
  backgroundFlow_dim(3)    =                  0.0E+0 ! Initial wind
  f_Coriolis_dim           =                  0.0E+0 ! Coriolis frequency
  corset                   =              'constant' ! 'constant' or 'periodic'
  z_tr_dim                 =                  1.2E+4 ! Tropopause height
  theta_tr_dim             =                  2.4E+2 ! Potential temperature in
                                                     ! troposphere
  gamma_t                  =                  0.0E+0 ! Lapse rate in troposphere
  gamma_s                  =                  0.0E+0 ! Lapse rate in
                                                     ! stratosphere
  tp_strato_dim            =                  2.0E+2 ! Temperature in
                                                     ! stratosphere for
                                                     ! 'HeldSuarez'
  tp_srf_trp_dim           =                  3.0E+2 ! Tropical surface
                                                     ! temperature for
                                                     ! 'HeldSuarez'
  tpdiffhor_tropo_dim      =                  5.0E+1 ! Temperature difference
                                                     ! between poles and
                                                     ! tropics for 'HeldSuarez'
  ptdiffvert_tropo_dim     =                  1.0E+1 ! Vertical potential
                                                     ! temperature difference
                                                     ! in troposphere for
                                                     ! 'HeldSuarez'
&end

&topographyList
  topography               =                       F ! Terrain-following
                                                     ! coordinates
  ipolTFC                  =                       2 ! Interpolation in the
                                                     ! transformation of w
  topographyTime           =                  0.0E+0 ! Topography growth time
  mountainHeight_dim       =                  5.0E+2 ! Maximum height
  mountainWidth_dim        =                  1.0E+6 ! Half width
  mountain_case            =                       2 ! Predefined topography
  height_factor            =                  1.0E+0 ! Ratio between large- and
                                                     ! small-scale wave
                                                     ! amplitudes
  width_factor             =                  1.0E+1 ! Ratio between large- and
                                                     ! small-scale wavelengths
  spectral_modes           =                       1 ! Number of spectral modes
  stretch_exponent         =                  1.0E+0 ! Exponent of vertical
                                                     ! grid stretching (1 for
                                                     ! no stretching)
&end

&boundaryList
  rhoFluxCorr              =                       F ! ...
  iceFluxCorr              =                       F ! ...
  uFluxCorr                =                       F ! ...
  vFluxCorr                =                       F ! ...
  wFluxCorr                =                       F ! ...
  thetaFluxCorr            =                       F ! ...
  nbCellCorr               =                       1 ! ...
  spongeLayer              =                       F ! General sponge layer
                                                     ! switch
  sponge_uv                =                       F ! Sponge layer for
                                                     ! horizontal wind if
                                                     ! unifiedSponge = .false.
  spongeHeight             =                  3.3E-1 ! Relative height of lower
                                                     ! sponge layer edge (scale
                                                     ! height for unifiedSponge
                                                     ! = .true. and spongeType
                                                     ! = 'exponential')
  spongeAlphaZ_dim         =                  2.0E-4 ! Maximum relaxation rate
                                                     ! for unifiedSponge =
                                                     ! .true.
  spongeAlphaZ_fac         =                  1.0E-2 ! Sponge layer factor for
                                                     ! unifiedSponge = .false.
  unifiedSponge            =                       F ! Unified sponge for both
                                                     ! time schemes, applied to
                                                     ! wind and density
  lateralSponge            =                       F ! Lateral sponge for
                                                     ! unifiedSponge = .true.
  spongeType               =            'polynomial' ! Sponge layer profile for
                                                     ! unifiedSponge = .true.
  spongeOrder              =                       1 ! Order of polynomial
                                                     ! sponge
  cosmoSteps               =                       1 ! Relative strength of
                                                     ! COSMO sponge
  relax_to_mean            =                       T ! Relax the wind to its
                                                     ! (terrain-following)
                                                     ! horizontal mean
                                                     ! (otherwise, relax to the
                                                     ! initial state) if
                                                     ! unifiedSponge == .true.
  relaxation_period        =                  0.0E+0 ! Period of an oscillation
                                                     ! superposed on the
                                                     ! background wind if
                                                     ! unifiedSponge == .true.
                                                     ! and relax_to_mean ==
                                                     ! .false. (0 for no
                                                     ! oscillation)
  relaxation_amplitude     =                  0.0E+0 ! Relative amplitude of an
                                                     ! oscillation superposed
                                                     ! on the background wind
                                                     ! if unifiedSponge ==
                                                     ! .true. and relax_to_mean
                                                     ! == .false.
  xBoundary                =              'periodic' ! Boundary conditions in x
                                                     ! ('periodic' only)
  yBoundary                =              'periodic' ! Boundary conditions in y
                                                     ! ('periodic' only)
  zBoundary                =            'solid_wall' ! Boundary conditions in z
                                                     ! ('periodic' or
                                                     ! 'solid_wall')
&end

&wkbList
  rayTracer                =                       F ! Ray-tracer switch
&end

&tracerList
  tracerSetup              =               'alpha_z' ! initial tracer
                                                     ! distribution
  include_trfrc_lo         =                       T ! leading-order GW tracer
                                                     ! forcing
  include_trfrc_no         =                       T ! next-order GW tracer
                                                     ! forcing
  include_trfrc_mix        =                       T ! diffusive tracer mixing
&end

&iceList
  inN                      =                       0 ! ...
  inQ                      =                       0 ! ...
  inQv                     =                       0 ! ...
  nVarIce                  =                       0 ! ...
  dt_ice                   =                  0.0E+0 ! ...
  no_ice_source            =                       F ! ...
&end
