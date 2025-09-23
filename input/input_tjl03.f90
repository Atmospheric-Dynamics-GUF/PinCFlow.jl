!------------------------------------------------------------------------------!
!                              PincFlow Namelists                              !
!------------------------------------------------------------------------------!

&domain
  sizeX                    =                      10 ! Cells in x
  sizeY                    =                       1 ! Cells in y
  sizeZ                    =                      10 ! Cells in z
  nbx                      =                       3 ! Halo/ghost cells in x
  nby                      =                       3 ! Halo/ghost cells in y
  nbz                      =                       3 ! Ghost cells in z
  lx_dim(0)                =                  0.0E+0 ! Minimum of x
  lx_dim(1)                =                  4.0E+4 ! Maximum of x
  ly_dim(0)                =                  0.0E+0 ! Minimum of y
  ly_dim(1)                =                  1.0E+4 ! Maximum of y
  lz_dim(0)                =                  0.0E+0 ! Minimum of z
  lz_dim(1)                =                  1.5E+4 ! Maximum of z
  nprocx                   =                1 ! Processors in x
  nprocy                   =                1 ! Processors in y
&end

&variables
  include_ice              =                       T ! Ice microphysics
                                                     ! parameterization
  include_tracer           =                       F ! Tracer equation
  include_testoutput       =                       T ! ...
&end

&outputList
  atmvarOut(1)             =                  'rhop' ! Atmospheric output
                                                     ! variable
  atmvarOut(2)             =                     'u' ! Atmospheric output
                                                     ! variable
  atmvarOut(3)             =                     'v' ! Atmospheric output
                                                     ! variable
  atmvarOut(4)             =                     'w' ! Atmospheric output
                                                     ! variable
  atmvarOut(5)             =                'thetap' ! Atmospheric output
                                                     ! variable
  atmvarOut(6)             =                   'pip' ! Atmospheric output
                                                    ! variable
  icevarOut(1)             =                     'n' ! Ice output variable
  icevarOut(2)             =                     'q' ! Ice output variable
  icevarOut(3)             =                    'qv' ! Ice output variable
  optvarOut(1)             =                     's' !
  optvarOut(2)             =                     'w' !
  optvarOut(3)             =                     't' !
  rayvarOut(1)             =                      '' ! Raytracer output variable
  icevarOut(1)             =                     'n' ! Ice output variable
  icevarOut(2)             =                     'q' ! Ice output variable
  icevarOut(3)             =                    'qv' ! Ice output variable
  saverayvols              =                       F ! Save ray volumes
  prepare_restart          =                       F ! Save everything needed
                                                     ! for a restart
  restart                  =                       F ! Restart the model from
                                                     ! state in previous
                                                     ! simulation
  iIn                      =                      -1 ! Restart at time step iIn
  runName                  =         'wavePacketIce' ! Run name for netCDF file
  outputType               =                  'time' ! 'timeStep' or 'time'
  nOutput                  =                       1 ! Output every nOutput
                                                     ! time steps for
                                                     ! outputType = 'timeStep'
  maxIter                  =                       1 ! Stop after maxIter time
                                                     ! steps for outputType =
                                                     ! 'timeStep'
  outputTimeDiff           =                  1.0E+0 ! Output every
                                                     ! outputTimeDiff seconds
                                                     ! for outputType = 'time'
  maxTime                  =                  12.0E+0 ! Stop after maxTime
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
  lambda_dim               =                  1.5E+3 ! Vertical wavelength
&end
&wavePacket
  wavePacketType           =                       2 ! 1 = Gaussian, 2 = Cosine
  wavePacketDim            =                       2 ! 1 = 1D, 2 = 2D, 3 = 3D
  lambdaX_dim              =                 10.0E+3 ! Wavelength in x (0.0 for
                                                     ! infinite wavelength)
  lambdaY_dim              =                       0 ! Wavelength in y (0.0 for
                                                     ! infinite wavelength)
  lambdaZ_dim              =                 -2.0E+3 ! Wavelength in z
  amplitudeFactor          =                  0.19   ! Normalized buoyancy
                                                     ! amplitude
  x0_dim                   =                  5.0E+3 ! Wave-packet center in x
  y0_dim                   =                  5.0E+3 ! Wave-packet center in y
  z0_dim                   =                  5.0E+3 ! Wave-packet center in z
  sigma_dim                =                  6.0E+3 ! Vertical width of
                                                     ! Gaussian wave packet
  sigma_hor_dim            =                  0.0E+0 ! Cosine distribution
                                                     ! width in x (0.0 for
                                                     ! infinite width)
  amp_mod_x                =                  1.0E+0 ! Fractional amplitude
                                                     ! modulation in x
  sigma_hor_yyy_dim        =                  0.0E+0 ! Cosine distribution
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
  MultipleWavePackets      =                       T ! superposition of WPs
&end

&MultipleWavePackets_List

  amplitudeFactor_sp       =             0.12, 0.12, ! Normalized buoyancy
                                                     ! amplitude
  omiSign_sp               =                   1, 1, !Frequency branch         
  lambdaX_dim_sp           =           10.e3, 10.e3, ! Wavelength in x (0.0 for
                                                     ! infinite wavelength)
  lambdaY_dim_sp           =                 0., 0., ! Wavelength in y (0.0 for
                                                     ! infinite wavelength) 
  lambdaZ_dim_sp           =            -2.e3, 2.e3, ! Wavelength in z 
  sigma_hor_dim_sp         =                 0., 0., !
  sigma_hor_yyy_dim_sp     =                 0., 0., !
  sigma_dim_sp             =             6.e3, 6.e3, !
  x0_dim_sp                =             5.e3, 5.e3, !
  y0_dim_sp                =             5.e3, 5.e3, !
  z0_dim_sp                =              7.e3, 9.e3 !

&end

&LagrangeRayTracing

  xrmin_dim = 0.0,         ! left bound of initial rays (in x direction) (m)
  xrmax_dim = 40.e3,        ! right bound of initial rays (in x dir.) (m)
  yrmin_dim = 0.0,         ! left bound of initial rays (in y direction) (m)
  yrmax_dim = 10.e3,        ! right bound of initial rays (in y dir.) (m)
  zrmin_dim = 0.e3,        ! bottom bound of initial rays (m)
  zrmax_dim = 15.e3,        ! top bound of initial rays (m)

  nrxl = 1,               ! no. of ray vol. init. within one hor. x column
  nryl = 1,               ! no. of ray vol. init. within one hor. y column
  nrzl = 1,               ! no. of ray vol. init. within one vert. layer

  fac_dk_init = 0.1,     ! init. width of total ray vol. in k space
                           ! (fraction of the initial wave number in x dir.)
  fac_dl_init = 0.1,     ! init. width of total ray vol. in l space
                           ! (fraction of the initial wave number in y dir.)
  fac_dm_init = 0.0016,     ! init. width of total ray vol. in m space
                           ! (fraction of the initial vert. wave number)

  nrk_init = 1,            ! no. of ray volumes initialized within dk
  nrl_init = 1,            ! no. of ray volumes initialized within dl
  nrm_init = 1,            ! no. of ray volumes initialized within dm

  nsmth_wkb = 2,           ! half (number -1) of cells f. smooth. wkb fluxes
  lsmth_wkb = .true.,      ! log. switch for smooth. wkb data (true/false)
  sm_filter = 2,
  
  lsaturation = .true.,    ! JaWi 16.12.16 (sat)
  alpha_sat = 1.0,         ! JaWi 16.12.16 (sat)

  case_wkb = 5,            ! 1/2: Gaussian/Cosine wave packet; 3: mountain
  nwm = 2
  amp_wkb = 0.19            ! amplitude of the wave packet (wrt saturation)

  wlrx_init = 10.e3,        ! initial lambda_x of the wave packet (m)
  wlry_init = 0.0          ! initial lambda_y of the wave packet (m)
                           ! (0 means infinity)
  wlrz_init = -2.e3,        ! initial lambda_z of the wave packet (m)
                           ! (0 means infinity)

  xr0_dim = 5.e3           ! center of the wave packet in hor. (x-dir.) (m)
  yr0_dim = 5.e3,          ! center of the wave packet in hor. (y-dir.) (m)
  zr0_dim = 5.e3,          ! center of the wave packet in vertical (m)

  sigwpx_dim = 0. ,      ! width of the wave packet in hor. (x-dir.) (m);
                           ! (0 means infinity)
  sigwpy_dim = 0.e0        ! width of the wave packet in hor. (y-dir.) (m);
                           ! (0 means infinity)
  sigwpz_dim = 6.e3,       ! width of the wave packet in vertical (m);

  branchr = 1,            ! frequency branch (dispersion relation)
  !presently not used:
  lindUinit = .false.,     ! ind. wind already at initial time (true/false)

  zmin_wkb_dim = 0.e0     ! minumum altitude (above the model bottom, in m)
                           ! for WKB wave-mean-flow interaction
                           ! (zmin_wkb > 0 can help preventing the
                           ! ray volumes being trapped by the self-induced
                           ! mean wind)

  nray_fac = 4            ! maximum factor (per wavenumber direction) by
                           ! which # of rays may increase in comparison to
                           ! initialization

  cons_merge = "en"        ! quantity to be conserved
                           ! ("wa" = wave action/ "en" = wave energy)
                           ! under ray-volume merging
!  nRayOutput =          10 ! Number of dominant ray
!                                                     ! volumes in output

&end

! raytracer + case_wkb=5

&case_wkb_5_list

  amp_wkb_sp = 0.12, 0.12,
!  amp_wkb_sp = 0.1, 0. ,!0.0144,
  branchr_sp = 1, 1,

  wlrx_init_sp = 10.e3, 10.e3,
  wlry_init_sp = 0., 0., 
  wlrz_init_sp = -2.e3, 2.e3, 

  sigwpx_dim_sp = 0., 0., 
  sigwpy_dim_sp = 0., 0., 
  sigwpz_dim_sp = 6.e3, 6.e3, 

  xr0_dim_sp = 5.e3, 5.e3,
  yr0_dim_sp = 5.e3, 5.e3,
  zr0_dim_sp = 7.e3, 9.e3

&end

&modelList
  model                    = 'pseudo_incompressible' ! Dynamic equations
  vert_theta               =                  9.0E+1 ! Rotation about x
  vert_alpha               =                  0.0E+0 ! Rotation about y
&end

&solverList
  cfl                      =                  5.0E-1 ! CFL number
  cfl_wave                 =                  5.0E-1 ! WKB CFL number
  dtMax_dim                =                  1.0E+0 ! Maximum time step
  tStepChoice              =                   'cfl' ! 'fix' or 'cfl'
  timeScheme               =           'LS_Will_RK3' ! 'LS_Will_RK3' or
                                                     ! 'semiimplicit'
  auxil_equ                =                       F ! Buoyancy equation
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
  tolPoisson               =                  1.0E-4 ! Abort criterion
  abs_tol                  =                  0.0E+0 ! Lower bound for tolerance
  tolCond                  =                 1.0E-23 ! Preconditioner tolerance
  maxIterPoisson           =                    1000 ! Maximum iterations
  preconditioner           =                   'yes' ! 'no' or 'yes'
  dtau                     =                  8.0E-1 ! Time parameter for
                                                     ! preconditioner
  maxIterADI               =                      10 ! Preconditioner iterations
  initialCleaning          =                       F ! Enforce initial non-
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
  f_Coriolis_dim           =                  1.0E-4 ! Coriolis frequency
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
  range_factor             =                  1.0E+1 ! Ratio between large and
                                                     ! small scales
  spectral_modes           =                       1 ! Number of spectral modes
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
&end

&boundaryList2
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
  inN                      =                       1 ! ...
  inQ                      =                       2 ! ...
  inQv                     =                       3 ! ...
  nVarIce                  =                       3 ! ...
  dt_ice                   =                  1.0E+0 ! ...
&end
