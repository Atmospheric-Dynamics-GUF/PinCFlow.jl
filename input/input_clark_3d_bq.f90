!------------------------------------------------
!   This is input file for pincFloit
!   The ending .f90 is chosen to facilitate
!   the f90 mode with Emacs
!------------------------------------------------


!------------------------------------------------
!                  Grid and Domain
!------------------------------------------------

&domain

  ! nb of global grid cells
  sizeX = 62,
  sizeY = 62,
  sizeZ = 82,
  ! nb. of ghost cells
  nbx = 2,
  nby = 2,
  nbz = 2,
  ! domain lenths in m
  lx_dim =    0.0, 37.2e3
  ly_dim =    0.0, 37.2e3
  lz_dim =    0.0, 82.0e2
  ! nb of processors in x and y direction must be set in the batch file
  nprocx = {nprocx},
  nprocy = {nprocy},

&end


!------------------------------------------------
!                Number of variables
!------------------------------------------------

&variables

  nVar = 8,              ! number of dependent variables
                         ! 7 = turbulent diffusivity from Dynamc Smagorinsky
                         ! 8 = heating term in the pressure solver
                         ! (due to GWs, only needed in WKB simulations)
  nOptVar = 4,           ! to switch on ice physics just set include_ice=.true.
                         ! and check that varIn, varOut and offset have the
                         ! right lengths
  include_ice = .false., ! include ice microphysics parametrization
                         ! automatically overwrites nVar, varOut, varIn and
                         ! offset by including additional dynamic variables
                         ! nIce, qIce and SIce

&end


!-----------------------------------------------------------------
!                          Model equations
!-----------------------------------------------------------------

&modelList

  model = "Boussinesq"            ! pseudo_incompressible / Boussinesq / WKB
  vert_theta = 90.0               ! deg    angle of rotation about y
  vert_alpha = 0.0                ! det    angle of rotation about z'

&end


!------------------------------------------------
!                   Pinc Solver
!------------------------------------------------

&solverList

  cfl = 0.5
  cfl_wave = 0.5                  ! passage rate of phase throuh a cell
  dtMax_dim = 1.e3                ! max time step in s
  tStepChoice = "cfl"             ! "fix" -> time step dtMax_dim is taken
                                  ! "cfl" -> stability criteria used
  timeScheme = "LS_Will_RK3"      ! LS_Will_RK3 -> Williamson / Euler /
                                  ! LS_TVD_RK3 / CL_TVD_RK3 / semiimplicit
  auxil_equ = .false.             ! auxiliary equation for the density
                                  ! fluctuations to be used in the explicit
                                  ! integration of the
                                  ! pseudo-incompressible system
  fluxType   = "upwind"           ! ILES / central / upwind
  reconstType = "MUSCL"           ! MUSCL / constant / SALD / ALDM
  limiterType1 = "MCVariant"      ! minmod / Cada / MCVariant
  fluctuationMode = .true.        ! use rho' as primary variable
  n_shap = 4                      ! (half) order of the shapiro filter
  shap_dts_fac = -1.              ! horizontal-Shapro-filter damping time
                                  ! scale (s < 0 means no filter)
  TurbScheme = .false.            ! Turbulence Schwme
  turb_dts = 5.e3                 ! (s) turbulent damping time scale for the
                                  ! smallest grid scales
  DySmaScheme = .false.           ! Dynamic Smagorinsky Scheme for the
                                  ! dynamic calculation of the turbulent
                                  ! damping time scale
  dtWave_on = .true.              ! .true. : include dtWave = pi/N to time
                                  ! step choice
  heatingONK14 = .false.          ! pseudo-incompressible dynamics with
                                  ! heating as in ONeill and Klein (2014)
  dens_relax = .false.            ! switch for replacement of relaxational
                                  ! heating by density relaxation

&end


!-------------------------------------------------
!                  Poisson solver
!-------------------------------------------------

&poissonSolverList
  tolcrit = "abs"                ! tolerance criterion:
                                 ! "abs" for absolute, faster convergence
                                 ! "rel" for relative tol., more demanding
  tolscale = .true.              ! scaling of the tolerance.
                                 ! If tolcrit = "abs" and tolscale = .true.,
                                 ! tolPoisson = tolPoisson*alpha,
                                 ! where alpha is dynamically calculated
                                 ! magnitude of gradients.
  tolPoisson = 1.0e-8            ! abort criterion
  tolCond = 1.e-23               ! tolerance value controlling the use of
                                 ! the preconditioner
  abs_tol = 0.                   ! it is unscaled abs. tol.,
                                 ! lower bound for tolerance.
  maxIterPoisson = 5000
  poissonSolverType = "bicgstab" ! "bicgstab" / "gcr" / "adi" / "hypre"
  storageType = "opr"            ! "csr" (compressed sparse row)
                                 !  "opr" (lin operator)
  preconditioner = "yes"          ! for operator-Solver: "no" / "yes"
  dtau = 4.0e-4                  ! time parameter for ADI (imperical value)
  maxIterADI = 2                 ! nb of iterations for ADI preconditioner
  initialCleaning = .true.       ! makes initial projection
  pressureScaling = .false.      ! .true. / .false. Scaling with PStrat
  useNAG = .false.               ! use NAG routine for TDMA algorithm
  correctMomentum = .true.       ! turn velocity projection on/off
  correctDivError = .false.      ! true -> subtract rho*div(u)

&end


!------------------------------------------------
!                    Atmosphere
!------------------------------------------------

&atmosphereList

  referenceQuantities = "Klein"  ! Klein / WKB / SI / general

  specifyReynolds = .false.            ! false -> give mu_viscous,
                                       ! true-> give ReInv
  ReInv = 0.                           ! inverse Reynolds number,
                                       ! ReInv = 0 -> inviscid flow
  mu_viscous_dim = 0.                  ! m^2/s
                                       ! kinematic viscosity: 0 for inviscid
                                       ! 1.5e-5 for z = 0 at bottom of
                                       ! atmosphere
                                       ! 1.5e-2 for z = 0 at 60km
  mu_conduct_dim = 0.                  ! m^2/s
                                       ! heat conductivity:
                                       ! 0 for non-diffusive
                                       ! 2 * mu_viscous_dim corresponds to
                                       ! Pr = 0.5
  background = "stratified_Boussinesq" ! const-N    -> set N_BruntVaisala_dim
                                       ! isothermal -> set Temp0_dim in K
                                       ! isentropic -> set theta0_dim in K
                                       ! uniform    -> constant density
                                       !               (Boussinesq)
                                       ! realistic  -> isentr. troposphere /
                                       !               isoth. stratosphere
                                       ! baroclinic  -> diff lapse rates in
                                       !                atmosphere and
                                       !                stratosphere
                                       ! diflapse  -> diff lapse rates in
                                       !              atmosphere and
                                       !              stratosphere,
                                       !              linear decline in
                                       !              temperature
                                       ! HeldSuarez -> reference atmosphere for
                                       !               Held & Suarez (1994) with
                                       !               decreasing temperature in
                                       !               tropoposphere and
                                       !               isothermal stratosphere
                                       !
                                       ! for realistic background...
  z_tr_dim = 12000.0                   ! m
                                       ! height of tropopause
                                       ! (need for baroclinic)
  theta_tr_dim = 300                   ! K
                                       ! const pot. temp. of troposphere
                                       ! (need for baroclinic)
  theta0_dim = 300                     ! K
                                       ! isentropic -> background pot. temp.
                                       ! const-N    -> ground pot. temp.
                                       ! uniform    -> background pot temp for
                                       !               Boussinesq
  Temp0_dim =  300                     ! K
                                       ! isothermal -> background temperature
  press0_dim =  100000.0               ! ground pressure (at z=0) in Pa:
                                       ! 101325.0 for z = 0 bottom of atmosphere
                                       ! 101.3250 for z = 0 at appr 60km
  N_BruntVaisala_dim = 0.01 !0.01            ! Brunt-Vaisala frequency for
                                       ! 1) "const-N" atmosphere in 1/s
                                       ! 2) "unifrom" Boussinesq
  backgroundFlow_dim =  4.0, 0.0, 0.0  ! m/s
                                       ! zonal background flow velocity u
  f_Coriolis_dim = 0.                  ! 1/s
                                       ! Coriolis parameter
  corset = "constant"                  ! constant/ periodic
  gamma_t = 0.000                      ! lapse rate in the troposphere
  gamma_s = 0.000                      ! lapse rate in the stratosphere
  tp_strato_dim = 200.                 ! stratosphere temperature (K)
                                       ! (used by HeldSuarez)
  tp_srf_trp_dim = 315.                ! tropical surface temperature (K)
                                       ! (used by HeldSuarez)
  tpdiffhor_tropo_dim = 60.            ! tropospheric temperature difference
                                       ! between poles and tropics (K)
                                       ! (used by HeldSuarez)
  ptdiffvert_tropo_dim = 10.           ! vertical potential-temperature
                                       ! difference in troposphere (K)
                                       ! (used by HeldSuarez)

&end


!-----------------------------------------------------------
!          Bottom topography
!-----------------------------------------------------------

&topographyList

  topography = .true.       ! switch for bottom topography
  testTFC = .false.         ! switch for TFC test
  spongeTFC = .true.        ! switch for unified sponge layer
  lateralSponge = .false.   ! switch for lateral sponge layers
  mountainHeight_dim = 1.e2 ! mountain height in m
  mountainWidth_dim = 3.e3  ! mountain half-width in m
  mountain_case = 4         ! shape of orography
                            ! 1 for cosine-shaped mountains
                            ! 2 for 3D cosine-shaped mountains (rotated)
                            ! 3 for bell-shaped mountain
                            ! 4 for 3D bell-shaped mountain
                            ! 5 for wave-packet-like mountain range
  range_fac = 10            ! factor by which mountain range is wider than
                            ! single mountains

&end


!-----------------------------------------------------------
!                          Boundary
!-----------------------------------------------------------

&boundaryList

  ! correction of solid wall boundary
  rhoFluxCorr = .false.       ! replace vertical mass flux by CDS at k=1,
                              ! k=nz-1
  uFluxCorr = .false.
  vFluxCorr = .false.
  wFluxCorr = .false.
  thetaFluxCorr = .false.     ! replace vertical theta flux by CDS at k=1,
                              ! k=nz-1
  nbCellCorr = 1
  ! sponge layer at upper boundary
  spongeLayer = .true.       ! sponge with relaxation to background
  spongeHeight = 0.5         ! relative height of sponge layer
  spongeAlphaZ_dim = 0.05 !1.4e-4 ! relaxation rate coeff in 1/s
  spongeAlphaZ_fac = 1.0

&end

&boundaryList2

  ! boundary types
  xBoundary = "periodic"   ! periodic
  yBoundary = "periodic"   ! periodic
  zBoundary = "solid_wall" ! periodic / solid_wall

&end


!------------------------------------------------
!                   Input / Output
!------------------------------------------------

&outputList

  outputType = "time"         ! timeStep / time
  nOutput = 1                 ! output every nOutput's time step
                              ! for outputType = "timeStep"
  maxIter = 10000             ! stop after maxIter time steps
  outputTimeDiff = 3.6e3      ! output every ... seconds
  maxTime = 86.4e3            ! stop after maxTime seconds
  dataFileName = ""           ! empty string "" -> dataFileName = testCase
  restartFile = "restart.ref" ! restart file in TEC360 format
  restart = .false.           ! true / false
  dimOut = .true., .false., .true.
                              ! 2D(x,z)-plot dimOut = 1,0,1, 3D with 1,1,1
  varOut = 1,1,1,1,1,1,1,0    ! 1 = output, 0 = no output
                              ! primary variables: rho,u,v,w,pi',theta',
                              ! dyn. Smagorinsky coeff.
                              ! if include_ice varOut must have length nVar+3
  varIn = 1,1,1,1,1,1,1,0     ! 1 = output, 0 = no output
                              ! data written into restart file pf_all_in.dat
                              ! ( = output file pf_all.dat from previous run)
                              ! primary variables: rho,u,v,w,pi',theta',
                              ! dyn. Smagorinsky coeff.
                              ! if include_ice varOut must have length nVar+3
  iIn = 97                    ! no. of record to be read from restart file
                              ! pf_all_in.dat
                              ! (first record in file has no. = 0)
  offset = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                              ! offset for primary variables
                              ! if include_ice offset must have length nVar+1,
                              ! alse nVar-2
  rhoOffset = .false.         ! subtract background
  ! optional variables
  optVarOut = 0,0,0, 0,0,0    ! 1 = output, 0 = no output for
                              ! 1) p-pBar in kPa
                              ! 2) buoyancy
                              ! 3) background pressure pBar in kPa
                              ! 4) background density rhoBar in kg/m^3
                              ! 5) div(Pu)
                              ! 6) stratification perturbation db/dz
  thetaOffset = .false.       ! subtract background
  ! WKB variables
  wkbVarOut = 0,0,0,0, 0,0,0,0, 0,0,0,0
                              ! 1) Wave action amplitude
                              ! 2) u00 in m/s
                              ! 3) th02
                              ! 4) pi02
                              ! 5) u10 in m/s
                              ! 6) w10 in m/s
                              ! 7) b11 in m/s^2
                              ! 8) pi12
                              ! 9) u21 in m/s
                              ! 10) w21 in m/s
                              ! 11) b22 in m/s^2
                              ! 12) pi23
  solutionTime = .true.       ! TECPLOT's "solution time" out yes/no
  solutionTimeUnit = "min"    ! "s", "min" or "h"
  showGhostCellsX = .false.   ! shows ghost cells in x
  showGhostCellsY = .false.   ! shows ghost cells in y
  showGhostCellsZ = .false.   ! shows ghost cells in z
  detailedinfo = .false.      ! provide info on the final state of Poisson
                              ! solver
  RHS_diagnostics = .true.    ! provide info about the RHS for Poisson problem

&end


!------------------------------------------------
!                 Parameter study
!------------------------------------------------

&parameterList

  parameterStudy = .false. ! .true. / .false.
  startParam = 1
  endParam = 10        stepParam = 1
  paramName = ""

&end


!------------------------------------------------
!           Debugging & Error handling
!------------------------------------------------

&debuggingList

  verbose = .false.
  dtMin_dim = 1.0e-5 ! stop program if dt < dtMin

&end


!------------------------------------------------
!                 WKB & Ray tracing
!------------------------------------------------

&wkbList

  rayTracer = .false. ! set up ray tracer
  nRayRatioX = 1      ! reduce amount of storage if raytracer is off
  nRayRatioY = 1
  nRayRatioZ = 1

&end


!------------------------------------------------
!  Testcases: initial and boundary conditions
!------------------------------------------------

! general
&testCaseList

  testCase = "mountainwave" ! Boussinesq: uniform_theta, wavePacket
                            ! agnesiMountain -> see topography
                            ! baroclinic_LC -> baroclinic life cycle with y-dep
                            ! tropopause
                            ! baroclinic_ID -> idealistic baroclinic life cycle
                            ! with const tropopause

&end

&monochromeWave

  lambdaZ_dim = 6000.0 ! m
                       ! vertical wave length

&end


!----------------
! Wave packets
!----------------

! Test cases: wavePacket, wavePacket_raytracer

&wavePacket

  wavePacketType = 1        ! 1 = Gaussian, 2 = Cosine
  wavePacketDim = 2         ! 1 = 1D, 2 = 2D, 3 = 3D
                            ! for a 2.5D Wave Packet use wavePacketDim = 2
  lambdaX_dim = 300000.0    ! wave length in x direction in m
                            ! lambdaX = 0.0 --> infinite wavelength
  lambdaY_dim = 0.0         ! wave length in y direction in m
                            ! lambday = 0.0 --> infinite wavelength
  lambdaZ_dim = -1000.0     ! vertical wave length in m
  amplitudeFactor = 0.0     ! normalilized buoyancy amplitude
  xCenter_dim = 450000.0    ! center of wave packet in x direction in m
  yCenter_dim = 2.e4        ! center of wave packet in y direction in m
  zCenter_dim = 30000.0     ! center of wave packet in z direction in m
  sigma_dim = 5000.0        ! vertical width of Gaussian wavepacket in m
  sigma_hor_dim = 1500000.0 ! cosine distribution width
                            ! (in x direction, 0 means infinity)
  sigma_hor_yyy_dim = 0.0   ! cosine distribution width
                            ! (in y direction, 0 means infinity)
  amp_mod_x = 1.e0          ! fractional amplitude of amplitude modulation
                            ! in x direction
                            ! (0 = no modulation, 1 = total modulation)
  amp_mod_y = 1.0           ! fractional amplitude of amplitude modulation
                            ! in y direction
                            ! (0 = no modulation, 1 = total modulation)
  L_cos_dim = 10000.0       ! half width of vertical cosine profile of GWP
  meanFlowX_dim = 0.0       ! mean flow in m/s / jet flow amplitude
  meanFlowZ_dim = 0.0       ! mean vertical flow m/s
  u0_jet_dim = 0.0          ! amplitude max of jet velocity
  L_jet_dim = 5000.0        ! half width vertical cosine profile of jet in 1m
  z0_jet_dim = 50000.0      ! center of jet stream
  omiSign = 1               ! frequency branch

&end


!-----------------------------------------------
!   WKB simulations wave packet or mountain wave
!-----------------------------------------------

&LagrangeRayTracing

  xrmin_dim = 0.0,      ! left bound of initial rays (in x direction) (m)
  xrmax_dim = 5.e6,     ! right bound of initial rays (in x dir.) (m)
  yrmin_dim = 0.0,      ! left bound of initial rays (in y direction) (m)
  yrmax_dim = 5.e4,     ! right bound of initial rays (in y dir.) (m)
  zrmin_dim = 3000.0,   ! bottom bound of initial rays (m)
  zrmax_dim = 70000.0,  ! top bound of initial rays (m)
  nrxl = 1,             ! no. of ray vol. init. within one hor. x column
  nryl = 1,             ! no. of ray vol. init. within one hor. y column
  nrzl = 1,             ! no. of ray vol. init. within one vert. layer
  fac_dk_init = 0.1,    ! init. width of total ray vol. in k space
                        ! (fraction of the initial wave number in x dir.)
  fac_dl_init = 0.1,    ! init. width of total ray vol. in l space
                        ! (fraction of the initial wave number in y dir.)
  fac_dm_init = 0.1,    ! init. width of total ray vol. in m space
                        ! (fraction of the initial vert. wave number)
  nrk_init = 2,         ! no. of ray volumes initialized within dk
  nrl_init = 1,         ! no. of ray volumes initialized within dl
  nrm_init = 2,         ! no. of ray volumes initialized within dm
  nsmth_wkb = 2,        ! half (number -1) of cells f. smooth. wkb fluxes
  lsmth_wkb = .true.,   ! log. switch for smooth. wkb data (true/false)
  sm_filter = 2,
  lsaturation = .true., ! JaWi 16.12.16 (sat)
  alpha_sat = 1.0,      ! JaWi 16.12.16 (sat)
  case_wkb = 3,         ! 1/2: Gaussian/Cosine wave packet; 3: mountain
  amp_wkb = 0.5         ! amplitude of the wave packet (wrt saturation)
  wlrx_init = 1.e5,     ! initial lambda_x of the wave packet (m)
  wlry_init = 0.0       ! initial lambda_y of the wave packet (m)
                        ! (0 means infinity)
  wlrz_init = 1.e3,     ! initial lambda_z of the wave packet (m)
                        ! (0 means infinity)
  xr0_dim = 2.5e6,      ! center of the wave packet in hor. (x-dir.) (m)
  yr0_dim = 2.5e4,      ! center of the wave packet in hor. (y-dir.) (m)
  zr0_dim = 3.e4,       ! center of the wave packet in vertical (m)
  sigwpx_dim = 1.e6,    ! width of the wave packet in hor. (x-dir.) (m);
                        ! (0 means infinity)
  sigwpy_dim = 0.e0     ! width of the wave packet in hor. (y-dir.) (m);
                        ! (0 means infinity)
  sigwpz_dim = 5.e3,    ! width of the wave packet in vertical (m);
  branchr = 1,          ! frequency branch (dispersion relation)
  ! presently not used
  lindUinit = .false.,  ! ind. wind already at initial time (true/false)
  ! oror_amp_dim = 50   ! orography amplititude height (m)
  oror_amp_dim = 5.e2   ! orography amplititude height (m)

  zmin_wkb_dim = 0.0    ! minumum altitude (above the model bottom, in m)
                        ! for WKB wave-mean-flow interaction
                        ! (zmin_wkb > 0 can help preventing the
                        ! ray volumes being trapped by the self-induced
                        ! mean wind)
  nray_fac = 5          ! maximum factor (per wavenumber direction) by
                        ! which # of rays may increase in comparison to
                        ! initialization
  cons_merge = "en"     ! quantity to be conserved
                        ! ("wa" = wave action/ "en" = wave energy)
                        ! under ray-volume merging

&end


!----------------------------------------------------
!   Wind relaxation for explicit or WKB mountain wave
!----------------------------------------------------

! test case initializes with zero flow that is rammped up gradually by a
! temporary wind relaxation.
! topography = .true. must be set separately (see above)!
! topography = .false. allows testing the temporary wind relaxation.

&mountainwavelist

  u_relax = 10.           ! [m/s] zonal wind to be attained by
                          ! temporary wind relexation
  t_relax = 5.e3          ! [s] total relaxation time
  t_ramp = 25.e2          ! [s] duration of ramping up/down the relaxation
  xextent_norelax = 1.e5  ! [m] zonal extent of
                          ! region without wind relaxation

&end


!----------------
!    Bubbles
!----------------

! Test cases:  ->  Robert_Bubble

&Robert_Bubble

  dTheta1_dim = 0.5    ! K pot temp offset
  a1_dim = 150.0       ! m radius of plateau
  sigma1_dim = 50.0    ! m Gaussian edge profile
  xCenter1_dim = 0.0   ! m zCenter1_dim = 300.0 !m
  zCenter1_dim = 300.0 ! m
  dTheta2_dim = -0.15  ! K pot temp offset
  a2_dim = 0.0         ! m radius of plateau
  sigma2_dim = 50.0    ! m Gaussian edge profile
  xCenter2_dim = 60.0  ! m
  zCenter2_dim = 640.0 ! m

&end

! Test cases: hotBubble, coldBubble, hotBubble3D

&bubble

  dTheta0_dim = 6.0    ! K
  xRadius_dim = 4000.0 ! m
  zRadius_dim = 2000.0 ! m
  xCenter_dim = 0.0    ! m
  zCenter_dim = 3000.0 ! m
  zExcentricity = 1.0

&end


!-----------------------------------------------------------------
!     baroclinic-wave life cycles
!
! either setup from Kuehnlein et al (2012): background = 'const-N'
! or setup of Held & Suarez (1994): background = 'HeldSuarez'
!-----------------------------------------------------------------

&baroclinic_LC

  ! for baroclinic test case

  ! used be Kuehnlein et al (2019):
  z_trpp0_dim = 8.e3            ! mean tropopause height (m)
  z_baro_dim = 2.4e4            ! alt. above which the atmo. is barotr. (m)
  deltht_dim = 3.e1             ! meridional potential-temp. contrast (K)
  thet0_dim  = 2.73e2           ! characteristic potential temperature (K)
  ntrp_dim = 1.e-2              ! Brunt-Vaisala frequency troposphere (1/s)
  nstr_dim = 2.45e-2            ! Brunt-Vaisala frequency stratosphere (1/s)
  kaptpp = 7.e1                 ! jet slope
  tau_relax = 864000.           ! relax. time, relax to the env. state if
                                ! RelaxHeating = 1, 2, [sec]
  tau_relax_low = 864000.       !432000.!345600.! 432000.!777600.432000.! !9d
  sigma_tau = 0.02              ! relative thickness of tau_relax(y) profile
  tau_jet = 7200.               ! time for a jet formation, used if
                                ! RelaxHeating = 2, [sec]
  ! used by Held & Suarez (1994):
  ta_hs_dim = 3.456e6           ! thermal-relaxation time scale outside
                                ! tropical boundary layer (s)
                                ! zero means infinty
  ts_hs_dim = 3.456e5           ! thermal-relaxation time scale in
                                ! tropical boundary layer (s)
                                ! zero means infinty
  tf_hs_dim = 8.64e4            ! boundary-layer Rayleigh-damping time
                                ! scale (s)
                                ! zero means infinty
  sigb_hs = 0.7                 ! sigma of boundary-layer top
  ! used both by Kuehnlein et al (2012) and Held & Suarez (1994):
  jwdth_dim = 2.4e6             !5.e6
                                ! half (!) jet width (m)
  add_ptptb = .true.            ! add local potential-temperature
                                ! perturbation
  ptptb_x_dim = 2.5e6           ! x coordinate of local
                                ! potential-temperature perturbation (m)
  ptptb_y_dim = 5.e6            !7.5e6 !10.e6
                                ! y coordinate of local
                                ! potential-temperature perturbation (m)
  ptptb_z_dim = 9.e3            ! z coordinate of local
                                ! potential-temperature perturbation (m)
  ptptb_dh_dim = 8.e5           !6.e5
                                ! horizontal width of local
                                ! potential-temperature perturbation  (m)
  ptptb_dz_dim = 2.e3           ! vertical width of local
                                ! potential-temperature perturbation  (m)
  ptptb_amp_dim = 3.0           !3.0
                                ! amplitute of local potential-temperature
                                ! perturbation (K)
  RelaxHeating = 1              !0 !1
                                ! "1" if sharp relax to environmental state
                                !     - heating sourse - RHS in divergence
                                !     constraint, 1/tau
                                ! "2" if smooth relax to env. state
                                !     - heating sourse - RHS in divergence
                                !     constraint, cos(...)/tau
                                ! "0" no relaxation
  Sponge_Rel_Type = "constant"  ! "constant" - in the sponge the
                                !              relaxation to the env state
                                !              is still with tau_relax
                                ! "linear" - in the sponge the relaxation
                                !          to the env state is linear
  Sponge_Rel_Bal_Type = "env"  ! the state to which the relaxation goes i
                                ! in the sponge
                                ! "hyd" - in the sponge the relaxation
                                !         goes to hydrostat. bal. state
                                ! "env" - in the sponge the relaxation
                                !         goes to geostr. bal. state
                                !  otherwise - in the sponge the
                                !              development goes freely for
                                !              u and \rho
  add_noise = .false.           ! add noise to the initial potential
                                ! temperature
  proc_noise = 0.01             ! amplitude of noise is
                                ! proc_noise*theta_tr_dim
  dTh_atm = 30.                 ! K, potential temperature deviation
                                ! between poles and equator

&end

&baroclinic_ID

  ! for idealsitic baroclinic test case

  init_2Dto3D = .false.                  ! initialize 3D state from 2D
                                         ! balanced state,
                                         ! zonally-uniform
  fileinitstate2D = 'pf_all_2D3D_in.dat' ! filename with 2D balanced state
  lastrecordnum = 88                     ! last record in 2D output called
                                         ! pf_all_in.dat
  bar_sigma_y = 0.10
  dTh_atm = 20.                          ! K,
                                         ! potential temperature deviation
  u_strength = 30.                       ! m/s: strength of jet
  tau_relax = 1296000.                   ! relaxation tim, erelax to the
                                         ! environm. state
                                         !  if RelaxHeating = 1, 2, [sec]
  tau_jet = 7200.                        ! time for a jet formation,
                                         ! used if RelaxHeating = 2, [sec]
  RelaxHeating = 0                       ! "1" if sharp relax to
                                         ! environmental state
                                         ! - heating sourse - RHS in
                                         ! divergence constraint, 1/tau
                                         ! Switched off, not tested well:
                                         ! "2" if smooth relax to
                                         ! environmental state
                                         ! - heating sourse - RHS in
                                         ! divergence constraint,
                                         ! cos(...)/tau
                                         ! "0" no relaxation
  init_bal = "geostr_id"                 ! "hydrost" if initialize with
                                         ! hydrostatically balanced state,
                                         ! "geostr_id"  if initialize with
                                         ! hydrostatically+geostrophically
                                         ! balanced state
  balance_eq = "PI_an"                   ! QG for quasigeostrophic thermal
                                         ! wind balance,
                                         ! PI for pseudoincompressible one
  Sponge_Rel_Type = "constant"           ! "constant" - in the sponge the
                                         ! relaxation to the env state is
                                         ! still with tau_relax
                                         ! "linear" - in the sponge the
                                         ! relaxation to the env state is
                                         ! linear
  Sponge_Rel_Bal_Type = "env"            ! the state to which the
                                         ! relaxation of u goes in the
                                         ! sponge
                                         ! "hyd" - in the sponge the
                                         ! relaxation goes to hydrostat.
                                         ! bal. state, u = 0
                                         ! "env" - in the sponge the
                                         ! relaxation goes to geostr. bal.
                                         ! state
                                         ! otherwise - in the sponge the
                                         ! development goes freely for u
                                         ! and \rho
                                         ! this functionality is switched
                                         ! off since the other options are
                                         ! not used
  add_noise = .true.                     ! add noise to the initial
                                         ! potential temperature
  proc_noise = 0.01                      ! amplitude of noise is
                                         ! proc_noise*dTh_atm*2
  output_theta_bgr = .false.             ! output environmental state
  output_rho_bgr = .true.                ! output environmental state
  output_br_vais_sq = .true.             ! output N^2 of environmental state
  output_heat = .true.

&end
