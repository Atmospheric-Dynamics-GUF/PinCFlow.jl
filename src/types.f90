module type_module

  !-----------------------------------------------------------------
  !    Definition of data types and variables accessible
  !    throughout the code! Use care when working on these
  !    data fields.
  !-----------------------------------------------------------------

  !-----------------------------------------------------------
  implicit none
  !-----------------------------------------------------------

  public
  ! all variables declared herein are public by default

  !-----------------------------------------------------------------
  !                     MPI & Domain composition
  !-----------------------------------------------------------------

  integer :: sizeX, sizeY, sizeZ
  integer :: nbx, nby, nbz
  integer :: nprocx, nprocy
  real, dimension(0:1) :: lx_dim, ly_dim, lz_dim ! dimensional domain

  namelist / domain / sizeX, sizeY, sizeZ, nbx, nby, nbz, lx_dim, ly_dim, &
      lz_dim, nprocx, nprocy

  integer :: sizeXX, sizeYY, sizeZZ
  integer :: nx1, in1, ny1, jn1
  integer :: iStart, jStart
  integer :: is, ie, js, je ! local start and end indices
  logical :: verboseMPI

  ! MPI include (parameters needed below)
  include 'mpif.h'

  ! MPI variables
  integer :: ierror
  integer, dimension(2) :: dims, coords
  logical, dimension(2) :: periods
  integer :: back, forw, right, left, rank
  integer :: idim, jdim, icoord, jcoord
  integer :: nbProc, comm
  logical :: master
  integer, dimension(mpi_status_size) :: sts_left, sts_right, sts_back, sts_forw

  integer, parameter :: root = 0

  !-----------------------------------------------------------------
  !                      (local) Grid & Domain
  !-----------------------------------------------------------------

  integer :: nx, ny, nz
  real, dimension(0:1) :: lx, ly, lz ! scaled domain

  real, dimension(:), allocatable :: x, y, z
  real :: dx, dy, dz
  integer :: nxx, nyy, nzz ! grid size inclusive ghost cells
  integer :: nxyz ! = nx*ny*nz

  !--------------------------------------------------------------
  !              for output global fields to single file
  !--------------------------------------------------------------

  real * 4, dimension(:, :), allocatable :: field_out, field_mst

  !-----------------------------------------------------------------
  ! for sponge: maximum damping rate in 1/dt
  !-----------------------------------------------------------------

  real :: alpspg

  !-----------------------------------------------------------------
  ! for
  ! (1) Rayleigh damping in land cells and
  ! (2) density-fluctuation relaxation in semi-implicit time stepping
  !-----------------------------------------------------------------

  real :: alprlx

  !-----------------------------------------------------------------
  !                          Variables
  !-----------------------------------------------------------------
  integer :: nVar, nOptVar, iVarT, nBscVar
  logical :: include_ice ! controls use of additional ice variables nAer,nIce,qIce and qv
  logical :: include_ice2 = .false., include_testoutput = .false.
  logical :: include_tracer = .false.
  namelist / variables / nVar, nOptVar, include_ice, include_ice2, &
      include_tracer, include_testoutput

  !-----------------------------------------------------------------
  !                    Input / Output variables
  !-----------------------------------------------------------------
  integer :: iOut ! output counter ; gagarina: moved from pinc

  character(len = 256) :: file_namelist
  character(len = 40) :: dataFileName
  character(len = 40) :: restartFile ! Tecplot file with restart data
  logical, dimension(3) :: dimOut ! (/1,0,1/) = 2D (x,z), (/1,1,1) = 3D
  integer, dimension(:), allocatable :: varOut ! 1 = output, 0 = no output

  ! achatzb
  ! data available in restart file
  integer, dimension(:), allocatable :: varIn ! 1 = available, 0 = not available

  ! record to be read from restart file (starting from 0!)
  integer :: iIn
  ! achatze

  real, dimension(:), allocatable :: offset ! offset for rho, u,v,w, p, nAer, nIce, qIce, qv
  integer, dimension(6) :: optVarOut ! 1 = output, 0 = no output
  integer, dimension(13) :: wkbVarOut !           --"--
  ! increase dimension of optVarOut for new optional variables

  ! subtract background yes/no
  logical :: thetaOffset
  logical :: rhoOffset

  character(len = 20) :: outputType ! "time" or "timeStep"
  integer :: nOutput ! output every nOutput's time step
  integer :: maxIter ! max nb. of time steps
  real :: outputTimeDiff ! output every ... seconds
  real :: maxTime ! max time in seconds
  logical :: solutionTime ! TECPLOT's "solution time"
  character(len = 3) :: solutionTimeUnit !
  logical :: restart ! reads restartFile in TECPLOT format
  logical :: showGhostCellsX ! plots include ghost cells in x-direction
  logical :: showGhostCellsY ! y-direction
  logical :: showGhostCellsZ ! z-direction

  logical :: detailedinfo
  logical :: RHS_diagnostics

  logical :: PVinversion

  namelist / outputList / dataFileName, restartFile, dimOut, varOut, varIn, &
      iIn, offset, optVarOut, wkbVarOut, outputType, nOutput, maxIter, &
      outputTimeDiff, maxTime, solutionTime, solutionTimeUnit, restart, &
      showGhostCellsX, showGhostCellsY, showGhostCellsZ, thetaOffset, &
      rhoOffset, detailedinfo, RHS_diagnostics, PVinversion
  !achatzb
  !achatze

  !-----------------------------------------------------------------
  !                         Parameter study
  !-----------------------------------------------------------------

  logical :: parameterStudy ! .true. / .false.
  integer :: startParam, endParam, stepParam
  character(len = 50) :: paramName
  namelist / parameterList / parameterStudy, startParam, endParam, stepParam, &
      paramName

  !-----------------------------------------------------------------
  !                    Debugging & Error handling
  !-----------------------------------------------------------------

  logical :: verbose
  real :: dtMin_dim
  namelist / debuggingList / verbose, dtMin_dim

  !-----------------------------------------------------------------
  !                           Test cases
  !-----------------------------------------------------------------

  ! general
  character(len = 50) :: testCase
  namelist / testCaseList / testCase

  !------------------------
  !  monochromatic wave
  !------------------------

  real :: lambda_dim ! wave length in (inclined) z direction

  namelist / monochromeWave / lambda_dim

  ! modified by Junhong Wei for 3DWP (20170828) *** starting line ***

  !-------------------------------------------------------
  !         Gravity wave packet: Wave resolving simulation
  !-------------------------------------------------------
  ! test cases:
  ! 1) wavePacket1D and wavePacket1D_raytracer
  ! 2) wavePacket2D and wavePacket2D_raytracer
  integer :: wavePacketDim
  integer :: wavePacketType
  real :: lambdaX_dim, lambdaY_dim, lambdaZ_dim
  real :: meanFlowX_dim, meanFlowZ_dim
  real :: amplitudeFactor
  real :: sigma_dim
  real :: sigma_hor_dim, sigma_hor_yyy_dim ! modified by Junhong Wei (20170214)
  real :: L_cos_dim
  integer :: omiSign
  real :: u0_jet_dim
  real :: z0_jet_dim
  real :: L_jet_dim
  real :: xCenter_dim, yCenter_dim, zCenter_dim
  ! achatzb
  real :: amp_mod_x, amp_mod_y
  logical :: inducedwind
  ! achatze

  namelist / wavePacket / wavePacketType, wavePacketDim, lambdaX_dim, &
      lambdaY_dim, lambdaZ_dim, meanFlowX_dim, meanFlowZ_dim, amplitudeFactor, &
      xCenter_dim, yCenter_dim, zCenter_dim, sigma_dim, sigma_hor_dim, &
      amp_mod_x, sigma_hor_yyy_dim, amp_mod_y, L_cos_dim, omiSign, u0_jet_dim, &
      z0_jet_dim, L_jet_dim, inducedwind
  ! achatzb
  ! achatze
  ! achatzb
  ! achatze

  !--------------------------------------------------------------
  ! Lagrangian ray tracer: Wave packet or mountain wave
  !--------------------------------------------------------------

  real :: xrmin_dim
  real :: xrmax_dim
  real :: yrmin_dim
  real :: yrmax_dim
  real :: zrmin_dim
  real :: zrmax_dim

  integer :: nrxl
  integer :: nryl
  integer :: nrzl

  real :: fac_dk_init
  real :: fac_dl_init
  real :: fac_dm_init

  integer :: nrk_init
  integer :: nrl_init
  integer :: nrm_init

  real :: dk_init
  real :: dl_init
  real :: dm_init

  integer :: nsmth_wkb
  logical :: lsmth_wkb
  integer :: sm_filter ! 1 = box filter/ 2 = Shapiro filter

  logical :: lsaturation ! JaWi 16.12.16 (sat)
  real :: alpha_sat ! JaWi 16.12.16 (sat)

  logical :: steady_state

  integer :: case_wkb !1/2: Gaussian/Cosine wave packet; 3: mountain

  real :: amp_wkb

  real :: wlrx_init
  real :: wlry_init ! added by Junhong Wei
  real :: wlrz_init

  real :: xr0_dim
  real :: yr0_dim
  real :: zr0_dim

  real :: sigwpx_dim
  real :: sigwpy_dim
  real :: sigwpz_dim

  integer :: branchr
  logical :: lindUinit

  real :: topographyTime_wkb ! FJFeb2023
  real :: mountainHeight_wkb_dim ! Jan Weinkaemmerer, 27.11.18
  real :: mountainWidth_wkb_dim ! FJJan2023
  integer :: mountain_case_wkb
  integer :: range_factor_wkb ! FJFeb2023

  ! Long number scaling of displacement (FJJan2023)
  logical :: blocking

  ! Number of wave modes (FJApr2023)
  integer :: nwm

  real :: zmin_wkb_dim, zmin_wkb

  integer :: nray_fac ! maximum factor by which the # of rays may increase
  ! compared to initialization

  character(len = 2) :: cons_merge ! quantity to be conserved
  ! ("wa" = wave action/
  !  "en" = wave energy)
  ! under ray-volume merging

  namelist / LagrangeRayTracing / xrmin_dim, xrmax_dim, yrmin_dim, yrmax_dim, &
      zrmin_dim, zrmax_dim, nrxl, nryl, nrzl, fac_dk_init, fac_dl_init, &
      fac_dm_init, nrk_init, nrl_init, nrm_init, nsmth_wkb, lsmth_wkb, &
      sm_filter, lsaturation, alpha_sat, steady_state, case_wkb, amp_wkb, &
      wlrx_init, wlry_init, wlrz_init, xr0_dim, yr0_dim, zr0_dim, sigwpx_dim, &
      sigwpy_dim, sigwpz_dim, branchr, lindUinit, topographyTime_wkb, &
      mountainHeight_wkb_dim, mountainWidth_wkb_dim, mountain_case_wkb, &
      range_factor_wkb, blocking, nwm, zmin_wkb_dim, nray_fac, cons_merge ! JaWi: new nml!
  ! Jan Weinkaemmerer, 27.11.18

  !------------------------------------------
  ! hotBubble, coldBubble, hotBubble3D
  !------------------------------------------

  real :: dTheta0_dim, xRadius_dim, zRadius_dim, rhoCenter_dim
  real :: zExcentricity
  namelist / bubble / dTheta0_dim, xRadius_dim, zRadius_dim, xCenter_dim, &
      zCenter_dim, zExcentricity, yCenter_dim, rhoCenter_dim

  ! hot and cold bubble by Robert
  real :: dTheta1_dim, a1_dim, sigma1_dim, xCenter1_dim, zCenter1_dim
  real :: dTheta2_dim, a2_dim, sigma2_dim, xCenter2_dim, zCenter2_dim
  namelist / robert_bubble / dTheta1_dim, a1_dim, sigma1_dim, xCenter1_dim, &
      zCenter1_dim, dTheta2_dim, a2_dim, sigma2_dim, xCenter2_dim, zCenter2_dim

  ! achatzb
  !-----------------------------------------------------------------
  !                  wind relaxation for mountain wave
  !-----------------------------------------------------------------

  ! zonal wind to be attained by temporary wind relexation
  real :: u_relax, v_relax, w_relax
  ! total relaxation time
  real :: t_relax
  ! duration of ramping up/down the relaxation
  real :: t_ramp
  ! zonal extent of region without wind relaxation
  real :: xextent_norelax
  ! TFC FJ
  ! Switch for background relaxation.
  logical :: wind_relaxation
  ! FJMar2023
  ! Surface layer depth.
  real :: surface_layer_depth
  namelist / mountainwavelist / u_relax, v_relax, w_relax, t_relax, t_ramp, &
      xextent_norelax, wind_relaxation, surface_layer_depth
  ! achatze

  !-----------------------------------------------------------------
  !                          Baroclinic life cycle: realistic
  !-----------------------------------------------------------------

  !UAB
  logical :: zero_initial_state ! zero_initial_state = .true. means that
  ! the equilibrium state is determined but
  ! that then the code is initialized with
  ! zero winds and and zero density and
  ! pressure fluctuations (to investigate
  ! the effect of the potential-temperature
  ! relaxation)
  !UAE
  real :: z_trpp0_dim ! mean tropopause height (m)
  real :: z_baro_dim ! alt. above which the atmo. is barotropic (m)
  real :: deltht_dim ! meridional potential-temp. contrast (K)
  real :: thet0_dim ! characteristic potential temperature (K)
  real :: ntrp_dim ! Brunt-Vaisala frequency troposphere (1/s)
  real :: nstr_dim ! Brunt-Vaisala frequency stratosphere (1/s)
  real :: jwdth_dim ! jet width (m)
  real :: kaptpp ! jet slope

  real :: ptptb_x_dim ! x coordinate of local potential-temperature
  ! perturbation (m)
  real :: ptptb_y_dim ! y coordinate of local potential-temperature
  ! perturbation (m)
  real :: ptptb_z_dim ! z coordinate of local potential-temperature
  ! perturbation (m)
  real :: ptptb_dh_dim ! horizontal width of local potential-temperature
  ! perturbation  (m)
  real :: ptptb_dz_dim ! vertical width of local potential-temperature
  ! perturbation  (m)
  real :: ptptb_amp_dim ! amplitute of local potential-temperature
  ! perturbation (K)
  logical :: add_ptptb ! switch local potential-temperature perturbation

  !UAB
  real :: ta_hs_dim ! thermal-relaxation time scale outside tropical
  ! boundary layer (s, used by HeldSuarez)
  real :: ts_hs_dim ! thermal-relaxation time scale in tropical
  ! boundary layer (s, used by HeldSuarez)
  real :: tf_hs_dim ! boundary-layer Rayleigh-damping time scale
  ! (s, used by HeldSuarez)
  real :: sigb_hs ! sigma of boundary-layer top (used by HeldSuarez)
  !UAB

  ! total relaxation time
  real :: t_relax_bar
  ! relaxation in massUpdate for the density
  real :: tau_relax, tau_relax_low, sigma_tau
  ! time scale for a jet formation
  real :: tau_jet

  ! environmental state
  real, dimension(:, :, :, :), allocatable :: var_env
  real, dimension(:, :, :), allocatable :: p_env_pp, the_env_pp, dens_env_pp, &
      u_env_pp, v_env_pp

  ! relaxation rates for Held & Suarez (1994)
  real, dimension(:, :), allocatable :: kt_hs, kv_hs
  real, dimension(:), allocatable :: kw_hs

  !UAB
  ! relaxation rates for Held & Suarez (1994)
  real, dimension(:, :), allocatable :: kr_sp, kr_sp_w
  !UAE

  real, dimension(:, :, :), allocatable :: kt_hs_tfc
  real, dimension(:, :, :), allocatable :: kr_sp_tfc, kr_sp_w_tfc

  ! switch of thermal relaxation in the divergence constraint
  integer :: RelaxHeating

  character(len = 20) :: fileinitstate2D

  ! type of relaxation in the sponge
  character(len = 20) :: Sponge_Rel_Type, Sponge_Rel_Bal_Type

  ! add noise
  logical :: add_noise, init_2Dto3D

  real :: proc_noise, dTh_atm

  namelist / baroclinic_LC / zero_initial_state, z_trpp0_dim, z_baro_dim, &
      deltht_dim, thet0_dim, ntrp_dim, nstr_dim, jwdth_dim, kaptpp, &
      t_relax_bar, add_ptptb, ptptb_x_dim, ptptb_y_dim, ptptb_z_dim, &
      ptptb_dh_dim, ptptb_dz_dim, ptptb_amp_dim, dTh_atm, add_noise, &
      proc_noise, tau_relax, tau_relax_low, sigma_tau, tau_jet, RelaxHeating, &
      Sponge_Rel_Type, Sponge_Rel_Bal_Type, init_2Dto3D, fileinitstate2D, &
      ta_hs_dim, ts_hs_dim, tf_hs_dim, sigb_hs
  !UAC & z_trpp0_dim, z_baro_dim, deltht_dim, thet0_dim, &
  !UAE
  !UAB
  !UAE

  !-----------------------------------------------------------------
  !                          Baroclinic life cycle: idealized
  !-----------------------------------------------------------------

  ! define the jet parameters
  real :: bar_sigma_y, u_strength
  integer :: lastrecordnum
  real :: height_eps_z

  real :: z_bl ! height of boundary layer

  ! type of initial balance: hydrostatic or geostrophic+hydrostatic
  character(len = 20) :: init_bal

  logical :: output_theta_bgr, output_rho_bgr ! output backgr. env. state
  logical :: output_br_vais_sq ! output N^2 of the background state
  logical :: output_heat ! output heating
  character(len = 20) :: balance_eq ! equations used for initial balance

  namelist / baroclinic_ID / bar_sigma_y, u_strength, dTh_atm, add_noise, &
      proc_noise, init_2Dto3D, tau_relax, tau_jet, RelaxHeating, init_bal, &
      Sponge_Rel_Type, Sponge_Rel_Bal_Type, lastrecordnum, fileinitstate2D, &
      height_eps_z, output_theta_bgr, output_br_vais_sq, output_heat, z_bl, &
      balance_eq, output_rho_bgr

  !-----------------------------------------------------------------
  !                          Model equations
  !-----------------------------------------------------------------

  ! vertical direction
  real, dimension(3) :: vertical
  real :: vert_theta, vert_alpha

  character(len = 25) :: model

  namelist / modelList / model, vert_theta, vert_alpha

  !-----------------------------------------------------------------
  !                               Solver
  !-----------------------------------------------------------------

  real :: cfl
  real :: cfl_wave
  real :: dtMax_dim
  real :: turb_dts
  integer :: n_shap ! (half) order of the Shapiro filter
  !UAC
  !real :: shap_dts_dim                     ! horizontal-Shapiro-filter
  !                                         ! damping time scale (s)
  !                                         ! < 0 means no filter
  real :: shap_dts_fac ! horizontal-Shapiro-filter
  ! damping time scale
  ! (in units of the time step)
  ! < 0 means no filter
  !UAE

  character(len = 20) :: tStepChoice ! "cfl", "fix"
  character(len = 20) :: timeScheme ! LS_Will_RK3 / Euler
  character(len = 20) :: timeSchemeType ! lowStorage / classical
  character(len = 20) :: fluxType ! ILES / central / upwind
  character(len = 20) :: reconstType ! ALDM / constant / SALD / MUSCL
  character(len = 20) :: musclType ! muscl1 / muscl2
  character(len = 20) :: limiterType1 ! minmod / ...
  logical :: fluctuationMode ! split rho = rhoStrat + rho'
  logical :: auxil_equ ! auxiliary equation for the
  ! density fluctuations to be
  ! used in the explicit
  ! integration of the
  ! pseudo-incompressible system
  logical :: TurbScheme ! turbulence Scheme
  logical :: DySmaScheme ! Dynamic Smagorinsky Scheme
  logical :: dtWave_on ! gagarina:
  ! time step for gravity wave
  ! resolution

  logical :: heatingONK14 ! pseudo-incompressible dynamics
  ! with heating as in
  ! ONeill & Klein (2014)
  ! GBcorr
  logical :: heating ! heating on/off defined in
  ! subroutine setup (init.f90)
  ! based on heatingONK14,turbScheme.
  ! raytracer

  logical :: dens_relax ! switch for replacement of
  ! relaxational heating by
  ! density relaxation

  namelist / solverList / cfl, cfl_wave, dtMax_dim, tStepChoice, timeScheme, &
      auxil_equ, fluxType, reconstType, musclType, limiterType1, &
      fluctuationMode, TurbScheme, turb_dts, DySmaScheme, dtWave_on, &
      heatingONK14, dens_relax, shap_dts_fac, n_shap
  !UAC & dens_relax, shap_dts_dim, n_shap

  integer :: nStages
  logical :: updateMass ! transport of mass=var(1)  on/off
  logical :: predictMomentum ! transport of momentum=var(2-4) on/off
  logical :: updateTheta ! transport of theta=var(6) on/off
  logical :: updateIce ! transport of ice=var(nVar-2:nVar) on/off
  logical :: updateTracer ! transport of tracer=var(iVarT) on/off

  !-----------------------------------------------------------------
  !                          Poisson solver
  !-----------------------------------------------------------------

  real :: tolPoisson
  ! achatzb
  real :: tolCond, tolref, abs_tol, scaled_atol, alpha_tol, b_norm
  logical :: tolscale
  ! achatze
  integer :: maxIterPoisson
  character(len = 20) :: poissonSolverType
  character(len = 20) :: storageType
  character(len = 10) :: preconditioner
  character(len = 10) :: tolcrit
  real :: dtau
  integer :: maxIterADI
  logical :: initialCleaning
  logical :: pressureScaling
  logical :: useNAG
  logical :: correctMomentum ! false -> momentumCorrector off
  logical :: correctDivError ! true -> subtract rho*div(u)
  namelist / poissonSolverList / tolPoisson, abs_tol, tolCond, maxIterPoisson, &
      poissonSolverType, storageType, preconditioner, dtau, maxIterADI, &
      initialCleaning, pressureScaling, useNAG, correctMomentum, &
      correctDivError, tolcrit, tolscale
  integer :: nnz ! number of nonzeros

  ! hypre and bicgstab objects
  ! integer * 8 grid_hypre, stencil_e, stencil_i, A_hp_e, A_hp_i, b_hp_e, &
  !     b_hp_i, x_hp_e, x_hp_i, solver_hp_e, solver_hp_i
  ! real, dimension (:), allocatable :: values_e
  ! real, dimension (:), allocatable :: values_i
  !UAC real, dimension(:,:,:), allocatable :: ac_b, al_b,ar_b, ab_b,af_b, &
  real, dimension(:, :, :), allocatable :: ac_b, acv_b, ach_b, al_b, ar_b, &
      ab_b, af_b, ad_b, au_b, alb_b, alf_b, arb_b, arf_b
  !UAE

  ! TFC FJ
  real, dimension(:, :, :), allocatable :: aru_b, ard_b, alu_b, ald_b, afu_b, &
      afd_b, abu_b, abd_b, auu_b, add_b, aruu_b, ardd_b, aluu_b, aldd_b, &
      afuu_b, afdd_b, abuu_b, abdd_b

  ! real, dimension (:), allocatable :: bvalue_vector_hypre, xvalue_vector_hypre

  !achatzb
  !integer, dimension(:), allocatable :: stencil_indices_hypre
  ! integer, dimension (:), allocatable :: stencil_indices_e
  ! integer, dimension (:), allocatable :: stencil_indices_i
  !achatze

  !achatzb
  !integer, parameter :: nentries_hypre = 7
  ! integer, parameter :: ne_hypre_e = 7
  ! integer, parameter :: ne_hypre_i = 11
  !achatze

  !-----------------------------------------------------------------
  !                           Constants
  !-----------------------------------------------------------------

  real :: pi ! you know...
  real, parameter :: small = 1.0e-20 ! to devision by zero
  complex, parameter :: imag = (0.0, 1.0) ! imaginary unit

  !-----------------------------------------------------------------
  !                         Flux specification
  !-----------------------------------------------------------------

  ! ILES (ALDM) parameter
  real :: sigmaC
  real :: sigma0
  real :: sigmaX, sigmaY, sigmaZ

  !-----------------------------------------------------------------
  !                             Atmosphere
  !-----------------------------------------------------------------

  character(len = 10) :: referenceQuantities ! set of reference quantities
  integer, dimension(3) :: bvarOut ! 1 = output, 0 = no output
  ! for the hydrostatically
  ! balanced background state
  logical :: specifyReynolds ! choose Reynolds or mu_viscous
  real :: ReInv ! reciprocal Reynolds number
  real :: mu_viscous_dim ! kinematic viscosity
  real :: mu_conduct_dim ! heat conductivity
  character(len = 30) :: background ! isentropic / isothermal /
  ! const-N / diflapse / HeldSuarez
  real :: N_BruntVaisala_dim ! Brunt-Vaisala frequency in 1/s
  real :: theta0_dim ! isentr. backgr. pot. temp. in K
  real :: Temp0_dim ! isoth. backgr. temp. in K
  real :: press0_dim ! pressure at z=0 in Pa
  real, dimension(3) :: backgroundFlow_dim
  real :: f_Coriolis_dim ! Coriolis parameter
  character(len = 30) :: corset ! constant / periodic
  real :: z_tr_dim ! height of tropopause
  real :: theta_tr_dim ! const pot. temp. of troposphere
  real :: gamma_t, gamma_s ! lapse rates in trop and strat
  !UAB
  real :: tp_strato_dim ! stratosphere temperature (K)
  ! (used by HeldSuarez)
  real :: tp_srf_trp_dim ! tropical surface temperature (K)
  ! (used by HeldSuarez)
  real :: tpdiffhor_tropo_dim ! tropospheric temperature
  ! difference between poles and
  ! tropics (K)
  ! (used by HeldSuarez)
  real :: ptdiffvert_tropo_dim ! vertical potential-temperature
  ! difference in troposphere (K)
  ! (used by HeldSuarez)
  !UAE

  namelist / atmosphereList / referenceQuantities, specifyReynolds, ReInv, &
      mu_viscous_dim, mu_conduct_dim, background, N_BruntVaisala_dim, &
      theta0_dim, Temp0_dim, press0_dim, backgroundFlow_dim, f_Coriolis_dim, &
      corset, z_tr_dim, theta_tr_dim, gamma_t, gamma_s, bvarOut, &
      tp_strato_dim, tp_srf_trp_dim, tpdiffhor_tropo_dim, ptdiffvert_tropo_dim
  !UAB
  !UAE

  real, dimension(3) :: backgroundFlow
  real :: theta00, rho00, P00 ! background values for Boussinesq

  real :: mu_conduct

  !-----------------------------------------------------------------
  !                         Bottom topography
  !-----------------------------------------------------------------

  logical :: topography ! via k = 1

  !UAD logical, dimension(:,:,:), allocatable :: topography_mask

  ! topography_surface x-y-dependent mountain surface
  real, dimension(:, :), allocatable :: topography_surface

  ! FJFeb2023
  real, dimension(:, :), allocatable :: final_topography_surface

  !UAB
  ! vertical index of (velocity) reconstruction points just above the
  ! mountain surface
  integer, dimension(:, :, :), allocatable :: kbl_topo
  ! topography gradients in x- and y-direction
  ! (below the reconstruction points)
  real, dimension(:, :, :), allocatable :: dhdx, dhdy
  ! location of (velocity) interpolation point to be used in determining
  ! the velocity at the reconstruction point
  real, dimension(:, :, :), allocatable :: x_ip, y_ip, z_ip
  ! factors between tangential or normal velocity at the reconstruction
  ! and interpolation points
  real, dimension(:, :, :), allocatable :: velocity_reconst_t
  real, dimension(:, :, :), allocatable :: velocity_reconst_n
  !roughness length
  real :: z_0_dim, z_0
  !UAE

  ! TFC FJ
  integer :: ipolTFC
  logical :: freeSlipTFC
  logical :: testTFC

  ! FJFeb2023
  real :: topographyTime

  real :: mountainHeight_dim
  real :: mountainWidth_dim
  integer :: mountain_case
  integer :: range_factor

  ! TFC FJ
  namelist / topographyList / topography, ipolTFC, freeSlipTFC, testTFC, &
      topographyTime, mountainHeight_dim, mountainWidth_dim, mountain_case, &
      range_factor, z_0_dim
  !UAB
  !UAE

  !-----------------------------------------------------------------
  !                         Boundary
  !-----------------------------------------------------------------

  logical :: rhoFluxCorr, iceFluxCorr, uFluxCorr, vFluxCorr, wFluxCorr, &
      thetaFluxCorr
  integer :: nbCellCorr

  ! sponge layer
  logical :: spongeLayer, sponge_uv
  real :: spongeHeight
  integer :: kSponge
  real :: zSponge
  real :: spongeAlphaZ_dim, spongeAlphaZ_fac

  ! Unified and lateral sponge layers (FJJun2023)
  logical :: unifiedSponge, lateralSponge
  real, dimension(:, :, :), allocatable :: alphaUnifiedSponge
  real :: xSponge0, ySponge0, xSponge1, ySponge1
  real :: dxSponge, dySponge, dzSponge

  ! Vertical sponge layer
  character(len = 50) :: spongeType

  ! Order of polynomial sponge
  integer :: spongeOrder

  ! Damping time of COSMO sponge (in time steps)
  integer :: cosmoSteps

  ! gaga: backup, delete later
  character(len = 10) :: utopcond ! dudz=0 or default u=0
  character(len = 10) :: rhocond ! rho_prime=0 or default drhodz=0
  character(len = 10) :: thcond ! theta_prime=0 or default dthdz=0
  real, dimension(:, :), allocatable :: u_const ! constant wind for baroclinic life cycle experiments
  ! gaga
  ! boundary types
  !  real :: xBoundary
  !  character(len=15) :: xxBoundary
  !  character(len=15) :: yBoundary
  !  character(len=15) :: zBoundary

  namelist / boundaryList / rhoFluxCorr, iceFluxCorr, uFluxCorr, vFluxCorr, &
      wFluxCorr, thetaFluxCorr, nbCellCorr, spongeLayer, sponge_uv, &
      spongeHeight, spongeAlphaZ_dim, spongeAlphaZ_fac, unifiedSponge, &
      lateralSponge, spongeType, spongeOrder, cosmoSteps, utopcond, rhocond, &
      thcond
  !UAC & spongeLayer, spongeHeight, &
  !UAC & spongeAlphaZ_dim, utopcond, rhocond, thcond

  ! boundary types
  character(len = 15) :: xBoundary
  character(len = 15) :: yBoundary
  character(len = 15) :: zBoundary

  namelist / boundaryList2 / xBoundary, yBoundary, zBoundary

  !-----------------------------------------------------------------
  !                     WKB arrays and variables
  !-----------------------------------------------------------------

  type rayType
    real :: x ! ray position
    real :: y
    real :: z

    real :: k ! ray wave vector
    real :: l
    real :: m

    real :: omega ! intrinsic frequency

    real :: area_xk ! x-k area associated with each ray volume
    real :: area_yl ! y-l area associated with each ray volume
    real :: area_zm ! z-m area associated with each ray volume

    real :: dkray ! dm associated with each ray volume
    real :: dlray ! dm associated with each ray volume
    real :: dmray ! dm associated with each ray volume

    real :: dxray ! dz associated with each ray volume
    real :: dyray ! dz associated with each ray volume
    real :: dzray ! dz associated with each ray volume

    real :: dens ! wave-action density

    !SD
    ! initial phase ray volume
    real :: dphi ! initial phase

  end type rayType

  ! total number of active rays
  integer, dimension(:, :, :), allocatable :: nRay

  ! maximum total number of rays
  integer :: nray_max

  ! maximum total number of rays provided in work space
  integer :: nray_wrk

  ! namelist "rayTracer"
  logical :: rayTracer ! run ray tracer
  integer :: nRayRatioX, nRayRatioY, nRayRatioZ

  ! maximum group velocities
  real :: cgx_max, cgy_max, cgz_max

  namelist / wkbList / rayTracer, nRayRatioX, nRayRatioY, nRayRatioZ

  !----------------------------------------------------------------
  !                           Tracer
  !----------------------------------------------------------------
  character(len = 20) :: tracerSetup
  logical :: include_gw_tracer_forcing, include_tracer_mixing, &
      include_env_tracer_forcing, tracerdifference, include_prime

  namelist / tracerList / tracerSetup, include_gw_tracer_forcing, &
      include_tracer_mixing, tracerdifference, include_prime, &
      include_env_tracer_forcing

	type waveAmp
		! components of the vectors given in (10.358)-(10.360) in
		! Achatz, Atmospheric Dynamics (2022)
		! (10.358): nowamp (next-order wave amplitudes)
		! (10.359): rhsamp (right-hand side of next-order wave amp. equations)
		! (10.360): lowamp (leading-order wave amplitudes)

		complex :: u  ! leading-order zonal wind amplitude uhat^(0) in lowamp
								  ! next-order zonal wind amp. uhat^(1) in nowamp
								  ! rhs of equation for uhat^(1) given in (10.352)
		complex :: v  ! wave as for u, but meridional wind
		complex :: w  ! vertical wind component
								  ! rhs equation given in (10.350)
		complex :: b	! buoyancy component as given in the book, so
									! divided by N !!
									! rhs equation given in (10.337)
		complex :: pi	! Exner-pressure component as given in book,
									! so multiplied by c_p/R * thetaStrat
									! rhs equation given in (10.327)
	end type waveAmp
  !-----------------------------------------------------------------
  !                           Ice physics
  !-----------------------------------------------------------------

  character(len = 20) :: iceTestcase ! choose initial ice variable setup
  ! possible: "homogeneous_qv", "homogeneous_SIce"
  real :: init_SIce ! initial vapor saturation with respect to ice
  real :: init_nAer, init_qv ! initial values for aerosols and humidity
  real :: init_m_ice ! initially assumed average ice crystal mass
  real :: radius_solution, sigma_r ! needed for log-normal correlation
  real :: T_nuc ! initial nucleation temperature for 1D_ISSR simple flow case
  real :: p_nuc ! initial nucleation pressure for 1D_ISSR simple flow case
  ! this overwrites press0_dim
  real :: dt_ice ! length of the microphysical time step
  character(len = 10) :: NUC_approx_type ! nucleation approximation type
  ! possible: "Koop", "linFit", "threshold"
  real :: ISSR_top ! top of ISSR relative to mountain height in qv_relaxation testcase
  logical :: super_simplified ! use a super_simplified DEP and SED scheme
  logical :: kT_linFit, dv_exp2, cm_dryAir, mu_linFit ! switch approximations on/off
  logical :: sedimentation_on ! turn sedimentation terms on or off
  logical :: nucleation_on ! turn nucleation on or off
  logical :: evaporation_on ! turn evaporation on or off
  character(len = 10) :: awi_type ! possible: "const", "linFit", "quadFit", "exact"
  character(len = 10) :: SIce_threshold_type ! possible: "linFit", "quadFit", "exact"

  namelist / iceList / iceTestcase, init_SIce, init_nAer, init_qv, init_m_ice, &
      radius_solution, sigma_r, T_nuc, p_nuc, dt_ice, NUC_approx_type, &
      ISSR_top, super_simplified, kT_linFit, dv_exp2, cm_dryAir, mu_linFit, &
      sedimentation_on, nucleation_on, evaporation_on, awi_type, &
      SIce_threshold_type

  !SD: Ice physics 2
  type opt_rayType
    real :: Fsat !full saturation ration within ray volume
    real :: Msat !mean saturation ration within ray volume
    real :: w ! vertical velocity amplitude
  end type opt_rayType

  integer, dimension(:), allocatable :: iVarIce
  type(opt_rayType), dimension(:, :, :, :), allocatable :: opt_ray
  integer :: inN, inQ, inQv, nVarIce

  real :: J_nuc !Nucleation rate [#/kg/s]
  real :: B_nuc !Exponent nucleation function (dimension less)
  real, parameter :: S_c = 1.5 !Critical saturation ratio

  !  real, parameter :: Mole_mass_water = 18.01528e-3, Mole_mass_dryAir = 28.9644e-3
  real, parameter :: epsil0 = 0.62 !\approx Mole_mass_water/Mole_mass_dryAir
  real, parameter :: meanMassIce = 1.E-12 ! mean mass ice crystals [kg]
  real :: mRef ! reference mass
  real, parameter :: PsatIceRef = 1 !reference saturation pressure [Pa]
  real :: Dep ! deposition coefficient
  real, dimension(:, :, :), allocatable :: PsatIce
  real, parameter :: thetaRef_trp = 210. ! reference temperature in the tropopause region [K]
  real, parameter :: L_ice = 2.8E6 ! constant latent heat  ice [J/kg]
  real, parameter :: R_v = 461. ! specific gas constant for water vapor [J/kg/K]

  real :: dt_ice2 = 0.1 ! length of the microphysical time step [s]

  real :: thetaRefRatio !thetaRef/thetaRef_tropo_pause
  real :: L_hat, Li_hat, epsil0hat
  real :: alr_ndim ! non-dimensioanl adiabatic lapse rate
  real * 4, dimension(:, :, :, :), allocatable :: ofield

  logical :: no_ice_source ! set to true if only ice advection and no ice source
  ! are considered
  logical, parameter :: compare_raytracer = .True. ! modify Wavepacket simulation to facilitate comparison with raytracer
  namelist / iceList2 / inN, inQ, inQv, nVarIce, dt_ice2, no_ice_source

  real, allocatable :: var_ww(:, :, :) ! flux <w'w'> from RayTracer

  !include_testoutput
  integer :: iVarO, iVarS, iVarP, iVarAW

  !output ray volumes
  real * 4, dimension(:), allocatable :: nor_mst
  real * 4, dimension(:, :), allocatable :: vct_out, vct_prc, vct_mst

  integer, parameter :: NFR = 15 ! number of quantities stored for each ray volume
  ! \vec x, \vec k, \vec dx, \vec dk, omega, density,
  ! vertical vel.
  integer :: TTNOR ! total number of rays
  integer :: NoR_out ! counter for ray volumes output

  contains

  subroutine default_values

    !-----------------------------------------------------------------
    !                     Set default values
    !-----------------------------------------------------------------

    ! FJ: Parameters marked as deprecated do not appear in the code or comments
    ! and can be removed!

    ! Variables
    nVar = 8
    nOptVar = 4 ! deprecated
    include_ice = .false.
    include_ice2 = .false.
    include_tracer = .false.
    include_testoutput = .false.

    ! Input/Output
    dataFileName = "" ! deprecated
    restartFile = "" ! deprecated
    dimOut = [.false., .false., .false.]
    ! varOut = [1, 1, 1, 1, 1, 1, 1, 0] ! must be allocated
    ! varIn = [1, 1, 1, 1, 1, 1, 1, 0] ! must be allocated
    iIn = 1
    ! offset = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] ! must be allocated
    optVarOut = [0, 0, 0, 0, 0, 0]
    wkbVarOut = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    outputType = "time"
    nOutput = 1
    maxIter = 1
    outputTimeDiff = 3600.0
    maxTime = 3600.0
    solutionTime = .false. ! deprecated
    solutionTimeUnit = "s" ! deprecated
    restart = .false.
    showGhostCellsX = .false. ! deprecated
    showGhostCellsY = .false. ! deprecated
    showGhostCellsZ = .false. ! deprecated
    thetaOffset = .false.
    rhoOffset = .false.
    detailedinfo = .false.
    RHS_diagnostics = .false.
    PVinversion = .false. ! deprecated

    ! Parameter study
    parameterStudy = .false. ! deprecated
    startParam = 0; endParam = 0; stepParam = 0 ! deprecated
    paramName = "" ! deprecated

    ! Debugging
    verbose = .false.
    dtMin_dim = 1.0e-6

    ! Test cases
    testCase = "IGW"

    ! Monochromatic wave
    lambda_dim = 0.1 * (lz_dim(1) - lz_dim(0))

    ! Wave packet
    wavePacketType = 1
    wavePacketDim = 1
    lambdaX_dim = 0.1 * (lx_dim(1) - lx_dim(0))
    lambdaY_dim = 0.1 * (ly_dim(1) - ly_dim(0))
    lambdaZ_dim = 0.1 * (lz_dim(1) - lz_dim(0))
    meanFlowX_dim = 0.0 ! deprecated
    meanFlowZ_dim = 0.0 ! deprecated
    amplitudeFactor = 1.0
    xCenter_dim = 0.5 * (lx_dim(0) + lx_dim(1))
    yCenter_dim = 0.5 * (ly_dim(0) + ly_dim(1))
    zCenter_dim = 0.5 * (lz_dim(0) + lz_dim(1))
    sigma_hor_dim = 0.5 * (lx_dim(1) - lx_dim(0))
    sigma_hor_yyy_dim = 0.5 * (ly_dim(1) - ly_dim(0))
    sigma_dim = 0.5 * (lz_dim(1) - lz_dim(0))
    amp_mod_x = 1.0
    amp_mod_y = 1.0
    L_cos_dim = 0.5 * (lz_dim(1) - lz_dim(0))
    omiSign = 1
    u0_jet_dim = 0.0
    z0_jet_dim = 0.0
    L_jet_dim = 0.0
    inducedwind = .false.

    ! Lagrangian ray tracer
    xrmin_dim = lx_dim(0); xrmax_dim = lx_dim(1)
    yrmin_dim = ly_dim(0); yrmax_dim = ly_dim(1)
    zrmin_dim = lz_dim(0); zrmax_dim = lz_dim(1)
    nrxl = 1; nryl = 1; nrzl = 1
    fac_dk_init = 0.1; fac_dl_init = 0.1; fac_dm_init = 0.1
    nrk_init = 1; nrl_init = 1; nrm_init = 1
    nsmth_wkb = 2
    lsmth_wkb = .true.
    sm_filter = 2
    lsaturation = .true.
    alpha_sat = 1.0
    steady_state = .false.
    case_wkb = 3
    amp_wkb = 1.0
    wlrx_init = 0.1 * (lx_dim(1) - lx_dim(0))
    wlry_init = 0.1 * (ly_dim(1) - ly_dim(0))
    wlrz_init = 0.1 * (lz_dim(1) - lz_dim(0))
    xr0_dim = 0.5 * (lx_dim(0) + lx_dim(1))
    yr0_dim = 0.5 * (ly_dim(0) + ly_dim(1))
    zr0_dim = 0.5 * (lz_dim(0) + lz_dim(1))
    sigwpx_dim = 0.5 * (lx_dim(1) - lx_dim(0))
    sigwpy_dim = 0.5 * (ly_dim(1) - ly_dim(0))
    sigwpz_dim = 0.5 * (lz_dim(1) - lz_dim(0))
    branchr = 1
    lindUinit = .false. ! deprecated
    topographyTime_wkb = 0.0
    mountainHeight_wkb_dim = 0.1 * (lz_dim(1) - lz_dim(0))
    mountainWidth_wkb_dim = 0.1 * (lx_dim(1) - lx_dim(0))
    mountain_case_wkb = 1
    range_factor_wkb = 1
    blocking = .false.
    nwm = 1
    zmin_wkb_dim = 0.0
    nray_fac = 20
    cons_merge = "en"

    ! Bubbles
    dTheta0_dim = 1.0
    xRadius_dim = 0.5 * (lx_dim(1) - lx_dim(0))
    zRadius_dim = 0.5 * (lz_dim(1) - lz_dim(0))
    xCenter_dim = 0.5 * (lx_dim(0) + lx_dim(1))
    yCenter_dim = 0.5 * (ly_dim(0) + ly_dim(1))
    zCenter_dim = 0.5 * (lz_dim(0) + lz_dim(1))
    zExcentricity = 1.0
    rhoCenter_dim = 1000.0

    ! Robert bubble
    dTheta1_dim = 1.0
    a1_dim = 0.25 * min(lx_dim(1) - lx_dim(0), lz_dim(1) - lz_dim(0))
    sigma1_dim = 0.25 * min(lx_dim(1) - lx_dim(0), lz_dim(1) - lz_dim(0))
    xCenter1_dim = 0.25 * (lx_dim(0) + lx_dim(1))
    zCenter1_dim = 0.25 * (lz_dim(0) + lz_dim(1))
    dTheta2_dim = - 1.0
    a2_dim = 0.25 * min(lx_dim(1) - lx_dim(0), lz_dim(1) - lz_dim(0))
    sigma2_dim = 0.25 * min(lx_dim(1) - lx_dim(0), lz_dim(1) - lz_dim(0))
    xCenter2_dim = 0.75 * (lx_dim(0) + lx_dim(1))
    zCenter2_dim = 0.75 * (lz_dim(0) + lz_dim(1))

    ! Mountain wave
    u_relax = 10.0
    v_relax = 0.0
    w_relax = 0.0
    t_relax = 0.0
    t_ramp = 0.0
    xextent_norelax = 0.0
    wind_relaxation = .false.
    surface_layer_depth = 0.0

    ! Baroclinic life cycle (realistic)
    zero_initial_state = .false.
    z_trpp0_dim = 10000.0
    z_baro_dim = 20000.0
    deltht_dim = 30.0 ! deprecated
    thet0_dim = 300.0
    ntrp_dim = 0.01
    nstr_dim = 0.02
    jwdth_dim = 0.0
    kaptpp = 0.0
    t_relax_bar = 0.0 ! deprecated
    add_ptptb = .false.
    ptptb_x_dim = 0.5 * (lx_dim(0) + lx_dim(1))
    ptptb_y_dim = 0.5 * (ly_dim(0) + ly_dim(1))
    ptptb_z_dim = 0.5 * (lz_dim(0) + lz_dim(1))
    ptptb_dh_dim = 0.1 * min(lx_dim(1) - lx_dim(0), ly_dim(1) - ly_dim(0))
    ptptb_dz_dim = 0.1 * (lz_dim(1) - lz_dim(0))
    ptptb_amp_dim = 1.0
    dTh_atm = 30.0
    add_noise = .false.
    proc_noise = 0.1
    tau_relax = 0.0
    tau_relax_low = 0.0
    sigma_tau = 0.0
    tau_jet = 0.0
    RelaxHeating = 0
    Sponge_Rel_Type = "constant"
    Sponge_Rel_Bal_Type = "env"
    init_2Dto3D = .false.
    fileinitstate2D = ""
    ta_hs_dim = 0.0
    ts_hs_dim = 0.0
    tf_hs_dim = 0.0
    sigb_hs = 0.0

    ! Baroclinic life cycle (idealized)
    bar_sigma_y = 0.5 * (ly_dim(1) - ly_dim(0))
    u_strength = 30.0
    dTh_atm = 30.0
    add_noise = .false.
    proc_noise = 0.1
    init_2Dto3D = .false.
    tau_relax = 0.0
    tau_jet = 0.0
    RelaxHeating = 0
    init_bal = "geostr_id"
    Sponge_Rel_Type = "constant"
    Sponge_Rel_Bal_Type = "env"
    lastrecordnum = 0
    fileinitstate2D = ""
    height_eps_z = 0.0 ! deprecated
    output_theta_bgr = .false.
    output_br_vais_sq = .false.
    output_heat = .false.
    z_bl = 0.0 ! deprecated
    balance_eq = "PI"
    output_rho_bgr = .false.

    ! Model equations
    model = "pseudo_incompressible"
    vert_theta = 90.0
    vert_alpha = 0.0

    ! Solver
    cfl = 0.5
    cfl_wave = 0.5
    dtMax_dim = 1.0e3
    tStepChoice = "cfl"
    timeScheme = "LS_Will_RK3"
    auxil_equ = .false.
    fluxType = "upwind"
    reconstType = "MUSCL"
    musclType = "muscl1"
    limiterType1 = "MCVariant"
    fluctuationMode = .true.
    TurbScheme = .false.
    turb_dts = 5.0e3
    DySmaScheme = .false.
    dtWave_on = .true.
    heatingONK14 = .false.
    dens_relax = .false.
    shap_dts_fac = - 1.0
    n_shap = 1

    ! Poisson solver
    tolPoisson = 1.0e-8
    abs_tol = 0.0
    tolCond = 1.0e-23
    maxIterPoisson = 1000
    poissonSolverType = "bicgstab"
    storageType = "opr" ! deprecated
    preconditioner = "yes"
    dtau = 4.0e-4
    maxIterADI = 2
    initialCleaning = .true.
    pressureScaling = .false.
    useNAG = .false. ! deprecated
    correctMomentum = .true.
    correctDivError = .false.
    tolcrit = "abs"
    tolscale = .true. ! deprecated

    ! Atmosphere
    referenceQuantities = "Klein"
    specifyReynolds = .false.
    ReInv = 0.0
    mu_viscous_dim = 0.0
    mu_conduct_dim = 0.0
    background = "isothermal"
    N_BruntVaisala_dim = 0.01
    theta0_dim = 300.0
    Temp0_dim = 300.0
    press0_dim = 100000.0
    backgroundFlow_dim = [10.0, 0.0, 0.0]
    f_Coriolis_dim = 0.0
    corset = "constant"
    z_tr_dim = 15000.0
    theta_tr_dim = 300.0
    gamma_t = 0.0
    gamma_s = 0.0
    bvarOut = 0 ! deprecated
    tp_strato_dim = 200.0
    tp_srf_trp_dim = 300.0
    tpdiffhor_tropo_dim = 50.0
    ptdiffvert_tropo_dim = 10.0

    ! Topography
    topography = .false.
    ipolTFC = 2
    freeSlipTFC = .false.
    testTFC = .false.
    topographyTime = 0.0
    mountainHeight_dim = 0.1 * (lz_dim(1) - lz_dim(0))
    mountainWidth_dim = 0.1 * (lx_dim(1) - lx_dim(0))
    mountain_case = 1
    range_factor = 1
    z_0_dim = 0.0 ! deprecated

    ! Boundaries
    rhoFluxCorr = .false.
    iceFluxCorr = .false.
    uFluxCorr = .false.; vFluxCorr = .false.; wFluxCorr = .false.
    thetaFluxCorr = .false.
    nbCellCorr = 0
    spongeLayer = .false.
    sponge_uv = .false.
    spongeHeight = 0.33
    spongeAlphaZ_dim = 0.01
    spongeAlphaZ_fac = 0.01
    unifiedSponge = .false.
    lateralSponge = .false.
    spongeType = "polynomial"
    spongeOrder = 1
    cosmoSteps = 1
    utopcond = "" ! deprecated
    rhocond = "" ! deprecated
    thcond = "" ! deprecated
    xBoundary = "periodic"
    yBoundary = "periodic"
    zBoundary = "solid_wall"

    ! WKB
    rayTracer = .false.
    nRayRatioX = 1; nRayRatioY = 1; nRayRatioZ = 1 ! deprecated

    ! Tracer
    tracerSetup = "thin_layer"
    include_gw_tracer_forcing = .true.
    include_tracer_mixing = .true.
    tracerdifference = .true.
    include_prime = .true.
    include_env_tracer_forcing = .true.

    ! Ice
    iceTestcase = "ice_cloud"
    init_SIce = 1.0e5
    init_nAer = 1.0e8
    init_qv = 1.0e-4
    init_m_ice = 1.0e-16
    radius_solution = 1.0e-8
    sigma_r = 1.0
    T_nuc = 200.0
    p_nuc = 2.0e4
    dt_ice = 1.0e-3
    NUC_approx_type = "linFit"
    ISSR_top = 0.0
    super_simplified = .false.
    kT_linFit = .false.
    dv_exp2 = .false.
    cm_dryAir = .false.
    mu_linFit = .false.
    sedimentation_on = .false.
    nucleation_on = .false.
    evaporation_on = .false.
    awi_type = "const"
    SIce_threshold_type = "linFit"
    inN = 0
    inQ = 0
    inQv = 0
    nVarIce = 0
    dt_ice2 = 0.0
    no_ice_source = .false.

  end subroutine default_values

end module type_module
