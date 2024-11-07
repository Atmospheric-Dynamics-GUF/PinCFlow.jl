module type_module

  use mpi
  !-----------------------------------------------------------------
  !    Definition of data types and variables accessible
  !    throughout the code! Use care when working on these
  !    data fields. When adding a new namelist parameter,
  !    set its default value in the subroutine default_values
  !    and add an appropriate line in the subroutine
  !    write_namelists.
  !-----------------------------------------------------------------

  implicit none

  public

  !-----------------------------------------------------------------
  !                     MPI & Domain composition
  !-----------------------------------------------------------------

  integer :: sizeX, sizeY, sizeZ
  integer :: nbx, nby, nbz
  integer :: nprocx, nprocy
  real, dimension(0:1) :: lx_dim, ly_dim, lz_dim ! dimensional domain

  namelist / domain / sizeX, sizeY, sizeZ, nbx, nby, nbz, lx_dim, ly_dim, &
      &lz_dim, nprocx, nprocy

  integer :: sizeXX, sizeYY, sizeZZ
  integer :: nx1, in1, ny1, jn1
  integer :: iStart, jStart
  integer :: is, ie, js, je ! local start and end indices
  logical :: verboseMPI

  ! MPI include (parameters needed below)
  !include 'mpif.h'

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
  logical :: include_ice = .false., include_testoutput = .false.
  logical :: include_tracer = .false.
  namelist / variables / include_ice, include_tracer, include_testoutput

  type var_type
    real, dimension(:, :, :), allocatable :: rho ! Density
    real, dimension(:, :, :), allocatable :: u ! Zonal wind
    real, dimension(:, :, :), allocatable :: v ! Meridional wind
    real, dimension(:, :, :), allocatable :: w ! Vertical wind
    real, dimension(:, :, :), allocatable :: pi ! Exner pressure
    real, dimension(:, :, :), allocatable :: rhop ! Density fluctuations
    real, dimension(:, :, :), allocatable :: DSC ! Dyn. Smagorinsky coeff.
    real, dimension(:, :, :), allocatable :: GWH ! GW heating
    real, dimension(:, :, :), allocatable :: P ! Mass-weighted pot. temp.
    real, dimension(:, :, :), allocatable :: chi ! Tracer
    real, dimension(:, :, :, :), allocatable :: ICE ! Alternative ice variables
    real, dimension(:, :, :, :), allocatable :: OPT ! Optional variables
  end type var_type

  type flux_type
    real, dimension(:, :, :, :), allocatable :: rho ! Density fluxes
    real, dimension(:, :, :, :), allocatable :: u ! Zonal-wind fluxes
    real, dimension(:, :, :, :), allocatable :: v ! Meridional-wind fluxes
    real, dimension(:, :, :, :), allocatable :: w ! Vertical-wind fluxes
    real, dimension(:, :, :, :), allocatable :: theta ! Pot.-temp. fluxes
    real, dimension(:, :, :, :), allocatable :: rhop ! Dens-fluct. fluxes
    real, dimension(:, :, :, :), allocatable :: P ! M.-w.-pot.-temp. fluxes
    real, dimension(:, :, :, :), allocatable :: chi ! Tracer fluxes
    real, dimension(:, :, :, :, :), allocatable :: ICE ! Alternative ice fluxes
  end type flux_type

  !-----------------------------------------------------------------
  !                    Input / Output variables
  !-----------------------------------------------------------------
  integer :: iOut ! output counter ; gagarina: moved from pinc

  character(len = 256) :: file_namelist
  integer, parameter :: nVar = 20 ! Maximum number of I/O variables

  character(len = 10), dimension(1:20) :: atmvarOut
  character(len = 10), dimension(1:20) :: rayvarOut
  character(len = 10), dimension(1:20) :: icevarOut

  logical :: saverayvols
  logical :: prepare_restart
  logical :: restart

  integer :: iIn

  character(len = 40) :: runName

  character(len = 20) :: outputType ! "time" or "timeStep"
  integer :: nOutput ! output every nOutput's time step
  integer :: maxIter ! max nb. of time steps
  real :: outputTimeDiff ! output every ... seconds
  real :: maxTime ! max time in seconds

  logical :: detailedinfo
  logical :: RHS_diagnostics

  logical :: fancy_namelists

  namelist / outputList / atmvarOut, rayvarOut, icevarOut, saverayvols, &
      &prepare_restart, restart, iIn, runName, outputType, nOutput, maxIter, &
      &outputTimeDiff, maxTime, detailedinfo, RHS_diagnostics, fancy_namelists
  !achatzb
  !achatze

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
  real :: amplitudeFactor
  real :: sigma_dim
  real :: sigma_hor_dim, sigma_hor_yyy_dim ! modified by Junhong Wei (20170214)
  real :: L_cos_dim
  integer :: omiSign
  real :: u0_jet_dim
  real :: z0_jet_dim
  real :: L_jet_dim
  real :: x0_dim, y0_dim, z0_dim
  ! achatzb
  real :: amp_mod_x, amp_mod_y
  logical :: inducedwind
  ! achatze

  namelist / wavePacket / wavePacketType, wavePacketDim, lambdaX_dim, &
      &lambdaY_dim, lambdaZ_dim, amplitudeFactor, x0_dim, y0_dim, z0_dim, &
      &sigma_dim, sigma_hor_dim, amp_mod_x, sigma_hor_yyy_dim, amp_mod_y, &
      &L_cos_dim, omiSign, u0_jet_dim, z0_jet_dim, L_jet_dim, inducedwind
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

  ! Long number scaling of displacement (FJJan2023)
  logical :: blocking

  ! Number of wave modes (FJApr2023)
  integer :: nwm

  character(len = 10) :: launch_algorithm

  real :: zmin_wkb_dim, zmin_wkb

  integer :: nray_fac ! maximum factor by which the # of rays may increase
  ! compared to initialization

  character(len = 2) :: cons_merge ! quantity to be conserved
  ! ("wa" = wave action/
  !  "en" = wave energy)
  ! under ray-volume merging

  integer :: nRayOutput

  namelist / LagrangeRayTracing / xrmin_dim, xrmax_dim, yrmin_dim, yrmax_dim, &
      &zrmin_dim, zrmax_dim, nrxl, nryl, nrzl, fac_dk_init, fac_dl_init, &
      &fac_dm_init, nrk_init, nrl_init, nrm_init, nsmth_wkb, lsmth_wkb, &
      &sm_filter, lsaturation, alpha_sat, steady_state, case_wkb, amp_wkb, &
      &wlrx_init, wlry_init, wlrz_init, xr0_dim, yr0_dim, zr0_dim, sigwpx_dim, &
      &sigwpy_dim, sigwpz_dim, branchr, lindUinit, blocking, nwm, &
      &launch_algorithm, zmin_wkb_dim, nray_fac, cons_merge, nRayOutput ! JaWi: new nml!
  ! Jan Weinkaemmerer, 27.11.18

  !------------------------------------------
  ! hotBubble, coldBubble, hotBubble3D
  !------------------------------------------

  real :: dTheta0_dim, xRadius_dim, zRadius_dim, rhoCenter_dim
  real :: xCenter_dim, yCenter_dim, zCenter_dim
  real :: zExcentricity
  namelist / bubble / dTheta0_dim, xRadius_dim, zRadius_dim, xCenter_dim, &
      &yCenter_dim, zCenter_dim, zExcentricity, rhoCenter_dim

  ! hot and cold bubble by Robert
  real :: dTheta1_dim, a1_dim, sigma1_dim, xCenter1_dim, zCenter1_dim
  real :: dTheta2_dim, a2_dim, sigma2_dim, xCenter2_dim, zCenter2_dim
  namelist / robert_bubble / dTheta1_dim, a1_dim, sigma1_dim, xCenter1_dim, &
      &zCenter1_dim, dTheta2_dim, a2_dim, sigma2_dim, xCenter2_dim, zCenter2_dim

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
  real :: xextent_relax
  real :: yextent_relax
  ! TFC FJ
  ! Switch for background relaxation.
  logical :: wind_relaxation
  ! FJMar2023
  ! Surface layer depth.
  real :: surface_layer_depth
  namelist / mountainwavelist / u_relax, v_relax, w_relax, t_relax, t_ramp, &
      &xextent_relax, yextent_relax, wind_relaxation, surface_layer_depth
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

  real :: ta_hs_dim ! thermal-relaxation time scale outside tropical
  ! boundary layer (s, used by HeldSuarez)
  real :: ts_hs_dim ! thermal-relaxation time scale in tropical
  ! boundary layer (s, used by HeldSuarez)
  real :: tf_hs_dim ! boundary-layer Rayleigh-damping time scale
  ! (s, used by HeldSuarez)
  real :: sigb_hs ! sigma of boundary-layer top (used by HeldSuarez)

  ! relaxation in massUpdate for the density
  real :: tau_relax, tau_relax_low, sigma_tau
  ! time scale for a jet formation
  real :: tau_jet

  ! environmental state
  type(var_type) :: var_env
  real, dimension(:, :, :), allocatable :: p_env_pp, the_env_pp, dens_env_pp, &
      &u_env_pp, v_env_pp

  ! relaxation rates for Held & Suarez (1994)
  real, dimension(:, :), allocatable :: kt_hs, kv_hs
  real, dimension(:), allocatable :: kw_hs

  ! relaxation rates for Held & Suarez (1994)
  real, dimension(:, :), allocatable :: kr_sp, kr_sp_w

  real, dimension(:, :, :), allocatable :: kt_hs_tfc
  real, dimension(:, :, :), allocatable :: kr_sp_tfc, kr_sp_w_tfc

  character(len = 20) :: fileinitstate2D

  ! type of relaxation in the sponge
  character(len = 20) :: Sponge_Rel_Bal_Type

  ! add noise
  logical :: add_noise, init_2Dto3D

  real :: proc_noise, dTh_atm

  namelist / baroclinic_LC / zero_initial_state, z_trpp0_dim, z_baro_dim, &
      &thet0_dim, ntrp_dim, nstr_dim, jwdth_dim, kaptpp, add_ptptb, &
      &ptptb_x_dim, ptptb_y_dim, ptptb_z_dim, ptptb_dh_dim, ptptb_dz_dim, &
      &ptptb_amp_dim, add_noise, proc_noise, tau_relax, tau_relax_low, &
      &sigma_tau, tau_jet, Sponge_Rel_Bal_Type, ta_hs_dim, ts_hs_dim, &
      &tf_hs_dim, sigb_hs

  !-----------------------------------------------------------------
  !                          Baroclinic life cycle: idealized
  !-----------------------------------------------------------------

  ! define the jet parameters
  real :: bar_sigma_y, u_strength
  integer :: lastrecordnum

  ! type of initial balance: hydrostatic or geostrophic+hydrostatic
  character(len = 20) :: init_bal

  logical :: output_theta_bgr, output_rho_bgr ! output backgr. env. state
  logical :: output_br_vais_sq ! output N^2 of the background state
  logical :: output_heat ! output heating
  character(len = 20) :: balance_eq ! equations used for initial balance

  namelist / baroclinic_ID / bar_sigma_y, u_strength, dTh_atm, init_2Dto3D, &
      &init_bal, lastrecordnum, fileinitstate2D, output_theta_bgr, &
      &output_br_vais_sq, output_heat, balance_eq, output_rho_bgr

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
  !real :: shap_dts_dim                     ! horizontal-Shapiro-filter
  !                                         ! damping time scale (s)
  !                                         ! < 0 means no filter
  real :: shap_dts_fac ! horizontal-Shapiro-filter
  ! damping time scale
  ! (in units of the time step)
  ! < 0 means no filter

  character(len = 20) :: tStepChoice ! "cfl", "fix"
  character(len = 20) :: timeScheme ! LS_Will_RK3 / Euler
  character(len = 20) :: timeSchemeType ! lowStorage / classical
  character(len = 20) :: fluxType ! ILES / central / upwind
  character(len = 20) :: reconstType ! ALDM / constant / SALD / MUSCL
  character(len = 20) :: musclType ! muscl1 / muscl2
  character(len = 20) :: limiterType1 ! minmod / ...
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
      &auxil_equ, fluxType, reconstType, musclType, limiterType1, TurbScheme, &
      &turb_dts, DySmaScheme, dtWave_on, heatingONK14, dens_relax, &
      &shap_dts_fac, n_shap
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
  real :: tolCond, tolref, abs_tol, scaled_atol, alpha_tol, b_norm
  integer :: maxIterPoisson
  character(len = 20) :: poissonSolverType
  character(len = 10) :: preconditioner
  character(len = 10) :: tolcrit
  real :: dtau
  integer :: maxIterADI
  logical :: initialCleaning
  logical :: pressureScaling
  logical :: correctMomentum ! false -> momentumCorrector off
  logical :: correctDivError ! true -> subtract rho*div(u)
  namelist / poissonSolverList / tolPoisson, abs_tol, tolCond, maxIterPoisson, &
      &poissonSolverType, preconditioner, dtau, maxIterADI, initialCleaning, &
      &pressureScaling, correctMomentum, correctDivError, tolcrit
  integer :: nnz ! number of nonzeros

  ! hypre and bicgstab objects
  ! integer * 8 grid_hypre, stencil_e, stencil_i, A_hp_e, A_hp_i, b_hp_e, &
  !     b_hp_i, x_hp_e, x_hp_i, solver_hp_e, solver_hp_i
  ! real, dimension (:), allocatable :: values_e
  ! real, dimension (:), allocatable :: values_i
  !UAC real, dimension(:,:,:), allocatable :: ac_b, al_b,ar_b, ab_b,af_b, &
  real, dimension(:, :, :), allocatable :: ac_b, acv_b, ach_b, al_b, ar_b, &
      &ab_b, af_b, ad_b, au_b, alb_b, alf_b, arb_b, arf_b
  !UAE

  ! TFC FJ
  real, dimension(:, :, :), allocatable :: aru_b, ard_b, alu_b, ald_b, afu_b, &
      &afd_b, abu_b, abd_b, auu_b, add_b, aruu_b, ardd_b, aluu_b, aldd_b, &
      &afuu_b, afdd_b, abuu_b, abdd_b

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

  namelist / atmosphereList / referenceQuantities, specifyReynolds, ReInv, &
      &mu_viscous_dim, mu_conduct_dim, background, N_BruntVaisala_dim, &
      &theta0_dim, Temp0_dim, press0_dim, backgroundFlow_dim, f_Coriolis_dim, &
      &corset, z_tr_dim, theta_tr_dim, gamma_t, gamma_s, tp_strato_dim, &
      &tp_srf_trp_dim, tpdiffhor_tropo_dim, ptdiffvert_tropo_dim

  real, dimension(3) :: backgroundFlow
  real :: theta00, rho00, P00 ! background values for Boussinesq

  real :: mu_conduct

  !-----------------------------------------------------------------
  !                         Bottom topography
  !-----------------------------------------------------------------

  logical :: topography ! via k = 1

  ! Resolved topography
  real, dimension(:, :), allocatable :: topography_surface, &
      &final_topography_surface

  ! Unresolved topography
  real, dimension(:, :, :), allocatable :: k_spectrum, l_spectrum, &
      &topography_spectrum, final_topography_spectrum

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
  real :: z_0

  integer :: ipolTFC
  logical :: freeSlipTFC
  logical :: testTFC

  real :: topographyTime

  real :: mountainHeight_dim
  real :: mountainWidth_dim
  integer :: mountain_case
  real :: range_factor
  integer :: spectral_modes

  namelist / topographyList / topography, ipolTFC, freeSlipTFC, testTFC, &
      &topographyTime, mountainHeight_dim, mountainWidth_dim, mountain_case, &
      &range_factor, spectral_modes
  !UAB
  !UAE

  !-----------------------------------------------------------------
  !                         Boundary
  !-----------------------------------------------------------------

  logical :: rhoFluxCorr, iceFluxCorr, uFluxCorr, vFluxCorr, wFluxCorr, &
      &thetaFluxCorr
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

  logical :: relax_to_mean

  ! gaga: backup, delete later
  real, dimension(:, :), allocatable :: u_const ! constant wind for baroclinic life cycle experiments

  namelist / boundaryList / rhoFluxCorr, iceFluxCorr, uFluxCorr, vFluxCorr, &
      &wFluxCorr, thetaFluxCorr, nbCellCorr, spongeLayer, sponge_uv, &
      &spongeHeight, spongeAlphaZ_dim, spongeAlphaZ_fac, unifiedSponge, &
      &lateralSponge, spongeType, spongeOrder, cosmoSteps, relax_to_mean

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

  ! maximum group velocities
  real :: cgx_max, cgy_max, cgz_max

  namelist / wkbList / rayTracer

  !----------------------------------------------------------------
  !                           Tracer
  !----------------------------------------------------------------
  character(len = 20) :: tracerSetup
  logical :: include_trfrc_lo, include_trfrc_mix, include_trfrc_no

  namelist / tracerList / tracerSetup, include_trfrc_lo, include_trfrc_mix, &
      &include_trfrc_no

  type :: waveAmpCompType
    ! components of the vectors given in (10.358)-(10.360) in
    ! Achatz, Atmospheric Dynamics (2022)
    ! (10.358): nowamp (next-order wave amplitudes)
    ! (10.359): rhsamp (right-hand side of next-order wave amp. equations)
    ! (10.360): lowamp (leading-order wave amplitudes)

    complex :: u ! leading-order zonal wind amplitude uhat^(0) in lowamp
    ! next-order zonal wind amp. uhat^(1) in nowamp
    ! rhs of equation for uhat^(1) given in (10.352)
    complex :: v ! wave as for u, but meridional wind
    complex :: w ! vertical wind component
    ! rhs equation given in (10.350)
    complex :: b ! buoyancy component as given in the book, but
    ! not divided by N !!
    ! rhs equation given in (10.337)
    complex :: pi ! Exner-pressure component as given in book, but
    ! multiplied by c_p/R * theta0 !!
    ! rhs equation given in (10.327)
    complex :: chi ! tracer mixing ratio component

  end type waveAmpCompType
  type :: waveAmpType

    type(waveAmpCompType) :: lowamp, rhsamp, nowamp, lowamp_prevts

    real :: phase

    real :: phaserhs

  end type waveAmpType
  type :: tracerFluxType

    real :: uflx, vflx, wflx, total

  end type tracerFluxType
  type :: tracerForceType

    type(tracerFluxType) :: loforce, noforce, mixingGW

  end type tracerForceType
  !-----------------------------------------------------------------
  !                           Ice physics
  !-----------------------------------------------------------------
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

  real :: dt_ice = 0.1 ! length of the microphysical time step [s]

  real :: thetaRefRatio !thetaRef/thetaRef_tropo_pause
  real :: L_hat, Li_hat, epsil0hat
  real :: alr_ndim ! non-dimensioanl adiabatic lapse rate
  real * 4, dimension(:, :, :, :), allocatable :: ofield

  logical :: no_ice_source ! set to true if only ice advection and no ice source
  ! are considered
  logical, parameter :: compare_raytracer = .True. ! modify Wavepacket simulation to facilitate comparison with raytracer
  namelist / iceList / inN, inQ, inQv, nVarIce, dt_ice, no_ice_source

  real, allocatable :: var_ww(:, :, :) ! flux <w'w'> from RayTracer

  !include_testoutput
  integer :: iVarO, iVarS, iVarAW

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

    ! Variables
    include_ice = .false.
    include_tracer = .false.
    include_testoutput = .false.

    ! Input/Output
    atmvarOut = ""
    rayvarOut = ""
    icevarOut = ""
    saverayvols = .false.
    prepare_restart = .false.
    restart = .false.
    iIn = - 1
    runName = "runName"
    outputType = "time"
    nOutput = 1
    maxIter = 1
    outputTimeDiff = 3600.0
    maxTime = 3600.0
    detailedinfo = .false.
    RHS_diagnostics = .false.
    fancy_namelists = .true.

    ! Debugging
    verbose = .false.
    dtMin_dim = 1.0e-6

    ! Test cases
    testCase = "SkamarockKlemp94"

    ! Monochromatic wave
    lambda_dim = 0.1 * (lz_dim(1) - lz_dim(0))

    ! Wave packet
    wavePacketType = 1
    wavePacketDim = 1
    lambdaX_dim = 0.1 * (lx_dim(1) - lx_dim(0))
    lambdaY_dim = 0.1 * (ly_dim(1) - ly_dim(0))
    lambdaZ_dim = 0.1 * (lz_dim(1) - lz_dim(0))
    amplitudeFactor = 1.0
    x0_dim = 0.5 * (lx_dim(0) + lx_dim(1))
    y0_dim = 0.5 * (ly_dim(0) + ly_dim(1))
    z0_dim = 0.5 * (lz_dim(0) + lz_dim(1))
    sigma_dim = 0.5 * (lz_dim(1) - lz_dim(0))
    sigma_hor_dim = 0.5 * (lx_dim(1) - lx_dim(0))
    amp_mod_x = 1.0
    sigma_hor_yyy_dim = 0.5 * (ly_dim(1) - ly_dim(0))
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
    lindUinit = .false.
    blocking = .false.
    nwm = 1
    launch_algorithm = "clip"
    zmin_wkb_dim = 0.0
    nray_fac = 20
    cons_merge = "en"
    nRayOutput = 10

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
    xextent_relax = 0.0
    yextent_relax = 0.0
    wind_relaxation = .false.
    surface_layer_depth = 0.0

    ! Baroclinic life cycle (realistic)
    zero_initial_state = .false.
    z_trpp0_dim = 10000.0
    z_baro_dim = 20000.0
    thet0_dim = 300.0
    ntrp_dim = 0.01
    nstr_dim = 0.02
    jwdth_dim = 0.0
    kaptpp = 0.0
    add_ptptb = .false.
    ptptb_x_dim = 0.5 * (lx_dim(0) + lx_dim(1))
    ptptb_y_dim = 0.5 * (ly_dim(0) + ly_dim(1))
    ptptb_z_dim = 0.5 * (lz_dim(0) + lz_dim(1))
    ptptb_dh_dim = 0.1 * min(lx_dim(1) - lx_dim(0), ly_dim(1) - ly_dim(0))
    ptptb_dz_dim = 0.1 * (lz_dim(1) - lz_dim(0))
    ptptb_amp_dim = 1.0
    add_noise = .false.
    proc_noise = 0.1
    tau_relax = 0.0
    tau_relax_low = 0.0
    sigma_tau = 0.0
    tau_jet = 0.0
    Sponge_Rel_Bal_Type = "env"
    ta_hs_dim = 0.0
    ts_hs_dim = 0.0
    tf_hs_dim = 0.0
    sigb_hs = 0.0

    ! Baroclinic life cycle (idealized)
    bar_sigma_y = 0.5 * (ly_dim(1) - ly_dim(0))
    u_strength = 30.0
    dTh_atm = 30.0
    init_2Dto3D = .false.
    init_bal = "geostr_id"
    lastrecordnum = 0
    fileinitstate2D = ""
    output_theta_bgr = .false.
    output_br_vais_sq = .false.
    output_heat = .false.
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
    preconditioner = "yes"
    dtau = 4.0e-4
    maxIterADI = 2
    initialCleaning = .true.
    pressureScaling = .false.
    correctMomentum = .true.
    correctDivError = .false.
    tolcrit = "abs"

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
    backgroundFlow_dim = 0.0
    f_Coriolis_dim = 0.0
    corset = "constant"
    z_tr_dim = 15000.0
    theta_tr_dim = 300.0
    gamma_t = 0.0
    gamma_s = 0.0
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
    range_factor = 1.0
    spectral_modes = 1

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
    relax_to_mean = .true.
    xBoundary = "periodic"
    yBoundary = "periodic"
    zBoundary = "solid_wall"

    ! WKB
    rayTracer = .false.

    ! Tracer
    tracerSetup = "alpha_z"
    include_trfrc_lo = .true.
    include_trfrc_mix = .true.
    include_trfrc_no = .true.

    ! Ice
    inN = 0
    inQ = 0
    inQv = 0
    nVarIce = 0
    dt_ice = 0.0
    no_ice_source = .false.

  end subroutine default_values

  subroutine write_namelists

    ! Write all namelists.

    implicit none

    integer, parameter :: parameter_space = 24, value_space = 23, &
        &comment_space = 72 - parameter_space - value_space
    character(len = 50) :: character_format, integer_format, logical_format, &
        &comment_format
    character(len = 2) :: counter
    integer :: iVar, jVar

    ! Adjust atmVarOut.
    jVar = 0
    do iVar = 1, size(atmVarOut)
      if(atmVarOut(iVar) == "" .or. (iVar > 1 .and. any(atmVarOut(:iVar - 1) &
          &== atmVarOut(iVar)))) cycle
      jVar = jVar + 1
      atmVarOut(jVar) = atmVarOut(iVar)
    end do
    atmVarOut(jVar + 1:) = ""

    ! Adjust rayVarOut.
    jVar = 0
    do iVar = 1, size(rayVarOut)
      if(rayVarOut(iVar) == "" .or. (iVar > 1 .and. any(rayVarOut(:iVar - 1) &
          &== rayVarOut(iVar)))) cycle
      jVar = jVar + 1
      rayVarOut(jVar) = rayVarOut(iVar)
    end do
    rayVarOut(jVar + 1:) = ""

    ! Adjust iceVarOut.
    jVar = 0
    do iVar = 1, size(iceVarOut)
      if(iceVarOut(iVar) == "" .or. (iVar > 1 .and. any(iceVarOut(:iVar - 1) &
          &== iceVarOut(iVar)))) cycle
      jVar = jVar + 1
      iceVarOut(jVar) = iceVarOut(iVar)
    end do
    iceVarOut(jVar + 1:) = ""

    ! Write namelists in standard format.
    if(master .and. .not. fancy_namelists) then
      open(unit = 90, file = "namelists.txt", action = "write", form &
          &= "formatted", status = "replace")
      write(unit = 90, nml = domain)
      write(90, "(a)") ""
      write(unit = 90, nml = variables)
      write(90, "(a)") ""
      write(unit = 90, nml = outputList)
      write(90, "(a)") ""
      write(unit = 90, nml = debuggingList)
      write(90, "(a)") ""
      write(unit = 90, nml = testCaseList)
      write(90, "(a)") ""
      write(unit = 90, nml = monochromeWave)
      write(90, "(a)") ""
      write(unit = 90, nml = wavePacket)
      write(90, "(a)") ""
      write(unit = 90, nml = LagrangeRayTracing)
      write(90, "(a)") ""
      write(unit = 90, nml = bubble)
      write(90, "(a)") ""
      write(unit = 90, nml = robert_bubble)
      write(90, "(a)") ""
      write(unit = 90, nml = mountainwavelist)
      write(90, "(a)") ""
      write(unit = 90, nml = baroclinic_LC)
      write(90, "(a)") ""
      write(unit = 90, nml = baroclinic_ID)
      write(90, "(a)") ""
      write(unit = 90, nml = modelList)
      write(90, "(a)") ""
      write(unit = 90, nml = solverList)
      write(90, "(a)") ""
      write(unit = 90, nml = poissonSolverList)
      write(90, "(a)") ""
      write(unit = 90, nml = atmosphereList)
      write(90, "(a)") ""
      write(unit = 90, nml = topographyList)
      write(90, "(a)") ""
      write(unit = 90, nml = boundaryList)
      write(90, "(a)") ""
      write(unit = 90, nml = boundaryList2)
      write(90, "(a)") ""
      write(unit = 90, nml = wkbList)
      write(90, "(a)") ""
      write(unit = 90, nml = tracerList)
      write(90, "(a)") ""
      write(unit = 90, nml = iceList)
      close(90)
      return
    end if

    ! Check if comment space is too small.
    if(comment_space <= 0) stop "Comment space too small! Decrease parameter &
        &space and/or value space!"

    ! Define character, integer and logical formats.
    write(counter, "(i2)") parameter_space
    character_format = "(2x, a" // trim(adjustl(counter)) // ", a3, a"
    integer_format = "(2x, a" // trim(adjustl(counter)) // ", a3, i"
    logical_format = "(2x, a" // trim(adjustl(counter)) // ", a3, l"
    write(counter, "(i2)") value_space
    character_format = trim(adjustl(character_format)) &
        &// trim(adjustl(counter)) // ", a3, a)"
    integer_format = trim(adjustl(integer_format)) // trim(adjustl(counter)) &
        &// ", a3, a)"
    logical_format = trim(adjustl(logical_format)) // trim(adjustl(counter)) &
        &// ", a3, a)"

    ! Define comment format.
    write(counter, "(i2)") parameter_space + value_space + 5
    comment_format = "(" // trim(adjustl(counter)) // "x, a3, a)"

    ! Open info file.
    if(master) then
      open(unit = 90, file = "namelists.txt", action = "write", form &
          &= "formatted", status = "replace")

      ! Write the header.
      write(90, "(a)") "!" // repeat("-", 78) // "!"
      write(90, "(a)") "!" // repeat(" ", 30) // "PincFlow Namelists" &
          &// repeat(" ", 30) // "!"
      write(90, "(a)") "!" // repeat("-", 78) // "!"
      write(90, "(a)") ""

      ! Write domain namelist.
      write(90, "(a)") "&domain"
      call write_integer("sizeX", sizeX, "Cells in x")
      call write_integer("sizeY", sizeY, "Cells in y")
      call write_integer("sizeZ", sizeZ, "Cells in z")
      call write_integer("nbx", nbx, "Halo/ghost cells in x")
      call write_integer("nby", nby, "Halo/ghost cells in y")
      call write_integer("nbz", nbz, "Ghost cells in z")
      call write_float("lx_dim(0)", lx_dim(0), "Minimum of x")
      call write_float("lx_dim(1)", lx_dim(1), "Maximum of x")
      call write_float("ly_dim(0)", ly_dim(0), "Minimum of y")
      call write_float("ly_dim(1)", ly_dim(1), "Maximum of y")
      call write_float("lz_dim(0)", lz_dim(0), "Minimum of z")
      call write_float("lz_dim(1)", lz_dim(1), "Maximum of z")
      call write_integer("nprocx", nprocx, "Processors in x")
      call write_integer("nprocy", nprocy, "Processors in y")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write variables namelist.
      write(90, "(a)") "&variables"
      call write_logical("include_ice", include_ice, "Ice microphysics &
          &parameterization")
      call write_logical("include_tracer", include_tracer, "Tracer equation")
      call write_logical("include_testoutput", include_testoutput, "...")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write output namelist.
      write(90, "(a)") "&outputList"
      do iVar = 1, size(atmvarOut)
        if(iVar == 1 .or. atmvarOut(iVar) /= "") then
          write(counter, "(i2)") iVar
          call write_character("atmvarOut(" // trim(adjustl(counter)) // ")", &
              &atmvarOut(iVar), "Atmospheric output variable")
        end if
      end do
      do iVar = 1, size(rayvarOut)
        if(iVar == 1 .or. rayvarOut(iVar) /= "") then
          write(counter, "(i2)") iVar
          call write_character("rayvarOut(" // trim(adjustl(counter)) // ")", &
              &rayvarOut(iVar), "Raytracer output variable")
        end if
      end do
      do iVar = 1, size(icevarOut)
        if(iVar == 1 .or. icevarOut(iVar) /= "") then
          write(counter, "(i2)") iVar
          call write_character("icevarOut(" // trim(adjustl(counter)) // ")", &
              &icevarOut(iVar), "Ice output variable")
        end if
      end do
      call write_logical("saverayvols", saverayvols, "Save ray volumes")
      call write_logical("prepare_restart", prepare_restart, "Save everything &
          &needed for a restart")
      call write_logical("restart", restart, "Restart the model from state in &
          &previous simulation")
      call write_integer("iIn", iIn, "Restart at time step iIn")
      call write_character("runName", runName, "Run name for netCDF file")
      call write_character("outputType", outputType, "'timeStep' or 'time'")
      call write_integer("nOutput", nOutput, "Output every nOutput time steps &
          &for outputType = 'timeStep'")
      call write_integer("maxIter", maxIter, "Stop after maxIter time steps &
          &for outputType = 'timeStep'")
      call write_float("outputTimeDiff", outputTimeDiff, "Output every &
          &outputTimeDiff seconds for outputType = 'time'")
      call write_float("maxTime", maxTime, "Stop after maxTime seconds for &
          &outputType = 'time'")
      call write_logical("detailedinfo", detailedinfo, "Provide info on the &
          &final state of Poisson solver")
      call write_logical("RHS_diagnostics", RHS_diagnostics, "Provide info &
          &about the RHS of Poisson equation")
      call write_logical("fancy_namelists", fancy_namelists, "Write all &
          &namelists with comments")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write debugging namelist.
      write(90, "(a)") "&debuggingList"
      call write_logical("verbose", verbose, "...")
      call write_float("dtMin_dim", dtMin_dim, "Stop if dt < dtMin")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write testcase namelist.
      write(90, "(a)") "&testCaseList"
      call write_character("testCase", testCase, "Predefined test case")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write monochromatic wave namelist.
      write(90, "(a)") "&monochromeWave"
      call write_float("lambda_dim", lambda_dim, "Vertical wavelength")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write wave packet namelist.
      write(90, "(a)") "&wavePacket"
      call write_integer("wavePacketType", wavePacketType, "1 = Gaussian, 2 &
          &= Cosine")
      call write_integer("wavePacketDim", wavePacketDim, "1 = 1D, 2 = 2D, 3 &
          &= 3D")
      call write_float("lambdaX_dim", lambdaX_dim, "Wavelength in x (0.0 for &
          &infinite wavelength)")
      call write_float("lambdaY_dim", lambdaY_dim, "Wavelength in y (0.0 for &
          &infinite wavelength)")
      call write_float("lambdaZ_dim", lambdaZ_dim, "Wavelength in z")
      call write_float("amplitudeFactor", amplitudeFactor, "Normalized &
          &buoyancy amplitude")
      call write_float("x0_dim", x0_dim, "Wave-packet center in x")
      call write_float("y0_dim", y0_dim, "Wave-packet center in y")
      call write_float("z0_dim", z0_dim, "Wave-packet center in z")
      call write_float("sigma_dim", sigma_dim, "Vertical width of Gaussian &
          &wave packet")
      call write_float("sigma_hor_dim", sigma_hor_dim, "Cosine distribution &
          &width in x (0.0 for infinite width)")
      call write_float("amp_mod_x", amp_mod_x, "Fractional amplitude &
          &modulation in x")
      call write_float("sigma_hor_yyy_dim", sigma_hor_yyy_dim, "Cosine &
          &distribution width in y (0.0 for infinite width)")
      call write_float("amp_mod_y", amp_mod_y, "Fractional amplitude &
          &modulation in y")
      call write_float("L_cos_dim", L_cos_dim, "Half width of vertical cosine &
          &profile of wave packet")
      call write_integer("omiSign", omiSign, "Frequency branch")
      call write_float("u0_jet_dim", u0_jet_dim, "Maximum amplitude of jet")
      call write_float("z0_jet_dim", z0_jet_dim, "Center of jet")
      call write_float("L_jet_dim", L_jet_dim, "Half width of vertical cosine &
          &profile of jet")
      call write_logical("inducedwind", inducedwind, "...")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write Lagrange ray tracing namelist.
      write(90, "(a)") "&LagrangeRayTracing"
      call write_float("xrmin_dim", xrmin_dim, "Left bound of initial rays")
      call write_float("xrmax_dim", xrmax_dim, "Right bound of initial rays")
      call write_float("yrmin_dim", yrmin_dim, "Backward bound of initial rays")
      call write_float("yrmax_dim", yrmax_dim, "Forward bound of initial rays")
      call write_float("zrmin_dim", zrmin_dim, "Bottom bound of initial rays")
      call write_float("zrmax_dim", zrmax_dim, "Top bound of initial rays")
      call write_integer("nrxl", nrxl, "Initial ray volumes within dx")
      call write_integer("nryl", nryl, "Initial ray volumes within dy")
      call write_integer("nrzl", nrzl, "Initial ray volumes within dz")
      call write_float("fac_dk_init", fac_dk_init, "Initial fraction dk/kh")
      call write_float("fac_dl_init", fac_dl_init, "Initial fraction dl/kh")
      call write_float("fac_dm_init", fac_dm_init, "Initial fraction dm/m")
      call write_integer("nrk_init", nrk_init, "Initial ray volumes within dk")
      call write_integer("nrl_init", nrl_init, "Initial ray volumes within dl")
      call write_integer("nrm_init", nrm_init, "Initial ray volumes within dm")
      call write_integer("nsmth_wkb", nsmth_wkb, "Width of smoothing operator &
          &for mean flow tendencies")
      call write_logical("lsmth_wkb", lsmth_wkb, "Smoothing operator for mean &
          &flow tendencies")
      call write_integer("sm_filter", sm_filter, "1 = Box filter, 2 = Shapiro &
          &filter")
      call write_logical("lsaturation", lsaturation, "Switch for saturation &
          &scheme")
      call write_float("alpha_sat", alpha_sat, "Saturation threshold")
      call write_logical("steady_state", steady_state, "Steady-state mode")
      call write_integer("case_wkb", case_wkb, "1 = Gaussian wave packet, 2 &
          &= cosine wave packet, 3 = mountain wave")
      call write_float("amp_wkb", amp_wkb, "Relative amplitude of wave packet")
      call write_float("wlrx_init", wlrx_init, "Wavelength in x of wave packet")
      call write_float("wlry_init", wlry_init, "Wavelength in y of wave packet")
      call write_float("wlrz_init", wlrz_init, "Wavelength in z of wave packet")
      call write_float("xr0_dim", xr0_dim, "Wave packet center in x")
      call write_float("yr0_dim", yr0_dim, "Wave packet center in y")
      call write_float("zr0_dim", zr0_dim, "Wave packet center in z")
      call write_float("sigwpx_dim", sigwpx_dim, "Wave packet width in x (0.0 &
          &for infinite width)")
      call write_float("sigwpy_dim", sigwpy_dim, "Wave packet width in y (0.0 &
          &for infinite width)")
      call write_float("sigwpz_dim", sigwpz_dim, "Wave packet width in z (0.0 &
          &for infinite width)")
      call write_integer("branchr", branchr, "Frequency branch")
      call write_logical("lindUinit", lindUinit, "Induced wind at initial time")
      call write_logical("blocking", blocking, "Simple blocked-layer scheme")
      call write_integer("nwm", nwm, "Number of initial wave modes")
      call write_character("launch_algorithm", launch_algorithm, "Ray-volume &
          &launch algorithm")
      call write_float("zmin_wkb_dim", zmin_wkb_dim, "Minimum altitude for &
          &wave-mean-flow interaction")
      call write_integer("nray_fac", nray_fac, "Maximum multiplication factor &
          &per spectral dimension")
      call write_character("cons_merge", cons_merge, "Conserved quantity in &
          &ray-volume merging ('wa' = wave action, 'en' = wave energy)")
      call write_integer("nRayOutput", nRayOutput, "Number of dominant ray &
          &volumes in output")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write bubble namelist.
      write(90, "(a)") "&bubble"
      call write_float("dTheta0_dim", dTheta0_dim, "...")
      call write_float("xRadius_dim", xRadius_dim, "...")
      call write_float("zRadius_dim", zRadius_dim, "...")
      call write_float("xCenter_dim", xCenter_dim, "...")
      call write_float("yCenter_dim", yCenter_dim, "...")
      call write_float("zCenter_dim", zCenter_dim, "...")
      call write_float("zExcentricity", zExcentricity, "...")
      call write_float("rhoCenter_dim", rhoCenter_dim, "...")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write Robert bubble namelist.
      write(90, "(a)") "&robert_bubble"
      call write_float("dTheta1_dim", dTheta1_dim, "Potential temperature &
          &offset")
      call write_float("a1_dim", a1_dim, "Radius of plateau")
      call write_float("sigma1_dim", sigma1_dim, "Gaussian edge profile")
      call write_float("xCenter1_dim", xCenter1_dim, "...")
      call write_float("zCenter1_dim", zCenter1_dim, "...")
      call write_float("dTheta2_dim", dTheta2_dim, "...")
      call write_float("a2_dim", a2_dim, "...")
      call write_float("sigma2_dim", sigma2_dim, "...")
      call write_float("xCenter2_dim", xCenter2_dim, "...")
      call write_float("zCenter2_dim", zCenter2_dim, "...")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write mountain wave namelist.
      write(90, "(a)") "&mountainwavelist"
      call write_float("u_relax", u_relax, "Zonal relaxation wind")
      call write_float("v_relax", v_relax, "Meridional relaxation wind")
      call write_float("w_relax", w_relax, "Vertical relaxation wind")
      call write_float("t_relax", t_relax, "Relaxation time")
      call write_float("t_ramp", t_ramp, "Not used at the moment")
      call write_float("xextent_relax", xextent_relax, "Zonal extent of &
          &relaxation region")
      call write_float("yextent_relax", yextent_relax, "Meridional extent of &
          &relaxation region")
      call write_logical("wind_relaxation", wind_relaxation, "Relaxation &
          &switch")
      call write_float("surface_layer_depth", surface_layer_depth, &
          &"Surface-layer depth")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write baroclinic LC namelist.
      write(90, "(a)") "&baroclinic_LC"
      call write_logical("zero_initial_state", zero_initial_state, "...")
      call write_float("z_trpp0_dim", z_trpp0_dim, "Mean tropopause height")
      call write_float("z_baro_dim", z_baro_dim, "Altitude above which the &
          &atmosphere is barotropic")
      call write_float("thet0_dim", thet0_dim, "Characteristic potential &
          &temperature")
      call write_float("ntrp_dim", ntrp_dim, "Buoyancy frequency of the &
          &troposphere")
      call write_float("nstr_dim", nstr_dim, "Buoyancy frequency of the &
          &stratosphere")
      call write_float("jwdth_dim", jwdth_dim, "Jet half width")
      call write_float("kaptpp", kaptpp, "Jet slope")
      call write_logical("add_ptptb", add_ptptb, "Add local &
          &potential-temperature perturbation")
      call write_float("ptptb_x_dim", ptptb_x_dim, "Location in x of local &
          &potential-temperature perturbation")
      call write_float("ptptb_y_dim", ptptb_y_dim, "Location in y of local &
          &potential-temperature perturbation")
      call write_float("ptptb_z_dim", ptptb_z_dim, "Location in z of local &
          &potential-temperature perturbation")
      call write_float("ptptb_dh_dim", ptptb_dh_dim, "Horizontal width of &
          &local potential-temperature perturbation")
      call write_float("ptptb_dz_dim", ptptb_dz_dim, "Vertical width of local &
          &potential-temperature perturbation")
      call write_float("ptptb_amp_dim", ptptb_amp_dim, "Amplitude of local &
          &potential-temperature perturbation")
      call write_logical("add_noise", add_noise, "Add noise to the initial &
          &potential temperature")
      call write_float("proc_noise", proc_noise, "Relative amplitude of &
          &potential-temperature noise")
      call write_float("tau_relax", tau_relax, "...")
      call write_float("tau_relax_low", tau_relax_low, "...")
      call write_float("sigma_tau", sigma_tau, "Relative thickness of &
          &relaxation profile")
      call write_float("tau_jet", tau_jet, "Time for jet formation")
      call write_character("Sponge_Rel_Bal_Type", Sponge_Rel_Bal_Type, &
          &"Relaxation state: 'hyd' = hydrostatic balance, 'env' = geostrophic &
          &balance")
      call write_float("ta_hs_dim", ta_hs_dim, "Thermal-relaxation time scale &
          &outside tropical boundary layer (zero means infinity)")
      call write_float("ts_hs_dim", ts_hs_dim, "Thermal-relaxation time scale &
          &in tropical boundary layer (zero means infinity)")
      call write_float("tf_hs_dim", tf_hs_dim, "Boundary-layer &
          &Rayleigh-damping  time scale")
      call write_float("sigb_hs", sigb_hs, "Sigma of boundary-layer top")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write baroclinic ID namelist.
      write(90, "(a)") "&baroclinic_ID"
      call write_float("bar_sigma_y", bar_sigma_y, "...")
      call write_float("u_strength", u_strength, "Jet strength")
      call write_float("dTh_atm", dTh_atm, "Potential-temperature difference &
          &between poles and tropics")
      call write_logical("init_2Dto3D", init_2Dto3D, "Initialize 3D state from &
          &2D balanced state")
      call write_character("init_bal", init_bal, "...")
      call write_integer("lastrecordnum", lastrecordnum, "Last record in 2D &
          &output")
      call write_character("fileinitstate2D", fileinitstate2D, "File name with &
          &2D balanced state")
      call write_logical("output_theta_bgr", output_theta_bgr, "Output &
          &environmental potential temperature")
      call write_logical("output_br_vais_sq", output_br_vais_sq, "Output &
          &environmental squared buoyancy frequency")
      call write_logical("output_heat", output_heat, "...")
      call write_character("balance_eq", balance_eq, "...")
      call write_logical("output_rho_bgr", output_rho_bgr, "Output &
          &environmental density")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write model namelist.
      write(90, "(a)") "&modelList"
      call write_character("model", model, "Dynamic equations")
      call write_float("vert_theta", vert_theta, "Rotation about x")
      call write_float("vert_alpha", vert_alpha, "Rotation about y")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write solver namelist.
      write(90, "(a)") "&solverList"
      call write_float("cfl", cfl, "CFL number")
      call write_float("cfl_wave", cfl_wave, "WKB CFL number")
      call write_float("dtMax_dim", dtMax_dim, "Maximum time step")
      call write_character("tStepChoice", tStepChoice, "'fix' or 'cfl'")
      call write_character("timeScheme", timeScheme, "'LS_Will_RK3' or &
          &'semiimplicit'")
      call write_logical("auxil_equ", auxil_equ, "Buoyancy equation")
      call write_character("fluxType", fluxType, "'ILES', 'central' or &
          &'upwind'")
      call write_character("reconstType", reconstType, "'MUSCL', 'constant', &
          &'SALD' or 'ALDM'")
      call write_character("musclType", musclType, "'muscl1' or 'muscl2'")
      call write_character("limiterType1", limiterType1, "'minmod', &
          &'MCVariant' or 'Cada'")
      call write_logical("TurbScheme", TurbScheme, "Turbulence scheme")
      call write_float("turb_dts", turb_dts, "Turbulent damping time")
      call write_logical("DySmaScheme", DySmaScheme, "Dynamic Smagorinsky &
          &scheme")
      call write_logical("dtWave_on", dtWave_on, "Limit time step by inverse &
          &buoyancy frequency")
      call write_logical("heatingONK14", heatingONK14, "Heating as in O'Neill &
          &and Klein (2014)")
      call write_logical("dens_relax", dens_relax, "Heating by density &
          &relaxation")
      call write_float("shap_dts_fac", shap_dts_fac, "Shapiro-filter damping &
          &time")
      call write_integer("n_shap", n_shap, "Order of Shapiro filter")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write Poisson solver namelist.
      write(90, "(a)") "&poissonSolverList"
      call write_float("tolPoisson", tolPoisson, "Abort criterion")
      call write_float("abs_tol", abs_tol, "Lower bound for tolerance")
      call write_float("tolCond", tolCond, "Preconditioner tolerance")
      call write_integer("maxIterPoisson", maxIterPoisson, "Maximum iterations")
      call write_character("poissonSolverType", poissonSolverType, &
          &"'bicgstab', 'gcr', 'adi' or 'hypre'")
      call write_character("preconditioner", preconditioner, "'no' or 'yes'")
      call write_float("dtau", dtau, "Time parameter for preconditioner")
      call write_integer("maxIterADI", maxIterADI, "Preconditioner iterations")
      call write_logical("initialCleaning", initialCleaning, "Enforce initial &
          &non-divergence")
      call write_logical("pressureScaling", pressureScaling, "Scale by P")
      call write_logical("correctMomentum", correctMomentum, "Correct momentum &
          &so that divergence constraint is fulfilled")
      call write_logical("correctDivError", correctDivError, "Subtract &
          &divergence")
      call write_character("tolcrit", tolcrit, "'abs' or 'rel'")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write atmosphere namelist.
      write(90, "(a)") "&atmosphereList"
      call write_character("referenceQuantities", referenceQuantities, &
          &"'Klein', 'WKB', 'SI' or 'general'")
      call write_logical("specifyReynolds", specifyReynolds, "Use inverse &
          &Reynolds number")
      call write_float("ReInv", ReInv, "Inverse Reynolds number")
      call write_float("mu_viscous_dim", mu_viscous_dim, "Kinematic viscosity")
      call write_float("mu_conduct_dim", mu_conduct_dim, "Heat conductivity")
      call write_character("background", background, "'realistic', &
          &'isothermal', 'isentropic', 'const-N', 'diflapse' or 'HeldSuarez'")
      call write_float("N_BruntVaisala_dim", N_BruntVaisala_dim, "Buoyancy &
          &frequency for 'const-N'")
      call write_float("theta0_dim", theta0_dim, "Background potential &
          &temperature for 'isentropic'")
      call write_float("Temp0_dim", Temp0_dim, "Background temperature for &
          &'isothermal'")
      call write_float("press0_dim", press0_dim, "Ground pressure")
      do iVar = 1, 3
        write(counter, "(i2)") iVar
        call write_float("backgroundFlow_dim(" // trim(adjustl(counter)) &
            &// ")", backgroundFlow_dim(iVar), "Initial wind")
      end do
      call write_float("f_Coriolis_dim", f_Coriolis_dim, "Coriolis frequency")
      call write_character("corset", corset, "'constant' or 'periodic'")
      call write_float("z_tr_dim", z_tr_dim, "Tropopause height")
      call write_float("theta_tr_dim", theta_tr_dim, "Potential temperature in &
          &troposphere")
      call write_float("gamma_t", gamma_t, "Lapse rate in troposphere")
      call write_float("gamma_s", gamma_s, "Lapse rate in stratosphere")
      call write_float("tp_strato_dim", tp_strato_dim, "Temperature in &
          &stratosphere for 'HeldSuarez'")
      call write_float("tp_srf_trp_dim", tp_srf_trp_dim, "Tropical surface &
          &temperature for 'HeldSuarez'")
      call write_float("tpdiffhor_tropo_dim", tpdiffhor_tropo_dim, &
          &"Temperature difference between poles and tropics for 'HeldSuarez'")
      call write_float("ptdiffvert_tropo_dim", ptdiffvert_tropo_dim, "Vertical &
          &potential temperature difference in troposphere for 'HeldSuarez'")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write topography namelist.
      write(90, "(a)") "&topographyList"
      call write_logical("topography", topography, "Terrain-following &
          &coordinates")
      call write_integer("ipolTFC", ipolTFC, "Interpolation in the &
          &transformation of w")
      call write_logical("freeSlipTFC", freeSlipTFC, "Transformed free-slip &
          &condition")
      call write_logical("testTFC", testTFC, "Various TFC tests")
      call write_float("topographyTime", topographyTime, "Topography growth &
          &time")
      call write_float("mountainHeight_dim", mountainHeight_dim, "Maximum &
          &height")
      call write_float("mountainWidth_dim", mountainWidth_dim, "Half width")
      call write_integer("mountain_case", mountain_case, "Predefined &
          &topography")
      call write_float("range_factor", range_factor, "Ratio between large and &
          &small scales")
      call write_integer("spectral_modes", spectral_modes, "Number of spectral &
          &modes")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write boundary namelist.
      write(90, "(a)") "&boundaryList"
      call write_logical("rhoFluxCorr", rhoFluxCorr, "...")
      call write_logical("iceFluxCorr", iceFluxCorr, "...")
      call write_logical("uFluxCorr", uFluxCorr, "...")
      call write_logical("vFluxCorr", vFluxCorr, "...")
      call write_logical("wFluxCorr", wFluxCorr, "...")
      call write_logical("thetaFluxCorr", thetaFluxCorr, "...")
      call write_integer("nbCellCorr", nbCellCorr, "...")
      call write_logical("spongeLayer", spongeLayer, "General sponge layer &
          &switch")
      call write_logical("sponge_uv", sponge_uv, "Sponge layer for horizontal &
          &wind if unifiedSponge = .false.")
      call write_float("spongeHeight", spongeHeight, "Relative height of lower &
          &sponge layer edge (scale height for unifiedSponge = .true. and &
          &spongeType = 'exponential')")
      call write_float("spongeAlphaZ_dim", spongeAlphaZ_dim, "Maximum &
          &relaxation rate for unifiedSponge = .true.")
      call write_float("spongeAlphaZ_fac", spongeAlphaZ_fac, "Sponge layer &
          &factor for unifiedSponge = .false.")
      call write_logical("unifiedSponge", unifiedSponge, "Unified sponge for &
          &both time schemes, applied to wind and density")
      call write_logical("lateralSponge", lateralSponge, "Lateral sponge for &
          &unifiedSponge = .true.")
      call write_character("spongeType", spongeType, "Sponge layer profile for &
          &unifiedSponge = .true.")
      call write_integer("spongeOrder", spongeOrder, "Order of polynomial &
          &sponge")
      call write_integer("cosmoSteps", cosmoSteps, "Relative strength of COSMO &
          &sponge")
      call write_logical("relax_to_mean", relax_to_mean, "Relax the wind to &
          &its (terrain-following) horizontal mean (otherwise, relax to the &
          &initial state) if unifiedSponge == .true.")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write second boundary namelist.
      write(90, "(a)") "&boundaryList2"
      call write_character("xBoundary", xBoundary, "Boundary conditions in x &
          &('periodic' only)")
      call write_character("yBoundary", yBoundary, "Boundary conditions in y &
          &('periodic' only)")
      call write_character("zBoundary", zBoundary, "Boundary conditions in z &
          &('periodic' or 'solid_wall')")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write WKB namelist.
      write(90, "(a)") "&wkbList"
      call write_logical("rayTracer", rayTracer, "Ray-tracer switch")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write tracer namelist.
      write(90, "(a)") "&tracerList"
      call write_character("tracerSetup", tracerSetup, "initial tracer &
          &distribution")
      call write_logical("include_trfrc_lo", include_trfrc_lo, "leading-order &
          &GW tracer forcing")
      call write_logical("include_trfrc_no", include_trfrc_no, "next-order GW &
          &tracer forcing")
      call write_logical("include_trfrc_mix", include_trfrc_mix, "diffusive &
          &tracer mixing")
      write(90, "(a)") "&end"
      write(90, "(a)") ""

      ! Write ice namelist.
      write(90, "(a)") "&iceList"
      call write_integer("inN", inN, "...")
      call write_integer("inQ", inQ, "...")
      call write_integer("inQv", inQv, "...")
      call write_integer("nVarIce", nVarIce, "...")
      call write_float("dt_ice", dt_ice, "...")
      call write_logical("no_ice_source", no_ice_source, "...")
      write(90, "(a)") "&end"

      ! Close info file.
      close(90)
    end if

    contains

    subroutine write_character(parameter, value, comment)

      character(len = *), intent(in) :: parameter, value, comment
      character(len = len_trim(adjustl(parameter))) :: parameter_trimmed
      character(len = len_trim(adjustl(value))) :: value_trimmed
      character(len = len_trim(adjustl(comment))) :: comment_trimmed
      character(len = parameter_space) :: parameter_output
      character(len = value_space) :: value_output
      character(len = comment_space) :: comment_output
      integer :: position

      ! Trim input.
      parameter_trimmed = trim(adjustl(parameter))
      value_trimmed = trim(adjustl(value))
      comment_trimmed = trim(adjustl(comment))

      ! Check if parameter space and value space are large enough.
      if(len(parameter_trimmed) > parameter_space) stop "Parameter space too &
          &small!"
      if(len(value_trimmed) + 2 > value_space) stop "Value space too small!"

      ! Set parameter output and value output.
      parameter_output = adjustl(parameter_trimmed)
      value_output = repeat(" ", value_space - len(value_trimmed) - 2) // "'" &
          &// value_trimmed // "'"

      ! Set comment output.
      if(len(comment_trimmed) > comment_space) then
        position = max(index(comment_trimmed(:comment_space), " ", back &
            &= .true.), index(comment_trimmed(:comment_space), "-", back &
            &= .true.), index(comment_trimmed(:comment_space), "/", back &
            &= .true.))
        if(position == 0) stop "Hyphenation required!"
        comment_output = adjustl(comment_trimmed(:position))
      else
        comment_output = adjustl(comment_trimmed)
      end if

      ! Write namelist line.
      write(90, character_format) parameter_output, " = ", value_output, " ! &
          &", trim(comment_output)

      ! Continue comment.
      if(len(comment_trimmed) > comment_space) then
        call write_comment(comment_trimmed, position)
      end if

    end subroutine write_character

    subroutine write_integer(parameter, value, comment)

      character(len = *), intent(in) :: parameter, comment
      integer, intent(in) :: value
      character(len = len_trim(adjustl(parameter))) :: parameter_trimmed
      character(len = len_trim(adjustl(comment))) :: comment_trimmed
      character(len = parameter_space) :: parameter_output
      character(len = comment_space) :: comment_output
      integer :: position

      ! Trim input.
      parameter_trimmed = trim(adjustl(parameter))
      comment_trimmed = trim(adjustl(comment))

      ! Check if parameter space and value space are large enough.
      if(len(parameter_trimmed) > parameter_space) stop "Parameter space too &
          &small!"
      if(integer_size(value) > value_space) stop "Value space too small!"

      ! Set parameter output.
      parameter_output = adjustl(parameter_trimmed)

      ! Set comment output.
      if(len(comment_trimmed) > comment_space) then
        position = max(index(comment_trimmed(:comment_space), " ", back &
            &= .true.), index(comment_trimmed(:comment_space), "-", back &
            &= .true.), index(comment_trimmed(:comment_space), "/", back &
            &= .true.))
        if(position == 0) stop "Hyphenation required!"
        comment_output = adjustl(comment_trimmed(:position))
      else
        comment_output = adjustl(comment_trimmed)
      end if

      ! Write namelist line.
      write(90, integer_format) parameter_output, " = ", value, " ! ", &
          &trim(comment_output)

      ! Continue comment.
      if(len(comment_trimmed) > comment_space) then
        call write_comment(comment_trimmed, position)
      end if

    end subroutine write_integer

    subroutine write_float(parameter, value, comment)

      character(len = *), intent(in) :: parameter, comment
      real, intent(in) :: value
      character(len = len_trim(adjustl(parameter))) :: parameter_trimmed
      character(len = len_trim(adjustl(comment))) :: comment_trimmed
      character(len = parameter_space) :: parameter_output
      character(len = 50) :: float_format
      character(len = comment_space) :: comment_output
      character(len = 2) counter
      integer, dimension(1:2) :: digit_count
      integer :: position

      ! Trim input.
      parameter_trimmed = trim(adjustl(parameter))
      comment_trimmed = trim(adjustl(comment))

      ! Compute significant digits.
      digit_count = significant_digits(value)

      ! Check if parameter space and value space are large enough.
      if(len(parameter_trimmed) > parameter_space) stop "Parameter space too &
          &small!"
      if(sum(digit_count) + 5 > value_space) stop "Value space too small!"

      ! Set parameter output.
      parameter_output = adjustl(parameter_trimmed)

      ! Define float format.
      write(counter, "(i2)") parameter_space
      float_format = "(2x, a" // trim(adjustl(counter)) // ", a3, es"
      write(counter, "(i2)") value_space
      float_format = trim(adjustl(float_format)) // trim(adjustl(counter)) &
          &// "."
      write(counter, "(i2)") digit_count(1)
      float_format = trim(adjustl(float_format)) // trim(adjustl(counter)) &
          &// "e"
      write(counter, "(i2)") digit_count(2)
      float_format = trim(adjustl(float_format)) // trim(adjustl(counter)) &
          &// ", a3, a)"

      ! Set comment output.
      if(len(comment_trimmed) > comment_space) then
        position = max(index(comment_trimmed(:comment_space), " ", back &
            &= .true.), index(comment_trimmed(:comment_space), "-", back &
            &= .true.), index(comment_trimmed(:comment_space), "/", back &
            &= .true.))
        if(position == 0) stop "Hyphenation required!"
        comment_output = adjustl(comment_trimmed(:position))
      else
        comment_output = adjustl(comment_trimmed)
      end if

      ! Write namelist line.
      write(90, float_format) parameter_output, " = ", value, " ! ", &
          &trim(comment_output)

      ! Continue comment.
      if(len(comment_trimmed) > comment_space) then
        call write_comment(comment_trimmed, position)
      end if

    end subroutine write_float

    subroutine write_logical(parameter, value, comment)

      character(len = *), intent(in) :: parameter, comment
      logical, intent(in) :: value
      character(len = len_trim(adjustl(parameter))) :: parameter_trimmed
      character(len = len_trim(adjustl(comment))) :: comment_trimmed
      character(len = parameter_space) :: parameter_output
      character(len = comment_space) :: comment_output
      integer :: position

      ! Trim input.
      parameter_trimmed = trim(adjustl(parameter))
      comment_trimmed = trim(adjustl(comment))

      ! Check if parameter space is large enough.
      if(len(parameter_trimmed) > parameter_space) stop "Parameter space too &
          &small!"

      ! Set parameter output.
      parameter_output = adjustl(parameter_trimmed)

      ! Set comment output.
      if(len(comment_trimmed) > comment_space) then
        position = max(index(comment_trimmed(:comment_space), " ", back &
            &= .true.), index(comment_trimmed(:comment_space), "-", back &
            &= .true.), index(comment_trimmed(:comment_space), "/", back &
            &= .true.))
        if(position == 0) stop "Hyphenation required!"
        comment_output = adjustl(comment_trimmed(:position))
      else
        comment_output = adjustl(comment_trimmed)
      end if

      ! Write namelist line.
      write(90, logical_format) parameter_output, " = ", value, " ! ", &
          &trim(comment_output)

      ! Continue comment.
      if(len(comment_trimmed) > comment_space) then
        call write_comment(comment_trimmed, position)
      end if

    end subroutine write_logical

    subroutine write_comment(comment, position)

      character(len = *), intent(in) :: comment
      integer, intent(in) :: position
      character(len = comment_space) :: comment_output
      integer :: left, right, shift

      left = position + 1
      do while(len_trim(comment(left:)) > 0)
        comment_output = adjustl(comment(left:))
        shift = max(index(comment_output, " ", back = .true.), &
            &index(comment_output, "-", back = .true.), index(comment_output, &
            &"/", back = .true.))
        if(shift == 0) stop "Hyphenation required!"
        comment_output = adjustl(comment_output(:shift))
        write(90, comment_format) " ! ", trim(comment_output)
        right = left + shift - 1
        left = right + 1
      end do

    end subroutine write_comment

  end subroutine write_namelists

  function trim_float(input) result(output)

    real, intent(in) :: input
    integer, dimension(1:2) :: digit_count
    character(len = 50) :: float_format
    character(len = 2) :: counter
    character(len = max(0, int(sign(1.0, - input))) &
        &+ sum(significant_digits(input)) + 4) :: output

    digit_count = significant_digits(input)

    write(counter, "(i2)") len(output)
    float_format = "(es" // trim(adjustl(counter)) // "."
    write(counter, "(i2)") digit_count(1)
    float_format = trim(adjustl(float_format)) // trim(adjustl(counter)) // "e"
    write(counter, "(i2)") digit_count(2)
    float_format = trim(adjustl(float_format)) // counter // ")"

    write(output, float_format) input

  end function trim_float

  function trim_integer(input) result(output)

    integer, intent(in) :: input
    character(len = 20) :: number
    character(len = integer_size(input)) :: output

    write(number, "(i20)") input
    output = trim(adjustl(number))

  end function trim_integer

  function trim_logical(input) result(output)

    logical, intent(in) :: input
    character(len = 1) :: output

    write(output, "(l1)") input

  end function trim_logical

  pure function convert_case(input, choice) result(output)

    ! Convert either lower to upper or upper to lower case.

    implicit none

    character(len = *), intent(in) :: input, choice
    character(len = len(input)) :: output
    character(len = *), parameter :: lower = "abcdefghijklmnopqrstuvwxyz"
    character(len = *), parameter :: upper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    integer :: i, j

    output = input
    if(choice == "upper") then
      do i = 1, len(output)
        j = index(lower, output(i:i))
        if(j > 0) output(i:i) = upper(j:j)
      end do
    else if(choice == "lower") then
      do i = 1, len(output)
        j = index(upper, output(i:i))
        if(j > 0) output(i:i) = lower(j:j)
      end do
    end if

  end function

  pure function significant_digits(input) result(digit_count)

    ! Determine the significant digits of a 64-bit float (the 16th digit is
    ! rounded).

    implicit none

    real, intent(in) :: input
    character(len = 23) :: number
    integer, dimension(1:2) :: digit_count
    integer :: position

    write(number, "(es23.15e3)") input
    digit_count = 0
    do position = 21, 23
      if(number(position:position) /= "0") exit
      digit_count(2) = digit_count(2) + 1
    end do
    digit_count(2) = max(1, 3 - digit_count(2))
    do position = 18, 4, - 1
      if(number(position:position) /= "0") exit
      digit_count(1) = digit_count(1) + 1
    end do
    digit_count(1) = max(1, 15 - digit_count(1))

  end function significant_digits

  pure function integer_size(input) result(size)

    ! Determine the size of a 64-bit integer.

    implicit none

    integer, intent(in) :: input
    character(len = 20) :: number
    integer :: size

    write(number, "(i20)") input
    size = len_trim(adjustl(number))

  end function integer_size

  subroutine allocate_var_type(var)

    implicit none

    type(var_type), intent(inout) :: var
    integer :: allocstat

    ! Allocate basic variables.
    allocate(var%rho(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "Allocation of var_type component rho failed!"
    allocate(var%u(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "Allocation of var_type component u failed!"
    allocate(var%v(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "Allocation of var_type component v failed!"
    allocate(var%w(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "Allocation of var_type component w failed!"
    allocate(var%pi(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "Allocation of var_type component pi failed!"
    allocate(var%rhop(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "Allocation of var_type component rhop failed!"

    ! Allocate dynamic Smagorinsky coefficient.
    if(turbScheme) then
      allocate(var%DSC(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "Allocation of var_type component DSC failed!"
    end if

    ! Allocate gravity-wave heating.
    if(rayTracer) then
      allocate(var%GWH(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "Allocation of var_type component GWH failed!"
    end if

    ! Allocate mass-weighted potential temperature.
    if(model == "compressible") then
      allocate(var%P(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "Allocation of var_type component P failed!"
    end if

    ! Allocate tracer.
    if(include_tracer) then
      allocate(var%chi(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "Allocation of var_type component chi failed!"
    end if

    ! Allocate alternative ice variables.
    if(include_ice) then
      allocate(var%ICE(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, &
          &1:nVarIce), stat = allocstat)
      if(allocstat /= 0) stop "Allocation of var_type component ICE failed!"
    end if

    ! Allocate optional variables.
    if(include_testoutput) then
      allocate(var%OPT(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3), &
          &stat = allocstat)
      if(allocstat /= 0) stop "Allocation of var_type component OPT failed!"
    end if

  end subroutine allocate_var_type

  subroutine allocate_flux_type(flux)

    implicit none

    type(flux_type), intent(inout) :: flux
    integer :: allocstat

    ! Allocate basic fluxes.
    allocate(flux%rho(- 1:nx, - 1:ny, - 1:nz, 1:3), stat = allocstat)
    if(allocstat /= 0) stop "Allocation of flux_type component rho failed!"
    allocate(flux%u(- 1:nx, - 1:ny, - 1:nz, 1:3), stat = allocstat)
    if(allocstat /= 0) stop "Allocation of flux_type component u failed!"
    allocate(flux%v(- 1:nx, - 1:ny, - 1:nz, 1:3), stat = allocstat)
    if(allocstat /= 0) stop "Allocation of flux_type component v failed!"
    allocate(flux%w(- 1:nx, - 1:ny, - 1:nz, 1:3), stat = allocstat)
    if(allocstat /= 0) stop "Allocation of flux_type component w failed!"
    allocate(flux%theta(- 1:nx, - 1:ny, - 1:nz, 1:3), stat = allocstat)
    if(allocstat /= 0) stop "Allocation of flux_type component theta failed!"
    allocate(flux%rhop(- 1:nx, - 1:ny, - 1:nz, 1:3), stat = allocstat)
    if(allocstat /= 0) stop "Allocation of flux_type component rhop failed!"

    ! Allocate mass-weighted-potential-temperature fluxes.
    if(model == "compressible") then
      allocate(flux%P(- 1:nx, - 1:ny, - 1:nz, 1:3), stat = allocstat)
      if(allocstat /= 0) stop "Allocation of flux_type component P failed!"
    end if

    ! Allocate tracer fluxes.
    if(include_tracer) then
      allocate(flux%chi(- 1:nx, - 1:ny, - 1:nz, 1:3), stat = allocstat)
      if(allocstat /= 0) stop "Allocation of flux_type component chi failed!"
    end if

    ! Allocate alternative ice fluxes.
    if(include_ice) then
      allocate(flux%ICE(- 1:nx, - 1:ny, - 1:nz, 1:3, 1:nVarIce), stat &
          &= allocstat)
      if(allocstat /= 0) stop "Allocation of flux_type component ICE failed!"
    end if

  end subroutine allocate_flux_type

  subroutine reset_var_type(var)

    implicit none

    type(var_type), intent(inout) :: var

    ! Reset basic variables.
    var%rho = 0.0
    var%u = 0.0
    var%v = 0.0
    var%w = 0.0
    var%pi = 0.0
    var%rhop = 0.0

    ! Reset dynamic Smagorinsky coefficient.
    if(turbScheme) then
      var%DSC = 0.0
    end if

    ! Reset gravity-wave heating.
    if(rayTracer) then
      var%GWH = 0.0
    end if

    ! Reset mass-weighted potential temperature.
    if(model == "compressible") then
      var%P = 0.0
    end if

    ! Reset tracer.
    if(include_tracer) then
      var%chi = 0.0
    end if

    ! Reset ice variables.
    if(include_ice) then
      var%ICE = 0.0
    end if

    ! Reset optional variables.
    if(include_testoutput) then
      var%OPT = 0.0
    end if

  end subroutine reset_var_type

  subroutine reset_flux_type(flux)

    implicit none

    type(flux_type), intent(inout) :: flux

    ! Reset basic fluxes.
    flux%rho = 0.0
    flux%u = 0.0
    flux%v = 0.0
    flux%w = 0.0
    flux%theta = 0.0
    flux%rhop = 0.0

    ! Reset mass-weighted-potential-temperature fluxes.
    if(model == "compressible") then
      flux%P = 0.0
    end if

    ! Reset tracer fluxes.
    if(include_tracer) then
      flux%chi = 0.0
    end if

    ! Reset ice fluxes.
    if(include_ice) then
      flux%ICE = 0.0
    end if

  end subroutine reset_flux_type

  subroutine deallocate_var_type(var)

    implicit none

    type(var_type), intent(inout) :: var
    integer :: allocstat

    ! Deallocate basic variables.
    deallocate(var%rho, stat = allocstat)
    if(allocstat /= 0) stop "Deallocation of var_type component rho failed!"
    deallocate(var%u, stat = allocstat)
    if(allocstat /= 0) stop "Deallocation of var_type component u failed!"
    deallocate(var%v, stat = allocstat)
    if(allocstat /= 0) stop "Deallocation of var_type component v failed!"
    deallocate(var%w, stat = allocstat)
    if(allocstat /= 0) stop "Deallocation of var_type component w failed!"
    deallocate(var%pi, stat = allocstat)
    if(allocstat /= 0) stop "Deallocation of var_type component pi failed!"
    deallocate(var%rhop, stat = allocstat)
    if(allocstat /= 0) stop "Deallocation of var_type component rhop failed!"

    ! Deallocate dynamic Smagorinsky coefficient.
    if(turbScheme) then
      deallocate(var%DSC, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of var_type component DSC failed!"
    end if

    ! Deallocate gravity-wave heating.
    if(rayTracer) then
      deallocate(var%GWH, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of var_type component GWH failed!"
    end if

    ! Deallocate mass-weighted potential temperature.
    if(model == "compressible") then
      deallocate(var%P, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of var_type component P failed!"
    end if

    ! Deallocate tracer.
    if(include_tracer) then
      deallocate(var%chi, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of var_type component chi failed!"
    end if

    ! Deallocate ice variables.
    if(include_ice) then
      deallocate(var%ICE, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of var_type component ICE failed!"
    end if

    ! Deallocate optional variables.
    if(include_testoutput) then
      deallocate(var%OPT, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of var_type component OPT failed!"
    end if

  end subroutine deallocate_var_type

  subroutine deallocate_flux_type(flux)

    implicit none

    type(flux_type), intent(inout) :: flux
    integer :: allocstat

    ! Deallocate basic fluxes.
    deallocate(flux%rho, stat = allocstat)
    if(allocstat /= 0) stop "Deallocation of flux_type component rho failed!"
    deallocate(flux%u, stat = allocstat)
    if(allocstat /= 0) stop "Deallocation of flux_type component u failed!"
    deallocate(flux%v, stat = allocstat)
    if(allocstat /= 0) stop "Deallocation of flux_type component v failed!"
    deallocate(flux%w, stat = allocstat)
    if(allocstat /= 0) stop "Deallocation of flux_type component w failed!"
    deallocate(flux%theta, stat = allocstat)
    if(allocstat /= 0) stop "Deallocation of flux_type component theta failed!"
    deallocate(flux%rhop, stat = allocstat)
    if(allocstat /= 0) stop "Deallocation of flux_type component rhop failed!"

    ! Deallocate mass-weighted-potential-temperature flux.
    if(model == "compressible") then
      deallocate(flux%P, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of flux_type component P failed!"
    end if

    ! Deallocate tracer flux.
    if(include_tracer) then
      deallocate(flux%chi, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of flux_type component chi failed!"
    end if

    ! Deallocate ice fluxes.
    if(include_ice) then
      deallocate(flux%ICE, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of flux_type component ICE failed!"
    end if

  end subroutine deallocate_flux_type

end module type_module
