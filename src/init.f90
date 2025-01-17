module init_module

  use type_module
  use atmosphere_module
  use ice_module
  use mpi_module
  use boundary_module
  use sizeof_module
  use mpi
  use output_netCDF_module

  implicit none

  private ! private module

  public :: initialise
  public :: setup

  private :: cphase

  contains

  subroutine setup(var, var0, var1, varG, flux, flux0, force, source, dRho, &
      &dRhop, dMom, dTheta, dPStrat, drhoStrat, w_0, dIce, dTracer, &
      &tracerforce, dPot)

    !-----------------------------------------
    ! allocate var and flux / read the namelist
    !-----------------------------------------

    ! in/out variables
    type(var_type), intent(out) :: var, var0, var1, varG, source
    type(flux_type), intent(out) :: flux, flux0
    real, dimension(:, :, :, :), allocatable, intent(out) :: force
    type(tracerForceType), dimension(:, :, :), allocatable, intent(out) :: &
        &tracerforce ! tracer forcing in WKB
    real, dimension(:, :, :), allocatable :: dRho, dRhop ! RK-Update for rho
    real, dimension(:, :, :, :), allocatable :: dMom ! ...rhoU,rhoV,rhoW
    real, dimension(:, :, :), allocatable :: dTheta ! RK-Update for theta
    real, dimension(:, :, :, :), allocatable :: dIce ! RK-Update for nIce,qIce,qAer,qv
    real, dimension(:, :, :), allocatable :: dTracer ! RK-Update for rhoTracer
    real, dimension(:, :, :), allocatable :: dPot ! RK-Update for P
    real, dimension(:), allocatable :: dPStrat, drhoStrat ! RK-Update for P
    real, dimension(:), allocatable :: w_0 !! w_0 from ONK14

    integer :: allocstat
    integer ::i, j, k, iVar

    ! Set constants.
    pi = 4 * atan(1.0)

    ! Set default values.
    call default_values

    !-------------------------------------
    !         Read namelists
    !-------------------------------------

    ! Open input file.
    open(unit = 10, file = file_namelist, action = "read", form = "formatted", &
        &status = "old")

    ! Read variable namelist.
    rewind(unit = 10)
    read(unit = 10, nml = variables, end = 1)
    1 continue

    ! Read output namelist.
    rewind(unit = 10)
    read(unit = 10, nml = outputList, end = 2)
    2 continue

    ! Read debugging namelist.
    rewind(unit = 10)
    read(unit = 10, nml = debuggingList, end = 3)
    3 continue

    ! Read test case namelist.
    rewind(unit = 10)
    read(unit = 10, nml = testCaseList, end = 4)
    4 continue

    ! Read monochromatic wave namelist.
    rewind(unit = 10)
    read(unit = 10, nml = monochromeWave, end = 5)
    5 continue

    ! Read wave packet namelist.
    rewind(unit = 10)
    read(unit = 10, nml = wavePacket, end = 6)
    6 continue

    ! Read Lagrange ray tracing namelist.
    rewind(unit = 10)
    read(unit = 10, nml = LagrangeRayTracing, end = 7)
    7 continue

    ! Read bubble namelist.
    rewind(unit = 10)
    read(unit = 10, nml = bubble, end = 8)
    8 continue

    ! Read Robert bubble namelist.
    rewind(unit = 10)
    read(unit = 10, nml = robert_bubble, end = 9)
    9 continue

    ! Read mountain wave namelist.
    rewind(unit = 10)
    read(unit = 10, nml = mountainwavelist, end = 10)
    10 continue

    ! Read baroclinic LC namelist.
    rewind(unit = 10)
    read(unit = 10, nml = baroclinic_LC, end = 11)
    11 continue

    ! Read baroclinic ID namelist.
    rewind(unit = 10)
    read(unit = 10, nml = baroclinic_ID, end = 12)
    12 continue

    ! Read model equations namelist.
    rewind(unit = 10)
    read(unit = 10, nml = modelList, end = 13)
    13 continue

    ! Read solver namelist.
    rewind(unit = 10)
    read(unit = 10, nml = solverList, end = 14)
    14 continue

    ! Read Poisson solver namelist.
    rewind(unit = 10)
    read(unit = 10, nml = poissonSolverList, end = 15)
    15 continue

    ! Read atmosphere namelist.
    rewind(unit = 10)
    read(unit = 10, nml = atmosphereList, end = 16)
    16 continue

    ! Read topography namelist.
    rewind(unit = 10)
    read(unit = 10, nml = topographyList, end = 17)
    17 continue

    ! Read boundary namelist.
    rewind(unit = 10)
    read(unit = 10, nml = boundaryList, end = 18)
    18 continue

    ! Read second boundary namelist.
    rewind(unit = 10)
    read(unit = 10, nml = boundaryList2, end = 19)
    19 continue

    ! Read WKB namelist.
    rewind(unit = 10)
    read(unit = 10, nml = wkbList, end = 20)
    20 continue

    ! Read tracer namelist.
    rewind(unit = 10)
    read(unit = 10, nml = tracerList, end = 21)
    21 continue

    ! Read ice namelist.
    rewind(unit = 10)
    read(unit = 10, nml = iceList, end = 22)
    22 continue

    ! Close input file.
    close(unit = 10)

    !---------------------------------------------------
    ! allocate x,y,z - cell centered coordinate fields
    !---------------------------------------------------

    allocate(x(- nbx:sizeX + nbx), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate x"

    allocate(y(- nby:sizeY + nby), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate y"

    allocate(z(- nbz:sizeZ + nbz), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate z"

    !---------------------------------------------------
    ! allocate surface and fields for immersed boundary
    !---------------------------------------------------

    if(rayTracer .and. case_wkb == 3 .and. .not. topography) stop "Error: &
        &topography must be .true. for rayTracer == .true. and case_wkb == 3!"

    ! Allocate resolved topography.
    if(topography) then
      allocate(topography_surface(- nbx:nx + nbx, - nby:ny + nby), stat &
          &= allocstat)
      if(allocstat /= 0) stop "setup: could not allocate topography_surface"
      if(topographyTime > 0.0) then
        allocate(final_topography_surface(- nbx:nx + nbx, - nby:ny + nby), &
            &stat = allocstat)
        if(allocstat /= 0) stop "setup: could not allocate &
            &final_topography_surface"
      end if
    end if

    ! Allocate unresolved topography.
    if(rayTracer .and. case_wkb == 3) then
      allocate(k_spectrum(1:nx, 1:ny, 1:nwm), stat = allocstat)
      if(allocstat /= 0) stop "setup: could not allocate k_spectrum"
      allocate(l_spectrum(1:nx, 1:ny, 1:nwm), stat = allocstat)
      if(allocstat /= 0) stop "setup: could not allocate l_spectrum"
      allocate(topography_spectrum(1:nx, 1:ny, 1:nwm), stat = allocstat)
      if(allocstat /= 0) stop "setup: could not allocate topography_spectrum"
      if(topographyTime > 0.0) then
        allocate(final_topography_spectrum(1:nx, 1:ny, 1:nwm), stat = allocstat)
        if(allocstat /= 0) stop "setup: could not allocate &
            &final_topography_spectrum"
      end if
    end if

    allocate(kbl_topo(- nbx:nx + nbx, - nby:ny + nby, 3), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate kbl_topo"
    allocate(dhdx(- nbx:nx + nbx, - nby:ny + nby, 3), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate dhdx"
    allocate(dhdy(- nbx:nx + nbx, - nby:ny + nby, 3), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate dhdy"
    allocate(x_ip(- nbx:nx + nbx, - nby:ny + nby, 3), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate x_ip"
    allocate(y_ip(- nbx:nx + nbx, - nby:ny + nby, 3), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate y_ip"
    allocate(z_ip(- nbx:nx + nbx, - nby:ny + nby, 3), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate z_ip"
    allocate(velocity_reconst_t(- nbx:nx + nbx, - nby:ny + nby, 3), stat &
        &= allocstat)
    if(allocstat /= 0) then
      stop "init.f90: could not allocate velocity_reconst_t"
    end if
    allocate(velocity_reconst_n(- nbx:nx + nbx, - nby:ny + nby, 3), stat &
        &= allocstat)
    if(allocstat /= 0) then
      stop "init.f90: could not allocate velocity_reconst_n"
    end if

    !-------------------------------------
    !      allocate variable fields
    !-------------------------------------

    if(include_ice) then
      nVarIce = 3
      inN = 1
      inQ = 2
      inQv = 3

      allocate(iVarIce(nVarIce), stat = allocstat)
      if(allocstat /= 0) stop "Allocation of iVarIce failed!"

      iVarIce = [inN, inQ, inQv]
    end if

    call allocate_var_type(var)
    call reset_var_type(var)

    call allocate_var_type(var0)
    call reset_var_type(var0)

    call allocate_var_type(var1)
    call reset_var_type(var1)

    call allocate_var_type(varG)
    call reset_var_type(varG)

    call allocate_var_type(source)
    call reset_var_type(source)

    ! allocate dRho
    allocate(dRho(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dRho."

    ! allocate dRhop
    allocate(dRhop(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dRhop."

    ! allocate dMom
    allocate(dMom(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 3), stat &
        &= allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dMom."

    ! allocate dTheta
    allocate(dTheta(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dTheta."

    ! allocate dPot
    if(model == "compressible") then
      allocate(dPot(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "init.f90: Could not allocate dPot."
    end if

    if(include_ice) then
      allocate(dIce(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 3), stat &
          &= allocstat)
      if(allocstat /= 0) stop "init.f90: Could not allocate dIce."
    end if !include_ice

    if(include_tracer) then
      allocate(dTracer(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "init.f90: Could not allocate dTracer."
    end if

    ! allocate dPStrat
    allocate(dPStrat(- nbz:nz + nbz), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dPStrat."

    ! allocate drhoStrat
    allocate(drhoStrat(- nbz:nz + nbz), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate drhoStrat."

    ! allocate w_0
    allocate(w_0(- nbz:nz + nbz), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate w_0"

    call allocate_flux_type(flux)
    call reset_flux_type(flux)

    call allocate_flux_type(flux0)
    call reset_flux_type(flux0)

    ! allocate force
    allocate(force(0:nx + 1, 0:ny + 1, 0:nz + 1, 3), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate force"

    ! allocate tracerforce
    allocate(tracerforce(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate tracerforce"

    allocate(p_env_pp(1:nx, 1:ny, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate p_env_pp"

    allocate(the_env_pp(1:nx, 1:ny, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate the_env_pp"

    allocate(dens_env_pp(1:nx, 1:ny, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate dens_env_pp"

    allocate(u_env_pp(0:nx, 1:ny, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate u_env"

    allocate(v_env_pp(1:nx, 0:ny, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate v_env"

    allocate(u_const(0:nx, 1:ny), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate u_const"

    ! GBcorr: decide if heating is applied or not + error messages
    if(heatingONK14 .and. (background == 'realistic' .or. background &
        &== 'isothermal' .or. background == 'isentropic' .or. background &
        &== 'const-N')) then
      print *, "WARNING: you are using an idealized background &
          &(realistic/isothermal/isentropic/const-N) AND heating is ON &
          &(heatingONK14=.T.)!!! To switch off heating set heatingONK14=.F. in &
          &the namelist."
      heating = (heatingONK14 .and. (TurbScheme .or. rayTracer))
    else
      heating = (heatingONK14 .and. (TurbScheme .or. rayTracer))
    endif

    if(TestCase == "baroclinic_LC") then
      call allocate_var_type(var_env)
      call reset_var_type(var_env)

      if(background == "HeldSuarez") then
        allocate(kt_hs(1:ny, 0:nz + 1), stat = allocstat)
        if(allocstat /= 0) stop "init.f90: could not allocate kt_hs"

        allocate(kv_hs(0:ny + 1, 0:nz + 1), stat = allocstat)
        if(allocstat /= 0) stop "init.f90: could not allocate kv_hs"

        allocate(kw_hs(0:nz + 1), stat = allocstat)
        if(allocstat /= 0) stop "init.f90: could not allocate kw_hs"

        if(topography) then
          allocate(kt_hs_tfc(1:nx, 1:ny, 1:nz), stat = allocstat)
          if(allocstat /= 0) stop "init.f90: could not allocate kt_hs_tfc"
        end if
      end if
    end if

    allocate(kr_sp(0:ny + 1, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate kr_sp"

    allocate(kr_sp_w(0:ny + 1, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate kr_sp_w"

    if(topography) then
      allocate(kr_sp_tfc(0:(nx + 1), 0:(ny + 1), 0:(nz + 1)), stat = allocstat)
      if(allocstat /= 0) stop "init.f90: could not allocate kr_sp_tfc"

      allocate(kr_sp_w_tfc(0:(nx + 1), 0:(ny + 1), 0:(nz + 1)), stat &
          &= allocstat)
      if(allocstat /= 0) stop "init.f90: could not allocate kr_sp_w_tfc"
    end if

    ! Allocate damping coefficient for unified sponge.
    if(spongeLayer .and. unifiedSponge) then
      allocate(alphaUnifiedSponge(0:(nx + 1), 0:(ny + 1), 0:(nz + 1)), stat &
          &= allocstat)
      if(allocstat /= 0) stop "init.f90: could not allocate alphaUnifiedSponge"
    end if

    ! Safety switch for general setup.
    if(topography) then
      if(model == "WKB" .or. poissonSolverType /= "bicgstab" .or. reconstType &
          &/= "MUSCL" .or. musclType /= "muscl1" .or. fluxType /= "upwind" &
          &.or. heatingONK14 .or. pressureScaling .or. dens_relax) then
        stop "Terrain-following coordinates not implemented for chosen setup!"
      end if
    end if

    ! Safety switch for test cases.
    if(topography) then
      if(.not. any(testCase == (/character(50)::"sinus", "uniform_theta", &
          &"monochromeWave", "wavePacket", "mountainwave", "raytracer", &
          &"Robert_Bubble", "coldBubble", "smoothVortex", "atmosphereatrest", &
          &"SkamarockKlemp94", "hotBubble", "hotBubble2D", "hotBubble3D"/))) &
          &then
        stop "Terrain-following coordinates not implemented for chosen testCase"
      end if
    end if

    ! Safety switch for non-TFC setup.
    if(.not. topography) then
      if(testTFC) stop "Set testTFC = .false."
      if(freeSlipTFC) stop "Set freeSlipTFC = .false."
    end if

    ! Safety switch for unified sponge.
    if(spongeLayer .and. unifiedSponge) then
      if(testCase == "baroclinic_LC" .or. testCase == "baroclinic_ID") then
        stop "Unified sponge not ready for chosen setup!"
      end if
    end if

    ! Safety switch for halos in TFC
    if(topography) then
      if(.not. (nbx >= 3 .and. nby >= 3 .and. nbz >= 3)) then
        stop "Three halos / ghost cells are needed in TFC!"
      end if
    end if

    !---------------------------------------
    !        Model equation settings
    !---------------------------------------

    select case(model)

    case("anelastic")

      pressureScaling = .false.

    case("Boussinesq")

      ! Boussinesq model with bicgstab
      ! Bouyancy is predicted via the auxiliary equation. Density
      ! fluctuations are therefore stored in var%rhop(i, j, k), whereas
      ! var%rho(i, j, k) must remain zero!
      updateMass = .true. ! necessary because of auxiliary equation
      predictMomentum = .true.
      correctMomentum = .true.
      updateTheta = .false.
      if(include_ice) then
        updateIce = .true.
      end if
      auxil_equ = .true.
      poissonSolverType = "bicgstab"
      reconstType = "MUSCL"
      musclType = "muscl1"
      fluxType = "upwind"
      heatingONK14 = .false.
      TurbScheme = .false.
      rayTracer = .false.
      pressureScaling = .false.

      updateTracer = .true.

    case("pseudo_incompressible")

      updateMass = .true.
      predictMomentum = .true.
      correctMomentum = .true.
      updateTheta = .false.
      if(include_ice) then
        updateIce = .true.
      end if
      if(include_tracer) then
        updateTracer = .true.
      end if

      !overwrite unsuitable input settings
      if(zBoundary == "periodic") then
        print *, "WARNING: zBoundary periodic not possible.  Reset to &
            &solid_wall!"
        zBoundary = "solid_wall"
      end if

    case("compressible")

      updateMass = .true.
      predictMomentum = .true.
      correctMomentum = .true.
      updateTheta = .false.

      !overwrite unsuitable input settings
      if(zBoundary == "periodic") then
        print *, "WARNING: zBoundary periodic not possible.  Reset to &
            &solid_wall!"
        zBoundary = "solid_wall"
      end if

      ! Allow compressible model only in TFC with semi-implicit time stepping
      if(.not. topography .or. timeScheme /= "semiimplicit") stop "TFC and &
          &semi-implicit time stepping needed in compressible model!"

    case default
      print *, "model = ", model
      stop "initialize: Unknown model"
    end select

    ! Write all namelists.
    call write_namelists

  end subroutine setup

  ! --------------------------------------------------------------------

  subroutine initialise(var, flux)
    implicit none

    !------------------
    ! setup test cases
    !------------------

    ! in/out variables
    type(var_type), intent(inout):: var
    type(flux_type), intent(inout) :: flux

    ! local variables
    integer :: i, j, k
    integer :: ivr

    integer :: i0, j0, k_test, k_1, k_2 ! modified by Junhong Wei (20161121)

    ! greshoVortex
    real :: xu, yu, zu, xv, yv, xw, zw, delX, delY
    real :: x0, y0, r, uVortex, rhodl, pcoeff, pcoeff_r1
    real :: uPhi, p, rho0
    real :: pInf, u0, r0

    ! densityDisc
    real :: v0, w0, dens, rhoDisc, z0

    ! Robert's hot and cold bubble
    real :: dTheta1, a1, sigma1, xCenter1, zCenter1, rhoCenter
    real :: dTheta2, a2, sigma2, xCenter2, zCenter2

    ! hotBubble
    real :: x_dim, y_dim, z_dim, delZ
    real :: dTheta, dTheta0, theta, rho
    real :: dTheta_dim ! , dTheta0_dim

    ! wavepacket: all quantities are scaled
    real :: kk, mm, kTot ! vertical, zonal, total wave number
    real :: ll ! meridional wave number
    real :: omi, omi2 ! intrinsic frequency, squared
    real :: bAmp, uAmp, wAmp, pAmp ! ampl. for buoyancy, u, w, Exner pr.
    real :: delx2, delz2 ! squared distance from center
    real :: envel ! envelope of wave packet
    real :: Gauss ! Gaussian distribution value
    real :: sigma_z ! vert. width of Gaussian distribution
    real :: sigma_x ! hor. width Gaussian distr. (x dir.)
    real :: sigma_y ! hor. width Gaussian distr. (y dir.)
    real :: L_cos ! width of Cosine distribution
    real :: xCenter, yCenter ! center of distribution (hor.)
    real :: zCenter ! center of distribution (vert.)
    real :: phi ! local phase
    real :: u, v, w, b, t ! buoyancy, velocities
    real :: lambdaX, lambdaY ! hor. wave lengths
    real :: lambdaZ ! vert. wave lengths

    complex, dimension(:, :, :, :, :), allocatable :: Psi

    integer :: allocstat

    real :: u1, w1, b1, p1, v1
    real :: u2, w2, b2, p2, v2

    logical, parameter :: initWave2 = .false.

    ! monochromatic wave
    real :: th, f, f2
    real :: lambda
    real :: amp, cot
    real :: uRot, vRot, wRot

    ! testing with Stefan
    real :: xx

    ! random noise on background
    real, allocatable, dimension(:, :, :) :: randNoise
    real, parameter :: randAmp = 0.0

    ! jet stream
    real :: u_jet ! jet stream velocity
    real :: L_jet
    real :: z0_jet

    ! baroclinic life cycle
    real :: term_B1, term_B2, term_A1, term_A2, term_P1, term_P2
    real :: alpha_t, alpha_s, z_diff
    real :: F_1, F_2, F_a, dFadz, z_1, z_2, smooth_diff, L_z
    real :: T_c, term_a, term_b, DgDy, dbdy, dT0dy, dTcdy, dZdy, dady
    real :: dpdy, dpsidy, T0_s, DgDy_tr, Ly_end, Lx_end, bar_sigma_x
    real :: tanhyP1, tanhyP2, ztdif, signmtanh_P1, signmtanh_P2
    real :: p_dim, p_dim1, tmp_dim, pressure, hor_wind, the_dim, rhe_dim, the, &
        &ampl_rho
    real, dimension(1:nx, 1:ny) :: P_ref, dPrefdy, Z_trop, T0_t
    real, dimension(1:nx, 1:ny) :: gamma_tr, gamma_st
    real, dimension(1:nx, 1:ny) :: dpdy_trop, P_trop
    real, dimension(0:nz + 1) :: F_0, dF0dz, rho_n
    real :: dthdz, streamfunc, cp, dexndy, dexndz, ex_pr_pert, dstdz
    real :: theta_bar_0, noise_mag
    real, dimension(1:nx, 1:ny, 1:nz) :: noise
    real, dimension(0:ny + 1, 0:nz + 1) :: pi_pr_yz
    real, dimension(0:nx + 1, 0:nz + 1) :: pi_pr_xz
    real, dimension(1:nx, 1:ny, 1:nz - 1) :: br_vais_sq, balance, balance1, &
        &balance2, balance3, balance4
    character(len = 20) :: noise_var, defth

    ! TFC FJ
    ! Baroclinic life cycles in TFC.
    real :: pEdgeR, pEdgeL, pFEdgeR, pFEdgeL, pEdgeF, pEdgeB, pREdgeF, pREdgeB
    real :: piUEdgeR, piDEdgeR, piUEdgeL, piDEdgeL, piFUEdgeR, piFDEdgeR, &
        &piFUEdgeL, piFDEdgeL, piUEdgeF, piDEdgeF, piUEdgeB, piDEdgeB, &
        &piRUEdgeF, piRDEdgeF, piRUEdgeB, piRDEdgeB
    real :: piGrad
    integer, dimension(1:nx, 1:ny) :: k_2_tfc
    real, dimension(1:nx, 1:ny) :: z_2_tfc, theta_bar_0_tfc, alpha_t_tfc
    real, dimension(1:nx, 1:ny, 1:nz) :: F_0_tfc, pi_pr_xz_tfc, pi_pr_yz_tfc

    integer :: i00, j00 ! modified by Junhong Wei (20161121)

    real :: rho_int_m0, rho_int_00, rho_int_mp, rho_int_0p, rho_int_0m, &
        &rho_int_pm, rho_int_p0

    real :: ymin, ymax, z_trpp0, deltht, thet0, ntrp, nstr, jwdth
    real :: z_trpp, fstrpp, fstht, xloc, yloc, yjet0, yjet, zloc, z_baro
    real :: zeta_baro, buoy0, buoy1, rho1, rhop0, rhop1
    real :: rptb, thtptb, rhotot

    real :: ptptb_x, ptptb_y, ptptb_z, ptptb_dh, ptptb_dz, ptptb_amp

    real, dimension(1:ny) :: s2_strtd, c2_strtd, c4_strtd
    real :: yjets, yjetn, dy_hs, tempev, pistar, thetastar, sig_pr, facsig
    real :: ka_hs, ks_hs, kf_hs
    real, dimension(1:nx, 1:ny) :: tp_sponge, tp_sp1, tp_sp2 !FS
    real, dimension(0:ny + 1) :: f_Coriolis_y !FS
    real :: spongeDz

    integer :: switch, k_tropopause

    real :: indwindcoeff
    integer :: iwm

    integer :: kshalf
    real :: fchtms, zhtmsd, zhtmsu
    real :: bmax
    bmax = 0.

    ! open the namelist file
    open(unit = 10, file = file_namelist, action = "read", form = "formatted", &
        &status = "old")

    !---------------------------
    !      General settings
    !---------------------------

    select case(model)

    case("Boussinesq")

      ! var%rho must remain zero! Background fields are constant!
      var%rho(:, :, :) = 0.0

      vertical = (/0.0, 0.0, 1.0/)

    case("pseudo_incompressible", "compressible")

      ! set vertical always paralle to z-axis
      vertical = (/0.0, 0.0, 1.0/)

    case("WKB")
      !

    case default
      stop "initialize: unknown case model."
    end select

    !-----------------------
    !      MPI stuff
    !-----------------------
    i0 = is + nbx - 1 ! 0 index, replace i -> i + i0 in x and y fields
    j0 = js + nby - 1

    ! on default there is no initial ice, humidity or aerosols in the atmosphere
    if(include_ice) var%ICE(:, :, :, 1:nVarIce) = 0.0
    ! just for safety reasons
    if(include_tracer) var%chi(:, :, :) = 0.0
    !---------------------------------------------------------------

    select case(testCase)

      !-------------------------------------
      !             Boussinesq only
      !-------------------------------------

    case("sinus")

      ! TFC FJ
      ! Density fluctuations are stored in var(i, j, k, 6),
      ! var(i, j, k, 1) must remain zero!
      var%rho(i, j, k) = 0.0
      var%u(:, :, :) = 0.0
      var%v(:, :, :) = 0.0
      var%w(:, :, :) = 0.0

      do i = - 2, nx + 3
        xx = lx(0) + (real(i - 1) + 0.5) * dx
        theta = sin(pi * xx * lRef)
        ! Density fluctuations are stored in var(i, j, k, 6)!
        var%rhop(i, 1, 1) = - theta * rho00 / theta00
      end do

    case("uniform_theta")

      ! zero atmospheric background flow
      var%u(:, :, :) = backgroundFlow_dim(1) / uRef
      var%v(:, :, :) = backgroundFlow_dim(2) / uRef
      var%w(:, :, :) = backgroundFlow_dim(3) / uRef

      ! constant pressure variable pi'
      var%pi(:, :, :) = 0.0

      ! uniform pot temp perturbation
      ! TFC FJ
      ! Density fluctuations are stored in var(i, j, k, 6)!
      var%rhop(:, :, :) = 0.0
      ! var(:,:,:,6) = 0.0 / thetaRef

    case("monochromeWave")

      ! read test case input data
      rewind(unit = 10)
      read(unit = 10, nml = monochromeWave)

      ! gather quantities
      th = vert_theta * pi / 180.0
      f = f_Coriolis_dim * tRef
      f2 = f ** 2
      omi = - sqrt(N2 * cos(th) ** 2 + f2 * sin(th) ** 2)

      ! wave vector
      lambda = lambda_dim / lRef
      kTot = 2.0 * pi / lambda

      kk = kTot * cos(th)
      mm = kTot * sin(th)

      ! amplitudes in non-rotated frame
      amp = amplitudeFactor * (- omi / kk)
      cot = 1.0 / tan(th)

      do k = 1, nz
        j = 1
        do i = 1, nx

          ! phase
          if(topography) then
            ! TFC FJ
            phi = - kk * x(i) - mm * heightTFC(i, j, k)
          else
            phi = - kk * x(i) - mm * z(k)
          end if

          ! amplitudes
          u = amp * cos(phi)
          v = amp * f / omi * sin(phi)
          w = amp * cot * cos(phi)
          b = amp * N2 / omi * cot * sin(phi)

          ! project on rotated coordinate system

          wRot = dot_product((/u, v, w/), vertical)
          vRot = v
          uRot = sqrt(u ** 2 + w ** 2 - wRot ** 2)
          theta = Fr2 * theta00 * b

          ! assign to var field
          var%u(i, j, k) = uRot
          var%v(i, j, k) = vRot
          var%w(i, j, k) = wRot
          ! TFC FJ
          ! Density fluctuations are stored in var(i, j, k, 6)!
          var%rhop(i, j, k) = - theta * rho00 / theta00
          ! var(i,j,k,6) = theta

          ! TFC FJ
          ! Compute terrain-following vertical wind.
          if(topography) then
            var%w(i, j, k) = var%w(i, j, k) / jac(i, j, k) + met(i, j, k, 1, &
                &3) * var%u(i, j, k) + met(i, j, k, 2, 3) * var%v(i, j, k)
          end if

        end do
      end do

      !-----------------------------
      !  Interpolate to cell faces
      !-----------------------------

      ! copy u values to ghost cells
      var%u(0, :, :) = var%u(nx, :, :)
      var%u(nx + 1, :, :) = var%u(1, :, :)

      ! average zonal velocities to cell face...
      do i = 0, nx
        var%u(i, :, :) = 0.5 * (var%u(i, :, :) + var%u(i + 1, :, :))
      end do

      ! copy v values to ghost cells
      var%v(:, 0, :) = var%v(:, ny, :)
      var%v(:, ny + 1, :) = var%v(:, 1, :)

      ! average meridional velocities to cell face...
      do j = 0, ny
        var%v(:, j, :) = 0.5 * (var%v(:, j, :) + var%v(:, j + 1, :))
      end do

      ! copy w values to ghost cells
      var%w(:, :, 0) = var%w(:, :, nz)
      var%w(:, :, nz + 1) = var%w(:, :, 1)

      ! average vertical velocities to cell face...
      do k = 0, nz
        var%w(:, :, k) = 0.5 * (var%w(:, :, k) + var%w(:, :, k + 1))
      end do

      !----------------------------------------------------------
      !          Gravity Waves: wave resolving simulations or WKB
      !----------------------------------------------------------

    case('wavePacket') ! 1D/2D wave packet

      !SD
      allocate(Psi(0:nx + 1, 0:ny + 1, 0:nz + 1, 5, 0:2), stat = allocstat)
      if(allocstat /= 0) stop "init.f90: Could not allocate Psi."
      allocate(randNoise(0:nx + 1, 0:ny + 1, 0:nz + 1), stat = allocstat)
      if(allocstat /= 0) stop "init.f90: Could not allocate randNoise."

      !---------------------
      ! set up random noise
      !---------------------

      call random_number(randNoise)
      do k = 1, nz
        randNoise(:, :, k) = randNoise(:, :, k) * randAmp
      end do

      !--------------------
      ! set up jet stream
      !--------------------

      ! read test case input data
      rewind(unit = 10)
      read(unit = 10, nml = wavePacket)

      !SDJul2024
      !MultipleWavePackets = .true.
      !if ( master ) print*, 'set MultipleWavePackets = .true.'

      if(MultipleWavePackets) then

        !init
        var%rho = 0.
        var%u = 0.
        var%v = 0.
        var%w = 0.
        var%pi = 0.
        TWM = NWM_WP

        read(unit = 10, nml = MultipleWavePackets_list)

        if(inducedwind .or. initWave2) then
          print *, 'MultipleWavePackets not implemented for &
              &inducedwind/initWave2 == true'
          stop
        end if
      else
        TWM = 1
      end if

      u0 = u0_jet_dim / uRef ! amplitude of jet
      L_jet = L_jet_dim / lRef ! half width of cos profile
      z0_jet = z0_jet_dim / lRef ! center of jet

      if(topography) then
        ! TFC FJ
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              delz = (heightTFC(i, j, k) - z0_jet)

              ! Cosine
              if(abs(delz) .le. L_jet) then
                u_jet = 0.5 * u0 * (1.0 + cos(pi * delz / L_jet))
              else
                u_jet = 0.0
              end if

              var%u(i, j, k) = u_jet
            end do
          end do
        end do
      else
        do k = 1, nz
          delz = (z(k) - z0_jet)

          ! Cosine
          if(abs(delz) .le. L_jet) then
            u_jet = 0.5 * u0 * (1.0 + cos(pi * delz / L_jet))
          else
            u_jet = 0.0
          end if

          var%u(:, :, k) = u_jet
        end do
      end if

      do iwm = 1, TWM

        ! overwite input if superposition of wavepackets considered
        if(MultipleWavePackets) then

          lambdaX_dim = lambdaX_dim_sp(iwm)
          lambdaY_dim = lambdaY_dim_sp(iwm)
          lambdaZ_dim = lambdaZ_dim_sp(iwm)

          x0_dim = x0_dim_sp(iwm)
          y0_dim = y0_dim_sp(iwm)
          z0_dim = z0_dim_sp(iwm)

          sigma_hor_dim = sigma_hor_dim_sp(iwm)
          sigma_hor_yyy_dim = sigma_hor_yyy_dim_sp(iwm)
          sigma_dim = sigma_dim_sp(iwm)

          amplitudeFactor = amplitudeFactor_sp(iwm)
          omiSign = omiSign_sp(iwm)
        end if

        !--------------------
        !     set up GWP
        !--------------------
        call init_GWP(Psi, kk, mm, ll, indwindcoeff)

        do k = 0, (nz + 1)
          do j = 0, (ny + 1)
            do i = 0, (nx + 1)
              if(topography) then
                ! TFC FJ
                phi = kk * x(i + i0) + ll * y(j + j0) + mm * heightTFC(i, j, k)
              else
                phi = kk * x(i + i0) + mm * z(k) + ll * y(j + j0)
              end if

              ! wave 1
              u1 = real(Psi(i, j, k, 1, 1) * exp(phi * imag))
              w1 = real(Psi(i, j, k, 2, 1) * exp(phi * imag))
              b1 = real(Psi(i, j, k, 3, 1) * exp(phi * imag))
              p1 = real(Psi(i, j, k, 4, 1) * exp(phi * imag))
              v1 = real(Psi(i, j, k, 5, 1) * exp(phi * imag))

              ! wave 2
              if(initWave2) then
                stop 'ERROR: 2ndary wave not ready for 2D or 3D wave p.'
                u2 = real(Psi(i, j, k, 1, 2) * exp(2. * phi * imag))
                w2 = real(Psi(i, j, k, 2, 2) * exp(2. * phi * imag))
                b2 = real(Psi(i, j, k, 3, 2) * exp(2. * phi * imag))
                p2 = real(Psi(i, j, k, 4, 2) * exp(2. * phi * imag))
                v2 = real(Psi(i, j, k, 5, 2) * exp(2. * phi * imag))
              end if

              ! sum of wave 1 and 2
              if(initWave2) then
                stop 'ERROR: 2ndary wave not ready for 2D or 3D wave p.'
                b = b1 + b2
                u = u1 + u2
                w = w1 + w2
                p = p1 + p2
                v = v1 + v2
              else
                b = b1
                u = u1
                w = w1
                p = p1
                v = v1
              end if

              !SDDec24
              var%rho(i, j, k) = var%rho(i, j, k) + b ! store b at 1
              var%u(i, j, k) = var%u(i, j, k) + u
              var%v(i, j, k) = var%v(i, j, k) + v
              var%w(i, j, k) = var%w(i, j, k) + w
              var%pi(i, j, k) = var%pi(i, j, k) + p

              if(iwm .eq. TWM) then

                !reset rho
                b = var%rho(i, j, k)
                var%rho(i, j, k) = 0.

                ! additional vars
                if(topography) then
                  ! TFC FJ
                  rho = 1.0 / (1.0 + Fr2 * b) * rhoStratTFC(i, j, k)
                else
                  rho = 1. / (1. + Fr2 * b) * rhoStrat(k)
                end if
                theta = Fr2 * theta00 * b

                ! write to field
                select case(model)
                case("pseudo_incompressible", "compressible")

                  ! add random noise
                  rho = rho + randNoise(i, j, k)

                  ! subtract background for fluctuation mode
                  if(topography) then
                    ! TFC FJ
                    rho = rho - rhoStratTFC(i, j, k)
                  else
                    rho = rho - rhoStrat(k)
                  end if

                  ! write to field
                  var%rho(i, j, k) = rho

                case("Boussinesq")

                  ! Density fluctuations are stored in var%rhop(i, j, k),
                  ! var%rho(i, j, k) must remain zero!
                  if(topography) then
                    var%rhop(i, j, k) = rho - rhoStratTFC(i, j, k)
                  else
                    var%rhop(i, j, k) = rho - rhoStrat(k)
                  end if

                  ! var(i,j,k,6) = theta

                case default
                  stop "initialize: unknown case model"
                end select ! model

                if(include_tracer) then
                  ! chi = <chi> + chi'
                  ! where chi' = alphaTracer/N^2 * b'
                  ! from inserting WKB ansatz into linearized
                  ! equation for chi' and using polarization
                  ! relation
                  if(topography) then
                    stop 'init.f90: wavepacket tracer prime and topography not &
                        &implemented'
                  else
                    ! only set up for <chi>=alphaTracer*z
                    ! large-scale tracer distribution
                    if(tracerSetup == "alpha_z") then
                      var%chi(i, j, k) = alphaTracer / N2 * b
                    else
                      stop 'init.f90: unknown initial tracer with wavepacket &
                          &tracer prime'
                    end if
                  end if
                end if

                if(inducedwind) then
                  !stop "Error: induced wind currently not possible. Potential error in code."
                  var%u(i, j, k) = var%u(i, j, k) + indwindcoeff * b ** 2.
                end if

                ! TFC FJ
                ! Compute terrain-following vertical wind.
                if(topography) then
                  var%w(i, j, k) = var%w(i, j, k) / jac(i, j, k) + met(i, j, &
                      &k, 1, 3) * var%u(i, j, k) + met(i, j, k, 2, 3) &
                      &* var%v(i, j, k)
                end if

                !CHANGES
                var%opt(i, j, k, 3) = (1. / (1. + Fr2 * Psi(i, j, k, 3, 1)) &
                    &* rhoStrat(k) - rhoStrat(k)) * rhoRef

              end if ! iwm==TWM
            end do !i
          end do ! l modified by Junhong Wei for 3DWP (20170922)
        end do ! k
      end do ! iwm

      ! average zonal velocities to cell face...
      do i = 0, nx
        var%u(i, :, :) = 0.5 * (var%u(i, :, :) + var%u(i + 1, :, :))
      end do

      ! average meridional velocities to cell face...
      do j = 0, ny
        var%v(:, j, :) = 0.5 * (var%v(:, j, :) + var%v(:, j + 1, :))
      end do

      ! average vertical velocities to cell faces
      do k = 0, nz ! modified by Junhong Wei for 3DWP (20171204)
        var%w(:, :, k) = 0.5 * (var%w(:, :, k) + var%w(:, :, k + 1))
      end do

      select case(model)
      case("pseudo_incompressible", "compressible")

        var%w(:, :, 0) = 0.0 ! reset velocity at wall to zero
        var%w(:, :, nz) = 0.0 ! reset velocity at wall to zero

      case("Boussinesq")

        ! TFC FJ
        if(zBoundary == "solid_wall") then
          var%w(:, :, 0) = 0.0
          var%w(:, :, nz) = 0.0
        end if

      case default
        stop "initialize: unknown case model"
      end select
      !SD
      if(include_ice) call setup_ice(var)

      !---------------------------------------------------------------

    case('mountainwave')
      ! for wave resolving simulation of mountain waves:
      ! read parameters for temporary wind relaxation
      ! zero-wind initial state for montain-wave simulations

      rewind(unit = 10)
      read(unit = 10, nml = mountainwavelist)

      ! nondimensionalization

      u_relax = u_relax / uRef
      v_relax = v_relax / uRef
      w_relax = w_relax / uRef

      t_relax = t_relax / tRef
      ! t_ramp = t_ramp / tRef

      xextent_relax = xextent_relax / lRef
      yextent_relax = yextent_relax / lRef

      if(wind_relaxation) then
        var%u(:, :, :) = u_relax / uRef
        var%v(:, :, :) = v_relax / uRef
        var%w(:, :, :) = w_relax / uRef
      else
        ! TFC FJ
        ! Provide initialization with constant background wind.
        var%u(:, :, :) = backgroundFlow_dim(1) / uRef
        var%v(:, :, :) = backgroundFlow_dim(2) / uRef
        var%w(:, :, :) = backgroundFlow_dim(3) / uRef

        ! FJMar2023
        ! Set background wind to zero in surface layer.
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              if(heightTFC(i, j, k) < surface_layer_depth / lRef) then
                var%u(i, j, k) = 0.0
                var%v(i, j, k) = 0.0
                var%w(i, j, k) = 0.0
              end if
            end do
          end do
        end do
      end if

      ! density, potential temperature, and pressure

      do k = 0, (nz + 1)
        do j = 0, (ny + 1)
          do i = 0, (nx + 1)
            select case(model)
            case("pseudo_incompressible", "compressible")
              ! initialization density = background density
              ! subtract background for fluctuation mode

              ! write to field
              var%rho(i, j, k) = 0.0

            case("Boussinesq")
              ! initialization zero buoyancy fluctuations

              var%rhop(i, j, k) = 0.0

            case default
              stop "initialize: unknown case model"
            end select

            ! initialization zero pressure fluctuations

            var%pi(i, j, k) = 0.0
          end do
        end do
      end do

      ! raytracer + case_wkb=5
      if(case_wkb == 5) then
        read(unit = 10, nml = case_wkb_5_list)
        if(nwm .ne. NWM_WP) then
          print *, 'nwm /= NWM_WP !!'
          print *, 'change NWM, NWM_WP'
          stop
        end if
      end if
      if(include_ice) call setup_ice(var)

      !-----------------------------------------------------------------

    case('raytracer')
      ! WKB simulations: Wave packet or mountain waves
      ! for the full set up see routine setup_wkb

      if(.not. raytracer) stop 'raytracer not set correctly'

      ! read namelist for wkb ray tracer
      rewind(unit = 10)
      read(unit = 10, nml = LagrangeRayTracing)

      ! Stop for unsuitable configuration.
      if(sizeX > 1 .and. fac_dk_init == 0.0) stop "Error in initialise: &
          &fac_dk_init = 0 and sizeX > 1!"
      if(sizeY > 1 .and. fac_dl_init == 0.0) stop "Error in initialise: &
          &fac_dl_init = 0 and sizeY > 1!"
      if(sizeZ == 1 .or. fac_dm_init == 0.0) stop "Error in initialise: &
          &fac_dm_init = 0 or sizeZ = 1!"
      if(wlrx_init == 0.0 .and. wlry_init == 0.0) stop "Error in initialise: &
          &wlrx_init = 0 and wlry_init = 0!"

      ! Compute initial ray-volume extent in k.
      if(sizeX == 1) then
        dk_init = 0.0
      else if(wlry_init == 0.0) then
        dk_init = fac_dk_init * 2.0 * pi / wlrx_init
      else if(wlrx_init == 0.0) then
        dk_init = fac_dk_init * 2.0 * pi / wlry_init
      else
        dk_init = fac_dk_init * 2.0 * pi * sqrt(1.0 / wlrx_init ** 2.0 + 1.0 &
            &/ wlry_init ** 2.0)
      end if

      ! Compute initial ray-volume extent in l.
      if(sizeY == 1) then
        dl_init = 0.0
      else if(wlrx_init == 0.0) then
        dl_init = fac_dl_init * 2.0 * pi / wlry_init
      else if(wlry_init == 0.0) then
        dl_init = fac_dl_init * 2.0 * pi / wlrx_init
      else
        dl_init = fac_dl_init * 2.0 * pi * sqrt(1.0 / wlrx_init ** 2.0 + 1.0 &
            &/ wlry_init ** 2.0)
      end if

      zmin_wkb = zmin_wkb_dim / lRef

      ! in WKB mountain-wave case read parameters for wind relaxation

      if(case_wkb == 3) then
        rewind(unit = 10)
        read(unit = 10, nml = mountainwavelist)

        if(zBoundary /= "solid_wall") stop "Error in initialise: zBoundary &
            &must be 'solid_wall' for case_wkb = 3!"

        ! nondimensionalization

        u_relax = u_relax / uRef
        v_relax = v_relax / uRef
        w_relax = w_relax / uRef

        t_relax = t_relax / tRef
        ! t_ramp = t_ramp / tRef

        xextent_relax = xextent_relax / lRef
        yextent_relax = yextent_relax / lRef

        ! increase relaxation wind u_relax so that u = u_relax after the
        ! relaxation period (in x-independent case without topography)

        if(wind_relaxation) then
          var%u(:, :, :) = u_relax / uRef
          var%v(:, :, :) = v_relax / uRef
          var%w(:, :, :) = w_relax / uRef
        else
          ! Provide initialization with constant background wind.
          var%u(:, :, :) = backgroundFlow_dim(1) / uRef
          var%v(:, :, :) = backgroundFlow_dim(2) / uRef
          var%w(:, :, :) = backgroundFlow_dim(3) / uRef

          ! FJMar2023
          ! Set background wind to zero in surface layer.
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                if(heightTFC(i, j, k) < surface_layer_depth / lRef) then
                  var%u(i, j, k) = 0.0
                  var%v(i, j, k) = 0.0
                  var%w(i, j, k) = 0.0
                end if
              end do
            end do
          end do
        end if
      else
        ! zero wind
        var%u(:, :, :) = 0.0 ! u
        var%v(:, :, :) = 0.0 ! v
        var%w(:, :, :) = 0.0 ! w
      end if

      ! density, potential temperature, and pressure

      do k = 0, (nz + 1)
        do j = 0, (ny + 1)
          do i = 0, (nx + 1)
            select case(model)
            case("pseudo_incompressible", "compressible")
              ! initialization density = background density
              ! subtract background for fluctuation mode

              ! write to field
              var%rho(i, j, k) = 0.0

            case("Boussinesq")
              ! initialization zero buoyancy fluctuations

              var%rhop(i, j, k) = 0.0

            case default
              stop "initialize: unknown case model"
            end select

            ! initialization zero pressure fluctuations

            var%pi(i, j, k) = 0.0
          end do
        end do
      end do

      !SD
      if(include_ice) call setup_ice(var)

      !---------------------------------------------------
      !                Hot and cold bubbles
      !---------------------------------------------------

    case('Robert_Bubble')
      ! read test case input data
      rewind(unit = 10)
      read(unit = 10, nml = robert_bubble)

      ! atmospheric background flow
      var%u(:, :, :) = backgroundFlow_dim(1) / uRef
      var%v(:, :, :) = backgroundFlow_dim(2) / uRef
      var%w(:, :, :) = backgroundFlow_dim(3) / uRef

      ! constant pressure variable pi'
      var%pi(:, :, :) = 0.0

      ! non-dimensionalize input parameters
      dTheta1 = dTheta1_dim / thetaRef
      a1 = a1_dim / lRef
      sigma1 = sigma1_dim / lRef
      xCenter1 = xCenter1_dim / lRef
      zCenter1 = zCenter1_dim / lRef

      dTheta2 = dTheta2_dim / thetaRef
      a2 = a2_dim / lRef
      sigma2 = sigma2_dim / lRef
      xCenter2 = xCenter2_dim / lRef
      zCenter2 = zCenter2_dim / lRef

      ! potential temperature and density
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx

            ! init theta off set
            dTheta = 0.0

            ! hot bubble
            delX = (x(i) - xCenter1)
            if(topography) then
              ! TFC FJ
              delz = heightTFC(i, j, k) - zCenter1
            else
              delZ = (z(k) - zCenter1)
            end if
            r = sqrt(delX ** 2 + delZ ** 2) ! scaled radius

            if(r <= a1) then
              dTheta = dTheta1
            else
              Gauss = exp(- (r - a1) ** 2 / sigma1 ** 2)
              dTheta = dTheta1 * Gauss
            end if

            ! cold bubble
            delX = (x(i) - xCenter2)
            if(topography) then
              ! TFC FJ
              delz = heightTFC(i, j, k) - zCenter2
            else
              delZ = (z(k) - zCenter2)
            end if
            r = sqrt(delX ** 2 + delZ ** 2) ! scaled radius

            if(r <= a2) then
              dTheta = dTheta + dTheta2
            else
              Gauss = exp(- (r - a2) ** 2 / sigma2 ** 2)
              dTheta = dTheta + dTheta2 * Gauss
            end if

            ! total potential temperature
            if(topography) then
              ! TFC FJ
              theta = thetaStratTFC(i, j, k) + dTheta
            else
              theta = thetaStrat(k) + dTheta
            end if

            select case(model)

            case("pseudo_incompressible", "compressible")

              ! calc pseudo-incompressible density rho*
              if(referenceQuantities == "SI") then
                rho = p0 ** kappa / Rsp * Pstrat(k) / theta
              else
                if(topography) then
                  ! TFC FJ
                  rho = pStratTFC(i, j, k) / theta
                else
                  rho = Pstrat(k) / theta
                end if
              end if

              ! subtract background for fluctuation mode
              if(topography) then
                ! TFC FJ
                rho = rho - rhoStratTFC(i, j, k)
              else
                rho = rho - rhoStrat(k)
              end if

              var%rho(i, j, k) = rho

            case("Boussinesq")

              ! set pot Temp deviation
              ! TFC FJ
              ! Density fluctuations are stored in var%rhop(i, j, k)!
              if(topography) then
                rho = pStratTFC(i, j, k) / theta - rhoStratTFC(i, j, k)
              else
                rho = pStrat(k) / theta - rhoStrat(k)
              end if
              var%rhop(i, j, k) = rho

            case default
              stop "initialize: unknown model."
            end select

          end do
        end do
      end do

      !------------------------------------------------------------------

    case('coldBubble')

      ! read test case input data
      rewind(unit = 10)
      read(unit = 10, nml = bubble)

      if(referenceQuantities == "SI") then
        stop "initialize: SI units not allowed"
      end if

      ! zero start velocity
      var%u(:, :, :) = 0.0
      var%v(:, :, :) = 0.0
      var%w(:, :, :) = 0.0

      ! constant pressure variable pi'
      var%pi(:, :, :) = 0.0

      i00 = is + nbx - 1 ! 0 index,

      ! potential temperature and density
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            x_dim = x(i + i00) * lRef ! dimensional lenghts
            if(topography) then
              ! TFC FJ
              z_dim = heightTFC(i, j, k) * lRef
            else
              z_dim = z(k) * lRef
            end if

            delX = (x_dim - xCenter_dim) / xRadius_dim
            delZ = (z_dim - zCenter_dim) / zRadius_dim

            r = sqrt(delX ** 2 + delZ ** 2) ! scaled radius

            if(r <= 1.0) then ! inside bubble
              dTheta_dim = - 7.5 * (1. + cos(pi * r))
              != 0.5*dTheta0_dim * (1.0 + (cos(pi*r/2.0))**2)
              if(topography) then
                ! TFC FJ
                theta = thetaStratTFC(i, j, k) + dTheta_dim / thetaRef
              else
                theta = thetaStrat(k) + dTheta_dim / thetaRef
              end if

              ! calc pseudo-incompressible density rho*
              if(referenceQuantities == "SI") then
                rho = p0 ** kappa / Rsp * Pstrat(k) / theta - rhoStrat(k)
              else
                if(topography) then
                  ! TFC FJ
                  rho = pStratTFC(i, j, k) / theta - rhoStratTFC(i, j, k)
                else
                  rho = Pstrat(k) / theta - rhoStrat(k)
                end if
              end if

              select case(model)

              case("pseudo_incompressible", "compressible")

                var%rho(i, j, k) = rho

              case("Boussinesq")

                ! TFC FJ
                ! Density fluctuations are stored in var%rhop(i, j, k),
                ! var%rho(i, j, k) must remain zero!
                var%rho(i, j, k) = 0.0
                var%rhop(i, j, k) = rho

              case default
                stop "initialize: unknown model."
              end select
            else ! outside bubble
              ! keep background density
              var%rho(i, j, k) = 0.
            end if

          end do
        end do
      end do

    case('smoothVortex')

      ! read test case input data
      rewind(unit = 10)
      read(unit = 10, nml = bubble)

      if(referenceQuantities == "SI") then
        stop "initialize: SI units not allowed"
      end if

      ! TFC FJ
      if(model == "Boussinesq") then
        stop "Boussinesq model not ready for smoothVortex!"
      end if

      i00 = is + nbx - 1 ! 0 index,
      j00 = js + nby - 1

      var%rho(:, :, :) = rhoCenter_dim
      var%u(:, :, :) = backgroundFlow_dim(1)
      var%v(:, :, :) = backgroundFlow_dim(2)
      var%w(:, :, :) = 0.

      ! potential temperature and density
      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            x_dim = x(i + i00) * lRef ! + 0.5*dx*lRef
            y_dim = y(j + j00) * lRef

            delX = (x_dim - xCenter_dim)
            delY = (y_dim - yCenter_dim)

            r = sqrt(delX ** 2 + delY ** 2) / 0.4

            th = atan(delY / delX)

            if(r < 1.0) then
              var%u(i, j, k) = - 1024. * delY / r * (1. - r) ** 6 * r ** 6 &
                  &+ var%u(i, j, k)
            end if

          end do
        end do
      end do

      do k = 1, nz
        do j = 0, ny
          do i = 1, nx
            x_dim = x(i + i00) * lRef
            y_dim = y(j + j00) * lRef !+ 0.5*dy*lRef

            delX = (x_dim - xCenter_dim)
            delY = (y_dim - yCenter_dim)

            r = sqrt(delX ** 2 + delY ** 2) / 0.4

            th = atan(delY / delX)

            if(r < 1.0) then
              var%v(i, j, k) = 1024. * delX / r * (1. - r) ** 6 * r ** 6 &
                  &+ var%v(i, j, k)
            end if
          end do
        end do
      end do

      var%pi(:, :, :) = 0.

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            x_dim = x(i + i00) * lRef
            y_dim = y(j + j00) * lRef

            delX = (x_dim - xCenter_dim)
            delY = (y_dim - yCenter_dim)

            r = sqrt(delX ** 2 + delY ** 2) / 0.4

            if(r < 1.0) then ! inside bubble
              var%rho(i, j, k) = var%rho(i, j, k) + 0.5 * (1. - r ** 2) ** 6

              rhodl = (1. + 0.5 * (1. - r ** 2) ** 6) / rhoCenter_dim
              pcoeff = 1024. ** 2 * r ** 12 * (1. / 72. * r ** 24 - 6. / 35. &
                  &* r ** 23 + 15. / 17. * r ** 22 - 74. / 33. * r ** 21 + 57. &
                  &/ 32. * r ** 20 + 174. / 31. * r ** 19 - 269. / 15. * r &
                  &** 18 + 450. / 29. * r ** 17 + 153. / 8. * r ** 16. - 1564. &
                  &/ 27. * r ** 15 + 510. / 13. * r ** 14 + 204. / 5. * r &
                  &** 13. - 1. / 24. * (2210. - rhodl) * r ** 12 + 12. / 23. &
                  &* (85. - rhodl) * r ** 11 + (510. / 11. + 3. * rhodl) * r &
                  &** 10 - 4. / 21. * (391. + 55. * rhodl) * r ** 9 + 9. / 40. &
                  &* (119. + 110. * rhodl) * r ** 8 + 18. / 19. * (25. - 44. &
                  &* rhodl) * r ** 7 - 1. / 9. * (269. - 462 * rhodl) * r ** 6 &
                  &+ 6. / 17. * (29. - 132 * rhodl) * r ** 5 + 3. / 16. * (19. &
                  &+ 165. * rhodl) * r ** 4 - 2. / 15. * (37. + 110. * rhodl) &
                  &* r ** 3 + 3. / 7. * (5. + 11. * rhodl) * r ** 2 - 6. / 13. &
                  &* (1. + 2. * rhodl) * r + 1. / 24. * (1. + 2. * rhodl))

              p_dim = pcoeff / pRef * rhoCenter_dim * ((backgroundFlow_dim(1)) &
                  &** 2 + (backgroundFlow_dim(2)) ** 2)

              rhodl = 1. / rhoCenter_dim
              pcoeff_r1 = 1024. ** 2 * (1. / 72. - 6. / 35. + 15. / 17. - 74. &
                  &/ 33. + 57. / 32. + 174. / 31. - 269. / 15. + 450. / 29. &
                  &+ 153. / 8. - 1564. / 27. + 510. / 13. + 204. / 5. - 1. &
                  &/ 24. * (2210. - rhodl) + 12. / 23. * (85. - rhodl) + (510. &
                  &/ 11. + 3. * rhodl) - 4. / 21. * (391. + 55. * rhodl) + 9. &
                  &/ 40. * (119. + 110. * rhodl) + 18. / 19. * (25. - 44. &
                  &* rhodl) - 1. / 9. * (269. - 462 * rhodl) + 6. / 17. * (29. &
                  &- 132 * rhodl) + 3. / 16. * (19. + 165. * rhodl) - 2. / 15. &
                  &* (37. + 110. * rhodl) + 3. / 7. * (5. + 11. * rhodl) - 6. &
                  &/ 13. * (1. + 2. * rhodl) + 1. / 24. * (1. + 2. * rhodl))

              p_dim1 = pcoeff_r1 / pRef * rhoCenter_dim &
                  &* ((backgroundFlow_dim(1)) ** 2 + (backgroundFlow_dim(2)) &
                  &** 2)

              var%pi(i, j, k) = kappaInv * (p_dim ** kappa - p_dim1 ** kappa)

            end if

            if(topography) then
              ! TFC FJ
              var%rho(i, j, k) = var%rho(i, j, k) / rhoRef - rhoStratTFC(i, j, &
                  &k)
            else
              var%rho(i, j, k) = var%rho(i, j, k) / rhoRef - rhoStrat(k)
            end if
            var%u(i, j, k) = var%u(i, j, k) / uRef
            var%v(i, j, k) = var%v(i, j, k) / uRef

            rhodl = 1. / rhoCenter_dim
            pcoeff_r1 = 1024. ** 2 * (1. / 72. - 6. / 35. + 15. / 17. - 74. &
                &/ 33. + 57. / 32. + 174. / 31. - 269. / 15. + 450. / 29. &
                &+ 153. / 8. - 1564. / 27. + 510. / 13. + 204. / 5. - 1. / 24. &
                &* (2210. - rhodl) + 12. / 23. * (85. - rhodl) + (510. / 11. &
                &+ 3. * rhodl) - 4. / 21. * (391. + 55. * rhodl) + 9. / 40. &
                &* (119. + 110. * rhodl) + 18. / 19. * (25. - 44. * rhodl) &
                &- 1. / 9. * (269. - 462 * rhodl) + 6. / 17. * (29. - 132 &
                &* rhodl) + 3. / 16. * (19. + 165. * rhodl) - 2. / 15. * (37. &
                &+ 110. * rhodl) + 3. / 7. * (5. + 11. * rhodl) - 6. / 13. &
                &* (1. + 2. * rhodl) + 1. / 24. * (1. + 2. * rhodl))

            p_dim1 = pcoeff_r1 / pRef * rhoCenter_dim &
                &* ((backgroundFlow_dim(1)) ** 2 + (backgroundFlow_dim(2)) ** 2)

            var%pi(i, j, k) = var%pi(i, j, k) - kappaInv * p_dim1 ** kappa

          end do
        end do
      end do

    case('atmosphereatrest')

      var%rho(:, :, :) = 0.
      var%u(:, :, :) = 0.
      var%v(:, :, :) = 0.
      var%w(:, :, :) = 0.
      var%pi(:, :, :) = 0.

    case('SkamarockKlemp94')

      ! read test case input data
      rewind(unit = 10)
      read(unit = 10, nml = bubble)

      if(referenceQuantities == "SI") then
        stop "initialize: SI units not allowed"
      end if

      ! zero start velocity
      var%u(:, :, :) = backgroundFlow_dim(1) / uRef
      var%v(:, :, :) = backgroundFlow_dim(2) / uRef
      var%w(:, :, :) = backgroundFlow_dim(3) / uRef

      ! constant pressure variable pi'
      var%pi(:, :, :) = 0.0

      i00 = is + nbx - 1 ! 0 index,
      ! replace i -> i + i0 in x and y fields
      j00 = js + nby - 1

      ! potential temperature and density
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            x_dim = x(i + i00) * lRef ! dimensional lenghts
            if(topography) then
              ! TFC FJ
              z_dim = heightTFC(i, j, k) * lRef
            else
              z_dim = z(k) * lRef
            end if

            ! delX = (x_dim - xCenter_dim)*60. / lx_dim(1)
            ! delZ = (z_dim - zCenter_dim) / zRadius_dim

            ! r = sqrt(delX**2 + delZ**2)  ! scaled radius

            ! if( r<=1.0 ) then  ! inside bubble
            dTheta_dim = 0.01 * sin((pi * z_dim) / 10000.) / (1. + ((x_dim &
                &- xCenter_dim) * 60. / lx_dim(1)) ** 2)
            ! = 0.5*dTheta0_dim * (1.0 + (cos(pi*r/2.0))**2)

            if(topography) then
              ! TFC FJ
              theta = thetaStratTFC(i, j, k) + dTheta_dim / thetaRef
            else
              theta = thetaStrat(k) + dTheta_dim / thetaRef
            end if

            ! calc pseudo-incompressible density rho*
            if(referenceQuantities == "SI") then
              rho = p0 ** kappa / Rsp * Pstrat(k) / theta - rhoStrat(k)
            else
              if(topography) then
                ! TFC FJ
                rho = pStratTFC(i, j, k) / theta - rhoStratTFC(i, j, k)
              else
                rho = Pstrat(k) / theta - rhoStrat(k)
              end if
            end if

            select case(model)

            case("pseudo_incompressible", "compressible")

              var%rho(i, j, k) = rho

            case("Boussinesq")

              ! Density fluctuations are stored in var%rhop(i, j, k),
              ! var%rho(i, j, k) must remain zero! Boussinesq model
              var%rho(i, j, k) = 0.0
              var%rhop(i, j, k) = rho

            case default
              stop "initialize: unknown model."
            end select

          end do
        end do
      end do

      !-----------------------------------------------------------------
    case('hotBubble_heat')

      ! start velocity
      var%rho(:, :, :) = 0.
      var%u(:, :, :) = 0.0
      var%v(:, :, :) = 0.0
      var%w(:, :, :) = 0.0

      ! constant pressure variable pi'
      var%pi(:, :, :) = 0.0

    case('heatedLayer')

      ! start velocity
      var%rho(:, :, :) = 0.
      var%u(:, :, :) = 0.0
      var%v(:, :, :) = 0.0
      var%w(:, :, :) = 0.0

      ! constant pressure variable pi'
      var%pi(:, :, :) = 0.0

    case('hotBubble_heatedLayer')

      ! start velocity
      var%rho(:, :, :) = 0.
      var%u(:, :, :) = 0.0
      var%v(:, :, :) = 0.0
      var%w(:, :, :) = 0.0

      ! constant pressure variable pi'
      var%pi(:, :, :) = 0.0

    case('hotBubble')

      ! read test case input data
      rewind(unit = 10)
      read(unit = 10, nml = bubble)

      ! start velocity
      var%u(:, :, :) = backgroundFlow(1)
      var%v(:, :, :) = backgroundFlow(2)
      var%w(:, :, :) = backgroundFlow(3)

      ! constant pressure variable pi'
      var%pi(:, :, :) = 0.0

      ! zero potential temperature
      var%rhop(:, :, :) = 0.0

      ! potential temperature and density

      dTheta0 = dTheta0_dim / thetaRef

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            x_dim = x(i) * lRef ! dimensional lenghts
            if(topography) then
              ! TFC FJ
              z_dim = heightTFC(i, j, k) * lRef
            else
              z_dim = z(k) * lRef
            end if

            delX = (x_dim - xCenter_dim) / xRadius_dim
            delZ = (z_dim - zCenter_dim) / zRadius_dim
            delZ = zExcentricity * delZ

            r = sqrt(delX ** 2 + delZ ** 2) ! scaled radius

            if(r <= 1.0) then
              !--------------------
              !    inside bubble
              !--------------------

              dTheta = dTheta0 * (cos(pi * r / 2.0)) ** 2

              if(topography) then
                ! TFC FJ
                theta = thetaStratTFC(i, j, k) + dTheta
              else
                theta = thetaStrat(k) + dTheta
              end if

              select case(model)

              case("pseudo_incompressible", "compressible")

                ! calc pseudo-incompressible density rho*
                if(referenceQuantities == "SI") then
                  rho = p0 ** kappa / Rsp * Pstrat(k) / theta - rhoStrat(k)
                else
                  if(topography) then
                    ! TFC FJ
                    rho = pStratTFC(i, j, k) / theta - rhoStratTFC(i, j, k)
                  else
                    rho = Pstrat(k) / theta - rhoStrat(k)
                  end if
                end if

                var%rho(i, j, k) = rho

              case("Boussinesq")

                ! TFC FJ
                ! Density fluctuations are stored in var%rhop(i, j, k),
                ! var%rho(i, j, k) must remain zero!
                if(topography) then
                  rho = pStratTFC(i, j, k) / theta - rhoStratTFC(i, j, k)
                else
                  rho = pStrat(k) / theta - rhoStrat(k)
                end if
                var%rhop(i, j, k) = rho

              case default
                stop "initialize: unknown model."
              end select

            else
              !------------------------------------------
              !  outside bubble keep background density
              !------------------------------------------

              var%rho(i, j, k) = 0.0
            end if

          end do
        end do
      end do

      !------------------------------------------------------------------

    case('hotBubble2D')

      ! read test case input data
      rewind(unit = 10)
      read(unit = 10, nml = wavePacket)
      rewind(unit = 10)
      read(unit = 10, nml = bubble)

      if(referenceQuantities == "SI") then
        stop "initialize: SI units not allowed"
      end if

      ! start velocity
      var%u(:, :, :) = 0.0
      var%v(:, :, :) = 0.0
      var%w(:, :, :) = 0.0

      ! constant pressure variable pi'
      var%pi(:, :, :) = 0.0

      ! set local index
      i00 = is + nbx - 1

      ! potential temperature and density
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            x_dim = x(i + i00) * lRef ! dimensional lenghts
            if(topography) then
              ! TFC FJ
              z_dim = heightTFC(i, j, k) * lRef
            else
              z_dim = z(k) * lRef
            end if

            delX = (x_dim - xCenter_dim) / xRadius_dim
            delZ = (z_dim - zCenter_dim) / zRadius_dim

            r = sqrt(delX ** 2 + delZ ** 2) ! scaled radius

            if(r <= 1.0) then ! inside bubble

              dTheta_dim = dTheta0_dim * (cos(pi * r / 2.0)) ** 2

              if(topography) then
                ! TFC FJ
                theta = thetaStratTFC(i, j, k) + dTheta_dim / thetaRef
              else
                theta = thetaStrat(k) + dTheta_dim / thetaRef
              end if

              ! calc pseudo-incompressible density rho*
              if(referenceQuantities == "SI") then
                rho = p0 ** kappa / Rsp * Pstrat(k) / theta - rhoStrat(k)
              else
                if(topography) then
                  rho = pStratTFC(i, j, k) / theta - rhoStratTFC(i, j, k)
                else
                  rho = Pstrat(k) / theta - rhoStrat(k)
                end if
              end if

              select case(model)

              case("pseudo_incompressible", "compressible")

                var%rho(i, j, k) = rho

              case("Boussinesq")

                ! Density fluctuations are stored in var%rhop(i, j, k),
                ! var%rho(i, j, k) must remain zero!
                var%rho(i, j, k) = 0.0
                var%rhop(i, j, k) = rho

              case default
                stop "initialize: unknown model."
              end select

            else ! outside bubble

              select case(model)

              case("pseudo_incompressible", "compressible")

                var%rho(i, j, k) = 0.0

              case("Boussinesq")

                ! TFC FJ
                ! Density fluctuations are stored in var(i, j, k, 6),
                ! var(i, j, k, 1) must remain zero!
                var%rho(i, j, k) = 0.0
                var%rhop(i, j, k) = 0.0

              case default
                stop "initialize: unknown model."
              end select

            end if

          end do
        end do
      end do

      !-------------------------------------------------------------------

    case('hotBubble3D')

      ! read test case input data
      rewind(unit = 10)
      read(unit = 10, nml = wavePacket)
      rewind(unit = 10)
      read(unit = 10, nml = bubble)

      if(referenceQuantities == "SI") then
        stop "initialize: SI units not allowed"
      end if

      ! zero start velocity
      var%u(:, :, :) = 0.0
      var%v(:, :, :) = 0.0
      var%w(:, :, :) = 0.0

      ! constant pressure variable pi'
      var%pi(:, :, :) = 0.0

      ! set local index
      i00 = is + nbx - 1
      j00 = js + nby - 1

      ! potential temperature and density
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            x_dim = x(i + i00) * lRef ! dimensional lengh
            y_dim = y(j + j00) * lRef
            if(topography) then
              ! TFC FJ
              z_dim = heightTFC(i, j, k) * lRef
            else
              z_dim = z(k) * lRef
            end if

            delX = (x_dim - xCenter_dim) / xRadius_dim
            delY = (y_dim - xCenter_dim) / xRadius_dim
            delZ = (z_dim - zCenter_dim) / zRadius_dim

            r = sqrt(delX ** 2 + delY ** 2 + delZ ** 2) ! scaled radius

            if(r <= 1.0) then ! inside bubble

              dTheta_dim = dTheta0_dim * (cos(pi * r / 2.0)) ** 2

              if(topography) then
                ! TFC FJ
                theta = thetaStratTFC(i, j, k) + dTheta_dim / thetaRef
              else
                theta = thetaStrat(k) + dTheta_dim / thetaRef
              end if

              ! calc pseudo-incompressible density rho*
              if(referenceQuantities == "SI") then
                rho = p0 ** kappa / Rsp * Pstrat(k) / theta - rhoStrat(k)
              else
                if(topography) then
                  rho = pStratTFC(i, j, k) / theta - rhoStratTFC(i, j, k)
                else
                  rho = Pstrat(k) / theta - rhoStrat(k)
                end if
              end if

              select case(model)

              case("pseudo_incompressible", "compressible")

                var%rho(i, j, k) = rho

              case("Boussinesq")

                ! Density fluctuations are stored in var(i, j, k, 6),
                ! var(i, j, k, 1) must remain zero!
                var%rho(i, j, k) = 0.0
                var%rhop(i, j, k) = rho

              case default
                stop "initialize: unknown model."
              end select
            else ! outside bubble
              select case(model)

              case("pseudo_incompressible", "compressible")

                var%rho(i, j, k) = 0.0

              case("Boussinesq")

                ! Density fluctuations are stored in var%rhop(i, j, k),
                ! var%rho(i, j, k) must remain zero!
                var%rho(i, j, k) = 0.0
                var%rhop(i, j, k) = 0.0

              case default
                stop "initialize: unknown model."
              end select
            end if
          end do
        end do
      end do

      !----------------------------------------------------------------
      !               Baroclinic life cycle: realistic
      !
      ! either setup from Kuehnlein et al (2012): background = 'const-N'
      ! or setup of Held & Suarez (1994): background = 'HeldSuarez'
      !----------------------------------------------------------------

    case('baroclinic_LC')

      if((background /= "const-N") .and. (background /= "HeldSuarez")) then
        stop 'ERROR: baroclinic_LC needs for background either const-N  or &
            &HeldSuarez'
      end if

      ! read test case input data
      rewind(unit = 10)
      read(unit = 10, nml = baroclinic_LC)
      !*****************************************************************
      ! initialize fields
      !*****************************************************************
      var%rho(:, :, :) = 0.0
      var%u(:, :, :) = 0.0
      var%v(:, :, :) = 0.0
      var%w(:, :, :) = 0.0
      var%pi(:, :, :) = 0.0
      var%rhop(:, :, :) = 0.0

      the_env_pp(:, :, :) = 0.0
      dens_env_pp(:, :, :) = 0.0
      p_env_pp(:, :, :) = 0.0

      u_env_pp(:, :, :) = 0.0
      v_env_pp(:, :, :) = 0.0

      ! non-dimensional quantities

      ymin = ly_dim(0) / lRef
      ymax = ly_dim(1) / lRef

      jwdth = jwdth_dim / lRef ! jet width

      if(background == "const-N") then
        z_trpp0 = z_trpp0_dim / lRef ! mean tropopause height

        z_baro = z_baro_dim / lRef ! alt. above which the atmosph. is
        ! barotropic

        deltht = 3.e1 / thetaRef ! merid. potential-temp. contrast

        thet0 = thet0_dim / thetaRef ! characteristic potential temperature

        ntrp = ntrp_dim * tRef ! Brunt-Vaisala frequency troposphere
        nstr = nstr_dim * tRef ! Brunt-Vaisala frequency stratosphere

      else if(background == "HeldSuarez") then
        yjets = 0.75 * ymin + 0.25 * ymax
        yjetn = 0.25 * ymin + 0.75 * ymax

        dy_hs = 2.0 * jwdth
      else
        stop 'ERROR: baroclinic_LC needs for background either const-N  or &
            &HeldSuarez'
      end if

      if(master .and. jwdth > 0.5 * (ymax - ymin)) then
        stop 'ERROR: jet width too large'
      end if

      if(sizeY <= 1) stop 'ERROR: Barocl LC expects sizeY > 1'

      ! set local index

      i00 = is + nbx - 1 ! 0 index,
      ! replace i -> i + i0 in x and y fields
      j00 = js + nby - 1

      if(background == "const-N") then
        ! potential-temperature field

        do j = 1, ny
          yloc = y(j + j00)

          ! tropopause height

          if(yloc < 0.5 * (ymax + ymin)) then
            yjet0 = 0.75 * ymin + 0.25 * ymax

            if(yloc - yjet0 < - jwdth) then
              fstrpp = - 1.e0
            else if(yloc - yjet0 >= - jwdth .and. yloc - yjet0 < jwdth) then
              fstrpp = sin(0.5 * pi * (yloc - yjet0) / jwdth)
            else if(yloc - yjet0 >= jwdth) then
              fstrpp = 1.e0
            end if
          else
            yjet0 = 0.25 * ymin + 0.75 * ymax

            if(yloc - yjet0 < - jwdth) then
              fstrpp = 1.e0
            else if(yloc - yjet0 >= - jwdth .and. yloc - yjet0 < jwdth) then
              fstrpp = - sin(0.5 * pi * (yloc - yjet0) / jwdth)
            else if(yloc - yjet0 >= jwdth) then
              fstrpp = - 1.e0
            end if
          end if

          z_trpp = z_trpp0 + g_ndim * deltht / (2.0 * thet0 * (nstr ** 2 &
              &- ntrp ** 2)) * fstrpp

          if(topography) then
            ! TFC FJ
            do i = 1, nx
              do k = 0, nz + 1
                zloc = heightTFC(i, j, k)

                if(zloc < z_trpp) then
                  ! below tropopause

                  if(yloc < 0.5 * (ymax + ymin)) then
                    yjet = yjet0 - kaptpp * zloc

                    if(yloc - yjet < - jwdth) then
                      fstht = - 1.0e0
                    else if(yloc - yjet >= - jwdth .and. yloc - yjet < jwdth) &
                        &then
                      fstht = sin(0.5 * pi * (yloc - yjet) / jwdth)
                    else if(yloc - yjet >= jwdth) then
                      fstht = 1.0e0
                    end if
                  else
                    yjet = yjet0 + kaptpp * zloc

                    if(yloc - yjet < - jwdth) then
                      fstht = 1.0e0
                    else if(yloc - yjet >= - jwdth .and. yloc - yjet < jwdth) &
                        &then
                      fstht = - sin(0.5 * pi * (yloc - yjet) / jwdth)
                    else if(yloc - yjet >= jwdth) then
                      fstht = - 1.0e0
                    end if
                  end if

                  the_env_pp(i, j, k) = thet0 * (1.e0 + ntrp ** 2.0 / g_ndim &
                      &* zloc) + 0.5 * deltht * fstht
                else
                  ! above tropopause

                  if(yloc < 0.5 * (ymax + ymin)) then
                    yjet = yjet0 - kaptpp * z_trpp

                    if(yloc - yjet < - jwdth) then
                      fstht = - 1.0e0
                    else if(yloc - yjet >= - jwdth .and. yloc - yjet < jwdth) &
                        &then
                      fstht = sin(0.5 * pi * (yloc - yjet) / jwdth)
                    else if(yloc - yjet >= jwdth) then
                      fstht = 1.0e0
                    end if
                  else
                    yjet = yjet0 + kaptpp * z_trpp

                    if(yloc - yjet < - jwdth) then
                      fstht = 1.0e0
                    else if(yloc - yjet >= - jwdth .and. yloc - yjet < jwdth) &
                        &then
                      fstht = - sin(0.5 * pi * (yloc - yjet) / jwdth)
                    else if(yloc - yjet >= jwdth) then
                      fstht = - 1.0e0
                    end if
                  end if

                  if(zloc < z_baro) then
                    zeta_baro = sin(0.5 * pi * (z_baro - zloc) / (z_baro &
                        &- z_trpp))
                  else
                    zeta_baro = 0.0
                  end if

                  the_env_pp(i, j, k) = thet0 * (1.0e0 + nstr ** 2.0 / g_ndim &
                      &* zloc) + zeta_baro * (- thet0 * (nstr ** 2.0 - ntrp &
                      &** 2.0) / g_ndim * z_trpp + 0.5 * deltht * fstht * (1.0 &
                      &+ 5.0 * (z_trpp - zloc) / z_trpp))
                end if

                dens_env_pp(i, j, k) = pStratTFC(i, j, k) / the_env_pp(i, j, k)

                var%rho(i, j, k) = dens_env_pp(i, j, k) - rhoStratTFC(i, j, k)
              end do
            end do
          else
            do k = 0, nz + 1
              zloc = z(k)

              if(zloc < z_trpp) then
                ! below tropopause

                if(yloc < 0.5 * (ymax + ymin)) then
                  yjet = yjet0 - kaptpp * zloc

                  if(yloc - yjet < - jwdth) then
                    fstht = - 1.e0
                  else if(yloc - yjet >= - jwdth .and. yloc - yjet < jwdth) then
                    fstht = sin(0.5 * pi * (yloc - yjet) / jwdth)
                  else if(yloc - yjet >= jwdth) then
                    fstht = 1.e0
                  end if
                else
                  yjet = yjet0 + kaptpp * zloc

                  if(yloc - yjet < - jwdth) then
                    fstht = 1.e0
                  else if(yloc - yjet >= - jwdth .and. yloc - yjet < jwdth) then
                    fstht = - sin(0.5 * pi * (yloc - yjet) / jwdth)
                  else if(yloc - yjet >= jwdth) then
                    fstht = - 1.e0
                  end if
                end if

                the_env_pp(:, j, k) = thet0 * (1.e0 + ntrp ** 2 / g_ndim &
                    &* zloc) + 0.5 * deltht * fstht
              else
                ! above tropopause

                if(yloc < 0.5 * (ymax + ymin)) then
                  yjet = yjet0 - kaptpp * z_trpp

                  if(yloc - yjet < - jwdth) then
                    fstht = - 1.e0
                  else if(yloc - yjet >= - jwdth .and. yloc - yjet < jwdth) then
                    fstht = sin(0.5 * pi * (yloc - yjet) / jwdth)
                  else if(yloc - yjet >= jwdth) then
                    fstht = 1.e0
                  end if
                else
                  yjet = yjet0 + kaptpp * z_trpp

                  if(yloc - yjet < - jwdth) then
                    fstht = 1.e0
                  else if(yloc - yjet >= - jwdth .and. yloc - yjet < jwdth) then
                    fstht = - sin(0.5 * pi * (yloc - yjet) / jwdth)
                  else if(yloc - yjet >= jwdth) then
                    fstht = - 1.e0
                  end if
                end if

                if(zloc < z_baro) then
                  zeta_baro = sin(0.5 * pi * (z_baro - zloc) / (z_baro &
                      &- z_trpp))
                else
                  zeta_baro = 0.0
                end if

                the_env_pp(:, j, k) = thet0 * (1.e0 + nstr ** 2 / g_ndim &
                    &* zloc) + zeta_baro * (- thet0 * (nstr ** 2 - ntrp ** 2) &
                    &/ g_ndim * z_trpp + 0.5 * deltht * fstht * (1.0 + 5.0 &
                    &* (z_trpp - zloc) / z_trpp))
              end if

              dens_env_pp(:, j, k) = Pstrat(k) / the_env_pp(:, j, k)

              var%rho(1:nx, j, k) = dens_env_pp(1:nx, j, k) - rhoStrat(k)
            end do
          end if
        end do

        ! Exner pressure just below the bottom
        ! (so that the bottom pressure vanishes)

        do j = 1, ny
          do i = 1, nx
            rhop0 = var%rho(i, j, 0)
            rhop1 = var%rho(i, j, 1)

            if(topography) then
              ! TFC FJ
              rho0 = var%rho(i, j, 0) + rhoStratTFC(i, j, 0)
              rho1 = var%rho(i, j, 1) + rhoStratTFC(i, j, 1)
            else
              rho0 = var%rho(i, j, 0) + rhoStrat(0)
              rho1 = var%rho(i, j, 1) + rhoStrat(1)
            end if

            buoy0 = - g_ndim * rhop0 / rho0
            buoy1 = - g_ndim * rhop1 / rho1

            rho = 0.5 * (rho0 + rho1)

            ! FJApr2024
            ! It seems like there is one factor 1/2 too much here...
            if(topography) then
              ! TFC FJ
              var%pi(i, j, 0) = 0.5 * heightTFC(i, j, 0) * Ma2 * kappa * rho &
                  &/ (0.5 * (pStratTFC(i, j, 0) + pStratTFC(i, j, 1))) * 0.5 &
                  &* (buoy0 + buoy1)
            else
              var%pi(i, j, 0) = - 0.25 * dz * Ma2 * kappa * rho &
                  &/ PstratTilde(0) * 0.5 * (buoy0 + buoy1)
            end if
          end do
        end do

        ! Exner pressure up to just above the lid

        do k = 0, nz
          do j = 1, ny
            do i = 1, nx
              rhop0 = var%rho(i, j, k)
              rhop1 = var%rho(i, j, k + 1)

              if(topography) then
                ! TFC FJ
                rho0 = var%rho(i, j, k) + rhoStratTFC(i, j, k)
                rho1 = var%rho(i, j, k + 1) + rhoStratTFC(i, j, k + 1)
              else
                rho0 = var%rho(i, j, k) + rhoStrat(k)
                rho1 = var%rho(i, j, k + 1) + rhoStrat(k + 1)
              end if

              buoy0 = - g_ndim * rhop0 / rho0
              buoy1 = - g_ndim * rhop1 / rho1

              rho = 0.5 * (rho0 + rho1)

              ! It seems like there is one factor 1/2 too much here...
              if(topography) then
                var%pi(i, j, k + 1) = var%pi(i, j, k) + 0.5 * 0.5 * (jac(i, j, &
                    &k) + jac(i, j, k + 1)) * dz * Ma2 * kappa * rho / (0.5 &
                    &* (pStratTFC(i, j, k) + pStratTFC(i, j, k + 1))) * 0.5 &
                    &* (buoy0 + buoy1)
              else
                var%pi(i, j, k + 1) = var%pi(i, j, k) + 0.5 * dz * Ma2 * kappa &
                    &* rho / PstratTilde(k) * 0.5 * (buoy0 + buoy1)
              end if
            end do
          end do
        end do
      else if(background == "HeldSuarez") then
        ! define stretched squared sine and cosine field
        ! define stretched quadrupled cosine field

        do j = 1, ny
          yloc = y(j + j00)

          if(corset == 'periodic') then

            s2_strtd(j) = cos(2. * pi * (yloc - ymin) / (ymax - ymin)) ** 2
            c2_strtd(j) = sin(2. * pi * (yloc - ymin) / (ymax - ymin)) ** 2
            c4_strtd(j) = sin(2. * pi * (yloc - ymin) / (ymax - ymin)) ** 4

          else if(corset == 'constant') then

            if(yloc < ymin) then
              stop 'ERROR: y < ymin'
            else if((yloc >= ymin) .and. (yloc < yjets - 0.5 * dy_hs)) then
              s2_strtd(j) = 1.0
              c2_strtd(j) = 0.0
              c4_strtd(j) = 0.0
            else if((yloc >= yjets - 0.5 * dy_hs) .and. (yloc < yjets + 0.5 &
                &* dy_hs)) then
              s2_strtd(j) = sin((yloc - (yjets + 0.5 * dy_hs)) / dy_hs * 0.5 &
                  &* pi) ** 2

              c2_strtd(j) = cos((yloc - (yjets + 0.5 * dy_hs)) / dy_hs * 0.5 &
                  &* pi) ** 2

              c4_strtd(j) = cos((yloc - (yjets + 0.5 * dy_hs)) / dy_hs * 0.5 &
                  &* pi) ** 4
            else if((yloc >= yjets + 0.5 * dy_hs) .and. (yloc < yjetn - 0.5 &
                &* dy_hs)) then
              s2_strtd(j) = 0.0
              c2_strtd(j) = 1.0
              c4_strtd(j) = 1.0
            else if((yloc >= yjetn - 0.5 * dy_hs) .and. (yloc < yjetn + 0.5 &
                &* dy_hs)) then
              s2_strtd(j) = sin((yloc - (yjetn - 0.5 * dy_hs)) / dy_hs * 0.5 &
                  &* pi) ** 2

              c2_strtd(j) = cos((yloc - (yjetn - 0.5 * dy_hs)) / dy_hs * 0.5 &
                  &* pi) ** 2

              c4_strtd(j) = cos((yloc - (yjetn - 0.5 * dy_hs)) / dy_hs * 0.5 &
                  &* pi) ** 4
            else if((yloc >= yjetn + 0.5 * dy_hs) .and. (yloc <= ymax)) then
              s2_strtd(j) = 1.0
              c2_strtd(j) = 0.0
              c4_strtd(j) = 0.0
            else if(yloc > ymax) then
              stop 'ERROR: y > ymax'
            end if

            ! no latitude dependence of the stratification
            c2_strtd(j) = 1.

          else

            stop 'ERROR: wrong corset'

          end if
        end do

        ! y- and z-dependent thermal relaxation rate
        ! z-dependent Rayleigh-damping-rate

        if(ta_hs_dim == 0.0) then
          ka_hs = 0.0
        else
          ka_hs = tRef / ta_hs_dim
        end if

        if(ts_hs_dim == 0.0) then
          ks_hs = 0.0
        else
          ks_hs = tRef / ts_hs_dim
        end if

        if(tf_hs_dim == 0.0) then
          kf_hs = 0.0
        else
          kf_hs = tRef / tf_hs_dim
        end if

        kv_hs = 0.
        kw_hs = 0.

        kr_sp = 0.0
        kr_sp_w = 0.0

        if(topography) then
          kr_sp_tfc = 0.0
          kr_sp_w_tfc = 0.0
          ! Held & Suarez (1994) boundary layer drag (wind)
          do k = 0, nz + 1
            do j = 0, ny + 1
              do i = 0, nx + 1
                sig_pr = piStratTFC(i, j, k) ** (1.0 / kappa)
                facsig = max(0.0, (sig_pr - sigb_hs) / (1.0 - sigb_hs))

                kr_sp_tfc(i, j, k) = kf_hs * facsig
                kr_sp_w_tfc(i, j, k) = 0.0
              end do
            end do
          end do
          ! Held & Suarez (1994) boundary layer drag (temperature)
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                sig_pr = piStratTFC(i, j, k) ** (1.0 / kappa)
                facsig = max(0.0, (sig_pr - sigb_hs) / (1.0 - sigb_hs))

                kt_hs_tfc(i, j, k) = ka_hs + (ks_hs - ka_hs) * c4_strtd(j) &
                    &* facsig
              end do
            end do
          end do
        else
          do k = 0, nz + 1
            sig_pr = pistrat(k) ** (1.0 / kappa)
            facsig = max(0.0, (sig_pr - sigb_hs) / (1.0 - sigb_hs))

            !Held + Suarez 1994: BL drag:
            kr_sp(:, k) = kf_hs * facsig
            kv_hs(:, k) = 0.
            kr_sp_w(:, k) = 0.

            do j = 1, ny
              ! Held Suarez 1994:
              kt_hs(j, k) = ka_hs + (ks_hs - ka_hs) * c4_strtd(j) * facsig

            end do
          end do
        end if

        do j = 1, ny
          ! (total) Exner pressure just above and below the surface so
          ! that it = 1 at the surface, by an Euler integration of
          ! hydrostatic equilibrium

          if(topography) then
            ! TFC FJ
            do i = 1, nx
              ! Set Exner pressure at the surface.
              tempev = max(tp_strato, tp_srf_trp - tpdiffhor_tropo &
                  &* s2_strtd(j))
              var%pi(i, j, 0) = 1.0 - heightTFC(i, j, 0) * kappa / tempev
              var%pi(i, j, 1) = 1.0 - heightTFC(i, j, 1) * kappa / tempev
              ! Set potential temperature at the surface.
              do k = 0, 1
                tempev = max(tp_strato, var%pi(i, j, k) * (tp_srf_trp &
                    &- tpdiffhor_tropo * s2_strtd(j) - ptdiffvert_tropo &
                    &/ kappa * log(var%pi(i, j, k)) * c2_strtd(j)))
                the_env_pp(i, j, k) = tempev / var%pi(i, j, k)
              end do
            end do
          else
            tempev = max(tp_strato, tp_srf_trp - tpdiffhor_tropo * s2_strtd(j))

            var%pi(1:nx, j, 0) = 1.0 + 0.5 * dz * kappa / tempev
            var%pi(1:nx, j, 1) = 1.0 - 0.5 * dz * kappa / tempev

            ! potential temperature just below and above the surface

            do k = 0, 1
              do i = 1, nx
                tempev = max(tp_strato, var%pi(i, j, k) * (tp_srf_trp &
                    &- tpdiffhor_tropo * s2_strtd(j) - ptdiffvert_tropo &
                    &/ kappa * log(var%pi(i, j, k)) * c2_strtd(j)))

                the_env_pp(i, j, k) = tempev / var%pi(i, j, k)
              end do
            end do
          end if

          ! for k > 1:
          ! Exner pressure and potential temperature by upward
          ! integration of hydrostatic equilibrium, using a trapezoidal
          ! leapfrog

          switch = 0

          do k = 2, nz + 1
            do i = 1, nx
              if(topography) then
                ! TFC FJ
                pistar = var%pi(i, j, k - 2) - 2.0 * dz * jac(i, j, k - 1) &
                    &* kappa / the_env_pp(i, j, k - 1)
              else
                pistar = var%pi(i, j, k - 2) - 2.0 * dz * kappa &
                    &/ the_env_pp(i, j, k - 1)
              end if

              tempev = max(tp_strato, pistar * (tp_srf_trp - tpdiffhor_tropo &
                  &* s2_strtd(j) - ptdiffvert_tropo / kappa * log(pistar) &
                  &* c2_strtd(j)))

              thetastar = tempev / pistar

              if(topography) then
                ! TFC FJ
                var%pi(i, j, k) = var%pi(i, j, k - 1) - 0.5 * dz * 0.5 &
                    &* (jac(i, j, k) + jac(i, j, k - 1)) * (kappa / thetastar &
                    &+ kappa / the_env_pp(i, j, k - 1))
              else
                var%pi(i, j, k) = var%pi(i, j, k - 1) - 0.5 * dz * (kappa &
                    &/ thetastar + kappa / the_env_pp(i, j, k - 1))
              end if

              tempev = max(tp_strato, var%pi(i, j, k) * (tp_srf_trp &
                  &- tpdiffhor_tropo * s2_strtd(j) - ptdiffvert_tropo / kappa &
                  &* log(var%pi(i, j, k)) * c2_strtd(j)))

              the_env_pp(i, j, k) = tempev / var%pi(i, j, k)
            end do
          end do

          ! reduction of the jets to zero
          ! (within the lower half of the sponge)
          ! by reduction of the pressure fluctuations

          if(spongeLayer) then
            if(topography) then
              do k = 0, nz + 1
                do i = 1, nx
                  if(heightTFC(i, j, k) >= zSponge) then
                    if(heightTFC(i, j, k) >= zSponge + 0.5 * (lz(1) &
                        &- zSponge)) then
                      var%pi(i, j, k) = piStratTFC(i, j, k)
                    else
                      var%pi(i, j, k) = piStratTFC(i, j, k) + cos(0.5 * pi &
                          &* (heightTFC(i, j, k) - zSponge) / (0.5 * (lz(1) &
                          &- zSponge))) ** 2.0 * (var%pi(i, j, k) &
                          &- piStratTFC(i, j, k))
                    end if
                    the_env_pp(i, j, k) = - kappa * 0.5 * (jac(i, j, k) &
                        &+ jac(i, j, k - 1)) * dz / (var%pi(i, j, k) &
                        &- var%pi(i, j, k - 1))
                  end if
                end do
              end do
            else
              kshalf = kSponge + int((nz - kSponge) / 2)

              do k = kSponge, kshalf
                do i = 1, nx
                  var%pi(i, j, k) = pistrat(k) + cos(0.5 * pi * (z(k) &
                      &- zSponge) / (z(kshalf) - zSponge)) ** 2 * (var%pi(i, &
                      &j, k) - pistrat(k))
                end do
              end do

              do k = kshalf + 1, nz + 1
                do i = 1, nx
                  var%pi(i, j, k) = pistrat(k)
                end do
              end do

              do k = kSponge, nz + 1
                do i = 1, nx
                  the_env_pp(i, j, k) = - kappa * dz / (var%pi(i, j, k) &
                      &- var%pi(i, j, k - 1))
                end do
              end do
            end if
          end if

          ! density
          if(topography) then
            do i = 1, nx
              do k = 0, nz + 1
                dens_env_pp(i, j, k) = pStratTFC(i, j, k) / the_env_pp(i, j, k)
                var%rho(i, j, k) = dens_env_pp(i, j, k) - rhoStratTFC(i, j, k)
              end do
            end do
          else
            do k = 0, nz + 1
              dens_env_pp(:, j, k) = Pstrat(k) / the_env_pp(:, j, k)

              var%rho(1:nx, j, k) = dens_env_pp(1:nx, j, k) - rhoStrat(k)
            end do
          end if
        end do

        ! in case of periodic Coriolis parameter choose horizontally
        ! homogeneous atmosphere at rest as equilibrium

        if(corset == 'periodic') then
          if(topography) then
            ! TFC FJ
            do i = 1, nx
              do j = 1, ny
                do k = 0, nz + 1
                  dens_env_pp(i, j, k) = rhoStratTFC(i, j, k)
                  var%rho(i, j, k) = 0.0
                  var%pi(i, j, k) = piStratTFC(i, j, k)
                end do
              end do
            end do
          else
            do k = 0, nz + 1
              dens_env_pp(:, :, k) = rhoStrat(k)

              var%rho(:, :, k) = 0.

              var%pi(:, :, k) = pistrat(k)
            end do
          end if
        end if

        ! subtract reference-atmosphere Exner pressure from the total

        if(topography) then
          ! TFC FJ
          do i = 1, nx
            do j = 1, ny
              do k = 0, nz + 1
                var%pi(i, j, k) = var%pi(i, j, k) - piStratTFC(i, j, k)
              end do
            end do
          end do
        else
          do k = 0, nz + 1
            var%pi(1:nx, 1:ny, k) = var%pi(1:nx, 1:ny, k) - pistrat(k)
          end do
        end if
      else
        stop 'ERROR: wrong background for baroclinic_LC'
      end if

      p_env_pp(1:nx, 1:ny, 0:nz + 1) = var%pi(1:nx, 1:ny, 0:nz + 1)

      call setHalos(var, "var")
      call setBoundary(var, flux, "var")

      ! determine horizontal wind from density and Exner-pressure
      ! fluctuations

      ! in case of periodic Coriolis parameter initilization zero wind
      ! (that is then also an equlibrium wind)

      if(corset == 'periodic') then

        var%u(:, :, :) = 0.

        u_env_pp = 0.
        v_env_pp = 0.

      else if(corset == 'constant') then

        do k = 0, nz + 1
          do j = 1, ny
            do i = 1, nx
              ! u-wind

              rho_int_0m = 0.5 * (var%rho(i, j - 1, k) + var%rho(i, j, k))
              rho_int_00 = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k))
              rho_int_pm = 0.5 * (var%rho(i + 1, j - 1, k) + var%rho(i + 1, j, &
                  &k))
              rho_int_p0 = 0.5 * (var%rho(i + 1, j, k) + var%rho(i + 1, j + 1, &
                  &k))

              rho = var%rho(i, j, k)

              if(topography) then
                ! TFC FJ
                rho_int_0m = rho_int_0m + 0.5 * (rhoStratTFC(i, j, k) &
                    &+ rhoStratTFC(i, j - 1, k))
                rho_int_00 = rho_int_00 + 0.5 * (rhoStratTFC(i, j, k) &
                    &+ rhoStratTFC(i, j + 1, k))
                rho_int_pm = rho_int_pm + 0.5 * (rhoStratTFC(i + 1, j, k) &
                    &+ rhoStratTFC(i + 1, j - 1, k))
                rho_int_p0 = rho_int_p0 + 0.5 * (rhoStratTFC(i + 1, j, k) &
                    &+ rhoStratTFC(i + 1, j + 1, k))
                rho = rho + rhoStratTFC(i, j, k)
              else
                rho_int_0m = rho_int_0m + rhoStrat(k)
                rho_int_00 = rho_int_00 + rhoStrat(k)
                rho_int_pm = rho_int_pm + rhoStrat(k)
                rho_int_p0 = rho_int_p0 + rhoStrat(k)

                rho = rho + rhoStrat(k)
              end if

              yloc = y(j + j00)
              ! latitude-dependent Coriolis parameter re-established
              f_Coriolis_y(j) = f_Coriolis_dim

              if(f_Coriolis_y(j) /= 0.0) then
                if(topography) then

                  ! Compute values at cell edges.
                  pEdgeF = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                      &+ pStratTFC(i, j + 1, k) / jac(i, j + 1, k))
                  pEdgeB = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                      &+ pStratTFC(i, j - 1, k) / jac(i, j - 1, k))
                  pREdgeF = 0.5 * (pStratTFC(i + 1, j, k) / jac(i + 1, j, k) &
                      &+ pStratTFC(i + 1, j + 1, k) / jac(i + 1, j + 1, k))
                  pREdgeB = 0.5 * (pStratTFC(i + 1, j, k) / jac(i + 1, j, k) &
                      &+ pStratTFC(i + 1, j - 1, k) / jac(i + 1, j - 1, k))
                  ! Compute pressure gradient component.
                  piUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                      &* var%pi(i, j, k + 1) + jac(i, j + 1, k + 1) * met(i, j &
                      &+ 1, k + 1, 2, 3) * var%pi(i, j + 1, k + 1))
                  piDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) &
                      &* var%pi(i, j, k - 1) + jac(i, j + 1, k - 1) * met(i, j &
                      &+ 1, k - 1, 2, 3) * var%pi(i, j + 1, k - 1))
                  piUEdgeB = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                      &* var%pi(i, j, k + 1) + jac(i, j - 1, k + 1) * met(i, j &
                      &- 1, k + 1, 2, 3) * var%pi(i, j - 1, k + 1))
                  piDEdgeB = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) &
                      &* var%pi(i, j, k - 1) + jac(i, j - 1, k - 1) * met(i, j &
                      &- 1, k - 1, 2, 3) * var%pi(i, j - 1, k - 1))
                  piRUEdgeF = 0.5 * (jac(i + 1, j, k + 1) * met(i + 1, j, k &
                      &+ 1, 2, 3) * var%pi(i + 1, j, k + 1) + jac(i + 1, j &
                      &+ 1, k + 1) * met(i + 1, j + 1, k + 1, 2, 3) * var%pi(i &
                      &+ 1, j + 1, k + 1))
                  piRDEdgeF = 0.5 * (jac(i + 1, j, k - 1) * met(i + 1, j, k &
                      &- 1, 2, 3) * var%pi(i + 1, j, k - 1) + jac(i + 1, j &
                      &+ 1, k - 1) * met(i + 1, j + 1, k - 1, 2, 3) * var%pi(i &
                      &+ 1, j + 1, k - 1))
                  piRUEdgeB = 0.5 * (jac(i + 1, j, k + 1) * met(i + 1, j, k &
                      &+ 1, 2, 3) * var%pi(i + 1, j, k + 1) + jac(i + 1, j &
                      &- 1, k + 1) * met(i + 1, j - 1, k + 1, 2, 3) * var%pi(i &
                      &+ 1, j - 1, k + 1))
                  piRDEdgeB = 0.5 * (jac(i + 1, j, k - 1) * met(i + 1, j, k &
                      &- 1, 2, 3) * var%pi(i + 1, j, k - 1) + jac(i + 1, j &
                      &- 1, k - 1) * met(i + 1, j - 1, k - 1, 2, 3) * var%pi(i &
                      &+ 1, j - 1, k - 1))
                  piGrad = 0.25 * (pEdgeF / rho_int_00 * ((jac(i, j + 1, k) &
                      &* var%pi(i, j + 1, k) - jac(i, j, k) * var%pi(i, j, k)) &
                      &/ dy + (piUEdgeF - piDEdgeF) * 0.5 / dz) + pEdgeB &
                      &/ rho_int_0m * ((jac(i, j, k) * var%pi(i, j, k) &
                      &- jac(i, j - 1, k) * var%pi(i, j - 1, k)) / dy &
                      &+ (piUEdgeB - piDEdgeB) * 0.5 / dz) + pREdgeF &
                      &/ rho_int_p0 * ((jac(i + 1, j + 1, k) * var%pi(i + 1, j &
                      &+ 1, k) - jac(i + 1, j, k) * var%pi(i + 1, j, k)) / dy &
                      &+ (piRUEdgeF - piRDEdgeF) * 0.5 / dz) + pREdgeB &
                      &/ rho_int_pm * ((jac(i + 1, j, k) * var%pi(i + 1, j, k) &
                      &- jac(i + 1, j - 1, k) * var%pi(i + 1, j - 1, k)) / dy &
                      &+ (piRUEdgeB - piRDEdgeB) * 0.5 / dz))
                  ! Compute zonal wind.
                  var%u(i, j, k) = - uRef / f_Coriolis_y(j) / lRef / (Ma2 &
                      &* kappa) * piGrad
                else
                  var%u(i, j, k) = - (uRef / f_Coriolis_y(j) / lRef) / (Ma2 &
                      &* kappa * dy) * Pstrat(k) * 0.25 * ((var%pi(i, j, k) &
                      &- var%pi(i, j - 1, k)) / rho_int_0m + (var%pi(i, j + 1, &
                      &k) - var%pi(i, j, k)) / rho_int_00 + (var%pi(i + 1, j, &
                      &k) - var%pi(i + 1, j - 1, k)) / rho_int_pm + (var%pi(i &
                      &+ 1, j + 1, k) - var%pi(i + 1, j, k)) / rho_int_p0)
                end if

                u_env_pp(i, j, k) = var%u(i, j, k)
                v_env_pp(i, j, k) = 0.
              end if

            end do
          end do
        end do

      else

        stop 'ERROR: wrong corset'

      end if

      ! density fluctuations

      if(timeScheme == "semiimplicit" .or. auxil_equ) then
        do k = - nbz, nz + nbz
          var%rhop(:, :, k) = var%rho(:, :, k)
        end do
      end if

      call setHalos(var, "var")
      call setBoundary(var, flux, "var")

      var_env = var ! store environmental state also in var_env

      !-----------------------------------------------------------
      ! add local potential-temperature perturbation
      !-----------------------------------------------------------

      if(add_ptptb) then
        ptptb_x = ptptb_x_dim / lRef
        ptptb_y = ptptb_y_dim / lRef
        ptptb_z = ptptb_z_dim / lRef

        ptptb_dh = ptptb_dh_dim / lRef
        ptptb_dz = ptptb_dz_dim / lRef

        ptptb_amp = ptptb_amp_dim / thetaRef

        do k = 1, nz
          zloc = z(k)

          do j = 1, ny
            yloc = y(j00 + j)

            do i = 1, nx
              xloc = x(i00 + i)

              ! TFC FJ
              if(topography) then
                zloc = heightTFC(i, j, k)
              end if

              if(ptptb_dh <= 0.) then
                rptb = sqrt(((zloc - ptptb_z) / ptptb_dz) ** 2)
              else

                rptb = sqrt(((xloc - ptptb_x) / ptptb_dh) ** 2 + ((yloc &
                    &- ptptb_y) / ptptb_dh) ** 2 + ((zloc - ptptb_z) &
                    &/ ptptb_dz) ** 2)
              end if

              if(rptb <= 1.0) then
                thtptb = ptptb_amp * cos(0.5 * pi * rptb) ** 2
              else
                thtptb = 0.0
              end if

              if(topography) then
                ! TFC FJ
                rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)
              else
                rho = var%rho(i, j, k) + rhoStrat(k)
              end if

              if(topography) then
                ! TFC FJ
                theta = pStratTFC(i, j, k) / rho + thtptb
              else
                theta = Pstrat(k) / rho + thtptb
              end if

              if(topography) then
                ! TFC FJ
                var%rho(i, j, k) = pStratTFC(i, j, k) / theta - rhoStratTFC(i, &
                    &j, k)
              else
                var%rho(i, j, k) = Pstrat(k) / theta - rhoStrat(k)
              end if

            end do
          end do
        end do

        ! add local PT perturbation on SH !FS
        ptptb_y = (- 1.) * ptptb_y
        do k = 1, nz
          zloc = z(k)

          do j = 1, ny
            yloc = y(j00 + j)

            do i = 1, nx
              xloc = x(i00 + i)

              ! TFC FJ
              if(topography) then
                zloc = heightTFC(i, j, k)
              end if

              rptb = sqrt(((xloc - ptptb_x) / ptptb_dh) ** 2 + ((yloc &
                  &- ptptb_y) / ptptb_dh) ** 2 + ((zloc - ptptb_z) / ptptb_dz) &
                  &** 2)

              if(rptb <= 1.0) then
                thtptb = ptptb_amp * cos(0.5 * pi * rptb) ** 2
              else
                thtptb = 0.0
              end if

              if(topography) then
                ! TFC FJ
                rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)
              else
                rho = var%rho(i, j, k) + rhoStrat(k) !Pstrat(k)
              end if

              if(topography) then
                ! TFC FJ
                theta = pStratTFC(i, j, k) / rho - thtptb
              else
                theta = Pstrat(k) / rho - thtptb
              end if

              if(topography) then
                ! TFC FJ
                var%rho(i, j, k) = pStratTFC(i, j, k) / theta - rhoStratTFC(i, &
                    &j, k)
              else
                var%rho(i, j, k) = Pstrat(k) / theta - rhoStrat(k)
              end if
            end do
          end do
        end do
      end if

      !------------------------------------------
      !        Adding the noise
      !------------------------------------------

      ! noise added to density fluctuations in a relative manner

      noise_var = "none"
      noise_mag = 1.0

      if(add_noise) then
        call Random_Seed()
        call noise_array(noise_mag, noise_var, noise)

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              !UAC noise to b emultiplied to full density in case of
              ! periodic Coriolis parameter
              if(corset == 'periodic') then
                if(topography) then
                  ! TFC FJ
                  var%rho(i, j, k) = rhoStratTFC(i, j, k) * noise(i, j, k)
                else
                  var%rho(i, j, k) = rhostrat(k) * noise(i, j, k)
                end if
              else if(corset == 'constant') then
                var%rho(i, j, k) = var%rho(i, j, k) * (1.0 + noise(i, j, k))
              else
                stop 'ERROR: wrong corset'
              end if
            end do
          end do
        end do
      end if

      !------------------------------------------------
      !               Baroclinic life cycle: idealistic
      !------------------------------------------------
    case('baroclinic_ID')

      ! read test case input data
      rewind(unit = 10)
      read(unit = 10, nml = baroclinic_ID)
      !*****************************************************************
      ! initialize fields
      !*****************************************************************
      var%rho(:, :, :) = 0.0
      var%u(:, :, :) = 0.0
      var%v(:, :, :) = 0.0
      var%w(:, :, :) = 0.0
      var%pi(:, :, :) = 0.0
      var%rhop(:, :, :) = 0.0

      the_env_pp(:, :, :) = 0.0

      dens_env_pp(:, :, :) = 0.0

      ! // TODO would be much nicer if everything was coded in
      ! non-dimensional units

      ! // TODO it is unfortunate that the below is so far most of the
      ! time applied with isothermal atmosphere
      ! would be better to have a constant lapse rate of 6K/km, e.g.,
      ! for the troposphere and isothermal above the tropopause

      Ly_end = ly_dim(1)

      ! to also enable test of 2D dynamics in x-z plane
      if(sizeX > 1 .and. sizeY == 1) then
        Lx_end = lx_dim(1)
        bar_sigma_x = bar_sigma_y
      end if

      L_z = z(nz) * lRef

      if(topography) then
        ! TFC FJ
        do i = 1, nx
          do j = 1, ny
            ! Find numerical tropopause.
            k_2_tfc(i, j) = int((levelTFC(i, j, z_tr_dim / lRef) - 0.5 * dz) &
                &/ dz) + 1
            z_2_tfc(i, j) = heightTFC(i, j, k_2_tfc(i, j)) * lRef
            theta_bar_0_tfc(i, j) = thetaStratTFC(i, j, k_2_tfc(i, j)) &
                &* thetaRef

            ! Set alpha_t_tfc.
            if(sizeX > 1 .and. sizeY == 1) then
              alpha_t_tfc(i, j) = u_strength * f_Coriolis_dim &
                  &* theta_bar_0_tfc(i, j) * Lx_end * bar_sigma_x / (dTh_atm &
                  &* g * z_2_tfc(i, j) * pi)
            else
              alpha_t_tfc(i, j) = u_strength * f_Coriolis_dim &
                  &* theta_bar_0_tfc(i, j) * Ly_end * bar_sigma_y / (dTh_atm &
                  &* g * z_2_tfc(i, j) * pi)
            end if
          end do
        end do
      else
        k_2 = int((z_tr_dim - 0.5 * dz * lRef) / (dz * lRef)) + 1
        z_2 = z(k_2) * lRef

        theta_bar_0 = thetaStrat(k_2) * thetaRef ! at the tropopause

        if(master) print *, 'theta_bar_0 = ', theta_bar_0

        ! to also enable test of 2D dynamics in x-z plane
        if(sizeX > 1 .and. sizeY == 1) then
          alpha_t = u_strength * f_Coriolis_dim * theta_bar_0 * Lx_end &
              &* bar_sigma_x / (dTh_atm * g * z_2 * pi)
        else
          alpha_t = u_strength * f_Coriolis_dim * theta_bar_0 * Ly_end &
              &* bar_sigma_y / (dTh_atm * g * z_2 * pi)
        end if
      end if

      ! TFC FJ
      if(master .and. .not. topography) then
        print *, "Barocl. test. instability measure alpha*dTh_atm/dz = ", &
            &alpha_t * dTh_atm / (dz * lRef)
        print *, "Barocl. test. magnitude of \theta alpha*dTh_atm = ", alpha_t &
            &* dTh_atm
      end if

      if(init_2Dto3D) then

        call read_netCDF(- 1, var)

        ! store env pot temp: local
        the_env_pp(1:nx, 1:ny, 1:nz) = var%rhop(1:nx, 1:ny, 1:nz)

        ! store env dens: local
        dens_env_pp(1:nx, 1:ny, 1:nz) = var%rho(1:nx, 1:ny, 1:nz)

        ! actually not used. left as backup for sponge.
        u_env_pp(0:nx, 1:ny, 1:nz) = var%u(0:nx, 1:ny, 1:nz)
      else
        !***************************************************************
        ! tropopause height: Z_trop,
        ! temperature at the ground: T0_t,
        ! lapse rates in the troposphere and stratosphere: gamma_tr,
        !                                                  gamma_st,
        ! pressure at a reference height z_ref: P_ref
        ! exact pressure at tropopause: P_trop
        !***************************************************************

        ! set local index
        i00 = is + nbx - 1 ! 0 index,
        ! replace i -> i + i0 in x and y fields
        j00 = js + nby - 1

        ! TFC FJ
        if(master .and. .not. topography) then
          print *, 'k_2 = ', k_2, 'Numerical tropopause is at z = ', z_2, ', &
              &while H_t = ', z_tr_dim
        end if

        cp = Rsp / kappa

        if(init_bal == "geostr_id") then
          ! determine the Exner-pressure fluctuations

          if(topography) then
            ! TFC FJ
            do i = 1, nx
              do j = 1, ny
                do k = 1, nz
                  ! Set F_a.
                  z_dim = heightTFC(i, j, k) * lRef
                  if(z_dim .lt. z_2_tfc(i, j)) then
                    F_0_tfc(i, j, k) = alpha_t_tfc(i, j) * (2.0 * z_dim &
                        &* sin(pi * z_dim / (2.0 * z_2_tfc(i, j))) - z_dim &
                        &** 2.0 / z_2_tfc(i, j))
                  else
                    F_0_tfc(i, j, k) = alpha_t_tfc(i, j) * (2.0 * (z_dim &
                        &- L_z) * sin(pi * (z_dim - L_z) / (2.0 * (z_2_tfc(i, &
                        &j) - L_z))) - (z_dim - L_z) ** 2.0 / (z_2_tfc(i, j) &
                        &- L_z)) * z_2_tfc(i, j) / (z_2_tfc(i, j) - L_z)
                  end if
                  F_a = F_0_tfc(i, j, k)

                  ! Determine Exner-pressure fluctuations.
                  if(sizeX > 1 .and. sizeY == 1) then
                    x_dim = x(i + i00) * lRef

                    if(x_dim < Lx_end * (0.25 - 0.5 * bar_sigma_x)) then
                      term_a = 1.0
                    else if(x_dim >= Lx_end * (0.25 - 0.5 * bar_sigma_x) .and. &
                        &x_dim < Lx_end * (0.25 + 0.5 * bar_sigma_x)) then
                      term_a = 0.5 * (1.0 - sin(pi * ((x_dim - 0.25 * Lx_end) &
                          &/ (Lx_end * bar_sigma_x))))
                    else if(x_dim >= Lx_end * (0.25 + 0.5 * bar_sigma_x) .and. &
                        &x_dim < Lx_end * (0.75 - 0.5 * bar_sigma_x)) then
                      term_a = 0.0
                    else if(x_dim >= Lx_end * (0.75 - 0.5 * bar_sigma_x) .and. &
                        &x_dim < Lx_end * (0.75 + 0.5 * bar_sigma_x)) then
                      term_a = 0.5 * (1.0 + sin(pi * ((x_dim - 0.75 * Lx_end) &
                          &/ (Lx_end * bar_sigma_x))))
                    else
                      term_a = 1.0
                    end if

                    term_b = (dTh_atm - 2.0 * dTh_atm * term_a) &
                        &/ theta_bar_0_tfc(i, j)

                    if(balance_eq == "QG") then
                      stop "ERROR: balance_eq == QG not provided"
                    else
                      streamfunc = g * F_a * term_b / f_Coriolis_dim
                      pi_pr_xz_tfc(i, j, k) = f_Coriolis_dim * streamfunc &
                          &/ (cp * thetaRef * thetaStratTFC(i, j, k))

                      if(balance_eq == "QG") then
                        stop "ERROR: balance_eq == QG not provided"
                      else
                        var%pi(i, j, k) = pi_pr_xz_tfc(i, j, k)
                      end if
                    end if
                  else
                    y_dim = y(j + j00) * lRef

                    if(y_dim < Ly_end * (0.25 - 0.5 * bar_sigma_y)) then
                      term_a = 1.0
                    else if(y_dim >= Ly_end * (0.25 - 0.5 * bar_sigma_y) .and. &
                        &y_dim < Ly_end * (0.25 + 0.5 * bar_sigma_y)) then
                      term_a = 0.5 * (1.0 - sin(pi * ((y_dim - 0.25 * Ly_end) &
                          &/ (Ly_end * bar_sigma_y))))
                    else if(y_dim >= Ly_end * (0.25 + 0.5 * bar_sigma_y) .and. &
                        &y_dim < Ly_end * (0.75 - 0.5 * bar_sigma_y)) then
                      term_a = 0.0
                    else if(y_dim >= Ly_end * (0.75 - 0.5 * bar_sigma_y) .and. &
                        &y_dim < Ly_end * (0.75 + 0.5 * bar_sigma_y)) then
                      term_a = 0.5 * (1.0 + sin(pi * ((y_dim - 0.75 * Ly_end) &
                          &/ (Ly_end * bar_sigma_y))))
                    else
                      term_a = 1.0
                    end if

                    term_b = (dTh_atm - 2.0 * dTh_atm * term_a) &
                        &/ theta_bar_0_tfc(i, j)

                    if(balance_eq == "QG") then
                      stop "ERROR: balance_eq == QG not provided"
                    else
                      streamfunc = g * F_a * term_b / f_Coriolis_dim

                      pi_pr_yz_tfc(i, j, k) = f_Coriolis_dim * streamfunc &
                          &/ (cp * thetaRef * thetaStratTFC(i, j, k))

                      if(balance_eq == "QG") then
                        stop "ERROR: balance_eq == QG not provided"
                      else
                        var%pi(i, j, k) = pi_pr_yz_tfc(i, j, k)
                      end if
                    end if
                  end if
                end do
              end do
            end do
          else
            do k = 0, nz + 1
              z_dim = z(k) * lRef

              if(z_dim .lt. z_2) then ! Troposphere
                F_0(k) = alpha_t * (2. * z_dim * sin(pi * z_dim / (2. * z_2)) &
                    &- z_dim ** 2 / z_2)
              else
                F_0(k) = alpha_t * (2. * (z_dim - L_z) * sin(pi * (z_dim &
                    &- L_z) / (2. * (z_2 - L_z))) - (z_dim - L_z) ** 2 / (z_2 &
                    &- L_z)) * z_2 / (z_2 - L_z)
              end if

              F_a = F_0(k)

              ! below also enables test of 2D dynamics in x-z plane

              if(sizeX > 1 .and. sizeY == 1) then
                do i = 1, nx
                  x_dim = x(i + i00) * lRef

                  if(x_dim < Lx_end * (0.25 - 0.5 * bar_sigma_x)) then
                    term_a = 1.0
                  else if(x_dim >= Lx_end * (0.25 - 0.5 * bar_sigma_x) .and. &
                      &x_dim < Lx_end * (0.25 + 0.5 * bar_sigma_x)) then
                    term_a = 0.5 * (1.0 - sin(pi * ((x_dim - 0.25 * Lx_end) &
                        &/ (Lx_end * bar_sigma_x))))
                  else if(x_dim >= Lx_end * (0.25 + 0.5 * bar_sigma_x) .and. &
                      &x_dim < Lx_end * (0.75 - 0.5 * bar_sigma_x)) then
                    term_a = 0.0
                  else if(x_dim >= Lx_end * (0.75 - 0.5 * bar_sigma_x) .and. &
                      &x_dim < Lx_end * (0.75 + 0.5 * bar_sigma_x)) then
                    term_a = 0.5 * (1.0 + sin(pi * ((x_dim - 0.75 * Lx_end) &
                        &/ (Lx_end * bar_sigma_x))))
                  else
                    term_a = 1.0
                  end if

                  term_b = (dTh_atm - 2. * dTh_atm * term_a) / theta_bar_0

                  if(balance_eq == 'QG') then
                    stop 'ERROR: balance_eq == QG not provided'
                  else
                    streamfunc = g * F_a * term_b / (f_Coriolis_dim)

                    pi_pr_xz(i, k) = f_Coriolis_dim * streamfunc / (cp &
                        &* thetaRef * thetaStrat(k))

                    if(balance_eq == 'QG') then
                      stop 'ERROR: balance_eq == QG not provided'
                    else
                      var%pi(i, :, k) = pi_pr_xz(i, k)
                    end if
                  end if
                end do
              else
                do j = 1, ny
                  y_dim = y(j + j00) * lRef

                  if(y_dim < Ly_end * (0.25 - 0.5 * bar_sigma_y)) then
                    term_a = 1.0
                  else if(y_dim >= Ly_end * (0.25 - 0.5 * bar_sigma_y) .and. &
                      &y_dim < Ly_end * (0.25 + 0.5 * bar_sigma_y)) then
                    term_a = 0.5 * (1.0 - sin(pi * ((y_dim - 0.25 * Ly_end) &
                        &/ (Ly_end * bar_sigma_y))))
                  else if(y_dim >= Ly_end * (0.25 + 0.5 * bar_sigma_y) .and. &
                      &y_dim < Ly_end * (0.75 - 0.5 * bar_sigma_y)) then
                    term_a = 0.0
                  else if(y_dim >= Ly_end * (0.75 - 0.5 * bar_sigma_y) .and. &
                      &y_dim < Ly_end * (0.75 + 0.5 * bar_sigma_y)) then
                    term_a = 0.5 * (1.0 + sin(pi * ((y_dim - 0.75 * Ly_end) &
                        &/ (Ly_end * bar_sigma_y))))
                  else
                    term_a = 1.0
                  end if

                  term_b = (dTh_atm - 2. * dTh_atm * term_a) / theta_bar_0

                  if(balance_eq == 'QG') then
                    stop 'ERROR: balance_eq == QG not provided'
                  else
                    streamfunc = g * F_a * term_b / (f_Coriolis_dim)

                    pi_pr_yz(j, k) = f_Coriolis_dim * streamfunc / (cp &
                        &* thetaRef * thetaStrat(k))

                    if(balance_eq == 'QG') then
                      stop 'ERROR: balance_eq == QG not provided'
                    else
                      var%pi(:, j, k) = pi_pr_yz(j, k)
                    end if
                  end if
                end do
              end if
            end do
          end if

          ! determine density from the Exner-pressure fluctuations

          if(topography) then
            ! TFC FJ
            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  var%rho(i, j, k) = - cp * thetaRef / (g * lRef * dz) * 0.25 &
                      &* ((pStratTFC(i, j, k) / jac(i, j, k) + pStratTFC(i, j, &
                      &k - 1) / jac(i, j, k - 1)) * (var%pi(i, j, k) &
                      &- var%pi(i, j, k - 1)) + (pStratTFC(i, j, k) / jac(i, &
                      &j, k) + pStratTFC(i, j, k + 1) / jac(i, j, k + 1)) &
                      &* (var%pi(i, j, k + 1) - var%pi(i, j, k)))
                end do
              end do
            end do
          else
            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  ! density
                  var%rho(i, j, k) = - cp * thetaRef / (g * lRef * dz) * 0.5 &
                      &* (PstratTilde(k - 1) * (var%pi(i, j, k) - var%pi(i, j, &
                      &k - 1)) + PstratTilde(k) * (var%pi(i, j, k + 1) &
                      &- var%pi(i, j, k)))

                end do
              end do
            end do
          end if

          call setHalos(var, "var")
          call setBoundary(var, flux, "var")

          ! determine horizontal wind from density and Exner-pressure
          ! fluctuations

          ! below also enables test of 2D dynamics in x-z plane

          if(sizeX > 1 .and. sizeY == 1) then
            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  ! v-wind

                  rho_int_m0 = 0.5 * (var%rho(i - 1, j, k) + var%rho(i, j, k)) &
                      &* rhoRef
                  rho_int_00 = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k)) &
                      &* rhoRef
                  rho_int_mp = 0.5 * (var%rho(i - 1, j + 1, k) + var%rho(i, j &
                      &+ 1, k)) * rhoRef
                  rho_int_0p = 0.5 * (var%rho(i, j + 1, k) + var%rho(i + 1, j &
                      &+ 1, k)) * rhoRef

                  rho = var%rho(i, j, k)

                  if(topography) then
                    ! TFC FJ
                    rho_int_m0 = rho_int_m0 + 0.5 * (rhoStratTFC(i, j, k) &
                        &+ rhoStratTFC(i - 1, j, k)) * rhoRef
                    rho_int_00 = rho_int_00 + 0.5 * (rhoStratTFC(i, j, k) &
                        &+ rhoStratTFC(i + 1, j, k)) * rhoRef
                    rho_int_mp = rho_int_mp + 0.5 * (rhoStratTFC(i, j + 1, k) &
                        &+ rhoStratTFC(i - 1, j + 1, k)) * rhoRef
                    rho_int_0p = rho_int_0p + 0.5 * (rhoStratTFC(i, j + 1, k) &
                        &+ rhoStratTFC(i + 1, j + 1, k)) * rhoRef
                    rho = rho + rhoStratTFC(i, j, k) * rhoRef
                  else
                    rho_int_m0 = rho_int_m0 + rhoStrat(k) * rhoRef
                    rho_int_00 = rho_int_00 + rhoStrat(k) * rhoRef
                    rho_int_mp = rho_int_mp + rhoStrat(k) * rhoRef
                    rho_int_0p = rho_int_0p + rhoStrat(k) * rhoRef

                    rho = rho + rhoStrat(k) * rhoRef
                  end if

                  if(topography) then
                    ! TFC FJ
                    ! Compute values at cell edges.
                    pEdgeR = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                        &+ pStratTFC(i + 1, j, k) / jac(i + 1, j, k))
                    pEdgeL = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                        &+ pStratTFC(i - 1, j, k) / jac(i - 1, j, k))
                    pFEdgeR = 0.5 * (pStratTFC(i, j + 1, k) / jac(i, j + 1, k) &
                        &+ pStratTFC(i + 1, j + 1, k) / jac(i + 1, j + 1, k))
                    pFEdgeL = 0.5 * (pStratTFC(i, j + 1, k) / jac(i, j + 1, k) &
                        &+ pStratTFC(i - 1, j + 1, k) / jac(i - 1, j + 1, k))
                    ! Compute pressure gradient components.
                    piUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, &
                        &3) * var%pi(i, j, k + 1) + jac(i + 1, j, k + 1) &
                        &* met(i + 1, j, k + 1, 1, 3) * var%pi(i + 1, j, k + 1))
                    piDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, &
                        &3) * var%pi(i, j, k - 1) + jac(i + 1, j, k - 1) &
                        &* met(i + 1, j, k - 1, 1, 3) * var%pi(i + 1, j, k - 1))
                    piUEdgeL = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, &
                        &3) * var%pi(i, j, k + 1) + jac(i - 1, j, k + 1) &
                        &* met(i - 1, j, k + 1, 1, 3) * var%pi(i - 1, j, k + 1))
                    piDEdgeL = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, &
                        &3) * var%pi(i, j, k - 1) + jac(i - 1, j, k - 1) &
                        &* met(i - 1, j, k - 1, 1, 3) * var%pi(i - 1, j, k - 1))
                    piFUEdgeR = 0.5 * (jac(i, j + 1, k + 1) * met(i, j + 1, k &
                        &+ 1, 1, 3) * var%pi(i, j + 1, k + 1) + jac(i + 1, j &
                        &+ 1, k + 1) * met(i + 1, j + 1, k + 1, 1, 3) &
                        &* var%pi(i + 1, j + 1, k + 1))
                    piFDEdgeR = 0.5 * (jac(i, j + 1, k - 1) * met(i, j + 1, k &
                        &- 1, 1, 3) * var%pi(i, j, k - 1) + jac(i + 1, j + 1, &
                        &k - 1) * met(i + 1, j + 1, k - 1, 1, 3) * var%pi(i &
                        &+ 1, j + 1, k - 1))
                    piFUEdgeL = 0.5 * (jac(i, j + 1, k + 1) * met(i, j + 1, k &
                        &+ 1, 1, 3) * var%pi(i, j + 1, k + 1) + jac(i - 1, j &
                        &+ 1, k + 1) * met(i - 1, j + 1, k + 1, 1, 3) &
                        &* var%pi(i - 1, j + 1, k + 1))
                    piFDEdgeL = 0.5 * (jac(i, j + 1, k - 1) * met(i, j + 1, k &
                        &- 1, 1, 3) * var%pi(i, j + 1, k - 1) + jac(i - 1, j &
                        &+ 1, k - 1) * met(i - 1, j + 1, k - 1, 1, 3) &
                        &* var%pi(i - 1, j + 1, k - 1))
                    piGrad = 0.25 * (pEdgeR / rho_int_00 * ((jac(i + 1, j, k) &
                        &* var%pi(i + 1, j, k) - jac(i, j, k) * var%pi(i, j, &
                        &k)) / dx + (piUEdgeR - piDEdgeR) * 0.5 / dz) + pEdgeL &
                        &/ rho_int_m0 * ((jac(i, j, k) * var%pi(i, j, k) &
                        &- jac(i - 1, j, k) * var%pi(i - 1, j, k)) / dx &
                        &+ (piUEdgeL - piDEdgeL) * 0.5 / dz) + pFEdgeR &
                        &/ rho_int_0p * ((jac(i + 1, j + 1, k) * var%pi(i + 1, &
                        &j + 1, k) - jac(i, j + 1, k) * var%pi(i, j + 1, k)) &
                        &/ dx + (piFUEdgeR - piFDEdgeR) * 0.5 / dz) + pFEdgeL &
                        &/ rho_int_mp * ((jac(i, j + 1, k) * var%pi(i, j + 1, &
                        &k) - jac(i - 1, j + 1, k) * var%pi(i - 1, j + 1, k)) &
                        &/ dx + (piFUEdgeL - piFDEdgeL) * 0.5 / dz))
                    ! Compute meridional wind.
                    var%v(i, j, k) = 1.0 / f_Coriolis_dim * cp * thetaRef &
                        &* rhoRef / (uRef * lRef) * piGrad
                  else
                    var%v(i, j, k) = 1.0 / f_Coriolis_dim * cp * thetaRef &
                        &* rhoRef / (uRef * lRef * dx) * Pstrat(k) * 0.25 &
                        &* ((var%pi(i, j, k) - var%pi(i - 1, j, k)) &
                        &/ rho_int_m0 + (var%pi(i + 1, j, k) - var%pi(i, j, &
                        &k)) / rho_int_00 + (var%pi(i, j + 1, k) - var%pi(i &
                        &- 1, j + 1, k)) / rho_int_mp + (var%pi(i + 1, j + 1, &
                        &k) - var%pi(i, j + 1, k)) / rho_int_0p)
                  end if

                  ! store env pot temp: local
                  if(topography) then
                    ! TFC FJ
                    the_env_pp(i, j, k) = pStratTFC(i, j, k) / rho * rhoRef
                  else
                    the_env_pp(i, j, k) = Pstrat(k) / rho * rhoRef
                  end if

                  ! store env dens: local
                  dens_env_pp(i, j, k) = rho / rhoRef

                  ! for sponge.
                  u_env_pp(i, j, k) = 0.
                  v_env_pp(i, j, k) = var%v(i, j, k)

                  if(timeScheme /= "semiimplicit" .and. .not. auxil_equ) then
                    ! Here it is pot temp (re-calculated below)
                    var%rhop(i, j, k) = the_env_pp(i, j, k)
                  end if
                end do
              end do
            end do

          else
            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  ! u-wind

                  rho_int_0m = 0.5 * (var%rho(i, j - 1, k) + var%rho(i, j, k)) &
                      &* rhoRef
                  rho_int_00 = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k)) &
                      &* rhoRef
                  rho_int_pm = 0.5 * (var%rho(i + 1, j - 1, k) + var%rho(i &
                      &+ 1, j, k)) * rhoRef
                  rho_int_p0 = 0.5 * (var%rho(i + 1, j, k) + var%rho(i + 1, j &
                      &+ 1, k)) * rhoRef

                  rho = var%rho(i, j, k) * rhoRef

                  if(topography) then
                    ! TFC FJ
                    rho_int_0m = rho_int_0m + 0.5 * (rhoStratTFC(i, j, k) &
                        &+ rhoStratTFC(i, j - 1, k)) * rhoRef
                    rho_int_00 = rho_int_00 + 0.5 * (rhoStratTFC(i, j, k) &
                        &+ rhoStratTFC(i, j + 1, k)) * rhoRef
                    rho_int_pm = rho_int_pm + 0.5 * (rhoStratTFC(i + 1, j, k) &
                        &+ rhoStratTFC(i + 1, j - 1, k)) * rhoRef
                    rho_int_p0 = rho_int_p0 + 0.5 * (rhoStratTFC(i + 1, j, k) &
                        &+ rhoStratTFC(i + 1, j + 1, k)) * rhoRef
                    rho = rho + rhoStratTFC(i, j, k) * rhoRef
                  else
                    rho_int_0m = rho_int_0m + rhoStrat(k) * rhoRef
                    rho_int_00 = rho_int_00 + rhoStrat(k) * rhoRef
                    rho_int_pm = rho_int_pm + rhoStrat(k) * rhoRef
                    rho_int_p0 = rho_int_p0 + rhoStrat(k) * rhoRef

                    rho = rho + rhoStrat(k) * rhoRef
                  end if

                  if(topography) then
                    ! TFC FJ
                    ! Compute values at cell edges.
                    pEdgeF = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                        &+ pStratTFC(i, j + 1, k) / jac(i, j + 1, k))
                    pEdgeB = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                        &+ pStratTFC(i, j - 1, k) / jac(i, j - 1, k))
                    pREdgeF = 0.5 * (pStratTFC(i + 1, j, k) / jac(i + 1, j, k) &
                        &+ pStratTFC(i + 1, j + 1, k) / jac(i + 1, j + 1, k))
                    pREdgeB = 0.5 * (pStratTFC(i + 1, j, k) / jac(i + 1, j, k) &
                        &+ pStratTFC(i + 1, j - 1, k) / jac(i + 1, j - 1, k))
                    ! Compute pressure gradient component.
                    piUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, &
                        &3) * var%pi(i, j, k + 1) + jac(i, j + 1, k + 1) &
                        &* met(i, j + 1, k + 1, 2, 3) * var%pi(i, j + 1, k + 1))
                    piDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, &
                        &3) * var%pi(i, j, k - 1) + jac(i, j + 1, k - 1) &
                        &* met(i, j + 1, k - 1, 2, 3) * var%pi(i, j + 1, k - 1))
                    piUEdgeB = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, &
                        &3) * var%pi(i, j, k + 1) + jac(i, j - 1, k + 1) &
                        &* met(i, j - 1, k + 1, 2, 3) * var%pi(i, j - 1, k + 1))
                    piDEdgeB = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, &
                        &3) * var%pi(i, j, k - 1) + jac(i, j - 1, k - 1) &
                        &* met(i, j - 1, k - 1, 2, 3) * var%pi(i, j - 1, k - 1))
                    piRUEdgeF = 0.5 * (jac(i + 1, j, k + 1) * met(i + 1, j, k &
                        &+ 1, 2, 3) * var%pi(i + 1, j, k + 1) + jac(i + 1, j &
                        &+ 1, k + 1) * met(i + 1, j + 1, k + 1, 2, 3) &
                        &* var%pi(i + 1, j + 1, k + 1))
                    piRDEdgeF = 0.5 * (jac(i + 1, j, k - 1) * met(i + 1, j, k &
                        &- 1, 2, 3) * var%pi(i + 1, j, k - 1) + jac(i + 1, j &
                        &+ 1, k - 1) * met(i + 1, j + 1, k - 1, 2, 3) &
                        &* var%pi(i + 1, j + 1, k - 1))
                    piRUEdgeB = 0.5 * (jac(i + 1, j, k + 1) * met(i + 1, j, k &
                        &+ 1, 2, 3) * var%pi(i + 1, j, k + 1) + jac(i + 1, j &
                        &- 1, k + 1) * met(i + 1, j - 1, k + 1, 2, 3) &
                        &* var%pi(i + 1, j - 1, k + 1))
                    piRDEdgeB = 0.5 * (jac(i + 1, j, k - 1) * met(i + 1, j, k &
                        &- 1, 2, 3) * var%pi(i + 1, j, k - 1) + jac(i + 1, j &
                        &- 1, k - 1) * met(i + 1, j - 1, k - 1, 2, 3) &
                        &* var%pi(i + 1, j - 1, k - 1))
                    piGrad = 0.25 * (pEdgeF / rho_int_00 * ((jac(i, j + 1, k) &
                        &* var%pi(i, j + 1, k) - jac(i, j, k) * var%pi(i, j, &
                        &k)) / dy + (piUEdgeF - piDEdgeF) * 0.5 / dz) + pEdgeB &
                        &/ rho_int_0m * ((jac(i, j, k) * var%pi(i, j, k) &
                        &- jac(i, j - 1, k) * var%pi(i, j - 1, k)) / dy &
                        &+ (piUEdgeB - piDEdgeB) * 0.5 / dz) + pREdgeF &
                        &/ rho_int_p0 * ((jac(i + 1, j + 1, k) * var%pi(i + 1, &
                        &j + 1, k) - jac(i + 1, j, k) * var%pi(i + 1, j, k)) &
                        &/ dy + (piRUEdgeF - piRDEdgeF) * 0.5 / dz) + pREdgeB &
                        &/ rho_int_pm * ((jac(i + 1, j, k) * var%pi(i + 1, j, &
                        &k) - jac(i + 1, j - 1, k) * var%pi(i + 1, j - 1, k)) &
                        &/ dy + (piRUEdgeB - piRDEdgeB) * 0.5 / dz))
                    ! Compute zonal wind.
                    var%u(i, j, k) = - 1.0 / f_Coriolis_dim * cp * thetaRef &
                        &* rhoRef / (uRef * lRef) * piGrad
                  else
                    var%u(i, j, k) = - 1.0 / f_Coriolis_dim * cp * thetaRef &
                        &* rhoRef / (uRef * lRef * dy) * Pstrat(k) * 0.25 &
                        &* ((var%pi(i, j, k) - var%pi(i, j - 1, k)) &
                        &/ rho_int_0m + (var%pi(i, j + 1, k) - var%pi(i, j, &
                        &k)) / rho_int_00 + (var%pi(i + 1, j, k) - var%pi(i &
                        &+ 1, j - 1, k)) / rho_int_pm + (var%pi(i + 1, j + 1, &
                        &k) - var%pi(i + 1, j, k)) / rho_int_p0)
                  end if

                  ! store env pot temp: local
                  if(topography) then
                    ! TFC FJ
                    the_env_pp(i, j, k) = pStratTFC(i, j, k) / rho * rhoRef
                  else
                    the_env_pp(i, j, k) = Pstrat(k) / rho * rhoRef
                  end if

                  ! store env dens: local
                  dens_env_pp(i, j, k) = rho / rhoRef

                  ! for sponge.
                  u_env_pp(i, j, k) = var%u(i, j, k)
                  v_env_pp(i, j, k) = 0.

                  if(timeScheme /= "semiimplicit" .and. .not. auxil_equ) then
                    ! Here it is pot temp (re-calculated below)
                    var%rhop(i, j, k) = the_env_pp(i, j, k)
                  end if
                end do
              end do
            end do

          end if

          ! // TODO this loop probably not needed?
          do j = 1, ny
            u_const(0:nx, j) = var%u(0:nx, j, nz)
          end do
        else if(init_bal == "hydrost") then
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                var%rho(i, j, k) = 0. ! density

                var%u(i, j, k) = 0.0 ! u
                var%v(i, j, k) = 0.0 ! v
                var%w(i, j, k) = 0.0 ! w
                var%pi(i, j, k) = 0.0 ! deviation of Exner pressure

                if(timeScheme /= "semiimplicit" .and. .not. auxil_equ) then
                  var%rhop(i, j, k) = thetaStrat(k) ! pot temp
                end if
              end do
            end do
          end do
        else
          stop "initialize: init_bal not def. for this model."
        end if
      end if

      do k = 1, nz - 1
        do j = 1, ny
          do i = 1, nx
            br_vais_sq(i, j, k) = g * (the_env_pp(i, j, k + 1) - the_env_pp(i, &
                &j, k)) / (dz * the_env_pp(i, j, k) * lRef)
            ! TFC FJ
            if(topography) then
              br_vais_sq(i, j, k) = br_vais_sq(i, j, k) / jac(i, j, k)
            end if
          end do
        end do
      end do

      !------------------------------------------
      !        Output of background \theta
      !------------------------------------------

      if(output_theta_bgr) then
        call output_background(the_env_pp, nz, 'theta_bgr.dat', thetaRef)
      end if

      if(output_rho_bgr) then
        call output_background(dens_env_pp, nz, 'rho_bgr.dat', rhoRef)
      end if

      !------------------------------------------
      !        Output of background \theta
      !------------------------------------------

      if(output_br_vais_sq) then
        call output_background(br_vais_sq, nz - 1, 'Br_Va_bgr.dat', 1.)
      end if

      if(topography) then
        ! TFC FJ
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              balance(i, j, k) = uRef * (var%u(i, j, k + 1) - var%u(i, j, k)) &
                  &/ (jac(i, j, k) * dz * lRef) + g / (jac(i, j, k) &
                  &* f_Coriolis_dim * lRef * thetaStratTFC(i, j, k)) &
                  &* ((jac(i, j + 1, k) * var%rhop(i, j + 1, k) - jac(i, j, k) &
                  &* var%rhop(i, j, k)) / dy + (jac(i, j, k + 1) * met(i, j, k &
                  &+ 1, 2, 3) * var%rhop(i, j, k + 1) - jac(i, j, k) * met(i, &
                  &j, k, 2, 3) * var%rhop(i, j, k)) / dz)

              balance1(i, j, k) = f_Coriolis_dim * uRef * var%u(i, j, k) &
                  &+ thetaRef * Rsp * var%rhop(i, j, k) + ((jac(i, j + 1, k) &
                  &* var%pi(i, j + 1, k) - jac(i, j, k) * var%pi(i, j, k)) &
                  &/ dy + (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                  &* var%pi(i, j, k + 1) - jac(i, j, k) * met(i, j, k, 2, 3) &
                  &* var%pi(i, j, k)) / dz) / (jac(i, j, k) * kappa * lRef)

              balance2(i, j, k) = thetaRef * Rsp * var%rhop(i, j, k) &
                  &* (var%pi(i, j, k + 1) - var%pi(i, j, k)) / (jac(i, j, k) &
                  &* kappa * dz * lRef) - g * (var%rhop(i, j, k) &
                  &- thetaStratTFC(i, j, k)) / thetaStratTFC(i, j, k)

              balance3(i, j, k) = g * (var%rhop(i, j, k) - thetaStratTFC(i, j, &
                  &k)) / thetaStratTFC(i, j, k)

              balance4(i, j, k) = thetaRef * Rsp * var%rhop(i, j, k) &
                  &* (var%pi(i, j, k + 1) - var%pi(i, j, k)) / (jac(i, j, k) &
                  &* kappa * dz * lRef)
            end do
          end do
        end do
      else
        do k = 1, nz - 1
          do j = 1, ny
            do i = 1, nx
              balance(i, j, k) = uRef * (var%u(i, j, k + 1) - var%u(i, j, k)) &
                  &/ (dz * lRef) + g * (var%rhop(i, j + 1, k) - var%rhop(i, j, &
                  &k)) / (f_Coriolis_dim * dy * lRef * thetaStrat(k))
            end do
          end do
        end do

        do k = 1, nz - 1
          do j = 1, ny
            do i = 1, nx
              balance1(i, j, k) = f_Coriolis_dim * uRef * var%u(i, j, k) &
                  &+ thetaRef * Rsp * var%rhop(i, j, k) * (var%pi(i, j + 1, k) &
                  &- var%pi(i, j, k)) / (kappa * dy * lRef)

              balance3(i, j, k) = g * (var%rhop(i, j, k) - thetaStrat(k)) &
                  &/ thetaStrat(k)

              balance4(i, j, k) = thetaRef * Rsp * var%rhop(i, j, k) &
                  &* (var%pi(i, j, k + 1) - var%pi(i, j, k)) / (kappa * dz &
                  &* lRef)

              balance2(i, j, k) = thetaRef * Rsp * var%rhop(i, j, k) &
                  &* (var%pi(i, j, k + 1) - var%pi(i, j, k)) / (kappa * dz &
                  &* lRef) - g * (var%rhop(i, j, k) - thetaStrat(k)) &
                  &/ thetaStrat(k)
            end do
          end do
        end do
      end if

      call output_background(balance, nz - 1, 'balance.dat', 1.)
      call output_background(balance1, nz - 1, 'balance1.dat', 1.)
      call output_background(balance2, nz - 1, 'balance2.dat', 1.)
      call output_background(balance3, nz - 1, 'balance3.dat', 1.)
      call output_background(balance4, nz - 1, 'balance4.dat', 1.)

      !------------------------------------------
      !        Adding the noise
      !------------------------------------------

      ! noise added to density fluctuations in a relative mainer

      noise_var = "none"
      noise_mag = 1.0

      if(add_noise) then
        call Random_Seed()
        call noise_array(noise_mag, noise_var, noise)

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              var%rho(i, j, k) = var%rho(i, j, k) * (1.0 + noise(i, j, k))
            end do
          end do
        end do
      end if

    case("densitySineTransport")

      rho0 = 1.0; u0 = 1.0; v0 = 0.0; w0 = 0.0

      ! density sine curve in x
      do i = - nbx, nx + nbx
        var%rho(i, :, :) = rho0 + 0.1 * sin(2 * pi * x(i))
      end do
      ! constant velocity vector in x-direction
      var%u(:, :, :) = u0
      var%v(:, :, :) = v0
      var%w(:, :, :) = w0

      ! constant Exner pressure
      var%pi(:, :, :) = p0

      ! ----------------------------------------------------------------

    case('densityDiscXYZ')
      ! center and radius of disc
      x0 = 0.5; y0 = 0.5; z0 = 0.5; r0 = 1. / 3.
      u0 = 1.0; v0 = 1.0; w0 = 1.0 ! constant backround velocity
      rhoDisc = 0.1 ! density of disc

      ! density
      var%rho(:, :, :) = 1.0
      do i = 1, nx
        do j = 1, ny
          do k = 1, nz
            delX = x(i) - x0
            delY = y(j) - y0
            delZ = z(k) - z0
            r = sqrt(delX ** 2 + delY ** 2 + delZ ** 2)
            if(r >= r0) var%rho(i, j, k) = var%rho(i, j, k) + rhoDisc
          end do
        end do
      end do

      ! horizontal velocity u
      var%u(:, :, :) = u0

      ! horizontal velocity v
      var%v(:, :, :) = v0

      ! vertical velocity w
      var%w(:, :, :) = w0
      var%w(:, :, - nbz:2) = 0.0 ! zero near walls
      var%w(:, :, nz - 2:nz + nbz) = 0.0

      ! pressure variable pi'
      var%pi(:, :, :) = 0.0

      ! ----------------------------------------------------------------

    case('densityDiscXZ')
      ! transport of a density disc with constant speed in x
      x0 = 0.5; z0 = 0.5; r0 = 1. / 3. ! center and radius of disc
      u0 = 1.0; v0 = 0.0; w0 = 0.0 ! constant backround velocity
      rhoDisc = 0.1 ! density of disc

      ! density
      var%rho(:, :, :) = 1.0
      do i = 1, nx
        do j = 1, ny
          do k = 1, nz
            delX = x(i) - x0
            delZ = z(k) - z0
            r = sqrt(delX ** 2 + delZ ** 2)
            if(r >= r0) var%rho(i, j, k) = var%rho(i, j, k) + rhoDisc
          end do
        end do
      end do

      ! horizontal velocity u
      var%u(:, :, :) = u0

      ! horizontal velocity v
      var%v(:, :, :) = v0

      ! vertical velocity w
      var%w(:, :, :) = w0
      if(w0 /= 0.0) then
        var%w(:, :, - nbz:2) = 0.0 ! zero near walls
        var%w(:, :, nz - 2:nz + nbz) = 0.0
      end if

      ! pressure variable pi'
      var%pi(:, :, :) = 0.0

      ! ----------------------------------------------------------------

    case('densityDiscXY')
      x0 = 0.5; y0 = 0.5; r0 = 1. / 3. ! center and radius of disc
      u0 = 1.0; v0 = 1.0; w0 = 0.0 ! constant backround velocity
      rhoDisc = 0.1 ! density of disc

      ! density
      var%rho(:, :, :) = 1.0
      do i = 1, nx
        do j = 1, ny
          do k = 1, nz
            delX = x(i) - x0
            delY = y(j) - y0
            r = sqrt(delX ** 2 + delY ** 2)
            if(r >= r0) var%rho(i, j, k) = var%rho(i, j, k) + rhoDisc
          end do
        end do
      end do

      ! horizontal velocity u
      var%u(:, :, :) = u0

      ! horizontal velocity v
      var%v(:, :, :) = v0

      ! vertical velocity w
      var%w(:, :, :) = w0

      ! pressure variable pi'
      var%pi(:, :, :) = 0.0

      ! -----------------------------------------------------------------

    case("greshoVortexXY")
      updateMass = .true.
      updateIce = .true.

      ! setup
      x0 = 0.5; y0 = 0.5; z0 = 0.5 ! center vortex
      uVortex = 1.0;
      r0 = 1. / 3.
      ! background
      pInf = 10.0
      rho0 = 1.0; u0 = 0.0; v0 = 0.0; w0 = 0.0

      ! density
      var%rho(:, :, :) = rho0

      ! velocities: background
      var%u(:, :, :) = u0; var%v(:, :, :) = v0; var%w(:, :, :) = w0

      ! vortex velocity u
      do i = 0, nx
        do j = 1, ny
          do k = 1, nz
            xu = x(i) + 0.5 * dx
            yu = y(j)
            delX = xu - x0
            delY = yu - y0
            r = sqrt(delX ** 2 + delY ** 2)
            uPhi = greshoVelocity(uVortex, r / r0)
            var%u(i, j, k) = var%u(i, j, k) - uPhi * delY / r
          end do
        end do
      end do

      ! vortex velocity v
      do i = 1, nx
        do j = 0, ny
          do k = 1, nz
            xv = x(i)
            yv = y(j) + 0.5 * dy
            delX = xv - x0
            delY = yv - y0
            r = sqrt(delX ** 2 + delY ** 2)
            uPhi = greshoVelocity(uVortex, r / r0)
            var%v(i, j, k) = var%v(i, j, k) + uPhi * delX / r
          end do
        end do
      end do

      ! pressure variable
      ! set pi'(inf) including ghost cells
      var%pi(:, :, :) = kappaInv * pInf ** kappa
      do i = 1, nx
        do j = 1, ny
          do k = 1, nz
            delX = x(i) - x0
            delY = y(j) - y0
            r = sqrt(delX ** 2 + delY ** 2)
            p = greshoPressure(pInf, uVortex, r / r0)

            ! pressure variable pi' (ref. Klein)
            var%pi(i, j, k) = kappaInv * (p ** kappa - pInf ** kappa)
          end do
        end do
      end do

      ! -----------------------------------------------------------------

    case("greshoVortexXZ")
      updateMass = .true.
      updateIce = .true.

      ! setup
      x0 = 0.0; z0 = 0.5 ! center vortex
      uVortex = 1.0;
      r0 = 1. / 3.
      ! background
      pInf = 10.0
      rho0 = 1.0; u0 = 1.0; v0 = 0.0; w0 = 0.0

      var%rho(:, :, :) = rho0

      ! velocities: background
      var%u(:, :, :) = u0; var%v(:, :, :) = v0; var%w(:, :, :) = w0

      ! vortex velocity u
      do i = 0, nx
        do j = 1, ny
          do k = 1, nz
            xu = x(i) + 0.5 * dx
            zu = z(k)
            delX = xu - x0
            delZ = zu - z0
            r = sqrt(delX ** 2 + delZ ** 2)
            uPhi = greshoVelocity(uVortex, r / r0)
            var%u(i, j, k) = var%u(i, j, k) - uPhi * delZ / r
          end do
        end do
      end do

      ! vortex velocity w
      do i = 1, nx
        do j = 1, ny
          do k = 0, nz
            xw = x(i)
            zw = z(k) + 0.5 * dz
            delX = xw - x0
            delZ = zw - z0
            r = sqrt(delX ** 2 + delZ ** 2)
            uPhi = greshoVelocity(uVortex, r / r0)
            var%w(i, j, k) = var%w(i, j, k) + uPhi * delX / r
          end do
        end do
      end do

      ! pressure variable
      ! set pi'(inf) including ghost cells
      var%pi(:, :, :) = kappaInv * pInf ** kappa
      do i = 1, nx
        do j = 1, ny
          do k = 1, nz
            delX = x(i) - x0
            delZ = z(k) - z0
            r = sqrt(delX ** 2 + delZ ** 2)
            p = greshoPressure(pInf, uVortex, r / r0)

            ! pressure variable pi' (ref. Klein)
            var%pi(i, j, k) = kappaInv * (p ** kappa - pInf ** kappa)
          end do
        end do
      end do

      !---------------------------------------------------------------------------

    case("projectionTest")
      ! test the projection step

      var%rho(:, :, :) = 1.0 ! constant density

      var%u(:, :, :) = 1.0 ! constant velocity in x

      var%v(:, :, :) = 0.0

      var%w(:, :, :) = 0.0

      ! pressure sine curve
      do i = 0, nx
        var%pi(i, :, :) = 1.0 + 0.01 * sin(2 * pi * x(i))
      end do

      ! -----------------------------------------------------------------

    case("matrixStructureTest")
      ! test the momentum transport at constant density and
      ! pressure along x.

      ! constant density
      var%rho(:, :, :) = 2.0

      ! velocity sine curve along z
      do k = 0, nz
        zw = z(k) + dz / 2.0
        var%w(:, :, k) = sin(2 * pi * zw)
      end do
      var%u(:, :, :) = 0.0
      var%v(:, :, :) = 0.0

      ! constant Exner pressure
      var%pi(:, :, :) = 1.0

      ! ----------------------------------------------------------------

    case("momentumTransportX")
      ! test the momentum transport at constant density and
      ! pressure along x.
      updateMass = .false.
      correctMomentum = .false.
      updateIce = .false.

      ! constant density
      var%rho(:, :, :) = 2.0

      ! velocity sine curve along x
      do i = 0, nx
        xu = x(i) + dx / 2.
        var%u(i, :, :) = 1.0 + 1.0e-5 * sin(2 * pi * xu)
      end do

      ! zero velocity along y and z
      var%v(:, :, :) = 0.0
      var%w(:, :, :) = 0.0

      ! constant Exner pressure
      var%pi(:, :, :) = 1.0

      ! ----------------------------------------------------------------

    case("momentumTransportY")
      ! test the momentum transport at constant density and
      ! pressure along y.
      updateMass = .false.
      correctMomentum = .false.
      updateIce = .false.

      ! constant density
      var%rho(:, :, :) = 2.0

      ! velocity sine curve along y
      do j = 0, ny
        yv = y(j) + dy / 2.0
        var%v(:, j, :) = 50.0 + sin(2 * pi * yv)
      end do

      ! zero velocity along x and z
      var%u(:, :, :) = 0.0
      var%w(:, :, :) = 0.0

      ! constant Exner pressure
      var%pi(:, :, :) = 1.0

      ! ----------------------------------------------------------------

    case("momentumFluxTest")
      ! constant density
      var%rho(:, :, :) = 2.0

      ! constant velocity fields
      var%u(:, :, :) = 1.0
      var%v(:, :, :) = 2.0
      var%w(:, :, :) = 3.0

      ! constant Exner pressure
      var%pi(:, :, :) = 1.0

      ! ----------------------------------------------------------------

    case("densityTransportX")
      ! density sine curve in x
      !       predictMomentum = .false.
      !       correctMomentum = .false.

      u0 = 1.0; v0 = 0.0; w0 = 0.0

      do i = 1, nx
        var%rho(i, :, :) = sin(2 * pi * x(i))
      end do
      ! constant velocity vector in x-direction
      var%u(:, :, :) = u0
      var%v(:, :, :) = v0
      var%w(:, :, :) = w0
      var%pi(:, :, :) = p0

      ! -----------------------------------------------------------------

    case("densityTransportY")
      !       predictMomentum = .false.
      !       correctMomentum = .false.
      ! density sine curve in y

      u0 = 0.0; v0 = 1.0; w0 = 0.0
      do j = 1, ny
        var%rho(:, j, :) = sin(2 * pi * (y(j) - 0.25))
      end do
      ! constant velocity vector in y-direction
      var%u(:, :, :) = u0
      var%v(:, :, :) = v0
      var%w(:, :, :) = w0

      ! constant Exner pressure
      var%pi(:, :, :) = p0

      ! ------------------------------------------------------------

    case("massFluxTest")
      ! constant density
      var%rho(:, :, :) = 1.5

      ! constant velocity field:
      var%u(:, :, :) = 1.0
      var%v(:, :, :) = 2.0
      var%w(:, :, :) = 3.0

      ! constant Exner pressure
      var%pi(:, :, :) = 1.0

    case("standard")

      do k = 1, nz
        var%rho(:, :, k) = z(k)
        var%u(:, :, k) = z(k)
        var%v(:, :, k) = z(k)
        var%w(:, :, k) = z(k)
        var%pi(:, :, k) = z(k)
      end do

    case("xparabel")
      ! rho parabolic in x
      do i = - nbx, nx + nbx
        var%rho(i, :, :) = (x(i) + dx / 2) ** 3 - (x(i) - dx / 2) ** 3
        var%rho(i, :, :) = 1.0 / 3.0 / dx * var%rho(i, :, :)
      end do

    case("yparabel")
      ! u-velocity parabolic in y
      do j = - nby, ny + nby
        var%u(:, j, :) = (y(j) + dy / 2) ** 3 - (y(j) - dy / 2) ** 3
        var%u(:, j, :) = 1.0 / 3.0 / dy * var%u(:, j, :)
      end do

    case("zparabel")
      ! u-velocity parabolic in z
      do k = - nbz, nz + nbz
        var%u(:, :, k) = (z(k) + dz / 2) ** 3 - (z(k) - dz / 2) ** 3
        var%u(:, :, k) = 1.0 / 3.0 / dz * var%u(:, :, k)
      end do

    case default

      print *, "init.f90/initialise: testCase = ", testCase
      stop "init.f90/initialise: This testCase is not valid. Stop."

    end select

    ! close input file pinc.f
    close(unit = 10)

    !----------------------------------
    !     Output system settings
    !----------------------------------

    if(master) then ! modified by Junhong Wei (20170216)

      print *, ""
      print *, "Initializing System: "
      print *, "  1) Reference quantities: "

      write(*, fmt = "(a25,f7.3,a)") "rhoRef = ", rhoRef, " kg/m^3"
      write(*, fmt = "(a25,f7.3,a)") "pRef = ", pRef * 1.e-3, " kPa"
      write(*, fmt = "(a25,f7.3,a)") "aRef,uRef = ", uRef, " m/s"
      write(*, fmt = "(a25,f7.3,a)") "lRef = ", lRef * 1.e-3, " km"
      write(*, fmt = "(a25,f7.3,a)") "tRef = ", tRef, " s"
      write(*, fmt = "(a25,f7.3,a)") "thetaRef = ", thetaRef, " K"
      write(*, fmt = "(a25,f7.3,a)") "FRef = ", FRef, " N/m^3"
      print *, ""

      print *, "  2) Non-dimensional numbers: "

      write(*, fmt = "(a25,es10.3)") "Ma = ", Ma
      write(*, fmt = "(a25,es10.3)") "Fr = ", Fr
      write(*, fmt = "(a25,es10.3)") "Re = ", Re
      print *, ""

      print *, "  3) Extreme values: "
      write(*, fmt = "(a25,es10.1,a,f5.1,a)") "PStrat = ", PStrat(nz) * pRef &
          &/ 1000.0, " kPa at z = ", z(nz) * lRef / 1000.0, " km"
      write(*, fmt = "(a25,es10.1,a,f5.1,a)") "rhoStrat = ", rhoStrat(nz) &
          &* rhoRef, " kg/m3 at z = ", z(nz) * lRef / 1000.0, " km"
      write(*, fmt = "(a25,es10.1,a,f5.1,a)") "thetaStrat = ", thetaStrat(nz) &
          &* thetaRef, " K at z = ", z(nz) * lRef / 1000.0, "km"
      print *, ""

      print *, "  4) Constants: "
      write(*, fmt = "(a25,f7.3,a)") "gamma = ", gamma, " "
      write(*, fmt = "(a25,f7.3,a)") "g = ", g, " m/s^2"
      write(*, fmt = "(a25,f7.3,a)") "R_sp = ", Rsp, " J/kg/K"
      write(*, fmt = "(a25,f7.3,a)") "f_Coriolis = ", f_Coriolis_dim, " 1/s"
      write(*, fmt = "(a25,f7.3,a)") "mu_viscous = ", mu_viscous_dim, " m^2/s"
      write(*, fmt = "(a25,f7.3,a)") "mu_conduct = ", mu_conduct_dim, " m^2/s"
      print *, ""

      print *, "  5) Background: "
      select case(background)
      case("isothermal")
        write(*, fmt = "(a25,a15)") "background = ", background
        write(*, fmt = "(a25,f7.3,a)") "T0 = ", Temp0_dim, " K"
        write(*, fmt = "(a25,f7.4,a)") "N = ", sqrt(N2) / tRef, " 1/s"
      case("isentropic")
        write(*, fmt = "(a25,a7)") "background = ", background
        write(*, fmt = "(a25,f7.3,a)") "theta_0 = ", theta0_dim, " K"
      case("const-N")
        write(*, fmt = "(a25,a7)") "background = ", background
        write(*, fmt = "(a25,f7.3,a)") "N = ", N_BruntVaisala_dim, " 1/s"
      case("uniform")
        write(*, fmt = "(a25,a15)") "background = ", background
        write(*, fmt = "(a25,f7.3,a)") "rho0 = ", rhoStrat(1) * rhoRef, " &
            &kg/m^3"
        write(*, fmt = "(a25,f7.1,a)") "theta0 = ", thetaStrat(1) * thetaRef, &
            &" K"
        write(*, fmt = "(a25,f7.2,a)") "N_Brunt-Vaisala = ", NN / tRef, " 1/s"
      end select
      write(*, fmt = "(a25,f7.1,a)") "u0 = ", backgroundFlow(1) * uRef, " m/s"
      write(*, fmt = "(a25,f7.1,a)") "v0 = ", backgroundFlow(2) * uRef, " m/s"
      write(*, fmt = "(a25,f7.1,a)") "w0 = ", backgroundFlow(3) * uRef, " m/s"
      print *, ""

      print *, "  6) Model equations: "
      write(*, fmt = "(a25,a15)") "model = ", model
      write(*, fmt = "(a25)") "inclination angle: "
      write(*, fmt = "(a25,f5.1,a)") "theta = ", vert_theta, " deg"
      write(*, fmt = "(a25,f5.1,a)") "alpha = ", vert_alpha, " deg"
      print *, ""

      print *, "  7) Domain: "
      !write(*,fmt="(a25,a,f6.1,a,f6.1,a)") "[xMin, xMax] = ", &
      !     & "[",lx(0)*lRef/1000.0,"km, ",lx(1)*lRef/1000.0,"km ]"
      print *, "[xMin, xMax] = ", "[", lx(0) * lRef / 1000.0, "km, ", lx(1) &
          &* lRef / 1000.0, "km ]"
      !write(*,fmt="(a25,a,f6.1,a,f6.1,a)") "[yMin, yMax] = ", &
      !     & "[",ly(0)*lRef/1000.0,"km, ",ly(1)*lRef/1000.0,"km ]"
      print *, "[yMin, yMax] = ", "[", ly(0) * lRef / 1000.0, "km, ", ly(1) &
          &* lRef / 1000.0, "km ]"
      !write(*,fmt="(a25,a,f6.1,a,f6.1,a)") "[zMin, zMax] = ", &
      !     & "[",lz(0)*lRef/1000.0,"km, ",lz(1)*lRef/1000.0,"km ]"
      print *, "[zMin, zMax] = ", "[", lz(0) * lRef / 1000.0, "km, ", lz(1) &
          &* lRef / 1000.0, "km ]"
      !write(*,fmt="(a25,i4,a,i4,a,i4)") "nx x ny x nz = ", &
      !     & nx," x ", ny, " x ", nz
      print *, "nx x ny x nz = ", nx, " x ", ny, " x ", nz
      print *, ""

      print *, "  8) Boundary conditions: "
      write(*, fmt = "(a25,a)") "xBoundary = ", trim(xBoundary)
      write(*, fmt = "(a25,a)") "yBoundary = ", trim(yBoundary)
      write(*, fmt = "(a25,a)") "zBoundary = ", trim(zBoundary)
      if(spongeLayer) then
        write(*, fmt = "(a25,a)") "sponge layer = ", "on"
        write(*, fmt = "(a25,f5.1,a)") "height = ", (lz(1) - lz(0)) &
            &* spongeHeight * lRef / 1000.0, " km"
        write(*, fmt = "(a25,es8.1,a)") "relaxation  = ", spongeAlphaZ_dim, " &
            &1/s"
        write(*, fmt = "(a25,es8.1,a)") "relaxation  = ", spongeAlphaZ_fac, " &
            &1/dt"
      else
        write(*, fmt = "(a25,a)") "sponge layer = ", "off"
      end if
      print *, ""

      print *, "  9) Poisson Solver: "
      write(*, fmt = "(a25,a)") "solver = ", poissonSolverType
      write(*, fmt = "(a25,es8.1)") "tolPoisson = ", tolPoisson
      write(*, fmt = "(a25,es8.1)") "tolCond = ", tolCond
      print *, ""

      print *, " 10) Topography: "
      if(topography) then
        write(*, fmt = "(a25,a)") "topography = ", "on"
        write(*, fmt = "(a25,f6.1,a)") "mountain height = ", &
            &mountainHeight_dim, " m"
        write(*, fmt = "(a25,es8.1,a)") "mountain width = ", &
            &mountainWidth_dim, " m"
      else
        write(*, fmt = "(a25,a)") "topography = ", "off"
      end if
      print *, ""

      print *, "11) Tracer: "
      if(include_tracer) then
        write(*, fmt = "(a25,a)") "tracer = ", "on"

        if(raytracer) then
          if(include_trfrc_lo) then
            write(*, fmt = "(a25,a)") "include_trfrc_lo = ", "on"
          else
            write(*, fmt = "(a25,a)") "include_trfrc_lo = ", "off"
          end if

          if(include_trfrc_mix) then
            write(*, fmt = "(a25,a)") "include_trfrc_mix = ", "on"
          else
            write(*, fmt = "(a25,a)") "include_trfrc_mix = ", "off"
          end if
        end if
      else
        write(*, fmt = "(a25,a)") "include_tracer = ", "off"
      end if

      print *, ""

    end if ! modified by Junhong Wei (20170216)

    !-------------------------------------------------------------------

    contains

    function discDensity(r)
      ! in/out
      real :: discDensity
      real, intent(in) :: r

      ! local variables
      real :: dens

      if(r >= 1.0) then ! outer region
        dens = 0.0
      else if(r < 0.5) then ! 0 <= r < r0/2  region
        dens = 4.0 * r ** 2
      else ! 1/2 <= r < 1
        dens = 4.0 * (1.0 - r) ** 2
      end if

      discDensity = dens

    end function discDensity

    !--------------------------------------------------------------------

    function greshoVelocity(u0, r) ! gresho vortex with normed radius r0=1
      ! in/out
      real :: greshoVelocity
      real, intent(in) :: u0, r

      ! local variables
      real :: u

      if(r >= 1.0) then ! outer region
        u = 0.0
      else if(r < 0.5) then ! 0 <= r < r0/2  region
        u = 2.0 * r
      else ! 1/2 <= r < 1
        u = 2.0 * (1.0 - r)
      end if

      greshoVelocity = u * u0

    end function greshoVelocity

    !--------------------------------------------------------------------

    function greshoPressure(pInf, u0, r)
      ! in/out
      real :: greshoPressure
      real, intent(in) :: pInf, u0, r

      ! local variables
      real :: p

      if(r >= 1.0) then ! outer region
        p = 0.0
      else if(r < 0.5) then ! 0 <= r < r0/2  region
        p = 4.0 * r ** 2
      else ! 1/2 <= r < 1
        p = 4.0 * (1.0 - r) ** 2
      end if

      greshoPressure = pInf - 0.5 * u0 ** 2 * p

    end function greshoPressure

    !---------------------------------------------------------------------

    subroutine init_GWP(Psi, kk, mm, ll, indwindcoeff)

      !------------------------------------------------
      !  calculate complex amplitudes for
      !    1) first harmonics,  leading order: Psi(:,:,:,0)
      !    2) second harmonics, leading order: Psi(:,:,:,1)
      !------------------------------------------------

      ! WARNING:
      ! 2nd harmonics probably not ready for 2D and 3D wave packets

      ! in/out variables
      ! wave amplitude
      complex, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, 5, 0:2), intent(out) :: &
          &Psi
      real, intent(out) :: kk, mm
      real, intent(out) :: ll
      real, intent(out) :: indwindcoeff

      ! local variables
      real :: A, rho
      real :: omi
      real :: kk2, mm2, omi2, kTot2
      complex :: u10, w10, b11, pi12
      integer :: iRay
      integer :: i, j, k
      integer :: ii, jj

      ! local amplitude values and derivatives
      complex :: u10_r, u10_c, u10_l, u10_t, u10_b
      complex :: w10_r, w10_c, w10_l, w10_t, w10_b
      complex :: theta11_r, theta11_c, theta11_l, theta11_t, theta11_b
      complex :: pi12_c
      complex :: pi0_t, pi0_c, pi0_b
      complex :: theta0_c
      complex :: du10_dx, du10_dz, dw10_dx, dw10_dz
      complex :: dpi0_dz, dtheta11_dx, dtheta11_dz
      complex :: Div
      complex :: Press, PressU, PressW
      complex :: d1u10, d1w10, d1theta11
      complex :: coeff, aux1
      complex :: M11, M12, M13, M14
      complex :: M21, M22, M23, M24
      complex :: M31, M32, M33, M34
      complex :: M41, M42, M43, M44
      complex, dimension(4) :: RHS, Sol
      complex, dimension(4, 4) :: M2, M2inv
      complex :: u21, w21, b22, pi23

      ! mean value calculation
      real :: rho0, rho0_t, rho0_b, d_dz, ypsi

      ! debugging stuff
      complex :: summe
      complex, dimension(4) :: term

      ! more debugging stuff
      real :: B11_pinc
      real :: D1TH11

      integer :: i0, j0 ! modified by Junhong Wei (20161201)

      complex :: tmp_var_3DWP

      real :: Ro_GWP, RoInv_GWP !FS

      if(f_Coriolis_dim /= 0.0) then !FS
        Ro_GWP = uRef / f_Coriolis_dim / lRef
        RoInv_GWP = 1.0 / Ro_GWP
      else
        Ro_GWP = 1.d40
        RoInv_GWP = 0.0
      end if

      !-----------------------
      !      MPI stuff
      !-----------------------
      i0 = is + nbx - 1 ! 0 index, replace i -> i + i0 in x and y fields
      j0 = js + nby - 1

      !------------------------
      !    Init data
      !-----------------------

      ! scale input data
      lambdaX = lambdaX_dim / lRef ! non-dim wave length in x dir.
      lambdaY = lambdaY_dim / lRef ! non-dim wave length in y dir.
      lambdaZ = lambdaZ_dim / lRef ! non-dim vert. wave length

      xCenter = x0_dim / lRef ! scaled position wave packtet x dir.
      yCenter = y0_dim / lRef ! scaled position wave packtet y dir.
      zCenter = z0_dim / lRef ! scaled position wave packtet z dir.

      sigma_x = sigma_hor_dim / lRef ! x width of Gaussian distribution
      sigma_y = sigma_hor_yyy_dim / lRef ! y width of Gaussian distribution
      sigma_z = sigma_dim / lRef ! vert. width of Gaussian distribution

      L_cos = L_cos_dim / lRef ! half length of cosine profile

      if(ABS(lambdaY_dim) /= 0.0) then
        ll = 2.0 * pi / lambdaY
      else
        ll = 0.0
      end if

      if(ABS(lambdaX_dim) /= 0.0) then
        kk = 2.0 * pi / lambdaX
      else
        kk = 0.0
      end if

      mm = 2.0 * pi / lambdaZ
      kk2 = kk ** 2
      mm2 = mm ** 2
      kTot2 = kk2 + mm2 + ll * ll
      kTot = sqrt(kTot2)

      ! intrinsic frequency
      omi = omiSign * sqrt(N2 * (kk * kk + ll * ll) + RoInv_GWP * RoInv_GWP &
          &* mm * mm) / kTot
      omi2 = omi ** 2

      ! amplitude coefficients for wave 1
      bAmp = amplitudeFactor * N2 / mm ! buoyancy
      uAmp = mm / kk * omi / N2 * bAmp
      wAmp = omi / N2 * bAmp
      pAmp = kappa * Ma2 * mm / kk ** 2 * omi2 / N2 * bAmp ! Exner pressure

      indwindcoeff = kk * omi * kTot2 / (N2 ** 2. * (kk * kk + ll * ll))
      close(10)

      !----------------------
      !  output of init data
      !----------------------

      if(master) then
        print *, "omi = ", omi / tRef
        print *, "mm = ", mm / lRef

        print *, "RoInv = ", RoInv_GWP / tRef ! modified by Junhong Wei
        print *, ""
        print *, "  0) Test case: "
        write(*, fmt = "(a25,a35)") "Test case  = ", "wave packet (full model)"
        write(*, fmt = "(a25,f10.1,a)") "lambda_x = ", lambdaX_dim, " m"
        write(*, fmt = "(a25,f10.1,a)") "lambda_z = ", lambdaZ_dim, " m"
        write(*, fmt = "(a25,f10.1a7)") "c_x  = ", omi / kk * uRef, " m/s"
        write(*, fmt = "(a25,f10.1,a7)") "c_z  = ", omi / mm * uRef, " m/s"
        write(*, fmt = "(a25,f10.1,a7)") "cg_x  = ", - NN * mm ** 2 / kTot &
            &** 3 * uRef, " m/s"
        write(*, fmt = "(a25,f10.1,a7)") "cg_z  = ", NN * mm * kk / kTot ** 3 &
            &* uRef, " m/s"
        write(*, fmt = "(a25,f10.1,a7)") "cg_z2  = ", - (NN ** 2. &
            &- (f_Coriolis_dim / lRef) ** 2.) * mm * (kk ** 2. + ll ** 2.) &
            &/ kTot2 ** 2. / omi * uRef, " m/s"
        write(*, fmt = "(a25,f10.1,a7)") "u_jet  = ", u0_jet_dim, " m/s"
        print *, ""
      end if ! modified by Junhong Wei (20170216)

      !---------------------------------------
      !        calc amplitude Psi_1^0
      !     (first harmonic, leading order)
      !---------------------------------------

      do k = 0, (nz + 1)
        do j = 0, (ny + 1)
          do i = 0, (nx + 1)
            ! profile: 1D and 2D

            if(wavePacketDim == 1) then
              delx = 0.0
            else
              delx = (x(i + i0) - xCenter)
            end if

            if(wavePacketDim == 3) then
              dely = (y(j + j0) - yCenter)
            else
              dely = 0.0
            end if

            if(topography) then
              ! TFC FJ
              delz = heightTFC(i, j, k) - zCenter
            else
              delz = (z(k) - zCenter)
            end if

            select case(wavePacketType)

            case(1)

              ! Gaussian
              ! cosine profile horizontally so that fields are zero
              ! at the horizontal boundaries
              ! in case of zero sigma in x or y direction use infinity

              !!$              if(sigma_x == 0.0) then
              !!$                envel = 1.0
              !!$              else if(abs(delx) < sigma_x) then
              !!$                envel = 1.0 - amp_mod_x + amp_mod_x * cos(delx * pi / (sigma_x &
              !!$                    * 2.0))
              !!$              else
              !!$                envel = 1.0 - amp_mod_x
              !!$              end if
              !!$
              !!$              if(sigma_y == 0.0) then
              !!$                envel = 1.0 * envel
              !!$              else if(abs(dely) < sigma_y) then
              !!$                envel = (1.0 - amp_mod_y + amp_mod_y * cos(dely * pi &
              !!$                    / (sigma_y * 2.0))) * envel
              !!$              else
              !!$                envel = envel * (1.0 - amp_mod_y)
              !!$              end if
              !!$
              !!$              envel = envel * exp(- (delz ** 2) / 2. / sigma_z ** 2)

              if(sigma_x == 0.0) then
                envel = 1.0
              else if(abs(delx) < sigma_x) then
                envel = 0.5 * (1.0 + cos(delx * pi / sigma_x))
              else
                envel = 0.
              end if

              if(sigma_y == 0.0) then
                envel = 1.0 * envel
              else if(abs(dely) < sigma_y) then
                envel = 0.5 * (1.0 + cos(dely * pi / sigma_y)) * envel
              else
                envel = 0.
              end if

              if(compare_raytracer) then
                envel = envel * exp(- (delz / sigma_z) ** 2)
              else
                envel = envel * exp(- (delz / sigma_z) ** 2 / 2.)
              end if

            case(2)

              ! IKAug2023 begin
              ! to match the raytracer cosine case
              if(sigma_x == 0.0) then
                envel = 1.0
              else if(abs(delx) < sigma_x) then
                envel = 0.5 * (1.0 + cos(delx * pi / sigma_x))
              else
                envel = 0.0
              end if
              if(sigma_y == 0.0) then
                envel = 1.0 * envel
              else if(abs(dely) < sigma_y) then
                envel = envel * 0.5 * (1.0 + cos(dely * pi / sigma_y))
              else
                envel = 0.0
              end if
              ! IKAug2023 end

              if(compare_raytracer) then
                ! Cosine
                if(abs(delz) .le. sigma_z) then
                  envel = 0.5 * (1.0 + cos(pi * delz / sigma_z)) * envel
                else
                  envel = 0.0
                end if
              else
                ! Cosine
                if(abs(delz) .le. L_cos) then
                  envel = envel * 0.5 * (1.0 + cos(pi * delz / L_cos))
                else
                  envel = 0.0
                end if
              end if
            case(4)
              envel = 0.0
              envel = 1. / (exp((z(k) - zCenter - sigma_z) / (lambdaZ)) + 1) &
                  &+ 1. / (exp(- (z(k) - zCenter + sigma_z) / (lambdaZ)) + 1) &
                  &- 1.
            case default
              stop "init.f90: unknown wavePacketType. Stop."

            end select

            if(compare_raytracer) then
              !SD
              !fix for wavepacket structure:
              !in raytracer cos-wavepacket in waveaction density A
              !here cos-wavepacket in buoyancy b with A \sim b^2
              envel = sqrt(envel)
            end if

            b11 = cmplx(envel * bAmp, 0.0)

            if(topography) then
              ! TFC FJ
              theta0 = thetaStratTFC(i, j, k)
            else
              theta0 = thetaStrat(k)
            end if

            tmp_var_3DWP = cmplx(0.0, (omi * omi - N2) / (mm * N2 * (omi * omi &
                &- RoInv_GWP * RoInv_GWP)))

            u10 = tmp_var_3DWP * cmplx(kk * omi, ll * RoInv_GWP) * b11

            w10 = cmplx(0.0, omi / N2) * b11

            pi12 = cmplx(0.0, kappa * Ma2 * (omi * omi - N2) / N2 / mm &
                &/ theta0) * b11

            Psi(i, j, k, :, 1) = (/u10, w10, b11, pi12, (tmp_var_3DWP &
                &* cmplx(ll * omi, - kk * RoInv_GWP) * b11)/)

          end do
        end do
      end do

      !---------------------------------------
      !        calc amplitude Psi_2^1
      !     (second harmonic, first order)
      !---------------------------------------

      ! WARNING:
      ! this part would probably still tp be adjusted fpor 2D or 3D wave p.

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx

            !-------------------------
            !       Set up RHS
            !-------------------------

            ! zonal velocities right, center, left, top, bottom
            u10_r = Psi(i + 1, j, k, 1, 1)
            u10_c = Psi(i, j, k, 1, 1)
            u10_l = Psi(i - 1, j, k, 1, 1)
            u10_t = Psi(i, j, k + 1, 1, 1)
            u10_b = Psi(i, j, k - 1, 1, 1)

            ! vertical velocities top, center, bottom
            w10_r = Psi(i + 1, j, k, 2, 1)
            w10_l = Psi(i - 1, j, k, 2, 1)
            w10_t = Psi(i, j, k + 1, 2, 1)
            w10_c = Psi(i, j, k, 2, 1)
            w10_b = Psi(i, j, k - 1, 2, 1)

            ! Buoyancy and pot. temp.
            if(topography) then
              ! TFC FJ
              theta11_r = Fr2 * thetaStratTFC(i + 1, j, k) * Psi(i + 1, j, k, &
                  &3, 1)
              theta11_c = Fr2 * thetaStratTFC(i, j, k) * Psi(i, j, k, 3, 1)
              theta11_l = Fr2 * thetaStratTFC(i - 1, j, k) * Psi(i - 1, j, k, &
                  &3, 1)
              theta11_t = Fr2 * thetaStratTFC(i, j, k + 1) * Psi(i, j, k + 1, &
                  &3, 1)
              theta11_b = Fr2 * thetaStratTFC(i, j, k - 1) * Psi(i, j, k - 1, &
                  &3, 1)
            else
              theta11_r = Fr2 * thetaStrat(k) * Psi(i + 1, j, k, 3, 1)
              theta11_c = Fr2 * thetaStrat(k) * Psi(i, j, k, 3, 1)
              theta11_l = Fr2 * thetaStrat(k) * Psi(i - 1, j, k, 3, 1)
              theta11_t = Fr2 * thetaStrat(k + 1) * Psi(i, j, k + 1, 3, 1)
              theta11_b = Fr2 * thetaStrat(k - 1) * Psi(i, j, k - 1, 3, 1)
            end if

            ! Second order Exner pressure
            pi12_c = Psi(i, j, k, 4, 1)

            ! Background Exner pressure and pot. temp.
            if(topography) then
              ! TFC FJ
              pi0_t = (pStratTFC(i, j, k + 1) / p0) ** gamma_1
              pi0_c = (pStratTFC(i, j, k) / p0) ** gamma_1
              pi0_b = (pStratTFC(i, j, k - 1) / p0) ** gamma_1
              theta0_c = thetaStratTFC(i, j, k)
            else
              pi0_t = (Pstrat(k + 1) / p0) ** gamma_1
              pi0_c = (Pstrat(k) / p0) ** gamma_1
              pi0_b = (Pstrat(k - 1) / p0) ** gamma_1
              theta0_c = thetaStrat(k)
            end if

            ! derivatives
            if(topography) then
              ! TFC FJ
              du10_dx = 0.5 * (jac(i + 1, j, k) * u10_r - jac(i - 1, j, k) &
                  &* u10_l) / dx / jac(i, j, k) + 0.5 * (jac(i, j, k + 1) &
                  &* met(i, j, k + 1, 1, 3) * u10_t - jac(i, j, k - 1) &
                  &* met(i, j, k - 1, 1, 3) * u10_b) / dz / jac(i, j, k)
              du10_dz = 0.5 * (u10_t - u10_b) / dz / jac(i, j, k)
              dw10_dx = 0.5 * (jac(i + 1, j, k) * w10_r - jac(i - 1, j, k) &
                  &* w10_l) / dx / jac(i, j, k) + 0.5 * (jac(i, j, k + 1) &
                  &* met(i, j, k + 1, 1, 3) * w10_t - jac(i, j, k - 1) &
                  &* met(i, j, k - 1, 1, 3) * w10_b) / dz / jac(i, j, k)
              dw10_dz = 0.5 * (w10_t - w10_b) / dz / jac(i, j, k)
              dpi0_dz = 0.5 * (pi0_t - pi0_b) / dz / jac(i, j, k)
              dtheta11_dx = 0.5 * (jac(i + 1, j, k) * theta11_r - jac(i - 1, &
                  &j, k) * theta11_l) / dx / jac(i, j, k) + 0.5 * (jac(i, j, k &
                  &+ 1) * met(i, j, k + 1, 1, 3) * theta11_t - jac(i, j, k &
                  &- 1) * met(i, j, k - 1, 1, 3) * theta11_b) / dz / jac(i, j, &
                  &k)
              dtheta11_dz = 0.5 * (theta11_t - theta11_b) / dz / jac(i, j, k)
            else
              du10_dx = 0.5 * (u10_r - u10_l) / dx
              du10_dz = 0.5 * (u10_t - u10_b) / dz
              dw10_dx = 0.5 * (w10_r - w10_l) / dx
              dw10_dz = 0.5 * (w10_t - w10_b) / dz
              dpi0_dz = 0.5 * (pi0_t - pi0_b) / dz
              dtheta11_dx = 0.5 * (theta11_r - theta11_l) / dx
              dtheta11_dz = 0.5 * (theta11_t - theta11_b) / dz
            end if

            ! divergence term -> Eq. (7.21)
            Div = - du10_dx - dw10_dz - (1. - kappa) / kappa * w10_c / pi0_c &
                &* dpi0_dz

            ! Pressure terms
            Press = 0.5 * kappaInv * MaInv2 * imag * theta11_c * pi12_c
            PressU = kk * Press
            PressW = mm * Press

            ! intermediate terms
            d1u10 = 0.5 * (u10_c * du10_dx + w10_c * du10_dz + Div * u10_c)
            d1w10 = 0.5 * (u10_c * dw10_dx + w10_c * dw10_dz + Div * w10_c)
            d1theta11 = 0.5 * (u10_c * dtheta11_dx + w10_c * dtheta11_dz + Div &
                &* theta11_c)

            RHS(1) = - d1u10 - PressU
            RHS(2) = - d1w10 - PressW
            RHS(3) = - FrInv2 / NN / theta0_c * d1theta11
            RHS(4) = (0.0, 0.0)

            !----------------------------------------------------
            !       Set up inverted system matrix M(2om,2k,2m)
            !----------------------------------------------------

            coeff = 1. / (4. * omi2 * kTot2 - kk2 * N2)
            aux1 = 4. * omi2 - N2

            M11 = 2. * imag * mm2 * omi; M12 = - 2. * imag * kk * mm * omi; &
                &M13 = kk * mm * NN; M14 = - 0.5 * imag * kk * aux1
            M21 = M12; M22 = 2. * imag * kk2 * omi; M23 = - kk2 * NN; M24 = &
                &- 2. * imag * omi2 * mm
            M31 = - M13; M32 = - M23; M33 = 2. * imag * omi * kTot2; M34 = &
                &- omi * mm * NN
            M41 = M14; M42 = M24; M43 = - M34; M44 = - 0.5 * imag * omi * aux1

            ! inverted matrix: check ok
            M2inv(1, :) = (/M11, M12, M13, M14/)
            M2inv(2, :) = (/M21, M22, M23, M24/)
            M2inv(3, :) = (/M31, M32, M33, M34/)
            M2inv(4, :) = (/M41, M42, M43, M44/)
            M2inv = coeff * M2inv

            !---------------------------------------
            !   Solve linear System -> save in Psi
            !---------------------------------------

            sol = matmul(M2inv, RHS)

            u21 = sol(1)
            w21 = sol(2)
            b22 = sol(3) * NN
            pi23 = sol(4) * kappa * Ma2 / theta0_c

            Psi(i, j, k, :, 2) = (/u21, w21, b22, pi23, (cmplx(0.0, 0.0) * &
                &b11)/)
          end do
        end do
      end do

    end subroutine init_GWP

    !-------------------------------------------------------------------

  end subroutine initialise

  !-------------------------------------------------------------------

  function cphase(c) result(phi) ! currently not being used

    !------------------------------------
    !  phase of complex number
    !------------------------------------

    ! in/out
    complex, intent(in) :: c
    real :: phi

    ! local vars
    real :: a, b

    a = real(c)
    b = aimag(c)

    ! first quadrant
    if(a >= 0 .and. b >= 0) then
      phi = atan(b / a)

      ! second quadrant
    else if(a < 0. .and. b >= 0.) then
      phi = pi - atan(- b / a)

      ! third quadrant
    else if(a < 0. .and. b < 0.) then
      phi = - pi + atan(b / a)

      ! fourth quadrant
    else if(a >= 0. .and. b < 0) then
      phi = - atan(- b / a)
    else
      stop "wkb.f90/cphase: case not included. Stop."
    end if

    phi = phi * 180. / pi

  end function cphase

  ! ----------------------------------------------------------------------
  subroutine noise_array(amplitude, ntovar, noise)

    ! local variables
    integer :: i, j, k
    integer :: allocstat, root
    real :: valRef, amplitude, perturbVal, sum_loc, sum_glob
    character(len = 20) :: ntovar
    real :: rand ! random number
    real :: Lx, Ly, Lz
    real, dimension(1:nx, 1:ny, 1:nz) :: noise

    !-------------------------------
    !  add noise
    !-------------------------------

    select case(ntovar)
    case("rho")
      valRef = rhoref
    case("vel")
      valRef = uref
    case("th")
      valRef = thetaRef
    case("none")
      valRef = 1.0
    case default
      stop "noise_array: Unknown variable"
    end select

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! calculate a random number between 0 and 1
          call Random_Number(rand)

          ! calc. perturbation scaled by amplitude
          perturbVal = 2.0 * (0.5 - rand) * proc_noise * amplitude

          ! non-dimensionalization ...

          if(referenceQuantities == "SI") then
            noise(i, j, k) = perturbVal
          else
            noise(i, j, k) = perturbVal / valRef
          end if
        enddo
      enddo
    enddo

  end subroutine noise_array

  !-------------------------------------------------------------------------

  subroutine output_background(th_bgr, max_nz, filename_bgr, ref)
    !-------------------------------
    !  writes background theta
    !-------------------------------

    integer max_nz
    ! argument fields
    real, dimension(1:nx, 1:ny, 1:max_nz), intent(in) :: th_bgr
    character(len = *) :: filename_bgr

    ! local and global output field and record nbs.
    real * 4, dimension(nx, ny) :: field_prc
    !   real*4,dimension(SizeX,SizeY) :: field_out
    integer irc_prc, irc_out

    ! local variables
    integer :: i, j, k

    integer :: i0, i1, j0, j1, k0, k1

    ! hotBubble output
    real :: rho, theta_dim, ref

    integer :: i_prc, i_mst, i_out, j_prc, j_mst, j_out

    if(master) then
      if((TestCase == "baroclinic_LC") .or. (TestCase == "baroclinic_ID")) then
        print *, ""
        print *, " Output of background into file ", filename_bgr
        print *, ""
      else
        print *, ""
        print *, " There is no environmental state for this testcase, EXIT "
        print *, ""
        stop
      end if
    end if

    !------------------------------
    !   prepare output file
    !------------------------------

    ! open output file

    !   open(40,file=dataFile,form="unformatted",access='direct',recl=nx*ny)
    if(master) then
      open(51, file = filename_bgr, form = "unformatted", access = 'direct', &
          &recl = SizeX * SizeY * sizeofreal4)
    end if

    !---------------------------------------
    !       dimensionalising and layerwise output
    !---------------------------------------

    irc_prc = 0
    do k = 1, max_nz
      ! dimensionalization
      do j = 1, ny
        do i = 1, nx
          select case(model)

          case("pseudo_incompressible", "compressible")

            if(referenceQuantities == "SI") then
              theta_dim = th_bgr(i, j, k)
            else
              theta_dim = th_bgr(i, j, k) * ref
            end if

            field_prc(i, j) = real(theta_dim, kind = 4)

          case("Boussinesq")
            stop "output_background: background undefined"

          case("WKB")
            stop "output_background: background undefined"
          case default
            stop "output_background: unknown model"
          end select ! model

        end do ! i
        call mpi_gather(field_prc(1, j), nx, mpi_real, field_mst(1, j), nx, &
            &mpi_real, 0, comm, ierror)
      end do ! j

      ! layerwise output
      irc_prc = irc_prc + 1
      call mpi_barrier(comm, ierror)
      if(master) then
        do j = 1, ny
          j_mst = j

          do j_prc = 1, nprocy
            j_out = ny * (j_prc - 1) + j

            do i_prc = 1, nprocx
              do i = 1, nx
                i_out = nx * (i_prc - 1) + i

                i_mst = nprocy * nx * (i_prc - 1) + (j_prc - 1) * nx + i

                field_out(i_out, j_out) = field_mst(i_mst, j_mst)
              end do
            end do
          end do
        end do

        write(51, rec = irc_prc) field_out
      end if
    end do ! k

    !------------------------------------
    !              close file
    !------------------------------------
    if(master) close(unit = 51)

  end subroutine output_background

end module init_module
