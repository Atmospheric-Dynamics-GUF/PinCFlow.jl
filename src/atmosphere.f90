module atmosphere_module

  use type_module
  use sizeof_module
  use mpi

  implicit none

  public ! all objects are known to programmes using this module

  !--------------
  !    public
  !--------------
  public :: init_atmosphere
  public :: terminate_atmosphere

  ! TFC
  public :: setHalosOfField2D
  public :: jac, met, chris, heightTFC, vertWindTFC, trafoTFC, stressTensTFC

  public :: update_topography

  ! compressible
  public :: add_JP_to_u

  real, dimension(:), allocatable :: PStrat, rhoStrat, thetaStrat, bvsStrat
  real, dimension(:), allocatable :: PStrat_0, rhoStrat_0
  real, dimension(:), allocatable :: rhoStrat_d, rhoStrat_s

  real, dimension(:), allocatable :: pistrat
  real, dimension(:), allocatable :: PStratTilde, rhoStratTilde, thetaStratTilde

  real, dimension(:), allocatable :: PStrat00, PStrat01, rhoStrat00, &
      &rhoStrat01, thetaStrat00, thetaStrat01, bvsStrat00, bvsStrat01, &
      &PStratTilde00, PStratTilde01, rhoStratTilde00, rhoStratTilde01, &
      &thetaStratTilde00, thetaStratTilde01

  real, dimension(:), allocatable :: Ro, RoInv

  ! reference quantites
  real :: rhoRef ! reference density
  real :: pRef ! reference pressure
  real :: uRef, aRef ! reference velocity / speed of sound
  real :: lRef ! reference length / pressure scale height h_sc
  real :: tRef ! ref time (time acoustic signal needs to pass domain)
  real :: thetaRef ! ref (potential) temperature   (Tref  = thetaRef)
  real :: Fref ! reference force

  ! natural constants
  real :: gamma ! adiabatic index
  real :: gammaInv ! 1/gamma
  real :: gamma_1 ! gamma-1
  real :: kappa ! (gamma-1)/gamma
  real :: kappaInv ! 1/kappa

  real, parameter :: g = 9.81 ! gravitational constant in m/s^2
  real, parameter :: Rsp = 287.0 ! spec. gas const. for dry air in J/kg/K
  real :: g_ndim ! nondimensional gravitational constant

  ! flow parameters
  real :: Re ! Reynolds number (calc. from input 1/Re)
  real :: Ma, MaInv2, Ma2 ! Mach number and 1/Ma^2, Ma^2
  real :: Fr, FrInv2, Fr2 ! Froude number Fr and 1/Fr^2, Fr^2
  real :: sig ! Ma^2/Fr^2

  ! pressure scale height
  real :: hp

  ! pot. temperature scale height
  real :: hTheta

  ! isentropic atmosphere
  real :: theta0
  real :: term

  ! stable atmosphere
  real :: N2 ! scaled square of Brunt-Vaisala frequency
  real :: NN ! scaled of Brunt-Vaisala frequency
  real :: coeff ! long coefficient

  ! Held-Suarez atmosphere
  real :: tp_strato ! stratosphere temperature
  real :: tp_srf_trp ! tropical surface temperature
  real :: tpdiffhor_tropo ! tropospheric temperature difference
  ! between poles and tropics
  real :: ptdiffvert_tropo ! vertical potential-temperature
  ! difference in troposphere

  ! isothermal
  real :: T0 ! scaled background temperature

  ! general
  real :: p0 ! scaled reference pressure at z=0

  real :: zk ! zk = z(k)

  real :: mountainHeight, mountainWidth, k_mountain
  real :: x_center, y_center

  ! 3D background fields.
  real, dimension(:, :, :), allocatable :: pStratTFC, thetaStratTFC, &
      &rhoStratTFC, bvsStratTFC, piStratTFC

  ! tracer variables
  real :: alphaTracer
  real, dimension(:, :, :), allocatable :: initialtracer, initialtracerrho

  contains

  subroutine init_atmosphere
    !----------------------------------------------------------------------
    ! allocate and initialize: coordinate system and background atmosphere
    ! set reference and thermodynamic quantities
    !----------------------------------------------------------------------

    ! local variables
    integer :: allocstat
    integer :: i, j, k
    real :: pBar, thetaBar
    real :: p1, p2, power ! exponents
    real :: eps
    real :: zmax
    real :: zk_half

    ! realistic atmosphere: isentropic troposphere / isothermal stratosphere
    real :: z_tr ! scaled height of troposphere
    real :: theta_tr ! scaled const. pot. temp. of troposphere
    real :: press_tr ! pressure at tropopause
    real :: P_tr ! pressure variable at tropopause
    real :: T_tr ! temperature at tropopause
    real :: delZ ! distance to tropopause

    ! for baroclinic case
    real :: T_c_b1, pow_t, pow_s, p_t_b ! tropopause quantities
    real :: T_bar, T_c_b, p_bar ! calculated quantities
    real :: tp_sponge, tp_sp1, tp_sp2

    real :: pistar, thetastar

    ! debugging
    integer, parameter :: errorlevel = 10 ! 0 -> no output

    integer :: i00, j00
    real :: yloc, ymax
    real, dimension(0:ny + 1) :: f_Coriolis_y

    real :: topos_u, topos_v, topos_w, x_sp_u, y_sp_u, z_sp_u, x_sp_v, y_sp_v, &
        &z_sp_v, x_sp_w, y_sp_w, z_sp_w, dRP_u, dIP_u, dRP_v, dIP_v, dRP_w, &
        &dIP_w

    real :: rhodl, pcoeff_r1, p_dim1

    ! allocate PStrat
    if(.not. allocated(pStrat)) then
      allocate(Pstrat(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate pStrat"
    end if

    ! allocate pStrat_0
    if(.not. allocated(pStrat_0)) then
      allocate(pStrat_0(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) then
        stop "atmosphere.f90: could not allocate pStrat_0"
      end if
    end if

    ! allocate pistrat
    if(.not. allocated(pistrat)) then
      allocate(pistrat(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate pistrat"
    end if

    ! allocate pStratTilde -> P at half levels
    if(.not. allocated(pStratTilde)) then
      allocate(PstratTilde(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not all. pStratTilde"
    end if

    ! allocate rhoStrat
    if(.not. allocated(rhoStrat)) then
      allocate(rhoStrat(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) then
        stop "atmosphere.f90: could not allocate rhoStrat"
      end if
    end if

    ! allocate rhoStrat_d
    if(.not. allocated(rhoStrat_d)) then
      allocate(rhoStrat_d(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) then
        stop "atmosphere.f90: could not allocate rhoStrat_d"
      end if
    end if

    if(.not. allocated(rhoStrat_s)) then
      allocate(rhoStrat_s(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) then
        stop "atmosphere.f90: could not allocate rhoStrat_s"
      end if
    end if

    ! allocate rhoStrat_0
    if(.not. allocated(rhoStrat_0)) then
      allocate(rhoStrat_0(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) then
        stop "atmosphere.f90: could not allocate rhoStrat_0"
      end if
    end if

    ! allocate rhoStratTilde
    if(.not. allocated(rhoStratTilde)) then
      allocate(rhoStratTilde(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could n. all. rhoStratTilde"
    end if

    ! allocate thetaStrat
    if(.not. allocated(thetaStrat)) then
      allocate(thetaStrat(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not all. thetaStrat"
    end if

    ! allocate squared Brunt-Vaisala frequency bvsStrat
    if(.not. allocated(bvsStrat)) then
      allocate(bvsStrat(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate bvsStrat"
    end if

    ! allocate thetaStratTilde
    if(.not. allocated(thetaStratTilde)) then
      allocate(thetaStratTilde(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: c. n. all. thetaStratTilde"
    end if

    ! allocate PStrat
    if(.not. allocated(pStrat00)) then
      allocate(Pstrat00(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate pStrat"
    end if

    ! allocate PStrat
    if(.not. allocated(pStrat01)) then
      allocate(Pstrat01(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate pStrat"
    end if

    ! allocate pStratTilde -> P at half levels
    if(.not. allocated(pStratTilde00)) then
      allocate(PstratTilde00(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not all. pStratTilde"
    end if

    ! allocate pStratTilde -> P at half levels
    if(.not. allocated(pStratTilde01)) then
      allocate(PstratTilde01(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not all. pStratTilde"
    end if

    ! allocate rhoStrat
    if(.not. allocated(rhoStrat00)) then
      allocate(rhoStrat00(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate rhoStrat"
    end if

    ! allocate rhoStratTilde
    if(.not. allocated(rhoStratTilde00)) then
      allocate(rhoStratTilde00(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could n. all. rhoStratTilde"
    end if

    ! allocate rhoStrat
    if(.not. allocated(rhoStrat01)) then
      allocate(rhoStrat01(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate rhoStrat"
    end if

    ! allocate rhoStratTilde
    if(.not. allocated(rhoStratTilde01)) then
      allocate(rhoStratTilde01(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could n. all. rhoStratTilde"
    end if

    if(.not. allocated(thetaStrat00)) then
      allocate(thetaStrat00(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not all. thetaStrat00"
    end if

    if(.not. allocated(thetaStratTilde00)) then
      allocate(thetaStratTilde00(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: c. n. all. thetaStratTilde"
    end if

    if(.not. allocated(thetaStrat01)) then
      allocate(thetaStrat01(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not all. thetaStrat"
    end if

    if(.not. allocated(thetaStratTilde01)) then
      allocate(thetaStratTilde01(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: c. n. all. thetaStratTilde"
    end if

    ! allocate squared Brunt-Vaisala frequency bvsStrat
    if(.not. allocated(bvsStrat00)) then
      allocate(bvsStrat00(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate bvsStrat"
    end if

    ! allocate squared Brunt-Vaisala frequency bvsStrat
    if(.not. allocated(bvsStrat01)) then
      allocate(bvsStrat01(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate bvsStrat"
    end if

    ! allocate Ro
    if(.not. allocated(Ro)) then
      allocate(Ro(0:ny + 1), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate Ro"
    end if

    ! allocate RoInv
    if(.not. allocated(RoInv)) then
      allocate(RoInv(0:ny + 1), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate RoInv"
    end if

    ! TFC
    ! Allocate 3D background fields.
    if(topography) then
      if(.not. allocated(pStratTFC)) then
        allocate(pStratTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
            &+ nbz)), stat = allocstat)
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate pStratTFC"
      end if

      if(.not. allocated(thetaStratTFC)) then
        allocate(thetaStratTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- &
            &1):(nz + 2)), stat = allocstat)
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate &
            &thetaStratTFC"
      end if

      if(.not. allocated(rhoStratTFC)) then
        allocate(rhoStratTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- 1):(nz &
            &+ 2)), stat = allocstat)
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate rhoStratTFC"
      end if

      if(.not. allocated(bvsStratTFC)) then
        allocate(bvsStratTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- 1):(nz &
            &+ 2)), stat = allocstat)
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate bvsStratTFC"
      end if

      if(.not. allocated(piStratTFC)) then
        allocate(piStratTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- 1):(nz &
            &+ 2)), stat = allocstat)
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate piStratTFC"
      end if
    end if

    ! save the initial tracer distributions in
    ! initialtracer and initialtracerrho = initialtracer*rho
    if(include_tracer) then
      if(.not. allocated(initialtracer)) then
        allocate(initialtracer(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate &
            &initialtracer"
      end if
      if(.not. allocated(initialtracerrho)) then
        allocate(initialtracerrho(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz &
            &+ nbz))
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate &
            &initialtracerrho"
      end if
    end if
    !----------------------------------
    !       auxiliary quantities
    !----------------------------------

    if(testCase == "agnesiMountain") then
      gamma = 1.0696864111498257 ! cf. Klein 2008
    else
      gamma = 1.4
    end if
    gamma_1 = gamma - 1.0
    kappa = (gamma - 1.0) / gamma ! = R/c_p = Gamma at Klein 2008
    kappaInv = 1.0 / kappa
    gammaInv = 1.0 / gamma

    !----------------------------------
    !      reference quantities &
    !      nondimensional numbers
    !----------------------------------

    select case(referenceQuantities)

    case("general")
      ! free references
      rhoRef = 1.184 ! in kg/m^3
      pRef = 101325.0 ! in Pa = kg/m/s^2
      lRef = pRef / rhoRef / g ! in m
      Ma = 0.1

      ! dependent references
      aRef = sqrt(pRef / rhoRef) ! in m/s
      uRef = Ma * aRef ! in m/s
      tRef = lRef / uRef ! in s
      thetaRef = aRef ** 2 / Rsp ! in K
      Fr = uRef / sqrt(g * lRef) !

    case("Klein")
      rhoRef = 1.184 ! in kg/m^3
      pRef = 101325.0 ! in Pa = kg/m/s^2
      !end if
      aRef = sqrt(pRef / rhoRef) ! in m/s
      uRef = aRef ! - "" -
      lRef = pRef / rhoRef / g ! in m
      tRef = lRef / aRef ! in s
      thetaRef = aRef ** 2 / Rsp ! in K
      FRef = rhoRef * uRef ** 2 / lRef ! in N/m^3 = reference force

      ! nondim numbers
      Ma = uRef / aRef ! Ma = 1
      Fr = uRef / sqrt(g * lRef) ! Fr = 1
      !                          ! sig = Ma^2/Fr^2 = 1

    case("WKB")
      eps = 0.1 ! asymptotic parameter ( Ma = Fr for eps = kappa! )
      rhoRef = 1.184 ! in kg/m^3
      pRef = 101325.0 ! in Pa = kg/m/s^2

      aRef = sqrt(pRef / rhoRef) ! in m/s
      thetaRef = aRef ** 2 / Rsp ! in K
      Ma = eps / sqrt(kappa)
      Fr = sqrt(eps)
      uRef = Ma * aRef ! in m/2
      lRef = eps * pRef / (kappa * rhoRef * g) ! in m
      tRef = lRef / uRef

    case("SI")
      stop "init_atmosphere: Problems w. Exner pr.. Use Klein's scaling."
      ! with this scaling all quantities are
      ! in SI units
      ! Note that in this case the thermodynamic
      ! equations have a different form
      ! while the Euler equations are unchanged
      rhoRef = 1.0 ! in kg/m^3
      pRef = 1.0 ! in Pa = kg/m/s^2
      lRef = 1.0 ! in m
      tRef = 1.0 ! in s
      uRef = 1.0 ! in m/s
      aRef = 1.0 ! in m/s
      thetaRef = 1.0 ! in K

      ! nondim numbers
      Ma = 1.0 ! normally 1/sqrt(kappa)
      !                          ! but we use Exner-pr./kappa like R. Klein
      Fr = 1.0 / sqrt(g) ! Fr = uRef/sqrt(g*lRef)

    case default
      print *, "referenceQuantities = ", referenceQuantities
      stop "init_atmosphere: unknown referenceQuantities. Stopping."
    end select

    ! auxiliary nondimensionals
    Fr2 = Fr ** 2
    Ma2 = Ma ** 2
    MaInv2 = 1.0 / Ma ** 2
    FrInv2 = 1.0 / Fr ** 2
    sig = Ma ** 2 / Fr ** 2

    ! Reynolds number
    if(.not. specifyReynolds) ReInv = mu_viscous_dim / uRef / lRef

    if(ReInv < 1.0e-20) then
      Re = 1.0e20 ! only used in timestep calculation routine
    else
      Re = 1.0 / ReInv
    end if

    ! Heat conduction
    mu_conduct = mu_conduct_dim / uRef / lRef + 1.0e-20

    ! scaled reference pressure at z = 0
    p0 = press0_dim / pRef

    ! scaled background flow
    backgroundFlow = backgroundFlow_dim / uRef

    ! nondimensional gravitational constant
    g_ndim = g / (uRef ** 2 / lRef)

    !----------------------------------
    !            setup domain
    !----------------------------------

    ! scale the domain by reference length lRef
    lx = lx_dim / lRef
    ly = ly_dim / lRef
    lz = lz_dim / lRef

    ! init cell size

    dx = (lx(1) - lx(0)) / real(sizeX)
    dy = (ly(1) - ly(0)) / real(sizeY)
    dz = (lz(1) - lz(0)) / real(sizeZ)

    ! init cell coordinates
    do i = - nbx, sizeX + nbx
      x(i) = lx(0) + real(i - 1) * dx + dx / 2.0
    end do

    do j = - nby, sizeY + nby
      y(j) = ly(0) + real(j - 1) * dy + dy / 2.0
    end do

    do k = - nbz, sizeZ + nbz
      z(k) = lz(0) + real(k - 1) * dz + dz / 2.0
    end do

    j00 = js + nby - 1 !FS
    if(TestCase == "baroclinic_LC") then
      ymax = ly_dim(1) / lRef
      do j = 0, ny + 1
        yloc = y(j + j00)
        f_Coriolis_y(j) = f_Coriolis_dim
        if(f_Coriolis_y(j) /= 0.0) then
          Ro(j) = uRef / f_Coriolis_y(j) / lRef
          RoInv(j) = 1.0 / Ro(j)
        else
          Ro(j) = 1.d40
          RoInv(j) = 0.0
        end if
      end do
    else
      if(f_Coriolis_dim /= 0.0) then
        Ro(0:ny + 1) = uRef / f_Coriolis_dim / lRef
        RoInv(0:ny + 1) = 1.0 / Ro(:)
      else
        Ro(0:ny + 1) = 1.d40
        RoInv(0:ny + 1) = 0.0
      end if
    end if

    ! initial large-scale tracer distribution
    ! X(z) = alpha*z
    ! for linear_increase_in_z
    if(include_tracer) then
      alphaTracer = lRef
    end if

    !----------------------------------
    !            setup topography
    !----------------------------------

    call setup_topography

    !---------------------------------------------
    !   Set up Sponge layer
    !---------------------------------------------

    if(spongeLayer) then
      if(unifiedSponge) then
        if(spongeType == "exponential") then
          kSponge = 1
        else
          kSponge = nz - ceiling(spongeHeight * real(nz))
        end if

        dzSponge = spongeHeight * (lz(1) - lz(0))
        zSponge = lz(1) - dzSponge

        if(lateralSponge) then
          dxSponge = 0.5 * spongeHeight * (lx(1) - lx(0))
          dySponge = 0.5 * spongeHeight * (ly(1) - ly(0))
          xSponge0 = lx(0) + dxSponge
          ySponge0 = ly(0) + dySponge
          xSponge1 = lx(1) - dxSponge
          ySponge1 = ly(1) - dySponge
        end if
      else
        kSponge = nz - ceiling(spongeHeight * real(nz))
        zSponge = lz(0) + (1.0 - spongeHeight) * (lz(1) - lz(0))
      end if
    end if

    select case(model)
      !--------------------------------------------------------------------
      !             Atmospheres for pseudo-incompressible
      !--------------------------------------------------------------------

    case("pseudo_incompressible", "WKB", "compressible")

      select case(background)

      case('smoothV')

        ! some stratified fields are still in use...
        rhodl = 1. / 0.5
        pcoeff_r1 = 1024. ** 2 * (1. / 72. - 6. / 35. + 15. / 17. - 74. / 33. &
            &+ 57. / 32. + 174. / 31. - 269. / 15. + 450. / 29. + 153. / 8. &
            &- 1564. / 27. + 510. / 13. + 204. / 5. - 1. / 24. * (2210. &
            &- rhodl) + 12. / 23. * (85. - rhodl) + (510. / 11. + 3. * rhodl) &
            &- 4. / 21. * (391. + 55. * rhodl) + 9. / 40. * (119. + 110. &
            &* rhodl) + 18. / 19. * (25. - 44. * rhodl) - 1. / 9. * (269. &
            &- 462 * rhodl) + 6. / 17. * (29. - 132 * rhodl) + 3. / 16. * (19. &
            &+ 165. * rhodl) - 2. / 15. * (37. + 110. * rhodl) + 3. / 7. * (5. &
            &+ 11. * rhodl) - 6. / 13. * (1. + 2. * rhodl) + 1. / 24. * (1. &
            &+ 2. * rhodl))

        p_dim1 = pcoeff_r1 / pRef * 0.5 * (1. ** 2 + 1. ** 2)

        NN = 0.0
        N2 = NN ** 2

        pistrat = kappaInv * p_dim1 ** kappa
        rho00 = 0.5 / rhoRef
        P00 = (kappaInv * p_dim1 ** kappa) ** ((1.0 - kappa) / kappa)
        theta00 = P00 / rho00

        rhostrat = rho00
        bvsstrat = 0.
        PStrat = P00
        thetaStrat = theta00
        PStratTilde = P00
        thetaStratTilde = theta00
        rhoStratTilde = rho00

        ! Define 3D background fields.
        if(topography) then
          ! Define pStratTFC.
          pStratTFC = p00
          ! Define thetaStratTFC.
          thetaStratTFC = theta00
          ! Define rhoStratTFC.
          rhoStratTFC = rho00
        end if

        !-----------------------------------------------------------
        !   Isentropic troposphere / isothermal stratosphere
        !-----------------------------------------------------------

      case('realistic')

        ! not implemented versions
        if(referenceQuantities == "SI") stop "atmosphere.f90: &
            &referenceQuantities = SI not implemented"

        ! quantities at tropopause
        z_tr = z_tr_dim / lRef
        theta_tr = theta_tr_dim / thetaRef
        ! presure
        press_tr = p0 * (1.0 - kappa * sig / theta_tr * z_tr) ** (1 / kappa)
        ! temperature
        T_tr = theta_tr * (press_tr / p0) ** kappa

        ! quantities at full levels
        do k = - 1, nz + 2
          zk = z(k)
          delZ = zk - z_tr

          if(zk <= z_tr) then ! isentropic troposphere

            thetaStrat(k) = theta_tr

            power = 1.0 / (gamma - 1.0)
            term = kappa * sig / theta_tr * zk
            if(term > 1.0) stop "init_atmosphere: negative term. Stop."

            Pstrat(k) = p0 * (1.0 - term) ** power
            rhoStrat(k) = PStrat(k) / thetaStrat(k)

          else ! isothermal stratosphere

            thetaStrat(k) = theta_tr * exp(kappa * sig / T_tr * delZ)
            Pstrat(k) = p0 ** kappa * press_tr ** (1 / gamma) * exp(- sig &
                &/ gamma / T_tr * delZ)
            rhoStrat(k) = PStrat(k) / thetaStrat(k)

          end if

        end do ! k loop

        ! quantities at half levels
        do k = - 1, nz + 2
          if(k == nz + 2) then
            zk_half = z(k) + dz / 2.
            delZ = zk_half - z_tr

            if(zk_half <= z_tr) then ! isentropic troposphere
              thetaStratTilde(k) = theta_tr

              power = 1.0 / (gamma - 1.0)
              term = kappa * sig / theta_tr * zk_half
              if(term > 1.0) stop "init_atmosphere: negative term. Stop."

              PstratTilde(k) = p0 * (1.0 - term) ** power
              rhoStratTilde(k) = PStratTilde(k) / thetaStratTilde(k)
            else ! isothermal stratosphere
              thetaStratTilde(k) = theta_tr * exp(kappa * sig / T_tr * delZ)
              PstratTilde(k) = p0 ** kappa * press_tr ** (1 / gamma) * exp(- &
                  &sig / gamma / T_tr * delZ)
              rhoStratTilde(k) = PStratTilde(k) / thetaStratTilde(k)
            end if
          else
            thetaStratTilde(k) = 0.5 * (thetaStrat(k) + thetaStrat(k + 1))
            PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k + 1))
            rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k + 1))
          end if
        end do ! k loop

        ! Define 3D background fields.
        if(topography) then
          do i = - nbx, nx + nbx
            do j = - nby, ny + nby
              do k = - 1, nz + 2
                if(heightTFC(i, j, k) <= z_tr) then
                  ! Isentropic troposphere.
                  power = 1.0 / (gamma - 1.0)
                  term = kappa * sig / theta_tr * heightTFC(i, j, k)
                  if(term > 1.0) then
                    stop "init_atmosphere: negative term. Stop"
                  end if
                  ! Define pStratTFC.
                  pStratTFC(i, j, k) = p0 * (1.0 - term) ** power
                  ! Define thetaStratTFC.
                  thetaStratTFC(i, j, k) = theta_tr
                  ! Define rhoStratTFC.
                  rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                      &j, k)
                else
                  ! Isothermal stratosphere.
                  ! Define pStratTFC.
                  pStratTFC(i, j, k) = p0 ** kappa * press_tr ** (1.0 / gamma) &
                      &* exp(- sig / gamma / T_tr * (heightTFC(i, j, k) - z_tr))
                  ! Define thetaStratTFC.
                  thetaStratTFC(i, j, k) = theta_tr * exp(kappa * sig / T_tr &
                      &* (heightTFC(i, j, k) - z_tr))
                  ! Define rhoStratTFC.
                  rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                      &j, k)
                end if
              end do
            end do
          end do
        end if

        !-----------------------------------------------------------
        !                    Isothermal atmosphere
        !-----------------------------------------------------------

      case('isothermal')

        if(referenceQuantities == "SI") then
          !------------------------------------
          !   original equations in SI units
          !------------------------------------
          stop "referenceQuantities = SI not possible currently."
          T0 = Temp0_dim / thetaRef ! T0 in K
          N2 = kappa * g ** 2 / Rsp / T0 ! isothermal Brunt-Vaisala fr.^2
          NN = sqrt(N2) !

          hp = Rsp * T0 / g ! pressure scale height
          hTheta = hp / kappa ! pot. temperature scale height

          do k = - 1, nz + 2
            PStrat(k) = p0 * exp(- z(k) / gamma / hp)
            thetaStrat(k) = T0 * exp(z(k) / hTheta)
            rhoStrat(k) = 1.0 / Rsp * PStrat(k) / thetaStrat(k)
          end do

        else
          !-----------------------------------------
          !    with reference quantities
          !-----------------------------------------

          T0 = Temp0_dim / thetaRef
          N2 = Ma ** 2 / Fr ** 4 * kappa / T0 ! isothermal Brunt-Vaisala fr.^2
          NN = sqrt(N2)

          do k = - 1, nz + 2
            PStrat(k) = p0 * exp(- sig * z(k) / gamma / T0)
            thetaStrat(k) = T0 * exp(kappa * sig / T0 * z(k))
            rhoStrat(k) = PStrat(k) / thetaStrat(k)
          end do

          ! rhoStrat and PStrat at half levels
          do k = - 1, nz + 2
            if(k == nz + 2) then
              zk_half = z(k) + 0.5 * dz ! half level height
              PStratTilde(k) = p0 * exp(- sig * zk_half / gamma / T0)
              thetaStratTilde(k) = T0 * exp(kappa * sig / T0 * zk_half)
              rhoStratTilde(k) = PStratTilde(k) / thetaStratTilde(k)
            else
              thetaStratTilde(k) = 0.5 * (thetaStrat(k) + thetaStrat(k + 1))
              PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k + 1))
              rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k + 1))
            end if
          end do

          ! Define 3D background fields.
          if(topography) then
            do i = - nbx, nx + nbx
              do j = - nby, ny + nby
                do k = - 1, nz + 2
                  ! Define pStratTFC.
                  pStratTFC(i, j, k) = p0 * exp(- sig * heightTFC(i, j, k) &
                      &/ gamma / T0)
                  ! Define thetaStratTFC.
                  thetaStratTFC(i, j, k) = T0 * exp(kappa * sig / T0 &
                      &* heightTFC(i, j, k))
                  ! Define rhoStratTFC.
                  rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                      &j, k)
                end do
              end do
            end do
          end if

        end if

        ! GBcorr
        bvsStrat = N2

        !-------------------------------------------------------------------
        !                        Isentropic atmosphere
        !-------------------------------------------------------------------

      case('isentropic')

        if(referenceQuantities == "SI") then
          !------------------------------------
          !   original equations in SI units
          !------------------------------------

          stop "init_atmosphere: fluctuationMode not implmented for SI!"

          NN = 0.0
          N2 = 0.0
          theta0 = theta0_dim / thetaRef ! theta0 in K
          power = 1.0 / (gamma - 1.0)
          do k = - 1, nz + 2
            term = kappa * g / Rsp / theta0 * z(k)
            if(term > 1.0) then
              stop "init_atmosphere: negative term with power.Stop."
            end if
            PStrat(k) = p0 * (1.0 - term) ** power
            thetaStrat(k) = theta0
            rhoStrat(k) = 1.0 / Rsp * PStrat(k) / thetaStrat(k)
          end do

        else
          !-----------------------------------------
          !    with reference quantities
          !-----------------------------------------
          thetaStrat = theta0_dim / thetaRef
          theta0 = theta0_dim / thetaRef
          NN = 0.0
          N2 = 0.0
          ! hydrostatic background pressure

          do k = - 1, nz + 2
            power = 1.0 / (gamma - 1.0)
            term = kappa * sig / theta0 * z(k)
            if(term > 1.0) then
              print *, "init_atmosphere: neg. term with power.Stopping."
              print *, "term = ", term
              print *, "kappa = ", kappa
              print *, "sig = ", sig
              print *, "theta0 = ", theta0
              print *, "z(nz) = ", z(k)
              print *, "z(nz)*l = ", z(k)
              print *, "lRef = ", lRef
              stop "stopping."
            end if
            Pstrat(k) = p0 * (1.0 - term) ** power

            rhoStrat(k) = PStrat(k) / theta0

          end do

          ! quantities at half levels
          do k = - 1, nz + 2
            if(k == nz + 2) then
              zk_half = z(k) + dz / 2.
              power = 1.0 / (gamma - 1.0)
              term = kappa * sig / theta0 * zk_half
              if(term > 1.0) then
                print *, "init_atmosphere: neg. term with power."
                print *, "term = ", term
                print *, "kappa = ", kappa
                print *, "sig = ", sig
                print *, "theta0 = ", theta0
                print *, "z(nz) = ", z(k)
                print *, "z(nz)*l = ", z(k)
                print *, "lRef = ", lRef
                stop "stopping."
              end if
              PstratTilde(k) = p0 * (1.0 - term) ** power
              rhoStratTilde(k) = PStratTilde(k) / theta0
            else
              PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k + 1))
              rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k + 1))
            end if
          end do

          ! Define 3D background fields.
          if(topography) then
            do i = - nbx, nx + nbx
              do j = - nby, ny + nby
                do k = - 1, nz + 2
                  power = 1.0 / (gamma - 1.0)
                  term = kappa * sig / theta0 * heightTFC(i, j, k)
                  if(term > 1.0) then
                    stop "init_atmosphere: negative term. Stop"
                  end if
                  ! Define pStratTFC.
                  pStratTFC(i, j, k) = p0 * (1.0 - term) ** power
                  ! Define thetaStratTFC.
                  thetaStratTFC(i, j, k) = theta0
                  ! Define rhoStratTFC.
                  rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                      &j, k)
                end do
              end do
            end do
          end if

        end if

        ! GBcorr
        bvsStrat = N2

        !----------------------------------------------------------------
        !                   Stable atmosphere with constant N
        !----------------------------------------------------------------

      case('const-N')

        if(referenceQuantities == "SI") then
          !------------------------------------
          !   original equations in SI units
          !------------------------------------
          stop "init_atmosphere: fluctuationMode not implmented for SI!"

          theta0 = theta0_dim / thetaRef ! theta0 at z=0 in K
          NN = N_BruntVaisala_dim * tRef
          N2 = NN ** 2 ! Brunt-Vaisala N0^2 in 1/s^2

          coeff = kappa * g ** 2 / (Rsp * N2 * theta0)
          power = 1. / (gamma - 1.)

          do k = - 1, nz + 2
            term = exp(- N2 / g * z(k))

            !if( term > 1.0 ) stop "init_atmosphere: root of a negative number."
            if(1.0 + coeff * (term - 1.0) < 0.0) stop "init_atmosphere: root &
                &of a negative number."

            PStrat(k) = p0 * (1.0 + coeff * (term - 1.0)) ** power
            thetaStrat(k) = theta0 * exp(N2 / g * z(k))
            rhoStrat(k) = 1.0 / Rsp * PStrat(k) / thetaStrat(k)
          end do

        else
          !-----------------------------------------
          !    with reference quantities
          !-----------------------------------------

          theta0 = theta0_dim / thetaRef ! theta0 at z=0
          N2 = (N_BruntVaisala_dim * tRef) ** 2 ! Brunt-Vaisala fr.^2
          NN = sqrt(N2)

          power = 1.0 / (gamma - 1.0)

          ! potential temperature and pressure Pstrat
          do k = - 1, nz + 2
            thetaStrat(k) = theta0 * exp(Fr ** 2 * N2 * z(k))

            term = exp(- Fr2 * N2 * z(k))

            if(1.0 + FrInv2 * kappa * sig / N2 / theta0 * (term - 1.0) < 0.0) &
                &then
              print *, "init_atmosphere: power of a neg. number."
              stop 'top of atmosphere too high for const-N'
            end if

            Pstrat(k) = p0 * (1.0 + FrInv2 * kappa * sig / N2 / theta0 * (term &
                &- 1.0)) ** power

            rhoStrat(k) = pStrat(k) / thetaStrat(k)

          end do

          ! rhoStrat at half levels
          do k = - 1, nz + 2
            if(k == nz + 2) then
              zk_half = z(k) + dz / 2.
              thetaStratTilde(k) = theta0 * exp(Fr ** 2 * N2 * zk_half)

              term = exp(- Fr2 * N2 * zk_half)

              if(term > 1.0) then
                stop "init_atmosphere: root of a negative number."
              end if

              if(1.0 + FrInv2 * kappa * sig / N2 / theta0 * (term - 1.0) &
                  &< 0.0) then
                print *, "init_atmosphere: power of a neg. number."
                stop 'top of atmosphere too high for const-N'
              end if

              PstratTilde(k) = p0 * (1.0 + FrInv2 * kappa * sig / N2 / theta0 &
                  &* (term - 1.0)) ** power
              rhoStratTilde(k) = PstratTilde(k) / thetaStratTilde(k)
            else
              thetaStratTilde(k) = 0.5 * (thetaStrat(k) + thetaStrat(k + 1))
              PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k + 1))
              rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k + 1))
            end if
          end do

          ! Define 3D background fields.
          if(topography) then
            do i = - nbx, nx + nbx
              do j = - nby, ny + nby
                do k = - 1, nz + 2
                  power = 1.0 / (gamma - 1.0)
                  term = exp(- Fr2 * N2 * heightTFC(i, j, k))
                  ! Define pStratTFC.
                  pStratTFC(i, j, k) = p0 * (1.0 + FrInv2 * kappa * sig / N2 &
                      &/ theta0 * (term - 1.0)) ** power
                  ! Define thetaStratTFC.
                  thetaStratTFC(i, j, k) = theta0 * exp(Fr ** 2.0 * N2 &
                      &* heightTFC(i, j, k))
                  ! Define rhoStratTFC.
                  rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                      &j, k)
                end do
              end do
            end do
          end if

        end if

        ! GBcorr
        bvsStrat = N2

        !-----------------------------------------------------------
        ! Setting for a baroclinic life cycle. Different lapse rates
        ! in troposphere and stratosphere that are also y-dependent.
        ! Change: Elena Gagarina, 15.03.2018
        !-----------------------------------------------------------

      case('diflapse')

        if(referenceQuantities == "SI") then
          stop "atmosphere.f90: referenceQuantities = SI not impl."
        end if

        ! quantities at tropopause

        if(gamma_t /= 0.) pow_t = g / (Rsp * gamma_t)
        if(gamma_s /= 0.) pow_s = g / (Rsp * gamma_s)

        if(gamma_t /= 0.) then
          p_t_b = press0_dim * (1. - gamma_t * z_tr_dim / Temp0_dim) ** pow_t
        else
          p_t_b = press0_dim * exp(- z_tr_dim / (Rsp * Temp0_dim / g))
        end if

        T_c_b = Temp0_dim - gamma_t * z_tr_dim

        ! quantities at full levels:

        do k = - 1, nz + 2
          zk = z(k) * lRef

          if(zk < z_tr_dim) then
            T_bar = Temp0_dim - gamma_t * zk

            if(gamma_t /= 0.) then
              p_bar = press0_dim * (1. - gamma_t * zk / Temp0_dim) ** pow_t
            else
              p_bar = press0_dim * exp(- zk / (Rsp * Temp0_dim / g))
            end if
          else
            T_bar = Temp0_dim - gamma_t * z_tr_dim - gamma_s * (zk - z_tr_dim)

            if(gamma_s /= 0.) then
              p_bar = p_t_b * (1. - gamma_s * (zk - z_tr_dim) / T_c_b) ** pow_s
            else
              p_bar = p_t_b * exp(- (zk - z_tr_dim) / (Rsp * T_c_b / g))
            end if
          endif

          thetaStrat(k) = T_bar * (press0_dim / p_bar) ** kappa / thetaRef

          rhoStrat(k) = p_bar / (Rsp * T_bar) / rhoRef

          Pstrat(k) = rhoStrat(k) * thetaStrat(k)
        enddo

        ! quantities at half levels

        ! with the exception of the uppermost ghost layer the
        ! values there are obtained by linear interpolation, in order
        ! to be consistent with the handling of the half levels in the
        ! semi-implicit time stepping and the corresponding pressure
        ! solver

        do k = - 1, nz + 2
          if(k == nz + 2) then
            zk_half = (z(k) + dz / 2.) * lRef

            if(zk_half < z_tr_dim) then
              T_bar = Temp0_dim - gamma_t * zk_half

              if(gamma_t /= 0.) then
                p_bar = press0_dim * (1. - gamma_t * zk_half / Temp0_dim) &
                    &** pow_t
              else
                p_bar = press0_dim * exp(- zk_half / (Rsp * Temp0_dim / g))
              end if
            else
              T_bar = Temp0_dim - gamma_t * z_tr_dim - gamma_s * (zk_half &
                  &- z_tr_dim)

              if(gamma_s /= 0.) then
                p_bar = p_t_b * (1. - gamma_s * (zk_half - z_tr_dim) / T_c_b) &
                    &** pow_s
              else
                p_bar = p_t_b * exp(- (zk_half - z_tr_dim) / (Rsp * T_c_b / g))
              end if
            endif

            thetaStratTilde(k) = T_bar * (press0_dim / p_bar) ** (kappa) &
                &/ thetaRef

            rhoStratTilde(k) = p_bar / (Rsp * T_bar) / rhoRef

            PstratTilde(k) = rhoStratTilde(k) * thetaStratTilde(k)
          else
            PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k + 1))

            rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k + 1))

            thetaStratTilde(k) = PStratTilde(k) / rhoStratTilde(k)
          end if
        enddo

        ! Define 3D background fields.
        if(topography) then
          do i = - nbx, nx + nbx
            do j = - nby, ny + nby
              do k = - 1, nz + 2
                zk = heightTFC(i, j, k) * lRef

                if(zk < z_tr_dim) then
                  T_bar = Temp0_dim - gamma_t * zk
                  if(gamma_t /= 0.0) then
                    p_bar = press0_dim * (1.0 - gamma_t * zk / Temp0_dim) &
                        &** pow_t
                  else
                    p_bar = press0_dim * exp(- zk / (Rsp * Temp0_dim / g))
                  end if
                else
                  T_bar = Temp0_dim - gamma_t * z_tr_dim - gamma_s * (zk &
                      &- z_tr_dim)
                  if(gamma_s /= 0.0) then
                    p_bar = p_t_b * (1.0 - gamma_s * (zk - z_tr_dim) / T_c_b) &
                        &** pow_s
                  else
                    p_bar = p_t_b * exp(- (zk - z_tr_dim) / (Rsp * T_c_b / g))
                  end if
                end if

                ! Define thetaStratTFC.
                thetaStratTFC(i, j, k) = T_bar * (press0_dim / p_bar) ** kappa &
                    &/ thetaRef
                ! Define rhoStratTFC.
                rhoStratTFC(i, j, k) = p_bar / (Rsp * T_bar) / rhoRef
                ! Define pStratTFC.
                pStratTFC(i, j, k) = rhoStratTFC(i, j, k) * thetaStratTFC(i, &
                    &j, k)
              end do
            end do
          end do
        end if

        !-----------------------------------------------------------
        ! Setting for an atmmosphere according to Held & Suarez (1994)
        !-----------------------------------------------------------

      case('HeldSuarez')

        if(referenceQuantities /= "Klein") then
          stop "atmosphere.f90: referenceQuantities = Klein expected!"
        end if

        ! non-dimensional parameters of the Held-Suarez reference
        ! atmosphere:

        ! stratosphere temperature
        tp_strato = tp_strato_dim / thetaRef
        ! tropical surface temperature
        tp_srf_trp = tp_srf_trp_dim / thetaRef
        ! tropospheric temperature difference between poles and tropics
        tpdiffhor_tropo = tpdiffhor_tropo_dim / thetaRef
        ! vertical potential-temperature difference in troposphere
        ptdiffvert_tropo = ptdiffvert_tropo_dim / thetaRef

        ! Exner pressure just above and below the surface so that it
        ! = 1 at the surface, by an Euler integration of hydrostatic
        ! equilibrium

        T_bar = max(tp_strato, tp_srf_trp - 0.5 * tpdiffhor_tropo)

        pistrat(0) = 1.0 + 0.5 * dz * kappa / T_bar
        pistrat(1) = 1.0 - 0.5 * dz * kappa / T_bar

        if(pistrat(1) <= 0.) then
          print *, 'ERROR: non-positive pistrat at k =', 1
          stop
        end if

        ! potential temperature just below and above the surface
        ! P and density determined as well

        do k = 0, 1
          T_bar = max(tp_strato, pistrat(k) * (tp_srf_trp - 0.5 &
              &* tpdiffhor_tropo - 0.5 * ptdiffvert_tropo / kappa &
              &* log(pistrat(k)))) !FS 0.5->1

          thetaStrat(k) = T_bar / pistrat(k)

          pStrat(k) = pistrat(k) ** ((1.0 - kappa) / kappa)

          rhoStrat(k) = pStrat(k) / thetaStrat(k)
        end do

        ! for k > 1:
        ! Exner pressure and potential temperature by upward
        ! integration of hydrostatic equilibrium, using a trapezoidal
        ! leapfrog
        ! P and density determined as well

        do k = 2, nz + 2 !nz-ceiling((0.4)*real(nz))!nz+2 !FS
          pistar = pistrat(k - 2) - 2.0 * dz * kappa / thetaStrat(k - 1)

          if(pistar <= 0.) then
            print *, 'ERROR: non-positive pistar at k =', k
            stop
          end if

          T_bar = max(tp_strato, pistar * (tp_srf_trp - 0.5 * tpdiffhor_tropo &
              &- 0.5 * ptdiffvert_tropo / kappa * log(pistar))) !FS 0.5->1

          thetastar = T_bar / pistar

          pistrat(k) = pistrat(k - 1) - 0.5 * dz * (kappa / thetastar + kappa &
              &/ thetaStrat(k - 1))

          if(pistrat(k) <= 0.) then
            print *, 'ERROR: negative non-positive pistrat at k =', k
            stop
          end if

          T_bar = max(tp_strato, pistrat(k) * (tp_srf_trp - 0.5 &
              &* tpdiffhor_tropo - 0.5 * ptdiffvert_tropo / kappa &
              &* log(pistrat(k)))) !FS 0.5->1

          thetaStrat(k) = T_bar / pistrat(k)

          pStrat(k) = pistrat(k) ** ((1.0 - kappa) / kappa)

          rhoStrat(k) = pStrat(k) / thetaStrat(k)

        end do

        thetaStrat(- 1) = thetaStrat(0)
        pStrat(- 1) = pStrat(0)
        rhoStrat(- 1) = rhoStrat(0)
        pistrat(- 1) = pistrat(0)

        ! quantities at half levels

        ! with the exception of the uppermost ghost layer the
        ! values there are obtained by linear interpolation, in order
        ! to be consistent with the handling of the half levels in the
        ! semi-implicit time stepping and the corresponding pressure
        ! solver

        do k = - 1, nz + 2
          if(k < nz + 2) then
            PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k + 1))

            rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k + 1))

            thetaStratTilde(k) = PStratTilde(k) / rhoStratTilde(k)
          else
            ! linear extrapolation at uppermost layer

            PStratTilde(k) = 2.0 * PStratTilde(k - 1) - PStratTilde(k - 2)

            rhoStratTilde(k) = 2.0 * rhoStratTilde(k - 1) - rhoStratTilde(k - 2)

            thetaStratTilde(k) = 2.0 * thetaStratTilde(k - 1) &
                &- thetaStratTilde(k - 2)
          end if
        enddo

        ! Define 3D background fields.
        ! This implementation does not work yet.
        if(topography) then
          do i = - nbx, nx + nbx
            do j = - nby, ny + nby
              ! Define Exner pressure and 3D background fields at
              ! the surface.
              do k = 0, 1
                T_bar = max(tp_strato, tp_srf_trp - 0.5 * tpdiffhor_tropo)
                ! Define Exner pressure.
                piStratTFC(i, j, k) = 1.0 - heightTFC(i, j, k) * kappa / T_bar
                if(k == 1 .and. piStratTFC(i, j, k) <= 0.0) then
                  stop "ERROR: non-positive piStratTFC at k = 1"
                end if
                T_bar = max(tp_strato, piStratTFC(i, j, k) * (tp_srf_trp - 0.5 &
                    &* tpdiffhor_tropo - 0.5 * ptdiffvert_tropo / kappa &
                    &* log(piStratTFC(i, j, k))))
                ! Define pStratTFC.
                pStratTFC(i, j, k) = piStratTFC(i, j, k) ** ((1.0 - kappa) &
                    &/ kappa)
                ! Define thetaStratTFC.
                thetaStratTFC(i, j, k) = T_bar / piStratTFC(i, j, k)
                ! Define rhoStratTFC.
                rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                    &j, k)
              end do

              ! Integrate upwards.
              do k = 2, nz + 2
                ! Leapfrog step.
                piStar = piStratTFC(i, j, k - 2) - 2.0 * dz * jac(i, j, k - 1) &
                    &* kappa / thetaStratTFC(i, j, k - 1)
                if(piStar <= 0.0) then
                  print *, "ERROR: non-positive piStar at k =", k
                  stop
                end if
                T_bar = max(tp_strato, piStar * (tp_srf_trp - 0.5 &
                    &* tpdiffhor_tropo - 0.5 * ptdiffvert_tropo / kappa &
                    &* log(piStar)))
                thetaStar = T_bar / piStar
                ! Trapezoidal step.
                piStratTFC(i, j, k) = piStratTFC(i, j, k - 1) - 0.5 * dz * 0.5 &
                    &* (jac(i, j, k) + jac(i, j, k - 1)) * (kappa / thetaStar &
                    &+ kappa / thetaStratTFC(i, j, k - 1))
                if(piStratTFC(i, j, k) <= 0.0) then
                  print *, "ERROR: non-positive piStratTFC at k =", k
                  stop
                end if
                T_bar = max(tp_strato, piStratTFC(i, j, k) * (tp_srf_trp - 0.5 &
                    &* tpdiffhor_tropo - 0.5 * ptdiffvert_tropo / kappa &
                    &* log(piStratTFC(i, j, k))))
                ! Define pStratTFC.
                pStratTFC(i, j, k) = piStratTFC(i, j, k) ** ((1.0 - kappa) &
                    &/ kappa)
                ! Define thetaStratTFC.
                thetaStratTFC(i, j, k) = T_bar / piStratTFC(i, j, k)
                ! Define rhoStratTFC.
                rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                    &j, k)
              end do
              ! Adjust at the lower boundary.
              piStratTFC(i, j, - 1) = piStratTFC(i, j, 0)
              pStratTFC(i, j, - 1) = pStratTFC(i, j, 0)
              thetaStratTFC(i, j, - 1) = thetaStratTFC(i, j, 0)
              rhoStratTFC(i, j, - 1) = rhoStratTFC(i, j, 0)
            end do
          end do
        end if

        !------------------------------------------------------------------

      case default
        print *, "background = ", trim(background)
        stop "atmosphere.f90/init_background: background not defined"
      end select

      ! non-dimensional squared Brunt-Vaisala frequency
      ! (this could be done a bit nicer)

      bvsStrat(- 1) = g_ndim / thetaStrat(0) * (thetaStrat(1) - thetaStrat(0)) &
          &/ dz

      bvsStrat(0) = g_ndim / thetaStrat(0) * (thetaStrat(1) - thetaStrat(0)) &
          &/ dz

      N2 = max(bvsStrat(- 1), bvsStrat(0))

      do k = 1, nz
        bvsStrat(k) = g_ndim / thetaStrat(k) * (thetaStrat(k + 1) &
            &- thetaStrat(k - 1)) / (2.0 * dz)

        N2 = max(N2, bvsStrat(k))
      end do

      bvsStrat(nz + 1) = g_ndim / thetaStrat(nz + 1) * (thetaStrat(nz + 1) &
          &- thetaStrat(nz)) / dz

      bvsStrat(nz + 2) = bvsStrat(nz + 1)

      N2 = max(N2, bvsStrat(nz + 1))

      if(N2 < 0.) then
        stop 'ERROR: N2 < 0'
      else
        NN = sqrt(N2)
      end if

      if(TestCase == "smoothVortex") then
        bvsstrat(:) = 0.
      end if

      ! Define bvsStratTFC.
      if(topography) then
        bvsStratTFC = 0.0
        do i = - nbx, nx + nbx
          do j = - nby, ny + nby
            ! Lower boundary.
            bvsStratTFC(i, j, - 1) = g_ndim / thetaStratTFC(i, j, 0) / jac(i, &
                &j, 0) * (thetaStratTFC(i, j, 1) - thetaStratTFC(i, j, 0)) / dz
            bvsStratTFC(i, j, 0) = bvsStratTFC(i, j, - 1)
            ! Between boundaries.
            do k = 1, nz
              bvsStratTFC(i, j, k) = g_ndim / thetaStratTFC(i, j, k) / jac(i, &
                  &j, k) * 0.5 * (thetaStratTFC(i, j, k + 1) &
                  &- thetaStratTFC(i, j, k - 1)) / dz
            end do
            ! Upper boundary.
            bvsStratTFC(i, j, nz + 1) = g_ndim / thetaStratTFC(i, j, nz + 1) &
                &/ jac(i, j, nz + 1) * (thetaStratTFC(i, j, nz + 1) &
                &- thetaStratTFC(i, j, nz)) / dz
            bvsStratTFC(i, j, nz + 2) = bvsStratTFC(i, j, nz + 1)
          end do
        end do
        if(testCase == "smoothVortex") then
          bvsStratTFC = 0.0
        end if
      end if

    case("Boussinesq")

      !-------------------------------------------------------------------
      !             Atmospheres for Boussinesq
      !-------------------------------------------------------------------

      select case(background)

      case('uniform_Boussinesq') ! uniform atmosphere for testing

        rho00 = 1.0
        theta00 = theta0_dim / thetaRef
        P00 = rho00 * theta00
        NN = 0.0
        N2 = NN ** 2

        ! some stratified fields are still in use...
        rhoStrat = rho00
        rhoStratTilde = rho00 ! background density at half levels
        thetaStrat = theta00
        PStrat = P00

        ! TFC FJ
        pStratTilde = p00
        thetaStratTilde = theta00

      case('stratified_Boussinesq')

        rho00 = 1.0
        theta00 = theta0_dim / thetaRef
        P00 = rho00 * theta00
        NN = N_BruntVaisala_dim * tRef
        N2 = NN ** 2

        ! some stratified fields are still in use...
        rhoStrat = rho00
        rhoStratTilde = rho00 ! background density at half levels
        thetaStrat = theta00
        PStrat = P00

        ! TFC FJ
        pStratTilde = p00
        thetaStratTilde = theta00

      case default
        print *, "background = ", trim(background)
        stop "atmosphere.f90/init_background: background not defined"
      end select

      bvsStrat = N2

      ! Background fields.
      if(topography) then
        pStratTFC = p00
        thetaStratTFC = theta00
        rhoStratTFC = rho00
        bvsStratTFC = N2
      end if

    case default
      print *, "model = ", model
      stop "init_atmosphere: unknown case model."
    end select

  end subroutine init_atmosphere

  !---------------------------------------------------------------------------

  subroutine terminate_atmosphere
    !---------------------------------
    ! deallocate background varaibles
    !---------------------------------

    ! local variables
    integer :: allocstat

    !---------------- deallocate variables -----------------------

    deallocate(Pstrat, stat = allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not deallocate pStrat"

    deallocate(pistrat, stat = allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not deallocate pistrat"

    deallocate(PstratTilde, stat = allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not deallocate pStratTilde"

    deallocate(rhoStrat, stat = allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not deallocate pStrat"

    deallocate(rhoStratTilde, stat = allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not deallocate pStrat"

    deallocate(thetaStrat, stat = allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not dealloc thetaStrat"

    deallocate(thetaStratTilde, stat = allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not dealloc thetaStratTilde"

    deallocate(Ro, stat = allocstat) !FS
    if(allocstat /= 0) stop "atmosphere.f90: could not dealloc Ro"

    deallocate(RoInv, stat = allocstat) !FS
    if(allocstat /= 0) stop "atmosphere.f90: could not dealloc RoInv"

    ! Deallocate 3D background fields.
    if(topography) then
      deallocate(pStratTFC, stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not dealloc pStratTFC"

      deallocate(thetaStratTFC, stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not dealloc thetaStratTFC"

      deallocate(rhoStratTFC, stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not dealloc rhoStratTFC"

      deallocate(bvsStratTFC, stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not dealloc bvsStratTFC"

      deallocate(piStratTFC, stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not dealloc piStratTFC"
    end if

  end subroutine terminate_atmosphere

  !---------------------------------------------------------------------------

  subroutine setHalosOfField2D(field)

    ! Subroutine needed for halos of topography.

    !-------------------------------
    !  set values in halo cells
    !-------------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby), intent(inout) :: field

    ! auxiliary fields
    real, dimension(nbx, - nby:ny + nby) :: xSliceLeft_send, xSliceRight_send
    real, dimension(nbx, - nby:ny + nby) :: xSliceLeft_recv, xSliceRight_recv
    real, dimension(- nbx:nx + nbx, nby) :: ySliceBack_send, ySliceForw_send
    real, dimension(- nbx:nx + nbx, nby) :: ySliceBack_recv, ySliceForw_recv

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: i, j, k
    integer :: i0, j0, k0

    !-----------------------------
    !     find neighbour procs
    !-----------------------------

    if(idim > 1) call mpi_cart_shift(comm, 0, 1, left, right, ierror)
    if(jdim > 1) call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

    !------------------------------
    !          x-direction
    !------------------------------

    if(idim > 1) then

      ! slice size
      sendcount = nbx * (ny + 2 * nby + 1)
      recvcount = sendcount

      ! read slice into contiguous array
      do i = 1, nbx
        xSliceLeft_send(i, :) = field(i, :)
        xSliceRight_send(i, :) = field(nx - nbx + i, :)
      end do

      ! left -> right
      source = left
      dest = right
      tag = 100

      i0 = 1; j0 = - nby

      call mpi_sendrecv(xSliceRight_send(i0, j0), sendcount, &
          &mpi_double_precision, dest, tag, xSliceLeft_recv(i0, j0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xSliceLeft_send(i0, j0), sendcount, &
          &mpi_double_precision, dest, tag, xSliceRight_recv(i0, j0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_right, ierror)

      ! write auxiliary slice to var field
      do i = 1, nbx
        ! right halos
        field(nx + i, :) = xSliceRight_recv(i, :)
        ! left halos
        field(- nbx + i, :) = xSliceLeft_recv(i, :)
      end do

    else

      do i = 1, nbx
        field(nx + i, :) = field(i, :)
        field(- i + 1, :) = field(nx - i + 1, :)
      end do

    end if

    !------------------------------
    !          y-direction
    !------------------------------

    if(jdim > 1) then

      ! slice size
      sendcount = nby * (nx + 2 * nbx + 1)
      recvcount = sendcount

      ! read slice into contiguous array
      do j = 1, nby
        ySliceBack_send(:, j) = field(:, j)
        ySliceForw_send(:, j) = field(:, ny - nby + j)
      end do

      ! back -> forw
      source = back
      dest = forw
      tag = 100

      i0 = - nbx; j0 = 1

      call mpi_sendrecv(ySliceForw_send(i0, j0), sendcount, &
          &mpi_double_precision, dest, tag, ySliceBack_recv(i0, j0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(ySliceBack_send(i0, j0), sendcount, &
          &mpi_double_precision, dest, tag, ySliceForw_recv(i0, j0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_forw, ierror)

      ! write auxiliary slice to var field
      do j = 1, nby
        ! right halos
        field(:, ny + j) = ySliceForw_recv(:, j)
        ! left halos
        field(:, - nby + j) = ySliceBack_recv(:, j)
      end do

    else

      do j = 1, nby
        field(:, ny + j) = field(:, j)
        field(:, - j + 1) = field(:, ny - j + 1)
      end do

    end if

  end subroutine setHalosOfField2D

  !---------------------------------------------------------------------------

  subroutine setup_topography

    real :: mountainHeight, mountainWidth, mountainWavenumber
    real :: x_center, y_center
    real :: kk, ll
    integer :: ix0, jy0
    integer :: ix, jy
    integer :: iwm

    if(lz(0) /= 0.0) stop "Error in setup_topography: lz(0) must be zero for &
        &TFC!"

    if(zBoundary == "periodic" .and. mountainHeight_dim /= 0.0) stop "Error in &
        &setup_topography: mountainHeight_dim must be zero for zBoundary &
        &= 'periodic'"

    mountainHeight = mountainHeight_dim / lRef
    mountainWidth = mountainWidth_dim / lRef
    mountainWavenumber = pi / mountainWidth

    x_center = 0.5 * (lx(1) + lx(0))
    y_center = 0.5 * (ly(1) + ly(0))

    ix0 = is + nbx - 1
    jy0 = js + nby - 1

    if(rayTracer .and. case_wkb == 3) then
      if(nwm < 1 .or. (mountain_case == 13 .and. nwm < spectral_modes)) stop &
          &"Error in setup_topography: nwm is too small!"

      if(mountain_case /= 0) then
        k_spectrum = 0.0
        l_spectrum = 0.0
        topography_spectrum = 0.0
        topography_surface = 0.0

        do jy = 1, ny
          do ix = 1, nx
            select case(mountain_case)
            case(1)
              ! 2D cosine mountains
              k_spectrum(ix, jy, 1) = mountainWavenumber
              topography_spectrum(ix, jy, 1) = 0.5 * mountainHeight
              topography_surface(ix, jy) = 0.5 * mountainHeight

            case(5)
              ! 2D cosine envelope and even background
              if(abs(x(ix + ix0) - x_center) <= mountainWidth * range_factor) &
                  &then
                k_spectrum(ix, jy, 1) = mountainWavenumber
                topography_spectrum(ix, jy, 1) = 0.25 * mountainHeight * (1.0 &
                    &+ cos(mountainWavenumber / range_factor * (x(ix + ix0) &
                    &- x_center)))
              end if
              topography_surface(ix, jy) = 0.5 * mountainHeight

            case(7)
              ! 2D Gaussian envelope and even background
              k_spectrum(ix, jy, 1) = mountainWavenumber
              topography_spectrum(ix, jy, 1) = 0.5 * mountainHeight * exp(- &
                  &(x(ix + ix0) - x_center) ** 2.0 / (mountainWidth &
                  &* range_factor) ** 2.0)
              topography_surface(ix, jy) = 0.5 * mountainHeight

            case(9)
              ! 2D cosine envelope and cosine background
              if(abs(x(ix + ix0) - x_center) <= mountainWidth * range_factor) &
                  &then
                k_spectrum(ix, jy, 1) = mountainWavenumber
                topography_spectrum(ix, jy, 1) = 0.25 * mountainHeight * (1.0 &
                    &+ cos(mountainWavenumber / range_factor * (x(ix + ix0) &
                    &- x_center)))
                topography_surface(ix, jy) = 0.25 * mountainHeight * (1.0 &
                    &+ cos(mountainWavenumber / range_factor * (x(ix + ix0) &
                    &- x_center)))
              end if

            case(11)
              ! 2D Gaussian envelope and Gaussian background
              k_spectrum(ix, jy, 1) = mountainWavenumber
              topography_spectrum(ix, jy, 1) = 0.5 * mountainHeight * exp(- &
                  &(x(ix + ix0) - x_center) ** 2.0 / (mountainWidth &
                  &* range_factor) ** 2.0)
              topography_surface(ix, jy) = 0.5 * mountainHeight * exp(- (x(ix &
                  &+ ix0) - x_center) ** 2.0 / (mountainWidth * range_factor) &
                  &** 2.0)

            case(13)
              ! 3D WKB topography
              if((x(ix + ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                  &** 2.0 <= (mountainWidth * range_factor) ** 2.0) then
                do iwm = 0, spectral_modes - 1
                  k_spectrum(ix, jy, iwm + 1) = mountainWavenumber * cos(pi &
                      &/ spectral_modes * iwm)
                  l_spectrum(ix, jy, iwm + 1) = mountainWavenumber * sin(pi &
                      &/ spectral_modes * iwm)
                  topography_spectrum(ix, jy, iwm + 1) = 0.25 * mountainHeight &
                      &* (1.0 + cos(mountainWavenumber / range_factor &
                      &* sqrt((x(ix + ix0) - x_center) ** 2.0 + (y(jy + jy0) &
                      &- y_center) ** 2.0))) / spectral_modes
                end do
                topography_surface(ix, jy) = 0.25 * mountainHeight * (1.0 &
                    &+ cos(mountainWavenumber / range_factor * sqrt((x(ix &
                    &+ ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                    &** 2.0)))
              end if

            case default
              stop "Error in setup_topography: Unknown mountain case!"
            end select
          end do
        end do
      else
        k_spectrum = k_spectrum * lRef
        l_spectrum = l_spectrum * lRef
        topography_spectrum = topography_spectrum / lRef
        topography_surface = topography_surface / lRef
      end if

      call setHalosOfField2D(topography_surface)

      if(topographyTime > 0.0) then
        final_topography_spectrum = topography_spectrum
        topography_spectrum = 0.0
        final_topography_surface = topography_surface
        topography_surface = 0.0
      end if
    else if(topography) then
      if(mountain_case /= 0) then
        topography_surface = 0.0
        do jy = 1, ny
          do ix = 1, nx
            select case(mountain_case)
            case(1)
              ! 2D cosine mountains
              topography_surface(ix, jy) = 0.5 * mountainHeight * (1.0 &
                  &+ cos(mountainWavenumber * (x(ix + ix0) - x_center)))

            case(2)
              ! 3D cosine mountains
              topography_surface(ix, jy) = 0.5 * mountainHeight * (1.0 &
                  &+ cos(mountainWavenumber * sqrt((x(ix + ix0) - x_center) &
                  &** 2.0 + (y(jy + jy0) - y_center) ** 2.0)))

            case(3)
              ! 2D isolated mountain
              topography_surface(ix, jy) = mountainHeight / (1.0 + (x(ix &
                  &+ ix0) - x_center) ** 2.0 / mountainWidth ** 2.0)

            case(4)
              ! 3D isolated mountain
              topography_surface(ix, jy) = mountainHeight / (1.0 + ((x(ix &
                  &+ ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                  &** 2.0) / mountainWidth ** 2.0)

            case(5)
              ! 2D cosine envelope and even background
              if(abs(x(ix + ix0) - x_center) <= mountainWidth * range_factor) &
                  &then
                topography_surface(ix, jy) = 0.5 * mountainHeight * (1.0 + 0.5 &
                    &* (1.0 + cos(mountainWavenumber / range_factor * (x(ix &
                    &+ ix0) - x_center))) * cos(mountainWavenumber * (x(ix &
                    &+ ix0) - x_center)))
              else
                topography_surface(ix, jy) = 0.5 * mountainHeight
              end if

            case(6)
              ! 3D cosine envelope and even background
              if((x(ix + ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                  &** 2.0 <= (mountainWidth * range_factor) ** 2.0) then
                topography_surface(ix, jy) = 0.5 * mountainHeight * (1.0 + 0.5 &
                    &* (1.0 + cos(mountainWavenumber / range_factor &
                    &* sqrt((x(ix + ix0) - x_center) ** 2.0 + (y(jy + jy0) &
                    &- y_center) ** 2.0))) * cos(mountainWavenumber &
                    &* sqrt((x(ix + ix0) - x_center) ** 2.0 + (y(jy + jy0) &
                    &- y_center) ** 2.0)))
              else
                topography_surface(ix, jy) = 0.5 * mountainHeight
              end if

            case(7)
              ! 2D Gaussian envelope and even background
              topography_surface(ix, jy) = 0.5 * mountainHeight * (1.0 + exp(- &
                  &(x(ix + ix0) - x_center) ** 2.0 / (mountainWidth &
                  &* range_factor) ** 2.0) * cos(mountainWavenumber * (x(ix &
                  &+ ix0) - x_center)))

            case(8)
              ! 3D Gaussian envelope and even background
              topography_surface(ix, jy) = 0.5 * mountainHeight * (1.0 + exp(- &
                  &((x(ix + ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                  &** 2.0) / (mountainWidth * range_factor) ** 2.0) &
                  &* cos(mountainWavenumber * sqrt((x(ix + ix0) - x_center) &
                  &** 2.0 + (y(jy + jy0) - y_center) ** 2.0)))

            case(9)
              ! 2D cosine envelope and cosine background
              if(abs(x(ix + ix0) - x_center) <= mountainWidth * range_factor) &
                  &then
                topography_surface(ix, jy) = 0.25 * mountainHeight * (1.0 &
                    &+ cos(mountainWavenumber / range_factor * (x(ix + ix0) &
                    &- x_center))) * (1.0 + cos(mountainWavenumber * (x(ix &
                    &+ ix0) - x_center)))
              end if

            case(10)
              ! 3D cosine envelope and cosine background
              if((x(ix + ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                  &** 2.0 <= (mountainWidth * range_factor) ** 2.0) then
                topography_surface(ix, jy) = 0.25 * mountainHeight * (1.0 &
                    &+ cos(mountainWavenumber / range_factor * sqrt((x(ix &
                    &+ ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                    &** 2.0))) * (1.0 + cos(mountainWavenumber * sqrt((x(ix &
                    &+ ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                    &** 2.0)))
              end if

            case(11)
              ! 2D Gaussian envelope and Gaussian background
              topography_surface(ix, jy) = 0.5 * mountainHeight * exp(- (x(ix &
                  &+ ix0) - x_center) ** 2.0 / (mountainWidth * range_factor) &
                  &** 2.0) * (1.0 + cos(mountainWavenumber * (x(ix + ix0) &
                  &- x_center)))

            case(12)
              ! 3D Gaussian envelope and Gaussian background
              topography_surface(ix, jy) = 0.5 * mountainHeight * exp(- ((x(ix &
                  &+ ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                  &** 2.0) / (mountainWidth * range_factor) ** 2.0) * (1.0 &
                  &+ cos(mountainWavenumber * sqrt((x(ix + ix0) - x_center) &
                  &** 2.0 + (y(jy + jy0) - y_center) ** 2.0)))

            case(13)
              ! 3D WKB topography
              if((x(ix + ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                  &** 2.0 <= (mountainWidth * range_factor) ** 2.0) then
                topography_surface(ix, jy) = 0.25 * mountainHeight * (1.0 &
                    &+ cos(mountainWavenumber / range_factor * sqrt((x(ix &
                    &+ ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                    &** 2.0)))
                do iwm = 0, spectral_modes - 1
                  kk = mountainWavenumber * cos(pi / spectral_modes * iwm)
                  ll = mountainWavenumber * sin(pi / spectral_modes * iwm)
                  topography_surface(ix, jy) = topography_surface(ix, jy) &
                      &+ 0.25 * mountainHeight * (1.0 + cos(mountainWavenumber &
                      &/ range_factor * sqrt((x(ix + ix0) - x_center) ** 2.0 &
                      &+ (y(jy + jy0) - y_center) ** 2.0))) * cos(kk * (x(ix &
                      &+ ix0) - x_center) + ll * (y(jy + jy0) - y_center)) &
                      &/ spectral_modes
                end do
              end if

            case default
              stop "Error in setup_topography: Unknown mountain case!"
            end select
          end do
        end do
      else
        topography_surface = topography_surface / lRef
      end if

      call setHalosOfField2D(topography_surface)

      if(topographyTime > 0.0) then
        final_topography_surface = topography_surface
        topography_surface = 0.0
      end if
    end if

  end subroutine setup_topography

  !---------------------------------------------------------------------------

  subroutine update_topography(time)

    real, intent(in) :: time

    real :: z_tr
    real :: theta_tr
    real :: press_tr
    real :: T_tr

    real :: power

    real :: pow_t, pow_s, p_t_b
    real :: T_bar, T_c_b, p_bar

    integer :: i, j, k

    if(topographyTime <= 0.0) return

    if(rayTracer .and. case_wkb == 3) then
      if(any(topography_spectrum /= final_topography_spectrum)) then
        if(time < topographyTime / tRef) then
          topography_spectrum = time / topographyTime * tRef &
              &* final_topography_spectrum
        else
          topography_spectrum = final_topography_spectrum
        end if
      end if
    end if

    if(any(topography_surface /= final_topography_surface)) then

      ! Update topography.
      if(time < topographyTime / tRef) then
        topography_surface = time / topographyTime * tRef &
            &* final_topography_surface
      else
        topography_surface = final_topography_surface
      end if

      ! Return if model is Boussinesq.
      if(model == "Boussinesq") return

      ! Update TFC background if necessary.
      select case(background)

      case("realistic")

        z_tr = z_tr_dim / lRef
        theta_tr = theta_tr_dim / thetaRef
        press_tr = p0 * (1.0 - kappa * sig / theta_tr * z_tr) ** (1 / kappa)
        T_tr = theta_tr * (press_tr / p0) ** kappa

        do i = - nbx, nx + nbx
          do j = - nby, ny + nby
            do k = - 1, nz + 2
              if(heightTFC(i, j, k) <= z_tr) then
                ! Isentropic troposphere.
                power = 1.0 / (gamma - 1.0)
                term = kappa * sig / theta_tr * heightTFC(i, j, k)
                if(term > 1.0) then
                  stop "init_atmosphere: negative term. Stop"
                end if
                ! Define pStratTFC.
                pStratTFC(i, j, k) = p0 * (1.0 - term) ** power
                ! Define thetaStratTFC.
                thetaStratTFC(i, j, k) = theta_tr
                ! Define rhoStratTFC.
                rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                    &j, k)
              else
                ! Isothermal stratosphere.
                ! Define pStratTFC.
                pStratTFC(i, j, k) = p0 ** kappa * press_tr ** (1.0 / gamma) &
                    &* exp(- sig / gamma / T_tr * (heightTFC(i, j, k) - z_tr))
                ! Define thetaStratTFC.
                thetaStratTFC(i, j, k) = theta_tr * exp(kappa * sig / T_tr &
                    &* (heightTFC(i, j, k) - z_tr))
                ! Define rhoStratTFC.
                rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                    &j, k)
              end if
            end do
          end do
        end do

      case("isothermal")

        T0 = Temp0_dim / thetaRef

        do i = - nbx, nx + nbx
          do j = - nby, ny + nby
            do k = - 1, nz + 2
              ! Define pStratTFC.
              pStratTFC(i, j, k) = p0 * exp(- sig * heightTFC(i, j, k) / gamma &
                  &/ T0)
              ! Define thetaStratTFC.
              thetaStratTFC(i, j, k) = T0 * exp(kappa * sig / T0 &
                  &* heightTFC(i, j, k))
              ! Define rhoStratTFC.
              rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, j, k)
            end do
          end do
        end do

      case("isentropic")

        theta0 = theta0_dim / thetaRef

        do i = - nbx, nx + nbx
          do j = - nby, ny + nby
            do k = - 1, nz + 2
              power = 1.0 / (gamma - 1.0)
              term = kappa * sig / theta0 * heightTFC(i, j, k)
              if(term > 1.0) then
                stop "init_atmosphere: negative term. Stop"
              end if
              ! Define pStratTFC.
              pStratTFC(i, j, k) = p0 * (1.0 - term) ** power
              ! Define thetaStratTFC.
              thetaStratTFC(i, j, k) = theta0
              ! Define rhoStratTFC.
              rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, j, k)
            end do
          end do
        end do

      case("const-N")

        theta0 = theta0_dim / thetaRef
        N2 = (N_BruntVaisala_dim * tRef) ** 2
        NN = sqrt(N2)

        do i = - nbx, nx + nbx
          do j = - nby, ny + nby
            do k = - 1, nz + 2
              power = 1.0 / (gamma - 1.0)
              term = exp(- Fr2 * N2 * heightTFC(i, j, k))
              ! Define pStratTFC.
              pStratTFC(i, j, k) = p0 * (1.0 + FrInv2 * kappa * sig / N2 &
                  &/ theta0 * (term - 1.0)) ** power
              ! Define thetaStratTFC.
              thetaStratTFC(i, j, k) = theta0 * exp(Fr ** 2.0 * N2 &
                  &* heightTFC(i, j, k))
              ! Define rhoStratTFC.
              rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, j, k)
            end do
          end do
        end do

      case("diflapse")

        if(gamma_t /= 0.) pow_t = g / (Rsp * gamma_t)
        if(gamma_s /= 0.) pow_s = g / (Rsp * gamma_s)

        if(gamma_t /= 0.) then
          p_t_b = press0_dim * (1. - gamma_t * z_tr_dim / Temp0_dim) ** pow_t
        else
          p_t_b = press0_dim * exp(- z_tr_dim / (Rsp * Temp0_dim / g))
        end if

        T_c_b = Temp0_dim - gamma_t * z_tr_dim

        do i = - nbx, nx + nbx
          do j = - nby, ny + nby
            do k = - 1, nz + 2
              zk = heightTFC(i, j, k) * lRef

              if(zk < z_tr_dim) then
                T_bar = Temp0_dim - gamma_t * zk
                if(gamma_t /= 0.0) then
                  p_bar = press0_dim * (1.0 - gamma_t * zk / Temp0_dim) ** pow_t
                else
                  p_bar = press0_dim * exp(- zk / (Rsp * Temp0_dim / g))
                end if
              else
                T_bar = Temp0_dim - gamma_t * z_tr_dim - gamma_s * (zk &
                    &- z_tr_dim)
                if(gamma_s /= 0.0) then
                  p_bar = p_t_b * (1.0 - gamma_s * (zk - z_tr_dim) / T_c_b) &
                      &** pow_s
                else
                  p_bar = p_t_b * exp(- (zk - z_tr_dim) / (Rsp * T_c_b / g))
                end if
              end if

              ! Define thetaStratTFC.
              thetaStratTFC(i, j, k) = T_bar * (press0_dim / p_bar) ** kappa &
                  &/ thetaRef
              ! Define rhoStratTFC.
              rhoStratTFC(i, j, k) = p_bar / (Rsp * T_bar) / rhoRef
              ! Define pStratTFC.
              pStratTFC(i, j, k) = rhoStratTFC(i, j, k) * thetaStratTFC(i, j, k)
            end do
          end do
        end do

      end select

      bvsStratTFC = 0.0
      do i = - nbx, nx + nbx
        do j = - nby, ny + nby
          ! Lower boundary.
          bvsStratTFC(i, j, - 1) = g_ndim / thetaStratTFC(i, j, 0) / jac(i, j, &
              &0) * (thetaStratTFC(i, j, 1) - thetaStratTFC(i, j, 0)) / dz
          bvsStratTFC(i, j, 0) = bvsStratTFC(i, j, - 1)
          ! Between boundaries.
          do k = 1, nz
            bvsStratTFC(i, j, k) = g_ndim / thetaStratTFC(i, j, k) / jac(i, j, &
                &k) * 0.5 * (thetaStratTFC(i, j, k + 1) - thetaStratTFC(i, j, &
                &k - 1)) / dz
          end do
          ! Upper boundary.
          bvsStratTFC(i, j, nz + 1) = g_ndim / thetaStratTFC(i, j, nz + 1) &
              &/ jac(i, j, nz + 1) * (thetaStratTFC(i, j, nz + 1) &
              &- thetaStratTFC(i, j, nz)) / dz
          bvsStratTFC(i, j, nz + 2) = bvsStratTFC(i, j, nz + 1)
        end do
      end do
      if(testCase == "smoothVortex") then
        bvsStratTFC = 0.0
      end if
    end if

  end subroutine update_topography

  !---------------------------------------------------------------------------

  subroutine add_JP_to_u(var, option)
    !--------------------------------
    ! in the compressible model the Euler steps are solved by updating Pu
    ! instead of u, so we write Pu in var
    !-------------------------------

    character(len = *), intent(in) :: option
    ! in/out variables
    type(var_type), intent(inout) :: var
    ! local variables
    integer :: i, j, k

    select case(option)

    case("forward")
      do k = - nbz, nz + nbz - 1
        do j = - nby, ny + nby - 1
          do i = - nbx, nx + nbx - 1
            var%w(i, j, k) = var%w(i, j, k) * 0.5 * (jac(i, j, k) * var%P(i, &
                &j, k) + jac(i, j, k + 1) * var%P(i, j, k + 1))

            var%v(i, j, k) = var%v(i, j, k) * 0.5 * (jac(i, j, k) * var%P(i, &
                &j, k) + jac(i, j + 1, k) * var%P(i, j + 1, k))

            var%u(i, j, k) = var%u(i, j, k) * 0.5 * (jac(i, j, k) * var%P(i, &
                &j, k) + jac(i + 1, j, k) * var%P(i + 1, j, k))
          end do
          var%w(nx + nbx, j, k) = 2 * var%w(nx + nbx - 2, j, k) - var%w(nx &
              &+ nbx - 1, j, k)
          var%v(nx + nbx, j, k) = 2 * var%v(nx + nbx - 2, j, k) - var%v(nx &
              &+ nbx - 1, j, k)
          var%u(nx + nbx, j, k) = 2 * var%u(nx + nbx - 2, j, k) - var%u(nx &
              &+ nbx - 1, j, k)
        end do
        var%w(:, ny + nby, k) = 2 * var%w(:, ny + nby - 2, k) - var%w(:, ny &
            &+ nby - 1, k)
        var%v(:, j, ny + nby) = 2 * var%v(:, ny + nby - 2, k) - var%v(:, ny &
            &+ nby - 1, k)
        var%u(:, ny + nby, k) = 2 * var%u(:, ny + nby - 2, k) - var%u(:, ny &
            &+ nby - 1, k)
      end do
      var%w(:, :, nz + nbz) = 2 * var%w(:, :, nz + nbz - 2) - var%w(:, :, nz &
          &+ nbz - 1)
      var%v(:, :, nz + nbz) = 2 * var%v(:, :, nz + nbz - 2) - var%v(:, :, nz &
          &+ nbz - 1)
      var%u(:, :, nz + nbz) = 2 * var%u(:, :, nz + nbz - 2) - var%u(:, :, nz &
          &+ nbz - 1)

    case("backward")
      do k = - nbz, nz + nbz - 1
        do j = - nby, ny + nby - 1
          do i = - nbx, nx + nbx - 1
            var%w(i, j, k) = var%w(i, j, k) / (0.5 * (jac(i, j, k) * var%P(i, &
                &j, k) + jac(i, j, k + 1) * var%P(i, j, k + 1)))

            var%v(i, j, k) = var%v(i, j, k) / (0.5 * (jac(i, j, k) * var%P(i, &
                &j, k) + jac(i, j + 1, k) * var%P(i, j + 1, k)))

            var%u(i, j, k) = var%u(i, j, k) / (0.5 * (jac(i, j, k) * var%P(i, &
                &j, k) + jac(i + 1, j, k) * var%P(i + 1, j, k)))
          end do
          var%w(nx + nbx, j, k) = 2 * var%w(nx + nbx - 2, j, k) - var%w(nx &
              &+ nbx - 1, j, k)
          var%v(nx + nbx, j, k) = 2 * var%v(nx + nbx - 2, j, k) - var%v(nx &
              &+ nbx - 1, j, k)
          var%u(nx + nbx, j, k) = 2 * var%u(nx + nbx - 2, j, k) - var%u(nx &
              &+ nbx - 1, j, k)
        end do
        var%w(:, ny + nby, k) = 2 * var%w(:, ny + nby - 2, k) - var%w(:, ny &
            &+ nby - 1, k)
        var%v(:, j, ny + nby) = 2 * var%v(:, ny + nby - 2, k) - var%v(:, ny &
            &+ nby - 1, k)
        var%u(:, ny + nby, k) = 2 * var%u(:, ny + nby - 2, k) - var%u(:, ny &
            &+ nby - 1, k)
      end do
      var%w(:, :, nz + nbz) = 2 * var%w(:, :, nz + nbz - 2) - var%w(:, :, nz &
          &+ nbz - 1)
      var%v(:, :, nz + nbz) = 2 * var%v(:, :, nz + nbz - 2) - var%v(:, :, nz &
          &+ nbz - 1)
      var%u(:, :, nz + nbz) = 2 * var%u(:, :, nz + nbz - 2) - var%u(:, :, nz &
          &+ nbz - 1)

    case default
      stop "add_JP_to_u unknown option."
    end select

  end subroutine

  !---------------------------------------------------------------------------

  function jac(i, j, k)
    ! Jacobian.

    real :: jac
    integer :: i, j, k

    jac = ((lz(1) - lz(0)) - topography_surface(i, j)) / (lz(1) - lz(0))

  end function jac

  function met(i, j, k, mu, nu)
    ! Metric tensor.

    real :: met
    integer :: i, j, k, mu, nu

    if((mu == 1 .and. nu == 3) .or. (mu == 3 .and. nu == 1)) then
      met = (topography_surface(i + 1, j) - topography_surface(i - 1, j)) &
          &/ (2.0 * dx) * (z(k) - (lz(1) - lz(0))) / ((lz(1) - lz(0)) &
          &- topography_surface(i, j))
    else if((mu == 2 .and. nu == 3) .or. (mu == 3 .and. nu == 2)) then
      met = (topography_surface(i, j + 1) - topography_surface(i, j - 1)) &
          &/ (2.0 * dy) * (z(k) - (lz(1) - lz(0))) / ((lz(1) - lz(0)) &
          &- topography_surface(i, j))
    else if(mu == 3 .and. nu == 3) then
      met = ((lz(1) - lz(0)) / ((lz(1) - lz(0)) - topography_surface(i, j))) &
          &** 2.0 + ((z(k) - (lz(1) - lz(0))) / ((lz(1) - lz(0)) &
          &- topography_surface(i, j))) ** 2.0 * (((topography_surface(i + 1, &
          &j) - topography_surface(i - 1, j)) / (2.0 * dx)) ** 2.0 &
          &+ ((topography_surface(i, j + 1) - topography_surface(i, j - 1)) &
          &/ (2.0 * dy)) ** 2.0)
    end if

  end function met

  function chris(i, j, k, mu, nu)
    ! Christophel tensor.

    real :: chris
    integer :: i, j, k, mu, nu

    if(mu == 1 .and. nu == 1) then
      chris = - (topography_surface(i - 1, j) - 2.0 * topography_surface(i, j) &
          &+ topography_surface(i + 1, j)) / (dx ** 2.0) * (z(k) - (lz(1) &
          &- lz(0))) / ((lz(1) - lz(0)) - topography_surface(i, j))
    else if((mu == 1 .and. nu == 2) .or. (mu == 2 .and. nu == 1)) then
      chris = - (topography_surface(i + 1, j + 1) - topography_surface(i - 1, &
          &j + 1) - topography_surface(i + 1, j - 1) + topography_surface(i &
          &- 1, j - 1)) / (4.0 * dx * dy) * (z(k) - (lz(1) - lz(0))) / ((lz(1) &
          &- lz(0)) - topography_surface(i, j))
    else if(mu == 2 .and. nu == 2) then
      chris = - (topography_surface(i, j - 1) - 2.0 * topography_surface(i, j) &
          &+ topography_surface(i, j + 1)) / (dy ** 2.0) * (z(k) - (lz(1) &
          &- lz(0))) / ((lz(1) - lz(0)) - topography_surface(i, j))
    else if((mu == 1 .and. nu == 3) .or. (mu == 3 .and. nu == 1)) then
      chris = - (topography_surface(i + 1, j) - topography_surface(i - 1, j)) &
          &/ (2.0 * dx) / ((lz(1) - lz(0)) - topography_surface(i, j))
    else if((mu == 2 .and. nu == 3) .or. (mu == 3 .and. nu == 2)) then
      chris = - (topography_surface(i, j + 1) - topography_surface(i, j - 1)) &
          &/ (2.0 * dy) / ((lz(1) - lz(0)) - topography_surface(i, j))
    end if

  end function chris

  function heightTFC(i, j, k)
    ! Transformation of vertical coordinate.

    real :: heightTFC
    integer :: i, j, k

    heightTFC = z(k) * ((lz(1) - lz(0)) - topography_surface(i, j)) / (lz(1) &
        &- lz(0)) + topography_surface(i, j)

  end function heightTFC

  function levelTFC(i, j, height)

    real :: levelTFC
    integer :: i, j
    real :: height

    levelTFC = (lz(1) - lz(0)) * (height - topography_surface(i, j)) / ((lz(1) &
        &- lz(0)) - topography_surface(i, j))

  end function levelTFC

  function vertWindTFC(i, j, k, var)
    ! Transformation of the vertical wind.

    type(var_type) :: var
    integer :: i, j, k

    real :: vertWindTFC
    real :: uEdgeR, uUEdgeR, uEdgeL, uUEdgeL, vEdgeF, vUEdgeF, vEdgeB, &
        &vUEdgeB, wEdgeU

    uEdgeR = var%u(i, j, k)
    uUEdgeR = var%u(i, j, k + 1)
    uEdgeL = var%u(i - 1, j, k)
    uUEdgeL = var%u(i - 1, j, k + 1)
    vEdgeF = var%v(i, j, k)
    vUEdgeF = var%v(i, j, k + 1)
    vEdgeB = var%v(i, j - 1, k)
    vUEdgeB = var%v(i, j - 1, k + 1)
    wEdgeU = var%w(i, j, k)

    vertWindTFC = trafoTFC(i, j, k, uEdgeR, uUEdgeR, uEdgeL, uUEdgeL, vEdgeF, &
        &vUEdgeF, vEdgeB, vUEdgeB, wEdgeU, "car")

  end function vertWindTFC

  function trafoTFC(i, j, k, uEdgeR, uUEdgeR, uEdgeL, uUEdgeL, vEdgeF, &
      &vUEdgeF, vEdgeB, vUEdgeB, wEdgeU, wind)

    integer :: i, j, k
    real :: uEdgeR, uUEdgeR, uEdgeL, uUEdgeL, vEdgeF, vUEdgeF, vEdgeB, &
        &vUEdgeB, wEdgeU
    character(len = 3) :: wind

    real :: trafoTFC
    real :: jacEdgeU
    real :: uC, uU, vC, vU
    real :: metEdgeR, metUEdgeR, metEdgeL, metUEdgeL, metEdgeF, metUEdgeF, &
        &metEdgeB, metUEdgeB
    real :: metEdgeRU, metEdgeLU, metEdgeFU, metEdgeBU
    real :: uEdgeRU, uEdgeLU, vEdgeFU, vEdgeBU

    select case(ipolTFC)

    case(1)

      ! Multiplication on rho grid

      jacEdgeU = 0.5 * (jac(i, j, k) + jac(i, j, k + 1))
      uC = 0.5 * (uEdgeR + uEdgeL)
      uU = 0.5 * (uUEdgeR + uUEdgeL)
      vC = 0.5 * (vEdgeF + vEdgeB)
      vU = 0.5 * (vUEdgeF + vUEdgeB)
      if(wind == "car") then
        trafoTFC = - 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * uC + jac(i, j, &
            &k + 1) * met(i, j, k + 1, 1, 3) * uU) - 0.5 * (jac(i, j, k) &
            &* met(i, j, k, 2, 3) * vC + jac(i, j, k + 1) * met(i, j, k + 1, &
            &2, 3) * vU) + jacEdgeU * wEdgeU
      else if(wind == "tfc") then
        trafoTFC = (0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * uC + jac(i, j, &
            &k + 1) * met(i, j, k + 1, 1, 3) * uU) + 0.5 * (jac(i, j, k) &
            &* met(i, j, k, 2, 3) * vC + jac(i, j, k + 1) * met(i, j, k + 1, &
            &2, 3) * vU) + wEdgeU) / jacEdgeU
      end if

    case(2)

      ! Multiplication on rho grid (inverted)

      jacEdgeU = 0.5 * (jac(i, j, k) + jac(i, j, k + 1))
      uC = 0.5 * (uEdgeR + uEdgeL)
      uU = 0.5 * (uUEdgeR + uUEdgeL)
      vC = 0.5 * (vEdgeF + vEdgeB)
      vU = 0.5 * (vUEdgeF + vUEdgeB)
      if(wind == "car") then
        trafoTFC = jacEdgeU * (- 0.5 * (met(i, j, k, 1, 3) * uC + met(i, j, k &
            &+ 1, 1, 3) * uU) - 0.5 * (met(i, j, k, 2, 3) * vC + met(i, j, k &
            &+ 1, 2, 3) * vU) + wEdgeU)
      else if(wind == "tfc") then
        trafoTFC = 0.5 * (met(i, j, k, 1, 3) * uC + met(i, j, k + 1, 1, 3) &
            &* uU) + 0.5 * (met(i, j, k, 2, 3) * vC + met(i, j, k + 1, 2, 3) &
            &* vU) + wEdgeU / jacEdgeU
      end if

    case(3)

      ! Multiplication on u grid

      jacEdgeU = 0.5 * (jac(i, j, k) + jac(i, j, k + 1))
      metEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) + jac(i + 1, j, k) &
          &* met(i + 1, j, k, 1, 3))
      metUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) + jac(i &
          &+ 1, j, k + 1) * met(i + 1, j, k + 1, 1, 3))
      metEdgeL = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) + jac(i - 1, j, k) &
          &* met(i - 1, j, k, 1, 3))
      metUEdgeL = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) + jac(i &
          &- 1, j, k + 1) * met(i - 1, j, k + 1, 1, 3))
      metEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) + jac(i, j + 1, k) &
          &* met(i, j + 1, k, 2, 3))
      metUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) + jac(i, j &
          &+ 1, k + 1) * met(i, j + 1, k + 1, 2, 3))
      metEdgeB = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) + jac(i, j - 1, k) &
          &* met(i, j - 1, k, 2, 3))
      metUEdgeB = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) + jac(i, j &
          &- 1, k + 1) * met(i, j - 1, k + 1, 2, 3))
      if(wind == "car") then
        trafoTFC = - 0.25 * (metEdgeR * uEdgeR + metUEdgeR * uUEdgeR &
            &+ metEdgeL * uEdgeL + metUEdgeL * uUEdgeL) - 0.25 * (metEdgeF &
            &* vEdgeF + metUEdgeF * vUEdgeF + metEdgeB * vEdgeB + metUEdgeB &
            &* vUEdgeB) + jacEdgeU * wEdgeU
      else if(wind == "tfc") then
        trafoTFC = (0.25 * (metEdgeR * uEdgeR + metUEdgeR * uUEdgeR + metEdgeL &
            &* uEdgeL + metUEdgeL * uUEdgeL) + 0.25 * (metEdgeF * vEdgeF &
            &+ metUEdgeF * vUEdgeF + metEdgeB * vEdgeB + metUEdgeB * vUEdgeB) &
            &+ wEdgeU) / jacEdgeU
      end if

    case(4)

      ! Multiplication on u grid (inverted)

      jacEdgeU = 0.5 * (jac(i, j, k) + jac(i, j, k + 1))
      metEdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))
      metUEdgeR = 0.5 * (met(i, j, k + 1, 1, 3) + met(i + 1, j, k + 1, 1, 3))
      metEdgeL = 0.5 * (met(i, j, k, 1, 3) + met(i - 1, j, k, 1, 3))
      metUEdgeL = 0.5 * (met(i, j, k + 1, 1, 3) + met(i - 1, j, k + 1, 1, 3))
      metEdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3))
      metUEdgeF = 0.5 * (met(i, j, k + 1, 2, 3) + met(i, j + 1, k + 1, 2, 3))
      metEdgeB = 0.5 * (met(i, j, k, 2, 3) + met(i, j - 1, k, 2, 3))
      metUEdgeB = 0.5 * (met(i, j, k + 1, 2, 3) + met(i, j - 1, k + 1, 2, 3))
      if(wind == "car") then
        trafoTFC = jacEdgeU * (- 0.25 * (metEdgeR * uEdgeR + metUEdgeR &
            &* uUEdgeR + metEdgeL * uEdgeL + metUEdgeL * uUEdgeL) - 0.25 &
            &* (metEdgeF * vEdgeF + metUEdgeF * vUEdgeF + metEdgeB * vEdgeB &
            &+ metUEdgeB * vUEdgeB) + wEdgeU)
      else if(wind == "tfc") then
        trafoTFC = 0.25 * (metEdgeR * uEdgeR + metUEdgeR * uUEdgeR + metEdgeL &
            &* uEdgeL + metUEdgeL * uUEdgeL) + 0.25 * (metEdgeF * vEdgeF &
            &+ metUEdgeF * vUEdgeF + metEdgeB * vEdgeB + metUEdgeB * vUEdgeB) &
            &+ wEdgeU / jacEdgeU
      end if

    case(5)

      ! Multiplication on uw grid

      jacEdgeU = 0.5 * (jac(i, j, k) + jac(i, j, k + 1))
      metEdgeRU = 0.25 * (jac(i, j, k) * met(i, j, k, 1, 3) + jac(i + 1, j, k) &
          &* met(i + 1, j, k, 1, 3) + jac(i, j, k + 1) * met(i, j, k + 1, 1, &
          &3) + jac(i + 1, j, k + 1) * met(i + 1, j, k + 1, 1, 3))
      metEdgeLU = 0.25 * (jac(i, j, k) * met(i, j, k, 1, 3) + jac(i - 1, j, k) &
          &* met(i - 1, j, k, 1, 3) + jac(i, j, k + 1) * met(i, j, k + 1, 1, &
          &3) + jac(i - 1, j, k + 1) * met(i - 1, j, k + 1, 1, 3))
      metEdgeFU = 0.25 * (jac(i, j, k) * met(i, j, k, 2, 3) + jac(i, j + 1, k) &
          &* met(i, j + 1, k, 2, 3) + jac(i, j, k + 1) * met(i, j, k + 1, 2, &
          &3) + jac(i, j + 1, k + 1) * met(i, j + 1, k + 1, 2, 3))
      metEdgeBU = 0.25 * (jac(i, j, k) * met(i, j, k, 2, 3) + jac(i, j - 1, k) &
          &* met(i, j - 1, k, 2, 3) + jac(i, j, k + 1) * met(i, j, k + 1, 2, &
          &3) + jac(i, j - 1, k + 1) * met(i, j - 1, k + 1, 2, 3))
      uEdgeRU = 0.5 * (uEdgeR + uUEdgeR)
      uEdgeLU = 0.5 * (uEdgeL + uUEdgeL)
      vEdgeFU = 0.5 * (vEdgeF + vUEdgeF)
      vEdgeBU = 0.5 * (vEdgeB + vUEdgeB)
      if(wind == "car") then
        trafoTFC = - 0.5 * (metEdgeRU * uEdgeRU - metEdgeLU * uEdgeLU) - 0.5 &
            &* (metEdgeFU * vEdgeFU - metEdgeBU * vEdgeBU) + jacEdgeU * wEdgeU
      else if(wind == "tfc") then
        trafoTFC = (0.5 * (metEdgeRU * uEdgeRU - metEdgeLU * uEdgeLU) + 0.5 &
            &* (metEdgeFU * vEdgeFU - metEdgeBU * vEdgeBU) + wEdgeU) / jacEdgeU
      end if

    case(6)

      ! Multiplication on uw grid (inverted)

      jacEdgeU = 0.5 * (jac(i, j, k) + jac(i, j, k + 1))
      metEdgeRU = 0.25 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3) + met(i, &
          &j, k + 1, 1, 3) + met(i + 1, j, k + 1, 1, 3))
      metEdgeLU = 0.25 * (met(i, j, k, 1, 3) + met(i - 1, j, k, 1, 3) + met(i, &
          &j, k + 1, 1, 3) + met(i - 1, j, k + 1, 1, 3))
      metEdgeFU = 0.25 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3) + met(i, &
          &j, k + 1, 2, 3) + met(i, j + 1, k + 1, 2, 3))
      metEdgeBU = 0.25 * (met(i, j, k, 2, 3) + met(i, j - 1, k, 2, 3) + met(i, &
          &j, k + 1, 2, 3) + met(i, j - 1, k + 1, 2, 3))
      uEdgeRU = 0.5 * (uEdgeR + uUEdgeR)
      uEdgeLU = 0.5 * (uEdgeL + uUEdgeL)
      vEdgeFU = 0.5 * (vEdgeF + vUEdgeF)
      vEdgeBU = 0.5 * (vEdgeB + vUEdgeB)
      if(wind == "car") then
        trafoTFC = jacEdgeU * (- 0.5 * (metEdgeRU * uEdgeRU - metEdgeLU &
            &* uEdgeLU) - 0.5 * (metEdgeFU * vEdgeFU - metEdgeBU * vEdgeBU) &
            &+ wEdgeU)
      else if(wind == "tfc") then
        trafoTFC = 0.5 * (metEdgeRU * uEdgeRU - metEdgeLU * uEdgeLU) + 0.5 &
            &* (metEdgeFU * vEdgeFU - metEdgeBU * vEdgeBU) + wEdgeU / jacEdgeU
      end if

    end select

  end function trafoTFC

  function stressTensTFC(i, j, k, mu, nu, var)
    ! Cartesian stress tensor.

    type(var_type) :: var
    integer :: i, j, k, mu, nu

    real :: stressTensTFC
    real :: jacEdgeR, jacEdgeL, jacEdgeF, jacEdgeB, jacEdgeU, jacEdgeD
    real :: met13EdgeU, met13EdgeD, met23EdgeU, met23EdgeD
    real :: uF, uB, uU, uD, vR, vL, vU, vD, wR, wL, wF, wB

    jacEdgeR = 0.5 * (jac(i, j, k) + jac(i + 1, j, k))
    jacEdgeL = 0.5 * (jac(i, j, k) + jac(i - 1, j, k))
    jacEdgeF = 0.5 * (jac(i, j, k) + jac(i, j + 1, k))
    jacEdgeB = 0.5 * (jac(i, j, k) + jac(i, j - 1, k))
    jacEdgeU = 0.5 * (jac(i, j, k) + jac(i, j, k + 1))
    jacEdgeD = 0.5 * (jac(i, j, k) + jac(i, j, k - 1))
    met13EdgeU = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) + jac(i, j, k + 1) &
        &* met(i, j, k + 1, 1, 3))
    met13EdgeD = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) + jac(i, j, k - 1) &
        &* met(i, j, k - 1, 1, 3))
    met23EdgeU = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) + jac(i, j, k + 1) &
        &* met(i, j, k + 1, 2, 3))
    met23EdgeD = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) + jac(i, j, k - 1) &
        &* met(i, j, k - 1, 2, 3))
    uF = 0.5 * (var%u(i, j + 1, k) + var%u(i - 1, j + 1, k))
    uB = 0.5 * (var%u(i, j - 1, k) + var%u(i - 1, j - 1, k))
    uU = 0.5 * (var%u(i, j, k + 1) + var%u(i - 1, j, k + 1))
    uD = 0.5 * (var%u(i, j, k - 1) + var%u(i - 1, j, k - 1))
    vR = 0.5 * (var%v(i + 1, j, k) + var%v(i + 1, j - 1, k))
    vL = 0.5 * (var%v(i - 1, j, k) + var%v(i - 1, j - 1, k))
    vU = 0.5 * (var%v(i, j, k + 1) + var%v(i, j - 1, k + 1))
    vD = 0.5 * (var%v(i, j, k - 1) + var%v(i, j - 1, k - 1))
    wR = 0.5 * (vertWindTFC(i + 1, j, k, var) + vertWindTFC(i + 1, j, k - 1, &
        &var))
    wL = 0.5 * (vertWindTFC(i - 1, j, k, var) + vertWindTFC(i - 1, j, k - 1, &
        &var))
    wF = 0.5 * (vertWindTFC(i, j + 1, k, var) + vertWindTFC(i, j + 1, k - 1, &
        &var))
    wB = 0.5 * (vertWindTFC(i, j - 1, k, var) + vertWindTFC(i, j - 1, k - 1, &
        &var))

    if(mu == 1 .and. nu == 1) then
      stressTensTFC = (2.0 * (jacEdgeR * var%u(i, j, k) - jacEdgeL * var%u(i &
          &- 1, j, k)) / dx + (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) * uU &
          &- jac(i, j, k - 1) * met(i, j, k - 1, 1, 3) * uD) / dz - 2.0 / 3.0 &
          &* ((jacEdgeR * var%u(i, j, k) - jacEdgeL * var%u(i - 1, j, k)) / dx &
          &+ (jacEdgeF * var%v(i, j, k) - jacEdgeB * var%v(i, j - 1, k)) / dy &
          &+ (jacEdgeU * var%w(i, j, k) - jacEdgeD * var%w(i, j, k - 1)) &
          &/ dz)) / jac(i, j, k)
    else if((mu == 1 .and. nu == 2) .or. (mu == 2 .and. nu == 1)) then
      stressTensTFC = (0.5 * (jac(i, j + 1, k) * uF - jac(i, j - 1, k) * uB) &
          &/ dy + 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) * uU &
          &- jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) * uD) / dz + 0.5 &
          &* (jac(i + 1, j, k) * vR - jac(i - 1, j, k) * vL) / dx + 0.5 &
          &* (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) * vU - jac(i, j, k &
          &- 1) * met(i, j, k - 1, 1, 3) * vD) / dz) / jac(i, j, k)
    else if((mu == 1 .and. nu == 3) .or. (mu == 3 .and. nu == 1)) then
      stressTensTFC = (0.5 * (uU - uD) / dz + 0.5 * (jac(i + 1, j, k) * wR &
          &- jac(i - 1, j, k) * wL) / dx + (met13EdgeU * vertWindTFC(i, j, k, &
          &var) - met13EdgeD * vertWindTFC(i, j, k - 1, var)) / dz) / jac(i, &
          &j, k)
    else if(mu == 2 .and. nu == 2) then
      stressTensTFC = (2.0 * (jacEdgeF * var%v(i, j, k) - jacEdgeB * var%v(i, &
          &j - 1, k)) / dy + (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) * vU &
          &- jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) * vD) / dz - 2.0 / 3.0 &
          &* ((jacEdgeR * var%u(i, j, k) - jacEdgeL * var%u(i - 1, j, k)) / dx &
          &+ (jacEdgeF * var%v(i, j, k) - jacEdgeB * var%v(i, j - 1, k)) / dy &
          &+ (jacEdgeU * var%w(i, j, k) - jacEdgeD * var%w(i, j, k - 1)) &
          &/ dz)) / jac(i, j, k)
    else if((mu == 2 .and. nu == 3) .or. (mu == 3 .and. nu == 2)) then
      stressTensTFC = (0.5 * (vU - vD) / dz + 0.5 * (jac(i, j + 1, k) * wF &
          &- jac(i, j - 1, k) * wB) / dy + (met23EdgeU * vertWindTFC(i, j, k, &
          &var) - met23EdgeD * vertWindTFC(i, j, k - 1, var)) / dz) / jac(i, &
          &j, k)
    else if(mu == 3 .and. nu == 3) then
      stressTensTFC = (2.0 * (vertWindTFC(i, j, k, var) - vertWindTFC(i, j, k &
          &- 1, var)) / dz - 2.0 / 3.0 * ((jacEdgeR * var%u(i, j, k) &
          &- jacEdgeL * var%u(i - 1, j, k)) / dx + (jacEdgeF * var%v(i, j, k) &
          &- jacEdgeB * var%v(i, j - 1, k)) / dy + (jacEdgeU * var%w(i, j, k) &
          &- jacEdgeD * var%w(i, j, k - 1)) / dz)) / jac(i, j, k)
    end if

  end function stressTensTFC

end module atmosphere_module
