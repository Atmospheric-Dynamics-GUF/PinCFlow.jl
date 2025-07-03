module atmosphere_module

  use type_module
  use sizeof_module
  use mpi_module
  use mpi

  implicit none

  public ! all objects are known to programmes using this module

  !--------------
  !    public
  !--------------
  public :: init_atmosphere
  public :: terminate_atmosphere

  ! TFC
  public :: vertWindTFC, trafoTFC, stressTensTFC
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

  ! Jacobian and metric tensor
  real, dimension(:, :, :), allocatable :: jac
  real, dimension(:, :, :, :, :), allocatable :: met

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

    ! Allocate Jacobian, metric tensor and 3D background fields.
    if(topography) then
      if(.not. allocated(jac)) then
        allocate(jac((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
            &+ nbz)), stat = allocstat)
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate jac"
      end if

      if(.not. allocated(met)) then
        allocate(met((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
            &+ nbz), 1:3, 1:3), stat = allocstat)
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate met"
      end if

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

    if(compute_cloudcover) then
      dxsc = dx / NSCX
      dysc = dy / NSCY
      dzsc = dz / NSCZ

      SizeX2 = sizeX * NSCX
      SizeY2 = sizeY * NSCY
      SizeZ2 = sizeZ * NSCZ

      allocate(x2(SizeX2), stat = allocstat)
      if(allocstat /= 0) stop "init.f90: could not allocate x2"

      allocate(y2(sizeY2), stat = allocstat)
      if(allocstat /= 0) stop "init.f90: could not allocate y2"

      allocate(z2(sizeZ2), stat = allocstat)
      if(allocstat /= 0) stop "init.f90: could not allocate z2"

      k = 0
      do i = 1, sizeX
        do j = 1, NSCX
          k = k + 1
          x2(k) = lx(0) + real(k - 1) * dxsc + dxsc / 2.0
        end do
      end do

      k = 0
      do j = 1, sizeY
        do i = 1, NSCY
          k = k + 1
          y2(k) = ly(0) + real(k - 1) * dysc + dysc / 2.0
        end do
      end do

      i = 0
      do k = 1, sizeZ
        do j = 1, NSCZ
          i = i + 1
          z2(i) = lz(0) + real(i - 1) * dzsc + dzsc / 2.0
        end do
      end do
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

    case("pseudo_incompressible", "compressible")

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

    case default
      print *, "model = ", model
      stop "init_atmosphere: unknown case model."
    end select

    ! Set background fields in TFC.
    if(topography) then
      call set_background_tfc
    end if

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

    ! Deallocate Jacobian, metric tensor and 3D background fields.
    if(topography) then
      deallocate(jac, stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not dealloc jac"

      deallocate(met, stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not dealloc met"

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

  subroutine setup_topography

    real :: mountainHeight, mountainWidth, mountainWavenumber
    real :: x_center, y_center
    real :: kk, ll
    integer :: ix0, jy0
    integer :: ix, jy, kz
    integer :: iwm

    if(.not. topography) return

    if(lz(0) /= 0.0) stop "Error in setup_topography: lz(0) must be zero for &
        &TFC!"

    if(zBoundary == "periodic") then
      if(mountainHeight_dim /= 0.0) stop "Error in setup_topography: &
          &mountainHeight_dim must be zero for zBoundary = 'periodic'"
      if(stretch_exponent /= 1.0) stop "Error in setup_topography: &
          &stretch_exponent must be one for zBoundary = 'periodic'"
    end if

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
              if(abs(x(ix + ix0) - x_center) <= mountainWidth * width_factor) &
                  &then
                k_spectrum(ix, jy, 1) = mountainWavenumber
                topography_spectrum(ix, jy, 1) = 0.25 * mountainHeight * (1.0 &
                    &+ cos(mountainWavenumber / width_factor * (x(ix + ix0) &
                    &- x_center)))
              end if
              topography_surface(ix, jy) = 0.5 * mountainHeight

            case(7)
              ! 2D Gaussian envelope and even background
              k_spectrum(ix, jy, 1) = mountainWavenumber
              topography_spectrum(ix, jy, 1) = 0.5 * mountainHeight * exp(- &
                  &(x(ix + ix0) - x_center) ** 2.0 / (mountainWidth &
                  &* width_factor) ** 2.0)
              topography_surface(ix, jy) = 0.5 * mountainHeight

            case(9)
              ! 2D cosine envelope and cosine background
              if(abs(x(ix + ix0) - x_center) <= mountainWidth * width_factor) &
                  &then
                k_spectrum(ix, jy, 1) = mountainWavenumber
                topography_spectrum(ix, jy, 1) = 0.25 * mountainHeight * (1.0 &
                    &+ cos(mountainWavenumber / width_factor * (x(ix + ix0) &
                    &- x_center)))
                topography_surface(ix, jy) = 0.25 * mountainHeight * (1.0 &
                    &+ cos(mountainWavenumber / width_factor * (x(ix + ix0) &
                    &- x_center)))
              end if

            case(11)
              ! 2D Gaussian envelope and Gaussian background
              k_spectrum(ix, jy, 1) = mountainWavenumber
              topography_spectrum(ix, jy, 1) = 0.5 * mountainHeight * exp(- &
                  &(x(ix + ix0) - x_center) ** 2.0 / (mountainWidth &
                  &* width_factor) ** 2.0)
              topography_surface(ix, jy) = 0.5 * mountainHeight * exp(- (x(ix &
                  &+ ix0) - x_center) ** 2.0 / (mountainWidth * width_factor) &
                  &** 2.0)

            case(13)
              ! 3D WKB topography
              if((x(ix + ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                  &** 2.0 <= (mountainWidth * width_factor) ** 2.0) then
                do iwm = 0, spectral_modes - 1
                  k_spectrum(ix, jy, iwm + 1) = mountainWavenumber * cos(pi &
                      &/ spectral_modes * iwm)
                  l_spectrum(ix, jy, iwm + 1) = mountainWavenumber * sin(pi &
                      &/ spectral_modes * iwm)
                  topography_spectrum(ix, jy, iwm + 1) = 0.5 * mountainHeight &
                      &* (1.0 + cos(mountainWavenumber / width_factor &
                      &* sqrt((x(ix + ix0) - x_center) ** 2.0 + (y(jy + jy0) &
                      &- y_center) ** 2.0))) / spectral_modes / (height_factor &
                      &+ 1.0)
                end do
                topography_surface(ix, jy) = 0.5 * mountainHeight * (1.0 &
                    &+ cos(mountainWavenumber / width_factor * sqrt((x(ix &
                    &+ ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                    &** 2.0))) * height_factor / (height_factor + 1.0)
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
        reference_topography_spectrum = topography_spectrum
        topography_spectrum = 0.0
        reference_topography_surface = topography_surface
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
              if(abs(x(ix + ix0) - x_center) <= mountainWidth * width_factor) &
                  &then
                topography_surface(ix, jy) = 0.5 * mountainHeight * (1.0 + 0.5 &
                    &* (1.0 + cos(mountainWavenumber / width_factor * (x(ix &
                    &+ ix0) - x_center))) * cos(mountainWavenumber * (x(ix &
                    &+ ix0) - x_center)))
              else
                topography_surface(ix, jy) = 0.5 * mountainHeight
              end if

            case(6)
              ! 3D cosine envelope and even background
              if((x(ix + ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                  &** 2.0 <= (mountainWidth * width_factor) ** 2.0) then
                topography_surface(ix, jy) = 0.5 * mountainHeight * (1.0 + 0.5 &
                    &* (1.0 + cos(mountainWavenumber / width_factor &
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
                  &* width_factor) ** 2.0) * cos(mountainWavenumber * (x(ix &
                  &+ ix0) - x_center)))

            case(8)
              ! 3D Gaussian envelope and even background
              topography_surface(ix, jy) = 0.5 * mountainHeight * (1.0 + exp(- &
                  &((x(ix + ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                  &** 2.0) / (mountainWidth * width_factor) ** 2.0) &
                  &* cos(mountainWavenumber * sqrt((x(ix + ix0) - x_center) &
                  &** 2.0 + (y(jy + jy0) - y_center) ** 2.0)))

            case(9)
              ! 2D cosine envelope and cosine background
              if(abs(x(ix + ix0) - x_center) <= mountainWidth * width_factor) &
                  &then
                topography_surface(ix, jy) = 0.25 * mountainHeight * (1.0 &
                    &+ cos(mountainWavenumber / width_factor * (x(ix + ix0) &
                    &- x_center))) * (1.0 + cos(mountainWavenumber * (x(ix &
                    &+ ix0) - x_center)))
              end if

            case(10)
              ! 3D cosine envelope and cosine background
              if((x(ix + ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                  &** 2.0 <= (mountainWidth * width_factor) ** 2.0) then
                topography_surface(ix, jy) = 0.25 * mountainHeight * (1.0 &
                    &+ cos(mountainWavenumber / width_factor * sqrt((x(ix &
                    &+ ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                    &** 2.0))) * (1.0 + cos(mountainWavenumber * sqrt((x(ix &
                    &+ ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                    &** 2.0)))
              end if

            case(11)
              ! 2D Gaussian envelope and Gaussian background
              topography_surface(ix, jy) = 0.5 * mountainHeight * exp(- (x(ix &
                  &+ ix0) - x_center) ** 2.0 / (mountainWidth * width_factor) &
                  &** 2.0) * (1.0 + cos(mountainWavenumber * (x(ix + ix0) &
                  &- x_center)))

            case(12)
              ! 3D Gaussian envelope and Gaussian background
              topography_surface(ix, jy) = 0.5 * mountainHeight * exp(- ((x(ix &
                  &+ ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                  &** 2.0) / (mountainWidth * width_factor) ** 2.0) * (1.0 &
                  &+ cos(mountainWavenumber * sqrt((x(ix + ix0) - x_center) &
                  &** 2.0 + (y(jy + jy0) - y_center) ** 2.0)))

            case(13)
              ! 3D WKB topography
              if((x(ix + ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                  &** 2.0 <= (mountainWidth * width_factor) ** 2.0) then
                topography_surface(ix, jy) = 0.5 * mountainHeight * (1.0 &
                    &+ cos(mountainWavenumber / width_factor * sqrt((x(ix &
                    &+ ix0) - x_center) ** 2.0 + (y(jy + jy0) - y_center) &
                    &** 2.0))) * height_factor / (height_factor + 1.0)
                do iwm = 0, spectral_modes - 1
                  kk = mountainWavenumber * cos(pi / spectral_modes * iwm)
                  ll = mountainWavenumber * sin(pi / spectral_modes * iwm)
                  topography_surface(ix, jy) = topography_surface(ix, jy) &
                      &+ 0.5 * mountainHeight * (1.0 + cos(mountainWavenumber &
                      &/ width_factor * sqrt((x(ix + ix0) - x_center) ** 2.0 &
                      &+ (y(jy + jy0) - y_center) ** 2.0))) * cos(kk * (x(ix &
                      &+ ix0) - x_center) + ll * (y(jy + jy0) - y_center)) &
                      &/ spectral_modes / (height_factor + 1.0)
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
        reference_topography_surface = topography_surface
        topography_surface = 0.0
      end if
    end if

    ! Compute the stretched vertical grid.
    do kz = - nbz, nz + nbz
      zTildeS(kz) = map(z(kz) + 0.5 * dz)
    end do
    do kz = - nbz + 1, nz + nbz
      zS(kz) = 0.5 * (zTildeS(kz) + zTildeS(kz - 1))
    end do
    zS(- nbz) = zTildeS(- nbz) - 0.5 * (zTildeS(nbz + 1) - zTildeS(nbz))

    ! Compute the vertical layers, Jacobian and metric tensor.
    call compute_grid_tfc

  end subroutine setup_topography

  !---------------------------------------------------------------------------

  subroutine update_topography(var, time, dt)

    type(var_type), intent(inout) :: var
    real, intent(in) :: time
    real, intent(in) :: dt

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: &
        &zTFCOld, zTildeTFCOld
    type(var_type) :: varOld

    integer :: ix, jy, kz

    real :: zc, zu, zd, phi
    real :: zwc, zwu, zwd, phiw
    integer :: kzu, kzd
    integer :: kzwu, kzwd

    ! Return for constant topography.
    if(topographyTime <= 0.0) return

    ! Allocate varOld.
    call allocate_var_type(varOld)

    ! Save the old grid and prognostic variables.
    zTildeTFCOld = zTildeTFC
    zTFCOld = zTFC
    varOld = var

    ! Grow the topography spectrum.
    if(rayTracer .and. case_wkb == 3) then
      topography_spectrum = min(1.0, time / topographyTime * tRef) &
          &* reference_topography_spectrum
    end if

    ! Grow the topography surface.
    topography_surface = min(1.0, time / topographyTime * tRef) &
        &* reference_topography_surface

    ! Update the grid.
    call compute_grid_tfc

    ! Return if the grid hasn't changed.
    if(all(zTFCOld == zTFC) .and. all(zTildeTFCOld == zTildeTFC)) return

    ! Update the background fields.
    call set_background_tfc

    ! Interpolate to the new grid.
    do kz = 1, nz
      do jy = 1, ny
        do ix = 1, nx
          zc = zTFC(ix, jy, kz)
          zwc = zTildeTFC(ix, jy, kz)

          ! Determine adjacent levels.
          kzu = minloc(abs(zTFCOld(ix, jy, :) - zc), dim = 1) - nbz - 1
          if(zTFCOld(ix, jy, kzu) < zc) kzu = kzu + 1
          if(kzu < - nbz + 1) kzu = - nbz + 1
          if(kzu > nz + nbz) kzu = nz + nbz
          kzd = kzu - 1
          zu = zTFCOld(ix, jy, kzu)
          zd = zTFCOld(ix, jy, kzd)

          ! Set interpolation factor.
          if(zc > zu) then
            phi = 0.0
          else if(zc > zd) then
            phi = (zu - zc) / (zu - zd)
          else
            phi = 1.0
          end if

          ! Determine adjacent half-levels.
          kzwu = minloc(abs(zTildeTFCOld(ix, jy, :) - zwc), dim = 1) - nbz - 1
          if(zTildeTFCOld(ix, jy, kzwu) < zwc) kzwu = kzwu + 1
          if(kzwu < - nbz + 1) kzwu = - nbz + 1
          if(kzwu > nz + nbz) kzwu = nz + nbz
          kzwd = kzwu - 1
          zwu = zTildeTFCOld(ix, jy, kzwu)
          zwd = zTildeTFCOld(ix, jy, kzwd)

          ! Set interpolation factor for half-levels.
          if(zwc > zwu) then
            phiw = 0.0
          else if(zwc > zwd) then
            phiw = (zwu - zwc) / (zwu - zwd)
          else
            phiw = 1.0
          end if

          ! Perform interpolation (not all of these are technically needed).
          var%rho(ix, jy, kz) = phi * varOld%rho(ix, jy, kzd) + (1.0 - phi) &
              &* varOld%rho(ix, jy, kzu)
          var%u(ix, jy, kz) = phi * varOld%u(ix, jy, kzd) + (1.0 - phi) &
              &* varOld%u(ix, jy, kzu)
          var%v(ix, jy, kz) = phi * varOld%v(ix, jy, kzd) + (1.0 - phi) &
              &* varOld%v(ix, jy, kzu)
          var%w(ix, jy, kz) = phiw * varOld%w(ix, jy, kzwd) + (1.0 - phiw) &
              &* varOld%w(ix, jy, kzwu)
          var%pi(ix, jy, kz) = phi * varOld%pi(ix, jy, kzd) + (1.0 - phi) &
              &* varOld%pi(ix, jy, kzu)
          var%rhop(ix, jy, kz) = phi * varOld%rhop(ix, jy, kzd) + (1.0 - phi) &
              &* varOld%rhop(ix, jy, kzu)
          if(turbScheme) then
            var%DSC(ix, jy, kz) = phi * varOld%DSC(ix, jy, kzd) + (1.0 - phi) &
                &* varOld%DSC(ix, jy, kzu)
          end if
          if(rayTracer) then
            var%GWH(ix, jy, kz) = phi * varOld%GWH(ix, jy, kzd) + (1.0 - phi) &
                &* varOld%GWH(ix, jy, kzu)
          end if
          if(model == "compressible") then
            var%P(ix, jy, kz) = phi * varOld%P(ix, jy, kzd) + (1.0 - phi) &
                &* varOld%P(ix, jy, kzu)
          end if
          if(include_tracer) then
            var%chi(ix, jy, kz) = phi * varOld%chi(ix, jy, kzd) + (1.0 - phi) &
                &* varOld%chi(ix, jy, kzu)
          end if
          if(include_ice) then
            var%ICE(ix, jy, kz, :) = phi * varOld%ICE(ix, jy, kzd, :) + (1.0 &
                &- phi) * varOld%ICE(ix, jy, kzu, :)
          end if
          if(include_testoutput) then
            var%OPT(ix, jy, kz, :) = phi * varOld%OPT(ix, jy, kzd, :) + (1.0 &
                &- phi) * varOld%OPT(ix, jy, kzu, :)
          end if
        end do
      end do
    end do

  end subroutine update_topography

  !---------------------------------------------------------------------------

  subroutine compute_grid_tfc

    ! Compute the vertical layers, Jacobian and metric tensor.

    integer :: ix, jy, kz

    ! Compute the vertical layers.
    do kz = - nbz, nz + nbz
      zTildeTFC(:, :, kz) = (lz(1) - topography_surface) / lz(1) * zTildeS(kz) &
          &+ topography_surface
      zTFC(:, :, kz) = (lz(1) - topography_surface) / lz(1) * zS(kz) &
          &+ topography_surface
    end do

    ! Compute the Jacobian.
    do kz = - nbz + 1, nz + nbz
      jac(:, :, kz) = (lz(1) - topography_surface) / lz(1) * (zTildeS(kz) &
          &- zTildeS(kz - 1)) / dz
    end do
    jac(:, :, - nbz) = jac(:, :, nbz + 1)

    ! Compute the metric tensor.
    met(:, :, :, 1, 1) = 1.0
    met(:, :, :, 1, 2) = 0.0
    do kz = - nbz + 1, nz + nbz
      do jy = 1, ny
        do ix = 1, nx
          met(ix, jy, kz, 1, 3) = (topography_surface(ix + 1, jy) &
              &- topography_surface(ix - 1, jy)) / (2.0 * dx) * (zS(kz) &
              &- lz(1)) / (lz(1) - topography_surface(ix, jy)) * dz &
              &/ (zTildeS(kz) - zTildeS(kz - 1))
        end do
      end do
      call setHalosOfField2D(met(:, :, kz, 1, 3))
    end do
    met(:, :, - nbz, 1, 3) = met(:, :, nbz + 1, 1, 3) * (zS(- nbz) - lz(1)) &
        &/ (zS(nbz + 1) - lz(1))
    met(:, :, :, 2, 1) = 0.0
    met(:, :, :, 2, 2) = 1.0
    do kz = - nbz + 1, nz + nbz
      do jy = 1, ny
        do ix = 1, nx
          met(ix, jy, kz, 2, 3) = (topography_surface(ix, jy + 1) &
              &- topography_surface(ix, jy - 1)) / (2.0 * dy) * (zS(kz) &
              &- lz(1)) / (lz(1) - topography_surface(ix, jy)) * dz &
              &/ (zTildeS(kz) - zTildeS(kz - 1))
        end do
      end do
      call setHalosOfField2D(met(:, :, kz, 2, 3))
    end do
    met(:, :, - nbz, 2, 3) = met(:, :, nbz + 1, 2, 3) * (zS(- nbz) - lz(1)) &
        &/ (zS(nbz + 1) - lz(1))
    met(:, :, :, 3, 1) = met(:, :, :, 1, 3)
    met(:, :, :, 3, 2) = met(:, :, :, 2, 3)
    do kz = - nbz + 1, nz + nbz
      do jy = 1, ny
        do ix = 1, nx
          met(ix, jy, kz, 3, 3) = ((lz(1) / (lz(1) - topography_surface(ix, &
              &jy))) ** 2.0 + ((zS(kz) - lz(1)) / (lz(1) &
              &- topography_surface(ix, jy))) ** 2.0 &
              &* (((topography_surface(ix + 1, jy) - topography_surface(ix &
              &- 1, jy)) / (2.0 * dx)) ** 2.0 + ((topography_surface(ix, jy &
              &+ 1) - topography_surface(ix, jy - 1)) / (2.0 * dy)) ** 2.0)) &
              &* (dz / (zTildeS(kz) - zTildeS(kz - 1))) ** 2.0
        end do
      end do
      call setHalosOfField2D(met(:, :, kz, 3, 3))
    end do
    do jy = 1, ny
      do ix = 1, nx
        met(ix, jy, - nbz, 3, 3) = ((lz(1) / (lz(1) - topography_surface(ix, &
            &jy))) ** 2.0 + ((zS(- nbz) - lz(1)) / (lz(1) &
            &- topography_surface(ix, jy))) ** 2.0 * (((topography_surface(ix &
            &+ 1, jy) - topography_surface(ix - 1, jy)) / (2.0 * dx)) ** 2.0 &
            &+ ((topography_surface(ix, jy + 1) - topography_surface(ix, jy &
            &- 1)) / (2.0 * dy)) ** 2.0)) * (dz / (zTildeS(nbz + 1) &
            &- zTildeS(nbz))) ** 2.0
      end do
    end do
    call setHalosOfField2D(met(:, :, - nbz, 3, 3))

  end subroutine compute_grid_tfc

  !---------------------------------------------------------------------------

  subroutine set_background_tfc

    ! Set all background fields in TFC.

    real :: z_tr
    real :: theta_tr
    real :: press_tr
    real :: T_tr

    real :: power

    real :: pow_t, pow_s, p_t_b
    real :: T_bar, T_c_b, p_bar

    real :: pistar, thetastar

    integer :: i, j, k

    select case(model)
    case("pseudo_incompressible", "compressible")

      select case(background)

      case("smoothV")

        pStratTFC = p00
        thetaStratTFC = theta00
        rhoStratTFC = rho00

      case("realistic")

        z_tr = z_tr_dim / lRef
        theta_tr = theta_tr_dim / thetaRef
        press_tr = p0 * (1.0 - kappa * sig / theta_tr * z_tr) ** (1 / kappa)
        T_tr = theta_tr * (press_tr / p0) ** kappa

        do i = - nbx, nx + nbx
          do j = - nby, ny + nby
            do k = - 1, nz + 2
              if(zTFC(i, j, k) <= z_tr) then
                ! Isentropic troposphere.
                power = 1.0 / (gamma - 1.0)
                term = kappa * sig / theta_tr * zTFC(i, j, k)
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
                    &* exp(- sig / gamma / T_tr * (zTFC(i, j, k) - z_tr))
                ! Define thetaStratTFC.
                thetaStratTFC(i, j, k) = theta_tr * exp(kappa * sig / T_tr &
                    &* (zTFC(i, j, k) - z_tr))
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
              pStratTFC(i, j, k) = p0 * exp(- sig * zTFC(i, j, k) / gamma / T0)
              ! Define thetaStratTFC.
              thetaStratTFC(i, j, k) = T0 * exp(kappa * sig / T0 * zTFC(i, j, &
                  &k))
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
              term = kappa * sig / theta0 * zTFC(i, j, k)
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
              term = exp(- Fr2 * N2 * zTFC(i, j, k))
              ! Define pStratTFC.
              pStratTFC(i, j, k) = p0 * (1.0 + FrInv2 * kappa * sig / N2 &
                  &/ theta0 * (term - 1.0)) ** power
              ! Define thetaStratTFC.
              thetaStratTFC(i, j, k) = theta0 * exp(Fr ** 2.0 * N2 * zTFC(i, &
                  &j, k))
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
              zk = zTFC(i, j, k) * lRef

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

      case("HeldSuarez")

        tp_strato = tp_strato_dim / thetaRef
        tp_srf_trp = tp_srf_trp_dim / thetaRef
        tpdiffhor_tropo = tpdiffhor_tropo_dim / thetaRef
        ptdiffvert_tropo = ptdiffvert_tropo_dim / thetaRef

        do i = - nbx, nx + nbx
          do j = - nby, ny + nby
            ! Define Exner pressure and 3D background fields at
            ! the surface.
            do k = 0, 1
              T_bar = max(tp_strato, tp_srf_trp - 0.5 * tpdiffhor_tropo)
              ! Define Exner pressure.
              piStratTFC(i, j, k) = 1.0 - zTFC(i, j, k) * kappa / T_bar
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
              rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, j, k)
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
              piStratTFC(i, j, k) = piStratTFC(i, j, k - 1) - dz * jac(i, j, &
                  &k) * jac(i, j, k - 1) / (jac(i, j, k) + jac(i, j, k - 1)) &
                  &* (kappa / thetaStar + kappa / thetaStratTFC(i, j, k - 1))
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
              rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, j, k)
            end do
            ! Adjust at the lower boundary.
            piStratTFC(i, j, - 1) = piStratTFC(i, j, 0)
            pStratTFC(i, j, - 1) = pStratTFC(i, j, 0)
            thetaStratTFC(i, j, - 1) = thetaStratTFC(i, j, 0)
            rhoStratTFC(i, j, - 1) = rhoStratTFC(i, j, 0)
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

    case("Boussinesq")

      pStratTFC = p00
      thetaStratTFC = theta00
      rhoStratTFC = rho00
      bvsStratTFC = N2

    end select

  end subroutine set_background_tfc

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
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            var%w(i, j, k) = var%w(i, j, k) * jac(i, j, k) * jac(i, j, k + 1) &
                &* (var%P(i, j, k) + var%P(i, j, k + 1)) / (jac(i, j, k) &
                &+ jac(i, j, k + 1))

            var%v(i, j, k) = var%v(i, j, k) * 0.5 * (jac(i, j, k) * var%P(i, &
                &j, k) + jac(i, j + 1, k) * var%P(i, j + 1, k))

            var%u(i, j, k) = var%u(i, j, k) * 0.5 * (jac(i, j, k) * var%P(i, &
                &j, k) + jac(i + 1, j, k) * var%P(i + 1, j, k))
          end do
        end do
      end do

    case("backward")
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            var%w(i, j, k) = var%w(i, j, k) / jac(i, j, k) / jac(i, j, k + 1) &
                &/ (var%P(i, j, k) + var%P(i, j, k + 1)) * (jac(i, j, k) &
                &+ jac(i, j, k + 1))

            var%v(i, j, k) = var%v(i, j, k) / (0.5 * (jac(i, j, k) * var%P(i, &
                &j, k) + jac(i, j + 1, k) * var%P(i, j + 1, k)))

            var%u(i, j, k) = var%u(i, j, k) / (0.5 * (jac(i, j, k) * var%P(i, &
                &j, k) + jac(i + 1, j, k) * var%P(i + 1, j, k)))
          end do
        end do
      end do

    case default
      stop "add_JP_to_u unknown option."
    end select

    call setHalosOfField(var%u)
    call setHalosOfField(var%v)
    call setHalosOfField(var%w)

    select case(zBoundary)
    case("periodic")
      var%w(:, :, 0) = var%w(:, :, nz)
      do k = 1, nbz
        var%u(:, :, nz + k) = var%u(:, :, k)
        var%u(:, :, - k + 1) = var%u(:, :, nz - k + 1)
        var%v(:, :, nz + k) = var%v(:, :, k)
        var%v(:, :, - k + 1) = var%v(:, :, nz - k + 1)
        var%w(:, :, nz + k) = var%w(:, :, k)
        var%w(:, :, - k) = var%w(:, :, nz - k)
      enddo
    case("solid_wall")
      var%w(:, :, 0) = 0.0
      var%w(:, :, nz) = 0.0
      do k = 1, nbz
        var%u(:, :, - k + 1) = var%u(:, :, k)
        var%u(:, :, nz + k) = var%u(:, :, nz - k + 1)
        var%v(:, :, - k + 1) = var%v(:, :, k)
        var%v(:, :, nz + k) = var%v(:, :, nz - k + 1)
        var%w(:, :, - k) = - var%w(:, :, k)
        var%w(:, :, nz + k) = - var%w(:, :, nz - k)
      end do
    case default
      stop "Error in add_JP_to_u: Unknown zBoundary!"
    end select

  end subroutine

  !---------------------------------------------------------------------------

  function map(level)
    ! Vertical grid stretching.

    real :: map
    real :: level

    if(level < 0) then
      map = - lz(1) * (- level / lz(1)) ** stretch_exponent
    else if(level > lz(1)) then
      map = 2 * lz(1) - lz(1) * ((2 * lz(1) - level) / lz(1)) &
          &** stretch_exponent
    else
      map = lz(1) * (level / lz(1)) ** stretch_exponent
    end if

  end function map

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

      jacEdgeU = 2.0 * jac(i, j, k) * jac(i, j, k + 1) / (jac(i, j, k) &
          &+ jac(i, j, k + 1))
      uC = 0.5 * (uEdgeR + uEdgeL)
      uU = 0.5 * (uUEdgeR + uUEdgeL)
      vC = 0.5 * (vEdgeF + vEdgeB)
      vU = 0.5 * (vUEdgeF + vUEdgeB)
      if(wind == "car") then
        trafoTFC = - jac(i, j, k) * jac(i, j, k + 1) * (met(i, j, k, 1, 3) &
            &* uC + met(i, j, k + 1, 1, 3) * uU + met(i, j, k, 2, 3) * vC &
            &+ met(i, j, k + 1, 2, 3) * vU) / (jac(i, j, k) + jac(i, j, k &
            &+ 1)) + jacEdgeU * wEdgeU
      else if(wind == "tfc") then
        trafoTFC = (jac(i, j, k) * jac(i, j, k + 1) * (met(i, j, k, 1, 3) * uC &
            &+ met(i, j, k + 1, 1, 3) * uU + met(i, j, k, 2, 3) * vC + met(i, &
            &j, k + 1, 2, 3) * vU) / (jac(i, j, k) + jac(i, j, k + 1)) &
            &+ wEdgeU) / jacEdgeU
      end if

    case(2)

      ! Multiplication on rho grid (inverted)

      jacEdgeU = 2.0 * jac(i, j, k) * jac(i, j, k + 1) / (jac(i, j, k) &
          &+ jac(i, j, k + 1))
      uC = 0.5 * (uEdgeR + uEdgeL)
      uU = 0.5 * (uUEdgeR + uUEdgeL)
      vC = 0.5 * (vEdgeF + vEdgeB)
      vU = 0.5 * (vUEdgeF + vUEdgeB)
      if(wind == "car") then
        trafoTFC = jacEdgeU * (- (jac(i, j, k + 1) * (met(i, j, k, 1, 3) * uC &
            &+ met(i, j, k, 2, 3) * vC) + jac(i, j, k) * (met(i, j, k + 1, 1, &
            &3) * uU + met(i, j, k + 1, 2, 3) * vU)) / (jac(i, j, k) + jac(i, &
            &j, k + 1)) + wEdgeU)
      else if(wind == "tfc") then
        trafoTFC = (jac(i, j, k + 1) * (met(i, j, k, 1, 3) * uC + met(i, j, k, &
            &2, 3) * vC) + jac(i, j, k) * (met(i, j, k + 1, 1, 3) * uU &
            &+ met(i, j, k + 1, 2, 3) * vU)) / (jac(i, j, k) + jac(i, j, k &
            &+ 1)) + wEdgeU / jacEdgeU
      end if

    case(3)

      ! Multiplication on u grid

      jacEdgeU = 2.0 * jac(i, j, k) * jac(i, j, k + 1) / (jac(i, j, k) &
          &+ jac(i, j, k + 1))
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
        trafoTFC = - 0.5 * (jac(i, j, k + 1) * (metEdgeR * uEdgeR + metEdgeL &
            &* uEdgeL + metEdgeF * vEdgeF + metEdgeB * vEdgeB) + jac(i, j, k) &
            &* (metUEdgeR * uUEdgeR + metUEdgeL * uUEdgeL + metUEdgeF &
            &* vUEdgeF + metUEdgeB * vUEdgeB)) / (jac(i, j, k) + jac(i, j, k &
            &+ 1)) + jacEdgeU * wEdgeU
      else if(wind == "tfc") then
        trafoTFC = (0.5 * (jac(i, j, k + 1) * (metEdgeR * uEdgeR + metEdgeL &
            &* uEdgeL + metEdgeF * vEdgeF + metEdgeB * vEdgeB) + jac(i, j, k) &
            &* (metUEdgeR * uUEdgeR + metUEdgeL * uUEdgeL + metUEdgeF &
            &* vUEdgeF + metUEdgeB * vUEdgeB)) / (jac(i, j, k) + jac(i, j, k &
            &+ 1)) + wEdgeU) / jacEdgeU
      end if

    case(4)

      ! Multiplication on u grid (inverted)

      jacEdgeU = 2.0 * jac(i, j, k) * jac(i, j, k + 1) / (jac(i, j, k) &
          &+ jac(i, j, k + 1))
      metEdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))
      metUEdgeR = 0.5 * (met(i, j, k + 1, 1, 3) + met(i + 1, j, k + 1, 1, 3))
      metEdgeL = 0.5 * (met(i, j, k, 1, 3) + met(i - 1, j, k, 1, 3))
      metUEdgeL = 0.5 * (met(i, j, k + 1, 1, 3) + met(i - 1, j, k + 1, 1, 3))
      metEdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3))
      metUEdgeF = 0.5 * (met(i, j, k + 1, 2, 3) + met(i, j + 1, k + 1, 2, 3))
      metEdgeB = 0.5 * (met(i, j, k, 2, 3) + met(i, j - 1, k, 2, 3))
      metUEdgeB = 0.5 * (met(i, j, k + 1, 2, 3) + met(i, j - 1, k + 1, 2, 3))
      if(wind == "car") then
        trafoTFC = jacEdgeU * (- 0.5 * (jac(i, j, k + 1) * (metEdgeR * uEdgeR &
            &+ metEdgeL * uEdgeL + metEdgeF * vEdgeF + metEdgeB * vEdgeB) &
            &+ jac(i, j, k) * (metUEdgeR * uUEdgeR + metUEdgeL * uUEdgeL &
            &+ metUEdgeF * vUEdgeF + metUEdgeB * vUEdgeB)) / (jac(i, j, k) &
            &+ jac(i, j, k + 1)) + wEdgeU)
      else if(wind == "tfc") then
        trafoTFC = 0.5 * (jac(i, j, k + 1) * (metEdgeR * uEdgeR + metEdgeL &
            &* uEdgeL + metEdgeF * vEdgeF + metEdgeB * vEdgeB) + jac(i, j, k) &
            &* (metUEdgeR * uUEdgeR + metUEdgeL * uUEdgeL + metUEdgeF &
            &* vUEdgeF + metUEdgeB * vUEdgeB)) / (jac(i, j, k) + jac(i, j, k &
            &+ 1)) + wEdgeU / jacEdgeU
      end if

    case(5)

      ! Multiplication on uw grid

      jacEdgeU = 2.0 * jac(i, j, k) * jac(i, j, k + 1) / (jac(i, j, k) &
          &+ jac(i, j, k + 1))
      metEdgeRU = 0.5 * ((jac(i, j, k + 1) + jac(i + 1, j, k + 1)) * (jac(i, &
          &j, k) * met(i, j, k, 1, 3) + jac(i + 1, j, k) * met(i + 1, j, k, 1, &
          &3)) + (jac(i, j, k) + jac(i + 1, j, k)) * (jac(i, j, k + 1) &
          &* met(i, j, k + 1, 1, 3) + jac(i + 1, j, k + 1) * met(i + 1, j, k &
          &+ 1, 1, 3))) / (jac(i, j, k) + jac(i + 1, j, k) + jac(i, j, k + 1) &
          &+ jac(i + 1, j, k + 1))
      metEdgeLU = 0.5 * ((jac(i, j, k + 1) + jac(i - 1, j, k + 1)) * (jac(i, &
          &j, k) * met(i, j, k, 1, 3) + jac(i - 1, j, k) * met(i - 1, j, k, 1, &
          &3)) + (jac(i, j, k) + jac(i - 1, j, k)) * (jac(i, j, k + 1) &
          &* met(i, j, k + 1, 1, 3) + jac(i - 1, j, k + 1) * met(i - 1, j, k &
          &+ 1, 1, 3))) / (jac(i, j, k) + jac(i - 1, j, k) + jac(i, j, k + 1) &
          &+ jac(i - 1, j, k + 1))
      metEdgeFU = 0.5 * ((jac(i, j, k + 1) + jac(i, j + 1, k + 1)) * (jac(i, &
          &j, k) * met(i, j, k, 2, 3) + jac(i, j + 1, k) * met(i, j + 1, k, 2, &
          &3)) + (jac(i, j, k) + jac(i, j + 1, k)) * (jac(i, j, k + 1) &
          &* met(i, j, k + 1, 2, 3) + jac(i, j + 1, k + 1) * met(i, j + 1, k &
          &+ 1, 2, 3))) / (jac(i, j, k) + jac(i, j + 1, k) + jac(i, j, k + 1) &
          &+ jac(i, j + 1, k + 1))
      metEdgeBU = 0.5 * ((jac(i, j, k + 1) + jac(i, j - 1, k + 1)) * (jac(i, &
          &j, k) * met(i, j, k, 2, 3) + jac(i, j - 1, k) * met(i, j - 1, k, 2, &
          &3)) + (jac(i, j, k) + jac(i, j - 1, k)) * (jac(i, j, k + 1) &
          &* met(i, j, k + 1, 2, 3) + jac(i, j - 1, k + 1) * met(i, j - 1, k &
          &+ 1, 2, 3))) / (jac(i, j, k) + jac(i, j - 1, k) + jac(i, j, k + 1) &
          &+ jac(i, j - 1, k + 1))
      uEdgeRU = ((jac(i, j, k + 1) + jac(i + 1, j, k + 1)) * uEdgeR + (jac(i, &
          &j, k) + jac(i + 1, j, k)) * uUEdgeR) / (jac(i, j, k) + jac(i + 1, &
          &j, k) + jac(i, j, k + 1) + jac(i + 1, j, k + 1))
      uEdgeLU = ((jac(i, j, k + 1) + jac(i - 1, j, k + 1)) * uEdgeL + (jac(i, &
          &j, k) + jac(i - 1, j, k)) * uUEdgeL) / (jac(i, j, k) + jac(i - 1, &
          &j, k) + jac(i, j, k + 1) + jac(i - 1, j, k + 1))
      vEdgeFU = ((jac(i, j, k + 1) + jac(i, j + 1, k + 1)) * vEdgeF + (jac(i, &
          &j, k) + jac(i, j + 1, k)) * vUEdgeF) / (jac(i, j, k) + jac(i, j &
          &+ 1, k) + jac(i, j, k + 1) + jac(i, j + 1, k + 1))
      vEdgeBU = ((jac(i, j, k + 1) + jac(i, j - 1, k + 1)) * vEdgeB + (jac(i, &
          &j, k) + jac(i, j - 1, k)) * vUEdgeB) / (jac(i, j, k) + jac(i, j &
          &- 1, k) + jac(i, j, k + 1) + jac(i, j - 1, k + 1))
      if(wind == "car") then
        trafoTFC = - 0.5 * (metEdgeRU * uEdgeRU + metEdgeLU * uEdgeLU) - 0.5 &
            &* (metEdgeFU * vEdgeFU + metEdgeBU * vEdgeBU) + jacEdgeU * wEdgeU
      else if(wind == "tfc") then
        trafoTFC = (0.5 * (metEdgeRU * uEdgeRU + metEdgeLU * uEdgeLU) + 0.5 &
            &* (metEdgeFU * vEdgeFU + metEdgeBU * vEdgeBU) + wEdgeU) / jacEdgeU
      end if

    case(6)

      ! Multiplication on uw grid (inverted)

      jacEdgeU = 2.0 * jac(i, j, k) * jac(i, j, k + 1) / (jac(i, j, k) &
          &+ jac(i, j, k + 1))
      metEdgeRU = 0.5 * ((jac(i, j, k + 1) + jac(i + 1, j, k + 1)) * (met(i, &
          &j, k, 1, 3) + met(i + 1, j, k, 1, 3)) + (jac(i, j, k) + jac(i + 1, &
          &j, k)) * (met(i, j, k + 1, 1, 3) + met(i + 1, j, k + 1, 1, 3))) &
          &/ (jac(i, j, k) + jac(i + 1, j, k) + jac(i, j, k + 1) + jac(i + 1, &
          &j, k + 1))
      metEdgeLU = 0.5 * ((jac(i, j, k + 1) + jac(i - 1, j, k + 1)) * (met(i, &
          &j, k, 1, 3) + met(i - 1, j, k, 1, 3)) + (jac(i, j, k) + jac(i - 1, &
          &j, k)) * (met(i, j, k + 1, 1, 3) + met(i - 1, j, k + 1, 1, 3))) &
          &/ (jac(i, j, k) + jac(i - 1, j, k) + jac(i, j, k + 1) + jac(i - 1, &
          &j, k + 1))
      metEdgeFU = 0.5 * ((jac(i, j, k + 1) + jac(i, j + 1, k + 1)) * (met(i, &
          &j, k, 2, 3) + met(i, j + 1, k, 2, 3)) + (jac(i, j, k) + jac(i, j &
          &+ 1, k)) * (met(i, j, k + 1, 2, 3) + met(i, j + 1, k + 1, 2, 3))) &
          &/ (jac(i, j, k) + jac(i, j + 1, k) + jac(i, j, k + 1) + jac(i, j &
          &+ 1, k + 1))
      metEdgeBU = 0.5 * ((jac(i, j, k + 1) + jac(i, j - 1, k + 1)) * (met(i, &
          &j, k, 2, 3) + met(i, j - 1, k, 2, 3)) + (jac(i, j, k) + jac(i, j &
          &- 1, k)) * (met(i, j, k + 1, 2, 3) + met(i, j - 1, k + 1, 2, 3))) &
          &/ (jac(i, j, k) + jac(i, j - 1, k) + jac(i, j, k + 1) + jac(i, j &
          &- 1, k + 1))
      uEdgeRU = ((jac(i, j, k + 1) + jac(i + 1, j, k + 1)) * uEdgeR + (jac(i, &
          &j, k) + jac(i + 1, j, k)) * uUEdgeR) / (jac(i, j, k) + jac(i + 1, &
          &j, k) + jac(i, j, k + 1) + jac(i + 1, j, k + 1))
      uEdgeLU = ((jac(i, j, k + 1) + jac(i - 1, j, k + 1)) * uEdgeL + (jac(i, &
          &j, k) + jac(i - 1, j, k)) * uUEdgeL) / (jac(i, j, k) + jac(i - 1, &
          &j, k) + jac(i, j, k + 1) + jac(i - 1, j, k + 1))
      vEdgeFU = ((jac(i, j, k + 1) + jac(i, j + 1, k + 1)) * vEdgeF + (jac(i, &
          &j, k) + jac(i, j + 1, k)) * vUEdgeF) / (jac(i, j, k) + jac(i, j &
          &+ 1, k) + jac(i, j, k + 1) + jac(i, j + 1, k + 1))
      vEdgeBU = ((jac(i, j, k + 1) + jac(i, j - 1, k + 1)) * vEdgeB + (jac(i, &
          &j, k) + jac(i, j - 1, k)) * vUEdgeB) / (jac(i, j, k) + jac(i, j &
          &- 1, k) + jac(i, j, k + 1) + jac(i, j - 1, k + 1))
      if(wind == "car") then
        trafoTFC = jacEdgeU * (- 0.5 * (metEdgeRU * uEdgeRU + metEdgeLU &
            &* uEdgeLU) - 0.5 * (metEdgeFU * vEdgeFU + metEdgeBU * vEdgeBU) &
            &+ wEdgeU)
      else if(wind == "tfc") then
        trafoTFC = 0.5 * (metEdgeRU * uEdgeRU + metEdgeLU * uEdgeLU) + 0.5 &
            &* (metEdgeFU * vEdgeFU + metEdgeBU * vEdgeBU) + wEdgeU / jacEdgeU
      end if

    end select

  end function trafoTFC

  function stressTensTFC(i, j, k, mu, nu, var)
    ! Cartesian stress tensor.

    type(var_type) :: var
    integer :: i, j, k, mu, nu

    real :: stressTensTFC
    real :: jacEdgeR, jacEdgeL, jacEdgeF, jacEdgeB, jacEdgeU, jacEdgeD
    real :: uF, uB, uU, uD, vR, vL, vU, vD, wR, wL, wF, wB

    jacEdgeR = 0.5 * (jac(i, j, k) + jac(i + 1, j, k))
    jacEdgeL = 0.5 * (jac(i, j, k) + jac(i - 1, j, k))
    jacEdgeF = 0.5 * (jac(i, j, k) + jac(i, j + 1, k))
    jacEdgeB = 0.5 * (jac(i, j, k) + jac(i, j - 1, k))
    jacEdgeU = 2.0 * jac(i, j, k) * jac(i, j, k + 1) / (jac(i, j, k) + jac(i, &
        &j, k + 1))
    jacEdgeD = 2.0 * jac(i, j, k) * jac(i, j, k - 1) / (jac(i, j, k) + jac(i, &
        &j, k - 1))
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
      stressTensTFC = 2.0 * (var%u(i, j, k) - var%u(i - 1, j, k)) / dx &
          &+ met(i, j, k, 1, 3) * (uU - uD) / dz - 2.0 / 3.0 * ((jacEdgeR &
          &* var%u(i, j, k) - jacEdgeL * var%u(i - 1, j, k)) / dx + (jacEdgeF &
          &* var%v(i, j, k) - jacEdgeB * var%v(i, j - 1, k)) / dy + (jacEdgeU &
          &* var%w(i, j, k) - jacEdgeD * var%w(i, j, k - 1)) / dz) / jac(i, j, &
          &k)
    else if((mu == 1 .and. nu == 2) .or. (mu == 2 .and. nu == 1)) then
      stressTensTFC = 0.5 * (uF - uB) / dy + 0.5 * met(i, j, k, 2, 3) * (uU &
          &- uD) / dz + 0.5 * (vR - vL) / dx + 0.5 * met(i, j, k, 1, 3) * (vU &
          &- vD) / dz
    else if((mu == 1 .and. nu == 3) .or. (mu == 3 .and. nu == 1)) then
      stressTensTFC = 0.5 * (uU - uD) / dz / jac(i, j, k) + 0.5 * (wR - wL) &
          &/ dx + met(i, j, k, 1, 3) * (vertWindTFC(i, j, k, var) &
          &- vertWindTFC(i, j, k - 1, var)) / dz
    else if(mu == 2 .and. nu == 2) then
      stressTensTFC = 2.0 * (var%v(i, j, k) - var%v(i, j - 1, k)) / dy &
          &+ met(i, j, k, 2, 3) * (vU - vD) / dz - 2.0 / 3.0 * ((jacEdgeR &
          &* var%u(i, j, k) - jacEdgeL * var%u(i - 1, j, k)) / dx + (jacEdgeF &
          &* var%v(i, j, k) - jacEdgeB * var%v(i, j - 1, k)) / dy + (jacEdgeU &
          &* var%w(i, j, k) - jacEdgeD * var%w(i, j, k - 1)) / dz) / jac(i, j, &
          &k)
    else if((mu == 2 .and. nu == 3) .or. (mu == 3 .and. nu == 2)) then
      stressTensTFC = 0.5 * (vU - vD) / dz / jac(i, j, k) + 0.5 * (wF - wB) &
          &/ dy + met(i, j, k, 2, 3) * (vertWindTFC(i, j, k, var) &
          &- vertWindTFC(i, j, k - 1, var)) / dz
    else if(mu == 3 .and. nu == 3) then
      stressTensTFC = 2.0 * (vertWindTFC(i, j, k, var) - vertWindTFC(i, j, k &
          &- 1, var)) / dz / jac(i, j, k) - 2.0 / 3.0 * ((jacEdgeR * var%u(i, &
          &j, k) - jacEdgeL * var%u(i - 1, j, k)) / dx + (jacEdgeF * var%v(i, &
          &j, k) - jacEdgeB * var%v(i, j - 1, k)) / dy + (jacEdgeU * var%w(i, &
          &j, k) - jacEdgeD * var%w(i, j, k - 1)) / dz) / jac(i, j, k)
    end if

  end function stressTensTFC

end module atmosphere_module
