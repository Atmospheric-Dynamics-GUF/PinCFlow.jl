module atmosphere_module

  use type_module
  use sizeof_module

  implicit none

  public ! all objects are known to programmes using this module

  !--------------
  !    public
  !--------------
  public :: init_atmosphere
  public :: terminate_atmosphere

  ! TFC FJ
  public :: setHalosOfField2D
  public :: jac, met, chris, heightTFC, vertWindTFC, trafoTFC, stressTensTFC

  ! FJFeb2023
  public :: update_topography

  real, dimension(:), allocatable :: PStrat, rhoStrat, thetaStrat, bvsStrat
  real, dimension(:), allocatable :: PStrat_0, rhoStrat_0
  real, dimension(:), allocatable :: rhoStrat_d, rhoStrat_s

  real, dimension(:), allocatable :: pistrat
  real, dimension(:), allocatable :: PStratTilde, rhoStratTilde, thetaStratTilde

  real, dimension(:), allocatable :: PStrat00, PStrat01, rhoStrat00, &
      rhoStrat01, thetaStrat00, thetaStrat01, bvsStrat00, bvsStrat01, &
      PStratTilde00, PStratTilde01, rhoStratTilde00, rhoStratTilde01, &
      thetaStratTilde00, thetaStratTilde01

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
  !FS real :: Ro, RoInv              ! Rossby number and its inverse

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

  !UAB
  ! Held-Suarez atmosphere
  real :: tp_strato ! stratosphere temperature
  real :: tp_srf_trp ! tropical surface temperature
  real :: tpdiffhor_tropo ! tropospheric temperature difference
  ! between poles and tropics
  real :: ptdiffvert_tropo ! vertical potential-temperature
  ! difference in troposphere
  !UAE

  ! isothermal
  real :: T0 ! scaled background temperature

  ! general
  real :: p0 ! scaled reference pressure at z=0

  real :: zk ! zk = z(k)

  ! achatzb
  real :: mountainHeight, mountainWidth, k_mountain
  real :: x_center, y_center
  ! achatze

  ! TFC FJ
  ! 3D background fields.
  real, dimension(:, :, :), allocatable :: pStratTFC, thetaStratTFC, &
      rhoStratTFC, bvsStratTFC, piStratTFC

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
    real :: gammaInv
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

    !gagarinab
    ! for baroclinic case
    real :: T_c_b1, pow_t, pow_s, p_t_b ! tropopause quantities
    real :: T_bar, T_c_b, p_bar ! calculated quantities
    real :: tp_sponge, tp_sp1, tp_sp2 !FS
    !gagarinae

    !UAB
    real :: pistar, thetastar
    !UAC

    ! debugging
    integer, parameter :: errorlevel = 10 ! 0 -> no output

    integer :: i00, j00
    real :: yloc, ymax
    real, dimension(0:ny + 1) :: f_Coriolis_y

    real :: topos_u, topos_v, topos_w, x_sp_u, y_sp_u, z_sp_u, x_sp_v, y_sp_v, &
        z_sp_v, x_sp_w, y_sp_w, z_sp_w, dRP_u, dIP_u, dRP_v, dIP_v, dRP_w, dIP_w

    real :: rhodl, pcoeff_r1, p_dim1

    ! allocate PStrat
    if(.not. allocated(pStrat)) then
      allocate(Pstrat(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate pStrat"
    end if

    !UAB 200413
    ! allocate pStrat_0
    if(.not. allocated(pStrat_0)) then
      allocate(pStrat_0(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) then
        stop "atmosphere.f90: could not allocate pStrat_0"
      end if
    end if
    !UAE 200413

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

    !UAB 200413
    ! allocate rhoStrat_0
    if(.not. allocated(rhoStrat_0)) then
      allocate(rhoStrat_0(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) then
        stop "atmosphere.f90: could not allocate rhoStrat_0"
      end if
    end if
    !UAE 200413

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
      !UAC allocate( bvsStrat(-1:nz+1),stat=allocstat)
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
      !UAC allocate( bvsStrat00(-1:nz+1),stat=allocstat)
      allocate(bvsStrat00(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate bvsStrat"
    end if

    ! allocate squared Brunt-Vaisala frequency bvsStrat
    if(.not. allocated(bvsStrat01)) then
      !UAC allocate( bvsStrat01(-1:nz+1),stat=allocstat)
      allocate(bvsStrat01(- 1:nz + 2), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate bvsStrat"
    end if

    ! allocate Ro !FS
    if(.not. allocated(Ro)) then
      allocate(Ro(0:ny + 1), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate Ro"
    end if

    ! allocate RoInv !FS
    if(.not. allocated(RoInv)) then
      allocate(RoInv(0:ny + 1), stat = allocstat)
      if(allocstat /= 0) stop "atmosphere.f90: could not allocate RoInv"
    end if

    ! TFC FJ
    ! Allocate 3D background fields.
    if(topography) then
      if(.not. allocated(pStratTFC)) then
        allocate(pStratTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- 1):(nz &
            + 2)), stat = allocstat)
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate pStratTFC"
      end if

      if(.not. allocated(thetaStratTFC)) then
        allocate(thetaStratTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- &
            1):(nz + 2)), stat = allocstat)
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate &
            thetaStratTFC"
      end if

      if(.not. allocated(rhoStratTFC)) then
        allocate(rhoStratTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- 1):(nz &
            + 2)), stat = allocstat)
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate rhoStratTFC"
      end if

      if(.not. allocated(bvsStratTFC)) then
        allocate(bvsStratTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- 1):(nz &
            + 2)), stat = allocstat)
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate bvsStratTFC"
      end if

      if(.not. allocated(piStratTFC)) then
        allocate(piStratTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- 1):(nz &
            + 2)), stat = allocstat)
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate piStratTFC"
      end if
    end if

    ! save the initial tracer distributions in
    ! initialtracer and initialtracerrho = initialtracer*rho
    if(include_tracer) then
      if(.not. allocated(initialtracer)) then
        allocate(initialtracer(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate &
            initialtracer"
      end if
      if(.not. allocated(initialtracerrho)) then
        allocate(initialtracerrho(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz &
            + nbz))
        if(allocstat /= 0) stop "atmosphere.f90: could not allocate &
            initialtracerrho"
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
      ! if (testCase == "smoothVortex")then
      !    rhoRef = 0.5!1.184             ! in kg/m^3
      !    pRef = 101625.!101325.0            ! in Pa = kg/m/s^2
      ! else
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

    ! Rossby number
    !   achatzb correction for zero Coriolis
    !   Ro = uRef/f_Coriolis_dim/lRef
    !   RoInv = 1.0/Ro
    ! FS: see now further below

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
    do i = - nbx, sizeX + nbx ! modified by Junhong Wei (20161104)
      x(i) = lx(0) + real(i - 1) * dx + dx / 2.0
    end do

    do j = - nby, sizeY + nby ! modified by Junhong Wei (20161104)
      y(j) = ly(0) + real(j - 1) * dy + dy / 2.0
    end do

    do k = - nbz, sizeZ + nbz ! modified by Junhong Wei (20161104)
      z(k) = lz(0) + real(k - 1) * dz + dz / 2.0
    end do

    j00 = js + nby - 1 !FS
    if(TestCase == "baroclinic_LC") then
      ymax = ly_dim(1) / lRef
      do j = 0, ny + 1
        yloc = y(j + j00)
        f_Coriolis_y(j) = f_Coriolis_dim !FSJuly2020*sin(4*atan(1.0)*yloc/ymax)
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
    ! TFC FJ
    ! Topography setup.
    if(topography) then
      if(lz(0) /= 0.0) stop "lz(0) must be zero for TFC!"

      if(raytracer .and. case_wkb == 3) then
        topography_surface = 0.0

        if(mountain_case_wkb == 0) then
          call read_topography
        else
          topography_surface = 0.5 * mountainHeight_wkb_dim / lRef
        end if

        call setHalosOfField2D(topography_surface)

        if(mountain_case_wkb /= 0) then
          call output_topography
        end if
      else
        mountainHeight = mountainHeight_dim / lRef
        mountainWidth = mountainWidth_dim / lRef

        x_center = 0.5 * (lx(1) + lx(0))
        y_center = 0.5 * (ly(1) + ly(0))

        k_mountain = pi / mountainWidth

        i00 = is + nbx - 1
        j00 = js + nby - 1

        topography_surface = 0.0

        do j = - nby, ny + nby
          do i = - nbx, nx + nbx
            if(mountain_case == 0) then ! custom topography
              exit
            else if(mountain_case == 1) then ! cosine mountains
              topography_surface(i, j) = 0.5 * mountainHeight * (1.0 &
                  + cos(k_mountain * (x(i + i00) - x_center)))
            else if(mountain_case == 2) then ! 3d cosine mountains
              topography_surface(i, j) = 0.5 * mountainHeight * (1.0 &
                  + cos(k_mountain / sqrt(2.0) * (x(i + i00) - x_center + y(j &
                  + j00 - y_center))))
            else if(mountain_case == 3) then ! bell mountain
              topography_surface(i, j) = mountainHeight / (1.0 + ((x(i + i00) &
                  - x_center) ** 2.0) / (mountainWidth ** 2.0))
            else if(mountain_case == 4) then ! 3d bell mountain
              topography_surface(i, j) = mountainHeight / (1.0 + ((x(i + i00) &
                  - x_center) ** 2.0 + (y(j + j00) - y_center) ** 2.0) &
                  / (mountainWidth ** 2.0)) ** 1.5
            else if(mountain_case == 5) then ! simple mountain range
              topography_surface(i, j) = mountainHeight * exp(- ((x(i + i00) &
                  - x_center) / mountainWidth) ** 2.0) * cos(5.0 * k_mountain &
                  * (x(i + i00) - x_center) / 4.0) ** 2.0
            else if(mountain_case == 6) then ! cosine mountain range
              if(abs(x(i + i00) - x_center) <= mountainWidth) then
                topography_surface(i, j) = 0.5 * mountainHeight * (1.0 + 0.5 &
                    * (1.0 + cos(k_mountain * (x(i + i00) - x_center))) &
                    * cos(range_factor * k_mountain * (x(i + i00) - x_center)))
              else
                topography_surface(i, j) = 0.5 * mountainHeight
              end if
            else if(mountain_case == 7) then ! 3d cosine mountain range
              if(abs(x(i + i00) - x_center) <= mountainWidth .and. abs(y(j &
                  + j00) - y_center) <= mountainWidth) then
                topography_surface(i, j) = 0.5 * mountainHeight * (1.0 + 0.25 &
                    * (1.0 + cos(k_mountain * (x(i + i00) - x_center))) * (1.0 &
                    + cos(k_mountain * (y(j + j00) - y_center))) &
                    * cos(range_factor * k_mountain * (x(i + i00) - x_center &
                    + y(j + j00) - y_center)))
              else
                topography_surface(i, j) = 0.5 * mountainHeight
              end if
            else if(mountain_case == 8) then ! gaussian mountain range
              topography_surface(i, j) = 0.5 * mountainHeight * (1.0 + exp(- &
                  ((x(i + i00) - x_center) / mountainWidth) ** 2.0) * cos(0.5 &
                  * range_factor * k_mountain * (x(i + i00) - x_center)))
            else if(mountain_case == 9) then ! 3d gaussian mountain range
              topography_surface(i, j) = 0.5 * mountainHeight * (1.0 + exp(- &
                  ((x(i + i00) - x_center) / mountainWidth) ** 2.0 - ((y(j &
                  + j00) - y_center) / mountainWidth) ** 2.0) * cos(0.5 &
                  * range_factor * k_mountain * (x(i + i00) - x_center + y(j &
                  + j00) - y_center)))
            else if(mountain_case == 10) then ! localized cosine mountain
              if(abs(x(i + i00) - x_center) <= 1.5 * mountainWidth) then
                topography_surface(i, j) = 0.5 * mountainHeight * (1.0 &
                    + cos(k_mountain * (x(i + i00) - x_center)))
              else
                topography_surface(i, j) = 0.5 * mountainHeight
              end if
            else
              if(master) stop "Mountain case not defined!"
            end if
            if(topography_surface(i, j) .lt. lz(0)) then
              print *, "Topography has negative values."
              print *, "lz_dim(0) should be set accordingly."
              stop
            end if
          end do
        end do

        ! Read topography data.
        if(mountain_case == 0) then
          call read_topography
        end if

        ! Ensure periodicity of topography.
        call setHalosOfField2D(topography_surface)

        ! Output topography data.
        if(mountain_case /= 0) then
          call output_topography
        end if
      end if

      ! FJFeb2023
      if(topographyTime > 0.0) then
        allocate(final_topography_surface(- nbx:nx + nbx, - nby:ny + nby), &
            stat = allocstat)
        if(allocstat /= 0) stop "init.f90: could not allocate &
            final_topography_surface"

        final_topography_surface = topography_surface
        topography_surface = 0.0
      end if

      ! kbl_topo(:,:,:) = nz + 2*nbz
      !
      ! do k = 1, nz
      !    do j = -nby, ny+nby-1
      !       do i = -nbx, nx+nbx-1
      !          ! vertical index of first u-point above the surface
      !          topos_u &
      !          = 0.5*(topography_surface(i,j) + topography_surface(i+1,j))
      !          if(z(k-1) <= topos_u .and. z(k) > topos_u) then
      !             kbl_topo(i,j,1) = k
      !          end if
      !
      !          ! vertical index of first v-point above the surface
      !          topos_v &
      !          = 0.5*(topography_surface(i,j) + topography_surface(i,j+1))
      !          if(z(k-1) <= topos_v .and. z(k) > topos_v) then
      !             kbl_topo(i,j,2) = k
      !          end if
      !
      !          ! vertical index of first w-point above the surface
      !          topos_w = topography_surface(i,j)
      !          if(z(k-1)+dz/2. <= topos_w .and. z(k)+dz/2. > topos_w) then
      !             kbl_topo(i,j,3) = k
      !          end if
      !       end do
      !    end do
      ! end do

      ! do j = -nby, ny+nby-1
      !    do i = -nbx, nx+nbx-1
      !       if (kbl_topo(i,j,1) < 1 .or. kbl_topo(i,j,1) > nz) then
      !       print*,'kbl_topo(',i,',',j,',1) =',kbl_topo(i,j,1), &
      !            & ' out of range'
      !       stop
      !       end if
      !
      !       if (kbl_topo(i,j,2) < 1 .or. kbl_topo(i,j,2) > nz) then
      !       print*,'kbl_topo(',i,',',j,',2) =',kbl_topo(i,j,2), &
      !            & ' out of range'
      !       stop
      !       end if
      !
      !       if (kbl_topo(i,j,3) < 1 .or. kbl_topo(i,j,3) > nz) then
      !       print*,'kbl_topo(',i,',',j,',3) =',kbl_topo(i,j,3), &
      !            & ' out of range'
      !       stop
      !       end if
      !    end do
      ! end do

      ! do j = 1, ny
      !    do i = 1, nx
      !       !--------------------------------------------------
      !       ! fields for reconstructing u just above the surface
      !       !--------------------------------------------------
      !
      !       k = kbl_topo(i,j,1)
      !
      !       ! gradient of topography below u-reconstruction point
      !
      !       dhdx(i,j,1) &
      !       = (topography_surface(i+1,j) - topography_surface(i,j))/dx
      !       dhdy(i,j,1) &
      !       = (  topography_surface(i,j+1) + topography_surface(i+1,j+1)&
      !          - topography_surface(i,j-1) - topography_surface(i+1,j-1)) &
      !         /(4.*dy)
      !
      !       ! height of topography below u-reconstruction point
      !
      !       topos_u &
      !       = 0.5*(topography_surface(i,j) + topography_surface(i+1,j))
      !
      !       ! coordinates of surface point connected by its surface normal
      !       ! with the u-reconstruction point
      !
      !       z_sp_u &
      !       = z(k) + (topos_u - z(k))/(1. + dhdx(i,j,1)**2 + dhdy(i,j,1)**2)
      !       x_sp_u = x(i00+i) + dx/2. - dhdx(i,j,1)*(z_sp_u - z(k))
      !       y_sp_u = y(j00+j) - dhdy(i,j,1)*(z_sp_u - z(k))
      !
      !       ! coordinates of free-atmosphere interpolation point connected
      !       ! by the surface normal to the u-reconstruction point
      !
      !       z_ip(i,j,1) = z(k+1)
      !
      !       if (abs(dhdx(i,j,1)*dz) > dx) then
      !          print*,'abs(dhdx(',i,',',j,',1)*dz) > dx'
      !          print*,'vertical grid to be refined!'
      !          stop
      !         else
      !          x_ip(i,j,1) = x(i00+i) + dx/2. - dhdx(i,j,1)*dz
      !       end if
      !
      !       if (abs(dhdy(i,j,1)*dz) > dy) then
      !          print*,'abs(dhdy(',i,',',j,',1)*dz) > dy'
      !          print*,'vertical grid to be refined!'
      !          stop
      !         else
      !          y_ip(i,j,1) = y(j00+j) - dhdy(i,j,1)*dz
      !       end if
      !
      !       ! distance of surface point to u-reconstruction and
      !       ! interpolation point
      !
      !       dRP_u = sqrt(  (x_sp_u - x(i00+i) - dx/2.)**2 &
      !                    + (y_sp_u - y(j00+j))**2 &
      !                    + (z_sp_u - z(k))**2)
      !
      !       dIP_u = sqrt(  (x_ip(i,j,1) - x_sp_u)**2 &
      !                    + (y_ip(i,j,1) - y_sp_u)**2 &
      !                    + (z_ip(i,j,1) - z_sp_u)**2)
      !
      !       ! factors for determining the normal and tangential velocity at
      !       ! the u-reconstruction point from the corresponding velocities
      !       ! at the interpolation point
      !
      !       if (z_0 > 0.) then
      !          if (dRP_u < z_0) then
      !             velocity_reconst_t(i,j,1) = 0.0
      !            else
      !             velocity_reconst_t(i,j,1) = log(dRP_u/z_0)/log(dIP_u/z_0)
      !          endif
      !         else
      !          velocity_reconst_t(i,j,1) = 1.0
      !       end if
      !
      !       velocity_reconst_n(i,j,1)=dRP_u/dIP_u
      !
      !       !--------------------------------------------------
      !       ! fields for reconstructing v just above the surface
      !       !--------------------------------------------------
      !
      !       k = kbl_topo(i,j,2)
      !
      !       ! gradient of topography below v-reconstruction point
      !
      !       dhdx(i,j,2) &
      !       = (  topography_surface(i+1,j) + topography_surface(i+1,j+1) &
      !          - topography_surface(i-1,j) - topography_surface(i-1,j+1)) &
      !          /(4.*dx)
      !       dhdy(i,j,2) &
      !       = (topography_surface(i,j+1) - topography_surface(i,j))/dy
      !
      !       ! height of topography below v-reconstruction point
      !
      !       topos_v &
      !       = 0.5*(topography_surface(i,j) + topography_surface(i,j+1))
      !
      !       ! coordinates of surface point connected by its surface normal
      !       ! with the v-reconstruction point
      !
      !       z_sp_v &
      !       = z(k) + (topos_v - z(k))/(1. + dhdx(i,j,2)**2 + dhdy(i,j,2)**2)
      !       x_sp_v = x(i00+i) - dhdx(i,j,2)*(z_sp_v - z(k))
      !       y_sp_v = y(j00+j) + dy/2. - dhdy(i,j,2)*(z_sp_v - z(k))
      !
      !       ! coordinates of free-atmosphere interpolation point connected
      !       ! by the surface normal to the v-reconstruction point
      !
      !       z_ip(i,j,2) = z(k+1)
      !
      !       if (abs(dhdx(i,j,2)*dz) > dx) then
      !          print*,'abs(dhdx(',i,',',j,',2)*dz) > dx'
      !          print*,'vertical grid to be refined!'
      !          stop
      !         else
      !          x_ip(i,j,2) = x(i00+i) - dhdx(i,j,2)*dz
      !       end if
      !
      !       if (abs(dhdy(i,j,2)*dz) > dy) then
      !          print*,'abs(dhdy(',i,',',j,',2)*dz) > dy'
      !          print*,'vertical grid to be refined!'
      !          stop
      !         else
      !          y_ip(i,j,2) = y(j00+j) + dy/2. - dhdy(i,j,2)*dz
      !       end if
      !
      !       ! distance of surface point to v-reconstruction and
      !       ! interpolation point
      !
      !       dRP_v = sqrt(  (x_sp_v - x(i00+i))**2 &
      !                    + (y_sp_v - y(j00+j) - dy/2.)**2 &
      !                    + (z_sp_v - z(k))**2)
      !
      !       dIP_v = sqrt(  (x_ip(i,j,2) - x_sp_v)**2 &
      !                    + (y_ip(i,j,2) - y_sp_v)**2 &
      !                    + (z_ip(i,j,2) - z_sp_v)**2)
      !
      !       ! factors for determining the normal and tangential velocity at
      !       ! the v-reconstruction point from the corresponding velocities
      !       ! at the interpolation point
      !
      !       if (z_0 > 0.) then
      !          if (dRP_v < z_0) then
      !             velocity_reconst_t(i,j,2) = 0.0
      !            else
      !             velocity_reconst_t(i,j,2) = log(dRP_v/z_0)/log(dIP_v/z_0)
      !          endif
      !         else
      !          velocity_reconst_t(i,j,2) = 1.0
      !       end if
      !
      !       velocity_reconst_n(i,j,2)=dRP_v/dIP_v
      !
      !       !--------------------------------------------------
      !       ! fields for reconstructing w just above the surface
      !       !--------------------------------------------------
      !
      !       k = kbl_topo(i,j,3)
      !
      !       ! gradient of topography below w-reconstruction point
      !
      !       dhdx(i,j,3) &
      !       = (topography_surface(i+1,j) - topography_surface(i-1,j))/(2.*dx)
      !       dhdy(i,j,3) &
      !       = (topography_surface(i,j+1) - topography_surface(i,j-1))/(2.*dy)
      !
      !       ! height of topography below w-reconstruction point
      !
      !       topos_w = topography_surface(i,j)
      !
      !       ! coordinates of surface point connected by its surface normal
      !       ! to the w-reconstruction point
      !
      !       z_sp_w &
      !       =   z(k) + dz/2. &
      !         + (topos_w - z(k) - dz/2.) &
      !           /(1. + dhdx(i,j,3)**2 + dhdy(i,j,3)**2)
      !       x_sp_w = x(i00+i) - dhdx(i,j,3)*(z_sp_w - z(k) - dz/2.)
      !       y_sp_w = y(j00+j) - dhdy(i,j,3)*(z_sp_w - z(k) - dz/2.)
      !
      !       ! coordinates of free-atmosphere interpolation point connected
      !       ! by the surface normal to the w-reconstruction point
      !
      !       z_ip(i,j,3) = z(k+1) + dz/2.
      !
      !       if (abs(dhdx(i,j,3)*dz) > dx) then
      !          print*,'abs(dhdx(',i,',',j,',3)*dz) > dx'
      !          print*,'vertical grid to be refined!'
      !          stop
      !         else
      !          x_ip(i,j,3) = x(i00+i) - dhdx(i,j,3)*dz
      !       end if
      !
      !       if (abs(dhdy(i,j,3)*dz) > dy) then
      !          print*,'abs(dhdy(',i,',',j,',3)*dz) > dy'
      !          print*,'vertical grid to be refined!'
      !          stop
      !         else
      !          y_ip(i,j,3) = y(j00+j) - dhdy(i,j,3)*dz
      !       end if
      !
      !       ! distance of surface point to w-reconstruction and
      !       ! interpolation point
      !
      !       dRP_w = sqrt(  (x_sp_w - x(i00+i))**2 + (y_sp_w - y(j00+j))**2 &
      !                    + (z_sp_w - z(k) - dz/2.)**2)
      !
      !       dIP_w = sqrt(  (x_ip(i,j,3) - x_sp_w)**2 &
      !                    + (y_ip(i,j,3) - y_sp_w)**2 &
      !                    + (z_ip(i,j,3) - z_sp_w)**2)
      !
      !       ! factors for determining the normal and tangential velocity at
      !       ! the w-reconstruction point from the corresponding velocities
      !       ! at the interpolation point
      !
      !       if (z_0 > 0.) then
      !          if (dRP_w < z_0) then
      !             velocity_reconst_t(i,j,3) = 0.0
      !            else
      !             velocity_reconst_t(i,j,3) = log(dRP_w/z_0)/log(dIP_w/z_0)
      !          endif
      !         else
      !          velocity_reconst_t(i,j,3) = 1.0
      !       end if
      !
      !       velocity_reconst_n(i,j,3)=dRP_w/dIP_w
      !    end do
      ! end do
    end if !topography

    !UAE

    !---------------------------------------------
    !   Set up Sponge layer
    !---------------------------------------------

    if(spongeLayer) then
      if(verticalSponge == "exponential") then
        if(timescheme == "semiimplicit") stop "ERROR: Wrong time scheme for &
            exponential sponge!"
        ksponge = 1
        zSponge = Rsp * Temp0_dim / g / lRef
      else
        kSponge = nz - ceiling(spongeHeight * real(nz))
        zSponge = lz(0) + (1.0 - spongeHeight) * (lz(1) - lz(0))
      end if

      ! Setup lateral sponge.
      if(lateralSponge) then
        xSponge0 = lx(0) + 0.5 * spongeHeight * (lx(1) - lx(0))
        ySponge0 = ly(0) + 0.5 * spongeHeight * (ly(1) - ly(0))
        xSponge1 = lx(1) - 0.5 * spongeHeight * (lx(1) - lx(0))
        ySponge1 = ly(1) - 0.5 * spongeHeight * (ly(1) - ly(0))
      end if
    end if

    select case(model)
      !--------------------------------------------------------------------
      !             Atmospheres for pseudo-incompressible
      !--------------------------------------------------------------------

    case("pseudo_incompressible", "WKB")

      select case(background)

      case('smoothV')

        ! some stratified fields are still in use...
        rhodl = 1. / 0.5
        pcoeff_r1 = 1024. ** 2 * (1. / 72. - 6. / 35. + 15. / 17. - 74. / 33. &
            + 57. / 32. + 174. / 31. - 269. / 15. + 450. / 29. + 153. / 8. &
            - 1564. / 27. + 510. / 13. + 204. / 5. - 1. / 24. * (2210. &
            - rhodl) + 12. / 23. * (85. - rhodl) + (510. / 11. + 3. * rhodl) &
            - 4. / 21. * (391. + 55. * rhodl) + 9. / 40. * (119. + 110. &
            * rhodl) + 18. / 19. * (25. - 44. * rhodl) - 1. / 9. * (269. - 462 &
            * rhodl) + 6. / 17. * (29. - 132 * rhodl) + 3. / 16. * (19. + 165. &
            * rhodl) - 2. / 15. * (37. + 110. * rhodl) + 3. / 7. * (5. + 11. &
            * rhodl) - 6. / 13. * (1. + 2. * rhodl) + 1. / 24. * (1. + 2. &
            * rhodl))

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

        ! TFC FJ
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
            referenceQuantities = SI not impl.."
        if(.not. fluctuationMode) stop "atmosphere.f90: only fluctuationMode &
            = TRUE impl.."

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
                / gamma / T_tr * delZ)
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
                  sig / gamma / T_tr * delZ)
              rhoStratTilde(k) = PStratTilde(k) / thetaStratTilde(k)
            end if
          else
            thetaStratTilde(k) = 0.5 * (thetaStrat(k) + thetaStrat(k + 1))
            PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k + 1))
            rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k + 1))
          end if
        end do ! k loop

        ! TFC FJ
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
                      j, k)
                else
                  ! Isothermal stratosphere.
                  ! Define pStratTFC.
                  pStratTFC(i, j, k) = p0 ** kappa * press_tr ** (1.0 / gamma) &
                      * exp(- sig / gamma / T_tr * (heightTFC(i, j, k) - z_tr))
                  ! Define thetaStratTFC.
                  thetaStratTFC(i, j, k) = theta_tr * exp(kappa * sig / T_tr &
                      * (heightTFC(i, j, k) - z_tr))
                  ! Define rhoStratTFC.
                  rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                      j, k)
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

          if(fluctuationMode) stop "init_atmosphere: fluctuationMode not impl. &
              for SI!"
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

          ! Removed inconsistency (FJJan2023).
          ! !xxx need rhoStratTilde(-1) in fluxes.f90, line 1927 \pm
          ! rhoStratTilde(- 1) = rhoStratTilde(0)
          ! rhoStratTilde(nz + 1) = rhoStratTilde(nz)

          ! TFC FJ
          ! Define 3D background fields.
          if(topography) then
            do i = - nbx, nx + nbx
              do j = - nby, ny + nby
                do k = - 1, nz + 2
                  ! Define pStratTFC.
                  pStratTFC(i, j, k) = p0 * exp(- sig * heightTFC(i, j, k) &
                      / gamma / T0)
                  ! Define thetaStratTFC.
                  thetaStratTFC(i, j, k) = T0 * exp(kappa * sig / T0 &
                      * heightTFC(i, j, k))
                  ! Define rhoStratTFC.
                  rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                      j, k)
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

        if(include_ice .and. (iceTestcase == "1D_ISSR")) then
          theta0_dim = T_nuc + 5.0
          press0_dim = p_nuc * (theta0_dim / T_nuc) ** kappaInv
          ! scaled reference pressure at z = 0
          p0 = press0_dim / pRef
        end if

        if(referenceQuantities == "SI") then
          !------------------------------------
          !   original equations in SI units
          !------------------------------------

          if(fluctuationMode) then
            stop "init_atmosphere: fluctuationMode not implmented for SI!"
          end if

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

          ! TFC FJ
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
                      j, k)
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

        if(include_ice .and. (iceTestcase == "1D_ISSR")) then
          theta0_dim = T_nuc + 5.0
          term = kappa * g ** 2 / (Rsp * N_BruntVaisala_dim ** 2)
          press0_dim = p_nuc / (1. + term / theta0_dim * ((theta0_dim - term) &
              / (T_nuc - term) - 1.)) ** kappaInv
          ! scaled reference pressure at z = 0
          p0 = press0_dim / pRef
        end if

        if(referenceQuantities == "SI") then
          !------------------------------------
          !   original equations in SI units
          !------------------------------------
          if(fluctuationMode) then
            stop "init_atmosphere: fluctuationMode not implmented for SI!"
          end if

          theta0 = theta0_dim / thetaRef ! theta0 at z=0 in K
          NN = N_BruntVaisala_dim * tRef
          N2 = NN ** 2 ! Brunt-Vaisala N0^2 in 1/s^2

          coeff = kappa * g ** 2 / (Rsp * N2 * theta0)
          power = 1. / (gamma - 1.)

          do k = - 1, nz + 2
            term = exp(- N2 / g * z(k))

            !if( term > 1.0 ) stop "init_atmosphere: root of a negative number."
            if(1.0 + coeff * (term - 1.0) < 0.0) stop "init_atmosphere: root &
                of a negative number."

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
                then
              print *, "init_atmosphere: power of a neg. number."
              stop 'top of atmosphere too high for const-N'
            end if

            Pstrat(k) = p0 * (1.0 + FrInv2 * kappa * sig / N2 / theta0 * (term &
                - 1.0)) ** power

            rhoStrat(k) = pStrat(k) / thetaStrat(k)

            !testb
            !print*,'Pstrat(',k,') =',Pstrat(k)
            !teste
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
                  < 0.0) then
                print *, "init_atmosphere: power of a neg. number."
                stop 'top of atmosphere too high for const-N'
              end if

              PstratTilde(k) = p0 * (1.0 + FrInv2 * kappa * sig / N2 / theta0 &
                  * (term - 1.0)) ** power
              rhoStratTilde(k) = PstratTilde(k) / thetaStratTilde(k)
            else
              thetaStratTilde(k) = 0.5 * (thetaStrat(k) + thetaStrat(k + 1))
              PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k + 1))
              rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k + 1))
            end if
          end do

          ! TFC FJ
          ! Define 3D background fields.
          if(topography) then
            do i = - nbx, nx + nbx
              do j = - nby, ny + nby
                do k = - 1, nz + 2
                  power = 1.0 / (gamma - 1.0)
                  term = exp(- Fr2 * N2 * heightTFC(i, j, k))
                  ! Define pStratTFC.
                  pStratTFC(i, j, k) = p0 * (1.0 + FrInv2 * kappa * sig / N2 &
                      / theta0 * (term - 1.0)) ** power
                  ! Define thetaStratTFC.
                  thetaStratTFC(i, j, k) = theta0 * exp(Fr ** 2.0 * N2 &
                      * heightTFC(i, j, k))
                  ! Define rhoStratTFC.
                  rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                      j, k)
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
                    ** pow_t
              else
                p_bar = press0_dim * exp(- zk_half / (Rsp * Temp0_dim / g))
              end if
            else
              T_bar = Temp0_dim - gamma_t * z_tr_dim - gamma_s * (zk_half &
                  - z_tr_dim)

              if(gamma_s /= 0.) then
                p_bar = p_t_b * (1. - gamma_s * (zk_half - z_tr_dim) / T_c_b) &
                    ** pow_s
              else
                p_bar = p_t_b * exp(- (zk_half - z_tr_dim) / (Rsp * T_c_b / g))
              end if
            endif

            thetaStratTilde(k) = T_bar * (press0_dim / p_bar) ** (kappa) &
                / thetaRef

            rhoStratTilde(k) = p_bar / (Rsp * T_bar) / rhoRef

            PstratTilde(k) = rhoStratTilde(k) * thetaStratTilde(k)
          else
            PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k + 1))

            rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k + 1))

            thetaStratTilde(k) = PStratTilde(k) / rhoStratTilde(k)
          end if
        enddo

        ! TFC FJ
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
                        ** pow_t
                  else
                    p_bar = press0_dim * exp(- zk / (Rsp * Temp0_dim / g))
                  end if
                else
                  T_bar = Temp0_dim - gamma_t * z_tr_dim - gamma_s * (zk &
                      - z_tr_dim)
                  if(gamma_s /= 0.0) then
                    p_bar = p_t_b * (1.0 - gamma_s * (zk - z_tr_dim) / T_c_b) &
                        ** pow_s
                  else
                    p_bar = p_t_b * exp(- (zk - z_tr_dim) / (Rsp * T_c_b / g))
                  end if
                end if

                ! Define thetaStratTFC.
                thetaStratTFC(i, j, k) = T_bar * (press0_dim / p_bar) ** kappa &
                    / thetaRef
                ! Define rhoStratTFC.
                rhoStratTFC(i, j, k) = p_bar / (Rsp * T_bar) / rhoRef
                ! Define pStratTFC.
                pStratTFC(i, j, k) = rhoStratTFC(i, j, k) * thetaStratTFC(i, &
                    j, k)
              end do
            end do
          end do
        end if

        !UAB
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

        !testb
        !print*,tp_strato, tp_srf_trp - 0.5*tpdiffhor_tropo, T_bar
        !teste

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
              * tpdiffhor_tropo - 0.5 * ptdiffvert_tropo / kappa &
              * log(pistrat(k)))) !FS 0.5->1

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
              - 0.5 * ptdiffvert_tropo / kappa * log(pistar))) !FS 0.5->1

          thetastar = T_bar / pistar

          pistrat(k) = pistrat(k - 1) - 0.5 * dz * (kappa / thetastar + kappa &
              / thetaStrat(k - 1))

          if(pistrat(k) <= 0.) then
            print *, 'ERROR: negative non-positive pistrat at k =', k
            stop
          end if

          T_bar = max(tp_strato, pistrat(k) * (tp_srf_trp - 0.5 &
              * tpdiffhor_tropo - 0.5 * ptdiffvert_tropo / kappa &
              * log(pistrat(k)))) !FS 0.5->1

          thetaStrat(k) = T_bar / pistrat(k)

          pStrat(k) = pistrat(k) ** ((1.0 - kappa) / kappa)

          rhoStrat(k) = pStrat(k) / thetaStrat(k)

          ! tp_sponge = T_bar

        end do

        !UAB
        thetaStrat(- 1) = thetaStrat(0)
        pStrat(- 1) = pStrat(0)
        rhoStrat(- 1) = rhoStrat(0)
        pistrat(- 1) = pistrat(0)
        !UAE

        !      ! close jets below sponge layer !FS
        !  do k = nz+1-ceiling((0.4)*real(nz)), nz+2 !FS
        !      pistar = pistrat(k-2) - 2.0*dz * kappa/thetaStrat(k-1)
        ! &
        !              T_bar &
        !              =  tp_sponge + pistar*(0.5*tpdiffhor_tropo &
        !                        + 1.*ptdiffvert_tropo/kappa * log(pistar))!FS 0.5->1

        !              thetastar = T_bar/pistar

        !              pistrat(k) &
        !              =   pistrat(k-1) &
        !                - 0.5*dz * (kappa/thetastar + kappa/thetaStrat(k-1))

        !              T_bar &
        !              = tp_sponge + pistrat(k)*(0.5*tpdiffhor_tropo &
        !                        + 1.*ptdiffvert_tropo/kappa * log(pistrat(k)))!FS 0.5->1

        !              thetaStrat(k) = T_bar/pistrat(k)

        !              pStrat(k) = pistrat(k)**((1.0 - kappa)/kappa)

        !              rhoStrat(k) = pStrat(k)/thetaStrat(k)

        !           end do

        ! do k = (nz - ceiling(0.25*real(nz)))+1,nz+2 !FS
        !    thetaStrat(k) = thetaStrat(nz - ceiling(0.25*real(nz)))
        !    pStrat(k) = pStrat(nz - ceiling(0.25*real(nz)))
        !    piStrat(k) = piStrat(nz - ceiling(0.25*real(nz)))
        !    rhoStrat(k) = rhoStrat(nz - ceiling(0.25*real(nz)))

        ! end do

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
                - thetaStratTilde(k - 2)
          end if
        enddo

        ! TFC FJ
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
                piStratTFC(i, j, k) = 1.0 + heightTFC(i, j, k) * kappa / T_bar
                if(k == 1 .and. piStratTFC(i, j, k) <= 0.0) then
                  stop "ERROR: non-positive piStratTFC at k = 1"
                end if
                T_bar = max(tp_strato, piStratTFC(i, j, k) * (tp_srf_trp - 0.5 &
                    * tpdiffhor_tropo - 0.5 * ptdiffvert_tropo / kappa &
                    * log(piStratTFC(i, j, k))))
                ! Define pStratTFC.
                pStratTFC(i, j, k) = piStratTFC(i, j, k) ** ((1.0 - kappa) &
                    / kappa)
                ! Define thetaStratTFC.
                thetaStratTFC(i, j, k) = T_bar / piStratTFC(i, j, k)
                ! Define rhoStratTFC.
                rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                    j, k)
              end do

              ! Integrate upwards.
              do k = 2, nz + 2
                ! Leapfrog step.
                piStar = piStratTFC(i, j, k - 2) - 2.0 * dz * jac(i, j, k - 1) &
                    * kappa / thetaStratTFC(i, j, k - 1)
                if(piStar <= 0.0) then
                  print *, "ERROR: non-positive piStar at k =", k
                  stop
                end if
                T_bar = max(tp_strato, piStar * (tp_srf_trp - 0.5 &
                    * tpdiffhor_tropo - 0.5 * ptdiffvert_tropo / kappa &
                    * log(piStar)))
                thetaStar = T_bar / piStar
                ! Trapezoidal step.
                piStratTFC(i, j, k) = piStratTFC(i, j, k - 1) - 0.5 * dz * 0.5 &
                    * (jac(i, j, k) + jac(i, j, k - 1)) * (kappa / thetaStar &
                    + kappa / thetaStratTFC(i, j, k - 1))
                if(piStratTFC(i, j, k) <= 0.0) then
                  print *, "ERROR: non-positive piStratTFC at k =", k
                  stop
                end if
                T_bar = max(tp_strato, piStratTFC(i, j, k) * (tp_srf_trp - 0.5 &
                    * tpdiffhor_tropo - 0.5 * ptdiffvert_tropo / kappa &
                    * log(piStratTFC(i, j, k))))
                ! Define pStratTFC.
                pStratTFC(i, j, k) = piStratTFC(i, j, k) ** ((1.0 - kappa) &
                    / kappa)
                ! Define thetaStratTFC.
                thetaStratTFC(i, j, k) = T_bar / piStratTFC(i, j, k)
                ! Define rhoStratTFC.
                rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                    j, k)
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
          / dz

      bvsStrat(0) = g_ndim / thetaStrat(0) * (thetaStrat(1) - thetaStrat(0)) &
          / dz

      !UAB
      N2 = max(bvsStrat(- 1), bvsStrat(0))
      !UAE

      do k = 1, nz
        bvsStrat(k) = g_ndim / thetaStrat(k) * (thetaStrat(k + 1) &
            - thetaStrat(k - 1)) / (2.0 * dz)

        !UAB
        N2 = max(N2, bvsStrat(k))
        !UAE
      end do

      bvsStrat(nz + 1) = g_ndim / thetaStrat(nz + 1) * (thetaStrat(nz + 1) &
          - thetaStrat(nz)) / dz

      !UAB
      bvsStrat(nz + 2) = bvsStrat(nz + 1)

      N2 = max(N2, bvsStrat(nz + 1))

      if(N2 < 0.) then
        stop 'ERROR: N2 < 0'
      else
        NN = sqrt(N2)
      end if
      !UAE

      if(TestCase == "smoothVortex") then
        bvsstrat(:) = 0.
      end if

      ! TFC FJ
      ! Define bvsStratTFC.
      if(topography) then
        bvsStratTFC = 0.0
        do i = - nbx, nx + nbx
          do j = - nby, ny + nby
            ! Lower boundary.
            bvsStratTFC(i, j, - 1) = g_ndim / thetaStratTFC(i, j, 0) / jac(i, &
                j, 0) * (thetaStratTFC(i, j, 1) - thetaStratTFC(i, j, 0)) / dz
            bvsStratTFC(i, j, 0) = bvsStratTFC(i, j, - 1)
            ! Between boundaries.
            do k = 1, nz
              bvsStratTFC(i, j, k) = g_ndim / thetaStratTFC(i, j, k) / jac(i, &
                  j, k) * 0.5 * (thetaStratTFC(i, j, k + 1) - thetaStratTFC(i, &
                  j, k - 1)) / dz
            end do
            ! Upper boundary.
            bvsStratTFC(i, j, nz + 1) = g_ndim / thetaStratTFC(i, j, nz + 1) &
                / jac(i, j, nz + 1) * (thetaStratTFC(i, j, nz + 1) &
                - thetaStratTFC(i, j, nz)) / dz
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

      ! TFC FJ
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

    !UAB
    deallocate(pistrat, stat = allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not deallocate pistrat"
    !UAE

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

    ! TFC FJ
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

    ! TFC FJ
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
          mpi_double_precision, dest, tag, xSliceLeft_recv(i0, j0), recvcount, &
          mpi_double_precision, source, mpi_any_tag, comm, sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xSliceLeft_send(i0, j0), sendcount, &
          mpi_double_precision, dest, tag, xSliceRight_recv(i0, j0), &
          recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          sts_right, ierror)

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
          mpi_double_precision, dest, tag, ySliceBack_recv(i0, j0), recvcount, &
          mpi_double_precision, source, mpi_any_tag, comm, sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(ySliceBack_send(i0, j0), sendcount, &
          mpi_double_precision, dest, tag, ySliceForw_recv(i0, j0), recvcount, &
          mpi_double_precision, source, mpi_any_tag, comm, sts_forw, ierror)

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

  subroutine output_topography

    real * 4, dimension(sizeX, sizeY) :: topography_out
    real * 4, dimension(sizeX * nprocy, ny) :: topography_mst
    real * 4, dimension(nx, ny) :: topography_prc
    integer :: i_out, i_mst, i_prc, j_out, j_mst, j_prc
    integer :: i, j, k

    ! Open file.
    if(master) then
      open(42, file = "topography.dat", form = "unformatted", access &
          = "direct", recl = sizeX * sizeY * sizeofreal4)
    end if

    do j = 1, ny
      ! Dimensionalize.
      do i = 1, nx
        topography_prc(i, j) = topography_surface(i, j) * lRef
      end do
      ! Distribute data over all processors.
      call mpi_gather(topography_prc(1, j), nx, mpi_real, topography_mst(1, &
          j), nx, mpi_real, 0, comm, ierror)
    end do

    call mpi_barrier(comm, ierror)

    ! Output layerwise.
    if(master) then
      do j = 1, ny
        j_mst = j
        do j_prc = 1, nprocy
          j_out = ny * (j_prc - 1) + j
          do i_prc = 1, nprocx
            do i = 1, nx
              i_out = nx * (i_prc - 1) + i
              i_mst = nprocy * nx * (i_prc - 1) + (j_prc - 1) * nx + i
              topography_out(i_out, j_out) = topography_mst(i_mst, j_mst)
            end do
          end do
        end do
      end do
      write(42, rec = 1) topography_out
    end if

    ! Close file.
    if(master) close(unit = 42)

  end subroutine output_topography

  !---------------------------------------------------------------------------

  subroutine read_topography

    real * 4, dimension(sizeX, sizeY) :: topography_out
    real * 4, dimension(sizeX * nprocy, ny) :: topography_mst
    real * 4, dimension(nx, ny) :: topography_prc
    integer :: i_out, i_mst, i_prc, j_out, j_mst, j_prc
    integer :: i, j, k

    ! Open file.
    if(master) then
      open(42, file = "topography.dat", form = "unformatted", access &
          = "direct", recl = sizeX * sizeY * sizeofreal4)
    end if

    ! Read data.
    if(master) then
      read(42, rec = 1) topography_out
      do j = 1, ny
        j_mst = j
        do j_prc = 1, nprocy
          j_out = ny * (j_prc - 1) + j
          do i_prc = 1, nprocx
            do i = 1, nx
              i_out = nx * (i_prc - 1) + i
              i_mst = nprocy * nx * (i_prc - 1) + (j_prc - 1) * nx + i
              topography_mst(i_mst, j_mst) = topography_out(i_out, j_out)
            end do
          end do
        end do
      end do
    end if

    call mpi_barrier(comm, ierror)

    do j = 1, ny
      ! Distribute data over all processors.
      call mpi_scatter(topography_mst(1, j), nx, mpi_real, topography_prc(1, &
          j), nx, mpi_real, 0, comm, ierror)
      do i = 1, nx
        ! Non-dimensionalize.
        topography_surface(i, j) = topography_prc(i, j) / lRef
      end do
    end do

    if(master) close(unit = 42)

  end subroutine read_topography

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

    if(any(topography_surface /= final_topography_surface)) then

      ! Update topography.
      if(time < topographyTime / tRef) then
        topography_surface = time / topographyTime * tRef &
            * final_topography_surface
      else
        topography_surface = final_topography_surface
      end if

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
                    j, k)
              else
                ! Isothermal stratosphere.
                ! Define pStratTFC.
                pStratTFC(i, j, k) = p0 ** kappa * press_tr ** (1.0 / gamma) &
                    * exp(- sig / gamma / T_tr * (heightTFC(i, j, k) - z_tr))
                ! Define thetaStratTFC.
                thetaStratTFC(i, j, k) = theta_tr * exp(kappa * sig / T_tr &
                    * (heightTFC(i, j, k) - z_tr))
                ! Define rhoStratTFC.
                rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                    j, k)
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
                  / T0)
              ! Define thetaStratTFC.
              thetaStratTFC(i, j, k) = T0 * exp(kappa * sig / T0 &
                  * heightTFC(i, j, k))
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
                  / theta0 * (term - 1.0)) ** power
              ! Define thetaStratTFC.
              thetaStratTFC(i, j, k) = theta0 * exp(Fr ** 2.0 * N2 &
                  * heightTFC(i, j, k))
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
                    - z_tr_dim)
                if(gamma_s /= 0.0) then
                  p_bar = p_t_b * (1.0 - gamma_s * (zk - z_tr_dim) / T_c_b) &
                      ** pow_s
                else
                  p_bar = p_t_b * exp(- (zk - z_tr_dim) / (Rsp * T_c_b / g))
                end if
              end if

              ! Define thetaStratTFC.
              thetaStratTFC(i, j, k) = T_bar * (press0_dim / p_bar) ** kappa &
                  / thetaRef
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
              0) * (thetaStratTFC(i, j, 1) - thetaStratTFC(i, j, 0)) / dz
          bvsStratTFC(i, j, 0) = bvsStratTFC(i, j, - 1)
          ! Between boundaries.
          do k = 1, nz
            bvsStratTFC(i, j, k) = g_ndim / thetaStratTFC(i, j, k) / jac(i, j, &
                k) * 0.5 * (thetaStratTFC(i, j, k + 1) - thetaStratTFC(i, j, k &
                - 1)) / dz
          end do
          ! Upper boundary.
          bvsStratTFC(i, j, nz + 1) = g_ndim / thetaStratTFC(i, j, nz + 1) &
              / jac(i, j, nz + 1) * (thetaStratTFC(i, j, nz + 1) &
              - thetaStratTFC(i, j, nz)) / dz
          bvsStratTFC(i, j, nz + 2) = bvsStratTFC(i, j, nz + 1)
        end do
      end do
      if(testCase == "smoothVortex") then
        bvsStratTFC = 0.0
      end if

    end if

  end subroutine update_topography

  !---------------------------------------------------------------------------

  ! TFC FJ
  ! Jacobian.

  function jac(i, j, k)

    real :: jac
    integer :: i, j, k

    jac = ((lz(1) - lz(0)) - topography_surface(i, j)) / (lz(1) - lz(0))

  end function jac

  ! TFC FJ
  ! Metric tensor.

  function met(i, j, k, mu, nu)

    real :: met
    integer :: i, j, k, mu, nu

    if((mu == 1 .and. nu == 3) .or. (mu == 3 .and. nu == 1)) then
      met = (topography_surface(i + 1, j) - topography_surface(i - 1, j)) &
          / (2.0 * dx) * (z(k) - (lz(1) - lz(0))) / ((lz(1) - lz(0)) &
          - topography_surface(i, j))
    else if((mu == 2 .and. nu == 3) .or. (mu == 3 .and. nu == 2)) then
      met = (topography_surface(i, j + 1) - topography_surface(i, j - 1)) &
          / (2.0 * dy) * (z(k) - (lz(1) - lz(0))) / ((lz(1) - lz(0)) &
          - topography_surface(i, j))
    else if(mu == 3 .and. nu == 3) then
      met = ((lz(1) - lz(0)) / ((lz(1) - lz(0)) - topography_surface(i, j))) &
          ** 2.0 + ((z(k) - (lz(1) - lz(0))) / ((lz(1) - lz(0)) &
          - topography_surface(i, j))) ** 2.0 * (((topography_surface(i + 1, &
          j) - topography_surface(i - 1, j)) / (2.0 * dx)) ** 2.0 &
          + ((topography_surface(i, j + 1) - topography_surface(i, j - 1)) &
          / (2.0 * dy)) ** 2.0)
    end if

  end function met

  ! TFC FJ
  ! Christophel tensor.

  function chris(i, j, k, mu, nu)

    real :: chris
    integer :: i, j, k, mu, nu

    if(mu == 1 .and. nu == 1) then
      chris = - (topography_surface(i - 1, j) - 2.0 * topography_surface(i, j) &
          + topography_surface(i + 1, j)) / (dx ** 2.0) * (z(k) - (lz(1) &
          - lz(0))) / ((lz(1) - lz(0)) - topography_surface(i, j))
    else if((mu == 1 .and. nu == 2) .or. (mu == 2 .and. nu == 1)) then
      chris = - (topography_surface(i + 1, j + 1) - topography_surface(i - 1, &
          j + 1) - topography_surface(i + 1, j - 1) + topography_surface(i &
          - 1, j - 1)) / (4.0 * dx * dy) * (z(k) - (lz(1) - lz(0))) / ((lz(1) &
          - lz(0)) - topography_surface(i, j))
    else if(mu == 2 .and. nu == 2) then
      chris = - (topography_surface(i, j - 1) - 2.0 * topography_surface(i, j) &
          + topography_surface(i, j + 1)) / (dy ** 2.0) * (z(k) - (lz(1) &
          - lz(0))) / ((lz(1) - lz(0)) - topography_surface(i, j))
    else if((mu == 1 .and. nu == 3) .or. (mu == 3 .and. nu == 1)) then
      chris = - (topography_surface(i + 1, j) - topography_surface(i - 1, j)) &
          / (2.0 * dx) / ((lz(1) - lz(0)) - topography_surface(i, j))
    else if((mu == 2 .and. nu == 3) .or. (mu == 3 .and. nu == 2)) then
      chris = - (topography_surface(i, j + 1) - topography_surface(i, j - 1)) &
          / (2.0 * dy) / ((lz(1) - lz(0)) - topography_surface(i, j))
    end if

  end function chris

  ! TFC FJ
  ! Transformation of vertical coordinate.

  function heightTFC(i, j, k)

    real :: heightTFC
    integer :: i, j, k

    heightTFC = z(k) * ((lz(1) - lz(0)) - topography_surface(i, j)) / (lz(1) &
        - lz(0)) + topography_surface(i, j)

  end function heightTFC

  ! TFC FJ
  ! Transformation of the vertical wind.

  function vertWindTFC(i, j, k, var)

    real, dimension((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        + nbz), nVar) :: var
    integer :: i, j, k

    real :: vertWindTFC
    real :: uEdgeR, uUEdgeR, uEdgeL, uUEdgeL, vEdgeF, vUEdgeF, vEdgeB, &
        vUEdgeB, wEdgeU

    uEdgeR = var(i, j, k, 2)
    uUEdgeR = var(i, j, k + 1, 2)
    uEdgeL = var(i - 1, j, k, 2)
    uUEdgeL = var(i - 1, j, k + 1, 2)
    vEdgeF = var(i, j, k, 3)
    vUEdgeF = var(i, j, k + 1, 3)
    vEdgeB = var(i, j - 1, k, 3)
    vUEdgeB = var(i, j - 1, k + 1, 3)
    wEdgeU = var(i, j, k, 4)

    vertWindTFC = trafoTFC(i, j, k, uEdgeR, uUEdgeR, uEdgeL, uUEdgeL, vEdgeF, &
        vUEdgeF, vEdgeB, vUEdgeB, wEdgeU, "car")

  end function vertWindTFC

  function trafoTFC(i, j, k, uEdgeR, uUEdgeR, uEdgeL, uUEdgeL, vEdgeF, &
      vUEdgeF, vEdgeB, vUEdgeB, wEdgeU, wind)

    integer :: i, j, k
    real :: uEdgeR, uUEdgeR, uEdgeL, uUEdgeL, vEdgeF, vUEdgeF, vEdgeB, &
        vUEdgeB, wEdgeU
    character(len = 3) :: wind

    real :: trafoTFC
    real :: jacEdgeU
    real :: uC, uU, vC, vU
    real :: metEdgeR, metUEdgeR, metEdgeL, metUEdgeL, metEdgeF, metUEdgeF, &
        metEdgeB, metUEdgeB
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
            k + 1) * met(i, j, k + 1, 1, 3) * uU) - 0.5 * (jac(i, j, k) &
            * met(i, j, k, 2, 3) * vC + jac(i, j, k + 1) * met(i, j, k + 1, 2, &
            3) * vU) + jacEdgeU * wEdgeU
      else if(wind == "tfc") then
        trafoTFC = (0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * uC + jac(i, j, &
            k + 1) * met(i, j, k + 1, 1, 3) * uU) + 0.5 * (jac(i, j, k) &
            * met(i, j, k, 2, 3) * vC + jac(i, j, k + 1) * met(i, j, k + 1, 2, &
            3) * vU) + wEdgeU) / jacEdgeU
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
            + 1, 1, 3) * uU) - 0.5 * (met(i, j, k, 2, 3) * vC + met(i, j, k &
            + 1, 2, 3) * vU) + wEdgeU)
      else if(wind == "tfc") then
        trafoTFC = 0.5 * (met(i, j, k, 1, 3) * uC + met(i, j, k + 1, 1, 3) &
            * uU) + 0.5 * (met(i, j, k, 2, 3) * vC + met(i, j, k + 1, 2, 3) &
            * vU) + wEdgeU / jacEdgeU
      end if

    case(3)

      ! Multiplication on u grid

      jacEdgeU = 0.5 * (jac(i, j, k) + jac(i, j, k + 1))
      metEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) + jac(i + 1, j, k) &
          * met(i + 1, j, k, 1, 3))
      metUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) + jac(i &
          + 1, j, k + 1) * met(i + 1, j, k + 1, 1, 3))
      metEdgeL = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) + jac(i - 1, j, k) &
          * met(i - 1, j, k, 1, 3))
      metUEdgeL = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) + jac(i &
          - 1, j, k + 1) * met(i - 1, j, k + 1, 1, 3))
      metEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) + jac(i, j + 1, k) &
          * met(i, j + 1, k, 2, 3))
      metUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) + jac(i, j &
          + 1, k + 1) * met(i, j + 1, k + 1, 2, 3))
      metEdgeB = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) + jac(i, j - 1, k) &
          * met(i, j - 1, k, 2, 3))
      metUEdgeB = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) + jac(i, j &
          - 1, k + 1) * met(i, j - 1, k + 1, 2, 3))
      if(wind == "car") then
        trafoTFC = - 0.25 * (metEdgeR * uEdgeR + metUEdgeR * uUEdgeR &
            + metEdgeL * uEdgeL + metUEdgeL * uUEdgeL) - 0.25 * (metEdgeF &
            * vEdgeF + metUEdgeF * vUEdgeF + metEdgeB * vEdgeB + metUEdgeB &
            * vUEdgeB) + jacEdgeU * wEdgeU
      else if(wind == "tfc") then
        trafoTFC = (0.25 * (metEdgeR * uEdgeR + metUEdgeR * uUEdgeR + metEdgeL &
            * uEdgeL + metUEdgeL * uUEdgeL) + 0.25 * (metEdgeF * vEdgeF &
            + metUEdgeF * vUEdgeF + metEdgeB * vEdgeB + metUEdgeB * vUEdgeB) &
            + wEdgeU) / jacEdgeU
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
            * uUEdgeR + metEdgeL * uEdgeL + metUEdgeL * uUEdgeL) - 0.25 &
            * (metEdgeF * vEdgeF + metUEdgeF * vUEdgeF + metEdgeB * vEdgeB &
            + metUEdgeB * vUEdgeB) + wEdgeU)
      else if(wind == "tfc") then
        trafoTFC = 0.25 * (metEdgeR * uEdgeR + metUEdgeR * uUEdgeR + metEdgeL &
            * uEdgeL + metUEdgeL * uUEdgeL) + 0.25 * (metEdgeF * vEdgeF &
            + metUEdgeF * vUEdgeF + metEdgeB * vEdgeB + metUEdgeB * vUEdgeB) &
            + wEdgeU / jacEdgeU
      end if

    case(5)

      ! Multiplication on uw grid

      jacEdgeU = 0.5 * (jac(i, j, k) + jac(i, j, k + 1))
      metEdgeRU = 0.25 * (jac(i, j, k) * met(i, j, k, 1, 3) + jac(i + 1, j, k) &
          * met(i + 1, j, k, 1, 3) + jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
          + jac(i + 1, j, k + 1) * met(i + 1, j, k + 1, 1, 3))
      metEdgeLU = 0.25 * (jac(i, j, k) * met(i, j, k, 1, 3) + jac(i - 1, j, k) &
          * met(i - 1, j, k, 1, 3) + jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
          + jac(i - 1, j, k + 1) * met(i - 1, j, k + 1, 1, 3))
      metEdgeFU = 0.25 * (jac(i, j, k) * met(i, j, k, 2, 3) + jac(i, j + 1, k) &
          * met(i, j + 1, k, 2, 3) + jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
          + jac(i, j + 1, k + 1) * met(i, j + 1, k + 1, 2, 3))
      metEdgeBU = 0.25 * (jac(i, j, k) * met(i, j, k, 2, 3) + jac(i, j - 1, k) &
          * met(i, j - 1, k, 2, 3) + jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
          + jac(i, j - 1, k + 1) * met(i, j - 1, k + 1, 2, 3))
      uEdgeRU = 0.5 * (uEdgeR + uUEdgeR)
      uEdgeLU = 0.5 * (uEdgeL + uUEdgeL)
      vEdgeFU = 0.5 * (vEdgeF + vUEdgeF)
      vEdgeBU = 0.5 * (vEdgeB + vUEdgeB)
      if(wind == "car") then
        trafoTFC = - 0.5 * (metEdgeRU * uEdgeRU - metEdgeLU * uEdgeLU) - 0.5 &
            * (metEdgeFU * vEdgeFU - metEdgeBU * vEdgeBU) + jacEdgeU * wEdgeU
      else if(wind == "tfc") then
        trafoTFC = (0.5 * (metEdgeRU * uEdgeRU - metEdgeLU * uEdgeLU) + 0.5 &
            * (metEdgeFU * vEdgeFU - metEdgeBU * vEdgeBU) + wEdgeU) / jacEdgeU
      end if

    case(6)

      ! Multiplication on uw grid (inverted)

      jacEdgeU = 0.5 * (jac(i, j, k) + jac(i, j, k + 1))
      metEdgeRU = 0.25 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3) + met(i, &
          j, k + 1, 1, 3) + met(i + 1, j, k + 1, 1, 3))
      metEdgeLU = 0.25 * (met(i, j, k, 1, 3) + met(i - 1, j, k, 1, 3) + met(i, &
          j, k + 1, 1, 3) + met(i - 1, j, k + 1, 1, 3))
      metEdgeFU = 0.25 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3) + met(i, &
          j, k + 1, 2, 3) + met(i, j + 1, k + 1, 2, 3))
      metEdgeBU = 0.25 * (met(i, j, k, 2, 3) + met(i, j - 1, k, 2, 3) + met(i, &
          j, k + 1, 2, 3) + met(i, j - 1, k + 1, 2, 3))
      uEdgeRU = 0.5 * (uEdgeR + uUEdgeR)
      uEdgeLU = 0.5 * (uEdgeL + uUEdgeL)
      vEdgeFU = 0.5 * (vEdgeF + vUEdgeF)
      vEdgeBU = 0.5 * (vEdgeB + vUEdgeB)
      if(wind == "car") then
        trafoTFC = jacEdgeU * (- 0.5 * (metEdgeRU * uEdgeRU - metEdgeLU &
            * uEdgeLU) - 0.5 * (metEdgeFU * vEdgeFU - metEdgeBU * vEdgeBU) &
            + wEdgeU)
      else if(wind == "tfc") then
        trafoTFC = 0.5 * (metEdgeRU * uEdgeRU - metEdgeLU * uEdgeLU) + 0.5 &
            * (metEdgeFU * vEdgeFU - metEdgeBU * vEdgeBU) + wEdgeU / jacEdgeU
      end if

    end select

  end function trafoTFC

  ! TFC FJ
  ! Cartesian stress tensor.

  function stressTensTFC(i, j, k, mu, nu, var)

    real, dimension((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        + nbz), nVar) :: var
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
        * met(i, j, k + 1, 1, 3))
    met13EdgeD = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) + jac(i, j, k - 1) &
        * met(i, j, k - 1, 1, 3))
    met23EdgeU = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) + jac(i, j, k + 1) &
        * met(i, j, k + 1, 2, 3))
    met23EdgeD = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) + jac(i, j, k - 1) &
        * met(i, j, k - 1, 2, 3))
    uF = 0.5 * (var(i, j + 1, k, 2) + var(i - 1, j + 1, k, 2))
    uB = 0.5 * (var(i, j - 1, k, 2) + var(i - 1, j - 1, k, 2))
    uU = 0.5 * (var(i, j, k + 1, 2) + var(i - 1, j, k + 1, 2))
    uD = 0.5 * (var(i, j, k - 1, 2) + var(i - 1, j, k - 1, 2))
    vR = 0.5 * (var(i + 1, j, k, 3) + var(i + 1, j - 1, k, 3))
    vL = 0.5 * (var(i - 1, j, k, 3) + var(i - 1, j - 1, k, 3))
    vU = 0.5 * (var(i, j, k + 1, 3) + var(i, j - 1, k + 1, 3))
    vD = 0.5 * (var(i, j, k - 1, 3) + var(i, j - 1, k - 1, 3))
    wR = 0.5 * (vertWindTFC(i + 1, j, k, var) + vertWindTFC(i + 1, j, k - 1, &
        var))
    wL = 0.5 * (vertWindTFC(i - 1, j, k, var) + vertWindTFC(i - 1, j, k - 1, &
        var))
    wF = 0.5 * (vertWindTFC(i, j + 1, k, var) + vertWindTFC(i, j + 1, k - 1, &
        var))
    wB = 0.5 * (vertWindTFC(i, j - 1, k, var) + vertWindTFC(i, j - 1, k - 1, &
        var))

    if(mu == 1 .and. nu == 1) then
      stressTensTFC = (2.0 * (jacEdgeR * var(i, j, k, 2) - jacEdgeL * var(i &
          - 1, j, k, 2)) / dx + (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
          * uU - jac(i, j, k - 1) * met(i, j, k - 1, 1, 3) * uD) / dz - 2.0 &
          / 3.0 * ((jacEdgeR * var(i, j, k, 2) - jacEdgeL * var(i - 1, j, k, &
          2)) / dx + (jacEdgeF * var(i, j, k, 3) - jacEdgeB * var(i, j - 1, k, &
          3)) / dy + (jacEdgeU * var(i, j, k, 4) - jacEdgeD * var(i, j, k - 1, &
          4)) / dz)) / jac(i, j, k)
    else if((mu == 1 .and. nu == 2) .or. (mu == 2 .and. nu == 1)) then
      stressTensTFC = (0.5 * (jac(i, j + 1, k) * uF - jac(i, j - 1, k) * uB) &
          / dy + 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) * uU &
          - jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) * uD) / dz + 0.5 &
          * (jac(i + 1, j, k) * vR - jac(i - 1, j, k) * vL) / dx + 0.5 &
          * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) * vU - jac(i, j, k - 1) &
          * met(i, j, k - 1, 1, 3) * vD) / dz) / jac(i, j, k)
    else if((mu == 1 .and. nu == 3) .or. (mu == 3 .and. nu == 1)) then
      stressTensTFC = (0.5 * (uU - uD) / dz + 0.5 * (jac(i + 1, j, k) * wR &
          - jac(i - 1, j, k) * wL) / dx + (met13EdgeU * vertWindTFC(i, j, k, &
          var) - met13EdgeD * vertWindTFC(i, j, k - 1, var)) / dz) / jac(i, j, &
          k)
    else if(mu == 2 .and. nu == 2) then
      stressTensTFC = (2.0 * (jacEdgeF * var(i, j, k, 3) - jacEdgeB * var(i, j &
          - 1, k, 3)) / dy + (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) * vU &
          - jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) * vD) / dz - 2.0 / 3.0 &
          * ((jacEdgeR * var(i, j, k, 2) - jacEdgeL * var(i - 1, j, k, 2)) &
          / dx + (jacEdgeF * var(i, j, k, 3) - jacEdgeB * var(i, j - 1, k, 3)) &
          / dy + (jacEdgeU * var(i, j, k, 4) - jacEdgeD * var(i, j, k - 1, 4)) &
          / dz)) / jac(i, j, k)
    else if((mu == 2 .and. nu == 3) .or. (mu == 3 .and. nu == 2)) then
      stressTensTFC = (0.5 * (vU - vD) / dz + 0.5 * (jac(i, j + 1, k) * wF &
          - jac(i, j - 1, k) * wB) / dy + (met23EdgeU * vertWindTFC(i, j, k, &
          var) - met23EdgeD * vertWindTFC(i, j, k - 1, var)) / dz) / jac(i, j, &
          k)
    else if(mu == 3 .and. nu == 3) then
      stressTensTFC = (2.0 * (vertWindTFC(i, j, k, var) - vertWindTFC(i, j, k &
          - 1, var)) / dz - 2.0 / 3.0 * ((jacEdgeR * var(i, j, k, 2) &
          - jacEdgeL * var(i - 1, j, k, 2)) / dx + (jacEdgeF * var(i, j, k, 3) &
          - jacEdgeB * var(i, j - 1, k, 3)) / dy + (jacEdgeU * var(i, j, k, 4) &
          - jacEdgeD * var(i, j, k - 1, 4)) / dz)) / jac(i, j, k)
    end if

  end function stressTensTFC

end module atmosphere_module
