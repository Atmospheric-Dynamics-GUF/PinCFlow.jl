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
  public :: jac, met, vertWindTFC, trafoTFC, stressTensTFC

  public :: update_topography

  real, dimension(:), allocatable :: PStrat, rhoStrat, thetaStrat, bvsStrat
  real, dimension(:), allocatable :: rhoStrat_d, rhoStrat_s

  real, dimension(:), allocatable :: pistrat
  real, dimension(:), allocatable :: PStratTilde, rhoStratTilde, thetaStratTilde

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

  ! stable atmosphere
  real :: N2 ! scaled square of Brunt-Vaisala frequency
  real :: NN ! scaled of Brunt-Vaisala frequency
  real :: coeff ! long coefficient

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
    end if

    !----------------------------------
    !       auxiliary quantities
    !----------------------------------

    gamma = 1.4
    gamma_1 = gamma - 1.0
    kappa = (gamma - 1.0) / gamma ! = R/c_p = Gamma at Klein 2008
    kappaInv = 1.0 / kappa
    gammaInv = 1.0 / gamma

    !----------------------------------
    !      reference quantities &
    !      nondimensional numbers
    !----------------------------------

    select case(referenceQuantities)

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

    if(f_Coriolis_dim /= 0.0) then
      Ro(0:ny + 1) = uRef / f_Coriolis_dim / lRef
      RoInv(0:ny + 1) = 1.0 / Ro(:)
    else
      Ro(0:ny + 1) = 1.d40
      RoInv(0:ny + 1) = 0.0
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

    case("pseudo_incompressible")

      select case(background)

      case('isothermal')

        
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
                  pStratTFC(i, j, k) = p0 * exp(- sig * zTFC(i, j, k) &
                      &/ gamma / T0)
                  ! Define thetaStratTFC.
                  thetaStratTFC(i, j, k) = T0 * exp(kappa * sig / T0 &
                      &* zTFC(i, j, k))
                  ! Define rhoStratTFC.
                  rhoStratTFC(i, j, k) = pStratTFC(i, j, k) / thetaStratTFC(i, &
                      &j, k)
                end do
              end do
            end do
          end if

        ! GBcorr
        bvsStrat = N2

  
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
    integer :: ix, jy, kz
    integer :: iwm

    if(.not. topography) return

    if(lz(0) /= 0.0) stop "Error in setup_topography: lz(0) must be zero for &
        &TFC!"

    mountainHeight = mountainHeight_dim / lRef
    mountainWidth = mountainWidth_dim / lRef
    mountainWavenumber = pi / mountainWidth

    x_center = 0.5 * (lx(1) + lx(0))
    y_center = 0.5 * (ly(1) + ly(0))

    ix0 = is + nbx - 1
    jy0 = js + nby - 1

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
                  &** 2.0))) * (1.0 + envelope_reduction)
              do iwm = 0, spectral_modes - 1
                kk = mountainWavenumber * cos(pi / spectral_modes * iwm)
                ll = mountainWavenumber * sin(pi / spectral_modes * iwm)
                topography_surface(ix, jy) = topography_surface(ix, jy) &
                    &+ 0.25 * mountainHeight * (1.0 + cos(mountainWavenumber &
                    &/ range_factor * sqrt((x(ix + ix0) - x_center) ** 2.0 &
                    &+ (y(jy + jy0) - y_center) ** 2.0))) * cos(kk * (x(ix &
                    &+ ix0) - x_center) + ll * (y(jy + jy0) - y_center)) &
                    &/ spectral_modes * (1.0 - envelope_reduction)
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

    ! Compute the stretched vertical grid.
    do kz = - nbz, nz + nbz
      zTildeS(kz) = map(z(kz) + 0.5 * dz)
    end do
    do kz = - nbz + 1, nz + nbz
      zS(kz) = 0.5 * (zTildeS(kz) + zTildeS(kz - 1))
    end do
    zS(- nbz) = zTildeS(- nbz) - 0.5 * (zTildeS(nbz + 1) - zTildeS(nbz))

    ! Compute the physical layers.
    do kz = - nbz, nz + nbz
      zTildeTFC(:, :, kz) = (lz(1) - topography_surface) / lz(1) * zTildeS(kz) &
          &+ topography_surface
      zTFC(:, :, kz) = (lz(1) - topography_surface) / lz(1) * zS(kz) &
          &+ topography_surface
    end do

  end subroutine setup_topography

  !---------------------------------------------------------------------------

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

  function jac(i, j, k)
    ! Jacobian.

    real :: jac
    integer :: i, j, k

    jac = (lz(1) - topography_surface(i, j)) / lz(1) * (zTildeS(k) - zTildeS(k &
        &- 1)) / dz

  end function jac

  function met(i, j, k, mu, nu)
    ! Metric tensor.

    real :: met
    integer :: i, j, k, mu, nu

    if((mu == 1 .and. nu == 3) .or. (mu == 3 .and. nu == 1)) then
      met = (topography_surface(i + 1, j) - topography_surface(i - 1, j)) &
          &/ (2.0 * dx) * (zS(k) - lz(1)) / (lz(1) - topography_surface(i, &
          &j)) * dz / (zTildeS(k) - zTildeS(k - 1))
    else if((mu == 2 .and. nu == 3) .or. (mu == 3 .and. nu == 2)) then
      met = (topography_surface(i, j + 1) - topography_surface(i, j - 1)) &
          &/ (2.0 * dy) * (zS(k) - lz(1)) / (lz(1) - topography_surface(i, &
          &j)) * dz / (zTildeS(k) - zTildeS(k - 1))
    else if(mu == 3 .and. nu == 3) then
      met = ((lz(1) / (lz(1) - topography_surface(i, j))) ** 2.0 + ((zS(k) &
          &- lz(1)) / (lz(1) - topography_surface(i, j))) ** 2.0 &
          &* (((topography_surface(i + 1, j) - topography_surface(i - 1, j)) &
          &/ (2.0 * dx)) ** 2.0 + ((topography_surface(i, j + 1) &
          &- topography_surface(i, j - 1)) / (2.0 * dy)) ** 2.0)) * (dz &
          &/ (zTildeS(k) - zTildeS(k - 1))) ** 2.0
    end if

  end function met

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
