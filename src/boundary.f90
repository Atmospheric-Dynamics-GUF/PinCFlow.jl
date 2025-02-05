module boundary_module

  use type_module
  use flux_module ! BC for rhoTilde, uTilde,...
  use atmosphere_module
  use mpi_module
  use mpi

  implicit none

  private

  !-------------------
  !  public routines
  !-------------------
  public :: setBoundary
  public :: setBoundary_x_periodic
  public :: setBoundary_y_periodic

  !--------------------
  !  private routines
  !--------------------
  private :: setBoundary_z_periodic
  private :: setBoundary_z_solidWall

  contains

  subroutine setBoundary(var, flux, option)
    !-------------------------------------------
    ! 1) set values in ghost cells according to
    !    boundary conditions
    ! 2) set fluxes explicitly to zero at walls
    !-------------------------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var
    type(flux_type), intent(inout) :: flux
    character(len = *), intent(in) :: option

    !------------------------------
    !          x-direction
    !------------------------------
    select case(xBoundary)

    case("periodic")

      if(idim > 1) then
        ! boundary conditions taken care of by setHalos
      else
        call setBoundary_x_periodic(var, flux, option)
      endif

    case default
      stop "setBoundary: unknown case xBoundary"
    end select

    !------------------------------
    !          y-direction
    !------------------------------
    select case(yBoundary)

    case("periodic")

      if(jdim > 1) then
        ! boundary conditions taken care of by setHalos
      else
        call setBoundary_y_periodic(var, flux, option)
      endif

    case default
      stop "setBoundary: unknown case yBoundary"
    end select

    !------------------------------
    !          z-direction
    !------------------------------
    select case(zBoundary)

    case("solid_wall")
      call setBoundary_z_solidWall(var, flux, option)

    case default
      stop "setBoundary: unknown case zBoundary"
    end select

  end subroutine setBoundary

  !-----------------------------------------------------------------------

  subroutine setBoundary_x_periodic(var, flux, option)
    !-------------------------------------------
    ! 1) set values in ghost cells according to
    !    boundary conditions
    ! 2) set fluxes explicitly to zero at walls
    !-------------------------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var
    type(flux_type), intent(inout) :: flux
    character(len = *), intent(in) :: option

    ! local variables
    integer :: i, j, k, iVar, ii

    real, dimension(- nby:ny + nby, - nbz:nz + nbz) :: uBound

    select case(option)

    case("var")
      !-----------------------------------
      !                var
      !-----------------------------------

      if(updateMass) then
        ! density -> iVar = 1
        do i = 1, nbx
          var%rho(nx + i, :, :) = var%rho(i, :, :)
          var%rho(- i + 1, :, :) = var%rho(nx - i + 1, :, :)
        end do

        if(verbose .and. master) then
          print *, "horizontalBoundary: x-horizontal BC for rho set."
        end if
      end if

      if(timeScheme == "semiimplicit" .or. auxil_equ) then
        ! density fluctuations -> iVar = 6
        do i = 1, nbx
          var%rhop(nx + i, :, :) = var%rhop(i, :, :)
          var%rhop(- i + 1, :, :) = var%rhop(nx - i + 1, :, :)
        end do

      end if

      if(predictMomentum) then
        ! velocity u (staggered along x) -> iVar = 2

        var%u(0, :, :) = var%u(nx, :, :)

        do i = 1, nbx
          ! velocity u (staggered along x) -> iVar = 2
          var%u(nx + i, :, :) = var%u(i, :, :)
          var%u(- i, :, :) = var%u(nx - i, :, :)
          ! velocity v -> iVar = 3                    ! ghost cells
          var%v(nx + i, :, :) = var%v(i, :, :) ! right
          var%v(- i + 1, :, :) = var%v(nx - i + 1, :, :) ! left
          ! velocity w -> iVar = 4
          var%w(nx + i, :, :) = var%w(i, :, :)
          var%w(- i + 1, :, :) = var%w(nx - i + 1, :, :)
        end do

      end if

      if(correctMomentum) then
        ! pressure Variable -> iVar = 5
        var%pi(nx + 1, :, :) = var%pi(1, :, :) ! right ghost cell
        var%pi(0, :, :) = var%pi(nx, :, :) ! left ghost cells

      end if

    case("varTilde")
      !-----------------------------------
      !             varTilde
      !-----------------------------------

      ! the following three boundary-condition calls can probably be
      ! removed
      ! probably only necessary for ALDM (that is not used any more
      ! anyway)

      if(updateMass) then
        ! reconstructed density needed in ghost cell i = nx+2
        rhoTilde(nx + 2, :, :, 1, 0) = rhoTilde(2, :, :, 1, 0)

        ! ...in ghost cell i = -1
        rhoTilde(- 1, :, :, 1, 1) = rhoTilde(nx - 1, :, :, 1, 1)

      end if

      if(timeScheme == "semiimplicit" .or. auxil_equ) then
        ! reconstructed density fluctuation needed in ghost cell i = nx+2
        rhopTilde(nx + 2, :, :, 1, 0) = rhopTilde(2, :, :, 1, 0)

        ! ...in ghost cell i = -1
        rhopTilde(- 1, :, :, 1, 1) = rhopTilde(nx - 1, :, :, 1, 1)

      end if


    case("flux")
      !-----------------------------------
      !              flux
      !-----------------------------------

      return

    case default
      stop "setBoundary_x: unknown option."
    end select

  end subroutine setBoundary_x_periodic

  !--------------------------------------------------------------------

  subroutine setBoundary_y_periodic(var, flux, option)
    !-------------------------------------------
    ! 1) set values in ghost cells according to
    !    boundary conditions
    ! 2) set fluxes explicitly to zero at walls
    !-------------------------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var
    type(flux_type), intent(inout) :: flux
    character(len = *), intent(in) :: option

    ! local variables
    integer :: i, j, k, iVar, ii

    real, dimension(- nbx:nx + nbx, - nbz:nz + nbz) :: vBound

    select case(option)

    case("var")
      !-----------------------------------
      !                var
      !-----------------------------------

      if(updateMass) then
        ! density -> iVar = 1
        do j = 1, nby
          var%rho(:, ny + j, :) = var%rho(:, j, :)
          var%rho(:, - j + 1, :) = var%rho(:, ny - j + 1, :)
        end do

      end if

      if(timeScheme == "semiimplicit" .or. auxil_equ) then
        ! density fluctuations -> iVar = 6
        do j = 1, nby
          var%rhop(:, ny + j, :) = var%rhop(:, j, :)
          var%rhop(:, - j + 1, :) = var%rhop(:, ny - j + 1, :)
        end do

      end if

      if(predictMomentum) then
        ! velocity v (staggared along y) -> iVar = 3

        var%v(:, 0, :) = var%v(:, ny, :)

        do j = 1, nby
          ! velocity u -> iVar = 2
          var%u(:, ny + j, :) = var%u(:, j, :)
          var%u(:, - j + 1, :) = var%u(:, ny - j + 1, :)
          ! velocity v (staggared along y) -> iVar = 3
          var%v(:, ny + j, :) = var%v(:, j, :)
          var%v(:, - j, :) = var%v(:, ny - j, :)
          ! velocity w -> iVar = 4
          var%w(:, ny + j, :) = var%w(:, j, :)
          var%w(:, - j + 1, :) = var%w(:, ny - j + 1, :)
        end do

      end if

      if(correctMomentum) then
        ! pressure variable -> iVar = 5
        var%pi(:, ny + 1, :) = var%pi(:, 1, :) ! forward ghost cell
        var%pi(:, 0, :) = var%pi(:, ny, :) ! backward

      end if

    case("varTilde")
      !-----------------------------------
      !              varTilde
      !-----------------------------------

      ! the following three boundary-condition calls can probably be
      ! removed
      ! probably only necessary for ALDM (that is not used any more
      ! anyway)

      if(updateMass) then
        ! reconstructed density needed in ghost cell j = ny+2
        rhoTilde(:, ny + 2, :, 2, 0) = rhoTilde(:, 2, :, 2, 0)

        ! ...in ghost cell j = -1
        rhoTilde(:, - 1, :, 2, 1) = rhoTilde(:, ny - 1, :, 2, 1)

      end if

      if(timeScheme == "semiimplicit" .or. auxil_equ) then
        ! reconstructed density fluctuations needed in ghost cell
        ! j = ny+2
        rhopTilde(:, ny + 2, :, 2, 0) = rhopTilde(:, 2, :, 2, 0)

        ! ...in ghost cell j = -1
        rhopTilde(:, - 1, :, 2, 1) = rhopTilde(:, ny - 1, :, 2, 1)

      end if

    case("flux")
      !-----------------------------------
      !              flux
      !-----------------------------------
      return

    case default
      stop "setBoundary_y: unknown option."
    end select

  end subroutine setBoundary_y_periodic

  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------

  subroutine setBoundary_z_solidWall(var, flux, option)
    !-------------------------------------------
    ! 1) set values in ghost cells according to
    !    boundary conditions
    ! 2) set fluxes explicitly to zero at walls
    !-------------------------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var
    type(flux_type), intent(inout) :: flux
    character(len = *), intent(in) :: option

    ! local variables
    integer :: i, j, k, iVar, ii
    real, dimension(- nby:ny + nby, - nbz:nz + nbz) :: uBound
    real, dimension(- nbx:nx + nbx, - nbz:nz + nbz) :: vBound

    ! rho flux correction
    real :: rhoU, rhoD, hRho, wSurf

    ! uvw flux correction
    real :: rhoEdge
    real :: hRhoU, hRhoV, hRhoW
    real :: uD, uU, wL, wR
    real :: vD, vU, wB, wF, wD, wU

    ! theta flux correction
    real :: thetaU, thetaD, hTheta

    ! uvw flux correction
    real :: thetaEdge
    real :: hThetaU, hThetaV, hThetaW

    select case(option)

    case("var")
      !-----------------------------------
      !                var
      !-----------------------------------

      if(updateMass) then
        ! reflect at boundary with change of sign
        ! rho -> iVar = 1

        ! in case of baroclinic life-cycle simulation only the
        ! deviations from the environmental state are reflected

        do k = 1, nbz

          var%rho(:, :, - k + 1) = - var%rho(:, :, k)
          var%rho(:, :, nz + k) = - var%rho(:, :, nz - k + 1)

        end do

        if(timeScheme == "semiimplicit" .or. auxil_equ) then
          ! vertical boundary condition for density fluctuations

          do k = 1, nbz
              var%rhop(:, :, - k + 1) = - var%rhop(:, :, k)
              var%rhop(:, :, nz + k) = - var%rhop(:, :, nz - k + 1)
            end do
          end if
      end if


      if(predictMomentum) then
        ! w -> set to zero at bound,
        !      reflect at bound with change of sign

        var%w(:, :, 0) = 0.0
        var%w(:, :, nz) = 0.0

        do k = 1, nbz
          var%w(:, :, - k) = - var%w(:, :, k)
          var%w(:, :, nz + k) = - var%w(:, :, nz - k)
        end do

        ! transverse velocities u,v:
        ! -> general case: reflect at bound. and change sign (no slip)
        !    life-cycle simulation: only deviations from environmental
        !                           state reflected,
        !                           sign not changed (free slip)
        ! FS changed sign back to no slip
          do k = 1, nbz
            var%u(:, :, - k + 1) = var%u(:, :, k)
            var%u(:, :, nz + k) = var%u(:, :, nz - k + 1)
            var%v(:, :, - k + 1) = var%v(:, :, k)
            var%v(:, :, nz + k) = var%v(:, :, nz - k + 1)
          end do
          ! No-slip condition.
          ! do k = 1, nbz
          !     var(:, :, - k + 1, 2) = - var(:, :, k, 2)
          !     var(:, :, nz + k, 2) = - var(:, :, nz - k + 1, 2)
          !     var(:, :, - k + 1, 3) = - var(:, :, k, 3)
          !     var(:, :, nz + k, 3) = - var(:, :, nz - k + 1, 3)
          ! end do
      end if

      if(correctMomentum) then
        ! set gradient at vertical boundary to 0
        ! life-cycle simulation: only gradient of deviations from
        !                        environmental pressure are set to 0
          var%pi(:, :, 0) = var%pi(:, :, 1)

          ! z = zMax
          var%pi(:, :, nz + 1) = var%pi(:, :, nz)
      end if

      ! additional call to setHalos in case baroclinic_LC in order to set
      ! the overlap betwen vertical and horizontal halos right

      !---------------------------------------------------------------

    case("varTilde ")
      !-----------------------------------
      !             varTilde
      !-----------------------------------

      return

    case("flux")
      !-----------------------------------
      !                flux
      !-----------------------------------

      if(updateMass) then
        ! set vertical fluxes at wall to 0

        ! these settings correspond to solid-wall condition (w = 0) at
        ! the bottom and top, whence advective and turbulent mass
        ! fluxes must vanish there
        ! a different b.c. seems appropriate if molecular diffusion is
        ! taken into account!

        ! density
        flux%rho(:, :, 0, 3) = 0.0
        flux%rho(:, :, nz, 3) = 0.0

        ! diffusive potential-temperature fluxes
        flux%theta(:, :, 0, 3) = 0.0
        flux%theta(:, :, nz, 3) = 0.0

        ! replace flux by CDS fluxes at upper / lower region
        if(rhoFluxCorr) stop 'ERROR: rhoFluxCorr = .false. expected'
      end if ! updateMass

      !if( updateIce ) then
      ! set vertical fluxes at wall to 0
      ! analog to mass flux modification above

      if(timeScheme == "semiimplicit" .or. auxil_equ) then
        ! density fluctuations
        flux%rhop(:, :, 0, 3) = 0.0
        flux%rhop(:, :, nz, 3) = 0.0
      end if

      if(predictMomentum) then
        ! momentum rho*u
        flux%u(:, :, 0, 3) = 0.0
        flux%u(:, :, nz, 3) = 0.0

        ! momentum rho*v
        flux%v(:, :, 0, 3) = 0.0
        flux%v(:, :, nz, 3) = 0.0

        ! momentum rho*w
        flux%w(:, :, - 1, 3) = 0.0
        flux%w(:, :, nz, 3) = 0.0

        ! replace flux by CDS fluxes at upper / lower region
      end if ! predictMomentum

    case default
      stop "setBoundary_z: unknown option."
    end select

  end subroutine setBoundary_z_solidWall

end module boundary_module
