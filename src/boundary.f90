module boundary_module

  use type_module
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
    integer :: i, j, k

    real, dimension(- nby:ny + nby, - nbz:nz + nbz) :: uBound

    select case(option)

    case("var")

      if(updateMass) then
        do i = 1, nbx
          var%rho(nx + i, :, :) = var%rho(i, :, :)
          var%rho(- i + 1, :, :) = var%rho(nx - i + 1, :, :)
        end do

        do i = 1, nbx
          var%rhop(nx + i, :, :) = var%rhop(i, :, :)
          var%rhop(- i + 1, :, :) = var%rhop(nx - i + 1, :, :)
        end do
      end if

      if(predictMomentum) then
        var%u(0, :, :) = var%u(nx, :, :)
        do i = 1, nbx
          var%u(nx + i, :, :) = var%u(i, :, :)
          var%u(- i, :, :) = var%u(nx - i, :, :)
          var%v(nx + i, :, :) = var%v(i, :, :)
          var%v(- i + 1, :, :) = var%v(nx - i + 1, :, :)
          var%w(nx + i, :, :) = var%w(i, :, :)
          var%w(- i + 1, :, :) = var%w(nx - i + 1, :, :)
        end do

      end if

      if(correctMomentum) then
        var%pi(nx + 1, :, :) = var%pi(1, :, :)
        var%pi(0, :, :) = var%pi(nx, :, :)

      end if

    case("flux")

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
    integer :: i, j, k

    select case(option)

    case("var")

      if(updateMass) then
        do j = 1, nby
          var%rho(:, ny + j, :) = var%rho(:, j, :)
          var%rho(:, - j + 1, :) = var%rho(:, ny - j + 1, :)
        end do

        do j = 1, nby
          var%rhop(:, ny + j, :) = var%rhop(:, j, :)
          var%rhop(:, - j + 1, :) = var%rhop(:, ny - j + 1, :)
        end do
      end if

      if(predictMomentum) then
        var%v(:, 0, :) = var%v(:, ny, :)
        do j = 1, nby
          var%u(:, ny + j, :) = var%u(:, j, :)
          var%u(:, - j + 1, :) = var%u(:, ny - j + 1, :)
          var%v(:, ny + j, :) = var%v(:, j, :)
          var%v(:, - j, :) = var%v(:, ny - j, :)
          var%w(:, ny + j, :) = var%w(:, j, :)
          var%w(:, - j + 1, :) = var%w(:, ny - j + 1, :)
        end do
      end if

      if(correctMomentum) then
        var%pi(:, ny + 1, :) = var%pi(:, 1, :)
        var%pi(:, 0, :) = var%pi(:, ny, :)

      end if

    case("flux")

      return

    case default
      stop "setBoundary_y: unknown option."
    end select

  end subroutine setBoundary_y_periodic

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
    integer :: i, j, k

    select case(option)

    case("var")

      if(updateMass) then
        ! reflect at boundary with change of sign
        do k = 1, nbz
          var%rho(:, :, - k + 1) = - var%rho(:, :, k)
          var%rho(:, :, nz + k) = - var%rho(:, :, nz - k + 1)
        end do

        ! vertical boundary condition for density fluctuations
        do k = 1, nbz
          var%rhop(:, :, - k + 1) = - var%rhop(:, :, k)
          var%rhop(:, :, nz + k) = - var%rhop(:, :, nz - k + 1)
        end do
      end if

      if(predictMomentum) then
        ! Vertical wind: zero at the boundary -> reflect with change of sign.
        var%w(:, :, 0) = 0.0
        var%w(:, :, nz) = 0.0
        do k = 1, nbz
          var%w(:, :, - k) = - var%w(:, :, k)
          var%w(:, :, nz + k) = - var%w(:, :, nz - k)
        end do

        ! Horizontal wind: free-slip -> reflect without change of sign.
        do k = 1, nbz
          var%u(:, :, - k + 1) = var%u(:, :, k)
          var%u(:, :, nz + k) = var%u(:, :, nz - k + 1)
          var%v(:, :, - k + 1) = var%v(:, :, k)
          var%v(:, :, nz + k) = var%v(:, :, nz - k + 1)
        end do
      end if

      if(correctMomentum) then
        ! Pressure: set gradient at vertical boundary to 0.
        var%pi(:, :, 0) = var%pi(:, :, 1)
        var%pi(:, :, nz + 1) = var%pi(:, :, nz)
      end if

    case("flux")

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
      end if

    case default
      stop "setBoundary_z: unknown option."
    end select

  end subroutine setBoundary_z_solidWall

end module boundary_module
