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

    case("periodic")
      call setBoundary_z_periodic(var, flux, option)

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

        if(verbose .and. master) then
          print *, "horizontalBoundary: x-horizontal BC for rhop set."
        end if
      end if

      if(model == "compressible") then
        ! set boundary for P mass weighted pot. temp.
        do i = 1, nbx
          var%P(nx + i, :, :) = var%P(i, :, :)
          var%P(- i + 1, :, :) = var%P(nx - i + 1, :, :)
        end do

        if(verbose .and. master) then
          print *, "horizontalBoundary: x-horizontal BC for P set."
        end if
      end if

      ! set boundaries if including tracer
      if(updateTracer) then
        do i = 1, nbx
          var%chi(nx + i, :, :) = var%chi(i, :, :)
          var%chi(- i + 1, :, :) = var%chi(nx - i + 1, :, :)
        end do

        if(verbose .and. master) then
          print *, "horizontalBoundary: x-horizontal BC for tracer set."
        end if
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

        if(verbose .and. master) then
          print *, "boundary.f90/horizontalBoundary:  x-horizontal BC for u, &
              &v, w set."
        end if
      end if

      if(correctMomentum) then
        ! pressure Variable -> iVar = 5
        var%pi(nx + 1, :, :) = var%pi(1, :, :) ! right ghost cell
        var%pi(0, :, :) = var%pi(nx, :, :) ! left ghost cells

        if(verbose .and. master) then
          print *, "boundary.f90/horizontalBoundary:  x-horizontal BC for p &
              &set."
        end if
      end if

    case("ice")
      ! ice variables
      do iVar = 1, nVarIce
        do i = 1, nbx
          var%ICE(nx + i, :, :, iVar) = var%ICE(i, :, :, iVar)
          var%ICE(- i + 1, :, :, iVar) = var%ICE(nx - i + 1, :, :, iVar)
        end do
      end do

      if(verbose .and. master) print *, "horizontalBoundary:  x-horizontal BC &
          &for ice variables set."

    case("iceTilde")
      ! reconstructed density needed in ghost cell i = nx+2
      nAerTilde(nx + 2, :, :, 1, 0) = nAerTilde(2, :, :, 1, 0)
      nIceTilde(nx + 2, :, :, 1, 0) = nIceTilde(2, :, :, 1, 0)
      qIceTilde(nx + 2, :, :, 1, 0) = qIceTilde(2, :, :, 1, 0)
      qvTilde(nx + 2, :, :, 1, 0) = qvTilde(2, :, :, 1, 0)

      ! ...in ghost cell i = -1
      nAerTilde(- 1, :, :, 1, 1) = nAerTilde(nx - 1, :, :, 1, 1)
      nIceTilde(- 1, :, :, 1, 1) = nIceTilde(nx - 1, :, :, 1, 1)
      qIceTilde(- 1, :, :, 1, 1) = qIceTilde(nx - 1, :, :, 1, 1)
      qvTilde(- 1, :, :, 1, 1) = qvTilde(nx - 1, :, :, 1, 1)

      if(verbose .and. master) print *, "horizontalBoundary:  x-horizontal BC &
          &for iceTilde variables set."

    case("iceFlux")
      ! nothing

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

        if(verbose .and. master) then
          print *, "horizontalBoundary: x-horizontal BC for rhoTilde set."
        end if
      end if

      if(timeScheme == "semiimplicit" .or. auxil_equ) then
        ! reconstructed density fluctuation needed in ghost cell i = nx+2
        rhopTilde(nx + 2, :, :, 1, 0) = rhopTilde(2, :, :, 1, 0)

        ! ...in ghost cell i = -1
        rhopTilde(- 1, :, :, 1, 1) = rhopTilde(nx - 1, :, :, 1, 1)

        if(verbose .and. master) then
          print *, "horizontalBoundary: x-horizontal BC for rhoTilde set."
        end if
      end if

      if(updateIce) then
        ! reconstructed density needed in ghost cell i = nx+2
        nAerTilde(nx + 2, :, :, 1, 0) = nAerTilde(2, :, :, 1, 0)
        nIceTilde(nx + 2, :, :, 1, 0) = nIceTilde(2, :, :, 1, 0)
        qIceTilde(nx + 2, :, :, 1, 0) = qIceTilde(2, :, :, 1, 0)
        qvTilde(nx + 2, :, :, 1, 0) = qvTilde(2, :, :, 1, 0)

        ! ...in ghost cell i = -1
        nAerTilde(- 1, :, :, 1, 1) = nAerTilde(nx - 1, :, :, 1, 1)
        nIceTilde(- 1, :, :, 1, 1) = nIceTilde(nx - 1, :, :, 1, 1)
        qIceTilde(- 1, :, :, 1, 1) = qIceTilde(nx - 1, :, :, 1, 1)
        qvTilde(- 1, :, :, 1, 1) = qvTilde(nx - 1, :, :, 1, 1)

        if(verbose .and. master) print *, "horizontalBoundary:  x-horizontal &
            &BC for iceTilde variables set."

      end if

      ! set boundaries if including tracer
      if(updateTracer) then
        tracerTilde(nx + 2, :, :, 1, 0) = tracerTilde(2, :, :, 1, 0)
        tracerTilde(- 1, :, :, 1, 1) = tracerTilde(nx - 1, :, :, 1, 1)

        if(verbose .and. master) then
          print *, "horizontalBoundary: x-horizontal BC for tracerTilde set."
        end if
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

        if(verbose .and. master) then
          print *, "horizontalBoundary: y-horizontal BC for rho set."
        end if
      end if

      if(timeScheme == "semiimplicit" .or. auxil_equ) then
        ! density fluctuations -> iVar = 6
        do j = 1, nby
          var%rhop(:, ny + j, :) = var%rhop(:, j, :)
          var%rhop(:, - j + 1, :) = var%rhop(:, ny - j + 1, :)
        end do

        if(verbose .and. master) then
          print *, "horizontalBoundary: y-horizontal BC for rho set."
        end if
      end if

      if(model == "compressible") then
        ! set boundary for P mass weighted pot. temp.
        do j = 1, nby
          var%P(:, ny + j, :) = var%P(:, j, :)
          var%P(:, - j + 1, :) = var%P(:, ny - j + 1, :)
        end do

        if(verbose .and. master) then
          print *, "horizontalBoundary: y-horizontal BC for P set."
        end if
      end if

      ! set boundaries if including tracer
      if(updateTracer) then
        do j = 1, nby
          var%chi(:, ny + j, :) = var%chi(:, j, :)
          var%chi(:, - j + 1, :) = var%chi(:, ny - j + 1, :)
        end do

        if(verbose .and. master) then
          print *, "horizontalBoundary: y-horizontal BC for tracer set."
        end if
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

        if(verbose .and. master) then
          print *, "boundary.f90/horizontalBoundary:  y-horizontal BC for u, &
              &v, w set."
        end if
      end if

      if(correctMomentum) then
        ! pressure variable -> iVar = 5
        var%pi(:, ny + 1, :) = var%pi(:, 1, :) ! forward ghost cell
        var%pi(:, 0, :) = var%pi(:, ny, :) ! backward

        if(verbose .and. master) then
          print *, "boundary.f90/horizontalBoundary:  y-horizontal BC for p &
              &set."
        end if
      end if

    case("ice")
      ! ice variables
      do iVar = 1, nVarIce
        do j = 1, nby
          var%ICE(:, ny + j, :, iVar) = var%ICE(:, j, :, iVar)
          var%ICE(:, - j + 1, :, iVar) = var%ICE(:, ny - j + 1, :, iVar)
        end do
      end do

      if(verbose .and. master) print *, "horizontalBoundary:  y-horizontal BC &
          &for ice variables set."

    case("iceTilde")
      ! see above, similar to rho
      nAerTilde(:, ny + 2, :, 2, 0) = nAerTilde(:, 2, :, 2, 0)
      nIceTilde(:, ny + 2, :, 2, 0) = nIceTilde(:, 2, :, 2, 0)
      qIceTilde(:, ny + 2, :, 2, 0) = qIceTilde(:, 2, :, 2, 0)
      qvTilde(:, ny + 2, :, 2, 0) = qvTilde(:, 2, :, 2, 0)

      ! see above, similar to rho
      nAerTilde(:, - 1, :, 2, 1) = nAerTilde(:, ny - 1, :, 2, 1)
      nIceTilde(:, - 1, :, 2, 1) = nIceTilde(:, ny - 1, :, 2, 1)
      qIceTilde(:, - 1, :, 2, 1) = qIceTilde(:, ny - 1, :, 2, 1)
      qvTilde(:, - 1, :, 2, 1) = qvTilde(:, ny - 1, :, 2, 1)

      if(verbose .and. master) print *, "horizontalBoundary:  y-horizontal BC &
          &for iceTilde variables set."

    case("iceFlux")
      ! nothing

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

        if(verbose .and. master) then
          print *, "horizontalBoundary: y-horizontal BC for rhoTilde set."
        end if
      end if

      if(timeScheme == "semiimplicit" .or. auxil_equ) then
        ! reconstructed density fluctuations needed in ghost cell
        ! j = ny+2
        rhopTilde(:, ny + 2, :, 2, 0) = rhopTilde(:, 2, :, 2, 0)

        ! ...in ghost cell j = -1
        rhopTilde(:, - 1, :, 2, 1) = rhopTilde(:, ny - 1, :, 2, 1)

        if(verbose .and. master) then
          print *, "horizontalBoundary: y-horizontal BC for rhoTilde set."
        end if
      end if

      if(updateIce) then
        ! see above, similar to rho
        nAerTilde(:, ny + 2, :, 2, 0) = nAerTilde(:, 2, :, 2, 0)
        nIceTilde(:, ny + 2, :, 2, 0) = nIceTilde(:, 2, :, 2, 0)
        qIceTilde(:, ny + 2, :, 2, 0) = qIceTilde(:, 2, :, 2, 0)
        qvTilde(:, ny + 2, :, 2, 0) = qvTilde(:, 2, :, 2, 0)

        ! see above, similar to rho
        nAerTilde(:, - 1, :, 2, 1) = nAerTilde(:, ny - 1, :, 2, 1)
        nIceTilde(:, - 1, :, 2, 1) = nIceTilde(:, ny - 1, :, 2, 1)
        qIceTilde(:, - 1, :, 2, 1) = qIceTilde(:, ny - 1, :, 2, 1)
        qvTilde(:, - 1, :, 2, 1) = qvTilde(:, ny - 1, :, 2, 1)

        if(verbose .and. master) print *, "horizontalBoundary:  y-horizontal &
            &BC for iceTilde variables set."

      end if

      ! set boundaries if including tracer
      if(updateTracer) then
        tracerTilde(:, ny + 2, :, 2, 0) = tracerTilde(:, 2, :, 2, 0)
        tracerTilde(:, - 1, :, 2, 1) = tracerTilde(:, ny - 1, :, 2, 1)

        if(verbose .and. master) then
          print *, "horizontalBoundary: y-horizontal BC for tracerTilde set."
        end if
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

  subroutine setBoundary_z_periodic(var, flux, option)
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

    real, dimension(- nbx:nx + nbx, - nby:ny + nby) :: wBound

    select case(option)

    case("var")
      !-----------------------------------
      !                var
      !-----------------------------------

      if(updateMass) then
        ! density -> iVar = 1
        do k = 1, nbz
          var%rho(:, :, nz + k) = var%rho(:, :, k)
          var%rho(:, :, - k + 1) = var%rho(:, :, nz - k + 1)
        end do

        if(verbose .and. master) then
          print *, "setBoundary_z_periodic: z-periodic BC for rho set."
        end if
      end if

      if(timeScheme == "semiimplicit" .or. auxil_equ) then
        ! density fluctuations -> iVar = 6
        do k = 1, nbz
          var%rhop(:, :, nz + k) = var%rhop(:, :, k)
          var%rhop(:, :, - k + 1) = var%rhop(:, :, nz - k + 1)
        end do

        if(verbose .and. master) then
          print *, "setBoundary_z_periodic: z-periodic BC for rhop set."
        end if
      end if

      ! set boundaries if including tracer
      if(updateTracer) then
        do k = 1, nbz
          var%chi(:, :, nz + k) = - var%chi(:, :, k)
          var%chi(:, :, - k + 1) = - var%chi(:, :, nz - k + 1)
        end do

        if(verbose .and. master) then
          print *, "setBoundary_z_periodic: z-periodic BC for tracer set."
        end if
      end if

      if(predictMomentum) then
        ! velocity w (staggared along z) -> iVar = 4

        var%w(:, :, 0) = var%w(:, :, nz)

        do k = 1, nbz
          ! velocity u -> iVar = 2
          var%u(:, :, nz + k) = var%u(:, :, k)
          var%u(:, :, - k + 1) = var%u(:, :, nz - k + 1)
          ! velocity v -> iVar = 3
          var%v(:, :, nz + k) = var%v(:, :, k)
          var%v(:, :, - k + 1) = var%v(:, :, nz - k + 1)
          ! velocity w (staggared along z) -> iVar = 4
          var%w(:, :, nz + k) = var%w(:, :, k)
          var%w(:, :, - k) = var%w(:, :, nz - k)
        end do

        if(verbose .and. master) then
          print *, "setBoundary_z_periodic:  z-periodic BC for u, v, w set."
        end if
      end if

      if(correctMomentum) then
        ! pressure variable -> iVar = 5
        var%pi(:, :, nz + 1) = var%pi(:, :, 1) ! forward ghost cell
        var%pi(:, :, 0) = var%pi(:, :, nz) ! backward

        if(verbose .and. master) then
          print *, "setBoundary_z_periodic:  z-periodic BC for p set."
        end if
      end if

    case("ice")
      ! ice variables
      do iVar = 1, nVarIce
        do k = 1, nbz
          var%ICE(:, :, nz + k, iVar) = var%ICE(:, :, k, iVar)
          var%ICE(:, :, - k + 1, iVar) = var%ICE(:, :, nz - k + 1, iVar)
        end do
      end do

      if(verbose .and. master) print *, "setBoundary_z_periodic:  z-periodic &
          &BC for ice variables set."

    case("iceTilde")
      ! see above, similar to rho
      nAerTilde(:, :, nz + 2, 3, 0) = nAerTilde(:, :, 2, 3, 0)
      nIceTilde(:, :, nz + 2, 3, 0) = nIceTilde(:, :, 2, 3, 0)
      qIceTilde(:, :, nz + 2, 3, 0) = qIceTilde(:, :, 2, 3, 0)
      qvTilde(:, :, nz + 2, 3, 0) = qvTilde(:, :, 2, 3, 0)

      ! see above, similar to rho
      nAerTilde(:, :, - 1, 3, 1) = nAerTilde(:, :, nz - 1, 3, 1)
      nIceTilde(:, :, - 1, 3, 1) = nIceTilde(:, :, nz - 1, 3, 1)
      qIceTilde(:, :, - 1, 3, 1) = qIceTilde(:, :, nz - 1, 3, 1)
      qvTilde(:, :, - 1, 3, 1) = qvTilde(:, :, nz - 1, 3, 1)

      if(verbose .and. master) print *, "setBoundary_z_periodic:   z-periodic &
          &BC for iceTilde variables set."

    case("iceFlux")
      ! nothing

    case("varTilde")
      !-----------------------------------
      !              varTilde
      !-----------------------------------

      ! the following three boundary-condition calls can probably be
      ! removed
      ! probably only necessary for ALDM (that is not used any more
      ! anyway)

      if(updateMass) then
        ! reconstructed density needed in ghost cell k = nz+2
        rhoTilde(:, :, nz + 2, 3, 0) = rhoTilde(:, :, 2, 3, 0)

        ! ...in ghost cell j = -1
        rhoTilde(:, :, - 1, 3, 1) = rhoTilde(:, :, nz - 1, 3, 1)

        if(verbose .and. master) then
          print *, "setBoundary_z_periodic:  z-periodic BC for rhoTilde set."
        end if
      end if

      if(timeScheme == "semiimplicit" .or. auxil_equ) then
        ! reconstructed density needed in ghost cell k = nz+2
        rhopTilde(:, :, nz + 2, 3, 0) = rhopTilde(:, :, 2, 3, 0)

        ! ...in ghost cell j = -1
        rhopTilde(:, :, - 1, 3, 1) = rhopTilde(:, :, nz - 1, 3, 1)

        if(verbose .and. master) then
          print *, "setBoundary_z_periodic:  z-periodic BC for rhoTilde set."
        end if
      end if

      if(updateIce) then
        ! see above, similar to rho
        nAerTilde(:, :, nz + 2, 3, 0) = nAerTilde(:, :, 2, 3, 0)
        nIceTilde(:, :, nz + 2, 3, 0) = nIceTilde(:, :, 2, 3, 0)
        qIceTilde(:, :, nz + 2, 3, 0) = qIceTilde(:, :, 2, 3, 0)
        qvTilde(:, :, nz + 2, 3, 0) = qvTilde(:, :, 2, 3, 0)

        ! see above, similar to rho
        nAerTilde(:, :, - 1, 3, 1) = nAerTilde(:, :, nz - 1, 3, 1)
        nIceTilde(:, :, - 1, 3, 1) = nIceTilde(:, :, nz - 1, 3, 1)
        qIceTilde(:, :, - 1, 3, 1) = qIceTilde(:, :, nz - 1, 3, 1)
        qvTilde(:, :, - 1, 3, 1) = qvTilde(:, :, nz - 1, 3, 1)

        if(verbose .and. master) print *, "setBoundary_z_periodic:   &
            &z-periodic BC for iceTilde variables set."

      end if

    case("flux")
      !-----------------------------------
      !              flux
      !-----------------------------------
      return

    case default
      stop "setBoundary_y: unknown option."
    end select

  end subroutine setBoundary_z_periodic

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

        if(testCase == "baroclinic_LC") then
          if(topography) then
            ! TFC FJ
            var%rho(1:nx, 1:ny, 0) = - var%rho(1:nx, 1:ny, 1) + dens_env_pp(:, &
                &:, 1) - rhoStratTFC(1:nx, 1:ny, 1) + dens_env_pp(:, :, 0) &
                &- rhoStratTFC(1:nx, 1:ny, 0)

            var%rho(1:nx, 1:ny, nz + 1) = - var%rho(1:nx, 1:ny, nz) &
                &+ dens_env_pp(:, :, nz) - rhoStratTFC(1:nx, 1:ny, nz) &
                &+ dens_env_pp(:, :, nz + 1) - rhoStratTFC(1:nx, 1:ny, nz + 1)
          else
            var%rho(1:nx, 1:ny, 0) = - var%rho(1:nx, 1:ny, 1) + dens_env_pp(:, &
                &:, 1) - rhoStrat(1) + dens_env_pp(:, :, 0) - rhoStrat(0)

            var%rho(1:nx, 1:ny, nz + 1) = - var%rho(1:nx, 1:ny, nz) &
                &+ dens_env_pp(:, :, nz) - rhoStrat(nz) + dens_env_pp(:, :, nz &
                &+ 1) - rhoStrat(nz + 1)
          end if
        else if(testCase == "smoothVortex") then
          var%rho(1:nx, 1:ny, 0) = var%rho(1:nx, 1:ny, 1)
          var%rho(1:nx, 1:ny, nz + 1) = var%rho(1:nx, 1:ny, nz)
        else
          do k = 1, nbz

            var%rho(:, :, - k + 1) = - var%rho(:, :, k)
            var%rho(:, :, nz + k) = - var%rho(:, :, nz - k + 1)

          end do
        end if

        if(timeScheme == "semiimplicit" .or. auxil_equ) then
          ! vertical boundary condition for density fluctuations

          if(testCase == "baroclinic_LC") then
            if(topography) then
              ! TFC FJ
              var%rhop(1:nx, 1:ny, 0) = - var%rhop(1:nx, 1:ny, 1) &
                  &+ dens_env_pp(:, :, 1) - rhoStratTFC(1:nx, 1:ny, 1) &
                  &+ dens_env_pp(:, :, 0) - rhoStratTFC(1:nx, 1:ny, 0)

              var%rhop(1:nx, 1:ny, nz + 1) = - var%rhop(1:nx, 1:ny, nz) &
                  &+ dens_env_pp(:, :, nz) - rhoStratTFC(1:nx, 1:ny, nz) &
                  &+ dens_env_pp(:, :, nz + 1) - rhoStratTFC(1:nx, 1:ny, nz + 1)
            else
              var%rhop(1:nx, 1:ny, 0) = - var%rhop(1:nx, 1:ny, 1) &
                  &+ dens_env_pp(:, :, 1) - rhoStrat(1) + dens_env_pp(:, :, 0) &
                  &- rhoStrat(0)

              var%rhop(1:nx, 1:ny, nz + 1) = - var%rhop(1:nx, 1:ny, nz) &
                  &+ dens_env_pp(:, :, nz) - rhoStrat(nz) + dens_env_pp(:, :, &
                  &nz + 1) - rhoStrat(nz + 1)
            end if
          else if(testCase == "smoothVortex") then
            var%rhop(1:nx, 1:ny, 0) = var%rhop(1:nx, 1:ny, 1)
            var%rhop(1:nx, 1:ny, nz + 1) = var%rhop(1:nx, 1:ny, nz)
          else
            do k = 1, nbz
              var%rhop(:, :, - k + 1) = - var%rhop(:, :, k)
              var%rhop(:, :, nz + k) = - var%rhop(:, :, nz - k + 1)
            end do
          end if
        end if
      end if

      if(model == "compressible") then
        ! set boundary for P mass weighted pot. temp.
        do k = 1, nbz
          var%P(:, :, - k + 1) = var%P(:, :, k)
          var%P(:, :, nz + k) = var%P(:, :, nz - k + 1)
        end do
      end if

      if(updateTracer) then
        do k = 1, nbz
          var%chi(:, :, - k + 1) = - var%chi(:, :, k)
          var%chi(:, :, nz + k) = - var%chi(:, :, nz - k + 1)
        end do
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

        ! if (testCase == "baroclinic_LC") then
        !    ! u
        !    var(1:nx,1:ny,0,2) &
        !    = var(1:nx,1:ny,1,2) - u_env_pp (:,:,1) + u_env_pp (:,:,0)
        !    var(1:nx,1:ny,nz+1,2) &
        !    = var(1:nx,1:ny,nz,2) - u_env_pp (:,:,nz) + u_env_pp (:,:,nz+1)

        !    ! v
        !    var(1:nx,1:ny,0,3) &
        !    = var(1:nx,1:ny,1,3) - v_env_pp (:,:,1) + v_env_pp (:,:,0)
        !    var(1:nx,1:ny,nz+1,3) &
        !    = var(1:nx,1:ny,nz,3) - v_env_pp (:,:,nz) + v_env_pp (:,:,nz+1)
        !   else
        if(testCase == "baroclinic_LC") then
          do k = 1, nbz
            ! u
            var%u(:, :, - k + 1) = - var%u(:, :, k) + var_env%u(:, :, k) &
                &+ var_env%u(:, :, - k + 1)
            var%u(:, :, nz + k) = - var%u(:, :, nz - k + 1) + var_env%u(:, :, &
                &nz - k + 1) + var_env%u(:, :, nz + k)

            ! v
            var%v(:, :, - k + 1) = - var%v(:, :, k)
            var%v(:, :, nz + k) = - var%v(:, :, nz - k + 1)
          end do
        else if(testCase == "smoothVortex") then
          var%u(1:nx, 1:ny, 0) = var%u(1:nx, 1:ny, 1)
          var%u(1:nx, 1:ny, nz + 1) = var%u(1:nx, 1:ny, nz)
          var%v(1:nx, 1:ny, 0) = var%v(1:nx, 1:ny, 1)
          var%v(1:nx, 1:ny, nz + 1) = var%v(1:nx, 1:ny, nz)
        else
          ! TFC FJ
          ! No-slip has been changed to free-slip!
          ! Free-slip condition.
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
      end if

      if(correctMomentum) then
        ! set gradient at vertical boundary to 0
        ! life-cycle simulation: only gradient of deviations from
        !                        environmental pressure are set to 0

        if(testCase == "baroclinic_LC") then
          ! z = 0
          var%pi(1:nx, 1:ny, 0) = var%pi(1:nx, 1:ny, 1) - p_env_pp(:, :, 1) &
              &+ p_env_pp(:, :, 0)

          ! z = zMax
          var%pi(1:nx, 1:ny, nz + 1) = var%pi(1:nx, 1:ny, nz) - p_env_pp(:, :, &
              &nz) + p_env_pp(:, :, nz + 1)
        else if(testCase == "smoothVortex") then
          var%pi(1:nx, 1:ny, 0) = var%pi(1:nx, 1:ny, 1)
          var%pi(1:nx, 1:ny, nz + 1) = var%pi(1:nx, 1:ny, nz)
        else
          ! z = 0
          var%pi(:, :, 0) = var%pi(:, :, 1)

          ! z = zMax
          var%pi(:, :, nz + 1) = var%pi(:, :, nz)
        end if
      end if

      ! additional call to setHalos in case baroclinic_LC in order to set
      ! the overlap betwen vertical and horizontal halos right

      if(testCase == "baroclinic_LC") call setHalos(var, "var")

    case("ice")
      ! reflect at boundary with change of sign
      ! at boundary var = 0
      do iVar = 1, nVarIce
        do k = 1, nbz
          var%ICE(:, :, - k + 1, iVar) = - var%ICE(:, :, k, iVar)
          var%ICE(:, :, nz + k, iVar) = - var%ICE(:, :, nz - k + 1, iVar)
        end do
      end do

    case("iceTilde")
      !nothing

    case("iceFlux")

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

        if(verbose .and. master) then
          print *, "boundary.f90/verticalBoundary: vertical BC for rho set"
        end if

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

      ! set vertical tracer fluxes at wall to 0
      ! if including tracer
      if(updateTracer) then
        flux%chi(:, :, 0, 3) = 0.0
        flux%chi(:, :, nz, 3) = 0.0
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

        if(verbose .and. master) then
          print *, "boundary.f90/verticalBoundary: vertical flux-BC for u,v,w &
              &set"
        end if

        ! replace flux by CDS fluxes at upper / lower region

        if(uFluxCorr) stop 'ERROR: uFluxCorr = .false. expected'

        if(vFluxCorr) stop 'ERROR: vFluxCorr = .false. expected'

        if(wFluxCorr) stop 'ERROR: wFluxCorr = .false. expected'
      end if ! predictMomentum

    case default
      stop "setBoundary_z: unknown option."
    end select

  end subroutine setBoundary_z_solidWall

end module boundary_module
