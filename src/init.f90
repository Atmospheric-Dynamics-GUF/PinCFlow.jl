module init_module

  use type_module
  use atmosphere_module
  use mpi_module
  use boundary_module
  use mpi
  use output_netCDF_module

  implicit none

  private

  public :: initialise
  public :: setup

  contains

  subroutine setup(var, var0, var1, flux, flux0, dRho, dRhop, dMom)

    !--------------------------------------------
    ! Allocate var and flux / read the namelists
    !--------------------------------------------

    ! in/out variables
    type(var_type), intent(out) :: var, var0, var1
    type(flux_type), intent(out) :: flux, flux0
    real, dimension(:, :, :), allocatable :: dRho, dRhop ! RK-Update for rho
    real, dimension(:, :, :, :), allocatable :: dMom ! ...rhoU,rhoV,rho

    integer :: allocstat
    integer :: i, j, k

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

    ! Allocate resolved topography.
    allocate(topography_surface(- nbx:nx + nbx, - nby:ny + nby), stat &
        &= allocstat)
    if(allocstat /= 0) stop "setup: could not allocate topography_surface"
    allocate(zTildeTFC(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "setup: could not allocate zTildeTFC"
    allocate(zTFC(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "setup: could not allocate zTFC"
    allocate(zTildeS(- nbz:nz + nbz), stat = allocstat)
    if(allocstat /= 0) stop "setup: could not allocate zTildeS"
    allocate(zS(- nbz:nz + nbz), stat = allocstat)
    if(allocstat /= 0) stop "setup: could not allocate zS"

    !-------------------------------------
    !      allocate variable fields
    !-------------------------------------

    call allocate_var_type(var)
    call reset_var_type(var)

    call allocate_var_type(var0)
    call reset_var_type(var0)

    call allocate_var_type(var1)
    call reset_var_type(var1)

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

    call allocate_flux_type(flux)
    call reset_flux_type(flux)

    call allocate_flux_type(flux0)
    call reset_flux_type(flux0)

    allocate(kr_sp_tfc(0:(nx + 1), 0:(ny + 1), 0:(nz + 1)), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate kr_sp_tfc"

    allocate(kr_sp_w_tfc(0:(nx + 1), 0:(ny + 1), 0:(nz + 1)), stat = allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate kr_sp_w_tfc"

    ! Allocate damping coefficient for unified sponge.
    if(spongeLayer .and. unifiedSponge) then
      allocate(alphaUnifiedSponge(0:(nx + 1), 0:(ny + 1), 0:(nz + 1)), stat &
          &= allocstat)
      if(allocstat /= 0) stop "init.f90: could not allocate alphaUnifiedSponge"
    end if

    ! Safety switch for halos in TFC
    if(nbx < 2 .or. nby < 2 .or. nbz < 2) then
      stop "Two halos / ghost cells are needed in each direction!"
    end if

    !---------------------------------------
    !        Model equation settings
    !---------------------------------------

    select case(model)

    case("pseudo_incompressible")

      updateMass = .true.
      predictMomentum = .true.
      correctMomentum = .true.

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
    type(var_type), intent(inout) :: var
    type(flux_type), intent(inout) :: flux

    ! local variables
    integer :: i, j, k
    integer :: allocstat

    select case(testCase)

    case('mountainwave')
      ! for wave resolving simulation of mountain waves:
      ! read parameters for temporary wind relaxation
      ! zero-wind initial state for montain-wave simulations

      ! nondimensionalization
      var%u(:, :, :) = backgroundFlow_dim(1) / uRef
      var%v(:, :, :) = backgroundFlow_dim(2) / uRef
      var%w(:, :, :) = backgroundFlow_dim(3) / uRef

      ! density and pressure
      do k = 0, (nz + 1)
        do j = 0, (ny + 1)
          do i = 0, (nx + 1)
            select case(model)
            case("pseudo_incompressible")

              ! initialization zero density fluctuations
              var%rho(i, j, k) = 0.0

            case default
              stop "initialize: unknown case model"
            end select

            ! initialization zero pressure fluctuations
            var%pi(i, j, k) = 0.0

          end do
        end do
      end do

    case default

      print *, "init.f90/initialise: testCase = ", testCase
      stop "init.f90/initialise: This testCase is not valid. Stop."

    end select

  end subroutine initialise

end module init_module
