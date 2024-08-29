module output_module

  use type_module
  use atmosphere_module
  use sizeof_module
  use ice2_sub_module, ONLY:output_ice
  implicit none

  private

  ! public subroutines
  public :: output_data
  public :: read_data
  public :: read_profile
  public :: output_wkb
  public :: output_field
  public :: output_profile
  public :: output_fluxes
  public :: output_background

  ! internal subroutines (listed for completeness)
  public :: init_output
  public :: terminate_output

  ! internal module variables (output variables)
  real, dimension(:, :, :), allocatable :: optVar

  contains

  subroutine output_data(iOut, var, iTime, time, cpuTime)

    !-------------------------------
    !  writes data to file pf_all.dat
    !-------------------------------

    ! output counter
    integer, intent(inout) :: iOut

    ! argument fields
    type(var_type), intent(in) :: var

    ! local and global output field and record nbs.
    real * 4, dimension(nx, ny) :: field_prc
    integer irc_prc, irc_out

    ! argument parameters
    integer, intent(in) :: iTime
    real, intent(in) :: time, cpuTime

    ! local variables
    integer :: i, j, k, iVar
    real :: time_dim

    real :: rhotracer
    real, dimension(1:nx, 1:ny, 1:nz) :: tracerdrho
    real, dimension(1:nx, 1:nz) :: meantracer

    ! needed for output to screen
    character(len = 20) :: fmt, form
    character(len = 40) :: cpuTimeChar

    ! CPU Time
    integer :: days, hours, mins, secs
    real :: timeVar

    ! hotBubble output
    real :: rho, theta_dim

    ! buoyancy
    real :: b, b_dim, theta

    integer :: i_prc, i_mst, i_out, j_prc, j_mst, j_out

    time_dim = time * tRef

    if(master) then ! modified by Junhong Wei (20161110)
      print *, ""
      print *, " Output into File"
      print *, ""
      write(*, fmt = "(a25,i15)") " at time step = ", iTime
      write(*, fmt = "(a25,f15.1,a8)") " at physical time = ", time_dim, " &
          &seconds"
    end if ! modified by Junhong Wei (20161110)

    !------------------------------
    !   prepare output file
    !------------------------------

    ! open output file

    if(master) then
      open(41, file = 'pf_all.dat', form = "unformatted", access = 'direct', &
          &recl = SizeX * SizeY * sizeofreal4)
    end if

    ! calc cpu-time in days/hours/minutes/seconds
    timeVar = cpuTime
    days = floor(timeVar / 86400.0)
    timeVar = timeVar - 86400.0 * days
    hours = floor(timeVar / 3600.0)
    timeVar = timeVar - 3600.0 * hours
    mins = floor(timeVar / 60.0)
    timeVar = timeVar - 60.0 * mins
    secs = int(timeVar)

    write(unit = cpuTimeChar, fmt = "(i2,a,2(i2.2,a),i2.2)") days, " ", hours, &
        &":", mins, ":", secs

    if(master) write(*, fmt = "(a25,a25)") "CPU time = ", cpuTimeChar

    !---------------------------------------
    !       dimensionalising and layerwise output
    !---------------------------------------

    irc_prc = 0
    do iVar = 1, nVar
      if(varOut(iVar) /= "") irc_prc = irc_prc + 1
    end do
    irc_prc = irc_prc * iOut * nz

    ! just for safety.
    if(include_tracer .and. include_ice) then
      stop "output.f90: tracer and ice not an option yet."
    end if

    do iVar = 1, nVar
      if(varOut(iVar) /= "") then
        do k = 1, nz
          ! dimensionalization

          do j = 1, ny
            do i = 1, nx
              select case(varOut(iVar))

              case("rho") ! density
                if(fluctuationMode) then
                  if(rhoOffset) then
                    field_prc(i, j) = var%rho(i, j, k) * rhoRef
                  else
                    if(topography) then
                      ! TFC FJ
                      ! Adjustment for 3D background field in
                      ! TFC.
                      field_prc(i, j) = (var%rho(i, j, k) + rhoStratTFC(i, j, &
                          &k)) * rhoRef
                    else
                      field_prc(i, j) = (var%rho(i, j, k) + rhoStrat(k)) &
                          &* rhoRef
                    end if
                  end if
                else
                  if(rhoOffset) then
                    field_prc(i, j) = (var%rho(i, j, k) - rhoStrat(k)) * rhoRef
                  else
                    field_prc(i, j) = var%rho(i, j, k) * rhoRef
                  end if
                end if

                ! average velocities to cell center

              case("u") ! u velocity
                field_prc(i, j) = (0.5 * (var%u(i, j, k) + var%u(i - 1, j, &
                    &k))) * uRef

              case("v") ! v velocity
                field_prc(i, j) = (0.5 * (var%v(i, j, k) + var%v(i, j - 1, &
                    &k))) * uRef

              case("w") ! w velocity
                field_prc(i, j) = (0.5 * (var%w(i, j, k) + var%w(i, j, k &
                    &- 1))) * uRef

              case("pi") ! Exner function pi'
                !(deviation from background)

                field_prc(i, j) = var%pi(i, j, k)

              case("theta") ! potential temperature theta'
                ! (deviation from background, Boussinesq)

                select case(model)

                case("pseudo_incompressible", "compressible")
                  if(fluctuationMode) then
                    if(topography) then
                      ! TFC FJ
                      rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)
                    else
                      rho = var%rho(i, j, k) + rhoStrat(k)
                    end if
                  else
                    rho = var%rho(i, j, k)
                  end if

                  if(topography) then
                    ! TFC FJ
                    theta_dim = pStratTFC(i, j, k) / rho * thetaRef
                  else
                    theta_dim = Pstrat(k) / rho * thetaRef
                  end if

                  if(thetaOffset) then
                    if(topography) then
                      theta_dim = theta_dim - thetaStratTFC(i, j, k) * thetaRef
                    else
                      theta_dim = theta_dim - thetaStrat(k) * thetaRef
                    end if
                  end if

                  field_prc(i, j) = theta_dim !thetaRef*(PStrat(k+1)/(var(i,j,k+1,1)+rhoStrat(k+1))-PStrat(k-1)/(var(i,j,k-1,1)+rhoStrat(k-1)))/(2.*dz*lRef)!theta_dim

                case("Boussinesq")
                  ! TFC FJ
                  ! Boussinesq: density fluctuations are stored in
                  ! var(:, :, :, 6)!
                  field_prc(i, j) = - var%rhop(i, j, k) * theta00 / rho00 &
                      &* thetaRef

                  ! field_prc(i,j) = var(i,j,k,iVar)*thetaRef

                case("WKB")

                case default
                  stop "Unknown model in output_data!"
                end select ! model

              case("DSC") ! dynamic Smagorinsky coefficient
                field_prc(i, j) = var%DSC(i, j, k) * uRef * lRef

              case("GWH") ! gravity-wave heating
                field_prc(i, j) = var%GWH(i, j, k) * rhoRef * thetaRef / tRef

              case("P") ! mass-weighted potential temperature
                field_prc(i, j) = pStratTFC(i, j, k) * rhoRef * thetaRef

              case("ICE-1") ! first ice variable
                field_prc(i, j) = real(var%ICE(i, j, k, 1) / (rhoRef * lRef &
                    &** 3), kind = 4)

              case("ICE-2") ! second ice variable
                field_prc(i, j) = real(var%ICE(i, j, k, 2) / (rhoRef * lRef &
                    &** 3), kind = 4)

              case("ICE-3") ! third ice variable
                field_prc(i, j) = real(var%ICE(i, j, k, 3), kind = 4)

              case("ICE-4") ! fourth ice variable
                field_prc(i, j) = real(var%ICE(i, j, k, 4), kind = 4)

              case("ICE2-N") ! first alternative ice variable
                call output_ice(i, j, k, inN, var, field_prc)

              case("ICE2-Q") ! second alternative ice variable
                call output_ice(i, j, k, inQ, var, field_prc)

              case("ICE2-QV") ! third alternative ice variable
                call output_ice(i, j, k, inQv, var, field_prc)

              case("OPT-1") ! first optional variable
                field_prc(i, j) = var%OPT(i, j, k, 1)
              case("OPT-2") ! second optional variable
                field_prc(i, j) = var%OPT(i, j, k, 2)
              case("OPT-3") ! third optional variable
                field_prc(i, j) = var%OPT(i, j, k, 3)

              case("chi") ! tracer

                if(fluctuationMode) then
                  if(topography) then
                    ! TFC FJ
                    ! Adjustment for 3D background field in
                    ! TFC.
                    rhotracer = (var%rho(i, j, k) + rhoStratTFC(i, j, k))
                  else
                    rhotracer = (var%rho(i, j, k) + rhoStrat(k))
                  end if
                  !rhotracer = var(i,j,k,iVar)
                else
                  rhotracer = var%rho(i, j, k)
                end if

                if(tracerdifference) then
                  field_prc(i, j) = var%chi(i, j, k) / rhotracer &
                      &- initialtracer(i, j, k)
                else
                  field_prc(i, j) = var%chi(i, j, k) / rhotracer
                end if
              case default
                stop "Error in output_data: unknown output variable!"
              end select ! iVar
              ! Subtract offset.
              field_prc(i, j) = field_prc(i, j) - offset(iVar)
            end do ! i
            call mpi_gather(field_prc(1, j), nx, mpi_real, field_mst(1, j), &
                &nx, mpi_real, 0, comm, ierror)
          end do ! j

          ! layerwise output
          irc_prc = irc_prc + 1
          !            write(40,rec=irc_prc) field_prc
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

            write(41, rec = irc_prc) field_out
          end if
        end do ! k
      end if
    end do ! iVar

    if(master) close(unit = 41)

    ! set counter
    iOut = iOut + 1

  end subroutine output_data

  !-------------------------------------------------------------------------

  subroutine read_data(iIn, var, time)

    !-------------------------------
    !  reads data from file pf_all_in.dat
    !-------------------------------

    ! input counter
    integer, intent(in) :: iIn

    ! argument fields
    type(var_type), intent(inout) :: var

    real, intent(out) :: time

    real :: rhotracer

    ! local and global output field and record nbs.
    real * 4, dimension(nx, ny) :: field_prc
    integer irc_prc, irc_out

    ! local variables
    integer :: i, j, k, iVar

    ! needed for output to screen
    character(len = 20) :: fmt, form

    ! hotBubble input
    real :: rho, theta_dim

    ! buoyancy
    real :: b, b_dim, theta

    integer :: i_prc, i_mst, i_out, j_prc, j_mst, j_out

    if(master) then
      print *, ""
      print *, " Input from File "
      print *, ""
      write(*, fmt = "(a25,i15)") " reading record no. ", iIn
    end if

    ! open input file

    if(master) then
      open(40, file = 'pf_all_in.dat', form = "unformatted", access &
          &= 'direct', recl = SizeX * SizeY * sizeofreal4)
      print *, "pf_all_in.dat opened"
    end if

    call reset_var_type(var)

    time = (iIn + 10) * outputTimeDiff / tRef

    !---------------------------------------
    !       layerwise input and non-dimensionalizing
    !---------------------------------------

    irc_prc = 0
    do iVar = 1, nVar
      if(varIn(iVar) /= "") irc_prc = irc_prc + 1
    end do
    irc_prc = irc_prc * iIn * nz

    do iVar = 1, nVar
      if(varIn(iVar) /= "") then
        do k = 1, nz
          ! read data layerwise

          irc_prc = irc_prc + 1

          if(master) then
            read(40, rec = irc_prc) field_out

            do j = 1, ny
              j_mst = j

              do j_prc = 1, nprocy
                j_out = ny * (j_prc - 1) + j

                do i_prc = 1, nprocx
                  do i = 1, nx
                    i_out = nx * (i_prc - 1) + i

                    i_mst = nprocy * nx * (i_prc - 1) + (j_prc - 1) * nx + i

                    field_mst(i_mst, j_mst) = field_out(i_out, j_out)
                  end do
                end do
              end do
            end do
          end if

          call mpi_barrier(comm, ierror)

          do j = 1, ny
            ! data distributed over all processors
            call mpi_scatter(field_mst(1, j), nx, mpi_real, field_prc(1, j), &
                &nx, mpi_real, 0, comm, ierror)
            do i = 1, nx
              ! Add offset.
              field_prc(i, j) = field_prc(i, j) + offset(iVar)

              select case(varIn(iVar))

              case("rho") ! density
                if(fluctuationMode) then
                  if(rhoOffset) then
                    var%rho(i, j, k) = field_prc(i, j) / rhoRef
                  else
                    if(topography) then
                      ! TFC FJ
                      ! Adjustment for 3D background field in TFC.
                      var%rho(i, j, k) = field_prc(i, j) / rhoRef &
                          &- rhoStratTFC(i, j, k)
                    else
                      var%rho(i, j, k) = field_prc(i, j) / rhoRef - rhoStrat(k)
                    end if
                  end if
                else
                  if(rhoOffset) then
                    var%rho(i, j, k) = field_prc(i, j) / rhoRef + rhoStrat(k)
                  else
                    var%rho(i, j, k) = field_prc(i, j) / rhoRef
                  end if
                end if

                ! interpolate velocities to cell faces

              case("u") ! u velocity
                var%u(i, j, k) = field_prc(i, j) / uRef

              case("v") ! v velocity
                var%v(i, j, k) = field_prc(i, j) / uRef

              case("w") ! w velocity
                var%w(i, j, k) = field_prc(i, j) / uRef

              case("pi") ! Exner function pi'
                !(deviation from background)
                var%pi(i, j, k) = field_prc(i, j)

              case("theta")
                select case(model)
                case("pseudo_incompressible", "compressible")
                  continue
                case("Boussinesq")
                  var%rhop(i, j, k) = - field_prc(i, j) / theta00 * rho00 &
                      &/ thetaRef
                case default
                  stop "Unknown model in read_data!"
                end select

              case("DSC") ! dynamic Smagorinsky coefficient
                var%DSC(i, j, k) = field_prc(i, j) / (uRef * lRef)

              case("GWH") ! gravity-wave heating
                var%GWH(i, j, k) = field_prc(i, j) / rhoRef / thetaRef * tRef

              case("P") ! mass-weighted potential temperature
                var%P(i, j, k) = field_prc(i, j) / rhoRef / thetaRef

              case("chi") ! tracer
                var%chi(i, j, k) = field_prc(i, j) + initialtracer(i, j, k)

              case default

                if(include_ice .or. include_ice2) then

                  stop "read_data: include_ice(2) not possible for restart"

                end if

                ! !--------------------------------------
                ! ! NEW: ice cases !    ! might be better to include reference units here
                ! if(iVar == nVar - 3) then ! aerosol particle number concentration nAer
                !   var(i, j, k, iVar) = field_prc(i, j) * rhoRef * lRef ** 3
                ! else if(iVar == nVar - 2) then ! ice particle number concentration nIce
                !   var(i, j, k, iVar) = field_prc(i, j) * rhoRef * lRef ** 3
                ! else if(iVar == nVar - 1) then ! ice particle mass concentration qIce
                !   var(i, j, k, iVar) = field_prc(i, j)
                ! else if(iVar == nVar) then ! water vapor mass concentration qv
                !   var(i, j, k, iVar) = field_prc(i, j)
                !   !---------------------------------
                ! else
                !   stop "tec360: unknown iVar"
                ! end if
              end select ! iVar
            end do ! i
          end do ! j
        end do ! k
      end if
    end do ! iVar

    ! Interpolate velocities to cell interfaces (halos are set later).
    var%u(1:nx, :, :) = 0.5 * (var%u(1:nx, :, :) + var%u(2:nx + 1, :, :))
    var%v(:, 1:ny, :) = 0.5 * (var%v(:, 1:ny, :) + var%v(:, 2:ny + 1, :))
    var%w(:, :, 1:nz) = 0.5 * (var%w(:, :, 1:nz) + var%w(:, :, 2:nz + 1))

    if(include_tracer) then
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            if(fluctuationMode) then
              if(topography) then
                rhotracer = (var%rho(i, j, k) + rhoStratTFC(i, j, k))
              else
                rhotracer = (var%rho(i, j, k) + rhoStrat(k))
              end if
            else
              rhotracer = var%rho(i, j, k)
            end if
            var%chi(i, j, k) = var%chi(i, j, k) * rhotracer
          end do
        end do
      end do
    end if

    !------------------------------------
    !              close file
    !------------------------------------

    if(master) close(unit = 40)

  end subroutine read_data

  !-------------------------------------------------------------------------

  subroutine output_wkb(iOut, ray, ray_var3D)

    !SD include opt_ray
    !*use type_module, ONLY : opt_ray
    !--------------------------------------------------
    !  WKB output to pf_wkb_mean.dat and pf_wkb_ray.dat
    !--------------------------------------------------

    ! output counter
    integer, intent(inout) :: iOut

    ! wkb arguments
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), &
        &intent(in) :: ray

    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, 1:13), intent(in) :: ray_var3D

    ! local variables
    integer :: i, j, k, iVar

    ! local and global output field and record nbs.
    real * 4, dimension(nx, ny) :: field_prc
    integer irc_prc, irc_out

    integer :: i_prc, i_mst, i_out, j_prc, j_mst, j_out

    !SD
    integer :: NoR, iRay, TNoR, cRay, NoR_prc
    real, dimension(NFR) :: vct

    real, dimension(nray_wrk) :: dens
    integer :: dim, jRay, iMax

    ! FJFeb2023
    ! set counter
    ! iOut = iOut - 1

    !------------------------------
    !   prepare output file
    !------------------------------

    ! open output file

    if(master) then
      open(40, file = 'pf_wkb_mean.dat', form = "unformatted", access &
          &= 'direct', recl = SizeX * SizeY * sizeofreal4)
    end if

    !---------------------------------------
    !       dimensionalising and layerwise output
    !---------------------------------------

    ! FJFeb2023
    irc_prc = 13 * (iOut - 1) * nz

    do iVar = 1, 13
      do k = 1, nz
        ! dimensionalization

        do j = 1, ny
          do i = 1, nx
            select case(iVar)

            case(1) ! (du/dt)_GW
              field_prc(i, j) = ray_var3D(i, j, k, 1) * uRef / tRef

            case(2) ! (dv/dt)_GW
              field_prc(i, j) = ray_var3D(i, j, k, 2) * uRef / tRef

            case(3) ! (dtheta/dt)_GW
              field_prc(i, j) = ray_var3D(i, j, k, 3) * thetaRef / tRef

            case(4) ! elastic term in x direction
              field_prc(i, j) = ray_var3D(i, j, k, 4) * uRef / tRef

              !case(5) ! elastic term in y direction
              !  field_prc(i,j) = ray_var3D(i,j,k,5) * uRef/tRef
            case(5) ! u'w' momentum flux , output chage by FDK
              field_prc(i, j) = ray_var3D(i, j, k, 5) * rhoRef * uRef ** 2

            case(6) ! GW energy
              field_prc(i, j) = ray_var3D(i, j, k, 6) * uRef ** 2 ! /tRef deleted by FDK; rhoRef *

            case(7) ! u'chi' (zonal gw tracer flux [m/s])
              if(include_tracer) then
                field_prc(i, j) = ray_var3D(i, j, k, 7) !* uRef
              else
                field_prc(i, j) = 0.0
              end if

            case(8) ! v'chi' (meridional gw tracer flux [m/s])
              if(include_tracer) then
                field_prc(i, j) = ray_var3D(i, j, k, 8) !* uRef
              else
                field_prc(i, j) = 0.0
              end if

            case(9) ! w'chi' (vertical gw tracer flux [m/s])
              if(include_tracer) then
                field_prc(i, j) = ray_var3D(i, j, k, 9) * uRef
              else
                field_prc(i, j) = 0.0
              end if

            case(10) ! w'chi' (next-order vertical gw tracer flux [m/s])
              ! currently leading-order buoyancy wave amplitude |bhat2^(2)|
              if(include_tracer) then
                field_prc(i, j) = ray_var3D(i, j, k, 10) !* uRef
              else
                field_prc(i, j) = 0.0
              end if

            case(11)
              ! tracer forcing (leading order gw tracer flux convergence [m^2/s])
              if(include_tracer) then
                field_prc(i, j) = ray_var3D(i, j, k, 11) / tRef
              else
                field_prc(i, j) = 0.0
              end if

            case(12)
              ! tracer forcing (next-order gw tracer flux convergence [m^2/s])
              if(include_tracer) then
                field_prc(i, j) = ray_var3D(i, j, k, 12) / tRef
              else
                field_prc(i, j) = 0.0
              end if

            case(13) ! tracer diffusion (diffusive mixing of tracer [1/s])
              if(include_tracer) then
                field_prc(i, j) = ray_var3D(i, j, k, 13) / tRef
              else
                field_prc(i, j) = 0.0
              end if

            case default
              stop "output_wkb: unkown iVar"
            end select ! iVar
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

          write(40, rec = irc_prc) field_out
        end if
      end do ! k
    end do ! iVar

    !------------------------------------
    !              close file
    !------------------------------------

    if(master) close(unit = 40)

    ! !------------------------------------
    ! !        Number of ray volumes
    ! !------------------------------------

    ! if(master) then
    !   open(40, file = 'pf_wkb_nray.dat', form = "unformatted", access &
    !       = 'direct', recl = SizeX * SizeY * sizeofreal4)
    ! end if

    ! irc_prc = (iOut - 1) * nz

    ! do k = 1, nz
    !   do j = 1, ny
    !     do i = 1, nx
    !       field_prc(i, j) = real(nRay(i, j, k), kind = 4)
    !     end do
    !     call mpi_gather(field_prc(1, j), nx, mpi_real, field_mst(1, j), nx, &
    !         mpi_real, 0, comm, ierror)
    !   end do

    !   irc_prc = irc_prc + 1
    !   call mpi_barrier(comm, ierror)
    !   if(master) then
    !     do j = 1, ny
    !       j_mst = j
    !       do j_prc = 1, nprocy
    !         j_out = ny * (j_prc - 1) + j
    !         do i_prc = 1, nprocx
    !           do i = 1, nx
    !             i_out = nx * (i_prc - 1) + i
    !             i_mst = nprocy * nx * (i_prc - 1) + (j_prc - 1) * nx + i
    !             field_out(i_out, j_out) = field_mst(i_mst, j_mst)
    !           end do
    !         end do
    !       end do
    !     end do

    !     write(40, rec = irc_prc) field_out
    !   end if
    ! end do

    !------------------------------------
    !        Ray volume properties
    !------------------------------------

    dim = 1
    if(fac_dk_init /= 0.0) dim = dim + 1
    if(fac_dl_init /= 0.0) dim = dim + 1

    if(master) then
      open(40, file = 'pf_wkb_rays.dat', form = "unformatted", access &
          &= 'direct', recl = SizeX * SizeY * sizeofreal4)
    end if

    irc_prc = (iOut - 1) * 14 * nRayOutput * nz

    do iVar = 1, 14
      do iRay = 1, nRayOutput
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              dens(:) = abs(ray(:, i, j, k)%dens)
              do jRay = 1, iRay
                iMax = maxloc(dens, dim = 1)
                dens(iMax) = - 1.0
              end do
              select case(iVar)
              case(1)
                field_prc(i, j) = real(ray(iMax, i, j, k)%x * lRef, kind = 4)
              case(2)
                field_prc(i, j) = real(ray(iMax, i, j, k)%y * lRef, kind = 4)
              case(3)
                field_prc(i, j) = real(ray(iMax, i, j, k)%z * lRef, kind = 4)
              case(4)
                field_prc(i, j) = real(ray(iMax, i, j, k)%k / lRef, kind = 4)
              case(5)
                field_prc(i, j) = real(ray(iMax, i, j, k)%l / lRef, kind = 4)
              case(6)
                field_prc(i, j) = real(ray(iMax, i, j, k)%m / lRef, kind = 4)
              case(7)
                field_prc(i, j) = real(ray(iMax, i, j, k)%dxray * lRef, kind &
                    &= 4)
              case(8)
                field_prc(i, j) = real(ray(iMax, i, j, k)%dyray * lRef, kind &
                    &= 4)
              case(9)
                field_prc(i, j) = real(ray(iMax, i, j, k)%dzray * lRef, kind &
                    &= 4)
              case(10)
                field_prc(i, j) = real(ray(iMax, i, j, k)%dkray / lRef, kind &
                    &= 4)
              case(11)
                field_prc(i, j) = real(ray(iMax, i, j, k)%dlray / lRef, kind &
                    &= 4)
              case(12)
                field_prc(i, j) = real(ray(iMax, i, j, k)%dmray / lRef, kind &
                    &= 4)
              case(13)
                field_prc(i, j) = real(ray(iMax, i, j, k)%omega / tRef, kind &
                    &= 4)
              case(14)
                field_prc(i, j) = real(ray(iMax, i, j, k)%dens * rhoref * tRef &
                    &* uRef ** 2.0 / lRef ** dim, kind = 4)
              end select
            end do
            call mpi_gather(field_prc(1, j), nx, mpi_real, field_mst(1, j), &
                &nx, mpi_real, 0, comm, ierror)
          end do

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

            write(40, rec = irc_prc) field_out
          end if
        end do
      end do
    end do

    if(master) close(unit = 40)

    !--------------------------------------------------------------
    ! here, after the parallelization of MS-GWaM, should also be
    ! optional output for all ray-volume quantities (x,y,z,k,l,m,N)
    !--------------------------------------------------------------

    ! FJFeb2023
    ! set counter
    ! iOut = iOut + 1

    !--------------------------------------------------------------
    ! SD
    ! output ray volumes
    !--------------------------------------------------------------

    if(include_ice2) then

      TNOR = 12800 * 16 ! total number of rays for pir09
      NOR_prc = TNOR / nprocx / nprocy ! number of rays per processor

      if(master) then
        open(40, file = 'pf_norays.dat', form = "unformatted", access &
            &= 'direct', recl = sizeofreal4)
        open(44, file = 'pf_rays.dat', form = "unformatted", access &
            &= 'direct', recl = TNOR * NFR * sizeofreal4)
      end if

      !irc_prc = NOR * (iOut - 1)

      cRay = 0

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx

            NOR = nRay(i, j, k)

            if(NOR < 1) cycle

            do iRay = 1, NOR

              !count all rays
              cRay = cRay + 1

              ! vector x
              vct(1) = ray(iRay, i, j, k)%x * lRef
              vct(2) = ray(iRay, i, j, k)%y * lRef
              vct(3) = ray(iRay, i, j, k)%z * lRef

              ! vector k
              vct(4) = ray(iRay, i, j, k)%k / lRef
              vct(5) = ray(iRay, i, j, k)%l / lRef
              vct(6) = ray(iRay, i, j, k)%m / lRef

              ! vector dx
              vct(7) = ray(iRay, i, j, k)%dxray * lRef
              vct(8) = ray(iRay, i, j, k)%dyray * lRef
              vct(9) = ray(iRay, i, j, k)%dzray * lRef

              ! vector dk
              vct(10) = ray(iRay, i, j, k)%dkray / lRef
              vct(11) = ray(iRay, i, j, k)%dlray / lRef
              vct(12) = ray(iRay, i, j, k)%dmray / lRef

              vct(13) = ray(iRay, i, j, k)%omega / tRef
              vct(14) = ray(iRay, i, j, k)%dens

              ! vertical vel.
              vct(15) = opt_ray(iRay, i, j, k)%w

              vct_prc(cRay, 1:NFR) = real(vct, 4)

            end do !iRay
          end do ! i
        end do ! j
      end do ! k

      print *, 'PRC', is, x(is + 2) * lRef / 10 + 3

      do j = 1, NFR
        call mpi_gather(vct_prc(1, j), NOR_prc, mpi_real, vct_mst(1, j), &
            &NOR_prc, mpi_real, 0, comm, ierror)
      end do
      !!$
      !!$!     call mpi_gather(real(TNOR,4), 1, mpi_real, nor_mst(1), 1, &
      !!$!          mpi_real, 0, comm, ierror)
      !*        irc_prc = irc_prc + 1

      call mpi_barrier(comm, ierror)
      if(master) then

        print *, 'TTNOR', TNOR, NOR_prc, cRay, sum(NRay)

        !output fields
        write(40, rec = iOut) real(TNOR, 4)

        !*************************
        !TEST output
        !*************************
        !!$!*        write(*, *) 'iOut', iOut, TNOR
        !!$
        !!$!*        vct_mst(1, 1) = iOut*2 ! x-coord.
        !!$!*        vct_mst(1, 3) = iOut ! z-coord.
        !!$!*
        !!$!*        vct_mst(1, 7) = 1. ! dx
        !!$!*        vct_mst(1, 9) = .5 ! dz

        !*************************
        !END TEST
        !*************************

        !!$        do j = 1, TNOR
        !!$           write(44, rec=NoR_out + j) vct_mst(j, 1:NFR)
        !!$        end do
        !NoR_out = NoR_out + TNOR

        write(44, rec = NoR_out) ((vct_mst(j, i), i = 1, NFR), j = 1, TNOR)

        NoR_out = NoR_out + 1

        close(40)
        close(44)

        print *, 'output_wkb finished'

      end if ! master

    end if ! include_ice2

  end subroutine output_wkb

  !---------------------------------------------------------------------

  subroutine init_output
    !-----------------------
    ! allocate optVar field
    !-----------------------

    ! local variables
    integer :: allocstat

    ! divPu
    allocate(optVar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "output.f90/init_output: could not allocate &
        &optVar. Stop."

    ! for output global fields
    allocate(field_mst(sizeX * nprocy, ny), stat = allocstat)
    if(allocstat /= 0) stop "output.f90/init_output: could not allocate &
        &field_mst. Stop."
    allocate(field_out(sizeX, sizeY), stat = allocstat)
    if(allocstat /= 0) stop "output.f90/init_output: could not allocate &
        &field_out. Stop."
  end subroutine init_output
  !
  !--------------------------------------------------------------------------

  subroutine terminate_output
    !-------------------------
    ! deallocate optVar field and varIn, varOut and offset (edited by Niklas Ehlert, 20190128)
    !-------------------------

    ! local variables
    integer :: allocstat

    !---------------- deallocate variables -----------------------

    deallocate(optVar, stat = allocstat)
    if(allocstat /= 0) stop "output.f90/init_output: could not deallocate &
        &optVar. Stop."

  end subroutine terminate_output

  !------------------------------------------------------------------------
  ! output some singe field in time
  subroutine output_field(iOut, field, filename, scaling_coef)

    !-------------------------------
    !  writes data to file filename.dat
    !-------------------------------

    ! output counter
    integer, intent(in) :: iOut

    ! argument fields
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(in) :: field

    ! local and global output field and record nbs.
    real * 4, dimension(nx, ny) :: field_prc
    integer irc_prc, irc_out

    ! local variables
    integer :: i, j, k, iVar
    real :: time_dim

    ! needed for output to screen
    character(len = 20) :: fmt, form
    character(len = 40) :: cpuTimeChar

    ! CPU Time
    integer :: days, hours, mins, secs
    real :: timeVar

    ! hotBubble output
    real :: rho, theta_dim

    ! scale the variable
    real, intent(in) :: scaling_coef
    character(len = *) :: filename

    ! buoyancy
    real :: b, b_dim, theta

    integer :: i_prc, i_mst, i_out, j_prc, j_mst, j_out

    if(master) then ! modified by Junhong Wei (20161110)
      print *, ""
      print *, " Output into File "
      print *, filename

    end if ! modified by Junhong Wei (20161110)

    !------------------------------
    !   prepare output file
    !------------------------------

    ! open output file

    if(master) then
      open(21, file = filename, form = "unformatted", access = 'direct', recl &
          &= SizeX * SizeY * sizeofreal4)
    end if

    !---------------------------------------
    !       dimensionalising and layerwise output
    !---------------------------------------

    irc_prc = 1
    irc_prc = irc_prc * iOut * nz

    do k = 1, nz
      ! dimensionalization

      do j = 1, ny
        do i = 1, nx
          field_prc(i, j) = real(field(i, j, k) * scaling_coef, kind = 4)
        end do ! i
        call mpi_gather(field_prc(1, j), nx, mpi_real, field_mst(1, j), nx, &
            &mpi_real, 0, comm, ierror)
      end do ! j

      ! layerwise output

      irc_prc = irc_prc + 1
      !            write(40,rec=irc_prc) field_prc
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

        write(21, rec = irc_prc) field_out
      end if
    end do ! k

    !------------------------------------
    !              close file
    !------------------------------------

    if(master) close(unit = 21)

  end subroutine output_field

  !--------------------------------------------------------------------------
  ! output some singe field in time
  subroutine output_profile(iOut, field, filename)

    !-------------------------------
    !  writes data to file filename.dat
    !-------------------------------

    ! output counter
    integer, intent(inout) :: iOut

    ! argument fields
    !UAC real,dimension(-nbz:nz+nbz),intent(in) :: field
    real, dimension(- 1:nz + 2), intent(in) :: field

    ! local and global output field and record nbs.
    real * 4, dimension(nx, ny) :: field_prc
    integer irc_prc, irc_out

    ! local variables
    integer :: i, j, k, iVar
    real :: time_dim

    ! needed for output to screen
    character(len = 20) :: fmt, form
    character(len = 40) :: cpuTimeChar

    ! CPU Time
    integer :: days, hours, mins, secs
    real :: timeVar

    ! hotBubble output
    real :: rho, theta_dim

    character(len = *) :: filename

    ! buoyancy
    real :: b, b_dim, theta

    integer :: i_prc, i_mst, i_out, j_prc, j_mst, j_out

    !  iOut = iOut - 1

    ! FJFeb2023
    ! if(master) then ! modified by Junhong Wei (20161110)
    !   print *, ""
    !   print *, " Output into File "
    !   print *, filename

    ! end if ! modified by Junhong Wei (20161110)

    !------------------------------
    !   prepare output file
    !------------------------------

    ! open output file

    if(master) then
      open(21, file = filename, form = "unformatted", access = 'direct', recl &
          &= SizeX * SizeY * sizeofreal4)
    end if

    !---------------------------------------
    !       dimensionalising and layerwise output
    !---------------------------------------

    ! FJFeb2023
    ! irc_prc = 1
    ! irc_prc = irc_prc * iOut * (nz)
    irc_prc = (iOut - 1) * nz

    do k = 1, nz
      ! dimensionalization

      !testb
      ! if(master) print *, k, field(k)
      !teste

      do j = 1, ny
        do i = 1, nx

          !UAC
          !field_prc(i,j) = real(field(k)  &
          !            & , kind=4)
          field_prc(i, j) = field(k)
          !UAE
        end do ! i
        call mpi_gather(field_prc(1, j), nx, mpi_real, field_mst(1, j), nx, &
            &mpi_real, 0, comm, ierror)
      end do ! j

      ! layerwise output

      irc_prc = irc_prc + 1
      !            write(40,rec=irc_prc) field_prc
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

        write(21, rec = irc_prc) field_out
      end if
    end do ! k

    !------------------------------------
    !              close file
    !------------------------------------

    if(master) close(unit = 21)

    ! iOut = iOut + 1

  end subroutine output_profile

  ! output some singe field in time !FS for restart
  subroutine read_profile(iIn, field, filename)

    !-------------------------------
    !  writes data to file filename.dat
    !-------------------------------

    ! output counter
    integer, intent(in) :: iIn

    ! argument fields
    real, dimension(- nbz:nz + nbz), intent(out) :: field

    ! local and global output field and record nbs.
    real * 4, dimension(nx, ny) :: field_prc
    integer irc_prc, irc_out

    ! local variables
    integer :: i, j, k, iVar
    real :: time_dim

    ! needed for output to screen
    character(len = 20) :: fmt, form
    character(len = 40) :: cpuTimeChar

    ! CPU Time
    integer :: days, hours, mins, secs
    real :: timeVar

    ! hotBubble output
    real :: rho, theta_dim

    character(len = *) :: filename

    ! buoyancy
    real :: b, b_dim, theta

    integer :: i_prc, i_mst, i_out, j_prc, j_mst, j_out

    !  iOut = iOut - 1

    if(master) then
      print *, ""
      print *, " Input from File "
      print *, ""
      write(*, fmt = "(a25,i15)") " reading record no. ", iIn
    end if

    !------------------------------
    !   prepare output file
    !------------------------------

    ! open output file

    if(master) then
      open(21, file = filename, form = "unformatted", access = 'direct', recl &
          &= SizeX * SizeY * sizeofreal4)
      !      print*,"pf_all_in.dat opened"
    end if

    field = 0.
    !---------------------------------------
    !       dimensionalising and layerwise output
    !---------------------------------------

    irc_prc = 1
    irc_prc = irc_prc * iIn * (nz + 3)

    do k = - 1, nz + 2
      ! read data layerwise

      irc_prc = irc_prc + 1
      if(master) then

        read(21, rec = irc_prc) field_out
        do j = 1, ny
          j_mst = j

          do j_prc = 1, nprocy
            j_out = ny * (j_prc - 1) + j

            do i_prc = 1, nprocx
              do i = 1, nx
                i_out = nx * (i_prc - 1) + i

                i_mst = nprocy * nx * (i_prc - 1) + (j_prc - 1) * nx + i

                field_mst(i_mst, j_mst) = field_out(i_out, j_out)
              end do
            end do
          end do
        end do

      end if
      call mpi_barrier(comm, ierror)

      ! dimensionalization

      do j = 1, ny
        ! data distributed over all processors

        call mpi_scatter(field_mst(1, j), nx, mpi_real, field_prc(1, j), nx, &
            &mpi_real, 0, comm, ierror)
        do i = 1, nx

          field(k) = field_prc(i, j)
        end do ! i
      end do ! j
    end do ! k

    !------------------------------------
    !              close file
    !------------------------------------

    if(master) close(unit = 21)

    ! iOut = iOut + 1

  end subroutine read_profile

  subroutine output_fluxes(iOut, var, flux, iTime, time, cpuTime)

    !-------------------------------
    !  writes data to file pf_all.dat
    !-------------------------------

    ! output counter
    integer, intent(inout) :: iOut

    ! argument fields
    type(var_type), intent(in) :: var

    type(flux_type), intent(in) :: flux

    ! local and global output field and record nbs.
    real * 4, dimension(nx, ny) :: field_prc
    integer irc_prc, irc_out

    ! argument parameters
    integer, intent(in) :: iTime
    real, intent(in) :: time, cpuTime

    ! local variables
    integer :: i, j, k, iVar
    real :: time_dim

    ! needed for output to screen
    character(len = 20) :: fmt, form
    character(len = 40) :: cpuTimeChar

    ! CPU Time
    integer :: days, hours, mins, secs
    real :: timeVar

    ! hotBubble output
    real :: rho, theta_dim

    ! buoyancy
    real :: b, b_dim, theta

    integer :: i_prc, i_mst, i_out, j_prc, j_mst, j_out

    real, dimension(0:nz + 1) :: sum_local, sum_global !\bar rho

    sum_local = 0.
    sum_global = 0.
    if(topography) then
      do k = 0, nz + 1
        sum_local(k) = sum(var%rho(1:nx, 1:ny, k) + rhoStratTFC(1:nx, 1:ny, k))
      end do
    else
      do k = 0, nz + 1
        sum_local(k) = sum(var%rho(1:nx, 1:ny, k) + rhoStrat(k))
      end do
    end if

    ! global sum and average

    call mpi_allreduce(sum_local(0), sum_global(0), nz + 1 - 0 + 1, &
        &mpi_double_precision, mpi_sum, comm, ierror)
    sum_global = sum_global / (sizeX * sizeY)

    !start output routine

    time_dim = time * tRef

    if(master) then ! modified by Junhong Wei (20161110)
      print *, ""
      print *, " Output into File "
      print *, ""
      write(*, fmt = "(a25,i15)") " at time step = ", iTime
      write(*, fmt = "(a25,f15.1,a8)") " at physical time = ", time_dim, " &
          &seconds"
    end if ! modified by Junhong Wei (20161110)

    !------------------------------
    !   prepare output file
    !------------------------------

    ! open output file

    if(master) then
      open(43, file = 'fluxes.dat', form = "unformatted", access = 'direct', &
          &recl = SizeX * SizeY * sizeofreal4)
    end if

    ! calc cpu-time in days/hours/minutes/seconds
    timeVar = cpuTime
    days = floor(timeVar / 86400.0)
    timeVar = timeVar - 86400.0 * days
    hours = floor(timeVar / 3600.0)
    timeVar = timeVar - 3600.0 * hours
    mins = floor(timeVar / 60.0)
    timeVar = timeVar - 60.0 * mins
    secs = int(timeVar)

    write(unit = cpuTimeChar, fmt = "(i2,a,2(i2.2,a),i2.2)") days, " ", hours, &
        &":", mins, ":", secs

    if(master) write(*, fmt = "(a25,a25)") "CPU time = ", cpuTimeChar

    !---------------------------------------
    !       dimensionalising and layerwise output
    !---------------------------------------

    ! irc_prc = 7 * iOut * nz
    irc_prc = 10 * (iOut - 1) * nz

    do iVar = 1, 10
      !if(varOut(iVar) == 1) then
      do k = 1, nz
        ! dimensionalization

        do j = 1, ny
          do i = 1, nx
            select case(iVar)

            case(1) ! f rho v
              if(topography) then
                field_prc(i, j) = rhoRef * uRef * f_Coriolis_dim * 0.5 &
                    &* (var%rho(i, j, k) + rhoStratTFC(i, j, k) + var%rho(i &
                    &+ 1, j, k) + rhoStratTFC(i + 1, j, k)) * 0.25 * (var%v(i, &
                    &j, k) + var%v(i + 1, j, k) + var%v(i, j + 1, k) + var%v(i &
                    &+ 1, j + 1, k))
              else
                field_prc(i, j) = rhoRef * uRef * f_Coriolis_dim * (0.5 &
                    &* (var%rho(i, j, k) + var%rho(i + 1, j, k)) &
                    &+ rhoStrat(k)) * (0.25 * (var%v(i, j, k) + var%v(i + 1, &
                    &j, k) + var%v(i, j + 1, k) + var%v(i + 1, j + 1, k)))
              end if

            case(2) ! d/dy (rho v u)
              field_prc(i, j) = rhoRef * uRef * uRef * (flux%u(i, j, k, 2))

            case(3) ! d/dz (rho w u)
              field_prc(i, j) = rhoRef * uRef * uRef * (flux%u(i, j, k, 3))

            case(4) ! Fx
              if(topography) then
                field_prc(i, j) = rhoRef * uRef * kr_sp_tfc(i, j, k) / tRef &
                    &* 0.5 * (var%rho(i, j, k) + rhoStratTFC(i, j, k) &
                    &+ var%rho(i + 1, j, k) + rhoStratTFC(i + 1, j, k)) &
                    &* (var%u(i, j, k) - var_env%u(i, j, k))
              else
                field_prc(i, j) = rhoRef * uRef * (kv_hs(j, k) + kr_sp(j, k)) &
                    &/ tRef * (0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k)) &
                    &+ rhoStrat(k)) * (var%u(i, j, k) - var_env%u(i, j, k))
              end if

            case(5) ! f v
              field_prc(i, j) = uRef * f_Coriolis_dim * (0.25 * (var%v(i, j, &
                  &k) + var%v(i + 1, j, k) + var%v(i, j + 1, k) + var%v(i + 1, &
                  &j + 1, k)))

            case(6) ! rhoStrat v u
              field_prc(i, j) = 0.5 * (var%v(i + 1, j, k) + var%v(i, j, k)) &
                  &* 0.5 * (var%u(i, j + 1, k) + var%u(i, j, k)) * uRef * uRef &
                  &* sum_global(k) * rhoRef

            case(7) ! rhoStrat w u
              if(topography) then
                field_prc(i, j) = 0.5 * (sum_global(k) + sum_global(k + 1)) &
                    &* 0.5 * (vertWindTFC(i + 1, j, k, var) + vertWindTFC(i, &
                    &j, k, var)) * 0.5 * (var%u(i, j, k + 1) + var%u(i, j, k)) &
                    &* uRef * uRef * rhoRef
              else
                field_prc(i, j) = 0.5 * (sum_global(k) + sum_global(k + 1)) &
                    &* 0.5 * (var%w(i + 1, j, k) + var%w(i, j, k)) * 0.5 &
                    &* (var%u(i, j, k + 1) + var%u(i, j, k)) * uRef * uRef &
                    &* rhoRef
              end if

            case(8) ! zonal tracer flux (rho u chi)
              if(include_tracer) then
                field_prc(i, j) = rhoRef * uRef * flux%chi(i, j, k, 1)
              else
                field_prc(i, j) = 0.0
              end if

            case(9) ! meridional tracer flux (rho v chi)
              if(include_tracer) then
                field_prc(i, j) = rhoRef * uRef * flux%chi(i, j, k, 2)
              else
                field_prc(i, j) = 0.0
              end if

            case(10) ! vertical tracer flux (rho w chi)
              if(include_tracer) then
                field_prc(i, j) = rhoRef * uRef * flux%chi(i, j, k, 3)
              else
                field_prc(i, j) = 0.0
              end if

            case default
            end select ! iVar
          end do ! i
          call mpi_gather(field_prc(1, j), nx, mpi_real, field_mst(1, j), nx, &
              &mpi_real, 0, comm, ierror)
        end do ! j

        ! layerwise output

        irc_prc = irc_prc + 1
        !            write(40,rec=irc_prc) field_prc
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

          write(43, rec = irc_prc) field_out
        end if
      end do ! k
      !end if
    end do ! iVar

    !------------------------------------
    !              close file
    !------------------------------------

    if(master) close(unit = 43)

    ! set counter
    !iOut = iOut + 1

  end subroutine output_fluxes

  !-----------------------------------------------------------------------------

  subroutine output_background(iOut)

    integer, intent(in) :: iOut
    character(len = 20) :: filename
    real * 4, dimension(nx, ny) :: field_prc
    real * 4 :: field
    integer :: i_out, i_mst, i_prc, j_out, j_mst, j_prc
    integer :: irc_prc
    integer :: i, j, k
    integer :: iVar

    do iVar = 1, 5
      ! Set filename.
      if(iVar == 1) then
        filename = "pStrat.dat"
      else if(iVar == 2) then
        filename = "thetaStrat.dat"
      else if(iVar == 3) then
        filename = "rhoStrat.dat"
      else if(iVar == 4) then
        filename = "bvsStrat.dat"
      end if

      irc_prc = (iOut - 1) * nz

      if(topography) then
        ! Open file.
        if(master) then
          open(42, file = filename, form = "unformatted", access = "direct", &
              &recl = sizeX * sizeY * sizeofreal4)
        end if

        do k = 1, nz
          do j = 1, ny
            ! Dimensionalize.
            do i = 1, nx
              if(iVar == 1) then
                field_prc(i, j) = pStratTFC(i, j, k) * rhoRef * thetaRef
              else if(iVar == 2) then
                field_prc(i, j) = thetaStratTFC(i, j, k) * thetaRef
              else if(iVar == 3) then
                field_prc(i, j) = rhoStratTFC(i, j, k) * rhoRef
              else if(iVar == 4) then
                field_prc(i, j) = bvsStratTFC(i, j, k) / tRef ** 2.0
              end if
            end do

            ! Distribute data over all processors.
            call mpi_gather(field_prc(1, j), nx, mpi_real, field_mst(1, j), &
                &nx, mpi_real, 0, comm, ierror)
          end do

          irc_prc = irc_prc + 1

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
                    field_out(i_out, j_out) = field_mst(i_mst, j_mst)
                  end do
                end do
              end do
            end do
            write(42, rec = irc_prc) field_out
          end if
        end do

        ! Close file.
        if(master) close(unit = 42)
      else
        if(master) then
          ! Open file.
          open(42, file = filename, form = "unformatted", access = "direct", &
              &recl = sizeofreal4)

          do k = 1, nz
            ! Dimensionalize.
            if(iVar == 1) then
              field = pStrat(k) * rhoRef * thetaRef
            else if(iVar == 2) then
              field = thetaStrat(k) * thetaRef
            else if(iVar == 3) then
              field = rhoStrat(k) * rhoRef
            else if(iVar == 4) then
              field = bvsStrat(k) / tRef ** 2.0
            end if

            irc_prc = irc_prc + 1

            ! Output layerwise.
            write(42, rec = irc_prc) field
          end do

          ! Close file.
          close(unit = 42)
        end if
      end if
    end do

  end subroutine output_background

end module output_module
