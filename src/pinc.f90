program pinc_prog

  !----------------------------------------------------------------------------!
  !                               PinCFlow                                     !
  !----------------------------------------------------------------------------!
  ! Pseudo-incompressible flow solver developed by Rieper et al. (2013),       !
  ! Wilhelm et al. (2018), Wei et al. (2019), Schmid et al. (2021), Jochum et  !
  ! al. (2025) and many others. This is a condensed variant of the main code.  !
  !----------------------------------------------------------------------------!

  use type_module
  use mpi_module
  use timeScheme_module
  use init_module
  use atmosphere_module
  use boundary_module
  use flux_module
  use update_module
  use poisson_module
  use finish_module
  use bicgstab_tools_module
  use mpi
  use output_netCDF_module

  implicit none

  integer :: iTime
  real :: time, dt

  ! CPU Time
  integer :: rate, startTimeCount, timeCount
  real :: cpuTime

  ! MPI stuff
  logical :: error_flag

  ! Variable fields
  type(var_type) :: var, var0, var1

  ! RK tendencies
  real, dimension(:, :, :), allocatable :: dRho, dRhop
  real, dimension(:, :, :, :), allocatable :: dMom

  ! Flux fields
  type(flux_type) :: flux, flux0

  ! Output per timeStep
  logical :: output
  real :: nextOutputTime ! scaled time for next output

  ! Indices
  integer :: i, j, k
  integer :: ix, jy, kz

  ! error handling and controlling
  logical :: errFlagBiCG, errFlagTStep
  integer :: nIterBicg, nTotalBicg
  real :: nAverageBicg
  character(len = 8) :: systemDate
  character(len = 10) :: fancySystemDate
  character(len = 10) :: systemTime
  character(len = 5) :: fancySystemTime

  real :: tolpoisson_s

  integer :: allocstat

  integer :: i00, j00
  real :: spongeAlphaZ, spongeAlphaY, spongeAlphaX
  real :: height

  !-------------------------------------------------
  !                    Set up
  !-------------------------------------------------

  file_namelist = 'input.f90'
  if(command_argument_count() /= 0) call get_command_argument(1, file_namelist)

  ! init counter and time
  iOut = 0
  iTime = 0; time = 0.0; cpuTime = 0.0

  call init_mpi(error_flag)

  if(error_flag) then
    print *, "error in the init_mpi"
    goto 666
  end if

  if(master) then
    call date_and_time(date = systemDate, time = systemTime)

    fancySystemDate = systemDate(1:4) // "-" // systemDate(5:6) // "-" &
        &// systemDate(7:8)
    fancySystemTime = systemTime(1:2) // ":" // systemTime(3:4)

    print "(a)", ""
    print "(a)", repeat("-", 80)
    print "(36x, a)", "PincFlow"
    print "(12x, a)", "developed by Rieper et al (2013) and Schmid et al (2021)"
    print "(28x, a)", "modified by many others"
    print "(a)", repeat("-", 80)
    print "(a)", ""
    print "(a)", "Date: " // fancySystemDate
    print "(a)", "Time: " // fancySystemTime
    print "(a)", ""
    print "(a)", "Virtual topology: [idim, jdim] = [" // trim_integer(idim) &
        &// ", " // trim_integer(jdim) // "]"
    print "(a)", ""

    call system_clock(count_rate = rate)
    call system_clock(count = startTimeCount)
  end if

  !-------------------------------------------------
  !                    Set up
  !-------------------------------------------------

  ! init counter and time
  nTotalBicg = 0

  ! 1) allocate variables
  ! 2) read input.f90
  call setup(var, var0, var1, flux, flux0, dRho, dRhop, dMom)

  ! Read topography.
  if(mountain_case == 0) then
    call read_topography_netCDF(iIn)
  end if

  call init_atmosphere ! set atmospheric background state

  call initialise(var, flux) ! set initial conditions

  call SetUpBiCGStab ! Set BiCGStab arrays

  call init_fluxes ! allocate reconstructed variables
  call init_timeScheme ! define Runge-Kutta parameters
  call init_poisson ! allocate dp

  !---------------------------------------------
  !        Initial divergence cleaning
  !---------------------------------------------

  if(initialCleaning) then
    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    call reset_flux_type(flux)
    dMom = 0.

    tolpoisson_s = tolPoisson
    tolPoisson = 1.e-8

    call Corrector(var, flux, 1.0, errFlagBicg, nIterBicg, "expl", 1., 1.)

    if(errFlagBicg) stop

    tolPoisson = tolpoisson_s
  end if

  !-------------------------------------------------
  !              Read initial data
  !-------------------------------------------------

  if(restart) then
    if(master) then
      print *, "reading restart files"
    end if

    call read_netCDF(iIn, var, time = time)

    if(maxTime < time * tRef) stop "restart error: maxTime < current time"

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")
  end if

  !------------------------------------------
  !              Initial output
  !------------------------------------------

  ! create netCDF file pincflow_data_out.nc
  call create_netCDF

  ! write initial state to netCDF file
  call write_netCDF(iOut, iTime, time, cpuTime, var)

  output = .false.
  nextOutputTime = time * tRef + outputTimeDiff

  !-----------------------------------------------------
  !                        Time loop
  !-----------------------------------------------------

  if(master) then
    print *, "...starting time loop"
  end if

  if(outputType == "time") maxIter = 2 ** 30

  time_loop: do iTime = 1, maxIter

    if(master) then
      print *, ""
      print "(a)", repeat("-", 80)
      write(*, fmt = "(a25,i15)") " Time step = ", iTime
      write(*, fmt = "(a25,f15.1,a8)") " Time = ", time * tRef, "seconds"
      print "(a)", repeat("-", 80)
    end if

    !----------------------------------
    !         Calc time step
    !----------------------------------

    call timestep(var, dt, errFlagTstep)

    ! abort if time step too small
    if(dt * tRef < dtMin_dim) then
      if(master) then
        print *, " TimeStep routine: "
        write(*, fmt = "(a25,es15.4,a8)") "dt < dtMin!!! dtMin = ", dtMin_dim, &
            &"seconds"
        write(*, fmt = "(a25,es15.4,a8)") "dt * tRef = ", dt * tRef, "seconds"

        call system_clock(count = timeCount)
        cpuTime = (timeCount - startTimeCount) / real(rate)
      end if

      ! final output
      call write_netCDF(iOut, iTime, time, cpuTime, var)

      call close_netCDF

      call mpi_barrier(comm, ierror)
      call mpi_finalize(ierror)
      print *, "pinc.f90: time step too small. dt = ", dt * tRef
      stop
    end if

    ! correct dt to hit desired output time for outputType 'time'
    if(outputType == 'time') then
      if((time + dt) * tRef + dtMin_dim > nextOutputTime) then
        dt = nextOutputTime / tRef - time
        output = .true.
        if(master) then
          write(*, fmt = "(a25,es15.4,a8)") "dt for output = ", dt * tRef, &
              &"seconds"
        end if
      end if
    end if

    ! advance time
    time = time + dt

    !-----------------------------------------------------------------
    !                         Sponge layer
    !-----------------------------------------------------------------

    if(spongeLayer) then
      if(unifiedSponge) then
        kr_sp_tfc = 0.0
        kr_sp_w_tfc = 0.0

        spongeAlphaZ = spongeAlphaZ_dim * tRef

        if(lateralSponge) then
          i00 = is + nbx - 1
          j00 = js + nby - 1
          if(sizeX > 1 .and. sizeY > 1) then
            spongeAlphaZ = spongeAlphaZ / 3.0
            spongeAlphaX = spongeAlphaZ
            spongeAlphaY = spongeAlphaZ
          else if(sizeX > 1) then
            spongeAlphaZ = spongeAlphaZ / 2.0
            spongeAlphaX = spongeAlphaZ
            spongeAlphaY = 0.0
          else if(sizeY > 1) then
            spongeAlphaZ = spongeAlphaZ / 2.0
            spongeAlphaX = 0.0
            spongeAlphaY = spongeAlphaZ
          end if
        end if

        alphaUnifiedSponge = 0.0

        do k = 1, nz
          do j = 0, ny + 1
            do i = 0, nx + 1
              height = zTFC(i, j, k)

              if(spongeType == "exponential") then
                if(sizeZ > 1) then
                  alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, k) &
                      &+ spongeAlphaZ * exp((height - lz(1)) / dzSponge)
                end if
                if(lateralSponge) then
                  if(sizeX > 1) then
                    if(x(i00 + i) <= 0.5 * (lx(0) + lx(1))) then
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + spongeAlphaX * exp((lx(0) - x(i00 + i)) &
                          &/ dxSponge)
                    else
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + spongeAlphaX * exp((x(i00 + i) - lx(1)) &
                          &/ dxSponge)
                    end if
                  end if
                  if(sizeY > 1) then
                    if(y(j00 + j) <= 0.5 * (ly(0) + ly(1))) then
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + spongeAlphaY * exp((ly(0) - y(j00 + j)) &
                          &/ dySponge)
                    else
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + spongeAlphaY * exp((y(j00 + j) - ly(1)) &
                          &/ dySponge)
                    end if
                  end if
                end if

              else if(spongeType == "cosmo") then
                if(sizeZ > 1) then
                  if(height >= zSponge) then
                    alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, k) &
                        &+ 0.5 / cosmoSteps / dt * (1.0 - cos(pi * (height &
                        &- zSponge) / dzSponge))
                  end if
                end if
                if(lateralSponge) then
                  if(sizeX > 1) then
                    if(x(i00 + i) <= xSponge0) then
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + 0.5 / cosmoSteps / dt * (1.0 - cos(pi &
                          &* (xSponge0 - x(i00 + i)) / dxSponge))
                    else if(x(i00 + i) >= xSponge1) then
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + 0.5 / cosmoSteps / dt * (1.0 - cos(pi * (x(i00 &
                          &+ i) - xSponge1) / dxSponge))
                    end if
                  end if
                  if(sizeY > 1) then
                    if(y(j00 + j) <= ySponge0) then
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + 0.5 / cosmoSteps / dt * (1.0 - cos(pi &
                          &* (ySponge0 - y(j00 + j)) / dySponge))
                    else if(y(j00 + j) >= ySponge1) then
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + 0.5 / cosmoSteps / dt * (1.0 - cos(pi * (y(j00 &
                          &+ j) - ySponge1) / dySponge))
                    end if
                  end if
                end if

              else if(spongeType == "polynomial") then
                if(sizeZ > 1) then
                  if(height >= zSponge) then
                    alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, k) &
                        &+ spongeAlphaZ * ((height - zSponge) / dzSponge) &
                        &** spongeOrder
                  end if
                end if
                if(lateralSponge) then
                  if(sizeX > 1) then
                    if(x(i00 + i) <= xSponge0) then
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + spongeAlphaX * ((xSponge0 - x(i00 + i)) &
                          &/ dxSponge) ** spongeOrder
                    else if(x(i00 + i) >= xSponge1) then
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + spongeAlphaX * ((x(i00 + i) - xSponge1) &
                          &/ dxSponge) ** spongeOrder
                    end if
                  end if
                  if(sizeY > 1) then
                    if(y(j00 + j) <= ySponge0) then
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + spongeAlphaY * ((ySponge0 - y(j00 + j)) &
                          &/ dySponge) ** spongeOrder
                    else if(y(j00 + j) >= ySponge1) then
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + spongeAlphaY * ((y(j00 + j) - ySponge1) &
                          &/ dySponge) ** spongeOrder
                    end if
                  end if
                end if

              else if(spongeType == "sinusoidal") then
                if(sizeZ > 1) then
                  if(height >= zSponge) then
                    alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, k) &
                        &+ spongeAlphaZ * sin(0.5 * pi * (height - zSponge) &
                        &/ dzSponge) ** 2.0
                  end if
                end if
                if(lateralSponge) then
                  if(sizeX > 1) then
                    if(x(i00 + i) <= xSponge0) then
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + spongeAlphaX * sin(0.5 * pi * (xSponge0 &
                          &- x(i00 + i)) / dxSponge) ** 2.0
                    else if(x(i00 + i) >= xSponge1) then
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + spongeAlphaX * sin(0.5 * pi * (x(i00 + i) &
                          &- xSponge1) / dxSponge) ** 2.0
                    end if
                  end if
                  if(sizeY > 1) then
                    if(y(j00 + j) <= ySponge0) then
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + spongeAlphaY * sin(0.5 * pi * (ySponge0 &
                          &- y(j00 + j)) / dySponge) ** 2.0
                    else if(y(j00 + j) >= ySponge1) then
                      alphaUnifiedSponge(i, j, k) = alphaUnifiedSponge(i, j, &
                          &k) + spongeAlphaY * sin(0.5 * pi * (y(j00 + j) &
                          &- ySponge1) / dySponge) ** 2.0
                    end if
                  end if
                end if
              end if
            end do
          end do
        end do

        alphaUnifiedSponge(:, :, 0) = alphaUnifiedSponge(:, :, 1)
        alphaUnifiedSponge(:, :, nz + 1) = alphaUnifiedSponge(:, :, nz)
      else
        alpspg = spongeAlphaZ_fac / dt

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              if(zTFC(i, j, k) >= zSponge) then
                kr_sp_tfc(i, j, k) = alpspg * sin(0.5 * pi * (zTFC(i, j, k) &
                    &- zSponge) / (lz(1) - zSponge)) ** 2.0
                kr_sp_w_tfc(i, j, k) = kr_sp_tfc(i, j, k) / jac(i, j, k)
              end if
            end do
          end do
        end do

        kr_sp_tfc(:, :, 0) = kr_sp_tfc(:, :, 1)
        kr_sp_tfc(:, :, nz + 1) = kr_sp_tfc(:, :, nz)
        kr_sp_w_tfc(:, :, 0) = kr_sp_w_tfc(:, :, 1)
        kr_sp_w_tfc(:, :, nz + 1) = kr_sp_w_tfc(:, :, nz)
      end if
    end if

    !---------------------------------------------------------------
    !                   Semi-implicit time scheme
    !---------------------------------------------------------------

    ! Synchronization of density fluctuations
    select case(model)
    case("pseudo_incompressible")
      var%rhop(:, :, :) = var%rho(:, :, :)
    case default
    end select

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    ! put initial state into var0 in order to save the advecting
    ! velocities
    var0 = var

    ! (1) explicit integration of convective and
    !     viscous-diffusive/turbulent fluxes over half a time step,
    !     with the advection velocity kept constant
    !     \psi^# = \psi^n + A^{dt/2} (\psi^n, v^n)

    if(master) print *, 'beginning a semi-implicit time step'
    if(master) print *, '(1) explicit integration lhs over dt/2'

    do RKstage = 1, nStages
      call setHalos(var, "var")
      call setBoundary(var, flux, "var")

      ! Reconstruction
      call reconstruction(var, "rho")
      call reconstruction(var, "rhop")
      call reconstruction(var, "uvw")

      ! Fluxes
      call massFlux(var0, var, flux)
      call momentumFlux(var0, var, flux)

      call setBoundary(var, flux, "flux")

      ! store initial flux
      if(RKstage == 1) flux0 = flux

      ! RK step for density and density fluctuations

      rhoOld = var%rho(:, :, :) ! rhoOld for momentum predictor

      call massUpdate(var, flux, 0.5 * dt, dRho, RKstage, "rho", "tot", &
          &"expl", 1.)

      call applyUnifiedSponge(var, stepFrac(RKstage) * 0.5 * dt, time, "rho")

      call massUpdate(var, flux, 0.5 * dt, dRhop, RKstage, "rhop", "lhs", &
          &"expl", 1.)

      call applyUnifiedSponge(var, stepFrac(RKstage) * 0.5 * dt, time, "rhop")

      ! RK step for momentum

      call setHalos(var, "var")
      call setBoundary(var, flux, "var")

      call momentumPredictor(var, flux, 0.5 * dt, dMom, RKstage, "lhs", &
          &"expl", 1.)

      call applyUnifiedSponge(var, stepFrac(RKstage) * 0.5 * dt, time, "uvw")

    end do

    ! (2) implicit integration of the linear right-hand sides of the
    !     equations for density fluctuations and momentum over half a
    !     time step, under consideration of the divergence constraint
    !     \psi^{n+1/2} = \psi^# + dt/2 Q(\psi^{n+1/2})

    if(master) print *, '(2) implicit integration rhs over dt/2'

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    ! use initial flux for update of reference atmosphere and w0
    flux = flux0

    ! uStar and vStar are needed for update of density fluctuations,
    ! therefore w is stored instead of rhop

    wOldTFC = var%w(:, :, :)

    ! update winds (uStar, vStar, wStar)

    call momentumPredictor(var, flux, 0.5 * dt, dMom, RKstage, "rhs", "impl", &
        &1.)

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    ! update density fluctuations (rhopStar)

    call massUpdate(var, flux, 0.5 * dt, dRhop, RKstage, "rhop", "rhs", &
        &"impl", 1.)

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    ! Correct momentum and density fluctuations
    call Corrector(var, flux, 0.5 * dt, errFlagBicg, nIterBicg, "impl", 1., 1.)

    if(errFlagBicg) then
      call write_netCDF(iOut, iTime, time, cpuTime, var)
      call close_netCDF

      if(master) then
        print *, 'output last state into record', iOut
      end if
      stop
    end if

    nTotalBicg = nTotalBicg + nIterBicg

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    ! put new state into var1 in order to save the advection velocities
    var1 = var

    ! (3) explicit integration of the linear right-hand sides of the
    !     equations for density fluctuations and momentum over half a
    !     time step, under consideration of the divergence constraint
    !     \psi^\ast = \psi^n + dt/2 Q(\psi^n)

    if(master) print *, '(3) explicit integration rhs over dt/2'

    ! (3) uses updated pressure field and (5) adjusts pressure over half a
    ! time step!
    var%rho(:, :, :) = var0%rho(:, :, :)
    var%u(:, :, :) = var0%u(:, :, :)
    var%v(:, :, :) = var0%v(:, :, :)
    var%w(:, :, :) = var0%w(:, :, :)
    var%rhop(:, :, :) = var0%rhop(:, :, :)

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    rhopOld = var%rhop(:, :, :) ! rhopOld for momentum predictor

    ! update density fluctuations (rhopStar)
    call massUpdate(var, flux, 0.5 * dt, dRhop, RKstage, "rhop", "rhs", &
        &"expl", 1.)

    ! update winds (uStar, vStar, wStar)

    call momentumPredictor(var, flux, 0.5 * dt, dMom, RKstage, "rhs", "expl", &
        &1.)

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    ! (4) explicit integration of convective and
    !     viscous-diffusive/turbulent fluxes over a full time step,
    !     with the advection velocity kept constant
    !     \psi^{\ast\ast} = \psi^\ast + A^dt (\psi^\ast, v^{n+1/2})

    if(master) print *, '(4) explicit integration lhs over dt'

    var0 = var1

    do RKstage = 1, nStages
      call setHalos(var, "var")
      call setBoundary(var, flux, "var")

      ! Reconstruction
      call reconstruction(var, "rho")
      call reconstruction(var, "rhop")
      call reconstruction(var, "uvw")

      ! Fluxes
      call massFlux(var0, var, flux)
      call momentumFlux(var0, var, flux)

      call setBoundary(var, flux, "flux")

      ! RK step for density and density fluctuations

      rhoOld = var%rho(:, :, :) ! rhoOld for momentum predictor

      call massUpdate(var, flux, dt, dRho, RKstage, "rho", "tot", "expl", 1.)

      call applyUnifiedSponge(var, stepFrac(RKstage) * dt, time, "rho")

      call massUpdate(var, flux, dt, dRhop, RKstage, "rhop", "lhs", "expl", 1.)

      call applyUnifiedSponge(var, stepFrac(RKstage) * dt, time, "rhop")

      ! RK step for momentum
      call setHalos(var, "var")
      call setBoundary(var, flux, "var")

      call momentumPredictor(var, flux, dt, dMom, RKstage, "lhs", "expl", 1.)

      call applyUnifiedSponge(var, stepFrac(RKstage) * dt, time, "uvw")
    end do

    ! (5) implicit integration of the linear right-hand sides of the
    !     equations for density fluctuations and momentum over half a
    !     time step, under consideration of the divergence constraint
    !     \psi^{n+1} = \psi^{\ast\ast} + dt/2 Q(\psi^{n+1})

    if(master) print *, '(5) implicit integration rhs over dt/2'

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    ! use initial flux for update of reference atmosphere and w0
    flux = flux0

    ! uStar and vStar are needed for update of density fluctuations,
    ! therefore w is stored instead of rhop
    wOldTFC = var%w(:, :, :)

    ! update winds (uStar, vStar, wStar)

    call momentumPredictor(var, flux, 0.5 * dt, dMom, RKstage, "rhs", "impl", &
        &2.)

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    ! update density fluctuations (rhopStar)

    call massUpdate(var, flux, 0.5 * dt, dRhop, RKstage, "rhop", "rhs", &
        &"impl", 2.)

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    ! (3) uses updated pressure field and (5) adjusts pressure over half a
    ! time step!
    call Corrector(var, flux, 0.5 * dt, errFlagBicg, nIterBicg, "impl", 2., 1.) ! pressure update over dt/2

    if(errFlagBicg) then
      call write_netCDF(iOut, iTime, time, cpuTime, var)
      call close_netCDF

      if(master) then
        print *, 'output last state into record', iOut
      end if
      stop
    end if

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    nTotalBicg = nTotalBicg + nIterBicg

    if(master) print *, 'semi-implicit time step done'

    !--------------------------------------------------------------
    !                           Output
    !--------------------------------------------------------------

    select case(outputType)
    case('time')
      if(output) then
        if(master) then
          call system_clock(count = timeCount)
          cpuTime = (timeCount - startTimeCount) / real(rate)
        end if

        call write_netCDF(iOut, iTime, time, cpuTime, var)
        output = .false.
        nextOutputTime = nextOutputTime + outputTimeDiff
        if(nextOutputTime >= maxTime) nextOutputTime = maxTime
      end if
    case('timeStep')
      if(modulo(iTime, nOutput) == 0) then
        if(master) then
          call system_clock(count = timeCount)
          cpuTime = (timeCount - startTimeCount) / real(rate)
        end if

        call write_netCDF(iOut, iTime, time, cpuTime, var)
      end if
    case default
      stop "main: unknown outputType"
    end select

    !-------------------------------------------
    !              Abort criteria
    !-------------------------------------------

    if(outputType == "time") then
      if(time * tRef >= maxTime) then
        if(master) then
          nAverageBicg = real(nTotalBicg) / real(iTime) / 2.0

          call system_clock(count = timeCount)
          cpuTime = (timeCount - startTimeCount) / real(rate)

          print *, ""
          print "(a)", repeat("-", 80)
          print *, " Pinc: Resume"
          print *, ""
          write(*, fmt = "(a31,es15.1)") "average Poisson iterations = ", &
              &nAverageBicg
          print "(a)", repeat("-", 80)
          print *, ""
        end if

        exit ! leave time_loop
      end if
    end if

  end do time_loop

  !-------------------------------------------
  !      Final output for timeStep
  !-------------------------------------------

  call close_netCDF

  if(master) then
    if(outputType == "timeStep") then
      nAverageBicg = real(nTotalBicg) / real(iTime - 1) / nStages

      call system_clock(count = timeCount)

      cpuTime = (timeCount - startTimeCount) / real(rate)

      print *, ""
      print "(a)", repeat("-", 80)
      print *, " Pinc: Resume"
      print *, ""
      write(*, fmt = "(a31,f15.1)") "average Poisson iterations = ", &
          &nAverageBicg
      if(preconditioner /= "no") write(*, fmt = "(a25,es15.1)") "ADI: dtau &
          &= ", dtau
      print "(a)", repeat("-", 80)
      print *, ""
    end if

  end if

  !-------------------------------------
  !       Deallocate variables
  !-------------------------------------

  call terminate_fluxes
  call terminate_poisson
  call terminate(var, var0, var1, flux, flux0, dRho, dRhop, dMom)
  call terminate_atmosphere

  if(master) then
    ! do nothing
    ! call CleanUpBiCGSTab ! Clean Up BiCGSTAB arrays
  else
    call CleanUpBiCGSTab ! Clean Up BiCGSTAB arrays
  end if

  666 if(master) then
    write(*, "(a)") ""
    write(*, "(a)") repeat("-", 80)
    write(*, "(a)") repeat(" ", 32) // "PincFlow finished" // repeat(" ", 33)
    write(*, "(a)") repeat("-", 80)
  end if

  call mpi_finalize(ierror)

  contains

  subroutine test(var, flux, variables)

    type(var_type), intent(in) :: var
    type(flux_type), intent(in) :: flux
    character(len = *) :: variables
    real :: local_sum, global_sum

    select case(variables)
    case("predictands")
      local_sum = sum(var%rho * rhoRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(rho) = ", global_sum

      local_sum = sum(var%rhop * rhoRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(rhop) = ", global_sum

      local_sum = sum(var%u * uRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(u) = ", global_sum

      local_sum = sum(var%v * uRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(v) = ", global_sum

      local_sum = sum(var%w * uRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(w) = ", global_sum

      local_sum = sum(var%pi)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(pip) = ", global_sum
    case("reconstructions")
      local_sum = sum(rhotilde * rhoRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(rhotilde) = ", global_sum

      local_sum = sum(rhoptilde * rhoRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(rhoptilde) = ", global_sum

      local_sum = sum(utilde * uRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(utilde) = ", global_sum

      local_sum = sum(vtilde * uRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(vtilde) = ", global_sum

      local_sum = sum(wtilde * uRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(wtilde) = ", global_sum
    case("fluxes")
      local_sum = sum(flux%rho * rhoRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(phirho) = ", global_sum

      local_sum = sum(flux%rhop * rhoRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(phirhop) = ", global_sum

      local_sum = sum(flux%u * uRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(phiu) = ", global_sum

      local_sum = sum(flux%v * uRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(phiv) = ", global_sum

      local_sum = sum(flux%w * uRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(phiw) = ", global_sum
    case("atmosphere")
      local_sum = sum(pStratTFC * rhoRef * thetaRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(pstrattfc) = ", global_sum

      local_sum = sum(thetaStratTFC * thetaRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(thetastrattfc) = ", global_sum

      local_sum = sum(rhoStratTFC * rhoRef)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(rhostrattfc) = ", global_sum

      local_sum = sum(bvsStratTFC / tRef ** 2)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(bvsstrattfc) = ", global_sum
    case("grid")
      if(master) then
        print *, "lx = ", lx
        print *, "ly = ", ly
        print *, "lz = ", lz
        print *, "dx = ", dx
        print *, "dy = ", dy
        print *, "dz = ", dz
        print *, "sum(x) = ", sum(x)
        print *, "sum(y) = ", sum(y)
        print *, "sum(z) = ", sum(z)
        print *, "sum(zs) = ", sum(zs)
        print *, "sum(ztildes) = ", sum(ztildes)
      end if

      local_sum = sum(topography_surface)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(topography_surface) = ", global_sum

      local_sum = sum(jac)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(jac) = ", global_sum

      local_sum = sum(met(:, :, :, 1, 3))
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(met(:, :, :, 1, 3)) = ", global_sum

      local_sum = sum(met(:, :, :, 2, 3))
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(met(:, :, :, 2, 3)) = ", global_sum

      local_sum = sum(met(:, :, :, 3, 3))
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(met(:, :, :, 3, 3)) = ", global_sum

      local_sum = sum(ztfc)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(ztfc) = ", global_sum

      local_sum = sum(ztildetfc)
      call mpi_allreduce(local_sum, global_sum, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      if(master) print *, "sum(ztildetfc) = ", global_sum
    end select

    stop

  end subroutine test

end program pinc_prog
