program pinc_prog

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !% calculates pseudo-incompressible model     %
  !% by Durran with ILES by Adams and Hickel    %
  !% cylflow (c) 2009 Felix Rieper              %
  !% latest revision: January 2010              %
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  use sizeof_module
  use bicgstab_tools_module
  use mpi
  use output_netCDF_module

  !----------------------------------------------------
  implicit none
  !---------------------------------------------------

  integer :: iTime
  real :: time, dt

  ! CPU Time
  integer :: rate, startTimeCount, timeCount
  real :: cpuTime

  ! MPI stuff
  logical :: error_flag
  real :: dt_local

  ! fields
  type(var_type) :: var, var0, var1, varG, source
  ! varG contains balance parts of current time step

  real, dimension(:, :, :), allocatable :: dRho, dRhop ! RK-Update for rho
  real, dimension(:, :, :, :), allocatable :: dMom ! RK for rhoU,rhoV,rhoW
  real, dimension(:, :, :), allocatable :: dTheta ! RK-Update for theta

  type(flux_type) :: flux, flux0

  ! output per timeStep
  logical :: output
  real :: nextOutputTime ! scaled time for next output

  ! general
  integer :: i, j, k, l
  integer :: ix, jy, kz

  ! restart
  ! logical :: scale

  ! classical RK
  real :: dt_Poisson

  ! error handling and controlling
  logical :: errFlagBiCG, errFlagTStep
  integer :: nIterBicg, nTotalBicg
  real :: nAverageBicg
  character(len = 8) :: systemDate
  character(len = 10) :: fancySystemDate
  character(len = 10) :: systemTime
  character(len = 5) :: fancySystemTime
  real :: a2

  real :: tolpoisson_s

  integer :: iVar

  real, dimension(:), allocatable :: alpbls
  integer :: allocstat

  integer :: j00
  real :: ymax, ymin, yloc, bla, wmax, wmax_loc

  ! TFC FJ
  integer :: i00
  real :: spongeAlphaZ, spongeAlphaY, spongeAlphaX

  real :: height

  !-------------------------------------------------
  !                    Set up
  !-------------------------------------------------

  file_namelist = 'input.f90' ! default
  if(command_argument_count() /= 0) call get_command_argument(1, file_namelist) ! take the given one

  ! init counter and time
  iOut = 0
  iTime = 0; time = 0.0; cpuTime = 0.0

  call getsize

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
  call setup(var, var0, var1, varG, flux, flux0, force, source, dRho, dRhop, &
      &dMom, dTheta)

  ! Read topography.
  if((topography) .and. mountain_case &
      &== 0) then
    call read_topography_netCDF(iIn)
  end if

  call init_atmosphere ! set atmospheric background state

  call initialise(var, flux) ! set initial conditions

  call SetUpBiCGStab ! Set BiCGStab arrays

  call init_fluxes ! allocate tilde variables
  call init_update
  call init_timeScheme ! define Runge-Kutta parameters
  call init_poisson ! allocate dp

  if(zero_initial_state) then
    call reset_var_type(var)
  end if

  !---------------------------------------------
  !        Initial divergence cleaning
  !---------------------------------------------

  if(initialCleaning) then
    call setHalos(var, "var")
    if(include_tracer) call setHalos(var, "tracer")
    call setBoundary(var, flux, "var")

    ! 1) allocate variables
    ! 2) read the namelist
    call reset_flux_type(flux)
    dMom = 0.

    tolpoisson_s = tolPoisson
    tolPoisson = 1.e-8

    call Corrector(var, flux, dMom, 1.0, errFlagBicg, nIterBicg, 1, "expl", &
        &1., 1.)
  
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
    if(include_tracer) call setHalos(var, "tracer")
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

  if(spongeLayer) then
    allocate(alpbls(0:ny + 1), stat = allocstat)
    if(allocstat /= 0) stop "pinc.f90: could not allocate alpbls"
  end if

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
        write(*, fmt = "(a25,i2.2)") "cycle paramLoop. iPara = ", iParam

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

    time = time + dt

    !-----------------------------------------------------------------
    ! relaxation rate for
    ! (1) Rayleigh damping in land cells and
    ! (2) density-fluctuation relaxation in semi-implicit time stepping
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! for sponge:
    !-----------------------------------------------------------------
    if(spongeLayer) then
      if(unifiedSponge) then
        kr_sp = 0.0
        kr_sp_w = 0.0
        if(topography) then
          kr_sp_tfc = 0.0
          kr_sp_w_tfc = 0.0
        end if

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
              if(topography) then
                height = zTFC(i, j, k)
              else
                height = z(k)
              end if

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
        ! allocate( alpbls(0:ny+1),stat=allocstat)
        ! if(allocstat /= 0) stop "pinc.f90: could not allocate alpbls"
        ! maximum damping rate

        alpspg = spongeAlphaZ_fac / dt

        ! sponge-layer relaxation

        if(topography) then
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                if(zTFC(i, j, k) >= zSponge) then
                  kr_sp_tfc(i, j, k) = alpspg * sin(0.5 * pi * (zTFC(i, &
                      &j, k) - zSponge) / (lz(1) - zSponge)) ** 2.0
                  kr_sp_w_tfc(i, j, k) = kr_sp_tfc(i, j, k) / jac(i, j, k)
                end if
              end do
            end do
          end do

          kr_sp_tfc(:, :, 0) = kr_sp_tfc(:, :, 1)
          kr_sp_tfc(:, :, nz + 1) = kr_sp_tfc(:, :, nz)
          kr_sp_w_tfc(:, :, 0) = kr_sp_w_tfc(:, :, 1)
          kr_sp_w_tfc(:, :, nz + 1) = kr_sp_w_tfc(:, :, nz)
        else
          alpbls(:) = 0. !kr_sp(:,kSponge-1)

          j00 = js + nby - 1
          ymax = ly_dim(1) / lRef
          ymin = ly_dim(0) / lRef

          do k = kSponge, nz
            kr_sp(:, k) = alpbls(:) + (alpspg - alpbls(:)) * sin(0.5 * pi &
                &* (z(k) - zSponge) / (z(nz) - zSponge)) ** 2
            kr_sp_w(:, k) = kr_sp(:, k)
          end do

          kr_sp(:, nz + 1) = kr_sp(:, nz)
          kr_sp_w(ny + 1, :) = kr_sp_w(ny, :)
          kr_sp_w(0, :) = kr_sp_w(1, :)
        end if
      end if
    end if



    !---------------------------------------------------------------
    !          Runge-Kutta stages or semi-implicit time step
    !---------------------------------------------------------------

    if(timeScheme == "semiimplicit") then
      ! just for safety ...

      ! initialize zero volume force
      ! (to be filled by ray tracer and wind relaxation at the
      ! horizontal boundaries)

      ! set density fluctuation (synchronization step)

      ! Boussinesq: density fluctuations are only stored in var%rhop!
      select case(model)
      case("pseudo_incompressible")
        if(topography) then
          ! Stationary background in TFC.
          var%rhop(:, :, :) = var%rho(:, :, :)
        end if
      case default
      end select

      ! turbulence scheme:
      ! either prescribed damping time scale for smallest spatial scales
      ! or dynamic Smagorinsky scheme
      ! diffusion coefficient (normalized by squared grid length scale)
      ! stored in var (...,7)

      ! put initial state into var0 in order to save the advecting
      ! velocities

      var0 = var

      call setHalos(var0, "var")
      call setBoundary(var0, flux, "var")

      ! (1) explicit integration of convective and
      !     viscous-diffusive/turbulent fluxes over half a time step,
      !     with the advection velocity kept constant
      !     \psi^# = \psi^n + A^{dt/2} (\psi^n, v^n)

      if(master) print *, 'beginning a semi-implicit time step'
      if(master) print *, '(1) explicit integration lhs over dt/2'

      do RKstage = 1, nStages
        ! Reconstruction

        call setHalos(var, "var")
        call setBoundary(var, flux, "var")

        call reconstruction(var, "rho")
        call reconstruction(var, "rhop")
        call reconstruction(var, "uvw")

        call setHalos(var, "varTilde")
        call setBoundary(var, flux, "varTilde")

        ! Fluxes and Forces

        call massFlux(var0, var, flux, "lin")
        call momentumFlux(var0, var, flux, "lin")

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

        call momentumPredictor(var, flux, force, 0.5 * dt, dMom, RKstage, &
            &"lhs", "expl", 1.)

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

      if(topography) then
        ! uStar and vStar are needed for update of density fluctuations,
        ! therefore w is stored instead of rhop

        wOldTFC = var%w(:, :, :)

        ! update winds (uStar, vStar, wStar)

        call momentumPredictor(var, flux, force, 0.5 * dt, dMom, RKstage, &
            &"rhs", "impl", 1.)

        call setHalos(var, "var")
        call setBoundary(var, flux, "var")

        ! update density fluctuations (rhopStar)

        call massUpdate(var, flux, 0.5 * dt, dRhop, RKstage, "rhop", "rhs", &
            &"impl", 1.)

        call setHalos(var, "var")
        call setBoundary(var, flux, "var")
      end if


      ! corrector: rhopStar, uStar, vStar, wStar
      !            -> new rhop, u, v, w
      call Corrector(var, flux, dMom, 0.5 * dt, errFlagBicg, nIterBicg, &
          &RKstage, "impl", 1., 1.)

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
      !
      !     could also be replaced by an implicit time step

      if(master) print *, '(3) explicit integration rhs over dt/2'

      ! (3) uses updated pressure field and (5) adjusts pressure over half a
      ! time step!
      ! var = var0
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

      call momentumPredictor(var, flux, force, 0.5 * dt, dMom, RKstage, "rhs", &
          &"expl", 1.)

      ! corrector: uStar, vStar, wStar
      !            -> new u, v, w

      call setHalos(var, "var")
      call setBoundary(var, flux, "var")

      ! (4) explicit integration of convective and
      !     viscous-diffusive/turbulent fluxes over a full time step,
      !     with the advection velocity kept constant
      !     \psi^{\ast\ast} = \psi^\ast + A^dt (\psi^\ast, v^{n+1/2})

      if(master) print *, '(4) explicit integration lhs over dt'

      var0 = var1

      call setHalos(var0, "var")
      call setBoundary(var0, flux, "var")

      do RKstage = 1, nStages
        ! Reconstruction

        call setHalos(var, "var")
        call setBoundary(var, flux, "var")

        call reconstruction(var, "rho")
        call reconstruction(var, "rhop")
        call reconstruction(var, "uvw")

        call setHalos(var, "varTilde")
        call setBoundary(var, flux, "varTilde")

        ! Fluxes and Forces

        call massFlux(var0, var, flux, "lin")
        call momentumFlux(var0, var, flux, "lin")

        call setBoundary(var, flux, "flux")

        ! RK step for density and density fluctuations

        rhoOld = var%rho(:, :, :) ! rhoOld for momentum predictor

        call massUpdate(var, flux, dt, dRho, RKstage, "rho", "tot", "expl", 1.)

        call applyUnifiedSponge(var, stepFrac(RKstage) * dt, time, "rho")

        call massUpdate(var, flux, dt, dRhop, RKstage, "rhop", "lhs", "expl", &
            &1.)

        call applyUnifiedSponge(var, stepFrac(RKstage) * dt, time, "rhop")

        ! RK step for momentum
        call setHalos(var, "var")
        call setBoundary(var, flux, "var")

        call momentumPredictor(var, flux, force, dt, dMom, RKstage, "lhs", &
            &"expl", 1.)

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

      if(topography) then
        ! TFC FJ
        ! uStar and vStar are needed for update of density fluctuations,
        ! therefore w is stored instead of rhop

        wOldTFC = var%w(:, :, :)

        ! update winds (uStar, vStar, wStar)

        call momentumPredictor(var, flux, force, 0.5 * dt, dMom, RKstage, &
            &"rhs", "impl", 2.)

        call setHalos(var, "var")
        call setBoundary(var, flux, "var")

        ! update density fluctuations (rhopStar)

        call massUpdate(var, flux, 0.5 * dt, dRhop, RKstage, "rhop", "rhs", &
            &"impl", 2.)

        call setHalos(var, "var")
        call setBoundary(var, flux, "var")
      
      end if

      ! corrector: rhopStar, uStar, vStar, wStar
      !            -> new rhop, u, v, w

      call setHalos(var, "var")
      call setBoundary(var, flux, "var")

      ! (3) uses updated pressure field and (5) adjusts pressure over half a
      ! time step!
      ! call Corrector ( var, flux, dMom, 0.5*dt, errFlagBicg, nIterBicg, &
      !                & RKstage, "impl", 2.,2.) ! pressure update over dt
      call Corrector(var, flux, dMom, 0.5 * dt, errFlagBicg, nIterBicg, &
          &RKstage, "impl", 2., 1.) ! pressure update over dt/2

      if(errFlagBicg) then
        call write_netCDF(iOut, iTime, time, cpuTime, var)
        call close_netCDF

        if(master) then
          print *, 'output last state into record', iOut
        end if
        stop
      end if

      ! for safety
      call setHalos(var, "var")
      call setBoundary(var, flux, "var")

      nTotalBicg = nTotalBicg + nIterBicg

      if(master) print *, 'semi-implicit time step done'

    
    end if ! timeScheme

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
          if(timeScheme == 'semiimplicit') then
            nAverageBicg = real(nTotalBicg) / real(iTime) / 2.0
          end if

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

  10 if(master) then
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
  call terminate(var, var0, var1, varG, source, flux, flux0, force, dRho, &
      &dRhop, dMom, dTheta, dIce, dTracer, tracerforce, dPot)
  call terminate_atmosphere

  !CHANGES : cleaning master leads to MPI error ?
  if(master) then
    ! do nothing
    !call CleanUpBiCGSTab ! Clean Up BiCGSTAB arrays
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

end program pinc_prog
