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
  use wkb_module
  use xweno_module
  use atmosphere_module
  use boundary_module
  use flux_module
  use update_module
  use poisson_module
  use finish_module
  use ice_module
  use sizeof_module
  use bicgstab_tools_module
  use tracer_module
  use mpi
  use output_netCDF_module

  !----------------------------------------------------
  implicit none
  !---------------------------------------------------

  integer :: iTime, Ice_RKstage
  integer :: ice_time_steps
  real :: time, dt

  ! CPU Time
  integer :: rate, startTimeCount, timeCount
  real :: cpuTime

  ! MPI stuff
  logical :: error_flag
  real :: dt_local

  ! fields
  type(var_type) :: var, var0, var1, varG, source
  ! varG contains balance parts of current time step !FS October 2020

  real, dimension(:, :, :), allocatable :: dRho, dRhop ! RK-Update for rho
  real, dimension(:, :, :, :), allocatable :: dMom ! RK for rhoU,rhoV,rhoW
  real, dimension(:, :, :), allocatable :: dTheta ! RK-Update for theta
  real, dimension(:, :, :, :), allocatable :: dIce ! RK-Update for nAer,nIce,qIce,qv
  real, dimension(:, :, :), allocatable :: dTracer
  real, dimension(:, :, :), allocatable :: dPot !RK-Update for P
  real, dimension(:, :, :), allocatable :: dPhase

  real, dimension(:), allocatable :: dPStrat, drhoStrat !RK update for P
  real, dimension(:), allocatable :: w_0

  type(flux_type) :: flux, flux0

  !--------------------
  !   WKB variables
  !--------------------
  type(rayType), dimension(:, :, :, :), allocatable :: ray
  real, dimension(:, :, :, :), allocatable :: ray_var3D
  real, dimension(:, :, :), allocatable :: diffusioncoeff
  type(waveAmpType), dimension(:, :, :), allocatable :: waveAmplitudes

  ! topography via force field
  real, dimension(:, :, :, :), allocatable :: force ! volume forces
  type(tracerForceType), dimension(:, :, :), allocatable :: tracerforce ! tracer forcing

  ! output per timeStep
  logical :: output
  real :: nextOutputTime ! scaled time for next output

  ! general
  integer :: i, j, k, l
  integer :: ix, jy, kz
  integer :: iRay
  integer :: nrlc

  ! restart
  logical :: scale

  ! parameter study
  integer :: iParam

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
  real :: fc_shap, shap_dts

  integer :: iVar

  integer :: heating_switch !FS

  !UAB
  real, dimension(:), allocatable :: alpbls
  real, dimension(:), allocatable :: sum_local, sum_global
  integer :: allocstat
  !UAE
  integer :: j00
  real :: ymax, ymin, yloc, bla, wmax, wmax_loc

  ! TFC FJ
  integer :: i00
  real :: spongeAlphaZ, spongeAlphaY, spongeAlphaX

  real :: height

  !SD
  integer :: n_step_ice, ii
  real :: dtt_ice
  real :: uTime, qTime

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
      &dMom, dTheta, dPStrat, drhoStrat, w_0, dIce, dTracer, tracerforce, dPot)

  ! Read topography.
  if((topography .or. (rayTracer .and. case_wkb == 3)) .and. mountain_case &
      &== 0) then
    call read_topography_netCDF(iIn)
  end if

  call init_atmosphere ! set atmospheric background state

  call initialise(var, flux) ! set initial conditions

  ! put PstratTFC into var
  if(model == "compressible") then
    var%P(:, :, :) = pStratTFC(:, :, :)
  end if

  if(include_tracer) call setup_tracer(var)

  if(.not. topography .and. model /= "Boussinesq") then

    ! determine difference between reference-atmosphere density and
    ! horizontal-mean density

    allocate(sum_local(1:nz), stat = allocstat)
    if(allocstat /= 0) stop "pinc.f90: could not allocate sum_local"
    allocate(sum_global(1:nz), stat = allocstat)
    if(allocstat /= 0) stop "pinc.f90: could not allocate sum_global"

    sum_local = 0.
    sum_global = 0.

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          sum_local(k) = sum_local(k) + var%rho(i, j, k) + rhoStrat(k)
        end do
      end do
    end do
    call mpi_allreduce(sum_local(1), sum_global(1), nz, mpi_double_precision, &
        &mpi_sum, comm, ierror)
    sum_global = sum_global / (sizeX * sizeY)

    rhoStrat_s = 0.

    do k = 1, nz
      rhoStrat_s(k) = rhoStrat(k) - sum_global(k)
    end do

  end if

  if(rayTracer) var%GWH = 0.0 ! Heating due to GWs in the rotating atmosphere

  if(poissonSolverType == 'bicgstab') then
    call SetUpBiCGStab ! Set BiCGStab arrays
  else
    stop 'ERROR: only BiCGStab ready to be used'
  end if

  call init_xweno ! set ILES parameters
  call init_fluxes ! allocate tilde variables
  call init_update
  call init_timeScheme ! define Runge-Kutta parameters
  call init_poisson ! allocate dp

  ice_time_steps = 1

  if(zero_initial_state) then
    call reset_var_type(var)
  end if

  ! TFC FJ
  ! TFC tests.
  if(topography .and. testTFC) then
    do k = 1, nz
      kr_sp_tfc(:, :, k) = k
      kr_sp_w_tfc(:, :, k) = k
    end do
    pStrat_0 = pStrat
    call random_number(var%rho)
    call random_number(var%u)
    call random_number(var%v)
    call random_number(var%w)
    call random_number(var%pi)
    call random_number(var%rhop)
    call random_number(flux%rho)
    call random_number(flux%u)
    call random_number(flux%v)
    call random_number(flux%w)
    call random_number(flux%rhop)
    call random_number(force)
    if(timeScheme == "semiimplicit") then
      ! Linear operator test
      call setHalos(var, "var")
      call setBoundary(var, flux, "var")
      call linearOperatorTestTFC(var, 1.0, "impl", 1.0)
      ! Corrector step test
      call setHalos(var, "var")
      call setBoundary(var, flux, "var")
      call correctorStepTestTFC(var, dMom, "impl")
      ! Momentum predictor test
      call setHalos(var, "var")
      call setBoundary(var, flux, "var")
      call setBoundary(var, flux, "flux")
      call momentumPredictorTestTFC(var, flux, force, dMom, "impl")
      ! Mass update test
      call setHalos(var, "var")
      call setBoundary(var, flux, "var")
      call setBoundary(var, flux, "flux")
      call massUpdateTestTFC(var, flux, dRho, "impl")
    else
      ! Linear operator test
      call setHalos(var, "var")
      call setBoundary(var, flux, "var")
      call linearOperatorTestTFC(var, 1.0, "expl", 1.0)
      ! Corrector step test
      call setHalos(var, "var")
      call setBoundary(var, flux, "var")
      call correctorStepTestTFC(var, dMom, "expl")
      ! Momentum predictor test
      call setHalos(var, "var")
      call setBoundary(var, flux, "var")
      call setBoundary(var, flux, "flux")
      call momentumPredictorTestTFC(var, flux, force, dMom, "expl")
      ! Mass update test
      call setHalos(var, "var")
      call setBoundary(var, flux, "var")
      call setBoundary(var, flux, "flux")
      call massUpdateTestTFC(var, flux, dRho, "expl")
    end if
    ! Reconstruction test
    call setHalos(var, "var")
    call setBoundary(var, flux, "var")
    call reconstructionTestTFC(var)
    ! Momentum flux test
    call setHalos(var, "var")
    call setBoundary(var, flux, "var")
    call reconstruction(var, "uvw")
    call setHalos(var, "varTilde")
    call setBoundary(var, flux, "varTilde")
    call momentumFlux(var, var, flux, "nln", PStrat, PStratTilde)
    call setBoundary(var, flux, "flux")
    call setHalos(var, "var")
    call setBoundary(var, flux, "var")
    call momentumFluxTestTFC(var, flux, RKStage)
    ! Mass flux test
    call setHalos(var, "varTilde")
    call setBoundary(var, flux, "varTilde")
    call massFluxTestTFC(var, flux)
    if(master) then
      stop "TFC tests completed!"
    end if
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
    if(testCase == 'smoothVortex') then
      tolPoisson = 1.e-12 !FS
    else
      tolPoisson = 1.e-8
    end if

    if(model == "compressible") then
      call add_JP_to_u(var, "forward")
    end if

    call Corrector(var, flux, dMom, 1.0, errFlagBicg, nIterBicg, 1, "expl", &
        &1., 1.)

    if(model == "compressible") then
      call add_JP_to_u(var, "backward")
    end if

    if(errFlagBicg) stop

    tolPoisson = tolpoisson_s
  end if

  !-------------------------------------------------------------------
  ! store initial reference atmosphere
  !-------------------------------------------------------------------

  rhoStrat_0 = rhoStrat
  pStrat_0 = pStrat

  !---------------------------------------------
  !               Init ray tracer
  !---------------------------------------------

  if(rayTracer) then

    ! allocate and initialize ray fields
    call setup_wkb(ray, ray_var3D, var, diffusioncoeff, waveAmplitudes, dPhase)

    if(include_ice) then
      uTime = 0 !set initial time
      !Init ofield to zero except omega, phi (2,3)
      ofield(:, :, :, 1) = 0. !p_i(0)
      ofield(:, :, :, 4:6) = 0.
      call calc_ice(ray, var)
    end if
  end if

  !-------------------------------------------------
  !              Read initial data
  !-------------------------------------------------

  if(restart) then

    if(master) then
      print *, "reading restart files"
    end if

    if(rayTracer) then
      call read_netCDF(iIn, var, ray, time)
    else
      call read_netCDF(iIn, var, time = time)
    end if

    if(maxTime < time * tRef) stop "restart error: maxTime < current time"

    if(rayTracer) then
      do kz = 0, nz + 1
        do jy = 1, ny
          do ix = 1, nx
            nrlc = 0
            do iRay = 1, nray_max
              if(ray(iRay, ix, jy, kz)%dens == 0.0) cycle

              ray(iRay, ix, jy, kz)%area_xk = ray(iRay, ix, jy, kz)%dxray &
                  &* ray(iRay, ix, jy, kz)%dkray
              ray(iRay, ix, jy, kz)%area_yl = ray(iRay, ix, jy, kz)%dyray &
                  &* ray(iRay, ix, jy, kz)%dlray
              ray(iRay, ix, jy, kz)%area_zm = ray(iRay, ix, jy, kz)%dzray &
                  &* ray(iRay, ix, jy, kz)%dmray

              nrlc = nrlc + 1
              ray(nrlc, ix, jy, kz) = ray(iRay, ix, jy, kz)
            end do
            nRay(ix, jy, kz) = nrlc
          end do
        end do
      end do
    end if

    if(include_tracer) then
      if(topography) then
        do ix = 1, nx
          do jy = 1, ny
            do kz = 1, nz
              var%chi(ix, jy, kz) = var%chi(ix, jy, kz) * (var%rho(ix, jy, kz) &
                  &+ rhoStratTFC(ix, jy, kz))
            end do
          end do
        end do
      else
        do kz = 1, nz
          var%chi(:, :, kz) = var%chi(:, :, kz) * (var%rho(:, :, kz) &
              &+ rhoStrat(kz))
        end do
      end if
    end if

    call setHalos(var, "var")
    if(include_tracer) call setHalos(var, "tracer")
    call setBoundary(var, flux, "var")

    if(model == "compressible") then
      pStratTFC = var%P
      call bvsUpdate(bvsStratTFC, var)
    end if
  end if

  !------------------------------------------
  !              Initial output
  !------------------------------------------

  ! create netCDF file pincflow_data_out.nc
  call create_netCDF

  ! write initial state to netCDF file
  if(rayTracer) then
    call write_netCDF(iOut, iTime, time, cpuTime, var, waveAmplitudes &
        &= waveAmplitudes)
  else
    call write_netCDF(iOut, iTime, time, cpuTime, var)
  end if

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

    if(TestCase == "hotBubble_heat" .or. TestCase == "hotBubble_heatedLayer") &
        &then
      if(time * tRef > 250.) then
        heatingONK14 = .false.
      end if
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
      if(rayTracer) then
        if(include_tracer) then
          call write_netCDF(iOut, iTime, time, cpuTime, var, ray_var3D &
              &= ray_var3D, tracerforce = tracerforce, waveAmplitudes &
              &= waveAmplitudes, ray = ray)
        else
          call write_netCDF(iOut, iTime, time, cpuTime, var, ray_var3D &
              &= ray_var3D, waveAmplitudes = waveAmplitudes, ray = ray)
        end if
      else
        call write_netCDF(iOut, iTime, time, cpuTime, var)
      end if

      call close_netCDF

      call mpi_barrier(comm, ierror)
      call mpi_finalize(ierror)
      print *, "pinc.f90: time step too small. dt = ", dt * tRef
      stop
    end if

    if(timeSchemeType == "classical" .and. iTime == 1) dt = 0.5 * dt

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

    ! Update the topography.
    if(topography .and. topographyTime > 0.0) then
      call update_topography(time)
    end if

    ! Check for static instability.
    ! if(topography) then
    !   do k = 2, nz
    !     do j = 1, ny
    !       do i = 1, nx
    !         if(pStratTFC(i, j, k) / (var(i, j, k, 1) + rhoStratTFC(i, j, k)) &
    !             < pStratTFC(i, j, k - 1) / (var(i, j, k - 1, 1) &
    !             + rhoStratTFC(i, j, k - 1))) then
    !           print *, "Static instability at z =", heightTFC(i, j, k), "m"
    !         end if
    !       end do
    !     end do
    !   end do
    ! end if

    !-----------------------------------------------------------------
    ! relaxation rate for
    ! (1) Rayleigh damping in land cells and
    ! (2) density-fluctuation relaxation in semi-implicit time stepping
    !-----------------------------------------------------------------

    alprlx = 0. !0.5/dt

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
                height = heightTFC(i, j, k)
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
                if(heightTFC(i, j, k) >= zSponge) then
                  kr_sp_tfc(i, j, k) = alpspg * sin(0.5 * pi * (heightTFC(i, &
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

    !------------------------------------
    !    Call wave saturation scheme
    !------------------------------------

    if(rayTracer) then
      if(lsaturation) then
        call saturation_3D(ray, dt, diffusioncoeff)
      end if
    end if

    !---------------------------------------------------------------
    !          Runge-Kutta stages or semi-implicit time step
    !---------------------------------------------------------------

    if(timeScheme == "semiimplicit") then
      ! just for safety ...

      if(correctDivError) then
        print *, 'ERROR: correction divergence error not  implemented properly'
        stop
      end if

      if(updateTheta) then
        print *, 'ERROR: semiimplicit time stepping does not allow  &
            &updateTheta = .true.'
        stop
      end if

      ! initialize zero volume force
      ! (to be filled by ray tracer and wind relaxation at the
      ! horizontal boundaries)

      force = 0.0

      ! Lagrangian WKB model (position-wavenumber space method)

      if(rayTracer) then
        do RKstage = 1, nStages
          call transport_rayvol(var, ray, dt, RKstage, time)

          if(RKstage == nStages) then
            call split_rayvol(ray)
            call shift_rayvol(ray)
            call merge_rayvol(ray)
            call boundary_rayvol(ray)

            call calc_meanFlow_effect(ray, var, force, ray_var3D)

            if(include_tracer) then
              call calc_tracerforce(ray, var, ray_var3D, tracerforce, &
                  &waveAmplitudes, dt)
            end if

            if(include_ice) then
              call calc_ice(ray, var)
            end if
          end if
        end do
      end if

      ! wind relaxation at horizontal boundaries

      if((testCase == "mountainwave") .or. (raytracer .and. case_wkb == 3) &
          &.or. (topography .and. topographyTime > 0.0)) then
        call volumeForce(var, time, force)
      end if

      ! set density fluctuation (synchronization step)

      ! TFC FJ
      ! Boussinesq: density fluctuations are only stored in var(:, :, :, 6)!
      select case(model)
      case("pseudo_incompressible")
        if(topography) then
          ! TFC FJ
          ! Stationary background in TFC.
          var%rhop(:, :, :) = var%rho(:, :, :)
        else
          do kz = - 1, nz + 2
            var%rhop(:, :, kz) = var%rho(:, :, kz) + rhoStrat(kz) &
                &- rhoStrat_0(kz) * PStrat(kz) / PStrat_0(kz)
          end do
        end if
      case("compressible")
        do kz = - 1, nz + 2
          var%rhop(:, :, kz) = var%rho(:, :, kz) + rhoStratTFC(:, :, kz) &
              &- pStratTFC(:, :, kz) / thetaStratTFC(:, :, kz)
        end do
      case default
      end select

      ! turbulence scheme:
      ! either prescribed damping time scale for smallest spatial scales
      ! or dynamic Smagorinsky scheme
      ! diffusion coefficient (normalized by squared grid length scale)
      ! stored in var (...,7)

      if(TurbScheme) then
        if(DySmaScheme) then
          call CoefDySma_update(var)

          ! limit Smagorinsky coefficient so that the damping time
          ! scale for the 2dx-wave is shorter than a time step
          ! var(:, :, :, 7) = min(var(:, :, :, 7), 1.e0 / (dt * pi ** 2))
        else
          var%DSC(:, :, :) = tRef / turb_dts
        end if
      end if

      ! put initial state into var0 in order to save the advecting
      ! velocities

      var0 = var

      !for update of P and rhoStrat in (5)

      PStrat00 = PStrat
      heating_switch = 0

      PStratTilde00 = PStratTilde

      call setHalos(var0, "var")
      if(include_tracer) call setHalos(var0, "tracer")
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
        if(include_tracer) call setHalos(var, "tracer")
        call setBoundary(var, flux, "var")

        call reconstruction(var, "rho")
        call reconstruction(var, "rhop")
        call reconstruction(var, "uvw")
        if(include_tracer) then
          call reconstruction(var, "tracer")
        end if

        call setHalos(var, "varTilde")
        call setBoundary(var, flux, "varTilde")

        ! Fluxes and Forces

        call massFlux(var0, var, flux, "lin", PStrat00, PStratTilde00)
        call momentumFlux(var0, var, flux, "lin", PStrat00, PStratTilde00)
        if(include_tracer) then
          call tracerFlux(var0, var, flux, "lin", PStrat00, PStratTilde00)
        end if

        call setBoundary(var, flux, "flux")

        ! store initial flux
        if(RKstage == 1) flux0 = flux

        ! RK step for density and density fluctuations

        rhoOld = var%rho(:, :, :) ! rhoOld for momentum predictor

        call massUpdate(var, flux, 0.5 * dt, dRho, RKstage, "rho", "tot", &
            &"expl", 1.)

        call applyUnifiedSponge(var, stepFrac(RKstage) * 0.5 * dt, "rho")

        call massUpdate(var, flux, 0.5 * dt, dRhop, RKstage, "rhop", "lhs", &
            &"expl", 1.)

        call applyUnifiedSponge(var, stepFrac(RKstage) * 0.5 * dt, "rhop")

        if(model == "compressible") then
          call massUpdate(var, flux, 0.5 * dt, dPot, RKstage, "P", "tot", &
              &"expl", 1.)
          call applyUnifiedSponge(var, stepFrac(RKstage) * 0.5 * dt, "P")
        end if

        if(include_tracer) then
          call tracerUpdate(var, flux, tracerforce, 0.5 * dt, dTracer, RKstage)
        end if

        ! RK step for momentum

        call setHalos(var, "var")
        if(include_tracer) call setHalos(var, "tracer")
        call setBoundary(var, flux, "var")

        call momentumPredictor(var, flux, force, 0.5 * dt, dMom, RKstage, &
            &"lhs", "expl", 1.)

        call applyUnifiedSponge(var, stepFrac(RKstage) * 0.5 * dt, "uvw")

      end do

      if(model == "compressible") then
        pStratTFC(:, :, :) = var%P(:, :, :)

        if(TurbScheme .or. rayTracer) then
          call piUpdate(var0, pinew, 0.5 * dt, "expl_heating", flux0)
          var%pi(:, :, :) = pinew(:, :, :)
        end if
        call applyUnifiedSponge(var, 0.5 * dt, "pi")
      end if

      ! (2) implicit integration of the linear right-hand sides of the
      !     equations for density fluctuations and momentum over half a
      !     time step, under consideration of the divergence constraint
      !     \psi^{n+1/2} = \psi^# + dt/2 Q(\psi^{n+1/2})

      if(master) print *, '(2) implicit integration rhs over dt/2'

      ! SK: Save JPu in var instead of u for updates
      if(model == "compressible") then
        call add_JP_to_u(var, "forward")
        call bvsUpdate(bvsStratTFC, var) ! Update N^2
      end if

      call setHalos(var, "var")
      if(include_tracer) call setHalos(var, "tracer")
      call setBoundary(var, flux, "var")

      ! use initial flux for update of reference atmosphere and w0
      flux = flux0

      if(heatingONK14 .or. TurbScheme .or. rayTracer) then
        if(model == 'Boussinesq') then
          print *, "main:ONeill+Klein2014 heating only for  &
              &pseudo-incompressible dyn."
          stop
        end if

        RKstage = 1
        dPStrat = 0.
        drhoStrat = 0.
        call BGstate_update(var, flux, 0.5 * dt, RKstage, dPStrat, drhoStrat, &
            &"impl", heating_switch)
      end if

      if(topography) then
        ! uStar and vStar are needed for update of density fluctuations,
        ! therefore w is stored instead of rhop

        wOldTFC = var%w(:, :, :)

        ! update winds (uStar, vStar, wStar)

        call momentumPredictor(var, flux, force, 0.5 * dt, dMom, RKstage, &
            &"rhs", "impl", 1.)

        call setHalos(var, "var")
        if(include_tracer) call setHalos(var, "tracer")
        call setBoundary(var, flux, "var")

        ! update density fluctuations (rhopStar)

        call massUpdate(var, flux, 0.5 * dt, dRhop, RKstage, "rhop", "rhs", &
            &"impl", 1.)

        call setHalos(var, "var")
        if(include_tracer) call setHalos(var, "tracer")
        call setBoundary(var, flux, "var")
      else
        rhopOld = var%rhop(:, :, :) ! rhopOld for momentum predictor

        ! update density fluctuations (rhopStar)

        call massUpdate(var, flux, 0.5 * dt, dRhop, RKstage, "rhop", "rhs", &
            &"impl", 1.)

        ! update winds (uStar, vStar, wStar)

        call momentumPredictor(var, flux, force, 0.5 * dt, dMom, RKstage, &
            &"rhs", "impl", 1.)

        call setHalos(var, "var")
        if(include_tracer) call setHalos(var, "tracer")
        call setBoundary(var, flux, "var")
      end if

      ! Shapiro filter
      if(shap_dts_fac > 0.) then
        ! smoothing of the fields in order to limit grid-point noise

        shap_dts = shap_dts_fac / tRef !FS shap_dts_fac * dt

        fc_shap = min(1.0, 0.5 * dt / shap_dts)

        call smooth_hor_shapiro(fc_shap, n_shap, flux, var, 0.5 * dt)
      end if

      ! corrector: rhopStar, uStar, vStar, wStar
      !            -> new rhop, u, v, w
      call Corrector(var, flux, dMom, 0.5 * dt, errFlagBicg, nIterBicg, &
          &RKstage, "impl", 1., 1.)

      if(errFlagBicg) then
        if(rayTracer) then
          if(include_tracer) then
            call write_netCDF(iOut, iTime, time, cpuTime, var, ray_var3D &
                &= ray_var3D, tracerforce = tracerforce, waveAmplitudes &
                &= waveAmplitudes, ray = ray)
          else
            call write_netCDF(iOut, iTime, time, cpuTime, var, ray_var3D &
                &= ray_var3D, waveAmplitudes = waveAmplitudes, ray = ray)
          end if
        else
          call write_netCDF(iOut, iTime, time, cpuTime, var)
        end if

        call close_netCDF

        if(master) then
          print *, 'output last state into record', iOut
        end if
        stop
      end if

      nTotalBicg = nTotalBicg + nIterBicg

      ! SK: remove JP from JPu
      if(model == "compressible") then
        call add_JP_to_u(var, "backward")
      end if

      call setHalos(var, "var")
      if(include_tracer) call setHalos(var, "tracer")
      call setBoundary(var, flux, "var")

      ! put new state into var1 in order to save the advection velocities

      var1 = var

      PStrat01 = PStrat
      PStratTilde01 = PStratTilde

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
      if(turbScheme) var%DSC(:, :, :) = var0%DSC(:, :, :)

      if(model == "compressible") then
        var%pi(:, :, :) = var0%pi(:, :, :) ! reset also pressure
        var%P(:, :, :) = var0%P(:, :, :)
        pStratTFC(:, :, :) = var0%P(:, :, :)
        call bvsUpdate(bvsStratTFC, var) ! Update N^2
        call add_JP_to_u(var, "forward") ! integrating with JPu not u on the right-hand side, so save JPu in var
      end if

      if(include_tracer) then
        var%chi(:, :, :) = var0%chi(:, :, :)
      end if

      PStrat = PStrat00
      PStratTilde = PStratTilde00

      call setHalos(var, "var")
      if(include_tracer) call setHalos(var, "tracer")
      call setBoundary(var, flux, "var")

      rhopOld = var%rhop(:, :, :) ! rhopOld for momentum predictor

      ! update density fluctuations (rhopStar)
      call massUpdate(var, flux, 0.5 * dt, dRhop, RKstage, "rhop", "rhs", &
          &"expl", 1.)

      ! update winds (uStar, vStar, wStar)

      if(model == "compressible") then
        ! Update pi' explicitly
        call piUpdate(var, pinew, 0.5 * dt, "expl", flux) ! pinew is updated field
      end if

      call momentumPredictor(var, flux, force, 0.5 * dt, dMom, RKstage, "rhs", &
          &"expl", 1.)

      ! SK: Save u in var
      if(model == "compressible") then
        call add_JP_to_u(var, "backward")
        var%pi(:, :, :) = pinew(:, :, :) ! update pi' with pinew from the piUpdate
      end if
      ! Shapiro filter

      if(shap_dts_fac > 0.) then
        ! smoothing of the fields in order to limit grid-point noise

        shap_dts = shap_dts_fac / tRef !FS shap_dts_fac * dt

        fc_shap = min(1.0, 0.5 * dt / shap_dts)

        call smooth_hor_shapiro(fc_shap, n_shap, flux, var, 0.5 * dt)
      end if

      ! corrector: uStar, vStar, wStar
      !            -> new u, v, w

      call setHalos(var, "var")
      if(include_tracer) call setHalos(var, "tracer")
      call setBoundary(var, flux, "var")

      ! (4) explicit integration of convective and
      !     viscous-diffusive/turbulent fluxes over a full time step,
      !     with the advection velocity kept constant
      !     \psi^{\ast\ast} = \psi^\ast + A^dt (\psi^\ast, v^{n+1/2})

      if(master) print *, '(4) explicit integration lhs over dt'

      var0 = var1

      if(model == "compressible") then
        var1 = var
        pStratTFC(:, :, :) = var0%P(:, :, :) ! set pstrattfc from (2)
      end if

      !FS
      PStrat = PStrat01
      PStratTilde = PStratTilde01

      call setHalos(var0, "var")
      if(include_tracer) call setHalos(var0, "tracer")
      call setBoundary(var0, flux, "var")

      do RKstage = 1, nStages
        ! Reconstruction

        call setHalos(var, "var")
        if(include_tracer) call setHalos(var, "tracer")
        call setBoundary(var, flux, "var")

        call reconstruction(var, "rho")
        call reconstruction(var, "rhop")
        call reconstruction(var, "uvw")
        if(include_tracer) then
          call reconstruction(var, "tracer")
        end if

        call setHalos(var, "varTilde")
        call setBoundary(var, flux, "varTilde")

        ! Fluxes and Forces

        call massFlux(var0, var, flux, "lin", PStrat01, PStratTilde01)
        call momentumFlux(var0, var, flux, "lin", PStrat01, PStratTilde01)

        if(include_tracer) then
          call tracerFlux(var0, var, flux, "lin", PStrat01, PStratTilde01)
        end if

        call setBoundary(var, flux, "flux")

        ! RK step for density and density fluctuations

        rhoOld = var%rho(:, :, :) ! rhoOld for momentum predictor

        call massUpdate(var, flux, dt, dRho, RKstage, "rho", "tot", "expl", 1.)

        call applyUnifiedSponge(var, stepFrac(RKstage) * dt, "rho")

        call massUpdate(var, flux, dt, dRhop, RKstage, "rhop", "lhs", "expl", &
            &1.)

        call applyUnifiedSponge(var, stepFrac(RKstage) * dt, "rhop")

        if(model == "compressible") then
          call massUpdate(var, flux, dt, dPot, RKstage, "P", "tot", "expl", 1.)
          call applyUnifiedSponge(var, stepFrac(RKstage) * dt, "P")
        end if

        if(include_tracer) then
          call tracerUpdate(var, flux, tracerforce, dt, dTracer, RKstage)
        end if

        ! RK step for momentum
        call setHalos(var, "var")
        if(include_tracer) call setHalos(var, "tracer")
        call setBoundary(var, flux, "var")

        call momentumPredictor(var, flux, force, dt, dMom, RKstage, "lhs", &
            &"expl", 1.)

        call applyUnifiedSponge(var, stepFrac(RKstage) * dt, "uvw")

        !SD
        !if(include_ice) call integrate_ice_advection(var, var0, flux, "lin", &
        !    source, dt, dIce, RKstage, PStrat01, PStratTilde01)

      end do

      if(model == "compressible") then
        pStratTFC = var%P(:, :, :)

        if(TurbScheme .or. rayTracer) then
          call piUpdate(var1, pinew, dt, "expl_heating", flux0) ! use explicitly updated pi'
          var%pi(:, :, :) = pinew(:, :, :)
        end if
        call applyUnifiedSponge(var, dt, "pi")
      end if

      ! (5) implicit integration of the linear right-hand sides of the
      !     equations for density fluctuations and momentum over half a
      !     time step, under consideration of the divergence constraint
      !     \psi^{n+1} = \psi^{\ast\ast} + dt/2 Q(\psi^{n+1})

      if(master) print *, '(5) implicit integration rhs over dt/2'

      ! SK: Save JPu in var instead of u for updates
      if(model == "compressible") then
        call add_JP_to_u(var, "forward")
        call bvsUpdate(bvsStratTFC, var) ! Update N^2
      end if

      call setHalos(var, "var")
      if(include_tracer) call setHalos(var, "tracer")
      call setBoundary(var, flux, "var")

      ! use initial flux for update of reference atmosphere and w0
      flux = flux0

      if(heatingONK14 .or. TurbScheme .or. rayTracer) then
        if(model == 'Boussinesq') then
          print *, "main:ONeill+Klein2014 heating only for  &
              &pseudo-incompressible dyn."
          stop
        end if

        RKstage = 1
        dPStrat = 0.
        drhoStrat = 0.
        heating_switch = 1
        call BGstate_update(var, flux, dt, RKstage, dPStrat, drhoStrat, &
            &"impl", heating_switch) !FS 0.5*dt -> dt
      end if

      if(topography) then
        ! TFC FJ
        ! uStar and vStar are needed for update of density fluctuations,
        ! therefore w is stored instead of rhop

        wOldTFC = var%w(:, :, :)

        ! update winds (uStar, vStar, wStar)

        call momentumPredictor(var, flux, force, 0.5 * dt, dMom, RKstage, &
            &"rhs", "impl", 2.)

        call setHalos(var, "var")
        if(include_tracer) call setHalos(var, "tracer")
        call setBoundary(var, flux, "var")

        ! update density fluctuations (rhopStar)

        call massUpdate(var, flux, 0.5 * dt, dRhop, RKstage, "rhop", "rhs", &
            &"impl", 2.)

        call setHalos(var, "var")
        if(include_tracer) call setHalos(var, "tracer")
        call setBoundary(var, flux, "var")
      else
        rhopOld = var%rhop(:, :, :) ! rhopOld for momentum predictor

        ! update density fluctuations (rhopStar)

        call massUpdate(var, flux, 0.5 * dt, dRhop, RKstage, "rhop", "rhs", &
            &"impl", 2.)

        call setHalos(var, "var")
        if(include_tracer) call setHalos(var, "tracer")
        call setBoundary(var, flux, "var")

        call momentumPredictor(var, flux, force, 0.5 * dt, dMom, RKstage, &
            &"rhs", "impl", 2.)
      end if

      ! Shapiro filter

      if(shap_dts_fac > 0.) then
        ! smoothing of the fields in order to limit grid-point noise

        shap_dts = shap_dts_fac / tRef !FS shap_dts_fac * dt

        fc_shap = min(1.0, 0.5 * dt / shap_dts)

        call smooth_hor_shapiro(fc_shap, n_shap, flux, var, 0.5 * dt)
      end if

      ! corrector: rhopStar, uStar, vStar, wStar
      !            -> new rhop, u, v, w

      call setHalos(var, "var")
      if(include_tracer) call setHalos(var, "tracer")
      call setBoundary(var, flux, "var")

      ! (3) uses updated pressure field and (5) adjusts pressure over half a
      ! time step!
      ! call Corrector ( var, flux, dMom, 0.5*dt, errFlagBicg, nIterBicg, &
      !                & RKstage, "impl", 2.,2.) ! pressure update over dt
      call Corrector(var, flux, dMom, 0.5 * dt, errFlagBicg, nIterBicg, &
          &RKstage, "impl", 2., 1.) ! pressure update over dt/2

      if(model == "compressible") then
        call add_JP_to_u(var, "backward")
      end if

      if(errFlagBicg) then
        if(rayTracer) then
          if(include_tracer) then
            call write_netCDF(iOut, iTime, time, cpuTime, var, ray_var3D &
                &= ray_var3D, tracerforce = tracerforce, waveAmplitudes &
                &= waveAmplitudes, ray = ray)
          else
            call write_netCDF(iOut, iTime, time, cpuTime, var, ray_var3D &
                &= ray_var3D, waveAmplitudes = waveAmplitudes, ray = ray)
          end if
        else
          call write_netCDF(iOut, iTime, time, cpuTime, var)
        end if

        call close_netCDF

        if(master) then
          print *, 'output last state into record', iOut
        end if
        stop
      end if

      ! for safety
      call setHalos(var, "var")
      if(include_tracer) call setHalos(var, "tracer")
      call setBoundary(var, flux, "var")

      nTotalBicg = nTotalBicg + nIterBicg

      if(master) print *, 'semi-implicit time step done'

    else ! (timeScheme /= "semiimplicit") explicit time stepping

      Runge_Kutta_Loop: do RKstage = 1, nStages

        ! turbulence scheme:
        ! either prescribed damping time scale for smallest spatial
        ! scales
        ! or dynamic Smagorinsky scheme
        ! diffusion coefficient (normali. by squared grid length scale)
        ! stored in var (...,7)

        if(TurbScheme) then
          if(DySmaScheme) then
            call CoefDySma_update(var)
          else
            var%DSC(:, :, :) = tRef / turb_dts
          end if
        end if

        if(timeSchemeType == "classical") then
          if(RKstage == 1) then
            var0 = var
          end if
        end if

        ! initialize density fluctuations for the integration
        if(auxil_equ .or. heatingONK14 .or. TurbScheme .or. rayTracer) then
          if(model /= "pseudo_incompressible" .and. model /= "Boussinesq") then
            stop "Auxiliary equation only ready for pseudo-incompressible and &
                &Boussinesq"
          end if

          alprlx = 0.

          if(model == "pseudo_incompressible") then
            if(topography) then
              ! Stationary background in TFC.
              var%rhop(:, :, :) = var%rho(:, :, :)
            else
              do kz = - 1, nz + 2
                var%rhop(:, :, kz) = var%rho(:, :, kz) + rhoStrat(kz) &
                    &- rhoStrat_0(kz) * PStrat(kz) / PStrat_0(kz)
              end do
            end if
          end if

          var0 = var

        end if

        ! initialize zero volume force

        force = 0.0

        ! Lag Ray tracer (position-wavenumber space method)

        if(rayTracer) then
          call transport_rayvol(var, ray, dt, RKstage, time)

          if(RKstage == nStages) then
            call split_rayvol(ray)
            call shift_rayvol(ray)
            call merge_rayvol(ray)
            call boundary_rayvol(ray)

            call calc_meanFlow_effect(ray, var, force, ray_var3D)

            if(include_tracer) then
              call calc_tracerforce(ray, var, ray_var3D, tracerforce, &
                  &waveAmplitudes, dt)
            end if

            if(include_ice) then
              call calc_ice(ray, var)
            end if
          end if
        end if

        ! Reconstruction

        call setHalos(var, "var")
        if(include_tracer) call setHalos(var, "tracer")
        call setBoundary(var, flux, "var")

        if(updateMass .or. (testcase == "nIce_w_test")) call &
            &reconstruction(var, "rho")
        if(updateMass) call reconstruction(var, "rho")
        if(updateMass .and. auxil_equ) then
          call reconstruction(var, "rhop")
        end if

        if(updateTheta) call reconstruction(var, "theta")
        if(predictMomentum .or. (testcase == "nIce_w_test")) call &
            &reconstruction(var, "uvw")
        if((include_tracer) .and. (updateTracer)) call reconstruction(var, &
            &"tracer")

        call setHalos(var, "varTilde")
        call setBoundary(var, flux, "varTilde")

        ! Fluxes and Forces

        if(updateMass) then
          call massFlux(var, var, flux, "nln", PStrat, PStratTilde)
          if(correctDivError) then
            print *, 'ERROR: correction divergence error not  implemented &
                &properly'
            stop
          end if
        end if

        if(updateTracer) then
          call tracerFlux(var, var, flux, "nln", PStrat, PStratTilde)
        end if

        if(updateTheta) then
          call thetaFlux(var, flux)
          call thetaSource(var, source)
        end if

        if(predictMomentum) then
          call momentumFlux(var, var, flux, "nln", PStrat, PStratTilde)
          call volumeForce(var, time, force)
        end if

        call setBoundary(var, flux, "flux")

        ! implementation of heating ONeill and Klein 2014

        if(heatingONK14 .or. TurbScheme .or. rayTracer) then

          if(model == 'Boussinesq') then
            print *, "main:ONeill+Klein2014 heating only for  &
                &pseudo-incompressible dyn."
            stop
          end if

          if(RKstage == 1) dPStrat = 0. ! init q
          if(RKstage == 1) drhoStrat = 0.

          call BGstate_update(var, flux, dt, RKstage, dPStrat, drhoStrat, &
              &"expl", heating_switch)

        end if

        if(updateMass) then
          ! rho_new

          if(RKstage == 1) dRho = 0. ! init q

          rhoOld = var%rho(:, :, :) ! rhoOld for momentum predictor
          call massUpdate(var, flux, dt, dRho, RKstage, "rho", "tot", "expl", &
              &1.)
          if(testCase /= 'baroclinic_LC') then
            call applyUnifiedSponge(var, stepFrac(RKstage) * dt, "rho")
          end if

          if(auxil_equ) then
            if(RKstage == 1) dRhop = 0. ! init q

            rhopOld = var%rhop(:, :, :) ! rhopOld for momentum predictor
            call massUpdate(var, flux, dt, dRhop, RKstage, "rhop", "tot", &
                &"expl", 1.)
            if(testCase /= 'baroclinic_LC') then
              call applyUnifiedSponge(var, stepFrac(RKstage) * dt, "rhop")
            end if
          end if

        else
          if(iTime == 1 .and. RKstage == 1 .and. master) print *, "main: &
              &MassUpdate off!"
        end if

        if(updateTracer) then
          if(RKstage == 1) dTracer = 0.0
          call tracerUpdate(var, flux, tracerforce, dt, dTracer, RKstage)
        end if

        if(updateTheta) then
          ! theta_new

          call setHalos(var, "var")
          call setBoundary(var, flux, "var")

          if(RKstage == 1) dTheta = 0. ! init q

          call thetaUpdate(var, var0, flux, source, dt, dTheta, RKstage)
        else
          if(iTime == 1 .and. RKstage == 1 .and. master) print *, "main: &
              &ThetaUpdate off!"
        end if

        if(predictMomentum) then
          ! predictor: uStar

          call setHalos(var, "var")
          call setBoundary(var, flux, "var")

          if(RKstage == 1) dMom = 0. ! init q

          call momentumPredictor(var, flux, force, dt, dMom, RKstage, "tot", &
              &"expl", 1.)
          if(testCase /= 'baroclinic_LC') then
            call applyUnifiedSponge(var, stepFrac(RKstage) * dt, "uvw")
          end if
        else
          if(iTime == 1 .and. RKstage == 1 .and. master) print *, "main: &
              &MomentumUpdate off!"
        end if

        if(shap_dts_fac > 0.) then
          ! smoothing of the fields in order to limit grid-point noise

          select case(timeSchemeType)
          case("lowStorage")
            dt_Poisson = betaRK(RKstage) * dt
          case("classical")
            dt_Poisson = rk(3, RKstage) * dt
          case default
            stop "thetaUpdate: unknown case timeSchemeType"
          end select

          shap_dts = shap_dts_fac / tRef

          fc_shap = min(1.0, dt_Poisson / shap_dts)

          call smooth_hor_shapiro(fc_shap, n_shap, flux, var, dt_Poisson)
        end if

        if(correctMomentum) then
          ! corrector: dp, du -> u_new, p_new

          call setHalos(var, "var")
          call setBoundary(var, flux, "var")

          select case(timeSchemeType)
          case("lowStorage")
            dt_Poisson = betaRK(RKstage) * dt
          case("classical")
            dt_Poisson = rk(3, RKstage) * dt
          case default
            stop "thetaUpdate: unknown case timeSchemeType"
          end select

          call Corrector(var, flux, dMom, dt_Poisson, errFlagBicg, nIterBicg, &
              &RKstage, "expl", 1., 1.)

          if(errFlagBicg) then
            if(rayTracer) then
              if(include_tracer) then
                call write_netCDF(iOut, iTime, time, cpuTime, var, ray_var3D &
                    &= ray_var3D, tracerforce = tracerforce, waveAmplitudes &
                    &= waveAmplitudes, ray = ray)
              else
                call write_netCDF(iOut, iTime, time, cpuTime, var, ray_var3D &
                    &= ray_var3D, waveAmplitudes = waveAmplitudes, ray = ray)
              end if
            else
              call write_netCDF(iOut, iTime, time, cpuTime, var)
            end if

            call close_netCDF

            if(master) then
              print *, 'output last state into record', iOut
            end if
            stop
          end if

          nTotalBicg = nTotalBicg + nIterBicg

          ! error handling
          if(errFlagBiCG) then
            print *, " Momentum Corrector error"
            write(*, fmt = "(a25,i2.2)") "maxIterPoisson! iPara=", iParam
            nAverageBicg = real(nTotalBicg) / real(iTime) / 3.0

            call system_clock(count = timeCount)
            cpuTime = (timeCount - startTimeCount) / real(rate)

            if(rayTracer) then
              if(include_tracer) then
                call write_netCDF(iOut, iTime, time, cpuTime, var, ray_var3D &
                    &= ray_var3D, tracerforce = tracerforce, waveAmplitudes &
                    &= waveAmplitudes, ray = ray)
              else
                call write_netCDF(iOut, iTime, time, cpuTime, var, ray_var3D &
                    &= ray_var3D, waveAmplitudes = waveAmplitudes, ray = ray)
              end if
            else
              call write_netCDF(iOut, iTime, time, cpuTime, var)
            end if

            call close_netCDF

            go to 10 ! dealloc fields
          end if
        else
          if(iTime == 1 .and. RKstage == 1) print *, "main: MomentumCorrector &
              &off!"
        end if

      end do Runge_Kutta_Loop
    end if ! timeScheme

    ! integrate ice physics/advection with fixed dry dynamics
    if(include_ice) then

      n_step_ice = ceiling(dt * tRef / dt_ice)
      dtt_ice = dt / n_step_ice

      do ii = 1, n_step_ice
        do RKstage = 1, nStages
          call integrate_ice(var, var0, flux, "nln", source, dtt_ice, dIce, &
              &RKstage, PStrat, PStratTilde, 'BOT', uTime, qTime)
        end do ! RKstage
      end do !ii
      if(master) then
        print *, 'maxval', maxval(var%ICE(:, :, :, inN)), maxval(var%ICE(:, :, &
            &:, inQ)), maxval(var%ICE(:, :, :, inQv))
      end if
    end if !include_ice

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

        if(rayTracer) then
          if(include_tracer) then
            call write_netCDF(iOut, iTime, time, cpuTime, var, ray_var3D &
                &= ray_var3D, tracerforce = tracerforce, waveAmplitudes &
                &= waveAmplitudes, ray = ray)
          else
            call write_netCDF(iOut, iTime, time, cpuTime, var, ray_var3D &
                &= ray_var3D, waveAmplitudes = waveAmplitudes, ray = ray)
          end if
        else
          call write_netCDF(iOut, iTime, time, cpuTime, var)
        end if

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

        if(rayTracer) then
          if(include_tracer) then
            call write_netCDF(iOut, iTime, time, cpuTime, var, ray_var3D &
                &= ray_var3D, tracerforce = tracerforce, waveAmplitudes &
                &= waveAmplitudes, ray = ray)
          else
            call write_netCDF(iOut, iTime, time, cpuTime, var, ray_var3D &
                &= ray_var3D, waveAmplitudes = waveAmplitudes, ray = ray)
          end if
        else
          call write_netCDF(iOut, iTime, time, cpuTime, var)
        end if

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
          else
            nAverageBicg = real(nTotalBicg) / real(iTime) / 3.0
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

  if(poissonSolverType == 'bicgstab') then

    !CHANGES : cleaning master leads to MPI error ?
    if(master) then
      ! do nothing
      !call CleanUpBiCGSTab ! Clean Up BiCGSTAB arrays
    else
      call CleanUpBiCGSTab ! Clean Up BiCGSTAB arrays
    end if

  else
    stop 'ERROR: BICGSTAB expected as Poisson solver'
  end if

  666 if(master) then
    write(*, "(a)") ""
    write(*, "(a)") repeat("-", 80)
    write(*, "(a)") repeat(" ", 32) // "PincFlow finished" // repeat(" ", 33)
    write(*, "(a)") repeat("-", 80)
  end if

  call mpi_finalize(ierror)

end program pinc_prog
