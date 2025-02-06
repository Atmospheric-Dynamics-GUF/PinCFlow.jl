module output_netCDF_module

  use type_module
  use atmosphere_module
  use netcdf
  use mpi

  implicit none

  public

  ! netCDF parameters
  character(len = *), parameter :: filenameOut = 'pincflow_data_out.nc'
  character(len = *), parameter :: filenameIn = 'pincflow_data_in.nc'
  integer :: mode_flagNC ! mode_flag for creating the netCDF file
  integer :: ncid ! id of the netCDF file
  integer :: atmvarid ! group ids

  ! start and count values for nf90_put_var(...)
  integer :: startxNC, startyNC, countxNC, countyNC
  integer :: starttNC, counttNC
  integer, dimension(4) :: startNC, countNC

  ! ids of the dimensions
  integer :: dimidx, dimidy, dimidz, dimidt ! x, y, z, t

  ! variable ids of the axis
  integer :: ncidx, ncidy, ncidz, ncidt ! x, y, z, t

  ! variable ids of the topography
  integer :: ncidhm ! large-scale topography

  ! variable ids in the group atmvar
  integer :: atmvarid_u, atmvarid_v, atmvarid_w ! u, v, w wind
  integer :: atmvarid_wTFC ! terrain-following w
  integer :: atmvarid_us, atmvarid_vs, atmvarid_ws ! u, v, w on staggered
  integer :: atmvarid_wsTFC ! terrain-following ws
  integer :: atmvarid_b ! buoyancy
  integer :: atmvarid_rhobar, atmvarid_p, atmvarid_n2, atmvarid_thetabar
  ! rhoStratTFC, pStratTFC, bvsStratTFC, thetaStratTFC
  integer :: atmvarid_rhop ! density fluctuations rho'
  integer :: atmvarid_thetap ! potential temperature fluctuations theta'
  integer :: atmvarid_pip ! Exner-pressure fluctuations pi'

  contains

  subroutine create_netCDF

    implicit none

    ! date, time, and zone when file was created
    character(8) :: dateNC
    character(10) :: timeNC
    character(5) :: zoneNC

    call date_and_time(date = dateNC, time = timeNC, zone = zoneNC)

    mode_flagNC = IOR(nf90_netcdf4, nf90_mpiio)
    mode_flagNC = IOR(mode_flagNC, nf90_noclobber)

    ! create netCDF file
    call handle_err(nf90_create(filenameOut, mode_flagNC, ncid, comm = comm, &
        &info = mpi_info_null))

    ! set global attributes (file information)
    call handle_err(nf90_put_att(ncid, nf90_global, 'Title', 'PinCFlow data'))
    call handle_err(nf90_put_att(ncid, nf90_global, 'Run', runName))
    call handle_err(nf90_put_att(ncid, nf90_global, 'Institution', 'Institute &
        &for Atmospheric and Environmental Sciences, Goethe University &
        &Frankfurt, Germany'))
    call handle_err(nf90_put_att(ncid, nf90_global, 'Date', dateNC(1:4) // '-' &
        &// dateNC(5:6) // '-' // dateNC(7:8)))
    call handle_err(nf90_put_att(ncid, nf90_global, 'Time', timeNC(1:2) // ':' &
        &// timeNC(3:4) // ':' // timeNC(5:6)))
    call handle_err(nf90_put_att(ncid, nf90_global, 'Zone', zoneNC(1:3)))

    ! define the x-, y-, z-, and t-dimensions
    call handle_err(nf90_def_dim(ncid, 'x', sizeX, dimidx))
    call handle_err(nf90_def_dim(ncid, 'y', sizeY, dimidy))
    call handle_err(nf90_def_dim(ncid, 'z', sizeZ, dimidz))
    call handle_err(nf90_def_dim(ncid, 't', nf90_unlimited, dimidt))

    !--------------------------
    ! define the grid variables
    !--------------------------

    ! x axis
    call handle_err(nf90_def_var(ncid, 'x', nf90_float, [dimidx], ncidx))
    call handle_err(nf90_var_par_access(ncid, ncidx, nf90_collective))
    call handle_err(nf90_put_att(ncid, ncidx, 'units', 'm'))

    ! y axis
    call handle_err(nf90_def_var(ncid, 'y', nf90_float, [dimidy], ncidy))
    call handle_err(nf90_var_par_access(ncid, ncidy, nf90_collective))
    call handle_err(nf90_put_att(ncid, ncidy, 'units', 'm'))

    ! z axis
    call handle_err(nf90_def_var(ncid, 'z', nf90_float, [dimidx, dimidy, &
        &dimidz, dimidt], ncidz))
    call handle_err(nf90_var_par_access(ncid, ncidz, nf90_collective))
    call handle_err(nf90_put_att(ncid, ncidz, 'units', 'm'))

    ! t axis
    call handle_err(nf90_def_var(ncid, 't', nf90_float, [dimidt], ncidt))
    call handle_err(nf90_var_par_access(ncid, ncidt, nf90_collective))
    call handle_err(nf90_put_att(ncid, ncidt, 'units', 's'))

    ! large-scale topography
    call handle_err(nf90_def_var(ncid, 'hm', nf90_float, [dimidx, dimidy, &
        &dimidt], ncidhm))
    call handle_err(nf90_var_par_access(ncid, ncidhm, nf90_collective))
    call handle_err(nf90_put_att(ncid, ncidhm, 'units', 'm'))

    !---------------------------------------------------
    ! create the groups (atmvar)
    !---------------------------------------------------

    call handle_err(nf90_def_grp(ncid, 'atmvar', atmvarid))

    !-------------------------------------
    ! define the variables in atmvar group
    !-------------------------------------

    ! potential-temperature background
    call handle_err(nf90_def_var(atmvarid, 'thetabar', nf90_float, [dimidx, &
        &dimidy, dimidz, dimidt], atmvarid_thetabar))
    call handle_err(nf90_var_par_access(atmvarid, atmvarid_thetabar, &
        &nf90_collective))
    call handle_err(nf90_put_att(atmvarid, atmvarid_thetabar, 'long_name', &
        &'potential-temperature background'))
    call handle_err(nf90_put_att(atmvarid, atmvarid_thetabar, 'units', 'K'))

    ! squared buoyancy frequency
    call handle_err(nf90_def_var(atmvarid, 'N2', nf90_float, [dimidx, dimidy, &
        &dimidz, dimidt], atmvarid_n2))
    call handle_err(nf90_var_par_access(atmvarid, atmvarid_n2, nf90_collective))
    call handle_err(nf90_put_att(atmvarid, atmvarid_n2, 'long_name', 'squared &
        &buoyancy frequency'))
    call handle_err(nf90_put_att(atmvarid, atmvarid_n2, 'units', '1/s**2'))

    ! mass-weighted potential temperature
    call handle_err(nf90_def_var(atmvarid, 'P', nf90_float, [dimidx, dimidy, &
        &dimidz, dimidt], atmvarid_p))
    call handle_err(nf90_var_par_access(atmvarid, atmvarid_p, nf90_collective))
    call handle_err(nf90_put_att(atmvarid, atmvarid_p, 'long_name', &
        &'mass-weighted potential temperature'))
    call handle_err(nf90_put_att(atmvarid, atmvarid_p, 'units', 'K*kg/m**3'))

    ! density background
    call handle_err(nf90_def_var(atmvarid, 'rhobar', nf90_float, [dimidx, &
        &dimidy, dimidz, dimidt], atmvarid_rhobar))
    call handle_err(nf90_var_par_access(atmvarid, atmvarid_rhobar, &
        &nf90_collective))
    call handle_err(nf90_put_att(atmvarid, atmvarid_rhobar, 'long_name', &
        &'background density'))
    call handle_err(nf90_put_att(atmvarid, atmvarid_rhobar, 'units', 'kg/m**3'))

    ! density fluctuations
    if(prepare_restart .or. any(atmvarOut == 'rhop')) then
      call handle_err(nf90_def_var(atmvarid, 'rhop', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_rhop))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_rhop, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_rhop, 'long_name', &
          &'density fluctuation'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_rhop, 'units', 'kg/m**3'))
    end if

    ! zonal wind
    if(any(atmvarOut == 'u')) then
      call handle_err(nf90_def_var(atmvarid, 'u', nf90_float, [dimidx, dimidy, &
          &dimidz, dimidt], atmvarid_u))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_u, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_u, 'long_name', 'zonal &
          &wind'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_u, 'units', 'm/s'))
    end if
    if(prepare_restart .or. any(atmvarOut == 'us')) then
      call handle_err(nf90_def_var(atmvarid, 'us', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_us))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_us, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_us, 'long_name', 'zonal &
          &wind staggered'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_us, 'units', 'm/s'))
    end if

    ! meridional wind
    if(any(atmvarOut == 'v')) then
      call handle_err(nf90_def_var(atmvarid, 'v', nf90_float, [dimidx, dimidy, &
          &dimidz, dimidt], atmvarid_v))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_v, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_v, 'long_name', &
          &'meridional wind'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_v, 'units', 'm/s'))
    end if
    if(prepare_restart .or. any(atmvarOut == 'vs')) then
      call handle_err(nf90_def_var(atmvarid, 'vs', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_vs))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_vs, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_vs, 'long_name', &
          &'meridional wind staggered'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_vs, 'units', 'm/s'))
    end if

    ! vertical wind
    if(any(atmvarOut == 'w')) then
      call handle_err(nf90_def_var(atmvarid, 'w', nf90_float, [dimidx, dimidy, &
          &dimidz, dimidt], atmvarid_w))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_w, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_w, 'long_name', &
          &'vertical wind'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_w, 'units', 'm/s'))
    end if
    if(any(atmvarOut == 'ws')) then
      call handle_err(nf90_def_var(atmvarid, 'ws', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_ws))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_ws, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_ws, 'long_name', &
          &'vertical wind staggered'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_ws, 'units', 'm/s'))
    end if

    ! terrain-following vertical wind
    if(any(atmvarOut == 'wTFC')) then
      call handle_err(nf90_def_var(atmvarid, 'wTFC', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_wTFC))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_wTFC, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_wTFC, 'long_name', &
          &'terrain-following vertical wind'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_wTFC, 'units', 'm/s'))
    end if
    if(prepare_restart .or. any(atmvarOut == 'wsTFC')) then
      call handle_err(nf90_def_var(atmvarid, 'wsTFC', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_wsTFC))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_wsTFC, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_wsTFC, 'long_name', &
          &'terrain-following vertical wind staggered'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_wsTFC, 'units', 'm/s'))
    end if

    ! buoyancy
    if(any(atmvarOut == 'b')) then
      call handle_err(nf90_def_var(atmvarid, 'b', nf90_float, [dimidx, dimidy, &
          &dimidz, dimidt], atmvarid_b))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_b, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_b, 'long_name', &
          &'buoyancy'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_b, 'units', 'm/s**2'))
    end if

    ! potential-temperature fluctuations
    if(any(atmvarOut == 'thetap')) then
      call handle_err(nf90_def_var(atmvarid, 'thetap', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_thetap))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_thetap, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_thetap, 'long_name', &
          &'potential temperature fluctuations'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_thetap, 'units', 'K'))
    end if

    ! Exner-pressure fluctuations.
    if(prepare_restart .or. any(atmvarOut == 'pip')) then
      call handle_err(nf90_def_var(atmvarid, 'pip', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_pip))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_pip, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_pip, 'long_name', &
          &'Exner-pressure fluctuations'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_pip, 'units', 'none'))
    end if

  end subroutine create_netCDF

  !-------------------------------------------------------------------------

  subroutine write_netCDF_background(iOutput)

    ! write the x-, y-, and z-axis to the netCDF file
    ! write background fields pStratTFC, thetaStratTFC, rhoStratTFC, bvsStratTFC
    ! to the netCDF file

    implicit none

    integer, intent(in) :: iOutput

    integer :: ix, jy, kz

    if(sizeX == 1) then
      startxNC = 1
      countxNC = 1
    else
      startxNC = is + nbx
      countxNC = nx
    end if

    if(sizeY == 1) then
      startyNC = 1
      countyNC = 1
    else
      startyNC = js + nby
      countyNC = ny
    end if

    startNC = (/startxNC, startyNC, 1, iOutput + 1/)
    countNC = (/countxNC, countyNC, nz, 1/)

    ! x-axis
    if(iOutput == 0) then
      if(sizeX == 1) then
        call handle_err(nf90_put_var(ncid, ncidx, x(1:sizeX) * lRef, start &
            &= [1]))
      else
        call handle_err(nf90_put_var(ncid, ncidx, x(1:sizeX) * lRef, start &
            &= [1], count = [sizeX]))
      end if

      ! y-axis
      if(sizeY == 1) then
        call handle_err(nf90_put_var(ncid, ncidy, y(1:sizeY) * lRef, start &
            &= [1]))
      else
        call handle_err(nf90_put_var(ncid, ncidy, y(1:sizeY) * lRef, start &
            &= [1], count = [sizeY]))
      end if

      ! initial time
      call handle_err(nf90_put_var(ncid, ncidt, 0.0, start = [1]))
    end if

    ! z-axis
    call handle_err(nf90_put_var(ncid, ncidz, zTFC(1:nx, 1:ny, 1:sizeZ) &
        &* lRef, start = startNC, count = countNC))

    ! mass-weighted potential temperature
    call handle_err(nf90_put_var(atmvarid, atmvarid_p, pStratTFC(1:nx, 1:ny, &
        &1:sizeZ) * rhoRef * thetaRef, start = startNC, count = countNC), &
        &'save P')

    ! potential-temperature background
    call handle_err(nf90_put_var(atmvarid, atmvarid_thetabar, &
        &thetaStratTFC(1:nx, 1:ny, 1:sizeZ) * thetaRef, start = startNC, count &
        &= countNC), 'save thetabar')

    ! density background
    call handle_err(nf90_put_var(atmvarid, atmvarid_rhobar, rhoStratTFC(1:nx, &
        &1:ny, 1:sizeZ) * rhoRef, start = startNC, count = countNC), 'save &
        &rhobar')

    ! squared buoyancy frequency
    call handle_err(nf90_put_var(atmvarid, atmvarid_n2, bvsStratTFC(1:nx, &
        &1:ny, 1:sizeZ) / tRef ** 2.0, start = startNC, count = countNC), &
        &'save N2')

  end subroutine write_netCDF_background

  !-------------------------------------------------------------------------

  subroutine write_netCDF(iOutput, iTime, time, cpuTime, var)

    implicit none

    integer, intent(inout) :: iOutput
    integer, intent(in) :: iTime
    real, intent(in) :: time
    real, intent(in) :: cpuTime

    type(var_type), intent(in) :: var

    real :: time_dim

    real :: cpuTimeLoc
    integer :: days, hours, minutes, seconds

    real, dimension(1:nx, 1:ny, 1:nz) :: rho, rhop
    real, dimension(1:nx, 1:ny, 1:nz) :: buoyancy
    real, dimension(1:nx, 1:ny, 1:nz) :: thetap
    real, dimension(1:nx, 1:ny, 0:nz) :: wOut

    integer :: ix, jy, kz

    time_dim = time * tRef

    cpuTimeLoc = cpuTime
    days = floor(cpuTimeLoc / 86400.0)
    cpuTimeLoc = cpuTimeLoc - days * 86400.0
    hours = floor(cpuTimeLoc / 3600.0)
    cpuTimeLoc = cpuTimeLoc - hours * 3600.0
    minutes = floor(cpuTimeLoc / 60.0)
    cpuTimeLoc = cpuTimeLoc - minutes * 60.0
    seconds = floor(cpuTimeLoc)

    if(master) then
      print *, "------------------------------------------------"
      print *, "Output into file pincflow_data_out.nc"
      print *, " "
      write(*, fmt = "(a25,i15)") " at time step = ", iTime
      write(*, fmt = "(a25,f15.1,a8)") " at physical time = ", time_dim, " &
          &seconds"
      write(*, fmt = "(a25,4x,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)") "CPU time = ", &
          &days, "-", hours, ":", minutes, ":", seconds
      print *, "------------------------------------------------"
    end if

    if(sizeX == 1) then
      startxNC = 1
      countxNC = 1
    else
      startxNC = is + nbx
      countxNC = nx
    end if

    if(sizeY == 1) then
      startyNC = 1
      countyNC = 1
    else
      startyNC = js + nby
      countyNC = ny
    end if

    startNC = (/startxNC, startyNC, 1, iOutput + 1/)
    countNC = (/countxNC, countyNC, nz, 1/)

    call write_netCDF_background(iOutput)

    rhop = var%rho(1:nx, 1:ny, 1:nz)

    rho = rhop + rhoStratTFC(1:nx, 1:ny, 1:nz)
    buoyancy = g * (pStratTFC(1:nx, 1:ny, 1:nz) / rho - thetaStratTFC(1:nx, &
        &1:ny, 1:nz)) / thetaStratTFC(1:nx, 1:ny, 1:nz)
    thetap = buoyancy / g * thetaStratTFC(1:nx, 1:ny, 1:nz)

    do ix = 1, nx
      do jy = 1, ny
        do kz = 0, nz
          wOut(ix, jy, kz) = vertWindTFC(ix, jy, kz, var)
        end do
      end do
    end do

    ! save time
    call handle_err(nf90_put_var(ncid, ncidt, time_dim, start = [iOutput + 1]))

    !--------------------------------
    ! save variables in atmvar group
    !--------------------------------

    ! density fluctuations
    if(prepare_restart .or. any(atmvarOut == 'rhop')) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_rhop, rhop(1:nx, 1:ny, &
          &1:nz) * rhoRef, start = startNC, count = countNC), 'save rhop')
    end if

    ! zonal wind
    if(any(atmvarOut == 'u')) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_u, (var%u(0:nx - 1, &
          &1:ny, 1:nz) + var%u(1:nx, 1:ny, 1:nz)) * uRef / 2., start &
          &= startNC, count = countNC), 'save u')
    end if
    if(prepare_restart .or. any(atmvarOut == 'us')) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_us, var%u(1:nx, 1:ny, &
          &1:nz) * uRef, start = startNC, count = countNC), 'save us')
    end if

    ! meridional wind
    if(any(atmvarOut == 'v')) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_v, (var%v(1:nx, 0:ny &
          &- 1, 1:nz) + var%v(1:nx, 1:ny, 1:nz)) * uRef / 2., start = startNC, &
          &count = countNC), 'save v')
    end if
    if(prepare_restart .or. any(atmvarOut == 'vs')) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_vs, var%v(1:nx, 1:ny, &
          &1:nz) * uRef, start = startNC, count = countNC), 'save vs')
    end if

    ! vertical wind
    if(any(atmvarOut == 'w')) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_w, (wOut(1:nx, 1:ny, &
          &0:nz - 1) + wOut(1:nx, 1:ny, 1:nz)) * uRef / 2., start = startNC, &
          &count = countNC), 'save w')
    end if
    if(any(atmvarOut == 'ws')) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_ws, wOut(1:nx, 1:ny, &
          &1:nz) * uRef, start = startNC, count = countNC), 'save ws')
    end if
    if(any(atmvarOut == 'wTFC')) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_wTFC, (var%w(1:nx, 1:ny, &
          &0:nz - 1) + var%w(1:nx, 1:ny, 1:nz)) * uRef / 2., start = startNC, &
          &count = countNC), 'save wTFC')
    end if
    if(prepare_restart .or. any(atmvarOut == 'wsTFC')) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_wsTFC, var%w(1:nx, 1:ny, &
          &1:nz) * uRef, start = startNC, count = countNC), 'save wsTFC')
    end if

    ! buoyancy
    if(any(atmvarOut == 'b')) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_b, buoyancy(1:nx, 1:ny, &
          &1:nz), start = startNC, count = countNC), 'save b')
    end if

    ! potential-temperature fluctuations
    if(any(atmvarOut == 'thetap')) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_thetap, thetap(1:nx, &
          &1:ny, 1:nz) * thetaRef, start = startNC, count = countNC), 'save &
          &thetap')
    end if

    ! Exner-pressure fluctuations
    if(prepare_restart .or. any(atmvarOut == 'pip')) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_pip, var%pi(1:nx, 1:ny, &
          &1:nz), start = startNC, count = countNC), 'save pi')
    end if

    call write_topography_netCDF(iOutput)

    iOutput = iOutput + 1

  end subroutine write_netCDF

  !-----------------------------------------------

  subroutine write_topography_netCDF(iOutput)

    implicit none

    integer, intent(in) :: iOutput

    if(sizeX == 1) then
      startxNC = 1
      countxNC = 1
    else
      startxNC = is + nbx
      countxNC = nx
    end if

    if(sizeY == 1) then
      startyNC = 1
      countyNC = 1
    else
      startyNC = js + nby
      countyNC = ny
    end if

    starttNC = iOutput + 1
    counttNC = 1

    ! Write large-scale topography.
    call handle_err(nf90_put_var(ncid, ncidhm, topography_surface(1:nx, 1:ny) &
        &* lRef, start = [startxNC, startyNC, starttNC], count = [countxNC, &
        &countyNC, counttNC]))

  end subroutine write_topography_netCDF

  !-------------------------------------------------------------------------

  subroutine close_netCDF

    implicit none

    call handle_err(nf90_close(ncid))

  end subroutine close_netCDF

  !-------------------------------------------------------------------------

  subroutine read_netCDF(timeStart, var, time)

    implicit none

    integer, intent(in) :: timeStart

    type(var_type), intent(inout) :: var

    real, dimension(1:nx, 1:ny, 1:nz) :: rhop

    real, intent(out), optional :: time

    integer :: Ntimesteps

    ! Open file.
    call handle_err(nf90_open(filenameIn, ior(nf90_nowrite, nf90_mpiio), ncid, &
        &comm = comm, info = mpi_info_null))

    ! Get group IDs.
    call handle_err(nf90_inq_ncid(ncid, "atmvar", atmvarid))

    ! Get last record.
    call handle_err(nf90_inq_dimid(ncid, 't', dimidt))
    call handle_err(nf90_inquire_dimension(ncid, dimidt, len = Ntimesteps))

    if(sizeX == 1) then
      startxNC = 1
      countxNC = 1
    else
      startxNC = is + nbx
      countxNC = nx
    end if

    if(sizeY == 1) then
      startyNC = 1
      countyNC = 1
    else
      startyNC = js + nby
      countyNC = ny
    end if

    if(timeStart >= 0) then
      startNC = (/startxNC, startyNC, 1, timeStart + 1/)
    else
      startNC = (/startxNC, startyNC, 1, Ntimesteps + timeStart + 1/)
    end if
    countNC = (/countxNC, countyNC, nz, 1/)

    ! Read time.
    call handle_err(nf90_inq_varid(ncid, 't', ncidt))
    if(timeStart >= 0) then
      call handle_err(nf90_get_var(ncid, ncidt, time, start = [timeStart + 1]))
    else
      call handle_err(nf90_get_var(ncid, ncidt, time, start = [Ntimesteps &
          &+ timeStart + 1]))
    end if

    if(master) then
      print *, "Reading data at time = ", time, 's'
    end if
    time = time / tRef

    ! Read density fluctuations.
    call handle_err(nf90_inq_varid(atmvarid, 'rhop', atmvarid_rhop))
    call handle_err(nf90_var_par_access(atmvarid, atmvarid_rhop, &
        &nf90_collective))
    call handle_err(nf90_get_var(atmvarid, atmvarid_rhop, rhop(1:nx, 1:ny, &
        &1:nz), start = startNC, count = countNC))
    var%rho(1:nx, 1:ny, 1:nz) = rhop(1:nx, 1:ny, 1:nz) / rhoRef

    ! Read zonal wind.
    call handle_err(nf90_inq_varid(atmvarid, 'us', atmvarid_us))
    call handle_err(nf90_var_par_access(atmvarid, atmvarid_us, nf90_collective))
    call handle_err(nf90_get_var(atmvarid, atmvarid_us, var%u(1:nx, 1:ny, &
        &1:nz), start = startNC, count = countNC))
    var%u(1:nx, 1:ny, 1:nz) = var%u(1:nx, 1:ny, 1:nz) / uRef

    ! Read meridional wind.
    call handle_err(nf90_inq_varid(atmvarid, 'vs', atmvarid_vs))
    call handle_err(nf90_var_par_access(atmvarid, atmvarid_vs, nf90_collective))
    call handle_err(nf90_get_var(atmvarid, atmvarid_vs, var%v(1:nx, 1:ny, &
        &1:nz), start = startNC, count = countNC))
    var%v(1:nx, 1:ny, 1:nz) = var%v(1:nx, 1:ny, 1:nz) / uRef

    ! Read vertical wind.
    call handle_err(nf90_inq_varid(atmvarid, 'wsTFC', atmvarid_wsTFC))
    call handle_err(nf90_var_par_access(atmvarid, atmvarid_wsTFC, &
        &nf90_collective))
    call handle_err(nf90_get_var(atmvarid, atmvarid_wsTFC, var%w(1:nx, 1:ny, &
        &1:nz), start = startNC, count = countNC))
    var%w(1:nx, 1:ny, 1:nz) = var%w(1:nx, 1:ny, 1:nz) / uRef

    ! Read Exner-pressure fluctuations.
    call handle_err(nf90_inq_varid(atmvarid, 'pip', atmvarid_pip))
    call handle_err(nf90_var_par_access(atmvarid, atmvarid_pip, &
        &nf90_collective))
    call handle_err(nf90_get_var(atmvarid, atmvarid_pip, var%pi(1:nx, 1:ny, &
        &1:nz), start = startNC, count = countNC))

    if(master) then
      print *, "reading file complete."
    end if

    call handle_err(nf90_close(ncid))

  end subroutine read_netCDF

  !-------------------------------------------------------------------------

  subroutine read_topography_netCDF(timeStart)

    implicit none

    integer, intent(in) :: timeStart
    integer :: nt

    ! Open file.
    call handle_err(nf90_open(filenameIn, ior(nf90_nowrite, nf90_mpiio), ncid, &
        &comm = comm, info = mpi_info_null))

    ! Get last record.
    call handle_err(nf90_inq_dimid(ncid, 't', dimidt))
    call handle_err(nf90_inquire_dimension(ncid, dimidt, len = nt))

    if(sizeX == 1) then
      startxNC = 1
      countxNC = 1
    else
      startxNC = is + nbx
      countxNC = nx
    end if

    if(sizeY == 1) then
      startyNC = 1
      countyNC = 1
    else
      startyNC = js + nby
      countyNC = ny
    end if

    if(timeStart >= 0) then
      starttNC = timeStart + 1
      counttNC = 1
    else
      starttNC = nt + timeStart + 1
      counttNC = 1
    endif

    ! Read large-scale topography (non-dimensionalization is done in
    ! setup_topography).
    call handle_err(nf90_inq_varid(ncid, 'hm', ncidhm))
    call handle_err(nf90_var_par_access(ncid, ncidhm, nf90_collective))
    call handle_err(nf90_get_var(ncid, ncidhm, topography_surface(1:nx, 1:ny), &
        &start = [startxNC, startyNC, starttNC], count = [countxNC, countyNC, &
        &counttNC]))

    call handle_err(nf90_close(ncid))

  end subroutine read_topography_netCDF

  !-------------------------------------------------------------------------

  subroutine handle_err(errcode, operation)

    implicit none
    integer, intent(in) :: errcode

    character(len = *), intent(in), optional :: operation

    if(errcode /= nf90_noerr) then
      print *, "Error encountered during ", operation
      print *, nf90_strerror(errcode)
      stop 2
    endif
  end subroutine handle_err

end module output_netCDF_module
