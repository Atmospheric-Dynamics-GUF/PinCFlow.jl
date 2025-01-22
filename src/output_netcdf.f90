module output_netCDF_module

  use type_module
  use atmosphere_module
  use sizeof_module
  use netcdf
  use mpi

  implicit none

  public

  ! netCDF parameters
  character(len = *), parameter :: filenameOut = 'pincflow_data_out.nc'
  character(len = *), parameter :: filenameIn = 'pincflow_data_in.nc'
  integer :: mode_flagNC ! mode_flag for creating the netCDF file
  integer :: ncid ! id of the netCDF file
  integer :: atmvarid, rayvarid, icevarid ! group ids
  integer :: rayvolid ! group id for ray volumes

  ! start and count values for nf90_put_var(...)
  integer :: startxNC, startyNC, countxNC, countyNC
  integer :: starthNC, starttNC, counthNC, counttNC
  integer, dimension(4) :: startNC, countNC
  integer, dimension(5) :: startNCray, countNCray

  ! dimensions of the grid
  integer :: dim

  ! ids of the dimensions
  integer :: dimidx, dimidy, dimidz, dimidt ! x, y, z, t
  integer :: dimidr, dimidzRay ! nray_max and ray z dimension id
  integer :: dimidh ! dimension id for topography spectrum

  ! variable ids of the axis
  integer :: ncidx, ncidy, ncidz, ncidt ! x, y, z, t

  ! variable ids of the topography
  integer :: ncidhm ! large-scale topography
  integer :: ncidhw, ncidkh, ncidlh ! small-scale topography

  ! variable ids in the group atmvar
  integer :: atmvarid_u, atmvarid_v, atmvarid_w ! u, v, w wind
  integer :: atmvarid_wTFC ! terrain-following w
  integer :: atmvarid_us, atmvarid_vs, atmvarid_ws ! u, v, w on staggered
  integer :: atmvarid_wsTFC ! terrain-following ws
  integer :: atmvarid_b ! buoyancy
  integer :: atmvarid_rhobar, atmvarid_p, atmvarid_n2, atmvarid_thetabar
  ! rhoStrat, pStrat, bvsStrat, thetaStrat
  integer :: atmvarid_rhop ! density fluctuations rho'
  integer :: atmvarid_thetap ! potential temperature fluctuations theta'
  integer :: atmvarid_tmr, atmvarid_tmrd ! tracer and tracer diff.
  integer :: atmvarid_pip ! Exner-pressure fluctuations pi'
  integer :: atmvarid_dsc ! dynamic Smagorinsky coefficient

  ! variable ids in the group rayvar
  integer :: rayvarid_dudt ! zonal wind-tendency
  integer :: rayvarid_dvdt ! meridional-wind tendency
  integer :: rayvarid_dthetadt ! potential-temperature tendency
  integer :: rayvarid_uw ! zonal-vertical momentum flux
  integer :: rayvarid_vw ! meridional-vertical momentum flux
  integer :: rayvarid_e ! gravity wave energy
  integer :: rayvarid_lsphase ! large-scale wave phase
  integer :: rayvarid_lobhat ! leading-order buoyancy wave amplitude
  integer :: rayvarid_nobhat ! next-order buoyancy wave amplitude
  integer :: rayvarid_tfrclot ! leading-order tracer flux convergence
  integer :: rayvarid_tfrclou ! leading-order zonal tracer flux
  integer :: rayvarid_tfrcnow ! next-order vertical tracer flux
  integer :: rayvarid_tfrcnot ! next-order tracer flux convergence
  integer :: rayvarid_gwh ! gravity-wave heating

  ! variable ids in the group rayvol
  integer :: rayvolid_x, rayvolid_y, rayvolid_z ! ray position
  integer :: rayvolid_k, rayvolid_l, rayvolid_m ! ray wave vector
  integer :: rayvolid_omega ! intrinsic frequency
  integer :: rayvolid_dkray, rayvolid_dlray, rayvolid_dmray
  integer :: rayvolid_dxray, rayvolid_dyray, rayvolid_dzray
  integer :: rayvolid_dens
  integer :: rayvolid_dphi

  contains

  subroutine create_netCDF

    implicit none

    ! date, time, and zone when file was created
    character(8) :: dateNC
    character(10) :: timeNC
    character(5) :: zoneNC

    call date_and_time(date = dateNC, time = timeNC, zone = zoneNC)

    dim = 1
    if(sizeX > 1) dim = dim + 1
    if(sizeY > 1) dim = dim + 1

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
    if(rayTracer) then
      call handle_err(nf90_put_att(ncid, nf90_global, 'Comment', 'Gravity &
          &waves are parameterized using MS-GWaM'))
    end if
    if(topography) then
      call handle_err(nf90_put_att(ncid, nf90_global, 'Comment', 'The grid is &
          &terrain-following'))
    end if

    ! define the x-, y-, z-, and t-dimensions
    call handle_err(nf90_def_dim(ncid, 'x', sizeX, dimidx))
    call handle_err(nf90_def_dim(ncid, 'y', sizeY, dimidy))
    call handle_err(nf90_def_dim(ncid, 'z', sizeZ, dimidz))
    call handle_err(nf90_def_dim(ncid, 't', nf90_unlimited, dimidt))

    if(rayTracer .and. (prepare_restart .or. saverayvols)) then
      call handle_err(nf90_def_dim(ncid, 'nray_max', nray_max, dimidr))
      call handle_err(nf90_def_dim(ncid, 'zRay', sizeZ + 2, dimidzRay))
    end if

    if(rayTracer .and. case_wkb == 3) then
      call handle_err(nf90_def_dim(ncid, 'h', nwm, dimidh))
    end if

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
    if(topography) then
      call handle_err(nf90_def_var(ncid, 'z', nf90_float, [dimidx, dimidy, &
          &dimidz, dimidt], ncidz))
    else
      call handle_err(nf90_def_var(ncid, 'z', nf90_float, [dimidz], ncidz))
    end if
    call handle_err(nf90_var_par_access(ncid, ncidz, nf90_collective))
    call handle_err(nf90_put_att(ncid, ncidz, 'units', 'm'))

    ! t axis
    call handle_err(nf90_def_var(ncid, 't', nf90_float, [dimidt], ncidt))
    call handle_err(nf90_var_par_access(ncid, ncidt, nf90_collective))
    call handle_err(nf90_put_att(ncid, ncidt, 'units', 's'))

    ! large-scale topography
    if(topography) then
      call handle_err(nf90_def_var(ncid, 'hm', nf90_float, [dimidx, dimidy, &
          &dimidt], ncidhm))
      call handle_err(nf90_var_par_access(ncid, ncidhm, nf90_collective))
      call handle_err(nf90_put_att(ncid, ncidhm, 'units', 'm'))
    end if

    ! small-scale topography
    if(rayTracer .and. case_wkb == 3) then
      ! zonal wavenumbers
      call handle_err(nf90_def_var(ncid, 'kh', nf90_float, [dimidx, dimidy, &
          &dimidh, dimidt], ncidkh))
      call handle_err(nf90_var_par_access(ncid, ncidkh, nf90_collective))
      call handle_err(nf90_put_att(ncid, ncidkh, 'long_name', &
          &'zonal wavenumbers of the small-scale topography spectrum'))
      call handle_err(nf90_put_att(ncid, ncidkh, 'units', '1/m'))

      ! meridional wavenumbers
      call handle_err(nf90_def_var(ncid, 'lh', nf90_float, [dimidx, dimidy, &
          &dimidh, dimidt], ncidlh))
      call handle_err(nf90_var_par_access(ncid, ncidlh, nf90_collective))
      call handle_err(nf90_put_att(ncid, ncidlh, 'long_name', &
          &'meridional wavenumbers of the small-scale topography spectrum'))
      call handle_err(nf90_put_att(ncid, ncidlh, 'units', '1/m'))

      ! wave amplitudes
      call handle_err(nf90_def_var(ncid, 'hw', nf90_float, [dimidx, dimidy, &
          &dimidh, dimidt], ncidhw))
      call handle_err(nf90_var_par_access(ncid, ncidhw, nf90_collective))
      call handle_err(nf90_put_att(ncid, ncidhw, 'long_name', &
          &'amplitudes of the small-scale topography spectrum'))
      call handle_err(nf90_put_att(ncid, ncidhw, 'units', 'm'))
    end if

    !---------------------------------------------------
    ! create the groups (atmvar, rayvar, icevar, rayvol)
    !---------------------------------------------------

    call handle_err(nf90_def_grp(ncid, 'atmvar', atmvarid))
    if(rayTracer) then
      call handle_err(nf90_def_grp(ncid, 'rayvar', rayvarid))
      if(prepare_restart .or. saverayvols) then
        call handle_err(nf90_def_grp(ncid, 'rayvol', rayvolid))
      end if
    end if
    if(include_ice) then
      call handle_err(nf90_def_grp(ncid, 'icevar', icevarid))
    end if

    !-------------------------------------
    ! define the variables in atmvar group
    !-------------------------------------

    ! potential-temperature background
    if(topography) then
      call handle_err(nf90_def_var(atmvarid, 'thetabar', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_thetabar))
    else
      call handle_err(nf90_def_var(atmvarid, 'thetabar', nf90_float, [dimidz], &
          &atmvarid_thetabar))
    end if
    call handle_err(nf90_var_par_access(atmvarid, atmvarid_thetabar, &
        &nf90_collective))
    call handle_err(nf90_put_att(atmvarid, atmvarid_thetabar, 'long_name', &
        &'potential-temperature background'))
    call handle_err(nf90_put_att(atmvarid, atmvarid_thetabar, 'units', 'K'))

    ! squared buoyancy frequency
    if(topography) then
      call handle_err(nf90_def_var(atmvarid, 'N2', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_n2))
    else
      call handle_err(nf90_def_var(atmvarid, 'N2', nf90_float, [dimidz], &
          &atmvarid_n2))
    end if
    call handle_err(nf90_var_par_access(atmvarid, atmvarid_n2, nf90_collective))
    call handle_err(nf90_put_att(atmvarid, atmvarid_n2, 'long_name', 'squared &
        &buoyancy frequency'))
    call handle_err(nf90_put_att(atmvarid, atmvarid_n2, 'units', '1/s**2'))

    ! mass-weighted potential temperature
    if(topography) then
      call handle_err(nf90_def_var(atmvarid, 'P', nf90_float, [dimidx, dimidy, &
          &dimidz, dimidt], atmvarid_p))
    else
      call handle_err(nf90_def_var(atmvarid, 'P', nf90_float, [dimidz, &
          &dimidt], atmvarid_p))
    end if
    call handle_err(nf90_var_par_access(atmvarid, atmvarid_p, nf90_collective))
    call handle_err(nf90_put_att(atmvarid, atmvarid_p, 'long_name', &
        &'mass-weighted potential temperature'))
    call handle_err(nf90_put_att(atmvarid, atmvarid_p, 'units', 'K*kg/m**3'))

    ! density background
    if(topography) then
      call handle_err(nf90_def_var(atmvarid, 'rhobar', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_rhobar))
    else
      call handle_err(nf90_def_var(atmvarid, 'rhobar', nf90_float, [dimidz], &
          &atmvarid_rhobar))
    end if
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
    if((.not. topography .and. prepare_restart) .or. any(atmvarOut == 'ws')) &
        &then
      call handle_err(nf90_def_var(atmvarid, 'ws', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_ws))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_ws, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_ws, 'long_name', &
          &'vertical wind staggered'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_ws, 'units', 'm/s'))
    end if

    ! terrain-following vertical wind
    if(topography) then
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

    ! dynamic Smagorinsky coefficient
    if(turbScheme .and. any(atmvarOut == 'dsc')) then
      call handle_err(nf90_def_var(atmvarid, 'dsc', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_dsc))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_dsc, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_dsc, 'long_name', &
          &'Dynamic Smagorinsky coefficient'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_dsc, 'units', 'm**2/s'))
    end if

    ! tracer mixing ratio (difference)
    if(include_tracer .and. any(atmvarOut == 'tmr')) then
      call handle_err(nf90_def_var(atmvarid, 'tmr', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_tmr))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_tmr, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_tmr, 'long_name', &
          &'tracer mixing ratio'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_tmr, 'units', 'none'))
    end if
    if(include_tracer .and. (prepare_restart .or. any(atmvarOut == 'tmrd'))) &
        &then
      call handle_err(nf90_def_var(atmvarid, 'tmrd', nf90_float, [dimidx, &
          &dimidy, dimidz, dimidt], atmvarid_tmrd))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_tmrd, &
          &nf90_collective))
      call handle_err(nf90_put_att(atmvarid, atmvarid_tmrd, 'long_name', &
          &'tracer mixing ratio difference'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_tmrd, 'units', 'none'))
      call handle_err(nf90_put_att(atmvarid, atmvarid_tmrd, 'description', &
          &'Change in tracer mixing ratio from initial distribution'))
    end if

    !-------------------------------------
    ! define the variables in rayvar group
    !-------------------------------------

    if(rayTracer) then
      if(include_tracer) then
        ! leading-order zonal tracer flux
        if(any(rayvarOut == 'tfrclou')) then
          call handle_err(nf90_def_var(rayvarid, 'tfrclou', nf90_float, &
              &[dimidx, dimidy, dimidz, dimidt], rayvarid_tfrclou))
          call handle_err(nf90_var_par_access(rayvarid, rayvarid_tfrclou, &
              &nf90_collective))
          call handle_err(nf90_put_att(rayvarid, rayvarid_tfrclou, &
              &'long_name', 'leading-order zonal tracer flux'))
          call handle_err(nf90_put_att(rayvarid, rayvarid_tfrclou, 'units', &
              &'1/s'))
        end if

        ! leading-order tracer flux convergence
        if(any(rayvarOut == 'tfrclot')) then
          call handle_err(nf90_def_var(rayvarid, 'tfrclot', nf90_float, &
              &[dimidx, dimidy, dimidz, dimidt], rayvarid_tfrclot))
          call handle_err(nf90_var_par_access(rayvarid, rayvarid_tfrclot, &
              &nf90_collective))
          call handle_err(nf90_put_att(rayvarid, rayvarid_tfrclot, &
              &'long_name', 'leading-order tracer flux convergence'))
          call handle_err(nf90_put_att(rayvarid, rayvarid_tfrclot, 'units', &
              &'1/s'))
        end if

        ! next-order tracer flux convergence
        if(any(rayvarOut == 'tfrcnot')) then
          call handle_err(nf90_def_var(rayvarid, 'tfrcnot', nf90_float, &
              &[dimidx, dimidy, dimidz, dimidt], rayvarid_tfrcnot))
          call handle_err(nf90_var_par_access(rayvarid, rayvarid_tfrcnot, &
              &nf90_collective))
          call handle_err(nf90_put_att(rayvarid, rayvarid_tfrcnot, &
              &'long_name', 'next-order tracer flux convergence'))
          call handle_err(nf90_put_att(rayvarid, rayvarid_tfrcnot, 'units', &
              &'1/s'))
        end if

        ! next-order vertical tracer flux
        if(any(rayvarOut == 'tfrcnow')) then
          call handle_err(nf90_def_var(rayvarid, 'tfrcnow', nf90_float, &
              &[dimidx, dimidy, dimidz, dimidt], rayvarid_tfrcnow))
          call handle_err(nf90_var_par_access(rayvarid, rayvarid_tfrcnow, &
              &nf90_collective))
          call handle_err(nf90_put_att(rayvarid, rayvarid_tfrcnow, &
              &'long_name', 'next-order vertical tracer flux'))
          call handle_err(nf90_put_att(rayvarid, rayvarid_tfrcnow, 'units', &
              &'m/s'))
        end if
      end if

      ! zonal-wind tendency
      if(any(rayvarOut == 'dudt')) then
        call handle_err(nf90_def_var(rayvarid, 'dudt', nf90_float, [dimidx, &
            &dimidy, dimidz, dimidt], rayvarid_dudt))
        call handle_err(nf90_var_par_access(rayvarid, rayvarid_dudt, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvarid, rayvarid_dudt, 'long_name', &
            &'zonal-wind tendency'))
        call handle_err(nf90_put_att(rayvarid, rayvarid_dudt, 'units', &
            &'m/s**2'))
      end if

      ! meridional-wind tendency
      if(any(rayvarOut == 'dvdt')) then
        call handle_err(nf90_def_var(rayvarid, 'dvdt', nf90_float, [dimidx, &
            &dimidy, dimidz, dimidt], rayvarid_dvdt))
        call handle_err(nf90_var_par_access(rayvarid, rayvarid_dvdt, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvarid, rayvarid_dvdt, 'long_name', &
            &'meridional-wind tendency'))
        call handle_err(nf90_put_att(rayvarid, rayvarid_dvdt, 'units', &
            &'m/s**2'))
      end if

      ! potential-temperature tendency
      if(any(rayvarOut == 'dthetadt')) then
        call handle_err(nf90_def_var(rayvarid, 'dthetadt', nf90_float, &
            &[dimidx, dimidy, dimidz, dimidt], rayvarid_dthetadt))
        call handle_err(nf90_var_par_access(rayvarid, rayvarid_dthetadt, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvarid, rayvarid_dthetadt, 'long_name', &
            &'potential-temperature tendency'))
        call handle_err(nf90_put_att(rayvarid, rayvarid_dthetadt, 'units', &
            &'K/s'))
      end if

      ! zonal-vertical momentum flux
      if(any(rayvarOut == 'uw')) then
        call handle_err(nf90_def_var(rayvarid, 'uw', nf90_float, [dimidx, &
            &dimidy, dimidz, dimidt], rayvarid_uw))
        call handle_err(nf90_var_par_access(rayvarid, rayvarid_uw, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvarid, rayvarid_uw, 'long_name', &
            &'zonal-vertical momentum flux'))
        call handle_err(nf90_put_att(rayvarid, rayvarid_uw, 'units', &
            &'m**2/s**2'))
      end if

      ! meridional-vertical momentum flux
      if(any(rayvarOut == 'vw')) then
        call handle_err(nf90_def_var(rayvarid, 'vw', nf90_float, [dimidx, &
            &dimidy, dimidz, dimidt], rayvarid_vw))
        call handle_err(nf90_var_par_access(rayvarid, rayvarid_vw, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvarid, rayvarid_vw, 'long_name', &
            &'meridional-vertical momentum flux'))
        call handle_err(nf90_put_att(rayvarid, rayvarid_vw, 'units', &
            &'m**2/s**2'))
      end if

      ! gravity-wave energy density
      if(any(rayvarOut == 'E')) then
        call handle_err(nf90_def_var(rayvarid, 'E', nf90_float, [dimidx, &
            &dimidy, dimidz, dimidt], rayvarid_e))
        call handle_err(nf90_var_par_access(rayvarid, rayvarid_e, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvarid, rayvarid_e, 'long_name', &
            &'gravity-wave energy density'))
        call handle_err(nf90_put_att(rayvarid, rayvarid_e, 'units', &
            &'kg/m/s**2'))
      end if

      ! large-scale wave phase
      if(any(rayvarOut == 'lsphase')) then
        call handle_err(nf90_def_var(rayvarid, 'lsphase', nf90_float, [dimidx, &
            &dimidy, dimidz, dimidt], rayvarid_lsphase))
        call handle_err(nf90_var_par_access(rayvarid, rayvarid_lsphase, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvarid, rayvarid_lsphase, 'long_name', &
            &'large-scale wave phase'))
      end if

      ! leading-order buoyancy amplitude
      if(any(rayvarOut == 'lobhat')) then
        call handle_err(nf90_def_var(rayvarid, 'lobhat', nf90_float, [dimidx, &
            &dimidy, dimidz, dimidt], rayvarid_lobhat))
        call handle_err(nf90_var_par_access(rayvarid, rayvarid_lobhat, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvarid, rayvarid_lobhat, 'long_name', &
            &'abs. magn. leading-order buoyancy wave amplitude'))
        call handle_err(nf90_put_att(rayvarid, rayvarid_lobhat, 'units', &
            &'m/s**2'))
      end if

      ! next-order buoyancy amplitude
      if(any(rayvarOut == 'nobhat')) then
        call handle_err(nf90_def_var(rayvarid, 'nobhat', nf90_float, [dimidx, &
            &dimidy, dimidz, dimidt], rayvarid_nobhat))
        call handle_err(nf90_var_par_access(rayvarid, rayvarid_nobhat, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvarid, rayvarid_nobhat, 'long_name', &
            &'abs. magn. next-order buoyancy wave amplitude'))
        call handle_err(nf90_put_att(rayvarid, rayvarid_nobhat, 'units', &
            &'m/s**2'))
      end if

      ! gravity-wave heating
      if(any(rayvarOut == 'gwh')) then
        call handle_err(nf90_def_var(rayvarid, 'gwh', nf90_float, [dimidx, &
            &dimidy, dimidz, dimidt], rayvarid_gwh))
        call handle_err(nf90_var_par_access(rayvarid, rayvarid_gwh, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvarid, rayvarid_gwh, 'long_name', &
            &'gravity-wave heating'))
        call handle_err(nf90_put_att(rayvarid, rayvarid_gwh, 'units', &
            &'K*kg/m**3/s'))
      end if

      if(prepare_restart .or. saverayvols) then
        ! x-position of ray volume
        call handle_err(nf90_def_var(rayvolid, 'ray_x', nf90_float, [dimidr, &
            &dimidx, dimidy, dimidzRay, dimidt], rayvolid_x))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_x, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_x, 'long_name', 'ray &
            &volume position in x-direction'))
        call handle_err(nf90_put_att(rayvolid, rayvolid_x, 'units', 'm'))

        ! y-position of ray volume
        call handle_err(nf90_def_var(rayvolid, 'ray_y', nf90_float, [dimidr, &
            &dimidx, dimidy, dimidzRay, dimidt], rayvolid_y))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_y, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_y, 'long_name', 'ray &
            &volume position in y-direction'))
        call handle_err(nf90_put_att(rayvolid, rayvolid_y, 'units', 'm'))

        ! z-position of ray volume
        call handle_err(nf90_def_var(rayvolid, 'ray_z', nf90_float, [dimidr, &
            &dimidx, dimidy, dimidzRay, dimidt], rayvolid_z))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_z, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_z, 'long_name', 'ray &
            &volume position in z-direction'))
        call handle_err(nf90_put_att(rayvolid, rayvolid_z, 'units', 'm'))

        ! k wavenumber of ray volume
        call handle_err(nf90_def_var(rayvolid, 'ray_k', nf90_float, [dimidr, &
            &dimidx, dimidy, dimidzRay, dimidt], rayvolid_k))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_k, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_k, 'long_name', 'ray &
            &volume wavenumber k'))
        call handle_err(nf90_put_att(rayvolid, rayvolid_k, 'units', '1/m'))

        ! l wavenumber of ray volume
        call handle_err(nf90_def_var(rayvolid, 'ray_l', nf90_float, [dimidr, &
            &dimidx, dimidy, dimidzRay, dimidt], rayvolid_l))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_l, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_l, 'long_name', 'ray &
            &volume wavenumber l'))
        call handle_err(nf90_put_att(rayvolid, rayvolid_l, 'units', '1/m'))

        ! m wavenumber of ray volume
        call handle_err(nf90_def_var(rayvolid, 'ray_m', nf90_float, [dimidr, &
            &dimidx, dimidy, dimidzRay, dimidt], rayvolid_m))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_m, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_m, 'long_name', 'ray &
            &volume wavenumber m'))
        call handle_err(nf90_put_att(rayvolid, rayvolid_m, 'units', '1/m'))

        ! intrinsic frequency of ray volume
        call handle_err(nf90_def_var(rayvolid, 'ray_omega', nf90_float, &
            &[dimidr, dimidx, dimidy, dimidzRay, dimidt], rayvolid_omega))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_omega, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_omega, 'long_name', &
            &'ray volume intrinsic frequency'))
        call handle_err(nf90_put_att(rayvolid, rayvolid_omega, 'units', '1/s'))

        ! ray volume extent in k
        call handle_err(nf90_def_var(rayvolid, 'ray_dkray', nf90_float, &
            &[dimidr, dimidx, dimidy, dimidzRay, dimidt], rayvolid_dkray))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_dkray, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dkray, 'long_name', &
            &'ray volume extent in k'))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dkray, 'units', '1/m'))

        ! ray volume extent in l
        call handle_err(nf90_def_var(rayvolid, 'ray_dlray', nf90_float, &
            &[dimidr, dimidx, dimidy, dimidzRay, dimidt], rayvolid_dlray))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_dlray, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dlray, 'long_name', &
            &'ray volume extent in l'))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dlray, 'units', '1/m'))

        ! ray volume extent in m
        call handle_err(nf90_def_var(rayvolid, 'ray_dmray', nf90_float, &
            &[dimidr, dimidx, dimidy, dimidzRay, dimidt], rayvolid_dmray))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_dmray, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dmray, 'long_name', &
            &'ray volume extent in m'))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dmray, 'units', '1/m'))

        ! ray volume extent in x
        call handle_err(nf90_def_var(rayvolid, 'ray_dxray', nf90_float, &
            &[dimidr, dimidx, dimidy, dimidzRay, dimidt], rayvolid_dxray))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_dxray, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dxray, 'long_name', &
            &'ray volume extent in x'))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dxray, 'units', 'm'))

        ! ray volume extent in y
        call handle_err(nf90_def_var(rayvolid, 'ray_dyray', nf90_float, &
            &[dimidr, dimidx, dimidy, dimidzRay, dimidt], rayvolid_dyray))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_dyray, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dyray, 'long_name', &
            &'ray volume extent in y'))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dyray, 'units', 'm'))

        ! ray volume extent in z
        call handle_err(nf90_def_var(rayvolid, 'ray_dzray', nf90_float, &
            &[dimidr, dimidx, dimidy, dimidzRay, dimidt], rayvolid_dzray))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_dzray, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dzray, 'long_name', &
            &'ray volume extent in z'))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dzray, 'units', 'm'))

        ! ray volume wave-action density
        call handle_err(nf90_def_var(rayvolid, 'ray_dens', nf90_float, &
            &[dimidr, dimidx, dimidy, dimidzRay, dimidt], rayvolid_dens))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_dens, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dens, 'long_name', &
            &'ray volume phase-space wave-action density'))
        select case(dim)
        case(1)
          call handle_err(nf90_put_att(rayvolid, rayvolid_dens, 'units', &
              &'kg/s'))
        case(2)
          call handle_err(nf90_put_att(rayvolid, rayvolid_dens, 'units', &
              &'kg*m/s'))
        case(3)
          call handle_err(nf90_put_att(rayvolid, rayvolid_dens, 'units', &
              &'kg*m^2/s'))
        case default
          stop "Error in create_netCDF: dim > 3!"
        end select

        ! ray volume initial phase
        call handle_err(nf90_def_var(rayvolid, 'ray_dphi', nf90_float, &
            &[dimidr, dimidx, dimidy, dimidzRay, dimidt], rayvolid_dphi))
        call handle_err(nf90_var_par_access(rayvolid, rayvolid_dphi, &
            &nf90_collective))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dphi, 'long_name', &
            &'ray volume initial phase'))
        call handle_err(nf90_put_att(rayvolid, rayvolid_dphi, 'units', 'none'))
      end if
    end if

  end subroutine create_netCDF

  !-------------------------------------------------------------------------

  subroutine write_netCDF_background(iOutput)

    ! write the x-, y-, and z-axis to the netCDF file
    ! write background fields pStrat, thetaStrat, rhoStrat, bvsStrat
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
    if(topography) then
      call handle_err(nf90_put_var(ncid, ncidz, zTFC(1:nx, 1:ny, 1:sizeZ) &
          &* lRef, start = startNC, count = countNC))
    else
      if(iOutput == 0) then
        call handle_err(nf90_put_var(ncid, ncidz, z(1:sizeZ) * lRef, start &
            &= [1], count = [sizeZ]))
      end if
    end if

    ! background is constant in Boussinesq mode
    if(model == "Boussinesq") return

    ! mass-weighted potential temperature
    if(topography) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_p, pStratTFC(1:nx, 1:ny, &
          &1:sizeZ) * rhoRef * thetaRef, start = startNC, count = countNC), &
          &'save P')
    else
      call handle_err(nf90_put_var(atmvarid, atmvarid_p, pStrat(1:sizeZ) &
          &* rhoRef * thetaRef, start = startNC(3:4), count = countNC(3:4)), &
          &'save P')
    end if

    ! potential-temperature background
    if(topography) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_thetabar, &
          &thetaStratTFC(1:nx, 1:ny, 1:sizeZ) * thetaRef, start = startNC, &
          &count = countNC), 'save thetabar')
    else
      if(iOutput == 0) then
        call handle_err(nf90_put_var(atmvarid, atmvarid_thetabar, &
            &thetaStrat(1:sizeZ) * thetaRef, start = [1], count = [sizeZ]), &
            &'save thetabar')
      end if
    end if

    ! density background
    if(topography) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_rhobar, &
          &rhoStratTFC(1:nx, 1:ny, 1:sizeZ) * rhoRef, start = startNC, count &
          &= countNC), 'save rhobar')
    else
      if(iOutput == 0) then
        call handle_err(nf90_put_var(atmvarid, atmvarid_rhobar, &
            &rhoStrat(1:sizeZ) * rhoRef, start = [1], count = [sizeZ]), 'save &
            &rhobar')
      end if
    end if

    ! squared buoyancy frequency
    if(topography) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_n2, bvsStratTFC(1:nx, &
          &1:ny, 1:sizeZ) / tRef ** 2.0, start = startNC, count = countNC), &
          &'save N2')
    else
      if(iOutput == 0) then
        call handle_err(nf90_put_var(atmvarid, atmvarid_n2, bvsStrat(1:sizeZ) &
            &/ tRef ** 2.0, start = [1], count = [sizeZ]), 'save N2')
      end if
    end if

  end subroutine write_netCDF_background

  !-------------------------------------------------------------------------

  subroutine write_netCDF(iOutput, iTime, time, cpuTime, var, ray_var3D, &
      &tracerforce, waveAmplitudes, ray)

    implicit none

    integer, intent(inout) :: iOutput
    integer, intent(in) :: iTime
    real, intent(in) :: time
    real, intent(in) :: cpuTime

    type(var_type), intent(in) :: var

    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, 1:6), intent(in), optional &
        &:: ray_var3D
    type(tracerForceType), dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz &
        &+ nbz), intent(in), optional :: tracerforce
    type(waveAmpType), dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz &
        &+ nbz), intent(inout), optional :: waveAmplitudes

    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(in), optional :: ray

    real :: time_dim

    real :: cpuTimeLoc
    integer :: days, hours, minutes, seconds

    real, dimension(1:nx, 1:ny, 1:nz) :: rho, rhop
    real, dimension(1:nx, 1:ny, 1:nz) :: buoyancy
    real, dimension(1:nx, 1:ny, 1:nz) :: thetap
    real, dimension(1:nx, 1:ny, 0:nz) :: wOut

    ! for saving derived type variables and avoiding temporary arrays
    real, dimension(1:nx, 1:ny, 1:nz) :: tmparray

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

    if(model == "Boussinesq") then
      rhop = var%rhop(1:nx, 1:ny, 1:nz)
    else
      rhop = var%rho(1:nx, 1:ny, 1:nz)
    end if

    if(topography) then
      rho = rhop + rhoStratTFC(1:nx, 1:ny, 1:nz)
      buoyancy = g * (pStratTFC(1:nx, 1:ny, 1:nz) / rho - thetaStratTFC(1:nx, &
          &1:ny, 1:nz)) / thetaStratTFC(1:nx, 1:ny, 1:nz)
      thetap = buoyancy / g * thetaStratTFC(1:nx, 1:ny, 1:nz)
    else
      do kz = 1, nz
        rho(:, :, kz) = rhop(:, :, kz) + rhoStrat(kz)
        buoyancy(:, :, kz) = g * (pStrat(kz) / rho(:, :, kz) - thetaStrat(kz)) &
            &/ thetaStrat(kz)
        thetap(:, :, kz) = buoyancy(:, :, kz) / g * thetaStrat(kz)
      end do
    end if

    if(topography) then
      do ix = 1, nx
        do jy = 1, ny
          do kz = 0, nz
            wOut(ix, jy, kz) = vertWindTFC(ix, jy, kz, var)
          end do
        end do
      end do
    end if

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
    if(topography) then
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
        call handle_err(nf90_put_var(atmvarid, atmvarid_wTFC, (var%w(1:nx, &
            &1:ny, 0:nz - 1) + var%w(1:nx, 1:ny, 1:nz)) * uRef / 2., start &
            &= startNC, count = countNC), 'save wTFC')
      end if
      if(prepare_restart .or. any(atmvarOut == 'wsTFC')) then
        call handle_err(nf90_put_var(atmvarid, atmvarid_wsTFC, var%w(1:nx, &
            &1:ny, 1:nz) * uRef, start = startNC, count = countNC), 'save &
            &wsTFC')
      end if
    else
      if(any(atmvarOut == 'w')) then
        call handle_err(nf90_put_var(atmvarid, atmvarid_w, (var%w(1:nx, 1:ny, &
            &0:nz - 1) + var%w(1:nx, 1:ny, 1:nz)) * uRef / 2., start &
            &= startNC, count = countNC), 'save w')
      end if
      if(prepare_restart .or. any(atmvarOut == 'ws')) then
        call handle_err(nf90_put_var(atmvarid, atmvarid_ws, var%w(1:nx, 1:ny, &
            &1:nz) * uRef, start = startNC, count = countNC), 'save ws')
      end if
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

    ! tracer mixing ratio
    if(include_tracer .and. any(atmvarOut == 'tmr')) then
      call handle_err(nf90_put_var(atmvarid, atmvarid_tmr, var%chi(1:nx, 1:ny, &
          &1:nz) / rho, start = startNC, count = countNC))
    end if

    ! tracer mixing ratio difference
    if(include_tracer .and. (prepare_restart .or. any(atmvarOut == 'tmrd'))) &
        &then
      call handle_err(nf90_put_var(atmvarid, atmvarid_tmrd, var%chi(1:nx, &
          &1:ny, 1:nz) / rho - initialtracer(1:nx, 1:ny, 1:nz), start &
          &= startNC, count = countNC))
    end if

    ! dynamic Smagorinsky coefficient
    if(TurbScheme .and. any(atmvarOut == 'dsc')) then
      call handle_Err(nf90_put_var(atmvarid, atmvarid_dsc, var%DSC(1:nx, 1:ny, &
          &1:nz) * uRef * lRef, start = startNC, count = countNC))
    end if

    !--------------------------------
    ! save variables in rayvar group
    !--------------------------------

    ! gravity-wave heating
    if(rayTracer .and. any(rayvarOut == 'gwh')) then
      call handle_err(nf90_put_var(rayvarid, rayvarid_gwh, var%GWH(1:nx, 1:ny, &
          &1:nz) * rhoRef * thetaRef / tRef, start = startNC, count = countNC))
    end if

    if(present(ray_var3D)) then
      ! zonal-wind tendency
      if(any(rayvarOut == 'dudt')) then
        call handle_err(nf90_put_var(rayvarid, rayvarid_dudt, ray_var3D(1:nx, &
            &1:ny, 1:nz, 1) * uRef / tRef, start = startNC, count = countNC))
      end if

      ! meridional-wind tendency
      if(any(rayvarOut == 'dvdt')) then
        call handle_err(nf90_put_var(rayvarid, rayvarid_dvdt, ray_var3D(1:nx, &
            &1:ny, 1:nz, 2) * uRef / tRef, start = startNC, count = countNC))
      end if

      ! potential-temperature tendency
      if(any(rayvarOut == 'dthetadt')) then
        call handle_err(nf90_put_var(rayvarid, rayvarid_dthetadt, &
            &ray_var3D(1:nx, 1:ny, 1:nz, 3) * thetaRef / tRef, start &
            &= startNC, count = countNC))
      end if

      ! zonal-vertical momentum flux
      if(any(rayvarOut == 'uw')) then
        call handle_err(nf90_put_var(rayvarid, rayvarid_uw, ray_var3D(1:nx, &
            &1:ny, 1:nz, 4) * uRef ** 2, start = startNC, count = countNC))
      end if

      ! meridional-vertical momentum flux
      if(any(rayvarOut == 'vw')) then
        call handle_err(nf90_put_var(rayvarid, rayvarid_vw, ray_var3D(1:nx, &
            &1:ny, 1:nz, 5) * uRef ** 2, start = startNC, count = countNC))
      end if

      ! gravity-wave energy density
      if(any(rayvarOut == 'E')) then
        call handle_err(nf90_put_var(rayvarid, rayvarid_e, ray_var3D(1:nx, &
            &1:ny, 1:nz, 6) * uRef ** 2, start = startNC, count = countNC))
      end if
    end if

    if(present(waveAmplitudes)) then
      ! large-scale wave phase
      if(any(rayvarOut == 'lsphase')) then
        tmparray = waveAmplitudes(1:nx, 1:ny, 1:nz)%phase
        call handle_err(nf90_put_var(rayvarid, rayvarid_lsphase, tmparray, &
            &start = startNC, count = countNC))
      end if

      ! leading-order buoyancy amplitude
      if(any(rayvarOut == 'lobhat')) then
        call handle_err(nf90_put_var(rayvarid, rayvarid_lobhat, &
            &abs(waveAmplitudes(1:nx, 1:ny, 1:nz)%lowamp%b * g), start &
            &= startNC, count = countNC))
      end if

      ! next-order buoyancy amplitude
      if(any(rayvarOut == 'nobhat')) then
        call handle_err(nf90_put_var(rayvarid, rayvarid_nobhat, &
            &abs(waveAmplitudes(1:nx, 1:ny, 1:nz)%nowamp%b * g), start &
            &= startNC, count = countNC))
      end if
    end if

    if(present(tracerforce)) then
      if(any(rayvarOut == 'tfrcnow')) then
        tmparray = tracerforce(1:nx, 1:ny, 1:nz)%noforce%wflx * uRef
        call handle_err(nf90_put_var(rayvarid, rayvarid_tfrcnow, tmparray, &
            &start = startNC, count = countNC))
      end if
      if(any(rayvarOut == 'tfrclou')) then
        tmparray = tracerforce(1:nx, 1:ny, 1:nz)%loforce%uflx * uRef
        call handle_err(nf90_put_var(rayvarid, rayvarid_tfrclou, tmparray, &
            &start = startNC, count = countNC))
      end if
      if(any(rayvarOut == 'tfrclot')) then
        tmparray = tracerforce(1:nx, 1:ny, 1:nz)%loforce%total / tRef
        call handle_err(nf90_put_var(rayvarid, rayvarid_tfrclot, tmparray, &
            &start = startNC, count = countNC))
      end if
      if(any(rayvarOut == 'tfrcnot')) then
        tmparray = tracerforce(1:nx, 1:ny, 1:nz)%noforce%total / tRef
        call handle_err(nf90_put_var(rayvarid, rayvarid_tfrcnot, tmparray, &
            &start = startNC, count = countNC))
      end if
    end if

    if(rayTracer .and. present(ray) .and. (prepare_restart .or. saverayvols)) &
        &then
      call write_rayvolumes(iOutput, ray)
    end if

    if(topography .or. (rayTracer .and. case_wkb == 3)) then
      call write_topography_netCDF(iOutput)
    end if

    iOutput = iOutput + 1

  end subroutine write_netCDF

  !-------------------------------------------------------------------------

  subroutine write_rayvolumes(iOutput, ray)

    implicit none

    integer, intent(in) :: iOutput

    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(in) :: ray

    real, dimension(nray_max, 1:nx, 1:ny, 0:nz + 1) :: tmparray

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

    startNCray = (/1, startxNC, startyNC, 1, iOutput + 1/)
    countNCray = (/nray_max, countxNC, countyNC, nz + 2, 1/)

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%x
    call handle_err(nf90_put_var(rayvolid, rayvolid_x, tmparray * lRef, start &
        &= startNCray, count = countNCray))

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%y
    call handle_err(nf90_put_var(rayvolid, rayvolid_y, tmparray * lRef, start &
        &= startNCray, count = countNCray))

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%z
    call handle_err(nf90_put_var(rayvolid, rayvolid_z, tmparray * lRef, start &
        &= startNCray, count = countNCray))

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%k
    call handle_err(nf90_put_var(rayvolid, rayvolid_k, tmparray / lRef, start &
        &= startNCray, count = countNCray))

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%l
    call handle_err(nf90_put_var(rayvolid, rayvolid_l, tmparray / lRef, start &
        &= startNCray, count = countNCray))

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%m
    call handle_err(nf90_put_var(rayvolid, rayvolid_m, tmparray / lRef, start &
        &= startNCray, count = countNCray))

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%omega
    call handle_err(nf90_put_var(rayvolid, rayvolid_omega, tmparray / tRef, &
        &start = startNCray, count = countNCray))

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dkray
    call handle_err(nf90_put_var(rayvolid, rayvolid_dkray, tmparray / lRef, &
        &start = startNCray, count = countNCray))

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dlray
    call handle_err(nf90_put_var(rayvolid, rayvolid_dlray, tmparray / lRef, &
        &start = startNCray, count = countNCray))

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dmray
    call handle_err(nf90_put_var(rayvolid, rayvolid_dmray, tmparray / lRef, &
        &start = startNCray, count = countNCray))

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dxray
    call handle_err(nf90_put_var(rayvolid, rayvolid_dxray, tmparray * lRef, &
        &start = startNCray, count = countNCray))

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dyray
    call handle_err(nf90_put_var(rayvolid, rayvolid_dyray, tmparray * lRef, &
        &start = startNCray, count = countNCray))

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dzray
    call handle_err(nf90_put_var(rayvolid, rayvolid_dzray, tmparray * lRef, &
        &start = startNCray, count = countNCray))

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dens
    call handle_err(nf90_put_var(rayvolid, rayvolid_dens, tmparray * rhoRef &
        &* uRef ** 2 * tRef * lRef ** dim, start = startNCray, count &
        &= countNCray))

    tmparray = ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dphi
    call handle_err(nf90_put_var(rayvolid, rayvolid_dphi, tmparray, start &
        &= startNCray, count = countNCray))

  end subroutine write_rayvolumes

  !-------------------------------------------------------------------------

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

    starthNC = 1
    counthNC = nwm

    starttNC = iOutput + 1
    counttNC = 1

    ! Write large-scale topography.
    if(topography) then
      call handle_err(nf90_put_var(ncid, ncidhm, topography_surface(1:nx, &
          &1:ny) * lRef, start = [startxNC, startyNC, starttNC], count &
          &= [countxNC, countyNC, counttNC]))
    end if

    ! Write small-scale topography.
    if(rayTracer .and. case_wkb == 3) then
      ! Write zonal wavenumbers.
      call handle_err(nf90_put_var(ncid, ncidkh, k_spectrum(1:nx, 1:ny, 1:nwm) &
          &/ lRef, start = [startxNC, startyNC, starthNC, starttNC], count &
          &= [countxNC, countyNC, counthNC, counttNC]))

      ! Write meridional wavenumbers.
      call handle_err(nf90_put_var(ncid, ncidlh, l_spectrum(1:nx, 1:ny, 1:nwm) &
          &/ lRef, start = [startxNC, startyNC, starthNC, starttNC], count &
          &= [countxNC, countyNC, counthNC, counttNC]))

      ! Write wave amplitudes.
      call handle_err(nf90_put_var(ncid, ncidhw, topography_spectrum(1:nx, &
          &1:ny, 1:nwm) * lRef, start = [startxNC, startyNC, starthNC, &
          &starttNC], count = [countxNC, countyNC, counthNC, counttNC]))
    end if

  end subroutine write_topography_netCDF

  !-------------------------------------------------------------------------

  subroutine close_netCDF

    implicit none

    call handle_err(nf90_close(ncid))

  end subroutine close_netCDF

  !-------------------------------------------------------------------------

  subroutine read_netCDF(timeStart, var, ray, time)

    implicit none

    integer, intent(in) :: timeStart

    type(var_type), intent(inout) :: var

    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(inout), optional :: ray

    real, dimension(1:nx, 1:ny, 1:nz) :: rhop

    real, dimension(nray_max, 1:nx, 1:ny, 0:nz + 1) :: tmparray

    real, intent(out), optional :: time

    integer :: Ntimesteps

    dim = 1
    if(sizeX > 1) dim = dim + 1
    if(sizeY > 1) dim = dim + 1

    ! Open file.
    call handle_err(nf90_open(filenameIn, ior(nf90_nowrite, nf90_mpiio), ncid, &
        &comm = comm, info = mpi_info_null))

    ! Get group IDs.
    call handle_err(nf90_inq_ncid(ncid, "atmvar", atmvarid))
    if(rayTracer) then
      call handle_err(nf90_inq_ncid(ncid, "rayvol", rayvolid), 'read rayvolid')
    end if

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
    if(model == "Boussinesq") then
      var%rhop(1:nx, 1:ny, 1:nz) = rhop(1:nx, 1:ny, 1:nz) / rhoRef
    else
      var%rho(1:nx, 1:ny, 1:nz) = rhop(1:nx, 1:ny, 1:nz) / rhoRef
    end if

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
    if(topography) then
      call handle_err(nf90_inq_varid(atmvarid, 'wsTFC', atmvarid_wsTFC))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_wsTFC, &
          &nf90_collective))
      call handle_err(nf90_get_var(atmvarid, atmvarid_wsTFC, var%w(1:nx, 1:ny, &
          &1:nz), start = startNC, count = countNC))
    else
      call handle_err(nf90_inq_varid(atmvarid, 'ws', atmvarid_ws))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_ws, &
          &nf90_collective))
      call handle_err(nf90_get_var(atmvarid, atmvarid_ws, var%w(1:nx, 1:ny, &
          &1:nz), start = startNC, count = countNC))
    end if
    var%w(1:nx, 1:ny, 1:nz) = var%w(1:nx, 1:ny, 1:nz) / uRef

    ! Read Exner-pressure fluctuations.
    call handle_err(nf90_inq_varid(atmvarid, 'pip', atmvarid_pip))
    call handle_err(nf90_var_par_access(atmvarid, atmvarid_pip, &
        &nf90_collective))
    call handle_err(nf90_get_var(atmvarid, atmvarid_pip, var%pi(1:nx, 1:ny, &
        &1:nz), start = startNC, count = countNC))

    ! Read tracer-mixing-ratio difference.
    if(include_tracer) then
      call handle_err(nf90_inq_varid(atmvarid, 'tmrd', atmvarid_tmrd))
      call handle_err(nf90_var_par_access(atmvarid, atmvarid_tmrd, &
          &nf90_collective))
      call handle_err(nf90_get_var(atmvarid, atmvarid_tmrd, var%chi(1:nx, &
          &1:ny, 1:nz), start = startNC, count = countNC))
      var%chi(1:nx, 1:ny, 1:nz) = var%chi(1:nx, 1:ny, 1:nz) &
          &+ initialtracer(1:nx, 1:ny, 1:nz)
    end if

    ! Read mass-weighted potential temperature.
    call handle_err(nf90_inq_varid(atmvarid, 'P', atmvarid_p))
    call handle_err(nf90_var_par_access(atmvarid, atmvarid_p, nf90_collective))
    if(model == 'compressible') then
      call handle_err(nf90_get_var(atmvarid, atmvarid_p, var%P(1:nx, 1:ny, &
          &1:nz), start = startNC, count = countNC))
      var%P(1:nx, 1:ny, 1:nz) = var%P(1:nx, 1:ny, 1:nz) / rhoRef / thetaRef
    else
      call handle_err(nf90_get_var(atmvarid, atmvarid_p, pStrat(1:nz), start &
          &= startNC(3:4), count = countNC(3:4)))
      pStrat(1:nz) = pStrat(1:nz) / rhoRef / thetaRef
    end if

    if(rayTracer) then

      if(timeStart >= 0) then
        startNCray = (/1, startxNC, startyNC, 1, timeStart + 1/)
      else
        startNCray = (/1, startxNC, startyNC, 1, Ntimesteps + timeStart + 1/)
      end if
      countNCray = (/nray_max, countxNC, countyNC, nz + 2, 1/)

      call handle_err(nf90_inq_varid(rayvolid, 'ray_x', rayvolid_x))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_x, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_x, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_x')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%x = tmparray / lRef

      call handle_err(nf90_inq_varid(rayvolid, 'ray_y', rayvolid_y))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_y, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_y, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_y')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%y = tmparray / lRef

      call handle_err(nf90_inq_varid(rayvolid, 'ray_z', rayvolid_z))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_z, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_z, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_z')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%z = tmparray / lRef

      call handle_err(nf90_inq_varid(rayvolid, 'ray_k', rayvolid_k))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_k, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_k, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_k')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%k = tmparray * lRef

      call handle_err(nf90_inq_varid(rayvolid, 'ray_l', rayvolid_l))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_l, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_l, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_l')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%l = tmparray * lRef

      call handle_err(nf90_inq_varid(rayvolid, 'ray_m', rayvolid_m))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_m, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_m, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_m')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%m = tmparray * lRef

      call handle_err(nf90_inq_varid(rayvolid, 'ray_omega', rayvolid_omega))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_omega, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_omega, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_omega')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%omega = tmparray * tRef

      call handle_err(nf90_inq_varid(rayvolid, 'ray_dkray', rayvolid_dkray))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_dkray, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_dkray, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_dkray')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dkray = tmparray * lRef

      call handle_err(nf90_inq_varid(rayvolid, 'ray_dlray', rayvolid_dlray))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_dlray, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_dlray, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_dlray')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dlray = tmparray * lRef

      call handle_err(nf90_inq_varid(rayvolid, 'ray_dmray', rayvolid_dmray))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_dmray, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_dmray, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_dmray')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dmray = tmparray * lRef

      call handle_err(nf90_inq_varid(rayvolid, 'ray_dxray', rayvolid_dxray))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_dxray, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_dxray, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_dxray')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dxray = tmparray / lRef

      call handle_err(nf90_inq_varid(rayvolid, 'ray_dyray', rayvolid_dyray))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_dyray, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_dyray, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_dyray')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dyray = tmparray / lRef

      call handle_err(nf90_inq_varid(rayvolid, 'ray_dzray', rayvolid_dzray))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_dzray, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_dzray, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_dzray')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dzray = tmparray * lRef

      call handle_err(nf90_inq_varid(rayvolid, 'ray_dens', rayvolid_dens))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_dens, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_dens, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_dens')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dens = tmparray / rhoRef / uRef &
          &** 2 / tRef / lRef ** dim

      call handle_err(nf90_inq_varid(rayvolid, 'ray_dphi', rayvolid_dphi))
      call handle_err(nf90_var_par_access(rayvolid, rayvolid_dphi, &
          &nf90_collective))
      call handle_err(nf90_get_var(rayvolid, rayvolid_dphi, tmparray, start &
          &= startNCray, count = countNCray), 'reading ray_dphi')
      ray(1:nray_max, 1:nx, 1:ny, 0:nz + 1)%dphi = tmparray
    end if

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

    starthNC = 1
    counthNC = nwm

    if(timeStart >= 0) then
      starttNC = timeStart + 1
      counttNC = 1
    else
      starttNC = nt + timeStart + 1
      counttNC = 1
    endif

    ! Read large-scale topography (non-dimensionalization is done in
    ! setup_topography).
    if(topography) then
      call handle_err(nf90_inq_varid(ncid, 'hm', ncidhm))
      call handle_err(nf90_var_par_access(ncid, ncidhm, nf90_collective))
      call handle_err(nf90_get_var(ncid, ncidhm, topography_surface(1:nx, &
          &1:ny), start = [startxNC, startyNC, starttNC], count = [countxNC, &
          &countyNC, counttNC]))
    end if

    ! Read small-scale topography (non-dimensionalization is done in
    ! setup_topography).
    if(rayTracer .and. case_wkb == 3) then
      ! Read zonal wavenumbers.
      call handle_err(nf90_inq_varid(ncid, 'kh', ncidkh))
      call handle_err(nf90_var_par_access(ncid, ncidkh, nf90_collective))
      call handle_err(nf90_get_var(ncid, ncidkh, k_spectrum(1:nx, 1:ny, &
          &1:nwm), start = [startxNC, startyNC, starthNC, starttNC], count &
          &= [countxNC, countyNC, counthNC, counttNC]))

      ! Read meridional wavenumbers.
      call handle_err(nf90_inq_varid(ncid, 'lh', ncidlh))
      call handle_err(nf90_var_par_access(ncid, ncidlh, nf90_collective))
      call handle_err(nf90_get_var(ncid, ncidlh, l_spectrum(1:nx, 1:ny, &
          &1:nwm), start = [startxNC, startyNC, starthNC, starttNC], count &
          &= [countxNC, countyNC, counthNC, counttNC]))

      ! Read wave amplitudes.
      call handle_err(nf90_inq_varid(ncid, 'hw', ncidhw))
      call handle_err(nf90_var_par_access(ncid, ncidhw, nf90_collective))
      call handle_err(nf90_get_var(ncid, ncidhw, topography_spectrum(1:nx, &
          &1:ny, 1:nwm), start = [startxNC, startyNC, starthNC, starttNC], &
          &count = [countxNC, countyNC, counthNC, counttNC]))
    end if

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
