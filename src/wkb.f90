module wkb_module

  !-----------------------------------------------------------
  ! This module couples the WKB model (MS-GWaM) to PincFloit
  !-----------------------------------------------------------

  use type_module
  use timeScheme_module
  use atmosphere_module
  use muscl_module
  use ice_module, ONLY:Psat_ice
  use mpi

  implicit none

  private ! all module variables are internal to the module

  !----------------------
  !   public routines
  !----------------------
  public :: setup_wkb
  public :: transport_rayvol
  public :: split_rayvol
  public :: shift_rayvol
  public :: merge_rayvol
  public :: boundary_rayvol
  public :: calc_meanFlow_effect
  public :: calc_tracerforce

  public :: saturation_3D

  public :: meanflow
  public :: stratification

  public :: smooth_wkb_box
  public :: smooth_wkb_shapiro

  public :: setboundary_wkb
  public :: setboundary_hor_wkb
  public :: setboundary_vrt_wkb
  public :: setboundary_x_periodic_wkb
  public :: setboundary_y_periodic_wkb
  public :: setboundary_z_periodic_wkb
  public :: setboundary_z_solidwall_wkb

  public :: setboundary_frc_wkb
  public :: setboundary_frc_hor_wkb
  public :: setboundary_frc_vrt_wkb
  public :: setboundary_frc_x_periodic_wkb
  public :: setboundary_frc_y_periodic_wkb
  public :: setboundary_frc_z_periodic_wkb
  public :: setboundary_frc_z_solidwall_wkb
  public :: setBoundary_waveAmp

  public :: setboundary_wkb_cmplx
  public :: setboundary_hor_wkb_cmplx
  public :: setboundary_vrt_wkb_cmplx
  public :: setboundary_x_periodic_wkb_cmplx
  public :: setboundary_y_periodic_wkb_cmplx
  public :: setboundary_z_periodic_wkb_cmplx
  public :: setboundary_z_solidwall_wkb_cmplx

  public :: calc_ice

  !------------------------------
  !   private module variables
  !------------------------------

  real, dimension(:, :, :, :, :), allocatable :: dxRay, dkRay
  real, dimension(:, :, :, :, :), allocatable :: ddxRay

  integer, dimension(:), allocatable :: ix2_sfc, jy2_sfc, kz2_sfc, ik_sfc, &
      &jl_sfc, km_sfc
  integer, dimension(:, :, :), allocatable :: ir_sfc

  ! FJApr2023
  integer, dimension(:), allocatable :: iwm_sfc

  real, dimension(:, :, :, :), allocatable :: dpRay

  integer :: iRay ! index of ray v. within a
  ! cell
  integer :: ixrv, jyrv, kzrv ! cell indices of ray volume
  integer :: nxRay_wrk, nyRay_wrk, nzRay_wrk ! nb. of rays along x,y,z in
  ! work space
  integer :: nxRay, nyRay, nzRay ! max. nb. of rays along x,y,z
  ! (beyond which ray volumes
  ! are merged)
  integer :: i_sfc, n_sfc ! indices for surface ray v.

  logical, parameter :: debugging = .false.

  real, dimension(:, :, :), allocatable :: zTFC, zTildeTFC

  contains

  !-----------------------------------------------------------------------
  subroutine calc_tracerforce(ray, var, ray_var3D, tracerforce, &
      &waveAmplitudes, dt)

    ! calculate the GW tracer fluxes

    implicit none

    ! in/out variables
    type(var_type), intent(in) :: var

    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(in) :: ray

    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, 1:6), intent(in) :: ray_var3D

    real, intent(in) :: dt

    type(tracerForceType), dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz &
        &+ nbz), intent(inout) :: tracerforce

    ! leading-order and next-order wave amplitudes lowamp and nowamp, resp.
    ! rhs to calculate nowamp in rhsamp
    type(waveAmpType), dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz &
        &+ nbz), intent(inout) :: waveAmplitudes

    ! leading-order fluxes
    complex, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1) :: louchi, lovchi, lowchi
    real :: rho

    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1) :: Kd
    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1) :: chiLS

    complex, dimension(5) :: rhsvector
    complex, dimension(5) :: novector

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: &
        &waveactiondensity

    ! tracer fluxes
    real, allocatable :: var_utracer(:, :, :) ! zonal
    real, allocatable :: var_vtracer(:, :, :) ! meridional
    real, allocatable :: var_wtracer(:, :, :) ! vertical

    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1) :: tmparray
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: tmparrayL
    complex, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1) :: tmparray_cmplx

    integer :: ix, jy, kz, im
    real :: dxi, dyi, dzi
    integer :: ixmin, ixmax
    integer :: jymin, jymax
    integer :: kzmin, kzmax
    integer :: ixr, jyr, kzr
    integer :: ix0, jy0

    real :: rhotot

    real :: wnrk
    real :: wnrl
    real :: wnrm
    real :: wnrh
    real :: dwnrk, dwnrl, dwnrm
    real :: f_cor_nd

    real :: cgirx, cgiry, cgirz
    real :: omir

    real :: xr, yr, zr
    real :: dxr, dyr, dzr

    real :: fcpspx, fcpspy, fcpspz
    real :: wadr
    real :: omegaMax, kMax, lMax, mMax
    real :: phi

    real :: NNR
    logical :: apply

    real :: tracerfluxcoeff, dchidx, dchidy, dchidz, rhotracerp, rhotracerm, &
        &dutracer, dvtracer, rhotracern

    real :: ff, kk, ll, mm, omega, kh
    real, dimension(3) :: cgroup ! x-, y-, z-component of group velocity
    complex :: uhattilde, vhattilde, whattilde, pihattilde

    ! allocate tracer fluxes
    allocate(var_utracer(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))
    allocate(var_vtracer(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))
    allocate(var_wtracer(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))

    var_utracer = 0.
    var_vtracer = 0.
    var_wtracer = 0.

    if(topography) then
      stop "wkb with TCF and tracer not implemented yet"
    end if

    ! declare wavepacket variables (monochromatic WP)
    ff = f_Coriolis_dim * tRef ! non dimensional Coriolis param.
    f_cor_nd = f_Coriolis_dim * tRef
    if(wlrx_init == 0.) then
      kk = 0.
    else
      kk = 2. * pi / (wlrx_init / lRef)
    end if
    if(wlry_init == 0.) then
      ll = 0.
    else
      ll = 2. * pi / (wlry_init / lRef)
    end if
    mm = 2. * pi / (wlrz_init / lRef)
    kh = sqrt(kk ** 2. + ll ** 2.)
    omega = branchr * sqrt((NN ** 2. * kh ** 2. + ff ** 2. * mm ** 2.) / (kh &
        &** 2. + mm ** 2.))
    cgroup(1) = NN ** 2. * kk * mm ** 2. / omega / (kh ** 2. + mm ** 2.) ** 2.
    cgroup(2) = NN ** 2. * ll * mm ** 2. / omega / (kh ** 2. + mm ** 2.) ** 2.
    cgroup(3) = - NN ** 2. * kh ** 2. * mm / omega / (kh ** 2. + mm ** 2.) ** 2.

    ! currently only for no rotation
    uhattilde = cmplx(0., (omega ** 2. - NN ** 2.) / (mm * NN ** 2. * omega &
        &** 2.) * kk * omega)
    vhattilde = cmplx(0., (omega ** 2. - NN ** 2.) / (mm * NN ** 2. * omega &
        &** 2.) * ll * omega)
    whattilde = cmplx(0., omega / NN ** 2.)
    pihattilde = cmplx(0., (omega ** 2. - NN ** 2.) / (mm * NN ** 2.))

    waveactiondensity = 0.

    chiLS = var%chi(0:nx + 1, 0:ny + 1, 0:nz + 1)

    omegaMax = 0.
    kMax = 0.
    lMax = 0.
    mMax = 0.

    do kzrv = 0, nz + 1
      do jyrv = 0, ny + 1
        ! loop including ghost cells in order to get all fluxes
        ! affecting a cell
        ! (assuming that ray volumes are not wider in y than dy)

        do ixrv = 0, nx + 1
          ! loop including ghost cells in order to get all fluxes
          ! affecting a cell
          ! (assuming that ray volumes are not wider in x than dx)

          if(nRay(ixrv, jyrv, kzrv) < 1) cycle

          do iRay = 1, nRay(ixrv, jyrv, kzrv)
            ! skip counting ray volumes with zero wave-action density

            if(ray(iRay, ixrv, jyrv, kzrv)%dens == 0.0) cycle

            xr = ray(iRay, ixrv, jyrv, kzrv)%x
            yr = ray(iRay, ixrv, jyrv, kzrv)%y
            zr = ray(iRay, ixrv, jyrv, kzrv)%z

            dxr = ray(iRay, ixrv, jyrv, kzrv)%dxray
            dyr = ray(iRay, ixrv, jyrv, kzrv)%dyray
            dzr = ray(iRay, ixrv, jyrv, kzrv)%dzray

            ! vertical boundary conditions:
            ! zBoundary = 'periodic': implement periodicity
            ! zBoundary = 'solid_wall': skip counting ray volumes
            !                           that have completely left the
            !                           model domain

            apply = .false.
            if(.not. topography .and. zr < lz(0)) then
              apply = .true.
            else if(topography) then
              if(zr < topography_surface(ixrv, jyrv)) then
                apply = .true.
              end if
            end if
            if(apply) then
              select case(zBoundary)
              case("periodic")
                zr = lz(1) + mod(zr - lz(0), lz(1) - lz(0))
              case("solid_wall")
                if((.not. topography .and. zr + 0.5 * dzr < lz(0)) .or. &
                    &(topography .and. zr + 0.5 * dzr &
                    &< topography_surface(ixrv, jyrv))) cycle
              case default
                stop "calc_meanflow_effect: unknown case zBoundary"
              end select
            elseif(zr > lz(1)) then
              select case(zBoundary)
              case("periodic")
                zr = lz(0) + mod(zr - lz(1), lz(1) - lz(0))
              case("solid_wall")
                if(zr - 0.5 * dzr > lz(1)) cycle
              case default
                stop "calc_meanflow_effect: unknown case zBoundary"
              end select
            end if

            ! implement horizontal boundary conditions for ray-volume
            ! positions

            if(sizeX > 1) then
              if(xBoundary /= "periodic") then
                print *, 'ERROR in calc_meanflow_effect:  boundary conditions &
                    &in x must be periodic'
                stop
              end if

              ! for leftmost cpu make sure that xr in
              ! ghost cell to the left is between x(0) - dx
              ! and x(0)

              if(ixrv == 0 .and. is + nbx == 1) then
                if(xr > lx(1) - dx .and. xr < lx(1)) then
                  xr = xr - lx(1) + lx(0)
                elseif(xr > lx(0) - dx .and. xr < lx(0)) then
                  xr = xr
                else
                  print *, 'ERROR in calc_meanflow_effect:'
                  print *, 'xr =', xr
                  print *, 'but r.v. is in ghost cell to the  left of the &
                      &leftmost cpu so that  one should have either'
                  print *, lx(1) - dx, '= lx(1) - dx < xr < lx(1) =', lx(1), ' &
                      &or'
                  print *, lx(0) - dx, '= lx(0) - dx < xr < lx(0) =', lx(0)
                  stop
                end if
              end if

              ! for rightmost cpu make sure that xr in
              ! ghost cell to the right is between x(0) + L_x
              ! and x(0) + L_x + dx

              if(ixrv == nx + 1 .and. is + nbx + nx == sizeX) then
                if(xr > lx(0) .and. xr < lx(0) + dx) then
                  xr = xr + lx(1) - lx(0)
                elseif(xr > lx(1) .and. xr < lx(1) + dx) then
                  xr = xr
                else
                  print *, 'ERROR in calc_meanflow_effect:'
                  print *, 'xr =', xr
                  print *, 'but r.v. is in ghost cell to the  right of the &
                      &rightmost cpu so that  one should have either'
                  print *, lx(0), '= lx(0) < xr < lx(0) + dx =', lx(0) + dx, ' &
                      &or'
                  print *, lx(1), '= lx(1) < xr < lx(1) + dx =', lx(1) + dx
                  stop
                end if
              end if
            end if

            if(sizeY > 1) then
              if(yBoundary /= "periodic") then
                print *, 'ERROR in calc_meanflow_effect:  boundary conditions &
                    &in y must be periodic'
                stop
              end if

              ! for first cpu in y direct. make sure that yr in
              ! ghost cell in front is between y(0) - dy
              ! and y(0)

              if(jyrv == 0 .and. js + nby == 1) then
                if(yr > ly(1) - dy .and. yr < ly(1)) then
                  yr = yr - ly(1) + ly(0)
                elseif(yr > ly(0) - dy .and. yr < ly(0)) then
                  yr = yr
                else
                  print *, 'ERROR in calc_meanflow_effect:'
                  print *, 'yr =', yr
                  print *, 'but r.v. is in ghost cell in front of the first &
                      &cpu in y dir. so that  one should have either'
                  print *, ly(1) - dy, '= ly(1) - dy < yr < ly(1) =', ly(1), ' &
                      &or'
                  print *, ly(0) - dy, '= ly(0) - dy < yr < ly(0) =', ly(0)
                  stop
                end if
              end if

              ! for last cpu in y direction make sure that yr
              ! in ghost cell behind is between
              ! y(0) + L_y and y(0) + L_y + dy

              if(jyrv == ny + 1 .and. js + nby + ny == sizeY) then
                if(yr > ly(0) .and. yr < ly(0) + dy) then
                  yr = yr + ly(1) - ly(0)
                elseif(yr > ly(1) .and. yr < ly(1) + dy) then
                  yr = yr
                else
                  print *, 'ERROR in calc_meanflow_effect:'
                  print *, 'yr =', yr
                  print *, 'but r.v. is in ghost cell behind the last cpu in y &
                      &dir. so that  one should have either'
                  print *, ly(0), '= ly(0) < yr < ly(0) + dy =', ly(0) + dy, ' &
                      &or'
                  print *, ly(1), '= ly(1) < yr < ly(1) + dy =', ly(1) + dy
                  stop
                end if
              end if
            end if

            wnrk = ray(iRay, ixrv, jyrv, kzrv)%k
            wnrl = ray(iRay, ixrv, jyrv, kzrv)%l
            wnrm = ray(iRay, ixrv, jyrv, kzrv)%m

            dwnrk = ray(iRay, ixrv, jyrv, kzrv)%dkray
            dwnrl = ray(iRay, ixrv, jyrv, kzrv)%dlray
            dwnrm = ray(iRay, ixrv, jyrv, kzrv)%dmray

            wnrh = sqrt(wnrk ** 2 + wnrl ** 2)

            apply = .false.
            if(.not. topography .and. zr < lz(0) - dz) then
              apply = .true.
            else if(topography) then
              if(zr < topography_surface(ixrv, jyrv) - jac(ixrv, jyrv, 0) &
                  &* dz) then
                apply = .true.
              end if
            end if
            if(apply) then
              print *, 'ERROR IN calc_meanflow_effect: RAY VOLUME', iRay, 'in &
                  &cell', ixrv, jyrv, kzrv, 'TOO LOW'
              stop
            end if

            call stratification(zr, 1, NNr)

            omir = branchr * sqrt(NNr * wnrh ** 2 + ff ** 2 * wnrm ** 2) &
                &/ sqrt(wnrh ** 2 + wnrm ** 2)

            if(abs(omir) > abs(omegaMax)) then
              omegaMax = omir
              kMax = wnrk
              lMax = wnrl
              mMax = wnrm
            end if

            cgirx = wnrk * (NNr - omir ** 2) / (omir * (wnrh ** 2 + wnrm ** 2))
            cgiry = wnrl * (NNr - omir ** 2) / (omir * (wnrh ** 2 + wnrm ** 2))

            cgirz = - wnrm * (omir ** 2 - ff ** 2) / (omir * (wnrh ** 2 + wnrm &
                &** 2))

            ! indices of range of cells touched by a ray volume

            if(sizeX > 1) then
              ! last x-index leftmost of cpu
              ix0 = is + nbx - 1

              ixmin = floor((xr - dxr * 0.5 - lx(0)) / dx) + 1 - ix0
              ixmax = floor((xr + dxr * 0.5 - lx(0)) / dx) + 1 - ix0

              if(ixmin > nx + 1) then
                print *, 'ixmin =', ixmin, '> nx+1 = ', nx + 1
                print *, 'ixrv = ', ixrv
                print *, 'xr =', xr
                print *, 'dxr =', dxr
                print *, 'dx =', dx
                print *, 'lx(0) =', lx(0)
                print *, 'lx(1) =', lx(1)
                print *, 'ix0 =', ix0
                print *, 'floor((xr - dxr*0.5 - lx(0)) / dx) + 1 =', floor((xr &
                    &- dxr * 0.5 - lx(0)) / dx) + 1
                stop
              else
                ! no fluxes calculated for the ghost cells
                ! (that are taken care of by the boundary-condition
                ! routines)

                ixmin = max(ixmin, 1)
              end if

              if(ixmax < 0) then
                print *, 'ixmax =', ixmax, '< 0'
                stop
              else
                ! no fluxes calculated for the ghost cells
                ! (that are taken care of by the boundary-condition
                ! routines)

                ixmax = min(ixmax, nx)
              end if
            else
              ixmin = 1
              ixmax = 1
            end if

            if(sizeY > 1) then
              ! last y-index in front of cpu
              jy0 = js + nby - 1

              jymin = floor((yr - dyr * 0.5 - ly(0)) / dy) + 1 - jy0
              jymax = floor((yr + dyr * 0.5 - ly(0)) / dy) + 1 - jy0

              if(jymin > ny + 1) then
                print *, 'jymin =', jymin, '> ny+1 = ', ny + 1
                stop
              else
                ! no fluxes calculated for the ghost cells
                ! (that are taken care of by the boundary-condition
                ! routines)

                jymin = max(jymin, 1)
              end if

              if(jymax < 0) then
                print *, 'jymax =', jymax, '< 0'
                stop
              else
                ! no fluxes calculated for the ghost cells
                ! (that are taken care of by the boundary-condition
                ! routines)

                jymax = min(jymax, ny)
              end if
            else
              jymin = 1
              jymax = 1
            end if

            kzmin = max(1, floor((zr - dzr * 0.5 - lz(0)) / dz) + 1)
            kzmax = min(nz, floor((zr + dzr * 0.5 - lz(0)) / dz) + 1)

            ! calculate momentum-flux / energy / elastic-term
            ! contribution from each ray volume

            do kz = kzmin, kzmax
              dzi = (min((zr + dzr * 0.5), lz(0) + kz * dz) - max((zr - dzr &
                  &* 0.5), lz(0) + (kz - 1) * dz))

              fcpspz = dwnrm * dzi / dz

              do jy = jymin, jymax
                if(sizeY > 1) then
                  dyi = (min((yr + dyr * 0.5), ly(0) + (jy + jy0) * dy) &
                      &- max((yr - dyr * 0.5), ly(0) + (jy + jy0 - 1) * dy))

                  fcpspy = dwnrl * dyi / dy
                else
                  fcpspy = 1.0
                end if

                do ix = ixmin, ixmax
                  if(sizeX > 1) then
                    dxi = (min((xr + dxr * 0.5), lx(0) + (ix + ix0) * dx) &
                        &- max((xr - dxr * 0.5), lx(0) + (ix + ix0 - 1) * dx))

                    fcpspx = dwnrk * dxi / dx
                  else
                    fcpspx = 1.0
                  end if

                  wadr = fcpspx * fcpspy * fcpspz * ray(iRay, ixrv, jyrv, &
                      &kzrv)%dens

                  waveactiondensity(ix, jy, kz) = waveactiondensity(ix, jy, &
                      &kz) + wadr

                  if(include_tracer) then
                    ! leading order gravity wave tracer fluxes
                    ! calculation of equations 2.91 - 2.93 in IK masters thesis

                    if(f_cor_nd /= 0.0) then

                      var_utracer(ix, jy, kz) = var_utracer(ix, jy, kz) &
                          &+ leading_order_tracer_flux(f_cor_nd, omir, wnrk, &
                          &wnrl, wnrm, wadr, 'x', ix, jy, kz, var)
                      var_vtracer(ix, jy, kz) = var_vtracer(ix, jy, kz) &
                          &+ leading_order_tracer_flux(f_cor_nd, omir, wnrk, &
                          &wnrl, wnrm, wadr, 'y', ix, jy, kz, var)
                      var_vtracer(ix, jy, kz) = var_vtracer(ix, jy, kz) &
                          &+ leading_order_tracer_flux(f_cor_nd, omir, wnrk, &
                          &wnrl, wnrm, wadr, 'z', ix, jy, kz, var)

                    end if ! f_cor_nd /= 0.0

                  end if

                end do
              end do
            end do
          end do
        end do
      end do
    end do

    call setboundary_wkb(var_utracer)
    call setboundary_wkb(var_vtracer)
    call setboundary_wkb(var_wtracer)

    if(lsmth_wkb) then
      if(sizeY == 1) then
        if(sizeX > 1) then
          if(sm_filter == 1) then
            call smooth_wkb_box(var_utracer, nsmth_wkb, 101)
            call smooth_wkb_box(var_vtracer, nsmth_wkb, 101)
            call smooth_wkb_box(var_wtracer, nsmth_wkb, 101)
          elseif(sm_filter == 2) then
            call smooth_wkb_shapiro(var_utracer, nsmth_wkb, 101)
            call smooth_wkb_shapiro(var_vtracer, nsmth_wkb, 101)
            call smooth_wkb_shapiro(var_wtracer, nsmth_wkb, 101)
          else
            stop 'WRONG sm_filter'
          end if
        else
          stop 'SMOOTHING JUST IN Z NOT YET IMPLEMENTED'
        endif
      elseif(sizeX == 1) then
        if(sm_filter == 1) then
          call smooth_wkb_box(var_utracer, nsmth_wkb, 11)
          call smooth_wkb_box(var_vtracer, nsmth_wkb, 11)
          call smooth_wkb_box(var_wtracer, nsmth_wkb, 11)
        elseif(sm_filter == 2) then
          call smooth_wkb_shapiro(var_utracer, nsmth_wkb, 11)
          call smooth_wkb_shapiro(var_vtracer, nsmth_wkb, 11)
          call smooth_wkb_shapiro(var_wtracer, nsmth_wkb, 11)
        else
          stop 'WRONG sm_filter'
        end if
      elseif(sizeX > 1) then
        if(sm_filter == 1) then
          call smooth_wkb_box(var_utracer, nsmth_wkb, 111)
          call smooth_wkb_box(var_vtracer, nsmth_wkb, 111)
          call smooth_wkb_box(var_wtracer, nsmth_wkb, 111)
        elseif(sm_filter == 2) then
          call smooth_wkb_shapiro(var_utracer, nsmth_wkb, 111)
          call smooth_wkb_shapiro(var_vtracer, nsmth_wkb, 111)
          call smooth_wkb_shapiro(var_wtracer, nsmth_wkb, 111)
        else
          stop 'WRONG sm_filter'
        end if
      endif
    endif

    do kz = 1, nz
      do jy = 1, ny
        do ix = 1, nx
          if(sizeX > 1) then

            ! d<u'chi'>/dx (for tracer flux convergence)
            dutracer = (var_utracer(ix + 1, jy, kz) - var_utracer(ix - 1, jy, &
                &kz)) / (2.0 * dx)
          else
            dutracer = 0.0
          end if

          if(sizeY > 1) then
            ! d<v'chi'>/dy (for tracer flux convergence)
            dvtracer = (var_vtracer(ix, jy + 1, kz) - var_vtracer(ix, jy - 1, &
                &kz)) / (2.0 * dy)
          else
          end if

          ! tracer flux convergence (leading order gw tracer fluxes)
          tracerforce(ix, jy, kz)%loforce%total = dutracer + dvtracer &
              &+ (var_wtracer(ix, jy, kz + 1) - var_wtracer(ix, jy, kz - 1)) &
              &/ (2.0 * dz)

          ! save diffusive mixing
          tracerforce(ix, jy, kz)%mixingGW%total = 0.
        end do
      end do
    end do

    tracerforce%loforce%uflx = var_utracer
    tracerforce%loforce%vflx = var_vtracer
    tracerforce%loforce%wflx = var_wtracer

    tmparray = tracerforce(0:nx + 1, 0:ny + 1, 0:nz + 1)%loforce%total
    call setboundary_frc_wkb(tmparray)
    tracerforce(0:nx + 1, 0:ny + 1, 0:nz + 1)%loforce%total = tmparray
    tmparray = tracerforce(0:nx + 1, 0:ny + 1, 0:nz + 1)%noforce%total
    call setboundary_frc_wkb(tmparray)
    tracerforce(0:nx + 1, 0:ny + 1, 0:nz + 1)%noforce%total = tmparray
    tmparray = tracerforce(0:nx + 1, 0:ny + 1, 0:nz + 1)%mixingGW%total
    call setboundary_frc_wkb(tmparray)
    tracerforce(0:nx + 1, 0:ny + 1, 0:nz + 1)%mixingGW%total = tmparray

  end subroutine calc_tracerforce

  !-----------------------------------------------------------------------

  subroutine calc_meanFlow_effect(ray, var, force, ray_var3D)

    ! supplemements cell-centered volume forces by WKB force
    ! as well as the heating by entropy-flux convergence

    implicit none

    ! in/out variables
    type(var_type), intent(inout) :: var

    ! volume forcing
    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, 3), intent(inout) :: force

    ! IKJuly2023 changed from 1:6 to 1:11 for tracer fluxes
    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, 1:6), intent(inout) :: &
        &ray_var3D

    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(in) :: ray

    ! real, intent(in) :: time

    ! local variables
    real :: F

    real, allocatable :: var_uu(:, :, :)
    real, allocatable :: var_uv(:, :, :)
    real, allocatable :: var_uw(:, :, :)

    real, allocatable :: var_vv(:, :, :)
    real, allocatable :: var_vw(:, :, :)

    real, allocatable :: var_ETx(:, :, :)
    real, allocatable :: var_ETy(:, :, :)

    real, allocatable :: var_ut(:, :, :)
    real, allocatable :: var_vt(:, :, :)

    real, allocatable :: var_E(:, :, :)

    real, allocatable :: var_drudt(:, :, :)
    real, allocatable :: var_drvdt(:, :, :)
    real, allocatable :: var_drtdt(:, :, :)

    real :: dxi, dyi, dzi

    integer :: ixmin, ixmax
    integer :: jymin, jymax
    integer :: kzmin, kzmax
    integer :: ixr, jyr, kzr
    integer :: ix, jy, kz, im

    integer :: ix0, jy0

    real :: rhotot

    real :: wnrk
    real :: wnrl
    real :: wnrm
    real :: wnrh
    real :: dwnrk, dwnrl, dwnrm
    real :: f_cor_nd

    real :: cgirx, cgiry, cgirz
    real :: omir

    real :: xr, yr, zr
    real :: dxr, dyr, dzr

    real :: fcpspx, fcpspy, fcpspz
    real :: wadr

    real :: NNR
    logical :: apply

    allocate(var_uu(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))
    allocate(var_uv(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))
    allocate(var_uw(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))

    allocate(var_vv(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))
    allocate(var_vw(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))

    allocate(var_ETx(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))
    allocate(var_ETy(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))

    allocate(var_ut(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))
    allocate(var_vt(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))

    allocate(var_E(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))

    allocate(var_drudt(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))
    allocate(var_drvdt(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))
    allocate(var_drtdt(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))

    ! ! Only allow mean flow impact after ray volumes have distributed (FJJul2023)
    ! if(case_wkb == 3 .and. time < topographyTime_wkb / tRef) return

    var_uu = 0.0
    var_uv = 0.0
    var_uw = 0.0

    var_vv = 0.0
    var_vw = 0.0

    var_ETx = 0.0
    var_ETy = 0.0

    var_ut = 0.0
    var_vt = 0.0

    var_E = 0.0

    var_drudt = 0.0
    var_drvdt = 0.0
    var_drtdt = 0.0

    ray_var3D = 0.0

    f_cor_nd = f_Coriolis_dim * tRef ! non dimensional Coriolis param.

    do kzrv = 0, nz + 1
      do jyrv = 0, ny + 1
        ! loop including ghost cells in order to get all fluxes
        ! affecting a cell
        ! (assuming that ray volumes are not wider in y than dy)

        do ixrv = 0, nx + 1
          ! loop including ghost cells in order to get all fluxes
          ! affecting a cell
          ! (assuming that ray volumes are not wider in x than dx)

          if(nRay(ixrv, jyrv, kzrv) < 1) cycle

          do iRay = 1, nRay(ixrv, jyrv, kzrv)
            ! skip counting ray volumes with zero wave-action density

            if(ray(iRay, ixrv, jyrv, kzrv)%dens == 0.0) cycle

            xr = ray(iRay, ixrv, jyrv, kzrv)%x
            yr = ray(iRay, ixrv, jyrv, kzrv)%y
            zr = ray(iRay, ixrv, jyrv, kzrv)%z

            dxr = ray(iRay, ixrv, jyrv, kzrv)%dxray
            dyr = ray(iRay, ixrv, jyrv, kzrv)%dyray
            dzr = ray(iRay, ixrv, jyrv, kzrv)%dzray

            ! vertical boundary conditions:
            ! zBoundary = 'periodic': implement periodicity
            ! zBoundary = 'solid_wall': skip counting ray volumes
            !                           that have completely left the
            !                           model domain

            ! FJApr2023
            !!$            if((.not. topography .and. zr < lz(0)) .or. (topography .and. zr &
            !!$                &< topography_surface(ixrv, jyrv))) then
            apply = .false.
            if(.not. topography .and. zr < lz(0)) then
              apply = .true.
            else if(topography) then
              if(zr < topography_surface(ixrv, jyrv)) then
                apply = .true.
              end if
            end if
            if(apply) then
              select case(zBoundary)
              case("periodic")
                zr = lz(1) + mod(zr - lz(0), lz(1) - lz(0))
              case("solid_wall")
                if((.not. topography .and. zr + 0.5 * dzr < lz(0)) .or. &
                    &(topography .and. zr + 0.5 * dzr &
                    &< topography_surface(ixrv, jyrv))) cycle
              case default
                stop "calc_meanflow_effect: unknown case zBoundary"
              end select
            elseif(zr > lz(1)) then
              select case(zBoundary)
              case("periodic")
                zr = lz(0) + mod(zr - lz(1), lz(1) - lz(0))
              case("solid_wall")
                if(zr - 0.5 * dzr > lz(1)) cycle
              case default
                stop "calc_meanflow_effect: unknown case zBoundary"
              end select
            end if

            ! implement horizontal boundary conditions for ray-volume
            ! positions

            if(sizeX > 1) then
              if(xBoundary /= "periodic") then
                print *, 'ERROR in calc_meanflow_effect:  boundary conditions &
                    &in x must be periodic'
                stop
              end if

              ! for leftmost cpu make sure that xr in
              ! ghost cell to the left is between x(0) - dx
              ! and x(0)

              if(ixrv == 0 .and. is + nbx == 1) then
                if(xr > lx(1) - dx .and. xr < lx(1)) then
                  xr = xr - lx(1) + lx(0)
                elseif(xr > lx(0) - dx .and. xr < lx(0)) then
                  xr = xr
                else
                  print *, 'ERROR in calc_meanflow_effect:'
                  print *, 'xr =', xr
                  print *, 'but r.v. is in ghost cell to the  left of the &
                      &leftmost cpu so that  one should have either'
                  print *, lx(1) - dx, '= lx(1) - dx < xr < lx(1) =', lx(1), ' &
                      &or'
                  print *, lx(0) - dx, '= lx(0) - dx < xr < lx(0) =', lx(0)
                  stop
                end if
              end if

              ! for rightmost cpu make sure that xr in
              ! ghost cell to the right is between x(0) + L_x
              ! and x(0) + L_x + dx

              if(ixrv == nx + 1 .and. is + nbx + nx == sizeX) then
                if(xr > lx(0) .and. xr < lx(0) + dx) then
                  xr = xr + lx(1) - lx(0)
                elseif(xr > lx(1) .and. xr < lx(1) + dx) then
                  xr = xr
                else
                  print *, 'ERROR in calc_meanflow_effect:'
                  print *, 'xr =', xr
                  print *, 'but r.v. is in ghost cell to the  right of the &
                      &rightmost cpu so that  one should have either'
                  print *, lx(0), '= lx(0) < xr < lx(0) + dx =', lx(0) + dx, ' &
                      &or'
                  print *, lx(1), '= lx(1) < xr < lx(1) + dx =', lx(1) + dx
                  stop
                end if
              end if
            end if

            if(sizeY > 1) then
              if(yBoundary /= "periodic") then
                print *, 'ERROR in calc_meanflow_effect:  boundary conditions &
                    &in y must be periodic'
                stop
              end if

              ! for first cpu in y direct. make sure that yr in
              ! ghost cell in front is between y(0) - dy
              ! and y(0)

              if(jyrv == 0 .and. js + nby == 1) then
                if(yr > ly(1) - dy .and. yr < ly(1)) then
                  yr = yr - ly(1) + ly(0)
                elseif(yr > ly(0) - dy .and. yr < ly(0)) then
                  yr = yr
                else
                  print *, 'ERROR in calc_meanflow_effect:'
                  print *, 'yr =', yr
                  print *, 'but r.v. is in ghost cell in front of the first &
                      &cpu in y dir. so that  one should have either'
                  print *, ly(1) - dy, '= ly(1) - dy < yr < ly(1) =', ly(1), ' &
                      &or'
                  print *, ly(0) - dy, '= ly(0) - dy < yr < ly(0) =', ly(0)
                  stop
                end if
              end if

              ! for last cpu in y direction make sure that yr
              ! in ghost cell behind is between
              ! y(0) + L_y and y(0) + L_y + dy

              if(jyrv == ny + 1 .and. js + nby + ny == sizeY) then
                if(yr > ly(0) .and. yr < ly(0) + dy) then
                  yr = yr + ly(1) - ly(0)
                elseif(yr > ly(1) .and. yr < ly(1) + dy) then
                  yr = yr
                else
                  print *, 'ERROR in calc_meanflow_effect:'
                  print *, 'yr =', yr
                  print *, 'but r.v. is in ghost cell behind the last cpu in y &
                      &dir. so that  one should have either'
                  print *, ly(0), '= ly(0) < yr < ly(0) + dy =', ly(0) + dy, ' &
                      &or'
                  print *, ly(1), '= ly(1) < yr < ly(1) + dy =', ly(1) + dy
                  stop
                end if
              end if
            end if

            wnrk = ray(iRay, ixrv, jyrv, kzrv)%k
            wnrl = ray(iRay, ixrv, jyrv, kzrv)%l
            wnrm = ray(iRay, ixrv, jyrv, kzrv)%m

            dwnrk = ray(iRay, ixrv, jyrv, kzrv)%dkray
            dwnrl = ray(iRay, ixrv, jyrv, kzrv)%dlray
            dwnrm = ray(iRay, ixrv, jyrv, kzrv)%dmray

            wnrh = sqrt(wnrk ** 2 + wnrl ** 2)

            !if((.not. topography .and. zr < lz(0) - dz) .or. (topography .and. &
            !    &zr < topography_surface(ixrv, jyrv) - jac(ixrv, jyrv, 0) &
            !    &* dz)) then
            apply = .false.
            if(.not. topography .and. zr < lz(0) - dz) then
              apply = .true.
            else if(topography) then
              if(zr < topography_surface(ixrv, jyrv) - jac(ixrv, jyrv, 0) &
                  &* dz) then
                apply = .true.
              end if
            end if
            if(apply) then
              print *, 'ERROR IN calc_meanflow_effect: RAY VOLUME', iRay, 'in &
                  &cell', ixrv, jyrv, kzrv, 'TOO LOW'
              stop
            end if

            call stratification(zr, 1, NNr)

            omir = branchr * sqrt(NNr * wnrh ** 2 + f_cor_nd ** 2 * wnrm ** 2) &
                &/ sqrt(wnrh ** 2 + wnrm ** 2)

            cgirx = wnrk * (NNr - omir ** 2) / (omir * (wnrh ** 2 + wnrm ** 2))
            cgiry = wnrl * (NNr - omir ** 2) / (omir * (wnrh ** 2 + wnrm ** 2))

            cgirz = - wnrm * (omir ** 2 - f_cor_nd ** 2) / (omir * (wnrh ** 2 &
                &+ wnrm ** 2))

            ! indices of range of cells touched by a ray volume

            if(sizeX > 1) then
              ! last x-index leftmost of cpu
              ix0 = is + nbx - 1

              ixmin = floor((xr - dxr * 0.5 - lx(0)) / dx) + 1 - ix0
              ixmax = floor((xr + dxr * 0.5 - lx(0)) / dx) + 1 - ix0

              if(ixmin > nx + 1) then
                print *, 'ixmin =', ixmin, '> nx+1 = ', nx + 1
                print *, 'ixrv = ', ixrv
                print *, 'xr =', xr
                print *, 'dxr =', dxr
                print *, 'dx =', dx
                print *, 'lx(0) =', lx(0)
                print *, 'lx(1) =', lx(1)
                print *, 'ix0 =', ix0
                print *, 'floor((xr - dxr*0.5 - lx(0)) / dx) + 1 =', floor((xr &
                    &- dxr * 0.5 - lx(0)) / dx) + 1
                stop
              else
                ! no fluxes calculated for the ghost cells
                ! (that are taken care of by the boundary-condition
                ! routines)

                ixmin = max(ixmin, 1)
              end if

              if(ixmax < 0) then
                print *, 'ixmax =', ixmax, '< 0'
                stop
              else
                ! no fluxes calculated for the ghost cells
                ! (that are taken care of by the boundary-condition
                ! routines)

                ixmax = min(ixmax, nx)
              end if
            else
              ixmin = 1
              ixmax = 1
            end if

            if(sizeY > 1) then
              ! last y-index in front of cpu
              jy0 = js + nby - 1

              jymin = floor((yr - dyr * 0.5 - ly(0)) / dy) + 1 - jy0
              jymax = floor((yr + dyr * 0.5 - ly(0)) / dy) + 1 - jy0

              if(jymin > ny + 1) then
                print *, 'jymin =', jymin, '> ny+1 = ', ny + 1
                stop
              else
                ! no fluxes calculated for the ghost cells
                ! (that are taken care of by the boundary-condition
                ! routines)

                jymin = max(jymin, 1)
              end if

              if(jymax < 0) then
                print *, 'jymax =', jymax, '< 0'
                stop
              else
                ! no fluxes calculated for the ghost cells
                ! (that are taken care of by the boundary-condition
                ! routines)

                jymax = min(jymax, ny)
              end if
            else
              jymin = 1
              jymax = 1
            end if

            ! FJMay2023
            if(topography) then

              do ix = ixmin, ixmax
                if(sizeX > 1) then
                  dxi = (min((xr + dxr * 0.5), lx(0) + (ix + ix0) * dx) &
                      &- max((xr - dxr * 0.5), lx(0) + (ix + ix0 - 1) * dx))

                  fcpspx = dwnrk * dxi / dx
                else
                  fcpspx = 1.0
                end if

                do jy = jymin, jymax
                  if(sizeY > 1) then
                    dyi = (min((yr + dyr * 0.5), ly(0) + (jy + jy0) * dy) &
                        &- max((yr - dyr * 0.5), ly(0) + (jy + jy0 - 1) * dy))

                    fcpspy = dwnrl * dyi / dy
                  else
                    fcpspy = 1.0
                  end if

                  kzmin = max(1, floor((levelTFC(ix, jy, zr - dzr * 0.5) &
                      &- lz(0)) / dz) + 1)
                  kzmax = min(nz, floor((levelTFC(ix, jy, zr + dzr * 0.5) &
                      &- lz(0)) / dz) + 1)

                  do kz = kzmin, kzmax
                    dzi = (min((zr + dzr * 0.5), zTildeTFC(ix, jy, kz)) &
                        &- max((zr - dzr * 0.5), zTildeTFC(ix, jy, kz - 1)))

                    fcpspz = dwnrm * dzi / jac(ix, jy, kz) / dz

                    wadr = fcpspx * fcpspy * fcpspz * ray(iRay, ixrv, jyrv, &
                        &kzrv)%dens

                    if(sizeX > 1) then
                      if(f_cor_nd /= 0.0) then
                        var_uu(ix, jy, kz) = var_uu(ix, jy, kz) + wadr * (wnrk &
                            &* cgirx - (wnrk * cgirx + wnrl * cgiry) / (1.0 &
                            &- (omir / f_cor_nd) ** 2))
                      else
                        var_uu(ix, jy, kz) = var_uu(ix, jy, kz) + wadr * wnrk &
                            &* cgirx
                      end if
                    end if

                    if(sizeX > 1 .or. sizeY > 1) then
                      var_uv(ix, jy, kz) = var_uv(ix, jy, kz) + wadr * cgirx &
                          &* wnrl
                    end if

                    var_uw(ix, jy, kz) = var_uw(ix, jy, kz) + wadr * wnrk &
                        &* cgirz / (1.0 - (f_cor_nd / omir) ** 2)

                    if(sizeY > 1) then
                      if(f_cor_nd /= 0.0) then
                        var_vv(ix, jy, kz) = var_vv(ix, jy, kz) + wadr * (wnrl &
                            &* cgiry - (wnrk * cgirx + wnrl * cgiry) / (1.0 &
                            &- (omir / f_cor_nd) ** 2))
                      else
                        var_vv(ix, jy, kz) = var_vv(ix, jy, kz) + wadr * wnrl &
                            &* cgiry
                      end if
                    end if

                    var_vw(ix, jy, kz) = var_vw(ix, jy, kz) + wadr * wnrl &
                        &* cgirz / (1.0 - (f_cor_nd / omir) ** 2)

                    if(f_cor_nd /= 0.0) then
                      var_ETx(ix, jy, kz) = var_ETx(ix, jy, kz) + wadr &
                          &* f_cor_nd ** 2 * NNr * wnrk * wnrm &
                          &/ (rhoStratTFC(ix, jy, kz) * g_ndim * omir * (wnrh &
                          &** 2 + wnrm ** 2))

                      var_ETy(ix, jy, kz) = var_ETy(ix, jy, kz) + wadr &
                          &* f_cor_nd ** 2 * NNr * wnrl * wnrm &
                          &/ (rhoStratTFC(ix, jy, kz) * g_ndim * omir * (wnrh &
                          &** 2 + wnrm ** 2))
                    end if

                    var_E(ix, jy, kz) = var_E(ix, jy, kz) + wadr * omir
                  end do
                end do
              end do
            else
              kzmin = max(1, floor((zr - dzr * 0.5 - lz(0)) / dz) + 1)
              kzmax = min(nz, floor((zr + dzr * 0.5 - lz(0)) / dz) + 1)

              ! calculate momentum-flux / energy / elastic-term
              ! contribution from each ray volume

              do kz = kzmin, kzmax
                dzi = (min((zr + dzr * 0.5), lz(0) + kz * dz) - max((zr - dzr &
                    &* 0.5), lz(0) + (kz - 1) * dz))

                fcpspz = dwnrm * dzi / dz

                do jy = jymin, jymax
                  if(sizeY > 1) then
                    dyi = (min((yr + dyr * 0.5), ly(0) + (jy + jy0) * dy) &
                        &- max((yr - dyr * 0.5), ly(0) + (jy + jy0 - 1) * dy))

                    fcpspy = dwnrl * dyi / dy
                  else
                    fcpspy = 1.0
                  end if

                  do ix = ixmin, ixmax
                    if(sizeX > 1) then
                      dxi = (min((xr + dxr * 0.5), lx(0) + (ix + ix0) * dx) &
                          &- max((xr - dxr * 0.5), lx(0) + (ix + ix0 - 1) * dx))

                      fcpspx = dwnrk * dxi / dx
                    else
                      fcpspx = 1.0
                    end if

                    wadr = fcpspx * fcpspy * fcpspz * ray(iRay, ixrv, jyrv, &
                        &kzrv)%dens

                    if(sizeX > 1) then
                      if(f_cor_nd /= 0.0) then
                        var_uu(ix, jy, kz) = var_uu(ix, jy, kz) + wadr * (wnrk &
                            &* cgirx - (wnrk * cgirx + wnrl * cgiry) / (1.0 &
                            &- (omir / f_cor_nd) ** 2))
                      else
                        var_uu(ix, jy, kz) = var_uu(ix, jy, kz) + wadr * wnrk &
                            &* cgirx
                      end if
                    end if

                    if(sizeX > 1 .or. sizeY > 1) then
                      var_uv(ix, jy, kz) = var_uv(ix, jy, kz) + wadr * cgirx &
                          &* wnrl
                    end if

                    if(steady_state) then
                      var_uw(ix, jy, kz) = var_uw(ix, jy, kz) + wadr * wnrk &
                          &* cgirz
                    else
                      var_uw(ix, jy, kz) = var_uw(ix, jy, kz) + wadr * wnrk &
                          &* cgirz / (1.0 - (f_cor_nd / omir) ** 2)
                    end if

                    if(sizeY > 1) then
                      if(f_cor_nd /= 0.0) then
                        var_vv(ix, jy, kz) = var_vv(ix, jy, kz) + wadr * (wnrl &
                            &* cgiry - (wnrk * cgirx + wnrl * cgiry) / (1.0 &
                            &- (omir / f_cor_nd) ** 2))
                      else
                        var_vv(ix, jy, kz) = var_vv(ix, jy, kz) + wadr * wnrl &
                            &* cgiry
                      end if
                    end if

                    if(steady_state) then
                      var_vw(ix, jy, kz) = var_vw(ix, jy, kz) + wadr * wnrl &
                          &* cgirz
                    else
                      var_vw(ix, jy, kz) = var_vw(ix, jy, kz) + wadr * wnrl &
                          &* cgirz / (1.0 - (f_cor_nd / omir) ** 2)
                    end if

                    if(f_cor_nd /= 0.0) then
                      var_ETx(ix, jy, kz) = var_ETx(ix, jy, kz) + wadr &
                          &* f_cor_nd ** 2 * NNr * wnrk * wnrm / (rhoStrat(kz) &
                          &* g_ndim * omir * (wnrh ** 2 + wnrm ** 2))

                      var_ETy(ix, jy, kz) = var_ETy(ix, jy, kz) + wadr &
                          &* f_cor_nd ** 2 * NNr * wnrl * wnrm / (rhoStrat(kz) &
                          &* g_ndim * omir * (wnrh ** 2 + wnrm ** 2))
                    end if

                    var_E(ix, jy, kz) = var_E(ix, jy, kz) + wadr * omir

                  end do
                end do
              end do
            end if
          end do
        end do
      end do
    end do

    ! for output of u'w', E_w, u, v, w:
    do kz = 1, nz
      do jy = 1, ny
        do ix = 1, nx
          ray_var3D(ix, jy, kz, 4) = var_uw(ix, jy, kz)
          ray_var3D(ix, jy, kz, 5) = var_vw(ix, jy, kz)
          ray_var3D(ix, jy, kz, 6) = var_E(ix, jy, kz)
        end do
      end do
    end do

    ! horizontal entropy fluxes
    if(f_cor_nd /= 0.0) then
      if(topography) then
        ! FJApr2023
        do kz = 1, nz
          do jy = 1, ny
            do ix = 1, nx
              var_ut(ix, jy, kz) = thetaStratTFC(ix, jy, kz) / f_cor_nd &
                  &* var_ETy(ix, jy, kz)
              var_vt(ix, jy, kz) = - thetaStratTFC(ix, jy, kz) / f_cor_nd &
                  &* var_ETx(ix, jy, kz)
            end do
          end do
        end do
      else
        do kz = 1, nz
          var_ut(:, :, kz) = thetaStrat(kz) / f_cor_nd * var_ETy(:, :, kz)
          var_vt(:, :, kz) = - thetaStrat(kz) / f_cor_nd * var_ETx(:, :, kz)
        end do
      end if
    end if

    if(steady_state) then
      var_uu = 0.0
      var_uv = 0.0
      var_vv = 0.0
      var_ETx = 0.0
      var_ETy = 0.0
      var_ut = 0.0
      var_vt = 0.0
    end if

    ! set boundary conditions for all fluxes

    ! vertical boundary conditions also set for the horizontal fluxes in
    ! order to be prepared for the vertical smooting

    call setboundary_wkb(var_uu)
    call setboundary_wkb(var_uv)
    call setboundary_wkb(var_uw)

    call setboundary_wkb(var_vv)
    call setboundary_wkb(var_vw)

    call setboundary_wkb(var_ETx)
    call setboundary_wkb(var_ETy)

    call setboundary_wkb(var_ut)
    call setboundary_wkb(var_vt)

    call setboundary_wkb(var_E)

    ! wave impact on horizontal momentum

    do kz = 1, nz
      do jy = 1, ny
        do ix = 1, nx
          select case(model)
          case("Boussinesq")
            rhotot = rho00

            ! Note by Junhong Wei (20170814): There may be
            ! something wrong with the case of Boussinesq here.
            ! The case of pseudo_incompressible should be fine.
          case("pseudo_incompressible", "compressible")
            rhotot = var%rho(ix, jy, kz)
            if(topography) then
              rhotot = rhotot + rhoStratTFC(ix, jy, kz)
            else
              rhotot = rhotot + rhoStrat(kz)
            end if
          case default
            stop "volumeForce: unknown case model."
          end select

          ! forcing in x direction

          if(topography) then
            ! FJJul2023
            var_drudt(ix, jy, kz) = - rhotot / rhoStratTFC(ix, jy, kz) &
                &/ jac(ix, jy, kz) * (var_uw(ix, jy, kz + 1) - var_uw(ix, jy, &
                &kz - 1)) / (2.0 * dz)

            if(sizeX > 1) then
              var_drudt(ix, jy, kz) = var_drudt(ix, jy, kz) - rhotot &
                  &/ rhoStratTFC(ix, jy, kz) / jac(ix, jy, kz) * ((jac(ix + 1, &
                  &jy, kz) * var_uu(ix + 1, jy, kz) - jac(ix - 1, jy, kz) &
                  &* var_uu(ix - 1, jy, kz)) / (2.0 * dx) + (jac(ix, jy, kz &
                  &+ 1) * met(ix, jy, kz + 1, 1, 3) * var_uu(ix, jy, kz + 1) &
                  &- jac(ix, jy, kz - 1) * met(ix, jy, kz - 1, 1, 3) &
                  &* var_uu(ix, jy, kz - 1)) / (2.0 * dz))
            end if

            if(sizeY > 1) then
              var_drudt(ix, jy, kz) = var_drudt(ix, jy, kz) - rhotot &
                  &/ rhoStratTFC(ix, jy, kz) / jac(ix, jy, kz) * ((jac(ix, jy &
                  &+ 1, kz) * var_uv(ix, jy + 1, kz) - jac(ix, jy - 1, kz) &
                  &* var_uv(ix, jy - 1, kz)) / (2.0 * dy) + (jac(ix, jy, kz &
                  &+ 1) * met(ix, jy, kz + 1, 2, 3) * var_uv(ix, jy, kz + 1) &
                  &- jac(ix, jy, kz - 1) * met(ix, jy, kz - 1, 2, 3) &
                  &* var_uv(ix, jy, kz - 1)) / (2.0 * dz))
            end if
          else
            var_drudt(ix, jy, kz) = - rhotot / rhoStrat(kz) * (var_uw(ix, jy, &
                &kz + 1) - var_uw(ix, jy, kz - 1)) / (2.0 * dz)

            if(sizeX > 1) then
              var_drudt(ix, jy, kz) = var_drudt(ix, jy, kz) - rhotot &
                  &/ rhoStrat(kz) * (var_uu(ix + 1, jy, kz) - var_uu(ix - 1, &
                  &jy, kz)) / (2.0 * dx)
            end if

            if(sizeY > 1) then
              var_drudt(ix, jy, kz) = var_drudt(ix, jy, kz) - rhotot &
                  &/ rhoStrat(kz) * (var_uv(ix, jy + 1, kz) - var_uv(ix, jy &
                  &- 1, kz)) / (2.0 * dy)
            end if
          end if

          var_drudt(ix, jy, kz) = var_drudt(ix, jy, kz) + rhotot * var_ETx(ix, &
              &jy, kz)

          ! forcing in y direction

          select case(model)
          case("Boussinesq")
            rhotot = rho00

            ! Note by Junhong Wei (20170814): There may be
            ! something wrong with the case of Boussinesq here.
            ! The case of pseudo_incompressible should be fine.
          case("pseudo_incompressible", "compressible")
            ! rhotot = 0.5 * (var(ix, jy, kz, 1) + var(ix, jy + 1, kz, 1))
            rhotot = var%rho(ix, jy, kz)
            if(topography) then
              ! FJJul2023
              ! rhotot = rhotot + 0.5 * (rhoStratTFC(ix, jy, kz) &
              !     + rhoStratTFC(ix, jy + 1, kz))
              rhotot = rhotot + rhoStratTFC(ix, jy, kz)
            else
              rhotot = rhotot + rhoStrat(kz)
            end if
          case default
            stop "volumeForce: unknown case model."
          end select

          if(topography) then
            ! FJJul2023
            var_drvdt(ix, jy, kz) = - rhotot / rhoStratTFC(ix, jy, kz) &
                &/ jac(ix, jy, kz) * (var_vw(ix, jy, kz + 1) - var_vw(ix, jy, &
                &kz - 1)) / (2.0 * dz)

            if(sizeX > 1) then
              var_drvdt(ix, jy, kz) = var_drvdt(ix, jy, kz) - rhotot &
                  &/ rhoStratTFC(ix, jy, kz) / jac(ix, jy, kz) * ((jac(ix + 1, &
                  &jy, kz) * var_uv(ix + 1, jy, kz) - jac(ix - 1, jy, kz) &
                  &* var_uv(ix - 1, jy, kz)) / (2.0 * dx) + (jac(ix, jy, kz &
                  &+ 1) * met(ix, jy, kz + 1, 1, 3) * var_uv(ix, jy, kz + 1) &
                  &- jac(ix, jy, kz - 1) * met(ix, jy, kz - 1, 1, 3) &
                  &* var_uv(ix, jy, kz - 1)) / (2.0 * dz))
            end if

            if(sizeY > 1) then
              var_drvdt(ix, jy, kz) = var_drvdt(ix, jy, kz) - rhotot &
                  &/ rhoStratTFC(ix, jy, kz) / jac(ix, jy, kz) * ((jac(ix, jy &
                  &+ 1, kz) * var_vv(ix, jy + 1, kz) - jac(ix, jy - 1, kz) &
                  &* var_vv(ix, jy - 1, kz)) / (2.0 * dy) + (jac(ix, jy, kz &
                  &+ 1) * met(ix, jy, kz + 1, 2, 3) * var_vv(ix, jy, kz + 1) &
                  &- jac(ix, jy, kz - 1) * met(ix, jy, kz - 1, 2, 3) &
                  &* var_vv(ix, jy, kz - 1)) / (2.0 * dz))
            end if
          else
            var_drvdt(ix, jy, kz) = - rhotot / rhoStrat(kz) * (var_vw(ix, jy, &
                &kz + 1) - var_vw(ix, jy, kz - 1)) / (2.0 * dz)

            if(sizeX > 1) then
              var_drvdt(ix, jy, kz) = var_drvdt(ix, jy, kz) - rhotot &
                  &/ rhoStrat(kz) * (var_uv(ix + 1, jy, kz) - var_uv(ix - 1, &
                  &jy, kz)) / (2.0 * dx)
            end if

            if(sizeY > 1) then
              var_drvdt(ix, jy, kz) = var_drvdt(ix, jy, kz) - rhotot &
                  &/ rhoStrat(kz) * (var_vv(ix, jy + 1, kz) - var_vv(ix, jy &
                  &- 1, kz)) / (2.0 * dy)
            end if
          end if

          var_drvdt(ix, jy, kz) = var_drvdt(ix, jy, kz) + rhotot * var_ETy(ix, &
              &jy, kz)
        end do
      end do
    end do

    ! negative wave-induced heating by GWs

    if((f_cor_nd /= 0.0) .and. (sizeX > 1 .or. sizeY > 1)) then
      do kz = 1, nz
        do jy = 1, ny
          do ix = 1, nx
            select case(model)
            case("Boussinesq")
              rhotot = rho00

              ! Note by Junhong Wei (20170814): There may be
              ! something wrong with the case of Boussinesq here.
              ! The case of pseudo_incompressible should be fine.
            case("pseudo_incompressible", "compressible")
              if(topography) then
                ! FJApr2023
                rhotot = var%rho(ix, jy, kz) + rhoStratTFC(ix, jy, kz)
              else
                rhotot = var%rho(ix, jy, kz) + rhoStrat(kz)
              end if
            case default
              stop "volumeForce: unknown case model."
            end select

            if(topography) then
              ! FJJul2023
              if(sizeX > 1) then
                var_drtdt(ix, jy, kz) = rhotot / jac(ix, jy, kz) * ((jac(ix &
                    &+ 1, jy, kz) * var_ut(ix + 1, jy, kz) - jac(ix - 1, jy, &
                    &kz) * var_ut(ix - 1, jy, kz)) / (2.0 * dx) + (jac(ix, jy, &
                    &kz + 1) * met(ix, jy, kz + 1, 1, 3) * var_ut(ix, jy, kz &
                    &+ 1) - jac(ix, jy, kz - 1) * met(ix, jy, kz - 1, 1, 3) &
                    &* var_ut(ix, jy, kz - 1)) / (2.0 * dz))
              end if

              if(sizeY > 1) then
                var_drtdt(ix, jy, kz) = var_drtdt(ix, jy, kz) + rhotot &
                    &/ jac(ix, jy, kz) * ((jac(ix, jy + 1, kz) * var_vt(ix, jy &
                    &+ 1, kz) - jac(ix, jy - 1, kz) * var_vt(ix, jy - 1, kz)) &
                    &/ (2.0 * dy) + (jac(ix, jy, kz + 1) * met(ix, jy, kz + 1, &
                    &2, 3) * var_vt(ix, jy, kz + 1) - jac(ix, jy, kz - 1) &
                    &* met(ix, jy, kz - 1, 2, 3) * var_vt(ix, jy, kz - 1)) &
                    &/ (2.0 * dz))
              end if
            else
              if(sizeX > 1) then
                var_drtdt(ix, jy, kz) = rhotot * (var_ut(ix + 1, jy, kz) &
                    &- var_ut(ix - 1, jy, kz)) / (2.0 * dx)
              end if

              if(sizeY > 1) then
                var_drtdt(ix, jy, kz) = var_drtdt(ix, jy, kz) + rhotot &
                    &* (var_vt(ix, jy + 1, kz) - var_vt(ix, jy - 1, kz)) &
                    &/ (2.0 * dy)
              end if
            end if
          end do
        end do
      end do
    end if

    ! boundary conditions (also as preparation for the smoothing)

    call setboundary_wkb(var_drudt)
    call setboundary_wkb(var_drvdt)
    call setboundary_wkb(var_drtdt)

    ! running average for the wave impacts

    if(lsmth_wkb) then
      if(sizeY == 1) then
        if(sizeX > 1) then
          if(sm_filter == 1) then
            call smooth_wkb_box(var_drudt, nsmth_wkb, 101)
            call smooth_wkb_box(var_drvdt, nsmth_wkb, 101)
            call smooth_wkb_box(var_drtdt, nsmth_wkb, 101)
          elseif(sm_filter == 2) then
            call smooth_wkb_shapiro(var_drudt, nsmth_wkb, 101)
            call smooth_wkb_shapiro(var_drvdt, nsmth_wkb, 101)
            call smooth_wkb_shapiro(var_drtdt, nsmth_wkb, 101)
          else
            stop 'WRONG sm_filter'
          end if
        else
          stop 'SMOOTHING JUST IN Z NOT YET IMPLEMENTED'
        endif
      elseif(sizeX == 1) then
        if(sm_filter == 1) then
          call smooth_wkb_box(var_drudt, nsmth_wkb, 11)
          call smooth_wkb_box(var_drvdt, nsmth_wkb, 11)
          call smooth_wkb_box(var_drtdt, nsmth_wkb, 11)
        elseif(sm_filter == 2) then
          call smooth_wkb_shapiro(var_drudt, nsmth_wkb, 11)
          call smooth_wkb_shapiro(var_drvdt, nsmth_wkb, 11)
          call smooth_wkb_shapiro(var_drtdt, nsmth_wkb, 11)
        else
          stop 'WRONG sm_filter'
        end if
      elseif(sizeX > 1) then
        if(sm_filter == 1) then
          call smooth_wkb_box(var_drudt, nsmth_wkb, 111)
          call smooth_wkb_box(var_drvdt, nsmth_wkb, 111)
          call smooth_wkb_box(var_drtdt, nsmth_wkb, 111)
        elseif(sm_filter == 2) then
          call smooth_wkb_shapiro(var_drudt, nsmth_wkb, 111)
          call smooth_wkb_shapiro(var_drvdt, nsmth_wkb, 111)
          call smooth_wkb_shapiro(var_drtdt, nsmth_wkb, 111)
        else
          stop 'WRONG sm_filter'
        end if
      endif
    endif

    ! add wave impact to tendencies ...

    do kz = 1, nz
      ! only allow wave impact on mean flow above lz(0) + zmin_wkb
      ! FJApr2023
      ! if(z(kz) < lz(0) + zmin_wkb) cycle
      if(.not. topography .and. z(kz) < lz(0) + zmin_wkb) cycle

      do jy = 1, ny
        do ix = 1, nx
          ! FJApr2023
          !*if(topography .and. zTFC(ix, jy, kz) < lz(0) + zmin_wkb) cycle
          ! SD
          if(topography) then
            if(zTFC(ix, jy, kz) < lz(0) + zmin_wkb) cycle
          end if
          select case(model)
          case("Boussinesq")
            rhotot = rho00

            ! Note by Junhong Wei (20170814): There may be
            ! something wrong with the case of Boussinesq here.
            ! The case of pseudo_incompressible should be fine.
          case("pseudo_incompressible", "compressible")
            ! rhotot = 0.5 * (var(ix, jy, kz, 1) + var(ix + 1, jy, kz, 1))
            rhotot = var%rho(ix, jy, kz)
            if(topography) then
              ! FJJul2023
              ! rhotot = rhotot + 0.5 * (rhoStratTFC(ix, jy, kz) &
              !     + rhoStratTFC(ix + 1, jy, kz))
              rhotot = rhotot + rhoStratTFC(ix, jy, kz)
            else
              rhotot = rhotot + rhoStrat(kz)
            end if
          case default
            stop "volumeForce: unknown case model."
          end select

          ! forcing in x direction

          force(ix, jy, kz, 1) = force(ix, jy, kz, 1) + var_drudt(ix, jy, kz)

          ! for output of mean-flow acceleration in x direction by GWs
          ray_var3D(ix, jy, kz, 1) = var_drudt(ix, jy, kz) / rhotot

          ! forcing in y direction

          select case(model)
          case("Boussinesq")
            rhotot = rho00

            ! Note by Junhong Wei (20170814): There may be
            ! something wrong with the case of Boussinesq here.
            ! The case of pseudo_incompressible should be fine.
          case("pseudo_incompressible", "compressible")
            ! rhotot = 0.5 * (var(ix, jy, kz, 1) + var(ix, jy + 1, kz, 1))
            rhotot = var%rho(ix, jy, kz)
            if(topography) then
              rhotot = rhotot + rhoStratTFC(ix, jy, kz)
            else
              rhotot = rhotot + rhoStrat(kz)
            end if
          case default
            stop "volumeForce: unknown case model."
          end select

          force(ix, jy, kz, 2) = force(ix, jy, kz, 2) + var_drvdt(ix, jy, kz)

          ! for output of mean-flow acceleration in y direction by GWs
          ray_var3D(ix, jy, kz, 2) = var_drvdt(ix, jy, kz) / rhotot

          ! Add forcing on terrain-following wind (FJJul2023).
          if(topography) then
            force(ix, jy, kz, 3) = force(ix, jy, kz, 3) + met(ix, jy, kz, 1, &
                &3) * var_drudt(ix, jy, kz) + met(ix, jy, kz, 2, 3) &
                &* var_drvdt(ix, jy, kz)
          end if
        end do
      end do
    end do

    call setboundary_frc_wkb(force(:, :, :, 1))
    call setboundary_frc_wkb(force(:, :, :, 2))

    var%GWH(:, :, :) = 0.0

    if((f_cor_nd /= 0.0) .and. (sizeX > 1 .or. sizeY > 1)) then
      do kz = 1, nz
        ! only allow wave impact on mean flow above lz(0) + zmin_wkb
        if(z(kz) < lz(0) + zmin_wkb) cycle

        do jy = 1, ny
          do ix = 1, nx
            select case(model)
            case("Boussinesq")
              rhotot = rho00

              ! Note by Junhong Wei (20170814): There may be
              ! something wrong with the case of Boussinesq here.
              ! The case of pseudo_incompressible should be fine.
            case("pseudo_incompressible", "compressible")
              if(topography) then
                ! FJApr2023
                rhotot = var%rho(ix, jy, kz) + rhoStratTFC(ix, jy, kz)
              else
                rhotot = var%rho(ix, jy, kz) + rhoStrat(kz)
              end if
            case default
              stop "volumeForce: unknown case model."
            end select

            var%GWH(ix, jy, kz) = var_drtdt(ix, jy, kz)

            ! for output of mean-flow potential-temperature tendency
            ! by GWs
            ray_var3D(ix, jy, kz, 3) = - var_drtdt(ix, jy, kz) / rhotot
          end do
        end do
      end do
    end if

  end subroutine calc_meanFlow_effect

  !---------------------------------------------------------------------

  subroutine setup_wkb(ray, ray_var3D, var, diffusioncoeff, waveAmplitudes, &
      &dPhase, ray_varIce, ray_cloud)

    !------------------------------------------------
    ! allocate ray field
    ! initialize position, wave vector and frequency
    ! for rays
    !------------------------------------------------

    implicit none

    ! argument list

    type(rayType), dimension(:, :, :, :), allocatable, intent(out) :: ray
    real, dimension(:, :, :, :), allocatable, intent(out) :: ray_var3D
    type(var_type), intent(inout) :: var

    ! turbulent eddy diffusivity. Needed for the mixing of tracer
    real, dimension(:, :, :), allocatable, intent(out) :: diffusioncoeff
    real, dimension(:, :, :), allocatable, intent(out) :: dPhase
    type(ice_rayType), dimension(:, :, :), allocatable, intent(out) :: &
        &ray_varIce
    type(ice_rayType2), dimension(:, :, :, :, :), allocatable, intent(out) :: &
        &ray_cloud

    type(waveAmpType), dimension(:, :, :), allocatable, intent(out) :: &
        &waveAmplitudes

    ! local variables
    integer :: allocstat

    ! FJApr2023
    real, allocatable :: omi_notop(:, :, :), omi_sfc(:, :, :), wnk_sfc(:, :, &
        &:), wnl_sfc(:, :, :), wnm_sfc(:, :, :)
    real, allocatable :: fld_amp(:, :, :, :)

    integer :: ix, jy, kz, im
    integer :: ix2, jy2, kz2
    integer :: ik, jl, km
    integer :: ixmin, ixmax
    integer :: jymin, jymax
    integer :: kzmin, kzmax
    integer :: kz2min

    integer :: ix0, jy0

    real :: wnrh
    real :: NN_nd
    real :: wnrh_init, wnrk_init, wnrl_init, wnrm_init
    real :: wnrk, wnrl, wnrm
    real :: cgirx, cgiry, cgirz
    real :: omir
    real :: NNr
    real :: f_cor_nd
    real :: sigwpx, sigwpy, sigwpz
    real :: xr0, yr0, zr0
    real :: xrmin, xrmax
    real :: yrmin, yrmax
    real :: zrmin, zrmax
    real :: displm
    real :: dk_ini_nd, dl_ini_nd, dm_ini_nd
    real :: wnk_0, wnl_0, wnm_0
    real :: pspvol
    real :: uxr, vyr, wzr
    real :: xr, yr, zr

    integer :: nrsuml, nr_sum

    ! Long number (FJJan2023)
    real :: long

    ! Wave mode (FJApr2023)
    integer :: iwm

    !SD
    real :: dphi

    ix0 = is + nbx - 1
    jy0 = js + nby - 1

    ! some non-dimensionalizations

    f_cor_nd = f_Coriolis_dim * tRef

    xr0 = xr0_dim / lRef
    yr0 = yr0_dim / lRef
    zr0 = zr0_dim / lRef

    sigwpx = sigwpx_dim / lRef
    sigwpy = sigwpy_dim / lRef
    sigwpz = sigwpz_dim / lRef

    xrmin = xrmin_dim / lRef
    xrmax = xrmax_dim / lRef

    yrmin = yrmin_dim / lRef
    yrmax = yrmax_dim / lRef

    zrmin = zrmin_dim / lRef
    zrmax = zrmax_dim / lRef

    !SDJul2024
    if(case_wkb == 5) then

      xr0_sp = xr0_dim_sp / lRef
      yr0_sp = yr0_dim_sp / lRef
      zr0_sp = zr0_dim_sp / lRef

      sigwpx_sp = sigwpx_dim_sp / lRef
      sigwpy_sp = sigwpy_dim_sp / lRef
      sigwpz_sp = sigwpz_dim_sp / lRef

    end if
    !-------------------------------------------
    ! compute maximum number of ray volumes  ...
    !-------------------------------------------

    ! factor z-m space

    if(zrmin < lz(0) .or. zrmax > lz(1)) then
      print *, 'zrmin too small or zrmax too large! --> exit'
      stop
    endif

    if(case_wkb == 3) then
      kzmin = 0
      kzmax = 0

      nzRay = nray_fac * nrzl * nrm_init
    else
      if(topography) then
        kzmin = 1
        kzmax = sizeZ
      else
        kzmin = max(1, int(floor((zrmin - lz(0)) / dz)) + 1)
        kzmax = min(sizeZ, int(floor((zrmax - lz(0)) / dz)) + 1)
      end if

      nzRay = nray_fac * nrzl * nrm_init
    end if

    ! factor x-k space
    if(case_wkb == 3) then ! topography
      ixmin = 1
      ixmax = nx
    else
      if(xrmin < lx(0) .or. xrmax > lx(1)) then
        print *, 'xrmin too small or xrmax too large! --> exit'
        stop
      endif

      ixmin = max(1, int(floor((xrmin - lx(0)) / dx)) + 1 - ix0)
      ixmax = min(nx, int(floor((xrmax - lx(0)) / dx)) + 1 - ix0)

      ! if the cpu domain is outside of the range where r.v. are to be
      ! generated one gets ixmin > ixmax. In this case do not do anything
      ! below
    end if

    if(sizeX == 1) then
      nxRay = 1
    else
      nxRay = nray_fac * nrxl * nrk_init
    end if

    ! factor y-l space
    if(case_wkb == 3) then
      jymin = 1
      jymax = ny
    else
      jymin = max(1, int(floor((yrmin - ly(0)) / dy)) + 1)
      jymax = min(ny, int(floor((yrmax - ly(0)) / dy)) + 1)

      ! if the cpu domain is outside of the range where r.v. are to be
      ! generated one gets jymin > jymax. In this case do not do anything
      ! below
    end if

    if(sizeY == 1) then
      nyRay = 1
    else
      nyRay = nray_fac * nryl * nrl_init
    end if

    ! maximum # of r.v. allowed in a cell before r.v. are merged
    ! FJApr2023
    nray_max = nxRay * nyRay * nzRay * nwm

    ! work-space size per wavenumber direction chosen to be double the
    ! maximum number of r.v. allowed before they are merged
    ! (should this turn out to be too small, change to triple, quadruple,
    ! etc)

    if(nxRay > 1) then
      nxRay_wrk = 2 * nxRay
    else
      nxRay_wrk = 1
    end if

    if(nyRay > 1) then
      nyRay_wrk = 2 * nyRay
    else
      nyRay_wrk = 1
    end if

    if(nzRay > 1) then
      nzRay_wrk = 2 * nzRay
    else
      nzRay_wrk = 1
    end if

    nray_wrk = nxRay_wrk * nyRay_wrk * nzRay_wrk

    ! Limit ray volume output.
    if(nRayOutput > nray_wrk) stop "Error: nRayOutput > nray_wrk"

    ! FJApr2023
    n_sfc = nwm
    if(nxRay > 1) n_sfc = n_sfc * nxRay / nray_fac
    if(nyRay > 1) n_sfc = n_sfc * nyRay / nray_fac
    if(nzRay > 1) n_sfc = n_sfc * nzRay / nray_fac

    ! field of ray volumes
    allocate(ray(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate ray"

    ! # of ray volumes per cell
    allocate(nRay(0:nx + 1, 0:ny + 1, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate nRay"

    nRay = 0

    if(case_wkb == 3) then
      ! pointers to surface ray volumes

      allocate(ir_sfc(n_sfc, 0:nx + 1, 0:ny + 1), stat = allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate ir_sfc"

      allocate(ix2_sfc(n_sfc), stat = allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate ix2_sfc"

      allocate(jy2_sfc(n_sfc), stat = allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate jy2_sfc"

      allocate(kz2_sfc(n_sfc), stat = allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate kz2_sfc"

      allocate(ik_sfc(n_sfc), stat = allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate ik_sfc"

      allocate(jl_sfc(n_sfc), stat = allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate jl_sfc"

      allocate(km_sfc(n_sfc), stat = allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate km_sfc"

      allocate(iwm_sfc(n_sfc), stat = allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate iwm_sfc"
    end if

    ! position displacement increment
    allocate(dxRay(3, nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate dxRay"

    ! wave vector increment
    allocate(dkRay(3, nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate dkRay"

    ! ray-volume extent increment
    allocate(ddxRay(3, nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), stat &
        &= allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate ddxRay"

    if(include_ice) then
      ! phase increment
      allocate(dpRay(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), stat &
          &= allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate dpRay"

      if(raytracer) then
        allocate(ray_varIce(0:nx + 1, 0:ny + 1, 0:nz + 1), stat = allocstat)
        if(allocstat /= 0) stop "setup_wkb: could not allocate ray_varIce"

        !if ( compute_cloudcover ) then
        !should be only allocated if cloud cover used/ requires to rewrite subroutines
        allocate(ray_cloud(0:nx + 1, 0:ny + 1, 0:nz + 1, NSCX, NSCY), stat &
            &= allocstat)
        if(allocstat /= 0) stop "setup_wkb: could not allocate ray_cloud"

        allocate(field_mst_cld(sizeX * NSCX * nprocy, ny * NSCY), stat &
            &= allocstat)
        allocate(field_out_cld(sizeX * NSCX, sizeY * NSCY), stat = allocstat)
        !end if

        !has to be placed probably somewhere else
        if(compute_cloudcover) then

          dxsc = dx / nscx
          dysc = dy / nscy

          do kz = 1, nz
            do jy = 1, ny
              do ix = 1, nx
                ray_cloud(ix, jy, kz, :, :)%Ni = var%ICE(ix, jy, kz, inN)
                ray_cloud(ix, jy, kz, :, :)%Qi = var%ICE(ix, jy, kz, inQ)
                ray_cloud(ix, jy, kz, :, :)%Qv = var%ICE(ix, jy, kz, inQv)
              end do
            end do
          end do
        end if
      end if
    end if

    ! fields for data WKB output
    allocate(ray_var3D(0:nx + 1, 0:ny + 1, 0:nz + 1, 1:6), stat = allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate ray_var3D"

    allocate(diffusioncoeff(0:nx + 1, 0:ny + 1, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate diffusioncoeff"

    allocate(dPhase(0:nx + 1, 0:ny + 1, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate dPhase"

    allocate(waveAmplitudes(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &stat = allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate waveAmplitudes"

    ! needed for initialization of ray volumes:
    if(case_wkb == 3) then
      ! FJApr2023
      allocate(omi_sfc(1:nx, 1:ny, 1:nwm))
      allocate(wnk_sfc(1:nx, 1:ny, 1:nwm))
      allocate(wnl_sfc(1:nx, 1:ny, 1:nwm))
      allocate(wnm_sfc(1:nx, 1:ny, 1:nwm))
    else
      ! FJApr2023
      allocate(omi_notop(1:nx, 1:ny, 1:sizeZ))
    end if

    ! FJApr2023
    allocate(fld_amp(1:nx, 1:ny, 0:sizeZ, 1:nwm)) ! 3D wave action field

    ! FJApr2023
    ! Store TFC levels for interpolations. Note that ray volumes beyond the
    ! vertical boundaries are interpolated at these boundaries.
    if(topography) then
      allocate(zTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz + nbz)))
      allocate(zTildeTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
          &+ nbz)))
      do ix = - nbx, nx + nbx
        do jy = - nby, ny + nby
          do kz = - nbz, nz + nbz
            zTFC(ix, jy, kz) = heightTFC(ix, jy, kz)
            zTildeTFC(ix, jy, kz) = heightTFC(ix, jy, kz) + 0.5 * jac(ix, jy, &
                &kz) * dz
          end do
        end do
      end do
    end if

    if(steady_state .and. case_wkb /= 3) stop "Steady state is implemented for &
        &case_wkb == 3 only!"

    ! non-dimensional wave numbers

    if(wlrx_init /= 0.0) then
      wnrk_init = 2.0 * pi / wlrx_init * lRef
    else
      wnrk_init = 0.0
    end if

    if(wlry_init /= 0.0) then
      wnrl_init = 2.0 * pi / wlry_init * lRef
    else
      wnrl_init = 0.0
    end if

    wnrh_init = sqrt(wnrk_init ** 2 + wnrl_init ** 2)

    wnrm_init = 2.0 * pi / wlrz_init * lRef

    if(case_wkb == 5) then

      do iwm = 1, NWM_WP

        if(wlrx_init_sp(iwm) /= 0.0) then
          wnrk_init_sp(iwm) = 2.0 * pi / wlrx_init_sp(iwm) * lRef
        else
          wnrk_init_sp(iwm) = 0.0
        end if

        if(wlry_init_sp(iwm) /= 0.0) then
          wnrl_init_sp(iwm) = 2.0 * pi / wlry_init_sp(iwm) * lRef
        else
          wnrl_init_sp(iwm) = 0.0
        end if

        wnrm_init_sp(iwm) = 2.0 * pi / wlrz_init_sp(iwm) * lRef
      end do

    end if

    ! achatzc:
    ! a slight inconsistency below is that the stratification
    ! is calculated at the cell centers, while it is later on
    ! interpolated to the ray-volume positions
    ! this is no issue in the isothermal case

    if(case_wkb == 3) then
      ! intrinsic frequency and horizontal wave number mountain wave

      kz = 0
      if(topography) then
        do jy = 1, ny
          do ix = 1, nx
            ! Local squared buoyancy frequency
            call stratification(zTFC(ix, jy, 1), 1, NN_nd)
            do iwm = 1, nwm
              ! Wavenumbers
              wnrk_init = k_spectrum(ix, jy, iwm)
              wnrl_init = l_spectrum(ix, jy, iwm)
              wnrh_init = sqrt(wnrk_init ** 2.0 + wnrl_init ** 2.0)
              wnrm_init = 0.0

              ! Intrinsic frequency
              omi_sfc(ix, jy, iwm) = - 0.5 * (var%u(ix, jy, 1) + var%u(ix - 1, &
                  &jy, 1)) * wnrk_init - 0.5 * (var%v(ix, jy, 1) + var%v(ix, &
                  &jy - 1, 1)) * wnrl_init

              ! Frequency branch
              if(omi_sfc(ix, jy, iwm) * branchr >= 0.0) then
                wnk_sfc(ix, jy, iwm) = wnrk_init
                wnl_sfc(ix, jy, iwm) = wnrl_init
              else
                omi_sfc(ix, jy, iwm) = - omi_sfc(ix, jy, iwm)

                wnk_sfc(ix, jy, iwm) = - wnrk_init
                wnl_sfc(ix, jy, iwm) = - wnrl_init
              end if

              ! Wave action density and vertical wavenumber
              if(abs(omi_sfc(ix, jy, iwm)) <= f_cor_nd) then
                fld_amp(ix, jy, kz, iwm) = 0.0
                wnrm = 0.0
              else if(abs(omi_sfc(ix, jy, iwm)) < sqrt(NN_nd)) then
                wnrm = - branchr * sqrt(wnrh_init ** 2 * (NN_nd - omi_sfc(ix, &
                    &jy, iwm) ** 2) / (omi_sfc(ix, jy, iwm) ** 2 - f_cor_nd &
                    &** 2))

                ! Displacement
                displm = abs(topography_spectrum(ix, jy, iwm))

                ! Long number scaling
                if(blocking) then
                  ! Compute Long number.
                  long = sqrt(NN_nd / (0.25 * (var%u(ix, jy, 1) + var%u(ix &
                      &- 1, jy, 1)) ** 2.0 + 0.25 * (var%v(ix, jy, 1) &
                      &+ var%v(ix, jy - 1, 1)) ** 2.0)) &
                      &* sum(abs(topography_spectrum(ix, jy, :)))
                  ! Apply scaling.
                  displm = displm * wave_amplitude_reduction(long)
                end if

                ! Surface wave-action density
                fld_amp(ix, jy, kz, iwm) = 0.5 * rhoStratTFC(ix, jy, 1) &
                    &* displm ** 2 * omi_sfc(ix, jy, iwm) * (wnrh_init ** 2 &
                    &+ wnrm ** 2) / wnrh_init ** 2
              else
                fld_amp(ix, jy, kz, iwm) = 0.0
                wnrm = 0.0
              end if
              wnm_sfc(ix, jy, iwm) = wnrm

            end do
          end do
        end do
      else
        call stratification(z(1), 1, NN_nd)
        do jy = 1, ny
          do ix = 1, nx
            ! FJApr2023
            ! Loop over all wave modes.
            do iwm = 1, nwm
              ! FJApr2023
              wnrk_init = k_spectrum(ix, jy, iwm)
              wnrl_init = l_spectrum(ix, jy, iwm)
              wnrh_init = sqrt(wnrk_init ** 2.0 + wnrl_init ** 2.0)
              wnrm_init = 0.0

              ! FJApr2023
              ! omi_sfc(ix, jy) = - var(ix, jy, 1, 2) * wnrk_init &
              !     - var(ix, jy, 1, 3) * wnrl_init
              omi_sfc(ix, jy, iwm) = - 0.5 * (var%u(ix, jy, 1) + var%u(ix - 1, &
                  &jy, 1)) * wnrk_init - 0.5 * (var%v(ix, jy, 1) + var%v(ix, &
                  &jy - 1, 1)) * wnrl_init

              ! choose correct sign of horizontal wavenumbers in order to
              ! be on the correct frequency branch

              ! FJApr2023
              ! if(omi_sfc(ix, jy) * branchr >= 0.0) then
              !   wnk_sfc(ix, jy) = wnrk_init
              !   wnl_sfc(ix, jy) = wnrl_init
              ! else
              !   omi_sfc(ix, jy) = - omi_sfc(ix, jy)

              !   wnk_sfc(ix, jy) = - wnrk_init
              !   wnl_sfc(ix, jy) = - wnrl_init
              ! end if
              if(omi_sfc(ix, jy, iwm) * branchr >= 0.0) then
                wnk_sfc(ix, jy, iwm) = wnrk_init
                wnl_sfc(ix, jy, iwm) = wnrl_init
              else
                omi_sfc(ix, jy, iwm) = - omi_sfc(ix, jy, iwm)

                wnk_sfc(ix, jy, iwm) = - wnrk_init
                wnl_sfc(ix, jy, iwm) = - wnrl_init
              end if
              ! end do
              ! end do

              ! FJMar2023
              ! local squared Brunt-Vaisala frequency
              ! call stratification(z(0), 1, NN_nd)

              ! vertical wave number and wave-action density to be distributed
              ! over the ray volumes

              ! do jy = 1, ny
              ! do ix = 1, nx
              ! FJApr2023
              ! fld_amp(ix, jy, 0) = 0.0
              fld_amp(ix, jy, kz, iwm) = 0.0
              wnrm = 0.0

              ! FJJan2023
              ! if ((sigwpx == 0.0 .or. abs(x(ix + ix0) - xr0) < sigwpx) .and. &
              !     (sigwpy == 0.0 .or. abs(y(jy + jy0) - yr0) < sigwpy)) then
              ! FJApr2023
              ! if(abs(omi_sfc(ix, jy)) <= f_cor_nd) then
              !   fld_amp(ix, jy, 0) = 0.0
              !   wnrm = 0.0
              ! elseif(abs(omi_sfc(ix, jy)) < sqrt(NN_nd)) then
              !   wnrm = - branchr * sqrt(wnrh_init ** 2 * (NN_nd - omi_sfc(ix, jy) &
              !       ** 2) / (omi_sfc(ix, jy) ** 2 - f_cor_nd ** 2))
              if(abs(omi_sfc(ix, jy, iwm)) <= f_cor_nd) then
                fld_amp(ix, jy, kz, iwm) = 0.0
                wnrm = 0.0
              elseif(abs(omi_sfc(ix, jy, iwm)) < sqrt(NN_nd)) then
                wnrm = - branchr * sqrt(wnrh_init ** 2 * (NN_nd - omi_sfc(ix, &
                    &jy, iwm) ** 2) / (omi_sfc(ix, jy, iwm) ** 2 - f_cor_nd &
                    &** 2))

                ! Displacement (factor two accounts for other frequency branch)
                displm = abs(topography_spectrum(ix, jy, iwm))

                ! FJJan2023
                ! displacement
                ! displm = mountainHeight_wkb_dim / lRef
                ! if (sigwpx > 0.0) then
                !   displm = displm * 0.5 * (1.0 + cos(pi * (x(ix + ix0) - xr0) &
                !       / sigwpx))
                ! end if

                ! if (sigwpy > 0.0) then
                !   displm = displm * 0.5 * (1.0 + cos(pi * (y(jy + jy0) - yr0) &
                !       / sigwpy))
                ! end if

                ! Long number scaling
                if(blocking) then
                  ! Compute Long number.
                  long = sqrt(NN_nd / (0.25 * (var%u(ix, jy, 1) + var%u(ix &
                      &- 1, jy, 1)) ** 2.0 + 0.25 * (var%v(ix, jy, 1) &
                      &+ var%v(ix, jy - 1, 1)) ** 2.0)) &
                      &* sum(abs(topography_spectrum(ix, jy, :)))
                  ! Apply scaling.
                  displm = displm * wave_amplitude_reduction(long)
                end if

                ! surface wave-action density
                ! FJApr2023
                !   fld_amp(ix, jy, 0) = 0.5 * rhoStrat(0) * displm ** 2 &
                !       * omi_sfc(ix, jy) * (wnrh_init ** 2 + wnrm ** 2) &
                !       / wnrh_init ** 2
                ! else
                !   fld_amp(ix, jy, 0) = 0.0
                !   wnrm = 0.0
                ! end if
                fld_amp(ix, jy, kz, iwm) = 0.5 * rhoStrat(1) * displm ** 2 &
                    &* omi_sfc(ix, jy, iwm) * (wnrh_init ** 2 + wnrm ** 2) &
                    &/ wnrh_init ** 2
              else
                fld_amp(ix, jy, kz, iwm) = 0.0
                wnrm = 0.0
              end if
              ! end if

              ! FJApr2023
              ! wnm_sfc(ix, jy) = wnrm
              wnm_sfc(ix, jy, iwm) = wnrm
            end do
          end do
        end do
      end if
    else ! case_wkb /= 3
      if(topography) then

        if(case_wkb == 5) then
          print *, 'case_wkb == 5 + topography NOT implemented yet in setup_wkb'
          stop
        end if

        do kz = 1, sizeZ
          do jy = 1, ny
            do ix = 1, nx
              ! local squared Brunt-Vaisala frequency
              call stratification(zTFC(ix, jy, kz), 1, NN_nd)

              ! intrinsic frequency
              omi_notop(ix, jy, kz) = branchr * sqrt((NN_nd * wnrh_init ** 2 &
                  &+ f_cor_nd ** 2 * wnrm_init ** 2) / (wnrh_init ** 2 &
                  &+ wnrm_init ** 2))

              ! wave-action density
              fld_amp(ix, jy, kz, :) = (amp_wkb / wnrm_init) ** 2 * (wnrh_init &
                  &** 2 + wnrm_init ** 2) / (2.0 * wnrh_init ** 2) &
                  &* omi_notop(ix, jy, kz) * rhoStratTFC(ix, jy, kz)

              if(case_wkb == 1) then
                fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * exp(- &
                    &((zTFC(ix, jy, kz) - zr0) / sigwpz) ** 2)

                if(sigwpx_dim > 0.0) then
                  fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * exp(- &
                      &((x(ix + ix0) - xr0) / sigwpx) ** 2)
                end if

                if(sigwpy_dim > 0.0) then
                  fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * exp(- &
                      &((y(jy + jy0) - yr0) / sigwpy) ** 2)
                end if
              elseif(case_wkb == 2) then
                if(abs(zTFC(ix, jy, kz) - zr0) < sigwpz) then
                  fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * 0.5 * (1.0 &
                      &+ cos(pi * (zTFC(ix, jy, kz) - zr0) / sigwpz))

                  if(sigwpx > 0.0) then
                    if(abs(x(ix + ix0) - xr0) < sigwpx) then
                      fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * 0.5 &
                          &* (1.0 + cos(pi * (x(ix + ix0) - xr0) / sigwpx))
                    else
                      fld_amp(ix, jy, kz, :) = 0.0
                    end if
                  end if

                  if(sigwpy > 0.0) then
                    if(abs(y(jy + jy0) - yr0) < sigwpy) then
                      fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * 0.5 &
                          &* (1.0 + cos(pi * (y(jy + jy0) - yr0) / sigwpy))
                    else
                      fld_amp(ix, jy, kz, :) = 0.0
                    end if
                  end if
                else
                  fld_amp(ix, jy, kz, :) = 0.0
                end if
              elseif(case_wkb == 4) then
                fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * exp(- &
                    &(zTFC(ix, jy, kz) - zr0) ** 2. / sigwpz ** 2.)

                if(sigwpx > 0.0) then
                  if(abs(x(ix + ix0) - xr0) < sigwpx) then
                    fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * cos(pi &
                        &* (x(ix + ix0) - xr0) / (2. * sigwpx)) ** 2.
                  else
                    fld_amp(ix, jy, kz, :) = 0.0
                  end if
                end if

                if(sigwpy > 0.0) then
                  if(abs(y(jy + jy0) - yr0) < sigwpy) then
                    fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * cos(pi &
                        &* (y(jy + jy0) - yr0) / (2. * sigwpy)) ** 2.
                  else
                    fld_amp(ix, jy, kz, :) = 0.0
                  end if
                end if

              end if ! case_wkb
            end do ! ix
          end do ! jy
        end do ! kz
      else
        do kz = 1, sizeZ
          ! local squared Brunt-Vaisala frequency
          call stratification(z(kz), 1, NN_nd)
          ! wave-action density
          do jy = 1, ny
            do ix = 1, nx

              if(case_wkb == 5) then

                do iwm = 1, NWM_WP

                  branchr = branchr_sp(iwm)
                  amp_wkb = amp_wkb_sp(iwm)

                  wnrk_init = wnrk_init_sp(iwm)
                  wnrl_init = wnrl_init_sp(iwm)
                  wnrm_init = wnrm_init_sp(iwm)

                  wnrh_init = sqrt(wnrk_init ** 2.0 + wnrl_init ** 2.0)

                  sigwpz = sigwpz_sp(iwm)
                  sigwpx = sigwpx_sp(iwm)
                  sigwpy = sigwpy_sp(iwm)

                  zr0 = zr0_sp(iwm)
                  xr0 = xr0_sp(iwm)
                  yr0 = yr0_sp(iwm)

                  ! intrinsic frequency
                  omi_notop(ix, jy, kz) = branchr * sqrt((NN_nd * wnrh_init &
                      &** 2 + f_cor_nd ** 2 * wnrm_init ** 2) / (wnrh_init &
                      &** 2 + wnrm_init ** 2))

                  fld_amp(ix, jy, kz, iwm) = (amp_wkb / wnrm_init) ** 2 &
                      &* (wnrh_init ** 2 + wnrm_init ** 2) / (2.0 * wnrh_init &
                      &** 2) * omi_notop(ix, jy, kz) * rhoStrat(kz)

                  if(abs(z(kz) - zr0) < sigwpz) then
                    fld_amp(ix, jy, kz, iwm) = fld_amp(ix, jy, kz, iwm) * 0.5 &
                        &* (1.0 + cos(pi * (z(kz) - zr0) / sigwpz))

                    if(sigwpx > 0.0) then
                      if(abs(x(ix + ix0) - xr0) < sigwpx) then
                        fld_amp(ix, jy, kz, iwm) = fld_amp(ix, jy, kz, iwm) &
                            &* 0.5 * (1.0 + cos(pi * (x(ix + ix0) - xr0) &
                            &/ sigwpx))
                      else
                        fld_amp(ix, jy, kz, iwm) = 0.0
                      end if
                    end if

                    if(sigwpy > 0.0) then
                      if(abs(y(jy + jy0) - yr0) < sigwpy) then
                        fld_amp(ix, jy, kz, iwm) = fld_amp(ix, jy, kz, iwm) &
                            &* 0.5 * (1.0 + cos(pi * (y(jy + jy0) - yr0) &
                            &/ sigwpy))
                      else
                        fld_amp(ix, jy, kz, iwm) = 0.0
                      end if
                    end if

                  else
                    fld_amp(ix, jy, kz, iwm) = 0.0
                  end if

                end do ! iwm

              else ! case_wkb /= 5
                ! intrinsic frequency
                omi_notop(ix, jy, kz) = branchr * sqrt((NN_nd * wnrh_init ** 2 &
                    &+ f_cor_nd ** 2 * wnrm_init ** 2) / (wnrh_init ** 2 &
                    &+ wnrm_init ** 2))

                fld_amp(ix, jy, kz, :) = (amp_wkb / wnrm_init) ** 2 &
                    &* (wnrh_init ** 2 + wnrm_init ** 2) / (2.0 * wnrh_init &
                    &** 2) * omi_notop(ix, jy, kz) * rhoStrat(kz)

                !TEST****
                !var%OPT(ix, jy, kz, 3) = amp_wkb !(amp_wkb / wnrm_init) ** 2 &
                !* (wnrh_init ** 2 + wnrm_init ** 2) / (2.0 * wnrh_init &
                !** 2) !* omi_notop(ix, jy, kz) * rhoStrat(kz) !sum(fld_amp(ix, jy, kz, :))

                if(case_wkb == 1) then
                  fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * exp(- &
                      &((z(kz) - zr0) / sigwpz) ** 2)
                  if(compare_raytracer) then
                    if(sigwpx > 0.0) then
                      if(abs(x(ix + ix0) - xr0) < sigwpx) then
                        fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * 0.5 &
                            &* (1.0 + cos(pi * (x(ix + ix0) - xr0) / sigwpx))
                      else
                        fld_amp(ix, jy, kz, :) = 0.0
                      end if
                    end if

                    if(sigwpy > 0.0) then
                      if(abs(y(jy + jy0) - yr0) < sigwpy) then
                        fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * 0.5 &
                            &* (1.0 + cos(pi * (y(jy + jy0) - yr0) / sigwpy))
                      else
                        fld_amp(ix, jy, kz, :) = 0.0
                      end if
                    end if
                  else
                    if(sigwpx_dim > 0.0) then
                      fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * exp(- &
                          &((x(ix + ix0) - xr0) / sigwpx) ** 2)
                    end if

                    if(sigwpy_dim > 0.0) then
                      fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * exp(- &
                          &((y(jy + jy0) - yr0) / sigwpy) ** 2)
                    end if
                  end if
                elseif(case_wkb == 2) then

                  if(abs(z(kz) - zr0) < sigwpz) then
                    fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * 0.5 &
                        &* (1.0 + cos(pi * (z(kz) - zr0) / sigwpz))

                    if(sigwpx > 0.0) then
                      if(abs(x(ix + ix0) - xr0) < sigwpx) then
                        fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * 0.5 &
                            &* (1.0 + cos(pi * (x(ix + ix0) - xr0) / sigwpx))
                      else
                        fld_amp(ix, jy, kz, :) = 0.0
                      end if
                    end if

                    if(sigwpy > 0.0) then
                      if(abs(y(jy + jy0) - yr0) < sigwpy) then
                        fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * 0.5 &
                            &* (1.0 + cos(pi * (y(jy + jy0) - yr0) / sigwpy))
                      else
                        fld_amp(ix, jy, kz, :) = 0.0
                      end if
                    end if

                  else
                    fld_amp(ix, jy, kz, :) = 0.0
                  end if

                elseif(case_wkb == 4) then ! to match the wavepacket case gaussian
                  fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * exp(- &
                      &(z(kz) - zr0) ** 2. / sigwpz ** 2.)

                  if(sigwpx > 0.0) then
                    if(abs(x(ix + ix0) - xr0) < sigwpx) then
                      fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * cos(pi &
                          &* (x(ix + ix0) - xr0) / (2. * sigwpx)) ** 2.
                    else
                      fld_amp(ix, jy, kz, :) = 0.0
                    end if
                  end if

                  if(sigwpy > 0.0) then
                    if(abs(y(jy + jy0) - yr0) < sigwpy) then
                      fld_amp(ix, jy, kz, :) = fld_amp(ix, jy, kz, :) * cos(pi &
                          &* (y(jy + jy0) - yr0) / (2. * sigwpy)) ** 2.
                    else
                      fld_amp(ix, jy, kz, :) = 0.0
                    end if
                  end if
                end if ! case_wkb

              end if !case_wkb /=5
            end do ! ix
          end do ! jy
        end do ! kz

      end if ! topography
    end if ! case_wkb == 3

    ! initialize with wave-induced zonal wind
    ! currently only available for 1D
    if(lindUinit) then
      if((sigwpx > 0.0) .or. (sigwpy > 0.0) .or. (wlry_init /= 0.0)) then
        stop "setup_wkb: lindUinit only possible for 1D wavepacket envelope in &
            &z and wavelength in y = 0."
      end if

      do kz = 1, (nz)
        do im = 1, nwm
          var%u(:, :, kz) = var%u(:, :, kz) + wnrk_init * fld_amp(:, :, kz, &
              &im) / rhoStrat(kz)
        end do
      end do
    end if

    if((case_wkb .ne. 3) .and. include_tracer .and. (f_cor_nd == 0.)) then
      waveAmplitudes%phase = 0.
      do kz = 1, nz
        do jy = 1, ny
          do ix = 1, nx
            ! calculating bhat^(2) from the wave-action density
            waveAmplitudes(ix, jy, kz)%lowamp%b = cmplx(sqrt(2. * NN_nd ** 4. &
                &/ rhoStrat(kz) * wnrh_init ** 2. * fld_amp(ix, jy, kz, 1) &
                &/ (omi_notop(ix, jy, kz) * (wnrh_init ** 2. + wnrm_init &
                &** 2.))), 0.)
            ! calculating remaining leading-order wave amplitudes using the
            ! polarization relations
            waveAmplitudes(ix, jy, kz)%lowamp%u = cmplx(0., 1.) &
                &* waveAmplitudes(ix, jy, kz)%lowamp%b * wnrk_init &
                &* (omi_notop(ix, jy, kz) ** 2. - NN_nd ** 2.) / omi_notop(ix, &
                &jy, kz) / wnrm_init / NN_nd ** 2.
            waveAmplitudes(ix, jy, kz)%lowamp%v = cmplx(0., 1.) &
                &* waveAmplitudes(ix, jy, kz)%lowamp%b * wnrl_init &
                &* (omi_notop(ix, jy, kz) ** 2. - NN_nd ** 2.) / omi_notop(ix, &
                &jy, kz) / wnrm_init / NN_nd ** 2.

            waveAmplitudes(ix, jy, kz)%lowamp%w = cmplx(0., 1.) &
                &* waveAmplitudes(ix, jy, kz)%lowamp%b * omi_notop(ix, jy, kz) &
                &/ NN_nd ** 2.
            waveAmplitudes(ix, jy, kz)%lowamp%pi = cmplx(0., 1.) &
                &* waveAmplitudes(ix, jy, kz)%lowamp%b * (omi_notop(ix, jy, &
                &kz) ** 2. - NN_nd ** 2.) / NN_nd ** 2. / wnrm_init

            waveAmplitudes(ix, jy, kz)%lowamp%chi = 0.
            if(sizeX > 1) then
              waveAmplitudes(ix, jy, kz)%lowamp%chi = waveAmplitudes(ix, jy, &
                  &kz)%lowamp%u * (initialtracer(ix + 1, jy, kz) &
                  &- initialtracer(ix - 1, jy, kz)) / (2. * dx)
            end if
            if(sizeY > 1) then
              waveAmplitudes(ix, jy, kz)%lowamp%chi = waveAmplitudes(ix, jy, &
                  &kz)%lowamp%chi + waveAmplitudes(ix, jy, kz)%lowamp%v &
                  &* (initialtracer(ix, jy + 1, kz) - initialtracer(ix, jy &
                  &- 1, kz)) / (2. * dy)
            end if
            waveAmplitudes(ix, jy, kz)%lowamp%chi = waveAmplitudes(ix, jy, &
                &kz)%lowamp%chi + waveAmplitudes(ix, jy, kz)%lowamp%w &
                &* (initialtracer(ix, jy, kz + 1) - initialtracer(ix, jy, kz &
                &- 1)) / (2. * dz)
            waveAmplitudes(ix, jy, kz)%lowamp%chi = - waveAmplitudes(ix, jy, &
                &kz)%lowamp%chi * cmplx(0., 1.) / omi_notop(ix, jy, kz)
          end do
        end do
      end do
      call setBoundary_waveAmp(waveAmplitudes)
    else
      waveAmplitudes%lowamp%u = 0.
      waveAmplitudes%lowamp%v = 0.
      waveAmplitudes%lowamp%w = 0.
      waveAmplitudes%lowamp%b = 0.
      waveAmplitudes%lowamp%pi = 0.
      waveAmplitudes%lowamp%chi = 0.
    end if

    cgx_max = 0.0
    cgy_max = 0.0
    cgz_max = 0.0

    ! in mountain-wave case initialization of only one layer of ray
    ! volumes just below the bottom surface:

    ! if(case_wkb == 3) then
    !   kz2min = nrzl
    ! else
    !   kz2min = 1
    ! end if

    ! nondimensional wave-number widths to be filled by ray volumes

    dk_ini_nd = dk_init * lRef
    dl_ini_nd = dl_init * lRef
    dm_ini_nd = dm_init * lRef

    if(ixmin <= ixmax .and. jymin <= jymax) then
      !  in x-k subspace, loop over all spatial cells with ray volumes
      do ix = ixmin, ixmax
        ! likewise for y-l subspace
        do jy = jymin, jymax
          ! likewise for z-m subspace
          do kz = kzmin, kzmax
            iRay = 0
            i_sfc = 0

            ! FJApr2023
            ! if(topography .and. (zTFC(ix, jy, kz) < zrmin &
            !     .or. zTFC(ix, jy, kz) > zrmax)) cycle

            ! in x-k subspace, loop over all r.v. within one spatial cell
            do ix2 = 1, nrxl
              ! in x-k subspace, loop over all r.v. within the
              ! wave-number extent to be filled with r.v.
              do ik = 1, nrk_init

                ! likewise for y-l subspace
                do jy2 = 1, nryl
                  do jl = 1, nrl_init

                    ! likewise for z-m subspace
                    ! do kz2 = kz2min, nrzl
                    do kz2 = 1, nrzl
                      do km = 1, nrm_init

                        ! FJApr2023
                        ! Loop over all wave modes.
                        do iwm = 1, nwm
                          if(case_wkb == 3) then
                            ! pointers for surface ray volumes

                            i_sfc = i_sfc + 1

                            ix2_sfc(i_sfc) = ix2
                            jy2_sfc(i_sfc) = jy2
                            kz2_sfc(i_sfc) = kz2

                            ik_sfc(i_sfc) = ik
                            jl_sfc(i_sfc) = jl
                            km_sfc(i_sfc) = km

                            iwm_sfc(i_sfc) = iwm

                            ! only add ray volumes with non-zero
                            ! wave-action density
                            ! (thus excluding intrinsic frequencies
                            ! outside of the allowed range)
                            ! excluded cases indicated by negative
                            ! ray-volume index

                            if(fld_amp(ix, jy, kz, iwm) == 0.0) then
                              ir_sfc(i_sfc, ix, jy) = - 1
                              cycle
                            else
                              iRay = iRay + 1
                              ir_sfc(i_sfc, ix, jy) = iRay
                            end if
                          else
                            iRay = iRay + 1
                          endif

                          ! ray-volume positions

                          ray(iRay, ix, jy, kz)%x = (x(ix + ix0) - 0.5 * dx &
                              &+ (ix2 - 0.5) * dx / nrxl)

                          ray(iRay, ix, jy, kz)%y = (y(jy + jy0) - 0.5 * dy &
                              &+ (jy2 - 0.5) * dy / nryl)

                          if(topography) then
                            ! FJApr2023
                            ray(iRay, ix, jy, kz)%z = (zTFC(ix, jy, kz) - 0.5 &
                                &* jac(ix, jy, kz) * dz + (kz2 - 0.5) &
                                &* jac(ix, jy, kz) * dz / nrzl)
                          else
                            ray(iRay, ix, jy, kz)%z = (z(kz) - 0.5 * dz + (kz2 &
                                &- 0.5) * dz / nrzl)
                          end if

                          xr = ray(iRay, ix, jy, kz)%x
                          yr = ray(iRay, ix, jy, kz)%y
                          zr = ray(iRay, ix, jy, kz)%z

                          ! local squared Brunt_Vaisala frequency

                          if(zr < lz(0) - dz) then
                            print *, 'ERROR IN setup_wkb: RAY VOLUME', iRay, &
                                &'at', ix, jy, kz, 'TOO LOW'
                            stop
                          end if

                          call stratification(zr, 1, NNr)

                          ! ray-volume spatial extensions

                          ray(iRay, ix, jy, kz)%dxray = dx / nrxl
                          ray(iRay, ix, jy, kz)%dyray = dy / nryl
                          if(topography) then
                            ! FJApr2023
                            ray(iRay, ix, jy, kz)%dzray = jac(ix, jy, kz) * dz &
                                &/ nrzl
                          else
                            ray(iRay, ix, jy, kz)%dzray = dz / nrzl
                          end if

                          ! ray-volume wave numbers

                          if(case_wkb == 3) then
                            wnk_0 = wnk_sfc(ix, jy, iwm)
                            wnl_0 = wnl_sfc(ix, jy, iwm)
                            wnm_0 = wnm_sfc(ix, jy, iwm)
                          elseif(case_wkb == 5) then
                            wnk_0 = wnrk_init_sp(iwm)
                            wnl_0 = wnrl_init_sp(iwm)
                            wnm_0 = wnrm_init_sp(iwm)
                          else
                            wnk_0 = wnrk_init
                            wnl_0 = wnrl_init
                            wnm_0 = wnrm_init
                          end if

                          ! FJApr2023
                          ! Ensure correct wavenumber extents.
                          if(case_wkb == 3 .and. sizeX > 1) then
                            dk_ini_nd = fac_dk_init * sqrt(wnk_0 ** 2.0 &
                                &+ wnl_0 ** 2.0)
                          end if

                          ray(iRay, ix, jy, kz)%k = (wnk_0 - 0.5 * dk_ini_nd &
                              &+ (real(ik) - 0.5) * dk_ini_nd / nrk_init)

                          ! FJApr2023
                          ! Ensure correct wavenumber extents.
                          if(case_wkb == 3 .and. sizeY > 1) then
                            dl_ini_nd = fac_dl_init * sqrt(wnk_0 ** 2.0 &
                                &+ wnl_0 ** 2.0)
                          end if

                          ray(iRay, ix, jy, kz)%l = (wnl_0 - 0.5 * dl_ini_nd &
                              &+ (real(jl) - 0.5) * dl_ini_nd / nrl_init)

                          if(wnm_0 == 0.0) then
                            stop "Error in setup_wkb: wnm_0 = 0!"
                          else
                            dm_ini_nd = fac_dm_init * abs(wnm_0)
                          end if

                          ray(iRay, ix, jy, kz)%m = (wnm_0 - 0.5 * dm_ini_nd &
                              &+ (real(km) - 0.5) * dm_ini_nd / nrm_init)

                          ! ray-volume wave-number extents

                          ray(iRay, ix, jy, kz)%dkray = dk_ini_nd / nrk_init

                          ray(iRay, ix, jy, kz)%dlray = dl_ini_nd / nrl_init

                          ray(iRay, ix, jy, kz)%dmray = dm_ini_nd / nrm_init

                          ! ray-volume phase-space volume

                          ray(iRay, ix, jy, kz)%area_xk = ray(iRay, ix, jy, &
                              &kz)%dxray * ray(iRay, ix, jy, kz)%dkray

                          ray(iRay, ix, jy, kz)%area_yl = ray(iRay, ix, jy, &
                              &kz)%dyray * ray(iRay, ix, jy, kz)%dlray

                          ray(iRay, ix, jy, kz)%area_zm = ray(iRay, ix, jy, &
                              &kz)%dzray * ray(iRay, ix, jy, kz)%dmray

                          pspvol = dm_ini_nd

                          if(sizeX > 1) then
                            pspvol = pspvol * dk_ini_nd
                          end if

                          if(sizeY > 1) then
                            pspvol = pspvol * dl_ini_nd
                          end if

                          ! phase-space wave-action density

                          if(kz == sizeZ) then
                            ray(iRay, ix, jy, kz)%dens = 0.0
                          else
                            ray(iRay, ix, jy, kz)%dens = fld_amp(ix, jy, kz, &
                                &iwm) / pspvol
                          endif

                          ! intrinsic frequency

                          if(case_wkb == 3) then
                            ray(iRay, ix, jy, kz)%omega = omi_sfc(ix, jy, iwm)
                          else
                            ray(iRay, ix, jy, kz)%omega = omi_notop(ix, jy, kz)
                          end if

                          ! intrinsic group velocities and maximum
                          ! group velocities

                          call meanflow(xr, yr, zr, var, 1, uxr)
                          call meanflow(xr, yr, zr, var, 2, vyr)
                          call meanflow(xr, yr, zr, var, 3, wzr)

                          wnrk = ray(iRay, ix, jy, kz)%k
                          wnrl = ray(iRay, ix, jy, kz)%l
                          wnrm = ray(iRay, ix, jy, kz)%m

                          wnrh = sqrt(wnrk ** 2 + wnrl ** 2)

                          omir = ray(iRay, ix, jy, kz)%omega

                          cgirx = wnrk * (NNr - omir ** 2) / (omir * (wnrh &
                              &** 2 + wnrm ** 2))

                          if(abs(uxr + cgirx) > abs(cgx_max)) then
                            cgx_max = abs(uxr + cgirx)
                          end if

                          cgiry = wnrl * (NNr - omir ** 2) / (omir * (wnrh &
                              &** 2 + wnrm ** 2))

                          if(abs(vyr + cgiry) > abs(cgy_max)) then
                            cgy_max = abs(vyr + cgiry)
                          end if

                          cgirz = - wnrm * (omir ** 2 - f_cor_nd ** 2) / (omir &
                              &* (wnrh ** 2 + wnrm ** 2))

                          if(abs(wzr + cgirz) > abs(cgz_max)) then
                            cgz_max = abs(wzr + cgirz)
                          end if

                          ! SD
                          if(include_ice) then
                            dphi = wnrk * xr + wnrl * yr + wnrm * zr
                            ray(iRay, ix, jy, kz)%dphi = dphi
                          end if

                        end do ! iwm
                      end do ! km
                    end do ! kz2
                  end do ! jl
                end do ! jy2

              end do ! ik
            end do ! ix2

            if(iRay > nray_wrk) then
              print *, 'ERROR at ix,jy,kz =', ix, jy, kz
              print *, 'iRay =', iRay, '> nray_wrk =', nray_wrk
              stop
            end if

            nRay(ix, jy, kz) = iRay

            if(nRay(ix, jy, kz) > nray_wrk) then
              print *, 'ERROR at ix,jy,kz =', ix, jy, kz
              print *, 'nRay =', nRay(ix, jy, kz), '> nray_wrk =', nray_wrk
              stop
            end if

            if(case_wkb == 3) then
              if(i_sfc /= n_sfc) then
                print *, 'ERROR at ix,jy,kz =', ix, jy, kz
                print *, 'i_sfc =', i_sfc, '/= n_sfc =', n_sfc
                stop
              end if
            end if
          end do ! kz
        end do ! jy
      end do ! ix

    end if

    !SD
    !double check if output rays needed
    if(include_ice) then

      allocate(nor_mst(nprocx * nprocy), stat = allocstat)
      if(allocstat /= 0) stop "wkb.f90/setup_wkb: could not allocate nor_mst. &
          &Stop."
      allocate(vct_prc(nray_wrk * nx * ny * nz, NFR), stat = allocstat)
      if(allocstat /= 0) stop "wkb.f90/setup_wkb: could not allocate vct_prc. &
          &Stop."
      allocate(vct_mst(nray_wrk * sizeX * sizeY * sizeZ, NFR), stat = allocstat)
      if(allocstat /= 0) stop "wkb.f90/setup_wkb: could not allocate vct_mst. &
          &Stop."
      !allocate(vct_out(nray_wrk*sizeX*sizeY*sizeZ, NFR), stat = allocstat)
      !if(allocstat /= 0) stop "wkb.f90/setup_wkb: could not allocate &
      !     vec_out. Stop."

      !Saturation field ect.
      allocate(opt_ray(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), stat &
          &= allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate opt_ray"

      ! init values
      nor_mst = 0.
      vct_prc = 0.
      vct_mst = 0.
      NoR_out = 1.

      if(case_wkb == 3) stop 'RayTracer+Ice+Topography not supported'
      if(topography) stop 'RayTracer+Ice+Topography not supported'

    end if
    !-------------------------------
    !       Feedback on Screen
    !-------------------------------

    nrsuml = sum(nRay(1:nx, 1:ny, 0:nz))

    call mpi_reduce(nrsuml, nr_sum, 1, mpi_integer, mpi_sum, root, comm, ierror)
    call mpi_bcast(nr_sum, 1, mpi_integer, root, comm, ierror)

    if(master) then
      ! Feedback to user
      print *, ""
      print *, " 11) Ray tracer: "

      if(rayTracer) then
        write(*, fmt = "(a25,a)") "ray tracer = ", "on"

        write(*, fmt = "(a25,i7)") "nb of rays = ", nr_sum
        write(*, fmt = "(a25,f7.1)") "nb rays / grid cell  = ", real(nr_sum) &
            &/ (sizeX * sizeY * (sizeZ + 1))
        write(*, fmt = "(a31,i7)") "max allowed nb of rays / cell = ", nray_max
      else
        write(*, fmt = "(a25,a)") "ray tracer = ", "off"
      end if
    end if

  end subroutine setup_wkb

  !------------------------------------------------------------------------

  subroutine stratification(zlc, strtpe, str)

    !---------------------------------------------------------------
    ! interpolation squared Brunt-Vasisala frequency or its vertical
    ! derivative to a specified vertical position
    ! strtpe = 1: N^2
    !          2: dN^2/dz
    !---------------------------------------------------------------

    real, intent(in) :: zlc
    integer, intent(in) :: strtpe
    real, intent(out) :: str

    integer :: kzu, kzd

    real :: zu, zd
    real :: strd, stru

    real :: factor

    ! interpolate NN = g/theta_b * d theta_b/dz and dNN/dz
    ! using grid distribution
    !
    !       theta_b, NN      z = z(k+1)
    !
    ! ------- dNN/dz ------- z = z(k+1/2) = lz(0) +     dz + (k-1)*dz
    !
    !       theta_b, NN      z = z(k)     = lz(0) + 0.5*dz + (k-1)*dz
    !
    ! ------- dNN/dz ------- z = z(k-1/2)   z |
    !                                         |
    !       theta_b, NN      z = z(k-1)       |
    !

    if(strtpe == 1) then
      if(topography) then
        ! Locate the two closest levels.
        kzd = max(- 1, floor((levelTFC(0, 0, zlc) - 0.5 * dz - lz(0)) / dz) + 1)
        kzu = kzd + 1

        if(kzu > nz + 1) then
          kzu = nz + 1
          kzd = nz
        end if

        ! Assign the values.
        zd = zTFC(0, 0, kzd)
        zu = zTFC(0, 0, kzu)
        strd = bvsStratTFC(0, 0, kzd)
        stru = bvsStratTFC(0, 0, kzu)
      else
        kzd = max(- 1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 1)

        kzu = kzd + 1

        if(kzu > nz + 1) then
          kzu = nz + 1
          kzd = nz
        end if

        zu = z(kzu)
        zd = z(kzd)

        strd = bvsStrat(kzd)
        stru = bvsStrat(kzu)
      end if
    elseif(strtpe == 2) then
      if(topography) then
        ! Locate the two closest levels.
        kzd = max(- 1, floor((levelTFC(0, 0, zlc) - lz(0)) / dz))
        kzu = kzd + 1

        if(kzu + 1 > nz + 1) then
          kzu = nz
          kzd = nz - 1
        end if

        ! Assign the values.
        zd = zTildeTFC(0, 0, kzd)
        zu = zTildeTFC(0, 0, kzu)
        strd = (bvsStratTFC(0, 0, kzd + 1) - bvsStratTFC(0, 0, kzd)) * 2.0 &
            &/ (jac(0, 0, kzd) + jac(0, 0, kzd + 1)) / dz
        stru = (bvsStratTFC(0, 0, kzu + 1) - bvsStratTFC(0, 0, kzu)) * 2.0 &
            &/ (jac(0, 0, kzu) + jac(0, 0, kzu + 1)) / dz
      else
        kzd = max(- 1, floor((zlc - lz(0)) / dz))

        kzu = kzd + 1

        if(kzu + 1 > nz + 1) then
          kzu = nz
          kzd = nz - 1
        end if

        zu = z(kzu) + 0.5 * dz
        zd = z(kzd) + 0.5 * dz

        strd = (bvsStrat(kzd + 1) - bvsStrat(kzd)) / dz
        stru = (bvsStrat(kzu + 1) - bvsStrat(kzu)) / dz
      end if
    else
      print *, 'ERROR: UNKNOWN strtpe =', strtpe, 'IN STRATIFICATION'
      stop
    end if

    ! interpolation in z

    if(zu < zd) then
      print *, 'ERROR IN STRATIFICATION: zu =', zu, '< zd =', zd
      stop
    elseif(zu == zd) then
      factor = 0.0
    elseif(zlc > zu) then
      factor = 0.0
    elseif(zlc > zd) then
      factor = (zu - zlc) / (zu - zd)
    else
      factor = 1.0
    end if

    str = factor * strd + (1.0 - factor) * stru

    ! testb
    if(strtpe == 1) then
      if(str <= 0.0) then
        print *, 'NN =', str, '<= 0.0 at zlc =', zlc
        print *, 'kzu =', kzu
        print *, 'kzd =', kzd
        print *, 'zu =', zu
        print *, 'dz =', dz
        print *, 'factor =', factor
        print *, 'strd =', strd
        print *, 'stru =', stru
        stop
      end if
    end if
    ! teste

  end subroutine stratification

  !------------------------------------------------------------------------
  subroutine tracerderivative(position, direction, position2, position3, var, &
      &dchidxyz)

    ! calculate the large-scale tracer derivative wrt x, y, and z

    real, intent(in) :: position ! at which the derivative of chi should be calculated
    integer, intent(in) :: direction ! 1: x, 2: y, 3: z
    real, intent(in) :: position2, position3 ! indices for the remaining two directions
    type(var_type), intent(in) :: var
    real, intent(out) :: dchidxyz

    integer :: kzu, kzd, jyf, jyb, ixr, ixl
    real :: xlc, ylc, zlc
    real :: zu, zd, yf, yb, xr, xl
    real :: tracu, tracd, tracl, tracr, tracf, tracb
    real :: rhodp, rhodm, rhoup, rhoum
    real :: rhofp, rhofm, rhobp, rhobm
    real :: rholp, rholm, rhorp, rhorm
    integer :: ixx, jyy, kzz

    real :: factor

    if(direction == 3) then
      xlc = position2
      ylc = position3

      ! find closest index on x-axis
      ixl = max(- 1, floor((xlc - lx(0)) / dx))
      ixr = ixl + 1

      if(ixr + 1 > nx + 1) then
        ixr = nx
        ixl = nx - 1
      end if

      if(abs(x(ixr) - xlc) > abs(xlc - x(ixl))) then
        ixx = ixl
      else
        ixx = ixr
      end if

      ! find closest index on y-axis
      jyb = max(- 1, floor((ylc - ly(0)) / dy))
      jyf = jyb + 1

      if(jyf + 1 > ny + 1) then
        jyf = ny
        jyb = ny - 1
      end if

      if(abs(y(jyf) - ylc) > abs(ylc - y(jyb))) then
        jyy = jyb
      else
        jyy = jyf
      end if

      zlc = position

      kzd = max(- 1, floor((zlc - lz(0)) / dz))
      kzu = kzd + 1

      if(kzu + 1 > nz + 1) then
        kzu = nz
        kzd = nz - 1
      end if

      zu = z(kzu) + 0.5 * dz
      zd = z(kzd) + 0.5 * dz

      rhodp = var%rho(ixx, jyy, kzd + 1) + rhoStrat(kzd + 1)
      rhodm = var%rho(ixx, jyy, kzd) + rhoStrat(kzd)

      rhoup = var%rho(ixx, jyy, kzu + 1) + rhoStrat(kzu + 1)
      rhoum = var%rho(ixx, jyy, kzu) + rhoStrat(kzu)

      tracd = (var%chi(ixx, jyy, kzd + 1) / rhodp - var%chi(ixx, jyy, kzd) &
          &/ rhodm) / dz
      tracu = (var%chi(ixx, jyy, kzu + 1) / rhoup - var%chi(ixx, jyy, kzu) &
          &/ rhoum) / dz

      if(zu < zd) then
        print *, 'ERROR IN TRACERDERIVATIVE: zu =', zu, '< zd =', zd
        stop
      elseif(zu == zd) then
        factor = 0.0
      elseif(zlc > zu) then
        factor = 0.0
      elseif(zlc > zd) then
        factor = (zu - zlc) / dz
      else
        factor = 1.0
      end if

      dchidxyz = factor * tracd + (1.0 - factor) * tracu

    elseif(direction == 2) then
      if(sizeY == 1) then
        dchidxyz = 0.0
      else
        xlc = position2
        ylc = position
        zlc = position3

        ! find closest index on x-axis
        ixl = max(- 1, floor((xlc - lx(0)) / dx))
        ixr = ixl + 1

        if(ixr + 1 > nx + 1) then
          ixr = nx
          ixl = nx - 1
        end if

        if(abs(x(ixr) - xlc) > abs(xlc - x(ixl))) then
          ixx = ixl
        else
          ixx = ixr
        end if

        ! find closest index on z-axis
        kzd = max(- 1, floor((zlc - lz(0)) / dz))
        kzu = kzd + 1

        if(kzu + 1 > nz + 1) then
          kzu = nz
          kzd = nz - 1
        end if

        if(abs(z(kzu) - zlc) > abs(zlc - z(kzd))) then
          kzz = kzd
        else
          kzz = kzu
        end if

        jyb = max(- 1, floor((ylc - ly(0)) / dy))
        jyf = jyb + 1

        if(jyf + 1 > ny + 1) then
          jyf = ny
          jyb = ny - 1
        end if

        yf = y(jyf) + 0.5 * dy
        yb = y(jyb) + 0.5 * dy

        rhobp = var%rho(ixx, jyb + 1, kzz) + rhoStrat(kzz)
        rhobm = var%rho(ixx, jyb, kzz) + rhoStrat(kzz)

        rhofp = var%rho(ixx, jyf + 1, kzz) + rhoStrat(kzz)
        rhofm = var%rho(ixx, jyf, kzz) + rhoStrat(kzz)

        tracb = (var%chi(ixx, jyb + 1, kzz) / rhobp - var%chi(ixx, jyb, kzz) &
            &/ rhobm) / dy
        tracf = (var%chi(ixx, jyf + 1, kzz) / rhofp - var%chi(ixx, jyf, kzz) &
            &/ rhofm) / dy

        if(yf < yb) then
          print *, 'ERROR IN TRACERDERIVATIVE: yf =', yf, '< yb =', yb
          stop
        elseif(yf == yb) then
          factor = 0.0
        elseif(ylc > yf) then
          factor = 0.0
        elseif(ylc > yb) then
          factor = (yf - ylc) / dy
        else
          factor = 1.0
        end if

        dchidxyz = factor * tracb + (1.0 - factor) * tracf
      end if

    elseif(direction == 1) then
      if(sizeX == 1) then
        dchidxyz = 0.0
      else
        xlc = position
        ylc = position2
        zlc = position3

        ! find closest index on y-axis
        jyb = max(- 1, floor((ylc - ly(0)) / dy))
        jyf = jyb + 1

        if(jyf + 1 > ny + 1) then
          jyf = ny
          jyb = ny - 1
        end if

        if(abs(y(jyf) - ylc) > abs(ylc - y(jyb))) then
          jyy = jyb
        else
          jyy = jyf
        end if

        ! find closest index on z-axis
        kzd = max(- 1, floor((zlc - lz(0)) / dz))
        kzu = kzd + 1

        if(kzu + 1 > nz + 1) then
          kzu = nz
          kzd = nz - 1
        end if

        if(abs(z(kzu) - zlc) > abs(zlc - z(kzd))) then
          kzz = kzd
        else
          kzz = kzu
        end if

        ixl = max(- 1, floor((xlc - lx(0)) / dx))
        ixr = ixl + 1

        if(ixr + 1 > nx + 1) then
          ixr = nx
          ixl = nx - 1
        end if

        xr = x(ixr) + 0.5 * dx
        xl = x(ixl) + 0.5 * dx

        rholp = var%rho(ixl + 1, jyy, kzz) + rhoStrat(kzz)
        rholm = var%rho(ixl, jyy, kzz) + rhoStrat(kzz)

        rhorp = var%rho(ixr + 1, jyy, kzz) + rhoStrat(kzz)
        rhorm = var%rho(ixr, jyy, kzz) + rhoStrat(kzz)

        tracl = (var%chi(ixl + 1, jyy, kzz) / rholp - var%chi(ixl, jyy, kzz) &
            &/ rholm) / dx
        tracr = (var%chi(ixr + 1, jyy, kzz) / rhorp - var%chi(ixr, jyy, kzz) &
            &/ rhorm) / dx

        if(xr < xl) then
          print *, 'ERROR IN TRACERDERIVATIVE: xr =', xr, '< xl =', xl
          stop
        elseif(xr == xl) then
          factor = 0.0
        elseif(xlc > xr) then
          factor = 0.0
        elseif(xlc > xl) then
          factor = (xr - xlc) / dx
        else
          factor = 1.0
        end if

        dchidxyz = factor * tracl + (1.0 - factor) * tracr
      end if
    else
      print *, "ERROR IN TRACERDERIVATIVE: direction must be 1, 2, or 3."
      stop
    end if

  end subroutine tracerderivative
  !------------------------------------------------------------------------

  subroutine meanflow(xlc_in, ylc_in, zlc_in, var, flwtpe, flw)

    !----------------------------------------------------------------
    ! interpolation mean flow or its spatial derivatives to specified
    ! local position
    ! flwtpe = 1 : u
    !          11: dudx
    !          12: dudy
    !          13: dudz
    !          2 : v
    !          21: dvdx
    !          22: dvdy
    !          23: dvdz
    !          3 : w
    !          31: dwdx
    !          32: dwdy
    !          33: dwdz
    !----------------------------------------------------------------

    real, intent(in) :: xlc_in, ylc_in, zlc_in
    type(var_type), intent(in) :: var
    integer, intent(in) :: flwtpe
    real, intent(out) :: flw

    integer :: kzu, kzd, jyf, jyb, ixr, ixl

    ! FJApr2023
    integer :: kzlbd, kzlbu, kzlfd, kzlfu, kzrbd, kzrbu, kzrfd, kzrfu
    real :: zlbd, zlbu, zlfd, zlfu, zrbd, zrbu, zrfd, zrfu
    real :: zbd, zbu, zfd, zfu

    real :: xlc, ylc, zlc
    real :: zu, zd, yf, yb, xr, xl
    real :: flwlbd, flwlbu, flwlfd, flwlfu, flwrbd, flwrbu, flwrfd, flwrfu
    real :: flwbd, flwbu, flwfd, flwfu
    real :: flwd, flwu

    real :: factor

    integer :: ix0, jy0

    ix0 = is + nbx - 1
    jy0 = js + nby - 1

    xlc = xlc_in
    ylc = ylc_in
    zlc = zlc_in

    ! vertical boundary conditions:
    ! zBoundary = 'periodic': implement periodicity
    ! zBoundary = 'solid_wall': noting done

    if(zlc < lz(0)) then
      select case(zBoundary)
      case("periodic")
        zlc = lz(1) + mod(zlc - lz(0), lz(1) - lz(0))
      case("solid_wall")
      case default
        stop "calc_meanflow_effect: unknown case zBoundary"
      end select
    elseif(zlc > lz(1)) then
      select case(zBoundary)
      case("periodic")
        zlc = lz(0) + mod(zlc - lz(1), lz(1) - lz(0))
      case("solid_wall")
      case default
        stop "calc_meanflow_effect: unknown case zBoundary"
      end select
    end if

    if(flwtpe == 1) then
      ! interpolate u using staggered-grid distribution

      !       |               |
      ! u(i-1,.,k+1)      u(i,.,k+1)  z = z(k+1)
      !       |               |
      !       |               |
      ! ----------------------------  z = z(k+1/2) = z(k) + 0.5*dz
      !       |               |         = lz(0) +  dz + (k-1)*dz
      !       |               |
      !   u(i-1,.,k)      u(i,.,k)    z = z(k)
      !       |               |         = lz(0) +  0.5*dz + (k-1)*dz
      !       |               |
      ! ----------------------------  z = z(k-1/2)
      !       |               |
      !       |               |
      ! u(i-1,.,k-1)      u(i,.,k-1)  z = z(k-1)
      !       |               |
      !    x(i-1/2)  x(i)  x(i+1/2)
      !
      !    x(i)     = lx(0) + 0.5*dx + (i-1)*dx
      !    x(i+1/2) = lx(0) +     dx + (i-1)*dx
      !
      !       |               |
      ! u(i-1,j+1,.)      u(i,j+1,.)  y = y(j+1)
      !       |               |
      !       |               |
      ! ----------------------------  y = y(j+1/2) = y(j) + 0.5*dy
      !       |               |         = ly(0) +  dy + (j-1)*dy
      !       |               |
      !   u(i-1,j,.)      u(i,j,.)    y = y(j)
      !       |               |         = ly(0) +  0.5*dy + (j-1)*dy
      !       |               |
      ! ----------------------------  y = y(j-1/2)
      !       |               |
      !       |               |
      ! u(i-1,j-1,.)      u(i,j-1,.)  y = y(j-1)
      !       |               |
      !    x(i-1/2)  x(i)  x(i+1/2)
      !

      ! for variable at full levels z(k) = lz(0) +  0.5*dz + (k-1)*dz:
      ! (levels below the model bottom are replaced by the first level
      ! in the model domain)

      if(topography) then
        ! FJApr2023
        ! Locate the closest points in zonal direction.
        if(sizeX == 1) then
          ixl = 1
          ixr = 1
        else
          ixl = floor((xlc - lx(0)) / dx) - ix0
          if(ixl < - nbx) then
            print *, "ERROR IN MEANFLOW: ixl =", ixl, "< - nbx =", - nbx
            stop
          end if
          ixr = ixl + 1
          if(ixr > nx + nbx) then
            print *, "ERROR IN MEANFLOW: ixr =", ixr, "> nx + nbx =", nx + nbx
            stop
          end if
        end if
        xr = x(ixr + ix0) + 0.5 * dx
        xl = x(ixl + ix0) + 0.5 * dx

        ! Locate the closest points in meridional direction.
        if(sizeY == 1) then
          jyb = 1
          jyf = 1
        else
          jyb = floor((ylc - 0.5 * dy - ly(0)) / dy) + 1 - jy0
          if(jyb < - nby) then
            print *, "ERROR IN MEANFLOW: jyb =", jyb, "< - nby =", - nby
            stop
          end if
          jyf = jyb + 1
          if(jyf > ny + nby) then
            print *, "ERROR IN MEANFLOW: jyf =", jyf, "> ny + nby =", ny + nby
            stop
          end if
        end if
        yf = y(jyf + jy0)
        yb = y(jyb + jy0)

        ! Locate the closest points in vertical direction.

        kzlbd = max(1, floor((levelTFC(ixl, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzlbu = max(1, floor((levelTFC(ixl, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzlbd > nz) then
          kzlbu = nz
          kzlbd = nz
        end if
        zlbd = zTFC(ixl, jyb, kzlbd)
        zlbu = zTFC(ixl, jyb, kzlbu)

        kzlfd = max(1, floor((levelTFC(ixl, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzlfu = max(1, floor((levelTFC(ixl, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzlfd > nz) then
          kzlfu = nz
          kzlfd = nz
        end if
        zlfd = zTFC(ixl, jyf, kzlfd)
        zlfu = zTFC(ixl, jyf, kzlfu)

        kzrbd = max(1, floor((levelTFC(ixr, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzrbu = max(1, floor((levelTFC(ixr, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzrbd > nz) then
          kzrbu = nz
          kzrbd = nz
        end if
        zrbd = zTFC(ixr, jyb, kzrbd)
        zrbu = zTFC(ixr, jyb, kzrbu)

        kzrfd = max(1, floor((levelTFC(ixr, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzrfu = max(1, floor((levelTFC(ixr, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzrfd > nz) then
          kzrfu = nz
          kzrfd = nz
        end if
        zrfd = zTFC(ixr, jyf, kzrfd)
        zrfu = zTFC(ixr, jyf, kzrfu)

        ! Assign the values.

        flwlbd = var%u(ixl, jyb, kzlbd)
        flwlbu = var%u(ixl, jyb, kzlbu)

        flwlfd = var%u(ixl, jyf, kzlfd)
        flwlfu = var%u(ixl, jyf, kzlfu)

        flwrbd = var%u(ixr, jyb, kzrbd)
        flwrbu = var%u(ixr, jyb, kzrbu)

        flwrfd = var%u(ixr, jyf, kzrfd)
        flwrfu = var%u(ixr, jyf, kzrfu)
      else
        ! index for lower level used for linear interpolation
        kzd = max(1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 1)

        ! index for upper level used for linear interpolation
        kzu = max(1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 2)

        if(kzd > nz) then
          kzu = nz
          kzd = nz
        end if

        ! full levels z(k)
        zu = z(kzu)
        zd = z(kzd)

        if(sizeY == 1) then
          jyb = 1
          jyf = 1
        else
          ! for variable at full levels y(j) = ly(0) +  0.5*dy + (j-1)*dy:

          ! index for backward level used for linear interpolation
          jyb = floor((ylc - 0.5 * dy - ly(0)) / dy) + 1 - jy0

          if(jyb < - nby) then
            print *, 'ERROR IN MEANFLOW: jyb =', jyb, '< -nby =', - nby
            stop
          end if

          ! index for forward level used for linear interpolation
          jyf = jyb + 1

          if(jyf > ny + nby) then
            print *, 'ERROR IN MEANFLOW: jyf =', jyf, '> ny + nby =', ny + nby
            stop
          end if
        end if

        ! full levels y(j)
        yf = y(jyf + jy0)
        yb = y(jyb + jy0)

        if(sizeX == 1) then
          ixl = 1
          ixr = 1
        else
          ! for variable at intermediate levels
          ! x(i+1/2) = x(i) + 0.5*dx = lx(0) + i*dx:

          ! index for leftmost level used for linear interpolation
          ixl = floor((xlc - lx(0)) / dx) - ix0

          if(ixl < - nbx) then
            print *, 'ERROR IN MEANFLOW: ixl =', ixl, '< -nbx =', - nbx
            ! testb
            print *, 'variable = u'
            print *, 'xlc =', xlc
            print *, 'lx(0) =', lx(0)
            print *, 'lx(1) =', lx(1)
            print *, 'dx =', dx
            print *, 'lRef =', lRef
            ! teste
            stop
          end if

          ! index for rightmost level used for linear interpolation
          ixr = ixl + 1

          if(ixr > nx + nbx) then
            print *, 'ERROR IN MEANFLOW: ixr =', ixr, '> nx + nbx =', nx + nbx
            stop
          end if
        end if

        ! intermediate levels x(i+1/2) = x(i) + 0.5*dx
        xr = x(ixr + ix0) + 0.5 * dx
        xl = x(ixl + ix0) + 0.5 * dx

        ! values of var. at the eight corners of the interpolation region
        ! var = u(i+1/2,j,k)
        flwlbd = var%u(ixl, jyb, kzd)
        flwlbu = var%u(ixl, jyb, kzu)

        flwlfd = var%u(ixl, jyf, kzd)
        flwlfu = var%u(ixl, jyf, kzu)

        flwrbd = var%u(ixr, jyb, kzd)
        flwrbu = var%u(ixr, jyb, kzu)

        flwrfd = var%u(ixr, jyf, kzd)
        flwrfu = var%u(ixr, jyf, kzu)
      end if
    elseif(flwtpe == 2) then
      ! interpolate v using staggered-grid distribution

      ! for variable at full levels z(k) = lz(0) +  0.5*dz + (k-1)*dz:
      ! (levels below the model bottom are replaced by the first level
      ! in the model domain)

      if(topography) then
        ! FJApr2023
        ! Locate the closest points in zonal direction.
        if(sizeX == 1) then
          ixl = 1
          ixr = 1
        else
          ixl = floor((xlc - 0.5 * dx - lx(0)) / dx) + 1 - ix0
          if(ixl < - nbx) then
            print *, "ERROR IN MEANFLOW: ixl =", ixl, "< - nbx =", - nbx
            stop
          end if
          ixr = ixl + 1
          if(ixr > nx + nbx) then
            print *, "ERROR IN MEANFLOW: ixr =", ixr, "> nx + nbx =", nx + nbx
            stop
          end if
        end if
        xr = x(ixr + ix0)
        xl = x(ixl + ix0)

        ! Locate the closest points in meridional direction.
        if(sizeY == 1) then
          jyb = 1
          jyf = 1
        else
          jyb = floor((ylc - ly(0)) / dy) - jy0
          if(jyb < - nby) then
            print *, "ERROR IN MEANFLOW: jyb =", jyb, "< - nby =", - nby
            stop
          end if
          jyf = jyb + 1
          if(jyf > ny + nby) then
            print *, "ERROR IN MEANFLOW: jyf =", jyf, "> ny + nby =", ny + nby
            stop
          end if
        end if
        yf = y(jyf + jy0) + 0.5 * dy
        yb = y(jyb + jy0) + 0.5 * dy

        ! Locate the closest points in vertical direction.

        kzlbd = max(1, floor((levelTFC(ixl, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzlbu = max(1, floor((levelTFC(ixl, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzlbd > nz) then
          kzlbu = nz
          kzlbd = nz
        end if
        zlbd = zTFC(ixl, jyb, kzlbd)
        zlbu = zTFC(ixl, jyb, kzlbu)

        kzlfd = max(1, floor((levelTFC(ixl, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzlfu = max(1, floor((levelTFC(ixl, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzlfd > nz) then
          kzlfu = nz
          kzlfd = nz
        end if
        zlfd = zTFC(ixl, jyf, kzlfd)
        zlfu = zTFC(ixl, jyf, kzlfu)

        kzrbd = max(1, floor((levelTFC(ixr, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzrbu = max(1, floor((levelTFC(ixr, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzrbd > nz) then
          kzrbu = nz
          kzrbd = nz
        end if
        zrbd = zTFC(ixr, jyb, kzrbd)
        zrbu = zTFC(ixr, jyb, kzrbu)

        kzrfd = max(1, floor((levelTFC(ixr, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzrfu = max(1, floor((levelTFC(ixr, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzrfd > nz) then
          kzrfu = nz
          kzrfd = nz
        end if
        zrfd = zTFC(ixr, jyf, kzrfd)
        zrfu = zTFC(ixr, jyf, kzrfu)

        ! Assign the values.

        flwlbd = var%v(ixl, jyb, kzlbd)
        flwlbu = var%v(ixl, jyb, kzlbu)

        flwlfd = var%v(ixl, jyf, kzlfd)
        flwlfu = var%v(ixl, jyf, kzlfu)

        flwrbd = var%v(ixr, jyb, kzrbd)
        flwrbu = var%v(ixr, jyb, kzrbu)

        flwrfd = var%v(ixr, jyf, kzrfd)
        flwrfu = var%v(ixr, jyf, kzrfu)
      else
        ! index for lower level used for linear interpolation
        kzd = max(1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 1)

        ! index for upper level used for linear interpolation
        kzu = max(1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 2)

        if(kzd > nz) then
          kzu = nz
          kzd = nz
        end if

        ! full levels z(k)
        zu = z(kzu)
        zd = z(kzd)

        if(sizeY == 1) then
          jyb = 1
          jyf = 1
        else
          ! for variable at intermediate levels
          ! y(j+1/2) = y(j) + 0.5*dy = ly(0) +  j*dy:

          ! index for backward level used for linear interpolation
          jyb = floor((ylc - ly(0)) / dy) - jy0

          if(jyb < - nby) then
            print *, 'ERROR IN MEANFLOW: jyb =', jyb, '< -nby =', - nby
            stop
          end if

          ! index for forward level used for linear interpolation
          jyf = jyb + 1

          if(jyf > ny + nby) then
            print *, 'ERROR IN MEANFLOW: jyf =', jyf, '> ny + nby =', ny + nby
            stop
          end if
        end if

        ! intermediate levels y(j+1/2) = y(j) + 0.5*dy
        yf = y(jyf + jy0) + 0.5 * dy
        yb = y(jyb + jy0) + 0.5 * dy

        if(sizeX == 1) then
          ixl = 1
          ixr = 1
        else
          ! for variable at full levels x(i) = lx(0) + 0.5*dx + (i-1)*dx

          ! index for leftmost level used for linear interpolation
          ixl = floor((xlc - 0.5 * dx - lx(0)) / dx) + 1 - ix0

          if(ixl < - nbx) then
            print *, 'ERROR IN MEANFLOW: ixl =', ixl, '< -nbx =', - nbx
            ! testb
            print *, 'variable = v'
            print *, 'xlc =', xlc
            print *, 'lx(0) =', lx(0)
            print *, 'lx(1) =', lx(1)
            print *, 'dx =', dx
            print *, 'lRef =', lRef
            ! teste
            stop
          end if

          ! index for rightmost level used for linear interpolation
          ixr = ixl + 1

          if(ixr > nx + nbx) then
            print *, 'ERROR IN MEANFLOW: ixr =', ixr, '> nx + nbx =', nx + nbx
            stop
          end if
        end if

        ! full levels x(i)
        xr = x(ixr + ix0)
        xl = x(ixl + ix0)

        ! values of var. at the eight corners of the interpolation region
        ! var = v(i,j+1/2,k)
        flwlbd = var%v(ixl, jyb, kzd)
        flwlbu = var%v(ixl, jyb, kzu)

        flwlfd = var%v(ixl, jyf, kzd)
        flwlfu = var%v(ixl, jyf, kzu)

        flwrbd = var%v(ixr, jyb, kzd)
        flwrbu = var%v(ixr, jyb, kzu)

        flwrfd = var%v(ixr, jyf, kzd)
        flwrfu = var%v(ixr, jyf, kzu)
      end if
    elseif(flwtpe == 3) then
      ! interpolate w using staggered-grid distribution

      ! for variable at intermediate levels
      ! z(k+1/2) = lz(0) + k*dz:

      if(topography) then
        ! FJApr2023
        ! Locate the closest points in zonal direction.
        if(sizeX == 1) then
          ixl = 1
          ixr = 1
        else
          ixl = floor((xlc - 0.5 * dx - lx(0)) / dx) + 1 - ix0
          if(ixl < - nbx) then
            print *, "ERROR IN MEANFLOW: ixl =", ixl, "< - nbx =", - nbx
            stop
          end if
          ixr = ixl + 1
          if(ixr > nx + nbx) then
            print *, "ERROR IN MEANFLOW: ixr =", ixr, "> nx + nbx =", nx + nbx
            stop
          end if
        end if
        xr = x(ixr + ix0)
        xl = x(ixl + ix0)

        ! Locate the closest points in meridional direction.
        if(sizeY == 1) then
          jyb = 1
          jyf = 1
        else
          jyb = floor((ylc - 0.5 * dy - ly(0)) / dy) + 1 - jy0
          if(jyb < - nby) then
            print *, "ERROR IN MEANFLOW: jyb =", jyb, "< - nby =", - nby
            stop
          end if
          jyf = jyb + 1
          if(jyf > ny + nby) then
            print *, "ERROR IN MEANFLOW: jyf =", jyf, "> ny + nby =", ny + nby
            stop
          end if
        end if
        yf = y(jyf + jy0)
        yb = y(jyb + jy0)

        ! Locate the closest points in vertical direction.

        kzlbd = max(- nbz, floor((levelTFC(ixl, jyb, zlc) - lz(0)) / dz))
        kzlbu = kzlbd + 1
        if(kzlbd > nz) then
          kzlbu = nz
          kzlbd = nz
        end if
        zlbd = zTildeTFC(ixl, jyb, kzlbd)
        zlbu = zTildeTFC(ixl, jyb, kzlbu)

        kzlfd = max(- nbz, floor((levelTFC(ixl, jyf, zlc) - lz(0)) / dz))
        kzlfu = kzlfd + 1
        if(kzlfd > nz) then
          kzlfu = nz
          kzlfd = nz
        end if
        zlfd = zTildeTFC(ixl, jyf, kzlfd)
        zlfu = zTildeTFC(ixl, jyf, kzlfu)

        kzrbd = max(- nbz, floor((levelTFC(ixr, jyb, zlc) - lz(0)) / dz))
        kzrbu = kzrbd + 1
        if(kzrbd > nz) then
          kzrbu = nz
          kzrbd = nz
        end if
        zrbd = zTildeTFC(ixr, jyb, kzrbd)
        zrbu = zTildeTFC(ixr, jyb, kzrbu)

        kzrfd = max(- nbz, floor((levelTFC(ixr, jyf, zlc) - lz(0)) / dz))
        kzrfu = kzrfd + 1
        if(kzrfd > nz) then
          kzrfu = nz
          kzrfd = nz
        end if
        zrfd = zTildeTFC(ixr, jyf, kzrfd)
        zrfu = zTildeTFC(ixr, jyf, kzrfu)

        ! Assign the values.

        if(zlbu < topography_surface(ixl, jyb)) then
          flwlbd = 0.0
          flwlbu = 0.0
        else if(zlbd < topography_surface(ixl, jyb)) then
          flwlbd = 0.0
          flwlbu = vertWindTFC(ixl, jyb, kzlbu, var)
        else
          flwlbd = vertWindTFC(ixl, jyb, kzlbd, var)
          flwlbu = vertWindTFC(ixl, jyb, kzlbu, var)
        end if

        if(zlfu < topography_surface(ixl, jyf)) then
          flwlfd = 0.0
          flwlfu = 0.0
        else if(zlfd < topography_surface(ixl, jyf)) then
          flwlfd = 0.0
          flwlfu = vertWindTFC(ixl, jyf, kzlfu, var)
        else
          flwlfd = vertWindTFC(ixl, jyf, kzlfd, var)
          flwlfu = vertWindTFC(ixl, jyf, kzlfu, var)
        end if

        if(zrbu < topography_surface(ixr, jyb)) then
          flwrbd = 0.0
          flwrbu = 0.0
        else if(zrbd < topography_surface(ixr, jyb)) then
          flwrbd = 0.0
          flwrbu = vertWindTFC(ixr, jyb, kzrbu, var)
        else
          flwrbd = vertWindTFC(ixr, jyb, kzrbd, var)
          flwrbu = vertWindTFC(ixr, jyb, kzrbu, var)
        end if

        if(zrfu < topography_surface(ixr, jyf)) then
          flwrfd = 0.0
          flwrfu = 0.0
        else if(zrfd < topography_surface(ixr, jyf)) then
          flwrfd = 0.0
          flwrfu = vertWindTFC(ixr, jyf, kzrfu, var)
        else
          flwrfd = vertWindTFC(ixr, jyf, kzrfd, var)
          flwrfu = vertWindTFC(ixr, jyf, kzrfu, var)
        end if
      else
        ! index for lower level used for linear interpolation
        kzd = max(- nbz, floor((zlc - lz(0)) / dz))

        ! index for upper level used for linear interpolation
        kzu = kzd + 1

        if(kzd > nz) then
          kzu = nz
          kzd = nz
        end if

        ! intermediate levels z(k+1/2) = z(k) + 0.5*dz
        zu = z(kzu) + 0.5 * dz
        zd = z(kzd) + 0.5 * dz

        if(sizeY == 1) then
          jyb = 1
          jyf = 1
        else
          ! for variable at full levels y(j) = ly(0) +  0.5*dy + (j-1)*dy:

          ! index for backward level used for linear interpolation
          jyb = floor((ylc - 0.5 * dy - ly(0)) / dy) + 1 - jy0

          if(jyb < - nby) then
            print *, 'ERROR IN MEANFLOW: jyb =', jyb, '< -nby =', - nby
            stop
          end if

          ! index for forward level used for linear interpolation
          jyf = jyb + 1

          if(jyf > ny + nby) then
            print *, 'ERROR IN MEANFLOW: jyf =', jyf, '> ny + nby =', ny + nby
            stop
          end if
        end if

        ! full levels y(j)
        yf = y(jyf + jy0)
        yb = y(jyb + jy0)

        if(sizeX == 1) then
          ixl = 1
          ixr = 1
        else
          ! for variable at full levels x(i) = lx(0) + 0.5*dx + (i-1)*dx

          ! index for leftmost level used for linear interpolation
          ixl = floor((xlc - 0.5 * dx - lx(0)) / dx) + 1 - ix0

          if(ixl < - nbx) then
            print *, 'ERROR IN MEANFLOW: ixl =', ixl, '< -nbx =', - nbx
            ! testb
            print *, 'variable = w'
            print *, 'xlc =', xlc
            print *, 'lx(0) =', lx(0)
            print *, 'lx(1) =', lx(1)
            print *, 'dx =', dx
            print *, 'lRef =', lRef
            ! teste
            stop
          end if

          ! index for rightmost level used for linear interpolation
          ixr = ixl + 1

          if(ixr > nx + nbx) then
            print *, 'ERROR IN MEANFLOW: ixr =', ixr, '> nx + nbx =', nx + nbx
            stop
          end if
        end if

        ! full levels x(i)
        xr = x(ixr + ix0)
        xl = x(ixl + ix0)

        ! values of var. at the eight corners of the interpolation region
        ! var = w(i,j,k+1/2)
        ! at levels below the model bottom use w = 0
        if(zu < lz(0)) then
          flwlbd = 0.0
          flwlbu = 0.0

          flwlfd = 0.0
          flwlfu = 0.0

          flwrbd = 0.0
          flwrbu = 0.0

          flwrfd = 0.0
          flwrfu = 0.0
        elseif(zd < lz(0)) then
          flwlbd = 0.0
          flwlbu = var%w(ixl, jyb, kzu)

          flwlfd = 0.0
          flwlfu = var%w(ixl, jyf, kzu)

          flwrbd = 0.0
          flwrbu = var%w(ixr, jyb, kzu)

          flwrfd = 0.0
          flwrfu = var%w(ixr, jyf, kzu)
        else
          flwlbd = var%w(ixl, jyb, kzd)
          flwlbu = var%w(ixl, jyb, kzu)

          flwlfd = var%w(ixl, jyf, kzd)
          flwlfu = var%w(ixl, jyf, kzu)

          flwrbd = var%w(ixr, jyb, kzd)
          flwrbu = var%w(ixr, jyb, kzu)

          flwrfd = var%w(ixr, jyf, kzd)
          flwrfu = var%w(ixr, jyf, kzu)
        end if
      end if
    elseif(flwtpe == 11) then
      ! interpolate du/dx using staggered-grid distribution
      ! du/dx (i,j,k)) = (u(i+1/2,j,k) - u(i-1/2,j,k))/dx, hence

      if(topography) then
        ! FJApr2023
        ! Locate the closest points in zonal direction.
        if(sizeX == 1) then
          flw = 0.0
          return
        else
          ixl = floor((xlc - 0.5 * dx - lx(0)) / dx) + 1 - ix0
          if(ixl - 1 < - nbx) then
            print *, "ERROR IN MEANFLOW: ixl - 1 =", ixl - 1, "< - nbx =", - nbx
            stop
          end if
          ixr = ixl + 1
          if(ixr > nx + nbx) then
            print *, "ERROR IN MEANFLOW: ixr =", ixr, "> nx + nbx =", nx + nbx
            stop
          end if
        end if
        xr = x(ixr + ix0)
        xl = x(ixl + ix0)

        ! Locate the closest points in meridional direction.
        if(sizeY == 1) then
          jyb = 1
          jyf = 1
        else
          jyb = floor((ylc - 0.5 * dy - ly(0)) / dy) + 1 - jy0
          if(jyb < - nby) then
            print *, "ERROR IN MEANFLOW: jyb =", jyb, "< - nby =", - nby
            stop
          end if
          jyf = jyb + 1
          if(jyf > ny + nby) then
            print *, "ERROR IN MEANFLOW: jyf =", jyf, "> ny + nby =", ny + nby
            stop
          end if
        end if
        yf = y(jyf + jy0)
        yb = y(jyb + jy0)

        ! Locate the closest points in vertical direction.

        kzlbd = max(1, floor((levelTFC(ixl, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzlbu = max(1, floor((levelTFC(ixl, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzlbd > nz) then
          kzlbu = nz
          kzlbd = nz
        end if
        zlbd = zTFC(ixl, jyb, kzlbd)
        zlbu = zTFC(ixl, jyb, kzlbu)

        kzlfd = max(1, floor((levelTFC(ixl, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzlfu = max(1, floor((levelTFC(ixl, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzlfd > nz) then
          kzlfu = nz
          kzlfd = nz
        end if
        zlfd = zTFC(ixl, jyf, kzlfd)
        zlfu = zTFC(ixl, jyf, kzlfu)

        kzrbd = max(1, floor((levelTFC(ixr, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzrbu = max(1, floor((levelTFC(ixr, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzrbd > nz) then
          kzrbu = nz
          kzrbd = nz
        end if
        zrbd = zTFC(ixr, jyb, kzrbd)
        zrbu = zTFC(ixr, jyb, kzrbu)

        kzrfd = max(1, floor((levelTFC(ixr, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzrfu = max(1, floor((levelTFC(ixr, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzrfd > nz) then
          kzrfu = nz
          kzrfd = nz
        end if
        zrfd = zTFC(ixr, jyf, kzrfd)
        zrfu = zTFC(ixr, jyf, kzrfu)

        ! Assign the values.

        flwlbd = (0.5 * ((jac(ixl, jyb, kzlbd) + jac(ixl + 1, jyb, kzlbd)) &
            &* var%u(ixl, jyb, kzlbd) - (jac(ixl, jyb, kzlbd) + jac(ixl - 1, &
            &jyb, kzlbd)) * var%u(ixl - 1, jyb, kzlbd)) / dx + 0.25 &
            &* (jac(ixl, jyb, kzlbd + 1) * met(ixl, jyb, kzlbd + 1, 1, 3) &
            &* (var%u(ixl, jyb, kzlbd + 1) + var%u(ixl - 1, jyb, kzlbd + 1)) &
            &- jac(ixl, jyb, kzlbd - 1) * met(ixl, jyb, kzlbd - 1, 1, 3) &
            &* (var%u(ixl, jyb, kzlbd - 1) + var%u(ixl - 1, jyb, kzlbd - 1))) &
            &/ dz) / jac(ixl, jyb, kzlbd)
        flwlbu = (0.5 * ((jac(ixl, jyb, kzlbu) + jac(ixl + 1, jyb, kzlbu)) &
            &* var%u(ixl, jyb, kzlbu) - (jac(ixl, jyb, kzlbu) + jac(ixl - 1, &
            &jyb, kzlbu)) * var%u(ixl - 1, jyb, kzlbu)) / dx + 0.25 &
            &* (jac(ixl, jyb, kzlbu + 1) * met(ixl, jyb, kzlbu + 1, 1, 3) &
            &* (var%u(ixl, jyb, kzlbu + 1) + var%u(ixl - 1, jyb, kzlbu + 1)) &
            &- jac(ixl, jyb, kzlbu - 1) * met(ixl, jyb, kzlbu - 1, 1, 3) &
            &* (var%u(ixl, jyb, kzlbu - 1) + var%u(ixl - 1, jyb, kzlbu - 1))) &
            &/ dz) / jac(ixl, jyb, kzlbu)

        flwlfd = (0.5 * ((jac(ixl, jyf, kzlfd) + jac(ixl + 1, jyf, kzlfd)) &
            &* var%u(ixl, jyf, kzlfd) - (jac(ixl, jyf, kzlfd) + jac(ixl - 1, &
            &jyf, kzlfd)) * var%u(ixl - 1, jyf, kzlfd)) / dx + 0.25 &
            &* (jac(ixl, jyf, kzlfd + 1) * met(ixl, jyf, kzlfd + 1, 1, 3) &
            &* (var%u(ixl, jyf, kzlfd + 1) + var%u(ixl - 1, jyf, kzlfd + 1)) &
            &- jac(ixl, jyf, kzlfd - 1) * met(ixl, jyf, kzlfd - 1, 1, 3) &
            &* (var%u(ixl, jyf, kzlfd - 1) + var%u(ixl - 1, jyf, kzlfd - 1))) &
            &/ dz) / jac(ixl, jyf, kzlfd)
        flwlfu = (0.5 * ((jac(ixl, jyf, kzlfu) + jac(ixl + 1, jyf, kzlfu)) &
            &* var%u(ixl, jyf, kzlfu) - (jac(ixl, jyf, kzlfu) + jac(ixl - 1, &
            &jyf, kzlfu)) * var%u(ixl - 1, jyf, kzlfu)) / dx + 0.25 &
            &* (jac(ixl, jyf, kzlfu + 1) * met(ixl, jyf, kzlfu + 1, 1, 3) &
            &* (var%u(ixl, jyf, kzlfu + 1) + var%u(ixl - 1, jyf, kzlfu + 1)) &
            &- jac(ixl, jyf, kzlfu - 1) * met(ixl, jyf, kzlfu - 1, 1, 3) &
            &* (var%u(ixl, jyf, kzlfu - 1) + var%u(ixl - 1, jyf, kzlfu - 1))) &
            &/ dz) / jac(ixl, jyf, kzlfu)

        flwrbd = (0.5 * ((jac(ixr, jyb, kzrbd) + jac(ixr + 1, jyb, kzrbd)) &
            &* var%u(ixr, jyb, kzrbd) - (jac(ixr, jyb, kzrbd) + jac(ixr - 1, &
            &jyb, kzrbd)) * var%u(ixr - 1, jyb, kzrbd)) / dx + 0.25 &
            &* (jac(ixr, jyb, kzrbd + 1) * met(ixr, jyb, kzrbd + 1, 1, 3) &
            &* (var%u(ixr, jyb, kzrbd + 1) + var%u(ixr - 1, jyb, kzrbd + 1)) &
            &- jac(ixr, jyb, kzrbd - 1) * met(ixr, jyb, kzrbd - 1, 1, 3) &
            &* (var%u(ixr, jyb, kzrbd - 1) + var%u(ixr - 1, jyb, kzrbd - 1))) &
            &/ dz) / jac(ixr, jyb, kzrbd)
        flwrbu = (0.5 * ((jac(ixr, jyb, kzrbu) + jac(ixr + 1, jyb, kzrbu)) &
            &* var%u(ixr, jyb, kzrbu) - (jac(ixr, jyb, kzrbu) + jac(ixr - 1, &
            &jyb, kzrbu)) * var%u(ixr - 1, jyb, kzrbu)) / dx + 0.25 &
            &* (jac(ixr, jyb, kzrbu + 1) * met(ixr, jyb, kzrbu + 1, 1, 3) &
            &* (var%u(ixr, jyb, kzrbu + 1) + var%u(ixr - 1, jyb, kzrbu + 1)) &
            &- jac(ixr, jyb, kzrbu - 1) * met(ixr, jyb, kzrbu - 1, 1, 3) &
            &* (var%u(ixr, jyb, kzrbu - 1) + var%u(ixr - 1, jyb, kzrbu - 1))) &
            &/ dz) / jac(ixr, jyb, kzrbu)

        flwrfd = (0.5 * ((jac(ixr, jyf, kzrfd) + jac(ixr + 1, jyf, kzrfd)) &
            &* var%u(ixr, jyf, kzrfd) - (jac(ixr, jyf, kzrfd) + jac(ixr - 1, &
            &jyf, kzrfd)) * var%u(ixr - 1, jyf, kzrfd)) / dx + 0.25 &
            &* (jac(ixr, jyf, kzrfd + 1) * met(ixr, jyf, kzrfd + 1, 1, 3) &
            &* (var%u(ixr, jyf, kzrfd + 1) + var%u(ixr - 1, jyf, kzrfd + 1)) &
            &- jac(ixr, jyf, kzrfd - 1) * met(ixr, jyf, kzrfd - 1, 1, 3) &
            &* (var%u(ixr, jyf, kzrfd - 1) + var%u(ixr - 1, jyf, kzrfd - 1))) &
            &/ dz) / jac(ixr, jyf, kzrfd)
        flwrfu = (0.5 * ((jac(ixr, jyf, kzrfu) + jac(ixr + 1, jyf, kzrfu)) &
            &* var%u(ixr, jyf, kzrfu) - (jac(ixr, jyf, kzrfu) + jac(ixr - 1, &
            &jyf, kzrfu)) * var%u(ixr - 1, jyf, kzrfu)) / dx + 0.25 &
            &* (jac(ixr, jyf, kzrfu + 1) * met(ixr, jyf, kzrfu + 1, 1, 3) &
            &* (var%u(ixr, jyf, kzrfu + 1) + var%u(ixr - 1, jyf, kzrfu + 1)) &
            &- jac(ixr, jyf, kzrfu - 1) * met(ixr, jyf, kzrfu - 1, 1, 3) &
            &* (var%u(ixr, jyf, kzrfu - 1) + var%u(ixr - 1, jyf, kzrfu - 1))) &
            &/ dz) / jac(ixr, jyf, kzrfu)
      else
        if(sizeX == 1) then
          ! no derivative if there is no x dependence
          flw = 0.0

          return
        else
          ! for variable at full levels x(i) = lx(0) + 0.5*dx + (i-1)*dx

          ! index for leftmost level used for linear interpolation
          ixl = floor((xlc - 0.5 * dx - lx(0)) / dx) + 1 - ix0

          if(ixl - 1 < - nbx) then
            print *, 'ERROR IN MEANFLOW: ixl - 1 =', ixl - 1, '< -nbx =', - nbx
            ! testb
            print *, 'variable = du/dx'
            print *, 'xlc =', xlc
            print *, 'lx(0) =', lx(0)
            print *, 'lx(1) =', lx(1)
            print *, 'dx =', dx
            print *, 'lRef =', lRef
            ! teste
            stop
          end if

          ! index for rightmost level used for linear interpolation
          ixr = ixl + 1

          if(ixr > nx + nbx) then
            print *, 'ERROR IN MEANFLOW: ixr =', ixr, '> nx + nbx =', nx + nbx
            stop
          end if
        end if

        ! full levels x(i)
        xr = x(ixr + ix0)
        xl = x(ixl + ix0)

        ! for variable at full levels z(k) = lz(0) +  0.5*dz + (k-1)*dz:
        ! (levels below the model bottom are replaced by the first level
        ! in the model domain)

        ! index for lower level used for linear interpolation
        kzd = max(1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 1)

        ! index for upper level used for linear interpolation
        kzu = max(1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 2)

        if(kzd > nz) then
          kzu = nz
          kzd = nz
        end if

        ! full levels z(k)
        zu = z(kzu)
        zd = z(kzd)

        if(sizeY == 1) then
          jyb = 1
          jyf = 1
        else
          ! for variable at full levels y(j) = ly(0) +  0.5*dy + (j-1)*dy:

          ! index for backward level used for linear interpolation
          jyb = floor((ylc - 0.5 * dy - ly(0)) / dy) + 1 - jy0

          if(jyb < - nby) then
            print *, 'ERROR IN MEANFLOW: jyb =', jyb, '< -nby =', - nby
            stop
          end if

          ! index for forward level used for linear interpolation
          jyf = jyb + 1

          if(jyf > ny + nby) then
            print *, 'MEANFLOW: jyf =', jyf, '> ny + nby =', ny + nby
            stop
          end if
        end if

        ! full levels y(j)
        yf = y(jyf + jy0)
        yb = y(jyb + jy0)

        ! values of var. at the eight corners of the interpolation region
        ! du/dx (i,j,k)) = (u(i+1/2,j,k) - u(i-1/2,j,k))/dx, hence
        flwlbd = (var%u(ixl, jyb, kzd) - var%u(ixl - 1, jyb, kzd)) / dx
        flwlbu = (var%u(ixl, jyb, kzu) - var%u(ixl - 1, jyb, kzu)) / dx

        flwlfd = (var%u(ixl, jyf, kzd) - var%u(ixl - 1, jyf, kzd)) / dx
        flwlfu = (var%u(ixl, jyf, kzu) - var%u(ixl - 1, jyf, kzu)) / dx

        flwrbd = (var%u(ixr, jyb, kzd) - var%u(ixr - 1, jyb, kzd)) / dx
        flwrbu = (var%u(ixr, jyb, kzu) - var%u(ixr - 1, jyb, kzu)) / dx

        flwrfd = (var%u(ixr, jyf, kzd) - var%u(ixr - 1, jyf, kzd)) / dx
        flwrfu = (var%u(ixr, jyf, kzu) - var%u(ixr - 1, jyf, kzu)) / dx
      end if
    elseif(flwtpe == 12) then
      ! interpolate du/dy using staggered-grid distribution
      ! du/dy (i+1/2,j+1/2,k)) = (u(i+1/2,j+1,k) - u(i+1/2,j,k))/dy, hence

      if(topography) then
        ! FJApr2023
        ! Locate the closest points in zonal direction.
        if(sizeX == 1) then
          ixl = 1
          ixr = 1
        else
          ixl = floor((xlc - lx(0)) / dx) - ix0
          if(ixl < - nbx) then
            print *, "ERROR IN MEANFLOW: ixl =", ixl, "< - nbx =", - nbx
            stop
          end if
          ixr = ixl + 1
          if(ixr > nx + nbx) then
            print *, "ERROR IN MEANFLOW: ixr =", ixr, "> nx + nbx =", nx + nbx
            stop
          end if
        end if
        xr = x(ixr + ix0) + 0.5 * dx
        xl = x(ixl + ix0) + 0.5 * dx

        ! Locate the closest points in meridional direction.
        if(sizeY == 1) then
          flw = 0.0
          return
        else
          jyb = floor((ylc - ly(0)) / dy) - jy0
          if(jyb < - nby) then
            print *, "ERROR IN MEANFLOW: jyb =", jyb, "< - nby =", - nby
            stop
          end if
          jyf = jyb + 1
          if(jyf + 1 > ny + nby) then
            print *, "ERROR IN MEANFLOW: jyf + 1 =", jyf + 1, "> ny + nby =", &
                &ny + nby
            stop
          end if
        end if
        yf = y(jyf + jy0) + 0.5 * dy
        yb = y(jyb + jy0) + 0.5 * dy

        ! Locate the closest points in vertical direction.

        kzlbd = max(1, floor((levelTFC(ixl, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzlbu = max(1, floor((levelTFC(ixl, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzlbd > nz) then
          kzlbu = nz
          kzlbd = nz
        end if
        zlbd = zTFC(ixl, jyb, kzlbd)
        zlbu = zTFC(ixl, jyb, kzlbu)

        kzlfd = max(1, floor((levelTFC(ixl, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzlfu = max(1, floor((levelTFC(ixl, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzlfd > nz) then
          kzlfu = nz
          kzlfd = nz
        end if
        zlfd = zTFC(ixl, jyf, kzlfd)
        zlfu = zTFC(ixl, jyf, kzlfu)

        kzrbd = max(1, floor((levelTFC(ixr, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzrbu = max(1, floor((levelTFC(ixr, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzrbd > nz) then
          kzrbu = nz
          kzrbd = nz
        end if
        zrbd = zTFC(ixr, jyb, kzrbd)
        zrbu = zTFC(ixr, jyb, kzrbu)

        kzrfd = max(1, floor((levelTFC(ixr, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzrfu = max(1, floor((levelTFC(ixr, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzrfd > nz) then
          kzrfu = nz
          kzrfd = nz
        end if
        zrfd = zTFC(ixr, jyf, kzrfd)
        zrfu = zTFC(ixr, jyf, kzrfu)

        ! Assign the values.

        flwlbd = (0.5 * ((jac(ixl, jyb + 1, kzlbd) + jac(ixl + 1, jyb + 1, &
            &kzlbd)) * var%u(ixl, jyb + 1, kzlbd) - (jac(ixl, jyb, kzlbd) &
            &+ jac(ixl + 1, jyb, kzlbd)) * var%u(ixl, jyb, kzlbd)) / dy &
            &+ 0.0625 * ((jac(ixl, jyb, kzlbd + 1) * met(ixl, jyb, kzlbd + 1, &
            &2, 3) + jac(ixl + 1, jyb, kzlbd + 1) * met(ixl + 1, jyb, kzlbd &
            &+ 1, 2, 3) + jac(ixl, jyb + 1, kzlbd + 1) * met(ixl, jyb + 1, &
            &kzlbd + 1, 2, 3) + jac(ixl + 1, jyb + 1, kzlbd + 1) * met(ixl &
            &+ 1, jyb + 1, kzlbd + 1, 2, 3)) * (var%u(ixl, jyb, kzlbd + 1) &
            &+ var%u(ixl, jyb + 1, kzlbd + 1)) - (jac(ixl, jyb, kzlbd - 1) &
            &* met(ixl, jyb, kzlbd - 1, 2, 3) + jac(ixl + 1, jyb, kzlbd - 1) &
            &* met(ixl + 1, jyb, kzlbd - 1, 2, 3) + jac(ixl, jyb + 1, kzlbd &
            &- 1) * met(ixl, jyb + 1, kzlbd - 1, 2, 3) + jac(ixl + 1, jyb + 1, &
            &kzlbd - 1) * met(ixl + 1, jyb + 1, kzlbd - 1, 2, 3)) &
            &* (var%u(ixl, jyb, kzlbd - 1) + var%u(ixl, jyb + 1, kzlbd - 1))) &
            &/ dz) * 4.0 / (jac(ixl, jyb, kzlbd) + jac(ixl + 1, jyb, kzlbd) &
            &+ jac(ixl, jyb + 1, kzlbd) + jac(ixl + 1, jyb + 1, kzlbd))
        flwlbu = (0.5 * ((jac(ixl, jyb + 1, kzlbu) + jac(ixl + 1, jyb + 1, &
            &kzlbu)) * var%u(ixl, jyb + 1, kzlbu) - (jac(ixl, jyb, kzlbu) &
            &+ jac(ixl + 1, jyb, kzlbu)) * var%u(ixl, jyb, kzlbu)) / dy &
            &+ 0.0625 * ((jac(ixl, jyb, kzlbu + 1) * met(ixl, jyb, kzlbu + 1, &
            &2, 3) + jac(ixl + 1, jyb, kzlbu + 1) * met(ixl + 1, jyb, kzlbu &
            &+ 1, 2, 3) + jac(ixl, jyb + 1, kzlbu + 1) * met(ixl, jyb + 1, &
            &kzlbu + 1, 2, 3) + jac(ixl + 1, jyb + 1, kzlbu + 1) * met(ixl &
            &+ 1, jyb + 1, kzlbu + 1, 2, 3)) * (var%u(ixl, jyb, kzlbu + 1) &
            &+ var%u(ixl, jyb + 1, kzlbu + 1)) - (jac(ixl, jyb, kzlbu - 1) &
            &* met(ixl, jyb, kzlbu - 1, 2, 3) + jac(ixl + 1, jyb, kzlbu - 1) &
            &* met(ixl + 1, jyb, kzlbu - 1, 2, 3) + jac(ixl, jyb + 1, kzlbu &
            &- 1) * met(ixl, jyb + 1, kzlbu - 1, 2, 3) + jac(ixl + 1, jyb + 1, &
            &kzlbu - 1) * met(ixl + 1, jyb + 1, kzlbu - 1, 2, 3)) &
            &* (var%u(ixl, jyb, kzlbu - 1) + var%u(ixl, jyb + 1, kzlbu - 1))) &
            &/ dz) * 4.0 / (jac(ixl, jyb, kzlbu) + jac(ixl + 1, jyb, kzlbu) &
            &+ jac(ixl, jyb + 1, kzlbu) + jac(ixl + 1, jyb + 1, kzlbu))

        flwlfd = (0.5 * ((jac(ixl, jyf + 1, kzlfd) + jac(ixl + 1, jyf + 1, &
            &kzlfd)) * var%u(ixl, jyf + 1, kzlfd) - (jac(ixl, jyf, kzlfd) &
            &+ jac(ixl + 1, jyf, kzlfd)) * var%u(ixl, jyf, kzlfd)) / dy &
            &+ 0.0625 * ((jac(ixl, jyf, kzlfd + 1) * met(ixl, jyf, kzlfd + 1, &
            &2, 3) + jac(ixl + 1, jyf, kzlfd + 1) * met(ixl + 1, jyf, kzlfd &
            &+ 1, 2, 3) + jac(ixl, jyf + 1, kzlfd + 1) * met(ixl, jyf + 1, &
            &kzlfd + 1, 2, 3) + jac(ixl + 1, jyf + 1, kzlfd + 1) * met(ixl &
            &+ 1, jyf + 1, kzlfd + 1, 2, 3)) * (var%u(ixl, jyf, kzlfd + 1) &
            &+ var%u(ixl, jyf + 1, kzlfd + 1)) - (jac(ixl, jyf, kzlfd - 1) &
            &* met(ixl, jyf, kzlfd - 1, 2, 3) + jac(ixl + 1, jyf, kzlfd - 1) &
            &* met(ixl + 1, jyf, kzlfd - 1, 2, 3) + jac(ixl, jyf + 1, kzlfd &
            &- 1) * met(ixl, jyf + 1, kzlfd - 1, 2, 3) + jac(ixl + 1, jyf + 1, &
            &kzlfd - 1) * met(ixl + 1, jyf + 1, kzlfd - 1, 2, 3)) &
            &* (var%u(ixl, jyf, kzlfd - 1) + var%u(ixl, jyf + 1, kzlfd - 1))) &
            &/ dz) * 4.0 / (jac(ixl, jyf, kzlfd) + jac(ixl + 1, jyf, kzlfd) &
            &+ jac(ixl, jyf + 1, kzlfd) + jac(ixl + 1, jyf + 1, kzlfd))
        flwlfu = (0.5 * ((jac(ixl, jyf + 1, kzlfu) + jac(ixl + 1, jyf + 1, &
            &kzlfu)) * var%u(ixl, jyf + 1, kzlfu) - (jac(ixl, jyf, kzlfu) &
            &+ jac(ixl + 1, jyf, kzlfu)) * var%u(ixl, jyf, kzlfu)) / dy &
            &+ 0.0625 * ((jac(ixl, jyf, kzlfu + 1) * met(ixl, jyf, kzlfu + 1, &
            &2, 3) + jac(ixl + 1, jyf, kzlfu + 1) * met(ixl + 1, jyf, kzlfu &
            &+ 1, 2, 3) + jac(ixl, jyf + 1, kzlfu + 1) * met(ixl, jyf + 1, &
            &kzlfu + 1, 2, 3) + jac(ixl + 1, jyf + 1, kzlfu + 1) * met(ixl &
            &+ 1, jyf + 1, kzlfu + 1, 2, 3)) * (var%u(ixl, jyf, kzlfu + 1) &
            &+ var%u(ixl, jyf + 1, kzlfu + 1)) - (jac(ixl, jyf, kzlfu - 1) &
            &* met(ixl, jyf, kzlfu - 1, 2, 3) + jac(ixl + 1, jyf, kzlfu - 1) &
            &* met(ixl + 1, jyf, kzlfu - 1, 2, 3) + jac(ixl, jyf + 1, kzlfu &
            &- 1) * met(ixl, jyf + 1, kzlfu - 1, 2, 3) + jac(ixl + 1, jyf + 1, &
            &kzlfu - 1) * met(ixl + 1, jyf + 1, kzlfu - 1, 2, 3)) &
            &* (var%u(ixl, jyf, kzlfu - 1) + var%u(ixl, jyf + 1, kzlfu - 1))) &
            &/ dz) * 4.0 / (jac(ixl, jyf, kzlfu) + jac(ixl + 1, jyf, kzlfu) &
            &+ jac(ixl, jyf + 1, kzlfu) + jac(ixl + 1, jyf + 1, kzlfu))

        flwrbd = (0.5 * ((jac(ixr, jyb + 1, kzrbd) + jac(ixr + 1, jyb + 1, &
            &kzrbd)) * var%u(ixr, jyb + 1, kzrbd) - (jac(ixr, jyb, kzrbd) &
            &+ jac(ixr + 1, jyb, kzrbd)) * var%u(ixr, jyb, kzrbd)) / dy &
            &+ 0.0625 * ((jac(ixr, jyb, kzrbd + 1) * met(ixr, jyb, kzrbd + 1, &
            &2, 3) + jac(ixr + 1, jyb, kzrbd + 1) * met(ixr + 1, jyb, kzrbd &
            &+ 1, 2, 3) + jac(ixr, jyb + 1, kzrbd + 1) * met(ixr, jyb + 1, &
            &kzrbd + 1, 2, 3) + jac(ixr + 1, jyb + 1, kzrbd + 1) * met(ixr &
            &+ 1, jyb + 1, kzrbd + 1, 2, 3)) * (var%u(ixr, jyb, kzrbd + 1) &
            &+ var%u(ixr, jyb + 1, kzrbd + 1)) - (jac(ixr, jyb, kzrbd - 1) &
            &* met(ixr, jyb, kzrbd - 1, 2, 3) + jac(ixr + 1, jyb, kzrbd - 1) &
            &* met(ixr + 1, jyb, kzrbd - 1, 2, 3) + jac(ixr, jyb + 1, kzrbd &
            &- 1) * met(ixr, jyb + 1, kzrbd - 1, 2, 3) + jac(ixr + 1, jyb + 1, &
            &kzrbd - 1) * met(ixr + 1, jyb + 1, kzrbd - 1, 2, 3)) &
            &* (var%u(ixr, jyb, kzrbd - 1) + var%u(ixr, jyb + 1, kzrbd - 1))) &
            &/ dz) * 4.0 / (jac(ixr, jyb, kzrbd) + jac(ixr + 1, jyb, kzrbd) &
            &+ jac(ixr, jyb + 1, kzrbd) + jac(ixr + 1, jyb + 1, kzrbd))
        flwrbu = (0.5 * ((jac(ixr, jyb + 1, kzrbu) + jac(ixr + 1, jyb + 1, &
            &kzrbu)) * var%u(ixr, jyb + 1, kzrbu) - (jac(ixr, jyb, kzrbu) &
            &+ jac(ixr + 1, jyb, kzrbu)) * var%u(ixr, jyb, kzrbu)) / dy &
            &+ 0.0625 * ((jac(ixr, jyb, kzrbu + 1) * met(ixr, jyb, kzrbu + 1, &
            &2, 3) + jac(ixr + 1, jyb, kzrbu + 1) * met(ixr + 1, jyb, kzrbu &
            &+ 1, 2, 3) + jac(ixr, jyb + 1, kzrbu + 1) * met(ixr, jyb + 1, &
            &kzrbu + 1, 2, 3) + jac(ixr + 1, jyb + 1, kzrbu + 1) * met(ixr &
            &+ 1, jyb + 1, kzrbu + 1, 2, 3)) * (var%u(ixr, jyb, kzrbu + 1) &
            &+ var%u(ixr, jyb + 1, kzrbu + 1)) - (jac(ixr, jyb, kzrbu - 1) &
            &* met(ixr, jyb, kzrbu - 1, 2, 3) + jac(ixr + 1, jyb, kzrbu - 1) &
            &* met(ixr + 1, jyb, kzrbu - 1, 2, 3) + jac(ixr, jyb + 1, kzrbu &
            &- 1) * met(ixr, jyb + 1, kzrbu - 1, 2, 3) + jac(ixr + 1, jyb + 1, &
            &kzrbu - 1) * met(ixr + 1, jyb + 1, kzrbu - 1, 2, 3)) &
            &* (var%u(ixr, jyb, kzrbu - 1) + var%u(ixr, jyb + 1, kzrbu - 1))) &
            &/ dz) * 4.0 / (jac(ixr, jyb, kzrbu) + jac(ixr + 1, jyb, kzrbu) &
            &+ jac(ixr, jyb + 1, kzrbu) + jac(ixr + 1, jyb + 1, kzrbu))

        flwrfd = (0.5 * ((jac(ixr, jyf + 1, kzrfd) + jac(ixr + 1, jyf + 1, &
            &kzrfd)) * var%u(ixr, jyf + 1, kzrfd) - (jac(ixr, jyf, kzrfd) &
            &+ jac(ixr + 1, jyf, kzrfd)) * var%u(ixr, jyf, kzrfd)) / dy &
            &+ 0.0625 * ((jac(ixr, jyf, kzrfd + 1) * met(ixr, jyf, kzrfd + 1, &
            &2, 3) + jac(ixr + 1, jyf, kzrfd + 1) * met(ixr + 1, jyf, kzrfd &
            &+ 1, 2, 3) + jac(ixr, jyf + 1, kzrfd + 1) * met(ixr, jyf + 1, &
            &kzrfd + 1, 2, 3) + jac(ixr + 1, jyf + 1, kzrfd + 1) * met(ixr &
            &+ 1, jyf + 1, kzrfd + 1, 2, 3)) * (var%u(ixr, jyf, kzrfd + 1) &
            &+ var%u(ixr, jyf + 1, kzrfd + 1)) - (jac(ixr, jyf, kzrfd - 1) &
            &* met(ixr, jyf, kzrfd - 1, 2, 3) + jac(ixr + 1, jyf, kzrfd - 1) &
            &* met(ixr + 1, jyf, kzrfd - 1, 2, 3) + jac(ixr, jyf + 1, kzrfd &
            &- 1) * met(ixr, jyf + 1, kzrfd - 1, 2, 3) + jac(ixr + 1, jyf + 1, &
            &kzrfd - 1) * met(ixr + 1, jyf + 1, kzrfd - 1, 2, 3)) &
            &* (var%u(ixr, jyf, kzrfd - 1) + var%u(ixr, jyf + 1, kzrfd - 1))) &
            &/ dz) * 4.0 / (jac(ixr, jyf, kzrfd) + jac(ixr + 1, jyf, kzrfd) &
            &+ jac(ixr, jyf + 1, kzrfd) + jac(ixr + 1, jyf + 1, kzrfd))
        flwrfu = (0.5 * ((jac(ixr, jyf + 1, kzrfu) + jac(ixr + 1, jyf + 1, &
            &kzrfu)) * var%u(ixr, jyf + 1, kzrfu) - (jac(ixr, jyf, kzrfu) &
            &+ jac(ixr + 1, jyf, kzrfu)) * var%u(ixr, jyf, kzrfu)) / dy &
            &+ 0.0625 * ((jac(ixr, jyf, kzrfu + 1) * met(ixr, jyf, kzrfu + 1, &
            &2, 3) + jac(ixr + 1, jyf, kzrfu + 1) * met(ixr + 1, jyf, kzrfu &
            &+ 1, 2, 3) + jac(ixr, jyf + 1, kzrfu + 1) * met(ixr, jyf + 1, &
            &kzrfu + 1, 2, 3) + jac(ixr + 1, jyf + 1, kzrfu + 1) * met(ixr &
            &+ 1, jyf + 1, kzrfu + 1, 2, 3)) * (var%u(ixr, jyf, kzrfu + 1) &
            &+ var%u(ixr, jyf + 1, kzrfu + 1)) - (jac(ixr, jyf, kzrfu - 1) &
            &* met(ixr, jyf, kzrfu - 1, 2, 3) + jac(ixr + 1, jyf, kzrfu - 1) &
            &* met(ixr + 1, jyf, kzrfu - 1, 2, 3) + jac(ixr, jyf + 1, kzrfu &
            &- 1) * met(ixr, jyf + 1, kzrfu - 1, 2, 3) + jac(ixr + 1, jyf + 1, &
            &kzrfu - 1) * met(ixr + 1, jyf + 1, kzrfu - 1, 2, 3)) &
            &* (var%u(ixr, jyf, kzrfu - 1) + var%u(ixr, jyf + 1, kzrfu - 1))) &
            &/ dz) * 4.0 / (jac(ixr, jyf, kzrfu) + jac(ixr + 1, jyf, kzrfu) &
            &+ jac(ixr, jyf + 1, kzrfu) + jac(ixr + 1, jyf + 1, kzrfu))
      else
        if(sizeY == 1) then
          !no derivative if there is no y dependence

          flw = 0.0

          return
        else
          ! for variable at intermediate levels
          ! y(j+1/2) = y(j) + 0.5*dy = ly(0) +  j*dy:

          ! index for backward level used for linear interpolation
          jyb = floor((ylc - ly(0)) / dy) - jy0

          if(jyb < - nby) then
            print *, 'ERROR IN MEANFLOW: jyb =', jyb, '< -nby =', - nby
            stop
          end if

          ! index for forward level used for linear interpolation
          jyf = jyb + 1

          if(jyf + 1 > ny + nby) then
            print *, 'MEANFLOW: jyf + 1 =', jyf + 1, '> ny + nby =', ny + nby
            stop
          end if
        end if

        ! intermediate levels y(j+1/2) = y(j) + 0.5*dy
        yf = y(jyf + jy0) + 0.5 * dy
        yb = y(jyb + jy0) + 0.5 * dy

        ! for variable at full levels z(k) = lz(0) +  0.5*dz + (k-1)*dz:
        ! (levels below the model bottom are replaced by the first level
        ! in the model domain)

        ! index for lower level used for linear interpolation
        kzd = max(1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 1)

        ! index for upper level used for linear interpolation
        kzu = max(1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 2)

        if(kzd > nz) then
          kzu = nz
          kzd = nz
        end if

        ! full levels z(k)
        zu = z(kzu)
        zd = z(kzd)

        if(sizeX == 1) then
          ixl = 1
          ixr = 1
        else
          ! for variable at intermediate levels
          ! x(i+1/2) = x(i) + 0.5*dx = lx(0) + i*dx:

          ! index for leftmost level used for linear interpolation
          ixl = floor((xlc - lx(0)) / dx) - ix0

          if(ixl < - nbx) then
            print *, 'ERROR IN MEANFLOW: ixl =', ixl, '< -nbx =', - nbx
            ! testb
            print *, 'variable = du/dy'
            print *, 'xlc =', xlc
            print *, 'lx(0) =', lx(0)
            print *, 'lx(1) =', lx(1)
            print *, 'dx =', dx
            print *, 'lRef =', lRef
            ! teste
            stop
          end if

          ! index for rightmost level used for linear interpolation
          ixr = ixl + 1

          if(ixr > nx + nbx) then
            print *, 'ERROR IN MEANFLOW: ixr =', ixr, '> nx + nbx =', nx + nbx
            stop
          end if
        end if

        ! intermediate levels x(i+1/2) = x(i) + 0.5*dx
        xr = x(ixr + ix0) + 0.5 * dx
        xl = x(ixl + ix0) + 0.5 * dx

        ! values of var. at the eight corners of the interpolation region
        ! du/dy (i+1/2,j+1/2,k)) = (u(i+1/2,j+1,k) - u(i+1/2,j,k))/dy, hence
        flwlbd = (var%u(ixl, jyb + 1, kzd) - var%u(ixl, jyb, kzd)) / dy
        flwlbu = (var%u(ixl, jyb + 1, kzu) - var%u(ixl, jyb, kzu)) / dy

        flwlfd = (var%u(ixl, jyf + 1, kzd) - var%u(ixl, jyf, kzd)) / dy
        flwlfu = (var%u(ixl, jyf + 1, kzu) - var%u(ixl, jyf, kzu)) / dy

        flwrbd = (var%u(ixr, jyb + 1, kzd) - var%u(ixr, jyb, kzd)) / dy
        flwrbu = (var%u(ixr, jyb + 1, kzu) - var%u(ixr, jyb, kzu)) / dy

        flwrfd = (var%u(ixr, jyf + 1, kzd) - var%u(ixr, jyf, kzd)) / dy
        flwrfu = (var%u(ixr, jyf + 1, kzu) - var%u(ixr, jyf, kzu)) / dy
      end if
    elseif(flwtpe == 13) then
      ! interpolate du/dz using staggered-grid distribution
      ! du/dz (i+1/2,j,k+1/2)) = (u(i+1/2,j,k+1) - u(i+1/2,j,k))/dz, hence

      if(topography) then
        ! FJApr2023
        ! Locate the closest points in zonal direction.
        if(sizeX == 1) then
          ixl = 1
          ixr = 1
        else
          ixl = floor((xlc - lx(0)) / dx) - ix0
          if(ixl < - nbx) then
            print *, "ERROR IN MEANFLOW: ixl =", ixl, "< - nbx =", - nbx
            stop
          end if
          ixr = ixl + 1
          if(ixr > nx + nbx) then
            print *, "ERROR IN MEANFLOW: ixr =", ixr, "> nx + nbx =", nx + nbx
            stop
          end if
        end if
        xr = x(ixr + ix0) + 0.5 * dx
        xl = x(ixl + ix0) + 0.5 * dx

        ! Locate the closest points in meridional direction.
        if(sizeY == 1) then
          jyb = 1
          jyf = 1
        else
          jyb = floor((ylc - 0.5 * dy - ly(0)) / dy) + 1 - jy0
          if(jyb < - nby) then
            print *, "ERROR IN MEANFLOW: jyb =", jyb, "< - nby =", - nby
            stop
          end if
          jyf = jyb + 1
          if(jyf > ny + nby) then
            print *, "ERROR IN MEANFLOW: jyf =", jyf, "> ny + nby =", ny + nby
            stop
          end if
        end if
        yf = y(jyf + jy0)
        yb = y(jyb + jy0)

        ! Locate the closest points in vertical direction.

        kzlbd = max(- nbz, floor((levelTFC(ixl, jyb, zlc) - lz(0)) / dz))
        kzlbu = kzlbd + 1
        if(kzlbd > nz) then
          kzlbu = nz + 1
          kzlbd = nz + 1
        end if
        zlbd = zTildeTFC(ixl, jyb, kzlbd)
        zlbu = zTildeTFC(ixl, jyb, kzlbu)

        kzlfd = max(- nbz, floor((levelTFC(ixl, jyf, zlc) - lz(0)) / dz))
        kzlfu = kzlfd + 1
        if(kzlfd > nz) then
          kzlfu = nz + 1
          kzlfd = nz + 1
        end if
        zlfd = zTildeTFC(ixl, jyf, kzlfd)
        zlfu = zTildeTFC(ixl, jyf, kzlfu)

        kzrbd = max(- nbz, floor((levelTFC(ixr, jyb, zlc) - lz(0)) / dz))
        kzrbu = kzrbd + 1
        if(kzrbd > nz) then
          kzrbu = nz + 1
          kzrbd = nz + 1
        end if
        zrbd = zTildeTFC(ixr, jyb, kzrbd)
        zrbu = zTildeTFC(ixr, jyb, kzrbu)

        kzrfd = max(- nbz, floor((levelTFC(ixr, jyf, zlc) - lz(0)) / dz))
        kzrfu = kzrfd + 1
        if(kzrfd > nz) then
          kzrfu = nz + 1
          kzrfd = nz + 1
        end if
        zrfd = zTildeTFC(ixr, jyf, kzrfd)
        zrfu = zTildeTFC(ixr, jyf, kzrfu)

        ! Assign the values.

        if(zlbu < topography_surface(ixl, jyb)) then
          flwlbd = 0.0
          flwlbu = 0.0
        else if(zlbd < topography_surface(ixl, jyb)) then
          flwlbd = 0.0
          flwlbu = (var%u(ixl, jyb, kzlbu + 1) - var%u(ixl, jyb, kzlbu)) / dz &
              &* 4.0 / (jac(ixl, jyb, kzlbu) + jac(ixl + 1, jyb, kzlbu) &
              &+ jac(ixl, jyb, kzlbu + 1) + jac(ixl + 1, jyb, kzlbu + 1))
        else
          if(zlbu < lz(1)) then
            flwlbd = (var%u(ixl, jyb, kzlbd + 1) - var%u(ixl, jyb, kzlbd)) &
                &/ dz * 4.0 / (jac(ixl, jyb, kzlbd) + jac(ixl + 1, jyb, kzlbd) &
                &+ jac(ixl, jyb, kzlbd + 1) + jac(ixl + 1, jyb, kzlbd + 1))
            flwlbu = (var%u(ixl, jyb, kzlbu + 1) - var%u(ixl, jyb, kzlbu)) &
                &/ dz * 4.0 / (jac(ixl, jyb, kzlbu) + jac(ixl + 1, jyb, kzlbu) &
                &+ jac(ixl, jyb, kzlbu + 1) + jac(ixl + 1, jyb, kzlbu + 1))
          else if(zlbd < lz(1)) then
            flwlbd = (var%u(ixl, jyb, kzlbd + 1) - var%u(ixl, jyb, kzlbd)) &
                &/ dz * 4.0 / (jac(ixl, jyb, kzlbd) + jac(ixl + 1, jyb, kzlbd) &
                &+ jac(ixl, jyb, kzlbd + 1) + jac(ixl + 1, jyb, kzlbd + 1))
            flwlbu = 0.0
          else
            flwlbd = 0.0
            flwlbu = 0.0
          end if
        end if

        if(zlfu < topography_surface(ixl, jyf)) then
          flwlfd = 0.0
          flwlfu = 0.0
        else if(zlfd < topography_surface(ixl, jyf)) then
          flwlfd = 0.0
          flwlfu = (var%u(ixl, jyf, kzlfu + 1) - var%u(ixl, jyf, kzlfu)) / dz &
              &* 4.0 / (jac(ixl, jyf, kzlfu) + jac(ixl + 1, jyf, kzlfu) &
              &+ jac(ixl, jyf, kzlfu + 1) + jac(ixl + 1, jyf, kzlfu + 1))
        else
          if(zlfu < lz(1)) then
            flwlfd = (var%u(ixl, jyf, kzlfd + 1) - var%u(ixl, jyf, kzlfd)) &
                &/ dz * 4.0 / (jac(ixl, jyf, kzlfd) + jac(ixl + 1, jyf, kzlfd) &
                &+ jac(ixl, jyf, kzlfd + 1) + jac(ixl + 1, jyf, kzlfd + 1))
            flwlfu = (var%u(ixl, jyf, kzlfu + 1) - var%u(ixl, jyf, kzlfu)) &
                &/ dz * 4.0 / (jac(ixl, jyf, kzlfu) + jac(ixl + 1, jyf, kzlfu) &
                &+ jac(ixl, jyf, kzlfu + 1) + jac(ixl + 1, jyf, kzlfu + 1))
          else if(zlfd < lz(1)) then
            flwlfd = (var%u(ixl, jyf, kzlfd + 1) - var%u(ixl, jyf, kzlfd)) &
                &/ dz * 4.0 / (jac(ixl, jyf, kzlfd) + jac(ixl + 1, jyf, kzlfd) &
                &+ jac(ixl, jyf, kzlfd + 1) + jac(ixl + 1, jyf, kzlfd + 1))
            flwlfu = 0.0
          else
            flwlbd = 0.0
            flwlbu = 0.0
          end if
        end if

        if(zrbu < topography_surface(ixr, jyb)) then
          flwrbd = 0.0
          flwrbu = 0.0
        else if(zrbd < topography_surface(ixr, jyb)) then
          flwrbd = 0.0
          flwrbu = (var%u(ixr, jyb, kzrbu + 1) - var%u(ixr, jyb, kzrbu)) / dz &
              &* 4.0 / (jac(ixr, jyb, kzrbu) + jac(ixr + 1, jyb, kzrbu) &
              &+ jac(ixr, jyb, kzrbu + 1) + jac(ixr + 1, jyb, kzrbu + 1))
        else
          if(zrbu < lz(1)) then
            flwrbd = (var%u(ixr, jyb, kzrbd + 1) - var%u(ixr, jyb, kzrbd)) &
                &/ dz * 4.0 / (jac(ixr, jyb, kzrbd) + jac(ixr + 1, jyb, kzrbd) &
                &+ jac(ixr, jyb, kzrbd + 1) + jac(ixr + 1, jyb, kzrbd + 1))
            flwrbu = (var%u(ixr, jyb, kzrbu + 1) - var%u(ixr, jyb, kzrbu)) &
                &/ dz * 4.0 / (jac(ixr, jyb, kzrbu) + jac(ixr + 1, jyb, kzrbu) &
                &+ jac(ixr, jyb, kzrbu + 1) + jac(ixr + 1, jyb, kzrbu + 1))
          else if(zrbd < lz(1)) then
            flwrbd = (var%u(ixr, jyb, kzrbd + 1) - var%u(ixr, jyb, kzrbd)) &
                &/ dz * 4.0 / (jac(ixr, jyb, kzrbd) + jac(ixr + 1, jyb, kzrbd) &
                &+ jac(ixr, jyb, kzrbd + 1) + jac(ixr + 1, jyb, kzrbd + 1))
            flwrbu = 0.0
          else
            flwrbd = 0.0
            flwrbu = 0.0
          end if
        end if

        if(zrfu < topography_surface(ixr, jyf)) then
          flwrfd = 0.0
          flwrfu = 0.0
        else if(zrfd < topography_surface(ixr, jyf)) then
          flwrfd = 0.0
          flwrfu = (var%u(ixr, jyf, kzrfu + 1) - var%u(ixr, jyf, kzrfu)) / dz &
              &* 4.0 / (jac(ixr, jyf, kzrfu) + jac(ixr + 1, jyf, kzrfu) &
              &+ jac(ixr, jyf, kzrfu + 1) + jac(ixr + 1, jyf, kzrfu + 1))
        else
          if(zrfu < lz(1)) then
            flwrfd = (var%u(ixr, jyf, kzrfd + 1) - var%u(ixr, jyf, kzrfd)) &
                &/ dz * 4.0 / (jac(ixr, jyf, kzrfd) + jac(ixr + 1, jyf, kzrfd) &
                &+ jac(ixr, jyf, kzrfd + 1) + jac(ixr + 1, jyf, kzrfd + 1))
            flwrfu = (var%u(ixr, jyf, kzrfu + 1) - var%u(ixr, jyf, kzrfu)) &
                &/ dz * 4.0 / (jac(ixr, jyf, kzrfu) + jac(ixr + 1, jyf, kzrfu) &
                &+ jac(ixr, jyf, kzrfu + 1) + jac(ixr + 1, jyf, kzrfu + 1))
          else if(zrfd < lz(1)) then
            flwrfd = (var%u(ixr, jyf, kzrfd + 1) - var%u(ixr, jyf, kzrfd)) &
                &/ dz * 4.0 / (jac(ixr, jyf, kzrfd) + jac(ixr + 1, jyf, kzrfd) &
                &+ jac(ixr, jyf, kzrfd + 1) + jac(ixr + 1, jyf, kzrfd + 1))
            flwrfu = 0.0
          else
            flwrfd = 0.0
            flwrfu = 0.0
          end if
        end if
      else
        ! for variable at intermediate levels
        ! z(k+1/2) = lz(0) + k*dz:

        ! index for lower level used for linear interpolation
        kzd = max(- nbz, floor((zlc - lz(0)) / dz))

        ! index for upper level used for linear interpolation
        kzu = kzd + 1

        ! if the ray volume is above the model domain make sure that the
        ! interpolation below leads to a zero gradient
        if(kzd > nz) then
          kzu = nz + 1
          kzd = nz + 1
        end if

        ! intermediate levels z(k+1/2) = z(k) + 0.5*dz
        zu = z(kzu) + 0.5 * dz
        zd = z(kzd) + 0.5 * dz

        if(sizeY == 1) then
          jyb = 1
          jyf = 1
        else
          ! for variable at full levels y(j) = ly(0) +  0.5*dy + (j-1)*dy:

          ! index for backward level used for linear interpolation
          jyb = floor((ylc - 0.5 * dy - ly(0)) / dy) + 1 - jy0

          if(jyb < - nby) then
            print *, 'ERROR IN MEANFLOW: jyb =', jyb, '< -nby =', - nby
            stop
          end if

          ! index for forward level used for linear interpolation
          jyf = jyb + 1

          if(jyf > ny + nby) then
            print *, 'ERROR IN MEANFLOW: jyf =', jyf, '> ny + nby =', ny + nby
            stop
          end if
        end if

        ! full levels y(j)
        yf = y(jyf + jy0)
        yb = y(jyb + jy0)

        if(sizeX == 1) then
          ixl = 1
          ixr = 1
        else
          ! for variable at intermediate levels
          ! x(i+1/2) = x(i) + 0.5*dx = lx(0) + i*dx:

          ! index for leftmost level used for linear interpolation
          ixl = floor((xlc - lx(0)) / dx) - ix0

          if(ixl < - nbx) then
            print *, 'ERROR IN MEANFLOW: ixl =', ixl, '< -nbx =', - nbx
            ! testb
            print *, 'variable = du/dz'
            print *, 'xlc =', xlc
            print *, 'lx(0) =', lx(0)
            print *, 'lx(1) =', lx(1)
            print *, 'dx =', dx
            print *, 'lRef =', lRef
            ! teste
            stop
          end if

          ! index for rightmost level used for linear interpolation
          ixr = ixl + 1

          if(ixr > nx + nbx) then
            print *, 'ERROR IN MEANFLOW: ixr =', ixr, '> nx + nbx =', nx + nbx
            stop
          end if
        end if

        ! intermediate levels x(i+1/2) = x(i) + 0.5*dx
        xr = x(ixr + ix0) + 0.5 * dx
        xl = x(ixl + ix0) + 0.5 * dx

        ! values of var. at the eight corners of the interpolation region
        ! du/dz (i+1/2,j,k+1/2)) = (u(i+1/2,j,k+1) - u(i+1/2,j,k))/dz, hence
        ! (using 0 at levels outside of the model domain)
        if(zu < lz(0)) then
          flwlbd = 0.0
          flwlbu = 0.0

          flwlfd = 0.0
          flwlfu = 0.0

          flwrbd = 0.0
          flwrbu = 0.0

          flwrfd = 0.0
          flwrfu = 0.0
        elseif(zd < lz(0)) then
          flwlbd = 0.0
          flwlbu = (var%u(ixl, jyb, kzu + 1) - var%u(ixl, jyb, kzu)) / dz

          flwlfd = 0.0
          flwlfu = (var%u(ixl, jyf, kzu + 1) - var%u(ixl, jyf, kzu)) / dz

          flwrbd = 0.0
          flwrbu = (var%u(ixr, jyb, kzu + 1) - var%u(ixr, jyb, kzu)) / dz

          flwrfd = 0.0
          flwrfu = (var%u(ixr, jyf, kzu + 1) - var%u(ixr, jyf, kzu)) / dz
        else
          if(zu < lz(1)) then
            flwlbd = (var%u(ixl, jyb, kzd + 1) - var%u(ixl, jyb, kzd)) / dz
            flwlbu = (var%u(ixl, jyb, kzu + 1) - var%u(ixl, jyb, kzu)) / dz

            flwlfd = (var%u(ixl, jyf, kzd + 1) - var%u(ixl, jyf, kzd)) / dz
            flwlfu = (var%u(ixl, jyf, kzu + 1) - var%u(ixl, jyf, kzu)) / dz

            flwrbd = (var%u(ixr, jyb, kzd + 1) - var%u(ixr, jyb, kzd)) / dz
            flwrbu = (var%u(ixr, jyb, kzu + 1) - var%u(ixr, jyb, kzu)) / dz

            flwrfd = (var%u(ixr, jyf, kzd + 1) - var%u(ixr, jyf, kzd)) / dz
            flwrfu = (var%u(ixr, jyf, kzu + 1) - var%u(ixr, jyf, kzu)) / dz
          elseif(zd < lz(1)) then
            flwlbd = (var%u(ixl, jyb, kzd + 1) - var%u(ixl, jyb, kzd)) / dz
            flwlbu = 0.0

            flwlfd = (var%u(ixl, jyf, kzd + 1) - var%u(ixl, jyf, kzd)) / dz
            flwlfu = 0.0

            flwrbd = (var%u(ixr, jyb, kzd + 1) - var%u(ixr, jyb, kzd)) / dz
            flwrbu = 0.0

            flwrfd = (var%u(ixr, jyf, kzd + 1) - var%u(ixr, jyf, kzd)) / dz
            flwrfu = 0.0
          else
            flwlbd = 0.0
            flwlbu = 0.0

            flwlfd = 0.0
            flwlfu = 0.0

            flwrbd = 0.0
            flwrbu = 0.0

            flwrfd = 0.0
            flwrfu = 0.0
          end if
        end if
      end if
    elseif(flwtpe == 21) then
      ! interpolate dv/dx using staggered-grid distribution
      ! dv/dx(i+1/2,j+1/2,k)) = (v(i+1,j+1/2,k) - v(i,j+1/2,k))/dx, hence

      if(topography) then
        ! FJApr2023
        ! Locate the closest points in zonal direction.
        if(sizeX == 1) then
          flw = 0.0
          return
        else
          ixl = floor((xlc - lx(0)) / dx) - ix0
          if(ixl < - nbx) then
            print *, "ERROR IN MEANFLOW: ixl =", ixl, "< - nbx =", - nbx
            stop
          end if
          ixr = ixl + 1
          if(ixr + 1 > nx + nbx) then
            print *, "ERROR IN MEANFLOW: ixr + 1 =", ixr + 1, "> nx + nbx =", &
                &nx + nbx
            stop
          end if
        end if
        xr = x(ixr + ix0) + 0.5 * dx
        xl = x(ixl + ix0) + 0.5 * dx

        ! Locate the closest points in meridional direction.
        if(sizeY == 1) then
          jyb = 1
          jyf = 1
        else
          jyb = floor((ylc - ly(0)) / dy) - jy0
          if(jyb < - nby) then
            print *, "ERROR IN MEANFLOW: jyb =", jyb, "< - nby =", - nby
            stop
          end if
          jyf = jyb + 1
          if(jyf + 1 > ny + nby) then
            print *, "ERROR IN MEANFLOW: jyf + 1 =", jyf + 1, "> ny + nby =", &
                &ny + nby
            stop
          end if
        end if
        yf = y(jyf + jy0) + 0.5 * dy
        yb = y(jyb + jy0) + 0.5 * dy

        ! Locate the closest points in vertical direction.

        kzlbd = max(1, floor((levelTFC(ixl, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzlbu = max(1, floor((levelTFC(ixl, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzlbd > nz) then
          kzlbu = nz
          kzlbd = nz
        end if
        zlbd = zTFC(ixl, jyb, kzlbd)
        zlbu = zTFC(ixl, jyb, kzlbu)

        kzlfd = max(1, floor((levelTFC(ixl, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzlfu = max(1, floor((levelTFC(ixl, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzlfd > nz) then
          kzlfu = nz
          kzlfd = nz
        end if
        zlfd = zTFC(ixl, jyf, kzlfd)
        zlfu = zTFC(ixl, jyf, kzlfu)

        kzrbd = max(1, floor((levelTFC(ixr, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzrbu = max(1, floor((levelTFC(ixr, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzrbd > nz) then
          kzrbu = nz
          kzrbd = nz
        end if
        zrbd = zTFC(ixr, jyb, kzrbd)
        zrbu = zTFC(ixr, jyb, kzrbu)

        kzrfd = max(1, floor((levelTFC(ixr, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzrfu = max(1, floor((levelTFC(ixr, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzrfd > nz) then
          kzrfu = nz
          kzrfd = nz
        end if
        zrfd = zTFC(ixr, jyf, kzrfd)
        zrfu = zTFC(ixr, jyf, kzrfu)

        ! Assign the values.

        flwlbd = (0.5 * ((jac(ixl + 1, jyb, kzlbd) + jac(ixl + 1, jyb + 1, &
            &kzlbd)) * var%v(ixl + 1, jyb, kzlbd) - (jac(ixl, jyb, kzlbd) &
            &+ jac(ixl, jyb + 1, kzlbd)) * var%v(ixl, jyb, kzlbd)) / dx &
            &+ 0.0625 * ((jac(ixl, jyb, kzlbd + 1) * met(ixl, jyb, kzlbd + 1, &
            &1, 3) + jac(ixl + 1, jyb, kzlbd + 1) * met(ixl + 1, jyb, kzlbd &
            &+ 1, 1, 3) + jac(ixl, jyb + 1, kzlbd + 1) * met(ixl, jyb + 1, &
            &kzlbd + 1, 1, 3) + jac(ixl + 1, jyb + 1, kzlbd + 1) * met(ixl &
            &+ 1, jyb + 1, kzlbd + 1, 1, 3)) * (var%v(ixl, jyb, kzlbd + 1) &
            &+ var%v(ixl + 1, jyb, kzlbd + 1)) - (jac(ixl, jyb, kzlbd - 1) &
            &* met(ixl, jyb, kzlbd - 1, 1, 3) + jac(ixl + 1, jyb, kzlbd - 1) &
            &* met(ixl + 1, jyb, kzlbd - 1, 1, 3) + jac(ixl, jyb + 1, kzlbd &
            &- 1) * met(ixl, jyb + 1, kzlbd - 1, 1, 3) + jac(ixl + 1, jyb + 1, &
            &kzlbd - 1) * met(ixl + 1, jyb + 1, kzlbd - 1, 1, 3)) &
            &* (var%v(ixl, jyb, kzlbd - 1) + var%v(ixl + 1, jyb, kzlbd - 1))) &
            &/ dz) * 4.0 / (jac(ixl, jyb, kzlbd) + jac(ixl + 1, jyb, kzlbd) &
            &+ jac(ixl, jyb + 1, kzlbd) + jac(ixl + 1, jyb + 1, kzlbd))
        flwlbu = (0.5 * ((jac(ixl + 1, jyb, kzlbu) + jac(ixl + 1, jyb + 1, &
            &kzlbu)) * var%v(ixl + 1, jyb, kzlbu) - (jac(ixl, jyb, kzlbu) &
            &+ jac(ixl, jyb + 1, kzlbu)) * var%v(ixl, jyb, kzlbu)) / dx &
            &+ 0.0625 * ((jac(ixl, jyb, kzlbu + 1) * met(ixl, jyb, kzlbu + 1, &
            &1, 3) + jac(ixl + 1, jyb, kzlbu + 1) * met(ixl + 1, jyb, kzlbu &
            &+ 1, 1, 3) + jac(ixl, jyb + 1, kzlbu + 1) * met(ixl, jyb + 1, &
            &kzlbu + 1, 1, 3) + jac(ixl + 1, jyb + 1, kzlbu + 1) * met(ixl &
            &+ 1, jyb + 1, kzlbu + 1, 1, 3)) * (var%v(ixl, jyb, kzlbu + 1) &
            &+ var%v(ixl + 1, jyb, kzlbu + 1)) - (jac(ixl, jyb, kzlbu - 1) &
            &* met(ixl, jyb, kzlbu - 1, 1, 3) + jac(ixl + 1, jyb, kzlbu - 1) &
            &* met(ixl + 1, jyb, kzlbu - 1, 1, 3) + jac(ixl, jyb + 1, kzlbu &
            &- 1) * met(ixl, jyb + 1, kzlbu - 1, 1, 3) + jac(ixl + 1, jyb + 1, &
            &kzlbu - 1) * met(ixl + 1, jyb + 1, kzlbu - 1, 1, 3)) &
            &* (var%v(ixl, jyb, kzlbu - 1) + var%v(ixl + 1, jyb, kzlbu - 1))) &
            &/ dz) * 4.0 / (jac(ixl, jyb, kzlbu) + jac(ixl + 1, jyb, kzlbu) &
            &+ jac(ixl, jyb + 1, kzlbu) + jac(ixl + 1, jyb + 1, kzlbu))

        flwlfd = (0.5 * ((jac(ixl + 1, jyf, kzlfd) + jac(ixl + 1, jyf + 1, &
            &kzlfd)) * var%v(ixl + 1, jyf, kzlfd) - (jac(ixl, jyf, kzlfd) &
            &+ jac(ixl, jyf + 1, kzlfd)) * var%v(ixl, jyf, kzlfd)) / dx &
            &+ 0.0625 * ((jac(ixl, jyf, kzlfd + 1) * met(ixl, jyf, kzlfd + 1, &
            &1, 3) + jac(ixl + 1, jyf, kzlfd + 1) * met(ixl + 1, jyf, kzlfd &
            &+ 1, 1, 3) + jac(ixl, jyf + 1, kzlfd + 1) * met(ixl, jyf + 1, &
            &kzlfd + 1, 1, 3) + jac(ixl + 1, jyf + 1, kzlfd + 1) * met(ixl &
            &+ 1, jyf + 1, kzlfd + 1, 1, 3)) * (var%v(ixl, jyf, kzlfd + 1) &
            &+ var%v(ixl + 1, jyf, kzlfd + 1)) - (jac(ixl, jyf, kzlfd - 1) &
            &* met(ixl, jyf, kzlfd - 1, 1, 3) + jac(ixl + 1, jyf, kzlfd - 1) &
            &* met(ixl + 1, jyf, kzlfd - 1, 1, 3) + jac(ixl, jyf + 1, kzlfd &
            &- 1) * met(ixl, jyf + 1, kzlfd - 1, 1, 3) + jac(ixl + 1, jyf + 1, &
            &kzlfd - 1) * met(ixl + 1, jyf + 1, kzlfd - 1, 1, 3)) &
            &* (var%v(ixl, jyf, kzlfd - 1) + var%v(ixl + 1, jyf, kzlfd - 1))) &
            &/ dz) * 4.0 / (jac(ixl, jyf, kzlfd) + jac(ixl + 1, jyf, kzlfd) &
            &+ jac(ixl, jyf + 1, kzlfd) + jac(ixl + 1, jyf + 1, kzlfd))
        flwlfu = (0.5 * ((jac(ixl + 1, jyf, kzlfu) + jac(ixl + 1, jyf + 1, &
            &kzlfu)) * var%v(ixl + 1, jyf, kzlfu) - (jac(ixl, jyf, kzlfu) &
            &+ jac(ixl, jyf + 1, kzlfu)) * var%v(ixl, jyf, kzlfu)) / dx &
            &+ 0.0625 * ((jac(ixl, jyf, kzlfu + 1) * met(ixl, jyf, kzlfu + 1, &
            &1, 3) + jac(ixl + 1, jyf, kzlfu + 1) * met(ixl + 1, jyf, kzlfu &
            &+ 1, 1, 3) + jac(ixl, jyf + 1, kzlfu + 1) * met(ixl, jyf + 1, &
            &kzlfu + 1, 1, 3) + jac(ixl + 1, jyf + 1, kzlfu + 1) * met(ixl &
            &+ 1, jyf + 1, kzlfu + 1, 1, 3)) * (var%v(ixl, jyf, kzlfu + 1) &
            &+ var%v(ixl + 1, jyf, kzlfu + 1)) - (jac(ixl, jyf, kzlfu - 1) &
            &* met(ixl, jyf, kzlfu - 1, 1, 3) + jac(ixl + 1, jyf, kzlfu - 1) &
            &* met(ixl + 1, jyf, kzlfu - 1, 1, 3) + jac(ixl, jyf + 1, kzlfu &
            &- 1) * met(ixl, jyf + 1, kzlfu - 1, 1, 3) + jac(ixl + 1, jyf + 1, &
            &kzlfu - 1) * met(ixl + 1, jyf + 1, kzlfu - 1, 1, 3)) &
            &* (var%v(ixl, jyf, kzlfu - 1) + var%v(ixl + 1, jyf, kzlfu - 1))) &
            &/ dz) * 4.0 / (jac(ixl, jyf, kzlfu) + jac(ixl + 1, jyf, kzlfu) &
            &+ jac(ixl, jyf + 1, kzlfu) + jac(ixl + 1, jyf + 1, kzlfu))

        flwrbd = (0.5 * ((jac(ixr + 1, jyb, kzrbd) + jac(ixr + 1, jyb + 1, &
            &kzrbd)) * var%v(ixr + 1, jyb, kzrbd) - (jac(ixr, jyb, kzrbd) &
            &+ jac(ixr, jyb + 1, kzrbd)) * var%v(ixr, jyb, kzrbd)) / dx &
            &+ 0.0625 * ((jac(ixr, jyb, kzrbd + 1) * met(ixr, jyb, kzrbd + 1, &
            &1, 3) + jac(ixr + 1, jyb, kzrbd + 1) * met(ixr + 1, jyb, kzrbd &
            &+ 1, 1, 3) + jac(ixr, jyb + 1, kzrbd + 1) * met(ixr, jyb + 1, &
            &kzrbd + 1, 1, 3) + jac(ixr + 1, jyb + 1, kzrbd + 1) * met(ixr &
            &+ 1, jyb + 1, kzrbd + 1, 1, 3)) * (var%v(ixr, jyb, kzrbd + 1) &
            &+ var%v(ixr + 1, jyb, kzrbd + 1)) - (jac(ixr, jyb, kzrbd - 1) &
            &* met(ixr, jyb, kzrbd - 1, 1, 3) + jac(ixr + 1, jyb, kzrbd - 1) &
            &* met(ixr + 1, jyb, kzrbd - 1, 1, 3) + jac(ixr, jyb + 1, kzrbd &
            &- 1) * met(ixr, jyb + 1, kzrbd - 1, 1, 3) + jac(ixr + 1, jyb + 1, &
            &kzrbd - 1) * met(ixr + 1, jyb + 1, kzrbd - 1, 1, 3)) &
            &* (var%v(ixr, jyb, kzrbd - 1) + var%v(ixr + 1, jyb, kzrbd - 1))) &
            &/ dz) * 4.0 / (jac(ixr, jyb, kzrbd) + jac(ixr + 1, jyb, kzrbd) &
            &+ jac(ixr, jyb + 1, kzrbd) + jac(ixr + 1, jyb + 1, kzrbd))
        flwrbu = (0.5 * ((jac(ixr + 1, jyb, kzrbu) + jac(ixr + 1, jyb + 1, &
            &kzrbu)) * var%v(ixr + 1, jyb, kzrbu) - (jac(ixr, jyb, kzrbu) &
            &+ jac(ixr, jyb + 1, kzrbu)) * var%v(ixr, jyb, kzrbu)) / dx &
            &+ 0.0625 * ((jac(ixr, jyb, kzrbu + 1) * met(ixr, jyb, kzrbu + 1, &
            &1, 3) + jac(ixr + 1, jyb, kzrbu + 1) * met(ixr + 1, jyb, kzrbu &
            &+ 1, 1, 3) + jac(ixr, jyb + 1, kzrbu + 1) * met(ixr, jyb + 1, &
            &kzrbu + 1, 1, 3) + jac(ixr + 1, jyb + 1, kzrbu + 1) * met(ixr &
            &+ 1, jyb + 1, kzrbu + 1, 1, 3)) * (var%v(ixr, jyb, kzrbu + 1) &
            &+ var%v(ixr + 1, jyb, kzrbu + 1)) - (jac(ixr, jyb, kzrbu - 1) &
            &* met(ixr, jyb, kzrbu - 1, 1, 3) + jac(ixr + 1, jyb, kzrbu - 1) &
            &* met(ixr + 1, jyb, kzrbu - 1, 1, 3) + jac(ixr, jyb + 1, kzrbu &
            &- 1) * met(ixr, jyb + 1, kzrbu - 1, 1, 3) + jac(ixr + 1, jyb + 1, &
            &kzrbu - 1) * met(ixr + 1, jyb + 1, kzrbu - 1, 1, 3)) &
            &* (var%v(ixr, jyb, kzrbu - 1) + var%v(ixr + 1, jyb, kzrbu - 1))) &
            &/ dz) * 4.0 / (jac(ixr, jyb, kzrbu) + jac(ixr + 1, jyb, kzrbu) &
            &+ jac(ixr, jyb + 1, kzrbu) + jac(ixr + 1, jyb + 1, kzrbu))

        flwrfd = (0.5 * ((jac(ixr + 1, jyf, kzrfd) + jac(ixr + 1, jyf + 1, &
            &kzrfd)) * var%v(ixr + 1, jyf, kzrfd) - (jac(ixr, jyf, kzrfd) &
            &+ jac(ixr, jyf + 1, kzrfd)) * var%v(ixr, jyf, kzrfd)) / dx &
            &+ 0.0625 * ((jac(ixr, jyf, kzrfd + 1) * met(ixr, jyf, kzrfd + 1, &
            &1, 3) + jac(ixr + 1, jyf, kzrfd + 1) * met(ixr + 1, jyf, kzrfd &
            &+ 1, 1, 3) + jac(ixr, jyf + 1, kzrfd + 1) * met(ixr, jyf + 1, &
            &kzrfd + 1, 1, 3) + jac(ixr + 1, jyf + 1, kzrfd + 1) * met(ixr &
            &+ 1, jyf + 1, kzrfd + 1, 1, 3)) * (var%v(ixr, jyf, kzrfd + 1) &
            &+ var%v(ixr + 1, jyf, kzrfd + 1)) - (jac(ixr, jyf, kzrfd - 1) &
            &* met(ixr, jyf, kzrfd - 1, 1, 3) + jac(ixr + 1, jyf, kzrfd - 1) &
            &* met(ixr + 1, jyf, kzrfd - 1, 1, 3) + jac(ixr, jyf + 1, kzrfd &
            &- 1) * met(ixr, jyf + 1, kzrfd - 1, 1, 3) + jac(ixr + 1, jyf + 1, &
            &kzrfd - 1) * met(ixr + 1, jyf + 1, kzrfd - 1, 1, 3)) &
            &* (var%v(ixr, jyf, kzrfd - 1) + var%v(ixr + 1, jyf, kzrfd - 1))) &
            &/ dz) * 4.0 / (jac(ixr, jyf, kzrfd) + jac(ixr + 1, jyf, kzrfd) &
            &+ jac(ixr, jyf + 1, kzrfd) + jac(ixr + 1, jyf + 1, kzrfd))
        flwrfu = (0.5 * ((jac(ixr + 1, jyf, kzrfu) + jac(ixr + 1, jyf + 1, &
            &kzrfu)) * var%v(ixr + 1, jyf, kzrfu) - (jac(ixr, jyf, kzrfu) &
            &+ jac(ixr, jyf + 1, kzrfu)) * var%v(ixr, jyf, kzrfu)) / dx &
            &+ 0.0625 * ((jac(ixr, jyf, kzrfu + 1) * met(ixr, jyf, kzrfu + 1, &
            &1, 3) + jac(ixr + 1, jyf, kzrfu + 1) * met(ixr + 1, jyf, kzrfu &
            &+ 1, 1, 3) + jac(ixr, jyf + 1, kzrfu + 1) * met(ixr, jyf + 1, &
            &kzrfu + 1, 1, 3) + jac(ixr + 1, jyf + 1, kzrfu + 1) * met(ixr &
            &+ 1, jyf + 1, kzrfu + 1, 1, 3)) * (var%v(ixr, jyf, kzrfu + 1) &
            &+ var%v(ixr + 1, jyf, kzrfu + 1)) - (jac(ixr, jyf, kzrfu - 1) &
            &* met(ixr, jyf, kzrfu - 1, 1, 3) + jac(ixr + 1, jyf, kzrfu - 1) &
            &* met(ixr + 1, jyf, kzrfu - 1, 1, 3) + jac(ixr, jyf + 1, kzrfu &
            &- 1) * met(ixr, jyf + 1, kzrfu - 1, 1, 3) + jac(ixr + 1, jyf + 1, &
            &kzrfu - 1) * met(ixr + 1, jyf + 1, kzrfu - 1, 1, 3)) &
            &* (var%v(ixr, jyf, kzrfu - 1) + var%v(ixr + 1, jyf, kzrfu - 1))) &
            &/ dz) * 4.0 / (jac(ixr, jyf, kzrfu) + jac(ixr + 1, jyf, kzrfu) &
            &+ jac(ixr, jyf + 1, kzrfu) + jac(ixr + 1, jyf + 1, kzrfu))
      else
        if(sizeX == 1) then
          ! no derivative if there is no x dependence
          flw = 0.0

          return
        else
          ! for variable at intermediate levels
          ! x(i+1/2) = x(i) + 0.5*dx = lx(0) + i*dx:

          ! index for leftmost level used for linear interpolation
          ixl = floor((xlc - lx(0)) / dx) - ix0

          if(ixl < - nbx) then
            print *, 'ERROR IN MEANFLOW: ixl =', ixl, '< -nbx =', - nbx
            ! testb
            print *, 'variable = dv/dx'
            print *, 'xlc =', xlc
            print *, 'lx(0) =', lx(0)
            print *, 'lx(1) =', lx(1)
            print *, 'dx =', dx
            print *, 'lRef =', lRef
            ! teste
            stop
          end if

          ! index for rightmost level used for linear interpolation
          ixr = ixl + 1

          if(ixr + 1 > nx + nbx) then
            print *, 'MEANFLOW: ixr + 1 =', ixr + 1, '> nx + nbx =', nx + nbx
            stop
          end if
        end if

        ! intermediate levels x(i+1/2) = x(i) + 0.5*dx
        xr = x(ixr + ix0) + 0.5 * dx
        xl = x(ixl + ix0) + 0.5 * dx

        ! for variable at full levels z(k) = lz(0) +  0.5*dz + (k-1)*dz:
        ! (levels below the model bottom are replaced by the first level
        ! in the model domain)

        ! index for lower level used for linear interpolation
        kzd = max(1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 1)

        ! index for upper level used for linear interpolation
        kzu = max(1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 2)

        if(kzd > nz) then
          kzu = nz
          kzd = nz
        end if

        ! full levels z(k)
        zu = z(kzu)
        zd = z(kzd)

        if(sizeY == 1) then
          jyb = 1
          jyf = 1
        else
          ! for variable at intermediate levels
          ! y(j+1/2) = y(j) + 0.5*dy = ly(0) +  j*dy:

          ! index for backward level used for linear interpolation
          jyb = floor((ylc - ly(0)) / dy) - jy0

          if(jyb < - nby) then
            print *, 'ERROR IN MEANFLOW: jyb =', jyb, '< -nby =', - nby
            stop
          end if

          ! index for forward level used for linear interpolation
          jyf = jyb + 1

          if(jyf + 1 > ny + nby) then
            print *, 'MEANFLOW: jyf + 1 =', jyf + 1, '> ny + nby =', ny + nby
            stop
          end if
        end if

        ! intermediate levels y(j+1/2) = y(j) + 0.5*dy
        yf = y(jyf + jy0) + 0.5 * dy
        yb = y(jyb + jy0) + 0.5 * dy

        ! values of var. at the eight corners of the interpolation region
        ! dv/dx (i+1/2,j+1/2,k)) = (v(i+1,j+1/2,k) - v(i,j+1/2,k))/dx, hence
        flwlbd = (var%v(ixl + 1, jyb, kzd) - var%v(ixl, jyb, kzd)) / dx
        flwlbu = (var%v(ixl + 1, jyb, kzu) - var%v(ixl, jyb, kzu)) / dx

        flwlfd = (var%v(ixl + 1, jyf, kzd) - var%v(ixl, jyf, kzd)) / dx
        flwlfu = (var%v(ixl + 1, jyf, kzu) - var%v(ixl, jyf, kzu)) / dx

        flwrbd = (var%v(ixr + 1, jyb, kzd) - var%v(ixr, jyb, kzd)) / dx
        flwrbu = (var%v(ixr + 1, jyb, kzu) - var%v(ixr, jyb, kzu)) / dx

        flwrfd = (var%v(ixr + 1, jyf, kzd) - var%v(ixr, jyf, kzd)) / dx
        flwrfu = (var%v(ixr + 1, jyf, kzu) - var%v(ixr, jyf, kzu)) / dx
      end if
    elseif(flwtpe == 22) then
      ! interpolate dv/dy using staggered-grid distribution
      ! dv/dy (i,j,k)) = (v(i,j+1/2,k) - v(i,j-1/2,k))/dy, hence

      if(topography) then
        ! FJApr2023
        ! Locate the closest points in zonal direction.
        if(sizeX == 1) then
          ixl = 1
          ixr = 1
        else
          ixl = floor((xlc - 0.5 * dx - lx(0)) / dx) + 1 - ix0
          if(ixl < - nbx) then
            print *, "ERROR IN MEANFLOW: ixl =", ixl, "< - nbx =", - nbx
            stop
          end if
          ixr = ixl + 1
          if(ixr > nx + nbx) then
            print *, "ERROR IN MEANFLOW: ixr =", ixr, "> nx + nbx =", nx + nbx
            stop
          end if
        end if
        xr = x(ixr + ix0)
        xl = x(ixl + ix0)

        ! Locate the closest points in meridional direction.
        if(sizeY == 1) then
          flw = 0.0
          return
        else
          jyb = floor((ylc - 0.5 * dy - ly(0)) / dy) + 1 - jy0
          if(jyb - 1 < - nby) then
            print *, "ERROR IN MEANFLOW: jyb - 1 =", jyb - 1, "< - nby =", - nby
            stop
          end if
          jyf = jyb + 1
          if(jyf > ny + nby) then
            print *, "ERROR IN MEANFLOW: jyf =", jyf, "> ny + nby =", ny + nby
            stop
          end if
        end if
        yf = y(jyf + jy0)
        yb = y(jyb + jy0)

        ! Locate the closest points in vertical direction.

        kzlbd = max(1, floor((levelTFC(ixl, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzlbu = max(1, floor((levelTFC(ixl, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzlbd > nz) then
          kzlbu = nz
          kzlbd = nz
        end if
        zlbd = zTFC(ixl, jyb, kzlbd)
        zlbu = zTFC(ixl, jyb, kzlbu)

        kzlfd = max(1, floor((levelTFC(ixl, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzlfu = max(1, floor((levelTFC(ixl, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzlfd > nz) then
          kzlfu = nz
          kzlfd = nz
        end if
        zlfd = zTFC(ixl, jyf, kzlfd)
        zlfu = zTFC(ixl, jyf, kzlfu)

        kzrbd = max(1, floor((levelTFC(ixr, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzrbu = max(1, floor((levelTFC(ixr, jyb, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzrbd > nz) then
          kzrbu = nz
          kzrbd = nz
        end if
        zrbd = zTFC(ixr, jyb, kzrbd)
        zrbu = zTFC(ixr, jyb, kzrbu)

        kzrfd = max(1, floor((levelTFC(ixr, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 1)
        kzrfu = max(1, floor((levelTFC(ixr, jyf, zlc) - 0.5 * dz - lz(0)) &
            &/ dz) + 2)
        if(kzrfd > nz) then
          kzrfu = nz
          kzrfd = nz
        end if
        zrfd = zTFC(ixr, jyf, kzrfd)
        zrfu = zTFC(ixr, jyf, kzrfu)

        ! Assign the values.

        flwlbd = (0.5 * ((jac(ixl, jyb, kzlbd) + jac(ixl, jyb + 1, kzlbd)) &
            &* var%v(ixl, jyb, kzlbd) - (jac(ixl, jyb, kzlbd) + jac(ixl, jyb &
            &- 1, kzlbd)) * var%v(ixl, jyb - 1, kzlbd)) / dy + 0.25 &
            &* (jac(ixl, jyb, kzlbd + 1) * met(ixl, jyb, kzlbd + 1, 2, 3) &
            &* (var%v(ixl, jyb, kzlbd + 1) + var%v(ixl, jyb - 1, kzlbd + 1)) &
            &- jac(ixl, jyb, kzlbd - 1) * met(ixl, jyb, kzlbd - 1, 2, 3) &
            &* (var%v(ixl, jyb, kzlbd - 1) + var%v(ixl, jyb - 1, kzlbd - 1))) &
            &/ dz) / jac(ixl, jyb, kzlbd)
        flwlbu = (0.5 * ((jac(ixl, jyb, kzlbu) + jac(ixl, jyb + 1, kzlbu)) &
            &* var%v(ixl, jyb, kzlbu) - (jac(ixl, jyb, kzlbu) + jac(ixl, jyb &
            &- 1, kzlbu)) * var%v(ixl, jyb - 1, kzlbu)) / dy + 0.25 &
            &* (jac(ixl, jyb, kzlbu + 1) * met(ixl, jyb, kzlbu + 1, 2, 3) &
            &* (var%v(ixl, jyb, kzlbu + 1) + var%v(ixl, jyb - 1, kzlbu + 1)) &
            &- jac(ixl, jyb, kzlbu - 1) * met(ixl, jyb, kzlbu - 1, 2, 3) &
            &* (var%v(ixl, jyb, kzlbu - 1) + var%v(ixl, jyb - 1, kzlbu - 1))) &
            &/ dz) / jac(ixl, jyb, kzlbu)

        flwlfd = (0.5 * ((jac(ixl, jyf, kzlfd) + jac(ixl, jyf + 1, kzlfd)) &
            &* var%v(ixl, jyf, kzlfd) - (jac(ixl, jyf, kzlfd) + jac(ixl, jyf &
            &- 1, kzlfd)) * var%v(ixl, jyf - 1, kzlfd)) / dy + 0.25 &
            &* (jac(ixl, jyf, kzlfd + 1) * met(ixl, jyf, kzlfd + 1, 2, 3) &
            &* (var%v(ixl, jyf, kzlfd + 1) + var%v(ixl, jyf - 1, kzlfd + 1)) &
            &- jac(ixl, jyf, kzlfd - 1) * met(ixl, jyf, kzlfd - 1, 2, 3) &
            &* (var%v(ixl, jyf, kzlfd - 1) + var%v(ixl, jyf - 1, kzlfd - 1))) &
            &/ dz) / jac(ixl, jyf, kzlfd)
        flwlfu = (0.5 * ((jac(ixl, jyf, kzlfu) + jac(ixl, jyf + 1, kzlfu)) &
            &* var%v(ixl, jyf, kzlfu) - (jac(ixl, jyf, kzlfu) + jac(ixl, jyf &
            &- 1, kzlfu)) * var%v(ixl, jyf - 1, kzlfu)) / dy + 0.25 &
            &* (jac(ixl, jyf, kzlfu + 1) * met(ixl, jyf, kzlfu + 1, 2, 3) &
            &* (var%v(ixl, jyf, kzlfu + 1) + var%v(ixl, jyf - 1, kzlfu + 1)) &
            &- jac(ixl, jyf, kzlfu - 1) * met(ixl, jyf, kzlfu - 1, 2, 3) &
            &* (var%v(ixl, jyf, kzlfu - 1) + var%v(ixl, jyf - 1, kzlfu - 1))) &
            &/ dz) / jac(ixl, jyf, kzlfu)

        flwrbd = (0.5 * ((jac(ixr, jyb, kzrbd) + jac(ixr, jyb + 1, kzrbd)) &
            &* var%v(ixr, jyb, kzrbd) - (jac(ixr, jyb, kzrbd) + jac(ixr, jyb &
            &- 1, kzrbd)) * var%v(ixr, jyb - 1, kzrbd)) / dy + 0.25 &
            &* (jac(ixr, jyb, kzrbd + 1) * met(ixr, jyb, kzrbd + 1, 2, 3) &
            &* (var%v(ixr, jyb, kzrbd + 1) + var%v(ixr, jyb - 1, kzrbd + 1)) &
            &- jac(ixr, jyb, kzrbd - 1) * met(ixr, jyb, kzrbd - 1, 2, 3) &
            &* (var%v(ixr, jyb, kzrbd - 1) + var%v(ixr, jyb - 1, kzrbd - 1))) &
            &/ dz) / jac(ixr, jyb, kzrbd)
        flwrbu = (0.5 * ((jac(ixr, jyb, kzrbu) + jac(ixr, jyb + 1, kzrbu)) &
            &* var%v(ixr, jyb, kzrbu) - (jac(ixr, jyb, kzrbu) + jac(ixr, jyb &
            &- 1, kzrbu)) * var%v(ixr, jyb - 1, kzrbu)) / dy + 0.25 &
            &* (jac(ixr, jyb, kzrbu + 1) * met(ixr, jyb, kzrbu + 1, 2, 3) &
            &* (var%v(ixr, jyb, kzrbu + 1) + var%v(ixr, jyb - 1, kzrbu + 1)) &
            &- jac(ixr, jyb, kzrbu - 1) * met(ixr, jyb, kzrbu - 1, 2, 3) &
            &* (var%v(ixr, jyb, kzrbu - 1) + var%v(ixr, jyb - 1, kzrbu - 1))) &
            &/ dz) / jac(ixr, jyb, kzrbu)

        flwrfd = (0.5 * ((jac(ixr, jyf, kzrfd) + jac(ixr, jyf + 1, kzrfd)) &
            &* var%v(ixr, jyf, kzrfd) - (jac(ixr, jyf, kzrfd) + jac(ixr, jyf &
            &- 1, kzrfd)) * var%v(ixr, jyf - 1, kzrfd)) / dy + 0.25 &
            &* (jac(ixr, jyf, kzrfd + 1) * met(ixr, jyf, kzrfd + 1, 2, 3) &
            &* (var%v(ixr, jyf, kzrfd + 1) + var%v(ixr, jyf - 1, kzrfd + 1)) &
            &- jac(ixr, jyf, kzrfd - 1) * met(ixr, jyf, kzrfd - 1, 2, 3) &
            &* (var%v(ixr, jyf, kzrfd - 1) + var%v(ixr, jyf - 1, kzrfd - 1))) &
            &/ dz) / jac(ixr, jyf, kzrfd)
        flwrfu = (0.5 * ((jac(ixr, jyf, kzrfu) + jac(ixr, jyf + 1, kzrfu)) &
            &* var%v(ixr, jyf, kzrfu) - (jac(ixr, jyf, kzrfu) + jac(ixr, jyf &
            &- 1, kzrfu)) * var%v(ixr, jyf - 1, kzrfu)) / dy + 0.25 &
            &* (jac(ixr, jyf, kzrfu + 1) * met(ixr, jyf, kzrfu + 1, 2, 3) &
            &* (var%v(ixr, jyf, kzrfu + 1) + var%v(ixr, jyf - 1, kzrfu + 1)) &
            &- jac(ixr, jyf, kzrfu - 1) * met(ixr, jyf, kzrfu - 1, 2, 3) &
            &* (var%v(ixr, jyf, kzrfu - 1) + var%v(ixr, jyf - 1, kzrfu - 1))) &
            &/ dz) / jac(ixr, jyf, kzrfu)
      else
        if(sizeY == 1) then
          ! no derivative if there is no y dependence
          flw = 0.0

          return
        else
          ! for variable at full levels y(j) = ly(0) +  0.5*dy + (j-1)*dy:

          ! index for backward level used for linear interpolation
          jyb = floor((ylc - 0.5 * dy - ly(0)) / dy) + 1 - jy0

          if(jyb - 1 < - nby) then
            print *, 'MEANFLOW: jyb - 1 =', jyb - 1, '< -nby =', - nby
            stop
          end if

          ! index for forward level used for linear interpolation
          jyf = jyb + 1

          if(jyf > ny + nby) then
            print *, 'MEANFLOW: jyf =', jyf, '> ny + nby =', ny + nby
            stop
          end if
        end if

        ! full levels y(j)
        yf = y(jyf + jy0)
        yb = y(jyb + jy0)

        ! for variable at full levels z(k) = lz(0) +  0.5*dz + (k-1)*dz:
        ! (levels below the model bottom are replaced by the first level
        ! in the model domain)

        ! index for lower level used for linear interpolation
        kzd = max(1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 1)

        ! index for upper level used for linear interpolation
        kzu = max(1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 2)

        if(kzd > nz) then
          kzu = nz
          kzd = nz
        end if

        ! full levels z(k)
        zu = z(kzu)
        zd = z(kzd)

        if(sizeX == 1) then
          ixl = 1
          ixr = 1
        else
          ! for variable at full levels x(i) = lx(0) + 0.5*dx + (i-1)*dx

          ! index for leftmost level used for linear interpolation
          ixl = floor((xlc - 0.5 * dx - lx(0)) / dx) + 1 - ix0

          if(ixl < - nbx) then
            print *, 'ERROR IN MEANFLOW: ixl =', ixl, '< -nbx =', - nbx
            ! testb
            print *, 'variable = dv/dy'
            print *, 'xlc =', xlc
            print *, 'lx(0) =', lx(0)
            print *, 'lx(1) =', lx(1)
            print *, 'dx =', dx
            print *, 'lRef =', lRef
            ! teste
            stop
          end if

          ! index for rightmost level used for linear interpolation
          ixr = ixl + 1

          if(ixr > nx + nbx) then
            print *, 'ERROR IN MEANFLOW: ixr =', ixr, '> nx + nbx =', nx + nbx
            stop
          end if
        end if

        ! full levels x(i)
        xr = x(ixr + ix0)
        xl = x(ixl + ix0)

        ! values of var. at the eight corners of the interpolation region
        ! dv/dy (i,j,k)) = (v(i,j+1/2,k) - v(i,j-1/2,k))/dy, hence
        flwlbd = (var%v(ixl, jyb, kzd) - var%v(ixl, jyb - 1, kzd)) / dy
        flwlbu = (var%v(ixl, jyb, kzu) - var%v(ixl, jyb - 1, kzu)) / dy

        flwlfd = (var%v(ixl, jyf, kzd) - var%v(ixl, jyf - 1, kzd)) / dy
        flwlfu = (var%v(ixl, jyf, kzu) - var%v(ixl, jyf - 1, kzu)) / dy

        flwrbd = (var%v(ixr, jyb, kzd) - var%v(ixr, jyb - 1, kzd)) / dy
        flwrbu = (var%v(ixr, jyb, kzu) - var%v(ixr, jyb - 1, kzu)) / dy

        flwrfd = (var%v(ixr, jyf, kzd) - var%v(ixr, jyf - 1, kzd)) / dy
        flwrfu = (var%v(ixr, jyf, kzu) - var%v(ixr, jyf - 1, kzu)) / dy
      end if
    elseif(flwtpe == 23) then
      ! interpolate dv/dz using staggered-grid distribution
      ! dv/dz (i,j+1/2,k+1/2)) = (v(i,j+1/2,k+1) - v(i,j+1/2,k))/dz, hence

      if(topography) then
        ! FJApr2023
        ! Locate the closest points in zonal direction.
        if(sizeX == 1) then
          ixl = 1
          ixr = 1
        else
          ixl = floor((xlc - 0.5 * dx - lx(0)) / dx) + 1 - ix0
          if(ixl < - nbx) then
            print *, "ERROR IN MEANFLOW: ixl =", ixl, "< - nbx =", - nbx
            stop
          end if
          ixr = ixl + 1
          if(ixr > nx + nbx) then
            print *, "ERROR IN MEANFLOW: ixr =", ixr, "> nx + nbx =", nx + nbx
            stop
          end if
        end if
        xr = x(ixr + ix0)
        xl = x(ixl + ix0)

        ! Locate the closest points in meridional direction.
        if(sizeY == 1) then
          jyb = 1
          jyf = 1
        else
          jyb = floor((ylc - ly(0)) / dy) - jy0
          if(jyb < - nby) then
            print *, "ERROR IN MEANFLOW: jyb =", jyb, "< - nby =", - nby
            stop
          end if
          jyf = jyb + 1
          if(jyf > ny + nby) then
            print *, "ERROR IN MEANFLOW: jyf =", jyf, "> ny + nby =", ny + nby
            stop
          end if
        end if
        yf = y(jyf + jy0) + 0.5 * dy
        yb = y(jyb + jy0) + 0.5 * dy

        ! Locate the closest points in vertical direction.

        kzlbd = max(- nbz, floor((levelTFC(ixl, jyb, zlc) - lz(0)) / dz))
        kzlbu = kzlbd + 1
        if(kzlbd > nz) then
          kzlbu = nz + 1
          kzlbd = nz + 1
        end if
        zlbd = zTildeTFC(ixl, jyb, kzlbd)
        zlbu = zTildeTFC(ixl, jyb, kzlbu)

        kzlfd = max(- nbz, floor((levelTFC(ixl, jyf, zlc) - lz(0)) / dz))
        kzlfu = kzlfd + 1
        if(kzlfd > nz) then
          kzlfu = nz + 1
          kzlfd = nz + 1
        end if
        zlfd = zTildeTFC(ixl, jyf, kzlfd)
        zlfu = zTildeTFC(ixl, jyf, kzlfu)

        kzrbd = max(- nbz, floor((levelTFC(ixr, jyb, zlc) - lz(0)) / dz))
        kzrbu = kzrbd + 1
        if(kzrbd > nz) then
          kzrbu = nz + 1
          kzrbd = nz + 1
        end if
        zrbd = zTildeTFC(ixr, jyb, kzrbd)
        zrbu = zTildeTFC(ixr, jyb, kzrbu)

        kzrfd = max(- nbz, floor((levelTFC(ixr, jyf, zlc) - lz(0)) / dz))
        kzrfu = kzrfd + 1
        if(kzrfd > nz) then
          kzrfu = nz + 1
          kzrfd = nz + 1
        end if
        zrfd = zTildeTFC(ixr, jyf, kzrfd)
        zrfu = zTildeTFC(ixr, jyf, kzrfu)

        ! Assign the values.

        if(zlbu < topography_surface(ixl, jyb)) then
          flwlbd = 0.0
          flwlbu = 0.0
        else if(zlbd < topography_surface(ixl, jyb)) then
          flwlbd = 0.0
          flwlbu = (var%v(ixl, jyb, kzlbu + 1) - var%v(ixl, jyb, kzlbu)) / dz &
              &* 4.0 / (jac(ixl, jyb, kzlbu) + jac(ixl, jyb + 1, kzlbu) &
              &+ jac(ixl, jyb, kzlbu + 1) + jac(ixl, jyb + 1, kzlbu + 1))
        else
          if(zlbu < lz(1)) then
            flwlbd = (var%v(ixl, jyb, kzlbd + 1) - var%v(ixl, jyb, kzlbd)) &
                &/ dz * 4.0 / (jac(ixl, jyb, kzlbd) + jac(ixl, jyb + 1, kzlbd) &
                &+ jac(ixl, jyb, kzlbd + 1) + jac(ixl, jyb + 1, kzlbd + 1))
            flwlbu = (var%v(ixl, jyb, kzlbu + 1) - var%v(ixl, jyb, kzlbu)) &
                &/ dz * 4.0 / (jac(ixl, jyb, kzlbu) + jac(ixl, jyb + 1, kzlbu) &
                &+ jac(ixl, jyb, kzlbu + 1) + jac(ixl, jyb + 1, kzlbu + 1))
          else if(zlbd < lz(1)) then
            flwlbd = (var%v(ixl, jyb, kzlbd + 1) - var%v(ixl, jyb, kzlbd)) &
                &/ dz * 4.0 / (jac(ixl, jyb, kzlbd) + jac(ixl, jyb + 1, kzlbd) &
                &+ jac(ixl, jyb, kzlbd + 1) + jac(ixl, jyb + 1, kzlbd + 1))
            flwlbu = 0.0
          else
            flwlbd = 0.0
            flwlbu = 0.0
          end if
        end if

        if(zlfu < topography_surface(ixl, jyf)) then
          flwlfd = 0.0
          flwlfu = 0.0
        else if(zlfd < topography_surface(ixl, jyf)) then
          flwlfd = 0.0
          flwlfu = (var%v(ixl, jyf, kzlfu + 1) - var%v(ixl, jyf, kzlfu)) / dz &
              &* 4.0 / (jac(ixl, jyf, kzlfu) + jac(ixl, jyf + 1, kzlfu) &
              &+ jac(ixl, jyf, kzlfu + 1) + jac(ixl, jyf + 1, kzlfu + 1))
        else
          if(zlfu < lz(1)) then
            flwlfd = (var%v(ixl, jyf, kzlfd + 1) - var%v(ixl, jyf, kzlfd)) &
                &/ dz * 4.0 / (jac(ixl, jyf, kzlfd) + jac(ixl, jyf + 1, kzlfd) &
                &+ jac(ixl, jyf, kzlfd + 1) + jac(ixl, jyf + 1, kzlfd + 1))
            flwlfu = (var%v(ixl, jyf, kzlfu + 1) - var%v(ixl, jyf, kzlfu)) &
                &/ dz * 4.0 / (jac(ixl, jyf, kzlfu) + jac(ixl, jyf + 1, kzlfu) &
                &+ jac(ixl, jyf, kzlfu + 1) + jac(ixl, jyf + 1, kzlfu + 1))
          else if(zlfd < lz(1)) then
            flwlfd = (var%v(ixl, jyf, kzlfd + 1) - var%v(ixl, jyf, kzlfd)) &
                &/ dz * 4.0 / (jac(ixl, jyf, kzlfd) + jac(ixl, jyf + 1, kzlfd) &
                &+ jac(ixl, jyf, kzlfd + 1) + jac(ixl, jyf + 1, kzlfd + 1))
            flwlfu = 0.0
          else
            flwlfd = 0.0
            flwlfu = 0.0
          end if
        end if

        if(zrbu < topography_surface(ixr, jyb)) then
          flwrbd = 0.0
          flwrbu = 0.0
        else if(zrbd < topography_surface(ixr, jyb)) then
          flwrbd = 0.0
          flwrbu = (var%v(ixr, jyb, kzrbu + 1) - var%v(ixr, jyb, kzrbu)) / dz &
              &* 4.0 / (jac(ixr, jyb, kzrbu) + jac(ixr, jyb + 1, kzrbu) &
              &+ jac(ixr, jyb, kzrbu + 1) + jac(ixr, jyb + 1, kzrbu + 1))
        else
          if(zrbu < lz(1)) then
            flwrbd = (var%v(ixr, jyb, kzrbd + 1) - var%v(ixr, jyb, kzrbd)) &
                &/ dz * 4.0 / (jac(ixr, jyb, kzrbd) + jac(ixr, jyb + 1, kzrbd) &
                &+ jac(ixr, jyb, kzrbd + 1) + jac(ixr, jyb + 1, kzrbd + 1))
            flwrbu = (var%v(ixr, jyb, kzrbu + 1) - var%v(ixr, jyb, kzrbu)) &
                &/ dz * 4.0 / (jac(ixr, jyb, kzrbu) + jac(ixr, jyb + 1, kzrbu) &
                &+ jac(ixr, jyb, kzrbu + 1) + jac(ixr, jyb + 1, kzrbu + 1))
          else if(zrbd < lz(1)) then
            flwrbd = (var%v(ixr, jyb, kzrbd + 1) - var%v(ixr, jyb, kzrbd)) &
                &/ dz * 4.0 / (jac(ixr, jyb, kzrbd) + jac(ixr, jyb + 1, kzrbd) &
                &+ jac(ixr, jyb, kzrbd + 1) + jac(ixr, jyb + 1, kzrbd + 1))
            flwrbu = 0.0
          else
            flwrbd = 0.0
            flwrbu = 0.0
          end if
        end if

        if(zrfu < topography_surface(ixr, jyf)) then
          flwrfd = 0.0
          flwrfu = 0.0
        else if(zrfd < topography_surface(ixr, jyf)) then
          flwrfd = 0.0
          flwrfu = (var%v(ixr, jyf, kzrfu + 1) - var%v(ixr, jyf, kzrfu)) / dz &
              &* 4.0 / (jac(ixr, jyf, kzrfu) + jac(ixr, jyf + 1, kzrfu) &
              &+ jac(ixr, jyf, kzrfu + 1) + jac(ixr, jyf + 1, kzrfu + 1))
        else
          if(zrfu < lz(1)) then
            flwrfd = (var%v(ixr, jyf, kzrfd + 1) - var%v(ixr, jyf, kzrfd)) &
                &/ dz * 4.0 / (jac(ixr, jyf, kzrfd) + jac(ixr, jyf + 1, kzrfd) &
                &+ jac(ixr, jyf, kzrfd + 1) + jac(ixr, jyf + 1, kzrfd + 1))
            flwrfu = (var%v(ixr, jyf, kzrfu + 1) - var%v(ixr, jyf, kzrfu)) &
                &/ dz * 4.0 / (jac(ixr, jyf, kzrfu) + jac(ixr, jyf + 1, kzrfu) &
                &+ jac(ixr, jyf, kzrfu + 1) + jac(ixr, jyf + 1, kzrfu + 1))
          else if(zrfd < lz(1)) then
            flwrfd = (var%v(ixr, jyf, kzrfd + 1) - var%v(ixr, jyf, kzrfd)) &
                &/ dz * 4.0 / (jac(ixr, jyf, kzrfd) + jac(ixr, jyf + 1, kzrfd) &
                &+ jac(ixr, jyf, kzrfd + 1) + jac(ixr, jyf + 1, kzrfd + 1))
            flwrfu = 0.0
          else
            flwrfd = 0.0
            flwrfu = 0.0
          end if
        end if
      else
        ! for variable at intermediate levels
        ! z(k+1/2) = lz(0) + k*dz:

        ! index for lower level used for linear interpolation
        kzd = max(- nbz, floor((zlc - lz(0)) / dz))

        ! index for upper level used for linear interpolation
        kzu = kzd + 1

        ! if the ray volume is above the model domain make sure that the
        ! interpolation below leads to a zero gradient
        if(kzd > nz) then
          kzu = nz + 1
          kzd = nz + 1
        end if

        ! intermediate levels z(k+1/2) = z(k) + 0.5*dz
        zu = z(kzu) + 0.5 * dz
        zd = z(kzd) + 0.5 * dz

        if(sizeY == 1) then
          jyb = 1
          jyf = 1
        else
          ! for variable at intermediate levels
          ! y(j+1/2) = y(j) + 0.5*dy = ly(0) +  j*dy:

          ! index for backward level used for linear interpolation
          jyb = floor((ylc - ly(0)) / dy) - jy0

          if(jyb < - nby) then
            print *, 'ERROR IN MEANFLOW: jyb =', jyb, '< -nby =', - nby
            stop
          end if

          ! index for forward level used for linear interpolation
          jyf = jyb + 1

          if(jyf > ny + nby) then
            print *, 'ERROR IN MEANFLOW: jyf =', jyf, '> ny + nby =', ny + nby
            stop
          end if
        end if

        ! intermediate levels y(j+1/2) = y(j) + 0.5*dy
        yf = y(jyf + jy0) + 0.5 * dy
        yb = y(jyb + jy0) + 0.5 * dy

        if(sizeX == 1) then
          ixl = 1
          ixr = 1
        else
          ! for variable at full levels x(i) = lx(0) + 0.5*dx + (i-1)*dx

          ! index for leftmost level used for linear interpolation
          ixl = floor((xlc - 0.5 * dx - lx(0)) / dx) + 1 - ix0

          if(ixl < - nbx) then
            print *, 'ERROR IN MEANFLOW: ixl =', ixl, '< -nbx =', - nbx
            ! testb
            print *, 'variable = dv/dz'
            print *, 'xlc =', xlc
            print *, 'lx(0) =', lx(0)
            print *, 'lx(1) =', lx(1)
            print *, 'dx =', dx
            print *, 'lRef =', lRef
            ! teste
            stop
          end if

          ! index for rightmost level used for linear interpolation
          ixr = ixl + 1

          if(ixr > nx + nbx) then
            print *, 'ERROR IN MEANFLOW: ixr =', ixr, '> nx + nbx =', nx + nbx
            stop
          end if
        end if

        ! full levels x(i)
        xr = x(ixr + ix0)
        xl = x(ixl + ix0)

        ! values of var. at the eight corners of the interpolation region
        ! dv/dz (i,j+1/2,k+1/2)) = (v(i,j+1/2,k+1) - v(i,j+1/2,k))/dz, hence
        ! (using 0 at levels below the model bottom)
        if(zu < lz(0)) then
          flwlbd = 0.0
          flwlbu = 0.0

          flwlfd = 0.0
          flwlfu = 0.0

          flwrbd = 0.0
          flwrbu = 0.0

          flwrfd = 0.0
          flwrfu = 0.0
        elseif(zd < lz(0)) then
          flwlbd = 0.0
          flwlbu = (var%v(ixl, jyb, kzu + 1) - var%v(ixl, jyb, kzu)) / dz

          flwlfd = 0.0
          flwlfu = (var%v(ixl, jyf, kzu + 1) - var%v(ixl, jyf, kzu)) / dz

          flwrbd = 0.0
          flwrbu = (var%v(ixr, jyb, kzu + 1) - var%v(ixr, jyb, kzu)) / dz

          flwrfd = 0.0
          flwrfu = (var%v(ixr, jyf, kzu + 1) - var%v(ixr, jyf, kzu)) / dz
        else
          if(zu < lz(1)) then
            flwlbd = (var%v(ixl, jyb, kzd + 1) - var%v(ixl, jyb, kzd)) / dz
            flwlbu = (var%v(ixl, jyb, kzu + 1) - var%v(ixl, jyb, kzu)) / dz

            flwlfd = (var%v(ixl, jyf, kzd + 1) - var%v(ixl, jyf, kzd)) / dz
            flwlfu = (var%v(ixl, jyf, kzu + 1) - var%v(ixl, jyf, kzu)) / dz

            flwrbd = (var%v(ixr, jyb, kzd + 1) - var%v(ixr, jyb, kzd)) / dz
            flwrbu = (var%v(ixr, jyb, kzu + 1) - var%v(ixr, jyb, kzu)) / dz

            flwrfd = (var%v(ixr, jyf, kzd + 1) - var%v(ixr, jyf, kzd)) / dz
            flwrfu = (var%v(ixr, jyf, kzu + 1) - var%v(ixr, jyf, kzu)) / dz
          elseif(zd < lz(1)) then
            flwlbd = (var%v(ixl, jyb, kzd + 1) - var%v(ixl, jyb, kzd)) / dz
            flwlbu = 0.0

            flwlfd = (var%v(ixl, jyf, kzd + 1) - var%v(ixl, jyf, kzd)) / dz
            flwlfu = 0.0

            flwrbd = (var%v(ixr, jyb, kzd + 1) - var%v(ixr, jyb, kzd)) / dz
            flwrbu = 0.0

            flwrfd = (var%v(ixr, jyf, kzd + 1) - var%v(ixr, jyf, kzd)) / dz
            flwrfu = 0.0
          else
            flwlbd = 0.0
            flwlbu = 0.0

            flwlfd = 0.0
            flwlfu = 0.0

            flwrbd = 0.0
            flwrbu = 0.0

            flwrfd = 0.0
            flwrfu = 0.0
          end if
        end if
      end if
    else
      print *, 'ERROR: UNKNOWN flwtpe =', flwtpe, 'IN MEANFLOW'
      stop
    end if

    ! interpolation in x

    if(sizeX == 1) then
      flwbd = flwlbd
      flwbu = flwlbu

      flwfd = flwlfd
      flwfu = flwlfu

      ! FJMay2023
      if(topography) then
        zbd = zlbd
        zbu = zlbu

        zfd = zlfd
        zfu = zlfu
      end if
    else
      if(xr < xl) then
        print *, 'ERROR IN MEANFLOW: xr =', xr, '< xl =', xl
        stop
      elseif(xr == xl) then
        factor = 0.0
      elseif(xlc > xr) then
        factor = 0.0
      elseif(xlc > xl) then
        factor = (xr - xlc) / dx
      else
        factor = 1.0
      end if

      flwbd = factor * flwlbd + (1.0 - factor) * flwrbd
      flwbu = factor * flwlbu + (1.0 - factor) * flwrbu

      flwfd = factor * flwlfd + (1.0 - factor) * flwrfd
      flwfu = factor * flwlfu + (1.0 - factor) * flwrfu

      ! FJMay2023
      if(topography) then
        zbd = factor * zlbd + (1.0 - factor) * zrbd
        zbu = factor * zlbu + (1.0 - factor) * zrbu

        zfd = factor * zlfd + (1.0 - factor) * zrfd
        zfu = factor * zlfu + (1.0 - factor) * zrfu
      end if
    end if

    ! interpolation in y

    if(sizeY == 1) then
      flwd = flwbd
      flwu = flwbu

      ! FJMay2023
      if(topography) then
        zd = zbd
        zu = zbu
      end if
    else
      if(yf < yb) then
        print *, 'ERROR IN MEANFLOW: yf =', yf, '< yb =', yb
        stop
      elseif(yf == yb) then
        factor = 0.0
      elseif(ylc > yf) then
        factor = 0.0
      elseif(ylc > yb) then
        factor = (yf - ylc) / dy
      else
        factor = 1.0
      end if

      flwd = factor * flwbd + (1.0 - factor) * flwfd
      flwu = factor * flwbu + (1.0 - factor) * flwfu

      ! FJMay2023
      if(topography) then
        zd = factor * zbd + (1.0 - factor) * zfd
        zu = factor * zbu + (1.0 - factor) * zfu
      end if
    end if

    ! interpolation in z

    if(zu < zd) then
      print *, 'ERROR IN MEANFLOW: zu =', zu, '< zd =', zd
      stop
    elseif(zu == zd) then
      factor = 0.0
    elseif(zlc > zu) then
      factor = 0.0
    elseif(zlc > zd) then
      ! FJApr2023
      ! factor = (zu - zlc) / dz
      factor = (zu - zlc) / (zu - zd)
    else
      factor = 1.0
    end if

    flw = factor * flwd + (1.0 - factor) * flwu

    return

  end subroutine meanflow

  !----------------------------------------------------------------------

  subroutine split_rayvol(ray)

    !-----------------------------------------------------------------
    ! splits ray volumes with a spatial extension (in any of the three
    ! directions) larger than the corresponding cell extension
    !-----------------------------------------------------------------

    implicit none

    ! argument list
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(inout) :: ray

    integer :: nrlc
    integer :: ix, jy, kz

    real :: xr, yr, zr
    real :: dxr, dyr, dzr
    real :: axk, ayl, azm
    real :: dphi, deltaPhi
    integer :: nrvtt0, nrvtt1, nrvloc

    if(steady_state) return

    ! total number of ray volumes before splitting

    nrvloc = sum(nRay(1:nx, 1:ny, 1:nz))

    ! testb
    ! print*,'before splitting nrvloc =',nrvloc
    ! teste

    call mpi_reduce(nrvloc, nrvtt0, 1, mpi_integer, mpi_sum, root, comm, ierror)
    call mpi_bcast(nrvtt0, 1, mpi_integer, root, comm, ierror)

    do kz = 1, nz
      do jy = 1, ny
        do ix = 1, nx
          if(nRay(ix, jy, kz) < 1) cycle

          nrlc = nRay(ix, jy, kz)

          !splitting in x direction

          if(sizeX > 1) then
            do iRay = 1, nRay(ix, jy, kz)
              xr = ray(iRay, ix, jy, kz)%x
              dxr = ray(iRay, ix, jy, kz)%dxray

              if(dxr > dx) then
                nrlc = nrlc + 1

                if(nrlc > nray_wrk) then
                  print *, 'r.v. getting too many by splitting in x  direction'
                  stop
                end if

                xr = ray(iRay, ix, jy, kz)%x

                axk = ray(iRay, ix, jy, kz)%area_xk

                ray(iRay, ix, jy, kz)%dxray = 0.5 * dxr

                ray(iRay, ix, jy, kz)%area_xk = 0.5 * axk

                ray(nrlc, ix, jy, kz) = ray(iRay, ix, jy, kz)

                ray(iRay, ix, jy, kz)%x = xr - 0.25 * dxr
                ray(nrlc, ix, jy, kz)%x = xr + 0.25 * dxr

                !SDJul2024
                if(update_phase) then
                  dphi = ray(iRay, ix, jy, kz)%dphi
                  deltaPhi = 0.25 * dxr * ray(iRay, ix, jy, kz)%k
                  ray(iRay, ix, jy, kz)%dphi = dphi - deltaPhi
                  ray(nrlc, ix, jy, kz)%dphi = dphi + deltaPhi
                end if
              end if
            end do

            if(nrlc > nRay(ix, jy, kz)) then
              nRay(ix, jy, kz) = nrlc

              if(nRay(ix, jy, kz) > nray_wrk) then
                print *, 'ERROR at ix,jy,kz =', ix, jy, kz
                print *, 'nRay =', nRay(ix, jy, kz), '> nray_wrk =', nray_wrk
                stop
              end if
            end if
          end if

          !splitting in y direction

          if(sizeY > 1) then
            do iRay = 1, nRay(ix, jy, kz)
              yr = ray(iRay, ix, jy, kz)%y
              dyr = ray(iRay, ix, jy, kz)%dyray

              if(dyr > dy) then
                nrlc = nrlc + 1

                if(nrlc > nray_wrk) then
                  print *, 'r.v. getting too many by splitting in y  direction'
                  stop
                end if

                yr = ray(iRay, ix, jy, kz)%y

                ayl = ray(iRay, ix, jy, kz)%area_yl

                ray(iRay, ix, jy, kz)%dyray = 0.5 * dyr

                ray(iRay, ix, jy, kz)%area_yl = 0.5 * ayl

                ray(nrlc, ix, jy, kz) = ray(iRay, ix, jy, kz)

                ray(iRay, ix, jy, kz)%y = yr - 0.25 * dyr
                ray(nrlc, ix, jy, kz)%y = yr + 0.25 * dyr

                !SDJul2024
                if(update_phase) then
                  dphi = ray(iRay, ix, jy, kz)%dphi
                  deltaPhi = 0.25 * dyr * ray(iRay, ix, jy, kz)%l
                  ray(iRay, ix, jy, kz)%dphi = dphi - deltaPhi
                  ray(nrlc, ix, jy, kz)%dphi = dphi + deltaPhi
                end if
              end if
            end do

            if(nrlc > nRay(ix, jy, kz)) then
              nRay(ix, jy, kz) = nrlc

              if(nRay(ix, jy, kz) > nray_wrk) then
                print *, 'ERROR at ix,jy,kz =', ix, jy, kz
                print *, 'nRay =', nRay(ix, jy, kz), '> nray_wrk =', nray_wrk
                stop
              end if
            end if
          end if

          !splitting in x direction

          do iRay = 1, nRay(ix, jy, kz)
            zr = ray(iRay, ix, jy, kz)%z
            dzr = ray(iRay, ix, jy, kz)%dzray

            if(dzr > dz) then
              nrlc = nrlc + 1

              if(nrlc > nray_wrk) then
                print *, 'r.v. getting too many by splitting in z  direction'
                stop
              end if

              zr = ray(iRay, ix, jy, kz)%z

              azm = ray(iRay, ix, jy, kz)%area_zm

              ray(iRay, ix, jy, kz)%dzray = 0.5 * dzr

              ray(iRay, ix, jy, kz)%area_zm = 0.5 * azm

              ray(nrlc, ix, jy, kz) = ray(iRay, ix, jy, kz)

              ray(iRay, ix, jy, kz)%z = zr - 0.25 * dzr
              ray(nrlc, ix, jy, kz)%z = zr + 0.25 * dzr

              if(update_phase) then
                dphi = ray(iRay, ix, jy, kz)%dphi
                deltaPhi = 0.25 * dzr * ray(iRay, ix, jy, kz)%m
                ray(iRay, ix, jy, kz)%dphi = dphi - deltaPhi
                ray(nrlc, ix, jy, kz)%dphi = dphi + deltaPhi
              end if
            end if
          end do

          if(nrlc > nRay(ix, jy, kz)) then
            nRay(ix, jy, kz) = nrlc

            if(nRay(ix, jy, kz) > nray_wrk) then
              print *, 'ERROR at ix,jy,kz =', ix, jy, kz
              print *, 'nRay =', nRay(ix, jy, kz), '> nray_wrk =', nray_wrk
              stop
            end if
          end if
        end do
      end do
    end do

    ! total number of ray volumes after splitting

    nrvloc = sum(nRay(1:nx, 1:ny, 1:nz))

    ! testb
    ! print*,'after splitting nrvloc =',nrvloc
    ! teste

    call mpi_reduce(nrvloc, nrvtt1, 1, mpi_integer, mpi_sum, root, comm, ierror)
    call mpi_bcast(nrvtt1, 1, mpi_integer, root, comm, ierror)

    if(master .and. nrvtt1 > nrvtt0) then
      print *, 'after splitting nray =', nrvtt1
    end if

    return

  end subroutine split_rayvol

  !----------------------------------------------------------------------

  subroutine shift_rayvol(ray)

    !-----------------------------------------------------------------
    ! shifts ray volumes to appropriate cell
    !-----------------------------------------------------------------

    implicit none

    ! argument list
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(inout) :: ray

    integer :: ix, jy, kz
    integer :: nrlc
    integer :: jRay
    integer :: nsl, nsr, nsb, nsf, nsd, nsu
    integer :: ish, jsh, ksh, nsh
    integer, dimension(:, :, :), allocatable :: nshl, nshr, nshb, nshf, nshd, &
        &nshu
    integer, dimension(:, :, :, :), allocatable :: irsl, irsr, irsb, irsf, &
        &irsd, irsu

    integer :: irsh(nray_wrk)

    logical :: lplace

    real :: xr, yr, zr

    integer :: ix0, jy0

    logical :: outside

    if(steady_state) return

    allocate(nshl(0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(nshr(0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(nshb(0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(nshf(0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(nshd(0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(nshu(0:nx + 1, 0:ny + 1, 0:nz + 1))

    allocate(irsl(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(irsr(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(irsb(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(irsf(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(irsd(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1))
    allocate(irsu(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1))

    ix0 = is + nbx - 1
    jy0 = js + nby - 1

    !------------------------
    ! shifting in x direction
    !------------------------

    if(sizeX > 1) then
      nshl = 0
      nshr = 0

      irsl = 0
      irsr = 0

      ! periodic boundary conditions in x

      call setboundary_rayvol_x(ray)

      ! move ray volumes to appropriate cell:
      ! check whether r.v. from neighboring cells have propagated into
      ! the cell and re-associate them

      do kz = 0, nz + 1
        do jy = 0, ny + 1
          do ix = 1, nx
            ! # of r.v. in the cell
            nrlc = nRay(ix, jy, kz)

            ! # r.v. taken from the cell to the left
            nsl = 0
            ! # r.v. taken from the cell to the right
            nsr = 0

            !  take over ray volumes from cell to the left

            if(nRay(ix - 1, jy, kz) > 0) then
              do iRay = 1, nRay(ix - 1, jy, kz)
                xr = ray(iRay, ix - 1, jy, kz)%x

                if(xr > x(ix - 1 + ix0) + 0.5 * dx) then
                  nrlc = nrlc + 1

                  nsl = nsl + 1

                  ray(nrlc, ix, jy, kz) = ray(iRay, ix - 1, jy, kz)

                  irsl(nsl, ix, jy, kz) = iRay
                end if
              end do

              ! # r.v. taken from the cell to the left
              nshl(ix, jy, kz) = nsl
            end if

            ! ray volumes from cell to the right

            if(nRay(ix + 1, jy, kz) > 0) then
              do iRay = 1, nRay(ix + 1, jy, kz)
                xr = ray(iRay, ix + 1, jy, kz)%x

                if(xr < x(ix + 1 + ix0) - 0.5 * dx) then
                  nrlc = nrlc + 1

                  nsr = nsr + 1

                  ray(nrlc, ix, jy, kz) = ray(iRay, ix + 1, jy, kz)

                  irsr(nsr, ix, jy, kz) = iRay
                end if
              end do

              ! # r.v. taken from the cell to the right
              nshr(ix, jy, kz) = nsr
            end if

            ! new # of r.v. in the cell
            if(nrlc > nray_wrk) then
              print *, 'at ix,jy,kz =', ix, jy, kz, 'nrlc > nray_wrk'
              stop
            else
              nRay(ix, jy, kz) = nrlc
            end if
          end do
        end do
      end do

      ! periodic boundary conditions in x for the index and number arrays
      ! filled above

      call setboundary_nshift_x(nshl, nshr)
      call setboundary_irshift_x(irsl, irsr)

      ! remove ray volumes from cells they have left

      do kz = 0, nz + 1
        do jy = 0, ny + 1
          do ix = 1, nx
            ! ordered list of ray indices for r.v. that have left a
            ! cell, stored in array irsh

            ! # of r.v. that have been transferred away
            nsh = 0

            ! r.v. that have been shifted to the right

            if(nshl(ix + 1, jy, kz) > 0) then
              do ish = 1, nshl(ix + 1, jy, kz)
                iRay = irsl(ish, ix + 1, jy, kz)

                ! list of r.v. shifted to the right already ordered
                ! by r.v. index
                irsh(ish) = iRay
              end do

              nsh = nshl(ix + 1, jy, kz)
            end if

            ! r.v. that have been shifted to the left

            if(nshr(ix - 1, jy, kz) > 0) then
              do ish = 1, nshr(ix - 1, jy, kz)
                iRay = irsr(ish, ix - 1, jy, kz)

                ! list of r.v. shifted to the right already ordered
                ! by r.v. index, but they have to be taken up into
                ! the complete list of all shifted r.v. so that
                ! they are  all ordered correctly by cell index
                if(nsh == 0) then
                  irsh(1) = iRay
                elseif(iRay < irsh(1)) then
                  do jsh = nsh, 1, - 1
                    irsh(jsh + 1) = irsh(jsh)
                  end do

                  irsh(1) = iRay
                elseif(iRay > irsh(nsh)) then
                  irsh(nsh + 1) = iRay
                else
                  !testb
                  if(nsh < 2) then
                    print *, 'ERROR: nsh < 2 in x-shifting where  this must &
                        &not be the case'
                    print *, 'ix,jy,kz =', ix, jy, kz
                    print *, 'nshl(ix+1,jy,kz) =', nshl(ix + 1, jy, kz)
                    print *, 'nshr(ix-1,jy,kz) =', nshr(ix - 1, jy, kz)
                    print *, 'ish =', ish
                    print *, 'iRay =', iRay
                    print *, 'nsh =', nsh
                    print *, 'irsh(1) =', irsh(1)
                    print *, 'irsh(nsh) =', irsh(nsh)
                    print *, 'xr =', ray(iRay, ix, jy, kz)%x
                    print *, 'x(ix0+ix) - 0.5*dx', x(ix0 + ix) - 0.5 * dx
                    print *, 'x(ix0+ix) + 0.5*dx', x(ix0 + ix) + 0.5 * dx
                    stop
                  end if
                  !teste

                  lplace = .false.

                  do jsh = 1, nsh - 1
                    if(lplace) cycle

                    if(iRay < irsh(jsh + 1)) then
                      do ksh = nsh, jsh + 1, - 1
                        irsh(ksh + 1) = irsh(ksh)
                      end do

                      irsh(jsh + 1) = iRay

                      lplace = .true.
                    end if
                  end do
                end if

                nsh = nsh + 1
              end do
            end if

            if(nsh /= nshl(ix + 1, jy, kz) + nshr(ix - 1, jy, kz)) then
              print *, 'at ix,jy,kz =', ix, jy, kz, 'nsh /= nshl + nshr'
              stop
            end if

            ! remove shifted ray volumes

            if(nsh > 0) then
              ! # of r.v. in the cell
              nrlc = nRay(ix, jy, kz)

              do ish = nsh, 1, - 1
                iRay = irsh(ish)

                if(iRay < nrlc) then
                  do jRay = iRay + 1, nrlc
                    !testb
                    if(jRay < 0) then
                      print *, 'ish =', ish
                      print *, 'nsh =', nsh
                      print *, 'iRay =', iRay
                      print *, 'nrlc =', nrlc
                      print *, 'jRay =', jRay
                      print *, 'ix,jy,kz =', ix, jy, kz
                      stop
                    end if
                    !teste
                    ray(jRay - 1, ix, jy, kz) = ray(jRay, ix, jy, kz)
                  end do
                end if

                nrlc = nrlc - 1
              end do

              ! new # of r.v. in the cell
              if(nrlc /= nRay(ix, jy, kz) - nshl(ix + 1, jy, kz) - nshr(ix &
                  &- 1, jy, kz)) then
                print *, 'at ix,jy,kz =', ix, jy, kz, 'nrlc &
                    &/= nRay(ix,jy,kz)   - nshl(ix+1,jy,kz) - nshr(ix-1,jy,kz)'
                stop
              else
                nRay(ix, jy, kz) = nrlc
              end if
            end if
          end do ! ix
        end do ! jy
      end do ! kz

      ! periodic boundary conditions in x

      call setboundary_rayvol_x(ray)
    end if

    !------------------------
    ! shifting in y direction
    !------------------------

    if(sizeY > 1) then
      ! procedure completely analogous to shifting in x ...

      nshb = 0
      nshf = 0

      irsb = 0
      irsf = 0

      ! periodic boundary conditions in y

      call setboundary_rayvol_y(ray)

      ! move ray volumes to appropriate cell

      do kz = 0, nz + 1
        do jy = 1, ny
          do ix = 0, nx + 1
            nrlc = nRay(ix, jy, kz)

            nsb = 0
            nsf = 0

            ! ray volumes from cell behind

            if(nRay(ix, jy - 1, kz) > 0) then
              do iRay = 1, nRay(ix, jy - 1, kz)
                yr = ray(iRay, ix, jy - 1, kz)%y

                if(yr > y(jy - 1 + jy0) + 0.5 * dy) then
                  nrlc = nrlc + 1

                  nsb = nsb + 1

                  ray(nrlc, ix, jy, kz) = ray(iRay, ix, jy - 1, kz)

                  irsb(nsb, ix, jy, kz) = iRay
                end if
              end do

              nshb(ix, jy, kz) = nsb
            end if

            ! ray volumes from cell in front

            if(nRay(ix, jy + 1, kz) > 0) then
              do iRay = 1, nRay(ix, jy + 1, kz)
                yr = ray(iRay, ix, jy + 1, kz)%y

                if(yr < y(jy + 1 + jy0) - 0.5 * dy) then
                  nrlc = nrlc + 1

                  nsf = nsf + 1

                  ray(nrlc, ix, jy, kz) = ray(iRay, ix, jy + 1, kz)

                  irsf(nsf, ix, jy, kz) = iRay
                end if
              end do

              nshf(ix, jy, kz) = nsf
            end if

            if(nrlc > nray_wrk) then
              print *, 'at ix,jy,kz =', ix, jy, kz, 'nrlc > nray_wrk'
              stop
            else
              nRay(ix, jy, kz) = nrlc
            end if
          end do
        end do
      end do

      ! periodic boundary conditions in y

      call setboundary_nshift_y(nshb, nshf)
      call setboundary_irshift_y(irsb, irsf)

      ! remove ray volumes from cells they have left

      do kz = 0, nz + 1
        do jy = 1, ny
          do ix = 0, nx + 1
            ! ordered list of ray indices for r.v. that have left the
            ! cell

            nsh = 0

            ! r.v. that have been shifted forward

            if(nshb(ix, jy + 1, kz) > 0) then
              do ish = 1, nshb(ix, jy + 1, kz)
                iRay = irsb(ish, ix, jy + 1, kz)
                irsh(ish) = iRay
              end do

              nsh = nshb(ix, jy + 1, kz)
            end if

            ! r.v. that have been shifted backward

            if(nshf(ix, jy - 1, kz) > 0) then
              do ish = 1, nshf(ix, jy - 1, kz)
                iRay = irsf(ish, ix, jy - 1, kz)

                if(nsh == 0) then
                  irsh(1) = iRay
                elseif(iRay < irsh(1)) then
                  do jsh = nsh, 1, - 1
                    irsh(jsh + 1) = irsh(jsh)
                  end do

                  irsh(1) = iRay
                elseif(iRay > irsh(nsh)) then
                  irsh(nsh + 1) = iRay
                else
                  !testb
                  if(nsh < 2) then
                    print *, 'ERROR: nsh < 2 in y-shifting where  this must &
                        &not be the case'
                    stop
                  end if
                  !teste

                  lplace = .false.

                  do jsh = 1, nsh - 1
                    if(lplace) cycle

                    if(iRay < irsh(jsh + 1)) then
                      do ksh = nsh, jsh + 1, - 1
                        irsh(ksh + 1) = irsh(ksh)
                      end do

                      irsh(jsh + 1) = iRay

                      lplace = .true.
                    end if
                  end do
                end if

                nsh = nsh + 1
              end do
            end if

            if(nsh /= nshb(ix, jy + 1, kz) + nshf(ix, jy - 1, kz)) then
              print *, 'at ix,jy,kz =', ix, jy, kz, 'nsh /= nshb + nshf'
              stop
            end if

            ! remove shifted ray volumes

            if(nsh > 0) then
              nrlc = nRay(ix, jy, kz)

              do ish = nsh, 1, - 1
                iRay = irsh(ish)

                if(iRay < nrlc) then
                  do jRay = iRay + 1, nrlc
                    ray(jRay - 1, ix, jy, kz) = ray(jRay, ix, jy, kz)
                  end do
                end if

                nrlc = nrlc - 1
              end do

              if(nrlc /= nRay(ix, jy, kz) - nshb(ix, jy + 1, kz) - nshf(ix, jy &
                  &- 1, kz)) then
                print *, 'at ix,jy,kz =', ix, jy, kz, 'nrlc &
                    &/= nRay(ix,jy,kz)   - nshb(ix,jy+1,kz) - nshf(ix,jy-1,kz)'
                stop
              else
                nRay(ix, jy, kz) = nrlc
              end if
            end if
          end do ! ix
        end do ! jy
      end do ! kz

      ! periodic boundary conditions in y

      call setboundary_rayvol_y(ray)
    end if

    !------------------------
    ! shifting in z direction
    !------------------------

    if(sizeZ > 1) then
      nshd = 0
      nshu = 0

      irsd = 0
      irsu = 0

      ! boundary conditions in z

      call setboundary_rayvol_z(ray)

      ! move ray volumes to appropriate cell

      do kz = 1, nz
        do jy = 0, ny + 1
          do ix = 0, nx + 1
            nrlc = nRay(ix, jy, kz)

            nsd = 0
            nsu = 0

            ! ray volumes from cell below

            if(nRay(ix, jy, kz - 1) > 0) then
              do iRay = 1, nRay(ix, jy, kz - 1)
                zr = ray(iRay, ix, jy, kz - 1)%z

                if(.not. topography .and. zr > z(kz - 1) + 0.5 * dz) then
                  outside = .true.
                else if(topography) then
                  if(zr > zTFC(ix, jy, kz - 1) + 0.5 * jac(ix, jy, kz - 1) &
                      &* dz) then
                    outside = .true.
                  else
                    outside = .false.
                  end if
                else
                  outside = .false.
                end if
                if(outside) then
                  nrlc = nrlc + 1

                  nsd = nsd + 1

                  ray(nrlc, ix, jy, kz) = ray(iRay, ix, jy, kz - 1)

                  irsd(nsd, ix, jy, kz) = iRay
                end if
              end do

              nshd(ix, jy, kz) = nsd
            end if

            ! ray volumes from cell above

            if(nRay(ix, jy, kz + 1) > 0) then
              do iRay = 1, nRay(ix, jy, kz + 1)
                zr = ray(iRay, ix, jy, kz + 1)%z

                if(.not. topography .and. zr < z(kz + 1) - 0.5 * dz) then
                  outside = .true.
                else if(topography) then
                  if(zr < zTFC(ix, jy, kz + 1) - 0.5 * jac(ix, jy, kz + 1) &
                      &* dz) then
                    outside = .true.
                  else
                    outside = .false.
                  end if
                else
                  outside = .false.
                end if
                if(outside) then
                  nrlc = nrlc + 1

                  nsu = nsu + 1

                  ray(nrlc, ix, jy, kz) = ray(iRay, ix, jy, kz + 1)

                  irsu(nsu, ix, jy, kz) = iRay
                end if
              end do

              nshu(ix, jy, kz) = nsu
            end if

            if(nrlc > nray_wrk) then
              print *, 'at ix,jy,kz =', ix, jy, kz, 'nrlc > nray_wrk'
              stop
            else
              nRay(ix, jy, kz) = nrlc
            end if
          end do
        end do
      end do

      ! periodic boundary conditions in z for the index and number arrays
      ! filled above

      if(zBoundary == "periodic") then
        nshd(:, :, 0) = nshd(:, :, nz)
        nshd(:, :, nz + 1) = nshd(:, :, 1)
        nshu(:, :, 0) = nshu(:, :, nz)
        nshu(:, :, nz + 1) = nshu(:, :, 1)
      end if

      ! remove ray volumes from cells they have left

      do kz = 1, nz
        do jy = 0, ny + 1
          do ix = 0, nx + 1
            ! ordered list of ray indices for r.v. that have left the
            ! cell

            nsh = 0

            ! r.v. that have been shifted upward

            if(nshd(ix, jy, kz + 1) > 0) then
              do ish = 1, nshd(ix, jy, kz + 1)
                iRay = irsd(ish, ix, jy, kz + 1)
                irsh(ish) = iRay
              end do

              nsh = nshd(ix, jy, kz + 1)
            end if

            ! r.v. that have been shifted downward

            if(nshu(ix, jy, kz - 1) > 0) then
              do ish = 1, nshu(ix, jy, kz - 1)
                iRay = irsu(ish, ix, jy, kz - 1)

                if(nsh == 0) then
                  irsh(1) = iRay
                elseif(iRay < irsh(1)) then
                  do jsh = nsh, 1, - 1
                    irsh(jsh + 1) = irsh(jsh)
                  end do

                  irsh(1) = iRay
                elseif(iRay > irsh(nsh)) then
                  irsh(nsh + 1) = iRay
                else
                  !testb
                  if(nsh < 2) then
                    print *, 'ERROR: nsh < 2 in z-shifting where  this must &
                        &not be the case'
                    stop
                  end if
                  !teste

                  lplace = .false.

                  do jsh = 1, nsh - 1
                    if(lplace) cycle

                    if(iRay < irsh(jsh + 1)) then
                      do ksh = nsh, jsh + 1, - 1
                        irsh(ksh + 1) = irsh(ksh)
                      end do

                      irsh(jsh + 1) = iRay

                      lplace = .true.
                    end if
                  end do
                end if

                nsh = nsh + 1
              end do
            end if

            if(nsh /= nshd(ix, jy, kz + 1) + nshu(ix, jy, kz - 1)) then
              print *, 'at ix,jy,kz =', ix, jy, kz, 'nsh /= nshd + nshu'
              stop
            end if

            ! remove shifted ray volumes

            if(nsh > 0) then
              nrlc = nRay(ix, jy, kz)

              do ish = nsh, 1, - 1
                iRay = irsh(ish)

                if(iRay < nrlc) then
                  do jRay = iRay + 1, nrlc
                    ray(jRay - 1, ix, jy, kz) = ray(jRay, ix, jy, kz)
                  end do
                end if

                nrlc = nrlc - 1
              end do

              if(nrlc /= nRay(ix, jy, kz) - nshd(ix, jy, kz + 1) - nshu(ix, &
                  &jy, kz - 1)) then
                print *, 'at ix,jy,kz =', ix, jy, kz, 'nrlc &
                    &/= nRay(ix,jy,kz)   - nshd(ix,jy,kz+1) - nshu(ix,jy,kz-1)'
                stop
              else
                nRay(ix, jy, kz) = nrlc
              end if
            end if
          end do ! ix
        end do ! jy
      end do ! kz

      ! boundary conditions in z

      call setboundary_rayvol_z(ray)
    end if

    ! testb
    do kz = 0, nz + 1
      do jy = 0, ny + 1
        do ix = 0, nx + 1
          if(nRay(ix, jy, kz) > 0) then
            do iRay = 1, nRay(ix, jy, kz)
              if(sizeX > 1) then
                xr = ray(iRay, ix, jy, kz)%x

                if(xr < x(ix + ix0) - 0.5 * dx) then
                  print *, 'ERROR in shift_rayvol:'
                  print *, 'xr =', xr, '< x(ix+ix0) - 0.5*dx =', x(ix + ix0) &
                      &- 0.5 * dx
                  print *, 'ix =', ix
                  print *, 'ix0 =', ix0
                  print *, 'x(ix+ix0) =', x(ix + ix0)
                  print *, 'dx =', dx
                  print *, 'iRay,ix,jy,kz =', iRay, ix, jy, kz
                  stop
                end if

                if(xr > x(ix + ix0) + 0.5 * dx) then
                  print *, 'ERROR in shift_rayvol:'
                  print *, 'xr =', xr, '> x(ix+ix0) + 0.5*dx =', x(ix + ix0) &
                      &- 0.5 * dx
                  print *, 'ix =', ix
                  print *, 'ix0 =', ix0
                  print *, 'x(ix+ix0) =', x(ix + ix0)
                  print *, 'dx =', dx
                  print *, 'iRay,ix,jy,kz =', iRay, ix, jy, kz
                  stop
                end if
              end if

              if(sizeY > 1) then
                yr = ray(iRay, ix, jy, kz)%y

                if(yr < y(jy + jy0) - 0.5 * dy) then
                  print *, 'ERROR in shift_rayvol:'
                  print *, 'yr =', yr, '< y(jy+jy0) - 0.5*dy =', y(jy + jy0) &
                      &- 0.5 * dy
                  print *, 'jy =', jy
                  print *, 'jy0 =', jy0
                  print *, 'y(jy+jy0) =', y(jy + jy0)
                  print *, 'dy =', dy
                  print *, 'iRay,ix,jy,kz =', iRay, ix, jy, kz
                  stop
                end if

                if(yr > y(jy + jy0) + 0.5 * dy) then
                  print *, 'ERROR in shift_rayvol:'
                  print *, 'yr =', yr, '> y(jy+jy0) + 0.5*dy =', y(jy + jy0) &
                      &- 0.5 * dy
                  print *, 'jy =', jy
                  print *, 'jy0 =', jy0
                  print *, 'y(jy+jy0) =', y(jy + jy0)
                  print *, 'dy =', dy
                  print *, 'iRay,ix,jy,kz =', iRay, ix, jy, kz
                  stop
                end if
              end if

              zr = ray(iRay, ix, jy, kz)%z

              if(topography) then
                ! FJApr2023
                if(zr < zTFC(ix, jy, kz) - 0.5 * jac(ix, jy, kz) * dz) then
                  print *, "ERROR in shift_rayvol:"
                  print *, "zr =", zr, "< zTFC(ix, jy, kz) - 0.5 * jac(ix, jy, &
                      &kz) * dz =", zTFC(ix, jy, kz) - 0.5 * jac(ix, jy, kz) &
                      &* dz
                  print *, "iRay, ix, jy, kz =", iRay, ix, jy, kz
                  stop
                end if

                if(zr > zTFC(ix, jy, kz) + 0.5 * jac(ix, jy, kz) * dz) then
                  print *, "ERROR in shift_rayvol:"
                  print *, "zr =", zr, "> zTFC(ix, jy, kz) + 0.5 * jac(ix, jy, &
                      &kz) * dz =", zTFC(ix, jy, kz) + 0.5 * jac(ix, jy, kz) &
                      &* dz
                  print *, "iRay, ix, jy, kz =", iRay, ix, jy, kz
                  stop
                end if
              else
                if(zr < z(kz) - 0.5 * dz) then
                  print *, 'ERROR in shift_rayvol:'
                  print *, 'zr =', zr, '< z(kz) - 0.5*dz =', z(kz) - 0.5 * dz
                  print *, 'kz =', kz
                  print *, 'z(kz) =', z(kz)
                  print *, 'dz =', dz
                  print *, 'iRay,ix,jy,kz =', iRay, ix, jy, kz
                  stop
                end if

                if(zr > z(kz) + 0.5 * dz) then
                  print *, 'ERROR in shift_rayvol:'
                  print *, 'zr =', zr, '> z(kz) + 0.5*dz =', z(kz) - 0.5 * dz
                  print *, 'kz =', kz
                  print *, 'z(kz) =', z(kz)
                  print *, 'dz =', dz
                  print *, 'iRay,ix,jy,kz =', iRay, ix, jy, kz
                  stop
                end if
              end if
            end do
          end if
        end do
      end do
    end do
    ! teste

    return

  end subroutine shift_rayvol

  !----------------------------------------------------------------------

  subroutine merge_rayvol(ray)

    !-----------------------------------------------------------------
    ! merges ray volumes in cell if their number exceeds nray_max
    ! use logarithmic spacing in wave-number space
    !-----------------------------------------------------------------

    implicit none

    ! argument list
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(inout) :: ray

    integer :: ix, jy, kz

    integer, dimension(:), allocatable :: nr_merge

    integer :: ir_k, ir_l, ir_m

    integer :: jRay

    ! testb
    ! logical :: lmrglc
    ! teste

    real :: wnrk_min, wnrk_max, wnrl_min, wnrl_max, wnrm_min, wnrm_max

    real :: wnrk_min_p, wnrk_max_p, wnrl_min_p, wnrl_max_p, wnrm_min_p, &
        &wnrm_max_p
    real :: wnrk_min_n, wnrk_max_n, wnrl_min_n, wnrl_max_n, wnrm_min_n, &
        &wnrm_max_n
    real :: wnrk, wnrl, wnrm, wnrh
    real :: dwnrk, dwnrl, dwnrm

    real :: wnrk_1, wnrk_2, wnrl_1, wnrl_2, wnrm_1, wnrm_2

    real :: xr, yr, zr
    real :: dxr, dyr, dzr

    real :: wdr

    real :: axk, ayl, azm

    real :: dwnrk_mg, dwnrl_mg, dwnrm_mg

    real :: dwnrk_mg_p, dwnrl_mg_p, dwnrm_mg_p
    real :: dwnrk_mg_n, dwnrl_mg_n, dwnrm_mg_n

    real, dimension(:), allocatable :: xrmnmg, xrmxmg, yrmnmg, yrmxmg, zrmnmg, &
        &zrmxmg

    real, dimension(:), allocatable :: krmnmg, krmxmg, lrmnmg, lrmxmg, mrmnmg, &
        &mrmxmg

    real, dimension(:), allocatable :: wadrmg

    real :: fcpspx, fcpspy, fcpspz

    real :: omir, NNr

    real :: f_cor_nd

    real :: wa_old, wa_new

    real :: en_old, en_new
    real :: endn_old, endn_new

    real :: wnrt

    integer :: nrvtt0, nrvtt1, nrvloc

    if(steady_state) return

    allocate(nr_merge(nray_max))

    allocate(xrmnmg(nray_max))
    allocate(xrmxmg(nray_max))
    allocate(yrmnmg(nray_max))
    allocate(yrmxmg(nray_max))
    allocate(zrmnmg(nray_max))
    allocate(zrmxmg(nray_max))

    allocate(krmnmg(nray_max))
    allocate(krmxmg(nray_max))
    allocate(lrmnmg(nray_max))
    allocate(lrmxmg(nray_max))
    allocate(mrmnmg(nray_max))
    allocate(mrmxmg(nray_max))

    allocate(wadrmg(nray_max))

    if(sizeX > 1 .and. mod(nxRay, 2) /= 0) stop 'ERROR: nxRay must be even!'
    if(sizeY > 1 .and. mod(nyRay, 2) /= 0) stop 'ERROR: nyRay must be even!'
    if(sizeZ > 1 .and. mod(nzRay, 2) /= 0) stop 'ERROR: nzRay must be even!'

    ! total number of ray volumes before merging

    nrvloc = sum(nRay(1:nx, 1:ny, 1:nz))

    ! testb
    ! print*,'before merging nrvloc =',nrvloc
    ! teste

    call mpi_reduce(nrvloc, nrvtt0, 1, mpi_integer, mpi_sum, root, comm, ierror)
    call mpi_bcast(nrvtt0, 1, mpi_integer, root, comm, ierror)

    f_cor_nd = f_Coriolis_dim * tRef

    do kz = 1, nz
      do jy = 1, ny
        do ix = 1, nx
          ! testb
          ! lmrglc = .false.
          ! teste

          if(nRay(ix, jy, kz) <= nray_max) cycle

          !testb
          !lmrglc = .true.

          !if (lmrglc) then
          !   print*,'merging at ix,jy,kz =',ix,jy,kz

          !   wa_old = 0.
          !   wa_new = 0.

          !   en_old = 0.
          !   en_new = 0.

          !   endn_old = 0.
          !   endn_new = 0.
          !end if
          !teste

          ! minimum and mximum wave numbers of r.v. to be merged
          ! discriminate between negative and positive wavenumber range
          ! a central interval is to capture all wavenumbers close to
          ! zero

          wnrk_min_p = 0.
          wnrk_max_p = 0.
          wnrl_min_p = 0.
          wnrl_max_p = 0.
          wnrm_min_p = 0.
          wnrm_max_p = 0.

          wnrk_min_n = 0.
          wnrk_max_n = 0.
          wnrl_min_n = 0.
          wnrl_max_n = 0.
          wnrm_min_n = 0.
          wnrm_max_n = 0.

          do iRay = 1, nRay(ix, jy, kz)
            wnrk = ray(iRay, ix, jy, kz)%k
            dwnrk = abs(ray(iRay, ix, jy, kz)%dkray)

            wnrl = ray(iRay, ix, jy, kz)%l
            dwnrl = abs(ray(iRay, ix, jy, kz)%dlray)

            wnrm = ray(iRay, ix, jy, kz)%m
            dwnrm = abs(ray(iRay, ix, jy, kz)%dmray)

            if(sizeX > 1) then
              ! subdivide the interval of wavenumbers to be processed
              ! into up to three ranges:
              ! (1) if there are negative wavenumbers wnrk < 0, an
              ! interval -wnrk_max_n <= wnrk <= -wnrk_min_n
              ! (2) if there are zero wavnumbers wnrk = 0, a slot
              ! for those
              ! (3) if there are positive wavenumbers wnrk > 0, an
              ! interval wnrk_min_p <= wnrk <= wnrk_max_p

              if(wnrk > 0.0) then
                if(wnrk_min_p == 0.0) then
                  wnrk_min_p = wnrk
                else
                  wnrk_min_p = min(wnrk_min_p, wnrk)
                end if

                wnrk_max_p = max(wnrk_max_p, wnrk)
              elseif(wnrk < 0.0) then
                if(wnrk_min_n == 0.0) then
                  wnrk_min_n = - wnrk
                else
                  wnrk_min_n = min(wnrk_min_n, - wnrk)
                end if

                wnrk_max_n = max(wnrk_max_n, - wnrk)
              end if

              !wnrk_1 = wnrk - 0.5*dwnrk
              !wnrk_2 = wnrk + 0.5*dwnrk

              !if (wnrk_1 > 0.0) then
              !   ! both k-edges of the r.v. positive

              !   if (wnrk_min_p > 0.0) then
              !      wnrk_min_p = min(wnrk_min_p, wnrk_1)
              !     else
              !      wnrk_min_p = wnrk_1
              !   end if

              !   wnrk_max_p = max(wnrk_max_p,wnrk_2)
              !  elseif (wnrk_2 < 0.0) then
              !   ! both k-edges of the r.v. negative

              !   if (wnrk_min_n > 0.0) then
              !      wnrk_min_n = min(wnrk_min_n, -wnrk_2)
              !     else
              !      wnrk_min_n = -wnrk_2
              !   end if

              !   wnrk_max_n = max(wnrk_max_n,-wnrk_1)
              !  else
              !   ! k-edges with opposite sign

              !   !wnrk_max_n = max(wnrk_max_n,-wnrk_1)
              !   !wnrk_max_p = max(wnrk_max_p,wnrk_2)

              !   if (wnrk_max_n == 0.0) then
              !      wnrk_max_n = -wnrk_1
              !     elseif (wnrk_min_n == 0.0) then
              !      wnrk_min_n = min(wnrk_max_n, -wnrk_1)
              !      wnrk_max_n = max(wnrk_max_n, -wnrk_1)
              !     else
              !      wnrk_min_n = min(wnrk_min_n, -wnrk_1)
              !      wnrk_max_n = max(wnrk_max_n, -wnrk_1)
              !   end if

              !   if (wnrk_max_p == 0.0) then
              !      wnrk_max_p = wnrk_2
              !     elseif (wnrk_min_p == 0.0) then
              !      wnrk_min_p = min(wnrk_max_p, wnrk_2)
              !      wnrk_max_p = max(wnrk_max_p, wnrk_2)
              !     else
              !      wnrk_min_p = min(wnrk_min_p, wnrk_2)
              !      wnrk_max_p = max(wnrk_max_p, wnrk_2)
              !   end if
              !end if
            end if

            if(sizeY > 1) then
              if(wnrl > 0.0) then
                if(wnrl_min_p == 0.0) then
                  wnrl_min_p = wnrl
                else
                  wnrl_min_p = min(wnrl_min_p, wnrl)
                end if

                wnrl_max_p = max(wnrl_max_p, wnrl)
              elseif(wnrl < 0.0) then
                if(wnrl_min_n == 0.0) then
                  wnrl_min_n = - wnrl
                else
                  wnrl_min_n = min(wnrl_min_n, - wnrl)
                end if

                wnrl_max_n = max(wnrl_max_n, - wnrl)
              end if

              !wnrl_1 = wnrl - 0.5*dwnrl
              !wnrl_2 = wnrl + 0.5*dwnrl

              !if (wnrl_1 > 0.0) then
              !   ! both l-edges of the r.v. positive

              !   if (wnrl_min_p > 0.0) then
              !      wnrl_min_p = min(wnrl_min_p, wnrl_1)
              !     else
              !      wnrl_min_p = wnrl_1
              !   end if

              !   wnrl_max_p = max(wnrl_max_p,wnrl_2)
              !  elseif (wnrl_2 < 0.0) then
              !   ! both l-edges of the r.v. negative

              !   if (wnrl_min_n > 0.0) then
              !      wnrl_min_n = min(wnrl_min_n, -wnrl_2)
              !     else
              !      wnrl_min_n = -wnrl_2
              !   end if

              !   wnrl_max_n = max(wnrl_max_n,-wnrl_1)
              !  else
              !   ! l-edges with opposite sign

              !   ! wnrl_max_n = max(wnrl_max_n,-wnrl_1)
              !   ! wnrl_max_p = max(wnrl_max_p,wnrl_2)

              !   if (wnrl_max_n == 0.0) then
              !      wnrl_max_n = -wnrl_1
              !     elseif (wnrl_min_n == 0.0) then
              !      wnrl_min_n = min(wnrl_max_n, -wnrl_1)
              !      wnrl_max_n = max(wnrl_max_n, -wnrl_1)
              !     else
              !      wnrl_min_n = min(wnrl_min_n, -wnrl_1)
              !      wnrl_max_n = max(wnrl_max_n, -wnrl_1)
              !   end if

              !   if (wnrl_max_p == 0.0) then
              !      wnrl_max_p = wnrl_2
              !     elseif (wnrl_min_p == 0.0) then
              !      wnrl_min_p = min(wnrl_max_p, wnrl_2)
              !      wnrl_max_p = max(wnrl_max_p, wnrl_2)
              !     else
              !      wnrl_min_p = min(wnrl_min_p, wnrl_2)
              !      wnrl_max_p = max(wnrl_max_p, wnrl_2)
              !   end if
              !end if
            end if

            if(wnrm > 0.0) then
              if(wnrm_min_p == 0.0) then
                wnrm_min_p = wnrm
              else
                wnrm_min_p = min(wnrm_min_p, wnrm)
              end if

              wnrm_max_p = max(wnrm_max_p, wnrm)
            elseif(wnrm < 0.0) then
              if(wnrm_min_n == 0.0) then
                wnrm_min_n = - wnrm
              else
                wnrm_min_n = min(wnrm_min_n, - wnrm)
              end if

              wnrm_max_n = max(wnrm_max_n, - wnrm)
            end if

            !wnrm_1 = wnrm - 0.5*dwnrm
            !wnrm_2 = wnrm + 0.5*dwnrm

            !if (wnrm_1 > 0.0) then
            !   ! both m-edges of the r.v. positive

            !   if (wnrm_min_p > 0.0) then
            !      wnrm_min_p = min(wnrm_min_p, wnrm_1)
            !     else
            !      wnrm_min_p = wnrm_1
            !   end if

            !   wnrm_max_p = max(wnrm_max_p,wnrm_2)
            !  elseif (wnrm_2 < 0.0) then
            !   ! both m-edges of the r.v. negative

            !   if (wnrm_min_n > 0.0) then
            !      wnrm_min_n = min(wnrm_min_n, -wnrm_2)
            !     else
            !      wnrm_min_n = -wnrm_2
            !   end if

            !   wnrm_max_n = max(wnrm_max_n,-wnrm_1)
            !  else
            !   ! m-edges with opposite sign

            !   ! wnrm_max_n = max(wnrm_max_n,-wnrm_1)
            !   ! wnrm_max_p = max(wnrm_max_p,wnrm_2)

            !   if (wnrm_max_n == 0.0) then
            !      wnrm_max_n = -wnrm_1
            !     elseif (wnrm_min_n == 0.0) then
            !      wnrm_min_n = min(wnrm_max_n, -wnrm_1)
            !      wnrm_max_n = max(wnrm_max_n, -wnrm_1)
            !     else
            !      wnrm_min_n = min(wnrm_min_n, -wnrm_1)
            !      wnrm_max_n = max(wnrm_max_n, -wnrm_1)
            !   end if

            !   if (wnrm_max_p == 0.0) then
            !      wnrm_max_p = wnrm_2
            !     elseif (wnrm_min_p == 0.0) then
            !      wnrm_min_p = min(wnrm_max_p, wnrm_2)
            !      wnrm_max_p = max(wnrm_max_p, wnrm_2)
            !     else
            !      wnrm_min_p = min(wnrm_min_p, wnrm_2)
            !      wnrm_max_p = max(wnrm_max_p, wnrm_2)
            !   end if
            !end if

            !testb
            !if (sizeX > 1) then
            !   fcpspx = ray(iRay,ix,jy,kz)%area_xk
            !  else
            !   fcpspx = 1.0
            !end if

            !if (sizeY > 1) then
            !   fcpspy = ray(iRay,ix,jy,kz)%area_yl
            !  else
            !   fcpspy = 1.0
            !end if

            !fcpspz = ray(iRay,ix,jy,kz)%area_zm

            !wa_old &
            != wa_old &
            !  + ray(iRay,ix,jy,kz)%dens * fcpspx * fcpspy * fcpspz

            !en_old &
            != en_old &
            !  + ray(iRay,ix,jy,kz)%dens * ray(iRay,ix,jy,kz)%omega &
            !    * fcpspx * fcpspy * fcpspz

            !endn_old &
            != endn_old &
            !  + ray(iRay,ix,jy,kz)%dens * ray(iRay,ix,jy,kz)%omega &
            !    * fcpspx * fcpspy * fcpspz / rhoStrat(kz)
            !teste
          end do

          if(sizeX > 1) then
            if(wnrk_min_n == 0.0 .and. wnrk_max_n == 0.0) then
              if(wnrk_min_p /= 0.0 .and. wnrk_max_p /= 0.0) then
                wnrk_min_n = wnrk_min_p
                wnrk_max_n = wnrk_max_p
              else
                ! all limits zero only applies if all wnrk = 0
                ! hence, just in order to provide some numbers ...
                wnrk_min_n = 1.0
                wnrk_max_n = 2.0
              end if
            end if

            if(wnrk_min_p == 0.0 .and. wnrk_max_p == 0.0) then
              if(wnrk_min_n /= 0.0 .and. wnrk_max_n /= 0.0) then
                wnrk_min_p = wnrk_min_n
                wnrk_max_p = wnrk_max_n
              else
                ! all limits zero only applies if all wnrk = 0
                ! hence, just in order to provide some numbers ...
                wnrk_min_p = 1.0
                wnrk_max_p = 2.0
              end if
            end if

            ! in order to prevent zero-width intervals ...

            if(wnrk_min_n == wnrk_max_n) then
              wnrk_min_n = 0.5 * wnrk_min_n
              wnrk_max_n = 2.0 * wnrk_max_n
            end if

            if(wnrk_min_p == wnrk_max_p) then
              wnrk_min_p = 0.5 * wnrk_min_p
              wnrk_max_p = 2.0 * wnrk_max_p
            end if

            !if (wnrk_min_n == 0.0) then
            !   if (wnrk_min_p == 0.0) then
            !      wnrk_min_n = wnrk_max_n/nxRay
            !     else
            !      wnrk_min_n = wnrk_min_p
            !    end if
            !end if

            !if (wnrk_max_n == 0.0) then
            !   if (wnrk_max_p == 0.0) then
            !      wnrk_max_n = wnrk_min_n*nxRay
            !     else
            !      wnrk_max_n = wnrk_max_p
            !    end if
            !end if

            !if (wnrk_min_p == 0.0) then
            !   if (wnrk_min_n == 0.0) then
            !      wnrk_min_p = wnrk_max_p/nxRay
            !     else
            !      wnrk_min_p = wnrk_min_n
            !    end if
            !end if

            !if (wnrk_max_p == 0.0) then
            !   if (wnrk_max_n == 0.0) then
            !      wnrk_max_p = wnrk_min_p*nxRay
            !     else
            !      wnrk_max_p = wnrk_max_n
            !    end if
            !end if
          end if

          if(sizeY > 1) then
            if(wnrl_min_n == 0.0 .and. wnrl_max_n == 0.0) then
              if(wnrl_min_p /= 0.0 .and. wnrl_max_p /= 0.0) then
                wnrl_min_n = wnrl_min_p
                wnrl_max_n = wnrl_max_p
              else
                ! all limits zero only applies if all wnrl = 0
                ! hence, just in order to provide some numbers ...
                wnrl_min_n = 1.0
                wnrl_max_n = 2.0
              end if
            end if

            if(wnrl_min_p == 0.0 .and. wnrl_max_p == 0.0) then
              if(wnrl_min_n /= 0.0 .and. wnrl_max_n /= 0.0) then
                wnrl_min_p = wnrl_min_n
                wnrl_max_p = wnrl_max_n
              else
                ! all limits zero only applies if all wnrl = 0
                ! hence, just in order to provide some numbers ...
                wnrl_min_p = 1.0
                wnrl_max_p = 2.0
              end if
            end if

            ! in order to prevent zero-width intervals ...

            if(wnrl_min_n == wnrl_max_n) then
              wnrl_min_n = 0.5 * wnrl_min_n
              wnrl_max_n = 2.0 * wnrl_max_n
            end if

            if(wnrl_min_p == wnrl_max_p) then
              wnrl_min_p = 0.5 * wnrl_min_p
              wnrl_max_p = 2.0 * wnrl_max_p
            end if

            !if (wnrl_min_n == 0.0) then
            !   if (wnrl_min_p == 0.0) then
            !      wnrl_min_n = wnrl_max_n/nyRay
            !     else
            !      wnrl_min_n = wnrl_min_p
            !    end if
            !end if

            !if (wnrl_max_n == 0.0) then
            !   if (wnrl_max_p == 0.0) then
            !      wnrl_max_n = wnrl_min_n*nyRay
            !     else
            !      wnrl_max_n = wnrl_max_p
            !    end if
            !end if

            !if (wnrl_min_p == 0.0) then
            !   if (wnrl_min_n == 0.0) then
            !      wnrl_min_p = wnrl_max_p/nyRay
            !     else
            !      wnrl_min_p = wnrl_min_n
            !    end if
            !end if

            !if (wnrl_max_p == 0.0) then
            !   if (wnrl_max_n == 0.0) then
            !      wnrl_max_p = wnrl_min_p*nyRay
            !     else
            !      wnrl_max_p = wnrl_max_n
            !    end if
            !end if
          end if

          if(wnrm_min_n == 0.0 .and. wnrm_max_n == 0.0) then
            if(wnrm_min_p /= 0.0 .and. wnrm_max_p /= 0.0) then
              wnrm_min_n = wnrm_min_p
              wnrm_max_n = wnrm_max_p
            else
              ! all limits zero only applies if all wnrm = 0
              ! hence, just in order to provide some numbers ...
              wnrm_min_n = 1.0
              wnrm_max_n = 2.0
            end if
          end if

          if(wnrm_min_p == 0.0 .and. wnrm_max_p == 0.0) then
            if(wnrm_min_n /= 0.0 .and. wnrm_max_n /= 0.0) then
              wnrm_min_p = wnrm_min_n
              wnrm_max_p = wnrm_max_n
            else
              ! all limits zero only applies if all wnrm = 0
              ! hence, just in order to provide some numbers ...
              wnrm_min_p = 1.0
              wnrm_max_p = 2.0
            end if
          end if

          ! in order to prevent zero-width intervals ...

          if(wnrm_min_n == wnrm_max_n) then
            wnrm_min_n = 0.5 * wnrm_min_n
            wnrm_max_n = 2.0 * wnrm_max_n
          end if

          if(wnrm_min_p == wnrm_max_p) then
            wnrm_min_p = 0.5 * wnrm_min_p
            wnrm_max_p = 2.0 * wnrm_max_p
          end if

          !if (wnrm_min_n == 0.0) then
          !   if (wnrm_min_p == 0.0) then
          !      wnrm_min_n = wnrm_max_n/nyRay
          !     else
          !      wnrm_min_n = wnrm_min_p
          !    end if
          !end if

          !if (wnrm_max_n == 0.0) then
          !   if (wnrm_max_p == 0.0) then
          !      wnrm_max_n = wnrm_min_n*nyRay
          !     else
          !      wnrm_max_n = wnrm_max_p
          !    end if
          !end if

          !if (wnrm_min_p == 0.0) then
          !   if (wnrm_min_n == 0.0) then
          !      wnrm_min_p = wnrm_max_p/nyRay
          !     else
          !      wnrm_min_p = wnrm_min_n
          !    end if
          !end if

          !if (wnrm_max_p == 0.0) then
          !   if (wnrm_max_n == 0.0) then
          !      wnrm_max_p = wnrm_min_p*nyRay
          !     else
          !      wnrm_max_p = wnrm_max_n
          !    end if
          !end if

          ! for each wavenumber direction and sign, the logarithmic
          ! width of the intervals within which the r.v. are to be
          ! merged

          ! each sign gets same amount of intervals,
          ! e.g. nxRay/2 - 1 for positive and negative k each
          ! a central interval for wavenumbers around zero is provided
          ! as well

          if(sizeX > 1) then
            if(wnrk_max_n == 0.0) stop 'ERROR: wnrk_max_n = 0'
            if(wnrk_min_n == 0.0) stop 'ERROR: wnrk_min_n = 0'
            if(wnrk_max_p == 0.0) stop 'ERROR: wnrk_max_p = 0'
            if(wnrk_min_p == 0.0) stop 'ERROR: wnrk_min_p = 0'
            if(nxRay < 3) stop 'ERROR: nxRay < 3'

            dwnrk_mg_n = log(wnrk_max_n / wnrk_min_n) / (nxRay / 2 - 1)
            dwnrk_mg_p = log(wnrk_max_p / wnrk_min_p) / (nxRay / 2 - 1)
          end if

          if(sizeY > 1) then
            if(wnrl_max_n == 0.0) stop 'ERROR: wnrl_max_n = 0'
            if(wnrl_min_n == 0.0) stop 'ERROR: wnrl_min_n = 0'
            if(wnrl_max_p == 0.0) stop 'ERROR: wnrl_max_p = 0'
            if(wnrl_min_p == 0.0) stop 'ERROR: wnrl_min_p = 0'
            if(nyRay < 3) stop 'ERROR: nyRay < 3'

            dwnrl_mg_n = log(wnrl_max_n / wnrl_min_n) / (nyRay / 2 - 1)
            dwnrl_mg_p = log(wnrl_max_p / wnrl_min_p) / (nyRay / 2 - 1)
          end if

          if(wnrm_max_n == 0.0) stop 'ERROR: wnrm_max_n = 0'
          if(wnrm_min_n == 0.0) stop 'ERROR: wnrm_min_n = 0'
          if(wnrm_max_p == 0.0) stop 'ERROR: wnrm_max_p = 0'
          if(wnrm_min_p == 0.0) stop 'ERROR: wnrm_min_p = 0'
          if(nzRay < 3) stop 'ERROR: nzRay < 3'

          dwnrm_mg_n = log(wnrm_max_n / wnrm_min_n) / (nzRay / 2 - 1)
          dwnrm_mg_p = log(wnrm_max_p / wnrm_min_p) / (nzRay / 2 - 1)

          ! generate merged r.v.

          nr_merge = 0

          do iRay = 1, nRay(ix, jy, kz)
            wnrk = ray(iRay, ix, jy, kz)%k
            dwnrk = ray(iRay, ix, jy, kz)%dkray

            wnrl = ray(iRay, ix, jy, kz)%l
            dwnrl = ray(iRay, ix, jy, kz)%dlray

            wnrm = ray(iRay, ix, jy, kz)%m
            dwnrm = ray(iRay, ix, jy, kz)%dmray

            xr = ray(iRay, ix, jy, kz)%x
            dxr = ray(iRay, ix, jy, kz)%dxray

            yr = ray(iRay, ix, jy, kz)%y
            dyr = ray(iRay, ix, jy, kz)%dyray

            zr = ray(iRay, ix, jy, kz)%z
            dzr = ray(iRay, ix, jy, kz)%dzray

            axk = ray(iRay, ix, jy, kz)%area_xk
            ayl = ray(iRay, ix, jy, kz)%area_yl
            azm = ray(iRay, ix, jy, kz)%area_zm

            wdr = ray(iRay, ix, jy, kz)%dens

            omir = ray(iRay, ix, jy, kz)%omega

            if(sizeX > 1) then
              fcpspx = axk

              ! nxRay - 1 intervals for k:
              ! indices 1 ... nxRay/2 - 1 for negative k
              ! index nxRay/2 for k = 0
              ! indices nxRay/2 + 1 ... nxRay - 1 for positive k

              if(wnrk < 0.0) then
                if(abs(log(- wnrk / wnrk_max_n) / dwnrk_mg_n) < 1.d-3) then
                  ir_k = nxRay / 2 - 1
                else
                  ir_k = int(log(- wnrk / wnrk_min_n) / dwnrk_mg_n) + 1
                end if
              elseif(wnrk == 0.0) then
                ir_k = nxRay / 2
              else
                if(abs(log(wnrk / wnrk_max_p) / dwnrk_mg_p) < 1.d-3) then
                  ir_k = nxRay - 1
                else
                  ir_k = int(log(wnrk / wnrk_min_p) / dwnrk_mg_p) + nxRay / 2 &
                      &+ 1
                end if
              end if

              if(ir_k < 1) then
                print *, 'ERROR in merge_rayvol: ir_k =', ir_k, '< 1'
                print *, 'wnrk = ', wnrk
                print *, 'wnrk_min_n = ', wnrk_min_n
                print *, 'wnrk_max_n = ', wnrk_max_n
                print *, 'wnrk_min_p = ', wnrk_min_p
                print *, 'wnrk_max_p = ', wnrk_max_p
                stop
              elseif(ir_k > nxRay - 1) then
                print *, 'ERROR in merge_rayvol: ir_k =', ir_k, '> nxRay - 1 &
                    &=', nxRay - 1
                print *, 'wnrk = ', wnrk
                print *, 'wnrk_min_n = ', wnrk_min_n
                print *, 'wnrk_max_n = ', wnrk_max_n
                print *, 'wnrk_min_p = ', wnrk_min_p
                print *, 'wnrk_max_p = ', wnrk_max_p
                stop
              end if
            else
              fcpspx = 1.0
              ir_k = 1
            end if

            if(sizeY > 1) then
              fcpspy = ayl

              if(wnrl < 0.0) then
                if(abs(log(- wnrl / wnrl_max_n) / dwnrl_mg_n) < 1.d-3) then
                  ir_l = nyRay / 2 - 1
                else
                  ir_l = int(log(- wnrl / wnrl_min_n) / dwnrl_mg_n) + 1
                end if
              elseif(wnrl == 0.0) then
                ir_l = nyRay / 2
              else
                if(abs(log(wnrl / wnrl_max_p) / dwnrl_mg_p) < 1.d-3) then
                  ir_l = nyRay - 1
                else
                  ir_l = int(log(wnrl / wnrl_min_p) / dwnrl_mg_p) + nyRay / 2 &
                      &+ 1
                end if
              end if

              if(ir_l < 1) then
                print *, 'ERROR in merge_rayvol: ir_l =', ir_l, '< 1'
                print *, 'wnrl = ', wnrl
                print *, 'wnrl_min_n = ', wnrl_min_n
                print *, 'wnrl_max_n = ', wnrl_max_n
                print *, 'wnrl_min_p = ', wnrl_min_p
                print *, 'wnrl_max_p = ', wnrl_max_p
                stop
              elseif(ir_l > nyRay - 1) then
                print *, 'ERROR in merge_rayvol: ir_l =', ir_l, '> nyRay - 1 &
                    &=', nyRay - 1
                print *, 'wnrl = ', wnrl
                print *, 'wnrl_min_n = ', wnrl_min_n
                print *, 'wnrl_max_n = ', wnrl_max_n
                print *, 'wnrl_min_p = ', wnrl_min_p
                print *, 'wnrl_max_p = ', wnrl_max_p
                stop
              end if
            else
              fcpspy = 1.0
              ir_l = 1
            end if

            fcpspz = azm

            if(wnrm < 0.0) then
              if(abs(log(- wnrm / wnrm_max_n) / dwnrm_mg_n) < 1.d-3) then
                ir_m = nzRay / 2 - 1
              else
                ir_m = int(log(- wnrm / wnrm_min_n) / dwnrm_mg_n) + 1
              end if
            elseif(wnrm == 0.0) then
              ir_m = nzRay / 2
            else
              if(abs(log(wnrm / wnrm_max_p) / dwnrm_mg_p) < 1.d-3) then
                ir_m = nzRay - 1
              else
                ir_m = int(log(wnrm / wnrm_min_p) / dwnrm_mg_p) + nzRay / 2 + 1
              end if
            end if

            if(ir_m < 1) then
              print *, 'ERROR in merge_rayvol: ir_m =', ir_m, '< 1'
              print *, 'wnrm = ', wnrm
              print *, 'wnrm_min_n = ', wnrm_min_n
              print *, 'wnrm_max_n = ', wnrm_max_n
              print *, 'wnrm_min_p = ', wnrm_min_p
              print *, 'wnrm_max_p = ', wnrm_max_p
              print *, 'wnrm_max_p - wnrm = ', wnrm_max_p - wnrm
              stop
            elseif(ir_m > nzRay - 1) then
              print *, 'ERROR in merge_rayvol: ir_m =', ir_m, '> nzRay - 1 =', &
                  &nzRay - 1
              print *, 'wnrm = ', wnrm
              print *, 'wnrm_min_n = ', wnrm_min_n
              print *, 'wnrm_max_n = ', wnrm_max_n
              print *, 'wnrm_min_p = ', wnrm_min_p
              print *, 'wnrm_max_p = ', wnrm_max_p
              print *, 'wnrm_max_p - wnrm = ', wnrm_max_p - wnrm
              stop
            end if

            !if (sizeX > 1) then
            !   if (sizeY > 1) then
            !      jRay &
            !      = (ir_m - 1) * (nyRay - 1)*(nxRay - 2) &
            !        + (ir_l - 1) * (nxRay -  1) + ir_k
            !     else
            !      jRay = (ir_m - 1) * (nxRay - 1) + ir_k
            !   end if
            !  else
            !   if (sizeY > 1) then
            !      jRay = (ir_m - 1) * (nyRay - 1) + ir_l
            !     else
            !      jRay = ir_m
            !   end if
            !end if

            if(sizeX > 1) then
              if(sizeY > 1) then
                jRay = (ir_m - 1) * (nyRay - 1) * (nxRay - 1) + (ir_l - 1) &
                    &* (nxRay - 1) + ir_k
              else
                jRay = (ir_m - 1) * (nxRay - 1) + ir_k
              end if
            else
              if(sizeY > 1) then
                jRay = (ir_m - 1) * (nyRay - 1) + ir_l
              else
                jRay = ir_m
              end if
            end if

            !testb
            !if (lmrglc) then
            !   print*,'iRay = ',iRay
            !   print*,'wnrk = ',wnrk
            !   print*,'wnrm = ',wnrm
            !   print*,'omir =',omir
            !end if
            !teste

            nr_merge(jRay) = nr_merge(jRay) + 1

            if(nr_merge(jRay) == 1) then
              xrmnmg(jRay) = xr - 0.5 * dxr
              xrmxmg(jRay) = xr + 0.5 * dxr

              yrmnmg(jRay) = yr - 0.5 * dyr
              yrmxmg(jRay) = yr + 0.5 * dyr

              zrmnmg(jRay) = zr - 0.5 * dzr
              zrmxmg(jRay) = zr + 0.5 * dzr

              krmnmg(jRay) = wnrk - 0.5 * dwnrk
              krmxmg(jRay) = wnrk + 0.5 * dwnrk

              lrmnmg(jRay) = wnrl - 0.5 * dwnrl
              lrmxmg(jRay) = wnrl + 0.5 * dwnrl

              mrmnmg(jRay) = wnrm - 0.5 * dwnrm
              mrmxmg(jRay) = wnrm + 0.5 * dwnrm

              if(cons_merge == "wa") then
                ! wave-action density after merging to be determined
                ! such that the wave action remains the same,
                ! hence ...

                wadrmg(jRay) = wdr * fcpspx * fcpspy * fcpspz
              elseif(cons_merge == "en") then
                ! wave-action density after merging to be determined
                ! such that the wave energy remains the same,
                ! hence ...

                wadrmg(jRay) = wdr * omir * fcpspx * fcpspy * fcpspz
              else
                stop 'wrong cons_merge in merge_rayvol'
              end if
            else
              xrmnmg(jRay) = min(xrmnmg(jRay), xr - 0.5 * dxr)
              xrmxmg(jRay) = max(xrmxmg(jRay), xr + 0.5 * dxr)

              yrmnmg(jRay) = min(yrmnmg(jRay), yr - 0.5 * dyr)
              yrmxmg(jRay) = max(yrmxmg(jRay), yr + 0.5 * dyr)

              zrmnmg(jRay) = min(zrmnmg(jRay), zr - 0.5 * dzr)
              zrmxmg(jRay) = max(zrmxmg(jRay), zr + 0.5 * dzr)

              krmnmg(jRay) = min(krmnmg(jRay), wnrk - 0.5 * dwnrk)
              krmxmg(jRay) = max(krmxmg(jRay), wnrk + 0.5 * dwnrk)

              lrmnmg(jRay) = min(lrmnmg(jRay), wnrl - 0.5 * dwnrl)
              lrmxmg(jRay) = max(lrmxmg(jRay), wnrl + 0.5 * dwnrl)

              mrmnmg(jRay) = min(mrmnmg(jRay), wnrm - 0.5 * dwnrm)
              mrmxmg(jRay) = max(mrmxmg(jRay), wnrm + 0.5 * dwnrm)

              if(cons_merge == "wa") then
                ! wave-action density after merging to be determined
                ! such that the wave action remains the same,
                ! hence ...

                wadrmg(jRay) = wadrmg(jRay) + wdr * fcpspx * fcpspy * fcpspz
              elseif(cons_merge == "en") then
                ! wave-action density after merging to be determined
                ! such that the wave energy remains the same,
                ! hence ...

                wadrmg(jRay) = wadrmg(jRay) + wdr * omir * fcpspx * fcpspy &
                    &* fcpspz
              else
                stop 'wrong cons_merge in merge_rayvol'
              end if
            end if
          end do

          ! replace old r.v. by the merged r.v.

          !testb
          ! print*,'merged ray volumes:'
          !teste

          iRay = 0

          do jRay = 1, nray_max
            if(nr_merge(jRay) < 1) cycle

            iRay = iRay + 1

            !testb
            ! if (lmrglc) then
            !    print*,'jRay = ',jRay
            !    print*,'iRay = ',iRay
            ! end if
            !teste

            ! position and width in physical space

            ray(iRay, ix, jy, kz)%x = 0.5 * (xrmxmg(jRay) + xrmnmg(jRay))
            ray(iRay, ix, jy, kz)%dxray = xrmxmg(jRay) - xrmnmg(jRay)

            ray(iRay, ix, jy, kz)%y = 0.5 * (yrmxmg(jRay) + yrmnmg(jRay))
            ray(iRay, ix, jy, kz)%dyray = yrmxmg(jRay) - yrmnmg(jRay)

            ray(iRay, ix, jy, kz)%z = 0.5 * (zrmxmg(jRay) + zrmnmg(jRay))
            ray(iRay, ix, jy, kz)%dzray = zrmxmg(jRay) - zrmnmg(jRay)

            ! position and width in wavenumber space

            !if (krmxmg(jRay) * krmnmg(jRay) > 0.0) then
            !   ray(iRay,ix,jy,kz)%k &
            !   = sqrt(krmxmg(jRay) * krmnmg(jRay))
            !  else
            !   ray(iRay,ix,jy,kz)%k &
            !   = 0.5*(krmxmg(jRay) + krmnmg(jRay))
            !end if

            ray(iRay, ix, jy, kz)%k = 0.5 * (krmxmg(jRay) + krmnmg(jRay))
            ray(iRay, ix, jy, kz)%dkray = krmxmg(jRay) - krmnmg(jRay)

            !if (lrmxmg(jRay) * lrmnmg(jRay) > 0.0) then
            !   ray(iRay,ix,jy,kz)%l &
            !   = sqrt(lrmxmg(jRay) * lrmnmg(jRay))
            !  else
            !   ray(iRay,ix,jy,kz)%l &
            !   = 0.5*(lrmxmg(jRay) + lrmnmg(jRay))
            !end if

            ray(iRay, ix, jy, kz)%l = 0.5 * (lrmxmg(jRay) + lrmnmg(jRay))
            ray(iRay, ix, jy, kz)%dlray = lrmxmg(jRay) - lrmnmg(jRay)

            !if (mrmxmg(jRay) * mrmnmg(jRay) > 0.0) then
            !   ray(iRay,ix,jy,kz)%m &
            !   = sqrt(mrmxmg(jRay) * mrmnmg(jRay))
            !  else
            !   ray(iRay,ix,jy,kz)%m &
            !   = 0.5*(mrmxmg(jRay) + mrmnmg(jRay))
            !end if

            ray(iRay, ix, jy, kz)%m = 0.5 * (mrmxmg(jRay) + mrmnmg(jRay))
            ray(iRay, ix, jy, kz)%dmray = mrmxmg(jRay) - mrmnmg(jRay)

            ! intrinsic frequency

            wnrk = ray(iRay, ix, jy, kz)%k
            wnrl = ray(iRay, ix, jy, kz)%l
            wnrm = ray(iRay, ix, jy, kz)%m

            wnrh = sqrt(wnrk ** 2 + wnrl ** 2)

            zr = ray(iRay, ix, jy, kz)%z

            call stratification(zr, 1, NNr)

            omir = branchr * sqrt(NNr * wnrh ** 2 + f_cor_nd ** 2 * wnrm ** 2) &
                &/ sqrt(wnrh ** 2 + wnrm ** 2)

            ray(iRay, ix, jy, kz)%omega = omir

            ! phase-space volumes

            ray(iRay, ix, jy, kz)%area_xk = ray(iRay, ix, jy, kz)%dxray &
                &* ray(iRay, ix, jy, kz)%dkray

            ray(iRay, ix, jy, kz)%area_yl = ray(iRay, ix, jy, kz)%dyray &
                &* ray(iRay, ix, jy, kz)%dlray

            ray(iRay, ix, jy, kz)%area_zm = ray(iRay, ix, jy, kz)%dzray &
                &* ray(iRay, ix, jy, kz)%dmray

            ! wave-action density

            if(sizeX > 1) then
              fcpspx = ray(iRay, ix, jy, kz)%area_xk
            else
              fcpspx = 1.0
            end if

            if(sizeY > 1) then
              fcpspy = ray(iRay, ix, jy, kz)%area_yl
            else
              fcpspy = 1.0
            end if

            fcpspz = ray(iRay, ix, jy, kz)%area_zm

            if(cons_merge == "wa") then
              ! wave-action density after merging to be determined
              ! such that the wave action remains the same,
              ! hence ...

              ray(iRay, ix, jy, kz)%dens = wadrmg(jRay) / (fcpspx * fcpspy &
                  &* fcpspz)
            elseif(cons_merge == "en") then
              ! wave-action density after merging to be determined
              ! such that the wave energy remains the same, hence ...

              ray(iRay, ix, jy, kz)%dens = wadrmg(jRay) / (omir * fcpspx &
                  &* fcpspy * fcpspz)
            else
              stop 'wrong cons_merge in merge_rayvol'
            end if

            !SDJul2024
            !set phase to zero after merging
            if(update_phase) ray(iRay, ix, jy, kz)%dphi = 0

            !testb
            !if (sizeX > 1) then
            !   fcpspx = ray(iRay,ix,jy,kz)%area_xk
            !  else
            !   fcpspx = 1.0
            !end if

            !if (sizeY > 1) then
            !   fcpspy = ray(iRay,ix,jy,kz)%area_yl
            !  else
            !   fcpspy = 1.0
            !end if

            !fcpspz = ray(iRay,ix,jy,kz)%area_zm

            !wa_new &
            != wa_new &
            !  + ray(iRay,ix,jy,kz)%dens * fcpspx * fcpspy * fcpspz

            !en_new &
            != en_new &
            !  + ray(iRay,ix,jy,kz)%dens * ray(iRay,ix,jy,kz)%omega &
            !    * fcpspx * fcpspy * fcpspz

            !endn_new &
            != endn_new &
            !  + ray(iRay,ix,jy,kz)%dens * ray(iRay,ix,jy,kz)%omega &
            !    * fcpspx * fcpspy * fcpspz / rhoStrat(kz)
            !teste

            !testb
            !if (lmrglc) then
            !   print*,'iRay = ',iRay
            !   print*,'wnrk = ',wnrk
            !   print*,'wnrm = ',wnrm
            !   print*,'omir =',omir
            !end if
            !teste
          end do

          !testb
          !if (lmrglc) then
          !   print*,' '
          !   print*,'after merging at ix, jy, kz =',ix, jy, kz
          !   print*,'old nRay =', nRay(ix,jy,kz)
          !   print*,'new nRay =', iRay
          !   print*,'while nxRay, nyRay, nzRay =', nxRay,nyRay,nzRay
          !   print*,'and hence nray_max =', nray_max
          !   print*,'old wave action =',wa_old
          !   print*,'new wave action =',wa_new
          !   print*,'old energy =',en_old
          !   print*,'new energy =',en_new
          !   print*,'old density-normalized energy =',endn_old
          !   print*,'new density-normalized energy =',endn_new
          !end if
          !teste

          if(iRay > nray_max) then
            print *, 'after merging at ix, jy, kz =', ix, jy, kz, 'nRay =', &
                &iRay, '> nray_max =', nray_max
            stop
          else
            nRay(ix, jy, kz) = iRay
          end if
        end do
      end do
    end do

    ! total number of ray volumes before after merging

    nrvloc = sum(nRay(1:nx, 1:ny, 1:nz))

    ! testb
    ! print*,'after merging nrvloc =',nrvloc
    ! teste

    call mpi_reduce(nrvloc, nrvtt1, 1, mpi_integer, mpi_sum, root, comm, ierror)
    call mpi_bcast(nrvtt1, 1, mpi_integer, root, comm, ierror)

    if(master .and. nrvtt1 < nrvtt0) then
      print *, 'after merging nray =', nrvtt1
    end if

    return

  end subroutine merge_rayvol

  !----------------------------------------------------------------------

  subroutine transport_rayvol(var, ray, dt, rkStage, time)

    !-------------------------------------------------------------
    ! integrates the eikonal equations for ray-volume position and
    ! wave number,
    ! as well as those for the ray-volume extents in position and
    ! wave-number space
    !-------------------------------------------------------------

    implicit none

    ! argument list
    type(var_type), intent(in) :: var
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(inout) :: ray

    real, intent(in) :: dt
    integer, intent(in) :: rkStage

    ! FJFeb2023
    real, intent(in) :: time

    integer :: ix, jy, kz
    integer :: ix2, jy2, kz2, ik, jl, km

    integer :: nskip

    integer :: nrlc

    integer :: ix0, jy0

    real :: F
    real :: dkdt, dldt, dmdt

    real :: wnrk, wnrl, wnrm
    real :: wnrh
    real :: omir, omir1, omir2

    real :: cgirx, cgiry, cgirz
    real :: cgirx1, cgiry1, cgirz1
    real :: cgirx2, cgiry2, cgirz2

    real :: cgrx, cgry, cgrz
    real :: cgrx1, cgry1, cgrz1
    real :: cgrx2, cgry2, cgrz2

    real :: uxr, vyr, wzr, NNr, dNNdzr
    real :: uxr1, vyr1, wzr1, NNr1
    real :: uxr2, vyr2, wzr2, NNr2

    real :: f_cor_nd

    real :: dk_ini_nd, dl_ini_nd, dm_ini_nd

    real :: pspvol

    real :: xr, yr, zr
    real :: xr1, yr1, zr1
    real :: xr2, yr2, zr2
    real :: dxr, dyr, dzr

    real :: dudxr, dudyr, dudzr
    real :: dvdxr, dvdyr, dvdzr

    real :: wnrk_init, wnrl_init, wnrm_init, wnrh_init

    real :: xr0, yr0, zr0
    real :: sigwpx, sigwpy, sigwpz

    real :: NN_nd

    real :: wnk_0, wnl_0, wnm_0

    real :: amp, displm

    real :: ddxdt, ddydt, ddzdt

    ! Long number (FJJan2023)
    real :: long

    ! Wave mode (FJApr2023)
    integer :: iwm

    ! Relaxation parameters (FJJun2023)
    real :: alphaSponge, betaSponge
    real :: spongeAlphaZ, spongeDz
    real :: spongeAlphaX, spongeAlphaY

    integer :: kz0
    real :: cgirz0

    real :: integral1, integral2, m2b2, m2b2k2, diffusion
    real :: dwnrk, dwnrl, dwnrm, dxi, dyi, dzi
    real :: facpsp

    ix0 = is + nbx - 1
    jy0 = js + nby - 1

    f_cor_nd = f_Coriolis_dim * tRef

    ! Update height.
    if(topography .and. topographyTime > 0.0) then
      do ix = - nbx, nx + nbx
        do jy = - nby, ny + nby
          do kz = - nbz, nz + nbz
            zTFC(ix, jy, kz) = heightTFC(ix, jy, kz)
            zTildeTFC(ix, jy, kz) = heightTFC(ix, jy, kz) + 0.5 * jac(ix, jy, &
                &kz) * dz
          end do
        end do
      end do
    end if

    if(case_wkb == 3 .and. steady_state) then
      call orographic_source(var, ray, time, stepFrac(RKStage) * dt)
    end if

    if(steady_state) then
      do kz = 1, nz
        do jy = 1, ny
          do ix = 1, nx
            ! Set ray-volume count.
            nRay(ix, jy, kz) = nRay(ix, jy, kz - 1)

            ! Set up saturation computation.
            integral1 = 0.0
            integral2 = 0.0
            m2b2 = 0.0
            m2b2k2 = 0.0

            ! Loop over ray volumes.
            do iRay = 1, nRay(ix, jy, kz)

              ! Prepare ray volume.
              ray(iRay, ix, jy, kz) = ray(iRay, ix, jy, kz - 1)

              ! Skip modes with zero wave-action density.
              if(ray(iRay, ix, jy, kz - 1)%dens == 0.0) cycle

              ! Set vertical position (and extent).
              if(topography) then
                ray(iRay, ix, jy, kz)%z = ray(iRay, ix, jy, kz - 1)%z + 0.5 &
                    &* (jac(ix, jy, kz - 1) + jac(ix, jy, kz)) * dz
                ray(iRay, ix, jy, kz)%dzray = ray(iRay, ix, jy, kz - 1)%dzray &
                    &* jac(ix, jy, kz) / jac(ix, jy, kz - 1)
              else
                ray(iRay, ix, jy, kz)%z = ray(iRay, ix, jy, kz - 1)%z + dz
              end if

              ! Get horizontal wavenumbers.
              wnrk = ray(iRay, ix, jy, kz)%k
              wnrl = ray(iRay, ix, jy, kz)%l
              wnrh = sqrt(wnrk ** 2.0 + wnrl ** 2.0)

              ! Set reference level.
              kz0 = max(1, kz - 1)

              ! Compute vertical group velocity at the level below.
              call stratification(ray(iRay, ix, jy, kz0)%z, 1, NN_nd)
              omir = ray(iRay, ix, jy, kz0)%omega
              if(branchr * omir > f_cor_nd .and. branchr * omir < sqrt(NN_nd)) &
                  &then
                wnrm = ray(iRay, ix, jy, kz0)%m
                cgirz0 = wnrm * (f_cor_nd ** 2 - NN_nd) * wnrh ** 2 / omir &
                    &/ (wnrh ** 2 + wnrm ** 2) ** 2
              else
                ray(iRay, ix, jy, kz)%dens = 0.0
                cycle
              end if

              ! Compute local intrinsic frequency, vertical
              ! wavenumber and vertical group velocity.
              call stratification(ray(iRay, ix, jy, kz)%z, 1, NN_nd)
              omir = - 0.5 * (var%u(ix, jy, kz) + var%u(ix - 1, jy, kz)) &
                  &* wnrk - 0.5 * (var%v(ix, jy, kz) + var%v(ix, jy - 1, kz)) &
                  &* wnrl
              if(branchr * omir > f_cor_nd .and. branchr * omir < sqrt(NN_nd)) &
                  &then
                wnrm = - branchr * sqrt(wnrh ** 2 * (NN_nd - omir ** 2) &
                    &/ (omir ** 2 - f_cor_nd ** 2))
                cgirz = wnrm * (f_cor_nd ** 2 - NN_nd) * wnrh ** 2 / omir &
                    &/ (wnrh ** 2 + wnrm ** 2) ** 2
              else
                ray(iRay, ix, jy, kz)%dens = 0.0
                cycle
              end if

              ! Set local intrinsic frequency and vertical wavenumber.
              ray(iRay, ix, jy, kz)%omega = omir
              ray(iRay, ix, jy, kz)%m = wnrm

              ! Set local wave action density.
              if(spongeLayer .and. unifiedSponge) then
                xr = ray(iRay, ix, jy, kz)%x
                yr = ray(iRay, ix, jy, kz)%y
                zr = ray(iRay, ix, jy, kz)%z
                alphaSponge = 2.0 * interpolate_sponge(xr, yr, zr)
                if(topography) then
                  ray(iRay, ix, jy, kz)%dens = 1.0 / (1.0 + alphaSponge &
                      &/ cgirz * 0.5 * (jac(ix, jy, kz - 1) + jac(ix, jy, kz)) &
                      &* dz) * cgirz0 * ray(iRay, ix, jy, kz0)%dens / cgirz
                else
                  ray(iRay, ix, jy, kz)%dens = 1.0 / (1.0 + alphaSponge &
                      &/ cgirz * dz) * cgirz0 * ray(iRay, ix, jy, kz0)%dens &
                      &/ cgirz
                end if
              else
                ray(iRay, ix, jy, kz)%dens = cgirz0 * ray(iRay, ix, jy, &
                    &kz0)%dens / cgirz
              end if

              ! Get ray volume extents.
              dxr = ray(iRay, ix, jy, kz)%dxray
              dyr = ray(iRay, ix, jy, kz)%dyray
              dzr = ray(iRay, ix, jy, kz)%dzray
              dwnrk = ray(iRay, ix, jy, kz)%dkray
              dwnrl = ray(iRay, ix, jy, kz)%dlray
              dwnrm = ray(iRay, ix, jy, kz)%dmray

              ! Compute phase space factor.
              if(topography) then
                dzi = min(dzr, jac(ix, jy, kz) * dz)
                facpsp = dzi / jac(ix, jy, kz) / dz * dwnrm
              else
                dzi = min(dzr, dz)
                facpsp = dzi / dz * dwnrm
              end if
              if(sizeX > 1) then
                dxi = min(dxr, dx)
                facpsp = facpsp * dxi / dx * dwnrk
              end if
              if(sizeY > 1) then
                dyi = min(dyr, dy)
                facpsp = facpsp * dyi / dy * dwnrl
              end if

              ! Update saturation amplitude.
              integral1 = wnrh ** 2 * wnrm ** 2 / ((wnrh ** 2 + wnrm ** 2) &
                  &* omir) * facpsp
              if(topography) then
                m2b2 = m2b2 + 2.0 * NN_nd ** 2 / rhoStratTFC(ix, jy, kz) &
                    &* integral1 * ray(iRay, ix, jy, kz)%dens
              else
                m2b2 = m2b2 + 2.0 * NN_nd ** 2 / rhoStrat(kz) * integral1 &
                    &* ray(iRay, ix, jy, kz)%dens
              end if
              integral2 = wnrh ** 2 * wnrm ** 2 / omir * facpsp
              if(topography) then
                m2b2k2 = m2b2k2 + 2.0 * NN_nd ** 2 / rhoStratTFC(ix, jy, kz) &
                    &* integral2 * ray(iRay, ix, jy, kz)%dens * jac(ix, jy, &
                    &kz) * dz / cgirz
              else
                m2b2k2 = m2b2k2 + 2.0 * NN_nd ** 2 / rhoStrat(kz) * integral2 &
                    &* ray(iRay, ix, jy, kz)%dens * dz / cgirz
              end if
            end do
            ! Compute diffusion coefficient
            if(topography) then
              call stratification(zTFC(ix, jy, kz), 1, NN_nd)
            else
              call stratification(z(kz), 1, NN_nd)
            end if
            if(m2b2k2 == 0.0 .or. m2b2 < alpha_sat ** 2 * NN_nd ** 2) then
              diffusion = 0.0
            else
              diffusion = (m2b2 - alpha_sat ** 2 * NN_nd ** 2) / (2.0 * m2b2k2)
            end if
            ! Reduce wave action density.
            do iRay = 1, nRay(ix, jy, kz)
              if(ray(iRay, ix, jy, kz)%dens == 0.0) cycle
              wnrk = ray(iRay, ix, jy, kz)%k
              wnrl = ray(iRay, ix, jy, kz)%l
              wnrm = ray(iRay, ix, jy, kz)%m
              wnrh = sqrt(wnrk ** 2 + wnrl ** 2)
              call stratification(ray(iRay, ix, jy, kz)%z, 1, NN_nd)
              omir = ray(iRay, ix, jy, kz)%omega
              if(branchr * omir > f_cor_nd .and. branchr * omir < sqrt(NN_nd)) &
                  &then
                cgirz = wnrm * (f_cor_nd ** 2 - NN_nd) * wnrh ** 2 / omir &
                    &/ (wnrh ** 2 + wnrm ** 2) ** 2
              else
                ray(iRay, ix, jy, kz)%dens = 0.0
                cycle
              end if
              if(topography) then
                ray(iRay, ix, jy, kz)%dens = ray(iRay, ix, jy, kz)%dens &
                    &* max(0.0, 1.0 - 0.5 * (jac(ix, jy, kz - 1) + jac(ix, jy, &
                    &kz)) * dz / cgirz * 2.0 * diffusion * (wnrh ** 2 + wnrm &
                    &** 2))
              else
                ray(iRay, ix, jy, kz)%dens = ray(iRay, ix, jy, kz)%dens &
                    &* max(0.0, 1.0 - dz / cgirz * 2.0 * diffusion * (wnrh &
                    &** 2 + wnrm ** 2))
              end if
            end do
          end do
        end do
      end do
      if(sizeX > 1) call setboundary_rayvol_x(ray)
      if(sizeY > 1) call setboundary_rayvol_y(ray)
      return
    end if

    ! initialize RK-tendencies at first RK stage
    if(RKstage == 1) then
      dxRay = 0.0
      dkRay = 0.0
      ddxRay = 0.0
    end if

    !SD
    ! init RK phase tendency
    if(include_ice) then
      if(RKstage == 1) then
        dpRay = 0.0
      end if
    end if

    cgx_max = 0.0
    cgy_max = 0.0
    cgz_max = 0.0

    if(case_wkb == 3) then
      kz0 = 0
    else
      kz0 = 1
    end if

    do kz = kz0, nz
      do jy = 1, ny
        do ix = 1, nx
          nskip = 0

          if(nRay(ix, jy, kz) < 1) cycle

          ray_loop: do iRay = 1, nRay(ix, jy, kz)
            wnrk = ray(iRay, ix, jy, kz)%k
            wnrl = ray(iRay, ix, jy, kz)%l
            wnrm = ray(iRay, ix, jy, kz)%m

            wnrh = sqrt(wnrk ** 2 + wnrl ** 2)

            xr = ray(iRay, ix, jy, kz)%x
            dxr = ray(iRay, ix, jy, kz)%dxray
            xr1 = xr - 0.5 * dxr
            xr2 = xr + 0.5 * dxr

            yr = ray(iRay, ix, jy, kz)%y
            dyr = ray(iRay, ix, jy, kz)%dyray
            yr1 = yr - 0.5 * dyr
            yr2 = yr + 0.5 * dyr

            zr = ray(iRay, ix, jy, kz)%z
            dzr = ray(iRay, ix, jy, kz)%dzray
            zr1 = zr - 0.5 * dzr
            zr2 = zr + 0.5 * dzr

            !skip ray volumes that have left the domain
            if(case_wkb /= 3) then
              if(topography) then
                ! FJApr2023
                if(zr1 < topography_surface(ix, jy) - jac(ix, jy, 0) * dz) then
                  nskip = nskip + 1
                  cycle
                end if
              else
                if(zr1 < lz(0) - dz) then
                  nskip = nskip + 1
                  cycle
                end if
              end if
            end if

            call stratification(zr1, 1, NNr1)
            call stratification(zr, 1, NNr)
            call stratification(zr2, 1, NNr2)

            omir1 = branchr * sqrt(NNr1 * wnrh ** 2 + f_cor_nd ** 2 * wnrm &
                &** 2) / sqrt(wnrh ** 2 + wnrm ** 2)

            omir = branchr * sqrt(NNr * wnrh ** 2 + f_cor_nd ** 2 * wnrm ** 2) &
                &/ sqrt(wnrh ** 2 + wnrm ** 2)

            ! testb
            if(NNr2 <= 0.0) then
              print *, 'NNr2 =', NNr2, '<= 0.0 at'
              print *, 'zr2 =', zr2, 'from'
              print *, 'ray(iRay,ix,jy,kz)%z =', ray(iRay, ix, jy, kz)%z
              print *, 'ray(iRay,ix,jy,kz)%dzray =', ray(iRay, ix, jy, kz)%dzray
              print *, 'iRay,ix,jy,kz =', iRay, ix, jy, kz
              stop
            end if

            if(wnrh <= 0.0) then
              print *, 'wnrh =', wnrh, '<= 0.0 from'
              print *, 'wnrk =', wnrk
              print *, 'wnrl =', wnrl
              print *, 'iRay,ix,jy,kz =', iRay, ix, jy, kz
              stop
            end if
            ! teste

            omir2 = branchr * sqrt(NNr2 * wnrh ** 2 + f_cor_nd ** 2 * wnrm &
                &** 2) / sqrt(wnrh ** 2 + wnrm ** 2)

            ray(iRay, ix, jy, kz)%omega = omir

            ! intrinsic group velocities at the respective edges of
            ! the ray volumes

            if(sizeX > 1) then
              ! intrinsic group velocity in x direction not depending
              ! on x
              cgirx = wnrk * (NNr - omir ** 2) / (omir * (wnrh ** 2 + wnrm &
                  &** 2))
            end if

            if(sizeY > 1) then
              ! intrinsic group velocity in y direction not depending
              ! on y
              cgiry = wnrl * (NNr - omir ** 2) / (omir * (wnrh ** 2 + wnrm &
                  &** 2))
            end if

            ! intrinsic vertical group velocity depending on z
            ! (via the stratification)
            cgirz1 = - wnrm * (omir1 ** 2 - f_cor_nd ** 2) / (omir1 * (wnrh &
                &** 2 + wnrm ** 2))
            cgirz2 = - wnrm * (omir2 ** 2 - f_cor_nd ** 2) / (omir2 * (wnrh &
                &** 2 + wnrm ** 2))

            !--------------------------------
            !     displacement in x direction
            !--------------------------------

            ! RK update

            if(sizeX > 1 .and. kz > 0) then
              call meanflow(xr1, yr, zr, var, 1, uxr1)
              call meanflow(xr2, yr, zr, var, 1, uxr2)

              ! group velocity in x direction at the two edges in x
              cgrx1 = cgirx + uxr1
              cgrx2 = cgirx + uxr2

              ! group velocity in x direction for the carrier ray
              cgrx = 0.5 * (cgrx1 + cgrx2)

              F = cgrx !allow horizontal ray propagation
              dxRay(1, iRay, ix, jy, kz) = dt * F + alphaRK(rkStage) &
                  &* dxRay(1, iRay, ix, jy, kz)
              ray(iRay, ix, jy, kz)%x = ray(iRay, ix, jy, kz)%x &
                  &+ betaRK(RKstage) * dxRay(1, iRay, ix, jy, kz)

              ! update maximum group velocity in x direction
              cgx_max = max(cgx_max, abs(cgrx))
            end if

            !-----------------------------
            !     displacement in y direction
            !-----------------------------

            ! RK update

            if(sizeY > 1 .and. kz > 0) then
              call meanflow(xr, yr1, zr, var, 2, vyr1)
              call meanflow(xr, yr2, zr, var, 2, vyr2)

              ! group velocity in y direction at the two edges in y
              cgry1 = cgiry + vyr1
              cgry2 = cgiry + vyr2

              ! group velocity in y direction for the carrier ray
              cgry = 0.5 * (cgry1 + cgry2)

              F = cgry !allow horizontal ray propagation
              dxRay(2, iRay, ix, jy, kz) = dt * F + alphaRK(rkStage) &
                  &* dxRay(2, iRay, ix, jy, kz)
              ray(iRay, ix, jy, kz)%y = ray(iRay, ix, jy, kz)%y &
                  &+ betaRK(RKstage) * dxRay(2, iRay, ix, jy, kz)

              ! update maximum group velocity in y direction
              cgy_max = max(cgy_max, abs(cgry))
            end if

            !-----------------------------
            !     vertical displacement
            !-----------------------------

            ! RK update

            ! in line with the asymptotic results the vertcal wind is
            ! NOT added to the intrinsic vertical group velocity
            ! should one want to change this, one would also have to
            ! take the vertical-wind gradient into account in the
            ! prognostic equations for the wave number

            ! call meanflow(xr,yr,zr1,var,3,wzr1)
            ! call meanflow(xr,yr,zr2,var,3,wzr2)

            ! group velocity in z direction at the two edges in z
            cgrz1 = cgirz1 ! + wzr1
            cgrz2 = cgirz2 ! + wzr2

            ! group velocity in z direction for the carrier ray
            cgrz = 0.5 * (cgrz1 + cgrz2)

            F = cgrz
            dxRay(3, iRay, ix, jy, kz) = dt * F + alphaRK(rkStage) * dxRay(3, &
                &iRay, ix, jy, kz)
            ray(iRay, ix, jy, kz)%z = ray(iRay, ix, jy, kz)%z &
                &+ betaRK(RKstage) * dxRay(3, iRay, ix, jy, kz)

            ! update maximum group velocity in z direction
            cgz_max = max(cgz_max, abs(cgrz))

            !-------------------------------
            !    change of wavenumber
            !-------------------------------

            ! wave refraction only above lz(0) + zmin_wkb
            if(zr > lz(0) + zmin_wkb) then
              ! RK procedure

              call meanflow(xr, yr, zr, var, 11, dudxr)
              call meanflow(xr, yr, zr, var, 12, dudyr)
              call meanflow(xr, yr, zr, var, 13, dudzr)

              call meanflow(xr, yr, zr, var, 21, dvdxr)
              call meanflow(xr, yr, zr, var, 22, dvdyr)
              call meanflow(xr, yr, zr, var, 23, dvdzr)

              if(zr < lz(0) - dz) then
                print *, 'ERROR IN setup_wkb: LOWER EDGE OF RAY  VOLUME', &
                    &iRay, ix, jy, kz, 'TOO LOW'
                stop
              end if

              if(zr < lz(0) - dz) then
                print *, 'ERROR IN transport_rayvol: RAY VOLUME', iRay, ix, &
                    &jy, kz, 'TOO LOW'
                stop
              end if

              call stratification(zr, 2, dnndzr)

              dkdt = - dudxr * wnrk - dvdxr * wnrl
              dldt = - dudyr * wnrk - dvdyr * wnrl
              dmdt = - dudzr * wnrk - dvdzr * wnrl - wnrh ** 2 * dnndzr / (2.0 &
                  &* omir + (wnrh ** 2 + wnrm ** 2))

              dkRay(1, iRay, ix, jy, kz) = dt * dkdt + alphaRK(rkStage) &
                  &* dkRay(1, iRay, ix, jy, kz)

              dkRay(2, iRay, ix, jy, kz) = dt * dldt + alphaRK(rkStage) &
                  &* dkRay(2, iRay, ix, jy, kz)

              dkRay(3, iRay, ix, jy, kz) = dt * dmdt + alphaRK(rkStage) &
                  &* dkRay(3, iRay, ix, jy, kz)

              ray(iRay, ix, jy, kz)%k = ray(iRay, ix, jy, kz)%k &
                  &+ betaRK(rkStage) * dkRay(1, iRay, ix, jy, kz)

              ray(iRay, ix, jy, kz)%l = ray(iRay, ix, jy, kz)%l &
                  &+ betaRK(rkStage) * dkRay(2, iRay, ix, jy, kz)

              ray(iRay, ix, jy, kz)%m = ray(iRay, ix, jy, kz)%m &
                  &+ betaRK(rkStage) * dkRay(3, iRay, ix, jy, kz)

              !----------------------------------------------
              !    change of wave-number width of ray volumes
              !----------------------------------------------

              ! dk

              if(sizeX > 1 .and. kz > 0) then
                ddxdt = cgrx2 - cgrx1

                ddxRay(1, iRay, ix, jy, kz) = dt * ddxdt + alphaRK(rkStage) &
                    &* ddxRay(1, iRay, ix, jy, kz)

                ray(iRay, ix, jy, kz)%dxray = ray(iRay, ix, jy, kz)%dxray &
                    &+ betaRK(rkStage) * ddxRay(1, iRay, ix, jy, kz)

                if(ray(iRay, ix, jy, kz)%dxray <= 0.0) then
                  print *, 'dxray(', iRay, ix, jy, kz, ') <= 0.0  ==> time &
                      &step too large?'
                  ray(iRay, ix, jy, kz)%dxray = - ray(iRay, ix, jy, kz)%dxray
                end if

                ray(iRay, ix, jy, kz)%dkray = ray(iRay, ix, jy, kz)%area_xk &
                    &/ ray(iRay, ix, jy, kz)%dxray
              end if

              ! dl

              if(sizeY > 1 .and. kz > 0) then
                ddydt = cgry2 - cgry1

                ddxRay(2, iRay, ix, jy, kz) = dt * ddydt + alphaRK(rkStage) &
                    &* ddxRay(2, iRay, ix, jy, kz)

                ray(iRay, ix, jy, kz)%dyray = ray(iRay, ix, jy, kz)%dyray &
                    &+ betaRK(rkStage) * ddxRay(2, iRay, ix, jy, kz)

                if(ray(iRay, ix, jy, kz)%dyray <= 0.0) then
                  print *, 'dyray(', iRay, ix, jy, kz, ') <= 0.0  ==> time &
                      &step too large?'
                  ray(iRay, ix, jy, kz)%dyray = - ray(iRay, ix, jy, kz)%dyray
                end if

                ray(iRay, ix, jy, kz)%dlray = ray(iRay, ix, jy, kz)%area_yl &
                    &/ ray(iRay, ix, jy, kz)%dyray
              end if

              ! dm

              ddzdt = cgrz2 - cgrz1

              ddxRay(3, iRay, ix, jy, kz) = dt * ddzdt + alphaRK(rkStage) &
                  &* ddxRay(3, iRay, ix, jy, kz)

              ray(iRay, ix, jy, kz)%dzray = ray(iRay, ix, jy, kz)%dzray &
                  &+ betaRK(rkStage) * ddxRay(3, iRay, ix, jy, kz)

              if(ray(iRay, ix, jy, kz)%dzray <= 0.0) then
                print *, 'dzray(', iRay, ix, jy, kz, ') <= 0.0  ==> time step &
                    &too large?'
                ray(iRay, ix, jy, kz)%dzray = - ray(iRay, ix, jy, kz)%dzray
              end if

              ray(iRay, ix, jy, kz)%dmray = ray(iRay, ix, jy, kz)%area_zm &
                  &/ ray(iRay, ix, jy, kz)%dzray
            end if

            ! SD
            !-----------------------------------
            ! update phase
            !-----------------------------------

            if(include_ice .and. update_phase) then

              call meanflow(xr, yr, zr, var, 1, uxr)
              call meanflow(xr, yr, zr, var, 2, vyr)

              !NB intrinsic group vel. orthogonal to wavenumber vector
              dpRay(iRay, ix, jy, kz) = - dt * (omir + 2. * (uxr * wnrk + vyr &
                  &* wnrl)) + alphaRK(rkStage) * dpRay(iRay, ix, jy, kz)
              !if((ix .eq. 1) .and. (jy .eq. 1) .and. (kz .eq. 1) .and. (iRay &
              !    .eq. 1)) then
              !  print *, 'dpRay', rkStage, dpRay(iRay, ix, jy, kz)
              !end if
              ray(iRay, ix, jy, kz)%dphi = ray(iRay, ix, jy, kz)%dphi &
                  &+ betaRK(rkStage) * dpRay(iRay, ix, jy, kz)

            end if
            !-----------------------------------
            ! update of the intrinsic frequency
            !-----------------------------------

            wnrk = ray(iRay, ix, jy, kz)%k
            wnrl = ray(iRay, ix, jy, kz)%l
            wnrm = ray(iRay, ix, jy, kz)%m

            wnrh = sqrt(wnrk ** 2 + wnrl ** 2)

            zr = ray(iRay, ix, jy, kz)%z

            call stratification(zr, 1, NNr)

            omir = branchr * sqrt(NNr * wnrh ** 2 + f_cor_nd ** 2 * wnrm ** 2) &
                &/ sqrt(wnrh ** 2 + wnrm ** 2)

            ray(iRay, ix, jy, kz)%omega = omir
          end do ray_loop

          if(nskip > 0) then
            print *, nskip, 'r.v. skipped in transport_rayvol out of', &
                &nRay(ix, jy, kz)
          end if
        end do ! ix
      end do ! jy
    end do ! kz

    ! FJJun2023
    ! Sponge layer
    if(spongeLayer .and. unifiedSponge) then
      do kz = 1, nz
        do jy = 1, ny
          do ix = 1, nx
            do iRay = 1, nRay(ix, jy, kz)
              xr = ray(iRay, ix, jy, kz)%x
              yr = ray(iRay, ix, jy, kz)%y
              zr = ray(iRay, ix, jy, kz)%z
              alphaSponge = 2.0 * interpolate_sponge(xr, yr, zr)
              betaSponge = 1.0 / (1.0 + alphaSponge * stepFrac(RKStage) * dt)
              ray(iRay, ix, jy, kz)%dens = betaSponge * ray(iRay, ix, jy, &
                  &kz)%dens
            end do
          end do
        end do
      end do
    end if

    if(case_wkb == 3 .and. .not. steady_state) then
      call orographic_source(var, ray, time, stepFrac(RKStage) * dt)
    end if

  end subroutine transport_rayvol

  !---------------------------------------------------------------------

  subroutine boundary_rayvol(ray)

    ! Update ray volumes at the boundaries.

    implicit none

    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(inout) :: ray

    if(steady_state) return

    ! Apply boundary conditions.
    if(sizeX > 1) call setboundary_rayvol_x(ray)
    if(sizeY > 1) call setboundary_rayvol_y(ray)
    call setboundary_rayvol_z(ray)

  end subroutine boundary_rayvol

  !---------------------------------------------------------------------

  subroutine saturation_3D(ray, dt, diffusioncoeff)

    ! Original subroutine: saturation by G. Boeloeni (2016)
    ! ---------------------------------------------------
    ! account for wave saturation
    ! if the buoyancy amplitude square m^2*B^2 at a given vertical layer
    ! i exceeds the threshold for static instability N^4, then the wave
    ! action density of each contributing ray is reduced so that:
    ! m^2*B^2 = alpha^2 N^4

    ! Extension on quasi-2D ray propagation by J. Wilhelm (12/2016)
    ! ---------------------------------------------------
    ! take x-position of ray volumes into account, 2D-field of reduction of
    ! wave action density (in each gridbox the WAD is reduced).

    ! 3D and non-dimensionalization by U. Achatz
    ! ------------------------------------------

    ! turbulent eddy diffusivity K(x,y,z) given in output diffusioncoeff
    ! this is necessary for mixing of tracer
    !-------------------------------------------

    implicit none

    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(inout) :: ray

    ! time step
    real, intent(in) :: dt

    ! give eddy diffusity K as output
    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(inout) :: &
        &diffusioncoeff

    ! indices, etc.
    integer iRay, kzmax, kzmin

    integer :: ix, jy, kz

    integer :: ix0, jy0

    integer :: nrlc

    ! ray-volume wave-action density
    real densr
    ! ray-volume proportion on the x-, y-, and z-grid
    real dxi, dyi, dzi
    ! m^2 * B^2
    real mB2(0:nx + 1, 0:ny + 1, 0:sizeZ + 1)
    ! m^2 * B^2 * (kh^2+m^2)
    real mB2K2(0:nx + 1, 0:ny + 1, 0:sizeZ + 1)
    ! diffusivity coefficient
    real diffusion(0:nx + 1, 0:ny + 1, 0:sizeZ + 1), kappa
    ! variables for integrals over m
    real integral1, integral2

    real :: wnrk, wnrl, wnrm, wnrhs, NN_nd, dwnrk, dwnrl, dwnrm, omir

    real :: xr, yr, zr
    real :: dxr, dyr, dzr

    real :: facpsp
    logical :: apply

    if(steady_state) return

    ix0 = is + nbx - 1
    jy0 = js + nby - 1

    ! compute saturation amplitude (m^2*B^2) for each level

    mB2 = 0.0
    mB2K2 = 0.0

    do kzrv = 1, nz
      do jyrv = 1, ny
        do ixrv = 1, nx
          if(nRay(ixrv, jyrv, kzrv) < 1) cycle

          do iRay = 1, nRay(ixrv, jyrv, kzrv)
            ! skip counting ray volumes with zero wave-action density

            if(ray(iRay, ixrv, jyrv, kzrv)%dens == 0.0) cycle

            xr = ray(iRay, ixrv, jyrv, kzrv)%x
            yr = ray(iRay, ixrv, jyrv, kzrv)%y
            zr = ray(iRay, ixrv, jyrv, kzrv)%z

            dxr = ray(iRay, ixrv, jyrv, kzrv)%dxray
            dyr = ray(iRay, ixrv, jyrv, kzrv)%dyray
            dzr = ray(iRay, ixrv, jyrv, kzrv)%dzray

            ! implement horizontal boundary conditions for ray-volume
            ! positions

            if(sizeX > 1) then
              if(xr < lx(0)) then
                select case(xBoundary)
                case("periodic")
                  xr = lx(1) + mod(xr - lx(0), lx(1) - lx(0))
                case default
                  stop "saturation_3D: unknown case xBoundary"
                end select
              elseif(xr > lx(1)) then
                select case(xBoundary)
                case("periodic")
                  xr = lx(0) + mod(xr - lx(1), lx(1) - lx(0))
                case default
                  stop "saturation_3D: unknown case xBoundary"
                end select
              end if

              ix = floor((xr - lx(0)) / dx) + 1 - ix0
            else
              ix = 1
            end if

            if(sizeY > 1) then
              if(yr < ly(0)) then
                select case(yBoundary)
                case("periodic")
                  yr = ly(1) + mod(yr - ly(0), ly(1) - ly(0))
                case default
                  stop "saturation_3D: unknown case yBoundary"
                end select
              elseif(yr > ly(1)) then
                select case(yBoundary)
                case("periodic")
                  yr = ly(0) + mod(yr - ly(1), ly(1) - ly(0))
                case default
                  stop "saturation_3D: unknown case yBoundary"
                end select
              end if

              jy = floor((yr - ly(0)) / dy) + 1 - jy0
            else
              jy = 1
            end if

            ! vertical boundary conditions:
            ! zBoundary = 'periodic': implement periodicity
            ! zBoundary = 'solid_wall': skip counting ray volumes
            !                           that have completely left the
            !                           model domain
            apply = .false.
            if(.not. topography .and. zr < lz(0)) then
              apply = .true.
            else if(topography) then
              if(zr < topography_surface(ix, jy)) then
                apply = .true.
              end if
            end if
            if(apply) then
              select case(zBoundary)
              case("periodic")
                zr = lz(1) + mod(zr - lz(0), lz(1) - lz(0))
              case("solid_wall")
                cycle
              case default
                stop "saturation_3D: unknown case zBoundary"
              end select
            elseif(zr > lz(1)) then
              select case(zBoundary)
              case("periodic")
                zr = lz(0) + mod(zr - lz(1), lz(1) - lz(0))
              case("solid_wall")
                cycle
              case default
                stop "saturation_3D: unknown case zBoundary"
              end select
            end if

            ! Skip rays propagating out of the domain.
            ! Why half-levels?
            if(topography) then
              kz = floor((levelTFC(ix, jy, zr) - lz(0)) / dz) + 1
              if(kz < 1 .or. kz > sizeZ) cycle
            else
              kz = floor((zr - lz(0)) / dz) + 1
              if(kz < 1 .or. kz > sizeZ) cycle
            end if

            ! Compute stratification.
            call stratification(zr, 1, NN_nd)

            wnrk = ray(iRay, ixrv, jyrv, kzrv)%k
            wnrl = ray(iRay, ixrv, jyrv, kzrv)%l
            wnrm = ray(iRay, ixrv, jyrv, kzrv)%m

            wnrhs = wnrk ** 2 + wnrl ** 2

            dwnrk = ray(iRay, ixrv, jyrv, kzrv)%dkray
            dwnrl = ray(iRay, ixrv, jyrv, kzrv)%dlray
            dwnrm = ray(iRay, ixrv, jyrv, kzrv)%dmray

            omir = ray(iRay, ixrv, jyrv, kzrv)%omega

            densr = ray(iRay, ixrv, jyrv, kzrv)%dens

            ! spatial extension of ray to be taken into account
            if(topography) then
              dzi = min(dzr, jac(ix, jy, kz) * dz)
              facpsp = dzi / jac(ix, jy, kz) / dz * dwnrm
            else
              dzi = min(dzr, dz)
              facpsp = dzi / dz * dwnrm
            end if

            if(sizeX > 1) then
              dxi = min(dxr, dx)
              facpsp = facpsp * dxi / dx * dwnrk
            end if

            if(sizeY > 1) then
              dyi = min(dyr, dy)
              facpsp = facpsp * dyi / dy * dwnrl
            end if

            integral1 = wnrhs * wnrm ** 2 / ((wnrhs + wnrm ** 2) * omir) &
                &* facpsp

            if(topography) then
              mB2(ix, jy, kz) = mB2(ix, jy, kz) + 2.0 * NN_nd ** 2 &
                  &/ rhoStratTFC(ix, jy, kz) * densr * integral1
            else
              mB2(ix, jy, kz) = mB2(ix, jy, kz) + 2.0 * NN_nd ** 2 &
                  &/ rhoStrat(kz) * densr * integral1
            end if

            integral2 = wnrhs * wnrm ** 2 / omir * facpsp

            if(topography) then
              mB2K2(ix, jy, kz) = mB2K2(ix, jy, kz) + 2.0 * NN_nd ** 2 &
                  &/ rhoStratTFC(ix, jy, kz) * densr * integral2
            else
              mB2K2(ix, jy, kz) = mB2K2(ix, jy, kz) + 2.0 * NN_nd ** 2 &
                  &/ rhoStrat(kz) * densr * integral2
            end if
          end do
        end do ! ixrv
      end do ! jyrv
    end do ! kzrv

    ! loop for computing the diffusivity coefficient
    diffusioncoeff = 0.0
    if(topography) then
      do kz = 1, sizeZ
        do jy = 1, ny
          do ix = 1, nx
            call stratification(zTFC(ix, jy, kz), 1, NN_nd)
            if(mB2K2(ix, jy, kz) == 0.0 .or. mB2(ix, jy, kz) < alpha_sat ** 2 &
                &* NN_nd ** 2) then
              diffusion(ix, jy, kz) = 0.0
            else
              diffusion(ix, jy, kz) = (mB2(ix, jy, kz) - alpha_sat ** 2 &
                  &* NN_nd ** 2) / (2.0 * dt * mB2K2(ix, jy, kz))
            endif
            if(include_tracer) then
              diffusioncoeff(ix, jy, kz) = diffusion(ix, jy, kz)
            end if
          end do
        end do
      end do
    else
      do kz = 1, sizeZ
        call stratification(z(kz), 1, NN_nd)
        do jy = 1, ny
          do ix = 1, nx
            if(mB2K2(ix, jy, kz) == 0.0 .or. mB2(ix, jy, kz) < alpha_sat ** 2 &
                &* NN_nd ** 2) then
              diffusion(ix, jy, kz) = 0.0
            else
              diffusion(ix, jy, kz) = (mB2(ix, jy, kz) - alpha_sat ** 2 &
                  &* NN_nd ** 2) / (2.0 * dt * mB2K2(ix, jy, kz))
            endif
            if(include_tracer) then
              diffusioncoeff(ix, jy, kz) = diffusion(ix, jy, kz)
            end if
          end do
        end do
      end do
    end if

    ! loop for reducing wave action density
    ! if m^2*B^2 exceeds saturation threshold
    do kzrv = 1, nz
      do jyrv = 1, ny
        do ixrv = 1, nx
          if(nRay(ixrv, jyrv, kzrv) < 1) cycle

          do iRay = 1, nRay(ixrv, jyrv, kzrv)
            ! skip counting ray volumes with zero wave-action density

            if(ray(iRay, ixrv, jyrv, kzrv)%dens == 0.0) cycle

            xr = ray(iRay, ixrv, jyrv, kzrv)%x
            yr = ray(iRay, ixrv, jyrv, kzrv)%y
            zr = ray(iRay, ixrv, jyrv, kzrv)%z

            dxr = ray(iRay, ixrv, jyrv, kzrv)%dxray
            dyr = ray(iRay, ixrv, jyrv, kzrv)%dyray
            dzr = ray(iRay, ixrv, jyrv, kzrv)%dzray

            ! implement horizontal boundary conditions for ray-volume
            ! positions

            if(sizeX > 1) then
              if(xr < lx(0)) then
                select case(xBoundary)
                case("periodic")
                  xr = lx(1) + mod(xr - lx(0), lx(1) - lx(0))
                case default
                  stop "saturation_3D: unknown case xBoundary"
                end select
              elseif(xr > lx(1)) then
                select case(xBoundary)
                case("periodic")
                  xr = lx(0) + mod(xr - lx(1), lx(1) - lx(0))
                case default
                  stop "saturation_3D: unknown case xBoundary"
                end select
              end if

              ix = floor((xr - lx(0)) / dx) + 1 - ix0
            else
              ix = 1
            end if

            if(sizeY > 1) then
              if(yr < ly(0)) then
                select case(yBoundary)
                case("periodic")
                  yr = ly(1) + mod(yr - ly(0), ly(1) - ly(0))
                case default
                  stop "saturation_3D: unknown case yBoundary"
                end select
              elseif(yr > ly(1)) then
                select case(yBoundary)
                case("periodic")
                  yr = ly(0) + mod(yr - ly(1), ly(1) - ly(0))
                case default
                  stop "saturation_3D: unknown case yBoundary"
                end select
              end if

              jy = floor((yr - ly(0)) / dy) + 1 - jy0
            else
              jy = 1
            end if

            ! vertical boundary conditions:
            ! zBoundary = 'periodic': implement periodicity
            ! zBoundary = 'solid_wall': skip counting ray volumes
            !                           that have completely left the
            !                           model domain

            apply = .false.
            if(.not. topography .and. zr < lz(0)) then
              apply = .true.
            else if(topography) then
              if(zr < topography_surface(ix, jy)) then
                apply = .true.
              end if
            end if
            if(apply) then
              select case(zBoundary)
              case("periodic")
                zr = lz(1) + mod(zr - lz(0), lz(1) - lz(0))
              case("solid_wall")
                cycle
              case default
                stop "saturation_3D: unknown case zBoundary"
              end select
            elseif(zr > lz(1)) then
              select case(zBoundary)
              case("periodic")
                zr = lz(0) + mod(zr - lz(1), lz(1) - lz(0))
              case("solid_wall")
                cycle
              case default
                stop "saturation_3D: unknown case zBoundary"
              end select
            end if

            ! Skip rays propagating out of the domain.
            ! Why half-levels?
            if(topography) then
              kz = floor((levelTFC(ix, jy, zr) - lz(0)) / dz) + 1
              if(kz < 1 .or. kz > sizeZ) cycle
            else
              kz = floor((zr - lz(0)) / dz) + 1
              if(kz < 1 .or. kz > sizeZ) cycle
            end if

            wnrk = ray(iRay, ixrv, jyrv, kzrv)%k
            wnrl = ray(iRay, ixrv, jyrv, kzrv)%l
            wnrm = ray(iRay, ixrv, jyrv, kzrv)%m

            kappa = diffusion(ix, jy, kz)

            ! it can hanppen that the damping coefficient becomes
            ! negative,
            ! hence set it to zero in that case
            ray(iRay, ixrv, jyrv, kzrv)%dens = ray(iRay, ixrv, jyrv, &
                &kzrv)%dens * max(0.0, 1.0 - dt * 2.0 * kappa * (wnrk ** 2 &
                &+ wnrl ** 2 + wnrm ** 2))
          end do
        end do ! ixrv
      end do ! jyrv
    end do ! kzrv

    ! diagnostic to check the impact of the saturation parametrization on
    ! energy

    ! check m^2*B^2 again and print diagnostics

    mB2 = 0.0

    do kzrv = 1, nz
      do jyrv = 1, ny
        do ixrv = 1, nx
          if(nRay(ixrv, jyrv, kzrv) < 1) cycle

          do iRay = 1, nRay(ixrv, jyrv, kzrv)
            ! skip counting ray volumes with zero wave-action density

            if(ray(iRay, ixrv, jyrv, kzrv)%dens == 0.0) cycle

            xr = ray(iRay, ixrv, jyrv, kzrv)%x
            yr = ray(iRay, ixrv, jyrv, kzrv)%y
            zr = ray(iRay, ixrv, jyrv, kzrv)%z

            dxr = ray(iRay, ixrv, jyrv, kzrv)%dxray
            dyr = ray(iRay, ixrv, jyrv, kzrv)%dyray
            dzr = ray(iRay, ixrv, jyrv, kzrv)%dzray

            ! implement horizontal boundary conditions for ray-volume
            ! positions

            if(sizeX > 1) then
              if(xr < lx(0)) then
                select case(xBoundary)
                case("periodic")
                  xr = lx(1) + mod(xr - lx(0), lx(1) - lx(0))
                case default
                  stop "saturation_3D: unknown case xBoundary"
                end select
              elseif(xr > lx(1)) then
                select case(xBoundary)
                case("periodic")
                  xr = lx(0) + mod(xr - lx(1), lx(1) - lx(0))
                case default
                  stop "saturation_3D: unknown case xBoundary"
                end select
              end if

              ix = floor((xr - lx(0)) / dx) + 1 - ix0
            else
              ix = 1
            end if

            if(sizeY > 1) then
              if(yr < ly(0)) then
                select case(yBoundary)
                case("periodic")
                  yr = ly(1) + mod(yr - ly(0), ly(1) - ly(0))
                case default
                  stop "saturation_3D: unknown case yBoundary"
                end select
              elseif(yr > ly(1)) then
                select case(yBoundary)
                case("periodic")
                  yr = ly(0) + mod(yr - ly(1), ly(1) - ly(0))
                case default
                  stop "saturation_3D: unknown case yBoundary"
                end select
              end if

              jy = floor((yr - ly(0)) / dy) + 1 - jy0
            else
              jy = 1
            end if

            ! vertical boundary conditions:
            ! zBoundary = 'periodic': implement periodicity
            ! zBoundary = 'solid_wall': skip counting ray volumes
            !                           that have completely left the
            !                           model domain

            apply = .false.
            if(.not. topography .and. zr < lz(0)) then
              apply = .true.
            else if(topography) then
              if(zr < topography_surface(ix, jy)) then
                apply = .true.
              end if
            end if
            if(apply) then
              select case(zBoundary)
              case("periodic")
                zr = lz(1) + mod(zr - lz(0), lz(1) - lz(0))
              case("solid_wall")
                cycle
              case default
                stop "saturation_3D: unknown case zBoundary"
              end select
            elseif(zr > lz(1)) then
              select case(zBoundary)
              case("periodic")
                zr = lz(0) + mod(zr - lz(1), lz(1) - lz(0))
              case("solid_wall")
                cycle
              case default
                stop "saturation_3D: unknown case zBoundary"
              end select
            end if

            ! Skip rays propagating out of the domain.
            ! Why half-levels?
            if(topography) then
              kz = floor((levelTFC(ix, jy, zr) - lz(0)) / dz) + 1
              if(kz < 1 .or. kz > sizeZ) cycle
            else
              kz = floor((zr - lz(0)) / dz) + 1
              if(kz < 1 .or. kz > sizeZ) cycle
            end if

            ! Compute stratification.
            call stratification(zr, 1, NN_nd)

            wnrk = ray(iRay, ixrv, jyrv, kzrv)%k
            wnrl = ray(iRay, ixrv, jyrv, kzrv)%l
            wnrm = ray(iRay, ixrv, jyrv, kzrv)%m

            wnrhs = wnrk ** 2 + wnrl ** 2

            dwnrk = ray(iRay, ixrv, jyrv, kzrv)%dkray
            dwnrl = ray(iRay, ixrv, jyrv, kzrv)%dlray
            dwnrm = ray(iRay, ixrv, jyrv, kzrv)%dmray

            omir = ray(iRay, ixrv, jyrv, kzrv)%omega

            densr = ray(iRay, ixrv, jyrv, kzrv)%dens

            ! spatial extension of ray to be taken into account
            if(topography) then
              dzi = min(dzr, jac(ix, jy, kz) * dz)
              facpsp = dzi / jac(ix, jy, kz) / dz * dwnrm
            else
              dzi = min(dzr, dz)
              facpsp = dzi / dz * dwnrm
            end if

            if(sizeX > 1) then
              dxi = min(dxr, dx)
              facpsp = facpsp * dxi / dx * dwnrk
            end if

            if(sizeY > 1) then
              dyi = min(dyr, dy)
              facpsp = facpsp * dyi / dy * dwnrl
            end if

            integral1 = wnrhs * wnrm ** 2 / ((wnrhs + wnrm ** 2) * omir) &
                &* facpsp

            if(topography) then
              mB2(ix, jy, kz) = mB2(ix, jy, kz) + 2.0 * NN_nd ** 2 &
                  &/ rhoStratTFC(ix, jy, kz) * densr * integral1
            else
              mB2(ix, jy, kz) = mB2(ix, jy, kz) + 2.0 * NN_nd ** 2 &
                  &/ rhoStrat(kz) * densr * integral1
            end if
          end do
        end do ! ixrv
      end do ! jyrv
    end do ! kzrv

    if(topography) then
      do kz = 1, sizeZ
        do jy = 1, ny
          do ix = 1, nx
            call stratification(zTFC(ix, jy, kz), 1, NN_nd)
            if(mB2(ix, jy, kz) - alpha_sat ** 2 * NN_nd ** 2 > 1.d-3 &
                &* alpha_sat ** 2 * NN_nd ** 2) then
              print *, 'SATURATION VIOLATED AT ix, jy, kz =', ix, jy, kz
              print *, 'mB2(ix,jy,kz) =', mB2(ix, jy, kz)
              print *, 'alpha_sat**2 * NN_nd**2 = ', alpha_sat ** 2 * NN_nd ** 2
            endif
          end do
        end do
      end do
    else
      do kz = 1, sizeZ
        call stratification(z(kz), 1, NN_nd)
        do jy = 1, ny
          do ix = 1, nx
            if(mB2(ix, jy, kz) - alpha_sat ** 2 * NN_nd ** 2 > 1.d-3 &
                &* alpha_sat ** 2 * NN_nd ** 2) then
              print *, 'SATURATION VIOLATED AT ix, jy, kz =', ix, jy, kz
              print *, 'mB2(ix,jy,kz) =', mB2(ix, jy, kz)
              print *, 'alpha_sat**2 * NN_nd**2 = ', alpha_sat ** 2 * NN_nd ** 2
            endif
          end do
        end do
      end do
    end if

    ! Remove ray volumes with zero wave action.
    do kz = 0, nz + 1
      do jy = 0, ny + 1
        do ix = 0, nx + 1
          if(nRay(ix, jy, kz) <= 0) cycle
          nrlc = 0
          do iRay = 1, nRay(ix, jy, kz)
            if(ray(iRay, ix, jy, kz)%dens == 0.0) cycle
            nrlc = nrlc + 1
            ray(nrlc, ix, jy, kz) = ray(iRay, ix, jy, kz)
          end do
          nRay(ix, jy, kz) = nrlc
        end do
      end do
    end do

    return

  end subroutine saturation_3D

  !------------------------------------------------------------------------

  subroutine smooth_wkb_shapiro(flxwkb, nsmth, homog_dir)

    !--------------------------------------------------------------------
    !    local smoothing of WKB fluxes
    !    use Shapiro weighting up to nsmth = 4
    !    homog_dir = 111: local smoothing in all three spatial directions
    !                101: local smoothing in x and z direction
    !                 11: local smoothing in y and z direction
    !                100: global smooting in x direction
    !--------------------------------------------------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb
    integer, intent(in) :: nsmth

    integer, intent(in) :: homog_dir

    ! allocatable fields
    real, dimension(:, :, :), allocatable :: flxwkb_0, flxwkb_1

    integer :: allocstat
    integer :: i, j, k

    allocate(flxwkb_0(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "smooth_wkb_shapiro:alloc failed"

    allocate(flxwkb_1(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "smooth_wkb_shapiro:alloc failed"

    ! set the values for flxwkb_0

    flxwkb_0 = flxwkb

    ! start to do the smoothing

    select case(homog_dir)
    case(111)
      ! boundaries only filling min(n.,nb.) ghost cells, hence

      if(min(nx, nbx) < nsmth) stop 'min(nx,nbx) too small for smoothing'
      if(min(ny, nby) < nsmth) stop 'min(ny,nby) too small for smoothing'
      if(min(nz, nbz) < nsmth) stop 'min(nz,nbz) too small for smoothing'

      if(nsmth == 1) then
        ! smooth in x

        do k = - nbz, nz + nbz
          do j = - nby, ny + nby
            do i = 1, nx
              flxwkb_0(i, j, k) = (flxwkb(i - 1, j, k) + flxwkb(i + 1, j, k) &
                  &+ 2.0 * flxwkb(i, j, k)) / 4.0
            end do
          end do
        end do

        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, &
                  &k) + 2.0 * flxwkb_0(i, j, k)) / 4.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1) &
                  &+ 2.0 * flxwkb_1(i, j, k)) / 4.0
            end do
          end do
        end do
      elseif(nsmth == 2) then
        ! smooth in x

        do k = - nbz, nz + nbz
          do j = - nby, ny + nby
            do i = 1, nx
              flxwkb_0(i, j, k) = (- flxwkb(i - 2, j, k) - flxwkb(i + 2, j, k) &
                  &+ 4.0 * (flxwkb(i - 1, j, k) + flxwkb(i + 1, j, k)) + 10.0 &
                  &* flxwkb(i, j, k)) / 16.0
            end do
          end do
        end do

        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (- flxwkb_0(i, j - 2, k) - flxwkb_0(i, j &
                  &+ 2, k) + 4.0 * (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, &
                  &k)) + 10.0 * flxwkb_0(i, j, k)) / 16.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (- flxwkb_1(i, j, k - 2) - flxwkb_1(i, j, k &
                  &+ 2) + 4.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k &
                  &+ 1)) + 10.0 * flxwkb_1(i, j, k)) / 16.0
            end do
          end do
        end do
      elseif(nsmth == 3) then
        ! smooth in x

        do k = - nbz, nz + nbz
          do j = - nby, ny + nby
            do i = 1, nx
              flxwkb_0(i, j, k) = (flxwkb(i - 3, j, k) + flxwkb(i + 3, j, k) &
                  &- 6.0 * (flxwkb(i - 2, j, k) + flxwkb(i + 2, j, k)) + 15.0 &
                  &* (flxwkb(i - 1, j, k) + flxwkb(i + 1, j, k)) + 44.0 &
                  &* flxwkb(i, j, k)) / 64.0
            end do
          end do
        end do

        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (flxwkb_0(i, j - 3, k) + flxwkb_0(i, j + 3, &
                  &k) - 6.0 * (flxwkb_0(i, j - 2, k) + flxwkb_0(i, j + 2, k)) &
                  &+ 15.0 * (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, k)) &
                  &+ 44.0 * flxwkb_0(i, j, k)) / 64.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 3) + flxwkb_1(i, j, k + 3) &
                  &- 6.0 * (flxwkb_1(i, j, k - 2) + flxwkb_1(i, j, k + 2)) &
                  &+ 15.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  &+ 44.0 * flxwkb_1(i, j, k)) / 64.0
            end do
          end do
        end do
      elseif(nsmth == 4) then
        ! smooth in x

        do k = - nbz, nz + nbz
          do j = - nby, ny + nby
            do i = 1, nx
              flxwkb_0(i, j, k) = (- flxwkb(i - 4, j, k) - flxwkb(i + 4, j, k) &
                  &+ 8.0 * (flxwkb(i - 3, j, k) + flxwkb(i + 3, j, k)) - 28.0 &
                  &* (flxwkb(i - 2, j, k) + flxwkb(i + 2, j, k)) + 56.0 &
                  &* (flxwkb(i - 1, j, k) + flxwkb(i + 1, j, k)) + 186.0 &
                  &* flxwkb(i, j, k)) / 256.0
            end do
          end do
        end do

        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (- flxwkb_0(i, j - 4, k) - flxwkb_0(i, j &
                  &+ 4, k) + 8.0 * (flxwkb_0(i, j - 3, k) + flxwkb_0(i, j + 3, &
                  &k)) - 28.0 * (flxwkb_0(i, j - 2, k) + flxwkb_0(i, j + 2, &
                  &k)) + 56.0 * (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, &
                  &k)) + 186.0 * flxwkb_0(i, j, k)) / 256.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 4) + flxwkb_1(i, j, k + 4) &
                  &+ 8.0 * (flxwkb_1(i, j, k - 3) + flxwkb_1(i, j, k + 3)) &
                  &- 28.0 * (flxwkb_1(i, j, k - 2) + flxwkb_1(i, j, k + 2)) &
                  &+ 56.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  &+ 186.0 * flxwkb_1(i, j, k)) / 256.0
            end do
          end do
        end do
      else
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = sum(flxwkb_0((i - nsmth):(i + nsmth), (j &
                  &- nsmth):(j + nsmth), (k - nsmth):(k + nsmth))) / real((2 &
                  &* nsmth + 1) ** 3)
            end do
          end do
        end do
      end if
    case(101)
      ! boundaries only filling min(n.,nb.) ghost cells, hence

      if(min(nx, nbx) < nsmth) stop 'min(nx,nbx) too small for smoothing'
      if(min(nz, nbz) < nsmth) stop 'min(nz,nbz) too small for smoothing'

      if(nsmth == 1) then
        ! smooth in x

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (flxwkb_0(i - 1, j, k) + flxwkb_0(i + 1, j, &
                  &k) + 2.0 * flxwkb_0(i, j, k)) / 4.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1) &
                  &+ 2.0 * flxwkb_1(i, j, k)) / 4.0
            end do
          end do
        end do
      elseif(nsmth == 2) then
        ! smooth in x

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (- flxwkb_0(i - 2, j, k) - flxwkb_0(i + 2, &
                  &j, k) + 4.0 * (flxwkb_0(i - 1, j, k) + flxwkb_0(i + 1, j, &
                  &k)) + 10.0 * flxwkb_0(i, j, k)) / 16.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (- flxwkb_1(i, j, k - 2) - flxwkb_1(i, j, k &
                  &+ 2) + 4.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k &
                  &+ 1)) + 10.0 * flxwkb_1(i, j, k)) / 16.0
            end do
          end do
        end do
      elseif(nsmth == 3) then
        ! smooth in x

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (flxwkb_0(i - 3, j, k) + flxwkb_0(i + 3, j, &
                  &k) - 6.0 * (flxwkb_0(i - 2, j, k) + flxwkb_0(i + 2, j, k)) &
                  &+ 15.0 * (flxwkb_0(i - 1, j, k) + flxwkb_0(i + 1, j, k)) &
                  &+ 44.0 * flxwkb_0(i, j, k)) / 64.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 3) + flxwkb_1(i, j, k + 3) &
                  &- 6.0 * (flxwkb_1(i, j, k - 2) + flxwkb_1(i, j, k + 2)) &
                  &+ 15.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  &+ 44.0 * flxwkb_1(i, j, k)) / 64.0
            end do
          end do
        end do
      elseif(nsmth == 4) then
        ! smooth in x

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (- flxwkb_0(i - 4, j, k) - flxwkb_0(i + 4, &
                  &j, k) + 8.0 * (flxwkb_0(i - 3, j, k) + flxwkb_0(i + 3, j, &
                  &k)) - 28.0 * (flxwkb_0(i - 2, j, k) + flxwkb_0(i + 2, j, &
                  &k)) + 56.0 * (flxwkb_0(i - 1, j, k) + flxwkb_0(i + 1, j, &
                  &k)) + 186.0 * flxwkb_0(i, j, k)) / 256.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 4) + flxwkb_1(i, j, k + 4) &
                  &+ 8.0 * (flxwkb_1(i, j, k - 3) + flxwkb_1(i, j, k + 3)) &
                  &- 28.0 * (flxwkb_1(i, j, k - 2) + flxwkb_1(i, j, k + 2)) &
                  &+ 56.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  &+ 186.0 * flxwkb_1(i, j, k)) / 256.0
            end do
          end do
        end do
      else
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = sum(flxwkb_0((i - nsmth):(i + nsmth), j, (k &
                  &- nsmth):(k + nsmth))) / real((2 * nsmth + 1) ** 2)
            end do
          end do
        end do
      end if
    case(11)
      ! boundaries only filling min(n.,nb.) ghost cells, hence

      if(min(ny, nby) < nsmth) stop 'min(nx,nbx) too small for smoothing'
      if(min(nz, nbz) < nsmth) stop 'min(nz,nbz) too small for smoothing'

      if(nsmth == 1) then
        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, &
                  &k) + 2.0 * flxwkb_0(i, j, k)) / 4.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1) &
                  &+ 2.0 * flxwkb_1(i, j, k)) / 4.0
            end do
          end do
        end do
      elseif(nsmth == 2) then
        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (- flxwkb_0(i, j - 2, k) - flxwkb_0(i, j &
                  &+ 2, k) + 4.0 * (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, &
                  &k)) + 10.0 * flxwkb_0(i, j, k)) / 16.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (- flxwkb_1(i, j, k - 2) - flxwkb_1(i, j, k &
                  &+ 2) + 4.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k &
                  &+ 1)) + 10.0 * flxwkb_1(i, j, k)) / 16.0
            end do
          end do
        end do
      elseif(nsmth == 3) then
        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (flxwkb_0(i, j - 3, k) + flxwkb_0(i, j + 3, &
                  &k) - 6.0 * (flxwkb_0(i, j - 2, k) + flxwkb_0(i, j + 2, k)) &
                  &+ 15.0 * (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, k)) &
                  &+ 44.0 * flxwkb_0(i, j, k)) / 64.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 3) + flxwkb_1(i, j, k + 3) &
                  &- 6.0 * (flxwkb_1(i, j, k - 2) + flxwkb_1(i, j, k + 2)) &
                  &+ 15.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  &+ 44.0 * flxwkb_1(i, j, k)) / 64.0
            end do
          end do
        end do
      elseif(nsmth == 4) then
        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (- flxwkb_0(i, j - 4, k) - flxwkb_0(i, j &
                  &+ 4, k) + 8.0 * (flxwkb_0(i, j - 3, k) + flxwkb_0(i, j + 3, &
                  &k)) - 28.0 * (flxwkb_0(i, j - 2, k) + flxwkb_0(i, j + 2, &
                  &k)) + 56.0 * (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, &
                  &k)) + 186.0 * flxwkb_0(i, j, k)) / 256.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 4) + flxwkb_1(i, j, k + 4) &
                  &+ 8.0 * (flxwkb_1(i, j, k - 3) + flxwkb_1(i, j, k + 3)) &
                  &- 28.0 * (flxwkb_1(i, j, k - 2) + flxwkb_1(i, j, k + 2)) &
                  &+ 56.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  &+ 186.0 * flxwkb_1(i, j, k)) / 256.0
            end do
          end do
        end do
      else
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = sum(flxwkb_0(i, (j - nsmth):(j + nsmth), (k &
                  &- nsmth):(k + nsmth))) / real((2 * nsmth + 1) ** 2)
            end do
          end do
        end do
      end if
    case(100)
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            flxwkb(i, j, k) = sum(flxwkb_0(1:nx, j, k)) / real(nx)
          end do
        end do
      end do
    case default
      stop "unknown case homog_dir."
    end select

    call setboundary_wkb(flxwkb)

    ! deallocate local fields
    deallocate(flxwkb_0, stat = allocstat); if(allocstat /= 0) stop &
        &"smooth_wkb_shapiro:dealloc failed"
    deallocate(flxwkb_1, stat = allocstat); if(allocstat /= 0) stop &
        &"smooth_wkb_shapiro:dealloc failed"

    return

  end subroutine smooth_wkb_shapiro

  !------------------------------------------------------------------------

  subroutine smooth_wkb_box(flxwkb, nsmth, homog_dir)

    !--------------------------------------------------------------------
    !    local smoothing of WKB fluxes
    !    homog_dir = 111: local smoothing in all three spatial directions
    !                101: local smoothing in x and z direction
    !                 11: local smoothing in y and z direction
    !                100: global smooting in x direction
    !--------------------------------------------------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb
    integer, intent(in) :: nsmth

    integer, intent(in) :: homog_dir

    ! allocatable fields
    real, dimension(:, :, :), allocatable :: flxwkb_0

    integer :: allocstat
    integer :: i, j, k

    allocate(flxwkb_0(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "smooth_wkb:alloc failed"

    ! set the values for flxwkb_0

    flxwkb_0 = flxwkb

    ! start to do the smoothing

    select case(homog_dir)
    case(111)
      ! boundaries only filling min(n.,nb.) ghost cells, hence

      if(min(nx, nbx) < nsmth) stop 'min(nx,nbx) too small for smoothing'
      if(min(ny, nby) < nsmth) stop 'min(ny,nby) too small for smoothing'
      if(min(nz, nbz) < nsmth) stop 'min(nz,nbz) too small for smoothing'

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            flxwkb(i, j, k) = sum(flxwkb_0((i - nsmth):(i + nsmth), (j &
                &- nsmth):(j + nsmth), (k - nsmth):(k + nsmth))) / real((2 &
                &* nsmth + 1) ** 3)
          end do
        end do
      end do
    case(101)
      ! boundaries only filling min(n.,nb.) ghost cells, hence

      if(min(nx, nbx) < nsmth) stop 'min(nx,nbx) too small for smoothing'
      if(min(nz, nbz) < nsmth) stop 'min(nz,nbz) too small for smoothing'

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            flxwkb(i, j, k) = sum(flxwkb_0((i - nsmth):(i + nsmth), j, (k &
                &- nsmth):(k + nsmth))) / real((2 * nsmth + 1) ** 2)
          end do
        end do
      end do
    case(11)
      ! boundaries only filling min(n.,nb.) ghost cells, hence

      if(min(ny, nby) < nsmth) stop 'min(nx,nbx) too small for smoothing'
      if(min(nz, nbz) < nsmth) stop 'min(nz,nbz) too small for smoothing'

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            flxwkb(i, j, k) = sum(flxwkb_0(i, (j - nsmth):(j + nsmth), (k &
                &- nsmth):(k + nsmth))) / real((2 * nsmth + 1) ** 2)
          end do
        end do
      end do
    case(100)
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            flxwkb(i, j, k) = sum(flxwkb_0(1:nx, j, k)) / real(nx)
          end do
        end do
      end do
    case default
      stop "unknown case homog_dir."
    end select

    call setboundary_wkb(flxwkb)

    ! deallocate local fields
    deallocate(flxwkb_0, stat = allocstat); if(allocstat /= 0) stop &
        &"smooth_wkb:dealloc failed"

    return

  end subroutine smooth_wkb_box

  ! ---------------------------------------------------------------------

  subroutine setboundary_wkb(flxwkb)

    ! -------------------------------------------------------------------
    ! boundary conditions for WKB fluxes
    ! so far only periodic boundary conditions allowed in horizontal
    ! -------------------------------------------------------------------

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb

    call setboundary_hor_wkb(flxwkb)
    call setboundary_vrt_wkb(flxwkb)

    return

  end subroutine setboundary_wkb

  ! ---------------------------------------------------------------------

  subroutine setboundary_hor_wkb(flxwkb)

    ! -------------------------------------------------------------------
    ! horizontal boundary conditions for WKB fluxes
    ! so far only periodic boundary conditions allowed
    ! -------------------------------------------------------------------

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb

    select case(xBoundary)
    case("periodic")
      call setboundary_x_periodic_wkb(flxwkb)
    case default
      stop "setboundary_hor_wkb: unknown case xBoundary"
    end select

    select case(yBoundary)
    case("periodic")
      call setboundary_y_periodic_wkb(flxwkb)
    case default
      stop "setboundary_hor_wkb: unknown case xBoundary"
    end select

    return

  end subroutine setboundary_hor_wkb

  ! ----------------------------------------------------------------------

  subroutine setboundary_x_periodic_wkb(flxwkb)

    ! -------------------------------------------------------------------
    ! periodic boundary conditions in x for WKB fluxes
    ! -------------------------------------------------------------------

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: ix

    ! auxiliary fields
    real, dimension(nbx, - nby:ny + nby, - nbz:nz + nbz) :: xSliceLeft_send, &
        &xSliceRight_send
    real, dimension(nbx, - nby:ny + nby, - nbz:nz + nbz) :: xSliceLeft_recv, &
        &xSliceRight_recv

    if(idim > 1) then
      ! more than 1 cpu in x direction

      call mpi_cart_shift(comm, 0, 1, left, right, ierror)

      ! slice size
      sendcount = nbx * (ny + 2 * nby + 1) * (nz + 2 * nbz + 1)
      recvcount = sendcount

      ! read slice into contiguous array
      do ix = 1, nbx
        xSliceLeft_send(ix, :, :) = flxwkb(ix, :, :)
        xSliceRight_send(ix, :, :) = flxwkb(nx - nbx + ix, :, :)
      end do

      ! left -> right
      source = left
      dest = right
      tag = 100

      call mpi_sendrecv(xSliceRight_send(1, - nby, - nbz), sendcount, &
          &mpi_double_precision, dest, tag, xSliceLeft_recv(1, - nby, - nbz), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xSliceLeft_send(1, - nby, - nbz), sendcount, &
          &mpi_double_precision, dest, tag, xSliceRight_recv(1, - nby, - nbz), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_right, ierror)

      ! write auxiliary slice to flux field
      do ix = 1, nbx
        ! right halos
        flxwkb(nx + ix, :, :) = xSliceRight_recv(ix, :, :)

        ! left halos
        flxwkb(- nbx + ix, :, :) = xSliceLeft_recv(ix, :, :)
      end do
    else
      ! only 1 cpu in x direction

      do ix = 1, min(nx, nbx)
        flxwkb(nx + ix, :, :) = flxwkb(ix, :, :)
        flxwkb(- ix + 1, :, :) = flxwkb(nx - ix + 1, :, :)
      end do
    end if

  end subroutine setboundary_x_periodic_wkb

  ! ----------------------------------------------------------------------

  subroutine setboundary_y_periodic_wkb(flxwkb)

    ! -------------------------------------------------------------------
    ! periodic boundary conditions in y for WKB fluxes
    ! -------------------------------------------------------------------

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: jy

    ! auxiliary fields
    real, dimension(- nbx:nx + nbx, nby, - nbz:nz + nbz) :: ySliceBack_send, &
        &ySliceForw_send
    real, dimension(- nbx:nx + nbx, nby, - nbz:nz + nbz) :: ySliceBack_recv, &
        &ySliceForw_recv

    if(jdim > 1) then
      ! more than 1 cpu in y direction

      call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

      ! slice size
      sendcount = nby * (nx + 2 * nbx + 1) * (nz + 2 * nbz + 1)
      recvcount = sendcount

      ! read slice into contiguous array
      do jy = 1, nby
        ySliceBack_send(:, jy, :) = flxwkb(:, jy, :)
        ySliceForw_send(:, jy, :) = flxwkb(:, ny - nby + jy, :)
      end do

      ! back -> forw
      source = back
      dest = forw
      tag = 100

      call mpi_sendrecv(ySliceForw_send(- nbx, 1, - nbz), sendcount, &
          &mpi_double_precision, dest, tag, ySliceBack_recv(- nbx, 1, - nbz), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(ySliceBack_send(- nbx, 1, - nbz), sendcount, &
          &mpi_double_precision, dest, tag, ySliceForw_recv(- nbx, 1, - nbz), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_forw, ierror)

      ! write auxiliary slice to var field
      do jy = 1, nby
        ! right halos
        flxwkb(:, ny + jy, :) = ySliceForw_recv(:, jy, :)

        ! left halos
        flxwkb(:, - nby + jy, :) = ySliceBack_recv(:, jy, :)
      end do
    else
      ! only 1 cpu in y direction

      do jy = 1, min(ny, nby)
        flxwkb(:, ny + jy, :) = flxwkb(:, jy, :)
        flxwkb(:, - jy + 1, :) = flxwkb(:, ny - jy + 1, :)
      end do
    end if

  end subroutine setboundary_y_periodic_wkb

  ! ---------------------------------------------------------------------

  subroutine setboundary_vrt_wkb(flxwkb)

    ! -------------------------------------------------------------------
    ! vertical boundary conditions for WKB fluxes
    ! -------------------------------------------------------------------

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb

    select case(zBoundary)
    case("periodic")
      call setboundary_z_periodic_wkb(flxwkb)
    case("solid_wall")
      call setboundary_z_solidwall_wkb(flxwkb)
    case default
      stop "setboundary_hor_wkb: unknown case xBoundary"
    end select

    return

  end subroutine setboundary_vrt_wkb

  ! ----------------------------------------------------------------------

  subroutine setboundary_z_periodic_wkb(flxwkb)

    ! -------------------------------------------------------------------
    ! periodic boundary conditions in z for WKB fluxes
    ! -------------------------------------------------------------------

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb

    integer :: k

    do k = 1, min(nz, nbz)
      flxwkb(:, :, nz + k) = flxwkb(:, :, k)
      flxwkb(:, :, - k + 1) = flxwkb(:, :, nz - k + 1)
    end do

  end subroutine setboundary_z_periodic_wkb

  ! ----------------------------------------------------------------------

  subroutine setboundary_z_solidwall_wkb(flxwkb)

    ! -------------------------------------------------------------------
    ! solid-wall boundary conditions in z for WKB fluxes
    !
    ! at lower boundary constant fluxes assumed, for consistency with
    ! launching from the bottom
    ! -------------------------------------------------------------------

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb

    integer :: k

    do k = 1, min(nz, nbz)
      flxwkb(:, :, nz + k) = flxwkb(:, :, nz - k + 1)
      flxwkb(:, :, - k + 1) = flxwkb(:, :, 1)
      !      flxwkb(:,:,-k+1) = flxwkb(:,:,k)
    end do

  end subroutine setboundary_z_solidwall_wkb

  ! ---------------------------------------------------------------------

  subroutine setboundary_frc_wkb(frcwkb)

    ! -------------------------------------------------------------------
    ! boundary conditions for WKB force
    ! so far only periodic boundary conditions allowed in horizontal
    ! -------------------------------------------------------------------

    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(inout) :: frcwkb

    call setboundary_frc_hor_wkb(frcwkb)
    call setboundary_frc_vrt_wkb(frcwkb)

    return

  end subroutine setboundary_frc_wkb

  ! ---------------------------------------------------------------------

  subroutine setboundary_frc_hor_wkb(frcwkb)

    ! -------------------------------------------------------------------
    ! horizontal boundary conditions for WKB forces
    ! so far only periodic boundary conditions allowed
    ! -------------------------------------------------------------------

    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(inout) :: frcwkb

    select case(xBoundary)
    case("periodic")
      call setboundary_frc_x_periodic_wkb(frcwkb)
    case default
      stop "setboundary_frc_hor_wkb: unknown case xBoundary"
    end select

    select case(yBoundary)
    case("periodic")
      call setboundary_frc_y_periodic_wkb(frcwkb)
    case default
      stop "setboundary_frc_hor_wkb: unknown case xBoundary"
    end select

    return

  end subroutine setboundary_frc_hor_wkb

  ! ----------------------------------------------------------------------

  subroutine setboundary_frc_x_periodic_wkb(frcwkb)

    ! -------------------------------------------------------------------
    ! periodic boundary conditions in x for WKB fluxes
    ! -------------------------------------------------------------------

    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(inout) :: frcwkb

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: ix

    ! auxiliary fields
    real, dimension(1, 0:ny + 1, 0:nz + 1) :: xSliceLeft_send, xSliceRight_send
    real, dimension(1, 0:ny + 1, 0:nz + 1) :: xSliceLeft_recv, xSliceRight_recv

    if(idim > 1) then
      ! more than 1 cpu in x direction

      call mpi_cart_shift(comm, 0, 1, left, right, ierror)

      ! slice size
      sendcount = (ny + 2) * (nz + 2)
      recvcount = sendcount

      ! read slice into contiguous array
      xSliceLeft_send(1, :, :) = frcwkb(1, :, :)
      xSliceRight_send(1, :, :) = frcwkb(nx, :, :)

      ! left -> right
      source = left
      dest = right
      tag = 100

      call mpi_sendrecv(xSliceRight_send(1, 0, 0), sendcount, &
          &mpi_double_precision, dest, tag, xSliceLeft_recv(1, 0, 0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xSliceLeft_send(1, 0, 0), sendcount, &
          &mpi_double_precision, dest, tag, xSliceRight_recv(1, 0, 0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_right, ierror)

      ! write auxiliary slice to forcing field
      ! right halos
      frcwkb(nx + 1, :, :) = xSliceRight_recv(1, :, :)

      ! left halos
      frcwkb(0, :, :) = xSliceLeft_recv(1, :, :)
    else
      ! only 1 cpu in x direction

      frcwkb(nx + 1, :, :) = frcwkb(1, :, :)
      frcwkb(0, :, :) = frcwkb(nx, :, :)
    end if

  end subroutine setboundary_frc_x_periodic_wkb

  ! ----------------------------------------------------------------------

  subroutine setboundary_frc_y_periodic_wkb(frcwkb)

    ! -------------------------------------------------------------------
    ! periodic boundary conditions in y for WKB forces
    ! -------------------------------------------------------------------

    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(inout) :: frcwkb

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! auxiliary fields
    real, dimension(0:nx + 1, 1, 0:nz + 1) :: ySliceBack_send, ySliceForw_send
    real, dimension(0:nx + 1, 1, 0:nz + 1) :: ySliceBack_recv, ySliceForw_recv

    if(jdim > 1) then
      ! more than 1 cpu in y direction

      call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

      ! slice size
      sendcount = (nx + 2) * (nz + 2)
      recvcount = sendcount

      ! read slice into contiguous array
      ySliceBack_send(:, 1, :) = frcwkb(:, 1, :)
      ySliceForw_send(:, 1, :) = frcwkb(:, ny, :)

      ! back -> forw
      source = back
      dest = forw
      tag = 100

      call mpi_sendrecv(ySliceForw_send(0, 1, 0), sendcount, &
          &mpi_double_precision, dest, tag, ySliceBack_recv(0, 1, 0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(ySliceBack_send(0, 1, 0), sendcount, &
          &mpi_double_precision, dest, tag, ySliceForw_recv(0, 1, 0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_forw, ierror)

      ! write auxiliary slice to var field
      ! right halos
      frcwkb(:, ny + 1, :) = ySliceForw_recv(:, 1, :)

      ! left halos
      frcwkb(:, 0, :) = ySliceBack_recv(:, 1, :)
    else
      frcwkb(:, ny + 1, :) = frcwkb(:, 1, :)
      frcwkb(:, 0, :) = frcwkb(:, ny, :)
    end if

  end subroutine setboundary_frc_y_periodic_wkb

  ! ---------------------------------------------------------------------

  subroutine setboundary_frc_vrt_wkb(frcwkb)

    ! -------------------------------------------------------------------
    ! vertical boundary conditions for WKB forces
    ! -------------------------------------------------------------------

    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(inout) :: frcwkb

    select case(zBoundary)
    case("periodic")
      call setboundary_frc_z_periodic_wkb(frcwkb)
    case("solid_wall")
      call setboundary_frc_z_solidwall_wkb(frcwkb)
    case default
      stop "setboundary_frc_vrt_wkb: unknown case zBoundary"
    end select

    return

  end subroutine setboundary_frc_vrt_wkb

  ! ----------------------------------------------------------------------

  subroutine setboundary_frc_z_periodic_wkb(frcwkb)

    ! -------------------------------------------------------------------
    ! periodic boundary conditions in z for WKB fluxes
    ! -------------------------------------------------------------------

    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(inout) :: frcwkb

    frcwkb(:, :, nz + 1) = frcwkb(:, :, 1)
    frcwkb(:, :, 0) = frcwkb(:, :, nz)

  end subroutine setboundary_frc_z_periodic_wkb

  ! ----------------------------------------------------------------------

  subroutine setboundary_frc_z_solidwall_wkb(frcwkb)

    ! -------------------------------------------------------------------
    ! solid-wall boundary conditions in z for WKB forces
    ! -------------------------------------------------------------------

    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(inout) :: frcwkb

    frcwkb(:, :, nz + 1) = - frcwkb(:, :, nz)
    frcwkb(:, :, 0) = - frcwkb(:, :, 1)

  end subroutine setboundary_frc_z_solidwall_wkb

  ! ----------------------------------------------------------------------

  subroutine setboundary_rayvol_x(ray)

    ! -------------------------------------------------------------------
    ! (periodic) boundary conditions for ray volumes in x direction
    ! phase-space volume not transferred but calculated from the
    ! transferred edge lengths
    ! -------------------------------------------------------------------

    implicit none

    ! argument list
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(inout) :: ray

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: ix, irprop, nrmxrl, nrmxll, nrmaxl, nrmaxr, kz, jy, iRay
    integer :: ix0

    real :: xr, xrt

    ! auxiliary fields
    real, dimension(:, :, :, :), allocatable :: xSliceLeft_send, &
        &xSliceRight_send
    real, dimension(:, :, :, :), allocatable :: xSliceLeft_recv, &
        &xSliceRight_recv

    if(xBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_rayvol_x must be  &
          &periodic!'
    end if

    ix0 = is + nbx - 1

    ! boundary conditions for number of ray volumes per cell

    call setboundary_nray_x

    if(idim > 1) then
      ! more than 1 cpu in x direction

      call mpi_cart_shift(comm, 0, 1, left, right, ierror)

      ! find global maximum of number of ray volumes per cell
      ! (at left and right edges of the cpu domains)

      nrmxll = maxval(nRay(1, :, :))
      nrmxrl = maxval(nRay(nx, :, :))

      call mpi_reduce(nrmxll, nrmaxl, 1, mpi_integer, mpi_max, root, comm, &
          &ierror)
      call mpi_bcast(nrmaxl, 1, mpi_integer, root, comm, ierror)

      call mpi_reduce(nrmxrl, nrmaxr, 1, mpi_integer, mpi_max, root, comm, &
          &ierror)
      call mpi_bcast(nrmaxr, 1, mpi_integer, root, comm, ierror)

      allocate(xSliceLeft_send(nrmaxl, 1, 0:ny + 1, 0:nz + 1))
      allocate(xSliceRight_send(nrmaxr, 1, 0:ny + 1, 0:nz + 1))

      allocate(xSliceLeft_recv(nrmaxr, 1, 0:ny + 1, 0:nz + 1))
      allocate(xSliceRight_recv(nrmaxl, 1, 0:ny + 1, 0:nz + 1))

      do irprop = 1, 18
        ! read slice into contiguous array

        do kz = 0, nz + 1
          do jy = 0, ny + 1
            if(nRay(1, jy, kz) > 0) then
              do iRay = 1, nRay(1, jy, kz)
                if(irprop == 1) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%x
                elseif(irprop == 2) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%y
                elseif(irprop == 3) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%z
                elseif(irprop == 4) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%k
                elseif(irprop == 5) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%l
                elseif(irprop == 6) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%m
                elseif(irprop == 7) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%omega
                elseif(irprop == 8) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%dkray
                elseif(irprop == 9) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%dlray
                elseif(irprop == 10) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%dmray
                elseif(irprop == 11) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%dxray
                elseif(irprop == 12) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%dyray
                elseif(irprop == 13) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%dzray
                elseif(irprop == 14) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%dens
                elseif(irprop == 15) then
                  xSliceLeft_send(iRay, 1, jy, kz) = ray(iRay, 1, jy, kz)%dphi
                end if
              end do
            end if
          end do
        end do

        do kz = 0, nz + 1
          do jy = 0, ny + 1
            if(nRay(nx, jy, kz) > 0) then
              do iRay = 1, nRay(nx, jy, kz)
                if(irprop == 1) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, kz)%x
                elseif(irprop == 2) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, kz)%y
                elseif(irprop == 3) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, kz)%z
                elseif(irprop == 4) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, kz)%k
                elseif(irprop == 5) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, kz)%l
                elseif(irprop == 6) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, kz)%m
                elseif(irprop == 7) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, &
                      &kz)%omega
                elseif(irprop == 8) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, &
                      &kz)%dkray
                elseif(irprop == 9) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, &
                      &kz)%dlray
                elseif(irprop == 10) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, &
                      &kz)%dmray
                elseif(irprop == 11) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, &
                      &kz)%dxray
                elseif(irprop == 12) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, &
                      &kz)%dyray
                elseif(irprop == 13) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, &
                      &kz)%dzray
                elseif(irprop == 14) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, kz)%dens
                elseif(irprop == 15) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, kz)%dphi
                end if
              end do
            end if
          end do
        end do

        if(irprop < 16) then
          ! left -> right
          sendcount = nrmaxr * (ny + 2) * (nz + 2)
          recvcount = sendcount
          source = left
          dest = right
          tag = 100

          call mpi_sendrecv(xSliceRight_send(1, 1, 0, 0), sendcount, &
              &mpi_double_precision, dest, tag, xSliceLeft_recv(1, 1, 0, 0), &
              &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              &sts_left, ierror)

          ! right -> left
          sendcount = nrmaxl * (ny + 2) * (nz + 2)
          recvcount = sendcount
          source = right
          dest = left
          tag = 100

          call mpi_sendrecv(xSliceLeft_send(1, 1, 0, 0), sendcount, &
              &mpi_double_precision, dest, tag, xSliceRight_recv(1, 1, 0, 0), &
              &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              &sts_right, ierror)
        end if

        ! write auxiliary slice to ray field

        ! right halos

        do kz = 0, nz + 1
          do jy = 0, ny + 1
            if(nRay(nx + 1, jy, kz) > 0) then
              do iRay = 1, nRay(nx + 1, jy, kz)
                if(irprop == 1) then
                  ray(iRay, nx + 1, jy, kz)%x = xSliceRight_recv(iRay, 1, jy, &
                      &kz)
                elseif(irprop == 2) then
                  ray(iRay, nx + 1, jy, kz)%y = xSliceRight_recv(iRay, 1, jy, &
                      &kz)
                elseif(irprop == 3) then
                  ray(iRay, nx + 1, jy, kz)%z = xSliceRight_recv(iRay, 1, jy, &
                      &kz)
                elseif(irprop == 4) then
                  ray(iRay, nx + 1, jy, kz)%k = xSliceRight_recv(iRay, 1, jy, &
                      &kz)
                elseif(irprop == 5) then
                  ray(iRay, nx + 1, jy, kz)%l = xSliceRight_recv(iRay, 1, jy, &
                      &kz)
                elseif(irprop == 6) then
                  ray(iRay, nx + 1, jy, kz)%m = xSliceRight_recv(iRay, 1, jy, &
                      &kz)
                elseif(irprop == 7) then
                  ray(iRay, nx + 1, jy, kz)%omega = xSliceRight_recv(iRay, 1, &
                      &jy, kz)
                elseif(irprop == 8) then
                  ray(iRay, nx + 1, jy, kz)%dkray = xSliceRight_recv(iRay, 1, &
                      &jy, kz)
                elseif(irprop == 9) then
                  ray(iRay, nx + 1, jy, kz)%dlray = xSliceRight_recv(iRay, 1, &
                      &jy, kz)
                elseif(irprop == 10) then
                  ray(iRay, nx + 1, jy, kz)%dmray = xSliceRight_recv(iRay, 1, &
                      &jy, kz)
                elseif(irprop == 11) then
                  ray(iRay, nx + 1, jy, kz)%dxray = xSliceRight_recv(iRay, 1, &
                      &jy, kz)
                elseif(irprop == 12) then
                  ray(iRay, nx + 1, jy, kz)%dyray = xSliceRight_recv(iRay, 1, &
                      &jy, kz)
                elseif(irprop == 13) then
                  ray(iRay, nx + 1, jy, kz)%dzray = xSliceRight_recv(iRay, 1, &
                      &jy, kz)
                elseif(irprop == 14) then
                  ray(iRay, nx + 1, jy, kz)%dens = xSliceRight_recv(iRay, 1, &
                      &jy, kz)
                elseif(irprop == 15) then
                  ray(iRay, nx + 1, jy, kz)%dphi = xSliceRight_recv(iRay, 1, &
                      &jy, kz)
                end if
              end do
            end if

            if(irprop == 1) then
              ! rightmost cpu:
              ! for ray volumes in last cell in x direction and in the
              ! ghost cell after that (that might have been handed
              ! over from the leftmost cpu) make sure that xr is
              ! adjusted accordingly

              if(is + nbx + nx - 1 == sizeX) then
                do ix = nx, nx + 1
                  if(nRay(ix, jy, kz) > 0) then
                    do iRay = 1, nRay(ix, jy, kz)
                      xr = ray(iRay, ix, jy, kz)%x

                      xrt = xr + lx(1) - lx(0)

                      if(abs(xrt - x(ix + ix0)) < abs(xr - x(ix + ix0))) then
                        xr = xrt
                      end if

                      ray(iRay, ix, jy, kz)%x = xr
                    end do
                  end if
                end do
              end if
            elseif(irprop > 15) then
              if(nRay(nx + 1, jy, kz) > 0) then
                do iRay = 1, nRay(nx + 1, jy, kz)
                  if(irprop == 16) then
                    ray(iRay, nx + 1, jy, kz)%area_xk = ray(iRay, nx + 1, jy, &
                        &kz)%dxray * ray(iRay, nx + 1, jy, kz)%dkray
                  else if(irprop == 17) then
                    ray(iRay, nx + 1, jy, kz)%area_yl = ray(iRay, nx + 1, jy, &
                        &kz)%dyray * ray(iRay, nx + 1, jy, kz)%dlray
                  else if(irprop == 18) then
                    ray(iRay, nx + 1, jy, kz)%area_zm = ray(iRay, nx + 1, jy, &
                        &kz)%dzray * ray(iRay, nx + 1, jy, kz)%dmray
                  end if
                end do
              end if
            end if
          end do
        end do

        ! left halos

        do kz = 0, nz + 1
          do jy = 0, ny + 1
            if(nRay(0, jy, kz) > 0) then
              do iRay = 1, nRay(0, jy, kz)
                if(irprop == 1) then
                  ray(iRay, 0, jy, kz)%x = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 2) then
                  ray(iRay, 0, jy, kz)%y = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 3) then
                  ray(iRay, 0, jy, kz)%z = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 4) then
                  ray(iRay, 0, jy, kz)%k = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 5) then
                  ray(iRay, 0, jy, kz)%l = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 6) then
                  ray(iRay, 0, jy, kz)%m = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 7) then
                  ray(iRay, 0, jy, kz)%omega = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 8) then
                  ray(iRay, 0, jy, kz)%dkray = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 9) then
                  ray(iRay, 0, jy, kz)%dlray = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 10) then
                  ray(iRay, 0, jy, kz)%dmray = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 11) then
                  ray(iRay, 0, jy, kz)%dxray = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 12) then
                  ray(iRay, 0, jy, kz)%dyray = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 13) then
                  ray(iRay, 0, jy, kz)%dzray = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 14) then
                  ray(iRay, 0, jy, kz)%dens = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 15) then
                  ray(iRay, 0, jy, kz)%dphi = xSliceLeft_recv(iRay, 1, jy, kz)
                end if
              end do
            end if

            if(irprop == 1) then
              ! leftmost cpu:
              ! for ray volumes in first cell in x direction and in
              ! the ghost cell before that (that might have been
              ! handed over from the rightmost cpu) make sure that
              ! xr is adjusted accordingly

              if(is + nbx == 1) then
                do ix = 0, 1
                  if(nRay(ix, jy, kz) > 0) then
                    do iRay = 1, nRay(ix, jy, kz)
                      xr = ray(iRay, ix, jy, kz)%x

                      xrt = xr - lx(1) + lx(0)

                      if(abs(xrt - x(ix + ix0)) < abs(xr - x(ix + ix0))) then
                        xr = xrt
                      end if

                      ray(iRay, ix, jy, kz)%x = xr
                    end do
                  end if
                end do
              end if
            elseif(irprop > 15) then
              if(nRay(0, jy, kz) > 0) then
                do iRay = 1, nRay(0, jy, kz)
                  if(irprop == 16) then
                    ray(iRay, 0, jy, kz)%area_xk = ray(iRay, 0, jy, kz)%dxray &
                        &* ray(iRay, 0, jy, kz)%dkray
                  else if(irprop == 17) then
                    ray(iRay, 0, jy, kz)%area_yl = ray(iRay, 0, jy, kz)%dyray &
                        &* ray(iRay, 0, jy, kz)%dlray
                  else if(irprop == 18) then
                    ray(iRay, 0, jy, kz)%area_zm = ray(iRay, 0, jy, kz)%dzray &
                        &* ray(iRay, 0, jy, kz)%dmray
                  end if
                end do
              end if
            end if
          end do
        end do
      end do

      deallocate(xSliceLeft_send)
      deallocate(xSliceRight_send)

      deallocate(xSliceLeft_recv)
      deallocate(xSliceRight_recv)
    else
      ! only 1 cpu in x direction

      do kz = 0, nz + 1
        do jy = 0, ny + 1
          if(nRay(0, jy, kz) > 0) then
            do iRay = 1, nRay(0, jy, kz)
              ray(iRay, 0, jy, kz) = ray(iRay, nx, jy, kz)
            end do
          end if

          do ix = 0, 1
            if(nRay(ix, jy, kz) > 0) then
              do iRay = 1, nRay(ix, jy, kz)
                xr = ray(iRay, ix, jy, kz)%x

                xrt = xr - lx(1) + lx(0)

                if(abs(xrt - x(ix + ix0)) < abs(xr - x(ix + ix0))) then
                  xr = xrt
                end if

                ray(iRay, ix, jy, kz)%x = xr
              end do
            end if
          end do
        end do
      end do

      do kz = 0, nz + 1
        do jy = 0, ny + 1
          if(nRay(nx + 1, jy, kz) > 0) then
            do iRay = 1, nRay(nx + 1, jy, kz)
              ray(iRay, nx + 1, jy, kz) = ray(iRay, 1, jy, kz)
            end do
          end if

          do ix = nx, nx + 1
            if(nRay(ix, jy, kz) > 0) then
              do iRay = 1, nRay(ix, jy, kz)
                xr = ray(iRay, ix, jy, kz)%x

                xrt = xr + lx(1) - lx(0)

                if(abs(xrt - x(ix + ix0)) < abs(xr - x(ix + ix0))) then
                  xr = xrt
                end if

                ray(iRay, ix, jy, kz)%x = xr
              end do
            end if
          end do
        end do
      end do
    end if

  end subroutine setboundary_rayvol_x

  ! ----------------------------------------------------------------------

  subroutine setboundary_irshift_x(irsl, irsr)

    ! -------------------------------------------------------------------
    ! (periodic) boundary conditions in x direction for indices of
    ! ray volumes shifted to the left or right
    ! -------------------------------------------------------------------

    implicit none

    ! argument list
    integer, dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), intent(inout) &
        &:: irsl, irsr

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: ix, irprop, nrmxrl, nrmxll, nrmaxl, nrmaxr
    integer :: jy, kz, iRay

    ! auxiliary fields
    integer, dimension(:, :, :, :), allocatable :: xSliceLeft_send, &
        &xSliceRight_send
    integer, dimension(:, :, :, :), allocatable :: xSliceLeft_recv, &
        &xSliceRight_recv

    if(xBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_irshift_x must be  &
          &periodic!'
    end if

    ! boundary conditions for number of ray volumes per cell

    call setboundary_nray_x

    if(idim > 1) then
      ! more than 1 cpu in x direction

      call mpi_cart_shift(comm, 0, 1, left, right, ierror)

      ! find global maximum of number of ray volumes per cell
      ! (at left and right edges of the cpu domains)

      nrmxll = maxval(nRay(1, :, :))
      nrmxrl = maxval(nRay(nx, :, :))

      call mpi_reduce(nrmxll, nrmaxl, 1, mpi_integer, mpi_max, root, comm, &
          &ierror)
      call mpi_bcast(nrmaxl, 1, mpi_integer, root, comm, ierror)

      call mpi_reduce(nrmxrl, nrmaxr, 1, mpi_integer, mpi_max, root, comm, &
          &ierror)
      call mpi_bcast(nrmaxr, 1, mpi_integer, root, comm, ierror)

      allocate(xSliceLeft_send(nrmaxl, 1, 0:ny + 1, 0:nz + 1))
      allocate(xSliceRight_send(nrmaxr, 1, 0:ny + 1, 0:nz + 1))

      allocate(xSliceLeft_recv(nrmaxr, 1, 0:ny + 1, 0:nz + 1))
      allocate(xSliceRight_recv(nrmaxl, 1, 0:ny + 1, 0:nz + 1))

      do irprop = 1, 2
        ! read slice into contiguous array

        do kz = 0, nz + 1
          do jy = 0, ny + 1
            if(nRay(1, jy, kz) > 0) then
              do iRay = 1, nRay(1, jy, kz)
                if(irprop == 1) then
                  xSliceLeft_send(iRay, 1, jy, kz) = irsl(iRay, 1, jy, kz)
                elseif(irprop == 2) then
                  xSliceLeft_send(iRay, 1, jy, kz) = irsr(iRay, 1, jy, kz)
                end if
              end do
            end if
          end do
        end do

        do kz = 0, nz + 1
          do jy = 0, ny + 1
            if(nRay(nx, jy, kz) > 0) then
              do iRay = 1, nRay(nx, jy, kz)
                if(irprop == 1) then
                  xSliceRight_send(iRay, 1, jy, kz) = irsl(iRay, nx, jy, kz)
                elseif(irprop == 2) then
                  xSliceRight_send(iRay, 1, jy, kz) = irsr(iRay, nx, jy, kz)
                end if
              end do
            end if
          end do
        end do

        ! left -> right
        sendcount = nrmaxr * (ny + 2) * (nz + 2)
        recvcount = sendcount
        source = left
        dest = right
        tag = 100

        call mpi_sendrecv(xSliceRight_send(1, 1, 0, 0), sendcount, &
            &mpi_integer, dest, tag, xSliceLeft_recv(1, 1, 0, 0), recvcount, &
            &mpi_integer, source, mpi_any_tag, comm, sts_left, ierror)

        ! right -> left
        sendcount = nrmaxl * (ny + 2) * (nz + 2)
        recvcount = sendcount
        source = right
        dest = left
        tag = 100

        call mpi_sendrecv(xSliceLeft_send(1, 1, 0, 0), sendcount, mpi_integer, &
            &dest, tag, xSliceRight_recv(1, 1, 0, 0), recvcount, mpi_integer, &
            &source, mpi_any_tag, comm, sts_right, ierror)

        ! write auxiliary slice to ray field

        ! right halos
        do kz = 0, nz + 1
          do jy = 0, ny + 1
            if(nRay(nx + 1, jy, kz) > 0) then
              do iRay = 1, nRay(nx + 1, jy, kz)
                if(irprop == 1) then
                  irsl(iRay, nx + 1, jy, kz) = xSliceRight_recv(iRay, 1, jy, kz)
                elseif(irprop == 2) then
                  irsr(iRay, nx + 1, jy, kz) = xSliceRight_recv(iRay, 1, jy, kz)
                end if
              end do
            end if
          end do
        end do

        ! left halos

        do kz = 0, nz + 1
          do jy = 0, ny + 1
            if(nRay(0, jy, kz) > 0) then
              do iRay = 1, nRay(0, jy, kz)
                if(irprop == 1) then
                  irsl(iRay, 0, jy, kz) = xSliceLeft_recv(iRay, 1, jy, kz)
                elseif(irprop == 2) then
                  irsr(iRay, 0, jy, kz) = xSliceLeft_recv(iRay, 1, jy, kz)
                end if
              end do
            end if
          end do
        end do
      end do

      deallocate(xSliceLeft_send)
      deallocate(xSliceRight_send)

      deallocate(xSliceLeft_recv)
      deallocate(xSliceRight_recv)
    else
      ! only 1 cpu in x direction

      irsl(:, 0, :, :) = irsl(:, nx, :, :)
      irsl(:, nx + 1, :, :) = irsl(:, 1, :, :)

      irsr(:, 0, :, :) = irsr(:, nx, :, :)
      irsr(:, nx + 1, :, :) = irsr(:, 1, :, :)
    end if

  end subroutine setboundary_irshift_x

  ! ----------------------------------------------------------------------

  subroutine setboundary_nshift_x(nshl, nshr)

    ! -------------------------------------------------------------------
    ! (periodic) boundary conditions in x direction for number of
    ! ray volumes shifted to the left or right
    ! -------------------------------------------------------------------

    implicit none

    ! argument list
    integer, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(inout) :: nshl, &
        &nshr

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: ix, irprop

    ! auxiliary fields
    integer, dimension(1, 0:ny + 1, 0:nz + 1) :: xSliceLeft_send, &
        &xSliceRight_send
    integer, dimension(1, 0:ny + 1, 0:nz + 1) :: xSliceLeft_recv, &
        &xSliceRight_recv

    if(xBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_nshift_x must be  &
          &periodic!'
    end if

    if(idim > 1) then
      ! more than 1 cpu in x direction

      call mpi_cart_shift(comm, 0, 1, left, right, ierror)

      ! slice size
      sendcount = (ny + 2) * (nz + 2)
      recvcount = sendcount

      do irprop = 1, 2
        ! read slice into contiguous array

        if(irprop == 1) then
          xSliceLeft_send(1, :, :) = nshl(1, :, :)
          xSliceRight_send(1, :, :) = nshl(nx, :, :)
        elseif(irprop == 2) then
          xSliceLeft_send(1, :, :) = nshr(1, :, :)
          xSliceRight_send(1, :, :) = nshr(nx, :, :)
        end if

        ! left -> right
        source = left
        dest = right
        tag = 100

        call mpi_sendrecv(xSliceRight_send(1, 0, 0), sendcount, mpi_integer, &
            &dest, tag, xSliceLeft_recv(1, 0, 0), recvcount, mpi_integer, &
            &source, mpi_any_tag, comm, sts_left, ierror)

        ! right -> left
        source = right
        dest = left
        tag = 100

        call mpi_sendrecv(xSliceLeft_send(1, 0, 0), sendcount, mpi_integer, &
            &dest, tag, xSliceRight_recv(1, 0, 0), recvcount, mpi_integer, &
            &source, mpi_any_tag, comm, sts_right, ierror)

        ! write auxiliary slice to ray field

        if(irprop == 1) then
          ! right halos
          nshl(nx + 1, :, :) = xSliceRight_recv(1, :, :)

          ! left halos
          nshl(0, :, :) = xSliceLeft_recv(1, :, :)
        elseif(irprop == 2) then
          nshr(nx + 1, :, :) = xSliceRight_recv(1, :, :)
          nshr(0, :, :) = xSliceLeft_recv(1, :, :)
        end if
      end do
    else
      ! only 1 cpu in x direction

      nshl(0, :, :) = nshl(nx, :, :)
      nshl(nx + 1, :, :) = nshl(1, :, :)

      nshr(0, :, :) = nshr(nx, :, :)
      nshr(nx + 1, :, :) = nshr(1, :, :)
    end if

  end subroutine setboundary_nshift_x

  ! ----------------------------------------------------------------------

  subroutine setboundary_rayvol_y(ray)

    ! -------------------------------------------------------------------
    ! (periodic) boundary conditions for ray volumes in y direction
    ! phase-space volume not transferred but calculated from the
    ! transferred edge lengths
    ! -------------------------------------------------------------------

    implicit none

    ! argument list
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(inout) :: ray

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: jy, irprop, ix, kz, nrmxbl, nrmxfl, nrmaxb, nrmaxf
    integer :: jy0

    real :: yr, yrt

    ! auxiliary fields
    real, dimension(:, :, :, :), allocatable :: ySliceBack_send, ySliceForw_send
    real, dimension(:, :, :, :), allocatable :: ySliceBack_recv, ySliceForw_recv

    if(yBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_rayvol_y must be  &
          &periodic!'
    end if

    jy0 = js + nby - 1

    ! boundary conditions for number of ray volumes per cell

    call setboundary_nray_y

    if(jdim > 1) then
      ! more than 1 cpu in y direction

      call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

      ! find global maximum of number of ray volumes per cell
      ! (at forward and backward edges of the cpu domains)

      nrmxbl = maxval(nRay(:, 1, :))
      nrmxfl = maxval(nRay(:, ny, :))

      call mpi_reduce(nrmxbl, nrmaxb, 1, mpi_integer, mpi_max, root, comm, &
          &ierror)
      call mpi_bcast(nrmaxb, 1, mpi_integer, root, comm, ierror)

      call mpi_reduce(nrmxfl, nrmaxf, 1, mpi_integer, mpi_max, root, comm, &
          &ierror)
      call mpi_bcast(nrmaxf, 1, mpi_integer, root, comm, ierror)

      allocate(ySliceBack_send(nrmaxb, 0:nx + 1, 1, 0:nz + 1))
      allocate(ySliceForw_send(nrmaxf, 0:nx + 1, 1, 0:nz + 1))

      allocate(ySliceBack_recv(nrmaxf, 0:nx + 1, 1, 0:nz + 1))
      allocate(ySliceForw_recv(nrmaxb, 0:nx + 1, 1, 0:nz + 1))

      do irprop = 1, 18
        ! read slice into contiguous array

        do kz = 0, nz + 1
          do ix = 0, nx + 1
            if(nRay(ix, 1, kz) > 0) then
              do iRay = 1, nRay(ix, 1, kz)
                if(irprop == 1) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%x
                elseif(irprop == 2) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%y
                elseif(irprop == 3) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%z
                elseif(irprop == 4) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%k
                elseif(irprop == 5) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%l
                elseif(irprop == 6) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%m
                elseif(irprop == 7) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%omega
                elseif(irprop == 8) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%dkray
                elseif(irprop == 9) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%dlray
                elseif(irprop == 10) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%dmray
                elseif(irprop == 11) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%dxray
                elseif(irprop == 12) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%dyray
                elseif(irprop == 13) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%dzray
                elseif(irprop == 14) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%dens
                elseif(irprop == 15) then
                  ySliceBack_send(iRay, ix, 1, kz) = ray(iRay, ix, 1, kz)%dphi
                end if
              end do
            end if
          end do
        end do

        do kz = 0, nz + 1
          do ix = 0, nx + 1
            if(nRay(ix, ny, kz) > 0) then
              do iRay = 1, nRay(ix, ny, kz)
                if(irprop == 1) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%x
                elseif(irprop == 2) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%y
                elseif(irprop == 3) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%z
                elseif(irprop == 4) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%k
                elseif(irprop == 5) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%l
                elseif(irprop == 6) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%m
                elseif(irprop == 7) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%omega
                elseif(irprop == 8) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%dkray
                elseif(irprop == 9) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%dlray
                elseif(irprop == 10) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%dmray
                elseif(irprop == 11) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%dxray
                elseif(irprop == 12) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%dyray
                elseif(irprop == 13) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%dzray
                elseif(irprop == 14) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%dens
                elseif(irprop == 15) then
                  ySliceForw_send(iRay, ix, 1, kz) = ray(iRay, ix, ny, kz)%dphi
                end if
              end do
            end if
          end do
        end do

        if(irprop < 16) then
          ! back -> forw
          sendcount = nrmaxf * (nx + 2) * (nz + 2)
          recvcount = sendcount
          source = back
          dest = forw
          tag = 100

          call mpi_sendrecv(ySliceForw_send(1, 0, 1, 0), sendcount, &
              &mpi_double_precision, dest, tag, ySliceBack_recv(1, 0, 1, 0), &
              &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              &sts_back, ierror)

          ! forw -> back
          sendcount = nrmaxb * (nx + 2) * (nz + 2)
          recvcount = sendcount
          source = forw
          dest = back
          tag = 100

          call mpi_sendrecv(ySliceBack_send(1, 0, 1, 0), sendcount, &
              &mpi_double_precision, dest, tag, ySliceForw_recv(1, 0, 1, 0), &
              &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              &sts_forw, ierror)
        end if

        ! write auxiliary slice to ray field

        ! forw halos

        do kz = 0, nz + 1
          do ix = 0, nx + 1
            if(nRay(ix, ny + 1, kz) > 0) then
              do iRay = 1, nRay(ix, ny + 1, kz)
                if(irprop == 1) then
                  ray(iRay, ix, ny + 1, kz)%x = ySliceForw_recv(iRay, ix, 1, kz)
                elseif(irprop == 2) then
                  ray(iRay, ix, ny + 1, kz)%y = ySliceForw_recv(iRay, ix, 1, kz)
                elseif(irprop == 3) then
                  ray(iRay, ix, ny + 1, kz)%z = ySliceForw_recv(iRay, ix, 1, kz)
                elseif(irprop == 4) then
                  ray(iRay, ix, ny + 1, kz)%k = ySliceForw_recv(iRay, ix, 1, kz)
                elseif(irprop == 5) then
                  ray(iRay, ix, ny + 1, kz)%l = ySliceForw_recv(iRay, ix, 1, kz)
                elseif(irprop == 6) then
                  ray(iRay, ix, ny + 1, kz)%m = ySliceForw_recv(iRay, ix, 1, kz)
                elseif(irprop == 7) then
                  ray(iRay, ix, ny + 1, kz)%omega = ySliceForw_recv(iRay, ix, &
                      &1, kz)
                elseif(irprop == 8) then
                  ray(iRay, ix, ny + 1, kz)%dkray = ySliceForw_recv(iRay, ix, &
                      &1, kz)
                elseif(irprop == 9) then
                  ray(iRay, ix, ny + 1, kz)%dlray = ySliceForw_recv(iRay, ix, &
                      &1, kz)
                elseif(irprop == 10) then
                  ray(iRay, ix, ny + 1, kz)%dmray = ySliceForw_recv(iRay, ix, &
                      &1, kz)
                elseif(irprop == 11) then
                  ray(iRay, ix, ny + 1, kz)%dxray = ySliceForw_recv(iRay, ix, &
                      &1, kz)
                elseif(irprop == 12) then
                  ray(iRay, ix, ny + 1, kz)%dyray = ySliceForw_recv(iRay, ix, &
                      &1, kz)
                elseif(irprop == 13) then
                  ray(iRay, ix, ny + 1, kz)%dzray = ySliceForw_recv(iRay, ix, &
                      &1, kz)
                elseif(irprop == 14) then
                  ray(iRay, ix, ny + 1, kz)%dens = ySliceForw_recv(iRay, ix, &
                      &1, kz)
                elseif(irprop == 15) then
                  ray(iRay, ix, ny + 1, kz)%dphi = ySliceForw_recv(iRay, ix, &
                      &1, kz)
                end if
              end do
            end if

            if(irprop == 2) then
              ! backwardmost cpu:
              ! for ray volumes in last cell in y direction and in the
              ! ghost cell after that (that might have been handed
              ! over from the forwardmost cpu) make sure that yr is
              ! adjusted accordingly

              if(js + nby + ny - 1 == sizeY) then
                do jy = ny, ny + 1
                  if(nRay(ix, jy, kz) > 0) then
                    do iRay = 1, nRay(ix, jy, kz)
                      yr = ray(iRay, ix, jy, kz)%y

                      yrt = yr + ly(1) - ly(0)

                      if(abs(yrt - y(jy + jy0)) < abs(yr - y(jy + jy0))) then
                        yr = yrt
                      end if

                      ray(iRay, ix, jy, kz)%y = yr
                    end do
                  end if
                end do
              end if
            else if(irprop > 15) then
              if(nRay(ix, ny + 1, kz) > 0) then
                do iRay = 1, nRay(ix, ny + 1, kz)
                  if(irprop == 16) then
                    ray(iRay, ix, ny + 1, kz)%area_xk = ray(iRay, ix, ny + 1, &
                        &kz)%dxray * ray(iRay, ix, ny + 1, kz)%dkray
                  elseif(irprop == 17) then
                    ray(iRay, ix, ny + 1, kz)%area_yl = ray(iRay, ix, ny + 1, &
                        &kz)%dyray * ray(iRay, ix, ny + 1, kz)%dlray
                  elseif(irprop == 18) then
                    ray(iRay, ix, ny + 1, kz)%area_zm = ray(iRay, ix, ny + 1, &
                        &kz)%dzray * ray(iRay, ix, ny + 1, kz)%dmray
                  end if
                end do
              end if
            end if
          end do
        end do

        ! back halos

        do kz = 0, nz + 1
          do ix = 0, nx + 1
            if(nRay(ix, 0, kz) > 0) then
              do iRay = 1, nRay(ix, 0, kz)
                if(irprop == 1) then
                  ray(iRay, ix, 0, kz)%x = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 2) then
                  ray(iRay, ix, 0, kz)%y = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 3) then
                  ray(iRay, ix, 0, kz)%z = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 4) then
                  ray(iRay, ix, 0, kz)%k = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 5) then
                  ray(iRay, ix, 0, kz)%l = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 6) then
                  ray(iRay, ix, 0, kz)%m = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 7) then
                  ray(iRay, ix, 0, kz)%omega = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 8) then
                  ray(iRay, ix, 0, kz)%dkray = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 9) then
                  ray(iRay, ix, 0, kz)%dlray = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 10) then
                  ray(iRay, ix, 0, kz)%dmray = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 11) then
                  ray(iRay, ix, 0, kz)%dxray = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 12) then
                  ray(iRay, ix, 0, kz)%dyray = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 13) then
                  ray(iRay, ix, 0, kz)%dzray = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 14) then
                  ray(iRay, ix, 0, kz)%dens = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 15) then
                  ray(iRay, ix, 0, kz)%dphi = ySliceBack_recv(iRay, ix, 1, kz)
                end if
              end do
            end if

            if(irprop == 2) then
              ! forwardmost cpu:
              ! for ray volumes in first cell in y direction and in
              ! the ghost cell before that (that might have been
              ! handed over from the backwardmost cpu) make sure
              ! that yr is adjusted accordingly

              if(js + nby == 1) then
                do jy = 0, 1
                  if(nRay(ix, jy, kz) > 0) then
                    do iRay = 1, nRay(ix, jy, kz)
                      yr = ray(iRay, ix, jy, kz)%y

                      yrt = yr - ly(1) + ly(0)

                      if(abs(yrt - y(jy + jy0)) < abs(yr - y(jy + jy0))) then
                        yr = yrt
                      end if

                      ray(iRay, ix, jy, kz)%y = yr
                    end do
                  end if
                end do
              end if
            elseif(irprop > 15) then
              if(nRay(ix, 0, kz) > 0) then
                do iRay = 1, nRay(ix, 0, kz)
                  if(irprop == 16) then
                    ray(iRay, ix, 0, kz)%area_xk = ray(iRay, ix, 0, kz)%dxray &
                        &* ray(iRay, ix, 0, kz)%dkray
                  elseif(irprop == 17) then
                    ray(iRay, ix, 0, kz)%area_yl = ray(iRay, ix, 0, kz)%dyray &
                        &* ray(iRay, ix, 0, kz)%dlray
                  elseif(irprop == 18) then
                    ray(iRay, ix, 0, kz)%area_zm = ray(iRay, ix, 0, kz)%dzray &
                        &* ray(iRay, ix, 0, kz)%dmray
                  end if
                end do
              end if
            end if
          end do
        end do
      end do

      deallocate(ySliceBack_send)
      deallocate(ySliceForw_send)

      deallocate(ySliceBack_recv)
      deallocate(ySliceForw_recv)
    else
      ! only 1 cpu in y direction

      do kz = 0, nz + 1
        do ix = 0, nx + 1
          if(nRay(ix, 0, kz) > 0) then
            do iRay = 1, nRay(ix, 0, kz)
              ray(iRay, ix, 0, kz) = ray(iRay, ix, ny, kz)
            end do
          end if

          do jy = 0, 1
            if(nRay(ix, jy, kz) > 0) then
              do iRay = 1, nRay(ix, jy, kz)
                yr = ray(iRay, ix, jy, kz)%y

                yrt = yr - ly(1) + ly(0)

                if(abs(yrt - y(jy + jy0)) < abs(yr - y(jy + jy0))) then
                  yr = yrt
                end if

                ray(iRay, ix, jy, kz)%y = yr
              end do
            end if
          end do
        end do
      end do

      do kz = 0, nz + 1
        do ix = 0, nx + 1
          if(nRay(ix, ny + 1, kz) > 0) then
            do iRay = 1, nRay(ix, ny + 1, kz)
              ray(iRay, ix, ny + 1, kz) = ray(iRay, ix, 1, kz)
            end do
          end if

          do jy = ny, ny + 1
            if(nRay(ix, jy, kz) > 0) then
              do iRay = 1, nRay(ix, jy, kz)
                yr = ray(iRay, ix, jy, kz)%y

                yrt = yr + ly(1) - ly(0)

                if(abs(yrt - y(jy + jy0)) < abs(yr - y(jy + jy0))) then
                  yr = yrt
                end if

                ray(iRay, ix, jy, kz)%y = yr
              end do
            end if
          end do
        end do
      end do
    end if

  end subroutine setboundary_rayvol_y

  ! ----------------------------------------------------------------------

  subroutine setboundary_irshift_y(irsb, irsf)

    ! -------------------------------------------------------------------
    ! (periodic) boundary conditions in x direction for indices of
    ! ray volumes shifted backward or forward
    ! -------------------------------------------------------------------

    implicit none

    ! argument list
    integer, dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), intent(inout) &
        &:: irsb, irsf

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: jy, irprop, nrmxbl, nrmxfl, nrmaxb, nrmaxf
    integer :: ix, kz, iRay

    ! auxiliary fields
    integer, dimension(:, :, :, :), allocatable :: ySliceBack_send, &
        &ySliceForw_send
    integer, dimension(:, :, :, :), allocatable :: ySliceBack_recv, &
        &ySliceForw_recv

    if(yBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_irshift_y must be  &
          &periodic!'
    end if

    ! boundary conditions for number of ray volumes per cell

    call setboundary_nray_y

    if(jdim > 1) then
      ! more than 1 cpu in x direction

      call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

      ! find global maximum of number of ray volumes per cell
      ! (at forward and backward edges of the cpu domains)

      nrmxbl = maxval(nRay(:, 1, :))
      nrmxfl = maxval(nRay(:, ny, :))

      call mpi_reduce(nrmxbl, nrmaxb, 1, mpi_integer, mpi_max, root, comm, &
          &ierror)
      call mpi_bcast(nrmaxb, 1, mpi_integer, root, comm, ierror)

      call mpi_reduce(nrmxfl, nrmaxf, 1, mpi_integer, mpi_max, root, comm, &
          &ierror)
      call mpi_bcast(nrmaxf, 1, mpi_integer, root, comm, ierror)

      allocate(ySliceBack_send(nrmaxb, 0:nx + 1, 1, 0:nz + 1))
      allocate(ySliceForw_send(nrmaxf, 0:nx + 1, 1, 0:nz + 1))

      allocate(ySliceBack_recv(nrmaxf, 0:nx + 1, 1, 0:nz + 1))
      allocate(ySliceForw_recv(nrmaxb, 0:nx + 1, 1, 0:nz + 1))

      do irprop = 1, 2
        ! read slice into contiguous array

        do kz = 0, nz + 1
          do ix = 0, nx + 1
            if(nRay(ix, 1, kz) > 0) then
              do iRay = 1, nRay(ix, 1, kz)
                if(irprop == 1) then
                  ySliceBack_send(iRay, ix, 1, kz) = irsb(iRay, ix, 1, kz)
                elseif(irprop == 2) then
                  ySliceBack_send(iRay, ix, 1, kz) = irsf(iRay, ix, 1, kz)
                end if
              end do
            end if
          end do
        end do

        do kz = 0, nz + 1
          do ix = 0, nx + 1
            if(nRay(ix, ny, kz) > 0) then
              do iRay = 1, nRay(ix, ny, kz)
                if(irprop == 1) then
                  ySliceForw_send(iRay, ix, 1, kz) = irsb(iRay, ix, ny, kz)
                elseif(irprop == 2) then
                  ySliceForw_send(iRay, ix, 1, kz) = irsf(iRay, ix, ny, kz)
                end if
              end do
            end if
          end do
        end do

        ! back -> forw
        sendcount = nrmaxf * (nx + 2) * (nz + 2)
        recvcount = sendcount
        source = back
        dest = forw
        tag = 100

        call mpi_sendrecv(ySliceForw_send(1, 0, 1, 0), sendcount, mpi_integer, &
            &dest, tag, ySliceBack_recv(1, 0, 1, 0), recvcount, mpi_integer, &
            &source, mpi_any_tag, comm, sts_back, ierror)

        ! forw -> back
        sendcount = nrmaxb * (nx + 2) * (nz + 2)
        recvcount = sendcount
        source = forw
        dest = back
        tag = 100

        call mpi_sendrecv(ySliceBack_send(1, 0, 1, 0), sendcount, mpi_integer, &
            &dest, tag, ySliceForw_recv(1, 0, 1, 0), recvcount, mpi_integer, &
            &source, mpi_any_tag, comm, sts_forw, ierror)

        ! write auxiliary slice to ray field

        ! forw halos

        do kz = 0, nz + 1
          do ix = 0, nx + 1
            if(nRay(ix, ny + 1, kz) > 0) then
              do iRay = 1, nRay(ix, ny + 1, kz)
                if(irprop == 1) then
                  irsb(iRay, ix, ny + 1, kz) = ySliceForw_recv(iRay, ix, 1, kz)
                elseif(irprop == 2) then
                  irsf(iRay, ix, ny + 1, kz) = ySliceForw_recv(iRay, ix, 1, kz)
                end if
              end do
            end if
          end do
        end do

        ! back halos

        do kz = 0, nz + 1
          do ix = 0, nx + 1
            if(nRay(ix, 0, kz) > 0) then
              do iRay = 1, nRay(ix, 0, kz)
                if(irprop == 1) then
                  irsb(iRay, ix, 0, kz) = ySliceBack_recv(iRay, ix, 1, kz)
                elseif(irprop == 2) then
                  irsf(iRay, ix, 0, kz) = ySliceBack_recv(iRay, ix, 1, kz)
                end if
              end do
            end if
          end do
        end do
      end do

      deallocate(ySliceBack_send)
      deallocate(ySliceForw_send)

      deallocate(ySliceBack_recv)
      deallocate(ySliceForw_recv)
    else
      ! only 1 cpu in x direction

      irsb(:, :, 0, :) = irsb(:, :, ny, :)
      irsb(:, :, ny + 1, :) = irsb(:, :, 1, :)

      irsf(:, :, 0, :) = irsf(:, :, ny, :)
      irsf(:, :, ny + 1, :) = irsf(:, :, 1, :)
    end if

  end subroutine setboundary_irshift_y

  ! ----------------------------------------------------------------------

  subroutine setboundary_nshift_y(nshb, nshf)

    ! -------------------------------------------------------------------
    ! (periodic) boundary conditions in x direction for number of
    ! ray volumes shifted to the back or forw
    ! -------------------------------------------------------------------

    implicit none

    ! argument list
    integer, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(inout) :: nshb, &
        &nshf

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: jy, irprop

    ! auxiliary fields
    integer, dimension(0:nx + 1, 1, 0:nz + 1) :: ySliceBack_send, &
        &ySliceForw_send
    integer, dimension(0:nx + 1, 1, 0:nz + 1) :: ySliceBack_recv, &
        &ySliceForw_recv

    if(yBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_nshift_y must be  &
          &periodic!'
    end if

    if(jdim > 1) then
      ! more than 1 cpu in x direction

      call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

      ! slice size
      sendcount = (nx + 2) * (nz + 2)
      recvcount = sendcount

      do irprop = 1, 2
        ! read slice into contiguous array

        if(irprop == 1) then
          ySliceBack_send(:, 1, :) = nshb(:, 1, :)
          ySliceForw_send(:, 1, :) = nshb(:, ny, :)
        elseif(irprop == 2) then
          ySliceBack_send(:, 1, :) = nshf(:, 1, :)
          ySliceForw_send(:, 1, :) = nshf(:, ny, :)
        end if

        ! back -> forw
        source = back
        dest = forw
        tag = 100

        call mpi_sendrecv(ySliceForw_send(0, 1, 0), sendcount, mpi_integer, &
            &dest, tag, ySliceBack_recv(0, 1, 0), recvcount, mpi_integer, &
            &source, mpi_any_tag, comm, sts_back, ierror)

        ! forw -> back
        source = forw
        dest = back
        tag = 100

        call mpi_sendrecv(ySliceBack_send(0, 1, 0), sendcount, mpi_integer, &
            &dest, tag, ySliceForw_recv(0, 1, 0), recvcount, mpi_integer, &
            &source, mpi_any_tag, comm, sts_forw, ierror)

        ! write auxiliary slice to ray field

        if(irprop == 1) then
          ! forw halos
          nshb(:, ny + 1, :) = ySliceForw_recv(:, 1, :)

          ! back halos
          nshb(:, 0, :) = ySliceBack_recv(:, 1, :)
        elseif(irprop == 2) then
          nshf(:, ny + 1, :) = ySliceForw_recv(:, 1, :)
          nshf(:, 0, :) = ySliceBack_recv(:, 1, :)
        end if
      end do
    else
      ! only 1 cpu in x direction

      nshb(:, 0, :) = nshb(:, ny, :)
      nshb(:, ny + 1, :) = nshb(:, 1, :)

      nshf(:, 0, :) = nshf(:, ny, :)
      nshf(:, ny + 1, :) = nshf(:, 1, :)
    end if

  end subroutine setboundary_nshift_y

  ! ----------------------------------------------------------------------

  subroutine setboundary_rayvol_z(ray)

    ! Update ray volumes at the vertical boundaries.

    implicit none

    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(inout) :: ray

    real :: xr, yr, zr, zrt
    real :: dzr
    real :: wnrm

    integer :: ix, jy, kz
    integer :: ixrv, jyrv
    integer :: ix0, jy0
    integer :: nrlc

    ix0 = is + nbx - 1
    jy0 = js + nby - 1

    select case(zBoundary)
    case("periodic")

      nRay(:, :, 0) = nRay(:, :, nz)
      nRay(:, :, nz + 1) = nRay(:, :, 1)

      do jy = 0, ny + 1
        do ix = 0, nx + 1
          if(nRay(ix, jy, 0) > 0) then
            do iRay = 1, nRay(ix, jy, 0)
              ray(iRay, ix, jy, 0) = ray(iRay, ix, jy, nz)
            end do
          end if

          do kz = 0, 1
            if(nRay(ix, jy, kz) > 0) then
              do iRay = 1, nRay(ix, jy, kz)
                zr = ray(iRay, ix, jy, kz)%z
                zrt = zr - lz(1) + lz(0)

                if(topography) then
                  xr = ray(iRay, ix, jy, kz)%x
                  yr = ray(iRay, ix, jy, kz)%y
                  ixrv = nint((xr - lx(0)) / dx + 0.5) - ix0
                  jyrv = nint((yr - ly(0)) / dy + 0.5) - jy0
                  if(abs(zrt - zTFC(ixrv, jyrv, kz)) < abs(zr - zTFC(ixrv, &
                      &jyrv, kz))) then
                    zr = zrt
                  end if
                else
                  if(abs(zrt - z(kz)) < abs(zr - z(kz))) then
                    zr = zrt
                  end if
                end if

                ray(iRay, ix, jy, kz)%z = zr
              end do
            end if
          end do
        end do
      end do

      do jy = 0, ny + 1
        do ix = 0, nx + 1
          if(nRay(ix, jy, nz + 1) > 0) then
            do iRay = 1, nRay(ix, jy, nz + 1)
              ray(iRay, ix, jy, nz + 1) = ray(iRay, ix, jy, 1)
            end do
          end if

          do kz = nz, nz + 1
            if(nRay(ix, jy, kz) > 0) then
              do iRay = 1, nRay(ix, jy, kz)
                zr = ray(iRay, ix, jy, kz)%z
                zrt = zr + lz(1) - lz(0)

                if(topography) then
                  xr = ray(iRay, ix, jy, kz)%x
                  yr = ray(iRay, ix, jy, kz)%y
                  ixrv = nint((xr - lx(0)) / dx + 0.5) - ix0
                  jyrv = nint((yr - ly(0)) / dy + 0.5) - jy0
                  if(abs(zrt - zTFC(ixrv, jyrv, kz)) < abs(zr - zTFC(ixrv, &
                      &jyrv, kz))) then
                    zr = zrt
                  end if
                else
                  if(abs(zrt - z(kz)) < abs(zr - z(kz))) then
                    zr = zrt
                  end if
                end if

                ray(iRay, ix, jy, kz)%z = zr
              end do
            end if
          end do
        end do
      end do

    case("solid_wall")

      ! Reflect ray volumes at the lower boundary.

      do kz = 1, 2
        do jy = 0, ny + 1
          do ix = 0, nx + 1
            if(nRay(ix, jy, kz) <= 0) cycle
            do iRay = 1, nRay(ix, jy, kz)

              zr = ray(iRay, ix, jy, kz)%z
              dzr = ray(iRay, ix, jy, kz)%dzray
              wnrm = ray(iRay, ix, jy, kz)%m

              if(topography) then
                xr = ray(iRay, ix, jy, kz)%x
                yr = ray(iRay, ix, jy, kz)%y
                ixrv = nint((xr - lx(0)) / dx + 0.5) - ix0
                jyrv = nint((yr - ly(0)) / dy + 0.5) - jy0
                if(zr - 0.5 * dzr < topography_surface(ixrv, jyrv)) then
                  ray(iRay, ix, jy, kz)%z = 2.0 * topography_surface(ixrv, &
                      &jyrv) - zr + dzr
                  ray(iRay, ix, jy, kz)%m = - wnrm
                end if
              else
                if(zr - 0.5 * dzr < lz(0)) then
                  ray(iRay, ix, jy, kz)%z = 2.0 * lz(0) - zr + dzr
                  ray(iRay, ix, jy, kz)%m = - wnrm
                end if
              end if
            end do
          end do
        end do
      end do

      ! Cut ray volumes at the upper boundary.

      do kz = nz - 1, nz
        do jy = 0, ny + 1
          do ix = 0, nx + 1
            if(nRay(ix, jy, kz) <= 0) cycle
            nrlc = 0
            do iRay = 1, nRay(ix, jy, kz)

              zr = ray(iRay, ix, jy, kz)%z
              dzr = ray(iRay, ix, jy, kz)%dzray

              if(zr - 0.5 * dzr > lz(1)) cycle

              if(zr + 0.5 * dzr > lz(1)) then
                ray(iRay, ix, jy, kz)%dzray = lz(1) - zr + 0.5 * dzr
                ray(iRay, ix, jy, kz)%z = lz(1) - 0.5 * ray(iRay, ix, jy, &
                    &kz)%dzray
                ray(iRay, ix, jy, kz)%area_zm = ray(iRay, ix, jy, kz)%dzray &
                    &* ray(iRay, ix, jy, kz)%dmray
              end if

              nrlc = nrlc + 1
              ray(nrlc, ix, jy, kz) = ray(iRay, ix, jy, kz)
            end do
            nRay(ix, jy, kz) = nrlc
          end do
        end do
      end do

    case default
      stop "Error in setboundary_rayvol_z: Unknown case zBoundary!"
    end select

  end subroutine setboundary_rayvol_z

  ! ----------------------------------------------------------------------

  subroutine setboundary_nray_x

    ! -------------------------------------------------------------------
    ! (periodic) boundary conditions in x direction for number nRay of
    ! ray volumes per cell
    ! -------------------------------------------------------------------

    implicit none

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: ix, irprop

    ! auxiliary fields
    integer, dimension(1, 0:ny + 1, 0:nz + 1) :: xSliceLeft_send, &
        &xSliceRight_send
    integer, dimension(1, 0:ny + 1, 0:nz + 1) :: xSliceLeft_recv, &
        &xSliceRight_recv

    if(xBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_nray_x must be  periodic!'
    end if

    if(idim > 1) then
      ! more than 1 cpu in x direction

      call mpi_cart_shift(comm, 0, 1, left, right, ierror)

      ! slice size
      sendcount = (ny + 2) * (nz + 2)
      recvcount = sendcount

      ! read slice into contiguous array

      xSliceLeft_send(1, :, :) = nRay(1, :, :)
      xSliceRight_send(1, :, :) = nRay(nx, :, :)

      ! left -> right
      source = left
      dest = right
      tag = 100

      call mpi_sendrecv(xSliceRight_send(1, 0, 0), sendcount, mpi_integer, &
          &dest, tag, xSliceLeft_recv(1, 0, 0), recvcount, mpi_integer, &
          &source, mpi_any_tag, comm, sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xSliceLeft_send(1, 0, 0), sendcount, mpi_integer, &
          &dest, tag, xSliceRight_recv(1, 0, 0), recvcount, mpi_integer, &
          &source, mpi_any_tag, comm, sts_right, ierror)

      ! write auxiliary slice to ray field

      ! right halos
      nRay(nx + 1, :, :) = xSliceRight_recv(1, :, :)

      ! left halos
      nRay(0, :, :) = xSliceLeft_recv(1, :, :)
    else
      ! only 1 cpu in x direction

      nRay(0, :, :) = nRay(nx, :, :)
      nRay(nx + 1, :, :) = nRay(1, :, :)
    end if

  end subroutine setboundary_nray_x

  ! ----------------------------------------------------------------------

  subroutine setboundary_nray_y

    ! -------------------------------------------------------------------
    ! (periodic) boundary conditions in y direction for number nRay of
    ! ray volumes per cell
    ! -------------------------------------------------------------------

    implicit none

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: jy, irprop

    ! auxiliary fields
    integer, dimension(0:nx + 1, 1, 0:nz + 1) :: ySliceBack_send, &
        &ySliceForw_send
    integer, dimension(0:nx + 1, 1, 0:nz + 1) :: ySliceBack_recv, &
        &ySliceForw_recv

    if(yBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_nRay_y must be  periodic!'
    end if

    if(jdim > 1) then
      ! more than 1 cpu in x direction

      call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

      ! slice size
      sendcount = (nx + 2) * (nz + 2)
      recvcount = sendcount

      ! read slice into contiguous array

      ySliceBack_send(:, 1, :) = nRay(:, 1, :)
      ySliceForw_send(:, 1, :) = nRay(:, ny, :)

      ! back -> forw
      source = back
      dest = forw
      tag = 100

      call mpi_sendrecv(ySliceForw_send(0, 1, 0), sendcount, mpi_integer, &
          &dest, tag, ySliceBack_recv(0, 1, 0), recvcount, mpi_integer, &
          &source, mpi_any_tag, comm, sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(ySliceBack_send(0, 1, 0), sendcount, mpi_integer, &
          &dest, tag, ySliceForw_recv(0, 1, 0), recvcount, mpi_integer, &
          &source, mpi_any_tag, comm, sts_forw, ierror)

      ! write auxiliary slice to ray field

      ! forw halos
      nRay(:, ny + 1, :) = ySliceForw_recv(:, 1, :)

      ! back halos
      nRay(:, 0, :) = ySliceBack_recv(:, 1, :)
    else
      ! only 1 cpu in x direction

      nRay(:, 0, :) = nRay(:, ny, :)
      nRay(:, ny + 1, :) = nRay(:, 1, :)
    end if

  end subroutine setboundary_nray_y

  ! ----------------------------------------------------------------------

  subroutine orographic_source(var, ray, time, dt)

    ! Launch ray volumes according to the lower boundary condition for mountain
    ! waves.

    implicit none

    ! Declare variables.
    type(var_type), intent(in) :: var
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(inout) :: ray
    real, intent(in) :: time
    real, intent(in) :: dt
    integer :: ix, jy, kz
    integer :: ix2, jy2, kz2, ik, jl, km, iwm
    integer :: nrlc
    integer :: ix0, jy0
    real :: wnrk, wnrl, wnrm, wnrh
    real :: omir
    real :: f_cor_nd, NN_nd
    real :: zr, dzr
    real :: dk_ini_nd, dl_ini_nd, dm_ini_nd
    real :: pspvol
    real :: wadr, displm
    real :: long
    real :: cgrz

    ! Determine index offset.
    ix0 = is + nbx - 1
    jy0 = js + nby - 1

    ! Set Coriolis parameter.
    f_cor_nd = f_Coriolis_dim * tRef

    ! Compute local buoyancy frequency.
    if(.not. topography) then
      call stratification(z(1), 1, NN_nd)
    end if

    ! Set launch level.
    kz = 0

    ! Iterate over surface grid cells.
    do jy = 1, ny
      do ix = 1, nx

        ! Compute local buoyancy frequency.
        if(topography) then
          call stratification(zTFC(ix, jy, 1), 1, NN_nd)
        end if

        ! Iterate over surface ray volumes.
        do i_sfc = 1, n_sfc
          iRay = ir_sfc(i_sfc, ix, jy)

          ! Set surface indices.
          ix2 = ix2_sfc(i_sfc)
          jy2 = jy2_sfc(i_sfc)
          kz2 = kz2_sfc(i_sfc)
          ik = ik_sfc(i_sfc)
          jl = jl_sfc(i_sfc)
          km = km_sfc(i_sfc)
          iwm = iwm_sfc(i_sfc)

          ! Set wavenumbers.
          wnrk = k_spectrum(ix, jy, iwm)
          wnrl = l_spectrum(ix, jy, iwm)
          wnrh = sqrt(wnrk ** 2.0 + wnrl ** 2.0)

          ! Get vertical position and extent of old ray volume.
          if(iRay > 0) then
            zr = ray(iRay, ix, jy, kz)%z
            dzr = ray(iRay, ix, jy, kz)%dzray
          else
            zr = 0.0
            dzr = 0.0
          end if

          ! Compute intrinsic frequency from orographic wavenumbers.
          omir = - 0.5 * (var%u(ix, jy, 1) + var%u(ix - 1, jy, 1)) * wnrk &
              &- 0.5 * (var%v(ix, jy, 1) + var%v(ix, jy - 1, 1)) * wnrl

          ! Adjust the signs to be consistent with the chosen frequency branch.
          if(omir * branchr < 0.0) then
            omir = - omir
            wnrk = - wnrk
            wnrl = - wnrl
          end if

          ! Compute vertical wavenumber and wave-action density.
          if(omir ** 2 > f_cor_nd ** 2 .and. omir ** 2 < NN_nd) then

            ! Compute vertical wavenumber.
            wnrm = - branchr * sqrt(wnrh ** 2 * (NN_nd - omir ** 2) / (omir &
                &** 2 - f_cor_nd ** 2))

            ! Get orographic mode.
            displm = abs(topography_spectrum(ix, jy, iwm))

            ! Apply reduction due to blocked-layer formation.
            if(blocking) then
              long = sqrt(NN_nd / (0.25 * (var%u(ix, jy, 1) + var%u(ix - 1, &
                  &jy, 1)) ** 2.0 + 0.25 * (var%v(ix, jy, 1) + var%v(ix, jy &
                  &- 1, 1)) ** 2.0)) * sum(abs(topography_spectrum(ix, jy, :)))
              displm = displm * wave_amplitude_reduction(long)
            end if

            ! Compute wave-action density.
            if(topography) then
              wadr = 0.5 * rhoStratTFC(ix, jy, 1) * displm ** 2 * omir * (wnrh &
                  &** 2 + wnrm ** 2) / wnrh ** 2
            else
              wadr = 0.5 * rhoStrat(1) * displm ** 2 * omir * (wnrh ** 2 &
                  &+ wnrm ** 2) / wnrh ** 2
            end if

            ! Set to zero if something went wrong.
            if(wadr /= wadr .or. wnrm /= wnrm) then
              wadr = 0.0
              wnrm = 0.0
            end if

            ! Account for critical and reflecting levels.
          else
            wadr = 0.0
            wnrm = 0.0
          end if

          ! Check if a new ray volume is to be launched and clip the old one.
          ! Three cases are distinguished.
          ! (1) There is no ray volume with nonzero wave-action density. A new
          !     ray volume is launched.
          ! (2) There is a ray volume with nonzero wave-action density, which
          !     has partially passed the lower boundary. It is clipped and the
          !     part below the lower boundary discarded before a new ray volume
          !     is launched.
          ! (3) There is a ray volume with nonzero wave-action density, which
          !     has not yet crossed the lower boundary. It is replaced with a
          !     new one.
          if(steady_state) then
            if(iRay < 0) then
              nRay(ix, jy, kz) = nRay(ix, jy, kz) + 1
              iRay = nRay(ix, jy, kz)
              ir_sfc(i_sfc, ix, jy) = iRay
            end if

            if(wadr == 0.0) then
              ray(iRay, ix, jy, kz)%dens = 0.0
              cycle
            end if
          else
            if(wadr /= 0.0) then
              if(launch_algorithm == "scale") iRay = - 1

              ! Check for case (2).
              if(iRay > 0 .and. .not. topography .and. zr + 0.5 * dzr > z(kz) &
                  &+ 0.5 * dz) then

                ! Shift the old ray volume.
                nRay(ix, jy, kz + 1) = nRay(ix, jy, kz + 1) + 1
                nrlc = nRay(ix, jy, kz + 1)
                if(nrlc > nray_wrk) stop "Error in orographic_source: nrlc &
                    &> nray_wrk!"
                ray(nrlc, ix, jy, kz + 1) = ray(iRay, ix, jy, kz)

                ! Clip or extend the old ray volume.
                if(zr - 0.5 * dzr < z(kz) + 0.5 * dz .or. kz2 == 1) then
                  ray(nrlc, ix, jy, kz + 1)%dzray = zr + 0.5 * dzr - (z(kz) &
                      &+ 0.5 * dz)
                  ray(nrlc, ix, jy, kz + 1)%z = zr + 0.5 * dzr - 0.5 &
                      &* ray(nrlc, ix, jy, kz + 1)%dzray
                  ray(nrlc, ix, jy, kz + 1)%area_zm = ray(nrlc, ix, jy, kz &
                      &+ 1)%dzray * ray(nrlc, ix, jy, kz + 1)%dmray
                end if

                ! Check for case (2) in TFC.
              else if(iRay > 0 .and. topography .and. zr + 0.5 * dzr &
                  &> zTildeTFC(ix, jy, kz)) then

                ! Shift the old ray volume.
                nRay(ix, jy, kz + 1) = nRay(ix, jy, kz + 1) + 1
                nrlc = nRay(ix, jy, kz + 1)
                if(nrlc > nray_wrk) stop "Error in orographic_source: nrlc &
                    &> nray_wrk!"
                ray(nrlc, ix, jy, kz + 1) = ray(iRay, ix, jy, kz)

                ! Clip or extend the old ray volume.
                if(zr - 0.5 * dzr < zTildeTFC(ix, jy, kz) .or. kz2 == 1) then
                  ray(nrlc, ix, jy, kz + 1)%dzray = zr + 0.5 * dzr &
                      &- zTildeTFC(ix, jy, kz)
                  ray(nrlc, ix, jy, kz + 1)%z = zr + 0.5 * dzr - 0.5 &
                      &* ray(nrlc, ix, jy, kz + 1)%dzray
                  ray(nrlc, ix, jy, kz + 1)%area_zm = ray(nrlc, ix, jy, kz &
                      &+ 1)%dzray * ray(nrlc, ix, jy, kz + 1)%dmray
                end if

                ! Check for case (1).
              elseif(iRay < 0) then
                nRay(ix, jy, kz) = nRay(ix, jy, kz) + 1
                iRay = nRay(ix, jy, kz)
                if(iRay > nray_wrk) stop "Error in orographic_source: iRay &
                    &> nray_wrk!"
                ir_sfc(i_sfc, ix, jy) = iRay
              end if

              ! No ray volume is launched for zero wave-action density.
            else
              ir_sfc(i_sfc, ix, jy) = - 1
              cycle
            end if
          end if

          ! Scale the wave action density.
          if(.not. steady_state .and. launch_algorithm == "scale") then
            cgrz = wnrm * (f_cor_nd ** 2.0 - NN_nd) * wnrh ** 2.0 / omir &
                &/ (wnrh ** 2.0 + wnrm ** 2.0) ** 2.0
            if(topography) then
              wadr = wadr * dt * cgrz / jac(ix, jy, kz) / dz
            else
              wadr = wadr * dt * cgrz / dz
            end if
          end if

          ! Set physical ray-volume positions.
          ray(iRay, ix, jy, kz)%x = (x(ix + ix0) - 0.5 * dx + (ix2 - 0.5) * dx &
              &/ nrxl)
          ray(iRay, ix, jy, kz)%y = (y(jy + jy0) - 0.5 * dy + (jy2 - 0.5) * dy &
              &/ nryl)
          if(topography) then
            ray(iRay, ix, jy, kz)%z = (zTFC(ix, jy, kz) - 0.5 * jac(ix, jy, &
                &kz) * dz + (kz2 - 0.5) * jac(ix, jy, kz) * dz / nrzl)
          else
            ray(iRay, ix, jy, kz)%z = (z(kz) - 0.5 * dz + (kz2 - 0.5) * dz &
                &/ nrzl)
          end if

          ! Set physical ray-volume extent.
          ray(iRay, ix, jy, kz)%dxray = dx / nrxl
          ray(iRay, ix, jy, kz)%dyray = dy / nryl
          if(topography) then
            ray(iRay, ix, jy, kz)%dzray = jac(ix, jy, kz) * dz / nrzl
          else
            ray(iRay, ix, jy, kz)%dzray = dz / nrzl
          end if

          ! Compute spectral ray-volume extent.
          if(sizeX == 1) then
            dk_ini_nd = 0.0
          else
            dk_ini_nd = fac_dk_init * sqrt(wnrk ** 2.0 + wnrl ** 2.0)
          end if
          if(sizeY == 1) then
            dl_ini_nd = 0.0
          else
            dl_ini_nd = fac_dl_init * sqrt(wnrk ** 2.0 + wnrl ** 2.0)
          end if
          if(wnrm == 0.0) then
            stop 'Error in orographic_source: wnrm = 0!'
          else
            dm_ini_nd = fac_dm_init * abs(wnrm)
          end if

          ! Set spectral ray-volume position.
          ray(iRay, ix, jy, kz)%k = (wnrk - 0.5 * dk_ini_nd + (real(ik) - 0.5) &
              &* dk_ini_nd / nrk_init)
          ray(iRay, ix, jy, kz)%l = (wnrl - 0.5 * dl_ini_nd + (real(jl) - 0.5) &
              &* dl_ini_nd / nrl_init)
          ray(iRay, ix, jy, kz)%m = (wnrm - 0.5 * dm_ini_nd + (real(km) - 0.5) &
              &* dm_ini_nd / nrm_init)

          ! Set spectral ray-voume extent.
          ray(iRay, ix, jy, kz)%dkray = dk_ini_nd / nrk_init
          ray(iRay, ix, jy, kz)%dlray = dl_ini_nd / nrl_init
          ray(iRay, ix, jy, kz)%dmray = dm_ini_nd / nrm_init

          ! Set phase-space volume.
          ray(iRay, ix, jy, kz)%area_xk = ray(iRay, ix, jy, kz)%dxray &
              &* ray(iRay, ix, jy, kz)%dkray
          ray(iRay, ix, jy, kz)%area_yl = ray(iRay, ix, jy, kz)%dyray &
              &* ray(iRay, ix, jy, kz)%dlray
          ray(iRay, ix, jy, kz)%area_zm = ray(iRay, ix, jy, kz)%dzray &
              &* ray(iRay, ix, jy, kz)%dmray

          ! Compute spectral volume.
          pspvol = dm_ini_nd
          if(sizeX > 1) then
            pspvol = pspvol * dk_ini_nd
          end if
          if(sizeY > 1) then
            pspvol = pspvol * dl_ini_nd
          end if

          ! Set phase-space wave-action density.
          ray(iRay, ix, jy, kz)%dens = wadr / pspvol

          ! Set intrinsic frequency.
          ray(iRay, ix, jy, kz)%omega = omir
        end do
      end do
    end do

  end subroutine orographic_source

  ! ----------------------------------------------------------------------

  function interpolate_sponge(xlc, ylc, zlc) result(alpha)

    real :: xlc, ylc, zlc
    real :: alpha

    integer :: ix0, jy0

    integer :: ixl, ixr
    integer :: jyb, jyf
    integer :: kzlbd, kzlbu, kzlfd, kzlfu, kzrbd, kzrbu, kzrfd, kzrfu
    integer :: kzd, kzu

    real :: xl, xr
    real :: yb, yf
    real :: zlbd, zlbu, zlfd, zlfu, zrbd, zrbu, zrfd, zrfu
    real :: zbd, zbu, zfd, zfu
    real :: zd, zu

    real :: alphalbd, alphalbu, alphalfd, alphalfu, alpharbd, alpharbu, &
        &alpharfd, alpharfu
    real :: alphabd, alphabu, alphafd, alphafu
    real :: alphad, alphau

    real :: factor

    ix0 = is + nbx - 1
    jy0 = js + nby - 1

    ! Dermine closest points in horizontal direction.
    if(sizeX > 1) then
      ixl = floor((xlc - 0.5 * dx - lx(0)) / dx) + 1 - ix0
      ixr = ixl + 1
    else
      ixl = 1
      ixr = 1
    end if
    xl = x(ix0 + ixl)
    xr = x(ix0 + ixr)

    ! Determine closest points in meridional direction.
    if(sizeY > 1) then
      jyb = floor((ylc - 0.5 * dy - ly(0)) / dy) + 1 - jy0
      jyf = jyb + 1
    else
      jyb = 1
      jyf = 1
    end if
    yb = y(jy0 + jyb)
    yf = y(jy0 + jyf)

    ! Determine closest points in vertical direction and set interpolation
    ! values.
    if(topography) then
      kzlbd = max(1, floor((levelTFC(ixl, jyb, zlc) - 0.5 * dz - lz(0)) / dz) &
          &+ 1)
      kzlbu = max(1, floor((levelTFC(ixl, jyb, zlc) - 0.5 * dz - lz(0)) / dz) &
          &+ 2)
      if(kzlbd > nz) then
        kzlbd = nz
        kzlbu = nz
      end if
      zlbd = zTFC(ixl, jyb, kzlbd)
      zlbu = zTFC(ixl, jyb, kzlbu)
      alphalbd = alphaUnifiedSponge(ixl, jyb, kzlbd)
      alphalbu = alphaUnifiedSponge(ixl, jyb, kzlbu)
      kzlfd = max(1, floor((levelTFC(ixl, jyf, zlc) - 0.5 * dz - lz(0)) / dz) &
          &+ 1)
      kzlfu = max(1, floor((levelTFC(ixl, jyf, zlc) - 0.5 * dz - lz(0)) / dz) &
          &+ 2)
      if(kzlfd > nz) then
        kzlfd = nz
        kzlfu = nz
      end if
      zlfd = zTFC(ixl, jyf, kzlfd)
      zlfu = zTFC(ixl, jyf, kzlfu)
      alphalfd = alphaUnifiedSponge(ixl, jyf, kzlfd)
      alphalfu = alphaUnifiedSponge(ixl, jyf, kzlfu)
      kzrbd = max(1, floor((levelTFC(ixr, jyb, zlc) - 0.5 * dz - lz(0)) / dz) &
          &+ 1)
      kzrbu = max(1, floor((levelTFC(ixr, jyb, zlc) - 0.5 * dz - lz(0)) / dz) &
          &+ 2)
      if(kzrbd > nz) then
        kzrbd = nz
        kzrbu = nz
      end if
      zrbd = zTFC(ixr, jyb, kzrbd)
      zrbu = zTFC(ixr, jyb, kzrbu)
      alpharbd = alphaUnifiedSponge(ixr, jyb, kzrbd)
      alpharbu = alphaUnifiedSponge(ixr, jyb, kzrbu)
      kzrfd = max(1, floor((levelTFC(ixr, jyf, zlc) - 0.5 * dz - lz(0)) / dz) &
          &+ 1)
      kzrfu = max(1, floor((levelTFC(ixr, jyf, zlc) - 0.5 * dz - lz(0)) / dz) &
          &+ 2)
      if(kzrfd > nz) then
        kzrfd = nz
        kzrfu = nz
      end if
      zrfd = zTFC(ixr, jyf, kzrfd)
      zrfu = zTFC(ixr, jyf, kzrfu)
      alpharfd = alphaUnifiedSponge(ixr, jyf, kzrfd)
      alpharfu = alphaUnifiedSponge(ixr, jyf, kzrfu)
    else
      kzd = max(1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 1)
      kzu = max(1, floor((zlc - 0.5 * dz - lz(0)) / dz) + 2)
      if(kzd > nz) then
        kzd = nz
        kzu = nz
      end if
      zd = z(kzd)
      zu = z(kzu)
      alphalbd = alphaUnifiedSponge(ixl, jyb, kzd)
      alphalbu = alphaUnifiedSponge(ixl, jyb, kzu)
      alphalfd = alphaUnifiedSponge(ixl, jyf, kzd)
      alphalfu = alphaUnifiedSponge(ixl, jyf, kzu)
      alpharbd = alphaUnifiedSponge(ixr, jyb, kzd)
      alpharbu = alphaUnifiedSponge(ixr, jyb, kzu)
      alpharfd = alphaUnifiedSponge(ixr, jyf, kzd)
      alpharfu = alphaUnifiedSponge(ixr, jyf, kzu)
    end if

    ! Interpolate in x.
    if(sizeX > 1) then
      if(xr == xl) then
        factor = 0.0
      elseif(xlc > xr) then
        factor = 0.0
      elseif(xlc > xl) then
        factor = (xr - xlc) / (xr - xl)
      else
        factor = 1.0
      end if

      if(topography) then
        zbd = factor * zlbd + (1.0 - factor) * zrbd
        zbu = factor * zlbu + (1.0 - factor) * zrbu

        zfd = factor * zlfd + (1.0 - factor) * zrfd
        zfu = factor * zlfu + (1.0 - factor) * zrfu
      end if

      alphabd = factor * alphalbd + (1.0 - factor) * alpharbd
      alphabu = factor * alphalbu + (1.0 - factor) * alpharbu

      alphafd = factor * alphalfd + (1.0 - factor) * alpharfd
      alphafu = factor * alphalfu + (1.0 - factor) * alpharfu
    else
      if(topography) then
        zbd = zlbd
        zbu = zlbu

        zfd = zlfd
        zfu = zlfu
      end if

      alphabd = alphalbd
      alphabu = alphalbu

      alphafd = alphalfd
      alphafu = alphalfu
    end if

    ! Interpolate in y.
    if(sizeY > 1) then
      if(yf == yb) then
        factor = 0.0
      elseif(ylc > yf) then
        factor = 0.0
      elseif(ylc > yb) then
        factor = (yf - ylc) / (yf - yb)
      else
        factor = 1.0
      end if

      if(topography) then
        zd = factor * zbd + (1.0 - factor) * zfd
        zu = factor * zbu + (1.0 - factor) * zfu
      end if

      alphad = factor * alphabd + (1.0 - factor) * alphafd
      alphau = factor * alphabu + (1.0 - factor) * alphafu
    else
      if(topography) then
        zd = zbd
        zu = zbu
      end if

      alphad = alphabd
      alphau = alphabu
    end if

    ! Interpolate in z.
    if(zu == zd) then
      factor = 0.0
    elseif(zlc > zu) then
      factor = 0.0
    elseif(zlc > zd) then
      factor = (zu - zlc) / (zu - zd)
    else
      factor = 1.0
    end if

    alpha = factor * alphad + (1.0 - factor) * alphau

  end function interpolate_sponge

  ! ----------------------------------------------------------------------

  function wave_amplitude_reduction(long) result(factor)

    real :: long
    real :: factor
    real, parameter :: cc = 0.25

    if(long == 0.0) then
      factor = 1.0
    else
      factor = min(1.0, cc / long)
    end if

  end function wave_amplitude_reduction

  ! ----------------------------------------------------------------------

  subroutine calc_ice(ray, var, ray_varIce, ray_cloud)

    ! supplemements cell-centered volume forces by WKB force
    ! as well as the heating by entropy-flux convergence

    implicit none

    ! in/out variables
    !TEST***
    type(var_type), intent(inout) :: var
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, 0:nz + 1), &
        &intent(in) :: ray
    type(ice_rayType), dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(out) :: &
        &ray_varIce
    type(ice_rayType2), dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, NSCX, NSCY), &
        &intent(out) :: ray_cloud

    real :: dxi, dyi, dzi

    integer :: ixmin, ixmax
    integer :: jymin, jymax
    integer :: kzmin, kzmax
    integer :: ixr, jyr, kzr
    integer :: ix, jy, kz

    integer :: ix0, jy0

    real :: rhotot

    real :: wnrk
    real :: wnrl
    real :: wnrm
    real :: wnrh
    real :: dwnrk, dwnrl, dwnrm
    real :: f_cor_nd

    real :: omir

    real :: xr, yr, zr
    real :: dxr, dyr, dzr

    real :: fcpspx, fcpspy, fcpspz
    real :: NNR

    !SD
    real :: fcpswn
    real :: fcpsar, amprw, wadr
    real :: rho
    real :: thetaPrime, expPrime, wPrime
    real :: dphi
    complex :: w10, b11, theta11, pi12
    real :: theta0
    real :: dxx, dyy, dzz !compute phase a cell center
    !real :: dxo, dyo, dzo !overlap
    real :: dxoh, dyoh, dzoh !overlap
    real, parameter :: fxo = 1., fyo = fxo, fzo = fxo ! overlap = 1-fxo
    real, parameter :: fxoh = 1., fyoh = fxoh, fzoh = fxoh ! overlap = 1-fxoh
    real :: xrh, yrh, zrh, dxrh, dyrh, dzrh !
    real :: wnrk0
    real :: wnrl0
    real :: wnrm0
    integer :: crv
    real :: sigdex, sigdey, sigdez
    real :: is2pi3, vv, nn, dens_gauss
    !to be deleted
    !    real :: dotThetaprime, dotExpprime
    !    real :: pres, temp, theta, psi, exn_p, tPrime, pPrime, PiPrime,
    !    real ::  SIce, Qv, PiMean
    integer :: ii, jj
    real :: xsc, ysc
    logical :: apply

    !real, allocatable :: var_uu(:, :, :)
    real, allocatable :: wadr_sum(:, :, :)
    real, allocatable :: fcpsar_max(:, :, :)

    allocate(wadr_sum(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))
    allocate(fcpsar_max(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz))

    !Init
    wadr_sum = 0.
    fcpsar_max = 0.
    ray_varIce%wwp = 0.
    ray_varIce%epp = 0.
    ray_varIce%thp = 0.

    !TEST***
    var%OPT(:, :, :, 3) = 0.

    if(compute_cloudcover) then
      ray_cloud%wwp = 0.
      ray_cloud%epp = 0.
      ray_cloud%thp = 0.
    end if

    f_cor_nd = f_Coriolis_dim * tRef

    do kzrv = 0, nz
      do jyrv = 0, ny + 1
        ! loop including ghost cells in order to get all fluxes
        ! affecting a cell
        ! (assuming that ray volumes are not wider in y than dy)

        do ixrv = 0, nx + 1
          ! loop including ghost cells in order to get all fluxes
          ! affecting a cell
          ! (assuming that ray volumes are not wider in x than dx)

          if(nRay(ixrv, jyrv, kzrv) < 1) cycle

          crv = 0

          do iRay = 1, nRay(ixrv, jyrv, kzrv)

            ! skip counting ray volumes with zero wave-action density
            if(ray(iRay, ixrv, jyrv, kzrv)%dens == 0.0) cycle

            xr = ray(iRay, ixrv, jyrv, kzrv)%x
            yr = ray(iRay, ixrv, jyrv, kzrv)%y
            zr = ray(iRay, ixrv, jyrv, kzrv)%z

            dxr = ray(iRay, ixrv, jyrv, kzrv)%dxray
            dyr = ray(iRay, ixrv, jyrv, kzrv)%dyray
            dzr = ray(iRay, ixrv, jyrv, kzrv)%dzray

            ! vertical boundary conditions:
            ! zBoundary = 'periodic': implement periodicity
            ! zBoundary = 'solid_wall': skip counting ray volumes
            !                           that have completely left the
            !                           model domain

            ! FJApr2023
            apply = .false.
            if(.not. topography .and. zr < lz(0)) then
              apply = .true.
            else if(topography) then
              if(zr < topography_surface(ixrv, jyrv)) then
                apply = .true.
              end if
            end if
            if(apply) then
              select case(zBoundary)
              case("periodic")
                zr = lz(1) + mod(zr - lz(0), lz(1) - lz(0))
              case("solid_wall")
                if((.not. topography .and. zr + 0.5 * dzr < lz(0)) .or. &
                    &(topography .and. zr + 0.5 * dzr &
                    &< topography_surface(ixrv, jyrv))) cycle
              case default
                stop "calc_meanflow_effect: unknown case zBoundary"
              end select
            elseif(zr > lz(1)) then
              select case(zBoundary)
              case("periodic")
                zr = lz(0) + mod(zr - lz(1), lz(1) - lz(0))
              case("solid_wall")
                if(zr - 0.5 * dzr > lz(1)) cycle
              case default
                stop "calc_meanflow_effect: unknown case zBoundary"
              end select
            end if

            ! implement horizontal boundary conditions for ray-volume
            ! positions

            if(sizeX > 1) then
              if(xBoundary /= "periodic") then
                print *, 'ERROR in calc_meanflow_effect:  boundary conditions &
                    &in x must be periodic'
                stop
              end if

              ! for leftmost cpu make sure that xr in
              ! ghost cell to the left is between x(0) - dx
              ! and x(0)

              if(ixrv == 0 .and. is + nbx == 1) then
                if(xr > lx(1) - dx .and. xr < lx(1)) then
                  xr = xr - lx(1) + lx(0)
                elseif(xr > lx(0) - dx .and. xr < lx(0)) then
                  xr = xr
                else
                  print *, 'ERROR in calc_meanflow_effect:'
                  print *, 'xr =', xr
                  print *, 'but r.v. is in ghost cell to the  left of the &
                      &leftmost cpu so that  one should have either'
                  print *, lx(1) - dx, '= lx(1) - dx < xr < lx(1) =', lx(1), ' &
                      &or'
                  print *, lx(0) - dx, '= lx(0) - dx < xr < lx(0) =', lx(0)
                  stop
                end if
              end if

              ! for rightmost cpu make sure that xr in
              ! ghost cell to the right is between x(0) + L_x
              ! and x(0) + L_x + dx

              if(ixrv == nx + 1 .and. is + nbx + nx == sizeX) then
                if(xr > lx(0) .and. xr < lx(0) + dx) then
                  xr = xr + lx(1) - lx(0)
                elseif(xr > lx(1) .and. xr < lx(1) + dx) then
                  xr = xr
                else
                  print *, 'ERROR in calc_meanflow_effect:'
                  print *, 'xr =', xr
                  print *, 'but r.v. is in ghost cell to the  right of the &
                      &rightmost cpu so that  one should have either'
                  print *, lx(0), '= lx(0) < xr < lx(0) + dx =', lx(0) + dx, ' &
                      &or'
                  print *, lx(1), '= lx(1) < xr < lx(1) + dx =', lx(1) + dx
                  stop
                end if
              end if
            end if

            if(sizeY > 1) then
              if(yBoundary /= "periodic") then
                print *, 'ERROR in calc_meanflow_effect:  boundary conditions &
                    &in y must be periodic'
                stop
              end if

              ! for first cpu in y direct. make sure that yr in
              ! ghost cell in front is between y(0) - dy
              ! and y(0)

              if(jyrv == 0 .and. js + nby == 1) then
                if(yr > ly(1) - dy .and. yr < ly(1)) then
                  yr = yr - ly(1) + ly(0)
                elseif(yr > ly(0) - dy .and. yr < ly(0)) then
                  yr = yr
                else
                  print *, 'ERROR in calc_meanflow_effect:'
                  print *, 'yr =', yr
                  print *, 'but r.v. is in ghost cell in front of the first &
                      &cpu in y dir. so that  one should have either'
                  print *, ly(1) - dy, '= ly(1) - dy < yr < ly(1) =', ly(1), ' &
                      &or'
                  print *, ly(0) - dy, '= ly(0) - dy < yr < ly(0) =', ly(0)
                  stop
                end if
              end if

              ! for last cpu in y direction make sure that yr
              ! in ghost cell behind is between
              ! y(0) + L_y and y(0) + L_y + dy

              if(jyrv == ny + 1 .and. js + nby + ny == sizeY) then
                if(yr > ly(0) .and. yr < ly(0) + dy) then
                  yr = yr + ly(1) - ly(0)
                elseif(yr > ly(1) .and. yr < ly(1) + dy) then
                  yr = yr
                else
                  print *, 'ERROR in calc_meanflow_effect:'
                  print *, 'yr =', yr
                  print *, 'but r.v. is in ghost cell behind the last cpu in y &
                      &dir. so that  one should have either'
                  print *, ly(0), '= ly(0) < yr < ly(0) + dy =', ly(0) + dy, ' &
                      &or'
                  print *, ly(1), '= ly(1) < yr < ly(1) + dy =', ly(1) + dy
                  stop
                end if
              end if
            end if

            wnrk = ray(iRay, ixrv, jyrv, kzrv)%k
            wnrl = ray(iRay, ixrv, jyrv, kzrv)%l
            wnrm = ray(iRay, ixrv, jyrv, kzrv)%m

            dwnrk = ray(iRay, ixrv, jyrv, kzrv)%dkray
            dwnrl = ray(iRay, ixrv, jyrv, kzrv)%dlray
            dwnrm = ray(iRay, ixrv, jyrv, kzrv)%dmray

            wnrh = sqrt(wnrk ** 2 + wnrl ** 2)

            apply = .false.
            if((.not. topography .and. zr < lz(0) - dz)) then
              apply = .true.
            else if(topography) then
              if(zr < topography_surface(ixrv, jyrv) - jac(ixrv, jyrv, kzrv) &
                  &* dz) then
                apply = .true.
              end if
            end if
            if(apply) then
              !            if((.not. topography .and. zr < lz(0) - dz) .or. (topography .and. &
              !                zr < topography_surface(ixrv, jyrv) - jac(ixrv, jyrv, kzrv) &
              !                * dz)) then
              print *, 'ERROR IN calc_meanflow_effect: RAY VOLUME', iRay, 'in &
                  &cell', ixrv, jyrv, kzrv, 'TOO LOW'
              stop
            end if

            call stratification(zr, 1, NNr)

            omir = branchr * sqrt(NNr * wnrh ** 2 + f_cor_nd ** 2 * wnrm ** 2) &
                &/ sqrt(wnrh ** 2 + wnrm ** 2)

            ! indices of range of cells touched by a ray volume

            if(sizeX > 1) then
              ! last x-index leftmost of cpu
              ix0 = is + nbx - 1

              ! SD
              if(average_cell) then
                ixmin = floor((xr - dxr * 0.5 - lx(0)) / dx) + 1 - ix0
                ixmax = floor((xr + dxr * 0.5 - lx(0)) / dx) + 1 - ix0
              elseif(average_cell_3) then
                ixmin = ixrv - 1
                ixmax = ixrv + 1
              else
                ixmin = ixrv
                ixmax = ixrv
              end if

              if(ixmin > nx + 1) then
                print *, 'ixmin =', ixmin, '> nx+1 = ', nx + 1
                print *, 'ixrv = ', ixrv
                print *, 'xr =', xr
                print *, 'dxr =', dxr
                print *, 'dx =', dx
                print *, 'lx(0) =', lx(0)
                print *, 'lx(1) =', lx(1)
                print *, 'ix0 =', ix0
                print *, 'floor((xr - dxr*0.5 - lx(0)) / dx) + 1 =', floor((xr &
                    &- dxr * 0.5 - lx(0)) / dx) + 1
                stop
              else
                ! no fluxes calculated for the ghost cells
                ! (that are taken care of by the boundary-condition
                ! routines)

                ixmin = max(ixmin, 1)
              end if

              if(ixmax < 0) then
                print *, 'ixmax =', ixmax, '< 0'
                stop
              else
                ! no fluxes calculated for the ghost cells
                ! (that are taken care of by the boundary-condition
                ! routines)

                ixmax = min(ixmax, nx)
              end if
            else
              ixmin = 1
              ixmax = 1

              !SD
              ix0 = 0

            end if

            if(sizeY > 1) then
              ! last y-index in front of cpu
              jy0 = js + nby - 1

              !SD
              if(average_cell) then
                jymin = floor((yr - dyr * 0.5 - ly(0)) / dy) + 1 - jy0
                jymax = floor((yr + dyr * 0.5 - ly(0)) / dy) + 1 - jy0
              elseif(average_cell_3) then
                jymin = jyrv - 1
                jymax = jyrv + 1
              else
                jymin = jyrv
                jymax = jyrv
              end if

              if(jymin > ny + 1) then
                print *, 'jymin =', jymin, '> ny+1 = ', ny + 1
                stop
              else
                ! no fluxes calculated for the ghost cells
                ! (that are taken care of by the boundary-condition
                ! routines)

                jymin = max(jymin, 1)
              end if

              if(jymax < 0) then
                print *, 'jymax =', jymax, '< 0'
                stop
              else
                ! no fluxes calculated for the ghost cells
                ! (that are taken care of by the boundary-condition
                ! routines)

                jymax = min(jymax, ny)
              end if
            else
              jymin = 1
              jymax = 1

              !SD
              jy0 = 0
            end if

            if(topography) then

              print *, 'superposition of GW not implemented yet for topography &
                  &case'
              stop

              if(include_tracer) then
                stop "Tracer, topography, and ray tracer not possible."
              end if

              do ix = ixmin, ixmax
                if(sizeX > 1) then
                  dxi = (min((xr + dxr * 0.5), lx(0) + (ix + ix0) * dx) &
                      &- max((xr - dxr * 0.5), lx(0) + (ix + ix0 - 1) * dx))

                  fcpspx = dwnrk * dxi / dx

                  !SD
                  !fcpswn = fcpswn * dwnrk

                else
                  fcpspx = 1.0
                end if

                do jy = jymin, jymax
                  if(sizeY > 1) then
                    dyi = (min((yr + dyr * 0.5), ly(0) + (jy + jy0) * dy) &
                        &- max((yr - dyr * 0.5), ly(0) + (jy + jy0 - 1) * dy))

                    fcpspy = dwnrl * dyi / dy
                    !SD
                    !fcpswn = fcpswn * dwnrl
                  else
                    fcpspy = 1.0
                  end if

                  if(average_cell) then
                    ! Jacobian is height-independent!
                    kzmin = max(1, floor((zr - dzr * 0.5 &
                        &- topography_surface(ix, jy)) / jac(ix, jy, 0) / dz) &
                        &+ 1)
                    kzmax = min(nz, floor((zr + dzr * 0.5 &
                        &- topography_surface(ix, jy)) / jac(ix, jy, 0) / dz) &
                        &+ 1)
                  elseif(average_cell_3) then
                    print *, 'average_cell_3 not working with topography'
                    stop
                  else
                    kzmin = kzrv
                    kzmax = kzrv
                  end if

                  do kz = kzmin, kzmax
                    dzi = (min((zr + dzr * 0.5), topography_surface(ix, jy) &
                        &+ kz * jac(ix, jy, kz) * dz) - max((zr - dzr * 0.5), &
                        &topography_surface(ix, jy) + (kz - 1) * jac(ix, jy, &
                        &kz - 1) * dz))

                    fcpspz = dwnrm * dzi / jac(ix, jy, kz) / dz

                    !SD
                    !fcpswn = fcpswn * dwnrm

                    wadr = fcpspx * fcpspy * fcpspz * ray(iRay, ixrv, jyrv, &
                        &kzrv)%dens

                    wadr_sum(ix, jy, kz) = wadr_sum(ix, jy, kz) + wadr !integrated waveaction density
                    fcpsar = dzi / jac(ix, jy, kz) / dz

                    if(sizeX > 1) then
                      fcpsar = (dxi / dx) * fcpsar
                    end if
                    if(sizeY > 1) then
                      fcpsar = (dyi / dy) * fcpsar
                    end if
                    fcpsar = abs(fcpsar) !fraction volume of RV inside some cell

                    ! save fluctutations corresponding to max volume
                    if(fcpsar .gt. fcpsar_max(ix, jy, kz) .and. wadr .gt. 0) &
                        &then

                      ! here we should extend for terrain following coord.
                      ! new: wadr missing under the root
                      amprw = sqrt(abs(omir) * 2. * wnrh ** 2 / (wnrh ** 2 &
                          &+ wnrm ** 2) / rhoStratTFC(ix, jy, kz))

                      !interpolate phase at coarse cell center
                      dxx = x(ix + ix0) - xr
                      dyy = y(jy + jy0) - yr
                      dzz = zTFC(ix, jy, kz) - zr

                      !
                      !heightTFC(i, j, k)

                      dphi = (ray(iRay, ixrv, jyrv, kzrv)%dphi + wnrk * dxx &
                          &+ wnrl * dyy + wnrm * dzz)

                      !to be consisten with LES wavepacket simulation
                      !there amplitude of b11 is real with sign from vert. wavenumber
                      b11 = amprw / abs(omir / NNr) * sign(1., wnrm)
                      w10 = cmplx(0.0, omir / NNr) * b11 ! amplitude w

                      theta0 = thetaStratTFC(ix, jy, kz)
                      theta11 = Fr2 * theta0 * b11

                      pi12 = cmplx(0.0, kappa * Ma2 * (omir * omir - NNr) &
                          &/ NNr / wnrm / thetaStratTFC(ix, jy, kz)) * b11

                      thetaPrime = real(theta11 * exp(dphi * imag))
                      expPrime = real(pi12 * exp(dphi * imag))
                      wPrime = real(w10 * exp(dphi * imag))

                      !store fields
                      ray_varIce(ix, jy, kz)%wwp = wPrime
                      ray_varIce(ix, jy, kz)%epp = expPrime
                      ray_varIce(ix, jy, kz)%thp = thetaPrime
                      fcpsar_max(ix, jy, kz) = fcpsar !max area coverd by RV

                    end if

                  end do !ix
                end do !jy
              end do !kz

            else ! (no topography)

              !*********************
              ! local GW fluctuations
              !*********************

              if(average_cell) then
                kzmin = max(1, floor((zr - dzr * 0.5 - lz(0)) / dz) + 1)
                kzmax = min(nz, floor((zr + dzr * 0.5 - lz(0)) / dz) + 1)
              elseif(average_cell_3) then
                kzmin = kzrv - 1
                kzmax = kzrv + 1
              else
                kzmin = kzrv
                kzmax = kzrv
              end if

              if(reconstruct_gw_field == 1) then

                do kz = kzmin, kzmax
                  dzi = (min((zr + dzr * 0.5), lz(0) + kz * dz) - max((zr &
                      &- dzr * 0.5), lz(0) + (kz - 1) * dz))

                  fcpspz = dwnrm * dzi / dz

                  !SD
                  !fcpswn = dwnrm

                  do jy = jymin, jymax
                    if(sizeY > 1) then
                      dyi = (min((yr + dyr * 0.5), ly(0) + (jy + jy0) * dy) &
                          &- max((yr - dyr * 0.5), ly(0) + (jy + jy0 - 1) * dy))

                      fcpspy = dwnrl * dyi / dy

                      !SD
                      !fcpswn = fcpswn * dwnrl
                    else
                      fcpspy = 1.0
                    end if

                    do ix = ixmin, ixmax
                      if(sizeX > 1) then
                        dxi = (min((xr + dxr * 0.5), lx(0) + (ix + ix0) * dx) &
                            &- max((xr - dxr * 0.5), lx(0) + (ix + ix0 - 1) &
                            &* dx))

                        fcpspx = dwnrk * dxi / dx

                        !SD
                        !fcpswn = fcpswn * dwnrk
                        !?fcpsar = fcpsar * dxi /dx
                      else
                        fcpspx = 1.0
                      end if

                      wadr = fcpspx * fcpspy * fcpspz * ray(iRay, ixrv, jyrv, &
                          &kzrv)%dens

                      wadr_sum(ix, jy, kz) = wadr_sum(ix, jy, kz) + wadr !integrated waveaction density
                      fcpsar = dzi / dz
                      if(sizeX > 1) then
                        fcpsar = (dxi / dx) * fcpsar
                      end if
                      if(sizeY > 1) then
                        fcpsar = (dyi / dy) * fcpsar
                      end if
                      fcpsar = abs(fcpsar) !fraction volume of RV inside some cell

                      ! save fluctutations corresponding to max ray volume
                      if(fcpsar .gt. fcpsar_max(ix, jy, kz) .and. wadr .gt. 0) &
                          &then

                        ! here we should extend for terrain following coord.
                        ! new: wadr missing under the root
                        amprw = sqrt(abs(omir) * 2. * wnrh ** 2 / (wnrh ** 2 &
                            &+ wnrm ** 2) / rhoStrat(kz))

                        !interpolate phase at coarse cell center
                        dxx = x(ix + ix0) - xr
                        dyy = y(jy + jy0) - yr
                        dzz = z(kz) - zr

                        !
                        !heightTFC(i, j, k)

                        dphi = (ray(iRay, ixrv, jyrv, kzrv)%dphi + wnrk * dxx &
                            &+ wnrl * dyy + wnrm * dzz)

                        !to be consisten with LES wavepacket simulation
                        !there amplitude of b11 is real with sign from vert. wavenumber
                        b11 = amprw / abs(omir / NNr) * sign(1., wnrm)
                        w10 = cmplx(0.0, omir / NNr) * b11 ! amplitude w

                        theta0 = thetaStrat(kz)
                        theta11 = Fr2 * theta0 * b11

                        pi12 = cmplx(0.0, kappa * Ma2 * (omir * omir - NNr) &
                            &/ NNr / wnrm / thetaStrat(kz)) * b11

                        thetaPrime = real(theta11 * exp(dphi * imag))
                        expPrime = real(pi12 * exp(dphi * imag))
                        wPrime = real(w10 * exp(dphi * imag))

                        !store fields
                        ray_varIce(ix, jy, kz)%wwp = wPrime
                        ray_varIce(ix, jy, kz)%epp = expPrime
                        ray_varIce(ix, jy, kz)%thp = thetaPrime
                        fcpsar_max(ix, jy, kz) = fcpsar !max fraction volume in a cell

                        !TEST***
                        !compare wave action
                        var%OPT(ix, jy, kz, 3) = wPrime !dphi

                      end if

                    end do ! ix
                  end do !jy
                end do !kz

              elseif(reconstruct_gw_field == 2) then

                do kz = kzmin, kzmax

                  fcpswn = dwnrm

                  do jy = jymin, jymax

                    if(sizeY > 1) then
                      fcpswn = fcpswn * dwnrl
                    else
                      fcpspy = 1.0
                    end if

                    do ix = ixmin, ixmax

                      if(sizeX > 1) then
                        fcpswn = fcpswn * dwnrk
                      else
                        fcpspx = 1.0
                      end if

                      !interpolate phase at coarse cell center
                      dxx = x(ix + ix0) - xr
                      dyy = y(jy + jy0) - yr
                      dzz = z(kz) - zr

                      dphi = (ray(iRay, ixrv, jyrv, kzrv)%dphi + wnrk * dxx &
                          &+ wnrl * dyy + wnrm * dzz)

                      if(.not. compute_cloudcover) then
                        if(gauss_smoothing) then

                          ! standard deviation Gauss dist. set to width RV
                          ! (sigdeX is variance)
                          sigdex = dxr ** 2
                          sigdey = dyr ** 2
                          sigdez = dzr ** 2

                          nn = 1. / sqrt(2. * pi) ** 3 / sqrt(sigdex * sigdey &
                              &* sigdez)
                          vv = dxr * dyr * dzr

                          ! dens = F exp{}
                          dens_gauss = nn * vv * ray(iRay, ixrv, jyrv, &
                              &kzrv)%dens * exp(- (dxx ** 2 / sigdex + dyy &
                              &** 2 / sigdey + dzz ** 2 / sigdez) / 2.)
                          amprw = sqrt(abs(omir) * 2. * wnrh ** 2 / (wnrh ** 2 &
                              &+ wnrm ** 2) / rhoStrat(kz) * fcpswn &
                              &* dens_gauss)
                        else
                          print *, 'compute_cloudcover=false works only with &
                              &gauss_smoothing=true'
                          stop
                        end if
                        !to be consisten with LES wavepacket simulation
                        !there amplitude of b11 is real with sign from vert. wavenumber
                        b11 = amprw / abs(omir / NNr) * sign(1., wnrm)
                        w10 = cmplx(0.0, omir / NNr) * b11 ! amplitude w

                        theta0 = thetaStrat(kz)
                        theta11 = Fr2 * theta0 * b11

                        pi12 = cmplx(0.0, kappa * Ma2 * (omir * omir - NNr) &
                            &/ NNr / wnrm / thetaStrat(kz)) * b11

                        thetaPrime = real(theta11 * exp(dphi * imag))
                        expPrime = real(pi12 * exp(dphi * imag))
                        wPrime = real(w10 * exp(dphi * imag))

                        !superimpose fields
                        ray_varIce(ix, jy, kz)%wwp = ray_varIce(ix, jy, &
                            &kz)%wwp + wPrime
                        ray_varIce(ix, jy, kz)%epp = ray_varIce(ix, jy, &
                            &kz)%epp + expPrime
                        ray_varIce(ix, jy, kz)%thp = ray_varIce(ix, jy, &
                            &kz)%thp + thetaPrime
                      elseif(compute_cloudcover) then

                        ! standard deviation Gauss dist. set to width RV
                        ! (sigdeX is variance)
                        sigdex = dxr ** 2
                        sigdey = dyr ** 2
                        sigdez = dzr ** 2

                        nn = 1. / sqrt(2. * pi) ** 3 / sqrt(sigdex * sigdey &
                            &* sigdez)
                        vv = dxr * dyr * dzr

                        do ii = 1, NSCX
                          do jj = 1, NSCY

                            !subcell center
                            xsc = x(ix + ix0) - dx / 2. + (ii - 0.5) * dxsc
                            ysc = y(jy + jy0) - dy / 2. + (jj - 0.5) * dysc

                            !interpolate phase/amplitude RV at cell center
                            dxx = xsc - xr
                            dyy = ysc - yr
                            dzz = z(kz) - zr

                            !****************
                            !CHANGES
                            !****************

                            if(sizeX > 1) then
                              dxi = (min((xr + dxr * 0.5), xsc + dxsc * 0.5) &
                                  &- max((xr - dxr * 0.5), xsc - dxsc * 0.5))

                              fcpspx = dwnrk * dxi / dx
                            else
                              fcpspx = 1.0
                            end if
                            if(sizeY > 1) then
                              dyi = (min((yr + dyr * 0.5), ysc + dysc * 0.5) &
                                  &- max((yr - dyr * 0.5), ysc - dysc * 0.5))

                              fcpspy = dwnrl * dyi / dy
                            else
                              fcpspy = 1.0
                            end if

                            dzi = (min((zr + dzr * 0.5), z(kz) + dz * 0.5) &
                                &- max((zr - dzr * 0.5), z(kz) - dz * 0.5))

                            fcpspz = dwnrm * dzi / dz

                            fcpswn = fcpspz * fcpspy * fcpspx

                            !****************
                            !END CHANGES
                            !****************

                            if(gauss_smoothing) then

                              ! dens = N exp{-x**2 / (2 \sigma)}
                              dens_gauss = nn * vv * ray(iRay, ixrv, jyrv, &
                                  &kzrv)%dens * exp(- (dxx ** 2 / sigdex + dyy &
                                  &** 2 / sigdey + dzz ** 2 / sigdez) / 2.)
                              amprw = sqrt(abs(omir) * 2. * wnrh ** 2 / (wnrh &
                                  &** 2 + wnrm ** 2) / rhoStrat(kz) * fcpswn &
                                  &* dens_gauss)

                            else
                              !cell center inside RV
                              if(abs(dxx) .le. dxr / 2. .and. abs(dyy) .le. &
                                  &dyr / 2. .and. abs(dzz) .le. dzr / 2.) then

                                amprw = sqrt(abs(omir) * 2. * wnrh ** 2 &
                                    &/ (wnrh ** 2 + wnrm ** 2) / rhoStrat(kz) &
                                    &* fcpswn * ray(iRay, ixrv, jyrv, &
                                    &kzrv)%dens)
                              end if
                            end if

                            !to be consisten with LES wavepacket simulation
                            !there amplitude of b11 is real with sign from vert. wavenumber
                            b11 = amprw / abs(omir / NNr) * sign(1., wnrm)
                            w10 = cmplx(0.0, omir / NNr) * b11 ! amplitude w

                            theta0 = thetaStrat(kz)
                            theta11 = Fr2 * theta0 * b11

                            pi12 = cmplx(0.0, kappa * Ma2 * (omir * omir &
                                &- NNr) / NNr / wnrm / thetaStrat(kz)) * b11

                            ! phase at cell center
                            dphi = (ray(iRay, ixrv, jyrv, kzrv)%dphi + wnrk &
                                &* dxx + wnrl * dyy + wnrm * dzz)

                            thetaPrime = real(theta11 * exp(dphi * imag))
                            expPrime = real(pi12 * exp(dphi * imag))
                            wPrime = real(w10 * exp(dphi * imag))

                            print *, wPrime
                            stop

                            !superimpose fields
                            ray_cloud(ix, jy, kz, ii, jj)%wwp = ray_cloud(ix, &
                                &jy, kz, ii, jj)%wwp + wPrime
                            ray_cloud(ix, jy, kz, ii, jj)%epp = ray_cloud(ix, &
                                &jy, kz, ii, jj)%epp + expPrime
                            ray_cloud(ix, jy, kz, ii, jj)%thp = ray_cloud(ix, &
                                &jy, kz, ii, jj)%thp + thetaPrime
                          end do ! ii
                        end do ! jj
                      end if ! compute_cloudcover

                    end do ! ix
                  end do !jy
                end do !kz

              else
                print *, 'reconstruct_gw_field only 1 or 2'
                stop
                stop
              end if

            end if ! topography

          end do !iray
        end do !ixrv
      end do !jyrv
    end do !kzrv

    if(reconstruct_gw_field == 1) then
      !***************************************
      ! compute effective fluctuations
      !***************************************

      do kz = 1, nz
        do jy = 1, ny
          do ix = 1, nx

            wadr = wadr_sum(ix, jy, kz)

            if(wadr .ge. 0) then

              ray_varIce(ix, jy, kz)%wwp = ray_varIce(ix, jy, kz)%wwp &
                  &* sqrt(wadr)
              ray_varIce(ix, jy, kz)%epp = ray_varIce(ix, jy, kz)%epp &
                  &* sqrt(wadr)
              ray_varIce(ix, jy, kz)%thp = ray_varIce(ix, jy, kz)%thp &
                  &* sqrt(wadr)

            end if

          end do !ix
        end do !jy
      end do !kz

    end if ! reconstruct_gw_field == 1

  end subroutine calc_ice

  ! -----------------------------------------------------------------------------
  function leading_order_tracer_flux(f0nondim, omega, wnkk, wnll, wnmm, &
      &wadens, direction, ix, jy, kz, var) result(LOtracerflx)

    real :: f0nondim, omega, wnkk, wnll, wnmm, wadens
    integer :: ix, jy, kz
    character(len = 1) :: direction
    type(var_type), intent(in) :: var

    real :: LOtracerflx

    real :: tracerflxcoeff, dchidx, dchidy, dchidz

    tracerflxcoeff = - f0nondim / omega * wnmm / (wnkk ** 2. + wnll ** 2. &
        &+ wnmm ** 2.) * wadens

    call tracerderivative(x(ix), 1, y(jy), z(kz), var, dchidx)
    call tracerderivative(y(jy), 2, x(ix), z(kz), var, dchidy)
    call tracerderivative(z(kz), 3, x(ix), y(jy), var, dchidz)

    if(direction == 'x') then
      LOtracerflx = tracerflxcoeff * (wnmm * dchidy - wnll * dchidz)
    elseif(direction == 'y') then
      LOtracerflx = tracerflxcoeff * (wnkk * dchidz - wnmm * dchidx)
    elseif(direction == 'z') then
      LOtracerflx = tracerflxcoeff * (wnll * dchidx - wnkk * dchidy)
    else
      stop "wkb.f90: function leading_order_tracer_flx incorrect direction"
    end if

  end function leading_order_tracer_flux

  ! ---------------------------------------------------------------------

  subroutine setboundary_waveAmp(waveampfield)

    type(waveAmpType), dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz &
        &+ nbz), intent(inout) :: waveampfield

    integer :: k

    complex, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: &
        &arrayrepl

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: &
        &arrayreplre

    arrayreplre = waveampfield(:, :, :)%phase
    call setboundary_hor_wkb(arrayreplre)
    waveampfield(:, :, :)%phase = arrayreplre

    do k = 1, min(nz, nbz)
      waveampfield(:, :, nz + k)%phase = waveampfield(:, :, nz - k + 1)%phase
      waveampfield(:, :, - k + 1)%phase = waveampfield(:, :, 1)%phase
    end do

    arrayrepl = waveampfield(:, :, :)%lowamp%u
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%lowamp%u = arrayrepl
    arrayrepl = waveampfield(:, :, :)%lowamp%v
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%lowamp%v = arrayrepl
    arrayrepl = waveampfield(:, :, :)%lowamp%w
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%lowamp%w = arrayrepl
    arrayrepl = waveampfield(:, :, :)%lowamp%b
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%lowamp%b = arrayrepl
    arrayrepl = waveampfield(:, :, :)%lowamp%pi
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%lowamp%pi = arrayrepl
    arrayrepl = waveampfield(:, :, :)%lowamp%chi
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%lowamp%chi = arrayrepl

    arrayrepl = waveampfield(:, :, :)%lowamp_prevts%u
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%lowamp_prevts%u = arrayrepl
    arrayrepl = waveampfield(:, :, :)%lowamp_prevts%v
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%lowamp_prevts%v = arrayrepl
    arrayrepl = waveampfield(:, :, :)%lowamp_prevts%w
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%lowamp_prevts%w = arrayrepl
    arrayrepl = waveampfield(:, :, :)%lowamp_prevts%b
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%lowamp_prevts%b = arrayrepl
    arrayrepl = waveampfield(:, :, :)%lowamp_prevts%pi
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%lowamp_prevts%pi = arrayrepl
    arrayrepl = waveampfield(:, :, :)%lowamp_prevts%chi
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%lowamp_prevts%chi = arrayrepl

    arrayrepl = waveampfield(:, :, :)%rhsamp%u
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%rhsamp%u = arrayrepl
    arrayrepl = waveampfield(:, :, :)%rhsamp%v
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%rhsamp%v = arrayrepl
    arrayrepl = waveampfield(:, :, :)%rhsamp%w
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%rhsamp%w = arrayrepl
    arrayrepl = waveampfield(:, :, :)%rhsamp%b
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%rhsamp%b = arrayrepl
    arrayrepl = waveampfield(:, :, :)%rhsamp%pi
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%rhsamp%pi = arrayrepl
    arrayrepl = waveampfield(:, :, :)%rhsamp%chi
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%rhsamp%chi = arrayrepl

    arrayrepl = waveampfield(:, :, :)%nowamp%u
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%nowamp%u = arrayrepl
    arrayrepl = waveampfield(:, :, :)%nowamp%v
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%nowamp%v = arrayrepl
    arrayrepl = waveampfield(:, :, :)%nowamp%w
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%nowamp%w = arrayrepl
    arrayrepl = waveampfield(:, :, :)%nowamp%b
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%nowamp%b = arrayrepl
    arrayrepl = waveampfield(:, :, :)%nowamp%pi
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%nowamp%pi = arrayrepl
    arrayrepl = waveampfield(:, :, :)%nowamp%chi
    call setboundary_wkb_cmplx(arrayrepl)
    waveampfield(:, :, :)%nowamp%chi = arrayrepl

    return

  end subroutine setboundary_waveAmp

  ! ---------------------------------------------------------------------

  subroutine setboundary_wkb_cmplx(flxwkb)

    ! -------------------------------------------------------------------
    ! boundary conditions for WKB fluxes
    ! so far only periodic boundary conditions allowed in horizontal
    ! -------------------------------------------------------------------

    complex, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb

    call setboundary_hor_wkb_cmplx(flxwkb)
    call setboundary_vrt_wkb_cmplx(flxwkb)

    return

  end subroutine setboundary_wkb_cmplx

  ! ---------------------------------------------------------------------

  subroutine setboundary_hor_wkb_cmplx(flxwkb)

    ! -------------------------------------------------------------------
    ! horizontal boundary conditions for WKB fluxes
    ! so far only periodic boundary conditions allowed
    ! -------------------------------------------------------------------

    complex, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb

    select case(xBoundary)
    case("periodic")
      call setboundary_x_periodic_wkb_cmplx(flxwkb)
    case default
      stop "setboundary_hor_wkb: unknown case xBoundary"
    end select

    select case(yBoundary)
    case("periodic")
      call setboundary_y_periodic_wkb_cmplx(flxwkb)
    case default
      stop "setboundary_hor_wkb: unknown case xBoundary"
    end select

    return

  end subroutine setboundary_hor_wkb_cmplx

  ! ----------------------------------------------------------------------

  subroutine setboundary_x_periodic_wkb_cmplx(flxwkb)

    ! -------------------------------------------------------------------
    ! periodic boundary conditions in x for WKB fluxes
    ! -------------------------------------------------------------------

    complex, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: ix

    ! auxiliary fields
    complex, dimension(nbx, - nby:ny + nby, - nbz:nz + nbz) :: &
        &xSliceLeft_send, xSliceRight_send
    complex, dimension(nbx, - nby:ny + nby, - nbz:nz + nbz) :: &
        &xSliceLeft_recv, xSliceRight_recv

    if(idim > 1) then
      ! more than 1 cpu in x direction

      call mpi_cart_shift(comm, 0, 1, left, right, ierror)

      ! slice size
      sendcount = nbx * (ny + 2 * nby + 1) * (nz + 2 * nbz + 1)
      recvcount = sendcount

      ! read slice into contiguous array
      do ix = 1, nbx
        xSliceLeft_send(ix, :, :) = flxwkb(ix, :, :)
        xSliceRight_send(ix, :, :) = flxwkb(nx - nbx + ix, :, :)
      end do

      ! left -> right
      source = left
      dest = right
      tag = 100

      call mpi_sendrecv(xSliceRight_send(1, - nby, - nbz), sendcount, &
          &mpi_double_precision, dest, tag, xSliceLeft_recv(1, - nby, - nbz), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xSliceLeft_send(1, - nby, - nbz), sendcount, &
          &mpi_double_precision, dest, tag, xSliceRight_recv(1, - nby, - nbz), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_right, ierror)

      ! write auxiliary slice to flux field
      do ix = 1, nbx
        ! right halos
        flxwkb(nx + ix, :, :) = xSliceRight_recv(ix, :, :)

        ! left halos
        flxwkb(- nbx + ix, :, :) = xSliceLeft_recv(ix, :, :)
      end do
    else
      ! only 1 cpu in x direction

      do ix = 1, min(nx, nbx)
        flxwkb(nx + ix, :, :) = flxwkb(ix, :, :)
        flxwkb(- ix + 1, :, :) = flxwkb(nx - ix + 1, :, :)
      end do
    end if

  end subroutine setboundary_x_periodic_wkb_cmplx

  ! ----------------------------------------------------------------------

  subroutine setboundary_y_periodic_wkb_cmplx(flxwkb)

    ! -------------------------------------------------------------------
    ! periodic boundary conditions in y for WKB fluxes
    ! -------------------------------------------------------------------

    complex, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: jy

    ! auxiliary fields
    complex, dimension(- nbx:nx + nbx, nby, - nbz:nz + nbz) :: &
        &ySliceBack_send, ySliceForw_send
    complex, dimension(- nbx:nx + nbx, nby, - nbz:nz + nbz) :: &
        &ySliceBack_recv, ySliceForw_recv

    if(jdim > 1) then
      ! more than 1 cpu in y direction

      call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

      ! slice size
      sendcount = nby * (nx + 2 * nbx + 1) * (nz + 2 * nbz + 1)
      recvcount = sendcount

      ! read slice into contiguous array
      do jy = 1, nby
        ySliceBack_send(:, jy, :) = flxwkb(:, jy, :)
        ySliceForw_send(:, jy, :) = flxwkb(:, ny - nby + jy, :)
      end do

      ! back -> forw
      source = back
      dest = forw
      tag = 100

      call mpi_sendrecv(ySliceForw_send(- nbx, 1, - nbz), sendcount, &
          &mpi_double_precision, dest, tag, ySliceBack_recv(- nbx, 1, - nbz), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(ySliceBack_send(- nbx, 1, - nbz), sendcount, &
          &mpi_double_precision, dest, tag, ySliceForw_recv(- nbx, 1, - nbz), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_forw, ierror)

      ! write auxiliary slice to var field
      do jy = 1, nby
        ! right halos
        flxwkb(:, ny + jy, :) = ySliceForw_recv(:, jy, :)

        ! left halos
        flxwkb(:, - nby + jy, :) = ySliceBack_recv(:, jy, :)
      end do
    else
      ! only 1 cpu in y direction

      do jy = 1, min(ny, nby)
        flxwkb(:, ny + jy, :) = flxwkb(:, jy, :)
        flxwkb(:, - jy + 1, :) = flxwkb(:, ny - jy + 1, :)
      end do
    end if

  end subroutine setboundary_y_periodic_wkb_cmplx

  ! ---------------------------------------------------------------------

  subroutine setboundary_vrt_wkb_cmplx(flxwkb)

    ! -------------------------------------------------------------------
    ! vertical boundary conditions for WKB fluxes
    ! -------------------------------------------------------------------

    complex, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb

    select case(zBoundary)
    case("periodic")
      call setboundary_z_periodic_wkb_cmplx(flxwkb)
    case("solid_wall")
      call setboundary_z_solidwall_wkb_cmplx(flxwkb)
    case default
      stop "setboundary_hor_wkb: unknown case xBoundary"
    end select

    return

  end subroutine setboundary_vrt_wkb_cmplx

  ! ----------------------------------------------------------------------

  subroutine setboundary_z_periodic_wkb_cmplx(flxwkb)

    ! -------------------------------------------------------------------
    ! periodic boundary conditions in z for WKB fluxes
    ! -------------------------------------------------------------------

    complex, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb

    integer :: k

    do k = 1, min(nz, nbz)
      flxwkb(:, :, nz + k) = flxwkb(:, :, k)
      flxwkb(:, :, - k + 1) = flxwkb(:, :, nz - k + 1)
    end do

  end subroutine setboundary_z_periodic_wkb_cmplx

  ! ----------------------------------------------------------------------

  subroutine setboundary_z_solidwall_wkb_cmplx(flxwkb)

    complex, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: flxwkb

    integer :: k

    do k = 1, min(nz, nbz)
      flxwkb(:, :, nz + k) = - flxwkb(:, :, nz - k + 1)
      flxwkb(:, :, - k + 1) = - flxwkb(:, :, k)
    end do

  end subroutine setboundary_z_solidwall_wkb_cmplx

end module wkb_module
