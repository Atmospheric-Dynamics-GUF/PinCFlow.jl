module wkb_module

  !-----------------------------------------------------------
  ! This module couples the WKB model (MS-GWaM) to PincFloit
  !-----------------------------------------------------------

  use type_module
  use timeScheme_module
  use atmosphere_module
  use muscl_module

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

  ! FJFeb2023
  public :: setup_topography_wkb
  public :: update_topography_wkb

  !------------------------------
  !   private module variables
  !------------------------------

  real, dimension(:, :, :, :, :), allocatable :: dxRay, dkRay
  real, dimension(:, :, :, :, :), allocatable :: ddxRay

  integer, dimension(:), allocatable :: ix2_sfc, jy2_sfc, ik_sfc, jl_sfc, km_sfc
  integer, dimension(:, :, :), allocatable :: ir_sfc

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

  ! FJMar2023
  real, dimension(:, :), allocatable :: k_spectrum, l_spectrum, &
      topography_spectrum, final_topography_spectrum

  ! FJMar2023
  real, dimension(:, :), allocatable :: long_tree

  contains

  !-----------------------------------------------------------------------

  subroutine calc_meanFlow_effect(ray, var, force, ray_var3D)

    ! supplemements cell-centered volume forces by WKB force
    ! as well as the heating by entropy-flux convergence

    implicit none

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var

    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, 3), intent(inout) :: force

    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, 1:6), intent(inout) :: &
        ray_var3D

    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), &
        intent(inout) :: ray

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
    integer :: ix, jy, kz

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

            if(zr < lz(0)) then
              select case(zBoundary)
              case("periodic")
                zr = lz(1) + mod(zr - lz(0), lz(1) - lz(0))
              case("solid_wall")
                if(zr + 0.5 * dzr < lz(0)) cycle
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
                    in x must be periodic'
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
                      leftmost cpu so that  one should have either'
                  print *, lx(1) - dx, '= lx(1) - dx < xr < lx(1) =', lx(1), ' &
                      or'
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
                      rightmost cpu so that  one should have either'
                  print *, lx(0), '= lx(0) < xr < lx(0) + dx =', lx(0) + dx, ' &
                      or'
                  print *, lx(1), '= lx(1) < xr < lx(1) + dx =', lx(1) + dx
                  stop
                end if
              end if

              !if(xr < lx(0)) then
              !   select case (xBoundary)
              !      case ("periodic")
              !         xr = lx(1) + mod(xr - lx(0),lx(1) - lx(0))
              !      case default
              !         stop"calc_meanflow_effect: unknown case &
              !              & xBoundary"
              !   end select
              !  elseif (xr > lx(1)) then
              !   select case (xBoundary)
              !      case ("periodic")
              !         xr = lx(0) + mod(xr - lx(1),lx(1) - lx(0))
              !      case default
              !         stop"calc_meanflow_effect: unknown case &
              !              & xBoundary"
              !   end select
              !end if
            end if

            if(sizeY > 1) then
              if(yBoundary /= "periodic") then
                print *, 'ERROR in calc_meanflow_effect:  boundary conditions &
                    in y must be periodic'
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
                      cpu in y dir. so that  one should have either'
                  print *, ly(1) - dy, '= ly(1) - dy < yr < ly(1) =', ly(1), ' &
                      or'
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
                      dir. so that  one should have either'
                  print *, ly(0), '= ly(0) < yr < ly(0) + dy =', ly(0) + dy, ' &
                      or'
                  print *, ly(1), '= ly(1) < yr < ly(1) + dy =', ly(1) + dy
                  stop
                end if
              end if

              !if(yr < ly(0)) then
              !   select case (yBoundary)
              !      case ("periodic")
              !         yr = ly(1) + mod(yr - ly(0),ly(1) - ly(0))
              !      case default
              !         stop"calc_meanflow_effect: unknown case &
              !              & yBoundary"
              !   end select
              !  elseif (yr > ly(1)) then
              !   select case (yBoundary)
              !      case ("periodic")
              !         yr = ly(0) + mod(yr - ly(1),ly(1) - ly(0))
              !      case default
              !         stop"calc_meanflow_effect: unknown case &
              !              & yBoundary"
              !   end select
              !end if
            end if

            wnrk = ray(iRay, ixrv, jyrv, kzrv)%k
            wnrl = ray(iRay, ixrv, jyrv, kzrv)%l
            wnrm = ray(iRay, ixrv, jyrv, kzrv)%m

            dwnrk = ray(iRay, ixrv, jyrv, kzrv)%dkray
            dwnrl = ray(iRay, ixrv, jyrv, kzrv)%dlray
            dwnrm = ray(iRay, ixrv, jyrv, kzrv)%dmray

            wnrh = sqrt(wnrk ** 2 + wnrl ** 2)

            if(zr < lz(0) - dz) then
              print *, 'ERROR IN calc_meanflow_effect: RAY VOLUME', iRay, 'in &
                  cell', ixrv, jyrv, kzrv, 'TOO LOW'
              stop
            end if

            call stratification(zr, 1, NNr)

            omir = branchr * sqrt(NNr * wnrh ** 2 + f_cor_nd ** 2 * wnrm ** 2) &
                / sqrt(wnrh ** 2 + wnrm ** 2)

            cgirx = wnrk * (NNr - omir ** 2) / (omir * (wnrh ** 2 + wnrm ** 2))
            cgiry = wnrl * (NNr - omir ** 2) / (omir * (wnrh ** 2 + wnrm ** 2))

            cgirz = - wnrm * (omir ** 2 - f_cor_nd ** 2) / (omir * (wnrh ** 2 &
                + wnrm ** 2))

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
                    - dxr * 0.5 - lx(0)) / dx) + 1
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
                  * 0.5), lz(0) + (kz - 1) * dz))

              fcpspz = dwnrm * dzi / dz

              do jy = jymin, jymax
                if(sizeY > 1) then
                  dyi = (min((yr + dyr * 0.5), ly(0) + (jy + jy0) * dy) &
                      - max((yr - dyr * 0.5), ly(0) + (jy + jy0 - 1) * dy))

                  fcpspy = dwnrl * dyi / dy
                else
                  fcpspy = 1.0
                end if

                do ix = ixmin, ixmax
                  if(sizeX > 1) then
                    dxi = (min((xr + dxr * 0.5), lx(0) + (ix + ix0) * dx) &
                        - max((xr - dxr * 0.5), lx(0) + (ix + ix0 - 1) * dx))

                    fcpspx = dwnrk * dxi / dx
                  else
                    fcpspx = 1.0
                  end if

                  wadr = fcpspx * fcpspy * fcpspz * ray(iRay, ixrv, jyrv, &
                      kzrv)%dens

                  if(sizeX > 1) then
                    if(f_cor_nd /= 0.0) then
                      var_uu(ix, jy, kz) = var_uu(ix, jy, kz) + wadr * (wnrk &
                          * cgirx - (wnrk * cgirx + wnrl * cgiry) / (1.0 &
                          - (omir / f_cor_nd) ** 2))
                    else
                      var_uu(ix, jy, kz) = var_uu(ix, jy, kz) + wadr * wnrk &
                          * cgirx
                    end if
                  end if

                  if(sizeX > 1 .or. sizeY > 1) then
                    var_uv(ix, jy, kz) = var_uv(ix, jy, kz) + wadr * cgirx &
                        * wnrl
                  end if

                  var_uw(ix, jy, kz) = var_uw(ix, jy, kz) + wadr * wnrk &
                      * cgirz / (1.0 - (f_cor_nd / omir) ** 2)

                  if(sizeY > 1) then
                    if(f_cor_nd /= 0.0) then
                      var_vv(ix, jy, kz) = var_vv(ix, jy, kz) + wadr * (wnrl &
                          * cgiry - (wnrk * cgirx + wnrl * cgiry) / (1.0 &
                          - (omir / f_cor_nd) ** 2))
                    else
                      var_vv(ix, jy, kz) = var_vv(ix, jy, kz) + wadr * wnrl &
                          * cgiry
                    end if
                  end if

                  var_vw(ix, jy, kz) = var_vw(ix, jy, kz) + wadr * wnrl &
                      * cgirz / (1.0 - (f_cor_nd / omir) ** 2)

                  if(f_cor_nd /= 0.0) then
                    var_ETx(ix, jy, kz) = var_ETx(ix, jy, kz) + wadr &
                        * f_cor_nd ** 2 * NNr * wnrk * wnrm / (rhoStrat(kz) &
                        * g_ndim * omir * (wnrh ** 2 + wnrm ** 2))

                    var_ETy(ix, jy, kz) = var_ETy(ix, jy, kz) + wadr &
                        * f_cor_nd ** 2 * NNr * wnrl * wnrm / (rhoStrat(kz) &
                        * g_ndim * omir * (wnrh ** 2 + wnrm ** 2))
                  end if

                  var_E(ix, jy, kz) = var_E(ix, jy, kz) + wadr * omir
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    !! for output of u'w', E_w/rho, u, v, w:
    ! for output of u'w', E_w, u, v, w:
    do kz = 1, nz
      do jy = 1, ny
        do ix = 1, nx
          ray_var3D(ix, jy, kz, 4) = var_ETx(ix, jy, kz)
          !ray_var3D(ix,jy,kz,5) =  var_ETy(ix,jy,kz)
          ray_var3D(ix, jy, kz, 5) = var_uw(ix, jy, kz) ! output change by FDK
          ray_var3D(ix, jy, kz, 6) = var_E(ix, jy, kz)
        end do
      end do
    end do

    ! testb
    ! print*,'in calc_meanFlow_effect:'
    ! print*,'peak density-normalized energy density =',&
    !       & maxval(ray_var3D(:,:,:,6)) * uRef**2
    ! print*,'total density-normalized energy =',&
    !       & sum(ray_var3D(1:nx,1:ny,1:nz,6)) * uRef**2 * dx * dy * dz
    ! teste

    ! horizontal entropy fluxes

    if(f_cor_nd /= 0.0) then
      do kz = 1, nz
        var_ut(:, :, kz) = thetaStrat(kz) / f_cor_nd * var_ETy(:, :, kz)
        var_vt(:, :, kz) = - thetaStrat(kz) / f_cor_nd * var_ETx(:, :, kz)
      end do
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
          case("pseudo_incompressible")
            rhotot = 0.5 * (var(ix, jy, kz, 1) + var(ix + 1, jy, kz, 1))
            if(fluctuationMode) rhotot = rhotot + rhoStrat(kz)
          case default
            stop "volumeForce: unknown case model."
          end select

          ! forcing in x direction

          var_drudt(ix, jy, kz) = - rhotot / rhoStrat(kz) * (var_uw(ix, jy, kz &
              + 1) - var_uw(ix, jy, kz - 1)) / (2.0 * dz)

          if(sizeX > 1) then
            var_drudt(ix, jy, kz) = var_drudt(ix, jy, kz) - rhotot &
                / rhoStrat(kz) * (var_uu(ix + 1, jy, kz) - var_uu(ix - 1, jy, &
                kz)) / (2.0 * dx)
          end if

          if(sizeY > 1) then
            var_drudt(ix, jy, kz) = var_drudt(ix, jy, kz) - rhotot &
                / rhoStrat(kz) * (var_uv(ix, jy + 1, kz) - var_uv(ix, jy - 1, &
                kz)) / (2.0 * dy)
          end if

          var_drudt(ix, jy, kz) = var_drudt(ix, jy, kz) + rhotot * var_ETx(ix, &
              jy, kz)

          ! forcing in y direction

          select case(model)
          case("Boussinesq")
            rhotot = rho00

            ! Note by Junhong Wei (20170814): There may be
            ! something wrong with the case of Boussinesq here.
            ! The case of pseudo_incompressible should be fine.
          case("pseudo_incompressible")
            rhotot = 0.5 * (var(ix, jy, kz, 1) + var(ix, jy + 1, kz, 1))
            if(fluctuationMode) rhotot = rhotot + rhoStrat(kz)
          case default
            stop "volumeForce: unknown case model."
          end select

          var_drvdt(ix, jy, kz) = - rhotot / rhoStrat(kz) * (var_vw(ix, jy, kz &
              + 1) - var_vw(ix, jy, kz - 1)) / (2.0 * dz)

          if(sizeX > 1) then
            var_drvdt(ix, jy, kz) = var_drvdt(ix, jy, kz) - rhotot &
                / rhoStrat(kz) * (var_uv(ix + 1, jy, kz) - var_uv(ix - 1, jy, &
                kz)) / (2.0 * dx)
          end if

          if(sizeY > 1) then
            var_drvdt(ix, jy, kz) = var_drvdt(ix, jy, kz) - rhotot &
                / rhoStrat(kz) * (var_vv(ix, jy + 1, kz) - var_vv(ix, jy - 1, &
                kz)) / (2.0 * dy)
          end if

          var_drvdt(ix, jy, kz) = var_drvdt(ix, jy, kz) + rhotot * var_ETy(ix, &
              jy, kz)
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
            case("pseudo_incompressible")
              if(fluctuationMode) then
                rhotot = var(ix, jy, kz, 1) + rhoStrat(kz)
              else
                rhotot = var(ix, jy, kz, 1)
              end if
            case default
              stop "volumeForce: unknown case model."
            end select

            if(sizeX > 1) then
              var_drtdt(ix, jy, kz) = rhotot * (var_ut(ix + 1, jy, kz) &
                  - var_ut(ix - 1, jy, kz)) / (2.0 * dx)
            end if

            if(sizeY > 1) then
              var_drtdt(ix, jy, kz) = var_drtdt(ix, jy, kz) + rhotot &
                  * (var_vt(ix, jy + 1, kz) - var_vt(ix, jy - 1, kz)) / (2.0 &
                  * dy)
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
      if(z(kz) < lz(0) + zmin_wkb) cycle

      do jy = 1, ny
        do ix = 1, nx
          select case(model)
          case("Boussinesq")
            rhotot = rho00

            ! Note by Junhong Wei (20170814): There may be
            ! something wrong with the case of Boussinesq here.
            ! The case of pseudo_incompressible should be fine.
          case("pseudo_incompressible")
            rhotot = 0.5 * (var(ix, jy, kz, 1) + var(ix + 1, jy, kz, 1))
            if(fluctuationMode) rhotot = rhotot + rhoStrat(kz)
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
          case("pseudo_incompressible")
            rhotot = 0.5 * (var(ix, jy, kz, 1) + var(ix, jy + 1, kz, 1))
            if(fluctuationMode) rhotot = rhotot + rhoStrat(kz)
          case default
            stop "volumeForce: unknown case model."
          end select

          force(ix, jy, kz, 2) = force(ix, jy, kz, 2) + var_drvdt(ix, jy, kz)

          ! for output of mean-flow acceleration in y direction by GWs
          ray_var3D(ix, jy, kz, 2) = var_drvdt(ix, jy, kz) / rhotot
        end do
      end do
    end do

    call setboundary_frc_wkb(force(:, :, :, 1))
    call setboundary_frc_wkb(force(:, :, :, 2))

    var(:, :, :, 8) = 0.0

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
            case("pseudo_incompressible")
              if(fluctuationMode) then
                rhotot = var(ix, jy, kz, 1) + rhoStrat(kz)
              else
                rhotot = var(ix, jy, kz, 1)
              end if
            case default
              stop "volumeForce: unknown case model."
            end select

            var(ix, jy, kz, 8) = var_drtdt(ix, jy, kz)

            ! for output of mean-flow potential-temperature tendency
            ! by GWs
            ray_var3D(ix, jy, kz, 3) = - var_drtdt(ix, jy, kz) / rhotot
          end do
        end do
      end do
    end if

  end subroutine calc_meanFlow_effect

  !---------------------------------------------------------------------

  subroutine setup_wkb(ray, ray_var3D, var)

    !------------------------------------------------
    ! allocate ray field
    ! initialize position, wave vector and frequency
    ! for rays
    !------------------------------------------------

    implicit none

    ! argument list

    type(rayType), dimension(:, :, :, :), allocatable, intent(out) :: ray
    real, dimension(:, :, :, :), allocatable, intent(out) :: ray_var3D
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(in) :: var

    ! local variables
    integer :: allocstat

    real, allocatable :: omi_notop(:), omi_sfc(:, :), wnk_sfc(:, :), &
        wnl_sfc(:, :), wnm_sfc(:, :)
    real, allocatable :: fld_amp(:, :, :) ! 3D wave action distribution

    ! testb
    ! real, allocatable :: fld_ene(:,:,:)  ! 3D density-normalized
    !                                      ! energy-density field
    ! teste

    integer :: ix, jy, kz
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

    ! Mountain properties (FJJan2023)
    real :: mountainHeight_wkb, mountainWidth_wkb
    real :: x_center, y_center
    real :: k_mountain_wkb

    ! Long number (FJJan2023)
    real :: long

    ix0 = is + nbx - 1
    jy0 = js + nby - 1

    ! some non-dimensionalizations

    f_cor_nd = f_Coriolis_dim * tRef

    xr0 = xr0_dim / lRef
    yr0 = yr0_dim / lRef
    zr0 = zr0_dim / lRef

    sigwpx = sigwpx_dim / lRef
    sigwpy = sigwpx_dim / lRef
    sigwpz = sigwpz_dim / lRef

    xrmin = xrmin_dim / lRef
    xrmax = xrmax_dim / lRef

    yrmin = yrmin_dim / lRef
    yrmax = yrmax_dim / lRef

    zrmin = zrmin_dim / lRef
    zrmax = zrmax_dim / lRef

    ! FJJan2023
    mountainHeight_wkb = mountainHeight_wkb_dim / lRef
    mountainWidth_wkb = mountainWidth_wkb_dim / lRef

    ! FJJan2023
    x_center = 0.5 * (lx(1) + lx(0))
    y_center = 0.5 * (ly(1) + ly(0))

    ! FJJan2023
    k_mountain_wkb = pi / mountainWidth_wkb

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
      kzmin = max(1, int(floor((zrmin - lz(0)) / dz)) + 1)
      kzmax = min(sizeZ, int(floor((zrmax - lz(0)) / dz)) + 1)

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

    if(fac_dk_init == 0.0) then
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

    if(fac_dl_init == 0.0) then
      nyRay = 1
    else
      nyRay = nray_fac * nryl * nrl_init
    end if

    ! maximum # of r.v. allowed in a cell before r.v. are merged
    nray_max = nxRay * nyRay * nzRay

    ! work-space size per wavenumber direction chosen to be twice the
    ! maximum number of r.v. allowed before they are merged
    ! (should this turn out to be too small, change to thrice, quadruple,
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

    n_sfc = nrm_init
    if(nxRay > 1) n_sfc = n_sfc * nxRay / nray_fac
    if(nyRay > 1) n_sfc = n_sfc * nyRay / nray_fac

    ! field of ray volumes
    allocate(ray(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), stat = allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate ray"

    ! # of ray volumes per cell
    allocate(nRay(0:nx + 1, 0:ny + 1, - 1:nz + 2), stat = allocstat)
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

      allocate(ik_sfc(n_sfc), stat = allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate ik_sfc"

      allocate(jl_sfc(n_sfc), stat = allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate jl_sfc"

      allocate(km_sfc(n_sfc), stat = allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate km_sfc"

      ! FJMar2023
      allocate(k_spectrum(1:nx, 1:ny), stat = allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate k_spectrum"
      allocate(l_spectrum(1:nx, 1:ny), stat = allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate l_spectrum"
      allocate(topography_spectrum(1:nx, 1:ny), stat = allocstat)
      if(allocstat /= 0) stop "setup_wkb: could not allocate &
          topography_spectrum"
    end if

    ! position displacement increment
    allocate(dxRay(3, nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), stat &
        = allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate dxRay"

    ! wave vector increment
    allocate(dkRay(3, nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), stat &
        = allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate dkRay"

    ! ray-volume extent increment
    allocate(ddxRay(3, nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), stat &
        = allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate ddxRay"

    ! fields for data WKB output
    allocate(ray_var3D(0:nx + 1, 0:ny + 1, 0:nz + 1, 1:6), stat = allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate ray_var3D"

    ! needed for initialization of ray volumes:
    if(case_wkb == 3) then
      allocate(omi_sfc(1:nx, 1:ny))
      allocate(wnk_sfc(1:nx, 1:ny))
      allocate(wnl_sfc(1:nx, 1:ny))
      allocate(wnm_sfc(1:nx, 1:ny))
    else
      allocate(omi_notop(1:sizeZ))
    end if

    allocate(fld_amp(1:nx, 1:ny, 0:sizeZ)) ! 3D wave action field

    ! testb
    ! allocate (fld_ene(1:nx,1:ny,0:sizeZ)) ! 3D density-normalized
    !                                             ! energy-density field
    ! teste

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

    ! achatzc:
    ! a slight inconsistency below is that the stratification
    ! is calculated at the cell centers, while it is later on
    ! interpolated to the ray-volume positions
    ! this is no issue in the isothermal case

    if(case_wkb == 3) then
      ! intrinsic frequency and horizontal wave number mountain wave

      ! FJFeb2023
      call setup_topography_wkb

      ! FJMar2023
      call stratification(z(0), 1, NN_nd)

      do jy = 1, ny
        do ix = 1, nx
          ! FJMar2023
          wnrk_init = k_spectrum(ix, jy)
          wnrl_init = l_spectrum(ix, jy)
          wnrh_init = sqrt(wnrk_init ** 2.0 + wnrl_init ** 2.0)
          wnrm_init = 0.0

          omi_sfc(ix, jy) = - var(ix, jy, 1, 2) * wnrk_init - var(ix, jy, 1, &
              3) * wnrl_init

          ! choose correct sign of horizontal wavenumbers in order to
          ! be on the correct frequency branch

          if(omi_sfc(ix, jy) * branchr >= 0.0) then
            wnk_sfc(ix, jy) = wnrk_init
            wnl_sfc(ix, jy) = wnrl_init
          else
            omi_sfc(ix, jy) = - omi_sfc(ix, jy)

            wnk_sfc(ix, jy) = - wnrk_init
            wnl_sfc(ix, jy) = - wnrl_init
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
          fld_amp(ix, jy, 0) = 0.0
          wnrm = 0.0

          ! FJJan2023
          ! if ((sigwpx == 0.0 .or. abs(x(ix + ix0) - xr0) < sigwpx) .and. &
          !     (sigwpy == 0.0 .or. abs(y(jy + jy0) - yr0) < sigwpy)) then
          if(abs(omi_sfc(ix, jy)) <= f_cor_nd) then
            fld_amp(ix, jy, 0) = 0.0
            wnrm = 0.0
          elseif(abs(omi_sfc(ix, jy)) < sqrt(NN_nd)) then
            wnrm = - branchr * sqrt(wnrh_init ** 2 * (NN_nd - omi_sfc(ix, jy) &
                ** 2) / (omi_sfc(ix, jy) ** 2 - f_cor_nd ** 2))

            ! Displacement (FJFeb2023)
            displm = topography_spectrum(ix, jy)

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

            ! Long number scaling (FJJan2023)
            if(long_scaling) then
              ! Compute Long number.
              long = displm * sqrt(NN_nd / (var(ix, jy, 1, 2) ** 2.0 + var(ix, &
                  jy, 1, 3) ** 2.0))
              if(long > 0.0) then
                if(long_fit == 0) then
                  ! Apply scaling from decision tree model.
                  displm = displm * compute_long_scaling(long)
                else if(long_fit == 1) then
                  ! Apply exponential scaling.
                  displm = displm * along * (1.0 - exp(- blong / long))
                else if(long_fit == 2) then
                  ! Apply fractional scaling.
                  displm = displm * along / long / (1.0 + blong / long)
                else if(long_fit == 3) then
                  ! Apply linear scaling.
                  displm = min(displm, along * displm / long)
                end if
              end if
            end if

            ! surface wave-action density
            fld_amp(ix, jy, 0) = 0.5 * rhoStrat(0) * displm ** 2 * omi_sfc(ix, &
                jy) * (wnrh_init ** 2 + wnrm ** 2) / wnrh_init ** 2
          else
            fld_amp(ix, jy, 0) = 0.0
            wnrm = 0.0
          end if
          ! end if

          wnm_sfc(ix, jy) = wnrm
        end do
      end do
    else
      do kz = 1, sizeZ
        ! local squared Brunt-Vaisala frequency

        call stratification(z(kz), 1, NN_nd)

        ! intrinsic frequency

        omi_notop(kz) = branchr * sqrt((NN_nd * wnrh_init ** 2 + f_cor_nd ** 2 &
            * wnrm_init ** 2) / (wnrh_init ** 2 + wnrm_init ** 2))

        ! wave-action density

        do jy = 1, ny
          do ix = 1, nx
            fld_amp(ix, jy, kz) = (amp_wkb / wnrm_init) ** 2 * (wnrh_init ** 2 &
                + wnrm_init ** 2) / (2.0 * wnrh_init ** 2) * omi_notop(kz) &
                * rhoStrat(kz)

            if(case_wkb == 1) then
              fld_amp(ix, jy, kz) = fld_amp(ix, jy, kz) * exp(- ((z(kz) - zr0) &
                  / sigwpz) ** 2)

              if(sigwpx_dim > 0.0) then
                fld_amp(ix, jy, kz) = fld_amp(ix, jy, kz) * exp(- ((x(ix &
                    + ix0) - xr0) / sigwpx) ** 2)
              end if

              if(sigwpy_dim > 0.0) then
                fld_amp(ix, jy, kz) = fld_amp(ix, jy, kz) * exp(- ((y(jy &
                    + jy0) - yr0) / sigwpy) ** 2)
              end if
            elseif(case_wkb == 2) then
              if(abs(z(kz) - zr0) < sigwpz) then
                fld_amp(ix, jy, kz) = fld_amp(ix, jy, kz) * 0.5 * (1.0 &
                    + cos(pi * (z(kz) - zr0) / sigwpz))

                if(sigwpx > 0.0) then
                  if(abs(x(ix + ix0) - xr0) < sigwpx) then
                    fld_amp(ix, jy, kz) = fld_amp(ix, jy, kz) * 0.5 * (1.0 &
                        + cos(pi * (x(ix + ix0) - xr0) / sigwpx))
                  else
                    fld_amp(ix, jy, kz) = 0.0
                  end if
                end if

                if(sigwpy > 0.0) then
                  if(abs(y(jy + jy0) - yr0) < sigwpy) then
                    fld_amp(ix, jy, kz) = fld_amp(ix, jy, kz) * 0.5 * (1.0 &
                        + cos(pi * (y(jy + jy0) - yr0) / sigwpy))
                  else
                    fld_amp(ix, jy, kz) = 0.0
                  end if
                end if
              else
                fld_amp(ix, jy, kz) = 0.0
              end if
            end if
            ! testb
            ! fld_ene(ix,jy,kz) &
            ! = fld_amp(ix,jy,kz) * omi_notop(kz) / rhoStrat(kz)
            ! teste
          end do
        end do
      end do
    end if

    ! testb
    ! print*,'in setup_wkb:'
    ! print*,'peak density-normalized energy density =',&
    !       & maxval(fld_ene(:,:,:)) * uRef**2
    ! teste

    cgx_max = 0.0
    cgy_max = 0.0
    cgz_max = 0.0

    ! in mountain-wave case initialization of only one layer of ray
    ! volumes just below the bottom surface:

    if(case_wkb == 3) then
      kz2min = nrzl
    else
      kz2min = 1
    end if

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

            ! in x-k subspace, loop over all r.v. within one spatial cell
            do ix2 = 1, nrxl
              ! in x-k subspace, loop over all r.v. within the
              ! wave-number extent to be filled with r.v.
              do ik = 1, nrk_init

                ! likewise for y-l subspace
                do jy2 = 1, nryl
                  do jl = 1, nrl_init

                    ! likewise for z-m subspace
                    do kz2 = kz2min, nrzl
                      do km = 1, nrm_init
                        if(case_wkb == 3) then
                          ! pointers for surface ray volumes

                          i_sfc = i_sfc + 1

                          ix2_sfc(i_sfc) = ix2
                          jy2_sfc(i_sfc) = jy2

                          ik_sfc(i_sfc) = ik
                          jl_sfc(i_sfc) = jl
                          km_sfc(i_sfc) = km

                          ! only add ray volumes with non-zero
                          ! wave-action density
                          ! (thus excluding intrinsic frequencies
                          ! outside of the allowed range)
                          ! excluded cases indicated by negative
                          ! ray-volume index

                          if(fld_amp(ix, jy, kz) == 0.0) then
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
                            + (ix2 - 0.5) * dx / nrxl)

                        ray(iRay, ix, jy, kz)%y = (y(jy + jy0) - 0.5 * dy &
                            + (jy2 - 0.5) * dy / nryl)

                        ray(iRay, ix, jy, kz)%z = (z(kz) - 0.5 * dz + (kz2 &
                            - 0.5) * dz / nrzl)

                        xr = ray(iRay, ix, jy, kz)%x
                        yr = ray(iRay, ix, jy, kz)%y
                        zr = ray(iRay, ix, jy, kz)%z

                        ! local squared Brunt_Vaisala frequency

                        if(zr < lz(0) - dz) then
                          print *, 'ERROR IN setup_wkb: RAY VOLUME', iRay, &
                              'at', ix, jy, kz, 'TOO LOW'
                          stop
                        end if

                        call stratification(zr, 1, NNr)

                        ! ray-volume spatial extensions

                        ray(iRay, ix, jy, kz)%dxray = dx / nrxl
                        ray(iRay, ix, jy, kz)%dyray = dy / nryl
                        ray(iRay, ix, jy, kz)%dzray = dz / nrzl

                        ! ray-volume wave numbers

                        if(case_wkb == 3) then
                          wnk_0 = wnk_sfc(ix, jy)
                          wnl_0 = wnl_sfc(ix, jy)
                          wnm_0 = wnm_sfc(ix, jy)
                        else
                          wnk_0 = wnrk_init
                          wnl_0 = wnrl_init
                          wnm_0 = wnrm_init
                        end if

                        ray(iRay, ix, jy, kz)%k = (wnk_0 - 0.5 * dk_ini_nd &
                            + (real(ik) - 0.5) * dk_ini_nd / nrk_init)

                        ray(iRay, ix, jy, kz)%l = (wnl_0 - 0.5 * dl_ini_nd &
                            + (real(jl) - 0.5) * dl_ini_nd / nrl_init)

                        if(fac_dm_init == 0.0) then
                          stop 'ERROR: FAC_DM_INIT = 0.0'
                        else if(wnm_0 == 0.0) then
                          stop 'ERROR: WNM_0 = 0.0'
                        else
                          dm_ini_nd = fac_dm_init * abs(wnm_0)
                        end if

                        ray(iRay, ix, jy, kz)%m = (wnm_0 - 0.5 * dm_ini_nd &
                            + (real(km) - 0.5) * dm_ini_nd / nrm_init)

                        ! ray-volume wave-number extents

                        ray(iRay, ix, jy, kz)%dkray = dk_ini_nd / nrk_init

                        ray(iRay, ix, jy, kz)%dlray = dl_ini_nd / nrl_init

                        ray(iRay, ix, jy, kz)%dmray = dm_ini_nd / nrm_init

                        ! ray-volume phase-space volume

                        ray(iRay, ix, jy, kz)%area_xk = ray(iRay, ix, jy, &
                            kz)%dxray * ray(iRay, ix, jy, kz)%dkray

                        ray(iRay, ix, jy, kz)%area_yl = ray(iRay, ix, jy, &
                            kz)%dyray * ray(iRay, ix, jy, kz)%dlray

                        ray(iRay, ix, jy, kz)%area_zm = ray(iRay, ix, jy, &
                            kz)%dzray * ray(iRay, ix, jy, kz)%dmray

                        pspvol = dm_ini_nd

                        if(fac_dk_init /= 0.0) then
                          pspvol = pspvol * dk_ini_nd
                        end if

                        if(fac_dl_init /= 0.0) then
                          pspvol = pspvol * dl_ini_nd
                        end if

                        ! phase-space wave-action density

                        if(kz == sizeZ) then
                          ray(iRay, ix, jy, kz)%dens = 0.0
                        else
                          ray(iRay, ix, jy, kz)%dens = fld_amp(ix, jy, kz) &
                              / pspvol
                        endif

                        ! intrinsic frequency

                        if(case_wkb == 3) then
                          ray(iRay, ix, jy, kz)%omega = omi_sfc(ix, jy)
                        else
                          ray(iRay, ix, jy, kz)%omega = omi_notop(kz)
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

                        cgirx = wnrk * (NNr - omir ** 2) / (omir * (wnrh ** 2 &
                            + wnrm ** 2))

                        if(abs(uxr + cgirx) > abs(cgx_max)) then
                          cgx_max = abs(uxr + cgirx)
                        end if

                        cgiry = wnrl * (NNr - omir ** 2) / (omir * (wnrh ** 2 &
                            + wnrm ** 2))

                        if(abs(vyr + cgiry) > abs(cgy_max)) then
                          cgy_max = abs(vyr + cgiry)
                        end if

                        cgirz = - wnrm * (omir ** 2 - f_cor_nd ** 2) / (omir &
                            * (wnrh ** 2 + wnrm ** 2))

                        if(abs(wzr + cgirz) > abs(cgz_max)) then
                          cgz_max = abs(wzr + cgirz)
                        end if
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
            / (sizeX * sizeY * (sizeZ + 1))
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
    elseif(strtpe == 2) then
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
      factor = (zu - zlc) / dz
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
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(in) :: var
    integer, intent(in) :: flwtpe
    real, intent(out) :: flw

    integer :: kzu, kzd, jyf, jyb, ixr, ixl

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

    !! implement horizontal boundary conditions for ray-volume positions

    !if (sizeX > 1) then
    !   if(xlc < lx(0)) then
    !      select case (xBoundary)
    !         case ("periodic")
    !            xlc = lx(1) + mod(xlc - lx(0),lx(1) - lx(0))
    !            ! testb
    !              if (xlc < lx(0)) then
    !                 print*,'ERROR IN MEANFLOW: xlc =',xlc,'< lx(0) =', &
    !                         & lx(0),'with lx(1) =',lx(1)
    !                 print*,'lRef =',lRef
    !                 stop
    !              end if
    !            ! teste
    !         case default
    !            stop"calc_meanflow_effect: unknown case xBoundary"
    !      end select
    !     elseif (xlc > lx(1)) then
    !      select case (xBoundary)
    !         case ("periodic")
    !            xlc = lx(0) + mod(xlc - lx(1),lx(1) - lx(0))
    !            ! testb
    !              if (xlc > lx(1)) then
    !                 print*,'ERROR IN MEANFLOW: xlc =',xlc,'> lx(1) =', &
    !                         & lx(1),'with lx(0) =',lx(0)
    !                 print*,'lRef =',lRef
    !                 stop
    !              end if
    !            ! teste
    !         case default
    !            stop"calc_meanflow_effect: unknown case xBoundary"
    !      end select
    !   end if
    !end if

    !if (sizeY > 1) then
    !   if(ylc < ly(0)) then
    !      select case (yBoundary)
    !         case ("periodic")
    !            ylc = ly(1) + mod(ylc - ly(0),ly(1) - ly(0))
    !         case default
    !            stop"calc_meanflow_effect: unknown case yBoundary"
    !      end select
    !     elseif (ylc > ly(1)) then
    !      select case (yBoundary)
    !         case ("periodic")
    !            ylc = ly(0) + mod(ylc - ly(1),ly(1) - ly(0))
    !         case default
    !            stop"calc_meanflow_effect: unknown case yBoundary"
    !      end select
    !   end if
    !end if

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
      flwlbd = var(ixl, jyb, kzd, 2)
      flwlbu = var(ixl, jyb, kzu, 2)

      flwlfd = var(ixl, jyf, kzd, 2)
      flwlfu = var(ixl, jyf, kzu, 2)

      flwrbd = var(ixr, jyb, kzd, 2)
      flwrbu = var(ixr, jyb, kzu, 2)

      flwrfd = var(ixr, jyf, kzd, 2)
      flwrfu = var(ixr, jyf, kzu, 2)
    elseif(flwtpe == 2) then
      ! interpolate v using staggered-grid distribution

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
      flwlbd = var(ixl, jyb, kzd, 3)
      flwlbu = var(ixl, jyb, kzu, 3)

      flwlfd = var(ixl, jyf, kzd, 3)
      flwlfu = var(ixl, jyf, kzu, 3)

      flwrbd = var(ixr, jyb, kzd, 3)
      flwrbu = var(ixr, jyb, kzu, 3)

      flwrfd = var(ixr, jyf, kzd, 3)
      flwrfu = var(ixr, jyf, kzu, 3)
    elseif(flwtpe == 3) then
      ! interpolate w using staggered-grid distribution

      ! for variable at intermediate levels
      ! z(k+1/2) = lz(0) + k*dz:

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
        flwlbu = var(ixl, jyb, kzu, 4)

        flwlfd = 0.0
        flwlfu = var(ixl, jyf, kzu, 4)

        flwrbd = 0.0
        flwrbu = var(ixr, jyb, kzu, 4)

        flwrfd = 0.0
        flwrfu = var(ixr, jyf, kzu, 4)
      else
        flwlbd = var(ixl, jyb, kzd, 4)
        flwlbu = var(ixl, jyb, kzu, 4)

        flwlfd = var(ixl, jyf, kzd, 4)
        flwlfu = var(ixl, jyf, kzu, 4)

        flwrbd = var(ixr, jyb, kzd, 4)
        flwrbu = var(ixr, jyb, kzu, 4)

        flwrfd = var(ixr, jyf, kzd, 4)
        flwrfu = var(ixr, jyf, kzu, 4)
      end if
    elseif(flwtpe == 11) then
      ! interpolate du/dx using staggered-grid distribution
      ! du/dx (i,j,k)) = (u(i+1/2,j,k) - u(i-1/2,j,k))/dx, hence

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
      flwlbd = (var(ixl, jyb, kzd, 2) - var(ixl - 1, jyb, kzd, 2)) / dx
      flwlbu = (var(ixl, jyb, kzu, 2) - var(ixl - 1, jyb, kzu, 2)) / dx

      flwlfd = (var(ixl, jyf, kzd, 2) - var(ixl - 1, jyf, kzd, 2)) / dx
      flwlfu = (var(ixl, jyf, kzu, 2) - var(ixl - 1, jyf, kzu, 2)) / dx

      flwrbd = (var(ixr, jyb, kzd, 2) - var(ixr - 1, jyb, kzd, 2)) / dx
      flwrbu = (var(ixr, jyb, kzu, 2) - var(ixr - 1, jyb, kzu, 2)) / dx

      flwrfd = (var(ixr, jyf, kzd, 2) - var(ixr - 1, jyf, kzd, 2)) / dx
      flwrfu = (var(ixr, jyf, kzu, 2) - var(ixr - 1, jyf, kzu, 2)) / dx
    elseif(flwtpe == 12) then
      ! interpolate du/dy using staggered-grid distribution
      ! du/dy (i+1/2,j+1/2,k)) = (u(i+1/2,j+1,k) - u(i+1/2,j,k))/dy, hence

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
      flwlbd = (var(ixl, jyb + 1, kzd, 2) - var(ixl, jyb, kzd, 2)) / dy
      flwlbu = (var(ixl, jyb + 1, kzu, 2) - var(ixl, jyb, kzu, 2)) / dy

      flwlfd = (var(ixl, jyf + 1, kzd, 2) - var(ixl, jyf, kzd, 2)) / dy
      flwlfu = (var(ixl, jyf + 1, kzu, 2) - var(ixl, jyf, kzu, 2)) / dy

      flwrbd = (var(ixr, jyb + 1, kzd, 2) - var(ixr, jyb, kzd, 2)) / dy
      flwrbu = (var(ixr, jyb + 1, kzu, 2) - var(ixr, jyb, kzu, 2)) / dy

      flwrfd = (var(ixr, jyf + 1, kzd, 2) - var(ixr, jyf, kzd, 2)) / dy
      flwrfu = (var(ixr, jyf + 1, kzu, 2) - var(ixr, jyf, kzu, 2)) / dy
    elseif(flwtpe == 13) then
      ! interpolate du/dz using staggered-grid distribution
      ! du/dz (i+1/2,j,k+1/2)) = (u(i+1/2,j,k+1) - u(i+1/2,j,k))/dz, hence

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
        flwlbu = (var(ixl, jyb, kzu + 1, 2) - var(ixl, jyb, kzu, 2)) / dz

        flwlfd = 0.0
        flwlfu = (var(ixl, jyf, kzu + 1, 2) - var(ixl, jyf, kzu, 2)) / dz

        flwrbd = 0.0
        flwrbu = (var(ixr, jyb, kzu + 1, 2) - var(ixr, jyb, kzu, 2)) / dz

        flwrfd = 0.0
        flwrfu = (var(ixr, jyf, kzu + 1, 2) - var(ixr, jyf, kzu, 2)) / dz
      else
        if(zu < lz(1)) then
          flwlbd = (var(ixl, jyb, kzd + 1, 2) - var(ixl, jyb, kzd, 2)) / dz
          flwlbu = (var(ixl, jyb, kzu + 1, 2) - var(ixl, jyb, kzu, 2)) / dz

          flwlfd = (var(ixl, jyf, kzd + 1, 2) - var(ixl, jyf, kzd, 2)) / dz
          flwlfu = (var(ixl, jyf, kzu + 1, 2) - var(ixl, jyf, kzu, 2)) / dz

          flwrbd = (var(ixr, jyb, kzd + 1, 2) - var(ixr, jyb, kzd, 2)) / dz
          flwrbu = (var(ixr, jyb, kzu + 1, 2) - var(ixr, jyb, kzu, 2)) / dz

          flwrfd = (var(ixr, jyf, kzd + 1, 2) - var(ixr, jyf, kzd, 2)) / dz
          flwrfu = (var(ixr, jyf, kzu + 1, 2) - var(ixr, jyf, kzu, 2)) / dz
        elseif(zd < lz(1)) then
          flwlbd = (var(ixl, jyb, kzd + 1, 2) - var(ixl, jyb, kzd, 2)) / dz
          flwlbu = 0.0

          flwlfd = (var(ixl, jyf, kzd + 1, 2) - var(ixl, jyf, kzd, 2)) / dz
          flwlfu = 0.0

          flwrbd = (var(ixr, jyb, kzd + 1, 2) - var(ixr, jyb, kzd, 2)) / dz
          flwrbu = 0.0

          flwrfd = (var(ixr, jyf, kzd + 1, 2) - var(ixr, jyf, kzd, 2)) / dz
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
    elseif(flwtpe == 21) then
      ! interpolate dv/dx using staggered-grid distribution
      ! dv/dx(i+1/2,j+1/2,k)) = (v(i+1,j+1/2,k) - v(i,j+1/2,k))/dx, hence

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
      flwlbd = (var(ixl + 1, jyb, kzd, 3) - var(ixl, jyb, kzd, 3)) / dx
      flwlbu = (var(ixl + 1, jyb, kzu, 3) - var(ixl, jyb, kzu, 3)) / dx

      flwlfd = (var(ixl + 1, jyf, kzd, 3) - var(ixl, jyf, kzd, 3)) / dx
      flwlfu = (var(ixl + 1, jyf, kzu, 3) - var(ixl, jyf, kzu, 3)) / dx

      flwrbd = (var(ixr + 1, jyb, kzd, 3) - var(ixr, jyb, kzd, 3)) / dx
      flwrbu = (var(ixr + 1, jyb, kzu, 3) - var(ixr, jyb, kzu, 3)) / dx

      flwrfd = (var(ixr + 1, jyf, kzd, 3) - var(ixr, jyf, kzd, 3)) / dx
      flwrfu = (var(ixr + 1, jyf, kzu, 3) - var(ixr, jyf, kzu, 3)) / dx
    elseif(flwtpe == 22) then
      ! interpolate dv/dy using staggered-grid distribution
      ! dv/dy (i,j,k)) = (v(i,j+1/2,k) - v(i,j-1/2,k))/dy, hence

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
      flwlbd = (var(ixl, jyb, kzd, 3) - var(ixl, jyb - 1, kzd, 3)) / dy
      flwlbu = (var(ixl, jyb, kzu, 3) - var(ixl, jyb - 1, kzu, 3)) / dy

      flwlfd = (var(ixl, jyf, kzd, 3) - var(ixl, jyf - 1, kzd, 3)) / dy
      flwlfu = (var(ixl, jyf, kzu, 3) - var(ixl, jyf - 1, kzu, 3)) / dy

      flwrbd = (var(ixr, jyb, kzd, 3) - var(ixr, jyb - 1, kzd, 3)) / dy
      flwrbu = (var(ixr, jyb, kzu, 3) - var(ixr, jyb - 1, kzu, 3)) / dy

      flwrfd = (var(ixr, jyf, kzd, 3) - var(ixr, jyf - 1, kzd, 3)) / dy
      flwrfu = (var(ixr, jyf, kzu, 3) - var(ixr, jyf - 1, kzu, 3)) / dy
    elseif(flwtpe == 23) then
      ! interpolate dv/dz using staggered-grid distribution
      ! dv/dz (i,j+1/2,k+1/2)) = (v(i,j+1/2,k+1) - v(i,j+1/2,k))/dz, hence

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
        flwlbu = (var(ixl, jyb, kzu + 1, 3) - var(ixl, jyb, kzu, 3)) / dz

        flwlfd = 0.0
        flwlfu = (var(ixl, jyf, kzu + 1, 3) - var(ixl, jyf, kzu, 3)) / dz

        flwrbd = 0.0
        flwrbu = (var(ixr, jyb, kzu + 1, 3) - var(ixr, jyb, kzu, 3)) / dz

        flwrfd = 0.0
        flwrfu = (var(ixr, jyf, kzu + 1, 3) - var(ixr, jyf, kzu, 3)) / dz
      else
        if(zu < lz(1)) then
          flwlbd = (var(ixl, jyb, kzd + 1, 3) - var(ixl, jyb, kzd, 3)) / dz
          flwlbu = (var(ixl, jyb, kzu + 1, 3) - var(ixl, jyb, kzu, 3)) / dz

          flwlfd = (var(ixl, jyf, kzd + 1, 3) - var(ixl, jyf, kzd, 3)) / dz
          flwlfu = (var(ixl, jyf, kzu + 1, 3) - var(ixl, jyf, kzu, 3)) / dz

          flwrbd = (var(ixr, jyb, kzd + 1, 3) - var(ixr, jyb, kzd, 3)) / dz
          flwrbu = (var(ixr, jyb, kzu + 1, 3) - var(ixr, jyb, kzu, 3)) / dz

          flwrfd = (var(ixr, jyf, kzd + 1, 3) - var(ixr, jyf, kzd, 3)) / dz
          flwrfu = (var(ixr, jyf, kzu + 1, 3) - var(ixr, jyf, kzu, 3)) / dz
        elseif(zd < lz(1)) then
          flwlbd = (var(ixl, jyb, kzd + 1, 3) - var(ixl, jyb, kzd, 3)) / dz
          flwlbu = 0.0

          flwlfd = (var(ixl, jyf, kzd + 1, 3) - var(ixl, jyf, kzd, 3)) / dz
          flwlfu = 0.0

          flwrbd = (var(ixr, jyb, kzd + 1, 3) - var(ixr, jyb, kzd, 3)) / dz
          flwrbu = 0.0

          flwrfd = (var(ixr, jyf, kzd + 1, 3) - var(ixr, jyf, kzd, 3)) / dz
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
    end if

    ! interpolation in y

    if(sizeY == 1) then
      flwd = flwbd
      flwu = flwbu
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
      factor = (zu - zlc) / dz
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
    !
    ! only done if r.v. is completely above the ground
    ! (in order to avoid inconsisteny with relaunches)
    !-----------------------------------------------------------------

    implicit none

    ! argument list
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), &
        intent(inout) :: ray

    integer nrlc
    integer ix, jy, kz

    real :: xr, yr, zr
    real :: dxr, dyr, dzr
    real :: axk, ayl, azm

    integer :: nrvtt0, nrvtt1, nrvloc

    ! total number of ray volumes before splitting

    nrvloc = sum(nRay(1:nx, 1:ny, 0:nz))

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
              dxr = ray(iRay, ix, jy, kz)%dxray

              zr = ray(iRay, ix, jy, kz)%z

              dzr = ray(iRay, ix, jy, kz)%dzray

              if(dxr > dx .and. zr - 0.5 * dzr > lz(0)) then
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
              dyr = ray(iRay, ix, jy, kz)%dyray

              zr = ray(iRay, ix, jy, kz)%z

              dzr = ray(iRay, ix, jy, kz)%dzray

              if(dyr > dy .and. zr - 0.5 * dzr > lz(0)) then
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
            dzr = ray(iRay, ix, jy, kz)%dzray

            zr = ray(iRay, ix, jy, kz)%z

            if(dzr > dz .and. zr - 0.5 * dzr > lz(0)) then
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

    ! total number of ray volumes before after splitting

    nrvloc = sum(nRay(1:nx, 1:ny, 0:nz))

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
    !
    ! assumed boundary conditions: periodic in x and y
    !                              solid (absorbing) in z
    !-----------------------------------------------------------------

    implicit none

    ! argument list
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), &
        intent(inout) :: ray

    integer :: ix, jy, kz
    integer :: nrlc
    integer :: jRay
    integer :: nsl, nsr, nsb, nsf, nsd, nsu
    integer :: ish, jsh, ksh, nsh
    integer, dimension(:, :, :), allocatable :: nshl, nshr, nshb, nshf, nshd, &
        nshu
    integer, dimension(:, :, :, :), allocatable :: irsl, irsr, irsb, irsf, &
        irsd, irsu

    integer :: irsh(nray_wrk)

    logical :: lplace

    real :: xr, yr, zr

    integer :: ix0, jy0

    allocate(nshl(0:nx + 1, 0:ny + 1, - 1:nz + 2))
    allocate(nshr(0:nx + 1, 0:ny + 1, - 1:nz + 2))
    allocate(nshb(0:nx + 1, 0:ny + 1, - 1:nz + 2))
    allocate(nshf(0:nx + 1, 0:ny + 1, - 1:nz + 2))
    allocate(nshd(0:nx + 1, 0:ny + 1, - 1:nz + 2))
    allocate(nshu(0:nx + 1, 0:ny + 1, - 1:nz + 2))

    allocate(irsl(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2))
    allocate(irsr(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2))
    allocate(irsb(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2))
    allocate(irsf(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2))
    allocate(irsd(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2))
    allocate(irsu(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2))

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

      do kz = - 1, nz + 2
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

      do kz = - 1, nz + 2
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
                        not be the case'
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
                  - 1, jy, kz)) then
                print *, 'at ix,jy,kz =', ix, jy, kz, 'nrlc &
                    /= nRay(ix,jy,kz)   - nshl(ix+1,jy,kz) - nshr(ix-1,jy,kz)'
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

      do kz = - 1, nz + 2
        do jy = 1, ny
          do ix = 0, nx + 1
            nrlc = nRay(ix, jy, kz)

            nsb = 0
            nsf = 0

            ! ray volumes from cell behind

            if(nRay(ix, jy - 1, kz) > 0) then
              do iRay = 1, nRay(ix, jy - 1, kz)
                yr = ray(iRay, ix, jy - 1, kz)%y

                if(yr > y(jy - 1) + 0.5 * dy) then
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

                if(yr < y(jy + 1) - 0.5 * dy) then
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

      do kz = - 1, nz + 2
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
                        not be the case'
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
                  - 1, kz)) then
                print *, 'at ix,jy,kz =', ix, jy, kz, 'nrlc &
                    /= nRay(ix,jy,kz)   - nshb(ix,jy+1,kz) - nshf(ix,jy-1,kz)'
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

      ! boundary conditions in z:
      ! open boundary at model top: remove all rays above the domain

      do kz = nz + 1, nz + 2
        nRay(:, :, kz) = 0
      end do

      ! also remove all r.v. from layers below the first layer below the
      ! bottom

      if(nbz > 0) then
        do kz = - 1, - 1
          nRay(:, :, kz) = 0
        end do
      end if

      ! move ray volumes to appropriate cell

      do kz = - 1, nz + 1
        do jy = 0, ny + 1
          do ix = 0, nx + 1
            nrlc = nRay(ix, jy, kz)

            nsd = 0
            nsu = 0

            ! ray volumes from cell below
            ! only transfer into cells above the model bottom
            ! (so that, e.g., kz = 0 does not receive anything from
            ! below)

            if(kz > 0) then
              if(nRay(ix, jy, kz - 1) > 0) then
                do iRay = 1, nRay(ix, jy, kz - 1)
                  zr = ray(iRay, ix, jy, kz - 1)%z

                  if(zr > z(kz - 1) + 0.5 * dz) then
                    nrlc = nrlc + 1

                    ray(nrlc, ix, jy, kz) = ray(iRay, ix, jy, kz - 1)

                    nsd = nsd + 1

                    irsd(nsd, ix, jy, kz) = iRay
                  end if
                end do

                nshd(ix, jy, kz) = nsd
              end if
            end if

            ! ray volumes from cell above
            ! only transfer into cells below the uppermost layer
            ! within the model domain
            ! (so that, e.g., kz = nz does not receive anything from
            ! above)

            if(kz < nz + 1) then
              if(nRay(ix, jy, kz + 1) > 0) then
                do iRay = 1, nRay(ix, jy, kz + 1)
                  zr = ray(iRay, ix, jy, kz + 1)%z

                  if(zr < z(kz + 1) - 0.5 * dz) then
                    ! r.v. having propagated into layers below the
                    ! model bottom are tagged to be deleted
                    ! below, but they are not transferred
                    ! (a reflecting boundary condition would be
                    ! more physical ...)
                    if(kz > 0) then
                      nrlc = nrlc + 1
                      ray(nrlc, ix, jy, kz) = ray(iRay, ix, jy, kz + 1)
                    end if

                    nsu = nsu + 1

                    irsu(nsu, ix, jy, kz) = iRay
                  end if
                end do

                nshu(ix, jy, kz) = nsu
              end if
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

      ! remove ray volumes from cells they have left

      do kz = 0, nz
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
                        not be the case'
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
                  jy, kz - 1)) then
                print *, 'at ix,jy,kz =', ix, jy, kz, 'nrlc &
                    /= nRay(ix,jy,kz)   - nshd(ix,jy,kz+1) - nshu(ix,jy,kz-1)'
                stop
              else
                nRay(ix, jy, kz) = nrlc
              end if
            end if
          end do ! ix
        end do ! jy
      end do ! kz
    end if

    ! testb
    do kz = - 1, nz + 2
      do jy = 0, ny + 1
        do ix = 0, nx + 1
          if(nRay(ix, jy, kz) > 0) then
            do iRay = 1, nRay(ix, jy, kz)
              if(sizeX > 1) then
                xr = ray(iRay, ix, jy, kz)%x

                if(xr < x(ix + ix0) - 0.5 * dx) then
                  print *, 'ERROR in shift_rayvol:'
                  print *, 'xr =', xr, '< x(ix+ix0) - 0.5*dx =', x(ix + ix0) &
                      - 0.5 * dx
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
                      - 0.5 * dx
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
                      - 0.5 * dy
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
                      - 0.5 * dy
                  print *, 'jy =', jy
                  print *, 'jy0 =', jy0
                  print *, 'y(jy+jy0) =', y(jy + jy0)
                  print *, 'dy =', dy
                  print *, 'iRay,ix,jy,kz =', iRay, ix, jy, kz
                  stop
                end if
              end if

              zr = ray(iRay, ix, jy, kz)%z

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
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), &
        intent(inout) :: ray

    integer :: ix, jy, kz

    integer, dimension(:), allocatable :: nr_merge

    integer :: ir_k, ir_l, ir_m

    integer :: jRay

    ! testb
    ! logical :: lmrglc
    ! teste

    real :: wnrk_min, wnrk_max, wnrl_min, wnrl_max, wnrm_min, wnrm_max

    real :: wnrk_min_p, wnrk_max_p, wnrl_min_p, wnrl_max_p, wnrm_min_p, &
        wnrm_max_p
    real :: wnrk_min_n, wnrk_max_n, wnrl_min_n, wnrl_max_n, wnrm_min_n, &
        wnrm_max_n
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
        zrmxmg

    real, dimension(:), allocatable :: krmnmg, krmxmg, lrmnmg, lrmxmg, mrmnmg, &
        mrmxmg

    real, dimension(:), allocatable :: wadrmg

    real :: fcpspx, fcpspy, fcpspz

    real :: omir, NNr

    real :: f_cor_nd

    real :: wa_old, wa_new

    real :: en_old, en_new
    real :: endn_old, endn_new

    real :: wnrt

    integer :: nrvtt0, nrvtt1, nrvloc

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

    nrvloc = sum(nRay(1:nx, 1:ny, 0:nz))

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
                      + 1
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
                    =', nxRay - 1
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
                      + 1
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
                    =', nyRay - 1
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
                  nzRay - 1
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
                    * (nxRay - 1) + ir_k
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
                    * fcpspz
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
                / sqrt(wnrh ** 2 + wnrm ** 2)

            ray(iRay, ix, jy, kz)%omega = omir

            ! phase-space volumes

            ray(iRay, ix, jy, kz)%area_xk = ray(iRay, ix, jy, kz)%dxray &
                * ray(iRay, ix, jy, kz)%dkray

            ray(iRay, ix, jy, kz)%area_yl = ray(iRay, ix, jy, kz)%dyray &
                * ray(iRay, ix, jy, kz)%dlray

            ray(iRay, ix, jy, kz)%area_zm = ray(iRay, ix, jy, kz)%dzray &
                * ray(iRay, ix, jy, kz)%dmray

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
                  * fcpspz)
            elseif(cons_merge == "en") then
              ! wave-action density after merging to be determined
              ! such that the wave energy remains the same, hence ...

              ray(iRay, ix, jy, kz)%dens = wadrmg(jRay) / (omir * fcpspx &
                  * fcpspy * fcpspz)
            else
              stop 'wrong cons_merge in merge_rayvol'
            end if

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
                iRay, '> nray_max =', nray_max
            stop
          else
            nRay(ix, jy, kz) = iRay
          end if
        end do
      end do
    end do

    ! total number of ray volumes before after merging

    nrvloc = sum(nRay(1:nx, 1:ny, 0:nz))

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
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(in) :: var
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), &
        intent(inout) :: ray

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

    ix0 = is + nbx - 1
    jy0 = js + nby - 1

    f_cor_nd = f_Coriolis_dim * tRef

    ! relaunching lower-boundary ray volumes

    if(case_wkb == 3 .and. RKstage == 1) then
      ! testb
      ! print*,'max u(:,:,1) =',maxval(abs(var(:,:,1,2)))*uRef
      ! if (cgz_max /= 0.0) then
      !    print*,'dz/cgz_max =',dz/cgz_max * tRef / 86400.0,'days'
      !   else
      !    print*,'dz/cgz_max = infinity'
      ! end if
      ! print*,'nRay =',sum(nRay(1:nx,1:ny,0:nz))
      ! teste

      ! FJFeb2023
      ! if (wlrx_init /= 0.0) then
      !   wnrk_init = 2.0 * pi / wlrx_init * lRef
      ! else
      !   wnrk_init = 0.0
      ! end if

      ! if (wlry_init /= 0.0) then
      !   wnrl_init = 2.0 * pi / wlry_init * lRef
      ! else
      !   wnrl_init = 0.0
      ! end if

      ! wnrh_init = sqrt(wnrk_init ** 2 + wnrl_init ** 2)

      ! wnrm_init = 2.0 * pi / wlrz_init * lRef

      ! FJJan2023
      ! xr0 = xr0_dim / lRef
      ! yr0 = yr0_dim / lRef
      ! zr0 = zr0_dim / lRef

      ! FJJan2023
      ! sigwpx = sigwpx_dim / lRef
      ! sigwpy = sigwpx_dim / lRef
      ! sigwpz = sigwpz_dim / lRef

      ! FJFeb2023
      if(topographyTime_wkb > 0.0) then
        call update_topography_wkb(time)
      end if

      ! nondimensional wave-number widths to be filled by ray volumes

      dk_ini_nd = dk_init * lRef
      dl_ini_nd = dl_init * lRef
      dm_ini_nd = dm_init * lRef

      ! local squared Brunt-Vaisala frequency

      call stratification(z(0), 1, NN_nd)

      kz = 0

      do jy = 1, ny
        do ix = 1, nx
          do i_sfc = 1, n_sfc
            iRay = ir_sfc(i_sfc, ix, jy)

            ! FJMar2023
            wnrk_init = k_spectrum(ix, jy)
            wnrl_init = l_spectrum(ix, jy)
            wnrh_init = sqrt(wnrk_init ** 2.0 + wnrl_init ** 2.0)
            wnrm_init = 0.0

            ! three cases:
            ! (1) no ray volume launched previously (iRay < 0)
            !     => check whether new ray volume is to be launched now
            ! (2) top of previously launched r.v. has passed lower
            !     boundary
            !     => clip old ray volume and keep part above lower
            !        boundary
            !        check whether new ray volume is to be launched
            !        now
            ! (3) top of previously launched r.v. below the model
            !     bottom
            !     => set wave-action density to zero
            !        check whether new ray volume is to be launched now

            if(iRay > 0) then
              zr = ray(iRay, ix, jy, kz)%z
              dzr = ray(iRay, ix, jy, kz)%dzray
            else
              zr = 0.0
              dzr = 0.0
            end if

            ix2 = ix2_sfc(i_sfc)
            jy2 = jy2_sfc(i_sfc)
            kz2 = nrzl

            ik = ik_sfc(i_sfc)
            jl = jl_sfc(i_sfc)
            km = km_sfc(i_sfc)

            omir = - var(ix, jy, 1, 2) * wnrk_init - var(ix, jy, 1, 3) &
                * wnrl_init

            ! choose correct sign of horizontal wavenumbers in order to
            ! be on the correct frequency branch

            if(omir * branchr >= 0.0) then
              wnrk = wnrk_init
              wnrl = wnrl_init
            else
              omir = - omir

              wnrk = - wnrk_init
              wnrl = - wnrl_init
            end if

            amp = 0.0
            wnrm = 0.0

            ! FJJan2023
            ! if ((sigwpx == 0.0 .or. abs(x(ix + ix0) - xr0) < sigwpx) .and. &
            !     (sigwpy == 0.0 .or. abs(y(jy + jy0) - yr0) < sigwpy)) then
            if(abs(omir) <= f_cor_nd) then
              amp = 0.0
              wnrm = 0.0
            elseif(omir ** 2 < NN_nd) then
              wnrm = - branchr * sqrt(wnrh_init ** 2 * (NN_nd - omir ** 2) &
                  / (omir ** 2 - f_cor_nd ** 2))

              ! Displacement (FJFeb2023)
              displm = topography_spectrum(ix, jy)

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

              ! Long number scaling (FJJan2023)
              if(long_scaling) then
                ! Compute Long number.
                long = displm * sqrt(NN_nd / (var(ix, jy, 1, 2) ** 2.0 &
                    + var(ix, jy, 1, 3) ** 2.0))
                if(long > 0.0) then
                  if(long_fit == 0) then
                    ! Apply scaling from decision tree model.
                    displm = displm * compute_long_scaling(long)
                  else if(long_fit == 1) then
                    ! Apply exponential scaling.
                    displm = displm * along * (1.0 - exp(- blong / long))
                  else if(long_fit == 2) then
                    ! Apply fractional scaling.
                    displm = displm * along / long / (1.0 + blong / long)
                  else if(long_fit == 3) then
                    ! Apply linear scaling.
                    displm = min(displm, along * displm / long)
                  end if
                end if
              end if

              ! surface wave-action density
              amp = 0.5 * rhoStrat(0) * displm ** 2 * omir * (wnrh_init ** 2 &
                  + wnrm ** 2) / wnrh_init ** 2
            else
              amp = 0.0
              wnrm = 0.0
            end if
            ! end if

            ! only launch new ray volume if wave-action-density is
            ! non-zero

            if(amp /= 0.0) then
              ! for cases (1) and (2) increase the number of r.v.
              ! in case (3) the old r.v. is replaced by the new one
              if(iRay < 0 .or. (iRay > 0 .and. zr + 0.5 * dzr > lz(0))) then
                nRay(ix, jy, kz) = nRay(ix, jy, kz) + 1

                if(nRay(ix, jy, kz) > nray_wrk) then
                  print *, 'ERROR at ix,jy,kz =', ix, jy, kz
                  print *, 'nRay =', nRay(ix, jy, kz), '> nray_wrk =', nray_wrk
                  stop
                end if
              end if

              if(iRay > 0 .and. zr + 0.5 * dzr > lz(0)) then
                nrlc = nRay(ix, jy, kz)

                ! case (2):
                ! move old ray volume to last in the row
                ray(nrlc, ix, jy, kz) = ray(iRay, ix, jy, kz)

                ! clip it so that only the part above the lower
                ! boundary is kept

                if(ray(nrlc, ix, jy, kz)%z - 0.5 * ray(nrlc, ix, jy, kz)%dzray &
                    < lz(0)) then
                  ray(nrlc, ix, jy, kz)%dzray = ray(nrlc, ix, jy, kz)%z + 0.5 &
                      * ray(nrlc, ix, jy, kz)%dzray - lz(0)

                  ray(nrlc, ix, jy, kz)%z = lz(0) + 0.5 * ray(nrlc, ix, jy, &
                      kz)%dzray

                  ray(nrlc, ix, jy, kz)%area_zm = ray(nrlc, ix, jy, kz)%dzray &
                      * ray(nrlc, ix, jy, kz)%dmray
                end if
              elseif(iRay < 0) then
                ! case (1):

                iRay = nRay(ix, jy, kz)
                ir_sfc(i_sfc, ix, jy) = iRay
              end if
              ! case (3): keep same iRay and fill it with new r.v.
              ! ie: nothing to be done here!
            else
              ir_sfc(i_sfc, ix, jy) = - 1
              cycle
            end if

            ! ray-volume positions

            ray(iRay, ix, jy, kz)%x = (x(ix + ix0) - dx / 2.0 + (ix2 - 0.5) &
                * dx / nrxl)
            ray(iRay, ix, jy, kz)%y = (y(jy + jy0) - dy / 2.0 + (jy2 - 0.5) &
                * dy / nryl)
            ray(iRay, ix, jy, kz)%z = (z(kz) - dz / 2.0 + (kz2 - 0.5) * dz &
                / nrzl)

            ! ray-volume spatial extensions

            ray(iRay, ix, jy, kz)%dxray = dx / nrxl
            ray(iRay, ix, jy, kz)%dyray = dy / nryl
            ray(iRay, ix, jy, kz)%dzray = dz / nrzl

            ! ray-volume wave numbers

            wnk_0 = wnrk
            wnl_0 = wnrl
            wnm_0 = wnrm

            ray(iRay, ix, jy, kz)%k = (wnk_0 - 0.5 * dk_ini_nd + (real(ik) &
                - 0.5) * dk_ini_nd / nrk_init)

            ray(iRay, ix, jy, kz)%l = (wnl_0 - 0.5 * dl_ini_nd + (real(jl) &
                - 0.5) * dl_ini_nd / nrl_init)

            if(fac_dm_init == 0.0) then
              stop 'ERROR: FAC_DM_INIT = 0.0'
            else if(wnm_0 == 0.0) then
              stop 'ERROR: WNM_0 = 0.0'
            else
              dm_ini_nd = fac_dm_init * abs(wnm_0)
            end if

            ray(iRay, ix, jy, kz)%m = (wnm_0 - 0.5 * dm_ini_nd + (real(km) &
                - 0.5) * dm_ini_nd / nrm_init)

            ! ray-volume wave-number extents

            ray(iRay, ix, jy, kz)%dkray = dk_ini_nd / nrk_init
            ray(iRay, ix, jy, kz)%dlray = dl_ini_nd / nrl_init
            ray(iRay, ix, jy, kz)%dmray = dm_ini_nd / nrm_init

            ! ray-volume phase-space volume

            ray(iRay, ix, jy, kz)%area_xk = ray(iRay, ix, jy, kz)%dxray &
                * ray(iRay, ix, jy, kz)%dkray
            ray(iRay, ix, jy, kz)%area_yl = ray(iRay, ix, jy, kz)%dyray &
                * ray(iRay, ix, jy, kz)%dlray
            ray(iRay, ix, jy, kz)%area_zm = ray(iRay, ix, jy, kz)%dzray &
                * ray(iRay, ix, jy, kz)%dmray

            pspvol = dm_ini_nd

            if(fac_dk_init /= 0.0) then
              pspvol = pspvol * dk_ini_nd
            end if

            if(fac_dl_init /= 0.0) then
              pspvol = pspvol * dl_ini_nd
            end if

            ! phase-space wave-action density

            ray(iRay, ix, jy, kz)%dens = amp / pspvol

            ! intrinsic frequency

            ray(iRay, ix, jy, kz)%omega = omir
          end do ! i_sfc
        end do ! jy
      end do ! ix
    endif

    ! initialize RK-tendencies at first RK stage
    if(RKstage == 1) then
      dxRay = 0.0
      dkRay = 0.0
      ddxRay = 0.0
    end if

    cgx_max = 0.0
    cgy_max = 0.0
    cgz_max = 0.0

    do kz = 0, nz
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
            if(zr1 < lz(0) - dz) then
              nskip = nskip + 1
              cycle
            end if

            call stratification(zr1, 1, NNr1)
            call stratification(zr, 1, NNr)
            call stratification(zr2, 1, NNr2)

            omir1 = branchr * sqrt(NNr1 * wnrh ** 2 + f_cor_nd ** 2 * wnrm &
                ** 2) / sqrt(wnrh ** 2 + wnrm ** 2)

            omir = branchr * sqrt(NNr * wnrh ** 2 + f_cor_nd ** 2 * wnrm ** 2) &
                / sqrt(wnrh ** 2 + wnrm ** 2)

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
                ** 2) / sqrt(wnrh ** 2 + wnrm ** 2)

            ray(iRay, ix, jy, kz)%omega = omir

            ! intrinsic group velocities at the respective edges of
            ! the ray volumes

            if(sizeX > 1) then
              ! intrinsic group velocity in x direction not depending
              ! on x
              cgirx = wnrk * (NNr - omir ** 2) / (omir * (wnrh ** 2 + wnrm &
                  ** 2))
            end if

            if(sizeY > 1) then
              ! intrinsic group velocity in y direction not depending
              ! on y
              cgiry = wnrl * (NNr - omir ** 2) / (omir * (wnrh ** 2 + wnrm &
                  ** 2))
            end if

            ! intrinsic vertical group velocity depending on z
            ! (via the stratification)
            cgirz1 = - wnrm * (omir1 ** 2 - f_cor_nd ** 2) / (omir1 * (wnrh &
                ** 2 + wnrm ** 2))
            cgirz2 = - wnrm * (omir2 ** 2 - f_cor_nd ** 2) / (omir2 * (wnrh &
                ** 2 + wnrm ** 2))

            !--------------------------------
            !     displacement in x direction
            !--------------------------------

            ! RK update

            if(sizeX > 1) then
              call meanflow(xr1, yr, zr, var, 1, uxr1)
              call meanflow(xr2, yr, zr, var, 1, uxr2)

              ! group velocity in x direction at the two edges in x
              cgrx1 = cgirx + uxr1
              cgrx2 = cgirx + uxr2

              ! group velocity in x direction for the carrier ray
              cgrx = 0.5 * (cgrx1 + cgrx2)

              F = cgrx !allow horizontal ray propagation
              dxRay(1, iRay, ix, jy, kz) = dt * F + alpha(rkStage) * dxRay(1, &
                  iRay, ix, jy, kz)
              ray(iRay, ix, jy, kz)%x = ray(iRay, ix, jy, kz)%x &
                  + beta(RKstage) * dxRay(1, iRay, ix, jy, kz)

              ! update maximum group velocity in x direction
              cgx_max = max(cgx_max, abs(cgrx))
            end if

            !-----------------------------
            !     displacement in y direction
            !-----------------------------

            ! RK update

            if(sizeY > 1) then
              call meanflow(xr, yr1, zr, var, 2, vyr1)
              call meanflow(xr, yr2, zr, var, 2, vyr2)

              ! group velocity in y direction at the two edges in y
              cgry1 = cgiry + vyr1
              cgry2 = cgiry + vyr2

              ! group velocity in y direction for the carrier ray
              cgry = 0.5 * (cgry1 + cgry2)

              F = cgry !allow horizontal ray propagation
              dxRay(2, iRay, ix, jy, kz) = dt * F + alpha(rkStage) * dxRay(2, &
                  iRay, ix, jy, kz)
              ray(iRay, ix, jy, kz)%y = ray(iRay, ix, jy, kz)%y &
                  + beta(RKstage) * dxRay(2, iRay, ix, jy, kz)

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
            dxRay(3, iRay, ix, jy, kz) = dt * F + alpha(rkStage) * dxRay(3, &
                iRay, ix, jy, kz)
            ray(iRay, ix, jy, kz)%z = ray(iRay, ix, jy, kz)%z + beta(RKstage) &
                * dxRay(3, iRay, ix, jy, kz)

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
                    iRay, ix, jy, kz, 'TOO LOW'
                stop
              end if

              if(zr < lz(0) - dz) then
                print *, 'ERROR IN transport_rayvol: RAY VOLUME', iRay, ix, &
                    jy, kz, 'TOO LOW'
                stop
              end if

              call stratification(zr, 2, dnndzr)

              dkdt = - dudxr * wnrk - dvdxr * wnrl
              dldt = - dudyr * wnrk - dvdyr * wnrl
              dmdt = - dudzr * wnrk - dvdzr * wnrl - wnrh ** 2 * dnndzr / (2.0 &
                  * omir + (wnrh ** 2 + wnrm ** 2))

              dkRay(1, iRay, ix, jy, kz) = dt * dkdt + alpha(rkStage) &
                  * dkRay(1, iRay, ix, jy, kz)

              dkRay(2, iRay, ix, jy, kz) = dt * dldt + alpha(rkStage) &
                  * dkRay(2, iRay, ix, jy, kz)

              dkRay(3, iRay, ix, jy, kz) = dt * dmdt + alpha(rkStage) &
                  * dkRay(3, iRay, ix, jy, kz)

              ray(iRay, ix, jy, kz)%k = ray(iRay, ix, jy, kz)%k &
                  + beta(rkStage) * dkRay(1, iRay, ix, jy, kz)

              ray(iRay, ix, jy, kz)%l = ray(iRay, ix, jy, kz)%l &
                  + beta(rkStage) * dkRay(2, iRay, ix, jy, kz)

              ray(iRay, ix, jy, kz)%m = ray(iRay, ix, jy, kz)%m &
                  + beta(rkStage) * dkRay(3, iRay, ix, jy, kz)

              !----------------------------------------------
              !    change of wave-number width of ray volumes
              !----------------------------------------------

              ! dk

              if(sizeX > 1) then
                ddxdt = cgrx2 - cgrx1

                ddxRay(1, iRay, ix, jy, kz) = dt * ddxdt + alpha(rkStage) &
                    * ddxRay(1, iRay, ix, jy, kz)

                ray(iRay, ix, jy, kz)%dxray = ray(iRay, ix, jy, kz)%dxray &
                    + beta(rkStage) * ddxRay(1, iRay, ix, jy, kz)

                if(ray(iRay, ix, jy, kz)%dxray <= 0.0) then
                  print *, 'dxray(', iRay, ix, jy, kz, ') <= 0.0  ==> time &
                      step too large?'
                  ray(iRay, ix, jy, kz)%dxray = - ray(iRay, ix, jy, kz)%dxray
                end if

                ray(iRay, ix, jy, kz)%dkray = ray(iRay, ix, jy, kz)%area_xk &
                    / ray(iRay, ix, jy, kz)%dxray
              end if

              ! dl

              if(sizeY > 1) then
                ddydt = cgry2 - cgry1

                ddxRay(2, iRay, ix, jy, kz) = dt * ddydt + alpha(rkStage) &
                    * ddxRay(2, iRay, ix, jy, kz)

                ray(iRay, ix, jy, kz)%dyray = ray(iRay, ix, jy, kz)%dyray &
                    + beta(rkStage) * ddxRay(2, iRay, ix, jy, kz)

                if(ray(iRay, ix, jy, kz)%dyray <= 0.0) then
                  print *, 'dyray(', iRay, ix, jy, kz, ') <= 0.0  ==> time &
                      step too large?'
                  ray(iRay, ix, jy, kz)%dyray = - ray(iRay, ix, jy, kz)%dyray
                end if

                ray(iRay, ix, jy, kz)%dlray = ray(iRay, ix, jy, kz)%area_yl &
                    / ray(iRay, ix, jy, kz)%dyray
              end if

              ! dm

              ddzdt = cgrz2 - cgrz1

              ddxRay(3, iRay, ix, jy, kz) = dt * ddzdt + alpha(rkStage) &
                  * ddxRay(3, iRay, ix, jy, kz)

              ray(iRay, ix, jy, kz)%dzray = ray(iRay, ix, jy, kz)%dzray &
                  + beta(rkStage) * ddxRay(3, iRay, ix, jy, kz)

              if(ray(iRay, ix, jy, kz)%dzray <= 0.0) then
                print *, 'dzray(', iRay, ix, jy, kz, ') <= 0.0  ==> time step &
                    too large?'
                ray(iRay, ix, jy, kz)%dzray = - ray(iRay, ix, jy, kz)%dzray
              end if

              ray(iRay, ix, jy, kz)%dmray = ray(iRay, ix, jy, kz)%area_zm &
                  / ray(iRay, ix, jy, kz)%dzray
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
                / sqrt(wnrh ** 2 + wnrm ** 2)

            ray(iRay, ix, jy, kz)%omega = omir
          end do ray_loop

          if(nskip > 0) then
            print *, nskip, 'r.v. skipped in transport_rayvol out of', &
                nRay(ix, jy, kz)
          end if
        end do ! ix
      end do ! jy
    end do ! kz

  end subroutine transport_rayvol

  !---------------------------------------------------------------------

  subroutine boundary_rayvol(ray)

    !------------------------------------------------------------------
    ! re-positioning of ray volumes so that they respect the periodic
    ! boundary conditions
    !
    ! in the case of solid-wall boundary conditions in the vertical,
    ! wave-action densities of ray volumes are set to zero if they have
    ! left the domain
    !------------------------------------------------------------------

    implicit none

    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), &
        intent(inout) :: ray
    ! time step

    integer :: iRay
    integer :: ix, jy, kz

    real :: xr, yr, zr
    real :: dzr

    do kz = 1, nz
      do jy = 1, ny
        do ix = 1, nx
          do iRay = 1, nRay(ix, jy, kz)
            ! implement horizontal boundary conditions for ray-volume
            ! positions

            if(sizeX > 1) then
              xr = ray(iRay, ix, jy, kz)%x

              if(xr < lx(0)) then
                select case(xBoundary)
                case("periodic")
                  xr = lx(1) + mod(xr - lx(0), lx(1) - lx(0))
                case default
                  stop "transport_rayvol: unknown case xBoundary"
                end select

                ray(iRay, ix, jy, kz)%x = xr
              elseif(xr > lx(1)) then
                select case(xBoundary)
                case("periodic")
                  xr = lx(0) + mod(xr - lx(1), lx(1) - lx(0))
                case default
                  stop "transport_rayvol: unknown case xBoundary"
                end select

                ray(iRay, ix, jy, kz)%x = xr
              end if
            end if

            if(sizeY > 1) then
              yr = ray(iRay, ix, jy, kz)%y

              if(yr < ly(0)) then
                select case(yBoundary)
                case("periodic")
                  yr = ly(1) + mod(yr - ly(0), ly(1) - ly(0))
                case default
                  stop "transport_rayvol: unknown case yBoundary"
                end select

                ray(iRay, ix, jy, kz)%y = yr
              elseif(yr > ly(1)) then
                select case(yBoundary)
                case("periodic")
                  yr = ly(0) + mod(yr - ly(1), ly(1) - ly(0))
                case default
                  stop "transport_rayvol: unknown case yBoundary"
                end select

                ray(iRay, ix, jy, kz)%y = yr
              end if
            end if

            ! vertical boundary conditions:
            ! zBoundary = 'periodic': implement periodicity
            ! zBoundary = 'solid_wall': remove wave-action density
            !                           for ray volumes that have
            !                           completely left the model
            !                           domain

            zr = ray(iRay, ix, jy, kz)%z

            if(zr < lz(0)) then
              select case(zBoundary)
              case("periodic")
                zr = lz(1) + mod(zr - lz(0), lz(1) - lz(0))

                ray(iRay, ix, jy, kz)%z = zr
              case("solid_wall")
                dzr = ray(iRay, ix, jy, kz)%dzray

                if(zr + 0.5 * dzr < lz(0)) ray(iRay, ix, jy, kz)%dens = 0.0
              case default
                stop "transport_rayvol: unknown case zBoundary"
              end select
            elseif(zr > lz(1)) then
              select case(zBoundary)
              case("periodic")
                zr = lz(0) + mod(zr - lz(1), lz(1) - lz(0))

                ray(iRay, ix, jy, kz)%z = zr
              case("solid_wall")
                dzr = ray(iRay, ix, jy, kz)%dzray

                if(zr - 0.5 * dzr > lz(1)) ray(iRay, ix, jy, kz)%dens = 0.0
              case default
                stop "transport_rayvol: unknown case zBoundary"
              end select
            end if
          end do
        end do ! ix
      end do ! jy
    end do ! kz

  end subroutine boundary_rayvol

  !---------------------------------------------------------------------

  subroutine saturation_3D(ray, dt)

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

    implicit none

    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), &
        intent(inout) :: ray
    ! time step
    real, intent(in) :: dt

    ! indices, etc.
    integer iRay, kzmax, kzmin

    integer :: ix, jy, kz

    integer :: ix0, jy0

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

            ! vertical boundary conditions:
            ! zBoundary = 'periodic': implement periodicity
            ! zBoundary = 'solid_wall': skip counting ray volumes
            !                           that have completely left the
            !                           model domain

            if(zr < lz(0)) then
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

            kz = floor((zr - lz(0)) / dz) + 1

            !  extra skip counting rays propagating out of the domain
            if(kz < 1 .or. kz > sizeZ) cycle

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

            call stratification(z(kz), 1, NN_nd)

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
            dzi = min(dzr, dz)

            facpsp = dzi / dz * dwnrm

            if(fac_dk_init /= 0.0) then
              dxi = min(dxr, dx)

              facpsp = facpsp * dxi / dx * dwnrk
            end if

            if(fac_dl_init /= 0.0) then
              dyi = min(dyr, dy)

              facpsp = facpsp * dyi / dy * dwnrl
            end if

            integral1 = wnrhs * wnrm ** 2 / ((wnrhs + wnrm ** 2) * omir) &
                * facpsp

            mB2(ix, jy, kz) = mB2(ix, jy, kz) + 2.0 * NN_nd ** 2 &
                / rhoStrat(kz) * densr * integral1

            integral2 = wnrhs * wnrm ** 2 / omir * facpsp

            mB2K2(ix, jy, kz) = mB2K2(ix, jy, kz) + 2.0 * NN_nd ** 2 &
                / rhoStrat(kz) * densr * integral2
          end do
        end do ! ixrv
      end do ! jyrv
    end do ! kzrv

    ! loop for computing the diffusivity coefficient

    do kz = 1, sizeZ
      call stratification(z(kz), 1, NN_nd)

      do jy = 1, ny
        do ix = 1, nx
          if(mB2K2(ix, jy, kz) == 0.0 .or. mB2(ix, jy, kz) < alpha_sat ** 2 &
              * NN_nd ** 2) then
            diffusion(ix, jy, kz) = 0.0
          else
            diffusion(ix, jy, kz) = (mB2(ix, jy, kz) - alpha_sat ** 2 * NN_nd &
                ** 2) / (2.0 * dt * mB2K2(ix, jy, kz))
          endif
        end do
      end do
    end do

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

            ! vertical boundary conditions:
            ! zBoundary = 'periodic': implement periodicity
            ! zBoundary = 'solid_wall': skip counting ray volumes
            !                           that have completely left the
            !                           model domain

            if(zr < lz(0)) then
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

            kz = floor((zr - lz(0)) / dz) + 1

            !  extra skip counting rays propagating out of the domain
            if(kz < 1 .or. kz > sizeZ) cycle

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

            wnrk = ray(iRay, ixrv, jyrv, kzrv)%k
            wnrl = ray(iRay, ixrv, jyrv, kzrv)%l
            wnrm = ray(iRay, ixrv, jyrv, kzrv)%m

            kappa = diffusion(ix, jy, kz)

            ! it can hanppen that the damping coefficient becomes
            ! negative,
            ! hence set it to zero in that case
            ray(iRay, ixrv, jyrv, kzrv)%dens = ray(iRay, ixrv, jyrv, &
                kzrv)%dens * max(0.0, 1.0 - dt * 2.0 * kappa * (wnrk ** 2 &
                + wnrl ** 2 + wnrm ** 2))
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

            ! vertical boundary conditions:
            ! zBoundary = 'periodic': implement periodicity
            ! zBoundary = 'solid_wall': skip counting ray volumes
            !                           that have completely left the
            !                           model domain

            if(zr < lz(0)) then
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

            kz = floor((zr - lz(0)) / dz) + 1

            !  extra skip counting rays propagating out of the domain
            if(kz < 1 .or. kz > sizeZ) cycle

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

            call stratification(z(kz), 1, NN_nd)

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
            dzi = min(dzr, dz)

            facpsp = dzi / dz * dwnrm

            if(fac_dk_init /= 0.0) then
              dxi = min(dxr, dx)

              facpsp = facpsp * dxi / dx * dwnrk
            end if

            if(fac_dl_init /= 0.0) then
              dyi = min(dyr, dy)

              facpsp = facpsp * dyi / dy * dwnrl
            end if

            integral1 = wnrhs * wnrm ** 2 / ((wnrhs + wnrm ** 2) * omir) &
                * facpsp

            mB2(ix, jy, kz) = mB2(ix, jy, kz) + 2.0 * NN_nd ** 2 &
                / rhoStrat(kz) * densr * integral1
          end do
        end do ! ixrv
      end do ! jyrv
    end do ! kzrv

    do kz = 1, sizeZ
      call stratification(z(kz), 1, NN_nd)

      do jy = 1, ny
        do ix = 1, nx
          if(mB2(ix, jy, kz) - alpha_sat ** 2 * NN_nd ** 2 > 1.d-3 * alpha_sat &
              ** 2 * NN_nd ** 2) then
            print *, 'SATURATION VIOLATED AT ix, jy, kz =', ix, jy, kz
            print *, 'mB2(ix,jy,kz) =', mB2(ix, jy, kz)
            print *, 'alpha_sat**2 * NN_nd**2 = ', alpha_sat ** 2 * NN_nd ** 2
            ! stop
          endif
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
        intent(inout) :: flxwkb
    integer, intent(in) :: nsmth

    integer, intent(in) :: homog_dir

    ! allocatable fields
    real, dimension(:, :, :), allocatable :: flxwkb_0, flxwkb_1

    integer :: allocstat
    integer :: i, j, k

    allocate(flxwkb_0(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        = allocstat)
    if(allocstat /= 0) stop "smooth_wkb_shapiro:alloc failed"

    allocate(flxwkb_1(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        = allocstat)
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
                  + 2.0 * flxwkb(i, j, k)) / 4.0
            end do
          end do
        end do

        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, &
                  k) + 2.0 * flxwkb_0(i, j, k)) / 4.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1) &
                  + 2.0 * flxwkb_1(i, j, k)) / 4.0
            end do
          end do
        end do
      elseif(nsmth == 2) then
        ! smooth in x

        do k = - nbz, nz + nbz
          do j = - nby, ny + nby
            do i = 1, nx
              flxwkb_0(i, j, k) = (- flxwkb(i - 2, j, k) - flxwkb(i + 2, j, k) &
                  + 4.0 * (flxwkb(i - 1, j, k) + flxwkb(i + 1, j, k)) + 10.0 &
                  * flxwkb(i, j, k)) / 16.0
            end do
          end do
        end do

        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (- flxwkb_0(i, j - 2, k) - flxwkb_0(i, j &
                  + 2, k) + 4.0 * (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, &
                  k)) + 10.0 * flxwkb_0(i, j, k)) / 16.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (- flxwkb_1(i, j, k - 2) - flxwkb_1(i, j, k &
                  + 2) + 4.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  + 10.0 * flxwkb_1(i, j, k)) / 16.0
            end do
          end do
        end do
      elseif(nsmth == 3) then
        ! smooth in x

        do k = - nbz, nz + nbz
          do j = - nby, ny + nby
            do i = 1, nx
              flxwkb_0(i, j, k) = (flxwkb(i - 3, j, k) + flxwkb(i + 3, j, k) &
                  - 6.0 * (flxwkb(i - 2, j, k) + flxwkb(i + 2, j, k)) + 15.0 &
                  * (flxwkb(i - 1, j, k) + flxwkb(i + 1, j, k)) + 44.0 &
                  * flxwkb(i, j, k)) / 64.0
            end do
          end do
        end do

        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (flxwkb_0(i, j - 3, k) + flxwkb_0(i, j + 3, &
                  k) - 6.0 * (flxwkb_0(i, j - 2, k) + flxwkb_0(i, j + 2, k)) &
                  + 15.0 * (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, k)) &
                  + 44.0 * flxwkb_0(i, j, k)) / 64.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 3) + flxwkb_1(i, j, k + 3) &
                  - 6.0 * (flxwkb_1(i, j, k - 2) + flxwkb_1(i, j, k + 2)) &
                  + 15.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  + 44.0 * flxwkb_1(i, j, k)) / 64.0
            end do
          end do
        end do
      elseif(nsmth == 4) then
        ! smooth in x

        do k = - nbz, nz + nbz
          do j = - nby, ny + nby
            do i = 1, nx
              flxwkb_0(i, j, k) = (- flxwkb(i - 4, j, k) - flxwkb(i + 4, j, k) &
                  + 8.0 * (flxwkb(i - 3, j, k) + flxwkb(i + 3, j, k)) - 28.0 &
                  * (flxwkb(i - 2, j, k) + flxwkb(i + 2, j, k)) + 56.0 &
                  * (flxwkb(i - 1, j, k) + flxwkb(i + 1, j, k)) + 186.0 &
                  * flxwkb(i, j, k)) / 256.0
            end do
          end do
        end do

        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (- flxwkb_0(i, j - 4, k) - flxwkb_0(i, j &
                  + 4, k) + 8.0 * (flxwkb_0(i, j - 3, k) + flxwkb_0(i, j + 3, &
                  k)) - 28.0 * (flxwkb_0(i, j - 2, k) + flxwkb_0(i, j + 2, k)) &
                  + 56.0 * (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, k)) &
                  + 186.0 * flxwkb_0(i, j, k)) / 256.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 4) + flxwkb_1(i, j, k + 4) &
                  + 8.0 * (flxwkb_1(i, j, k - 3) + flxwkb_1(i, j, k + 3)) &
                  - 28.0 * (flxwkb_1(i, j, k - 2) + flxwkb_1(i, j, k + 2)) &
                  + 56.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  + 186.0 * flxwkb_1(i, j, k)) / 256.0
            end do
          end do
        end do
      else
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = sum(flxwkb_0((i - nsmth):(i + nsmth), (j &
                  - nsmth):(j + nsmth), (k - nsmth):(k + nsmth))) / real((2 &
                  * nsmth + 1) ** 3)
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
                  k) + 2.0 * flxwkb_0(i, j, k)) / 4.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1) &
                  + 2.0 * flxwkb_1(i, j, k)) / 4.0
            end do
          end do
        end do
      elseif(nsmth == 2) then
        ! smooth in x

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (- flxwkb_0(i - 2, j, k) - flxwkb_0(i + 2, &
                  j, k) + 4.0 * (flxwkb_0(i - 1, j, k) + flxwkb_0(i + 1, j, &
                  k)) + 10.0 * flxwkb_0(i, j, k)) / 16.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (- flxwkb_1(i, j, k - 2) - flxwkb_1(i, j, k &
                  + 2) + 4.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  + 10.0 * flxwkb_1(i, j, k)) / 16.0
            end do
          end do
        end do
      elseif(nsmth == 3) then
        ! smooth in x

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (flxwkb_0(i - 3, j, k) + flxwkb_0(i + 3, j, &
                  k) - 6.0 * (flxwkb_0(i - 2, j, k) + flxwkb_0(i + 2, j, k)) &
                  + 15.0 * (flxwkb_0(i - 1, j, k) + flxwkb_0(i + 1, j, k)) &
                  + 44.0 * flxwkb_0(i, j, k)) / 64.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 3) + flxwkb_1(i, j, k + 3) &
                  - 6.0 * (flxwkb_1(i, j, k - 2) + flxwkb_1(i, j, k + 2)) &
                  + 15.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  + 44.0 * flxwkb_1(i, j, k)) / 64.0
            end do
          end do
        end do
      elseif(nsmth == 4) then
        ! smooth in x

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (- flxwkb_0(i - 4, j, k) - flxwkb_0(i + 4, &
                  j, k) + 8.0 * (flxwkb_0(i - 3, j, k) + flxwkb_0(i + 3, j, &
                  k)) - 28.0 * (flxwkb_0(i - 2, j, k) + flxwkb_0(i + 2, j, k)) &
                  + 56.0 * (flxwkb_0(i - 1, j, k) + flxwkb_0(i + 1, j, k)) &
                  + 186.0 * flxwkb_0(i, j, k)) / 256.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 4) + flxwkb_1(i, j, k + 4) &
                  + 8.0 * (flxwkb_1(i, j, k - 3) + flxwkb_1(i, j, k + 3)) &
                  - 28.0 * (flxwkb_1(i, j, k - 2) + flxwkb_1(i, j, k + 2)) &
                  + 56.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  + 186.0 * flxwkb_1(i, j, k)) / 256.0
            end do
          end do
        end do
      else
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = sum(flxwkb_0((i - nsmth):(i + nsmth), j, (k &
                  - nsmth):(k + nsmth))) / real((2 * nsmth + 1) ** 2)
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
                  k) + 2.0 * flxwkb_0(i, j, k)) / 4.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1) &
                  + 2.0 * flxwkb_1(i, j, k)) / 4.0
            end do
          end do
        end do
      elseif(nsmth == 2) then
        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (- flxwkb_0(i, j - 2, k) - flxwkb_0(i, j &
                  + 2, k) + 4.0 * (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, &
                  k)) + 10.0 * flxwkb_0(i, j, k)) / 16.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (- flxwkb_1(i, j, k - 2) - flxwkb_1(i, j, k &
                  + 2) + 4.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  + 10.0 * flxwkb_1(i, j, k)) / 16.0
            end do
          end do
        end do
      elseif(nsmth == 3) then
        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (flxwkb_0(i, j - 3, k) + flxwkb_0(i, j + 3, &
                  k) - 6.0 * (flxwkb_0(i, j - 2, k) + flxwkb_0(i, j + 2, k)) &
                  + 15.0 * (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, k)) &
                  + 44.0 * flxwkb_0(i, j, k)) / 64.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 3) + flxwkb_1(i, j, k + 3) &
                  - 6.0 * (flxwkb_1(i, j, k - 2) + flxwkb_1(i, j, k + 2)) &
                  + 15.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  + 44.0 * flxwkb_1(i, j, k)) / 64.0
            end do
          end do
        end do
      elseif(nsmth == 4) then
        ! smooth in y

        do k = - nbz, nz + nbz
          do j = 1, ny
            do i = 1, nx
              flxwkb_1(i, j, k) = (- flxwkb_0(i, j - 4, k) - flxwkb_0(i, j &
                  + 4, k) + 8.0 * (flxwkb_0(i, j - 3, k) + flxwkb_0(i, j + 3, &
                  k)) - 28.0 * (flxwkb_0(i, j - 2, k) + flxwkb_0(i, j + 2, k)) &
                  + 56.0 * (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, k)) &
                  + 186.0 * flxwkb_0(i, j, k)) / 256.0
            end do
          end do
        end do

        ! smooth in z

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = (flxwkb_1(i, j, k - 4) + flxwkb_1(i, j, k + 4) &
                  + 8.0 * (flxwkb_1(i, j, k - 3) + flxwkb_1(i, j, k + 3)) &
                  - 28.0 * (flxwkb_1(i, j, k - 2) + flxwkb_1(i, j, k + 2)) &
                  + 56.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
                  + 186.0 * flxwkb_1(i, j, k)) / 256.0
            end do
          end do
        end do
      else
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              flxwkb(i, j, k) = sum(flxwkb_0(i, (j - nsmth):(j + nsmth), (k &
                  - nsmth):(k + nsmth))) / real((2 * nsmth + 1) ** 2)
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
        "smooth_wkb_shapiro:dealloc failed"
    deallocate(flxwkb_1, stat = allocstat); if(allocstat /= 0) stop &
        "smooth_wkb_shapiro:dealloc failed"

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
        intent(inout) :: flxwkb
    integer, intent(in) :: nsmth

    integer, intent(in) :: homog_dir

    ! allocatable fields
    real, dimension(:, :, :), allocatable :: flxwkb_0

    integer :: allocstat
    integer :: i, j, k

    allocate(flxwkb_0(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        = allocstat)
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
                - nsmth):(j + nsmth), (k - nsmth):(k + nsmth))) / real((2 &
                * nsmth + 1) ** 3)
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
                - nsmth):(k + nsmth))) / real((2 * nsmth + 1) ** 2)
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
                - nsmth):(k + nsmth))) / real((2 * nsmth + 1) ** 2)
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
        "smooth_wkb:dealloc failed"

    return

  end subroutine smooth_wkb_box

  ! ---------------------------------------------------------------------

  subroutine setboundary_wkb(flxwkb)

    ! -------------------------------------------------------------------
    ! boundary conditions for WKB fluxes
    ! so far only periodic boundary conditions allowed in horizontal
    ! -------------------------------------------------------------------

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        intent(inout) :: flxwkb

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
        intent(inout) :: flxwkb

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
        intent(inout) :: flxwkb

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: ix

    ! auxiliary fields
    real, dimension(nbx, - nby:ny + nby, - nbz:nz + nbz) :: xSliceLeft_send, &
        xSliceRight_send
    real, dimension(nbx, - nby:ny + nby, - nbz:nz + nbz) :: xSliceLeft_recv, &
        xSliceRight_recv

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
          mpi_double_precision, dest, tag, xSliceLeft_recv(1, - nby, - nbz), &
          recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xSliceLeft_send(1, - nby, - nbz), sendcount, &
          mpi_double_precision, dest, tag, xSliceRight_recv(1, - nby, - nbz), &
          recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          sts_right, ierror)

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
        intent(inout) :: flxwkb

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: jy

    ! auxiliary fields
    real, dimension(- nbx:nx + nbx, nby, - nbz:nz + nbz) :: ySliceBack_send, &
        ySliceForw_send
    real, dimension(- nbx:nx + nbx, nby, - nbz:nz + nbz) :: ySliceBack_recv, &
        ySliceForw_recv

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
          mpi_double_precision, dest, tag, ySliceBack_recv(- nbx, 1, - nbz), &
          recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(ySliceBack_send(- nbx, 1, - nbz), sendcount, &
          mpi_double_precision, dest, tag, ySliceForw_recv(- nbx, 1, - nbz), &
          recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          sts_forw, ierror)

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
        intent(inout) :: flxwkb

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
        intent(inout) :: flxwkb

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
        intent(inout) :: flxwkb

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
          mpi_double_precision, dest, tag, xSliceLeft_recv(1, 0, 0), &
          recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xSliceLeft_send(1, 0, 0), sendcount, &
          mpi_double_precision, dest, tag, xSliceRight_recv(1, 0, 0), &
          recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          sts_right, ierror)

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
          mpi_double_precision, dest, tag, ySliceBack_recv(0, 1, 0), &
          recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(ySliceBack_send(0, 1, 0), sendcount, &
          mpi_double_precision, dest, tag, ySliceForw_recv(0, 1, 0), &
          recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          sts_forw, ierror)

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
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), &
        intent(inout) :: ray

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: ix, irprop, nrmxrl, nrmxll, nrmaxl, nrmaxr, kz, jy, iRay
    integer :: ix0

    real :: xr, xrt

    ! auxiliary fields
    real, dimension(:, :, :, :), allocatable :: xSliceLeft_send, &
        xSliceRight_send
    real, dimension(:, :, :, :), allocatable :: xSliceLeft_recv, &
        xSliceRight_recv

    if(xBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_rayvol_x must be  &
          periodic!'
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
          ierror)
      call mpi_bcast(nrmaxl, 1, mpi_integer, root, comm, ierror)

      call mpi_reduce(nrmxrl, nrmaxr, 1, mpi_integer, mpi_max, root, comm, &
          ierror)
      call mpi_bcast(nrmaxr, 1, mpi_integer, root, comm, ierror)

      allocate(xSliceLeft_send(nrmaxl, 1, 0:ny + 1, - 1:nz + 2))
      allocate(xSliceRight_send(nrmaxr, 1, 0:ny + 1, - 1:nz + 2))

      allocate(xSliceLeft_recv(nrmaxr, 1, 0:ny + 1, - 1:nz + 2))
      allocate(xSliceRight_recv(nrmaxl, 1, 0:ny + 1, - 1:nz + 2))

      do irprop = 1, 17
        ! read slice into contiguous array

        do kz = - 1, nz + 2
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
                  !elseif (irprop == 15) then
                  ! xSliceLeft_send (iRay,1,jy,kz) &
                  ! = ray(iRay,1,jy,kz)%area_xk
                  !elseif (irprop == 16) then
                  ! xSliceLeft_send (iRay,1,jy,kz) &
                  ! = ray(iRay,1,jy,kz)%area_yl
                  !elseif (irprop == 17) then
                  ! xSliceLeft_send (iRay,1,jy,kz) &
                  ! = ray(iRay,1,jy,kz)%area_zm
                end if
              end do
            end if
          end do
        end do

        do kz = - 1, nz + 2
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
                      kz)%omega
                elseif(irprop == 8) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, &
                      kz)%dkray
                elseif(irprop == 9) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, &
                      kz)%dlray
                elseif(irprop == 10) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, &
                      kz)%dmray
                elseif(irprop == 11) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, &
                      kz)%dxray
                elseif(irprop == 12) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, &
                      kz)%dyray
                elseif(irprop == 13) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, &
                      kz)%dzray
                elseif(irprop == 14) then
                  xSliceRight_send(iRay, 1, jy, kz) = ray(iRay, nx, jy, kz)%dens
                  !elseif (irprop == 15) then
                  ! xSliceRight_send(iRay,1,jy,kz) &
                  ! = ray(iRay,nx,jy,kz)%area_xk
                  !elseif (irprop == 16) then
                  ! xSliceRight_send(iRay,1,jy,kz) &
                  ! = ray(iRay,nx,jy,kz)%area_yl
                  !elseif (irprop == 17) then
                  ! xSliceRight_send(iRay,1,jy,kz) &
                  ! = ray(iRay,nx,jy,kz)%area_zm
                end if
              end do
            end if
          end do
        end do

        if(irprop < 15) then
          ! left -> right
          sendcount = nrmaxr * (ny + 2) * (nz + 4)
          recvcount = sendcount
          source = left
          dest = right
          tag = 100

          call mpi_sendrecv(xSliceRight_send(1, 1, 0, - 1), sendcount, &
              mpi_double_precision, dest, tag, xSliceLeft_recv(1, 1, 0, - 1), &
              recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              sts_left, ierror)

          ! right -> left
          sendcount = nrmaxl * (ny + 2) * (nz + 4)
          recvcount = sendcount
          source = right
          dest = left
          tag = 100

          call mpi_sendrecv(xSliceLeft_send(1, 1, 0, - 1), sendcount, &
              mpi_double_precision, dest, tag, xSliceRight_recv(1, 1, 0, - 1), &
              recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              sts_right, ierror)
        end if

        ! write auxiliary slice to ray field

        ! right halos

        do kz = - 1, nz + 2
          do jy = 0, ny + 1
            if(nRay(nx + 1, jy, kz) > 0) then
              do iRay = 1, nRay(nx + 1, jy, kz)
                if(irprop == 1) then
                  ray(iRay, nx + 1, jy, kz)%x = xSliceRight_recv(iRay, 1, jy, &
                      kz)
                elseif(irprop == 2) then
                  ray(iRay, nx + 1, jy, kz)%y = xSliceRight_recv(iRay, 1, jy, &
                      kz)
                elseif(irprop == 3) then
                  ray(iRay, nx + 1, jy, kz)%z = xSliceRight_recv(iRay, 1, jy, &
                      kz)
                elseif(irprop == 4) then
                  ray(iRay, nx + 1, jy, kz)%k = xSliceRight_recv(iRay, 1, jy, &
                      kz)
                elseif(irprop == 5) then
                  ray(iRay, nx + 1, jy, kz)%l = xSliceRight_recv(iRay, 1, jy, &
                      kz)
                elseif(irprop == 6) then
                  ray(iRay, nx + 1, jy, kz)%m = xSliceRight_recv(iRay, 1, jy, &
                      kz)
                elseif(irprop == 7) then
                  ray(iRay, nx + 1, jy, kz)%omega = xSliceRight_recv(iRay, 1, &
                      jy, kz)
                elseif(irprop == 8) then
                  ray(iRay, nx + 1, jy, kz)%dkray = xSliceRight_recv(iRay, 1, &
                      jy, kz)
                elseif(irprop == 9) then
                  ray(iRay, nx + 1, jy, kz)%dlray = xSliceRight_recv(iRay, 1, &
                      jy, kz)
                elseif(irprop == 10) then
                  ray(iRay, nx + 1, jy, kz)%dmray = xSliceRight_recv(iRay, 1, &
                      jy, kz)
                elseif(irprop == 11) then
                  ray(iRay, nx + 1, jy, kz)%dxray = xSliceRight_recv(iRay, 1, &
                      jy, kz)
                elseif(irprop == 12) then
                  ray(iRay, nx + 1, jy, kz)%dyray = xSliceRight_recv(iRay, 1, &
                      jy, kz)
                elseif(irprop == 13) then
                  ray(iRay, nx + 1, jy, kz)%dzray = xSliceRight_recv(iRay, 1, &
                      jy, kz)
                elseif(irprop == 14) then
                  ray(iRay, nx + 1, jy, kz)%dens = xSliceRight_recv(iRay, 1, &
                      jy, kz)
                  !elseif (irprop == 15) then
                  ! ray(iRay,nx+1,jy,kz)%area_xk &
                  ! = xSliceRight_recv(iRay,1,jy,kz)
                  !elseif (irprop == 16) then
                  ! ray(iRay,nx+1,jy,kz)%area_yl &
                  ! = xSliceRight_recv(iRay,1,jy,kz)
                  !elseif (irprop == 17) then
                  ! ray(iRay,nx+1,jy,kz)%area_zm &
                  ! = xSliceRight_recv(iRay,1,jy,kz)
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
            elseif(irprop > 14) then
              if(nRay(nx + 1, jy, kz) > 0) then
                do iRay = 1, nRay(nx + 1, jy, kz)
                  if(irprop == 15) then
                    ray(iRay, nx + 1, jy, kz)%area_xk = ray(iRay, nx + 1, jy, &
                        kz)%dxray * ray(iRay, nx + 1, jy, kz)%dkray
                  else if(irprop == 16) then
                    ray(iRay, nx + 1, jy, kz)%area_yl = ray(iRay, nx + 1, jy, &
                        kz)%dyray * ray(iRay, nx + 1, jy, kz)%dlray
                  else if(irprop == 17) then
                    ray(iRay, nx + 1, jy, kz)%area_zm = ray(iRay, nx + 1, jy, &
                        kz)%dzray * ray(iRay, nx + 1, jy, kz)%dmray
                  end if
                end do
              end if
            end if
          end do
        end do

        ! left halos

        do kz = - 1, nz + 2
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
                  !elseif (irprop == 15) then
                  ! ray(iRay,0,jy,kz)%area_xk &
                  ! =  xSliceLeft_recv(iRay,1,jy,kz)
                  !elseif (irprop == 16) then
                  ! ray(iRay,0,jy,kz)%area_yl &
                  ! =  xSliceLeft_recv(iRay,1,jy,kz)
                  !elseif (irprop == 17) then
                  ! ray(iRay,0,jy,kz)%area_zm &
                  ! =  xSliceLeft_recv(iRay,1,jy,kz)
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
            elseif(irprop > 14) then
              if(nRay(0, jy, kz) > 0) then
                do iRay = 1, nRay(0, jy, kz)
                  if(irprop == 15) then
                    ray(iRay, 0, jy, kz)%area_xk = ray(iRay, 0, jy, kz)%dxray &
                        * ray(iRay, 0, jy, kz)%dkray
                  else if(irprop == 16) then
                    ray(iRay, 0, jy, kz)%area_yl = ray(iRay, 0, jy, kz)%dyray &
                        * ray(iRay, 0, jy, kz)%dlray
                  else if(irprop == 17) then
                    ray(iRay, 0, jy, kz)%area_zm = ray(iRay, 0, jy, kz)%dzray &
                        * ray(iRay, 0, jy, kz)%dmray
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

      do kz = - 1, nz + 2
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

      do kz = - 1, nz + 2
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
    integer, dimension(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), &
        intent(inout) :: irsl, irsr

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: ix, irprop, nrmxrl, nrmxll, nrmaxl, nrmaxr
    integer :: jy, kz, iRay

    ! auxiliary fields
    integer, dimension(:, :, :, :), allocatable :: xSliceLeft_send, &
        xSliceRight_send
    integer, dimension(:, :, :, :), allocatable :: xSliceLeft_recv, &
        xSliceRight_recv

    if(xBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_irshift_x must be  &
          periodic!'
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
          ierror)
      call mpi_bcast(nrmaxl, 1, mpi_integer, root, comm, ierror)

      call mpi_reduce(nrmxrl, nrmaxr, 1, mpi_integer, mpi_max, root, comm, &
          ierror)
      call mpi_bcast(nrmaxr, 1, mpi_integer, root, comm, ierror)

      allocate(xSliceLeft_send(nrmaxl, 1, 0:ny + 1, - 1:nz + 2))
      allocate(xSliceRight_send(nrmaxr, 1, 0:ny + 1, - 1:nz + 2))

      allocate(xSliceLeft_recv(nrmaxr, 1, 0:ny + 1, - 1:nz + 2))
      allocate(xSliceRight_recv(nrmaxl, 1, 0:ny + 1, - 1:nz + 2))

      do irprop = 1, 2
        ! read slice into contiguous array

        do kz = - 1, nz + 2
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

        do kz = - 1, nz + 2
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
        sendcount = nrmaxr * (ny + 2) * (nz + 4)
        recvcount = sendcount
        source = left
        dest = right
        tag = 100

        call mpi_sendrecv(xSliceRight_send(1, 1, 0, - 1), sendcount, &
            mpi_integer, dest, tag, xSliceLeft_recv(1, 1, 0, - 1), recvcount, &
            mpi_integer, source, mpi_any_tag, comm, sts_left, ierror)

        ! right -> left
        sendcount = nrmaxl * (ny + 2) * (nz + 4)
        recvcount = sendcount
        source = right
        dest = left
        tag = 100

        call mpi_sendrecv(xSliceLeft_send(1, 1, 0, - 1), sendcount, &
            mpi_integer, dest, tag, xSliceRight_recv(1, 1, 0, - 1), recvcount, &
            mpi_integer, source, mpi_any_tag, comm, sts_right, ierror)

        ! write auxiliary slice to ray field

        ! right halos
        do kz = - 1, nz + 2
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

        do kz = - 1, nz + 2
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
    integer, dimension(0:nx + 1, 0:ny + 1, - 1:nz + 2), intent(inout) :: nshl, &
        nshr

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: ix, irprop

    ! auxiliary fields
    integer, dimension(1, 0:ny + 1, - 1:nz + 2) :: xSliceLeft_send, &
        xSliceRight_send
    integer, dimension(1, 0:ny + 1, - 1:nz + 2) :: xSliceLeft_recv, &
        xSliceRight_recv

    if(xBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_nshift_x must be  &
          periodic!'
    end if

    if(idim > 1) then
      ! more than 1 cpu in x direction

      call mpi_cart_shift(comm, 0, 1, left, right, ierror)

      ! slice size
      sendcount = (ny + 2) * (nz + 4)
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

        call mpi_sendrecv(xSliceRight_send(1, 0, - 1), sendcount, mpi_integer, &
            dest, tag, xSliceLeft_recv(1, 0, - 1), recvcount, mpi_integer, &
            source, mpi_any_tag, comm, sts_left, ierror)

        ! right -> left
        source = right
        dest = left
        tag = 100

        call mpi_sendrecv(xSliceLeft_send(1, 0, - 1), sendcount, mpi_integer, &
            dest, tag, xSliceRight_recv(1, 0, - 1), recvcount, mpi_integer, &
            source, mpi_any_tag, comm, sts_right, ierror)

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
    type(rayType), dimension(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), &
        intent(inout) :: ray

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
          periodic!'
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
          ierror)
      call mpi_bcast(nrmaxb, 1, mpi_integer, root, comm, ierror)

      call mpi_reduce(nrmxfl, nrmaxf, 1, mpi_integer, mpi_max, root, comm, &
          ierror)
      call mpi_bcast(nrmaxf, 1, mpi_integer, root, comm, ierror)

      allocate(ySliceBack_send(nrmaxb, 0:nx + 1, 1, - 1:nz + 2))
      allocate(ySliceForw_send(nrmaxf, 0:nx + 1, 1, - 1:nz + 2))

      allocate(ySliceBack_recv(nrmaxf, 0:nx + 1, 1, - 1:nz + 2))
      allocate(ySliceForw_recv(nrmaxb, 0:nx + 1, 1, - 1:nz + 2))

      do irprop = 1, 17
        ! read slice into contiguous array

        do kz = - 1, nz + 2
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
                  !elseif (irprop == 15) then
                  ! ySliceBack_send(iRay,ix,1,kz) &
                  ! = ray(iRay,ix,1,kz)%area_xk
                  !elseif (irprop == 16) then
                  ! ySliceBack_send(iRay,ix,1,kz) &
                  ! = ray(iRay,ix,1,kz)%area_yl
                  !elseif (irprop == 17) then
                  ! ySliceBack_send(iRay,ix,1,kz) &
                  ! = ray(iRay,ix,1,kz)%area_zm
                end if
              end do
            end if
          end do
        end do

        do kz = - 1, nz + 2
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
                  !elseif (irprop == 15) then
                  ! ySliceForw_send(iRay,ix,1,kz) &
                  ! = ray(iRay,ix,ny,kz)%area_xk
                  !elseif (irprop == 16) then
                  ! ySliceForw_send(iRay,ix,1,kz) &
                  ! = ray(iRay,ix,ny,kz)%area_yl
                  !elseif (irprop == 17) then
                  ! ySliceForw_send(iRay,ix,1,kz) &
                  ! = ray(iRay,ix,ny,kz)%area_zm
                end if
              end do
            end if
          end do
        end do

        if(irprop < 15) then
          ! back -> forw
          sendcount = nrmaxf * (nx + 2) * (nz + 4)
          recvcount = sendcount
          source = back
          dest = forw
          tag = 100

          call mpi_sendrecv(ySliceForw_send(1, 0, 1, - 1), sendcount, &
              mpi_double_precision, dest, tag, ySliceBack_recv(1, 0, 1, - 1), &
              recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              sts_back, ierror)

          ! forw -> back
          sendcount = nrmaxb * (nx + 2) * (nz + 4)
          recvcount = sendcount
          source = forw
          dest = back
          tag = 100

          call mpi_sendrecv(ySliceBack_send(1, 0, 1, - 1), sendcount, &
              mpi_double_precision, dest, tag, ySliceForw_recv(1, 0, 1, - 1), &
              recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              sts_forw, ierror)
        end if

        ! write auxiliary slice to ray field

        ! forw halos

        do kz = - 1, nz + 2
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
                      1, kz)
                elseif(irprop == 8) then
                  ray(iRay, ix, ny + 1, kz)%dkray = ySliceForw_recv(iRay, ix, &
                      1, kz)
                elseif(irprop == 9) then
                  ray(iRay, ix, ny + 1, kz)%dlray = ySliceForw_recv(iRay, ix, &
                      1, kz)
                elseif(irprop == 10) then
                  ray(iRay, ix, ny + 1, kz)%dmray = ySliceForw_recv(iRay, ix, &
                      1, kz)
                elseif(irprop == 11) then
                  ray(iRay, ix, ny + 1, kz)%dxray = ySliceForw_recv(iRay, ix, &
                      1, kz)
                elseif(irprop == 12) then
                  ray(iRay, ix, ny + 1, kz)%dyray = ySliceForw_recv(iRay, ix, &
                      1, kz)
                elseif(irprop == 13) then
                  ray(iRay, ix, ny + 1, kz)%dzray = ySliceForw_recv(iRay, ix, &
                      1, kz)
                elseif(irprop == 14) then
                  ray(iRay, ix, ny + 1, kz)%dens = ySliceForw_recv(iRay, ix, &
                      1, kz)
                  !elseif (irprop == 15) then
                  ! ray(iRay,ix,ny+1,kz)%area_xk &
                  ! = ySliceForw_recv(iRay,ix,1,kz)
                  !elseif (irprop == 16) then
                  ! ray(iRay,ix,ny+1,kz)%area_yl &
                  ! = ySliceForw_recv(iRay,ix,1,kz)
                  !elseif (irprop == 17) then
                  ! ray(iRay,ix,ny+1,kz)%area_zm &
                  ! = ySliceForw_recv(iRay,ix,1,kz)
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
            else if(irprop > 14) then
              if(nRay(ix, ny + 1, kz) > 0) then
                do iRay = 1, nRay(ix, ny + 1, kz)
                  if(irprop == 15) then
                    ray(iRay, ix, ny + 1, kz)%area_xk = ray(iRay, ix, ny + 1, &
                        kz)%dxray * ray(iRay, ix, ny + 1, kz)%dkray
                  elseif(irprop == 16) then
                    ray(iRay, ix, ny + 1, kz)%area_yl = ray(iRay, ix, ny + 1, &
                        kz)%dyray * ray(iRay, ix, ny + 1, kz)%dlray
                  elseif(irprop == 17) then
                    ray(iRay, ix, ny + 1, kz)%area_zm = ray(iRay, ix, ny + 1, &
                        kz)%dzray * ray(iRay, ix, ny + 1, kz)%dmray
                  end if
                end do
              end if
            end if
          end do
        end do

        ! back halos

        do kz = - 1, nz + 2
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
                  !elseif (irprop == 15) then
                  ! ray(iRay,ix,0,kz)%area_xk &
                  ! = ySliceBack_recv(iRay,ix,1,kz)
                  !elseif (irprop == 16) then
                  ! ray(iRay,ix,0,kz)%area_yl &
                  ! = ySliceBack_recv(iRay,ix,1,kz)
                  !elseif (irprop == 17) then
                  ! ray(iRay,ix,0,kz)%area_zm &
                  ! = ySliceBack_recv(iRay,ix,1,kz)
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
            elseif(irprop > 14) then
              if(nRay(ix, 0, kz) > 0) then
                do iRay = 1, nRay(ix, 0, kz)
                  if(irprop == 15) then
                    ray(iRay, ix, 0, kz)%area_xk = ray(iRay, ix, 0, kz)%dxray &
                        * ray(iRay, ix, 0, kz)%dkray
                  elseif(irprop == 16) then
                    ray(iRay, ix, 0, kz)%area_yl = ray(iRay, ix, 0, kz)%dyray &
                        * ray(iRay, ix, 0, kz)%dlray
                  elseif(irprop == 17) then
                    ray(iRay, ix, 0, kz)%area_zm = ray(iRay, ix, 0, kz)%dzray &
                        * ray(iRay, ix, 0, kz)%dmray
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

      do kz = - 1, nz + 2
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

      do kz = - 1, nz + 2
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
    integer, dimension(nray_wrk, 0:nx + 1, 0:ny + 1, - 1:nz + 2), &
        intent(inout) :: irsb, irsf

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: jy, irprop, nrmxbl, nrmxfl, nrmaxb, nrmaxf
    integer :: ix, kz, iRay

    ! auxiliary fields
    integer, dimension(:, :, :, :), allocatable :: ySliceBack_send, &
        ySliceForw_send
    integer, dimension(:, :, :, :), allocatable :: ySliceBack_recv, &
        ySliceForw_recv

    if(yBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_irshift_y must be  &
          periodic!'
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
          ierror)
      call mpi_bcast(nrmaxb, 1, mpi_integer, root, comm, ierror)

      call mpi_reduce(nrmxfl, nrmaxf, 1, mpi_integer, mpi_max, root, comm, &
          ierror)
      call mpi_bcast(nrmaxf, 1, mpi_integer, root, comm, ierror)

      allocate(ySliceBack_send(nrmaxb, 0:nx + 1, 1, - 1:nz + 2))
      allocate(ySliceForw_send(nrmaxf, 0:nx + 1, 1, - 1:nz + 2))

      allocate(ySliceBack_recv(nrmaxf, 0:nx + 1, 1, - 1:nz + 2))
      allocate(ySliceForw_recv(nrmaxb, 0:nx + 1, 1, - 1:nz + 2))

      do irprop = 1, 2
        ! read slice into contiguous array

        do kz = - 1, nz + 2
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

        do kz = - 1, nz + 2
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
        sendcount = nrmaxf * (nx + 2) * (nz + 4)
        recvcount = sendcount
        source = back
        dest = forw
        tag = 100

        call mpi_sendrecv(ySliceForw_send(1, 0, 1, - 1), sendcount, &
            mpi_integer, dest, tag, ySliceBack_recv(1, 0, 1, - 1), recvcount, &
            mpi_integer, source, mpi_any_tag, comm, sts_back, ierror)

        ! forw -> back
        sendcount = nrmaxb * (nx + 2) * (nz + 4)
        recvcount = sendcount
        source = forw
        dest = back
        tag = 100

        call mpi_sendrecv(ySliceBack_send(1, 0, 1, - 1), sendcount, &
            mpi_integer, dest, tag, ySliceForw_recv(1, 0, 1, - 1), recvcount, &
            mpi_integer, source, mpi_any_tag, comm, sts_forw, ierror)

        ! write auxiliary slice to ray field

        ! forw halos

        do kz = - 1, nz + 2
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

        do kz = - 1, nz + 2
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
    integer, dimension(0:nx + 1, 0:ny + 1, - 1:nz + 2), intent(inout) :: nshb, &
        nshf

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: jy, irprop

    ! auxiliary fields
    integer, dimension(0:nx + 1, 1, - 1:nz + 2) :: ySliceBack_send, &
        ySliceForw_send
    integer, dimension(0:nx + 1, 1, - 1:nz + 2) :: ySliceBack_recv, &
        ySliceForw_recv

    if(yBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_nshift_y must be  &
          periodic!'
    end if

    if(jdim > 1) then
      ! more than 1 cpu in x direction

      call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

      ! slice size
      sendcount = (nx + 2) * (nz + 4)
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

        call mpi_sendrecv(ySliceForw_send(0, 1, - 1), sendcount, mpi_integer, &
            dest, tag, ySliceBack_recv(0, 1, - 1), recvcount, mpi_integer, &
            source, mpi_any_tag, comm, sts_back, ierror)

        ! forw -> back
        source = forw
        dest = back
        tag = 100

        call mpi_sendrecv(ySliceBack_send(0, 1, - 1), sendcount, mpi_integer, &
            dest, tag, ySliceForw_recv(0, 1, - 1), recvcount, mpi_integer, &
            source, mpi_any_tag, comm, sts_forw, ierror)

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
    integer, dimension(1, 0:ny + 1, - 1:nz + 2) :: xSliceLeft_send, &
        xSliceRight_send
    integer, dimension(1, 0:ny + 1, - 1:nz + 2) :: xSliceLeft_recv, &
        xSliceRight_recv

    if(xBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_nray_x must be  periodic!'
    end if

    if(idim > 1) then
      ! more than 1 cpu in x direction

      call mpi_cart_shift(comm, 0, 1, left, right, ierror)

      ! slice size
      sendcount = (ny + 2) * (nz + 4)
      recvcount = sendcount

      ! read slice into contiguous array

      xSliceLeft_send(1, :, :) = nRay(1, :, :)
      xSliceRight_send(1, :, :) = nRay(nx, :, :)

      ! left -> right
      source = left
      dest = right
      tag = 100

      call mpi_sendrecv(xSliceRight_send(1, 0, - 1), sendcount, mpi_integer, &
          dest, tag, xSliceLeft_recv(1, 0, - 1), recvcount, mpi_integer, &
          source, mpi_any_tag, comm, sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xSliceLeft_send(1, 0, - 1), sendcount, mpi_integer, &
          dest, tag, xSliceRight_recv(1, 0, - 1), recvcount, mpi_integer, &
          source, mpi_any_tag, comm, sts_right, ierror)

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

  end subroutine setboundary_nRay_x

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
    integer, dimension(0:nx + 1, 1, - 1:nz + 2) :: ySliceBack_send, &
        ySliceForw_send
    integer, dimension(0:nx + 1, 1, - 1:nz + 2) :: ySliceBack_recv, &
        ySliceForw_recv

    if(yBoundary /= "periodic") then
      stop 'ERROR: boundary conditions in setboundary_nRay_y must be  periodic!'
    end if

    if(jdim > 1) then
      ! more than 1 cpu in x direction

      call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

      ! slice size
      sendcount = (nx + 2) * (nz + 4)
      recvcount = sendcount

      ! read slice into contiguous array

      ySliceBack_send(:, 1, :) = nRay(:, 1, :)
      ySliceForw_send(:, 1, :) = nRay(:, ny, :)

      ! back -> forw
      source = back
      dest = forw
      tag = 100

      call mpi_sendrecv(ySliceForw_send(0, 1, - 1), sendcount, mpi_integer, &
          dest, tag, ySliceBack_recv(0, 1, - 1), recvcount, mpi_integer, &
          source, mpi_any_tag, comm, sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(ySliceBack_send(0, 1, - 1), sendcount, mpi_integer, &
          dest, tag, ySliceForw_recv(0, 1, - 1), recvcount, mpi_integer, &
          source, mpi_any_tag, comm, sts_forw, ierror)

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

  end subroutine setboundary_nRay_y

  ! ----------------------------------------------------------------------

  subroutine setup_topography_wkb

    integer :: ix, jy, kz
    integer :: ix0, jy0
    integer :: iostat, size, allocstat
    integer :: row, column

    real :: mountainHeight_wkb, mountainWidth_wkb
    real :: x_center, y_center
    real :: k_mountain_wkb

    ix0 = is + nbx - 1
    jy0 = js + nby - 1

    mountainHeight_wkb = mountainHeight_wkb_dim / lRef
    mountainWidth_wkb = mountainWidth_wkb_dim / lRef

    x_center = 0.5 * (lx(1) + lx(0))
    y_center = 0.5 * (ly(1) + ly(0))

    k_mountain_wkb = pi / mountainWidth_wkb

    topography_spectrum = 0.0
    k_spectrum = 0.0
    l_spectrum = 0.0
    do jy = 1, ny
      do ix = 1, nx
        if(mountain_case_wkb == 0) then
          exit
        else if(mountain_case_wkb == 1) then
          topography_spectrum(ix, jy) = 0.5 * mountainHeight_wkb
          k_spectrum(ix, jy) = k_mountain_wkb
          l_spectrum(ix, jy) = 0.0
        else if(mountain_case_wkb == 2) then
          topography_spectrum(ix, jy) = 0.5 * mountainHeight_wkb
          k_spectrum(ix, jy) = k_mountain_wkb / sqrt(2.0)
          l_spectrum(ix, jy) = k_mountain_wkb / sqrt(2.0)
        else if(mountain_case_wkb == 6) then
          if(abs(x(ix + ix0) - x_center) <= mountainWidth_wkb) then
            topography_spectrum(ix, jy) = 0.25 * mountainHeight_wkb * (1.0 &
                + cos(k_mountain_wkb * (x(ix + ix0) - x_center)))
            k_spectrum(ix, jy) = range_factor_wkb * k_mountain_wkb
            l_spectrum(ix, jy) = 0.0
          end if
        else if(mountain_case_wkb == 7) then
          if(abs(x(ix + ix0) - x_center) <= mountainWidth_wkb .and. abs(y(jy &
              + jy0) - y_center) <= mountainWidth_wkb) then
            topography_spectrum(ix, jy) = 0.125 * mountainHeight_wkb * (1.0 &
                + cos(k_mountain_wkb * (x(ix + ix0) - x_center))) * (1.0 &
                + cos(k_mountain_wkb * (y(jy + jy0) - y_center)))
            k_spectrum(ix, jy) = range_factor_wkb * k_mountain_wkb
            l_spectrum(ix, jy) = range_factor_wkb * k_mountain_wkb
          end if
        else if(mountain_case_wkb == 8) then
          topography_spectrum(ix, jy) = 0.5 * mountainHeight_wkb * exp(- &
              ((x(ix + ix0) - x_center) / mountainWidth_wkb) ** 2.0)
          k_spectrum(ix, jy) = 0.5 * range_factor_wkb * k_mountain_wkb
          l_spectrum(ix, jy) = 0.0
        else if(mountain_case_wkb == 9) then
          topography_spectrum(ix, jy) = 0.5 * mountainHeight_wkb * exp(- &
              ((x(ix + ix0) - x_center) / mountainWidth_wkb) ** 2.0 - ((y(jy &
              + jy0) - y_center) / mountainWidth_wkb) ** 2.0)
          k_spectrum(ix, jy) = 0.5 * range_factor_wkb * k_mountain_wkb
          l_spectrum(ix, jy) = 0.5 * range_factor_wkb * k_mountain_wkb
        else
          if(master) stop "Mountain case not defined!"
        end if
      end do
    end do

    if(mountain_case_wkb == 0) then
      call read_wkb_topography
    end if

    if(topographyTime_wkb > 0.0) then
      allocate(final_topography_spectrum(1:nx, 1:ny), stat = allocstat)
      if(allocstat /= 0) stop "setup_topography_wkb: could not allocate &
          final_topography_spectrum"

      final_topography_spectrum = topography_spectrum
      topography_spectrum = 0.0
    end if

    if(long_scaling) then
      if(master) then
        open(42, file = "long_scaling.txt")
        size = 0
        do
          read(42, *, iostat = iostat)
          if(iostat /= 0) exit
          size = size + 1
        end do
        close(42)
      end if
      call mpi_bcast(size, 1, mpi_integer, root, mpi_comm_world, ierror)

      allocate(long_tree(1:size, 1:4), stat = allocstat)
      if(allocstat /= 0) stop "setup_topography_wkb: could not allocate &
          long_tree"

      if(master) then
        open(42, file = "long_scaling.txt")
        do row = 1, size
          read(42, *) (long_tree(row, column), column = 1, 4)
        end do
        close(42)
      end if
      call mpi_bcast(long_tree, size * 4, mpi_double_precision, root, &
          mpi_comm_world, ierror)
    end if

  end subroutine setup_topography_wkb

  ! ----------------------------------------------------------------------

  subroutine update_topography_wkb(time)

    real, intent(in) :: time

    if(topographyTime_wkb <= 0.0) return

    if(any(topography_spectrum /= final_topography_spectrum)) then
      if(time < topographyTime_wkb / tRef) then
        topography_spectrum = time / topographyTime_wkb * tRef &
            * final_topography_spectrum
      else
        topography_spectrum = final_topography_spectrum
      end if
    end if

  end subroutine update_topography_wkb

  ! ----------------------------------------------------------------------

  subroutine read_wkb_topography

    real * 4, dimension(sizeX, sizeY) :: field_in
    real * 4, dimension(sizeX * nprocy, ny) :: field_mst
    real * 4, dimension(nx, ny) :: field_prc
    integer :: i_out, i_mst, i_prc, j_out, j_mst, j_prc
    integer :: irc_prc
    integer :: i, j, k

    ! Open file.
    if(master) then
      open(42, file = "wkb_topography.dat", form = "unformatted", access &
          = "direct", recl = sizeX * sizeY * sizeofreal4)
    end if

    irc_prc = 0

    ! Read data.
    do k = 1, 3
      irc_prc = irc_prc + 1
      if(master) then
        read(42, rec = irc_prc) field_in
        do j = 1, ny
          j_mst = j
          do j_prc = 1, nprocy
            j_out = ny * (j_prc - 1) + j
            do i_prc = 1, nprocx
              do i = 1, nx
                i_out = nx * (i_prc - 1) + i
                i_mst = nprocy * nx * (i_prc - 1) + (j_prc - 1) * nx + i
                field_mst(i_mst, j_mst) = field_in(i_out, j_out)
              end do
            end do
          end do
        end do
      end if

      call mpi_barrier(comm, ierror)

      do j = 1, ny
        ! Distribute data over all processors.
        call mpi_scatter(field_mst(1, j), nx, mpi_real, field_prc(1, j), nx, &
            mpi_real, 0, comm, ierror)
        do i = 1, nx
          ! Non-dimensionalize.
          if(k == 1) then
            topography_spectrum(i, j) = field_prc(i, j) / lRef
          else if(k == 2) then
            k_spectrum(i, j) = field_prc(i, j) * lRef
          else if(k == 3) then
            l_spectrum(i, j) = field_prc(i, j) * lRef
          end if
        end do
      end do
    end do

    if(master) close(unit = 42)

  end subroutine read_wkb_topography

  ! ----------------------------------------------------------------------

  function compute_long_scaling(long) result(scaling)

    real :: long
    real :: scaling

    real :: node
    logical :: leaf
    integer :: index

    node = 0.0
    leaf = .false.
    do while(.not. leaf)
      if(.not. any(long_tree(:, 1) == node)) then
        leaf = .true.
        cycle
      end if
      index = findloc(long_tree(:, 1), node, dim = 1)
      scaling = long_tree(index, 4)
      if(long <= long_tree(index, 3)) then
        node = 2.0 * node + 1.0
      else if(long >= long_tree(index, 3)) then
        node = 2.0 * node + 2.0
      end if
    end do

  end function compute_long_scaling

end module wkb_module