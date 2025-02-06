module flux_module

  use type_module
  use muscl_module
  use atmosphere_module
  use mpi

  implicit none

  private

  ! Public routines
  public :: reconstruction
  public :: massFlux
  public :: momentumFlux
  public :: init_fluxes
  public :: terminate_fluxes
  public :: setHalosOfField

  ! Internal module variables
  real, dimension(:, :, :), allocatable :: rhoBar, rhopBar, rhoOld, rhopOld
  real, dimension(:, :, :), allocatable :: uBar
  real, dimension(:, :, :), allocatable :: vBar
  real, dimension(:, :, :), allocatable :: wBar

  ! Needed for semi-implicit time scheme in TFC.
  real, dimension(:, :, :), allocatable :: uOldTFC, vOldTFC, wOldTFC

  ! Reconstructed variables (uTilde, vTilde and wTilde are the reconstructed
  ! momenta). The indices are i, j, k, dir and edge with dir=1|2|3 for
  ! reconstruction in x|y|z-direction and edge=0|1 for left|right cell edge.
  real, dimension(:, :, :, :, :), allocatable :: rhoTilde
  real, dimension(:, :, :, :, :), allocatable :: rhopTilde
  real, dimension(:, :, :, :, :), allocatable :: uTilde
  real, dimension(:, :, :, :, :), allocatable :: vTilde
  real, dimension(:, :, :, :, :), allocatable :: wTilde

  ! Public variables
  ! Needed for
  ! 1) BC correction
  ! 2) explicit boundary setting
  ! 3) update module
  public :: rhoTilde, rhopTilde
  public :: uTilde, vTilde, wTilde
  public :: rhoOld, rhopOld

  ! Needed for semi-implicit time scheme in TFC.
  public :: uOldTFC, vOldTFC, wOldTFC

  contains

  subroutine reconstruction(var, variable)

    ! MUSCL reconstruction of prognostic variables.

    ! I/O variables
    type(var_type), intent(in) :: var
    character(len = *), intent(in) :: variable

    ! Indices
    integer :: ix, jy, kz

    ! TFC variables
    real :: rhoEdgeR, rhoEdgeF, rhoEdgeU
    real :: pEdgeR, pEdgeF, pEdgeU

    select case(variable)

    case("rho")

      ! Compute \rho/P for reconstruction.
      rhoBar = 0.0
      do ix = - nbx, nx + nbx
        do jy = - nby, ny + nby
          do kz = 0, nz + 1
            if(pStratTFC(ix, jy, kz) == 0.0) then
              print *, "Error in reconstruction: pStratTFC(" &
                  &// trim_integer(ix) // "," // trim_integer(jy) // "," &
                  &// trim_integer(kz) // ") = 0"
              stop
            end if
            rhoBar(ix, jy, kz) = var%rho(ix, jy, kz) / pStratTFC(ix, jy, kz)
          end do
        end do
      end do

      call reconstruct_MUSCL(rhoBar, rhoTilde, nxx, nyy, nzz, limiterType1)

    case("rhop")

      ! Compute \rho'/P for reconstruction.
      rhopBar = 0.0
      do ix = - nbx, nx + nbx
        do jy = - nby, ny + nby
          do kz = 0, nz + 1
            if(pStratTFC(ix, jy, kz) == 0.0) then
              print *, "Error in reconstruction: pStratTFC(" &
                  &// trim_integer(ix) // "," // trim_integer(jy) // "," &
                  &// trim_integer(kz) // ") = 0"
              stop
            end if
            rhopBar(ix, jy, kz) = var%rhop(ix, jy, kz) / pStratTFC(ix, jy, kz)
          end do
        end do
      end do

      call reconstruct_MUSCL(rhopBar, rhopTilde, nxx, nyy, nzz, limiterType1)

    case("uvw")

      ! Compute \rho*u/P for reconstruction.
      do ix = - nbx, nx + nbx - 1
        do jy = - nby, ny + nby
          do kz = 0, nz + 1
            rhoEdgeR = 0.5 * (var%rho(ix, jy, kz) + var%rho(ix + 1, jy, kz) &
                &+ rhoStratTFC(ix, jy, kz) + rhoStratTFC(ix + 1, jy, kz))
            pEdgeR = 0.5 * (pStratTFC(ix, jy, kz) + pStratTFC(ix + 1, jy, kz))
            uBar(ix, jy, kz) = var%u(ix, jy, kz) * rhoEdgeR / pEdgeR
          end do
        end do
      end do

      ! Compute \rho*v/P for reconstruction.
      do ix = - nbx, nx + nbx
        do jy = - nby, ny + nby - 1
          do kz = 0, nz + 1
            rhoEdgeF = 0.5 * (var%rho(ix, jy, kz) + var%rho(ix, jy + 1, kz) &
                &+ rhoStratTFC(ix, jy, kz) + rhoStratTFC(ix, jy + 1, kz))
            pEdgeF = 0.5 * (pStratTFC(ix, jy, kz) + pStratTFC(ix, jy + 1, kz))
            vBar(ix, jy, kz) = var%v(ix, jy, kz) * rhoEdgeF / pEdgeF
          end do
        end do
      end do

      ! Compute \rho*w/P for reconstruction.
      wBar(:, :, 0:(nz + 1)) = var%w(:, :, 0:(nz + 1))
      do ix = 1, nx
        do jy = 1, ny
          do kz = 0, nz + 1
            wBar(ix, jy, kz) = vertWindTFC(ix, jy, kz, var)
          end do
        end do
      end do
      call setHalosOfField(wBar)
      do ix = - nbx, nx + nbx
        do jy = - nby, ny + nby
          do kz = 0, nz + 1
            rhoEdgeU = (jac(ix, jy, kz + 1) * (var%rho(ix, jy, kz) &
                &+ rhoStratTFC(ix, jy, kz)) + jac(ix, jy, kz) * (var%rho(ix, &
                &jy, kz + 1) + rhoStratTFC(ix, jy, kz + 1))) / (jac(ix, jy, &
                &kz) + jac(ix, jy, kz + 1))
            pEdgeU = (jac(ix, jy, kz + 1) * pStratTFC(ix, jy, kz) + jac(ix, &
                &jy, kz) * pStratTFC(ix, jy, kz + 1)) / (jac(ix, jy, kz) &
                &+ jac(ix, jy, kz + 1))
            wBar(ix, jy, kz) = wBar(ix, jy, kz) * rhoEdgeU / pEdgeU
          end do
        end do
      end do

      ! Reconstruct \rho*u/P, \rho*v/P and \rho*w/P.
      call reconstruct_MUSCL(uBar, uTilde, nxx, nyy, nzz, limiterType1)
      call reconstruct_MUSCL(vBar, vTilde, nxx, nyy, nzz, limiterType1)
      call reconstruct_MUSCL(wBar, wTilde, nxx, nyy, nzz, limiterType1)

    case default
      stop "Error in reconstruction: unknown case variable."
    end select

  end subroutine reconstruction

  !--------------------------------------------------------------------

  subroutine massFlux(vara, var, flux)

    !---------------------------------------------------------------------
    ! Computes the mass flux at all cell edges using reconstructed values.
    ! MUSCL assumes that the reconstructed densities are \rho/P.
    !---------------------------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: vara, var
    type(flux_type), intent(inout) :: flux

    integer :: i, j, k
    real :: rhoL, rhoR ! L=Left i-1/2, R=Right i+1/2
    real :: rhoB, rhoF ! B=Backward j-1/2, F=Forward j+1/2
    real :: rhoD, rhoU ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurf, vSurf, wSurf ! velocities at cell surface

    real :: rhoStratEdgeR, rhoStratEdgeF, rhoStratEdgeU
    real :: pEdgeR, pEdgeF, pEdgeU

    real :: fRho, gRho, hRho

    !-----------------------------------------
    !       Zonal rho fluxes in x: f
    !-----------------------------------------

    do k = 1, nz
      do j = 1, ny
        do i = 0, nx
          rhoStratEdgeR = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, j, &
              &k))
          pEdgeR = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i + 1, j, k))
          rhoR = rhoTilde(i + 1, j, k, 1, 0) + rhoStratEdgeR / pEdgeR
          rhoL = rhoTilde(i, j, k, 1, 1) + rhoStratEdgeR / pEdgeR

          pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, j, k) &
              &* pStratTFC(i + 1, j, k))
          uSurf = pEdgeR * vara%u(i, j, k)

          fRho = flux_muscl(uSurf, rhoL, rhoR)

          flux%rho(i, j, k, 1) = fRho
        end do
      end do
    end do

    !-----------------------------------------
    !    Meridional rho fluxes in y: g
    !-----------------------------------------

    do k = 1, nz
      do j = 0, ny
        do i = 1, nx
          rhoStratEdgeF = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, j + 1, &
              &k))
          pEdgeF = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j + 1, k))
          rhoF = rhoTilde(i, j + 1, k, 2, 0) + rhoStratEdgeF / pEdgeF
          rhoB = rhoTilde(i, j, k, 2, 1) + rhoStratEdgeF / pEdgeF

          pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j + 1, k) &
              &* pStratTFC(i, j + 1, k))
          vSurf = pEdgeF * vara%v(i, j, k)

          gRho = flux_muscl(vSurf, rhoB, rhoF)

          flux%rho(i, j, k, 2) = gRho
        end do
      end do
    end do

    !-----------------------------------------
    !      Vertical rho fluxes in z: h
    !-----------------------------------------

    do k = 0, nz
      do j = 1, ny
        do i = 1, nx
          rhoStratEdgeU = (jac(i, j, k + 1) * rhoStratTFC(i, j, k) + jac(i, j, &
              &k) * rhoStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k &
              &+ 1))
          pEdgeU = (jac(i, j, k + 1) * pStratTFC(i, j, k) + jac(i, j, k) &
              &* pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
          rhoU = rhoTilde(i, j, k + 1, 3, 0) + rhoStratEdgeU / pEdgeU
          rhoD = rhoTilde(i, j, k, 3, 1) + rhoStratEdgeU / pEdgeU

          pEdgeU = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
              &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
          wSurf = pEdgeU * vara%w(i, j, k)

          hRho = flux_muscl(wSurf, rhoD, rhoU)

          flux%rho(i, j, k, 3) = hRho
        end do
      end do
    end do

    ! --------------------------------------------
    !        Density-fluctuation fluxes
    ! --------------------------------------------

    ! Zonal rhop fluxes in x: f
    do k = 1, nz
      do j = 1, ny
        do i = 0, nx
          rhoR = rhopTilde(i + 1, j, k, 1, 0)
          rhoL = rhopTilde(i, j, k, 1, 1)

          pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, j, k) &
              &* pStratTFC(i + 1, j, k))
          uSurf = pEdgeR * vara%u(i, j, k)

          fRho = flux_muscl(uSurf, rhoL, rhoR)

          flux%rhop(i, j, k, 1) = fRho
        end do
      end do
    end do

    ! Meridional rhop fluxes in y: g
    do k = 1, nz
      do j = 0, ny
        do i = 1, nx
          rhoF = rhopTilde(i, j + 1, k, 2, 0)
          rhoB = rhopTilde(i, j, k, 2, 1)

          pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j + 1, k) &
              &* pStratTFC(i, j + 1, k))
          vSurf = pEdgeF * vara%v(i, j, k)

          gRho = flux_muscl(vSurf, rhoB, rhoF)

          flux%rhop(i, j, k, 2) = gRho
        end do
      end do
    end do

    ! Vertical rhop fluxes in z: h
    do k = 0, nz
      do j = 1, ny
        do i = 1, nx
          rhoU = rhopTilde(i, j, k + 1, 3, 0)
          rhoD = rhopTilde(i, j, k, 3, 1)

          pEdgeU = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
              &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
          wSurf = pEdgeU * vara%w(i, j, k)

          hRho = flux_muscl(wSurf, rhoD, rhoU)

          flux%rhop(i, j, k, 3) = hRho
        end do
      end do
    end do

  end subroutine massFlux

  !----------------------------------------------------------------------

  function flux_muscl(uSurf, phiUp, phiDown)

    !----------------------------
    !   upwind flux function
    !----------------------------

    ! in/out arguments
    real, intent(in) :: uSurf ! cell face value
    real, intent(in) :: phiUp, phiDown ! upwind, downwind values
    real :: flux_muscl

    if(uSurf > 0.0) then
      flux_muscl = uSurf * phiUp
    else
      flux_muscl = uSurf * phiDown
    end if

  end function flux_muscl

  !-----------------------------------------------------------------------

  subroutine momentumFlux(vara, var, flux)

    !----------------------------------------------------------------------
    ! Computes the momentum fluxes at the cell edges using reconstr. values.
    ! MUSCL assumes that the reconstructed momenta are \rho/P * \vec v.
    !----------------------------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: vara, var
    type(flux_type), intent(inout) :: flux

    ! local variables
    integer :: i, j, k

    ! uTilde at cell edges
    real :: uL, uR, vL, vR, wL, wR ! L=Left at i+1, R=Right at i
    real :: uB, uF, vB, vF, wB, wF ! B=Backward at j+1, F=Forward at j
    real :: uD, uU, vD, vU, wD, wU ! D=Downward at k+1, U=Upward at k

    ! mass-weigted potential temperature and stress tensor
    real :: pEdgeR, pREdgeR, pEdgeF, pREdgeF, pEdgeU, pREdgeU, pFEdgeR, &
        &pFEdgeF, pFedgeU, pUEdgeR, pUEdgeF, pUEdgeU
    real :: stressTens13, stressTens13R, stressTens13U, stressTens13RU, &
        &stressTens23, stressTens23F, stressTens23U, stressTens23FU

    ! upwinding
    real :: uSurf, vSurf, wSurf

    ! local flux variables
    real :: fRhoU, gRhoU, hRhoU ! rho*U momentum fluxes
    real :: fRhoV, gRhoV, hRhoV ! rho*V momentum fluxes
    real :: fRhoW, gRhoW, hRhoW ! rho*W momentum fluxes

    ! viscous fluxes
    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1) :: divU
    real :: fRhoU_visc, gRhoU_visc, hRhoU_visc
    real :: fRhoV_visc, gRhoV_visc, hRhoV_visc
    real :: fRhoW_visc, gRhoW_visc, hRhoW_visc
    real :: coef_v

    !------------------------------
    !       Fluxes for rho*u
    !------------------------------

    ! Flux fRhoU
    do k = 1, nz
      do j = 1, ny
        do i = - 1, nx
          ! The uTilde are the reconstructed specific momenta, divided by P.
          ! These are to be multiplied by the linearly interpolated velocities
          ! (times P) in order to obtain the desired momentum fluxes.

          uR = uTilde(i + 1, j, k, 1, 0)
          uL = uTilde(i, j, k, 1, 1)

          pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, j, k) &
              &* pStratTFC(i + 1, j, k))
          pREdgeR = 0.5 * (jac(i + 1, j, k) * pStratTFC(i + 1, j, k) + jac(i &
              &+ 2, j, k) * pStratTFC(i + 2, j, k))
          uSurf = 0.5 * (pEdgeR * vara%u(i, j, k) + pREdgeR * vara%u(i + 1, j, &
              &k))

          fRhoU = flux_muscl(uSurf, uL, uR)

          flux%u(i, j, k, 1) = fRhoU
        end do
      end do
    end do

    !----------------------------------------------------------

    !  Flux gRhoU
    do k = 1, nz
      do j = 0, ny
        do i = 0, nx
          ! The uTilde are the reconstructed specific momenta, divided by P.
          ! These are to be multiplied by the linearly interpolated velocities
          ! (times P) in order to obtain the desired momentum fluxes.

          uF = uTilde(i, j + 1, k, 2, 0)
          uB = uTilde(i, j, k, 2, 1)

          pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j + 1, k) &
              &* pStratTFC(i, j + 1, k))
          pREdgeF = 0.5 * (jac(i + 1, j, k) * pStratTFC(i + 1, j, k) + jac(i &
              &+ 1, j + 1, k) * pStratTFC(i + 1, j + 1, k))
          vSurf = 0.5 * (pEdgeF * vara%v(i, j, k) + pREdgeF * vara%v(i + 1, j, &
              &k))

          gRhoU = flux_muscl(vSurf, uB, uF)

          flux%u(i, j, k, 2) = gRhoU
        end do
      end do
    end do

    !---------------------------------------------------------------------

    ! Flux hRhoU
    do k = 0, nz
      do j = 1, ny
        do i = 0, nx
          ! The uTilde are the reconstructed specific momenta, divided by P.
          ! These are to be multiplied by the linearly interpolated velocities
          ! (times P) in order to obtain the desired momentum fluxes.

          uU = uTilde(i, j, k + 1, 3, 0)
          uD = uTilde(i, j, k, 3, 1)

          pEdgeU = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
              &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
          pREdgeU = jac(i + 1, j, k) * jac(i + 1, j, k + 1) * (pStratTFC(i &
              &+ 1, j, k) + pStratTFC(i + 1, j, k + 1)) / (jac(i + 1, j, k) &
              &+ jac(i + 1, j, k + 1))
          wSurf = 0.5 * (pEdgeU * vara%w(i, j, k) + pREdgeU * vara%w(i + 1, j, &
              &k))

          hRhoU = flux_muscl(wSurf, uD, uU)

          flux%u(i, j, k, 3) = hRhoU
        end do
      end do
    end do

    !------------------------------
    !      Fluxes for rho*v
    !------------------------------

    !  Flux fRhoV
    do k = 1, nz
      do j = 0, ny
        do i = 0, nx
          ! The vTilde are the reconstructed specific momenta, divided by P.
          ! These are to be multiplied by the linearly interpolated velocities
          ! (times P) in order to obtain the desired momentum fluxes.

          vR = vTilde(i + 1, j, k, 1, 0)
          vL = vTilde(i, j, k, 1, 1)

          pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, j, k) &
              &* pStratTFC(i + 1, j, k))
          pFEdgeR = 0.5 * (jac(i, j + 1, k) * pStratTFC(i, j + 1, k) + jac(i &
              &+ 1, j + 1, k) * pStratTFC(i + 1, j + 1, k))
          uSurf = 0.5 * (pEdgeR * vara%u(i, j, k) + pFEdgeR * vara%u(i, j + 1, &
              &k))

          fRhoV = flux_muscl(uSurf, vL, vR)

          flux%v(i, j, k, 1) = fRhoV
        end do
      end do
    end do

    !---------------------------------------------------------------------

    ! Flux gRhoV
    do k = 1, nz
      do j = - 1, ny
        do i = 1, nx
          ! The vTilde are the reconstructed specific momenta, divided by P.
          ! These are to be multiplied by the linearly interpolated velocities
          ! (times P) in order to obtain the desired momentum fluxes.

          vF = vTilde(i, j + 1, k, 2, 0)
          vB = vTilde(i, j, k, 2, 1)

          pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j + 1, k) &
              &* pStratTFC(i, j + 1, k))
          pFEdgeF = 0.5 * (jac(i, j + 1, k) * pStratTFC(i, j + 1, k) + jac(i, &
              &j + 2, k) * pStratTFC(i, j + 2, k))
          vSurf = 0.5 * (pEdgeF * vara%v(i, j, k) + pFEdgeF * vara%v(i, j + 1, &
              &k))

          gRhoV = flux_muscl(vSurf, vB, vF)

          flux%v(i, j, k, 2) = gRhoV
        end do
      end do
    end do

    !----------------------------------------------------------------------

    ! Flux hRhoV
    do k = 0, nz
      do j = 0, ny
        do i = 1, nx
          ! The vTilde are the reconstructed specific momenta, divided by P.
          ! These are to be multiplied by the linearly interpolated velocities
          ! (times P) in order to obtain the desired momentum fluxes.

          vU = vTilde(i, j, k + 1, 3, 0)
          vD = vTilde(i, j, k, 3, 1)

          pEdgeU = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
              &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
          pFEdgeU = jac(i, j + 1, k) * jac(i, j + 1, k + 1) * (pStratTFC(i, j &
              &+ 1, k) + pStratTFC(i, j + 1, k + 1)) / (jac(i, j + 1, k) &
              &+ jac(i, j + 1, k + 1))
          wSurf = 0.5 * (pEdgeU * vara%w(i, j, k) + pFEdgeU * vara%w(i, j + 1, &
              &k))

          hRhoV = flux_muscl(wSurf, vD, vU)

          flux%v(i, j, k, 3) = hRhoV
        end do
      end do
    end do

    !------------------------------
    !      Fluxes for rho*w
    !------------------------------

    ! Flux fRhoW
    do k = 0, nz
      do j = 1, ny
        do i = 0, nx
          ! The wTilde are the reconstructed specific momenta, divided by P.
          ! These are to be multiplied by the linearly interpolated velocities
          ! (times P) in order to obtain the desired momentum fluxes.

          wR = wTilde(i + 1, j, k, 1, 0)
          wL = wTilde(i, j, k, 1, 1)

          pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, j, k) &
              &* pStratTFC(i + 1, j, k))
          pUEdgeR = 0.5 * (jac(i, j, k + 1) * pStratTFC(i, j, k + 1) + jac(i &
              &+ 1, j, k + 1) * pStratTFC(i + 1, j, k + 1))
          uSurf = ((jac(i, j, k + 1) + jac(i + 1, j, k + 1)) * pEdgeR &
              &* vara%u(i, j, k) + (jac(i, j, k) + jac(i + 1, j, k)) * pUEdgeR &
              &* vara%u(i, j, k + 1)) / (jac(i, j, k) + jac(i + 1, j, k) &
              &+ jac(i, j, k + 1) + jac(i + 1, j, k + 1))

          fRhoW = flux_muscl(uSurf, wL, wR)

          flux%w(i, j, k, 1) = fRhoW
        end do
      end do
    end do

    !-------------------------------------------------------------------

    ! Flux gRhoW
    do k = 0, nz
      do j = 0, ny
        do i = 1, nx
          ! The wTilde are the reconstructed specific momenta, divided by P.
          ! These are to be multiplied by the linearly interpolated velocities
          ! (times P) in order to obtain the desired momentum fluxes.

          wF = wTilde(i, j + 1, k, 2, 0)
          wB = wTilde(i, j, k, 2, 1)

          pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j + 1, k) &
              &* pStratTFC(i, j + 1, k))
          pUEdgeF = 0.5 * (jac(i, j, k + 1) * pStratTFC(i, j, k + 1) + jac(i, &
              &j + 1, k + 1) * pStratTFC(i, j + 1, k + 1))
          vSurf = ((jac(i, j, k + 1) + jac(i, j + 1, k + 1)) * pEdgeF &
              &* vara%v(i, j, k) + (jac(i, j, k) + jac(i, j + 1, k)) * pUEdgeF &
              &* vara%v(i, j, k + 1)) / (jac(i, j, k) + jac(i, j + 1, k) &
              &+ jac(i, j, k + 1) + jac(i, j + 1, k + 1))

          gRhoW = flux_muscl(vSurf, wB, wF)

          flux%w(i, j, k, 2) = gRhoW
        end do
      end do
    end do

    !---------------------------------------------------------------------

    ! Flux hRhoW
    do k = - 1, nz
      do j = 1, ny
        do i = 1, nx
          ! The wTilde are the reconstructed specific momenta, divided by P.
          ! These are to be multiplied by the linearly interpolated velocities
          ! (times P) in order to obtain the desired momentum fluxes.

          wU = wTilde(i, j, k + 1, 3, 0)
          wD = wTilde(i, j, k, 3, 1)

          pEdgeU = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
              &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
          pUEdgeU = jac(i, j, k + 1) * jac(i, j, k + 2) * (pStratTFC(i, j, k &
              &+ 1) + pStratTFC(i, j, k + 2)) / (jac(i, j, k + 1) + jac(i, j, &
              &k + 2))
          wSurf = 0.5 * (pEdgeU * vara%w(i, j, k) + pUEdgeU * vara%w(i, j, k &
              &+ 1))

          hRhoW = flux_muscl(wSurf, wD, wU)

          flux%w(i, j, k, 3) = hRhoW
        end do
      end do
    end do

    !-------------------------------------------------------------------
    !                          Viscous Fluxes
    !-------------------------------------------------------------------

    if(ReInv == 0.0) return

    !------------------------------------
    !       Calculate divergence
    !------------------------------------

    do k = 1, nz
      do j = 0, ny + 1
        do i = 0, nx + 1
          uR = var%u(i, j, k)
          uL = var%u(i - 1, j, k)
          vF = var%v(i, j, k)
          vB = var%v(i, j - 1, k)
          wU = var%w(i, j, k)
          wD = var%w(i, j, k - 1)

          divU(i, j, k) = (uR - uL) / dx + (vF - vB) / dy + (wU - wD) / dz
        end do
      end do
    end do
    divU(:, :, 0) = divU(:, :, 1)
    divU(:, :, nz + 1) = divU(:, :, nz)

    !------------------------------
    !      Fluxes for rho*u
    !------------------------------

    ! Flux fRhoU
    do k = 1, nz
      do j = 1, ny
        do i = - 1, nx
          coef_v = ReInv * rhoStratTFC(i + 1, j, 1)

          fRhoU_visc = coef_v * jac(i + 1, j, k) * stressTensTFC(i + 1, j, k, &
              &1, 1, var)

          flux%u(i, j, k, 1) = flux%u(i, j, k, 1) - fRhoU_visc
        end do
      end do
    end do

    ! Flux gRhoU
    do k = 1, nz
      do j = 0, ny
        do i = 0, nx
          coef_v = ReInv * 0.25 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i + 1, &
              &j, 1) + rhoStratTFC(i, j + 1, 1) + rhoStratTFC(i + 1, j + 1, 1))

          gRhoU_visc = coef_v * 0.25 * (jac(i, j, k) * stressTensTFC(i, j, k, &
              &1, 2, var) + jac(i + 1, j, k) * stressTensTFC(i + 1, j, k, 1, &
              &2, var) + jac(i, j + 1, k) * stressTensTFC(i, j + 1, k, 1, 2, &
              &var) + jac(i + 1, j + 1, k) * stressTensTFC(i + 1, j + 1, k, 1, &
              &2, var))

          flux%u(i, j, k, 2) = flux%u(i, j, k, 2) - gRhoU_visc
        end do
      end do
    end do

    ! Flux hRhoU
    do k = 0, nz
      do j = 1, ny
        do i = 0, nx
          coef_v = ReInv * 0.5 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i + 1, j, &
              &1))

          stressTens13 = met(i, j, k, 1, 3) * stressTensTFC(i, j, k, 1, 1, &
              &var) + met(i, j, k, 2, 3) * stressTensTFC(i, j, k, 1, 2, var) &
              &+ stressTensTFC(i, j, k, 1, 3, var) / jac(i, j, k)
          stressTens13R = met(i + 1, j, k, 1, 3) * stressTensTFC(i + 1, j, k, &
              &1, 1, var) + met(i + 1, j, k, 2, 3) * stressTensTFC(i + 1, j, &
              &k, 1, 2, var) + stressTensTFC(i + 1, j, k, 1, 3, var) / jac(i &
              &+ 1, j, k)
          stressTens13U = met(i, j, k + 1, 1, 3) * stressTensTFC(i, j, k + 1, &
              &1, 1, var) + met(i, j, k + 1, 2, 3) * stressTensTFC(i, j, k &
              &+ 1, 1, 2, var) + stressTensTFC(i, j, k + 1, 1, 3, var) &
              &/ jac(i, j, k + 1)
          stressTens13RU = met(i + 1, j, k + 1, 1, 3) * stressTensTFC(i + 1, &
              &j, k + 1, 1, 1, var) + met(i + 1, j, k + 1, 2, 3) &
              &* stressTensTFC(i + 1, j, k + 1, 1, 2, var) + stressTensTFC(i &
              &+ 1, j, k + 1, 1, 3, var) / jac(i + 1, j, k + 1)
          hRhoU_visc = coef_v * 0.5 * (jac(i, j, k) * jac(i, j, k + 1) &
              &* (stressTens13 + stressTens13U) / (jac(i, j, k) + jac(i, j, k &
              &+ 1)) + jac(i + 1, j, k) * jac(i + 1, j, k + 1) &
              &* (stressTens13R + stressTens13RU) / (jac(i + 1, j, k) + jac(i &
              &+ 1, j, k + 1)))

          flux%u(i, j, k, 3) = flux%u(i, j, k, 3) - hRhoU_visc
        end do
      end do
    end do

    !------------------------------
    !      Fluxes for rho*v
    !------------------------------

    ! Flux fRhoV
    do k = 1, nz
      do j = 0, ny
        do i = 0, nx
          coef_v = ReInv * 0.25 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i + 1, &
              &j, 1) + rhoStratTFC(i, j + 1, 1) + rhoStratTFC(i + 1, j + 1, 1))

          fRhoV_visc = coef_v * 0.25 * (jac(i, j, k) * stressTensTFC(i, j, k, &
              &2, 1, var) + jac(i + 1, j, k) * stressTensTFC(i + 1, j, k, 2, &
              &1, var) + jac(i, j + 1, k) * stressTensTFC(i, j + 1, k, 2, 1, &
              &var) + jac(i + 1, j + 1, k) * stressTensTFC(i + 1, j + 1, k, 2, &
              &1, var))

          flux%v(i, j, k, 1) = flux%v(i, j, k, 1) - fRhoV_visc
        end do
      end do
    end do

    ! Flux gRhoV
    do k = 1, nz
      do j = - 1, ny
        do i = 1, nx
          coef_v = ReInv * rhoStratTFC(i, j + 1, 1)

          gRhoV_visc = coef_v * jac(i, j + 1, k) * stressTensTFC(i, j + 1, k, &
              &2, 2, var)

          flux%v(i, j, k, 2) = flux%v(i, j, k, 2) - gRhoV_visc
        end do
      end do
    end do

    ! Flux hRhoV
    do k = 0, nz
      do j = 0, ny
        do i = 1, nx
          coef_v = ReInv * 0.5 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i, j + 1, &
              &1))

          stressTens23 = met(i, j, k, 1, 3) * stressTensTFC(i, j, k, 2, 1, &
              &var) + met(i, j, k, 2, 3) * stressTensTFC(i, j, k, 2, 2, var) &
              &+ stressTensTFC(i, j, k, 2, 3, var) / jac(i, j, k)
          stressTens23F = met(i, j + 1, k, 1, 3) * stressTensTFC(i, j + 1, k, &
              &2, 1, var) + met(i, j + 1, k, 2, 3) * stressTensTFC(i, j + 1, &
              &k, 2, 2, var) + stressTensTFC(i, j + 1, k, 2, 3, var) / jac(i, &
              &j + 1, k)
          stressTens23U = met(i, j, k + 1, 1, 3) * stressTensTFC(i, j, k + 1, &
              &2, 1, var) + met(i, j, k + 1, 2, 3) * stressTensTFC(i, j, k &
              &+ 1, 2, 2, var) + stressTensTFC(i, j, k + 1, 2, 3, var) &
              &/ jac(i, j, k + 1)
          stressTens23FU = met(i, j + 1, k + 1, 1, 3) * stressTensTFC(i, j &
              &+ 1, k + 1, 2, 1, var) + met(i, j + 1, k + 1, 2, 3) &
              &* stressTensTFC(i, j + 1, k + 1, 2, 2, var) + stressTensTFC(i, &
              &j + 1, k + 1, 2, 3, var) / jac(i, j + 1, k + 1)
          hRhoV_visc = coef_v * 0.5 * (jac(i, j, k) * jac(i, j, k + 1) &
              &* (stressTens23 + stressTens23U) / (jac(i, j, k) + jac(i, j, k &
              &+ 1)) + jac(i, j + 1, k) * jac(i, j + 1, k + 1) &
              &* (stressTens23F + stressTens23FU) / (jac(i, j + 1, k) + jac(i, &
              &j + 1, k + 1)))

          flux%v(i, j, k, 3) = flux%v(i, j, k, 3) - hRhoV_visc
        end do
      end do
    end do

    !------------------------------
    !      Fluxes for rho*w
    !------------------------------

    ! Flux fRhoW
    do k = 0, nz
      do j = 1, ny
        do i = 0, nx
          coef_v = ReInv * 0.5 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i + 1, j, &
              &1))

          fRhoW_visc = coef_v * 0.5 * (jac(i, j, k) * jac(i, j, k + 1) &
              &* (stressTensTFC(i, j, k, 3, 1, var) + stressTensTFC(i, j, k &
              &+ 1, 3, 1, var)) / (jac(i, j, k) + jac(i, j, k + 1)) + jac(i &
              &+ 1, j, k) * jac(i + 1, j, k + 1) * (stressTensTFC(i + 1, j, k, &
              &3, 1, var) + stressTensTFC(i + 1, j, k + 1, 3, 1, var)) &
              &/ (jac(i + 1, j, k) + jac(i + 1, j, k + 1)))

          flux%w(i, j, k, 1) = flux%w(i, j, k, 1) - fRhoW_visc
        end do
      end do
    end do

    ! Flux gRhoW
    do k = 0, nz
      do j = 0, ny
        do i = 1, nx
          coef_v = ReInv * 0.5 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i, j + 1, &
              &1))

          gRhoW_visc = coef_v * 0.5 * (jac(i, j, k) * jac(i, j, k + 1) &
              &* (stressTensTFC(i, j, k, 3, 1, var) + stressTensTFC(i, j, k &
              &+ 1, 3, 1, var)) / (jac(i, j, k) + jac(i, j, k + 1)) + jac(i, j &
              &+ 1, k) * jac(i, j + 1, k + 1) * (stressTensTFC(i, j + 1, k, 3, &
              &1, var) + stressTensTFC(i, j + 1, k + 1, 3, 1, var)) / (jac(i, &
              &j + 1, k) + jac(i, j + 1, k + 1)))

          flux%w(i, j, k, 2) = flux%w(i, j, k, 2) - gRhoW_visc
        end do
      end do
    end do

    ! Flux hRhoW
    do k = - 1, nz
      do j = 1, ny
        do i = 1, nx
          coef_v = ReInv * rhoStratTFC(i, j, 1)

          hRhoW_visc = coef_v * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
              &* stressTensTFC(i, j, k + 1, 3, 1, var) + jac(i, j, k + 1) &
              &* met(i, j, k + 1, 2, 3) * stressTensTFC(i, j, k + 1, 3, 2, &
              &var) + stressTensTFC(i, j, k + 1, 3, 3, var))

          flux%w(i, j, k, 3) = flux%w(i, j, k, 3) - hRhoW_visc
        end do
      end do
    end do

  end subroutine momentumFlux

  !---------------------------------------------------------------------------

  subroutine init_fluxes

    !---------------------------------------------
    ! Allocate flux module variables
    !---------------------------------------------

    integer :: allocstat

    allocate(rhoBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhoBar"

    allocate(rhopBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhopBar"

    allocate(rhoOld(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "init_fluxes: alloc of rhoOld failed"

    allocate(rhopOld(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "init_fluxes: alloc of rhopOld failed"

    allocate(uBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "init_fluxes: could not allocate uBar"

    allocate(vBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate vBar"

    allocate(wBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate wBar"

    allocate(rhoTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
        &0:1), stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhoTilde"

    allocate(rhopTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
        &0:1), stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhopTilde"

    allocate(uTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, 0:1), &
        &stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate uTilde"

    allocate(vTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, 0:1), &
        &stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate vTilde"

    allocate(wTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, 0:1), &
        &stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate wTilde"

    allocate(uOldTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        &+ nbz)), stat = allocstat)
    if(allocstat /= 0) stop "init_fluxes: alloc of uOldTFC failed"

    allocate(vOldTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        &+ nbz)), stat = allocstat)
    if(allocstat /= 0) stop "init_fluxes: alloc of vOldTFC failed"

    allocate(wOldTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        &+ nbz)), stat = allocstat)
    if(allocstat /= 0) stop "init_fluxes: alloc of wOldTFC failed"

  end subroutine init_fluxes

  !---------------------------------------------------------------------------

  subroutine terminate_fluxes

    !-----------------------------------
    ! Deallocate flux module variables
    !-----------------------------------

    integer :: allocstat

    deallocate(rhoBar, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate rhoBar"

    deallocate(rhopBar, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate rhopBar"

    deallocate(rhoOld, stat = allocstat)
    if(allocstat /= 0) stop "terminate_fluxes: dealloc of rhoOld failed"

    deallocate(rhopOld, stat = allocstat)
    if(allocstat /= 0) stop "terminate_fluxes: dealloc of rhopOld failed"

    deallocate(uBar, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate uBar"

    deallocate(vBar, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate vBar"

    deallocate(wBar, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate wBar"

    deallocate(rhoTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate rhoTilde"

    deallocate(rhopTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate rhopTilde"

    deallocate(uTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate uTilde"

    deallocate(vTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate vTilde"

    deallocate(wTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate wTilde"

    deallocate(uOldTFC, stat = allocstat)
    if(allocstat /= 0) stop "terminate_fluxes: dealloc of uOldTFC failed"

    deallocate(vOldTFC, stat = allocstat)
    if(allocstat /= 0) stop "terminate_fluxes: dealloc of vOldTFC failed"

    deallocate(wOldTFC, stat = allocstat)
    if(allocstat /= 0) stop "terminate_fluxes: dealloc of wOldTFC failed"

  end subroutine terminate_fluxes

  ! ----------------------------------------------------

  subroutine setHalosOfField(field)

    !-------------------------------
    !  set values in halo cells
    !-------------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: field

    ! auxiliary fields
    real, dimension(nbx, - nby:ny + nby, - nbz:nz + nbz) :: xSliceLeft_send, &
        &xSliceRight_send
    real, dimension(nbx, - nby:ny + nby, - nbz:nz + nbz) :: xSliceLeft_recv, &
        &xSliceRight_recv

    real, dimension(- nbx:nx + nbx, nby, - nbz:nz + nbz) :: ySliceBack_send, &
        &ySliceForw_send
    real, dimension(- nbx:nx + nbx, nby, - nbz:nz + nbz) :: ySliceBack_recv, &
        &ySliceForw_recv

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: i, j, k
    integer :: i0, j0, k0

    !-----------------------------
    !     find neighbour procs
    !-----------------------------

    if(idim > 1) call mpi_cart_shift(comm, 0, 1, left, right, ierror)
    if(jdim > 1) call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

    !------------------------------
    !          x-direction
    !------------------------------

    if(idim > 1) then

      ! slice size
      sendcount = nbx * (ny + 2 * nby + 1) * (nz + 2 * nbz + 1)
      recvcount = sendcount

      ! read slice into contiguous array
      do i = 1, nbx
        xSliceLeft_send(i, :, :) = field(i, :, :)
        xSliceRight_send(i, :, :) = field(nx - nbx + i, :, :)
      end do

      ! left -> right
      source = left
      dest = right
      tag = 100

      i0 = 1; j0 = - nby; k0 = - nbz

      call mpi_sendrecv(xSliceRight_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, xSliceLeft_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xSliceLeft_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, xSliceRight_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_right, ierror)

      ! write auxiliary slice to var field
      do i = 1, nbx
        ! right halos
        field(nx + i, :, :) = xSliceRight_recv(i, :, :)
        ! left halos
        field(- nbx + i, :, :) = xSliceLeft_recv(i, :, :)
      end do

    else

      do i = 1, nbx
        field(nx + i, :, :) = field(i, :, :)
        field(- i + 1, :, :) = field(nx - i + 1, :, :)
      end do

    end if

    !------------------------------
    !          y-direction
    !------------------------------

    if(jdim > 1) then

      ! slice size
      sendcount = nby * (nx + 2 * nbx + 1) * (nz + 2 * nbz + 1)
      recvcount = sendcount

      ! read slice into contiguous array
      do j = 1, nby
        ySliceBack_send(:, j, :) = field(:, j, :)
        ySliceForw_send(:, j, :) = field(:, ny - nby + j, :)
      end do

      ! back -> forw
      source = back
      dest = forw
      tag = 100

      i0 = - nbx; j0 = 1; k0 = - nbz

      call mpi_sendrecv(ySliceForw_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, ySliceBack_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(ySliceBack_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, ySliceForw_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_forw, ierror)

      ! write auxiliary slice to var field
      do j = 1, nby
        ! right halos
        field(:, ny + j, :) = ySliceForw_recv(:, j, :)
        ! left halos
        field(:, - nby + j, :) = ySliceBack_recv(:, j, :)
      end do

    else

      do j = 1, nby
        field(:, ny + j, :) = field(:, j, :)
        field(:, - j + 1, :) = field(:, ny - j + 1, :)
      end do

    end if

  end subroutine setHalosOfField

end module flux_module
