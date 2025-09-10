module flux_module

  use type_module
  use muscl_module
  use atmosphere_module
  use sizeof_module
  use mpi

  implicit none

  private

  ! Public routines
  public :: reconstruction
  public :: massFlux
  public :: momentumFlux
  public :: volumeForce
  public :: init_fluxes
  public :: terminate_fluxes
  public :: tracerFlux
  public :: iceFlux

  ! TFC routines
  public :: setHalosOfField

  ! Internal module variables
  real, dimension(:, :, :), allocatable :: rhoBar, rhopBar, rhoOld, rhopOld
  real, dimension(:, :, :), allocatable :: uBar
  real, dimension(:, :, :), allocatable :: vBar
  real, dimension(:, :, :), allocatable :: wBar
  real, dimension(:, :, :), allocatable :: PBar
  real, dimension(:, :, :), allocatable :: nAerBar
  real, dimension(:, :, :), allocatable :: nIceBar
  real, dimension(:, :, :), allocatable :: qIceBar
  real, dimension(:, :, :), allocatable :: qvBar
  real, dimension(:, :, :), allocatable :: tracerBar
  real, dimension(:, :, :), allocatable :: IceBar

  ! Needed for semi-implicit time scheme in TFC.
  real, dimension(:, :, :), allocatable :: uOldTFC, vOldTFC, wOldTFC

  ! Needed for compressible explicit Euler.
  real, dimension(:, :, :), allocatable :: pinew

  ! Reconstructed variables (uTilde, vTilde and wTilde are the reconstructed
  ! momenta). The indices are i, j, k, dir and edge with dir=1|2|3 for
  ! reconstruction in x|y|z-direction and edge=0|1 for left|right cell edge.
  real, dimension(:, :, :, :, :), allocatable :: rhoTilde
  real, dimension(:, :, :, :, :), allocatable :: rhoTilde_mom
  real, dimension(:, :, :, :, :), allocatable :: rhopTilde
  real, dimension(:, :, :, :, :), allocatable :: uTilde
  real, dimension(:, :, :, :, :), allocatable :: vTilde
  real, dimension(:, :, :, :, :), allocatable :: wTilde
  real, dimension(:, :, :, :, :), allocatable :: PTilde
  real, dimension(:, :, :, :, :), allocatable :: nAerTilde
  real, dimension(:, :, :, :, :), allocatable :: nIceTilde
  real, dimension(:, :, :, :, :), allocatable :: qIceTilde
  real, dimension(:, :, :, :, :), allocatable :: qvTilde
  real, dimension(:, :, :, :, :), allocatable :: tracerTilde

  ! Public variables
  ! Needed for
  ! 1) BC correction
  ! 2) explicit boundary setting
  ! 3) update module
  public :: rhoTilde, rhopTilde, rhoTilde_mom, PTilde
  public :: uTilde, vTilde, wTilde
  public :: rhoOld, rhopOld
  public :: nIceTilde, qIceTilde, qvTilde, nAerTilde
  public :: tracerTilde

  ! TFC FJ
  ! Needed for semi-implicit time scheme in TFC.
  public :: uOldTFC, vOldTFC, wOldTFC

  ! SK Needed for compressible explicit Euler
  public :: pinew

  contains

  subroutine reconstruction(var, variable)

    ! MUSCL reconstruction of prognostic variables.

    ! I/O variables
    type(var_type), intent(in) :: var
    character(len = *), intent(in) :: variable

    ! Indices
    integer :: ix, jy, kz
    integer :: iVar

    ! TFC variables
    real :: rhoEdgeR, rhoEdgeF, rhoEdgeU
    real :: pEdgeR, pEdgeF, pEdgeU

    select case(variable)

    case("rho")

      ! Compute \rho/P for reconstruction.
      rhoBar = 0.0
      if(topography) then
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
      else
        do kz = 0, nz + 1
          if(Pstrat(kz) == 0.0) then
            print *, "Error in reconstruction: pStrat(" // trim_integer(kz) &
                &// ") = 0"
            stop
          end if
          rhoBar(:, :, kz) = (var%rho(:, :, kz)) / Pstrat(kz)
        end do
      end if
      call reconstruct_MUSCL(rhoBar, rhoTilde, nxx, nyy, nzz, limiterType1)

    case("rhop")

      ! Compute \rho'/P for reconstruction.
      rhopBar = 0.0
      if(topography) then
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
      else
        do kz = 0, nz + 1
          if(Pstrat(kz) == 0.0) then
            print *, "Error in reconstruction: pStrat(" // trim_integer(kz) &
                &// ") = 0"
            stop
          end if
          rhopBar(:, :, kz) = (var%rhop(:, :, kz)) / Pstrat(kz)
        end do
      end if
      call reconstruct_MUSCL(rhopBar, rhopTilde, nxx, nyy, nzz, limiterType1)

    case("uvw")

      ! Compute \rho*u/P for reconstruction.
      if(topography) then
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
      else
        do kz = 0, nz + 1
          do ix = - nbx, nx + nbx - 1
            uBar(ix, :, kz) = var%u(ix, :, kz) * (0.5 * (var%rho(ix, :, kz) &
                &+ var%rho(ix + 1, :, kz)) + rhoStrat(kz)) / Pstrat(kz)
          end do
        end do
      end if

      ! Compute \rho*v/P for reconstruction.
      if(topography) then
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
      else
        do kz = 0, nz + 1
          do jy = - nby, ny + nby - 1
            vBar(:, jy, kz) = var%v(:, jy, kz) * (0.5 * (var%rho(:, jy, kz) &
                &+ var%rho(:, jy + 1, kz)) + rhoStrat(kz)) / Pstrat(kz)
          end do
        end do
      end if

      ! Compute \rho*w/P for reconstruction.
      if(topography) then
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
      else
        do kz = 0, nz + 1
          wBar(:, :, kz) = var%w(:, :, kz) * (0.5 * (var%rho(:, :, kz) &
              &+ var%rho(:, :, kz + 1)) + rhoStratTilde(kz)) / PstratTilde(kz)
        end do
      end if

      ! Reconstruct \rho*u/P, \rho*v/P and \rho*w/P.
      call reconstruct_MUSCL(uBar, uTilde, nxx, nyy, nzz, limiterType1)
      call reconstruct_MUSCL(vBar, vTilde, nxx, nyy, nzz, limiterType1)
      call reconstruct_MUSCL(wBar, wTilde, nxx, nyy, nzz, limiterType1)

    case("P") ! reconstruct 1 = P/P
      PBar(:, :, :) = 1
      call reconstruct_MUSCL(PBar, PTilde, nxx, nyy, nzz, limiterType1)

    case("ice")

      do iVar = 1, nVarIce

        iceBar = 0.0
        if(topography) then
          ! Adjust reconstruction for 3D fields.
          do ix = - nbx, nx + nbx
            do jy = - nby, ny + nby
              do kz = 0, nz + 1
                if(pStratTFC(ix, jy, kz) == 0.0) then
                  print *, "Error in reconstruction: pStratTFC(" &
                      &// trim_integer(ix) // "," // trim_integer(jy) // "," &
                      &// trim_integer(kz) // ") = 0"
                  stop
                end if
                iceBar(ix, jy, kz) = var%ICE(ix, jy, kz, iVar) / pStratTFC(ix, &
                    &jy, kz)
              end do
            end do
          end do
        else
          do kz = 0, nz + 1
            if(Pstrat(kz) == 0.0) then
              print *, "Error in reconstruction: pStrat(" // trim_integer(kz) &
                  &// ") = 0"
              stop
            end if
            iceBar(:, :, kz) = (var%ICE(:, :, kz, iVar)) / Pstrat(kz)
          end do
        end if

        if(iVar .eq. inN) then

          call reconstruct_MUSCL(iceBar, nIceTilde, nxx, nyy, nzz, limiterType1)

        elseif(iVar .eq. inQ) then

          call reconstruct_MUSCL(iceBar, qIceTilde, nxx, nyy, nzz, limiterType1)

        elseif(iVar .eq. inQv) then

          call reconstruct_MUSCL(iceBar, qvTilde, nxx, nyy, nzz, limiterType1)

        end if
      end do

    case("tracer")
      ! calucate tracerBar

      tracerBar = 0.0
      if(topography) then
        do ix = - nbx, nx + nbx
          do jy = - nby, ny + nby
            do kz = 0, nz + 1
              if(pStratTFC(ix, jy, kz) == 0.0) then
                print *, "Error in reconstruction: pStratTFC(" &
                    &// trim_integer(ix) // "," // trim_integer(jy) // "," &
                    &// trim_integer(kz) // ") = 0"
                stop
              end if
              tracerBar(ix, jy, kz) = var%chi(ix, jy, kz) / pStratTFC(ix, jy, &
                  &kz)
            end do
          end do
        end do
      else
        do kz = 0, nz + 1
          if(Pstrat(kz) == 0.0) then
            print *, "Error in reconstruction: pStrat(" // trim_integer(kz) &
                &// ") = 0"
            stop
          end if

          tracerBar(:, :, kz) = var%chi(:, :, kz) / Pstrat(kz)
        end do
      end if

      call reconstruct_MUSCL(tracerBar, tracerTilde, nxx, nyy, nzz, &
          &limiterType1)

    case default
      stop "Error in reconstruction: unknown case variable."
    end select

  end subroutine reconstruction

  !--------------------------------------------------------------------

  subroutine massFlux(vara, var, flux, fluxmode, Pstrata, PStratTildea)
    !---------------------------------------------------------------------
    ! computes the mass flux at all cell edges using reconstructed values
    ! fluxmode = lin => linear flux, advecting velocities prescribed in
    !                   vara
    !            nln => nonlinear flux, advecting velocities from var
    !
    ! MUSCL assumes that the reconstructed densities are \rho/P
    !---------------------------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: vara, var
    character(len = *), intent(in) :: fluxmode

    type(flux_type), intent(inout) :: flux
    real, dimension(- 1:nz + 2), intent(in) :: PStrata, PStratTildea

    integer :: i, j, k, l
    real :: rhoL, rhoR, uL, uR ! L=Left i-1/2, R=Right i+1/2
    real :: rhoB, rhoF, vB, vF ! B=Backward j-1/2, F=Forward j+1/2
    real :: rhoD, rhoU, wD, wU ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurf, vSurf, wSurf ! velocities at cell surface

    real :: pL, pR, pB, pF, pD, pU
    real :: metD, metU

    ! TFC FJ
    real :: rhoStratEdgeR, rhoStratEdgeF, rhoStratEdgeU
    real :: pEdgeR, pEdgeF, pEdgeU

    real :: fRho, gRho, hRho

    ! avoid abs() for linerisation
    real :: delta
    real, parameter :: delta0 = 1.0e-6

    !   variables for the turbulence scheme
    real :: Pr_turb
    real :: coef_t, drho_dxi, dtht_dxi

    real :: delta_hs, delta_vs

    Pr_turb = 1.0 !0.5   !FS

    ! squared grid scales for the anisotropic turbulence scheme

    if(TurbScheme) then
      if(ny == 1 .and. nx == 1) then
        stop 'ERROR: turbulence scheme assumes either nx > 1 or ny > 1'
      else
        if(nx == 1) then
          delta_hs = dy ** 2 ! 2D problems in y and z
        else if(ny == 1) then
          delta_hs = dx ** 2 ! 2D problems in x and z
        else
          delta_hs = dx * dy ! 3D problems

          if(dx / dy > 10.) then
            print *, 'WARNING: dx/dy > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          elseif(dy / dx > 10.) then
            print *, 'WARNING: dy/dx > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          end if
        end if

        delta_vs = dz ** 2
      end if
    end if

    !-----------------------------------------
    !       Zonal rho fluxes in x: f
    !-----------------------------------------

    do k = 1, nz
      do j = 1, ny
        do i = 0, nx
          if(topography) then
            rhoStratEdgeR = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, &
                &j, k))
            pEdgeR = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i + 1, j, k))
            rhoR = rhoTilde(i + 1, j, k, 1, 0) + rhoStratEdgeR / pEdgeR
            rhoL = rhoTilde(i, j, k, 1, 1) + rhoStratEdgeR / pEdgeR
          else
            rhoR = rhoTilde(i + 1, j, k, 1, 0) + rhoStrat(k) / Pstrat(k)
            rhoL = rhoTilde(i, j, k, 1, 1) + rhoStrat(k) / Pstrat(k)
          end if

          if(topography) then
            pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, j, &
                &k) * pStratTFC(i + 1, j, k))
            if(fluxmode == "nln") then
              uSurf = pEdgeR * var%u(i, j, k)
            else if(fluxmode == "lin") then
              uSurf = pEdgeR * vara%u(i, j, k)
            else
              stop "Error in massFlux: Unknown fluxmode!"
            end if
          else
            if(fluxmode == "nln") then
              uSurf = var%u(i, j, k) * Pstrat(k)
            else if(fluxmode == "lin") then
              uSurf = vara%u(i, j, k) * Pstrata(k)
            else
              stop "Error in massFlux: Unknown fluxmode!"
            end if
          end if

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
          if(topography) then
            rhoStratEdgeF = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, j &
                &+ 1, k))
            pEdgeF = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j + 1, k))
            rhoF = rhoTilde(i, j + 1, k, 2, 0) + rhoStratEdgeF / pEdgeF
            rhoB = rhoTilde(i, j, k, 2, 1) + rhoStratEdgeF / pEdgeF
          else
            rhoF = rhoTilde(i, j + 1, k, 2, 0) + rhoStrat(k) / Pstrat(k)
            rhoB = rhoTilde(i, j, k, 2, 1) + rhoStrat(k) / Pstrat(k)
          end if

          if(topography) then
            pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j + 1, &
                &k) * pStratTFC(i, j + 1, k))
            if(fluxmode == "nln") then
              vSurf = pEdgeF * var%v(i, j, k)
            else if(fluxmode == "lin") then
              vSurf = pEdgeF * vara%v(i, j, k)
            else
              stop "Error in massFlux: Unknown fluxmode!"
            end if
          else
            if(fluxmode == "nln") then
              vSurf = var%v(i, j, k) * Pstrat(k)
            else if(fluxmode == "lin") then
              vSurf = vara%v(i, j, k) * Pstrata(k)
            else
              stop "Error in massFlux: Unknown fluxmode!"
            end if
          end if

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
          if(topography) then
            rhoStratEdgeU = (jac(i, j, k + 1) * rhoStratTFC(i, j, k) + jac(i, &
                &j, k) * rhoStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, &
                &k + 1))
            pEdgeU = (jac(i, j, k + 1) * pStratTFC(i, j, k) + jac(i, j, k) &
                &* pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            rhoU = rhoTilde(i, j, k + 1, 3, 0) + rhoStratEdgeU / pEdgeU
            rhoD = rhoTilde(i, j, k, 3, 1) + rhoStratEdgeU / pEdgeU
          else
            rhoU = rhoTilde(i, j, k + 1, 3, 0) + rhoStratTilde(k) &
                &/ PstratTilde(k)
            rhoD = rhoTilde(i, j, k, 3, 1) + rhoStratTilde(k) / PstratTilde(k)
          end if

          if(topography) then
            pEdgeU = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
                &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            if(fluxmode == "nln") then
              wSurf = pEdgeU * var%w(i, j, k)
            else if(fluxmode == "lin") then
              wSurf = pEdgeU * vara%w(i, j, k)
            else
              stop "Error in massFlux: Unknown fluxmode!"
            end if
          else
            if(fluxmode == "nln") then
              wSurf = var%w(i, j, k) * PstratTilde(k)
            else if(fluxmode == "lin") then
              wSurf = vara%w(i, j, k) * PstratTildea(k)
            else
              stop "Error in massFlux: Unknown fluxmode!"
            end if
          end if

          hRho = flux_muscl(wSurf, rhoD, rhoU)

          flux%rho(i, j, k, 3) = hRho
        end do
      end do
    end do

    ! --------------------------------------------
    !        Density-fluctuation fluxes
    ! --------------------------------------------

    if(timeScheme == "semiimplicit" .or. auxil_equ) then

      ! Zonal rhop fluxes in x: f
      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            rhoR = rhopTilde(i + 1, j, k, 1, 0)
            rhoL = rhopTilde(i, j, k, 1, 1)

            if(topography) then
              pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, &
                  &j, k) * pStratTFC(i + 1, j, k))
              if(fluxmode == "nln") then
                uSurf = pEdgeR * var%u(i, j, k)
              else if(fluxmode == "lin") then
                uSurf = pEdgeR * vara%u(i, j, k)
              else
                stop "Error in massFlux: Unknown fluxmode!"
              end if
            else
              if(fluxmode == "nln") then
                uSurf = var%u(i, j, k) * Pstrat(k)
              else if(fluxmode == "lin") then
                uSurf = vara%u(i, j, k) * Pstrata(k)
              else
                stop "Error in massFlux: Unknown fluxmode!"
              end if
            end if

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

            if(topography) then
              pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                  &+ 1, k) * pStratTFC(i, j + 1, k))
              if(fluxmode == "nln") then
                vSurf = pEdgeF * var%v(i, j, k)
              else if(fluxmode == "lin") then
                vSurf = pEdgeF * vara%v(i, j, k)
              else
                stop "Error in massFlux: Unknown fluxmode!"
              end if
            else
              if(fluxmode == "nln") then
                vSurf = var%v(i, j, k) * Pstrat(k)
              else if(fluxmode == "lin") then
                vSurf = vara%v(i, j, k) * Pstrata(k)
              else
                stop "Error in massFlux: Unknown fluxmode!"
              end if
            end if

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

            if(topography) then
              pEdgeU = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
                  &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
              if(fluxmode == "nln") then
                wSurf = pEdgeU * var%w(i, j, k)
              else if(fluxmode == "lin") then
                wSurf = pEdgeU * vara%w(i, j, k)
              else
                stop "Error in massFlux: Unknown fluxmode!"
              end if
            else
              if(fluxmode == "nln") then
                wSurf = var%w(i, j, k) * PstratTilde(k)
              else if(fluxmode == "lin") then
                wSurf = vara%w(i, j, k) * PstratTildea(k)
              else
                stop "Error in massFlux: Unknown fluxmode!"
              end if
            end if

            hRho = flux_muscl(wSurf, rhoD, rhoU)

            flux%rhop(i, j, k, 3) = hRho
          end do
        end do
      end do
    end if

    ! --------------------------------------------
    !  Mass-weighted potential-temperature fluxes
    !---------------------------------------------

    if(model == "compressible") then

      ! Zonal fluxes in x: f
      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            rhoR = 1.0
            rhoL = 1.0

            if(topography) then
              pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, &
                  &j, k) * pStratTFC(i + 1, j, k))
              if(fluxmode == "nln") then
                uSurf = pEdgeR * var%u(i, j, k)
              else if(fluxmode == "lin") then
                uSurf = pEdgeR * vara%u(i, j, k)
              else
                stop "Error in massFlux: Unknown fluxmode!"
              end if
            else
              if(fluxmode == "nln") then
                uSurf = var%u(i, j, k) * Pstrat(k)
              else if(fluxmode == "lin") then
                uSurf = vara%u(i, j, k) * Pstrata(k)
              else
                stop "Error in massFlux: Unknown fluxmode!"
              end if
            end if

            fRho = flux_muscl(uSurf, rhoL, rhoR)

            flux%P(i, j, k, 1) = fRho
          end do
        end do
      end do

      ! Meridional rhop fluxes in y: g
      do k = 1, nz
        do j = 0, ny
          do i = 1, nx
            rhoF = 1.0
            rhoB = 1.0

            if(topography) then
              pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                  &+ 1, k) * pStratTFC(i, j + 1, k))
              if(fluxmode == "nln") then
                vSurf = pEdgeF * var%v(i, j, k)
              else if(fluxmode == "lin") then
                vSurf = pEdgeF * vara%v(i, j, k)
              else
                stop "Error in massFlux: Unknown fluxmode!"
              end if
            else
              if(fluxmode == "nln") then
                vSurf = var%v(i, j, k) * Pstrat(k)
              else if(fluxmode == "lin") then
                vSurf = vara%v(i, j, k) * Pstrata(k)
              else
                stop "Error in massFlux: Unknown fluxmode!"
              end if
            end if

            gRho = flux_muscl(vSurf, rhoB, rhoF)

            flux%P(i, j, k, 2) = gRho
          end do
        end do
      end do

      ! Vertical rhop fluxes in z: h
      do k = 0, nz
        do j = 1, ny
          do i = 1, nx
            rhoU = 1.0
            rhoD = 1.0

            if(topography) then
              pEdgeU = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
                  &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
              if(fluxmode == "nln") then
                wSurf = pEdgeU * var%w(i, j, k)
              else if(fluxmode == "lin") then
                wSurf = pEdgeU * vara%w(i, j, k)
              else
                stop "Error in massFlux: Unknown fluxmode!"
              end if
            else
              if(fluxmode == "nln") then
                wSurf = var%w(i, j, k) * PstratTilde(k) !UA
              else if(fluxmode == "lin") then
                wSurf = vara%w(i, j, k) * PstratTildea(k) !UA
              else
                stop "Error in massFlux: Unknown fluxmode!"
              end if
            end if

            hRho = flux_muscl(wSurf, rhoD, rhoU)

            flux%P(i, j, k, 3) = hRho
          end do
        end do
      end do
    end if

    !--------------------------------------------------------
    !  Contributions from molecular and turbulent diffusion
    !  to the potential-temperature fluxes
    !--------------------------------------------------------

    if(mu_conduct == 0.0 .and. .not. TurbScheme) return

    ! Flux in x direction
    if(topography) then
      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            select case(model)
            case("pseudo_incompressible", "compressible")
              coef_t = mu_conduct * 0.5 * (rhoStratTFC(i, j, 1) &
                  &/ rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, j, 1) &
                  &/ rhoStratTFC(i + 1, j, k))
            case("Boussinesq")
              coef_t = mu_conduct
            case default
              stop "Error in massFlux: Unkown case model."
            end select

            if(TurbScheme) then
              coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i + 1, j, &
                  &k)) * delta_hs / Pr_turb
            end if

            rhoL = var%rho(i, j, k) + rhoStratTFC(i, j, k)
            rhoR = var%rho(i + 1, j, k) + rhoStratTFC(i + 1, j, k)
            rhoD = 0.5 * (var%rho(i, j, k - 1) + rhoStratTFC(i, j, k - 1) &
                &+ var%rho(i + 1, j, k - 1) + rhoStratTFC(i + 1, j, k - 1))
            rhoU = 0.5 * (var%rho(i, j, k + 1) + rhoStratTFC(i, j, k + 1) &
                &+ var%rho(i + 1, j, k + 1) + rhoStratTFC(i + 1, j, k + 1))

            pL = pStratTFC(i, j, k)
            pR = pStratTFC(i + 1, j, k)
            pD = 0.5 * (pStratTFC(i, j, k - 1) + pStratTFC(i + 1, j, k - 1))
            pU = 0.5 * (pStratTFC(i, j, k + 1) + pStratTFC(i + 1, j, k + 1))

            dtht_dxi = 0.5 * (jac(i, j, k) + jac(i + 1, j, k)) * (pR / rhoR &
                &- pL / rhoL) / dx + 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) &
                &+ jac(i + 1, j, k) * met(i + 1, j, k, 1, 3)) * (pU / rhoU &
                &- pD / rhoD) / (2.0 * dz)

            flux%theta(i, j, k, 1) = - coef_t * dtht_dxi
          end do
        end do
      end do
    else
      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            select case(model)
            case("pseudo_incompressible")
              coef_t = mu_conduct * rhoStrat(1) / rhoStrat(k)
            case("Boussinesq")
              coef_t = mu_conduct
            case default
              stop "Error in massFlux: Unkown case model."
            end select

            if(TurbScheme) then
              coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i + 1, j, &
                  &k)) * delta_hs / Pr_turb
            end if

            rhoL = var%rho(i, j, k) + rhoStrat(k)
            rhoR = var%rho(i + 1, j, k) + rhoStrat(k)

            dtht_dxi = (Pstrat(k) / rhoR - Pstrat(k) / rhoL) / dx

            flux%theta(i, j, k, 1) = - coef_t * dtht_dxi
          end do
        end do
      end do
    end if

    ! Flux in y direction
    if(topography) then
      do k = 1, nz
        do j = 0, ny
          do i = 1, nx
            select case(model)
            case("pseudo_incompressible", "compressible")
              coef_t = mu_conduct * 0.5 * (rhoStratTFC(i, j, 1) &
                  &/ rhoStratTFC(i, j, k) + rhoStratTFC(i, j + 1, 1) &
                  &/ rhoStratTFC(i, j + 1, k))
            case("Boussinesq")
              coef_t = mu_conduct
            case default
              stop "Error in massFlux: Unkown case model."
            end select

            if(TurbScheme) then
              coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i, j + 1, &
                  &k)) * delta_hs / Pr_turb
            end if

            rhoB = var%rho(i, j, k) + rhoStratTFC(i, j, k)
            rhoF = var%rho(i, j + 1, k) + rhoStratTFC(i, j + 1, k)
            rhoD = 0.5 * (var%rho(i, j, k - 1) + rhoStratTFC(i, j, k - 1) &
                &+ var%rho(i, j + 1, k - 1) + rhoStratTFC(i, j + 1, k - 1))
            rhoU = 0.5 * (var%rho(i, j, k + 1) + rhoStratTFC(i, j, k + 1) &
                &+ var%rho(i, j + 1, k + 1) + rhoStratTFC(i, j + 1, k + 1))

            pB = pStratTFC(i, j, k)
            pF = pStratTFC(i, j + 1, k)
            pD = 0.5 * (pStratTFC(i, j, k - 1) + pStratTFC(i, j + 1, k - 1))
            pU = 0.5 * (pStratTFC(i, j, k + 1) + pStratTFC(i, j + 1, k + 1))

            dtht_dxi = 0.5 * (jac(i, j, k) + jac(i, j + 1, k)) * (pF / rhoF &
                &- pB / rhoB) / dy + 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) &
                &+ jac(i, j + 1, k) * met(i, j + 1, k, 2, 3)) * (pU / rhoU &
                &- pD / rhoD) / (2.0 * dz)

            flux%theta(i, j, k, 2) = - coef_t * dtht_dxi
          end do
        end do
      end do
    else
      do k = 1, nz
        do j = 0, ny
          do i = 1, nx
            select case(model)
            case("pseudo_incompressible")
              coef_t = mu_conduct * rhoStrat(1) / rhoStrat(k)
            case("Boussinesq")
              coef_t = mu_conduct
            case default
              stop "Error in massFlux: Unkown case model."
            end select

            if(TurbScheme) then
              coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i, j + 1, &
                  &k)) * delta_hs / Pr_turb
            end if

            rhoF = var%rho(i, j + 1, k) + rhoStrat(k)
            rhoB = var%rho(i, j, k) + rhoStrat(k)

            dtht_dxi = (Pstrat(k) / rhoF - Pstrat(k) / rhoB) / dy

            flux%theta(i, j, k, 2) = - coef_t * dtht_dxi
          end do
        end do
      end do
    end if

    ! Flux in z direction
    if(topography) then
      do k = 0, nz
        do j = 1, ny
          do i = 1, nx
            select case(model)
            case("pseudo_incompressible", "compressible")
              coef_t = mu_conduct * rhoStratTFC(i, j, 1) * (jac(i, j, k + 1) &
                  &/ rhoStratTFC(i, j, k) + jac(i, j, k) / rhoStratTFC(i, j, k &
                  &+ 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            case("Boussinesq")
              coef_t = mu_conduct
            case default
              stop "Error in massFlux: Unkown case model."
            end select

            if(TurbScheme) then
              coef_t = coef_t + jac(i, j, k) * jac(i, j, k + 1) * (jac(i, j, &
                  &k) * var%DSC(i, j, k) + jac(i, j, k + 1) * var%DSC(i, j, k &
                  &+ 1)) / (jac(i, j, k) + jac(i, j, k + 1)) * delta_vs &
                  &/ Pr_turb
            end if

            rhoL = (jac(i - 1, j, k + 1) * (var%rho(i - 1, j, k) &
                &+ rhoStratTFC(i - 1, j, k)) + jac(i - 1, j, k) * (var%rho(i &
                &- 1, j, k + 1) + rhoStratTFC(i - 1, j, k + 1))) / (jac(i - 1, &
                &j, k) + jac(i - 1, j, k + 1))
            rhoR = (jac(i + 1, j, k + 1) * (var%rho(i + 1, j, k) &
                &+ rhoStratTFC(i + 1, j, k)) + jac(i + 1, j, k) * (var%rho(i &
                &+ 1, j, k + 1) + rhoStratTFC(i + 1, j, k + 1))) / (jac(i + 1, &
                &j, k) + jac(i + 1, j, k + 1))
            rhoB = (jac(i, j - 1, k + 1) * (var%rho(i, j - 1, k) &
                &+ rhoStratTFC(i, j - 1, k)) + jac(i, j - 1, k) * (var%rho(i, &
                &j - 1, k + 1) + rhoStratTFC(i, j - 1, k + 1))) / (jac(i, j &
                &- 1, k) + jac(i, j - 1, k + 1))
            rhoF = (jac(i, j + 1, k + 1) * (var%rho(i, j + 1, k) &
                &+ rhoStratTFC(i, j + 1, k)) + jac(i, j + 1, k) * (var%rho(i, &
                &j + 1, k + 1) + rhoStratTFC(i, j + 1, k + 1))) / (jac(i, j &
                &+ 1, k) + jac(i, j + 1, k + 1))
            rhoD = var%rho(i, j, k) + rhoStratTFC(i, j, k)
            rhoU = var%rho(i, j, k + 1) + rhoStratTFC(i, j, k + 1)

            pL = (jac(i - 1, j, k + 1) * pStratTFC(i - 1, j, k) + jac(i - 1, &
                &j, k) * pStratTFC(i - 1, j, k + 1)) / (jac(i - 1, j, k) &
                &+ jac(i - 1, j, k + 1))
            pR = (jac(i + 1, j, k + 1) * pStratTFC(i + 1, j, k) + jac(i + 1, &
                &j, k) * pStratTFC(i + 1, j, k + 1)) / (jac(i + 1, j, k) &
                &+ jac(i + 1, j, k + 1))
            pB = (jac(i, j - 1, k + 1) * pStratTFC(i, j - 1, k) + jac(i, j &
                &- 1, k) * pStratTFC(i, j - 1, k + 1)) / (jac(i, j - 1, k) &
                &+ jac(i, j - 1, k + 1))
            pF = (jac(i, j + 1, k + 1) * pStratTFC(i, j + 1, k) + jac(i, j &
                &+ 1, k) * pStratTFC(i, j + 1, k + 1)) / (jac(i, j + 1, k) &
                &+ jac(i, j + 1, k + 1))
            pD = pStratTFC(i, j, k)
            pU = pStratTFC(i, j, k + 1)

            dtht_dxi = jac(i, j, k) * jac(i, j, k + 1) * ((met(i, j, k, 1, 3) &
                &+ met(i, j, k + 1, 1, 3)) * (pR / rhoR - pL / rhoL) / (2.0 &
                &* dx) + (met(i, j, k, 2, 3) + met(i, j, k + 1, 2, 3)) * (pF &
                &/ rhoF - pB / rhoB) / (2.0 * dy) + (met(i, j, k, 3, 3) &
                &+ met(i, j, k + 1, 3, 3)) * (pU / rhoU - pD / rhoD) / dz) &
                &/ (jac(i, j, k) + jac(i, j, k + 1))

            flux%theta(i, j, k, 3) = - coef_t * dtht_dxi
          end do
        end do
      end do
    else
      do k = 0, nz
        do j = 1, ny
          do i = 1, nx
            select case(model)
            case("pseudo_incompressible")
              coef_t = mu_conduct * rhoStrat(1) / rhoStrat(k)
            case("Boussinesq")
              coef_t = mu_conduct
            case default
              stop "Error in massFlux: Unkown case model."
            end select

            if(TurbScheme) then
              coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i, j, k &
                  &+ 1)) * delta_vs / Pr_turb
            end if

            rhoU = var%rho(i, j, k + 1) + rhoStrat(k + 1)
            rhoD = var%rho(i, j, k) + rhoStrat(k)

            dtht_dxi = (Pstrat(k + 1) / rhoU - Pstrat(k) / rhoD) / dz

            flux%theta(i, j, k, 3) = - coef_t * dtht_dxi
          end do
        end do
      end do
    end if

    if(verbose) print *, "rhoFlux: rho fluxes fRho, gRho and fRho calculated"

  end subroutine massFlux

  !-----------------------------------------------------------------------

  subroutine tracerFlux(vara, var, flux, fluxmode, Pstrata, PStratTildea)
    ! calculate tracer fluxes
    ! zonal:      u*rho*chi = flux%chi(:, :, :, 1)
    ! meridional: v*rho*chi = flux%chi(:, :, :, 2)
    ! vertical:   w*rho*chi = flux%chi(:, :, :, 3)
    implicit none

    type(var_type), intent(in) :: vara, var
    character(len = *), intent(in) :: fluxmode

    type(flux_type), intent(inout) :: flux

    real, dimension(- 1:nz + 2), intent(in) :: PStrata, PStratTildea

    integer :: i, j, k, l
    real :: tracerL, tracerR, uL, uR ! L=Left i-1/2, R=Right i+1/2
    real :: tracerB, tracerF, vB, vF ! B=Backward j-1/2, F=Forward j+1/2
    real :: tracerD, tracerU, wD, wU ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurf, vSurf, wSurf ! velocities at cell surface
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: rho

    real :: rhoStratEdgeR, rhoStratEdgeF, rhoStratEdgeU
    real :: pEdgeR, pEdgeF, pEdgeU

    real :: fTracer, gTracer, hTracer

    if(TurbScheme) then
      stop "fluxes.f90: TurbScheme not implemented for tracer."
    end if

    !-----------------------------------------
    !       Zonal tracer fluxes in x: f
    !-----------------------------------------
    do k = 1, nz
      do j = 1, ny
        do i = 0, nx
          tracerR = tracerTilde(i + 1, j, k, 1, 0)
          tracerL = tracerTilde(i, j, k, 1, 1)

          if(topography) then
            pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, j, &
                &k) * pStratTFC(i + 1, j, k))
            if(fluxmode == "nln") then
              uSurf = pEdgeR * var%u(i, j, k)
            else if(fluxmode == "lin") then
              uSurf = pEdgeR * vara%u(i, j, k)
            else
              stop "fluxes.f90: wrong fluxmode"
            end if
          else
            if(fluxmode == "nln") then
              uSurf = var%u(i, j, k) * Pstrat(k)
            else if(fluxmode == "lin") then
              uSurf = vara%u(i, j, k) * Pstrata(k)
            else
              stop "fluxes.f90: wrong fluxmode"
            end if
          end if

          fTracer = flux_muscl(uSurf, tracerL, tracerR)

          flux%chi(i, j, k, 1) = fTracer
        end do
      end do
    end do

    !-----------------------------------
    ! Meridional tracer fluxes in y: g
    !-----------------------------------
    do k = 1, nz
      do j = 0, ny
        do i = 1, nx
          tracerF = tracerTilde(i, j + 1, k, 2, 0)
          tracerB = tracerTilde(i, j, k, 2, 1)

          if(topography) then
            pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j + 1, &
                &k) * pStratTFC(i, j + 1, k))
            if(fluxmode == "nln") then
              vSurf = pEdgeF * var%v(i, j, k)
            else if(fluxmode == "lin") then
              vSurf = pEdgeF * vara%v(i, j, k)
            else
              stop "fluxes.f90: wrong fluxmode"
            end if
          else
            if(fluxmode == "nln") then
              vSurf = var%v(i, j, k) * Pstrat(k)
            else if(fluxmode == "lin") then
              vSurf = vara%v(i, j, k) * Pstrata(k)
            else
              stop "fluxes.f90: wrong fluxmode"
            end if
          end if

          gTracer = flux_muscl(vSurf, tracerB, tracerF)

          flux%chi(i, j, k, 2) = gTracer
        end do
      end do
    end do

    !--------------------------------
    ! Vertical tracer fluxes in z: h
    !-------------------------------
    do k = 0, nz
      do j = 1, ny
        do i = 1, nx
          tracerU = tracerTilde(i, j, k + 1, 3, 0)
          tracerD = tracerTilde(i, j, k, 3, 1)

          if(topography) then
            pEdgeU = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
                &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            if(fluxmode == "nln") then
              wSurf = pEdgeU * var%w(i, j, k)
            else if(fluxmode == "lin") then
              wSurf = pEdgeU * vara%w(i, j, k)
            else
              stop "fluxes.f90: wrong fluxmode"
            end if
          else
            if(fluxmode == "nln") then
              wSurf = var%w(i, j, k) * PstratTilde(k)
            else if(fluxmode == "lin") then
              wSurf = vara%w(i, j, k) * PstratTildea(k)
            else
              stop "fluxes.f90: wrong fluxmode"
            end if
          end if

          hTracer = flux_muscl(wSurf, tracerD, tracerU)

          flux%chi(i, j, k, 3) = hTracer
        end do
      end do
    end do

  end subroutine tracerFlux

  !-----------------------------------------------------------------------!

  subroutine volumeForce(var, time, force)

    !------------------------------------------------------------
    ! supplememts volume forces on a grid cell by gravitation and
    ! Coriolis force
    !------------------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var

    ! volume forces
    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, 3), intent(inout) :: force

    ! local variables
    integer :: i, j, k, l

    ! gravitational forces
    real :: dRho, theta
    real, dimension(3) :: gForce

    ! Coriolis force
    real :: n1, n2, n3
    real :: u1, u2, u3
    real :: f1, f2, f3
    real :: rho

    ! zonal wind relaxation
    real :: ft_relax
    real :: time
    real, dimension(0:nx + 1) :: fx_relax
    real, dimension(0:ny + 1) :: fy_relax
    integer :: i0, j0

    if(timeScheme /= 'semiimplicit') then
      !--------------------------------------------
      !             Gravitational force
      !--------------------------------------------

      do k = 0, nz + 1
        do j = 0, ny + 1
          do i = 0, nx + 1

            select case(model)

            case("pseudo_incompressible", "Boussinesq", "compressible")

              if(auxil_equ) then
                dRho = var%rhop(i, j, k)
              else
                dRho = var%rho(i, j, k)
              end if

              gForce = - FrInv2 * dRho * vertical

            case default
              stop "volumeForce: unknown case model."
            end select

            ! Adjust gravitational force.
            if(topography) then
              gForce = gForce / jac(i, j, k)
            end if

            force(i, j, k, 1:3) = force(i, j, k, 1:3) + gForce

          end do
        end do
      end do

      !--------------------------------------------
      !             Coriolis force
      !--------------------------------------------
      ! Coriolis force defined in scalar volume cells
      ! -> needs interpolation

      do k = 0, nz + 1
        do j = 0, ny + 1
          do i = 0, nx + 1
            u1 = 0.5 * (var%u(i, j, k) + var%u(i - 1, j, k))
            if(testCase == "SkamarockKlemp94") then
              u1 = u1 - backgroundFlow_dim(1) / uRef
            end if

            u2 = 0.5 * (var%v(i, j, k) + var%v(i, j - 1, k))
            u3 = 0.5 * (var%w(i, j, k) + var%w(i, j, k - 1))

            select case(model)

            case("Boussinesq")
              rho = rho00

            case("pseudo_incompressible", "compressible")
              if(topography) then
                ! Adjust for 3D fields.
                rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)
              else
                rho = var%rho(i, j, k) + rhoStrat(k)
              end if

            case default
              stop "volumeForce: unknown case model."
            end select

            n1 = vertical(1)
            n2 = vertical(2)
            n3 = vertical(3)

            f1 = n2 * u3 - n3 * u2
            f2 = n3 * u1 - n1 * u3
            f3 = n1 * u2 - n2 * u1

            ! Coriolis force normally written in LHS with "+"
            ! gets now a "-" since force is assumed on the RHS

            force(i, j, k, 1:3) = force(i, j, k, 1:3) - (rho * RoInv(j) * &
                &(/f1, f2, f3/))

            ! Add vertical Coriolis force component in TFC.
            if(topography) then
              force(i, j, k, 1:3) = force(i, j, k, 1:3) + rho * RoInv(j) &
                  &* vertical * (met(i, j, k, 3, 1) * u2 - met(i, j, k, 3, 2) &
                  &* u1)
            end if

          end do
        end do
      end do

    end if ! not semiimplicit

    !--------------------------------------------
    !             topography growth
    !--------------------------------------------

    if(topography .and. topographyTime > 0.0) then
      if(any(topography_surface /= final_topography_surface)) then
        do k = 0, nz + 1
          do j = 0, ny + 1
            do i = 0, nx + 1
              ! Zonal part
              force(i, j, k, 3) = force(i, j, k, 3) + (var%rho(i, j, k) &
                  &+ rhoStratTFC(i, j, k)) * 0.5 * (var%u(i, j, k) + var%u(i &
                  &- 1, j, k)) * dMet13Dt(i, j, k)
              ! Meridional part
              force(i, j, k, 3) = force(i, j, k, 3) + (var%rho(i, j, k) &
                  &+ rhoStratTFC(i, j, k)) * 0.5 * (var%v(i, j, k) + var%v(i, &
                  &j - 1, k)) * dMet23Dt(i, j, k)
              ! Vertical part
              force(i, j, k, 3) = force(i, j, k, 3) + (var%rho(i, j, k) &
                  &+ rhoStratTFC(i, j, k)) * 0.5 * (vertWindTFC(i, j, k, var) &
                  &+ vertWindTFC(i, j, k - 1, var)) * dJacInvDt(i, j, k)
            end do
          end do
        end do
      end if
    end if

    !--------------------------------------------
    !               wind relaxation
    !--------------------------------------------

    ! TFC FJ
    if(.not. wind_relaxation) return

    i0 = is + nbx - 1
    j0 = js + nby - 1

    if((testCase == "mountainwave") .or. (raytracer .and. case_wkb == 3)) then

      if(xextent_relax > 0.0) then
        do i = 0, nx + 1
          if(x(i0 + i) < lx(0) + 0.5 * xextent_relax) then
            fx_relax(i) = cos((x(i0 + i) - lx(0)) * pi / xextent_relax)
          else if(x(i0 + i) < lx(1) - 0.5 * xextent_relax) then
            fx_relax(i) = 0.0
          else
            fx_relax(i) = cos((lx(1) - x(i0 + i)) * pi / xextent_relax)
          end if
        end do
      else
        do i = 0, nx + 1
          fx_relax(i) = 0.0
        end do
      end if

      if(yextent_relax > 0.0) then
        do j = 0, ny + 1
          if(y(j0 + j) < ly(0) + 0.5 * yextent_relax) then
            fy_relax(j) = cos((y(j0 + j) - ly(0)) * pi / yextent_relax)
          else if(y(j0 + j) < ly(1) - 0.5 * yextent_relax) then
            fy_relax(j) = 0.0
          else
            fy_relax(j) = cos((ly(1) - y(j0 + j)) * pi / yextent_relax)
          end if
        end do
      else
        do j = 0, ny + 1
          fy_relax(j) = 0.0
        end do
      end if

      do k = 0, nz + 1
        do j = 0, ny + 1
          do i = 0, nx + 1

            select case(model)

            case("Boussinesq")
              rho = rho00

            case("pseudo_incompressible", "compressible")

              if(topography) then
                rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)
              else
                rho = var%rho(i, j, k) + rhoStrat(k)
              end if

            case default
              stop "volumeForce: unknown case model."
            end select

            force(i, j, k, 1) = force(i, j, k, 1) - rho * (0.5 * (var%u(i, j, &
                &k) + var%u(i - 1, j, k)) - u_relax) * fx_relax(i) &
                &* fy_relax(j) / t_relax

            force(i, j, k, 2) = force(i, j, k, 2) - rho * (0.5 * (var%v(i, j, &
                &k) + var%v(i, j - 1, k)) - v_relax) * fx_relax(i) &
                &* fy_relax(j) / t_relax

            force(i, j, k, 3) = force(i, j, k, 3) - rho * (0.5 * (var%w(i, j, &
                &k) + var%w(i, j, k - 1)) - w_relax) * fx_relax(i) &
                &* fy_relax(j) / t_relax

          end do
        end do
      end do

    end if ! mountainwave

  end subroutine volumeForce

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

  subroutine momentumFlux(vara, var, flux, fluxmode, PStrata, PStratTildea)

    !----------------------------------------------------------------------
    ! computes the momentum fluxes at the cell edges using reconstr. values
    ! fluxmode = lin => linear flux, advecting velocities prescribed in
    !                   vara
    !            nln => nonlinear flux, advecting velocities from var
    !
    ! MUSCL assumes that the reconstructed momenta are \rho/P * \vec v
    !----------------------------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: vara, var
    character(len = *), intent(in) :: fluxmode

    type(flux_type), intent(inout) :: flux
    ! flux(i,j,k,dir,iFlux)
    ! dir = 1..3 > f,g,h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, fRhoV, fRhoW

    real, dimension(- 1:nz + 2), intent(in) :: PStrata, PStratTildea

    ! local variables
    integer :: i, j, k, l

    ! density at cell edge
    real :: rhoEdge

    ! uTilde at cell edges
    real :: uL, uR, vL, vR, wL, wR ! L=Left at i+1, R=Right at i
    real :: uB, uF, vB, vF, wB, wF ! B=Backward at j+1, F=Forward at j
    real :: uD, uU, vD, vU, wD, wU ! D=Downward at k+1, U=Upward at k

    ! advecting velocities
    real :: uL0, uR0, vL0, vR0, wL0, wR0
    real :: uB0, uF0, vB0, vF0, wB0, wF0
    real :: uD0, uU0, vD0, vU0, wD0, wU0

    ! cell averaged values at cell centres
    real :: uBarL, uBarR, vBarL, vBarR, wBarL, wBarR ! at i and i+1
    real :: uBarB, uBarF, vBarB, vBarF, wBarB, wBarF ! at j and j+1
    real :: uBarD, uBarU, vBarD, vBarU, wBarD, wBarU ! at k and k+1

    ! TFC FJ
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
    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1) :: divU ! div(u) field of visc.
    real :: div
    real :: du_dx, du_dy, du_dz ! partial derivatives
    real :: dv_dx, dv_dy, dv_dz
    real :: dw_dx, dw_dy, dw_dz
    real :: fRhoU_visc, gRhoU_visc, hRhoU_visc ! viscous mom. fluxes
    real :: fRhoV_visc, gRhoV_visc, hRhoV_visc
    real :: fRhoW_visc, gRhoW_visc, hRhoW_visc

    ! debugging
    real :: rhoEdge2

    ! avoid abs() for linerisation
    real :: delta

    !for basic-state density
    real :: rhos

    !achatzb for putting together molecular and turbulent vioscosity
    real :: coef_v
    !achatze

    real :: delta_hs, delta_vs

    ! squared grid scales for the anisotropic turbulence scheme

    if(TurbScheme) then
      if(ny == 1 .and. nx == 1) then
        stop 'ERROR: turbulence scheme assumes either nx > 1 or ny > 1'
      else
        if(nx == 1) then
          delta_hs = dy ** 2 ! 2D problems in y and z
        else if(ny == 1) then
          delta_hs = dx ** 2 ! 2D problems in x and z
        else
          delta_hs = dx * dy ! 3D problems

          if(dx / dy > 10.) then
            print *, 'WARNING: dx/dy > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          elseif(dy / dx > 10.) then
            print *, 'WARNING: dy/dx > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          end if
        end if

        delta_vs = dz ** 2
      end if
    end if

    if(verbose) print *, "fluxes.f90/momentumFlux: Entering subroutine..."

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

          if(topography) then
            pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, j, &
                &k) * pStratTFC(i + 1, j, k))
            pREdgeR = 0.5 * (jac(i + 1, j, k) * pStratTFC(i + 1, j, k) + jac(i &
                &+ 2, j, k) * pStratTFC(i + 2, j, k))
            if(fluxmode == "nln") then
              uSurf = 0.5 * (pEdgeR * var%u(i, j, k) + pREdgeR * var%u(i + 1, &
                  &j, k))
            else if(fluxmode == "lin") then
              uSurf = 0.5 * (pEdgeR * vara%u(i, j, k) + pREdgeR * vara%u(i &
                  &+ 1, j, k))
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          else
            if(fluxmode == "nln") then
              usurf = 0.5 * (var%u(i, j, k) + var%u(i + 1, j, k)) * Pstrat(k)
            else if(fluxmode == "lin") then
              usurf = 0.5 * (vara%u(i, j, k) + vara%u(i + 1, j, k)) * Pstrata(k)
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          end if

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

          if(topography) then
            pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j + 1, &
                &k) * pStratTFC(i, j + 1, k))
            pREdgeF = 0.5 * (jac(i + 1, j, k) * pStratTFC(i + 1, j, k) + jac(i &
                &+ 1, j + 1, k) * pStratTFC(i + 1, j + 1, k))
            if(fluxmode == "nln") then
              vSurf = 0.5 * (pEdgeF * var%v(i, j, k) + pREdgeF * var%v(i + 1, &
                  &j, k))
            else if(fluxmode == "lin") then
              vSurf = 0.5 * (pEdgeF * vara%v(i, j, k) + pREdgeF * vara%v(i &
                  &+ 1, j, k))
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          else
            if(fluxmode == "nln") then
              vsurf = 0.5 * (var%v(i, j, k) + var%v(i + 1, j, k)) * Pstrat(k)
            else if(fluxmode == "lin") then
              vsurf = 0.5 * (vara%v(i, j, k) + vara%v(i + 1, j, k)) * Pstrata(k)
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          end if

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

          if(topography) then
            pEdgeU = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
                &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            pREdgeU = jac(i + 1, j, k) * jac(i + 1, j, k + 1) * (pStratTFC(i &
                &+ 1, j, k) + pStratTFC(i + 1, j, k + 1)) / (jac(i + 1, j, k) &
                &+ jac(i + 1, j, k + 1))
            if(fluxmode == "nln") then
              wSurf = 0.5 * (pEdgeU * var%w(i, j, k) + pREdgeU * var%w(i + 1, &
                  &j, k))
            else if(fluxmode == "lin") then
              wSurf = 0.5 * (pEdgeU * vara%w(i, j, k) + pREdgeU * vara%w(i &
                  &+ 1, j, k))
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          else
            if(fluxmode == "nln") then
              wsurf = 0.5 * (var%w(i, j, k) + var%w(i + 1, j, k)) &
                  &* PstratTilde(k)
            else if(fluxmode == "lin") then
              wsurf = 0.5 * (vara%w(i, j, k) + vara%w(i + 1, j, k)) &
                  &* PstratTildea(k)
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          end if

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

          if(topography) then
            pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, j, &
                &k) * pStratTFC(i + 1, j, k))
            pFEdgeR = 0.5 * (jac(i, j + 1, k) * pStratTFC(i, j + 1, k) + jac(i &
                &+ 1, j + 1, k) * pStratTFC(i + 1, j + 1, k))
            if(fluxmode == "nln") then
              uSurf = 0.5 * (pEdgeR * var%u(i, j, k) + pFEdgeR * var%u(i, j &
                  &+ 1, k))
            else if(fluxmode == "lin") then
              uSurf = 0.5 * (pEdgeR * vara%u(i, j, k) + pFEdgeR * vara%u(i, j &
                  &+ 1, k))
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          else
            if(fluxmode == "nln") then
              usurf = 0.5 * (var%u(i, j, k) + var%u(i, j + 1, k)) * Pstrat(k)
            else if(fluxmode == "lin") then
              usurf = 0.5 * (vara%u(i, j, k) + vara%u(i, j + 1, k)) * Pstrata(k)
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          end if

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

          if(topography) then
            pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j + 1, &
                &k) * pStratTFC(i, j + 1, k))
            pFEdgeF = 0.5 * (jac(i, j + 1, k) * pStratTFC(i, j + 1, k) &
                &+ jac(i, j + 2, k) * pStratTFC(i, j + 2, k))
            if(fluxmode == "nln") then
              vSurf = 0.5 * (pEdgeF * var%v(i, j, k) + pFEdgeF * var%v(i, j &
                  &+ 1, k))
            else if(fluxmode == "lin") then
              vSurf = 0.5 * (pEdgeF * vara%v(i, j, k) + pFEdgeF * vara%v(i, j &
                  &+ 1, k))
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          else
            if(fluxmode == "nln") then
              vsurf = 0.5 * (var%v(i, j, k) + var%v(i, j + 1, k)) * Pstrat(k)
            else if(fluxmode == "lin") then
              vsurf = 0.5 * (vara%v(i, j, k) + vara%v(i, j + 1, k)) * Pstrata(k)
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          end if

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

          if(topography) then
            pEdgeU = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
                &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            pFEdgeU = jac(i, j + 1, k) * jac(i, j + 1, k + 1) * (pStratTFC(i, &
                &j + 1, k) + pStratTFC(i, j + 1, k + 1)) / (jac(i, j + 1, k) &
                &+ jac(i, j + 1, k + 1))
            if(fluxmode == "nln") then
              wSurf = 0.5 * (pEdgeU * var%w(i, j, k) + pFEdgeU * var%w(i, j &
                  &+ 1, k))
            else if(fluxmode == "lin") then
              wSurf = 0.5 * (pEdgeU * vara%w(i, j, k) + pFEdgeU * vara%w(i, j &
                  &+ 1, k))
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          else
            if(fluxmode == "nln") then
              wsurf = 0.5 * (var%w(i, j, k) + var%w(i, j + 1, k)) &
                  &* PstratTilde(k)
            else if(fluxmode == "lin") then
              wsurf = 0.5 * (vara%w(i, j, k) + vara%w(i, j + 1, k)) &
                  &* PstratTildea(k)
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          end if

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

          if(topography) then
            pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, j, &
                &k) * pStratTFC(i + 1, j, k))
            pUEdgeR = 0.5 * (jac(i, j, k + 1) * pStratTFC(i, j, k + 1) + jac(i &
                &+ 1, j, k + 1) * pStratTFC(i + 1, j, k + 1))
            if(fluxmode == "nln") then
              uSurf = ((jac(i, j, k + 1) + jac(i + 1, j, k + 1)) * pEdgeR &
                  &* var%u(i, j, k) + (jac(i, j, k) + jac(i + 1, j, k)) &
                  &* pUEdgeR * var%u(i, j, k + 1)) / (jac(i, j, k) + jac(i &
                  &+ 1, j, k) + jac(i, j, k + 1) + jac(i + 1, j, k + 1))
            else if(fluxmode == "lin") then
              uSurf = ((jac(i, j, k + 1) + jac(i + 1, j, k + 1)) * pEdgeR &
                  &* vara%u(i, j, k) + (jac(i, j, k) + jac(i + 1, j, k)) &
                  &* pUEdgeR * vara%u(i, j, k + 1)) / (jac(i, j, k) + jac(i &
                  &+ 1, j, k) + jac(i, j, k + 1) + jac(i + 1, j, k + 1))
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          else
            if(fluxmode == "nln") then
              usurf = 0.5 * (var%u(i, j, k) * Pstrat(k) + var%u(i, j, k + 1) &
                  &* Pstrat(k + 1))
            else if(fluxmode == "lin") then
              usurf = 0.5 * (vara%u(i, j, k) * Pstrata(k) + vara%u(i, j, k &
                  &+ 1) * Pstrata(k + 1))
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          end if

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

          if(topography) then
            pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j + 1, &
                &k) * pStratTFC(i, j + 1, k))
            pUEdgeF = 0.5 * (jac(i, j, k + 1) * pStratTFC(i, j, k + 1) &
                &+ jac(i, j + 1, k + 1) * pStratTFC(i, j + 1, k + 1))
            if(fluxmode == "nln") then
              vSurf = ((jac(i, j, k + 1) + jac(i, j + 1, k + 1)) * pEdgeF &
                  &* var%v(i, j, k) + (jac(i, j, k) + jac(i, j + 1, k)) &
                  &* pUEdgeF * var%v(i, j, k + 1)) / (jac(i, j, k) + jac(i, j &
                  &+ 1, k) + jac(i, j, k + 1) + jac(i, j + 1, k + 1))
            else if(fluxmode == "lin") then
              vSurf = ((jac(i, j, k + 1) + jac(i, j + 1, k + 1)) * pEdgeF &
                  &* vara%v(i, j, k) + (jac(i, j, k) + jac(i, j + 1, k)) &
                  &* pUEdgeF * vara%v(i, j, k + 1)) / (jac(i, j, k) + jac(i, j &
                  &+ 1, k) + jac(i, j, k + 1) + jac(i, j + 1, k + 1))
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          else
            if(fluxmode == "nln") then
              vsurf = 0.5 * (var%v(i, j, k) * Pstrat(k) + var%v(i, j, k + 1) &
                  &* Pstrat(k + 1))
            else if(fluxmode == "lin") then
              vsurf = 0.5 * (vara%v(i, j, k) * Pstrata(k) + vara%v(i, j, k &
                  &+ 1) * Pstrata(k + 1))
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          end if

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

          if(topography) then
            pEdgeU = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
                &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            pUEdgeU = jac(i, j, k + 1) * jac(i, j, k + 2) * (pStratTFC(i, j, k &
                &+ 1) + pStratTFC(i, j, k + 2)) / (jac(i, j, k + 1) + jac(i, &
                &j, k + 2))
            if(fluxmode == "nln") then
              wSurf = 0.5 * (pEdgeU * var%w(i, j, k) + pUEdgeU * var%w(i, j, k &
                  &+ 1))
            else if(fluxmode == "lin") then
              wSurf = 0.5 * (pEdgeU * vara%w(i, j, k) + pUEdgeU * vara%w(i, j, &
                  &k + 1))
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          else
            if(fluxmode == "nln") then
              wsurf = 0.5 * (var%w(i, j, k) * PstratTilde(k) + var%w(i, j, k &
                  &+ 1) * PstratTilde(k + 1))
            else if(fluxmode == "lin") then
              wsurf = 0.5 * (vara%w(i, j, k) * PstratTildea(k) + vara%w(i, j, &
                  &k + 1) * PstratTildea(k + 1))
            else
              stop "Error in momentumFlux: Unknown fluxmode!"
            end if
          end if

          hRhoW = flux_muscl(wSurf, wD, wU)

          flux%w(i, j, k, 3) = hRhoW
        end do
      end do
    end do

    !-------------------------------------------------------------------
    !                          Viscous Fluxes
    !-------------------------------------------------------------------

    if(ReInv == 0.0 .and. .not. turbScheme) return

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
          if(model == "Boussinesq") then
            coef_v = ReInv

            if(turbScheme) then
              coef_v = coef_v + delta_hs * rho00 * var%DSC(i + 1, j, k)
            end if
          else
            if(topography) then
              coef_v = ReInv * rhoStratTFC(i + 1, j, 1)
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(turbScheme) then
              if(topography) then
                coef_v = coef_v + delta_hs * (var%rho(i + 1, j, k) &
                    &+ rhoStratTFC(i + 1, j, k)) * var%DSC(i + 1, j, k)
              else
                coef_v = coef_v + delta_hs * (var%rho(i + 1, j, k) &
                    &+ rhoStrat(k)) * var%DSC(i + 1, j, k)
              end if
            end if
          end if

          if(topography) then
            fRhoU_visc = coef_v * jac(i + 1, j, k) * stressTensTFC(i + 1, j, &
                &k, 1, 1, var)
          else
            du_dx = (var%u(i + 1, j, k) - var%u(i, j, k)) / dx
            div = divU(i + 1, j, k)
            fRhoU_visc = coef_v * (du_dx + du_dx - 2. / 3. * div)
          end if

          flux%u(i, j, k, 1) = flux%u(i, j, k, 1) - fRhoU_visc
        end do
      end do
    end do

    ! Flux gRhoU
    do k = 1, nz
      do j = 0, ny
        do i = 0, nx
          if(model == "Boussinesq") then
            coef_v = ReInv

            if(turbScheme) then
              coef_v = coef_v + delta_hs * rho00 * 0.25 * (var%DSC(i, j, k) &
                  &+ var%DSC(i + 1, j, k) + var%DSC(i, j + 1, k) + var%DSC(i &
                  &+ 1, j + 1, k))
            end if
          else
            if(topography) then
              coef_v = ReInv * 0.25 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i &
                  &+ 1, j, 1) + rhoStratTFC(i, j + 1, 1) + rhoStratTFC(i + 1, &
                  &j + 1, 1))
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(turbScheme) then
              if(topography) then
                coef_v = coef_v + delta_hs * 0.25 * ((var%rho(i, j, k) &
                    &+ rhoStratTFC(i, j, k)) * var%DSC(i, j, k) + (var%rho(i &
                    &+ 1, j, k) + rhoStratTFC(i + 1, j, k)) * var%DSC(i + 1, &
                    &j, k) + (var%rho(i, j + 1, k) + rhoStratTFC(i, j + 1, k)) &
                    &* var%DSC(i, j + 1, k) + (var%rho(i + 1, j + 1, k) &
                    &+ rhoStratTFC(i + 1, j + 1, k)) * var%DSC(i + 1, j + 1, k))
              else
                coef_v = coef_v + delta_hs * 0.25 * ((var%rho(i, j, k) &
                    &+ rhoStrat(k)) * var%DSC(i, j, k) + (var%rho(i + 1, j, k) &
                    &+ rhoStrat(k)) * var%DSC(i + 1, j, k) + (var%rho(i, j &
                    &+ 1, k) + rhoStrat(k)) * var%DSC(i, j + 1, k) &
                    &+ (var%rho(i + 1, j + 1, k) + rhoStrat(k)) * var%DSC(i &
                    &+ 1, j + 1, k))
              end if
            end if
          end if

          if(topography) then
            gRhoU_visc = coef_v * 0.25 * (jac(i, j, k) * stressTensTFC(i, j, &
                &k, 1, 2, var) + jac(i + 1, j, k) * stressTensTFC(i + 1, j, k, &
                &1, 2, var) + jac(i, j + 1, k) * stressTensTFC(i, j + 1, k, 1, &
                &2, var) + jac(i + 1, j + 1, k) * stressTensTFC(i + 1, j + 1, &
                &k, 1, 2, var))
          else
            du_dy = (var%u(i, j + 1, k) - var%u(i, j, k)) / dy
            dv_dx = (var%v(i + 1, j, k) - var%v(i, j, k)) / dx
            gRhoU_visc = coef_v * (du_dy + dv_dx)
          end if

          flux%u(i, j, k, 2) = flux%u(i, j, k, 2) - gRhoU_visc
        end do
      end do
    end do

    ! Flux hRhoU
    do k = 0, nz
      do j = 1, ny
        do i = 0, nx
          if(model == "Boussinesq") then
            coef_v = ReInv

            if(turbScheme) then
              if(topography) then
                coef_v = coef_v + delta_vs * rho00 * 0.5 * (jac(i, j, k) &
                    &* jac(i, j, k + 1) * (jac(i, j, k) * var%DSC(i, j, k) &
                    &+ jac(i, j, k + 1) * var%DSC(i, j, k + 1)) / (jac(i, j, &
                    &k) + jac(i, j, k + 1)) + jac(i + 1, j, k) * jac(i + 1, j, &
                    &k + 1) * (jac(i + 1, j, k) * var%DSC(i + 1, j, k) + jac(i &
                    &+ 1, j, k + 1) * var%DSC(i + 1, j, k + 1)) / (jac(i + 1, &
                    &j, k) + jac(i + 1, j, k + 1)))
              else
                coef_v = coef_v + delta_vs * rho00 * 0.25 * (var%DSC(i, j, k) &
                    &+ var%DSC(i + 1, j, k) + var%DSC(i, j, k + 1) + var%DSC(i &
                    &+ 1, j, k + 1))
              end if
            end if
          else
            if(topography) then
              coef_v = ReInv * 0.5 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i &
                  &+ 1, j, 1))
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(turbScheme) then
              if(topography) then
                coef_v = coef_v + delta_vs * 0.5 * (jac(i, j, k) * jac(i, j, k &
                    &+ 1) * (jac(i, j, k) * (var%rho(i, j, k) + rhoStratTFC(i, &
                    &j, k)) * var%DSC(i, j, k) + jac(i, j, k + 1) &
                    &* (var%rho(i, j, k + 1) + rhoStratTFC(i, j, k + 1)) &
                    &* var%DSC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k &
                    &+ 1)) + jac(i + 1, j, k) * jac(i + 1, j, k + 1) * (jac(i &
                    &+ 1, j, k) * (var%rho(i + 1, j, k) + rhoStratTFC(i + 1, &
                    &j, k)) * var%DSC(i + 1, j, k) + jac(i + 1, j, k + 1) &
                    &* (var%rho(i + 1, j, k + 1) + rhoStratTFC(i + 1, j, k &
                    &+ 1)) * var%DSC(i + 1, j, k + 1)) / (jac(i + 1, j, k) &
                    &+ jac(i + 1, j, k + 1)))
              else
                coef_v = coef_v + delta_vs * 0.25 * ((var%rho(i, j, k) &
                    &+ rhoStrat(k)) * var%DSC(i, j, k) + (var%rho(i + 1, j, k) &
                    &+ rhoStrat(k)) * var%DSC(i + 1, j, k) + (var%rho(i, j, k &
                    &+ 1) + rhoStrat(k + 1)) * var%DSC(i, j, k + 1) &
                    &+ (var%rho(i + 1, j, k + 1) + rhoStrat(k + 1)) &
                    &* var%DSC(i + 1, j, k + 1))
              end if
            end if
          end if

          if(topography) then
            stressTens13 = met(i, j, k, 1, 3) * stressTensTFC(i, j, k, 1, 1, &
                &var) + met(i, j, k, 2, 3) * stressTensTFC(i, j, k, 1, 2, var) &
                &+ stressTensTFC(i, j, k, 1, 3, var) / jac(i, j, k)
            stressTens13R = met(i + 1, j, k, 1, 3) * stressTensTFC(i + 1, j, &
                &k, 1, 1, var) + met(i + 1, j, k, 2, 3) * stressTensTFC(i + 1, &
                &j, k, 1, 2, var) + stressTensTFC(i + 1, j, k, 1, 3, var) &
                &/ jac(i + 1, j, k)
            stressTens13U = met(i, j, k + 1, 1, 3) * stressTensTFC(i, j, k &
                &+ 1, 1, 1, var) + met(i, j, k + 1, 2, 3) * stressTensTFC(i, &
                &j, k + 1, 1, 2, var) + stressTensTFC(i, j, k + 1, 1, 3, var) &
                &/ jac(i, j, k + 1)
            stressTens13RU = met(i + 1, j, k + 1, 1, 3) * stressTensTFC(i + 1, &
                &j, k + 1, 1, 1, var) + met(i + 1, j, k + 1, 2, 3) &
                &* stressTensTFC(i + 1, j, k + 1, 1, 2, var) + stressTensTFC(i &
                &+ 1, j, k + 1, 1, 3, var) / jac(i + 1, j, k + 1)
            hRhoU_visc = coef_v * 0.5 * (jac(i, j, k) * jac(i, j, k + 1) &
                &* (stressTens13 + stressTens13U) / (jac(i, j, k) + jac(i, j, &
                &k + 1)) + jac(i + 1, j, k) * jac(i + 1, j, k + 1) &
                &* (stressTens13R + stressTens13RU) / (jac(i + 1, j, k) &
                &+ jac(i + 1, j, k + 1)))
          else
            du_dz = (var%u(i, j, k + 1) - var%u(i, j, k)) / dz
            dw_dx = (var%w(i + 1, j, k) - var%w(i, j, k)) / dx
            hRhoU_visc = coef_v * (du_dz + dw_dx)
          end if

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
          if(model == "Boussinesq") then
            coef_v = ReInv

            if(turbScheme) then
              coef_v = coef_v + delta_hs * rho00 * 0.25 * (var%DSC(i, j, k) &
                  &+ var%DSC(i + 1, j, k) + var%DSC(i, j + 1, k) + var%DSC(i &
                  &+ 1, j + 1, k))
            end if
          else
            if(topography) then
              coef_v = ReInv * 0.25 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i &
                  &+ 1, j, 1) + rhoStratTFC(i, j + 1, 1) + rhoStratTFC(i + 1, &
                  &j + 1, 1))
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(turbScheme) then
              if(topography) then
                coef_v = coef_v + delta_hs * 0.25 * ((var%rho(i, j, k) &
                    &+ rhoStratTFC(i, j, k)) * var%DSC(i, j, k) + (var%rho(i &
                    &+ 1, j, k) + rhoStratTFC(i + 1, j, k)) * var%DSC(i + 1, &
                    &j, k) + (var%rho(i, j + 1, k) + rhoStratTFC(i, j + 1, k)) &
                    &* var%DSC(i, j + 1, k) + (var%rho(i + 1, j + 1, k) &
                    &+ rhoStratTFC(i + 1, j + 1, k)) * var%DSC(i + 1, j + 1, k))
              else
                coef_v = coef_v + delta_hs * 0.25 * ((var%rho(i, j, k) &
                    &+ rhoStrat(k)) * var%DSC(i, j, k) + (var%rho(i + 1, j, k) &
                    &+ rhoStrat(k)) * var%DSC(i + 1, j, k) + (var%rho(i, j &
                    &+ 1, k) + rhoStrat(k)) * var%DSC(i, j + 1, k) &
                    &+ (var%rho(i + 1, j + 1, k) + rhoStrat(k)) * var%DSC(i &
                    &+ 1, j + 1, k))
              end if
            end if
          end if

          if(topography) then
            fRhoV_visc = coef_v * 0.25 * (jac(i, j, k) * stressTensTFC(i, j, &
                &k, 2, 1, var) + jac(i + 1, j, k) * stressTensTFC(i + 1, j, k, &
                &2, 1, var) + jac(i, j + 1, k) * stressTensTFC(i, j + 1, k, 2, &
                &1, var) + jac(i + 1, j + 1, k) * stressTensTFC(i + 1, j + 1, &
                &k, 2, 1, var))
          else
            dv_dx = (var%v(i + 1, j, k) - var%v(i, j, k)) / dx
            du_dy = (var%u(i, j + 1, k) - var%u(i, j, k)) / dy
            fRhoV_visc = coef_v * (dv_dx + du_dy)
          end if

          flux%v(i, j, k, 1) = flux%v(i, j, k, 1) - fRhoV_visc
        end do
      end do
    end do

    ! Flux gRhoV
    do k = 1, nz
      do j = - 1, ny
        do i = 1, nx
          if(model == "Boussinesq") then
            coef_v = ReInv

            if(turbScheme) then
              coef_v = coef_v + delta_hs * rho00 * var%DSC(i, j + 1, k)
            end if
          else
            if(topography) then
              coef_v = ReInv * rhoStratTFC(i, j + 1, 1)
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(turbScheme) then
              if(topography) then
                coef_v = coef_v + delta_hs * (var%rho(i, j + 1, k) &
                    &+ rhoStratTFC(i, j + 1, k)) * var%DSC(i, j + 1, k)
              else
                coef_v = coef_v + delta_hs * (var%rho(i, j + 1, k) &
                    &+ rhoStrat(k)) * var%DSC(i, j + 1, k)
              end if
            end if
          end if

          if(topography) then
            gRhoV_visc = coef_v * jac(i, j + 1, k) * stressTensTFC(i, j + 1, &
                &k, 2, 2, var)
          else
            dv_dy = (var%v(i, j + 1, k) - var%v(i, j, k)) / dy
            div = divU(i, j + 1, k)
            gRhoV_visc = coef_v * (dv_dy + dv_dy - 2. / 3. * div)
          end if

          flux%v(i, j, k, 2) = flux%v(i, j, k, 2) - gRhoV_visc
        end do
      end do
    end do

    ! Flux hRhoV
    do k = 0, nz
      do j = 0, ny
        do i = 1, nx
          if(model == "Boussinesq") then
            coef_v = ReInv

            if(turbScheme) then
              if(topography) then
                coef_v = coef_v + delta_vs * rho00 * 0.5 * (jac(i, j, k) &
                    &* jac(i, j, k + 1) * (jac(i, j, k) * var%DSC(i, j, k) &
                    &+ jac(i, j, k + 1) * var%DSC(i, j, k + 1)) / (jac(i, j, &
                    &k) + jac(i, j, k + 1)) + jac(i, j + 1, k) * jac(i, j + 1, &
                    &k + 1) * (jac(i, j + 1, k) * var%DSC(i, j + 1, k) &
                    &+ jac(i, j + 1, k + 1) * var%DSC(i, j + 1, k + 1)) &
                    &/ (jac(i, j + 1, k) + jac(i, j + 1, k + 1)))
              else
                coef_v = coef_v + delta_vs * rho00 * 0.25 * (var%DSC(i, j, k) &
                    &+ var%DSC(i, j + 1, k) + var%DSC(i, j, k + 1) &
                    &+ var%DSC(i, j + 1, k + 1))
              end if
            end if
          else
            if(topography) then
              coef_v = ReInv * 0.5 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i, j &
                  &+ 1, 1))
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(turbScheme) then
              if(topography) then
                coef_v = coef_v + delta_vs * 0.5 * (jac(i, j, k) * jac(i, j, k &
                    &+ 1) * (jac(i, j, k) * (var%rho(i, j, k) + rhoStratTFC(i, &
                    &j, k)) * var%DSC(i, j, k) + jac(i, j, k + 1) &
                    &* (var%rho(i, j, k + 1) + rhoStratTFC(i, j, k + 1)) &
                    &* var%DSC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k &
                    &+ 1)) + jac(i, j + 1, k) * jac(i, j + 1, k + 1) * (jac(i, &
                    &j + 1, k) * (var%rho(i, j + 1, k) + rhoStratTFC(i, j + 1, &
                    &k)) * var%DSC(i, j + 1, k) + jac(i, j + 1, k + 1) &
                    &* (var%rho(i, j + 1, k + 1) + rhoStratTFC(i, j + 1, k &
                    &+ 1)) * var%DSC(i, j + 1, k + 1)) / (jac(i, j + 1, k) &
                    &+ jac(i, j + 1, k + 1)))
              else
                coef_v = coef_v + delta_vs * 0.25 * ((var%rho(i, j, k) &
                    &+ rhoStrat(k)) * var%DSC(i, j, k) + (var%rho(i, j + 1, k) &
                    &+ rhoStrat(k)) * var%DSC(i, j + 1, k) + (var%rho(i, j, k &
                    &+ 1) + rhoStrat(k + 1)) * var%DSC(i, j, k + 1) &
                    &+ (var%rho(i, j + 1, k + 1) + rhoStrat(k + 1)) &
                    &* var%DSC(i, j + 1, k + 1))
              end if
            end if
          end if

          if(topography) then
            stressTens23 = met(i, j, k, 1, 3) * stressTensTFC(i, j, k, 2, 1, &
                &var) + met(i, j, k, 2, 3) * stressTensTFC(i, j, k, 2, 2, var) &
                &+ stressTensTFC(i, j, k, 2, 3, var) / jac(i, j, k)
            stressTens23F = met(i, j + 1, k, 1, 3) * stressTensTFC(i, j + 1, &
                &k, 2, 1, var) + met(i, j + 1, k, 2, 3) * stressTensTFC(i, j &
                &+ 1, k, 2, 2, var) + stressTensTFC(i, j + 1, k, 2, 3, var) &
                &/ jac(i, j + 1, k)
            stressTens23U = met(i, j, k + 1, 1, 3) * stressTensTFC(i, j, k &
                &+ 1, 2, 1, var) + met(i, j, k + 1, 2, 3) * stressTensTFC(i, &
                &j, k + 1, 2, 2, var) + stressTensTFC(i, j, k + 1, 2, 3, var) &
                &/ jac(i, j, k + 1)
            stressTens23FU = met(i, j + 1, k + 1, 1, 3) * stressTensTFC(i, j &
                &+ 1, k + 1, 2, 1, var) + met(i, j + 1, k + 1, 2, 3) &
                &* stressTensTFC(i, j + 1, k + 1, 2, 2, var) &
                &+ stressTensTFC(i, j + 1, k + 1, 2, 3, var) / jac(i, j + 1, k &
                &+ 1)
            hRhoV_visc = coef_v * 0.5 * (jac(i, j, k) * jac(i, j, k + 1) &
                &* (stressTens23 + stressTens23U) / (jac(i, j, k) + jac(i, j, &
                &k + 1)) + jac(i, j + 1, k) * jac(i, j + 1, k + 1) &
                &* (stressTens23F + stressTens23FU) / (jac(i, j + 1, k) &
                &+ jac(i, j + 1, k + 1)))
          else
            dv_dz = (var%v(i, j, k + 1) - var%v(i, j, k)) / dz
            dw_dy = (var%w(i, j + 1, k) - var%w(i, j, k)) / dy
            hRhoV_visc = coef_v * (dv_dz + dw_dy)
          end if

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
          if(model == "Boussinesq") then
            coef_v = ReInv

            if(turbScheme) then
              if(topography) then
                coef_v = coef_v + delta_vs * rho00 * 0.5 * (jac(i, j, k) &
                    &* jac(i, j, k + 1) * (jac(i, j, k) * var%DSC(i, j, k) &
                    &+ jac(i, j, k + 1) * var%DSC(i, j, k + 1)) / (jac(i, j, &
                    &k) + jac(i, j, k + 1)) + jac(i + 1, j, k) * jac(i + 1, j, &
                    &k + 1) * (jac(i + 1, j, k) * var%DSC(i + 1, j, k) + jac(i &
                    &+ 1, j, k + 1) * var%DSC(i + 1, j, k + 1)) / (jac(i + 1, &
                    &j, k) + jac(i + 1, j, k + 1)))
              else
                coef_v = coef_v + delta_vs * rho00 * 0.25 * (var%DSC(i, j, k) &
                    &+ var%DSC(i + 1, j, k) + var%DSC(i, j, k + 1) + var%DSC(i &
                    &+ 1, j, k + 1))
              end if
            end if
          else
            if(topography) then
              coef_v = ReInv * 0.5 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i &
                  &+ 1, j, 1))
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(turbScheme) then
              if(topography) then
                coef_v = coef_v + delta_vs * 0.5 * (jac(i, j, k) * jac(i, j, k &
                    &+ 1) * (jac(i, j, k) * (var%rho(i, j, k) + rhoStratTFC(i, &
                    &j, k)) * var%DSC(i, j, k) + jac(i, j, k + 1) &
                    &* (var%rho(i, j, k + 1) + rhoStratTFC(i, j, k + 1)) &
                    &* var%DSC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k &
                    &+ 1)) + jac(i + 1, j, k) * jac(i + 1, j, k + 1) * (jac(i &
                    &+ 1, j, k) * (var%rho(i + 1, j, k) + rhoStratTFC(i + 1, &
                    &j, k)) * var%DSC(i + 1, j, k) + jac(i + 1, j, k + 1) &
                    &* (var%rho(i + 1, j, k + 1) + rhoStratTFC(i + 1, j, k &
                    &+ 1)) * var%DSC(i + 1, j, k + 1)) / (jac(i + 1, j, k) &
                    &+ jac(i + 1, j, k + 1)))
              else
                coef_v = coef_v + delta_vs * 0.25 * ((var%rho(i, j, k) &
                    &+ rhoStrat(k)) * var%DSC(i, j, k) + (var%rho(i + 1, j, k) &
                    &+ rhoStrat(k)) * var%DSC(i + 1, j, k) + (var%rho(i, j, k &
                    &+ 1) + rhoStrat(k + 1)) * var%DSC(i, j, k + 1) &
                    &+ (var%rho(i + 1, j, k + 1) + rhoStrat(k + 1)) &
                    &* var%DSC(i + 1, j, k + 1))
              end if
            end if
          end if

          if(topography) then
            fRhoW_visc = coef_v * 0.5 * (jac(i, j, k) * jac(i, j, k + 1) &
                &* (stressTensTFC(i, j, k, 3, 1, var) + stressTensTFC(i, j, k &
                &+ 1, 3, 1, var)) / (jac(i, j, k) + jac(i, j, k + 1)) + jac(i &
                &+ 1, j, k) * jac(i + 1, j, k + 1) * (stressTensTFC(i + 1, j, &
                &k, 3, 1, var) + stressTensTFC(i + 1, j, k + 1, 3, 1, var)) &
                &/ (jac(i + 1, j, k) + jac(i + 1, j, k + 1)))
          else
            dw_dx = (var%w(i + 1, j, k) - var%w(i, j, k)) / dx
            du_dz = (var%u(i, j, k + 1) - var%u(i, j, k)) / dz
            fRhoW_visc = coef_v * (dw_dx + du_dz)
          end if

          flux%w(i, j, k, 1) = flux%w(i, j, k, 1) - fRhoW_visc
        end do
      end do
    end do

    ! Flux gRhoW
    do k = 0, nz
      do j = 0, ny
        do i = 1, nx
          if(model == "Boussinesq") then
            coef_v = ReInv

            if(turbScheme) then
              if(topography) then
                coef_v = coef_v + delta_vs * rho00 * 0.5 * (jac(i, j, k) &
                    &* jac(i, j, k + 1) * (jac(i, j, k) * var%DSC(i, j, k) &
                    &+ jac(i, j, k + 1) * var%DSC(i, j, k + 1)) / (jac(i, j, &
                    &k) + jac(i, j, k + 1)) + jac(i, j + 1, k) * jac(i, j + 1, &
                    &k + 1) * (jac(i, j + 1, k) * var%DSC(i, j + 1, k) &
                    &+ jac(i, j + 1, k + 1) * var%DSC(i, j + 1, k + 1)) &
                    &/ (jac(i, j + 1, k) + jac(i, j + 1, k + 1)))
              else
                coef_v = coef_v + delta_vs * rho00 * 0.25 * (var%DSC(i, j, k) &
                    &+ var%DSC(i, j + 1, k) + var%DSC(i, j, k + 1) &
                    &+ var%DSC(i, j + 1, k + 1))
              end if
            end if
          else
            if(topography) then
              coef_v = ReInv * 0.5 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i, j &
                  &+ 1, 1))
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(turbScheme) then
              if(topography) then
                coef_v = coef_v + delta_vs * 0.5 * (jac(i, j, k) * jac(i, j, k &
                    &+ 1) * (jac(i, j, k) * (var%rho(i, j, k) + rhoStratTFC(i, &
                    &j, k)) * var%DSC(i, j, k) + jac(i, j, k + 1) &
                    &* (var%rho(i, j, k + 1) + rhoStratTFC(i, j, k + 1)) &
                    &* var%DSC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k &
                    &+ 1)) + jac(i, j + 1, k) * jac(i, j + 1, k + 1) * (jac(i, &
                    &j + 1, k) * (var%rho(i, j + 1, k) + rhoStratTFC(i, j + 1, &
                    &k)) * var%DSC(i, j + 1, k) + jac(i, j + 1, k + 1) &
                    &* (var%rho(i, j + 1, k + 1) + rhoStratTFC(i, j + 1, k &
                    &+ 1)) * var%DSC(i, j + 1, k + 1)) / (jac(i, j + 1, k) &
                    &+ jac(i, j + 1, k + 1)))
              else
                coef_v = coef_v + delta_vs * 0.25 * ((var%rho(i, j, k) &
                    &+ rhoStrat(k)) * var%DSC(i, j, k) + (var%rho(i, j + 1, k) &
                    &+ rhoStrat(k)) * var%DSC(i, j + 1, k) + (var%rho(i, j, k &
                    &+ 1) + rhoStrat(k + 1)) * var%DSC(i, j, k + 1) &
                    &+ (var%rho(i, j + 1, k + 1) + rhoStrat(k + 1)) &
                    &* var%DSC(i, j + 1, k + 1))
              end if
            end if
          end if

          if(topography) then
            gRhoW_visc = coef_v * 0.5 * (jac(i, j, k) * jac(i, j, k + 1) &
                &* (stressTensTFC(i, j, k, 3, 1, var) + stressTensTFC(i, j, k &
                &+ 1, 3, 1, var)) / (jac(i, j, k) + jac(i, j, k + 1)) + jac(i, &
                &j + 1, k) * jac(i, j + 1, k + 1) * (stressTensTFC(i, j + 1, &
                &k, 3, 1, var) + stressTensTFC(i, j + 1, k + 1, 3, 1, var)) &
                &/ (jac(i, j + 1, k) + jac(i, j + 1, k + 1)))
          else
            dw_dy = (var%w(i, j + 1, k) - var%w(i, j, k)) / dy
            dv_dz = (var%u(i, j, k + 1) - var%u(i, j, k)) / dz
            gRhoW_visc = coef_v * (dw_dy + dv_dz)
          end if

          flux%w(i, j, k, 2) = flux%w(i, j, k, 2) - gRhoW_visc
        end do
      end do
    end do

    ! Flux hRhoW
    do k = - 1, nz
      do j = 1, ny
        do i = 1, nx
          if(model == "Boussinesq") then
            coef_v = ReInv

            if(turbScheme) then
              coef_v = coef_v + delta_vs * rho00 * jac(i, j, k + 1) ** 2.0 &
                  &* var%DSC(i, j, k + 1)
            end if
          else
            if(topography) then
              coef_v = ReInv * rhoStratTFC(i, j, 1)
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(turbScheme) then
              if(topography) then
                coef_v = coef_v + delta_vs * jac(i, j, k + 1) ** 2.0 &
                    &* (var%rho(i, j, k + 1) + rhoStratTFC(i, j, k + 1)) &
                    &* var%DSC(i, j, k + 1)
              else
                coef_v = coef_v + delta_vs * (var%rho(i, j, k + 1) &
                    &+ rhoStrat(k + 1)) * var%DSC(i, j, k + 1)
              end if
            end if
          end if

          if(topography) then
            hRhoW_visc = coef_v * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
                &* stressTensTFC(i, j, k + 1, 3, 1, var) + jac(i, j, k + 1) &
                &* met(i, j, k + 1, 2, 3) * stressTensTFC(i, j, k + 1, 3, 2, &
                &var) + stressTensTFC(i, j, k + 1, 3, 3, var))
          else
            dw_dz = (var%w(i, j, k + 1) - var%w(i, j, k)) / dz
            div = divU(i, j, k + 1)
            hRhoW_visc = coef_v * (dw_dz + dw_dz - 2. / 3. * div)
          end if

          flux%w(i, j, k, 3) = flux%w(i, j, k, 3) - hRhoW_visc
        end do
      end do
    end do

    if(verbose) print *, "fluxes.f90/momentumFlux:  momentum fluxes fRhoU, &
        &fRhoV, fRhoW,  gRhoU, gRhoV, gRhoW hRhoU, hRhoV, hRhoW calculated."

  end subroutine momentumFlux

  !---------------------------------------------------------------------------

  subroutine init_fluxes
    !---------------------------------------------
    ! 1) set parameter for central or ILES flux
    ! 2) allocate flux module variables
    !---------------------------------------------

    ! local variables
    integer :: allocstat

    ! module variables

    ! rhoBar
    allocate(rhoBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhoBar"

    ! rhopBar
    allocate(rhopBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhopBar"

    ! rhoOld
    allocate(rhoOld(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "init_fluxes: alloc of rhoOld failed"

    ! rhopOld
    allocate(rhopOld(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "init_fluxes: alloc of rhopOld failed"

    ! uBar
    allocate(uBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "init_fluxes: could not allocate uBar"

    ! vBar
    allocate(vBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate vBar"

    ! wBar
    allocate(wBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate wBar"

    ! rhoTilde
    allocate(rhoTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
        &0:1), stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhoTilde"

    !UAB
    ! rhoTilde_mom
    allocate(rhoTilde_mom(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
        &0:1), stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhoTilde"
    !UAE

    ! rhopTilde
    allocate(rhopTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
        &0:1), stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhopTilde"

    ! uTilde
    allocate(uTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, 0:1), &
        &stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate uTilde"

    ! vTilde
    allocate(vTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, 0:1), &
        &stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate vTilde"

    ! wTilde
    allocate(wTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, 0:1), &
        &stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate wTilde"

    if(model == "compressible") then
      allocate(PBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate PBar"

      allocate(PTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
          &0:1), stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate PTilde"

      allocate(pinew(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate pinew"

    end if

    !SD added if
    if(include_ice) then
      !  nIceTilde
      allocate(nIceTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
          &0:1), stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate nIceTilde"

      ! qIceTilde
      allocate(qIceTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
          &0:1), stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate qIceTilde"

      ! qvTilde
      allocate(qvTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
          &0:1), stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate qvTilde"

      ! nAerTilde
      allocate(nAerTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
          &0:1), stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate nAerTilde"
    end if

    !SD

    if(include_ice) then

      allocate(IceBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate IceBar"

    end if

    allocate(tracerBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate tracerBar"

    allocate(tracerTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
        &0:1), stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate tracerTilde"

    ! TFC FJ
    if(topography) then

      allocate(uOldTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
          &+ nbz)), stat = allocstat)
      if(allocstat /= 0) stop "init_fluxes: alloc of uOldTFC failed"

      allocate(vOldTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
          &+ nbz)), stat = allocstat)
      if(allocstat /= 0) stop "init_fluxes: alloc of vOldTFC failed"

      allocate(wOldTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
          &+ nbz)), stat = allocstat)
      if(allocstat /= 0) stop "init_fluxes: alloc of wOldTFC failed"

    end if

    if(verbose) print *, "init_fluxes: rhoBar, uBar, vBar, wBar, thetaBar,  &
        &rhoTilde, uTilde, vTilde, wTilde allocated."

  end subroutine init_fluxes

  !---------------------------------------------------------------------------

  subroutine terminate_fluxes
    !-----------------------------------
    ! deallocates flux module variables
    !-----------------------------------

    ! local variables
    integer :: allocstat

    !---------------- deallocate variables -----------------------

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

    !UAB
    deallocate(rhoTilde_mom, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate rhoTilde"
    !UAE

    deallocate(rhopTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate rhopTilde"

    deallocate(uTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate uTilde"

    deallocate(vTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate vTilde"

    deallocate(wTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate wTilde"

    ! SK compressible model
    if(model == "compressible") then
      deallocate(PBar, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate PBar"

      deallocate(PTilde, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate PTilde"

      deallocate(pinew, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate pinew"

    end if

    !SD
    if(include_ice) then

      deallocate(nIceTilde, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate nIceTilde"

      deallocate(qIceTilde, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate qIceTilde"

      deallocate(qvTilde, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate qvTilde"

      deallocate(nAerTilde, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate nAerTilde"

    end if

    if(include_ice) then

      deallocate(IceBar, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate IceBar"

    end if

    if(include_tracer) then
      deallocate(tracerBar, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate tracerBar"

      deallocate(tracerTilde, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate tracerTilde"
    end if

    ! TFC FJ
    if(topography) then

      deallocate(uOldTFC, stat = allocstat)
      if(allocstat /= 0) stop "terminate_fluxes: dealloc of uOldTFC failed"

      deallocate(vOldTFC, stat = allocstat)
      if(allocstat /= 0) stop "terminate_fluxes: dealloc of vOldTFC failed"

      deallocate(wOldTFC, stat = allocstat)
      if(allocstat /= 0) stop "terminate_fluxes: dealloc of wOldTFC failed"

    end if

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

  ! ----------------------------------------------------

  subroutine iceFlux(vara, var, flux, fluxmode, Pstrata, PStratTildea)
    !---------------------------------------------------------------------
    ! computes the tracer flux at all cell edges using reconstructed values
    ! fluxmode = lin => linear flux, advecting velocities prescribed in
    !                   vara
    !            nln => nonlinear flux, advecting velocities from var
    !
    ! MUSCL assumes that the reconstructed tracer is tracer/P
    !---------------------------------------------------------------------
    implicit none

    ! in/out variables
    type(var_type), intent(in) :: vara, var
    character(len = *), intent(in) :: fluxmode

    type(flux_type), intent(inout) :: flux
    real, dimension(- 1:nz + 2), intent(in) :: PStrata, PStratTildea

    integer :: i, j, k !,l
    !    real :: rhoL,rhoR, uL,uR       ! L=Left i-1/2, R=Right i+1/2
    !    real :: rhoB,rhoF, vB,vF       ! B=Backward j-1/2, F=Forward j+1/2
    !    real :: rhoD,rhoU, wD,wU       ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurf, vSurf, wSurf ! velocities at cell surface

    real :: iceL, iceR ! L=Left i-1/2, R=Right i+1/2
    real :: fIce, gIce, hIce

    ! TFC FJ
    real :: rhoStratEdgeR, rhoStratEdgeF, rhoStratEdgeU
    real :: pEdgeR, pEdgeF, pEdgeU

    !    real  :: fRho, gRho, hRho

    ! avoid abs() for linerisation
    real :: delta
    real, parameter :: delta0 = 1.0e-6

    !   variables for the turbulence scheme
    real :: Pr_turb
    real :: coef_t, drho_dxi, dtht_dxi

    real :: delta_hs, delta_vs

    integer :: iVar, ii

    Pr_turb = 1.0 !0.5   !FS

    ! squared grid scales for the anisotropic turbulence scheme

    if(TurbScheme) then
      if(ny == 1 .and. nx == 1) then
        stop 'ERROR: turbulence scheme assumes either nx > 1 or ny > 1'
      else
        if(nx == 1) then
          delta_hs = dy ** 2 ! 2D problems in y and z
        else if(ny == 1) then
          delta_hs = dx ** 2 ! 2D problems in x and z
        else
          delta_hs = dx * dy ! 3D problems

          if(dx / dy > 10.) then
            print *, 'WARNING: dx/dy > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          elseif(dy / dx > 10.) then
            print *, 'WARNING: dy/dx > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          end if
        end if

        delta_vs = dz ** 2
      end if
    end if

    do iVar = 1, nVarIce

      !-----------------------------------------
      !       Zonal rho fluxes in x: f
      !-----------------------------------------

      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            if(iVar .eq. inN) then

              iceR = nIceTilde(i + 1, j, k, 1, 0)
              iceL = nIceTilde(i, j, k, 1, 1)

            elseif(iVar .eq. inQ) then

              iceR = qIceTilde(i + 1, j, k, 1, 0)
              iceL = qIceTilde(i, j, k, 1, 1)

            elseif(iVar .eq. inQv) then

              iceR = qvTilde(i + 1, j, k, 1, 0)
              iceL = qvTilde(i, j, k, 1, 1)

            end if

            if(topography) then
              pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, &
                  &j, k) * pStratTFC(i + 1, j, k))
              if(fluxmode == "nln") then
                uSurf = pEdgeR * var%u(i, j, k)
              else if(fluxmode == "lin") then
                uSurf = pEdgeR * vara%u(i, j, k)
              else
                stop "ERROR: wrong fluxmode"
              end if
            else
              if(fluxmode == "nln") then
                uSurf = var%u(i, j, k) * Pstrat(k) !UA
              else if(fluxmode == "lin") then
                uSurf = vara%u(i, j, k) * Pstrata(k) !UA
              else
                stop 'ERROR: worng fluxmode'
              end if
            end if

            fIce = flux_muscl(uSurf, iceL, iceR)

            flux%ICE(i, j, k, 1, iVar) = fIce
          end do
        end do
      end do

      !-----------------------------------------
      !    Meridional rho fluxes in y: g
      !-----------------------------------------

      do k = 1, nz
        do j = 0, ny
          do i = 1, nx
            if(iVar .eq. inN) then

              iceR = nIceTilde(i, j + 1, k, 2, 0)
              iceL = nIceTilde(i, j, k, 2, 1)

            elseif(iVar .eq. inQ) then

              iceR = qIceTilde(i, j + 1, k, 2, 0)
              iceL = qIceTilde(i, j, k, 2, 1)

            elseif(iVar .eq. inQv) then

              iceR = qvTilde(i, j + 1, k, 2, 0)
              iceL = qvTilde(i, j, k, 2, 1)

            end if

            if(topography) then
              pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                  &+ 1, k) * pStratTFC(i, j + 1, k))
              if(fluxmode == "nln") then
                vSurf = pEdgeF * var%v(i, j, k)
              else if(fluxmode == "lin") then
                vSurf = pEdgeF * vara%v(i, j, k)
              else
                stop "ERROR: wrong fluxmode"
              end if
            else
              if(fluxmode == "nln") then
                vSurf = var%v(i, j, k) * Pstrat(k) !UA
              else if(fluxmode == "lin") then
                vSurf = vara%v(i, j, k) * Pstrata(k) !UA
              else
                stop 'ERROR: worng fluxmode'
              end if
            end if

            gIce = flux_muscl(vSurf, iceL, iceR)

            flux%ICE(i, j, k, 2, iVar) = gIce
          end do
        end do
      end do

      !-----------------------------------------
      !      Vertical rho fluxes in z: h
      !-----------------------------------------

      do k = 0, nz
        do j = 1, ny
          do i = 1, nx
            if(iVar .eq. inN) then

              iceR = nIceTilde(i, j, k + 1, 3, 0)
              iceL = nIceTilde(i, j, k, 3, 1)

            elseif(iVar .eq. inQ) then

              iceR = qIceTilde(i, j, k + 1, 3, 0)
              iceL = qIceTilde(i, j, k, 3, 1)

            elseif(iVar .eq. inQv) then

              iceR = qvTilde(i, j, k + 1, 3, 0)
              iceL = qvTilde(i, j, k, 3, 1)

            end if

            if(topography) then
              pEdgeU = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
                  &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
              if(fluxmode == "nln") then
                wSurf = pEdgeU * var%w(i, j, k)
              else if(fluxmode == "lin") then
                wSurf = pEdgeU * vara%w(i, j, k)
              else
                stop "ERROR: wrong fluxmode"
              end if
            else
              if(fluxmode == "nln") then
                wSurf = var%w(i, j, k) * PstratTilde(k) !UA
              else if(fluxmode == "lin") then
                wSurf = vara%w(i, j, k) * PstratTildea(k) !UA
              else
                stop 'ERROR: worng fluxmode'
              end if
            end if

            hIce = flux_muscl(wSurf, iceL, iceR)

            flux%ICE(i, j, k, 3, iVar) = hIce
          end do
        end do
      end do

    end do ! ii

  end subroutine iceFlux

end module flux_module
