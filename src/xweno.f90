module xweno_module

  ! This module is independent of the rest of the code
  ! I.E. it always assumes 2 additional cells on both ends
  ! of the array for interpolation

  use type_module
  use debug_module

  implicit none

  private
  ! all module variables are internal to the module

  public :: reconstruct_ALDM
  public :: reconstruct_SALD
  public :: reconstruct3D
  public :: reconstruct1D
  public :: reconstruct1D_old
  public :: init_xweno

  ! ------------- module variables  ---------------------
  real, dimension(1:3, 0:2, 0:2, - 1:1) :: alpha
  ! alpha(k,r,l,lambda*2)

  real, dimension(1:3, 0:2, - 1:1) :: gamma
  ! gamma(k,r,lambda)

  contains

  subroutine reconstruct_SALD(uBar, uTilde)
    !-----------------------------------------------------------------
    ! simplified ALDM (SALD) means: reconstruct a variable only once
    ! in each spatial direction onto the cell face
    !------------------------------------------------------------------

    ! in/out variables
    real, dimension(nxx, nyy, nzz), intent(in) :: uBar
    real, dimension(nxx, nyy, nzz, 1:3, 0:1), intent(out) :: uTilde
    ! uTilde(i,j,k,dir,Left/Right) with
    ! dir = 1,2,3 for x,y and z direction
    ! 0 = Left and 1 = Right

    ! local
    !     real, dimension(nxx,nyy,nzz)  :: phiBar
    real, dimension(nxx, nyy, nzz, - 1:1) :: phiTilde
    ! reconstruction to left edge (-1), cell centre (0), right edge (1)

    logical :: centered

    ! ------------- reconstruction to x +/- 1/2 --------------

    centered = .false.
    call reconstruct3D(uBar, "x", centered, phiTilde)
    uTilde(:, :, :, 1, :) = phiTilde(:, :, :, - 1:1:2)

    ! ------------- reconstruction to y +/- 1/2 --------------

    centered = .false.
    call reconstruct3D(uBar, "y", centered, phiTilde)
    uTilde(:, :, :, 2, :) = phiTilde(:, :, :, - 1:1:2)

    ! ------------- reconstruction to z +/- 1/2 --------------

    centered = .false.
    call reconstruct3D(uBar, "z", centered, phiTilde)
    uTilde(:, :, :, 3, :) = phiTilde(:, :, :, - 1:1:2)

  end subroutine reconstruct_SALD

  !---------------------------------------------------------------------------

  subroutine reconstruct_ALDM(uBar, uTilde)
    !-------------------------------------------------
    ! standard ALDM procedure:
    ! reconstruct variable uBar for all cell faces
    ! using 3 consequential reconstructions - 2 onto
    ! the cell center and the last to the wanted edge
    !-------------------------------------------------

    ! in/out variables
    real, dimension(nxx, nyy, nzz), intent(in) :: uBar
    real, dimension(nxx, nyy, nzz, 1:3, 0:1), intent(out) :: uTilde
    ! uTilde(i,j,k,dir,Left/Right) with
    ! dir = 1,2,3 for x,y and z direction
    ! 0 = Left and 1 = Right

    ! local
    real, dimension(nxx, nyy, nzz) :: phiBar
    real, dimension(nxx, nyy, nzz, - 1:1) :: phiTilde
    ! reconstruction to left edge (-1), cell centre (0), right edge (1)

    logical :: centered

    ! ------------- reconstruction to x +/- 1/2 --------------

    phiBar = uBar
    centered = .true.
    call reconstruct3D(phiBar, "y", centered, phiTilde)

    phiBar = phiTilde(:, :, :, 0)
    centered = .true.
    call reconstruct3D(phiBar, "z", centered, phiTilde)

    phiBar = phiTilde(:, :, :, 0)
    centered = .false.
    call reconstruct3D(phiBar, "x", centered, phiTilde)
    uTilde(:, :, :, 1, :) = phiTilde(:, :, :, - 1:1:2)

    ! ------------- reconstruction to y +/- 1/2 --------------

    phiBar = uBar
    centered = .true.
    call reconstruct3D(phiBar, "z", centered, phiTilde)

    phiBar = phiTilde(:, :, :, 0)
    centered = .true.
    call reconstruct3D(phiBar, "x", centered, phiTilde)

    phiBar = phiTilde(:, :, :, 0)
    centered = .false.
    call reconstruct3D(phiBar, "y", centered, phiTilde)
    uTilde(:, :, :, 2, :) = phiTilde(:, :, :, - 1:1:2)

    ! ------------- reconstruction to z +/- 1/2 --------------

    phiBar = uBar
    centered = .true.
    call reconstruct3D(phiBar, "x", centered, phiTilde)

    phiBar = phiTilde(:, :, :, 0)
    centered = .true.
    call reconstruct3D(phiBar, "y", centered, phiTilde)

    phiBar = phiTilde(:, :, :, 0)
    centered = .false.
    call reconstruct3D(phiBar, "z", centered, phiTilde)
    uTilde(:, :, :, 3, :) = phiTilde(:, :, :, - 1:1:2)

  end subroutine reconstruct_ALDM

  !---------------------------------------------------------------------------

  subroutine reconstruct3D(phiBar, direction, centered, phiTilde)
    !------------------------------------------------
    ! reconstruct 3D-field phiBar along "direction"
    ! onto the cell center if "centered = .true." or
    ! onto the cell edges if "centered = .false."
    !------------------------------------------------

    ! inout
    real, dimension(nxx, nyy, nzz), intent(in) :: phiBar
    character(len = 1), intent(in) :: direction
    logical, intent(in) :: centered
    real, dimension(nxx, nyy, nzz, - 1:1), intent(out) :: phiTilde

    ! local 1D fields
    integer :: i, j, k
    real, dimension(nxx) :: phiBarX
    real, dimension(nyy) :: phiBarY
    real, dimension(nzz) :: phiBarZ
    real, dimension(nxx, - 1:1) :: phiTildeX
    real, dimension(nyy, - 1:1) :: phiTildeY
    real, dimension(nzz, - 1:1) :: phiTildeZ

    select case(direction)

    case("x") ! interpolation in x
      do k = 1, nzz
        do j = 1, nyy
          phiBarX(:) = phiBar(:, j, k)
          call reconstruct1D(phiBarX, nxx, phiTildeX, centered)
          phiTilde(:, j, k, :) = phiTildeX(:, :)
        end do
      end do

    case("y") ! interpolation in y
      do k = 1, nzz
        do i = 1, nxx
          phiBarY(:) = phiBar(i, :, k)
          call reconstruct1D(phiBarY, nyy, phiTildeY, centered)
          phiTilde(i, :, k, :) = phiTildeY(:, :)
        end do
      end do

    case("z") ! interpolation in z
      do j = 1, nyy
        do i = 1, nxx
          phiBarZ(:) = phiBar(i, j, :)
          call reconstruct1D(phiBarZ, nzz, phiTildeZ, centered)
          phiTilde(i, j, :, :) = phiTildeZ(:, :)
        end do
      end do

    case default
      print *, "xweno.f90: direction not correctly speciefied. STOP"
      stop

    end select

  end subroutine reconstruct3D

  !---------------------------------------------------------------------------

  subroutine reconstruct1D(phiBar, n, phiTilde, centered)
    !--------------------------------------------------
    ! reconstruct 1D-field phiBar with dimenstion "n"
    ! onto cell center if "centered = .true." or
    ! onto cell faces  if "centered = .false."
    !--------------------------------------------------

    ! in/out variables
    integer, intent(in) :: n
    logical, intent(in) :: centered
    real, dimension(n), intent(in) :: phiBar
    real, dimension(n, - 1:1), intent(out) :: phiTilde

    ! local fields
    real, dimension(3:n, 1:3, 0:2, - 1:1) :: phiHat
    ! reconstructed data phiHat(i,k,r,lambda)
    ! k = 1 constant, k = 2 linear, k=3 quadratic
    ! r stencil: left, right, center shifted stencil
    ! lambda: reconstruction to center, left/right edge

    real, dimension(0:2, 3:n, 1:3, - 1:1) :: omega ! omega(r,i,k,lambda)
    ! dphi2(i) = ( phi(i+1) - phi(i) )^2
    real, dimension(1:n - 1) :: dPhi2

    real, dimension(3:n, 2:3) :: beta ! beta(i,k)   ...r is not necessary
    ! since beta(k,r,x(i)) = beta(k,r-1,x(i-1))

    real :: lSum, sSum, third
    integer :: i, k, r, l, s, lambda2
    integer :: l0, l1

    ! parameter
    real, parameter :: epsBeta = 1.0e-16

    ! improve computational efficiency
    ! beta(k,r,xi) = beta(k,r-1,xi-1)

    ! compute only interpolants at xi or at the edges xi +/- dx/2
    if(centered) then
      l0 = 0
      l1 = 0
    else
      l0 = - 1
      l1 = 1
    end if

    third = 1. / 3.
    ! init
    phiTilde = 10.0

    ! reconstruct phiHat for constant, lin. and quad. interpolation

    do lambda2 = l0, l1, 2
      i_loop: do i = 3, n - 2 ! apart from 2 ghost cells
        do k = 1, 3
          do r = 0, k - 1
            lSum = 0.0
            do l = 0, k - 1
              lSum = lSum + alpha(k, r, l, lambda2) * phiBar(i - r + l)
            end do
            phiHat(i, k, r, lambda2) = lSum
          end do
        end do
      end do i_loop
    end do

    ! calc jumps dPhi2
    do i = 1, n - 1
      dPhi2(i) = (phiBar(i + 1) - phiBar(i)) * (phiBar(i + 1) - phiBar(i))
    end do

    ! --------
    !  k = 1
    ! --------
    k = 1
    r = 0
    omega(r, :, k, :) = third

    ! --------
    !  k = 2
    ! --------
    k = 2

    ! beta - the smoothness measure
    do i = 3, n - 1
      ! beta(i,k) = 1./( epsBeta + dPhi2(i-1) )**2
      ! faster version:
      beta(i, k) = 1. / ((epsBeta + dPhi2(i - 1)) * (epsBeta + dPhi2(i - 1)))
      !       beta(i,k) = 1.0
    end do

    ! omega - the weight for ALDM/SALD
    do lambda2 = l0, l1, 2
      i_loop2: do i = 3, n - 2
        sSum = gamma(2, 0, lambda2) * beta(i + 1, k)
        sSum = sSum + gamma(2, 1, lambda2) * beta(i, k)
        r = 0
        omega(r, i, k, lambda2) = third * gamma(k, r, lambda2) * beta(i + 1, &
            &k) / sSum
        r = 1
        omega(r, i, k, lambda2) = third * gamma(k, r, lambda2) * beta(i, k) &
            &/ sSum
      end do i_loop2
    end do

    ! --------
    !  k = 3
    ! --------
    k = 3

    ! beta - the smoothness measure
    do i = 3, n
      beta(i, k) = 1. / (epsBeta + dPhi2(i - 2) + dPhi2(i - 1)) ** 2
      !       beta(i,k) = 1.0
    end do

    ! omega - the weight for ALDM/SALD
    do lambda2 = l0, l1, 2
      i_loop3: do i = 3, n - 2
        sSum = gamma(3, 0, lambda2) * beta(i + 2, k)
        sSum = sSum + gamma(3, 1, lambda2) * beta(i + 1, k)
        sSum = sSum + gamma(3, 2, lambda2) * beta(i, k)
        do r = 0, k - 1
          omega(r, i, k, lambda2) = third * gamma(k, r, lambda2) * beta(i + 2 &
              &- r, k) / sSum
        end do
      end do i_loop3
    end do

    ! ---------
    !  interpolation to left, middle, right of cell i ---------
    ! --------

    do lambda2 = l0, l1, 2
      i_loop4: do i = 3, n - 2
        lSum = 0.0
        do k = 1, 3
          do r = 0, k - 1
            lSum = lSum + omega(r, i, k, lambda2) * phiHat(i, k, r, lambda2)
          end do
        end do
        phiTilde(i, lambda2) = lSum
      end do i_loop4
    end do

  end subroutine reconstruct1D

  !---------------------------------------------------------------------------

  subroutine reconstruct1D_old(phiBar, n, phiTilde, centered)
    !------------------------------------------------
    ! old version of "reconstruct1D": closer to the
    ! literature but not optimized
    !------------------------------------------------

    ! in/out variables
    integer, intent(in) :: n
    logical, intent(in) :: centered
    real, dimension(n), intent(in) :: phiBar
    real, dimension(n, - 1:1), intent(out) :: phiTilde

    ! main
    real, dimension(1:3, 0:2, - 1:1) :: phiHat
    real, parameter :: epsBeta = 1.0e-16
    real, dimension(1:3, 0:2) :: beta
    real, dimension(1:3, 0:2, - 1:1) :: omega
    real :: lSum, sSum
    integer :: i, k, r, l, s, lambda2
    integer :: l0, l1

    ! improve computational efficiency
    ! beta(k,r,xi) = beta(k,r-1,xi-1)

    print *, "xweno.f90/reconstruct1D_old: phiBarX(n/2,n/2+1) = ", n / 2, n &
        &/ 2 + 1, phiBar(n / 2), phiBar(n / 2 + 1)

    ! compute only interpolants at xi or at the edges xi +/- dx/2
    if(centered) then
      l0 = 0
      l1 = 0
    else
      l0 = - 1
      l1 = 1
    end if

    ! init
    phiTilde = 10.0

    ! -------- constant, linear, quadratic interpolation terms -------------

    i_loop: do i = 3, n - 2 ! interpolation phiBar at position lambda2
      ! apart from 2 ghost cells
      do k = 1, 3
        do r = 0, k - 1
          do lambda2 = l0, l1, 2
            lSum = 0.0
            do l = 0, k - 1
              lSum = lSum + alpha(k, r, l, lambda2) * phiBar(i - r + l)
            end do
            phiHat(k, r, lambda2) = lSum
          end do
        end do
      end do

      ! beta - the smoothness measure

      do k = 1, 3
        do r = 0, k - 1
          lSum = epsBeta
          do l = - r, k - r - 2
            lSum = lSum + (phiBar(i + l + 1) - phiBar(i + l)) ** 2
          end do
          beta(k, r) = 1.0 / lSum ** 2
        end do
      end do

      ! omega - the weight for ALDM

      do k = 1, 3
        do r = 0, k - 1
          do lambda2 = l0, l1, 2
            sSum = 0.0
            do s = 0, k - 1
              sSum = sSum + beta(k, s) * gamma(k, s, lambda2)
            end do
            omega(k, r, lambda2) = 1.0 / 3.0 * beta(k, r) / sSum * gamma(k, r, &
                &lambda2)

          end do
        end do
      end do

      ! ---------- interpolation to left, middle, right of cell i ---------

      do lambda2 = l0, l1, 2
        lSum = 0.0
        do k = 1, 3 ! k = 1,K with K = 3 for ALDM/SALD: quadratic interpolation
          ! is highest order implemented
          do r = 0, k - 1 ! r = 0 (right stencil)
            lSum = lSum + phiHat(k, r, lambda2) * omega(k, r, lambda2)
          end do
        end do
        phiTilde(i, lambda2) = lSum
      end do

    end do i_loop

  end subroutine reconstruct1D_old

  !--------------------------------------------------------------------------

  subroutine init_xweno
    !-------------------------------------------------
    ! 1) compute weights alpha(), gamma()
    ! 2) set parameters sigmaX, sigmaY, sigmaZ, L_ijk
    !-------------------------------------------------

    ! variables for alpha
    real :: lambda
    real :: mSum, pSum, nProd
    integer :: k, r, l, m, n, p, lambda2

    ! variables for gammas
    real :: c1, c2, c3, c4

    ! variables for sigmas
    real :: L0, h0, L_ijk
    real, parameter :: s = 1. / 3.

    ! -------------- alpha(k,r,l,lambda) ----------------

    ! init alpha
    alpha = 1.0

    do lambda2 = - 1, 1, 1
      lambda = real(lambda2) / 2.0
      do k = 1, 3
        do r = 0, k - 1
          do l = 0, k - 1

            mSum = 0.0;
            do m = l + 1, k

              pSum = 0.0;
              do p = 0, k
                if(p == m) cycle

                nProd = 1.0;
                do n = 0, k
                  if(n == p .or. n == m) cycle
                  nProd = nProd * (lambda + r - n + 0.5);
                end do
                pSum = pSum + nProd;
              end do

              nProd = 1.0;
              do n = 0, k
                if(n == m) cycle
                nProd = nProd * (m - n);
              end do
              mSum = mSum + pSum / nProd;
            end do

            alpha(k, r, l, lambda2) = mSum;
          end do
        end do
      end do
    end do

    ! ------------------- gamma(k,r,lambda) -----------------

    c1 = 0.05003;
    c2 = 1.0;
    c3 = 0.01902;
    c4 = 0.0855;

    ! gamma(k,r+1,lambda*2+2) with k = 1,2,3,  r = 0,..k-1,
    ! lambda = -1/2, 0, +1/2
    k = 1
    r = 0
    lambda2 = - 1
    gamma(k, r, lambda2) = 1.0
    lambda2 = 0
    gamma(k, r, lambda2) = 1.0
    lambda2 = 1
    gamma(k, r, lambda2) = 1.0

    k = 2
    r = 0
    lambda2 = - 1
    gamma(k, r, lambda2) = 0.0
    lambda2 = 0
    gamma(k, r, lambda2) = 0.5
    lambda2 = 1
    gamma(k, r, lambda2) = 1.0

    r = 1
    lambda2 = - 1
    gamma(k, r, lambda2) = 1.0
    lambda2 = 0
    gamma(k, r, lambda2) = 0.5
    lambda2 = 1
    gamma(k, r, lambda2) = 0.0

    k = 3
    r = 0
    lambda2 = - 1
    gamma(k, r, lambda2) = 1.0 - c3 - c4
    lambda2 = 0
    gamma(k, r, lambda2) = 0.5 * (1.0 - c1)
    lambda2 = 1
    gamma(k, r, lambda2) = c3

    r = 1
    lambda2 = - 1
    gamma(k, r, lambda2) = c4
    lambda2 = 0
    gamma(k, r, lambda2) = c1
    lambda2 = 1
    gamma(k, r, lambda2) = c4

    r = 2
    lambda2 = - 1
    gamma(k, r, lambda2) = c3
    lambda2 = 0
    gamma(k, r, lambda2) = 0.5 * (1.0 - c1)
    lambda2 = 1
    gamma(k, r, lambda2) = 1.0 - c3 - c4

    ! ----------------- sigma ------------------------

    ! Info by S. Hickel, mail from 30.3.2011

    ! for momentum transport
    sigma0 = 0.06891
    !xxx
    !sigma0 = 0.0
    !print*,"WARNING: sigma0 = 0.0"
    sigmaX = sigma0
    sigmaY = sigma0
    sigmaZ = sigma0

    ! for mass transport
    sigmaC = 0.615
    !xxx
    !print*,"WARNING: sigmaC = 0.0"
    !sigmaC = 0.0

    !---------------------
    !  no longer used:
    !---------------------

    L_ijk = lx(1) - lx(0) ! "...is a local estimate of the current
    ! integral flow scale" (Hickel, 2008)
    ! larger L_ijk -> larger viscosity_

    ! original formula, cf. PhD thesis (2.41)
    ! L0/h0 = 32.0
    !sigmaX = sigma0 * (32.0 * dx/L_ijk)**(-s)
    !sigmaY = sigma0 * (32.0 * dy/L_ijk)**(-s)
    !sigmaZ = sigma0 * (32.0 * dz/L_ijk)**(-s)

  end subroutine init_xweno

end module xweno_module
