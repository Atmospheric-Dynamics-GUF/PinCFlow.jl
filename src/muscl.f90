module muscl_module
  !-----------------------------------------------------
  ! MUSCL reconstruction routines
  ! This module is independent of the rest of the code
  ! 1) it always assumes 1 additional cells on both ends
  ! of the array for interpolation
  !-----------------------------------------------------

  use type_module

  implicit none

  private

  public :: reconstruct_MUSCL
  public :: muscl_reconstruct1D_mcvariant
  public :: muscl_reconstruct3D

  contains

  subroutine reconstruct_MUSCL(uBar, uTilde, sizeX, sizeY, sizeZ, limiter)

    !-------------------------------------------------------------
    !  MUSCL reconstruction of a single 3D data field u
    !  of size uSizeX * uSizeY * sizeZ into data field uTilde with
    !  limiter of limterType - in ALL spatial directions
    !  Note: all fields have a single ghost cell at the boundary
    !-------------------------------------------------------------

    ! in/out
    integer, intent(in) :: sizeX, sizeY, sizeZ
    real, dimension(sizeX, sizeY, sizeZ), intent(in) :: uBar
    real, dimension(sizeX, sizeY, sizeZ, 1:3, 0:1), intent(out) :: uTilde
    ! 1:3 -> spation dimension to be reconstrueced, 1=x, 2=y, 3=z
    ! 0:1 -> 0 = left/right or backward/forward or down/up
    character(len = *), intent(in) :: limiter

    ! -------- reconstruction to all directions +/- 1/2 ----------

    call muscl_reconstruct3D(uBar, sizeX, sizeY, sizeZ, uTilde, limiter)

  end subroutine reconstruct_MUSCL

  !------------------------------------------------------------------

  subroutine muscl_reconstruct3D(u, sizeX, sizeY, sizeZ, uTilde, limiterType)

    !-------------------------------------------------------------
    !  MUSCL reconstruction of a single 3D data field u
    !  of size uSizeX * uSizeY * sizeZ into data field uTilde with
    !  limiter of limterType - in a SINGLE specified direction
    !  Note: all fields have a single ghost cell at the boundary
    !-------------------------------------------------------------

    ! in/out
    integer, intent(in) :: sizeX, sizeY, sizeZ
    real, dimension(sizeX, sizeY, sizeZ), intent(in) :: u
    real, dimension(sizeX, sizeY, sizeZ, 1:3, 0:1), intent(out) :: uTilde
    ! 1:3 -> spation dimension to be reconstrueced, 1=x, 2=y, 3=z
    ! 0:1 -> 0 = left/right or backward/forward or down/up
    character(len = *), intent(in) :: limiterType

    ! local fields
    real, dimension(sizeX) :: phiX
    real, dimension(sizeX, 0:1) :: phiTildeX
    real, dimension(sizeY) :: phiY
    real, dimension(sizeY, 0:1) :: phiTildeY
    real, dimension(sizeZ) :: phiZ
    real, dimension(sizeZ, 0:1) :: phiTildeZ

    ! local vars
    integer :: i, j, k

    select case(limiterType)

    case('MCVariant')

      ! reconstruction in x-direction
      do k = 2, sizeZ - 1
        do j = 2, sizeY - 1
          phiX = u(:, j, k)
          call muscl_reconstruct1D_mcvariant(phiX, sizeX, phiTildeX)
          uTilde(:, j, k, 1, :) = phiTildeX
        end do
      end do

      ! reconstruction in y-direction
      do k = 2, sizeZ - 1
        do i = 2, sizeX - 1
          phiY = u(i, :, k)
          call muscl_reconstruct1D_mcvariant(phiY, sizeY, phiTildeY)
          uTilde(i, :, k, 2, :) = phiTildeY
        end do
      end do

      ! reconstruction in z-direction
      do j = 2, sizeY - 1
        do i = 2, sizeX - 1
          phiZ = u(i, j, :)
          call muscl_reconstruct1D_mcvariant(phiZ, sizeZ, phiTildeZ)
          uTilde(i, j, :, 3, :) = phiTildeZ
        end do
      end do

    case default
      stop "muscl.f90/limit: unknown limiter type.Stop."
    end select

  end subroutine muscl_reconstruct3D

  !----------------------------------------------------------------------------

  subroutine muscl_reconstruct1D_mcvariant(phi, phiSize, phiTilde)

    !-------------------------------------
    !  MUSCL reconstruction of a single
    !  1D data field phi of size phiSize
    !  into data field phiTilde with
    !  limiter of limterType
    !-------------------------------------

    ! in/out
    integer, intent(in) :: phiSize
    real, dimension(phiSize), intent(in) :: phi
    real, dimension(phiSize, 0:1), intent(out) :: phiTilde

    ! local varibales
    real :: sigmaR, sigmaL, s
    real :: deltaR, deltaL, theta

    ! local vars
    integer :: i, j, k

    ! init phiTilde with
    phiTilde = 1000.0

    do i = 2, phiSize - 1
      deltaL = phi(i) - phi(i - 1)
      deltaR = phi(i + 1) - phi(i)

      if(deltaR == 0.) then
        phiTilde(i, 1) = phi(i)
        phiTilde(i, 0) = phi(i)
      else
        if(deltaL == 0.) then
          theta = deltaL / deltaR
          s = (2.0 + theta) / 3.0
          sigmaL = max(0.0, min(2 * theta, s, 2.0))

          phiTilde(i, 1) = phi(i) + 0.5 * sigmaL * deltaR
          phiTilde(i, 0) = phi(i)
        else
          theta = deltaL / deltaR

          s = (2.0 + theta) / 3.0
          sigmaL = max(0.0, min(2 * theta, s, 2.0))

          s = (2.0 + 1.0 / theta) / 3.0
          sigmaR = max(0.0, min(2 / theta, s, 2.0))

          phiTilde(i, 1) = phi(i) + 0.5 * sigmaL * deltaR
          phiTilde(i, 0) = phi(i) - 0.5 * sigmaR * deltaL
        end if
      end if
    end do

  end subroutine muscl_reconstruct1D_mcvariant

end module muscl_module
