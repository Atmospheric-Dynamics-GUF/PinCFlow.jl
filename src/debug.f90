module debug_module

  use type_module
  use atmosphere_module
  !  use flux_module ! for rhoTilde

  implicit none

  private ! private module

  public :: initSymm1D
  public :: checkSymm1D
  public :: checkSymmTilde1D

  public :: initSymmetry
  public :: checkSymmetry

  public :: initSymmFlux
  public :: checkSymmFlux

  public :: initSymmTilde
  public :: checkSymmTilde

  private :: identic

  contains

  subroutine initSymm1D(phi, nx)
    !--------------------------------------------
    ! create random but symmetric data field phi
    ! with symmetry axis at nx/2
    !---------------------------------------------

    ! in/out variables
    real, dimension(1:nx), intent(inout) :: phi
    integer, intent(in) :: nx

    ! local vars
    integer :: m

    ! warning
    if(modulo(nx, 2) .ne. 0) stop "debug.f90/intiSymm1D: nx not even"

    ! reset var
    phi = 0.0

    ! create random field
    call random_seed ! create new numbers...
    call random_number(phi)

    ! make symmetric
    do m = 1, nx / 2
      phi(nx / 2 + m) = phi(nx / 2 - m + 1)
    end do

    !    print*,"debug.f90/initSymmetry: Fehler einbauen in w"
    !    var(4,1,6,2) = 1.0

  end subroutine initSymm1D

  ! ----------------------------------------------------------

  subroutine checkSymm1D(phi, nx)
    !------------------------------------------------
    ! check symmetry of phi with symmetry axis nx/2
    !------------------------------------------------

    ! in/out variables
    real, dimension(nx), intent(in) :: phi
    integer, intent(in) :: nx

    ! local vars
    integer :: m
    real :: dphi, phiR, phiL
    logical :: symmetric = .false.

    ! warning
    if(modulo(nx, 2) .ne. 0) stop "debug.f90/checkSymm1D: nx not even"

    print *, "debug.f90/checkSymm1D: checking phi for symmetry"

    ! check symmetry
    do m = 1, nx / 2

      phiR = phi(nx / 2 + m)
      phiL = phi(nx / 2 - m + 1)
      dphi = phiR - phiL
      if(identic(phiL, phiR)) then
        symmetric = .true.
      else
        print *, "i_left, i_right,j, k, dphi: ", nx / 2 - m + 1, nx / 2 + m, &
            &dphi
        symmetric = .false.
      end if

      if(.not. symmetric) exit

    end do

    ! output
    if(symmetric) then
      print *, "debug.f90/checkSymmTilde1D: phi is symmetric"
    else
      print *, "debug.f90/checkSymmTilde1D: phi is not symmetric"
      stop
    end if

  end subroutine checkSymm1D

  ! ----------------------------------------------------------

  subroutine checkSymmTilde1D(phiTilde, nx)
    !---------------------------------------------------
    ! check symmetry of reconstruced variable phiTilde
    ! with symmetry axis nx/2
    !---------------------------------------------------

    ! in/out variables
    real, dimension(nx, - 1:1), intent(in) :: phiTilde
    integer, intent(in) :: nx

    ! local vars
    integer :: m, lambda
    real :: dphiTilde
    logical :: symmetric = .false.
    real :: phiR, phiL

    ! warning
    if(modulo(nx, 2) .ne. 0) stop "debug.f90/checkSymmTilde1D: nx not even"

    print *, "debug.f90/checkSymmTilde1D: checking phiTilde for symmetry"
    print *, "debug.f90/checkSymmTilde1D: monitoring phiTilde at nx/2 left and &
        &nx/2+1 right:", nx / 2, nx / 2 + 1, phiTilde(nx / 2, - 1), &
        &phiTilde(nx / 2 + 1, 1)

    ! check symmetric
    do m = 1, nx / 2
      do lambda = - 1, 1, 2

        phiR = phiTilde(nx / 2 + m, lambda)
        phiL = phiTilde(nx / 2 - m + 1, - lambda)
        dphiTilde = phiR - phiL

        if(identic(phiL, phiR)) then
          symmetric = .true.
        else
          print *, "i_left, i_right,j, k, dphi: ", nx / 2 - m + 1, nx / 2 + m, &
              &dphiTilde
          symmetric = .false.
        end if

      end do ! lambda

      if(.not. symmetric) exit

    end do ! m

    ! output
    if(symmetric) then
      print *, "debug.f90/checkSymmTilde1D: phiTilde is symmetric"
    else
      print *, "debug.f90/checkSymmTilde1D: phiTilde is not symmetric"; print &
          &*, ""
      return
    end if

  end subroutine checkSymmTilde1D

  ! ----------------------------------------------------------

  subroutine initSymmTilde(phiTilde)
    !---------------------------------------------
    ! initialize random symmetric field phiTilde
    ! with symmetry axis nx/2
    !---------------------------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, 0:1), &
        &intent(inout) :: phiTilde

    ! local vars
    integer :: m
    logical :: collocated

    collocated = .true.

    ! warning
    print *, "debug.f90/initSymmTilde: phiTilde is reset randomly! collocated &
        &= ", collocated

    ! reset var
    phiTilde = 0.0

    ! create random field
    ! call random_seed ! create new numbers...
    call random_number(phiTilde)

    if(collocated) then
      ! make symmetric
      do m = 1, nx / 2 + nbx
        phiTilde(nx / 2 + m, :, :, 1, 1) = phiTilde(nx / 2 - m + 1, :, :, 1, 0) ! rho: right = left
        phiTilde(nx / 2 + m, :, :, 1, 0) = phiTilde(nx / 2 - m + 1, :, :, 1, 1) ! rho: left = right
      end do
    else
      do m = 1, nx / 2 + nbx
        stop "not implemented"
        ! u: staggered and antisymmetric
        ! phiTilde(nx/2+m,:,:,2) = -phiTilde(nx/2-m,:,:,2)   ! uTilde
      end do
    end if

    !    print*,"debug.f90/initSymmetry: Fehler einbauen in w"
    !    phiTilde(4,1,6,2) = 1.0

  end subroutine initSymmTilde

  ! ----------------------------------------------------------

  subroutine initSymmetry(var)
    !------------------------------------------------------------
    ! create random symmetric field var with symmetry axis nx/2
    !------------------------------------------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var

    ! local vars
    integer :: m

    ! warning
    print *, "debug.f90/initSymmetry: var is reset randomly!!!"

    ! reset var
    call reset_var_type(var)

    ! create random field
    ! call random_seed ! create new numbers...
    call random_number(var%rho)
    call random_number(var%u)
    call random_number(var%v)
    call random_number(var%w)
    call random_number(var%pi)

    ! make symmetric
    do m = 1, nx / 2 + nbx
      var%rho(nx / 2 + m, :, :) = var%rho(nx / 2 - m + 1, :, :) ! rho
      var%v(nx / 2 + m, :, :) = var%v(nx / 2 - m + 1, :, :) ! v
      var%w(nx / 2 + m, :, :) = var%w(nx / 2 - m + 1, :, :) ! w
      var%pi(nx / 2 + m, :, :) = var%pi(nx / 2 - m + 1, :, :) ! pi'

      ! u: staggered and antisymmetric
      var%u(nx / 2 + m, :, :) = - var%u(nx / 2 - m, :, :) ! u
    end do

    !    print*,"debug.f90/initSymmetry: Fehler einbauen in w"
    !    var(4,1,6,2) = 1.0

  end subroutine initSymmetry

  ! ----------------------------------------------------------

  subroutine initSymmFlux(flux)
    !--------------------------------------------------------------
    ! ceate random symmetric flux variable with symmetry axis nx/2
    !--------------------------------------------------------------

    ! in/out variables
    type(flux_type), intent(inout) :: flux
    ! flux(i,j,k,dir,iFlux)
    ! dir = 1..3 > f,g,h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW

    ! local vars
    integer :: m

    ! warning
    print *, "debug.f90/initSymmFlux: flux is reset randomly!!!"

    ! reset var
    call reset_flux_type(flux)

    ! create random field
    ! call random_seed ! create new numbers...
    call random_number(flux%rho)
    call random_number(flux%u)
    call random_number(flux%v)
    call random_number(flux%w)

    ! make anti-symmetric
    do m = 1, nx / 2
      flux%rho(nx / 2 + m, :, :, 1) = - flux%rho(nx / 2 - m, :, :, 1) ! frho
      flux%v(nx / 2 + m, :, :, 1) = - flux%v(nx / 2 - m, :, :, 1) ! fv
      flux%w(nx / 2 + m, :, :, 1) = - flux%w(nx / 2 - m, :, :, 1) ! fw
    end do

    ! staggered u + symmetric because of rho*u^2
    do m = 0, nx / 2
      flux%u(nx / 2 + m, :, :, 1) = flux%u(nx / 2 - m - 1, :, :, 1) ! fu
    end do

    !    print*,"debug.f90/initSymmetry: Fehler einbauen in w"
    !    flux(4,1,6,1,2) = 1.0

  end subroutine initSymmFlux

  ! ----------------------------------------------------------

  subroutine checkSymmTilde(phiTilde, pos)
    !--------------------------------------------
    ! check phiTilde for symmetry with axis nx/2
    !--------------------------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, 0:1), &
        &intent(in) :: phiTilde
    integer, intent(in) :: pos

    ! local vars
    integer :: i, j, k, m, lambda
    logical :: collocated, symmetric
    real :: dphi, phiR, phiL

    collocated = .true.
    symmetric = .false. ! init asymmetric to avoid errors

    ! warning
    print *, "debug.f90/checkSymmTilde: checking phiTilde for symmetry at &
        &pos", pos
    print *, "collocated = ", collocated

    if(collocated) then
      ! check symmetry
      do k = 1, nz
        do j = 1, ny

          do m = 1, nx / 2 + 1 ! Note: symmetry only for 1. ghost cell assured
            do lambda = 0, 1
              phiR = phiTilde(nx / 2 + m, j, k, 1, 1 - lambda)
              phiL = phiTilde(nx / 2 - m + 1, j, k, 1, lambda)
              dphi = phiR - phiL ! rho: right = left?
              if(identic(phiL, phiR)) then
                symmetric = .true.
              else
                print *, "i_left, i_right,j,k, dphi, phiL, phiR: ", nx / 2 - m &
                    &+ 1, nx / 2 + m, j, k, dphi, phiL, phiR
                symmetric = .false.
              end if

            end do ! lambda
            if(.not. symmetric) exit
          end do !m
          if(.not. symmetric) exit
        end do !j
        if(.not. symmetric) exit
      end do !k

    else ! staggered
      do k = - nbz, nz + nbz
        do j = - nby, ny + nby

          do m = 1, nx / 2 + nbx ! Note: symmetry for all ghost cells assured
            stop "not implemented"
            ! u: staggered and antisymmetric
            ! phiTilde(nx/2+m,j,k,2) = -phiTilde(nx/2-m,j,k,2)   ! uTilde
          end do !m

        end do ! j
      end do ! k
    end if

    ! output
    if(symmetric) then
      print *, "debug.f90/checkSymmTilde: phiTilde symmetric at pos", pos
    else
      print *, "debug.f90/checkSymmTilde: phiTilde not symmetric at pos", pos
      stop
    end if

  end subroutine checkSymmTilde

  ! ----------------------------------------------------------

  subroutine checkSymmFlux(flux, pos)
    !-----------------------------------------
    ! check symmetry for flux with axsis nx/2
    !-----------------------------------------

    ! in/out variables
    type(flux_type), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux)
    ! dir = 1..3 > f,g,h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW

    integer, intent(in) :: pos

    ! local vars
    integer :: i, j, k, m
    real :: dfrho, dfu, dfv, dfw
    logical :: symmetric

    ! intit vars
    symmetric = .true.
    dfrho = 0.0; dfu = 0.0; dfv = 0.0; dfw = 0.0

    ! warning
    print *, "debug.f90/checkSymmFlux: checking f-flux symmetry at pos", pos

    ! check anti-symmetry
    do k = - 1, nz
      do j = - 1, ny

        do m = 1, nx / 2
          dfrho = flux%rho(nx / 2 + m, j, k, 1) + flux%rho(nx / 2 - m, j, k, 1) ! frho
          if(dfrho .ne. 0.0) then
            print *, "i_left, i_right,j,k, dfrho: ", nx / 2 - m, nx / 2 + m, &
                &j, k, dfrho
            symmetric = .false.
          end if

          dfv = flux%v(nx / 2 + m, j, k, 1) + flux%v(nx / 2 - m, j, k, 1) ! fv
          if(dfv .ne. 0.0) then
            print *, "i_left, i_right,j,k, dfv: ", nx / 2 - m, nx / 2 + m, j, &
                &k, dfv
            symmetric = .false.
          end if

          dfw = flux%w(nx / 2 + m, j, k, 1) + flux%w(nx / 2 - m, j, k, 1) ! fw
          if(dfw .ne. 0.0) then
            print *, "i_left, i_right,j,k, dfw: ", nx / 2 - m, nx / 2 + m, j, &
                &k, dfw
            symmetric = .false.
          end if

        end do

        ! staggered u: f is symmetric because rho*u^2 is symmetric
        do m = 0, nx / 2
          dfu = flux%u(nx / 2 + m, j, k, 1) - flux%u(nx / 2 - m - 1, j, k, 1) ! fu
          if(dfu .ne. 0.0) then
            print *, "i_left, i_right,j,k, dfu: ", nx / 2 - m - 1, nx / 2 + m, &
                &j, k, dfu
            symmetric = .false.
          end if
        end do

      end do
    end do

    ! output
    if(symmetric) then
      print *, "debug.f90/checkSymmFlux: f-fluxes symmetric at pos", pos
    else
      print *, "debug.f90/checkSymmFlux: f-fluxes not symmetric at pos", pos
      stop
    end if

  end subroutine checkSymmFlux

  ! ----------------------------------------------------------

  subroutine checkSymmetry(var, pos)
    !-----------------------------------------
    ! check symmetry for var with axsis nx/2
    !-----------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var
    integer, intent(in) :: pos

    ! local vars
    integer :: i, j, k, m
    real :: drho, du, dv, dw, dp
    real :: rhoR, rhoL, uR, uL, vR, vL, wR, wL, pR, pL
    logical :: symmetric

    ! intit vars
    symmetric = .false.

    ! warning
    print *, "debug.f90/checkSymmetry: checking var symmetry at pos", pos

    ! check symmetry
    do k = - nbz, nz + nbz
      do j = - nby, ny + nby
        do m = 1, nx / 2 + nbx

          ! rho
          rhoR = var%rho(nx / 2 + m, j, k)
          rhoL = var%rho(nx / 2 - m + 1, j, k)
          drho = rhoR - rhoL
          if(identic(rhoR, rhoL)) then
            symmetric = .true.
          else
            symmetric = .false.
            print *, "i_left, i_right,j,k, drho: ", nx / 2 - m + 1, nx / 2 &
                &+ m, j, k, drho
          end if

          ! v
          vR = var%v(nx / 2 + m, j, k)
          vL = var%v(nx / 2 - m + 1, j, k)
          dv = vR - vL
          if(identic(vR, vL)) then
            symmetric = .true.
          else
            print *, "i_left, i_right,j,k, dv: ", nx / 2 - m + 1, nx / 2 + m, &
                &j, k, dv
            symmetric = .false.
          end if

          ! w
          wR = var%w(nx / 2 + m, j, k)
          wL = var%w(nx / 2 - m + 1, j, k)
          dw = wR - wL
          if(identic(wR, wL)) then
            symmetric = .true.
          else
            print *, "i_left, i_right,j,k, dw: ", nx / 2 - m + 1, nx / 2 + m, &
                &j, k, dw
            symmetric = .false.
          end if

          ! pi'
          pR = var%pi(nx / 2 + m, j, k)
          pL = var%pi(nx / 2 - m + 1, j, k)
          dp = pR - pL
          if(identic(pL, pR)) then
            symmetric = .true.
          else
            print *, "i_left, i_right,j,k, dp: ", nx / 2 - m + 1, nx / 2 + m, &
                &j, k, dp
            symmetric = .false.
          end if

          ! u: staggered and anti-symmetric
          uR = var%u(nx / 2 + m, j, k)
          uL = var%u(nx / 2 - m, j, k)
          du = uR + uL ! note: uR = -uL for anti-symmetric vars
          if(identic(uR, - uL)) then
            symmetric = .true.
          else
            print *, "i_left, i_right,j,k, du: ", nx / 2 - m, nx / 2 + m, j, &
                &k, du
            symmetric = .false.
          end if

          if(.not. symmetric) exit
        end do ! m
        if(.not. symmetric) exit
      end do ! j
      if(.not. symmetric) exit
    end do ! k

    ! output
    if(symmetric) then
      print *, "debug.f90/checkSymmetrie: vars symmetric at pos", pos
    else
      print *, "debug.f90/checkSymmetrie: vars not symmetric at pos", pos
      stop
    end if

  end subroutine checkSymmetry

  ! -----------------------------------------------------------------------

  function identic(a, b)
    !-------------------------------
    ! check whether a and b deviate
    !-------------------------------

    logical :: identic
    real, intent(in) :: a, b

    ! local vars
    real, parameter :: tol = 1.0e-15
    real, parameter :: epsi = 1.0
    real :: diff, av

    av = (a + b) / 2.0 + epsi
    diff = abs(a - b) / av

    if(diff < tol) then
      identic = .true.
    else
      identic = .false.
    end if

  end function identic

end module debug_module
