module timeScheme_module

  use type_module

  implicit none

  public

  !------------------------
  !   public subroutines
  !------------------------
  public :: init_timeScheme

  !-------------------------------
  !    public module variables
  !------------------------------

  integer :: RKstage
  real, dimension(3), public :: alphaRK
  real, dimension(3), public :: betaRK
  real, dimension(3), public :: stepFrac
  real, dimension(3, 3), public :: rk

  contains

  subroutine init_timeScheme
    !----------------------------------------------------
    ! initialize parameters for Runge-Kutta time schemes
    !----------------------------------------------------

    ! local variables
    real :: z1, z2, z3, z4, z5, z6, c2 ! TVD-RK3

    ! -------------------------------
    !        setup time scheme
    ! -------------------------------

    select case(timeScheme)

    case("semiimplicit") ! ref. Durran (????), Benaccio & Klein (2019),
      ! and Achatz et al (????)
      if(verbose) then
        print *, "update.f90/init_update: timeScheme = semiimplicit"
      end if
      alphaRK = (/0.0, - 5. / 9., - 153. / 128./)
      betaRK = (/1. / 3., 15. / 16., 8. / 15./)
      stepFrac = (/1. / 3., 5. / 12., 1. / 4./)
      ! set phi = t, F(phi) = 1 in time scheme to calc stepFrac
      nStages = 3
      timeSchemeType = "lowStorage"

    
    case default
      stop "update.f90/init_update: unknown timeScheme."

    end select

  end subroutine init_timeScheme

end module timeScheme_module
