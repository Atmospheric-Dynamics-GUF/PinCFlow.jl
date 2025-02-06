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

  contains

  subroutine init_timeScheme

    !----------------------------------------------------
    ! Set parameters for Runge-Kutta time scheme
    !----------------------------------------------------

    alphaRK = (/0.0, - 5. / 9., - 153. / 128./)
    betaRK = (/1. / 3., 15. / 16., 8. / 15./)
    stepFrac = (/1. / 3., 5. / 12., 1. / 4./)
    nStages = 3

  end subroutine init_timeScheme

end module timeScheme_module
