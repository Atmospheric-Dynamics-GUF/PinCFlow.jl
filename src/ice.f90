module ice_module

  use type_module
  
  implicit none

  private 

  public :: SIce_crit
  public :: alpha_ice, gamma_ice, delta_ice
  public :: J0_ice, A_ice


contains

  real function SIce_crit(T)
    ! in/out variables
    real :: T
    ! parameter
    real :: s0,s1

   s1 = -3.4057422903191969E-3 
   s0 = 2.2525409204521596
   SIce_crit = s0 + s1*T

  end function SIce_crit

!----------------------------------------------

  real function alpha_ice(T)
    ! in/out variables
    real :: T

    alpha_ice = 1 ! to be corrected

  end function alpha_ice

!---------------------------------------------

  real function gamma_ice(T)
    ! in/out variables
    real :: T

    gamma_ice = 1 ! to be corrected

  end function gamma_ice

!---------------------------------------------
 
 real function delta_ice(T)
    ! in/out variables
    real :: T

    delta_ice = 1 ! to be corrected

  end function delta_ice

!---------------------------------------------

  real function J0_ice(T)
    ! in/out variables
    real :: T

    J0_ice = 37 ! to be corrected

  end function J0_ice

!---------------------------------------------

  real function A_ice(T)
    ! in/out variables
    real :: T

    A_ice = 337 ! to be corrected

  end function A_ice

!---------------------------------------------


end module ice_module
