module finish_module

  use type_module
  
  implicit none

  private 

  public :: terminate

contains

  
  subroutine terminate (var, var0, var1, dRho, dRhop, dMom, dTheta)
    !-------------------
    ! deallocate fields 
    !-------------------

    ! in/out variables
    real, dimension(:,:,:,:), allocatable :: var,var0,var1,dMom
    real, dimension(:,:,:), allocatable :: dRho, dRhop, dTheta
 
    ! argument list
    integer :: allocstat

!--------------- deallcoate grid -----------------------

    deallocate(x,stat=allocstat)
    if(allocstat /= 0) stop"finish.f90: could not deallocate x"

    deallocate(y,stat=allocstat)
    if(allocstat /= 0) stop"finish.f90: could not deallocate y"

    deallocate(z,stat=allocstat)
    if(allocstat /= 0) stop"finish.f90: could not deallocate z"


!---------------- deallocate variables -----------------------

    deallocate(var,stat=allocstat)
    if(allocstat /= 0) stop"finish.f90: could not deallocate var"

    deallocate(var0,stat=allocstat)
    if(allocstat /= 0) stop"finish.f90: could not deallocate var"

    deallocate(var1,stat=allocstat)
    if(allocstat /= 0) stop"finish.f90: could not deallocate var"

    deallocate(dRho,stat=allocstat)
    if(allocstat /= 0) stop"finish.f90: could not deallocate var"

    deallocate(dRhop,stat=allocstat)
    if(allocstat /= 0) stop"finish.f90: could not deallocate var"

    deallocate(dTheta,stat=allocstat)
    if(allocstat /= 0) stop"finish.f90: could not deallocate dTheta"

    deallocate(dMom,stat=allocstat)
    if(allocstat /= 0) stop"finish.f90: could not deallocate var"



  end subroutine terminate


end module finish_module
