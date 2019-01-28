module finish_module

  use type_module
  
  implicit none

  private 

  public :: terminate

contains

  
  subroutine terminate (var, var0, dRho, dMom, dTheta, dIce)
    !-------------------
    ! deallocate fields 
    !-------------------

    ! in/out variables
    real, dimension(:,:,:,:), allocatable :: var,var0,dMom, dIce
    real, dimension(:,:,:), allocatable :: dRho, dTheta
 
    ! argument list
    integer :: allocstat

!--------------- deallcoate grid -----------------------

    deallocate(x,stat=allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate x"

    deallocate(y,stat=allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate y"

    deallocate(z,stat=allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate z"


!---------------- deallocate variables -----------------------

    deallocate(var,stat=allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate var"

    deallocate(var0,stat=allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate var0"

    deallocate(dRho,stat=allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate dRho"

    deallocate(dTheta,stat=allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate dTheta"

    deallocate(dMom,stat=allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate dMom"

    if (include_ice) then
      deallocate(dIce,stat=allocstat)
      if(allocstat /= 0) stop "finish.f90: could not deallocate dIce"
    end if

  end subroutine terminate


end module finish_module
