module finish_module

  use type_module

  implicit none

  private

  public :: terminate

  contains

  subroutine terminate(var, var0, var1, flux, flux0, dRho, dRhop, dMom)

    !-------------------
    ! deallocate fields
    !-------------------

    ! in/out variables
    type(var_type) :: var, var0, var1
    type(flux_type) :: flux, flux0
    real, dimension(:, :, :, :), allocatable :: dMom
    real, dimension(:, :, :), allocatable :: dRho, dRhop

    integer :: allocstat

    !--------------- deallocate grid -----------------------

    deallocate(x, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate x"

    deallocate(y, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate y"

    deallocate(z, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate z"

    !---------------- deallocate variables -----------------------

    call deallocate_var_type(var)

    call deallocate_var_type(var0)

    call deallocate_var_type(var1)

    call deallocate_flux_type(flux)

    call deallocate_flux_type(flux0)

    deallocate(dRho, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate dRho"

    deallocate(dRhop, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate dRhop"

    deallocate(dMom, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate dMom"

  end subroutine terminate

end module finish_module
