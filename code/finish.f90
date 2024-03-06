module finish_module

  use type_module

  implicit none

  private

  public :: terminate

  contains

  subroutine terminate(var, var0, var1, flux, force, source, dRho, dRhop, &
      dMom, dTheta, dIce, dTracer, tracerforce)
    !-------------------
    ! deallocate fields
    !-------------------

    ! in/out variables
    real, dimension(:, :, :, :), allocatable :: var, var0, var1, force, &
        source, dMom, dIce, tracerforce
    real, dimension (:, :, :, :, :), allocatable :: flux
    real, dimension (:, :, :), allocatable :: dRho, dRhop, dTheta, dTracer

    ! argument list
    integer :: allocstat

    !--------------- deallcoate grid -----------------------

    deallocate(x, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate x"

    deallocate(y, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate y"

    deallocate(z, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate z"

    !---------------- deallocate variables -----------------------

    deallocate(var, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate var"

    deallocate(var0, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate var0"

    deallocate(flux, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate flux"

    deallocate(force, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate force"

    deallocate(source, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate source"

    deallocate(dRho, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate dRho"

    deallocate(var1, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate var1"

    deallocate(dRhop, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate dRhop"

    deallocate(dTheta, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate dTheta"

    deallocate(dMom, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate dMom"

    if(include_ice) then
      deallocate(dIce, stat = allocstat)
      if(allocstat /= 0) stop "finish.f90: could not deallocate dIce"
    end if

    if(include_tracer) then
      deallocate(dTracer, stat = allocstat)
      if(allocstat /= 0) stop "finish.f90: could not deallocate dTracer"
    end if

    deallocate (tracerforce, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate tracerforce"

  end subroutine terminate

end module finish_module
