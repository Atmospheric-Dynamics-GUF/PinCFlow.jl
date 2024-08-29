module finish_module

  use type_module

  implicit none

  private

  public :: terminate

  contains

  subroutine terminate(var, var0, var1, varG, source, flux, flux0, force, &
      &dRho, dRhop, dMom, dTheta, dIce, dTracer, tracerforce, dPot)
    !-------------------
    ! deallocate fields
    !-------------------

    ! in/out variables
    type(var_type) :: var, var0, var1, varG, source
    type(flux_type) :: flux, flux0
    real, dimension(:, :, :, :), allocatable :: force, dMom, dIce, tracerforce
    real, dimension(:, :, :), allocatable :: dRho, dRhop, dTheta, dTracer, dPot

    ! argument list
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

    call deallocate_var_type(varG)

    call deallocate_var_type(source)

    call deallocate_flux_type(flux)

    call deallocate_flux_type(flux0)

    deallocate(force, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate force"

    deallocate(dRho, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate dRho"

    deallocate(dRhop, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate dRhop"

    deallocate(dTheta, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate dTheta"

    deallocate(dMom, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate dMom"

    if(model == "compressible") then
      deallocate(dPot, stat = allocstat)
      if(allocstat /= 0) stop "finish.f90: could not deallocate dPot"
    end if

    if(include_ice) then
      deallocate(dIce, stat = allocstat)
      if(allocstat /= 0) stop "finish.f90: could not deallocate dIce"
    end if

    if(include_tracer) then
      deallocate(dTracer, stat = allocstat)
      if(allocstat /= 0) stop "finish.f90: could not deallocate dTracer"
    end if

    deallocate(tracerforce, stat = allocstat)
    if(allocstat /= 0) stop "finish.f90: could not deallocate tracerforce"

  end subroutine terminate

end module finish_module
