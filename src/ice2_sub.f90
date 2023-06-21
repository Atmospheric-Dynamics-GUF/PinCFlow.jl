module opt_field_mod

!!$                   if(include_ice2) then
!!$                      if ( any(iVar == ofIce%indVarPF(1:ofIce%numVar))) then
!!$                         !We should put dimensions
!!$                         !field_prc(i,j) = real(var(i,j,k,iVar),kind=4)
!!$                         
!!$                         !.and. any(iVar == iVarIce)
!!$                         call output_ice(i, j, k, iVar, var, field_prc)
!!$                      end if
!!$                   end if

  implicit none
  integer, parameter ::  MaxNumVarOptFld = 5

  type :: opt_field

     integer :: numVar
     logical :: output
     integer :: indVarPF(MaxNumVarOptFld)
     
   !contains
   !  procedure :: init_opt_field
  end type opt_field

  type(opt_field) :: ofTracer, ofIce, ofSomething  

  contains

  subroutine init_opt_field(ofld, nvar, out, numOfld, indOfldVarM1)

    class(opt_field), intent(inout) :: ofld
    integer :: nvar !number of variables for opt_field
                    !tracer : 1; ice2 : 3 (n,q,q_v) 
    logical :: out ! True if output required

    integer :: numOfld, indOfldVarM1 !auxillary counters
    integer :: i

    if(nvar .gt. MaxNumVarOptFld ) then
       print*, 'opt_field%numVar > MaxNumVarOptFld !'
       print*, 'stop'
       stop
    end if 

    ofld%numVar = nvar
    ofld%output = out

    ofld%indVarPF(:) = 0 !init
    do i = 1, nvar
       ofld%indVarPF(i) = indOfldVarM1 + i
    end do

     numOfld = numOfld + 1
     indOfldVarM1 = indOfldVarM1 + nvar    

    print*, 'Init opt_field'
    
  end subroutine init_opt_field

  subroutine print_opt_field(ofld)
    class(opt_field), intent(in) :: ofld
    integer :: i 

    print*, 'produce output ', ofld%output
    print*, 'number variables ', ofld%numVar

    do i = 1, ofld%numVar
       print*, 'index of variable ',i, ' is', &
            ofld%indVarPF(i)
    end do
  end subroutine print_opt_field

!  subroutine set_opt_field(ofTracer, ofIce, ofSomething)
  subroutine set_opt_field

  use type_module, ONLY : include_tracer, include_ice2, master, nVar
  !logical, parameter :: include_tracer=.true., include_ice2=.false.  
  integer :: numOfld, indOfldVarM1 !auxillary counters
  !class(opt_field), intent(inout) :: ofTracer, ofIce, ofSomething

  open(555, file='index_opt_field.py')

  numOfld = 0      ! first no optional fields
  indOfldVarM1 = 8 ! number of "basic" PincFlow fields
 
  if (include_tracer) then 
     call init_opt_field(ofTracer, 1, .true., numOfld, indOfldVarM1)
     nVar = nVar + ofTracer%numVar
     if (master) then
     !call print_opt_field(ofTracer)
     write(555,"(a,i2.1)") 'iTr = ',  ofTracer%indVarPF(1)  
     end if
  end if
 
  if (include_ice2) then
     call init_opt_field(ofIce, 3, .true., numOfld, indOfldVarM1)
     nVar = nVar + ofIce%numVar
     if (master) then
     !call print_opt_field(ofIce) 
     write(555,"(a,i2.1)") 'iIce1 = ',  ofIce%indVarPF(1)  
     write(555,"(a,i2.1)") 'iIce2 = ',  ofIce%indVarPF(2)  
     write(555,"(a,i2.1)") 'iIce3 = ',  ofIce%indVarPF(3)  
     end if
  end if

  close(555)
 
end subroutine set_opt_field

end module opt_field_mod

module ice2_sub_module

  contains

  subroutine output_ice(i, j, k, iVar, var, field_prc)

    use type_module, ONLY:nx, ny, nz, nVar, nbx, nby, nbz, include_ice2, &
        nVarIce, inN, inQ, inQv, master, model, fluctuationMode, topography, &
        thetaRefRatio, timeScheme, L_ice, R_v, S_c, Dep, pSatIceRef, epsil0, &
        epsil0hat, L_hat

    use mpi_module
    use atmosphere_module, ONLY:PStrat01, PStratTilde01, PStrat, rhoStrat, &
        piStrat, kappaInv, PStratTFC, piStratTfc, rhoStratTFC, gamma_1, p0, &
        mRef, uRef, kappa, Rsp, thetaRef, tRef, g, pRef
    use boundary_module !, ONLY : setHalos, setBoundary, reconstruction
    !use flux_module

    implicit none

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(in) :: var

    ! local and global output field and record nbs.
    real * 4, dimension(nx, ny), intent(out) :: field_prc

    real :: rho, pres, temp, theta, psi, Qv, SIce, NIce, exn_p, A
    integer, intent(in) :: i, j, k, iVar

    if(model == "pseudo_incompressible") then

      if(fluctuationMode) then

        if(topography) then

          rho = var(i, j, k, 1) + rhoStratTFC(i, j, k)

          theta = PstratTFC(i, j, k) / rho

          !if ( timeScheme == "semiimplicit" ) then
          !   print*, 'ice2Sources works only with explicit integration'
          !   stop
          !else
          !problems in \pi if heating switched on
          exn_p = var(i, j, k, 5) + (PstratTFC(i, j, k) / p0) ** gamma_1
          !end if

        else

          rho = var(i, j, k, 1) + rhoStrat(k)

          theta = Pstrat(k) / rho

          !problems in \pi if heating of background
          exn_p = var(i, j, k, 5) + (PStrat(k) / p0) ** gamma_1

        end if ! topography

      else

        print *, ' ice2Sources works only with fluctuationMode == T '
        stop

      end if ! fluctuationMode

      pres = p0 * exn_p ** kappaInv !kappaInv = c_p/R

      temp = theta * exn_p

      psi = exp(L_hat * (1. - 1. / (temp * thetaRefRatio))) !Psat_ice(temp)

      Qv = var(i, j, k, inQv) ! Q_v = \rho q_v

      SIce = temp * Qv / psi / epsil0hat !SatRatio(Qv, temp, psi)

    else

      print *, ' ice2Sources works only with model == pseudo_incompressible '
      stop

    end if ! pseudo_inc

    if(iVar .eq. inN) then

      field_prc(i, j) = real(var(i, j, k, iVar) / rho / mRef, kind = 4) ! N_v = \rho n

    elseif(iVar .eq. inQ) then

      field_prc(i, j) = real(SIce, kind = 4)

    elseif(iVar .eq. inQv) then

      field_prc(i, j) = real(var(i, j, k, iVar) / rho, kind = 4) !Q_v = \rho q_v

    else

      field_prc(i, j) = real(var(i, j, k, iVar), kind = 4)

    end if

    !asymptotic solution
    !A = var(i, j, k, 4)*uRef*L_ice*g*kappa/Rsp/thetaRef**2/R_v
    !var0(i,j,k, inQ) = 2*A*S_c / (S_c - 1.) / (Dep * temp * psi / pres /tRef) / mRef
    !var0(i,j,k, inQ) = 2*A*S_c / ((S_c - 1.)* Dep * temp * pRef/pSatIceRef/epsil0 /tRef) / mRef

  end subroutine output_ice

end module ice2_sub_module
