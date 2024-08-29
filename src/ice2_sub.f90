module ice2_sub_module

  contains

  subroutine output_ice(i, j, k, iVar, var, field_prc)

    use type_module, ONLY:nx, ny, nz, nVar, nbx, nby, nbz, include_ice2, &
        &nVarIce, inN, inQ, inQv, master, model, fluctuationMode, topography, &
        &thetaRefRatio, timeScheme, L_ice, R_v, S_c, Dep, pSatIceRef, epsil0, &
        &epsil0hat, L_hat, master, include_testoutput, var_type

    use mpi_module
    use atmosphere_module, ONLY:PStrat01, PStratTilde01, PStrat, rhoStrat, &
        &piStrat, kappaInv, PStratTFC, piStratTfc, rhoStratTFC, gamma_1, p0, &
        &mRef, uRef, kappa, Rsp, thetaRef, tRef, g, pRef
    use boundary_module !, ONLY : setHalos, setBoundary, reconstruction
    !use flux_module

    implicit none

    type(var_type), intent(in) :: var

    ! local and global output field and record nbs.
    real * 4, dimension(nx, ny), intent(out) :: field_prc

    real :: rho, pres, temp, theta, psi, Qv, SIce, NIce, exn_p, A
    integer, intent(in) :: i, j, k, iVar

    !some diagnostics
    if(master .and. i == 1 .and. j == 1 .and. k == 1 .and. iVar == inN) then
      print *, ""
      print "(a)", repeat("-", 80)
      print *, " subroutine output_ice:  "
      print "(a)", repeat("-", 80)
      print *, ""
      print *, "max value N in master =", maxval(real(var%ICE2(:, :, :, inN)))
      print *, "max value Q in master =", maxval(real(var%ICE2(:, :, :, inQ)))
      print *, "max value Qv in master =", maxval(real(var%ICE2(:, :, :, inQv)))
      print *, ""
    end if

    if(model == "pseudo_incompressible") then

      if(fluctuationMode) then

        if(topography) then

          rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)

          theta = PstratTFC(i, j, k) / rho

          !if ( timeScheme == "semiimplicit" ) then
          !   print*, 'ice2Sources works only with explicit integration'
          !   stop
          !else
          !problems in \pi if heating switched on
          exn_p = var%pi(i, j, k) + (PstratTFC(i, j, k) / p0) ** gamma_1
          !end if

        else

          rho = var%rho(i, j, k) + rhoStrat(k)

          theta = Pstrat(k) / rho

          !problems in \pi if heating of background
          exn_p = var%pi(i, j, k) + (PStrat(k) / p0) ** gamma_1

        end if ! topography

      else

        print *, ' ice2Sources works only with fluctuationMode == T '
        stop

      end if ! fluctuationMode

      pres = p0 * exn_p ** kappaInv !kappaInv = c_p/R

      temp = theta * exn_p

      psi = exp(L_hat * (1. - 1. / (temp * thetaRefRatio))) !Psat_ice(temp)

      Qv = var%ICE2(i, j, k, inQv) ! Q_v = \rho q_v

      SIce = temp * Qv / psi / epsil0hat !SatRatio(Qv, temp, psi)

    else

      print *, ' ice2Sources works only with model == pseudo_incompressible '
      stop

    end if ! pseudo_inc

    if(iVar .eq. inN) then

      field_prc(i, j) = real(var%ICE2(i, j, k, iVar) / rho / mRef, kind = 4) ! N_v = \rho n

    elseif(iVar .eq. inQ) then

      field_prc(i, j) = real(SIce, kind = 4)

    elseif(iVar .eq. inQv) then

      field_prc(i, j) = real(var%ICE2(i, j, k, iVar) / rho, kind = 4) !Q_v = \rho q_v

      !if (include_testoutput) then
      !  field_prc(i, j) = real(var(i, j, k, iVar), kind = 4)
      !end if

    else
      print *, 'wrong value for iVar in output_ice in ice2_sub.f90'
      !field_prc(i, j) = real(var(i, j, k, iVar), kind = 4)

    end if

    !asymptotic solution
    !A = var(i, j, k, 4)*uRef*L_ice*g*kappa/Rsp/thetaRef**2/R_v
    !var0(i,j,k, inQ) = 2*A*S_c / (S_c - 1.) / (Dep * temp * psi / pres /tRef) / mRef
    !var0(i,j,k, inQ) = 2*A*S_c / ((S_c - 1.)* Dep * temp * pRef/pSatIceRef/epsil0 /tRef) / mRef

  end subroutine output_ice

end module ice2_sub_module
