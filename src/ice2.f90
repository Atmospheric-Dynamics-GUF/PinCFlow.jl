module ice2_module

  contains

  subroutine setup_ice2(var)

    use type_module
    use atmosphere_module, ONLY:heightTFC, tRef, rhoRef, lRef, thetaRef, pRef, &
        PStrat, rhoStrat, piStrat, kappaInv, PStratTFC, piStratTfc, &
        rhoStratTFC, gamma_1, p0, g, Rsp

    implicit none

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var

    integer, parameter :: ic_ice = 5

    real :: dz_tr
    real :: rho, exn_p, pres, temp, theta, psi

    integer :: iVar, ii, k, hzn, dzn, i, j
    real :: n0, q0, qv0, S0
    integer :: irec

    !case 3
    real :: z0_issr, sig_issr, S_issr

    !case 5
    real :: presMean, exn_pMean, tempMean, thetaMean, psiMean

    !set variables
    J_nuc = 4.9E4 !nucleation rate
    B_nuc = 337. !nucleation exponent
    Dep = 4.3E-8 !C_0 coeff paper
    epsil0hat = epsil0 * PsatIceRef / pRef

    !nondimensionalize variables/ fields

    thetaRefRatio = thetaRef / thetaRef_trp

    L_hat = L_ice / R_v / thetaRef_trp

    Li_hat = L_ice / R_v / thetaRef

    if(referenceQuantities .eq. "Klein") then
      !non-dimensional adiabatic lapse rate: Gamma = g/c_p
      alr_ndim = g / (7. * Rsp / 2.) * lRef / thetaRef
    else
      print *, 'uRef not eq. aRef'
      print *, 'check non-dimensionalization alr_ndim and others'
      stop
    end if

    mRef = rhoRef * lRef ** 3 !reverence mass

    J_nuc = J_nuc * tRef * mRef

    Dep = Dep * thetaRef * tRef * meanMassIce ** (1. / 3.) * PsatIceRef / pRef &
        / mRef

    !call compare_reduced_model
    !if ( master ) stop

    ! init to zero
    do ii = 1, nVarIce
      iVar = iVarIce(ii)

      var(:, :, :, iVar) = 0.

    end do

    select case(ic_ice)

    case(1)

      !!$       if(master) then
      !!$          open(44,file='test_output.dat',form="unformatted",access='direct',&
      !!$               & recl=nx*ny)
      !!$       end if

      ! non-dimensional variables
      ! top-hat/kink distribution
      hzn = nz / 2
      dzn = nz / 8
      do k = 1, nz
        !if ( k .ge. hzn .and. k .le. hzn+dzn ) then

        do j = 1, ny
          do i = 1, nx

            if(topography) then

              rho = var(i, j, k, 1) + rhoStratTFC(i, j, k)

              theta = PstratTFC(i, j, k) / rho

              !problems in \pi
              !if heating switched on
              !if ( timeScheme == "semiimplicit" ) then
              !   print*, 'ice2Sources works only with explicit integration'
              !   stop
              !else
              exn_p = var(i, j, k, 5) + (PstratTFC(i, j, k) / p0) ** gamma_1
              !end if

            else

              rho = var(i, j, k, 1) + rhoStrat(k)

              !problems in \pi
              !if heating switched on
              theta = Pstrat(k) / rho

              exn_p = var(i, j, k, 5) + (PStrat(k) / p0) ** gamma_1

            end if ! topography

            pres = p0 * exn_p ** kappaInv !kappaInv = c_p/R

            temp = theta * exn_p

            psi = Psat_ice(temp)

            ! IC asymptotic model
            !n = 0.1 * 2.E6
            !S = 1.4
            !q = meanMassIce * n

            !dimensional IC for n, q_v, q
            n0 = 0.1 * 2.E6 ![kg**-1]
            S0 = 1.4 !1.48
            qv0 = epsil0hat * S0 * psi / pres ! [kg/kg]
            q0 = meanMassIce * n0

            !double kink profile
            !S0 = 1+0.5*tanh(5*(z(k)-8/lRef)/2)+ (1-0.5*tanh(5*(z(k)-10/lRef)/5)) -1.4

            var(:, :, k, inN) = rho * n0 * mRef !\hat N = \hat \rho \hat n
            var(:, :, k, inQv) = rho * qv0
            var(:, :, k, inQ) = rho * q0

          end do
        end do

        !end if !k in range

      end do !k

    case(2)
      !linear profile
      dz_tr = 0.1 * dz ! increment tracer

      if(topography) then
        if(master) then
          print *, 'NEW TRACER PROF'
        end if
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx

              rho = var(i, j, k, 1) + rhoStratTFC(i, j, k)
              ! N = \rho n
              var(:, :, k, inN) = dz_tr * heightTFC(i, j, k) * rho

            end do
          end do
        end do

      else

        do k = 1, nz

          ! N = \rho n
          !var(:,:,k, inN) = (k-1)*0.1
          var(:, :, k, inN) = dz_tr * (z(k) - z(1)) * (var(:, :, k, 1) &
              + rhoStrat(k))

        end do

      end if

    case(3)

      !ISSRegion

      !center ISSR
      z0_issr = 8.e3 ! [m]
      !vertical width ISSR (standard deviation of gaussian dist.)
      sig_issr = 1.e3 ! [m]
      !max value S in ISSR
      S_issr = 1.45

      !nondim.
      z0_issr = z0_issr / lRef
      sig_issr = sig_issr / lRef

      do k = 1, nz

        do j = 1, ny
          do i = 1, nx

            if(topography) then

              rho = var(i, j, k, 1) + rhoStratTFC(i, j, k)

              theta = PstratTFC(i, j, k) / rho

              !problems in \pi
              !if heating switched on
              !if ( timeScheme == "semiimplicit" ) then
              !   print*, 'ice2Sources works only with explicit integration'
              !   stop
              !else
              exn_p = var(i, j, k, 5) + (PstratTFC(i, j, k) / p0) ** gamma_1
              !end if

            else

              rho = var(i, j, k, 1) + rhoStrat(k)

              !problems in \pi
              !if heating switched on
              theta = Pstrat(k) / rho

              exn_p = var(i, j, k, 5) + (PStrat(k) / p0) ** gamma_1

            end if ! topography

            pres = p0 * exn_p ** kappaInv !kappaInv = c_p/R

            temp = theta * exn_p

            psi = Psat_ice(temp)

            ! IC asymptotic model
            !n = 0.1 * 2.E6
            !S = 1.4
            !q = meanMassIce * n

            !dimensional IC for n, q_v, q
            n0 = 0. !0.1 * 2.E6 ![kg**-1]
            S0 = S_issr * exp(- (z(k) - z0_issr) ** 2 / 2. / sig_issr ** 2)
            qv0 = epsil0hat * S0 * psi / pres ! [kg/kg]
            q0 = meanMassIce * n0

            var(i, j, k, inN) = rho * n0 * mRef !\hat N = \hat \rho \hat n
            var(i, j, k, inQv) = rho * qv0
            var(i, j, k, inQ) = rho * q0

          end do !i
        end do !j

      end do !k
      !end case 3

    case(4)

      !ISSRegion

      !center ISSR
      z0_issr = 8.e3 ! [m]
      !vertical width ISSR (standard deviation of gaussian dist.)
      sig_issr = 4.e3 ! [m]
      !max value S in ISSR
      S_issr = 1.45

      !nondim.
      z0_issr = z0_issr / lRef
      sig_issr = sig_issr / lRef

      do k = 1, nz

        do j = 1, ny
          do i = 1, nx

            if(topography) then

              rho = var(i, j, k, 1) + rhoStratTFC(i, j, k)

              theta = PstratTFC(i, j, k) / rho

              !problems in \pi
              !if heating switched on
              !if ( timeScheme == "semiimplicit" ) then
              !   print*, 'ice2Sources works only with explicit integration'
              !   stop
              !else
              exn_p = var(i, j, k, 5) + (PstratTFC(i, j, k) / p0) ** gamma_1
              !end if

            else

              rho = var(i, j, k, 1) + rhoStrat(k)

              !problems in \pi
              !if heating switched on
              theta = Pstrat(k) / rho

              exn_p = var(i, j, k, 5) + (PStrat(k) / p0) ** gamma_1

            end if ! topography

            pres = p0 * exn_p ** kappaInv !kappaInv = c_p/R

            temp = theta * exn_p

            psi = Psat_ice(temp)

            ! IC asymptotic model
            !n = 0.1 * 2.E6
            !S = 1.4
            !q = meanMassIce * n

            !dimensional IC for n, q_v, q
            n0 = 0. !0.1 * 2.E6 ![kg**-1]
            S0 = S_issr * exp(- (z(k) - z0_issr) ** 2 / 2. / sig_issr ** 2)
            qv0 = epsil0hat * S0 * psi / pres ! [kg/kg]
            q0 = meanMassIce * n0

            var(i, j, k, inN) = rho * n0 * mRef !\hat N = \hat \rho \hat n
            var(i, j, k, inQv) = rho * qv0
            var(i, j, k, inQ) = rho * q0

            if(include_testoutput) then

              var(i, j, k, iVarO) = S0
              var(i, j, k, iVarO + 1) = pres * pRef
              var(i, j, k, iVarO + 2) = pres / psi * pRef / PsatIceRef

              if(raytracer) then
                var(i, j, k, iVarO) = S0
                var(i, j, k, iVarO + 1) = 0.
                var(i, j, k, iVarO + 2) = 0.
              end if

            end if

          end do !i
        end do !j

      end do !k

    case(5)

      !ISSRegion

      !center ISSR
      z0_issr = 8.e3 ! [m]
      !vertical width ISSR (standard deviation of gaussian dist.)
      sig_issr = 4.e3 ! [m]
      !max value S in ISSR
      S_issr = 1.45

      !nondim.
      z0_issr = z0_issr / lRef
      sig_issr = sig_issr / lRef

      do k = 1, nz

        do j = 1, ny
          do i = 1, nx

            if(topography) then

              rho = var(i, j, k, 1) + rhoStratTFC(i, j, k)

              theta = PstratTFC(i, j, k) / rho

              !problems in \pi
              !if heating switched on
              !if ( timeScheme == "semiimplicit" ) then
              !   print*, 'ice2Sources works only with explicit integration'
              !   stop
              !else
              exn_p = var(i, j, k, 5) + (PstratTFC(i, j, k) / p0) ** gamma_1
              !end if

            else

              rho = var(i, j, k, 1) + rhoStrat(k)

              !problems in \pi
              !if heating switched on
              theta = Pstrat(k) / rho

              thetaMean = Pstrat(k) / rhoStrat(k)

              exn_p = var(i, j, k, 5) + (PStrat(k) / p0) ** gamma_1

              exn_pMean = (PStrat(k) / p0) ** gamma_1

            end if ! topography

            pres = p0 * exn_p ** kappaInv !kappaInv = c_p/R

            presMean = p0 * exn_pMean ** kappaInv !kappaInv = c_p/R

            temp = theta * exn_p

            tempMean = thetaMean * exn_pMean

            psi = Psat_ice(temp)

            psiMean = Psat_ice(tempMean)

            ! IC asymptotic model
            !n = 0.1 * 2.E6
            !S = 1.4
            !q = meanMassIce * n

            !dimensional IC for n, q_v, q
            n0 = 0. !0.1 * 2.E6 ![kg**-1]
            S0 = S_issr * exp(- (z(k) - z0_issr) ** 2 / 2. / sig_issr ** 2)

            !qv0 = epsil0hat * S0 * psi / pres ! [kg/kg]
            !CHANGES
            qv0 = epsil0hat * S0 * psiMean / presMean ! [kg/kg]

            q0 = meanMassIce * n0

            var(i, j, k, inN) = rho * n0 * mRef !\hat N = \hat \rho \hat n
            var(i, j, k, inQv) = rho * qv0
            var(i, j, k, inQ) = rho * q0

            if(include_testoutput) then

              var(i, j, k, iVarO) = S0
              var(i, j, k, iVarO + 1) = pres / psi - presMean / psiMean ! store IC \tilde p_i(0)
              !var(i, j, k, iVarO+1) = pres*pRef
              var(i, j, k, iVarO + 2) = pres / psi * pRef / PsatIceRef

              if(raytracer) then
                var(i, j, k, iVarO) = S0
                var(i, j, k, iVarO + 1) = 0.
                var(i, j, k, iVarO + 2) = 0.
              end if

            end if

          end do !i
        end do !j

      end do !k
      !end case(5)

    end select

  end subroutine setup_ice2

  subroutine ice2Sources(var, source)

    use type_module, ONLY:nx, ny, nz, nVar, nbx, nby, nbz, include_ice2, &
        nVarIce, inN, inQ, inQv, master, model, fluctuationMode, topography, &
        thetaRefRatio, timeScheme, rayTracer, L_ice, R_v, var_ww, iVarO, &
        no_ice_source, ofield, pSatIceRef
    use mpi_module
    use atmosphere_module, ONLY:PStrat01, PStratTilde01, PStrat, rhoStrat, &
        piStrat, kappaInv, PStratTFC, piStratTfc, rhoStratTFC, gamma_1, p0, g, &
        Rsp, pRef
    use boundary_module !, ONLY : setHalos, setBoundary, reconstruction
    use flux_module
    !   use update_module, ONLY : dIce, ice2Update

    implicit none

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var
    !    real, dimension(-1:nx,-1:ny,-1:nz, nVarIce), intent(out) :: source
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(out) :: source
    real :: rho, pres, temp, theta, psi, Qv, SIce, NIce, exn_p
    real :: amp_pi, w_gw, PiMean, PiPrime, delta_t !RayTracer coupling
    integer :: i, j, k
    real :: pi0

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

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

            psi = Psat_ice(temp)

            Qv = var(i, j, k, inQv) ! Q_v = \rho q_v

            NIce = var(i, j, k, inN) ! N_v = \rho n

            !CHANGES
            !SIce = SatRatio(Qv, temp, psi)
            pi0 = var(i, j, k, iVarO + 1)
            SIce = SatRatio_rm_ini_pi(Qv, temp, psi, rho, pi0)

            !output variables
            var(i, j, k, iVarO) = SIce
            !var(i, j, k, iVarO+1) = pres*pRef
            var(i, j, k, iVarO + 2) = pres / psi * pRef / PsatIceRef
          else

            print *, ' ice2Sources works only with model  &
                == pseudo_incompressible '
            stop

          end if ! pseudo_inc

          if(.not. no_ice_source) then
            source(i, j, k, inN) = nucleation_n(SIce, rho)
            source(i, j, k, inQv) = deposition_qv(SIce, NIce, temp, pres, psi)
            source(i, j, k, inQ) = - source(i, j, k, inQv)
          end if

        end do ! i
      end do ! j
    end do ! k

  end subroutine ice2Sources

  subroutine ice2Sources_raytr(var, source, uTime, pTime)

    use type_module, ONLY:nx, ny, nz, nVar, nbx, nby, nbz, include_ice2, &
        nVarIce, inN, inQ, inQv, master, model, fluctuationMode, topography, &
        thetaRefRatio, timeScheme, rayTracer, L_ice, R_v, var_ww, iVarO, &
        no_ice_source, ofield, PsatIceRef, Li_hat, alr_ndim

    use type_module, ONLY:PI ! TEST
    use atmosphere_module, ONLY:tRef !TEST

    use mpi_module
    use atmosphere_module, ONLY:PStrat01, PStratTilde01, PStrat, rhoStrat, &
        piStrat, kappaInv, PStratTFC, piStratTfc, rhoStratTFC, gamma_1, p0, g, &
        Rsp, uRef, pRef
    use boundary_module !, ONLY : setHalos, setBoundary, reconstruction
    use flux_module
    !   use update_module, ONLY : dIce, ice2Update

    implicit none

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var
    !    real, dimension(-1:nx,-1:ny,-1:nz, nVarIce), intent(out) :: source
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(out) :: source
    real :: rho, pres, temp, theta, psi, Qv, SIce, NIce, exn_p
    real :: amp_pi, w_gw, PiMean, PiPrime, delta_t !RayTracer coupling
    real :: uTime, pTime
    integer :: i, j, k
    real :: dphi, omg

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

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

            psi = Psat_ice(temp)

            Qv = var(i, j, k, inQv) ! Q_v = \rho q_v

            NIce = var(i, j, k, inN) ! N_v = \rho n

            SIce = SatRatio(Qv, temp, psi)

            !var(i, j, k, iVarO) = SIce

            ! modify saturation ratio S
            ! to account for GW tendency in p/p_si
            !
            ! NB: we assume constant GW vertical velocity during
            ! each pincflow time step, this is paticular critical
            ! for high-frequency GWs and large integration time steps

            !Pi = p/p_si
            !Pi = <Pi> + Pi'
            !<Pi> large scale field, Pi' GW fluctuations

            ! var_ww : get from MSGWam
            ! amp_pi =
            delta_t = uTime - pTime
            if(rayTracer) then

              PiMean = pres / psi

              !amp_pi = g * pres / temp /psi / (7.*Rsp/2.) * &
              ! (L_ice / R_v / temp - kappaInv)
              !PiPrime(t_0+dt) = PiPrime(t_0) + amp_pi * var(i, j, k, iVarO+1) * delta_t
              ! work with STD of w
              !PiPrime = ofield(i, j, k, 1) + amp_pi * ofield(i, j, k, 5) * delta_t

              ! work with monochromatic wave w = A_w cos(delta_phi + wt )
              ! \dot p_i = B cos(delta_phi-wt )
              ! p_i(t) = p_i(t0) + B/\omega [ sin(delta_phi + wt ) ]^t_t0
              !
              !compute phase
              !delta_phi = kx + ly + mz + phi_0
              !this should work only for initial time t=0
              !(and if ray volumes do not leave cells and omega=const)

              !nondimensional version
              omg = ofield(i, j, k, 2)
              dphi = ofield(i, j, k, 3)

              amp_pi = (Li_hat / temp - kappaInv) * alr_ndim * PiMean / temp &
                  * ofield(i, j, k, 6) / omg

              if(amp_pi .gt. 0.) then
                PiPrime = ofield(i, j, k, 1) + amp_pi * (cos(dphi - omg &
                    * uTime) - cos(dphi - omg * pTime))
              else
                PiPrime = 0.
              end if

              !save PiPrime(t_0+dt)
              ofield(i, j, k, 1) = PiPrime

              !CHNAGES
              !SIce = PiPrime / PiMean
              SIce = SIce * (1 + PiPrime / PiMean)

            end if

            !output variables
            var(i, j, k, iVarO) = SIce
            !STD vert. vel.
            !var(i, j, k, iVarO+1) = ofield(i, j, k, 5)*uRef !\tilde w

            !max GW vert. vel.
            if(amp_pi .gt. 0.) then
              var(i, j, k, iVarO + 1) = ofield(i, j, k, 6) * uRef * sin(dphi &
                  - omg * uTime) ! w_max
            else
              var(i, j, k, iVarO + 1) = 0.
            end if
            var(i, j, k, iVarO + 2) = ofield(i, j, k, 1) * pRef / PsatIceRef !p/p_si
            !var(i, j, k, iVarO+2) = ofield(i, j, k, 1)
          else

            print *, ' ice2Sources works only with model  &
                == pseudo_incompressible '
            stop

          end if ! pseudo_inc

          if(.not. no_ice_source) then
            ! CHANGES
            source(i, j, k, inN) = nucleation_n(SIce, rho)
            source(i, j, k, inQv) = deposition_qv(SIce, NIce, temp, pres, psi)
            source(i, j, k, inQ) = - source(i, j, k, inQv)
          end if

        end do ! i
      end do ! j
    end do ! k

  end subroutine ice2Sources_raytr

  function deposition_qv(S, N, T, p, p_si)
    ! compute deposition rate
    ! - D p_si/ p (S-1) T N

    use type_module, ONLY:Dep

    implicit none

    real, intent(in) :: S, T, N, p, p_si
    real :: deposition_qv

    deposition_qv = - Dep * (S - 1.) * T * N * p_si / p

  end function deposition_qv

  function source_S(S, N, T, p, p_si, rho)
    ! compute rhs of eq. for S in asymptotic model
    ! - D p_si/ p (S-1) T N + A S cos (w t)

    use type_module, ONLY:Dep, pSatIceRef, epsil0
    use atmosphere_module, ONLY:heightTFC, tRef, rhoRef, lRef, thetaRef, pRef

    implicit none

    real, intent(in) :: S, T, N, p, p_si, rho
    real :: source_S

    source_S = - Dep * (S - 1.) * T * N * pRef / pSatIceRef / epsil0 / rho

  end function source_S

  function nucleation_n(S, rho)
    ! compute nucleation rate
    ! \rho * J * exp[B(S-S_c)]

    use type_module, ONLY:J_nuc, B_nuc, S_c

    implicit none

    real, intent(in) :: S, rho
    real :: nucleation_n

    nucleation_n = rho * J_nuc * exp(B_nuc * (S - S_c))

  end function nucleation_n

  function SatRatio(Qv, T, p_si)
    ! compute saturation ratio
    ! S= p * q_v / p_si / epsil0 = R T \rho q_v / p_si / epsil0 =  R T Qv / p_si / epsil0
    use type_module, ONLY:epsil0hat
    implicit none

    real, intent(in) :: Qv, T, p_si
    real :: SatRatio

    SatRatio = Qv * T / p_si / epsil0hat

  end function SatRatio

  function SatRatio_rm_ini_pi(Qv, T, p_si, rho, pi0)
    ! compute saturation ratio, remove wave contribution \tilde p_i(t_0)=p/p_si(t_0)
    ! S= p * q_v / p_si / epsil0 = R T \rho q_v / p_si / epsil0 =  R T Qv / p_si / epsil0
    use type_module, ONLY:epsil0hat
    implicit none

    real, intent(in) :: Qv, T, p_si, rho, pi0
    real :: SatRatio_rm_ini_pi

    SatRatio_rm_ini_pi = Qv * (T / p_si - pi0 / rho) / epsil0hat

  end function SatRatio_Rm_Ini_Pi

  function Psat_ice(T)
    ! compute saturation pressure over ice
    ! NB: temperature nondimensionalized with thetaRef_tropopause
    ! bc. Psat_ice_ref = Psat_ice_ref (thetaRef_tropopause)
    use type_module, ONLY:L_hat, thetaRefRatio

    implicit none

    real, intent(in) :: T
    real :: Psat_ice

    Psat_ice = exp(L_hat * (1. - 1. / (T * thetaRefRatio)))

  end function Psat_Ice

  subroutine compare_reduced_model
    !compare tendency with reduced model from pyton script

    use type_module
    use atmosphere_module, ONLY:heightTFC, tRef, rhoRef, lRef, thetaRef, pRef, &
        PStrat01, PStratTilde01, PStrat, rhoStrat, piStrat, kappaInv, &
        PStratTFC, piStratTfc, rhoStratTFC, gamma_1

    implicit none
    real :: temp, pres, psi, rho, n, S, q, dn, ds, dq

    !set values form paper
    temp = 210 / thetaRef
    pres = 30000 / pRef
    psi = Psat_ice(temp) !non. dimen. sat. pressure over ice
    print *, 'PSatRef', psi, PSatIceRef, pRef, epsil0
    rho = 0.5 / rhoRef

    !initial conditions from Model in XX/kg
    n = 0.1 * 2.E6
    S = 1.4
    q = meanMassIce * n

    !transorm to \hat N , S, \hat Q
    !\hat N = \hat \rho * \hat n
    n = rho * (n * mRef)
    !\hat Q = \hat \rho * \hat q
    q = rho * q

    dn = nucleation_n(S, rho)
    ds = source_S(S, n, temp, pres, psi, rho)

    !compute mass specific tendencies
    dn = dn / rho
    !dq = dq/rho

    !dimensionalize tendencies
    dn = dn / mRef / tRef
    ds = ds / tRef

    print *, 'tendency n', dn
    print *, 'tendency S', ds
    !print*, 'tendency q', dq

    print *, 'REF: , 1.1336318519308867e-10, -0.00011651612903225811, &
        2.330322580645162e-07'

  end subroutine compare_reduced_model

  subroutine redim_fields(var, var0)

    use type_module, ONLY:nx, ny, nz, nVar, nbx, nby, nbz, include_ice2, &
        nVarIce, inN, inQ, inQv, master, model, fluctuationMode, topography, &
        thetaRefRatio, timeScheme, L_ice, R_v, S_c, Dep, pSatIceRef, epsil0
    use mpi_module
    use atmosphere_module, ONLY:PStrat01, PStratTilde01, PStrat, rhoStrat, &
        piStrat, kappaInv, PStratTFC, piStratTfc, rhoStratTFC, gamma_1, p0, &
        mRef, uRef, kappa, Rsp, thetaRef, tRef, g, pRef
    use boundary_module !, ONLY : setHalos, setBoundary, reconstruction
    use flux_module
    !   use update_module, ONLY : dIce, ice2Update

    implicit none

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(in) :: var
    !    real, dimension(-1:nx,-1:ny,-1:nz, nVarIce), intent(out) :: source
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(out) :: var0
    real :: rho, pres, temp, theta, psi, Qv, SIce, NIce, exn_p, A
    integer :: i, j, k

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

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
                exn_p = var(i, j, k, 5) + (PstratTFC(i, j, k) / p0) ** gamma_1 !end if

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

            psi = Psat_ice(temp)

            Qv = var(i, j, k, inQv) ! Q_v = \rho q_v

            NIce = var(i, j, k, inN) ! N_v = \rho n

            SIce = SatRatio(Qv, temp, psi)

          else

            print *, ' ice2Sources works only with model  &
                == pseudo_incompressible '
            stop

          end if ! pseudo_inc

          var0(i, j, k, inN) = NIce / rho / mRef

          !var0(i,j,k, inQ) = SIce

          var0(i, j, k, inQv) = Qv / rho

          !var0(i,j,k, inQ) = var(i, j, k, 4)*uRef !store w at 10

          !asymptotic solution
          A = var(i, j, k, 4) * uRef * L_ice * g * kappa / Rsp / thetaRef ** 2 &
              / R_v
          !var0(i,j,k, inQ) = 2*A*S_c / (S_c - 1.) / (Dep * temp * psi / pres /tRef) / mRef
          var0(i, j, k, inQ) = 2 * A * S_c / ((S_c - 1.) * Dep * temp * pRef &
              / pSatIceRef / epsil0 / tRef) / mRef

        end do
      end do
    end do

  end subroutine redim_fields

  subroutine integrate_ice2(var, var0, flux, fluxmode, source, dt, q, RKstage, &
      PS0, PSTilde0, update_type, uTime, qTime)

    !new routine to handle both advection and ice physics

    use type_module, ONLY:nx, ny, nz, nVar, nbx, nby, nbz, include_ice2, &
        master, nVarIce, heatingONK14, raytracer
    use mpi_module

    !use atmosphere_module, ONLY : PStrat01, PStratTilde01
    use boundary_module, ONLY:setBoundary
    use flux_module, ONLY:ice2Flux, reconstruction
    use update_module, ONLY:ice2Update_apb, timeUpdate

    implicit none

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var, source
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(in) :: var0
    real, dimension(- 1:nx, - 1:ny, - 1:nz, 3, nVar) :: flux
    real, dimension(- 1:nz + 2), intent(in) :: PS0, PSTilde0
    character(len = *), intent(in) :: fluxmode
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVarIce), &
        intent(inout) :: q ! RK update for ice fields
    real :: dt
    integer :: RKstage
    character(len = 3), intent(in) :: update_type
    real :: uTime, qTime, pTime

    if(update_type == 'ADV' .or. update_type == 'BOT') then

      call setHalos(var, "ice2")
      call setBoundary(var, flux, "ice2")

      call reconstruction(var, "ice2")

      !** ?? call setHalos( var, "iceTilde" )
      call setBoundary(var, flux, "iceTilde")

      ! Fluxes and Sources
      call ice2Flux(var0, var, flux, fluxmode, PS0, PSTilde0)

    end if

    if(update_type == 'PHY' .or. update_type == 'BOT') then

      if(heatingONK14) then
        print *, 'ice2sources does not work with heating'
        stop
      else
        if(raytracer) then
          pTime = uTime !save previous time
          call timeUpdate(uTime, dt, qTime, RKstage)
          call ice2Sources_raytr(var, source, uTime, pTime)
        else
          call ice2Sources(var, source)
        end if ! raytracer
      end if

    end if

    if(update_type == 'BOT') then
      call ice2Update_apb(var, flux, source, dt, q, RKstage, 'BOT')
    elseif(update_type == 'ADV') then
      call ice2Update_apb(var, flux, source, dt, q, RKstage, 'ADV')
    elseif(update_type == 'PHY') then
      call ice2Update_apb(var, flux, source, dt, q, RKstage, 'PHY')
    else
      print *, 'unknown update_type in integrate_ice2'
      stop
    end if

  end subroutine integrate_ice2

end module ice2_module
