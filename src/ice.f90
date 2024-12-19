module ice_module

  use mpi

  contains

  subroutine setup_ice(var)

    use type_module
    use atmosphere_module, ONLY:heightTFC, tRef, rhoRef, lRef, thetaRef, pRef, &
        &PStrat, rhoStrat, piStrat, kappaInv, PStratTFC, piStratTfc, &
        &rhoStratTFC, gamma_1, p0, g, Rsp

    implicit none

    type(var_type), intent(inout) :: var

    integer, parameter :: ic_ice = 5

    real :: dz_tr
    real :: rho, exn_p, pres, temp, theta, psi

    integer :: iVar, ii, k, hzn, dzn, i, j
    real :: n0, q0, qv0, S0
    integer :: irec

    !case 3
    real :: z0_issr, sig_issr, S_issr

    !case 5
    real :: presMean, exn_pMean, tempMean, thetaMean, psiMean, rhoMean

    !set variables
    J_nuc = 4.9E4 !nucleation rate
    B_nuc = 337. !nucleation exponent
    Dep0 = 4.3E-8 !C_0 coeff paper
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

    Dep = Dep0 * thetaRef * tRef * meanMassIce ** (1. / 3.) * PsatIceRef &
        &/ pRef / mRef

    DepS = Dep0 * thetaRef * tRef * meanMassIce ** (1. / 3.) / mRef

    !call compare_reduced_model
    !if ( master ) stop

    ! init to zero
    do iVar = 1, nVarIce
      var%ICE(:, :, :, iVar) = 0.
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

              rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)

              theta = PstratTFC(i, j, k) / rho

              !problems in \pi
              !if heating switched on
              !if ( timeScheme == "semiimplicit" ) then
              !   print*, 'iceSources works only with explicit integration'
              !   stop
              !else
              exn_p = var%pi(i, j, k) + (PstratTFC(i, j, k) / p0) ** gamma_1
              !end if

            else

              rho = var%rho(i, j, k) + rhoStrat(k)

              !problems in \pi
              !if heating switched on
              theta = Pstrat(k) / rho

              exn_p = var%pi(i, j, k) + (PStrat(k) / p0) ** gamma_1

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

            var%ICE(:, :, k, inN) = rho * n0 * mRef !\hat N = \hat \rho \hat n
            var%ICE(:, :, k, inQv) = rho * qv0
            var%ICE(:, :, k, inQ) = rho * q0

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

              rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)
              ! N = \rho n
              var%ICE(:, :, k, inN) = dz_tr * heightTFC(i, j, k) * rho

            end do
          end do
        end do

      else

        do k = 1, nz

          ! N = \rho n
          !var(:,:,k, inN) = (k-1)*0.1
          var%ICE(:, :, k, inN) = dz_tr * (z(k) - z(1)) * (var%rho(:, :, k) &
              &+ rhoStrat(k))

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

              rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)

              theta = PstratTFC(i, j, k) / rho

              !problems in \pi
              !if heating switched on
              !if ( timeScheme == "semiimplicit" ) then
              !   print*, 'iceSources works only with explicit integration'
              !   stop
              !else
              exn_p = var%pi(i, j, k) + (PstratTFC(i, j, k) / p0) ** gamma_1
              !end if

            else

              rho = var%rho(i, j, k) + rhoStrat(k)

              !problems in \pi
              !if heating switched on
              theta = Pstrat(k) / rho

              exn_p = var%pi(i, j, k) + (PStrat(k) / p0) ** gamma_1

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

            var%ICE(i, j, k, inN) = rho * n0 * mRef !\hat N = \hat \rho \hat n
            var%ICE(i, j, k, inQv) = rho * qv0
            var%ICE(i, j, k, inQ) = rho * q0

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

              rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)

              theta = PstratTFC(i, j, k) / rho

              !problems in \pi
              !if heating switched on
              !if ( timeScheme == "semiimplicit" ) then
              !   print*, 'iceSources works only with explicit integration'
              !   stop
              !else
              exn_p = var%pi(i, j, k) + (PstratTFC(i, j, k) / p0) ** gamma_1
              !end if

            else

              rho = var%rho(i, j, k) + rhoStrat(k)

              !problems in \pi
              !if heating switched on
              theta = Pstrat(k) / rho

              exn_p = var%pi(i, j, k) + (PStrat(k) / p0) ** gamma_1

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

            var%ICE(i, j, k, inN) = rho * n0 * mRef !\hat N = \hat \rho \hat n
            var%ICE(i, j, k, inQv) = rho * qv0
            var%ICE(i, j, k, inQ) = rho * q0

            if(include_testoutput) then

              var%OPT(i, j, k, 1) = S0
              var%OPT(i, j, k, 2) = pres * pRef
              var%OPT(i, j, k, 3) = pres / psi * pRef / PsatIceRef

              if(raytracer) then
                var%OPT(i, j, k, 1) = S0
                var%OPT(i, j, k, 2) = 0.
                var%OPT(i, j, k, 3) = 0.
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

      ! define flat topography
      if(topography) then
        topography_surface = 0.
      end if

      do k = 1, nz

        do j = 1, ny
          do i = 1, nx

            if(topography) then

              rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)

              theta = PstratTFC(i, j, k) / rho

              rhoMean = rhoStratTFC(i, j, k)

              thetaMean = PstratTFC(i, j, k) / rhoMean

              exn_p = var%pi(i, j, k) + (PstratTFC(i, j, k) / p0) ** gamma_1

              exn_pMean = (PstratTFC(i, j, k) / p0) ** gamma_1

            else

              rho = var%rho(i, j, k) + rhoStrat(k)
              rhoMean = rhoStrat(k)

              theta = Pstrat(k) / rho

              thetaMean = Pstrat(k) / rhoStrat(k)

              exn_p = var%pi(i, j, k) + (PStrat(k) / p0) ** gamma_1

              exn_pMean = (PStrat(k) / p0) ** gamma_1

            end if ! topography

            pres = p0 * exn_p ** kappaInv !kappaInv = c_p/R

            presMean = p0 * exn_pMean ** kappaInv !kappaInv = c_p/R

            temp = theta * exn_p

            tempMean = thetaMean * exn_pMean

            psi = Psat_ice(temp)

            psiMean = Psat_ice(tempMean)

            !dimensional IC for n, q_v, q
            n0 = 0. !0.1 * 2.E6 ![kg**-1]

            !use S0 to define q_v
            !to be consistent with RayTracer simulation
            !NB: S0 is not initial S_i
            !S0  = q_v(0)* \mean p/ \mean p_si
            !S_i = q_v(0)* p/p_si
            if(topography) then
              S0 = S_issr * exp(- (heightTFC(i, j, k) - z0_issr) ** 2 / 2. &
                  &/ sig_issr ** 2)
            else
              S0 = S_issr * exp(- (z(k) - z0_issr) ** 2 / 2. / sig_issr ** 2)
            end if

            !qv0 = epsil0hat * S0 * psi / pres ! [kg/kg]
            !CHANGES
            qv0 = epsil0hat * S0 * psiMean / presMean ! [kg/kg]

            q0 = meanMassIce * n0

            var%ICE(i, j, k, inN) = rhoMean * n0 * mRef !\hat N = \hat \rho \hat n
            var%ICE(i, j, k, inQv) = rhoMean * qv0
            var%ICE(i, j, k, inQ) = rhoMean * q0

            if(include_testoutput) then

              !var%OPT(i, j, k, 1) = S0
              !var%OPT(i, j, k, 2) = qv0 / epsil0hat * pres / psi ! store IC S

              if(raytracer) then
                !var(i, j, k, iVarO) = S0
                !var%OPT(i, j, k, 2) = 0.
                !var%OPT(i, j, k, 3) = 0.
              end if
            end if

          end do !i
        end do !j

      end do !k
      !end case(5)

    end select

  end subroutine setup_ice

  subroutine iceSources_raytr(var, source, ray_varIce, ray_cloud)

    use type_module

    use atmosphere_module, ONLY:tRef !TEST
    use mpi_module
    use atmosphere_module, ONLY:PStrat01, PStratTilde01, PStrat, rhoStrat, &
        &piStrat, kappaInv, PStratTFC, piStratTfc, rhoStratTFC, gamma_1, p0, &
        &g, Rsp, uRef, pRef
    use flux_module

    implicit none

    type(var_type), intent(inout) :: var
    !    real, dimension(-1:nx,-1:ny,-1:nz, nVarIce), intent(out) :: source
    type(var_type), intent(inout) :: source
    type(ice_rayType), dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(in) :: &
        &ray_varIce
    type(ice_rayType2), dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, nscx, nscy), &
        &intent(inout) :: ray_cloud
    real :: rho, pres, temp, theta, psi, Qv, SIce, NIce, exn_p
    real :: amp_pi, w_gw, PiMean, PiPrime, delta_t !RayTracer coupling
    !real :: uTime, pTime
    integer :: i, j, k
    real :: dphi, omg
    real :: thetaPrime, expPrime, tPrime, pPrime, wPrime, wPrimeInt
    real :: dotThetaprime, dotExpprime, dotTPrime, dotpPrime, dotPiPrime, FORgw
    real :: N_pre, N_post
    real :: rhoMean
    real :: pres_ls, temp_ls
    integer :: ii, jj

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

          if(model == "pseudo_incompressible") then

            if(topography) then

              rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)
              rhoMean = rhoStratTFC(i, j, k)

              theta = PstratTFC(i, j, k) / rho

              !problems in \pi if heating of background
              exn_p = var%pi(i, j, k) + (PstratTFC(i, j, k) / p0) ** gamma_1

            else
              rho = var%rho(i, j, k) + rhoStrat(k)
              rhoMean = rhoStrat(k)

              theta = Pstrat(k) / rho

              !problems in \pi if heating of background
              exn_p = var%pi(i, j, k) + (PStrat(k) / p0) ** gamma_1

            end if ! topography

            pres = p0 * exn_p ** kappaInv !kappaInv = c_p/R

            temp = theta * exn_p

            psi = Psat_ice(temp)

            Qv = var%ICE(i, j, k, inQv) ! Q_v = \rho q_v

            NIce = var%ICE(i, j, k, inN) ! N_v = \rho n

            if(raytracer) then

              if(compute_cloudcover) then

                !store large scale values
                pres_ls = pres
                temp_ls = temp

                do ii = 1, NSCX
                  do jj = 1, NSCY

                    expPrime = ray_cloud(i, j, k, ii, jj)%epp
                    thetaPrime = ray_cloud(i, j, k, ii, jj)%thp
                    !compute p', T'
                    !pPrime = PStrat(k) * expPrime
                    !HOWEVER to be consistent with other implementation here we use
                    pPrime = p0 * (exn_p + expPrime) ** kappaInv - pres_ls
                    tPrime = thetaPrime * exn_p + expPrime * theta &
                        &+ thetaPrime * expPrime
                    !add GW fluctuations to large-scale fields
                    pres = pres_ls + pPrime
                    temp = temp_ls + tPrime

                    psi = Psat_ice(temp)

                    Qv = ray_cloud(i, j, k, ii, jj)%Qv
                    NIce = ray_cloud(i, j, k, ii, jj)%Ni
                    SIce = SatRatio_new(Qv, pres, psi, rhoMean)

                    if(SIce .ge. S_c) then !only if critical S reached
                      ray_cloud(i, j, k, ii, jj)%tNi = nucleation_n(SIce, &
                          &rhoMean)
                    else
                      ray_cloud(i, j, k, ii, jj)%tNi = 0.
                    end if
                    !CHANGES
                    SIce = S_c

                    ray_cloud(i, j, k, ii, jj)%tQv = deposition_qv(SIce, NIce, &
                        &temp, pres, psi)
                    ray_cloud(i, j, k, ii, jj)%tQi = - ray_cloud(i, j, k, ii, &
                        &jj)%tQv

                  end do !jj
                end do !ii
              else ! /= compute_cloudcover

                ! modify saturation ratio S
                ! to account for GW fluctuations in p/p_si
                !
                ! NB: we assume constant GW vertical velocity during
                ! each pincflow time step, this is  critical
                ! for high-frequency GWs and large integration time steps
                ! alternatively: one should position computation of Pi' in
                ! calc_ice to icesources_raytr and advance phase
                !
                !Pi = p/p_si
                !Pi = <Pi> + Pi'
                !<Pi> large scale field, Pi' GW fluctuations

                expPrime = ray_varIce(i, j, k)%epp
                thetaPrime = ray_varIce(i, j, k)%thp
                !compute p', T'
                !pPrime = PStrat(k) * expPrime
                !HOWEVER to be consistent with other implementation here we use
                pPrime = p0 * (exn_p + expPrime) ** kappaInv - pres
                tPrime = thetaPrime * exn_p + expPrime * theta + thetaPrime &
                    &* expPrime
                !add GW fluctuations to large-scale fields
                pres = pres + pPrime
                temp = temp + tPrime
                psi = Psat_ice(temp)

                !Changes in rho due to GW not accounted for

                SIce = SatRatio_new(Qv, pres, psi, rhoMean)
                !CHANGES
                !*SIce = SatRatio(Qv, temp, psi)

              end if ! cloudcover
            else ! /= raytracer
              SIce = SatRatio_new(Qv, pres, psi, rhoMean)
            end if ! raytracer

            !output variables
            if(include_testoutput) then
              if(rayTracer .and. (.not. compute_cloudcover)) then
                var%OPT(i, j, k, 1) = SIce !full SIce in LES/ MeanSIce in RT
                var%OPT(i, j, k, 2) = ray_varIce(i, j, k)%wwp ! vertical vel.
                !var%OPT(i, j, k, 3) = 0. ! Full S
              else
                var%OPT(i, j, k, 1) = SIce !full SIce in LES/ MeanSIce in RT
                var%OPT(i, j, k, 2) = 0.
                var%OPT(i, j, k, 3) = 0. ! Full S
              end if
            end if ! include_testoutput

          else

            print *, ' iceSources works only with model  &
                &== pseudo_incompressible '
            stop

          end if ! pseudo_inc

          if(SIce .ge. S_c) then
            print *, 'S HITS S_c'
            print *, '**********'
          end if

          if(.not. no_ice_source) then
            if(parameterized_nucleation .and. raytracer .and. (.not. &
                &compute_cloudcover)) then

              if(SIce .ge. S_c) then !check if nucleation condition fulfilled
                wPrime = ray_varIce(i, j, k)%wwp !vertical vel.

                !new GW forcing term: expressed terms of Exner pressure tendency only and
                !include advection of background pressure by GW vertical vel.

                !***dotPiprime = pres/psi * pRef/PsatIceRef * ( kappaInv - Li_hat / temp ) &
                !     * wPrime * gamma_1 * (PStrat(k+1) - PStrat(k-1))/2./dz / PStrat(k) !ONLY if no topography
                !add extra  /exn_p*( dotExpPrime !+ wPrime * ( piStrat(k+1)-piStrat(k-1))/2./dz )

                dotPiprime = (kappaInv - Li_hat / temp) * wPrime * gamma_1 &
                    &* (PStrat(k + 1) - PStrat(k - 1)) / 2. / dz / PStrat(k) !ONLY if no topography

                FORgw = dotPiPrime * S_c

                n_pre = NIce / rho
                n_post = nIce_param_nuc(n_pre, SIce, FORgw, temp, pres, psi)

                var%ICE(i, j, k, inN) = n_post * rho ! N_ice = \rho n
                source%ICE(i, j, k, inN) = 0. ! \dot N_ice=0.

                if(include_testoutput) then
                  var%OPT(i, j, k, 3) = FORgw / tref
                end if

              end if ! S ge Sc

              source%ICE(i, j, k, inQv) = deposition_qv(SIce, NIce, temp, &
                  &pres, psi)
              source%ICE(i, j, k, inQ) = - source%ICE(i, j, k, inQv)

            elseif(.not. raytracer .or. (raytracer .and. (.not. &
                &compute_cloudcover))) then ! no parameterized_nucleation

              if(SIce .ge. S_c) then !only if critical S reached

                !CHANGES
                !*source%ICE(i, j, k, inN) = nucleation_n(SIce, rho)
                source%ICE(i, j, k, inN) = nucleation_n(SIce, rhoMean)
              else
                source%ICE(i, j, k, inN) = 0.
              end if
              !CHANGES
              SIce = S_c

              source%ICE(i, j, k, inQv) = deposition_qv(SIce, NIce, temp, &
                  &pres, psi)
              source%ICE(i, j, k, inQ) = - source%ICE(i, j, k, inQv)
            end if
          end if

        end do ! i
      end do ! j
    end do ! k

    !**************************
    !TEST set tendency to zero
    !**************************
    !!$    source%ICE = 0
    !!$
    !!$    do i = -nbx, nx
    !!$       var%ICE(i, :, :, :) = sin(2*PI*(i + (is + nbx - 1) -1)/SizeX)
    !!$    end do
    !!$
    !**************************

  end subroutine iceSources_raytr

  function NIce_param_nuc(N, S, dotS, T, p, p_si)
    ! compute nucleated number of ice crystals
    ! using the asymptotic solution

    use type_module, ONLY:S_c, DepS
    implicit none

    real, intent(in) :: dotS, S, T, N, p, p_si
    real :: NIce_param_nuc

    NIce_param_nuc = 2. * dotS / (DepS * (S_c - 1.) * T) - N

    if(NIce_param_nuc .gt. n) then
      NIce_param_nuc = NIce_param_nuc - N
    else
      NIce_param_nuc = N
    end if

  end function NIce_param_nuc

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

  function SatRatio_new(Qv, p, p_si, rhoMean)
    ! compute saturation ratio
    ! S= p * q_v / p_si / epsil0 = R T \rho q_v / p_si / epsil0 =  R T Qv / p_si / epsil0
    use type_module, ONLY:epsil0hat
    implicit none

    real, intent(in) :: Qv, p, p_si, rhoMean
    real :: SatRatio_new

    SatRatio_new = Qv / rhoMean * p / p_si / epsil0hat

  end function SatRatio_New

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

  subroutine integrate_ice(var, var0, flux, fluxmode, source, dt, q, RKstage, &
      &PS0, PSTilde0, update_type, uTime, qTime, ray_varIce, ray_cloud)

    !new routine to handle both advection and ice physics

    use mpi_module
    use type_module
    !use atmosphere_module, ONLY : PStrat01, PStratTilde01
    use boundary_module, ONLY:setBoundary
    use flux_module, ONLY:iceFlux, reconstruction
    use update_module, ONLY:iceUpdate_apb, timeUpdate

    implicit none
    type(var_type), intent(inout) :: var, source
    type(var_type), intent(in) :: var0
    type(flux_type) :: flux
    type(ice_rayType), dimension(0:nx + 1, 0:ny + 1, 0:nz + 1), intent(in) :: &
        &ray_varIce

    type(ice_rayType2), dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, nscx, nscy), &
        &intent(inout) :: ray_cloud

    real, dimension(- 1:nz + 2), intent(in) :: PS0, PSTilde0
    character(len = *), intent(in) :: fluxmode
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVarIce), &
        &intent(inout) :: q ! RK update for ice fields
    real :: dt
    integer :: RKstage
    character(len = 3), intent(in) :: update_type
    real :: uTime, qTime, pTime

    if(update_type == 'ADV' .or. update_type == 'BOT') then

      call setHalos(var, "ice")
      call setBoundary(var, flux, "ice")

      call reconstruction(var, "ice")

      !** ?? call setHalos( var, "iceTilde" )
      call setBoundary(var, flux, "iceTilde")

      ! Fluxes and Sources
      call iceFlux(var0, var, flux, fluxmode, PS0, PSTilde0)

    end if

    if(update_type == 'PHY' .or. update_type == 'BOT') then

      !pTime = uTime !save previous time
      !call timeUpdate(uTime, dt, qTime, RKstage)
      call iceSources_raytr(var, source, ray_varIce, ray_cloud)

    end if

    if(update_type == 'BOT') then
      call iceUpdate_apb(var, flux, source, dt, q, RKstage, 'BOT', ray_cloud)
    elseif(update_type == 'ADV') then
      call iceUpdate_apb(var, flux, source, dt, q, RKstage, 'ADV', ray_cloud)
    elseif(update_type == 'PHY') then
      call iceUpdate_apb(var, flux, source, dt, q, RKstage, 'PHY', ray_cloud)
    else
      print *, 'unknown update_type in integrate_ice'
      stop
    end if

  end subroutine integrate_ice

end module ice_module
