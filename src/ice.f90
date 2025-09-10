module ice_module

  use mpi

  INTERFACE OPERATOR(*)
    MODULE PROCEDURE sxc
  END INTERFACE OPERATOR(*)

  INTERFACE OPERATOR(*)
    MODULE PROCEDURE sxsc
  END INTERFACE OPERATOR(*)

  INTERFACE OPERATOR(/)
    MODULE PROCEDURE scds
  END INTERFACE OPERATOR(/)

  !INTERFACE OPERATOR(+)
  !  MODULE PROCEDURE cpc
  !END INTERFACE OPERATOR(+)

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE scpsc
  END INTERFACE OPERATOR(+)

  INTERFACE OPERATOR(-)
    MODULE PROCEDURE cmc
  END INTERFACE OPERATOR(-)

  INTERFACE OPERATOR(-)
    MODULE PROCEDURE scmsc
  END INTERFACE OPERATOR(-)

  CONTAINS
  function sxc(sc, cl) result(scl)
    use type_module, ONLY:ice_rayType2
    implicit none
    real, INTENT(IN) :: sc

    type(ice_rayType2), dimension(:, :, :), INTENT(IN) :: cl
    type(ice_rayType2), dimension(size(cl, dim = 1), size(cl, dim = 2), &
        &size(cl, dim = 3)) :: scl

    scl%wwp = sc * cl%wwp
    scl%epp = sc * cl%epp
    scl%thp = sc * cl%thp
    scl%si = sc * cl%si

  end function sxc

  function sxsc(sc, cl) result(scl)
    use type_module, ONLY:ice_rayType2
    implicit none
    real, INTENT(IN) :: sc

    type(ice_rayType2), INTENT(IN) :: cl
    type(ice_rayType2) :: scl

    scl%wwp = sc * cl%wwp
    scl%epp = sc * cl%epp
    scl%thp = sc * cl%thp
    scl%si = sc * cl%si

  end function sxsc

  function scds(cl, sc) result(scl)
    use type_module, ONLY:ice_rayType2
    implicit none
    real, INTENT(IN) :: sc

    type(ice_rayType2), INTENT(IN) :: cl
    type(ice_rayType2) :: scl

    scl%wwp = cl%wwp / sc
    scl%epp = cl%epp / sc
    scl%thp = cl%thp / sc
    scl%si = cl%si / sc

  end function scds

  function scmsc(cl1, cl2) result(sum)
    use type_module, ONLY:ice_rayType2
    implicit none

    type(ice_rayType2), INTENT(IN) :: cl1, cl2
    type(ice_rayType2) :: sum

    sum%wwp = cl1%wwp - cl2%wwp
    sum%epp = cl1%epp - cl2%epp
    sum%thp = cl1%thp - cl2%thp
    sum%si = cl1%si - cl2%si

  end function scmsc

  function scpsc(cl1, cl2) result(sum)
    use type_module, ONLY:ice_rayType2
    implicit none

    type(ice_rayType2), INTENT(IN) :: cl1, cl2
    type(ice_rayType2) :: sum

    sum%wwp = cl1%wwp + cl2%wwp
    sum%epp = cl1%epp + cl2%epp
    sum%thp = cl1%thp + cl2%thp
    sum%si = cl1%si + cl2%si

  end function scpsc

  function cmc(cl1, cl2) result(sum)
    use type_module, ONLY:ice_rayType2
    implicit none

    type(ice_rayType2), dimension(:, :, :), INTENT(IN) :: cl1, cl2
    type(ice_rayType2), dimension(size(cl1, dim = 1), size(cl1, dim = 2), &
        &size(cl1, dim = 3)) :: sum

    sum%wwp = cl1%wwp - cl2%wwp
    sum%epp = cl1%epp - cl2%epp
    sum%thp = cl1%thp - cl2%thp
    sum%si = cl1%si - cl2%si

  end function cmc

  function cpc(cl1, cl2) result(sum)
    use type_module, ONLY:ice_rayType2
    implicit none

    type(ice_rayType2), dimension(:, :, :), INTENT(IN) :: cl1, cl2
    type(ice_rayType2), dimension(size(cl1, dim = 1), size(cl1, dim = 2), &
        &size(cl1, dim = 3)) :: sum

    sum%wwp = cl1%wwp + cl2%wwp
    sum%epp = cl1%epp + cl2%epp
    sum%thp = cl1%thp + cl2%thp
    sum%si = cl1%si + cl2%si

  end function cpc

  subroutine setup_ice(var)

    use type_module
    use atmosphere_module, ONLY:tRef, rhoRef, lRef, thetaRef, pRef, PStrat, &
        &rhoStrat, piStrat, kappaInv, PStratTFC, piStratTfc, rhoStratTFC, &
        &gamma_1, p0, g, Rsp

    implicit none

    type(var_type), intent(inout) :: var

    integer, parameter :: ic_ice = 1

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

      !define upper/lower bounds of ISSR
      zMin_issr = z0_issr - sig_issr
      zMax_issr = z0_issr + sig_issr

      ! define flat topography
      !if(topography) then
      !  topography_surface = 0.
      !  print *, 'TOPOGRAPHY SURFACE SET TO ZERO in ice.f90 !'
      !**stop
      !end if

      !end case(1)

    case default
      print *, 'invalid ic_ice selected'
      stop
    end select

    S0 = 0.
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
            if(zTFC(i, j, k) .lt. zMin_issr .or. zTFC(i, j, k) .gt. zMax_issr) &
                &cycle
            S0 = S_issr * exp(- (zTFC(i, j, k) - z0_issr) ** 2 / 2. / sig_issr &
                &** 2)
          else
            if(z(k) .lt. zMin_issr .or. z(k) .gt. zMax_issr) cycle
            S0 = S_issr * exp(- (z(k) - z0_issr) ** 2 / 2. / sig_issr ** 2)
          end if

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
    type(ice_rayType2), dimension(1 - nbxSC:nxNSCX + nbxSC, 1 - nbySC:nyNSCY &
        &+ nbySC, 1 - nbzSC:nzNSCZ + nbzSC), intent(inout), optional :: &
        &ray_cloud
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
    integer :: ii, jj, kk, ii2, jj2, kk2

    do k = 1, nz

      if(.not. topography .and. (z(k) .lt. zMin_issr .or. z(k) .gt. &
          &zMax_issr)) cycle

      do j = 1, ny
        do i = 1, nx

          if(topography .and. (zTFC(i, j, k) .lt. zMin_issr .or. zTFC(i, j, k) &
              &.gt. zMax_issr)) cycle

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
                  ii2 = (i - 1) * NSCX + ii
                  do jj = 1, NSCY
                    jj2 = (j - 1) * NSCY + jj
                    do kk = 1, NSCZ
                      kk2 = (k - 1) * NSCZ + kk

                      expPrime = ray_cloud(ii2, jj2, kk2)%epp
                      thetaPrime = ray_cloud(ii2, jj2, kk2)%thp
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

                      Qv = ray_cloud(ii2, jj2, kk2)%Qv
                      NIce = ray_cloud(ii2, jj2, kk2)%Ni
                      SIce = SatRatio_new(Qv, pres, psi, rhoMean)

                      if(SIce .ge. S_c) then !only if critical S reached
                        !CHANGES
                        SIce = S_c
                        ray_cloud(ii2, jj2, kk2)%tNi = nucleation_n(SIce, &
                            &rhoMean)

                      else
                        ray_cloud(ii2, jj2, kk2)%tNi = 0.
                      end if

                      ray_cloud(ii2, jj2, kk2)%tQv = deposition_qv(SIce, NIce, &
                          &temp, pres, psi)
                      ray_cloud(ii2, jj2, kk2)%tQi = - ray_cloud(ii2, jj2, &
                          &kk2)%tQv

                      if(include_testoutput) then
                        ray_cloud(ii2, jj2, kk2)%si = SIce !full SIce in RT
                      end if
                    end do ! kk
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

                !                if(include_testoutput) then
                !                  var%OPT(i, j, k, 3) = FORgw / tref
                !                end if

              end if ! S ge Sc

              source%ICE(i, j, k, inQv) = deposition_qv(SIce, NIce, temp, &
                  &pres, psi)
              source%ICE(i, j, k, inQ) = - source%ICE(i, j, k, inQv)

            elseif(.not. raytracer .or. (raytracer .and. (.not. &
                &compute_cloudcover))) then ! no parameterized_nucleation

              if(SIce .ge. S_c) then !only if critical S reached
                !CHANGES
                SIce = S_c
                !*source%ICE(i, j, k, inN) = nucleation_n(SIce, rho)
                source%ICE(i, j, k, inN) = nucleation_n(SIce, rhoMean)

              else
                source%ICE(i, j, k, inN) = 0.
              end if

              source%ICE(i, j, k, inQv) = deposition_qv(SIce, NIce, temp, &
                  &pres, psi)
              source%ICE(i, j, k, inQ) = - source%ICE(i, j, k, inQv)

              if(include_testoutput) then
                var%OPT(i, j, k, 1) = SIce !full SIce in LES
              end if
            end if
          end if

        end do ! i
      end do ! j
    end do ! k

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
    use atmosphere_module, ONLY:tRef, rhoRef, lRef, thetaRef, pRef

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

    type(ice_rayType2), dimension(1 - nbxSC:nxNSCX + nbxSC, 1 - nbySC:nyNSCY &
        &+ nbySC, 1 - nbzSC:nzNSCZ + nbzSC), intent(inout), optional :: &
        &ray_cloud

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
      if(compute_cloudcover) then
        call iceSources_raytr(var, source, ray_varIce, ray_cloud = ray_cloud)
      else
        call iceSources_raytr(var, source, ray_varIce)
      end if

    end if

    if(update_type == 'BOT') then

      if(compute_cloudcover) then
        call iceUpdate_apb(var, flux, source, dt, q, RKstage, 'BOT', ray_cloud &
            &= ray_cloud)
      else
        call iceUpdate_apb(var, flux, source, dt, q, RKstage, 'BOT')
      end if

    elseif(update_type == 'ADV') then

      if(compute_cloudcover) then
        call iceUpdate_apb(var, flux, source, dt, q, RKstage, 'ADV', ray_cloud &
            &= ray_cloud)
      else
        call iceUpdate_apb(var, flux, source, dt, q, RKstage, 'ADV')
      end if

    elseif(update_type == 'PHY') then

      if(compute_cloudcover) then
        call iceUpdate_apb(var, flux, source, dt, q, RKstage, 'PHY', ray_cloud &
            &= ray_cloud)
      else
        call iceUpdate_apb(var, flux, source, dt, q, RKstage, 'PHY')
      end if

    else
      print *, 'unknown update_type in integrate_ice'
      stop
    end if

  end subroutine integrate_ice

  subroutine smooth_shapiro_sgs(fc_shap, ray_cloud)

    !--------------------------------------------------------------------
    !    local smoothing of density, winds, pressure,
    !    and density fluctuations
    !    order of horizontal shapiro filter given by 2*n_shap
    !    order of vertical shapiro filter is 2
    !    0 <= fcshap <= 1 is fraction of shapiro filter applied
    !
    !    the implementation of the boundary conditions is sub-optimal:
    !    (1) all elements of var are processed
    !        (although pressure and - in the explicit case - the density
    !         fluctuations are not filtered)
    !    (2) flux only transferred to the subroutine because it is an
    !        argument of setBoundary
    !        (it is not used, however)
    !    this could be done more efficiently
    !-------------------------------------------------------------------

    use type_module
    use mpi_module
    use boundary_module, ONLY:setBoundary, setBoundary_x_periodic, &
        &setBoundary_y_periodic

    implicit none
    ! in/out variables
    type(ice_rayType2), dimension(1 - nbxSC:nxNSCX + nbxSC, 1 - nbySC:nyNSCY &
        &+ nbySC, 1 - nbzSC:nzNSCZ + nbzSC), intent(inout) :: ray_cloud
    real, intent(in) :: fc_shap
    !integer, intent(in) :: n_shap

    ! allocatable fields
    type(ice_rayType2), dimension(:, :, :), allocatable :: field, var_l, &
        &flxwkb, flxwkb_0, flxwkb_1

    integer :: allocstat
    integer :: i, j, k
    integer :: nsmth
    integer :: iVar, ivmax
    integer :: i_lapl
    integer :: nz_max
    integer, parameter :: sm_order = 4

    allocate(field(1 - nbxSC:nxNSCX + nbxSC, 1 - nbySC:nyNSCY + nbySC, 1 &
        &- nbzSC:nzNSCZ + nbzSC), stat = allocstat)
    if(allocstat /= 0) stop "smooth_shapiro:alloc failed"

    allocate(var_l(1 - nbxSC:nxNSCX + nbxSC, 1 - nbySC:nyNSCY + nbySC, 1 &
        &- nbzSC:nzNSCZ + nbzSC), stat = allocstat)
    if(allocstat /= 0) stop "smooth_shapiro:alloc failed"

    ! make sure that boundary conditions are satisfied

    !*call setHalos(var, "var")

    !CHANGE
    !call setBoundary_subgridcell(ray_cloud, "sgs")

    !if(min(nxNSCX, nbxSC) < sm_order) stop 'min(nxNSCX,nbxSC) too small for smoothing'
    !if(min(nyNSCY, nbySC) < sm_order) stop 'min(nyNSCY,nbySC) too small for smoothing'
    !if(min(nzNSCZ, nbzSC) < nsmth) stop 'min(nzNSCZ,nbzSC) too small for smoothing'

    if(sm_order == 1) then
      var_l = ray_cloud

      !smoothing in x-direction

      ! if(sizeX > 1) then
      !   ! 2n-th-order x-derivative of all fields that are to be smoothed
      !   field = var_l
      !     do i = 1, nxNSCX
      !       var_l(i:i, 1:nyNSCY, 1:nzNSCZ) = (1./4.) * (field(i - 1:i - 1, 1:nyNSCY, 1:nzNSCZ) + field(i + 1:i+1, &
      !           &1:nyNSCY, 1:nzNSCZ) + 2.0 * field(i:i, 1:nyNSCY, 1:nzNSCZ))
      !     end do

      ! end if

      ! !smoothing in y-direction

      ! if(sizeY > 1) then
      !   ! 2n-th-order y-derivative of all fields that are to be smoothed

      !   field = var_l
      !     do j = 1, nyNSCY
      !       var_l(1:nxNSCX, j:j, 1:nzNSCZ) = (1./ 4.0)*(field(1:nxNSCX, j-1:j-1, 1:nzNSCZ) + &
      !         field(1:nxNSCX, j+1:j+1, 1:nzNSCZ) + 2.0 * field(1:nxNSCX, j:j, 1:nzNSCZ))
      !     end do
      ! end if

      ! ! filtering in z-direction

      ! ! 2nd-order z-derivative of all fields that are to be smoothed

      ! field = var_l

      ! do k = 1, nz
      !   var_l(1:nxNSCX, 1:nyNSCY, k:k) = ( 1.0/ 4.0) * (field(1:nxNSCX, 1:nyNSCY, k-1 :k - 1) + field(1:nxNSCX, 1:nyNSCY, k+1:k &
      !       &+ 1) + 2.0 * field(1:nxNSCX, 1:nyNSCY, k:k))
      ! end do

      ! apply filter

      ray_cloud(1:nxNSCX, 1:nyNSCY, 1:nzNSCZ - 1) = var_l(1:nxNSCX, 1:nyNSCY, &
          &1:nzNSCZ - 1)

      deallocate(field, stat = allocstat); if(allocstat /= 0) stop &
          &"smooth_shapiro:dealloc failed"
      deallocate(var_l, stat = allocstat); if(allocstat /= 0) stop &
          &"smooth_shapiro:dealloc failed"

    elseif(sm_order == 4) then

      allocate(flxwkb(- nbxSC:nxNSCX + nbxSC, - nbySC:nyNSCY + nbySC, &
          &- nbzSC:nzNSCZ + nbzSC), stat = allocstat)
      allocate(flxwkb_0(- nbxSC:nxNSCX + nbxSC, - nbySC:nyNSCY + nbySC, &
          &- nbzSC:nzNSCZ + nbzSC), stat = allocstat)
      allocate(flxwkb_1(- nbxSC:nxNSCX + nbxSC, - nbySC:nyNSCY + nbySC, &
          &- nbzSC:nzNSCZ + nbzSC), stat = allocstat)
      if(allocstat /= 0) stop "smooth_shapiro:alloc flxwkb_X failed"

      flxwkb(1:nxNSCX, 1:nyNSCY, 1:nzNSCZ) = ray_cloud(1:nxNSCX, 1:nyNSCY, &
          &1:nzNSCZ)

      !CHANGE
      call setBoundary_subgridcell(flxwkb(1 - nbxSC:nxNSCX + nbxSC, 1 &
          &- nbySC:nyNSCY + nbySC, 1 - nbzSC:nzNSCZ + nbzSC), "sgs")

      !check if nbX, nbY, nbZ are consistent with size of the field
      if(min(nxNSCX, nbxSC) >= sm_order) then
        do k = 1 - nbzSC, nzNSCZ + nbzSC
          do j = 1 - nbySC, nyNSCY + nbySC
            do i = 1, nxNSCX

              flxwkb_0(i, j, k)%wwp = ((- 1.0) * flxwkb(i - 4, j, k)%wwp &
                  &- flxwkb(i + 4, j, k)%wwp + 8.0 * (flxwkb(i - 3, j, k)%wwp &
                  &+ flxwkb(i + 3, j, k)%wwp) - 28.0 * (flxwkb(i - 2, j, &
                  &k)%wwp + flxwkb(i + 2, j, k)%wwp) + 56.0 * (flxwkb(i - 1, &
                  &j, k)%wwp + flxwkb(i + 1, j, k)%wwp) + 186.0 * flxwkb(i, j, &
                  &k)%wwp) / 256.0

              flxwkb_0(i, j, k)%epp = ((- 1.0) * flxwkb(i - 4, j, k)%epp &
                  &- flxwkb(i + 4, j, k)%epp + 8.0 * (flxwkb(i - 3, j, k)%epp &
                  &+ flxwkb(i + 3, j, k)%epp) - 28.0 * (flxwkb(i - 2, j, &
                  &k)%epp + flxwkb(i + 2, j, k)%epp) + 56.0 * (flxwkb(i - 1, &
                  &j, k)%epp + flxwkb(i + 1, j, k)%epp) + 186.0 * flxwkb(i, j, &
                  &k)%epp) / 256.0

              flxwkb_0(i, j, k)%thp = ((- 1.0) * flxwkb(i - 4, j, k)%thp &
                  &- flxwkb(i + 4, j, k)%thp + 8.0 * (flxwkb(i - 3, j, k)%thp &
                  &+ flxwkb(i + 3, j, k)%thp) - 28.0 * (flxwkb(i - 2, j, &
                  &k)%thp + flxwkb(i + 2, j, k)%thp) + 56.0 * (flxwkb(i - 1, &
                  &j, k)%thp + flxwkb(i + 1, j, k)%thp) + 186.0 * flxwkb(i, j, &
                  &k)%thp) / 256.0
            end do
          end do
        end do
      end if

      ! ! smooth in y
      ! if (min(nyNSCY, nbySC) >= sm_order ) then
      !   do k = 1 - nbzSC, nzNSCZ + nbzSC
      !     do j = 1, nyNSCY
      !       do i = 1, nxNSCX
      !         flxwkb_1(i, j, k) = ( (-1.0) * flxwkb_0(i, j - 4, k) - flxwkb_0(i, j &
      !             &+ 4, k) + 8.0 * (flxwkb_0(i, j - 3, k) + flxwkb_0(i, j + 3, &
      !             &k)) - 28.0 * (flxwkb_0(i, j - 2, k) + flxwkb_0(i, j + 2, &
      !             &k)) + 56.0 * (flxwkb_0(i, j - 1, k) + flxwkb_0(i, j + 1, &
      !             &k)) + 186.0 * flxwkb_0(i, j, k)) / 256.0
      !       end do
      !     end do
      !   end do
      ! endif

      flxwkb_1 = flxwkb_0

      ! smooth in z
      !if (min(nzNSCZ, nbzSC) >= sm_order ) then
      ! do k = 1, nzNSCZ
      !   do j = 1, nyNSCY
      !     do i = 1, nxNSCX
      !       flxwkb(i, j, k) = ((-1.0) *flxwkb_1(i, j, k - 4) - flxwkb_1(i, j, k + 4) &
      !           &+ 8.0 * (flxwkb_1(i, j, k - 3) + flxwkb_1(i, j, k + 3)) &
      !           &- 28.0 * (flxwkb_1(i, j, k - 2) + flxwkb_1(i, j, k + 2)) &
      !           &+ 56.0 * (flxwkb_1(i, j, k - 1) + flxwkb_1(i, j, k + 1)) &
      !           &+ 186.0 * flxwkb_1(i, j, k)) / 256.0
      !     end do
      !   end do
      ! end do
      !endif

      do k = 1, nzNSCZ
        do j = 1, nyNSCY
          do i = 1, nxNSCX

            flxwkb(i, j, k)%wwp = ((- 1.0) * flxwkb_1(i, j, k - 4)%wwp &
                &- flxwkb_1(i, j, k + 4)%wwp + 8.0 * (flxwkb_1(i, j, k &
                &- 3)%wwp + flxwkb_1(i, j, k + 3)%wwp) - 28.0 * (flxwkb_1(i, &
                &j, k - 2)%wwp + flxwkb_1(i, j, k + 2)%wwp) + 56.0 &
                &* (flxwkb_1(i, j, k - 1)%wwp + flxwkb_1(i, j, k + 1)%wwp) &
                &+ 186.0 * flxwkb_1(i, j, k)%wwp) / 256.0

            flxwkb(i, j, k)%epp = ((- 1.0) * flxwkb_1(i, j, k - 4)%epp &
                &- flxwkb_1(i, j, k + 4)%epp + 8.0 * (flxwkb_1(i, j, k &
                &- 3)%epp + flxwkb_1(i, j, k + 3)%epp) - 28.0 * (flxwkb_1(i, &
                &j, k - 2)%epp + flxwkb_1(i, j, k + 2)%epp) + 56.0 &
                &* (flxwkb_1(i, j, k - 1)%epp + flxwkb_1(i, j, k + 1)%epp) &
                &+ 186.0 * flxwkb_1(i, j, k)%epp) / 256.0

            flxwkb(i, j, k)%thp = ((- 1.0) * flxwkb_1(i, j, k - 4)%thp &
                &- flxwkb_1(i, j, k + 4)%thp + 8.0 * (flxwkb_1(i, j, k &
                &- 3)%thp + flxwkb_1(i, j, k + 3)%thp) - 28.0 * (flxwkb_1(i, &
                &j, k - 2)%thp + flxwkb_1(i, j, k + 2)%thp) + 56.0 &
                &* (flxwkb_1(i, j, k - 1)%thp + flxwkb_1(i, j, k + 1)%thp) &
                &+ 186.0 * flxwkb_1(i, j, k)%thp) / 256.0

            flxwkb(i, j, k)%si = ((- 1.0) * flxwkb_1(i, j, k - 4)%si &
                &- flxwkb_1(i, j, k + 4)%si + 8.0 * (flxwkb_1(i, j, k - 3)%si &
                &+ flxwkb_1(i, j, k + 3)%si) - 28.0 * (flxwkb_1(i, j, k &
                &- 2)%si + flxwkb_1(i, j, k + 2)%si) + 56.0 * (flxwkb_1(i, j, &
                &k - 1)%si + flxwkb_1(i, j, k + 1)%si) + 186.0 * flxwkb_1(i, &
                &j, k)%si) / 256.0
          end do
        end do
      end do

      ! print*, maxval(abs(ray_cloud%si))
      ! print*, maxval(abs(flxwkb%si))
      ! print*, ' *** '

      ! do k = 1, nzNSCZ
      !   do j = 1, nyNSCY
      !     do i = 1, nxNSCX

      !       flxwkb(i, j, k) = flxwkb(i, j, k)+flxwkb(i, j, k)
      !     end do
      !   end do
      ! end do

      ray_cloud(1:nxNSCX, 1:nyNSCY, 1:nzNSCZ) = flxwkb(1:nxNSCX, 1:nyNSCY, &
          &1:nzNSCZ)

      deallocate(flxwkb, stat = allocstat); if(allocstat /= 0) stop &
          &"smooth_shapiro:dealloc failed"
      deallocate(flxwkb_0, stat = allocstat); if(allocstat /= 0) stop &
          &"smooth_shapiro:dealloc failed"
      deallocate(flxwkb_1, stat = allocstat); if(allocstat /= 0) stop &
          &"smooth_shapiro:dealloc failed"
    else
      print *, 'smooth_shapiro_sgs: unknown smoothing order'
      stop
    end if
    ! boundary conditions again

    !*call setHalos(var, "var")

    !CHANGE
    !**call setBoundary_subgridcell(ray_cloud, "sgs")

    ! deallocate local fields

    return
  end subroutine smooth_shapiro_sgs

  subroutine setBoundary_x_periodic_subgridcell(ray_cloud, option)
    !-------------------------------------------
    ! 1) set values in ghost cells according to
    !    boundary conditions
    ! 2) set fluxes explicitly to zero at walls
    !-------------------------------------------
    use type_module
    implicit none
    ! in/out variables
    type(ice_rayType2), dimension(1 - nbxSC:nxNSCX + nbxSC, 1 - nbySC:nyNSCY &
        &+ nbySC, 1 - nbzSC:nzNSCZ + nbzSC), intent(inout) :: ray_cloud
    character(len = *), intent(in) :: option

    ! local variables
    integer :: i, j, k, iVar, ii

    select case(option)

    case("sgs")

      do i = 1, nbxsc

        ray_cloud(nxNSCX + i, :, :) = ray_cloud(i, :, :)
        ray_cloud(- i + 1, :, :) = ray_cloud(nxNSCX - i + 1, :, :)

      end do
    case default
      stop "setBoundary_x_subgridcell: unknown option."
    end select

  end subroutine setBoundary_x_periodic_subgridcell

  subroutine setBoundary_y_periodic_subgridcell(ray_cloud, option)
    !-------------------------------------------
    ! 1) set values in ghost cells according to
    !    boundary conditions
    ! 2) set fluxes explicitly to zero at walls
    !-------------------------------------------
    use type_module
    implicit none
    ! in/out variables
    type(ice_rayType2), dimension(1 - nbxSC:nxNSCX + nbxSC, 1 - nbySC:nyNSCY &
        &+ nbySC, 1 - nbzSC:nzNSCZ + nbzSC), intent(inout) :: ray_cloud
    character(len = *), intent(in) :: option

    ! local variables
    integer :: i, j, k, iVar, ii

    if(sizeY > 1) then
      select case(option)
      case("sgs")
        do j = 1, nbysc
          ray_cloud(:, nyNSCY + j, :) = ray_cloud(:, j, :)
          ray_cloud(:, - j + 1, :) = ray_cloud(:, nyNSCY - j + 1, :)
        end do
      case default
        stop "setBou    ndary_y: unknown option."
      end select
    end if

  end subroutine setBoundary_y_periodic_subgridcell

  subroutine setBoundary_subgridcell(ray_cloud, option)
    !-------------------------------------------
    ! 1) set values in ghost cells according to
    !    boundary conditions
    ! 2) set fluxes explicitly to zero at walls
    !-------------------------------------------

    use type_module
    implicit none
    ! in/out variables
    type(ice_rayType2), dimension(1 - nbxSC:nxNSCX + nbxSC, 1 - nbySC:nyNSCY &
        &+ nbySC, 1 - nbzSC:nzNSCZ + nbzSC), intent(inout) :: ray_cloud
    character(len = *), intent(in) :: option

    !------------------------------
    !          x-direction
    !------------------------------
    select case(xBoundary)

    case("periodic")

      if(idim > 1) then
        ! boundary conditions taken care of by setHalos
      else
        call setBoundary_x_periodic_subgridcell(ray_cloud, option)
      endif

    case default
      stop "setBoundary: unknown case xBoundary_subgridcell"
    end select

    !------------------------------
    !          y-direction
    !------------------------------
    select case(yBoundary)

    case("periodic")

      if(jdim > 1) then
        ! boundary conditions taken care of by setHalos
      else
        call setBoundary_y_periodic_subgridcell(ray_cloud, option)
      endif

    case default
      stop "setBoundary: unknown case yBoundary_subgridcell"
    end select

    !------------------------------
    !          z-direction
    !------------------------------
    select case(zBoundary)

    case("solid_wall")
      call setBoundary_z_solidWall_subgridscale(ray_cloud, option)

    case default
      stop "setBoundary: unknown case zBoundary_subgridcell"
    end select

  end subroutine setBoundary_subgridcell

  subroutine setBoundary_z_solidWall_subgridscale(ray_cloud, option)
    !-------------------------------------------
    ! 1) set values in ghost cells according to
    !    boundary conditions
    ! 2) set fluxes explicitly to zero at walls
    !-------------------------------------------
    use type_module
    ! in/out variables
    type(ice_rayType2), dimension(1 - nbxSC:nxNSCX + nbxSC, 1 - nbySC:nyNSCY &
        &+ nbySC, 1 - nbzSC:nzNSCZ + nbzSC), intent(inout) :: ray_cloud
    character(len = *), intent(in) :: option

    selectcase(option)

    case("sgs")
      ! reflect at boundary with change of sign
      ! at boundary var = 0
      do k = 1, nbzsc
        ray_cloud(:, :, nzNSCZ + k) = ray_cloud(:, :, k)
        ray_cloud(:, :, - k + 1:- k + 1) = (- 1.) * ray_cloud(:, :, nzNSCZ - k &
            &+ 1:nzNSCZ - k + 1)
      end do

    case default
      stop "setBoundary_z_subgridscale: unknown option."
    end select

  end subroutine setBoundary_z_solidWall_subgridscale

end module ice_module
