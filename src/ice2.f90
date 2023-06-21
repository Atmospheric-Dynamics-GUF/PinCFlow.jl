module ice2_module

  contains

  subroutine setup_ice2(var)

    use type_module
    use atmosphere_module, ONLY:heightTFC, tRef, rhoRef, lRef, thetaRef, pRef, &
        PStrat, rhoStrat, piStrat, kappaInv, PStratTFC, piStratTfc, &
        rhoStratTFC, gamma_1, p0

    implicit none

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var

    integer, parameter :: ic_ice = 3

    real :: dz_tr
    real :: rho, exn_p, pres, temp, theta, psi

    integer :: iVar, ii, k, hzn, dzn, i, j
    real :: n0, q0, qv0, S0
    integer :: irec

    !case 3
    real :: z0_issr, sig_issr, S_issr

    !set variables
    J_nuc = 4.9E4 !nucleation rate
    B_nuc = 337. !nucleation exponent
    Dep = 4.3E-8 !C_0 coeff paper
    epsil0hat = epsil0 * PsatIceRef / pRef

    !nondimensionalize variables/ fields

    thetaRefRatio = thetaRef / thetaRef_trp

    L_hat = L_ice / R_v / thetaRef_trp

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

            !save for output
            !!$                ofield(i, j, k, 1) = pres
            !!$                ofield(i, j, k, 2) = temp
            !!$                ofield(i, j, k, 3) = rho
            !!$                ofield(i, j, k, 4) = PstratTFC(i, j, k)
            !!$                ofield(i, j, k, 5) = exn_p
            !!$                ofield(i, j, k, 6) = psi

          end do
        end do

        !end if !k in range

      end do !k

      !output p, T, rho, theta, P, p_si

      !!$       irec = 1
      !!$       if (master) then
      !!$          do ii = 1, 6
      !!$
      !!$             do k = 1, nz
      !!$
      !!$                write(44,rec=irec) ((ofield(i, j, k, ii), i=1, nx), j = 1, ny)
      !!$                irec = irec + 1
      !!$
      !!$             end do
      !!$          end do
      !!$
      !!$          close(44)
      !!$       end if

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

            var(:, :, k, inN) = rho * n0 * mRef !\hat N = \hat \rho \hat n
            var(:, :, k, inQv) = rho * qv0
            var(:, :, k, inQ) = rho * q0

          end do
        end do

      end do !k

    end select

  end subroutine setup_ice2

  subroutine ice2Sources(var, source)

    use type_module, ONLY:nx, ny, nz, nVar, nbx, nby, nbz, include_ice2, &
        nVarIce, inN, inQ, inQv, master, model, fluctuationMode, topography, &
        thetaRefRatio, timeScheme
    use mpi_module
    use atmosphere_module, ONLY:PStrat01, PStratTilde01, PStrat, rhoStrat, &
        piStrat, kappaInv, PStratTFC, piStratTfc, rhoStratTFC, gamma_1, p0
    use boundary_module !, ONLY : setHalos, setBoundary, reconstruction
    use flux_module
    !   use update_module, ONLY : dIce, ice2Update

    implicit none

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(in) :: var
    !    real, dimension(-1:nx,-1:ny,-1:nz, nVarIce), intent(out) :: source
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(out) :: source
    real :: rho, pres, temp, theta, psi, Qv, SIce, NIce, exn_p
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

          else

            print *, ' ice2Sources works only with model &
                == pseudo_incompressible '
            stop

          end if ! pseudo_inc

          source(i, j, k, inN) = nucleation_n(SIce, rho)

          source(i, j, k, inQv) = deposition_qv(SIce, NIce, temp, pres, psi)

          source(i, j, k, inQ) = - source(i, j, k, inQv)

        end do
      end do
    end do

  end subroutine ice2Sources

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

    SatRatio = T * Qv / p_si / epsil0hat

  end function SatRatio

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
      PS0, PSTilde0)

    use type_module, ONLY:nx, ny, nz, nVar, nbx, nby, nbz, include_ice2, &
        master, nVarIce, heatingONK14
    use mpi_module

    !use atmosphere_module, ONLY : PStrat01, PStratTilde01
    use boundary_module, ONLY:setBoundary
    use flux_module, ONLY:ice2Flux, reconstruction
    use update_module, ONLY:ice2Update, ice2Update_source

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

    call setHalos(var, "ice2")
    call setBoundary(var, flux, "ice2")

    call reconstruction(var, "ice2")

    !** ?? call setHalos( var, "iceTilde" )
    call setBoundary(var, flux, "iceTilde")

    !print*, 'maxval flux 1', maxval(abs(flux(1:nx,1:ny,1:nz,3,9:10)))

    ! Fluxes and Sources
    call ice2Flux(var0, var, flux, fluxmode, PS0, PSTilde0)

    !print*, 'maxval flux 2', maxval(abs(flux(1:nx,1:ny,1:nz,3,9:10)))

    !ONLY works with Pstat, PstratTFC: no heating
    if(heatingONK14) then
      print *, 'ice2sources does not work with heating'
      stop
    else
      call ice2Sources(var, source)
    end if

    !call setBoundary (var, flux, "flux") !do nothing for "flux"?

    ! OLD: NO SOURCE
    !call ice2Update(var, flux, dt, q, RKstage, &
    !     "expl", 1.)

    call ice2Update_source(var, flux, source, dt, q, RKstage)

  end subroutine integrate_ice

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

            print *, ' ice2Sources works only with model &
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

end module ice2_module

!!$    subroutine test_print_out
!!$
!!$      !SD
!!$      initialCleaning = .false.
!!$      !SD
!!$      if (master) then
!!$         print*, '--1.0--'
!!$         print*,  var(1,1,2,2)
!!$         print*,  var(1,1,1,2)
!!$      end if
!!$      !fluctuationmode = .false.
!!$
!!$      !SD
!!$      if (master) then
!!$         print*, '--2.0--'
!!$         print*,  var(1,1,2,2)
!!$         print*,  var(1,1,1,2)
!!$         stop
!!$      end if
!!$
!!$    end subroutine test_print_out

subroutine integrate_ice_TEST(var, var0, flux, dt, dIce, RKstage)

  use type_module, ONLY:nx, ny, nz, nVar, nbx, nby, nbz, include_ice2
  use mpi_module
  !!$  use timeScheme_module
  !!$  use init_module
  !!$  use debug_module
  !!$  use wkb_module
  !!$  use output_module
  !!$  use xweno_module
  use atmosphere_module, ONLY:PStrat01, PStratTilde01
  use boundary_module !, ONLY : setHalos, setBoundary, reconstruction
  use flux_module
  use update_module, ONLY:dIce, ice2Update
  !!$  use ice_module

  implicit none

  real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
      intent(inout) :: var0, var
  real, dimension(- 1:nx, - 1:ny, - 1:nz, 3, nVar) :: flux
  real, dimension(:, :, :, :), allocatable :: dIce ! RK update for ice fields
  real :: dt
  integer :: RKstage

  !SD
  !!$           !TEST
  !!$           call random_number(var(:,:,:,1))
  !!$           var(:,:,:,9) = var(:,:,:,1)
  !!$           var(:,:,:,10) = var(:,:,:,1)
  !!$           var(:,:,:,11) = var(:,:,:,1)

  call setHalos(var, "ice2")
  call setBoundary(var, flux, "ice2")

  !TEST
  !!$           call setHalos( var, "var" )
  !!$           call setBoundary (var, flux, "var")
  !!$           !set z-Boundary for mass as for ice fields
  !!$           var(:,:,-3:0,1) =  var(:,:,-3:0,9)
  !!$           var(:,:,11:13,1) =  var(:,:,11:13,9)

  call reconstruction(var, "ice2")

  !TEST
  !!$           call reconstruction (var, "rho")
  !!$           if ( master ) then
  !!$           print*, sum(abs(nIceTilde(:,:,:,:,:)- rhoTilde(:,:,:,:,:)))
  !!$           print*, 'TEST 1'
  !!$           print*, sum(abs(var(:,:,:,10)- var(:,:,:,1)))
  !!$           print*, sum(abs(var(:,:,:,11)- var(:,:,:,1)))
  !!$           end if

  !** ?? call setHalos( var, "iceTilde" )
  call setBoundary(var, flux, "iceTilde")

  ! Fluxes and Sources

  call ice2Flux(var0, var, flux, "nln", PStrat01, PStratTilde01)
  !call setBoundary (var, flux, "flux") !do nothing for "flux"?
  call ice2Update(var, flux, dt, dIce, RKstage, "expl", 1.)

  !TEST
  !!$           call massFlux (var0,var,flux,"nln",PStrat01,PStratTilde01)
  !!$           !call setBoundary (var, flux, "flux") !do nothing for "flux"?
  !!$           call massUpdate(var, flux, dt, dRho, RKstage, &
  !!$           "rho", "tot", "expl", 1.)

  !!$           !test
  !!$           if ( master ) then
  !!$           is=1 !-3
  !!$           ie=3 !35
  !!$           js=1 !-3
  !!$           je=1  !7
  !!$           ks=1 !-3
  !!$           ke=1 !13
  !*print*, flux(is:ie,js:je,ks:ke, 1, 9)
  !*print*, 'QTilde', qIceTilde(is:ie,js:je,ks:ke, 2, :)
  !!$           print*, 'TEST 2'
  !!$           print*, sum(abs(var(is:ie,js:je,ks:ke, 10) - var(is:ie,js:je,ks:ke, 1)) )
  !!$           print*, sum(abs(var(is:ie,js:je,ks:ke, 9) - var(is:ie,js:je,ks:ke, 1)) )
  !!$           print*, 'Ice-Rho', sum(abs( var(:,:,:,11)- var(:,:,:,1)))
  !!$           print*, 'Diff flux', sum(abs( flux(:,:,:,:,2)- flux(:,:,:,:,1)))
  !!$           !print*, 'stop'
  !!$           !stop
  !!$           !end test
  !!$           end if ! master
  !!$           !print*, 'integrate ice', RKstage

  !!         do RKstage = 1, nStages
  !!$     ! Reconstruction
  !!$
  !!$     call massUpdate(var, flux, dt, dRhop, RKstage, &
  !!$          & "rhop", "lhs", "expl",1.)
  !!$
  !!$     if(topography .and. spongeTFC) then
  !!$        call set_spongeLayer(var, stepFrac(RKstage) * dt, "uvw")
  !!$     end if
  !!$
  !!$  end do

end subroutine integrate_ice_TEST

!----------------------------------------------------
!SD
!test ice_integration
!----------------------------------------------------

!!$  dt = 1 ! [s]
!!$  dt = dt/tRef
!!$  maxIter = 3600./tRef/dt !
!!$
!!$  !output
!!$
!!$  irec = 1
!!$  iOut = 1  !in units of dt
!!$
!!$  if(master) then
!!$     open(44,file='test_output.dat',form="unformatted",access='direct',&
!!$          & recl=nx*ny)
!!$
!!$     print*, ' **************************** '
!!$     print*, 'SET TOPOGRAPHY HEIGHT TO ZERO'
!!$     print*, ' **************************** '
!!$  end if
!!$
!!$  var(:,:,:,1) = 0. !dens fluctuations set to 0
!!$  !prescribe constant vertical velocity
!!$  var(:,:,:,2) = 0. !u=0
!!$  var(:,:,:,3) = 0. !v=0
!!$  !NB this is vertical velocity in TFC
!!$  do k = 1, nz
!!$  var(:,:,k,4) = .2 / uRef/ Pstrat(k) !w \ne 0
!!$  end do
!!$
!!$  var0(:,:,:,2:4) = var(:,:,:,2:4)
!!$
!!$  print*, 'Factor', N2/Tref**2, g**2*L_ice*kappa/Rsp/R_v/300**2
!!$
!!$  do iTime = 0, maxIter
!!$
!!$     !output
!!$     if ( mod(iTime, iOut) == 0) then
!!$
!!$        if (master) then
!!$
!!$           print*, 'iTime', iTime
!!$
!!$           !compute dimensional n, q, qv, S
!!$           !store in var0
!!$           call redim_fields(var, var0)
!!$
!!$
!!$           do ii = 1, 3
!!$              iVar = iVarIce(ii)
!!$              do k = 1, nz
!!$
!!$                 !OUTPUT var0
!!$                 !ofield(1:nx, 1:ny, k, ii) = real(var0(1:nx,1:ny,k, iVar),4)
!!$
!!$                 !ofield(1:nx, 1:ny, k, ii) = real(var(1:nx,1:ny,k, iVar),4)
!!$                 !OUTPUT FLUX
!!$                 ofield(1:nx, 1:ny, k, ii) = real(flux(1:nx,1:ny,k, 3, iVar),4)
!!$                 !ofield(1:nx, 1:ny, k, ii) = real(source(1:nx,1:ny,k, iVar),4)
!!$                 write(44,rec=irec) ((ofield(i, j, k, ii), i=1, nx), j = 1, ny)
!!$                 irec = irec + 1
!!$
!!$              end do !k
!!$
!!$                 !if (itime == 0) then
!!$                 !  print*, 'itime', itime, ii, maxval(ofield(:, :, :, ii)), minval(ofield(:, :, :, ii)), ofield(1,1,1,1)
!!$                 !end if
!!$
!!$           end do !ii
!!$
!!$           !save veritical velocity
!!$
!!$        end if ! master
!!$     end if ! mod
!!$
!!$     do RKstage = 1, nStages
!!$        call integrate_ice(var, var0, flux, "lin", source, dt, dIce, RKstage, PStrat, PStratTilde)
!!$     end do
!!$
!!$  end do !iTime
!!$
!!$  if (master) then
!!$     close(44)
!!$     print*, 'TEST FINISHED ', ' irec=', irec-1
!!$     stop
!!$  end if

!----------------------------------------------------
! end test ice_integration
!----------------------------------------------------
