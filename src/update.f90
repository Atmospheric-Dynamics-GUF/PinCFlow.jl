module update_module

  use type_module
  use timeScheme_module
  use atmosphere_module
  use flux_module
  use poisson_module
  use boundary_module
  use mpi_module
  use mpi

  implicit none

  private
  ! all module objects (variables, functions, subroutines)
  ! are internal to the module by default

  !------------------------
  !   public subroutines
  !------------------------
  public :: momentumPredictor
  public :: massUpdate
  public :: tracerUpdate
  public :: timestep
  public :: init_update
  public :: CoefDySma_update
  public :: Var3DSmthDySma
  public :: iceUpdate_apb, timeUpdate
  public :: setHaloAndBoundary

  public :: smooth_shapiro
  public :: smooth_hor_shapiro

  public :: BGstate_update

  public :: piUpdate ! Update of pi' in compressible model
  public :: bvsUpdate ! Update of N^2 in compressible model

  public :: applyUnifiedSponge

  !-------------------------------
  !    private module variables
  !------------------------------

  contains

  subroutine applyUnifiedSponge(var, dt, time, variable)
    !--------------------------------------
    ! relaxes the predicted solution to
    ! the background state
    !--------------------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var
    real, intent(in) :: dt
    real, intent(in) :: time
    character(len = *), intent(in) :: variable

    ! local variables
    integer :: i, j, k, iVar

    ! relaxation parameters
    real :: alpha, beta

    ! variables for rho
    real :: rho_old, rho_bg, rho_new
    real :: uOld, uBG, uNew
    real :: vOld, vBG, vNew
    real :: wOld, wBG, wNew
    real :: P_old, P_bg, P_new
    real :: pi_old, pi_bg, pi_new, dPdPi

    ! variables for ice
    real :: nAer_bg, nIce_bg, qIce_bg, qv_bg
    real :: T, p

    real, dimension(1:nz) :: sum_local, sum_global

    if(.not. spongeLayer .or. .not. unifiedSponge) return

    select case(variable)

    case("rho")

      ! No sponge is applied to rho in Boussinesq model.
      if(model /= "Boussinesq") then
        do k = 1, nz
          rho_bg = 0.0
          do j = 1, ny
            do i = 1, nx
              alpha = alphaUnifiedSponge(i, j, k)
              rho_old = var%rho(i, j, k)
              beta = 1.0 / (1.0 + alpha * dt)
              rho_new = (1.0 - beta) * rho_bg + beta * rho_old
              var%rho(i, j, k) = rho_new
            end do
          end do
        end do
      end if

    case("P")

      ! No sponge is applied to rho in Boussinesq model.
      if(model /= "Boussinesq") then
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              P_bg = rhoStratTFC(i, j, k) * pStratTFC(i, j, k) / (var%rho(i, &
                  &j, k) + rhoStratTFC(i, j, k))
              alpha = alphaUnifiedSponge(i, j, k)
              P_old = var%P(i, j, k)
              beta = 1.0 / (1.0 + alpha * dt)
              P_new = (1.0 - beta) * P_bg + beta * P_old
              var%P(i, j, k) = P_new
            end do
          end do
        end do
      end if

    case("pi")

      ! No sponge is applied to rho in Boussinesq model.
      if(model /= "Boussinesq") then
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              dPdPi = 1 / (gamma - 1) * (Rsp / pref) ** (1 - gamma) * var%P(i, &
                  &j, k) ** (2 - gamma)
              pi_bg = rhoStratTFC(i, j, k) * pStratTFC(i, j, k) / (var%rho(i, &
                  &j, k) + rhoStratTFC(i, j, k)) / dPdPi
              alpha = alphaUnifiedSponge(i, j, k)
              pi_old = var%pi(i, j, k)
              pi_new = pi_old - alpha * dt * (pStratTFC(i, j, k) / dPdPi &
                  &- pi_bg)
              var%pi(i, j, k) = pi_new
            end do
          end do
        end do
      end if

    case("rhop")

      rho_bg = 0.0
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            alpha = alphaUnifiedSponge(i, j, k)
            rho_old = var%rhop(i, j, k)
            beta = 1.0 / (1.0 + alpha * dt)
            rho_new = (1.0 - beta) * rho_bg + beta * rho_old
            var%rhop(i, j, k) = rho_new
          end do
        end do
      end do

    case("uvw")

      ! Sponge for zonal wind

      ! Determine relaxation wind.
      if(relax_to_mean) then
        ! Compute local horizontal sum.
        do k = 1, nz
          sum_local(k) = sum(var%u(1:nx, 1:ny, k))
        end do

        ! Compute global sum and average.
        call mpi_allreduce(sum_local(1), sum_global(1), nz, &
            &mpi_double_precision, mpi_sum, comm, ierror)
        sum_global = sum_global / (sizeX * sizeY)
      else
        uBG = backgroundFlow_dim(1) / uRef
        if(relaxation_period > 0.0) uBG = uBG * (1.0 + relaxation_amplitude &
            &* sin(2.0 * pi * time / relaxation_period * tRef))
      end if

      do k = 1, nz
        if(relax_to_mean) uBG = sum_global(k)
        do j = 1, ny
          do i = 1, nx
            alpha = 0.5 * (alphaUnifiedSponge(i, j, k) + alphaUnifiedSponge(i &
                &+ 1, j, k))
            uOld = var%u(i, j, k)
            beta = 1.0 / (1.0 + alpha * dt)
            uNew = (1.0 - beta) * uBG + beta * uOld
            var%u(i, j, k) = uNew
          end do
        end do
      end do

      ! Sponge for meridional wind

      ! Determine relaxation wind.
      if(relax_to_mean) then
        ! Compute local horizontal sum.
        do k = 1, nz
          sum_local(k) = sum(var%v(1:nx, 1:ny, k))
        end do

        ! Compute global sum and average.
        call mpi_allreduce(sum_local(1), sum_global(1), nz, &
            &mpi_double_precision, mpi_sum, comm, ierror)
        sum_global = sum_global / (sizeX * sizeY)
      else
        vBG = backgroundFlow_dim(2) / uRef
        if(relaxation_period > 0.0) vBG = vBG * (1.0 + relaxation_amplitude &
            &* sin(2.0 * pi * time / relaxation_period * tRef))
      end if

      do k = 1, nz
        if(relax_to_mean) vBG = sum_global(k)
        do j = 1, ny
          do i = 1, nx
            alpha = 0.5 * (alphaUnifiedSponge(i, j, k) + alphaUnifiedSponge(i, &
                &j + 1, k))
            vOld = var%v(i, j, k)
            beta = 1.0 / (1.0 + alpha * dt)
            vNew = (1.0 - beta) * vBG + beta * vOld
            var%v(i, j, k) = vNew
          end do
        end do
      end do

      ! Sponge for vertical wind

      ! Determine relaxation wind.
      if(relax_to_mean) then
        ! Compute local horizontal sum.
        do k = 1, nz
          sum_local(k) = sum(var%w(1:nx, 1:ny, k))
        end do

        ! Compute global sum and average.
        call mpi_allreduce(sum_local(1), sum_global(1), nz, &
            &mpi_double_precision, mpi_sum, comm, ierror)
        sum_global = sum_global / (sizeX * sizeY)
      else
        wBG = backgroundFlow_dim(3) / uRef
        if(relaxation_period > 0.0) wBG = wBG * (1.0 + relaxation_amplitude &
            &* sin(2.0 * pi * time / relaxation_period * tRef))
      end if

      do k = 1, nz
        if(relax_to_mean) wBG = sum_global(k)
        do j = 1, ny
          do i = 1, nx
            if(topography) then
              alpha = (jac(i, j, k + 1) * alphaUnifiedSponge(i, j, k) + jac(i, &
                  &j, k) * alphaUnifiedSponge(i, j, k + 1)) / (jac(i, j, k) &
                  &+ jac(i, j, k + 1))
            else
              alpha = 0.5 * (alphaUnifiedSponge(i, j, k) &
                  &+ alphaUnifiedSponge(i, j, k + 1))
            end if
            wOld = var%w(i, j, k)
            beta = 1.0 / (1.0 + alpha * dt)
            wNew = (1.0 - beta) * wBG + beta * wOld
            var%w(i, j, k) = wNew
          end do
        end do
      end do

    case default
      stop "applyUnifiedSponge: Unknown variable!"
    end select

  end subroutine applyUnifiedSponge

  !---------------------------------------------------------------------

  subroutine momentumPredictor(var, flux, force, dt, q, m, mmp_mod, int_mod, &
      &facray)
    !----------------------------------
    !  calculates the velocities u^*
    !----------------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var

    type(flux_type), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux)
    ! dir = 1..3 > f,g,h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW

    ! mmp_mod decides, which part of the momentum equation is to be used:
    ! tot => total momentum equation
    ! lhs => only advection and molecular and turbulent viscous fluxes on
    !        the left-hand side of the equation
    ! rhs => only pressure-gradient, Coriolis and gravitational force on
    !        the right-hand side of the equation

    ! int_mod discriminates between implicit and explicit time stepping:
    ! expl => explicit time stepping
    !         (always the case for integration of the lhs)
    !         RK sub step for integration of the lhs
    !         Euler step for the rhs of the momentum equations
    ! impl => implicit-time-step part without pressure-gradient term
    !         (only for rhs of the momentum equations)

    ! facray multiplies the Rayleigh-damping terms so that they are only
    ! handled in the implicit time stepping (sponge and immersed boundary)
    character(len = *), intent(in) :: mmp_mod, int_mod

    ! volume forces
    ! mmp_mod = tot =>
    ! 1) gravitational / buoyancy force (cell-centered)
    ! 2) Coriolis force (cell-centered)
    ! 3) WKB wave driving (cell-centered)
    ! mmp_mod = rhs =>
    ! 1) WKB wave driving (cell-centered)
    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, 3), intent(in) :: force

    !UAC real, intent(in) :: dt
    real, intent(in) :: dt, facray
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 3), &
        &intent(inout) :: q
    integer, intent(in) :: m

    !UAB
    logical :: spongeLayer_s, topography_s
    !UAE

    ! local variables
    real :: fL, fR, gB, gF, hD, hU
    ! flux Left/Right, Backward/Forward, Downward/Upward

    ! usave to keep the new u until v has been updated as well
    ! (for mmp_mod = rhs)
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: usave

    ! other stuff
    real :: rhoM_1, rhoM ! rho(m-1), rho(m)
    real :: drho_e
    real :: uM_1, vM_1, wM_1 ! u(m-1), v(m-1) and w(m-1)
    real :: momM_1, momM ! momentum at t(m-1) and t(m)
    real :: piR, piL, piF, piB, piU, piD
    real :: fluxDiff ! conv. and viscous flux contr.
    real :: piGrad ! pressure gradient
    real :: piGradx, piGrady ! horizontal pressure-gradient
    ! components
    real :: F ! update part for Runge-Kutta step
    real :: uAst, vAst, wAst ! predicted velocities u*, v* and w*
    real :: uhorx, vhory, wvert

    integer :: i, j, k
    integer :: i0, i1, j0, j1, k0, k1

    ! gravity
    real :: thetaEdge
    real :: volForce

    ! local interpolation values
    real :: pBarEdge ! stratified background pressure interpolated
    real :: thetaBar ! stratified pot. temp. interpolated
    real :: rhoEdge ! interpolated density to velocity edge

    ! TFC variables
    real :: jacEdgeR, jacEdgeF, jacEdgeU
    real :: pEdgeR, pEdgeF, pEdgeU
    real :: met13EdgeR, met23EdgeF, met13EdgeU, met23EdgeU, met33EdgeU
    real :: piEdgeR, piUEdgeR, piUUEdgeR, piDEdgeR, piDDEdgeR, piEdgeF, &
        &piUEdgeF, piUUEdgeF, piDEdgeF, piDDEdgeF, piREdgeU, piLEdgeU, &
        &piFEdgeU, piBEdgeU
    real :: rhoStratEdgeR, rhoStratEdgeF, rhoStratEdgeU
    real :: vC, vR, uC, uF, vU, uU
    real, dimension(0:1, 0:1) :: fluxDiffU, fluxDiffV
    real, dimension(1:nx, 1:ny, 1:nz) :: fluxDiffW
    ! real :: metEdgeR1, metEdgeR2, metEdgeL1, metEdgeL2, &
    !         metEdgeF1, metEdgeF2, metEdgeB1, metEdgeB2
    integer :: ll, mm

    ! classical RK3
    real :: rhoM_0, uM_0, vM_0, wM_0, momM_0

    ! non-dimensional Corilois parameter (= inverse Rossby number)
    !FS real :: f_cor_nd
    real, dimension(0:ny + 1) :: f_cor_nd

    !real :: rho, rhop, rhou, rhov, rhow, facu, facv, facw, facr, pstw, buoy
    real :: rho, rhop, rhou, rhov, rhow, facu, facv, facw, facr, pstw, pstw_0, &
        &buoy
    real :: rho10, rho01
    real :: rhov0m, rhov00, rhov1m, rhov10
    real :: rhou00, rhoum0, rhou01, rhoum1
    real :: rho000, rho001
    real :: volfcx, volfcy, volfcz
    real :: bvsstw

    ! SK compressible: JP on interfaces right, left, forward, backward, upward, downward at k and k+1
    real :: JPR, JPL, JPF, JPB, JPU, JPD, JPUR, JPUL, JPUF, JPUB

    real :: rho_e, rhou_e, rhov_e, rhow_e, pstw_e
    real :: rho10_e, rho01_e
    real :: rhov0m_e, rhov00_e, rhov1m_e, rhov10_e
    real :: rhou00_e, rhoum0_e, rhou01_e, rhoum1_e
    real :: rho000_e, rho001_e
    real :: piR_e, piL_e, piF_e, piB_e, piU_e, piD_e

    real :: f_cor_v

    ! integer :: i00,j00

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: heat

    real :: rhop_0, rhop_1
    real :: rho_p0, rho_p1

    real, dimension(- nbz:nz + nbz) :: w_0
    real, dimension(- nbz:nz + nbz) :: S_bar
    real :: heat_flc, heat0, heat1

    real :: ymax, yloc, ymin

    if(corset == 'periodic') then
      ymax = ly_dim(1) / lRef
      ymin = ly_dim(0) / lRef

      j0 = js + nby - 1

      do j = 0, ny + 1
        yloc = y(j + j0)

        f_cor_nd(j) = - 4. * pi / 8.64e4 * tRef * cos(2. * pi * (yloc - ymin) &
            &/ (ymax - ymin))
      end do
    else if(corset == 'constant') then
      f_cor_nd(0:ny + 1) = f_Coriolis_dim * tRef
    else
      stop 'ERROR: wrong corset'
    end if

    ! if(topography) then
    !    i00=is+nbx-1
    !    j00=js+nby-1
    ! end if

    !if (int_mod == "impl") then
    !   ! environmental heating
    !   call calculate_heating(var,flux,heat)

    !   ! heating by GW entropy-flux convergence
    !   if (raytracer) heat(:,:,:) = heat(:,:,:) + var(:,:,:,8)
    !end if

    if(correctDivError) then
      print *, 'ERROR: correction divergence error not allowed'
      stop
    end if

    ! init q
    if(m == 1) q = 0.

    if(mmp_mod == 'rhs') then
      if(int_mod == 'expl') then
        spongeLayer_s = spongeLayer
        spongeLayer = .false.
      else if(int_mod == 'impl') then
        kr_sp = kr_sp * facray
        kr_sp_w = kr_sp_w * facray
        alprlx = alprlx * facray
        if(topography) then
          kr_sp_tfc = kr_sp_tfc * facray
          kr_sp_w_tfc = kr_sp_w_tfc * facray
        end if
      end if
    end if

    ! -------------------------------------
    !            predict u -> u*
    ! -------------------------------------

    select case(xBoundary)

    case("solid_wall")
      i0 = 1
      i1 = nx - 1
    case("periodic")
      i0 = 0
      i1 = nx
    case default
      stop "momentumPredictor: unknown case xBoundary."
    end select

    if(mmp_mod == "tot" .or. mmp_mod == "lhs") then
      if(int_mod /= "expl") then
        stop 'ERROR: wrong int_mod for mmp_mod = tot or mmp_mod = lhs'
      end if

      do k = 1, nz
        do j = 1, ny
          do i = i0, i1

            !--- convective fluxes -> conv
            fR = flux%u(i, j, k, 1)
            fL = flux%u(i - 1, j, k, 1)
            gF = flux%u(i, j, k, 2)
            gB = flux%u(i, j - 1, k, 2)
            hU = flux%u(i, j, k, 3)
            hD = flux%u(i, j, k - 1, 3)
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz ! diverg.

            ! Adjust zonal momentum flux divergence.
            if(topography) then
              jacEdgeR = 0.5 * (jac(i, j, k) + jac(i + 1, j, k))
              fluxDiff = fluxDiff / jacEdgeR
            end if

            volForce = 0.

            if(mmp_mod == "tot") then
              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piR = var%pi(i + 1, j, k) - var_env%pi(i + 1, j, k)
                piL = var%pi(i, j, k) - var_env%pi(i, j, k)
              else
                piR = var%pi(i + 1, j, k)
                piL = var%pi(i, j, k)
              end if

              if(topography) then
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) - var_env%pi(:, :, :)
                end if
                ! Compute values at cell edges.
                pEdgeR = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i + 1, j, k))
                met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))
                ! Compute pressure gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  piUUEdgeR = 0.5 * (var%pi(i, j, k + 2) + var%pi(i + 1, j, k &
                      &+ 2))
                  piUEdgeR = 0.5 * (var%pi(i, j, k + 1) + var%pi(i + 1, j, k &
                      &+ 1))
                  piEdgeR = 0.5 * (var%pi(i, j, k) + var%pi(i + 1, j, k))
                  piGrad = kappaInv * MaInv2 * pEdgeR * ((piR - piL) / dx &
                      &+ met13EdgeR * (- piUUEdgeR + 4.0 * piUEdgeR - 3.0 &
                      &* piEdgeR) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeR = 0.5 * (var%pi(i, j, k - 2) + var%pi(i + 1, j, k &
                      &- 2))
                  piDEdgeR = 0.5 * (var%pi(i, j, k - 1) + var%pi(i + 1, j, k &
                      &- 1))
                  piEdgeR = 0.5 * (var%pi(i, j, k) + var%pi(i + 1, j, k))
                  piGrad = kappaInv * MaInv2 * pEdgeR * ((piR - piL) / dx &
                      &+ met13EdgeR * (piDDEdgeR - 4.0 * piDEdgeR + 3.0 &
                      &* piEdgeR) * 0.5 / dz)
                else
                  piUEdgeR = 0.5 * (var%pi(i, j, k + 1) + var%pi(i + 1, j, k &
                      &+ 1))
                  piDEdgeR = 0.5 * (var%pi(i, j, k - 1) + var%pi(i + 1, j, k &
                      &- 1))
                  piGrad = kappaInv * MaInv2 * pEdgeR * ((piR - piL) / dx &
                      &+ met13EdgeR * (piUEdgeR - piDEdgeR) * 0.5 / dz)
                end if
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) + var_env%pi(:, :, :)
                end if
              else
                piGrad = kappaInv * MaInv2 * Pstrat(k) * (piR - piL) / dx

                if(TestCase == "baroclinic_LC") then !FS
                  piGrad = piGrad + kappaInv * MaInv2 * (PStrat(k) &
                      &- pStrat_0(k)) * (var_env%pi(i + 1, j, k) &
                      &- var_env%pi(i, j, k)) / dx
                end if
              end if

              !---- volume forces
              volForce = 0.5 * (force(i, j, k, 1) + force(i + 1, j, k, 1))

              if(TestCase == "baroclinic_LC") then
                select case(model)
                case("pseudo_incompressible", "compressible")
                  rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i + 1, j, k))

                  if(topography) then
                    rhoM_1 = rhoM_1 + 0.5 * (rhoStratTFC(i, j, k) &
                        &+ rhoStratTFC(i + 1, j, k))
                  else
                    rhoM_1 = rhoM_1 + rhoStrat(k)
                  end if
                case("Boussinesq")
                  rhoM_1 = rho00
                case default
                  stop "momentumPredictor: unkown model."
                end select

                volforce = volforce - rhoM_1 * RoInv(j) * 0.25 * (var_env%v(i, &
                    &j - 1, k) + var_env%v(i + 1, j - 1, k) + var_env%v(i, j, &
                    &k) + var_env%v(i + 1, j, k))

              end if

            end if

            if(TestCase == "baroclinic_LC") then
              if(background == "HeldSuarez") then
                ! Rayleigh damping

                select case(model)
                case("pseudo_incompressible", "compressible")
                  rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i + 1, j, k))

                  if(topography) then
                    rhoM_1 = rhoM_1 + 0.5 * (rhoStratTFC(i, j, k) &
                        &+ rhoStratTFC(i + 1, j, k))
                  else
                    rhoM_1 = rhoM_1 + rhoStrat(k)
                  end if
                case("Boussinesq")
                  rhoM_1 = rho00
                case default
                  stop "momentumPredictor: unkown model."
                end select

                volForce = volForce - kv_hs(j, k) * rhoM_1 * (var%u(i, j, k) &
                    &- var_env%u(i, j, k))

              end if
            end if

            ! Explicit integration of Coriolis force in TFC.
            if(topography .and. mmp_mod == "lhs") then
              uOldTFC(i, j, k) = var%u(i, j, k)
              vC = 0.5 * (var%v(i, j, k) + var%v(i, j - 1, k))
              vR = 0.5 * (var%v(i + 1, j, k) + var%v(i + 1, j - 1, k))
              if(testCase == "baroclinic_LC") then
                vC = vC - 0.5 * (var_env%v(i, j, k) + var_env%v(i, j - 1, k))
                vR = vR - 0.5 * (var_env%v(i + 1, j, k) + var_env%v(i + 1, j &
                    &- 1, k))
              end if
              volForce = volForce + 0.5 * f_cor_nd(j) * ((rhoOld(i, j, k) &
                  &+ rhoStratTFC(i, j, k)) * vC + (rhoOld(i + 1, j, k) &
                  &+ rhoStratTFC(i + 1, j, k)) * vR)
            end if

            !--------------------
            !   d/dt ... = F(phi) (RHS of ODE)
            !--------------------
            ! fluxDiff -> convective and viscous fluxes
            ! piGrad   -> pressure gradient along x scaled with 1/Ma^2
            ! volForce -> Gravity, Coriolis
            if(mmp_mod == "tot") then
              F = - fluxDiff - piGrad + volForce
            else if(mmp_mod == "lhs") then
              F = - fluxDiff + volForce !200413
            else
              stop 'ERROR: wrong mmp_mod'
            end if

            ! interpolated density
            select case(model)

            case("pseudo_incompressible", "compressible")

              rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i + 1, j, k))
              rhoM = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k))

              if(topography) then
                rhoStratEdgeR = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i &
                    &+ 1, j, k))
                rhoM_1 = rhoM_1 + rhoStratEdgeR
                rhoM = rhoM + rhoStratEdgeR
              else
                rhoM_1 = rhoM_1 + rhoStrat(k)
                rhoM = rhoM + rhoStrat(k)
              end if
            case("Boussinesq")
              rhoM_1 = rho00
              rhoM = rho00
            case default
              stop "momentumPredictor: unkown case model."
            end select

            ! velocity and momentum at t(m-1)
            uM_1 = var%u(i, j, k)
            momM_1 = rhoM_1 * uM_1

            ! q(m-1) -> q(m)

            q(i, j, k, 1) = dt * F + alphaRK(m) * q(i, j, k, 1)

            ! rhoU(m-1) -> rhoU(m)
            momM = momM_1 + betaRK(m) * q(i, j, k, 1)

            ! calc u(m,*)
            uAst = momM / rhoM

            ! uAst -> var
            var%u(i, j, k) = uAst
          end do
        end do
      end do
    else if(mmp_mod == "rhs") then
      if(int_mod == "expl") then
        do k = 1, nz
          do j = 1, ny
            do i = i0, i1
              rhou = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k))
              if(topography) then
                rhoStratEdgeR = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i &
                    &+ 1, j, k))
                rhou = rhou + rhoStratEdgeR
              else
                rhou = rhou + rhoStrat(k)
              end if

              if(TestCase == "baroclinic_LC") then
                rhou_e = 0.5 * (var_env%rho(i, j, k) + var_env%rho(i + 1, j, k))
                if(topography) then
                  rhou_e = rhou_e + 0.5 * (rhoStratTFC(i, j, k) &
                      &+ rhoStratTFC(i + 1, j, k))
                else
                  rhou_e = rhou_e + rhoStrat_0(k)
                end if
              end if

              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piR = var%pi(i + 1, j, k) - var_env%pi(i + 1, j, k)
                piL = var%pi(i, j, k) - var_env%pi(i, j, k)

                piR_e = var_env%pi(i + 1, j, k)
                piL_e = var_env%pi(i, j, k)
              else
                piR = var%pi(i + 1, j, k)
                piL = var%pi(i, j, k)
              end if

              if(topography) then
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) - var_env%pi(:, :, :)
                end if
                ! Compute values at cell edges.
                pEdgeR = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i + 1, j, k))
                met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))
                ! Compute pressure gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  piUUEdgeR = 0.5 * (var%pi(i, j, k + 2) + var%pi(i + 1, j, k &
                      &+ 2))
                  piUEdgeR = 0.5 * (var%pi(i, j, k + 1) + var%pi(i + 1, j, k &
                      &+ 1))
                  piEdgeR = 0.5 * (var%pi(i, j, k) + var%pi(i + 1, j, k))
                  piGrad = kappaInv * MaInv2 * pEdgeR / rhou * ((piR - piL) &
                      &/ dx + met13EdgeR * (- piUUEdgeR + 4.0 * piUEdgeR - 3.0 &
                      &* piEdgeR) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeR = 0.5 * (var%pi(i, j, k - 2) + var%pi(i + 1, j, k &
                      &- 2))
                  piDEdgeR = 0.5 * (var%pi(i, j, k - 1) + var%pi(i + 1, j, k &
                      &- 1))
                  piEdgeR = 0.5 * (var%pi(i, j, k) + var%pi(i + 1, j, k))
                  piGrad = kappaInv * MaInv2 * pEdgeR / rhou * ((piR - piL) &
                      &/ dx + met13EdgeR * (piDDEdgeR - 4.0 * piDEdgeR + 3.0 &
                      &* piEdgeR) * 0.5 / dz)
                else
                  piUEdgeR = 0.5 * (var%pi(i, j, k + 1) + var%pi(i + 1, j, k &
                      &+ 1))
                  piDEdgeR = 0.5 * (var%pi(i, j, k - 1) + var%pi(i + 1, j, k &
                      &- 1))
                  piGrad = kappaInv * MaInv2 * pEdgeR / rhou * ((piR - piL) &
                      &/ dx + met13EdgeR * (piUEdgeR - piDEdgeR) * 0.5 / dz)
                end if
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) + var_env%pi(:, :, :)
                end if
              else
                piGrad = kappaInv * MaInv2 * Pstrat(k) / rhou * (piR - piL) / dx
              end if

              if(TestCase == "baroclinic_LC") then !FS
                if(topography) then
                  ! Compute values at cell edges.
                  pEdgeR = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i + 1, j, k))
                  met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, &
                      &3))
                  ! Compute pressure gradient component.
                  if(k == 1 .and. zBoundary == "solid_wall") then
                    piUUEdgeR = 0.5 * (var_env%pi(i, j, k + 2) + var_env%pi(i &
                        &+ 1, j, k + 2))
                    piUEdgeR = 0.5 * (var_env%pi(i, j, k + 1) + var_env%pi(i &
                        &+ 1, j, k + 1))
                    piEdgeR = 0.5 * (var_env%pi(i, j, k) + var_env%pi(i + 1, &
                        &j, k))
                    piGrad = piGrad + kappaInv * MaInv2 * (pEdgeR / rhou &
                        &- pEdgeR / rhou_e) * ((piR_e - piL_e) / dx &
                        &+ met13EdgeR * (- piUUEdgeR + 4.0 * piUEdgeR - 3.0 &
                        &* piEdgeR) * 0.5 / dz)
                  else if(k == nz .and. zBoundary == "solid_wall") then
                    piDDEdgeR = 0.5 * (var_env%pi(i, j, k - 2) + var_env%pi(i &
                        &+ 1, j, k - 2))
                    piDEdgeR = 0.5 * (var_env%pi(i, j, k - 1) + var_env%pi(i &
                        &+ 1, j, k - 1))
                    piEdgeR = 0.5 * (var_env%pi(i, j, k) + var_env%pi(i + 1, &
                        &j, k))
                    piGrad = piGrad + kappaInv * MaInv2 * (pEdgeR / rhou &
                        &- pEdgeR / rhou_e) * ((piR_e - piL_e) / dx &
                        &+ met13EdgeR * (piDDEdgeR - 4.0 * piDEdgeR + 3.0 &
                        &* piEdgeR) * 0.5 / dz)
                  else
                    piUEdgeR = 0.5 * (var_env%pi(i, j, k + 1) + var_env%pi(i &
                        &+ 1, j, k + 1))
                    piDEdgeR = 0.5 * (var_env%pi(i, j, k - 1) + var_env%pi(i &
                        &+ 1, j, k - 1))
                    piGrad = piGrad + kappaInv * MaInv2 * (pEdgeR / rhou &
                        &- pEdgeR / rhou_e) * ((piR_e - piL_e) / dx &
                        &+ met13EdgeR * (piUEdgeR - piDEdgeR) * 0.5 / dz)
                  end if
                else
                  piGrad = piGrad + kappaInv * MaInv2 * (Pstrat(k) / rhou &
                      &- Pstrat_0(k) / rhou_e) * (piR_e - piL_e) / dx
                end if
              end if

              ! gravity-wave forcing
              if(raytracer .or. testCase == "mountainwave") then
                volfcx = 0.5 * (force(i, j, k, 1) + force(i + 1, j, k, 1))
              else
                volfcx = 0.0
              end if

              ! ustar
              if(TestCase == "baroclinic_LC") then
                uhorx = var%u(i, j, k) - var_env%u(i, j, k)
              else
                uhorx = var%u(i, j, k)
              end if

              if(topography) then
                ! Coriolis force is integrated on LHS.
                if(model == "compressible") then ! Muliply with JP
                  JPR = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i + 1, j, &
                      &k) * var%P(i + 1, j, k))
                  uAst = uhorx + dt * (- piGrad + volfcx / rhou) * JPR
                else
                  uAst = uhorx + dt * (- piGrad + volfcx / rhou)
                end if
              else
                vhory = 0.25 * (var%v(i, j - 1, k) + var%v(i, j, k) + var%v(i &
                    &+ 1, j - 1, k) + var%v(i + 1, j, k))

                uAst = uhorx + dt * (f_cor_nd(j) * vhory - piGrad + volfcx &
                    &/ rhou)
              end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  uAst = uAst - dt * kv_hs(j, k) * uhorx
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                if(topography) then
                  uAst = uAst - dt * 0.5 * (kr_sp_tfc(i, j, k) + kr_sp_tfc(i &
                      &+ 1, j, k)) * uhorx
                else
                  uAst = uAst - dt * kr_sp(j, k) * uhorx
                end if
              end if

              usave(i, j, k) = uAst

              if(TestCase == "baroclinic_LC") then
                usave(i, j, k) = usave(i, j, k) + var_env%u(i, j, k)
              end if
            end do
          end do
        end do
      else if(int_mod == "impl") then
        do k = 1, nz
          do j = 1, ny
            do i = i0, i1
              rhou = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k))

              rhov0m = 0.5 * (var%rho(i, j, k) + var%rho(i, j - 1, k))
              rhov00 = 0.5 * (var%rho(i, j + 1, k) + var%rho(i, j, k))
              rhov1m = 0.5 * (var%rho(i + 1, j, k) + var%rho(i + 1, j - 1, k))
              rhov10 = 0.5 * (var%rho(i + 1, j + 1, k) + var%rho(i + 1, j, k))

              if(topography) then
                rhoStratEdgeR = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i &
                    &+ 1, j, k))
                rhou = rhou + rhoStratEdgeR
              else
                rhou = rhou + rhoStrat(k)
                rhov0m = rhov0m + rhoStrat(k)
                rhov00 = rhov00 + rhoStrat(k)
                rhov1m = rhov1m + rhoStrat(k)
                rhov10 = rhov10 + rhoStrat(k)
              end if

              if(TestCase == "baroclinic_LC") then
                rhou_e = 0.5 * (var_env%rho(i, j, k) + var_env%rho(i + 1, j, k))

                rhov0m_e = 0.5 * (var_env%rho(i, j, k) + var_env%rho(i, j - 1, &
                    &k))
                rhov00_e = 0.5 * (var_env%rho(i, j + 1, k) + var_env%rho(i, j, &
                    &k))
                rhov1m_e = 0.5 * (var_env%rho(i + 1, j, k) + var_env%rho(i &
                    &+ 1, j - 1, k))
                rhov10_e = 0.5 * (var_env%rho(i + 1, j + 1, k) + var_env%rho(i &
                    &+ 1, j, k))

                if(topography) then
                  rhou_e = rhou_e + 0.5 * (rhoStratTFC(i, j, k) &
                      &+ rhoStratTFC(i + 1, j, k))
                else
                  rhou_e = rhou_e + rhoStrat_0(k)
                  rhov0m_e = rhov0m_e + rhoStrat_0(k)
                  rhov00_e = rhov00_e + rhoStrat_0(k)
                  rhov1m_e = rhov1m_e + rhoStrat_0(k)
                  rhov10_e = rhov10_e + rhoStrat_0(k)
                end if
              end if

              !--- pressure gradient terms -> piGradx, piGrady
              if(TestCase == "baroclinic_LC") then
                piR = var%pi(i + 1, j, k) - var_env%pi(i + 1, j, k)
                piL = var%pi(i, j, k) - var_env%pi(i, j, k)

                piR_e = var_env%pi(i + 1, j, k)
                piL_e = var_env%pi(i, j, k)
              else
                piR = var%pi(i + 1, j, k)
                piL = var%pi(i, j, k)
              end if

              if(topography) then
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) - var_env%pi(:, :, :)
                end if
                ! Compute values at cell edges.
                pEdgeR = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i + 1, j, k))
                met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))
                ! Compute pressure gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  piUUEdgeR = 0.5 * (var%pi(i, j, k + 2) + var%pi(i + 1, j, k &
                      &+ 2))
                  piUEdgeR = 0.5 * (var%pi(i, j, k + 1) + var%pi(i + 1, j, k &
                      &+ 1))
                  piEdgeR = 0.5 * (var%pi(i, j, k) + var%pi(i + 1, j, k))
                  piGradX = kappaInv * MaInv2 * pEdgeR / rhou * ((piR - piL) &
                      &/ dx + met13EdgeR * (- piUUEdgeR + 4.0 * piUEdgeR - 3.0 &
                      &* piEdgeR) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeR = 0.5 * (var%pi(i, j, k - 2) + var%pi(i + 1, j, k &
                      &- 2))
                  piDEdgeR = 0.5 * (var%pi(i, j, k - 1) + var%pi(i + 1, j, k &
                      &- 1))
                  piEdgeR = 0.5 * (var%pi(i, j, k) + var%pi(i + 1, j, k))
                  piGradX = kappaInv * MaInv2 * pEdgeR / rhou * ((piR - piL) &
                      &/ dx + met13EdgeR * (piDDEdgeR - 4.0 * piDEdgeR + 3.0 &
                      &* piEdgeR) * 0.5 / dz)
                else
                  piUEdgeR = 0.5 * (var%pi(i, j, k + 1) + var%pi(i + 1, j, k &
                      &+ 1))
                  piDEdgeR = 0.5 * (var%pi(i, j, k - 1) + var%pi(i + 1, j, k &
                      &- 1))
                  piGradX = kappaInv * MaInv2 * pEdgeR / rhou * ((piR - piL) &
                      &/ dx + met13EdgeR * (piUEdgeR - piDEdgeR) * 0.5 / dz)
                end if
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) + var_env%pi(:, :, :)
                end if
              else
                piGradx = kappaInv * MaInv2 * Pstrat(k) / rhou * (piR - piL) &
                    &/ dx !FS
              end if

              if(TestCase == "baroclinic_LC") then !FS
                if(topography) then
                  ! Compute values at cell edges.
                  pEdgeR = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i + 1, j, k))
                  met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, &
                      &3))
                  ! Compute pressure gradient component.
                  if(k == 1 .and. zBoundary == "solid_wall") then
                    piUUEdgeR = 0.5 * (var_env%pi(i, j, k + 2) + var_env%pi(i &
                        &+ 1, j, k + 2))
                    piUEdgeR = 0.5 * (var_env%pi(i, j, k + 1) + var_env%pi(i &
                        &+ 1, j, k + 1))
                    piEdgeR = 0.5 * (var_env%pi(i, j, k) + var_env%pi(i + 1, &
                        &j, k))
                    piGradX = piGradX + kappaInv * MaInv2 * (pEdgeR / rhou &
                        &- pEdgeR / rhou_e) * ((piR_e - piL_e) / dx &
                        &+ met13EdgeR * (- piUUEdgeR + 4.0 * piUEdgeR - 3.0 &
                        &* piEdgeR) * 0.5 / dz)
                  else if(k == nz .and. zBoundary == "solid_wall") then
                    piDDEdgeR = 0.5 * (var_env%pi(i, j, k - 2) + var_env%pi(i &
                        &+ 1, j, k - 2))
                    piDEdgeR = 0.5 * (var_env%pi(i, j, k - 1) + var_env%pi(i &
                        &+ 1, j, k - 1))
                    piEdgeR = 0.5 * (var_env%pi(i, j, k) + var_env%pi(i + 1, &
                        &j, k))
                    piGradX = piGradX + kappaInv * MaInv2 * (pEdgeR / rhou &
                        &- pEdgeR / rhou_e) * ((piR_e - piL_e) / dx &
                        &+ met13EdgeR * (piDDEdgeR - 4.0 * piDEdgeR + 3.0 &
                        &* piEdgeR) * 0.5 / dz)
                  else
                    piUEdgeR = 0.5 * (var_env%pi(i, j, k + 1) + var_env%pi(i &
                        &+ 1, j, k + 1))
                    piDEdgeR = 0.5 * (var_env%pi(i, j, k - 1) + var_env%pi(i &
                        &+ 1, j, k - 1))
                    piGradX = piGradX + kappaInv * MaInv2 * (pEdgeR / rhou &
                        &- pEdgeR / rhou_e) * ((piR_e - piL_e) / dx &
                        &+ met13EdgeR * (piUEdgeR - piDEdgeR) * 0.5 / dz)
                  end if
                else
                  piGradx = piGradx + kappaInv * MaInv2 * (PStrat(k) / rhou &
                      &- pStrat_0(k) / rhou_e) * (piR_e - piL_e) / dx
                end if
              end if

              if(.not. topography) then
                if(TestCase == "baroclinic_LC") then
                  piGrady = kappaInv * MaInv2 * 0.25 * (Pstrat(k) / rhov0m &
                      &* (var%pi(i, j, k) - var%pi(i, j - 1, k) &
                      &- var_env%pi(i, j, k) + var_env%pi(i, j - 1, k)) / dy &
                      &+ Pstrat(k) / rhov00 * (var%pi(i, j + 1, k) - var%pi(i, &
                      &j, k) - var_env%pi(i, j + 1, k) + var_env%pi(i, j, k)) &
                      &/ dy + Pstrat(k) / rhov1m * (var%pi(i + 1, j, k) &
                      &- var%pi(i + 1, j - 1, k) - var_env%pi(i + 1, j, k) &
                      &+ var_env%pi(i + 1, j - 1, k)) / dy + Pstrat(k) &
                      &/ rhov10 * (var%pi(i + 1, j + 1, k) - var%pi(i + 1, j, &
                      &k) - var_env%pi(i + 1, j + 1, k) + var_env%pi(i + 1, j, &
                      &k)) / dy)

                  piGrady = piGrady + kappaInv * MaInv2 * 0.25 * ((Pstrat(k) &
                      &/ rhov0m - Pstrat_0(k) / rhov0m_e) * (var_env%pi(i, j, &
                      &k) - var_env%pi(i, j - 1, k)) / dy + (Pstrat(k) &
                      &/ rhov00 - Pstrat_0(k) / rhov00_e) * (var_env%pi(i, j &
                      &+ 1, k) - var_env%pi(i, j, k)) / dy + (Pstrat(k) &
                      &/ rhov1m - Pstrat_0(k) / rhov1m_e) * (var_env%pi(i + 1, &
                      &j, k) - var_env%pi(i + 1, j - 1, k)) / dy + (Pstrat(k) &
                      &/ rhov10 - Pstrat_0(k) / rhov10_e) * (var_env%pi(i + 1, &
                      &j + 1, k) - var_env%pi(i + 1, j, k)) / dy)
                else
                  piGrady = kappaInv * MaInv2 * 0.25 * (Pstrat(k) / rhov0m &
                      &* (var%pi(i, j, k) - var%pi(i, j - 1, k)) / dy &
                      &+ Pstrat(k) / rhov00 * (var%pi(i, j + 1, k) - var%pi(i, &
                      &j, k)) / dy + Pstrat(k) / rhov1m * (var%pi(i + 1, j, k) &
                      &- var%pi(i + 1, j - 1, k)) / dy + Pstrat(k) / rhov10 &
                      &* (var%pi(i + 1, j + 1, k) - var%pi(i + 1, j, k)) / dy)
                end if
              end if

              ! gravity-wave forcing
              if(raytracer .or. testCase == "mountainwave") then
                volfcx = 0.5 * (force(i, j, k, 1) + force(i + 1, j, k, 1))
                volfcy = 0.5 * (force(i, j, k, 2) + force(i, j + 1, k, 2))
              else
                volfcx = 0.0
                volfcy = 0.0
              end if

              ! ustar
              if(TestCase == "baroclinic_LC") then
                uhorx = var%u(i, j, k) - var_env%u(i, j, k)
              else
                uhorx = var%u(i, j, k)
              end if

              vhory = 0.25 * (var%v(i, j - 1, k) + var%v(i, j, k) + var%v(i &
                  &+ 1, j - 1, k) + var%v(i + 1, j, k))

              facu = 1.0

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  facu = facu + dt * kv_hs(j, k)
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                if(topography) then
                  facu = facu + dt * 0.5 * (kr_sp_tfc(i, j, k) + kr_sp_tfc(i &
                      &+ 1, j, k))
                else
                  facu = facu + dt * kr_sp(j, k)
                end if
              end if

              facv = facu

              if(topography) then
                ! Coriolis force is integrated on LHS.
                if(model == "compressible") then ! Muliply with JP
                  JPR = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i + 1, j, &
                      &k) * var%P(i + 1, j, k))
                  uAst = 1.0 / facu * (uhorx + dt * (- piGradX + volfcx &
                      &/ rhou) * JPR)
                else
                  uAst = 1.0 / facu * (uhorx + dt * (- piGradX + volfcx / rhou))
                end if
              else
                if(testCase == "SkamarockKlemp94") then
                  uAst = 1.0 / (facu * facv + (f_cor_nd(j) * dt) ** 2) * (facv &
                      &* (uhorx + dt * (volfcx / rhou - piGradx)) &
                      &+ f_cor_nd(j) * dt * (vhory + dt * (volfcy / rhou &
                      &- piGrady)) + f_cor_nd(j) ** 2 * dt ** 2 &
                      &* backgroundFlow_dim(1) / uRef)
                else
                  uAst = 1.0 / (facu * facv + (f_cor_nd(j) * dt) ** 2) * (facv &
                      &* (uhorx + dt * (volfcx / rhou - piGradx)) &
                      &+ f_cor_nd(j) * dt * (vhory + dt * (volfcy / rhou &
                      &- piGrady)))
                end if
              end if

              usave(i, j, k) = uAst

              if(TestCase == "baroclinic_LC") then
                usave(i, j, k) = usave(i, j, k) + var_env%u(i, j, k)
              end if
            end do
          end do
        end do
      else
        stop 'ERROR: unknown int_mod'
      end if
    else
      stop 'ERROR: unknown mmp_mod'
    end if

    ! -------------------------------------
    !            predict v -> v*
    ! -------------------------------------

    select case(yBoundary)

    case("solid_wall")
      j0 = 1
      j1 = ny - 1
    case("periodic")
      j0 = 0
      j1 = ny
    case default
      stop "momentumPredictor: unknown case yBoundary."
    end select

    if(mmp_mod == "tot" .or. mmp_mod == "lhs") then
      if(int_mod /= "expl") then
        stop 'ERROR: wrong int_mod for mmp_mod = tot or mmp_mod = lhs'
      end if

      do k = 1, nz
        do j = j0, j1
          do i = 1, nx

            !--- convective part -> conv
            fR = flux%v(i, j, k, 1)
            fL = flux%v(i - 1, j, k, 1)
            gF = flux%v(i, j, k, 2)
            gB = flux%v(i, j - 1, k, 2)
            hU = flux%v(i, j, k, 3)
            hD = flux%v(i, j, k - 1, 3)
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

            ! Adjust meridional momentum flux divergence.
            if(topography) then
              jacEdgeF = 0.5 * (jac(i, j, k) + jac(i, j + 1, k))
              fluxDiff = fluxDiff / jacEdgeF
            end if

            volForce = 0.

            if(mmp_mod == "tot") then
              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piF = var%pi(i, j + 1, k) - var_env%pi(i, j + 1, k)
                piB = var%pi(i, j, k) - var_env%pi(i, j, k)
              else
                piF = var%pi(i, j + 1, k)
                piB = var%pi(i, j, k)
              end if

              if(topography) then
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) - var_env%pi(:, :, :)
                end if
                ! Compute values at cell edges.
                pEdgeF = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j + 1, k))
                met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3))
                ! Compute pressure gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  piUUEdgeF = 0.5 * (var%pi(i, j, k + 2) + var%pi(i, j + 1, k &
                      &+ 2))
                  piUEdgeF = 0.5 * (var%pi(i, j, k + 1) + var%pi(i, j + 1, k &
                      &+ 1))
                  piEdgeF = 0.5 * (var%pi(i, j, k) + var%pi(i, j + 1, k))
                  piGrad = kappaInv * MaInv2 * pEdgeF * ((piF - piB) / dy &
                      &+ met23EdgeF * (- piUUEdgeF + 4.0 * piUEdgeF - 3.0 &
                      &* piEdgeF) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeF = 0.5 * (var%pi(i, j, k - 2) + var%pi(i, j + 1, k &
                      &- 2))
                  piDEdgeF = 0.5 * (var%pi(i, j, k - 1) + var%pi(i, j + 1, k &
                      &- 1))
                  piEdgeF = 0.5 * (var%pi(i, j, k) + var%pi(i, j + 1, k))
                  piGrad = kappaInv * MaInv2 * pEdgeF * ((piF - piB) / dy &
                      &+ met23EdgeF * (piDDEdgeF - 4.0 * piDEdgeF + 3.0 &
                      &* piEdgeF) * 0.5 / dz)
                else
                  piUEdgeF = 0.5 * (var%pi(i, j, k + 1) + var%pi(i, j + 1, k &
                      &+ 1))
                  piDEdgeF = 0.5 * (var%pi(i, j, k - 1) + var%pi(i, j + 1, k &
                      &- 1))
                  piGrad = kappaInv * MaInv2 * pEdgeF * ((piF - piB) / dy &
                      &+ met23EdgeF * (piUEdgeF - piDEdgeF) * 0.5 / dz)
                end if
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) + var_env%pi(:, :, :)
                end if
              else
                piGrad = kappaInv * MaInv2 * pStrat(k) * (piF - piB) / dy

                if(TestCase == "baroclinic_LC") then !FS
                  piGrad = piGrad + kappaInv * MaInv2 * (PStrat(k) &
                      &- pStrat_0(k)) * (var_env%pi(i, j + 1, k) &
                      &- var_env%pi(i, j, k)) / dy
                end if
              end if

              !---- volume forces
              volForce = 0.5 * (force(i, j, k, 2) + force(i, j + 1, k, 2))

              if(TestCase == "baroclinic_LC") then
                select case(model)
                case("pseudo_incompressible", "compressible")
                  rhoM = rhoOld(i, j, k)
                  rhoM_1 = rhoOld(i, j + 1, k)

                  if(topography) then
                    rhoM = rhoM + rhoStratTFC(i, j, k)
                    rhoM_1 = rhoM_1 + rhoStratTFC(i, j + 1, k)
                  else
                    rhoM = rhoM + rhoStrat(k)
                    rhoM_1 = rhoM_1 + rhoStrat(k)
                  end if
                case("Boussinesq")
                  rhoM = rho00
                  rhoM_1 = rho00
                case default
                  stop "momentumPredictor: unkown model."
                end select

                volforce = volforce + RoInv(j) * 0.25 * (rhoM * (var_env%u(i &
                    &- 1, j, k) + var_env%u(i, j, k)) + rhoM_1 * (var_env%u(i &
                    &- 1, j + 1, k) + var_env%u(i, j + 1, k)))

              end if

            end if

            if(TestCase == "baroclinic_LC") then
              if(background == "HeldSuarez") then
                ! Rayleigh damping

                select case(model)
                case("pseudo_incompressible", "compressible")
                  rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j + 1, k))

                  if(topography) then
                    rhoM_1 = rhoM_1 + 0.5 * (rhoStratTFC(i, j, k) &
                        &+ rhoStratTFC(i, j + 1, k))
                  else
                    rhoM_1 = rhoM_1 + rhoStrat(k)
                  end if
                case("Boussinesq")
                  rhoM_1 = rho00
                case default
                  stop "momentumPredictor: unkown model."
                end select

                volForce = volForce - 0.5 * (kv_hs(j, k) + kv_hs(j + 1, k)) &
                    &* rhoM_1 * (var%v(i, j, k) - var_env%v(i, j, k))

              end if
            end if

            ! Explicit integration of Coriolis force in TFC.
            if(topography .and. mmp_mod == "lhs") then
              vOldTFC(i, j, k) = var%v(i, j, k)
              uC = 0.5 * (uOldTFC(i, j, k) + uOldTFC(i - 1, j, k))
              uF = 0.5 * (uOldTFC(i, j + 1, k) + uOldTFC(i - 1, j + 1, k))
              if(testCase == "SkamarockKlemp94") then
                uC = uC - backgroundFlow_dim(1) / uRef
                uF = uF - backgroundFlow_dim(1) / uRef
              else if(testCase == "baroclinic_LC") then
                uC = uC - 0.5 * (var_env%u(i, j, k) + var_env%u(i - 1, j, k))
                uF = uF - 0.5 * (var_env%u(i, j + 1, k) + var_env%u(i - 1, j &
                    &+ 1, k))
              end if
              volForce = volForce - 0.5 * (f_cor_nd(j) * (rhoOld(i, j, k) &
                  &+ rhoStratTFC(i, j, k)) * uC + f_cor_nd(j + 1) * (rhoOld(i, &
                  &j + 1, k) + rhoStratTFC(i, j + 1, k)) * uF)
            end if

            !--------------------
            !   F(phi) = RHS
            !--------------------
            ! fluxDiff -> convective and viscous fluxes
            ! piGrad   -> pressure gradient along x
            ! volForce -> Gravity, Coriolis
            if(mmp_mod == "tot") then
              F = - fluxDiff - piGrad + volForce
            else if(mmp_mod == "lhs") then
              F = - fluxDiff + volForce !UA 200413
            else
              stop 'ERROR: wrong mmp_mod'
            end if

            ! interpolated density
            select case(model)

            case("pseudo_incompressible", "compressible")

              rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j + 1, k))
              rhoM = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k))

              if(topography) then
                rhoStratEdgeF = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, j &
                    &+ 1, k))
                rhoM_1 = rhoM_1 + rhoStratEdgeF
                rhoM = rhoM + rhoStratEdgeF
              else
                rhoM_1 = rhoM_1 + rhoStrat(k)
                rhoM = rhoM + rhoStrat(k)
              end if

            case("Boussinesq")
              rhoM_1 = rho00
              rhoM = rho00
            case default
              stop "momentumPredictor: unkown case model."
            end select

            ! velocity and momentum at t(m-1)
            vM_1 = var%v(i, j, k)
            momM_1 = rhoM_1 * vM_1

            ! q(m-1) -> q(m)
            q(i, j, k, 2) = dt * F + alphaRK(m) * q(i, j, k, 2)

            ! rhoV(m-1) -> rhoV(m)
            momM = momM_1 + betaRK(m) * q(i, j, k, 2)

            ! calc v(m,*)
            vAst = momM / rhoM

            ! vAst -> var
            var%v(i, j, k) = vAst
          end do
        end do
      end do
    else if(mmp_mod == "rhs") then
      if(int_mod == "expl") then
        do k = 1, nz
          do j = j0, j1
            do i = 1, nx
              rhov = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k))
              if(topography) then
                rhoStratEdgeF = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, j &
                    &+ 1, k))
                rhov = rhov + rhoStratEdgeF
              else
                rhov = rhov + rhoStrat(k)
              end if

              if(TestCase == "baroclinic_LC") then
                rhov_e = 0.5 * (var_env%rho(i, j, k) + var_env%rho(i, j + 1, k))

                if(topography) then
                  rhov_e = rhov_e + 0.5 * (rhoStratTFC(i, j, k) &
                      &+ rhoStratTFC(i, j + 1, k))
                else
                  rhov_e = rhov_e + rhoStrat_0(k)
                end if
              end if

              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piF = var%pi(i, j + 1, k) - var_env%pi(i, j + 1, k)
                piB = var%pi(i, j, k) - var_env%pi(i, j, k)

                if(TestCase == "baroclinic_LC") then !FS
                  piF_e = var_env%pi(i, j + 1, k)
                  piB_e = var_env%pi(i, j, k)
                end if
              else
                piF = var%pi(i, j + 1, k)
                piB = var%pi(i, j, k)
              end if

              if(topography) then
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) - var_env%pi(:, :, :)
                end if
                ! Compute values at cell edges.
                pEdgeF = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j + 1, k))
                met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3))
                ! Compute pressure gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  piUUEdgeF = 0.5 * (var%pi(i, j, k + 2) + var%pi(i, j + 1, k &
                      &+ 2))
                  piUEdgeF = 0.5 * (var%pi(i, j, k + 1) + var%pi(i, j + 1, k &
                      &+ 1))
                  piEdgeF = 0.5 * (var%pi(i, j, k) + var%pi(i, j + 1, k))
                  piGrad = kappaInv * MaInv2 * pEdgeF / rhov * ((piF - piB) &
                      &/ dy + met23EdgeF * (- piUUEdgeF + 4.0 * piUEdgeF - 3.0 &
                      &* piEdgeF) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeF = 0.5 * (var%pi(i, j, k - 2) + var%pi(i, j + 1, k &
                      &- 2))
                  piDEdgeF = 0.5 * (var%pi(i, j, k - 1) + var%pi(i, j + 1, k &
                      &- 1))
                  piEdgeF = 0.5 * (var%pi(i, j, k) + var%pi(i, j + 1, k))
                  piGrad = kappaInv * MaInv2 * pEdgeF / rhov * ((piF - piB) &
                      &/ dy + met23EdgeF * (piDDEdgeF - 4.0 * piDEdgeF + 3.0 &
                      &* piEdgeF) * 0.5 / dz)
                else
                  piUEdgeF = 0.5 * (var%pi(i, j, k + 1) + var%pi(i, j + 1, k &
                      &+ 1))
                  piDEdgeF = 0.5 * (var%pi(i, j, k - 1) + var%pi(i, j + 1, k &
                      &- 1))
                  piGrad = kappaInv * MaInv2 * pEdgeF / rhov * ((piF - piB) &
                      &/ dy + met23EdgeF * (piUEdgeF - piDEdgeF) * 0.5 / dz)
                end if
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) + var_env%pi(:, :, :)
                end if
              else
                piGrad = kappaInv * MaInv2 * Pstrat(k) / rhov * (piF - piB) / dy
              end if

              if(TestCase == "baroclinic_LC") then !FS
                if(topography) then
                  ! Compute values at cell edges.
                  pEdgeF = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j + 1, k))
                  met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, &
                      &3))
                  ! Compute pressure gradient component.
                  if(k == 1 .and. zBoundary == "solid_wall") then
                    piUUEdgeF = 0.5 * (var_env%pi(i, j, k + 2) + var_env%pi(i, &
                        &j + 1, k + 2))
                    piUEdgeF = 0.5 * (var_env%pi(i, j, k + 1) + var_env%pi(i, &
                        &j + 1, k + 1))
                    piEdgeF = 0.5 * (var_env%pi(i, j, k) + var_env%pi(i, j &
                        &+ 1, k))
                    piGrad = piGrad + kappaInv * MaInv2 * (pEdgeF / rhov &
                        &- pEdgeF / rhov_e) * ((piF_e - piB_e) / dy &
                        &+ met23EdgeF * (- piUUEdgeF + 4.0 * piUEdgeF - 3.0 &
                        &* piEdgeF) * 0.5 / dz)
                  else if(k == nz .and. zBoundary == "solid_wall") then
                    piDDEdgeF = 0.5 * (var_env%pi(i, j, k - 2) + var_env%pi(i, &
                        &j + 1, k - 2))
                    piDEdgeF = 0.5 * (var_env%pi(i, j, k - 1) + var_env%pi(i, &
                        &j + 1, k - 1))
                    piEdgeF = 0.5 * (var_env%pi(i, j, k) + var_env%pi(i, j &
                        &+ 1, k))
                    piGrad = piGrad + kappaInv * MaInv2 * (pEdgeF / rhov &
                        &- pEdgeF / rhov_e) * ((piF_e - piB_e) / dy &
                        &+ met23EdgeF * (piDDEdgeF - 4.0 * piDEdgeF + 3.0 &
                        &* piEdgeF) * 0.5 / dz)
                  else
                    piUEdgeF = 0.5 * (var_env%pi(i, j, k + 1) + var_env%pi(i, &
                        &j + 1, k + 1))
                    piDEdgeF = 0.5 * (var_env%pi(i, j, k - 1) + var_env%pi(i, &
                        &j + 1, k - 1))
                    piGrad = piGrad + kappaInv * MaInv2 * (pEdgeF / rhov &
                        &- pEdgeF / rhov_e) * ((piF_e - piB_e) / dy &
                        &+ met23EdgeF * (piUEdgeF - piDEdgeF) * 0.5 / dz)
                  end if
                else
                  piGrad = piGrad + kappaInv * MaInv2 * (Pstrat(k) / rhov &
                      &- Pstrat_0(k) / rhov_e) * (piF_e - piB_e) / dy
                end if
              end if

              ! gravity-wave forcing
              if(raytracer .or. testCase == "mountainwave") then
                volfcy = 0.5 * (force(i, j, k, 2) + force(i, j + 1, k, 2))
              else
                volfcy = 0.0
              end if

              ! vstar
              if(TestCase == "baroclinic_LC") then
                uhorx = 0.25 * (var%u(i - 1, j, k) + var%u(i - 1, j + 1, k) &
                    &- var_env%u(i - 1, j, k) - var_env%u(i - 1, j + 1, k) &
                    &+ var%u(i, j, k) + var%u(i, j + 1, k) - var_env%u(i, j, &
                    &k) - var_env%u(i, j + 1, k))
              else
                uhorx = 0.25 * (var%u(i - 1, j, k) + var%u(i - 1, j + 1, k) &
                    &+ var%u(i, j, k) + var%u(i, j + 1, k))
              end if

              vhory = var%v(i, j, k)

              f_cor_v = 0.5 * (f_cor_nd(j) + f_cor_nd(j + 1))

              if(topography) then
                ! Coriolis force is integrated on LHS.
                if(model == "compressible") then ! Muliply with JP
                  JPF = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i, j + 1, &
                      &k) * var%P(i, j + 1, k))
                  vAst = vhory + dt * (- piGrad + volfcy / rhov) * JPF
                else
                  vAst = vhory + dt * (- piGrad + volfcy / rhov)
                end if
              else
                if(testCase == "SkamarockKlemp94") then
                  vAst = vhory + dt * (- f_cor_v * (uhorx &
                      &- backgroundFlow_dim(1) / uRef) - piGrad + volfcy / rhov)
                else
                  vAst = vhory + dt * (- f_cor_v * uhorx - piGrad + volfcy &
                      &/ rhov)
                end if
              end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  vAst = vAst - dt * 0.5 * (kv_hs(j, k) + kv_hs(j + 1, k)) &
                      &* vhory
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                if(topography) then
                  vAst = vAst - dt * 0.5 * (kr_sp_tfc(i, j, k) + kr_sp_tfc(i, &
                      &j + 1, k)) * vhory
                else
                  vAst = vAst - dt * 0.5 * (kr_sp(j, k) + kr_sp(j + 1, k)) &
                      &* vhory
                end if
              end if

              var%v(i, j, k) = vAst
            end do
          end do
        end do
      else if(int_mod == "impl") then
        do k = 1, nz
          do j = j0, j1
            do i = 1, nx
              rhoum0 = 0.5 * (var%rho(i, j, k) + var%rho(i - 1, j, k))
              rhou00 = 0.5 * (var%rho(i + 1, j, k) + var%rho(i, j, k))
              rhoum1 = 0.5 * (var%rho(i, j + 1, k) + var%rho(i - 1, j + 1, k))
              rhou01 = 0.5 * (var%rho(i + 1, j + 1, k) + var%rho(i, j + 1, k))

              rhov = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k))

              if(topography) then
                rhoStratEdgeF = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, j &
                    &+ 1, k))
                rhov = rhov + rhoStratEdgeF
              else
                rhov = rhov + rhoStrat(k)
                rhoum0 = rhoum0 + rhoStrat(k)
                rhou00 = rhou00 + rhoStrat(k)
                rhoum1 = rhoum1 + rhoStrat(k)
                rhou01 = rhou01 + rhoStrat(k)
              end if

              if(TestCase == "baroclinic_LC") then
                rhoum0_e = 0.5 * (var_env%rho(i, j, k) + var_env%rho(i - 1, j, &
                    &k))
                rhou00_e = 0.5 * (var_env%rho(i + 1, j, k) + var_env%rho(i, j, &
                    &k))
                rhoum1_e = 0.5 * (var_env%rho(i, j + 1, k) + var_env%rho(i &
                    &- 1, j + 1, k))
                rhou01_e = 0.5 * (var_env%rho(i + 1, j + 1, k) &
                    &+ var_env%rho(i, j + 1, k))

                rhov_e = 0.5 * (var_env%rho(i, j, k) + var_env%rho(i, j + 1, k))

                if(topography) then
                  rhov_e = rhov_e + 0.5 * (rhoStratTFC(i, j, k) &
                      &+ rhoStratTFC(i, j + 1, k))
                else
                  rhov_e = rhov_e + rhoStrat_0(k)
                  rhoum0_e = rhoum0_e + rhoStrat_0(k)
                  rhou00_e = rhou00_e + rhoStrat_0(k)
                  rhoum1_e = rhoum1_e + rhoStrat_0(k)
                  rhou01_e = rhou01_e + rhoStrat_0(k)
                end if
              end if

              !--- pressure gradient terms -> piGradx, piGrady
              if(.not. topography) then
                if(TestCase == "baroclinic_LC") then
                  piGradx = kappaInv * MaInv2 * 0.25 * (Pstrat(k) / rhou00 &
                      &* (var%pi(i + 1, j, k) - var%pi(i, j, k) - var_env%pi(i &
                      &+ 1, j, k) + var_env%pi(i, j, k)) / dx + Pstrat(k) &
                      &/ rhoum0 * (var%pi(i, j, k) - var%pi(i - 1, j, k) &
                      &- var_env%pi(i, j, k) + var_env%pi(i - 1, j, k)) / dx &
                      &+ Pstrat(k) / rhou01 * (var%pi(i + 1, j + 1, k) &
                      &- var%pi(i, j + 1, k) - var_env%pi(i + 1, j + 1, k) &
                      &+ var_env%pi(i, j + 1, k)) / dx + Pstrat(k) / rhoum1 &
                      &* (var%pi(i, j + 1, k) - var%pi(i - 1, j + 1, k) &
                      &- var_env%pi(i, j + 1, k) + var_env%pi(i - 1, j + 1, &
                      &k)) / dx)

                  piGradx = piGradx + kappaInv * MaInv2 * 0.25 * ((Pstrat(k) &
                      &/ rhou00 - Pstrat_0(k) / rhou00_e) * (var_env%pi(i + 1, &
                      &j, k) - var_env%pi(i, j, k)) / dx + (Pstrat(k) / rhoum0 &
                      &- Pstrat_0(k) / rhoum0_e) * (var_env%pi(i, j, k) &
                      &- var_env%pi(i - 1, j, k)) / dx + (Pstrat(k) / rhou01 &
                      &- Pstrat_0(k) / rhou01_e) * (var_env%pi(i + 1, j + 1, &
                      &k) - var_env%pi(i, j + 1, k)) / dx + (Pstrat(k) &
                      &/ rhoum1 - Pstrat_0(k) / rhoum1_e) * (var_env%pi(i, j &
                      &+ 1, k) - var_env%pi(i - 1, j + 1, k)) / dx)
                else
                  piGradx = kappaInv * MaInv2 * 0.25 * (Pstrat(k) / rhou00 &
                      &* (var%pi(i + 1, j, k) - var%pi(i, j, k)) / dx &
                      &+ Pstrat(k) / rhoum0 * (var%pi(i, j, k) - var%pi(i - 1, &
                      &j, k)) / dx + Pstrat(k) / rhou01 * (var%pi(i + 1, j &
                      &+ 1, k) - var%pi(i, j + 1, k)) / dx + Pstrat(k) &
                      &/ rhoum1 * (var%pi(i, j + 1, k) - var%pi(i - 1, j + 1, &
                      &k)) / dx)
                end if
              end if

              if(TestCase == "baroclinic_LC") then
                piF = var%pi(i, j + 1, k) - var_env%pi(i, j + 1, k)
                piB = var%pi(i, j, k) - var_env%pi(i, j, k)

                piF_e = var_env%pi(i, j + 1, k)
                piB_e = var_env%pi(i, j, k)
              else
                piF = var%pi(i, j + 1, k)
                piB = var%pi(i, j, k)
              end if

              if(topography) then
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) - var_env%pi(:, :, :)
                end if
                ! Compute values at cell edges.
                pEdgeF = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j + 1, k))
                met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3))
                ! Compute pressure gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  piUUEdgeF = 0.5 * (var%pi(i, j, k + 2) + var%pi(i, j + 1, k &
                      &+ 2))
                  piUEdgeF = 0.5 * (var%pi(i, j, k + 1) + var%pi(i, j + 1, k &
                      &+ 1))
                  piEdgeF = 0.5 * (var%pi(i, j, k) + var%pi(i, j + 1, k))
                  piGradY = kappaInv * MaInv2 * pEdgeF / rhov * ((piF - piB) &
                      &/ dy + met23EdgeF * (- piUUEdgeF + 4.0 * piUEdgeF - 3.0 &
                      &* piEdgeF) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeF = 0.5 * (var%pi(i, j, k - 2) + var%pi(i, j + 1, k &
                      &- 2))
                  piDEdgeF = 0.5 * (var%pi(i, j, k - 1) + var%pi(i, j + 1, k &
                      &- 1))
                  piEdgeF = 0.5 * (var%pi(i, j, k) + var%pi(i, j + 1, k))
                  piGradY = kappaInv * MaInv2 * pEdgeF / rhov * ((piF - piB) &
                      &/ dy + met23EdgeF * (piDDEdgeF - 4.0 * piDEdgeF + 3.0 &
                      &* piEdgeF) * 0.5 / dz)
                else
                  piUEdgeF = 0.5 * (var%pi(i, j, k + 1) + var%pi(i, j + 1, k &
                      &+ 1))
                  piDEdgeF = 0.5 * (var%pi(i, j, k - 1) + var%pi(i, j + 1, k &
                      &- 1))
                  piGradY = kappaInv * MaInv2 * pEdgeF / rhov * ((piF - piB) &
                      &/ dy + met23EdgeF * (piUEdgeF - piDEdgeF) * 0.5 / dz)
                end if
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) + var_env%pi(:, :, :)
                end if
              else
                piGrady = kappaInv * MaInv2 * Pstrat(k) / rhov * (piF - piB) &
                    &/ dy
              end if

              if(TestCase == "baroclinic_LC") then !FS
                if(topography) then
                  ! Compute values at cell edges.
                  pEdgeF = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j + 1, k))
                  met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, &
                      &3))
                  ! Compute pressure gradient component.
                  if(k == 1 .and. zBoundary == "solid_wall") then
                    piUUEdgeF = 0.5 * (var_env%pi(i, j, k + 2) + var_env%pi(i, &
                        &j + 1, k + 2))
                    piUEdgeF = 0.5 * (var_env%pi(i, j, k + 1) + var_env%pi(i, &
                        &j + 1, k + 1))
                    piEdgeF = 0.5 * (var_env%pi(i, j, k) + var_env%pi(i, j &
                        &+ 1, k))
                    piGradY = piGradY + kappaInv * MaInv2 * (pEdgeF / rhov &
                        &- pEdgeF / rhov_e) * ((piF_e - piB_e) / dy &
                        &+ met23EdgeF * (- piUUEdgeF + 4.0 * piUEdgeF - 3.0 &
                        &* piEdgeF) * 0.5 / dz)
                  else if(k == nz .and. zBoundary == "solid_wall") then
                    piDDEdgeF = 0.5 * (var_env%pi(i, j, k - 2) + var_env%pi(i, &
                        &j + 1, k - 2))
                    piDEdgeF = 0.5 * (var_env%pi(i, j, k - 1) + var_env%pi(i, &
                        &j + 1, k - 1))
                    piEdgeF = 0.5 * (var_env%pi(i, j, k) + var_env%pi(i, j &
                        &+ 1, k))
                    piGradY = piGradY + kappaInv * MaInv2 * (pEdgeF / rhov &
                        &- pEdgeF / rhov_e) * ((piF_e - piB_e) / dy &
                        &+ met23EdgeF * (piDDEdgeF - 4.0 * piDEdgeF + 3.0 &
                        &* piEdgeF) * 0.5 / dz)
                  else
                    piUEdgeF = 0.5 * (var_env%pi(i, j, k + 1) + var_env%pi(i, &
                        &j + 1, k + 1))
                    piDEdgeF = 0.5 * (var_env%pi(i, j, k - 1) + var_env%pi(i, &
                        &j + 1, k - 1))
                    piGradY = piGradY + kappaInv * MaInv2 * (pEdgeF / rhov &
                        &- pEdgeF / rhov_e) * ((piF_e - piB_e) / dy &
                        &+ met23EdgeF * (piUEdgeF - piDEdgeF) * 0.5 / dz)
                  end if
                else
                  piGrady = piGrady + kappaInv * MaInv2 * (Pstrat(k) / rhov &
                      &- Pstrat_0(k) / rhov_e) * (piF_e - piB_e) / dy
                end if
              end if

              ! gravity-wave forcing
              if(raytracer .or. testCase == "mountainwave") then
                volfcx = 0.5 * (force(i, j, k, 1) + force(i + 1, j, k, 1))
                volfcy = 0.5 * (force(i, j, k, 2) + force(i, j + 1, k, 2))
              else
                volfcx = 0.0
                volfcy = 0.0
              end if

              ! vstar
              if(TestCase == "baroclinic_LC") then
                uhorx = 0.25 * (var%u(i - 1, j, k) + var%u(i - 1, j + 1, k) &
                    &- var_env%u(i - 1, j, k) - var_env%u(i - 1, j + 1, k) &
                    &+ var%u(i, j, k) + var%u(i, j + 1, k) - var_env%u(i, j, &
                    &k) - var_env%u(i, j + 1, k))
              else
                uhorx = 0.25 * (var%u(i - 1, j, k) + var%u(i - 1, j + 1, k) &
                    &+ var%u(i, j, k) + var%u(i, j + 1, k))
              end if

              vhory = var%v(i, j, k)

              facv = 1.0

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  facv = facv + dt * 0.5 * (kv_hs(j, k) + kv_hs(j + 1, k))
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                if(topography) then
                  facv = facv + dt * 0.5 * (kr_sp_tfc(i, j, k) + kr_sp_tfc(i, &
                      &j + 1, k))
                else
                  facv = facv + dt * 0.5 * (kr_sp(j, k) + kr_sp(j + 1, k))
                end if
              end if

              facu = facv

              if(topography) then
                ! Coriolis force is integrated on LHS.
                if(model == "compressible") then ! Muliply with JP
                  JPF = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i, j + 1, &
                      &k) * var%P(i, j + 1, k))
                  vAst = 1.0 / facv * (vhory + dt * (- piGradY + volfcy &
                      &/ rhov) * JPF)
                else
                  vAst = 1.0 / facv * (vhory + dt * (- piGradY + volfcy / rhov))
                end if
              else
                if(testCase == "SkamarockKlemp94") then
                  vAst = 1.0 / (facu * facv + (0.5 * (f_cor_nd(j) + f_cor_nd(j &
                      &+ 1)) * dt) ** 2) * (- 0.5 * (f_cor_nd(j) + f_cor_nd(j &
                      &+ 1)) * dt * ((uhorx - backgroundFlow_dim(1) / uRef) &
                      &+ dt * (volfcx / rhov - piGradx)) + facu * (vhory + dt &
                      &* (volfcy / rhov - piGrady)))
                else
                  vAst = 1.0 / (facu * facv + (0.5 * (f_cor_nd(j) + f_cor_nd(j &
                      &+ 1)) * dt) ** 2) * (- 0.5 * (f_cor_nd(j) + f_cor_nd(j &
                      &+ 1)) * dt * (uhorx + dt * (volfcx / rhov - piGradx)) &
                      &+ facu * (vhory + dt * (volfcy / rhov - piGrady)))
                end if
              end if

              var%v(i, j, k) = vAst
            end do
          end do
        end do
      else
        stop 'ERROR: unknown int_mod'
      end if

      ! now the new u can be put into the proper array
      var%u(:, :, :) = usave(:, :, :)
    else
      stop 'ERROR: unknown mmp_mod'
    end if

    !testb
    !write(42) var
    !stop
    !teste

    ! -------------------------------------
    !            predict w -> w*
    ! -------------------------------------

    select case(zBoundary)

    case("solid_wall")
      k0 = 1
      k1 = nz - 1
    case("periodic")
      k0 = 0
      k1 = nz
    case default
      stop "momentumPredictor: unknown case zBoundary."
    end select

    if(mmp_mod == "tot" .or. mmp_mod == "lhs") then
      if(int_mod /= "expl") then
        stop 'ERROR: wrong int_mod for mmp_mod = tot or mmp_mod = lhs'
      end if

      do k = k0, k1
        do j = 1, ny
          do i = 1, nx

            !--- convective part -> conv
            fR = flux%w(i, j, k, 1)
            fL = flux%w(i - 1, j, k, 1)
            gF = flux%w(i, j, k, 2)
            gB = flux%w(i, j - 1, k, 2)
            hU = flux%w(i, j, k, 3)
            hD = flux%w(i, j, k - 1, 3)
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

            ! Adjust vertical momentum flux divergence.
            if(topography) then
              ! Adjust Cartesian vertical momentum flux divergence.
              jacEdgeU = 2.0 * jac(i, j, k) * jac(i, j, k + 1) / (jac(i, j, k) &
                  &+ jac(i, j, k + 1))
              fluxDiff = fluxDiff / jacEdgeU
              ! Compute zonal momentum flux divergences.
              do ll = 0, 1
                do mm = 0, 1
                  fR = flux%u(i - ll, j, k + mm, 1)
                  fL = flux%u(i - 1 - ll, j, k + mm, 1)
                  gF = flux%u(i - ll, j, k + mm, 2)
                  gB = flux%u(i - ll, j - 1, k + mm, 2)
                  hU = flux%u(i - ll, j, k + mm, 3)
                  hD = flux%u(i - ll, j, k - 1 + mm, 3)
                  fluxDiffU(ll, mm) = (fR - fL) / dx + (gF - gB) / dy + (hU &
                      &- hD) / dz
                  jacEdgeR = 0.5 * (jac(i - ll, j, k + mm) + jac(i + 1 - ll, &
                      &j, k + mm))
                  fluxDiffU(ll, mm) = fluxDiffU(ll, mm) / jacEdgeR
                end do
              end do
              ! Compute meridional momentum flux divergences.
              do ll = 0, 1
                do mm = 0, 1
                  fR = flux%v(i, j - ll, k + mm, 1)
                  fL = flux%v(i - 1, j - ll, k + mm, 1)
                  gF = flux%v(i, j - ll, k + mm, 2)
                  gB = flux%v(i, j - 1 - ll, k + mm, 2)
                  hU = flux%v(i, j - ll, k + mm, 3)
                  hD = flux%v(i, j - ll, k - 1 + mm, 3)
                  fluxDiffV(ll, mm) = (fR - fL) / dx + (gF - gB) / dy + (hU &
                      &- hD) / dz
                  jacEdgeF = 0.5 * (jac(i, j - ll, k + mm) + jac(i, j + 1 &
                      &- ll, k + mm))
                  fluxDiffV(ll, mm) = fluxDiffV(ll, mm) / jacEdgeF
                end do
              end do
              ! Compute transformed vertical momentum flux divergence.
              fluxDiff = trafoTFC(i, j, k, fluxDiffU(0, 0), fluxDiffU(0, 1), &
                  &fluxDiffU(1, 0), fluxDiffU(1, 1), fluxDiffV(0, 0), &
                  &fluxDiffV(0, 1), fluxDiffV(1, 0), fluxDiffV(1, 1), &
                  &fluxDiff, "tfc")
            end if

            if(mmp_mod == "tot") then
              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piU = var%pi(i, j, k + 1) - var_env%pi(i, j, k + 1)
                piD = var%pi(i, j, k) - var_env%pi(i, j, k)
              else
                piU = var%pi(i, j, k + 1)
                piD = var%pi(i, j, k)
              end if

              if(topography) then
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) - var_env%pi(:, :, :)
                end if
                ! Compute values at cell edges.
                pEdgeU = (jac(i, j, k + 1) * pStratTFC(i, j, k) + jac(i, j, k) &
                    &* pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k &
                    &+ 1))
                met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) + jac(i, &
                    &j, k) * met(i, j, k + 1, 1, 3)) / (jac(i, j, k) + jac(i, &
                    &j, k + 1))
                met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) + jac(i, &
                    &j, k) * met(i, j, k + 1, 2, 3)) / (jac(i, j, k) + jac(i, &
                    &j, k + 1))
                met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) + jac(i, &
                    &j, k) * met(i, j, k + 1, 3, 3)) / (jac(i, j, k) + jac(i, &
                    &j, k + 1))
                piREdgeU = (jac(i + 1, j, k + 1) * var%pi(i + 1, j, k) + jac(i &
                    &+ 1, j, k) * var%pi(i + 1, j, k + 1)) / (jac(i + 1, j, k) &
                    &+ jac(i + 1, j, k + 1))
                piLEdgeU = (jac(i - 1, j, k + 1) * var%pi(i - 1, j, k) + jac(i &
                    &- 1, j, k) * var%pi(i - 1, j, k + 1)) / (jac(i - 1, j, k) &
                    &+ jac(i - 1, j, k + 1))
                piFEdgeU = (jac(i, j + 1, k + 1) * var%pi(i, j + 1, k) &
                    &+ jac(i, j + 1, k) * var%pi(i, j + 1, k + 1)) / (jac(i, j &
                    &+ 1, k) + jac(i, j + 1, k + 1))
                piBEdgeU = (jac(i, j - 1, k + 1) * var%pi(i, j - 1, k) &
                    &+ jac(i, j - 1, k) * var%pi(i, j - 1, k + 1)) / (jac(i, j &
                    &- 1, k) + jac(i, j - 1, k + 1))
                ! Compute pressure gradient component.
                piGrad = kappaInv * MaInv2 * pEdgeU * (met13EdgeU * (piREdgeU &
                    &- piLEdgeU) * 0.5 / dx + met23EdgeU * (piFEdgeU &
                    &- piBEdgeU) * 0.5 / dy + met33EdgeU * (var%pi(i, j, k &
                    &+ 1) - var%pi(i, j, k)) / dz)
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) + var_env%pi(:, :, :)
                end if
              else
                piGrad = 0.5 * kappaInv * MaInv2 * (Pstrat(k) + Pstrat(k + 1)) &
                    &* (piU - piD) / dz

                if(TestCase == "baroclinic_LC") then !FS
                  piGrad = piGrad + 0.5 * kappaInv * MaInv2 * (Pstrat(k) &
                      &+ Pstrat(k + 1) - pStrat_0(k) - pStrat_0(k + 1)) &
                      &* (var_env%pi(i, j, k + 1) - var_env%pi(i, j, k)) / dz
                end if
              end if

              !---- volume forces
              if(topography) then
                volForce = (jac(i, j, k + 1) * force(i, j, k, 3) + jac(i, j, &
                    &k) * force(i, j, k + 1, 3)) / (jac(i, j, k) + jac(i, j, k &
                    &+ 1))
              else
                volForce = 0.5 * (force(i, j, k, 3) + force(i, j, k + 1, 3))
              end if

              if(TestCase == "baroclinic_LC") then
                select case(model)
                case("pseudo_incompressible", "compressible")
                  if(topography) then
                    drho_e = (jac(i, j, k + 1) * var_env%rho(i, j, k) + jac(i, &
                        &j, k) * var_env%rho(i, j, k + 1)) / (jac(i, j, k) &
                        &+ jac(i, j, k + 1))
                  else
                    drho_e = 0.5 * (var_env%rho(i, j, k) + var_env%rho(i, j, k &
                        &+ 1))
                  end if

                case("Boussinesq")
                  stop 'ERROR: baroclinic LC not ready yet for  Boussinesq'
                case default
                  stop "momentumPredictor: unkown model."
                end select

                if(topography) then
                  volForce = volForce + FrInv2 * drho_e / (2.0 * jac(i, j, k) &
                      &* jac(i, j, k + 1) / (jac(i, j, k) + jac(i, j, k + 1)))
                else
                  volForce = volForce + FrInv2 * drho_e
                end if
              end if

            end if

            if(TestCase == "baroclinic_LC") then
              if(background == "HeldSuarez") then
                ! Rayleigh damping
                select case(model)
                case("pseudo_incompressible", "compressible")
                  if(topography) then
                    rhoM_1 = (jac(i, j, k + 1) * (rhoOld(i, j, k) &
                        &+ rhoStratTFC(i, j, k)) + jac(i, j, k) * (rhoOld(i, &
                        &j, k + 1) + rhoStratTFC(i, j, k + 1))) / (jac(i, j, &
                        &k) + jac(i, j, k + 1))
                  else
                    rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j, k + 1)) &
                        &+ rhoStratTilde(k)
                  end if
                case("Boussinesq")
                  rhoM_1 = rho00
                case default
                  stop "momentumPredictor: unkown model."
                end select

                volForce = volforce - 0.5 * (kw_hs(k) + kw_hs(k + 1)) * rhoM_1 &
                    &* var%w(i, j, k)
              end if
            end if

            ! Explicit integration of Coriolis force in TFC.
            if(topography .and. mmp_mod == "lhs") then
              vC = 0.5 * (vOldTFC(i, j, k) + vOldTFC(i, j - 1, k))
              vU = 0.5 * (vOldTFC(i, j, k + 1) + vOldTFC(i, j - 1, k + 1))
              uC = 0.5 * (uOldTFC(i, j, k) + uOldTFC(i - 1, j, k))
              uU = 0.5 * (uOldTFC(i, j, k + 1) + uOldTFC(i - 1, j, k + 1))
              if(testCase == "SkamarockKlemp94") then
                uC = uC - backgroundFlow_dim(1) / uRef
                uU = uU - backgroundFlow_dim(1) / uRef
              else if(testCase == "baroclinic_LC") then
                vC = vC - 0.5 * (var_env%v(i, j, k) + var_env%v(i, j - 1, k))
                vU = vU - 0.5 * (var_env%v(i, j, k + 1) + var_env%v(i, j - 1, &
                    &k + 1))
                uC = uC - 0.5 * (var_env%u(i, j, k) + var_env%u(i - 1, j, k))
                uU = vU - 0.5 * (var_env%u(i, j, k + 1) + var_env%u(i - 1, j, &
                    &k + 1))
              end if
              volForce = volForce + f_cor_nd(j) * (jac(i, j, k + 1) * met(i, &
                  &j, k, 1, 3) * (rhoOld(i, j, k) + rhoStratTFC(i, j, k)) * vC &
                  &+ jac(i, j, k) * met(i, j, k + 1, 1, 3) * (rhoOld(i, j, k &
                  &+ 1) + rhoStratTFC(i, j, k + 1)) * vU) / (jac(i, j, k) &
                  &+ jac(i, j, k + 1)) - f_cor_nd(j) * (jac(i, j, k + 1) &
                  &* met(i, j, k, 2, 3) * (rhoOld(i, j, k) + rhoStratTFC(i, j, &
                  &k)) * uC + jac(i, j, k) * met(i, j, k + 1, 2, 3) &
                  &* (rhoOld(i, j, k + 1) + rhoStratTFC(i, j, k + 1)) * uU) &
                  &/ (jac(i, j, k) + jac(i, j, k + 1))
            end if

            !--------------------
            !   F(phi) = RHS
            !--------------------
            ! fluxDiff -> convective and viscous fluxes
            ! piGrad   -> pressure gradient along x
            ! volForce -> Gravity, Coriolis
            if(mmp_mod == "tot") then
              F = - fluxDiff - piGrad + volForce
            else if(mmp_mod == "lhs" .and. topography) then
              F = - fluxDiff + volForce
            else if(mmp_mod == "lhs") then
              F = - fluxDiff
            else
              stop 'ERROR: wrong mmp_mod'
            end if

            ! interpolated densities
            select case(model)

            case("pseudo_incompressible", "compressible")
              if(topography) then
                rhoM_1 = (jac(i, j, k + 1) * rhoOld(i, j, k) + jac(i, j, k) &
                    &* rhoOld(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
                rhoM = (jac(i, j, k + 1) * var%rho(i, j, k) + jac(i, j, k) &
                    &* var%rho(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))

                rhoStratEdgeU = (jac(i, j, k + 1) * rhoStratTFC(i, j, k) &
                    &+ jac(i, j, k) * rhoStratTFC(i, j, k + 1)) / (jac(i, j, &
                    &k) + jac(i, j, k + 1))
                rhoM_1 = rhoM_1 + rhoStratEdgeU
                rhoM = rhoM + rhoStratEdgeU
              else
                rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j, k + 1))
                rhoM = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1))

                rhoM_1 = rhoM_1 + rhoStratTilde(k)
                rhoM = rhoM + rhoStratTilde(k)
              end if

            case("Boussinesq")
              rhoM_1 = rho00
              rhoM = rho00
            case default
              stop "momentumPredictor: unkown case model."
            end select

            ! velocity and momentum at t(m-1)
            wM_1 = var%w(i, j, k)
            momM_1 = rhoM_1 * wM_1

            ! q(m-1) -> q(m)
            q(i, j, k, 3) = dt * F + alphaRK(m) * q(i, j, k, 3)

            ! rhoW(m-1) -> rhoW(m)
            momM = momM_1 + betaRK(m) * q(i, j, k, 3)

            ! calc w(m,*)
            wAst = momM / rhoM

            ! wAst -> var
            var%w(i, j, k) = wAst
          end do
        end do
      end do
    else if(mmp_mod == "rhs") then
      if(int_mod == "expl") then
        do k = k0, k1
          pstw = 0.5 * (Pstrat(k) + Pstrat(k + 1))

          if(TestCase == "baroclinic_LC") then
            pstw_e = 0.5 * (Pstrat_0(k) + Pstrat_0(k + 1))
          end if

          do j = 1, ny
            do i = 1, nx
              rho000 = var%rho(i, j, k)
              rho001 = var%rho(i, j, k + 1)

              if(topography) then
                rhow = (jac(i, j, k + 1) * var%rho(i, j, k) + jac(i, j, k) &
                    &* var%rho(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))

                rhoStratEdgeU = (jac(i, j, k + 1) * rhoStratTFC(i, j, k) &
                    &+ jac(i, j, k) * rhoStratTFC(i, j, k + 1)) / (jac(i, j, &
                    &k) + jac(i, j, k + 1))
                rho000 = rho000 + rhoStratTFC(i, j, k)
                rho001 = rho001 + rhoStratTFC(i, j, k + 1)
                rhow = rhow + rhoStratEdgeU
              else
                rhow = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1))

                rho000 = rho000 + rhoStrat(k)
                rho001 = rho001 + rhoStrat(k + 1)

                rhow = rhow + rhoStratTilde(k)
              end if

              if(TestCase == "baroclinic_LC") then
                rho000_e = var_env%rho(i, j, k)
                rho001_e = var_env%rho(i, j, k + 1)

                if(topography) then
                  rhow_e = (jac(i, j, k + 1) * var_env%rho(i, j, k) + jac(i, &
                      &j, k) * var_env%rho(i, j, k + 1)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))

                  rho000_e = rho000_e + rhoStratTFC(i, j, k)
                  rho001_e = rho001_e + rhoStratTFC(i, j, k + 1)
                  rhow_e = rhow_e + (jac(i, j, k + 1) * rhoStratTFC(i, j, k) &
                      &+ jac(i, j, k) * rhoStratTFC(i, j, k + 1)) / (jac(i, j, &
                      &k) + jac(i, j, k + 1))
                else
                  rhow_e = 0.5 * (var_env%rho(i, j, k) + var_env%rho(i, j, k &
                      &+ 1))

                  rho000_e = rho000_e + rhoStrat_0(k)
                  rho001_e = rho001_e + rhoStrat_0(k + 1)

                  rhow_e = rhow_e + 0.5 * (rhoStrat_0(k) + rhoStrat_0(k + 1))
                end if
              end if

              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piU = var%pi(i, j, k + 1) - var_env%pi(i, j, k + 1)
                piD = var%pi(i, j, k) - var_env%pi(i, j, k)

                piU_e = var_env%pi(i, j, k + 1)
                piD_e = var_env%pi(i, j, k)
              else
                piU = var%pi(i, j, k + 1)
                piD = var%pi(i, j, k)
              end if

              if(topography) then
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) - var_env%pi(:, :, :)
                end if
                ! Compute values at cell edges.
                pEdgeU = (jac(i, j, k + 1) * pStratTFC(i, j, k) + jac(i, j, k) &
                    &* pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k &
                    &+ 1))
                met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) + jac(i, &
                    &j, k) * met(i, j, k + 1, 1, 3)) / (jac(i, j, k) + jac(i, &
                    &j, k + 1))
                met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) + jac(i, &
                    &j, k) * met(i, j, k + 1, 2, 3)) / (jac(i, j, k) + jac(i, &
                    &j, k + 1))
                met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) + jac(i, &
                    &j, k) * met(i, j, k + 1, 3, 3)) / (jac(i, j, k) + jac(i, &
                    &j, k + 1))
                piREdgeU = (jac(i + 1, j, k + 1) * var%pi(i + 1, j, k) + jac(i &
                    &+ 1, j, k) * var%pi(i + 1, j, k + 1)) / (jac(i + 1, j, k) &
                    &+ jac(i + 1, j, k + 1))
                piLEdgeU = (jac(i - 1, j, k + 1) * var%pi(i - 1, j, k) + jac(i &
                    &- 1, j, k) * var%pi(i - 1, j, k + 1)) / (jac(i - 1, j, k) &
                    &+ jac(i - 1, j, k + 1))
                piFEdgeU = (jac(i, j + 1, k + 1) * var%pi(i, j + 1, k) &
                    &+ jac(i, j + 1, k) * var%pi(i, j + 1, k + 1)) / (jac(i, j &
                    &+ 1, k) + jac(i, j + 1, k + 1))
                piBEdgeU = (jac(i, j - 1, k + 1) * var%pi(i, j - 1, k) &
                    &+ jac(i, j - 1, k) * var%pi(i, j - 1, k + 1)) / (jac(i, j &
                    &- 1, k) + jac(i, j - 1, k + 1))
                ! Compute pressure gradient component.
                piGrad = kappaInv * MaInv2 * pEdgeU / rhow * (met13EdgeU &
                    &* (piREdgeU - piLEdgeU) * 0.5 / dx + met23EdgeU &
                    &* (piFEdgeU - piBEdgeU) * 0.5 / dy + met33EdgeU &
                    &* (var%pi(i, j, k + 1) - var%pi(i, j, k)) / dz)
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) + var_env%pi(:, :, :)
                end if
              else
                piGrad = kappaInv * MaInv2 * pstw / rhow * (piU - piD) / dz
              end if

              if(TestCase == "baroclinic_LC") then !FS
                if(topography) then
                  ! Compute values at cell edges.
                  pEdgeU = (jac(i, j, k + 1) * pStratTFC(i, j, k) + jac(i, j, &
                      &k) * pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, &
                      &j, k + 1))
                  met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) + jac(i, &
                      &j, k) * met(i, j, k + 1, 1, 3)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))
                  met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) + jac(i, &
                      &j, k) * met(i, j, k + 1, 2, 3)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))
                  met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) + jac(i, &
                      &j, k) * met(i, j, k + 1, 3, 3)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))
                  piREdgeU = (jac(i + 1, j, k + 1) * var_env%pi(i + 1, j, k) &
                      &+ jac(i + 1, j, k) * var_env%pi(i + 1, j, k + 1)) &
                      &/ (jac(i + 1, j, k) + jac(i + 1, j, k + 1))
                  piLEdgeU = (jac(i - 1, j, k + 1) * var_env%pi(i - 1, j, k) &
                      &+ jac(i - 1, j, k) * var_env%pi(i - 1, j, k + 1)) &
                      &/ (jac(i - 1, j, k) + jac(i - 1, j, k + 1))
                  piFEdgeU = (jac(i, j + 1, k + 1) * var_env%pi(i, j + 1, k) &
                      &+ jac(i, j + 1, k) * var_env%pi(i, j + 1, k + 1)) &
                      &/ (jac(i, j + 1, k) + jac(i, j + 1, k + 1))
                  piBEdgeU = (jac(i, j - 1, k + 1) * var_env%pi(i, j - 1, k) &
                      &+ jac(i, j - 1, k) * var_env%pi(i, j - 1, k + 1)) &
                      &/ (jac(i, j - 1, k) + jac(i, j - 1, k + 1))
                  ! Compute pressure gradient component.
                  piGrad = piGrad + kappaInv * MaInv2 * (pEdgeU / rhow &
                      &- pEdgeU / rhow_e) * (met13EdgeU * (piREdgeU &
                      &- piLEdgeU) * 0.5 / dx + met23EdgeU * (piFEdgeU &
                      &- piBEdgeU) * 0.5 / dy + met33EdgeU * (var_env%pi(i, j, &
                      &k + 1) - var_env%pi(i, j, k)) / dz)
                else
                  piGrad = piGrad + kappaInv * MaInv2 * (pstw / rhow - pstw_e &
                      &/ rhow_e) * (piU_e - piD_e) / dz
                end if
              end if

              ! Gravity-wave forcing and metric terms due to topography growth.
              if(raytracer .or. testCase == "mountainwave") then
                volfcz = (jac(i, j, k + 1) * force(i, j, k, 3) + jac(i, j, k) &
                    &* force(i, j, k + 1, 3)) / (jac(i, j, k) + jac(i, j, k &
                    &+ 1))
              else
                volfcz = 0.0
              end if

              ! wstar
              wvert = var%w(i, j, k)

              if(TestCase == "baroclinic_LC") then
                if(topography) then
                  buoy = - g_ndim * (jac(i, j, k + 1) * (rhopOld(i, j, k) &
                      &/ rho000 / jac(i, j, k) - var_env%rhop(i, j, k) &
                      &/ rho000_e / jac(i, j, k)) + jac(i, j, k) * (rhopOld(i, &
                      &j, k + 1) / rho001 / jac(i, j, k + 1) - var_env%rhop(i, &
                      &j, k + 1) / rho001_e / jac(i, j, k + 1))) / (jac(i, j, &
                      &k) + jac(i, j, k + 1))
                else
                  buoy = - g_ndim * 0.5 * (rhopOld(i, j, k) / rho000 &
                      &- var_env%rhop(i, j, k) / rho000_e + rhopOld(i, j, k &
                      &+ 1) / rho001 - var_env%rhop(i, j, k + 1) / rho001_e)
                end if
              else
                if(topography) then
                  buoy = - g_ndim * (jac(i, j, k + 1) * rhopOld(i, j, k) &
                      &/ rho000 / jac(i, j, k) + jac(i, j, k) * rhopOld(i, j, &
                      &k + 1) / rho001 / jac(i, j, k + 1)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))
                else
                  buoy = - g_ndim * 0.5 * (rhopOld(i, j, k) / rho000 &
                      &+ rhopOld(i, j, k + 1) / rho001)
                end if
              end if

              if(model == "compressible") then ! Muliply with JP
                JPU = jac(i, j, k) * jac(i, j, k + 1) * (var%P(i, j, k) &
                    &+ var%P(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
                wAst = wvert + dt * (buoy - piGrad + volfcz / rhow) * JPU
              else
                wAst = wvert + dt * (buoy - piGrad + volfcz / rhow)
              end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping

                  wAst = wAst - dt * 0.5 * (kw_hs(k) + kw_hs(k + 1)) * wvert
                end if
              end if

              if(spongeLayer) then
                if(topography) then
                  wAst = wAst - dt * (jac(i, j, k + 1) * kr_sp_w_tfc(i, j, k) &
                      &+ jac(i, j, k) * kr_sp_w_tfc(i, j, k + 1)) / (jac(i, j, &
                      &k) + jac(i, j, k + 1)) * wvert
                else
                  wAst = wAst - dt * 0.5 * (kr_sp_w(j, k) + kr_sp_w(j, k + 1)) &
                      &* wvert
                end if
              end if

              var%w(i, j, k) = wAst
            end do
          end do
        end do
      else if(int_mod == "impl") then
        ! heating due to relaxation, entropy diffusion and GWs, its
        ! horizontal mean and the horizontal-mean vertical wind
        ! resulting from it

        if(model == "compressible") then
          heat = 0.
          S_bar = 0.
          w_0 = 0.
        elseif(heatingONK14 .or. TurbScheme .or. rayTracer) then
          !call heat_w0(var,flux,dt,heat,S_bar,w_0)
          call calculate_heating(var, flux, heat)
        else
          heat = 0.
          S_bar = 0.
          w_0 = 0.
        end if

        do k = k0, k1
          pstw = 0.5 * (Pstrat(k) + Pstrat(k + 1))
          pstw_0 = 0.5 * (Pstrat_0(k) + Pstrat_0(k + 1))

          if(TestCase == "baroclinic_LC") then
            pstw_e = 0.5 * (Pstrat_0(k) + Pstrat_0(k + 1))
          end if

          do j = 1, ny
            do i = 1, nx
              rho000 = var%rho(i, j, k)
              rho001 = var%rho(i, j, k + 1)

              if(topography) then
                rhow = (jac(i, j, k + 1) * var%rho(i, j, k) + jac(i, j, k) &
                    &* var%rho(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))

                rhoStratEdgeU = (jac(i, j, k + 1) * rhoStratTFC(i, j, k) &
                    &+ jac(i, j, k) * rhoStratTFC(i, j, k + 1)) / (jac(i, j, &
                    &k) + jac(i, j, k + 1))
                rho000 = rho000 + rhoStratTFC(i, j, k)
                rho001 = rho001 + rhoStratTFC(i, j, k + 1)
                rhow = rhow + rhoStratEdgeU
              else
                rhow = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1))

                rho000 = rho000 + rhoStrat(k)
                rho001 = rho001 + rhoStrat(k + 1)

                rhow = rhow + rhoStratTilde(k)
              end if

              if(TestCase == "baroclinic_LC") then
                rho000_e = var_env%rho(i, j, k)
                rho001_e = var_env%rho(i, j, k + 1)

                if(topography) then
                  rhow_e = (jac(i, j, k + 1) * var_env%rho(i, j, k) + jac(i, &
                      &j, k) * var_env%rho(i, j, k + 1)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))

                  rho000_e = rho000_e + rhoStratTFC(i, j, k)
                  rho001_e = rho001_e + rhoStratTFC(i, j, k + 1)
                  rhow_e = rhow_e + (jac(i, j, k + 1) * rhoStratTFC(i, j, k) &
                      &+ jac(i, j, k) * rhoStratTFC(i, j, k + 1)) / (jac(i, j, &
                      &k) + jac(i, j, k + 1))
                else
                  rhow_e = 0.5 * (var_env%rho(i, j, k) + var_env%rho(i, j, k &
                      &+ 1))

                  rho000_e = rho000_e + rhoStrat_0(k)
                  rho001_e = rho001_e + rhoStrat_0(k + 1)

                  rhow_e = rhow_e + 0.5 * (rhoStrat_0(k) + rhoStrat_0(k + 1))
                end if
              end if

              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piU = var%pi(i, j, k + 1) - var_env%pi(i, j, k + 1)
                piD = var%pi(i, j, k) - var_env%pi(i, j, k)

                piU_e = var_env%pi(i, j, k + 1)
                piD_e = var_env%pi(i, j, k)
              else
                piU = var%pi(i, j, k + 1)
                piD = var%pi(i, j, k)
              end if

              if(topography) then
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) - var_env%pi(:, :, :)
                end if
                ! Compute values at cell edges.
                pEdgeU = (jac(i, j, k + 1) * pStratTFC(i, j, k) + jac(i, j, k) &
                    &* pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k &
                    &+ 1))
                met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) + jac(i, &
                    &j, k) * met(i, j, k + 1, 1, 3)) / (jac(i, j, k) + jac(i, &
                    &j, k + 1))
                met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) + jac(i, &
                    &j, k) * met(i, j, k + 1, 2, 3)) / (jac(i, j, k) + jac(i, &
                    &j, k + 1))
                met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) + jac(i, &
                    &j, k) * met(i, j, k + 1, 3, 3)) / (jac(i, j, k) + jac(i, &
                    &j, k + 1))
                piREdgeU = (jac(i + 1, j, k + 1) * var%pi(i + 1, j, k) + jac(i &
                    &+ 1, j, k) * var%pi(i + 1, j, k + 1)) / (jac(i + 1, j, k) &
                    &+ jac(i + 1, j, k + 1))
                piLEdgeU = (jac(i - 1, j, k + 1) * var%pi(i - 1, j, k) + jac(i &
                    &- 1, j, k) * var%pi(i - 1, j, k + 1)) / (jac(i - 1, j, k) &
                    &+ jac(i - 1, j, k + 1))
                piFEdgeU = (jac(i, j + 1, k + 1) * var%pi(i, j + 1, k) &
                    &+ jac(i, j + 1, k) * var%pi(i, j + 1, k + 1)) / (jac(i, j &
                    &+ 1, k) + jac(i, j + 1, k + 1))
                piBEdgeU = (jac(i, j - 1, k + 1) * var%pi(i, j - 1, k) &
                    &+ jac(i, j - 1, k) * var%pi(i, j - 1, k + 1)) / (jac(i, j &
                    &- 1, k) + jac(i, j - 1, k + 1))
                ! Compute pressure gradient component.
                piGrad = kappaInv * MaInv2 * pEdgeU / rhow * (met13EdgeU &
                    &* (piREdgeU - piLEdgeU) * 0.5 / dx + met23EdgeU &
                    &* (piFEdgeU - piBEdgeU) * 0.5 / dy + met33EdgeU &
                    &* (var%pi(i, j, k + 1) - var%pi(i, j, k)) / dz)
                if(testCase == "baroclinic_LC") then
                  var%pi(:, :, :) = var%pi(:, :, :) + var_env%pi(:, :, :)
                end if
              else
                piGrad = kappaInv * MaInv2 * pstw / rhow * (piU - piD) / dz
              end if

              if(TestCase == "baroclinic_LC") then !FS
                if(topography) then
                  ! Compute values at cell edges.
                  pEdgeU = (jac(i, j, k + 1) * pStratTFC(i, j, k) + jac(i, j, &
                      &k) * pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, &
                      &j, k + 1))
                  met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) + jac(i, &
                      &j, k) * met(i, j, k + 1, 1, 3)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))
                  met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) + jac(i, &
                      &j, k) * met(i, j, k + 1, 2, 3)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))
                  met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) + jac(i, &
                      &j, k) * met(i, j, k + 1, 3, 3)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))
                  piREdgeU = (jac(i + 1, j, k + 1) * var_env%pi(i + 1, j, k) &
                      &+ jac(i + 1, j, k) * var_env%pi(i + 1, j, k + 1)) &
                      &/ (jac(i + 1, j, k) + jac(i + 1, j, k + 1))
                  piLEdgeU = (jac(i - 1, j, k + 1) * var_env%pi(i - 1, j, k) &
                      &+ jac(i - 1, j, k) * var_env%pi(i - 1, j, k + 1)) &
                      &/ (jac(i - 1, j, k) + jac(i - 1, j, k + 1))
                  piFEdgeU = (jac(i, j + 1, k + 1) * var_env%pi(i, j + 1, k) &
                      &+ jac(i, j + 1, k) * var_env%pi(i, j + 1, k + 1)) &
                      &/ (jac(i, j + 1, k) + jac(i, j + 1, k + 1))
                  piBEdgeU = (jac(i, j - 1, k + 1) * var_env%pi(i, j - 1, k) &
                      &+ jac(i, j - 1, k) * var_env%pi(i, j - 1, k + 1)) &
                      &/ (jac(i, j - 1, k) + jac(i, j - 1, k + 1))
                  ! Compute pressure gradient component.
                  piGrad = piGrad + kappaInv * MaInv2 * (pEdgeU / rhow &
                      &- pEdgeU / rhow_e) * (met13EdgeU * (piREdgeU &
                      &- piLEdgeU) * 0.5 / dx + met23EdgeU * (piFEdgeU &
                      &- piBEdgeU) * 0.5 / dy + met33EdgeU * (var_env%pi(i, j, &
                      &k + 1) - var_env%pi(i, j, k)) / dz)
                else
                  piGrad = piGrad + kappaInv * MaInv2 * (pstw / rhow - pstw_e &
                      &/ rhow_e) * (piU_e - piD_e) / dz
                end if
              end if

              ! Gravity-wave forcing and metric terms due to topography growth.
              if(raytracer .or. testCase == "mountainwave") then
                volfcz = (jac(i, j, k + 1) * force(i, j, k, 3) + jac(i, j, k) &
                    &* force(i, j, k + 1, 3)) / (jac(i, j, k) + jac(i, j, k &
                    &+ 1))
              else
                volfcz = 0.0
              end if

              ! wstar
              wvert = var%w(i, j, k)

              ! squared Brunt-Vaisala frequency averaged to half
              ! levels
              ! (could be done a bit nicer by determining this without
              ! averaging directly from the reference-atmosphere
              ! density)
              if(topography) then
                bvsstw = (jac(i, j, k + 1) * bvsStratTFC(i, j, k) + jac(i, j, &
                    &k) * bvsStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, &
                    &j, k + 1))
              else
                bvsstw = 0.5 * (bvsStrat(k) + bvsStrat(k + 1))
              end if

              facw = 1.0

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping

                  facw = facw + dt * 0.5 * (kw_hs(k) + kw_hs(k + 1))
                end if
              end if

              if(spongeLayer) then
                if(topography) then
                  facw = facw + dt * (jac(i, j, k + 1) * kr_sp_w_tfc(i, j, k) &
                      &+ jac(i, j, k) * kr_sp_w_tfc(i, j, k + 1)) / (jac(i, j, &
                      &k) + jac(i, j, k + 1))
                else
                  facw = facw + dt * 0.5 * (kr_sp_w(j, k) + kr_sp_w(j, k + 1))
                end if
              end if

              heat0 = heat(i, j, k)

              heat1 = heat(i, j, k + 1)

              ! Buoyancy is predicted after momentum in implicit steps.
              if(topography) then
                if(testCase == "baroclinic_LC") then
                  buoy = - g_ndim * (jac(i, j, k + 1) * (var%rhop(i, j, k) &
                      &/ rho000 / jac(i, j, k) - var_env%rhop(i, j, k) &
                      &/ rho000_e / jac(i, j, k)) + jac(i, j, k) &
                      &* (var%rhop(i, j, k + 1) / rho001 / jac(i, j, k + 1) &
                      &- var_env%rhop(i, j, k + 1) / rho001_e / jac(i, j, k &
                      &+ 1))) / (jac(i, j, k) + jac(i, j, k + 1))
                else
                  buoy = - g_ndim * (jac(i, j, k + 1) * var%rhop(i, j, k) &
                      &/ rho000 / jac(i, j, k) + jac(i, j, k) * var%rhop(i, j, &
                      &k + 1) / rho001 / jac(i, j, k + 1)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))
                end if
              end if

              if(topography) then
                if(model == "compressible") then ! Interpolate (U/JP)
                  ! Calculate JP on cell-interfaces.
                  JPR = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i + 1, j, &
                      &k) * var%P(i + 1, j, k)) ! right
                  JPL = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i - 1, j, &
                      &k) * var%P(i - 1, j, k)) ! left
                  JPF = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i, j + 1, &
                      &k) * var%P(i, j + 1, k)) ! forward
                  JPB = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i, j - 1, &
                      &k) * var%P(i, j - 1, k)) ! forward
                  JPUR = 0.5 * (jac(i, j, k + 1) * var%P(i, j, k + 1) + jac(i &
                      &+ 1, j, k + 1) * var%P(i + 1, j, k + 1)) ! right at k+1
                  JPUL = 0.5 * (jac(i, j, k + 1) * var%P(i, j, k + 1) + jac(i &
                      &- 1, j, k + 1) * var%P(i - 1, j, k + 1)) ! left at k+1
                  JPUF = 0.5 * (jac(i, j, k + 1) * var%P(i, j, k + 1) + jac(i, &
                      &j + 1, k + 1) * var%P(i, j + 1, k + 1)) ! forward at k+1
                  JPUB = 0.5 * (jac(i, j, k + 1) * var%P(i, j, k + 1) + jac(i, &
                      &j - 1, k + 1) * var%P(i, j - 1, k + 1)) ! forward at k+1

                  ! Calculate U/JP at k and k+1
                  uC = 0.5 * (var%u(i, j, k) / JPR + var%u(i - 1, j, k) / JPL)
                  uU = 0.5 * (var%u(i, j, k + 1) / JPUR + var%u(i - 1, j, k &
                      &+ 1) / JPUL)
                  vC = 0.5 * (var%v(i, j, k) / JPF + var%v(i, j - 1, k) / JPB)
                  vU = 0.5 * (var%v(i, j, k + 1) / JPUF + var%v(i, j - 1, k &
                      &+ 1) / JPUB)
                else ! Interpolate u
                  uC = 0.5 * (var%u(i, j, k) + var%u(i - 1, j, k))
                  uU = 0.5 * (var%u(i, j, k + 1) + var%u(i - 1, j, k + 1))
                  vC = 0.5 * (var%v(i, j, k) + var%v(i, j - 1, k))
                  vU = 0.5 * (var%v(i, j, k + 1) + var%v(i, j - 1, k + 1))
                end if

                if(model == "compressible") then ! Muliply with JP
                  JPU = jac(i, j, k) * jac(i, j, k + 1) * (var%P(i, j, k) &
                      &+ var%P(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
                  wAst = 1.0 / (facw + bvsstw * dt ** 2.0) * (wvert - dt &
                      &* piGrad * JPU + dt * buoy * JPU + dt * volfcz / rhow &
                      &* JPU + JPU * bvsstw * dt ** 2.0 * (jac(i, j, k + 1) &
                      &* (met(i, j, k, 1, 3) * uC + met(i, j, k, 2, 3) * vC) &
                      &+ jac(i, j, k) * (met(i, j, k + 1, 1, 3) * uU + met(i, &
                      &j, k + 1, 2, 3) * vU)) / (jac(i, j, k) + jac(i, j, k &
                      &+ 1)))
                else
                  wAst = 1.0 / (facw + rhoStratEdgeU / rhow * bvsstw * dt &
                      &** 2.0) * (wvert - dt * piGrad + dt * buoy + dt &
                      &* volfcz / rhow + rhoStratEdgeU / rhow * bvsstw * dt &
                      &** 2.0 * (jac(i, j, k + 1) * (met(i, j, k, 1, 3) * uC &
                      &+ met(i, j, k, 2, 3) * vC) + jac(i, j, k) * (met(i, j, &
                      &k + 1, 1, 3) * uU + met(i, j, k + 1, 2, 3) * vU)) &
                      &/ (jac(i, j, k) + jac(i, j, k + 1)))
                end if
              else
                if(TestCase == "baroclinic_LC") then
                  wAst = 1.0 / (facw + rhoStratTilde(k) / rhow * pstw / pstw_0 &
                      &* bvsstw * dt ** 2) * (wvert - dt * piGrad - dt &
                      &* g_ndim * 0.5 * (rhopOld(i, j, k) / rho000 &
                      &- var_env%rhop(i, j, k) / rho000_e + rhopOld(i, j, k &
                      &+ 1) / rho001 - var_env%rhop(i, j, k + 1) / rho001_e &
                      &+ dt * (rhoStrat(k) / Pstrat_0(k) * heat0 / rho000 &
                      &+ rhoStrat(k + 1) / Pstrat_0(k + 1) * heat1 / rho001)))
                else
                  wAst = 1.0 / (facw + rhoStratTilde(k) / rhow * pstw / pstw_0 &
                      &* bvsstw * dt ** 2) * (wvert - dt * piGrad - dt &
                      &* g_ndim * 0.5 * (rhopOld(i, j, k) / rho000 &
                      &+ rhopOld(i, j, k + 1) / rho001 + dt * (rhoStrat(k) &
                      &/ Pstrat_0(k) * heat0 / rho000 + rhoStrat(k + 1) &
                      &/ Pstrat_0(k + 1) * heat1 / rho001)))
                end if
              end if

              var%w(i, j, k) = wAst
            end do
          end do
        end do
      else
        stop 'ERROR: unknown int_mod'
      end if
    else
      stop 'ERROR: unknown mmp_mod'
    end if

    if(mmp_mod == 'rhs') then
      if(int_mod == 'expl') then
        spongeLayer = spongeLayer_s
      else if(int_mod == 'impl') then
        kr_sp = kr_sp / facray
        kr_sp_w = kr_sp_w / facray
        alprlx = alprlx / facray
        if(topography) then
          kr_sp_tfc = kr_sp_tfc / facray
          kr_sp_w_tfc = kr_sp_w_tfc / facray
        end if
      end if
    end if

  end subroutine momentumPredictor

  !--------------------------------------------------------------------------

  subroutine massUpdate(var, flux, dt, q, m, upd_var, upd_mod, int_mod, facray)
    !-----------------------------
    ! adds mass flux to cell mass
    !-----------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var

    ! upd_var decides what is to be propagated in time:
    ! rho => total density
    ! rhop => density fluctuations

    ! upd_mod decides which part of the equation is to be used:
    ! tot => total equation (always the case for the total density)
    ! lhs => only advection and molecular and turbulent diffusive fluxes
    !        on the left-hand side of the density-fluctuation equation
    ! rhs => only the right-hand side of the density-fluctuation equation

    ! int_mod discriminates implicit and explicit time stepping:
    ! expl => explicit time stepping
    !         (always the case for the total density)
    !         RK sub step for the total density
    !         Euler step for the rhs of the density-fluctuation equation
    ! impl => implicit-time-step part without pressure-gradient term
    !         (only for the density fluctuations, only for rhs)

    ! facray multiplies the Rayleigh-damping terms so that they are only
    ! handled in the implicit time stepping (sponge and immersed boundary)
    character(len = *), intent(in) :: upd_var, upd_mod, int_mod

    type(flux_type), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux)
    ! dir = 1..3 > f-, g- and h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW

    !UAC real, intent(in) :: dt
    real, intent(in) :: dt, facray
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: q

    integer, intent(in) :: m
    integer :: i00, j00

    ! local variables
    integer :: i, j, k, l
    real :: fL, fR ! flux Left/Right
    real :: gB, gF ! flux Backward/Forward
    real :: hD, hU ! flux Downward/Upward
    real :: fluxDiff ! convective part
    real :: F ! F(phi)

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: heat

    real :: buoy0, buoy, rho, rhow, rhowm, rhop, wvrt, facw, facr, pstw, &
        &pstwm, piU, piD, piGrad

    ! SK compressible: JP on interfaces right, left, forward, backward, upward, downward
    real :: JPR, JPL, JPF, JPB, JPU, JPD

    ! TFC variables
    real :: pEdgeU, pEdgeD
    real :: met13EdgeU, met23EdgeU, met33EdgeU, met13EdgeD, met23EdgeD, &
        &met33EdgeD
    real :: piREdgeU, piLEdgeU, piFEdgeU, piBEdgeU, piREdgeD, piLEdgeD, &
        &piFEdgeD, piBEdgeD
    real :: piGradZEdgeU, piGradZEdgeD

    real :: rho_p

    real, dimension(- nbz:nz + nbz) :: w_0
    real, dimension(- nbz:nz + nbz) :: S_bar
    real :: heat_flc

    real :: rho_e, pstw_e, pstwm_e, rhow_e, rhowm_e

    real, dimension(1:nz) :: sum_local, sum_global

    real, dimension(- nbz:nz + nbz) :: rhopw_bar

    real :: ymax, yloc

    ymax = ly_dim(1) / lRef

    if(correctDivError) then
      print *, 'ERROR: correction divergence error not allowed'
      stop
    end if

    ! Constant background density in Boussinesq model.
    if(model == "Boussinesq" .and. upd_var == "rho") then
      return
    end if

    ! init q
    if(m == 1) q = 0.

    if(upd_var == "rho") then
      if(upd_mod /= "tot" .and. upd_mod /= "lhs") then
        print *, 'ERROR: wrong upd_mod for upd_var = rho'
        stop
      end if

      if(int_mod /= "expl") stop 'ERROR: wrong int_mod for upd_var = rho'

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            fL = flux%rho(i - 1, j, k, 1) ! mass flux accros left cell edge
            fR = flux%rho(i, j, k, 1) ! right
            gB = flux%rho(i, j - 1, k, 2) ! backward
            gF = flux%rho(i, j, k, 2) ! forward
            hD = flux%rho(i, j, k - 1, 3) ! downward
            hU = flux%rho(i, j, k, 3) ! upward

            ! convective part
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

            ! Adjust mass flux divergence.
            if(topography) then
              fluxDiff = fluxDiff / jac(i, j, k)
            end if

            ! F(phi)
            F = - fluxDiff

            ! density relaxation
            if(dens_relax) then
              if(background /= "HeldSuarez") then
                stop 'ERROR: density relaxation only ready for background &
                    &= HeldSuarez'
              end if

              rho = var%rho(i, j, k) + rhoStrat(k)

              rho_e = Pstrat(k) / the_env_pp(i, j, k)

              F = F - kt_hs(j, k) * (rho - rho_e)
            end if

            ! update: q(m-1) -> q(m)
            q(i, j, k) = dt * F + alphaRK(m) * q(i, j, k)

            ! update density
            var%rho(i, j, k) = var%rho(i, j, k) + betaRK(m) * q(i, j, k)
          end do
        end do
      end do
    else if(upd_var == "rhop") then
      if(upd_mod == "tot") then
        if(int_mod /= "expl") then
          stop 'ERROR: wrong int_mod for upd_mod = tot'
        end if

        ! heating due to relaxation, entropy diffusion and GWs, its
        ! horizontal mean and the horizontal-mean vertical wind
        ! resulting from it

        if(heatingONK14 .or. TurbScheme .or. rayTracer) then
          !call heat_w0(var,flux,dt,heat,S_bar,w_0)
          call calculate_heating(var, flux, heat)
        else
          heat = 0.
          S_bar = 0.
          w_0 = 0.
        end if

        !! horizontal mean of the vertical density-fluctuation flux

        !rhopw_bar = 0.

        !sum_local = 0.
        !sum_global = 0.

        !do k = 1,nz-1
        !   do j = 1,ny
        !      do i = 1,nx
        !         sum_local(k) = sum_local(k) + flux(i,j,k,3,6)
        !      end do
        !   end do
        !end do
        !call mpi_allreduce(sum_local(1),sum_global(1),&
        !     nz-1+1,&
        !     mpi_double_precision,mpi_sum,comm,ierror)
        !sum_global = sum_global/(sizeX*sizeY)

        !rhopw_bar(1:nz-1) = sum_global(1:nz-1)

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              fL = flux%rhop(i - 1, j, k, 1) ! mass flux accros left cell edge
              fR = flux%rhop(i, j, k, 1) ! right
              gB = flux%rhop(i, j - 1, k, 2) ! backward
              gF = flux%rhop(i, j, k, 2) ! forward
              hD = flux%rhop(i, j, k - 1, 3) ! downward
              hU = flux%rhop(i, j, k, 3) ! upward

              ! convective part
              fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

              if(topography) then
                fluxDiff = fluxDiff / jac(i, j, k)
              end if

              rhop = var%rhop(i, j, k)

              rho = var%rho(i, j, k) + rhoStrat(k)

              if(topography) then
                wvrt = 0.5 * (vertWindTFC(i, j, k, var) + vertWindTFC(i, j, k &
                    &- 1, var))
              else
                wvrt = 0.5 * (var%w(i, j, k) + var%w(i, j, k - 1))
              end if

              heat_flc = heat(i, j, k)

              if(topography) then
                F = - fluxDiff + rhoStratTFC(i, j, k) / g_ndim &
                    &* bvsStratTFC(i, j, k) * wvrt
              else
                F = - fluxDiff + PStrat(k) / PStrat_0(k) * rhoStrat(k) &
                    &/ g_ndim * bvsStrat(k) * wvrt + rhoStrat(k) / Pstrat_0(k) &
                    &* heat_flc
              end if

              ! density relaxation
              if(dens_relax) then
                if(background /= "HeldSuarez") then
                  stop 'ERROR: density relaxation only ready for background &
                      &= HeldSuarez'
                end if

                rho = var%rho(i, j, k) + rhoStrat(k)

                rho_e = Pstrat_0(k) / the_env_pp(i, j, k)

                F = F - kt_hs(j, k) * (rho - rho_e)
              end if

              ! update: q(m-1) -> q(m)
              q(i, j, k) = dt * F + alphaRK(m) * q(i, j, k)

              ! update density
              var%rhop(i, j, k) = var%rhop(i, j, k) + betaRK(m) * q(i, j, k)
            end do
          end do
        end do
      else if(upd_mod == "lhs") then
        if(int_mod /= "expl") then
          stop 'ERROR: wrong int_mod for upd_mod = lhs'
        end if

        !! horizontal mean of the vertical density-fluctuation flux

        !rhopw_bar = 0.

        !sum_local = 0.
        !sum_global = 0.

        !do k = 1,nz-1
        !   do j = 1,ny
        !      do i = 1,nx
        !         sum_local(k) = sum_local(k) + flux(i,j,k,3,6)
        !      end do
        !   end do
        !end do
        !call mpi_allreduce(sum_local(1),sum_global(1),&
        !     nz-1+1,&
        !     mpi_double_precision,mpi_sum,comm,ierror)
        !sum_global = sum_global/(sizeX*sizeY)

        !rhopw_bar(1:nz-1) = sum_global(1:nz-1)

        ! heating due to relaxation, entropy diffusion and GWs
        if(model == "compressible" .and. (TurbScheme .or. rayTracer)) then
          call calculate_heating(var, flux, heat)
        else
          heat = 0.0
        end if

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              fL = flux%rhop(i - 1, j, k, 1) ! mass flux accros left cell edge
              fR = flux%rhop(i, j, k, 1) ! right
              gB = flux%rhop(i, j - 1, k, 2) ! backward
              gF = flux%rhop(i, j, k, 2) ! forward
              hD = flux%rhop(i, j, k - 1, 3) ! downward
              hU = flux%rhop(i, j, k, 3) ! upward

              ! convective part
              fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

              if(topography) then
                fluxDiff = fluxDiff / jac(i, j, k)
              end if

              ! F(phi)
              if(model == "compressible") then
                F = - fluxDiff + heat(i, j, k) / thetaStratTFC(i, j, k)
              else
                F = - fluxDiff
              end if

              ! density relaxation
              if(dens_relax) then
                if(background /= "HeldSuarez") then
                  stop 'ERROR: density relaxation only ready for background &
                      &= HeldSuarez'
                end if

                rho = var%rho(i, j, k) + rhoStrat(k)

                rho_e = Pstrat_0(k) / the_env_pp(i, j, k)

                F = F - kt_hs(j, k) * (rho - rho_e)
              end if

              ! update: q(m-1) -> q(m)
              q(i, j, k) = dt * F + alphaRK(m) * q(i, j, k)

              ! update density
              var%rhop(i, j, k) = var%rhop(i, j, k) + betaRK(m) * q(i, j, k)
            end do
          end do
        end do
      else if(upd_mod == "rhs") then
        ! calculate bstar ...

        ! heating due to relaxation, entropy diffusion and GWs, its
        ! horizontal mean and the horizontal-mean vertical wind
        ! resulting from it

        if(model == "compressible") then
          heat = 0.
          S_bar = 0.
          w_0 = 0.
        elseif(heatingONK14 .or. TurbScheme .or. rayTracer) then
          !call heat_w0(var,flux,dt,heat,S_bar,w_0)
          call calculate_heating(var, flux, heat)
        else
          heat = 0.
          S_bar = 0.
          w_0 = 0.
        end if

        if(int_mod == "impl") then
          kr_sp = kr_sp * facray
          kr_sp_w = kr_sp_w * facray
          alprlx = alprlx * facray
          if(topography) then
            kr_sp_tfc = kr_sp_tfc * facray
            kr_sp_w_tfc = kr_sp_w_tfc * facray
          end if

          do k = 1, nz
            pstw = 0.5 * (Pstrat(k) + Pstrat(k + 1))
            pstwm = 0.5 * (Pstrat(k - 1) + Pstrat(k))

            do j = 1, ny
              do i = 1, nx
                rho = var%rho(i, j, k)

                if(topography) then
                  rhow = (jac(i, j, k + 1) * var%rho(i, j, k) + jac(i, j, k) &
                      &* var%rho(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k &
                      &+ 1))
                  rhowm = (jac(i, j, k - 1) * var%rho(i, j, k) + jac(i, j, k) &
                      &* var%rho(i, j, k - 1)) / (jac(i, j, k) + jac(i, j, k &
                      &- 1))

                  rho = rho + rhoStratTFC(i, j, k)
                  rhow = rhow + (jac(i, j, k + 1) * rhoStratTFC(i, j, k) &
                      &+ jac(i, j, k) * rhoStratTFC(i, j, k + 1)) / (jac(i, j, &
                      &k) + jac(i, j, k + 1))
                  rhowm = rhowm + (jac(i, j, k - 1) * rhoStratTFC(i, j, k) &
                      &+ jac(i, j, k) * rhoStratTFC(i, j, k - 1)) / (jac(i, j, &
                      &k) + jac(i, j, k - 1))
                else
                  rhow = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1))
                  rhowm = 0.5 * (var%rho(i, j, k - 1) + var%rho(i, j, k))

                  rho = rho + rhoStrat(k)
                  rhow = rhow + rhoStratTilde(k)
                  rhowm = rhowm + rhoStratTilde(k - 1)
                end if

                if(topography) then
                  ! Momentum is predicted before buoyancy in implicit
                  ! steps.
                  if(model == "compressible") then
                    ! Calculate w/(JP).
                    JPU = jac(i, j, k) * jac(i, j, k + 1) * (var%P(i, j, k) &
                        &+ var%P(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k &
                        &+ 1)) !upward
                    JPD = jac(i, j, k) * jac(i, j, k - 1) * (var%P(i, j, k) &
                        &+ var%P(i, j, k - 1)) / (jac(i, j, k) + jac(i, j, k &
                        &- 1)) ! downward
                    wvrt = 0.5 * (wOldTFC(i, j, k) / JPU + wOldTFC(i, j, k &
                        &- 1) / JPD)
                  else
                    wvrt = 0.5 * (wOldTFC(i, j, k) + wOldTFC(i, j, k - 1))
                  end if
                else
                  wvrt = 0.5 * (var%w(i, j, k) + var%w(i, j, k - 1))
                end if

                heat_flc = heat(i, j, k)

                if(topography) then
                  if(testCase == "baroclinic_LC") then
                    var%pi(:, :, :) = var%pi(:, :, :) - var_env%pi(:, :, :)
                  end if
                  ! Compute P coefficients.
                  pEdgeU = (jac(i, j, k + 1) * pStratTFC(i, j, k) + jac(i, j, &
                      &k) * pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, &
                      &j, k + 1))
                  pEdgeD = (jac(i, j, k - 1) * pStratTFC(i, j, k) + jac(i, j, &
                      &k) * pStratTFC(i, j, k - 1)) / (jac(i, j, k) + jac(i, &
                      &j, k - 1))
                  ! Interpolate metric-tensor elements.
                  met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) + jac(i, &
                      &j, k) * met(i, j, k + 1, 1, 3)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))
                  met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) + jac(i, &
                      &j, k) * met(i, j, k + 1, 2, 3)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))
                  met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) + jac(i, &
                      &j, k) * met(i, j, k + 1, 3, 3)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))
                  met13EdgeD = (jac(i, j, k - 1) * met(i, j, k, 1, 3) + jac(i, &
                      &j, k) * met(i, j, k - 1, 1, 3)) / (jac(i, j, k) &
                      &+ jac(i, j, k - 1))
                  met23EdgeD = (jac(i, j, k - 1) * met(i, j, k, 2, 3) + jac(i, &
                      &j, k) * met(i, j, k - 1, 2, 3)) / (jac(i, j, k) &
                      &+ jac(i, j, k - 1))
                  met33EdgeD = (jac(i, j, k - 1) * met(i, j, k, 3, 3) + jac(i, &
                      &j, k) * met(i, j, k - 1, 3, 3)) / (jac(i, j, k) &
                      &+ jac(i, j, k - 1))
                  ! Interpolate pressure differences.
                  piREdgeU = (jac(i + 1, j, k + 1) * var%pi(i + 1, j, k) &
                      &+ jac(i + 1, j, k) * var%pi(i + 1, j, k + 1)) / (jac(i &
                      &+ 1, j, k) + jac(i + 1, j, k + 1))
                  piLEdgeU = (jac(i - 1, j, k + 1) * var%pi(i - 1, j, k) &
                      &+ jac(i - 1, j, k) * var%pi(i - 1, j, k + 1)) / (jac(i &
                      &- 1, j, k) + jac(i - 1, j, k + 1))
                  piREdgeD = (jac(i + 1, j, k - 1) * var%pi(i + 1, j, k) &
                      &+ jac(i + 1, j, k) * var%pi(i + 1, j, k - 1)) / (jac(i &
                      &+ 1, j, k) + jac(i + 1, j, k - 1))
                  piLEdgeD = (jac(i - 1, j, k - 1) * var%pi(i - 1, j, k) &
                      &+ jac(i - 1, j, k) * var%pi(i - 1, j, k - 1)) / (jac(i &
                      &- 1, j, k) + jac(i - 1, j, k - 1))
                  piFEdgeU = (jac(i, j + 1, k + 1) * var%pi(i, j + 1, k) &
                      &+ jac(i, j + 1, k) * var%pi(i, j + 1, k + 1)) / (jac(i, &
                      &j + 1, k) + jac(i, j + 1, k + 1))
                  piBEdgeU = (jac(i, j - 1, k + 1) * var%pi(i, j - 1, k) &
                      &+ jac(i, j - 1, k) * var%pi(i, j - 1, k + 1)) / (jac(i, &
                      &j - 1, k) + jac(i, j - 1, k + 1))
                  piFEdgeD = (jac(i, j + 1, k - 1) * var%pi(i, j + 1, k) &
                      &+ jac(i, j + 1, k) * var%pi(i, j + 1, k - 1)) / (jac(i, &
                      &j + 1, k) + jac(i, j + 1, k - 1))
                  piBEdgeD = (jac(i, j - 1, k - 1) * var%pi(i, j - 1, k) &
                      &+ jac(i, j - 1, k) * var%pi(i, j - 1, k - 1)) / (jac(i, &
                      &j - 1, k) + jac(i, j - 1, k - 1))
                  ! Compute pressure gradients.
                  piGradZEdgeU = kappaInv * MaInv2 * pEdgeU / rhow * (0.5 &
                      &* met13EdgeU * (piREdgeU - piLEdgeU) / dx + 0.5 &
                      &* met23EdgeU * (piFEdgeU - piBEdgeU) / dy + met33EdgeU &
                      &* (var%pi(i, j, k + 1) - var%pi(i, j, k)) / dz)
                  piGradZEdgeD = kappaInv * MaInv2 * pEdgeD / rhowm * (0.5 &
                      &* met13EdgeD * (piREdgeD - piLEdgeD) / dx + 0.5 &
                      &* met23EdgeD * (piFEdgeD - piBEdgeD) / dy + met33EdgeD &
                      &* (var%pi(i, j, k) - var%pi(i, j, k - 1)) / dz)
                  ! Adjust at boundaries.
                  if(k == 1 .and. zBoundary == "solid_wall") then
                    piGradZEdgeD = 0.0
                  else if(k == nz .and. zBoundary == "solid_wall") then
                    piGradZEdgeU = 0.0
                  end if
                  ! Interpolate.
                  piGrad = 0.5 * (piGradZEdgeU + piGradZEdgeD)
                  ! Adjust for baroclinic LC.
                  if(testCase == "baroclinic_LC") then
                    var%pi(:, :, :) = var%pi(:, :, :) + var_env%pi(:, :, :)
                    ! Interpolate densities.
                    rhow_e = (jac(i, j, k + 1) * (var_env%rho(i, j, k) &
                        &+ rhoStratTFC(i, j, k)) + jac(i, j, k) &
                        &* (var_env%rho(i, j, k + 1) + rhoStratTFC(i, j, k &
                        &+ 1))) / (jac(i, j, k) + jac(i, j, k + 1))
                    rhowm_e = (jac(i, j, k - 1) * (var_env%rho(i, j, k) &
                        &+ rhoStratTFC(i, j, k)) + jac(i, j, k) &
                        &* (var_env%rho(i, j, k - 1) + rhoStratTFC(i, j, k &
                        &- 1))) / (jac(i, j, k) + jac(i, j, k - 1))
                    ! Interpolate pressure differences.
                    piREdgeU = (jac(i + 1, j, k + 1) * var_env%pi(i + 1, j, k) &
                        &+ jac(i + 1, j, k) * var_env%pi(i + 1, j, k + 1)) &
                        &/ (jac(i + 1, j, k) + jac(i + 1, j, k + 1))
                    piLEdgeU = (jac(i - 1, j, k + 1) * var_env%pi(i - 1, j, k) &
                        &+ jac(i - 1, j, k) * var_env%pi(i - 1, j, k + 1)) &
                        &/ (jac(i - 1, j, k) + jac(i - 1, j, k + 1))
                    piREdgeD = (jac(i + 1, j, k - 1) * var_env%pi(i + 1, j, k) &
                        &+ jac(i + 1, j, k) * var_env%pi(i + 1, j, k - 1)) &
                        &/ (jac(i + 1, j, k) + jac(i + 1, j, k - 1))
                    piLEdgeD = (jac(i - 1, j, k - 1) * var_env%pi(i - 1, j, k) &
                        &+ jac(i - 1, j, k) * var_env%pi(i - 1, j, k - 1)) &
                        &/ (jac(i - 1, j, k) + jac(i - 1, j, k - 1))
                    piFEdgeU = (jac(i, j + 1, k + 1) * var_env%pi(i, j + 1, k) &
                        &+ jac(i, j + 1, k) * var_env%pi(i, j + 1, k + 1)) &
                        &/ (jac(i, j + 1, k) + jac(i, j + 1, k + 1))
                    piBEdgeU = (jac(i, j - 1, k + 1) * var_env%pi(i, j - 1, k) &
                        &+ jac(i, j - 1, k) * var_env%pi(i, j - 1, k + 1)) &
                        &/ (jac(i, j - 1, k) + jac(i, j - 1, k + 1))
                    piFEdgeD = (jac(i, j + 1, k - 1) * var_env%pi(i, j + 1, k) &
                        &+ jac(i, j + 1, k) * var_env%pi(i, j + 1, k - 1)) &
                        &/ (jac(i, j + 1, k) + jac(i, j + 1, k - 1))
                    piBEdgeD = (jac(i, j - 1, k - 1) * var_env%pi(i, j - 1, k) &
                        &+ jac(i, j - 1, k) * var_env%pi(i, j - 1, k - 1)) &
                        &/ (jac(i, j - 1, k) + jac(i, j - 1, k - 1))
                    ! Compute pressure gradients.
                    piGradZEdgeU = kappaInv * MaInv2 * (pEdgeU / rhow - pEdgeU &
                        &/ rhow_e) * (0.5 * met13EdgeU * (piREdgeU - piLEdgeU) &
                        &/ dx + 0.5 * met23EdgeU * (piFEdgeU - piBEdgeU) / dy &
                        &+ met33EdgeU * (var_env%pi(i, j, k + 1) &
                        &- var_env%pi(i, j, k)) / dz)
                    piGradZEdgeD = kappaInv * MaInv2 * (pEdgeD / rhowm &
                        &- pEdgeD / rhowm_e) * (0.5 * met13EdgeD * (piREdgeD &
                        &- piLEdgeD) / dx + 0.5 * met23EdgeD * (piFEdgeD &
                        &- piBEdgeD) / dy + met33EdgeD * (var_env%pi(i, j, k) &
                        &- var_env%pi(i, j, k - 1)) / dz)
                    ! Adjust at boundaries.
                    if(k == 1 .and. zBoundary == "solid_wall") then
                      piGradZEdgeD = 0.0
                    else if(k == nz .and. zBoundary == "solid_wall") then
                      piGradZEdgeU = 0.0
                    end if
                    ! Interpolate.
                    piGrad = piGrad + 0.5 * (piGradZEdgeU + piGradZEdgeD)
                  end if
                else
                  if(TestCase == "baroclinic_LC") then
                    piGrad = kappaInv * MaInv2 * 0.5 * (pstw / rhow &
                        &* (var%pi(i, j, k + 1) - var%pi(i, j, k) &
                        &- var_env%pi(i, j, k + 1) + var_env%pi(i, j, k)) / dz &
                        &+ pstwm / rhowm * (var%pi(i, j, k) - var%pi(i, j, k &
                        &- 1) - var_env%pi(i, j, k) + var_env%pi(i, j, k - 1)) &
                        &/ dz)

                    pstw_e = 0.5 * (pStrat_0(k) + pStrat_0(k + 1))
                    pstwm_e = 0.5 * (pStrat_0(k - 1) + pStrat_0(k))

                    rhow_e = 0.5 * (var_env%rho(i, j, k) + var_env%rho(i, j, k &
                        &+ 1))
                    rhowm_e = 0.5 * (var_env%rho(i, j, k - 1) + var_env%rho(i, &
                        &j, k))

                    rhow_e = rhow_e + 0.5 * (rhoStrat_0(k) + rhoStrat_0(k + 1))
                    rhowm_e = rhowm_e + 0.5 * (rhoStrat_0(k - 1) &
                        &+ rhoStrat_0(k))

                    piGrad = piGrad + kappaInv * MaInv2 * 0.5 * ((pstw / rhow &
                        &- pstw_e / rhow_e) * (var_env%pi(i, j, k + 1) &
                        &- var_env%pi(i, j, k)) / dz + (pstwm / rhowm &
                        &- pstwm_e / rhowm_e) * (var_env%pi(i, j, k) &
                        &- var_env%pi(i, j, k - 1)) / dz)
                  else
                    piGrad = kappaInv * MaInv2 * 0.5 * (pstw / rhow &
                        &* (var%pi(i, j, k + 1) - var%pi(i, j, k)) / dz &
                        &+ pstwm / rhowm * (var%pi(i, j, k) - var%pi(i, j, k &
                        &- 1)) / dz)
                  end if
                end if

                facw = 1.0

                if(TestCase == "baroclinic_LC") then
                  if(background == "HeldSuarez") then
                    ! Rayleigh damping

                    facw = facw + dt * kw_hs(k)
                  end if
                end if

                if(spongeLayer) then
                  if(topography) then
                    facw = facw + dt * kr_sp_w_tfc(i, j, k)
                  else
                    facw = facw + dt * kr_sp_w(j, k)
                  end if
                end if

                if(TestCase == "baroclinic_LC") then
                  rho_e = var_env%rho(i, j, k)

                  if(topography) then
                    rho_e = rho_e + rhoStratTFC(i, j, k)
                  else
                    rho_e = rho_e + rhoStrat_0(k)
                  end if

                  if(topography) then
                    ! Predict buoyancy.
                    buoy = - g_ndim * (var%rhop(i, j, k) / rho &
                        &- var_env%rhop(i, j, k) / rho_e)
                    buoy = 1.0 / (facw + rhoStratTFC(i, j, k) / rho &
                        &* bvsStratTFC(i, j, k) * dt ** 2.0) * (- &
                        &rhoStratTFC(i, j, k) / rho * bvsStratTFC(i, j, k) &
                        &* dt * jac(i, j, k) * (wvrt - dt * piGrad) + facw &
                        &* buoy + rhoStratTFC(i, j, k) / rho * bvsStratTFC(i, &
                        &j, k) * dt * jac(i, j, k) * facw * 0.5 * (met(i, j, &
                        &k, 1, 3) * (var%u(i, j, k) + var%u(i - 1, j, k)) &
                        &+ met(i, j, k, 2, 3) * (var%v(i, j, k) + var%v(i, j &
                        &- 1, k))))
                  else
                    buoy = 1.0 / (facw + rhoStrat(k) / rho * PStrat(k) &
                        &/ PStrat_0(k) * bvsStrat(k) * dt ** 2) * (- &
                        &rhoStrat(k) / rho * PStrat(k) / PStrat_0(k) &
                        &* bvsStrat(k) * dt * (wvrt - dt * piGrad) - facw &
                        &* g_ndim * (var%rhop(i, j, k) / rho - var_env%rhop(i, &
                        &j, k) / rho_e + dt / rho * rhoStrat(k) / Pstrat_0(k) &
                        &* heat_flc))
                  end if

                  buoy = buoy - g_ndim * var_env%rhop(i, j, k) / rho_e
                else
                  if(topography) then
                    ! Predict buoyancy.
                    buoy = - g_ndim * var%rhop(i, j, k) / rho
                    if(model == "compressible") then
                      JPR = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i + 1, &
                          &j, k) * var%P(i + 1, j, k)) ! right
                      JPL = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i - 1, &
                          &j, k) * var%P(i - 1, j, k)) ! left
                      JPF = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i, j &
                          &+ 1, k) * var%P(i, j + 1, k)) ! forward
                      JPB = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i, j &
                          &- 1, k) * var%P(i, j - 1, k)) ! backward
                      buoy = 1.0 / (facw + bvsStratTFC(i, j, k) * dt ** 2.0) &
                          &* (- bvsStratTFC(i, j, k) * dt * jac(i, j, k) &
                          &* (wvrt - dt * piGrad) + facw * buoy &
                          &+ bvsStratTFC(i, j, k) * dt * jac(i, j, k) * facw &
                          &* 0.5 * (met(i, j, k, 1, 3) * (var%u(i, j, k) / JPR &
                          &+ var%u(i - 1, j, k) / JPL) + met(i, j, k, 2, 3) &
                          &* (var%v(i, j, k) / JPF + var%v(i, j - 1, k) / JPB)))
                    else
                      buoy = 1.0 / (facw + rhoStratTFC(i, j, k) / rho &
                          &* bvsStratTFC(i, j, k) * dt ** 2.0) * (- &
                          &rhoStratTFC(i, j, k) / rho * bvsStratTFC(i, j, k) &
                          &* dt * jac(i, j, k) * (wvrt - dt * piGrad) + facw &
                          &* buoy + rhoStratTFC(i, j, k) / rho &
                          &* bvsStratTFC(i, j, k) * dt * jac(i, j, k) * facw &
                          &* 0.5 * (met(i, j, k, 1, 3) * (var%u(i, j, k) &
                          &+ var%u(i - 1, j, k)) + met(i, j, k, 2, 3) &
                          &* (var%v(i, j, k) + var%v(i, j - 1, k))))
                    end if
                  else
                    buoy = 1.0 / (facw + rhoStrat(k) / rho * PStrat(k) &
                        &/ PStrat_0(k) * bvsStrat(k) * dt ** 2) * (- &
                        &rhoStrat(k) / rho * PStrat(k) / PStrat_0(k) &
                        &* bvsStrat(k) * dt * (wvrt - dt * piGrad) - facw &
                        &* g_ndim / rho * (var%rhop(i, j, k) + dt &
                        &* rhoStrat(k) / Pstrat_0(k) * heat_flc))
                  end if
                end if

                var%rhop(i, j, k) = - buoy * rho / g_ndim
              end do
            end do
          end do

          kr_sp = kr_sp / facray
          kr_sp_w = kr_sp_w / facray
          alprlx = alprlx / facray
          if(topography) then
            kr_sp_tfc = kr_sp_tfc / facray
            kr_sp_w_tfc = kr_sp_w_tfc / facray
          end if

        else if(int_mod == "expl") then
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                rhop = var%rhop(i, j, k)

                rho = var%rho(i, j, k)
                if(topography) then
                  rho = rho + rhoStratTFC(i, j, k)
                else
                  rho = rho + rhoStrat(k)
                end if

                if(topography) then
                  if(model == "compressible") then
                    JPU = jac(i, j, k) * jac(i, j, k + 1) * (var%P(i, j, k) &
                        &+ var%P(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k &
                        &+ 1))
                    JPD = jac(i, j, k) * jac(i, j, k - 1) * (var%P(i, j, k) &
                        &+ var%P(i, j, k - 1)) / (jac(i, j, k) + jac(i, j, k &
                        &- 1))
                    wvrt = 0.5 * (vertWindTFC(i, j, k, var) / JPU &
                        &+ vertWindTFC(i, j, k - 1, var) / JPD)
                  else
                    wvrt = 0.5 * (vertWindTFC(i, j, k, var) + vertWindTFC(i, &
                        &j, k - 1, var))
                  end if
                else
                  wvrt = 0.5 * (var%w(i, j, k) + var%w(i, j, k - 1))
                end if

                !heat_flc= heat(i,j,k) - S_bar(k)
                heat_flc = heat(i, j, k)

                if(topography) then
                  buoy = - g_ndim * rhop / rho
                  if(model == "compressible") then
                    buoy = buoy - dt * bvsStratTFC(i, j, k) * wvrt
                  else
                    buoy = buoy - dt * rhoStratTFC(i, j, k) / rho &
                        &* bvsStratTFC(i, j, k) * wvrt
                  end if
                else
                  buoy = - g_ndim * rhop / rho - dt * (rhoStrat(k) / rho &
                      &* PStrat(k) / PStrat_0(k) * bvsStrat(k) * wvrt + g_ndim &
                      &/ rho * rhoStrat(k) / Pstrat_0(k) * heat_flc)
                end if

                var%rhop(i, j, k) = - buoy * rho / g_ndim
              end do
            end do
          end do
        else
          stop 'int_mod unknown'
        end if
      else
        stop 'upd_mod unknown'
      end if
    else if(upd_var == "P") then
      ! SK: update of mass-weighted potential temp. for compressible model
      if(upd_mod /= "tot") then
        stop 'upd_mod unknown'
      else
        if(int_mod /= "expl") then
          stop 'ERROR: wrong int_mod for upd_mod = lhs'
        end if

        if(TurbScheme .or. rayTracer) then
          call calculate_heating(var, flux, heat)
        else
          heat = 0.0
        end if

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              ! using JPu as the carrier flux
              fL = flux%P(i - 1, j, k, 1) ! mass flux accros left cell edge
              fR = flux%P(i, j, k, 1) ! right
              gB = flux%P(i, j - 1, k, 2) ! backward
              gF = flux%P(i, j, k, 2) ! forward
              hD = flux%P(i, j, k - 1, 3) ! downward
              hU = flux%P(i, j, k, 3) ! upward

              ! convective part
              fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

              if(topography) then
                fluxDiff = fluxDiff / jac(i, j, k)
              end if

              ! F(phi)
              F = - fluxDiff - heat(i, j, k)

              ! update: q(m-1) -> q(m)
              q(i, j, k) = dt * F + alphaRK(m) * q(i, j, k)

              ! update density
              var%P(i, j, k) = var%P(i, j, k) + betaRK(m) * q(i, j, k)
            end do
          end do
        end do
      end if
    else
      stop 'upd_var unknown'
    end if

  end subroutine massUpdate

  !-----------------------------------------------------------------------
  subroutine tracerUpdate(var, flux, tracerforce, dt, q, m)

    ! in/out variables
    type(var_type), intent(inout) :: var

    type(flux_type), intent(in) :: flux

    type(tracerForceType), dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz &
        &+ nbz), intent(in) :: tracerforce

    real, intent(in) :: dt
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: q

    integer, intent(in) :: m

    ! local variables
    integer :: i, j, k, l
    real :: fL, fR ! flux Left/Right
    real :: gB, gF ! flux Backward/Forward
    real :: hD, hU ! flux Downward/Upward
    real :: fluxDiff ! convective part
    real :: F ! F(phi)
    real :: rho
    real :: forcetracer ! additional rhs terms from gw parameterization

    if(correctDivError) then
      print *, 'ERROR: correction divergence error not allowed'
      stop
    end if

    ! init q
    if(m == 1) q = 0.

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          fL = flux%chi(i - 1, j, k, 1) ! mass flux accros left cell edge
          fR = flux%chi(i, j, k, 1) ! right
          gB = flux%chi(i, j - 1, k, 2) ! backward
          gF = flux%chi(i, j, k, 2) ! forward
          hD = flux%chi(i, j, k - 1, 3) ! downward
          hU = flux%chi(i, j, k, 3) ! upward

          if(topography) then
            rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)
          else
            rho = var%rho(i, j, k) + rhoStrat(k)
          end if

          ! convective part, advection due to wind
          fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

          if(topography) then
            fluxDiff = fluxDiff / jac(i, j, k)
          end if

          ! F(phi)
          F = - fluxDiff
          if(rayTracer) then
            ! additional rhs terms due to gw and turbulence impact
            ! saved in forcetracer
            forcetracer = 0.0

            ! include leading order gw tracer flux convergence
            if(include_trfrc_lo) then
              forcetracer = forcetracer + tracerforce(i, j, k)%loforce%total
            end if

            ! include next-order gw tracer flux convergence
            if(include_trfrc_no) then
              forcetracer = forcetracer + rho * tracerforce(i, j, &
                  &k)%noforce%total
            end if

            ! include diffusive mixing of tracer
            if(include_trfrc_mix) then
              forcetracer = forcetracer - tracerforce(i, j, k)%mixingGW%total
            end if

            F = F - forcetracer ! rho *
          end if

          if(dens_relax) then
            stop "update.f90: dens_relax not implemented in tracerUpdate"
          end if

          ! update: q(m-1) -> q(m)
          q(i, j, k) = dt * F + alphaRK(m) * q(i, j, k)

          ! update density
          var%chi(i, j, k) = var%chi(i, j, k) + betaRK(m) * q(i, j, k)

        end do
      end do
    end do

  end subroutine tracerUpdate

  !-----------------------------------------------------------------------

  subroutine timeUpdate(time, dt, q, m)

    implicit none
    ! in/out variables
    real, intent(inout) :: time
    real, intent(inout) :: q
    real, intent(in) :: dt
    integer, intent(in) :: m

    ! init q
    !if(m == 1) q = 0. !alphaRK(1) allways 0 ?!

    ! update: q(m-1) -> q(m)
    q = dt + alphaRK(m) * q

    ! update time
    time = time + betaRK(m) * q

  end subroutine timeUpdate

  !!$  subroutine iceUpdate_apb(var, flux, source, dt, q, m, update_type, ray_cloud)
  !!$
  !!$    ! in/out variables
  !!$    type(var_type), intent(inout) :: var
  !!$    type(flux_type), intent(in) :: flux
  !!$    type(var_type), intent(in) :: source
  !!$    real, intent(in) :: dt
  !!$    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVarIce), &
  !!$         &intent(inout) :: q
  !!$    type(ice_rayType2), dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, nscx, nscy), &
  !!$         &intent(inout) :: ray_cloud
  !!$    integer, intent(in) :: m
  !!$    character(len = 3), intent(in) :: update_type
  !!$
  !!$    ! local variables
  !!$    integer :: i, j, k, l
  !!$    real :: fL, fR ! flux Left/Right
  !!$    real :: gB, gF ! flux Backward/Forward
  !!$    real :: hD, hU ! flux Downward/Upward
  !!$    real :: fluxDiff ! convective part
  !!$    real :: F ! F(phi)
  !!$
  !!$    integer :: ii, iVar
  !!$
  !!$    if(compute_cloudcover) then
  !!$
  !!$       ! init q
  !!$       if(m == 1) then
  !!$          ray_cloud%qNi = 0.
  !!$          ray_cloud%qQi = 0.
  !!$          ray_cloud%qQv = 0.
  !!$       end if
  !!$
  !!$       ! update: q(m-1) -> q(m)
  !!$       ray_cloud%qNi = dt * ray_cloud%tNi + alphaRK(m) * ray_cloud%qNi
  !!$       ray_cloud%qQi = dt * ray_cloud%tQi + alphaRK(m) * ray_cloud%qQi
  !!$       ray_cloud%qQv = dt * ray_cloud%tQv + alphaRK(m) * ray_cloud%qQv
  !!$
  !!$       ! update fields
  !!$       ray_cloud%Ni = ray_cloud%Ni + betaRK(m) * ray_cloud%qNi
  !!$       ray_cloud%Qi = ray_cloud%Qi + betaRK(m) * ray_cloud%qQi
  !!$       ray_cloud%Qv = ray_cloud%Qv + betaRK(m) * ray_cloud%qQv
  !!$
  !!$    else
  !!$
  !!$       ! init q
  !!$       if(m == 1) q = 0.
  !!$
  !!$       do iVar = 1, nVarIce
  !!$
  !!$          do k = 1, nz
  !!$             do j = 1, ny
  !!$                do i = 1, nx
  !!$
  !!$                   if(update_type .eq. 'ADV' .or. update_type .eq. 'BOT') then
  !!$
  !!$                      fL = flux%ICE(i - 1, j, k, 1, iVar) ! flux accros left cell edge
  !!$                      fR = flux%ICE(i, j, k, 1, iVar) ! right
  !!$                      gB = flux%ICE(i, j - 1, k, 2, iVar) ! backward
  !!$                      gF = flux%ICE(i, j, k, 2, iVar) ! forward
  !!$                      hD = flux%ICE(i, j, k - 1, 3, iVar) ! downward
  !!$                      hU = flux%ICE(i, j, k, 3, iVar) ! upward
  !!$
  !!$                      ! convective part
  !!$                      fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz
  !!$
  !!$                      ! TFC FJ
  !!$                      ! Adjust mass flux divergence.
  !!$                      if(topography) then
  !!$                         fluxDiff = fluxDiff / jac(i, j, k)
  !!$                      end if
  !!$
  !!$                      if(update_type .eq. 'BOT') then
  !!$                         ! F(phi)
  !!$                         F = - fluxDiff + source%ICE(i, j, k, iVar)
  !!$                      elseif(update_type .eq. 'ADV') then
  !!$                         F = - fluxDiff
  !!$                      elseif(update_type .eq. 'PHY') then
  !!$                         F = source%ICE(i, j, k, iVar)
  !!$                      else
  !!$                         print *, 'wrong update_type in iceUpdate_apb'
  !!$                         stop
  !!$                      end if
  !!$
  !!$                      ! update: q(m-1) -> q(m)
  !!$                      q(i, j, k, iVar) = dt * F + alphaRK(m) * q(i, j, k, iVar)
  !!$
  !!$                      ! update fields
  !!$                      var%ICE(i, j, k, iVar) = var%ICE(i, j, k, iVar) + betaRK(m) &
  !!$                           &* q(i, j, k, iVar)
  !!$
  !!$                   end do !i
  !!$                end do !j
  !!$             end do !k
  !!$
  !!$          end do !ii
  !!$       end if !compute cloudcover
  !!$     end subroutine iceUpdate_apb
  subroutine iceUpdate_apb(var, flux, source, dt, q, m, update_type, ray_cloud)

    ! in/out variables
    type(var_type), intent(inout) :: var
    type(flux_type), intent(in) :: flux
    type(var_type), intent(in) :: source
    real, intent(in) :: dt
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVarIce), &
        &intent(inout) :: q
    type(ice_rayType2), dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, nscx, nscy), &
        &intent(inout) :: ray_cloud
    integer, intent(in) :: m
    character(len = 3), intent(in) :: update_type

    ! local variables
    integer :: i, j, k, l
    real :: fL, fR ! flux Left/Right
    real :: gB, gF ! flux Backward/Forward
    real :: hD, hU ! flux Downward/Upward
    real :: fluxDiff ! convective part
    real :: F !
    integer :: ii, iVar

    if(compute_cloudcover) then

      ! init q
      if(m == 1) then
        ray_cloud%qNi = 0.
        ray_cloud%qQi = 0.
        ray_cloud%qQv = 0.
      end if

      ! update: q(m-1) -> q(m)
      ray_cloud%qNi = dt * ray_cloud%tNi + alphaRK(m) * ray_cloud%qNi
      ray_cloud%qQi = dt * ray_cloud%tQi + alphaRK(m) * ray_cloud%qQi
      ray_cloud%qQv = dt * ray_cloud%tQv + alphaRK(m) * ray_cloud%qQv

      ! update fields
      ray_cloud%Ni = ray_cloud%Ni + betaRK(m) * ray_cloud%qNi
      ray_cloud%Qi = ray_cloud%Qi + betaRK(m) * ray_cloud%qQi
      ray_cloud%Qv = ray_cloud%Qv + betaRK(m) * ray_cloud%qQv

    else

      ! init q
      if(m == 1) q = 0.

      do iVar = 1, nVarIce

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx

              if(update_type .eq. 'ADV' .or. update_type .eq. 'BOT') then

                fL = flux%ICE(i - 1, j, k, 1, iVar) ! flux accros left cell edge
                fR = flux%ICE(i, j, k, 1, iVar) ! right
                gB = flux%ICE(i, j - 1, k, 2, iVar) ! backward
                gF = flux%ICE(i, j, k, 2, iVar) ! forward
                hD = flux%ICE(i, j, k - 1, 3, iVar) ! downward
                hU = flux%ICE(i, j, k, 3, iVar) ! upward

                ! convective part
                fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

                ! TFC FJ
                ! Adjust mass flux divergence.
                if(topography) then
                  fluxDiff = fluxDiff / jac(i, j, k)
                end if

              end if

              if(update_type .eq. 'BOT') then
                ! F(phi)
                F = - fluxDiff + source%ICE(i, j, k, iVar)
              elseif(update_type .eq. 'ADV') then
                F = - fluxDiff
              elseif(update_type .eq. 'PHY') then
                F = source%ICE(i, j, k, iVar)
              else
                print *, 'wrong update_type in iceUpdate_apb'
                stop
              end if

              ! update: q(m-1) -> q(m)
              q(i, j, k, iVar) = dt * F + alphaRK(m) * q(i, j, k, iVar)

              ! update fields
              var%ICE(i, j, k, iVar) = var%ICE(i, j, k, iVar) + betaRK(m) &
                  &* q(i, j, k, iVar)

            end do !i
          end do !j
        end do !k

      end do !ii
    end if !compute cloudcover
  end subroutine iceUpdate_apb

  subroutine timestep(var, dt, errFlag)
    !---------------------------------------------
    ! compute time step from stability criteria:
    ! 1) CFL criterion for advection
    ! 2) von Neumann cirterion for dissipation
    ! 3) set maximum time step
    ! 4) gravity-wave group condition: c_g_x * dt < dx, ...dy,dz
    !---------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var
    real, intent(out) :: dt
    logical, intent(out) :: errFlag

    ! locals
    real :: uMax, vMax, wMax
    real :: dtConv, dtVisc, dtCond, dtWKB
    real :: dtConv_loc, dtVisc_loc, dtCond_loc, dtWKB_loc
    real :: dtMax
    real :: dtWave

    ! local integer
    integer :: i, j, k

    errFlag = .false.

    !-------------------------------------------
    !              Fixed time step
    !-------------------------------------------

    if(tStepChoice == "fix") then

      dt = dtMax_dim / tRef
      errFlag = .false.

      if(master) then
        write(*, fmt = "(a25,es15.1,a8)") "dt = dtFix = ", dt * tRef, "seconds"
      end if

    else

      !-------------------------------------------
      !           Variable time step
      !-------------------------------------------

      select case(model)

      case("Boussinesq", "pseudo_incompressible", "compressible")

        !----------------------------
        !   Full model time step
        !----------------------------

        !----------------------
        !     CFL condition
        !----------------------
        uMax = maxval(abs(var%u(1:nx, 1:ny, 1:nz))) + small
        vMax = maxval(abs(var%v(1:nx, 1:ny, 1:nz))) + small
        wMax = maxval(abs(var%w(1:nx, 1:ny, 1:nz))) + small

        dtConv_loc = cfl * min(dx / uMax, dy / vMax, dz / wMax)

        if(topography) then
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                dtConv_loc = min(dtConv_loc, cfl * jac(i, j, k) * dz &
                    &/ (abs(0.5 * (vertWindTFC(i, j, k, var) + vertWindTFC(i, &
                    &j, k - 1, var))) + small))
              end do
            end do
          end do
        end if

        ! find global minimum

        call mpi_reduce(dtConv_loc, dtConv, 1, mpi_double_precision, mpi_min, &
            &root, comm, ierror)

        call mpi_bcast(dtConv, 1, mpi_double_precision, root, comm, ierror)

        !---------------------------
        !   von Neumann condition
        !----------------------------

        dtVisc = 0.5 * min(dx ** 2, dy ** 2, dz ** 2) * Re
        dtCond = 0.5 * min(dx ** 2, dy ** 2, dz ** 2) / mu_conduct

        if(topography) then
          dtVisc_loc = dtVisc
          dtCond_loc = dtCond
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                dtVisc_loc = min(dtVisc_loc, 0.5 * (jac(i, j, k) * dz) ** 2.0 &
                    &* Re)
                dtCond_loc = min(dtCond_loc, 0.5 * (jac(i, j, k) * dz) ** 2.0 &
                    &/ mu_conduct)
              end do
            end do
          end do
          call mpi_reduce(dtVisc_loc, dtVisc, 1, mpi_double_precision, &
              &mpi_min, root, comm, ierror)
          call mpi_reduce(dtCond_loc, dtCond, 1, mpi_double_precision, &
              &mpi_min, root, comm, ierror)
          call mpi_bcast(dtVisc, 1, mpi_double_precision, root, comm, ierror)
          call mpi_bcast(dtCond, 1, mpi_double_precision, root, comm, ierror)
        end if

        !----------------------------
        !    Maximal time step
        !----------------------------
        dtMax = dtMax_dim / tRef

        !------------------------------------
        !    Gravity wave time period
        !------------------------------------

        dtWave = 1. / (NN + small)

        !------------------------------------
        !     WKB "CFL" criterion
        !------------------------------------

        if(raytracer) then
          dtWKB_loc = dz / (cgz_max + small)

          if(topography) then
            do k = 0, nz
              do j = 1, ny
                do i = 1, nx
                  dtWKB_loc = min(dtWKB_loc, jac(i, j, k) * dz &
                      &/ (cgz_max_tfc(i, j, k) + small))
                end do
              end do
            end do
          end if

          if(sizeX > 1) dtWKB_loc = min(dtWKB_loc, dx / (cgx_max + small))
          if(sizeY > 1) dtWKB_loc = min(dtWKB_loc, dy / (cgy_max + small))

          dtWKB_loc = cfl_wave * dtWKB_loc

          ! find global minimum

          call mpi_reduce(dtWKB_loc, dtWKB, 1, mpi_double_precision, mpi_min, &
              &root, comm, ierror)

          call mpi_bcast(dtWKB, 1, mpi_double_precision, root, comm, ierror)
        end if

        !-------------------------------
        !        Make your choice
        !-------------------------------

        if(dtWave_on .and. timeScheme /= 'semiimplicit') then
          dt = min(dtVisc, dtCond, dtConv, dtMax, dtWave)
        else
          dt = min(dtVisc, dtCond, dtConv, dtMax)
        end if

        if(raytracer) dt = min(dt, dtWKB)

        !-----------------------------------------
        !     Inform on time step restrictions
        !-----------------------------------------

        if(master) then

          write(*, fmt = "(a25,es15.1,a8)") "dtVisc =", dtVisc * tRef, "seconds"
          write(*, fmt = "(a25,es15.1,a8)") "dtCond =", dtCond * tRef, "seconds"
          write(*, fmt = "(a25,es15.1,a8)") "dtConv =", dtConv * tRef, "seconds"
          write(*, fmt = "(a25,es15.1,a8)") "dtMax =", dtMax * tRef, "seconds"
          write(*, fmt = "(a25,es15.1,a8)") "dtWave =", dtWave * tRef, "seconds"
          if(raytracer) then
            write(*, fmt = "(a25,es15.1,a8)") "dtWKB =", dtWKB * tRef, "seconds"
          end if
          print *, ""

          if(dt == dtMax) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtMax = ", dt * tRef, &
                &"seconds"
          else if(dt == dtConv) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtConv = ", dt * tRef, &
                &"seconds"
          else if(dt == dtVisc) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtVisc = ", dt * tRef, &
                &"seconds"
          else if(dt == dtCond) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtCond = ", dt * tRef, &
                &"seconds"
          else if(dt == dtWave) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtWave = ", dt * tRef, &
                &"seconds"
          else if(dt == dtWKB) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtWKB =", dt * tRef, &
                &"seconds"
          else
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = ????? = ", dt * tRef, &
                &"seconds"
          end if
          print *, ""

        end if

      case default
        stop "timestep: unknown case model."
      end select ! WKB / full model

    end if

    if(dt * tRef < dtMin_dim) errFlag = .true.

  end subroutine timestep

  !-------------------------------------------------------------------------

  subroutine init_update
    !--------------------------------------
    ! allocate variables for update module
    !--------------------------------------

    ! local variables
    integer :: allocstat

    ! number of equations
    nxyz = nx * ny * nz

    ! number of non-zeros in the equations system
    nnz = 7 * nxyz - 2 * nx * ny
    ! xxx: verify this value when changing boundary conditions

  end subroutine init_update

  !---------------------------------------------------------------------

  subroutine CoefDySma_update(var)
    !--------------------------------------
    ! calculate the Coefficient for Dynamic Smagorinsky Scheme
    !--------------------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var

    ! more variables
    real :: delta_hs, delta_vs
    real :: uL, uR, uB, uF, uD, uU
    real :: vL, vR, vB, vF, vD, vU
    real :: wL, wR, wB, wF, wD, wU
    real :: du_dx, du_dy, du_dz
    real :: dv_dx, dv_dy, dv_dz
    real :: dw_dx, dw_dy, dw_dz

    ! allocatable fields
    real, dimension(:, :, :, :, :), allocatable :: Sij, Lij, Mij
    real, dimension(:, :, :), allocatable :: S_norm

    real, dimension(:, :, :, :, :), allocatable :: uiuj_smth, S_Sij_smth, &
        &Sij_smth
    real, dimension(:, :, :, :), allocatable :: ui_smth
    real, dimension(:, :, :), allocatable :: Sn_smth
    real, dimension(:, :, :), allocatable :: LijMij_smth, MijMij_smth
    real, dimension(:, :, :), allocatable :: CS2_DySma

    integer :: allocstat
    integer :: i, j, k
    integer :: iw, jw

    integer :: smth_npts1_DySma, smth_npts2_DySma
    parameter(smth_npts1_DySma = 1) ! revised by JW (20160824)
    parameter(smth_npts2_DySma = 2) ! revised by JW (20160824)

    ! Allocate local fields
    allocate(Sij(1:nx, 1:ny, 1:nz, 1:3, 1:3), stat = allocstat)
    if(allocstat /= 0) stop "CoefDySma_update:alloc failed"

    allocate(Lij(1:nx, 1:ny, 1:nz, 1:3, 1:3), stat = allocstat)
    if(allocstat /= 0) stop "CoefDySma_update:alloc failed"

    allocate(Mij(1:nx, 1:ny, 1:nz, 1:3, 1:3), stat = allocstat)
    if(allocstat /= 0) stop "CoefDySma_update:alloc failed"

    allocate(S_norm(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "CoefDySma_update:alloc failed"

    allocate(uiuj_smth(1:nx, 1:ny, 1:nz, 1:3, 1:3), stat = allocstat)
    if(allocstat /= 0) stop "CoefDySma_update:alloc failed"

    allocate(S_Sij_smth(1:nx, 1:ny, 1:nz, 1:3, 1:3), stat = allocstat)
    if(allocstat /= 0) stop "CoefDySma_update:alloc failed"

    allocate(Sij_smth(1:nx, 1:ny, 1:nz, 1:3, 1:3), stat = allocstat)
    if(allocstat /= 0) stop "CoefDySma_update:alloc failed"

    allocate(ui_smth(1:nx, 1:ny, 1:nz, 1:3), stat = allocstat)
    if(allocstat /= 0) stop "CoefDySma_update:alloc failed"

    allocate(Sn_smth(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "CoefDySma_update:alloc failed"

    allocate(LijMij_smth(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "CoefDySma_update:alloc failed"

    allocate(MijMij_smth(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "CoefDySma_update:alloc failed"

    allocate(CS2_DySma(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "CoefDySma_update:alloc failed"

    ! calculate delta

    if(TurbScheme) then
      if(ny == 1 .and. nx == 1) then
        stop 'ERROR: turbulence assumes either nx > 1 or ny > 1'
      else
        if(nx == 1) then
          delta_hs = dy ** 2 ! 2D problems in y and z
        else if(ny == 1) then
          delta_hs = dx ** 2 ! 2D problems in x and z
        else
          delta_hs = dx * dy ! 3D problems

          if(dx / dy > 10.) then
            print *, 'WARNING: dx/dy > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          elseif(dy / dx > 10.) then
            print *, 'WARNING: dy/dx > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          end if
        end if

        delta_vs = dz ** 2
      end if
    end if

    ! calculate S_ij

    !---------------------------------
    !         Loop over field
    !---------------------------------

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

          uL = var%u(i - 1, j, k)
          uR = var%u(i, j, k)

          ! UA not sure whether this averaging is necessary.
          ! Compare, e.g. the viscous fluxes. There it is not done.

          ! replied by JW (20160824): For now, there is no revision on
          ! this part, since it is a relatively minor issue at the moment.

          uB = 0.5 * (var%u(i - 1, j - 1, k) + var%u(i, j - 1, k))
          uF = 0.5 * (var%u(i - 1, j + 1, k) + var%u(i, j + 1, k))
          uD = 0.5 * (var%u(i - 1, j, k - 1) + var%u(i, j, k - 1))
          uU = 0.5 * (var%u(i - 1, j, k + 1) + var%u(i, j, k + 1))

          vL = 0.5 * (var%v(i - 1, j, k) + var%v(i - 1, j - 1, k))
          vR = 0.5 * (var%v(i + 1, j, k) + var%v(i + 1, j - 1, k))
          vB = var%v(i, j - 1, k)
          vF = var%v(i, j, k)
          vD = 0.5 * (var%v(i, j, k - 1) + var%v(i, j - 1, k - 1))
          vU = 0.5 * (var%v(i, j, k + 1) + var%v(i, j - 1, k + 1))

          if(topography) then
            wL = 0.5 * (vertWindTFC(i - 1, j, k - 1, var) + vertWindTFC(i - 1, &
                &j, k, var))
            wR = 0.5 * (vertWindTFC(i + 1, j, k - 1, var) + vertWindTFC(i + 1, &
                &j, k, var))
            wB = 0.5 * (vertWindTFC(i, j - 1, k - 1, var) + vertWindTFC(i, j &
                &- 1, k, var))
            wF = 0.5 * (vertWindTFC(i, j + 1, k - 1, var) + vertWindTFC(i, j &
                &+ 1, k, var))
            wD = vertWindTFC(i, j, k - 1, var)
            wU = vertWindTFC(i, j, k, var)
          else
            wL = 0.5 * (var%w(i - 1, j, k - 1) + var%w(i - 1, j, k))
            wR = 0.5 * (var%w(i + 1, j, k - 1) + var%w(i + 1, j, k))
            wB = 0.5 * (var%w(i, j - 1, k - 1) + var%w(i, j - 1, k))
            wF = 0.5 * (var%w(i, j + 1, k - 1) + var%w(i, j + 1, k))
            wD = var%w(i, j, k - 1)
            wU = var%w(i, j, k)
          end if

          if(topography) then
            du_dx = (uR - uL) / dx + met(i, j, k, 1, 3) * (uU - uD) / (2.0 * dz)
            du_dy = (uF - uB) / (2.0 * dy) + met(i, j, k, 2, 3) * (uU - uD) &
                &/ (2.0 * dz)
            du_dz = (uU - uD) / (2.0 * dz) / jac(i, j, k)

            dv_dx = (vR - vL) / (2.0 * dx) + met(i, j, k, 1, 3) * (vU - vD) &
                &/ (2.0 * dz)
            dv_dy = (vF - vB) / dy + met(i, j, k, 2, 3) * (vU - vD) / (2.0 * dz)
            dv_dz = (vU - vD) / (2.0 * dz) / jac(i, j, k)

            dw_dx = (wR - wL) / (2.0 * dx) + met(i, j, k, 1, 3) * (wU - wD) / dz
            dw_dy = (wF - wB) / (2.0 * dy) + met(i, j, k, 2, 3) * (wU - wD) / dz
            dw_dz = (wU - wD) / dz / jac(i, j, k)
          else
            du_dx = (uR - uL) / dx

            !UA in case without averaging, no factor 2!

            ! replied by JW (20160824): For now, there is no revision on
            ! this part, since it is a relatively minor issue at the moment.

            du_dy = (uF - uB) / (2.0 * dy)
            du_dz = (uU - uD) / (2.0 * dz)

            dv_dx = (vR - vL) / (2.0 * dx)
            dv_dy = (vF - vB) / dy
            dv_dz = (vU - vD) / (2.0 * dz)

            dw_dx = (wR - wL) / (2.0 * dx)
            dw_dy = (wF - wB) / (2.0 * dy)
            dw_dz = (wU - wD) / dz
          end if

          Sij(i, j, k, 1, 1) = 0.5 * (du_dx + du_dx)
          Sij(i, j, k, 1, 2) = 0.5 * (du_dy + dv_dx)
          Sij(i, j, k, 1, 3) = 0.5 * (du_dz + dw_dx)
          Sij(i, j, k, 2, 1) = 0.5 * (dv_dx + du_dy)
          Sij(i, j, k, 2, 2) = 0.5 * (dv_dy + dv_dy)
          Sij(i, j, k, 2, 3) = 0.5 * (dv_dz + dw_dy)
          Sij(i, j, k, 3, 1) = 0.5 * (dw_dx + du_dz)
          Sij(i, j, k, 3, 2) = 0.5 * (dw_dy + dv_dz)
          Sij(i, j, k, 3, 3) = 0.5 * (dw_dz + dw_dz)

          Sij(i, j, k, 1, 1) = Sij(i, j, k, 1, 1) - (du_dx + dv_dy + dw_dz) &
              &/ 3.0
          Sij(i, j, k, 2, 2) = Sij(i, j, k, 2, 2) - (du_dx + dv_dy + dw_dz) &
              &/ 3.0
          Sij(i, j, k, 3, 3) = Sij(i, j, k, 3, 3) - (du_dx + dv_dy + dw_dz) &
              &/ 3.0

          S_norm(i, j, k) = 0.0

          do jw = 1, 3
            do iw = 1, 3
              S_norm(i, j, k) = S_norm(i, j, k) + Sij(i, j, k, iw, jw) &
                  &* Sij(i, j, k, iw, jw)
            end do
          end do

          S_norm(i, j, k) = S_norm(i, j, k) * 2.0
          S_norm(i, j, k) = sqrt(S_norm(i, j, k))
        end do
      end do
    end do

    !---------------------------------
    !         Loop over field
    !         Prepare the data (before smoothing at the 1st step)
    !---------------------------------

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ui_smth(i, j, k, 1) = 0.5 * (var%u(i, j, k) + var%u(i - 1, j, k))
          ui_smth(i, j, k, 2) = 0.5 * (var%v(i, j, k) + var%v(i, j - 1, k))
          if(topography) then
            ui_smth(i, j, k, 3) = 0.5 * (vertWindTFC(i, j, k, var) &
                &+ vertWindTFC(i, j, k - 1, var))
          else
            ui_smth(i, j, k, 3) = 0.5 * (var%w(i, j, k) + var%w(i, j, k - 1))
          end if

          Sn_smth(i, j, k) = S_norm(i, j, k)

          do jw = 1, 3
            do iw = 1, 3
              uiuj_smth(i, j, k, iw, jw) = ui_smth(i, j, k, iw) * ui_smth(i, &
                  &j, k, jw)

              Sij_smth(i, j, k, iw, jw) = Sij(i, j, k, iw, jw)

              S_Sij_smth(i, j, k, iw, jw) = S_norm(i, j, k) * Sij(i, j, k, iw, &
                  &jw)
            end do
          end do
        end do
      end do
    end do

    !---------------------------------
    !         Smoothing at the 1st step
    !---------------------------------

    if(ny .eq. 1) then
      do iw = 1, 3
        call Var3DSmthDySma(ui_smth(1:nx, 1:ny, 1:nz, iw), smth_npts1_DySma, &
            &"XZ_local_smth")
      end do

      call Var3DSmthDySma(Sn_smth(1:nx, 1:ny, 1:nz), smth_npts1_DySma, &
          &"XZ_local_smth")

      do jw = 1, 3
        do iw = 1, 3
          call Var3DSmthDySma(uiuj_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              &smth_npts1_DySma, "XZ_local_smth")
          call Var3DSmthDySma(Sij_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              &smth_npts1_DySma, "XZ_local_smth")
          call Var3DSmthDySma(S_Sij_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              &smth_npts1_DySma, "XZ_local_smth")
        end do
      end do
    elseif(nx .eq. 1) then
      do iw = 1, 3
        call Var3DSmthDySma(ui_smth(1:nx, 1:ny, 1:nz, iw), smth_npts1_DySma, &
            &"YZ_local_smth")
      end do

      call Var3DSmthDySma(Sn_smth(1:nx, 1:ny, 1:nz), smth_npts1_DySma, &
          &"YZ_local_smth")

      do jw = 1, 3
        do iw = 1, 3
          call Var3DSmthDySma(uiuj_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              &smth_npts1_DySma, "YZ_local_smth")
          call Var3DSmthDySma(Sij_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              &smth_npts1_DySma, "YZ_local_smth")
          call Var3DSmthDySma(S_Sij_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              &smth_npts1_DySma, "YZ_local_smth")
        end do
      end do
    else
      do iw = 1, 3
        call Var3DSmthDySma(ui_smth(1:nx, 1:ny, 1:nz, iw), smth_npts1_DySma, &
            &"XYZ_local_smth")
      end do

      call Var3DSmthDySma(Sn_smth(1:nx, 1:ny, 1:nz), smth_npts1_DySma, &
          &"XYZ_local_smth")

      do jw = 1, 3
        do iw = 1, 3
          call Var3DSmthDySma(uiuj_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              &smth_npts1_DySma, "XYZ_local_smth")
          call Var3DSmthDySma(Sij_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              &smth_npts1_DySma, "XYZ_local_smth")
          call Var3DSmthDySma(S_Sij_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              &smth_npts1_DySma, "XYZ_local_smth")
        end do
      end do
    end if

    !---------------------------------
    !         Loop over field
    !         (Smoothing at the 2nd step)
    !---------------------------------

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          LijMij_smth(i, j, k) = 0.0
          MijMij_smth(i, j, k) = 0.0

          do jw = 1, 3
            do iw = 1, 3
              Lij(i, j, k, iw, jw) = uiuj_smth(i, j, k, iw, jw) - ui_smth(i, &
                  &j, k, iw) * ui_smth(i, j, k, jw)

              Mij(i, j, k, iw, jw) = S_Sij_smth(i, j, k, iw, jw) - (2.0 &
                  &* smth_npts1_DySma + 1.0) ** 2 * Sn_smth(i, j, k) &
                  &* Sij_smth(i, j, k, iw, jw)

              ! allow for grid anisotropy

              if(iw == 3 .or. jw == 3) then
                if(topography) then
                  Mij(i, j, k, iw, jw) = Mij(i, j, k, iw, jw) * jac(i, j, k) &
                      &** 2.0 * delta_vs
                else
                  Mij(i, j, k, iw, jw) = Mij(i, j, k, iw, jw) * delta_vs
                end if
              else
                Mij(i, j, k, iw, jw) = Mij(i, j, k, iw, jw) * delta_hs
              end if

              LijMij_smth(i, j, k) = LijMij_smth(i, j, k) + Lij(i, j, k, iw, &
                  &jw) * Mij(i, j, k, iw, jw)

              MijMij_smth(i, j, k) = MijMij_smth(i, j, k) + Mij(i, j, k, iw, &
                  &jw) * Mij(i, j, k, iw, jw)
            end do
          end do
        end do
      end do
    end do

    if(ny .eq. 1) then
      call Var3DSmthDySma(LijMij_smth(1:nx, 1:ny, 1:nz), smth_npts2_DySma, &
          &"XZ_local_smth")
      call Var3DSmthDySma(MijMij_smth(1:nx, 1:ny, 1:nz), smth_npts2_DySma, &
          &"XZ_local_smth")
    elseif(nx .eq. 1) then
      call Var3DSmthDySma(LijMij_smth(1:nx, 1:ny, 1:nz), smth_npts2_DySma, &
          &"YZ_local_smth")
      call Var3DSmthDySma(MijMij_smth(1:nx, 1:ny, 1:nz), smth_npts2_DySma, &
          &"YZ_local_smth")
    else
      call Var3DSmthDySma(LijMij_smth(1:nx, 1:ny, 1:nz), smth_npts2_DySma, &
          &"XYZ_local_smth")
      call Var3DSmthDySma(MijMij_smth(1:nx, 1:ny, 1:nz), smth_npts2_DySma, &
          &"XYZ_local_smth")
    end if

    !---------------------------------
    !         Get the final results
    !---------------------------------

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

          if(MijMij_smth(i, j, k) /= 0.) then
            CS2_DySma(i, j, k) = 0.5 * LijMij_smth(i, j, k) / MijMij_smth(i, &
                &j, k)
          else
            CS2_DySma(i, j, k) = 0.
          end if

          if(CS2_DySma(i, j, k) < 0.0) then
            CS2_DySma(i, j, k) = 0.0
          end if

          var%DSC(i, j, k) = CS2_DySma(i, j, k) * S_norm(i, j, k)
        end do
      end do
    end do

    ! *** set values for the ghost celss ***
    call setHaloAndBoundary(var%DSC(:, :, :), nbx, nby, nbz)

    ! deallocate local fields
    deallocate(Sij, stat = allocstat); if(allocstat /= 0) stop &
        &"update.f90:dealloc failed"
    deallocate(Lij, stat = allocstat); if(allocstat /= 0) stop &
        &"update.f90:dealloc failed"
    deallocate(Mij, stat = allocstat); if(allocstat /= 0) stop &
        &"update.f90:dealloc failed"
    deallocate(S_norm, stat = allocstat); if(allocstat /= 0) stop &
        &"update.f90:dealloc failed"
    deallocate(uiuj_smth, stat = allocstat); if(allocstat /= 0) stop &
        &"update.f90:dealloc failed"
    deallocate(S_Sij_smth, stat = allocstat); if(allocstat /= 0) stop &
        &"update.f90:dealloc failed"
    deallocate(Sij_smth, stat = allocstat); if(allocstat /= 0) stop &
        &"update.f90:dealloc failed"
    deallocate(ui_smth, stat = allocstat); if(allocstat /= 0) stop &
        &"update.f90:dealloc failed"
    deallocate(Sn_smth, stat = allocstat); if(allocstat /= 0) stop &
        &"update.f90:dealloc failed"
    deallocate(LijMij_smth, stat = allocstat); if(allocstat /= 0) stop &
        &"update.f90:dealloc failed"
    deallocate(MijMij_smth, stat = allocstat); if(allocstat /= 0) stop &
        &"update.f90:dealloc failed"
    deallocate(CS2_DySma, stat = allocstat); if(allocstat /= 0) stop &
        &"update.f90:dealloc failed"

    return

  end subroutine CoefDySma_update

  subroutine Var3DSmthDySma(var3D_DySma, nsmth_DySma, homog_dir_DySma)
    !--------------------------------------
    ! calculate the Coefficient for Dynamic Smagorinsky Scheme
    !--------------------------------------

    ! in/out variables
    real, dimension(1:nx, 1:ny, 1:nz), intent(inout) :: var3D_DySma
    integer, intent(in) :: nsmth_DySma

    character(len = *), intent(in) :: homog_dir_DySma

    ! allocatable fields
    real, dimension(:, :, :), allocatable :: var3D_DySma_Extend

    integer :: allocstat
    integer :: i, j, k
    integer :: kmin, kmax, nsmthv
    integer :: nsmthall, ismth, jsmth, ksmth
    integer :: i0, j0

    allocate(var3D_DySma_Extend((0 - nsmth_DySma):(nx + nsmth_DySma), (0 &
        &- nsmth_DySma):(ny + nsmth_DySma), (0 - nsmth_DySma):(nz &
        &+ nsmth_DySma)), stat = allocstat)
    if(allocstat /= 0) stop "Var3DSmthDySma:alloc failed"

    ! set the values for var3D_DySma_Extend

    var3D_DySma_Extend(1:nx, 1:ny, 1:nz) = var3D_DySma(1:nx, 1:ny, 1:nz)

    call setHaloAndBoundary(var3D_DySma_Extend(:, :, :), nsmth_DySma, &
        &nsmth_DySma, nsmth_DySma)

    ! start to do the smoothing

    i0 = is + nbx - 1
    j0 = js + nby - 1

    select case(homog_dir_DySma)

    case("XYZ_local_smth")

      if(xBoundary /= "periodic") stop "DYNAMIC SMAGORINSKY NOT READY FOR &
          &NON-PERIODIC BOUNDARY CONDITIONS IN X"

      if(yBoundary /= "periodic") stop "DYNAMIC SMAGORINSKY NOT READY FOR &
          &NON-PERIODIC BOUNDARY CONDITIONS IN Y"

      !---------------------------------
      !         Loop over field
      !---------------------------------

      if(nz /= sizeZ) stop " DYNAMIC SMAGORINSKY NOT READY FOR MPI IN Z"

      do k = 1, nz
        ! correct handling of solid and periodic boundaries in z

        if(zBoundary == "solid_wall") then
          kmin = max(1, k - nsmth_DySma)
          kmax = min(nz, k + nsmth_DySma)
        else if(zBoundary == "periodic") then
          kmin = k - nsmth_DySma
          kmax = k + nsmth_DySma
        else
          stop "vertical smoothing: unknown case zBoundary."
        end if

        nsmthv = kmax - kmin + 1

        do j = 1, ny
          do i = 1, nx
            ! here and below removed the conditioning of the
            ! averaging on the topography

            !UAC averaging dyn. Smag. coeff. only over atmosphere cells

            ! if(topography) then
            !   nsmthall=0
            !   var3D_DySma(i,j,k)=0.0

            !   do ksmth=kmin,kmax
            !      do jsmth=j-nsmth_DySma,j+nsmth_DySma
            !         do ismth=i-nsmth_DySma,i+nsmth_DySma
            !            if(.not.&
            !               topography_mask(i0+ismth,j0+jsmth,ksmth)) &
            !               then
            !               var3D_DySma(i,j,k)&
            !               =var3D_DySma(i,j,k)&
            !                +var3D_DySma_Extend(ismth,jsmth,ksmth)

            !               nsmthall=nsmthall+1
            !            end if
            !         end do
            !       end do
            !    end do

            !   if(nsmthall > 0) then
            !      var3D_DySma(i,j,k)=var3D_DySma(i,j,k)/nsmthall
            !   end if
            !  else
            !   var3D_DySma(i,j,k)&
            !   =&
            !   sum(&
            !   var3D_DySma_Extend( &
            !   (i-nsmth_DySma):(i+nsmth_DySma), &
            !   (j-nsmth_DySma):(j+nsmth_DySma), &
            !   kmin:kmax &
            !   )&
            !   )&
            !   /( ((2*nsmth_DySma + 1)**2) * nsmthv)
            !end if
            var3D_DySma(i, j, k) = sum(var3D_DySma_Extend((i - nsmth_DySma):(i &
                &+ nsmth_DySma), (j - nsmth_DySma):(j + nsmth_DySma), &
                &kmin:kmax)) / (((2 * nsmth_DySma + 1) ** 2) * nsmthv)
            !UAE
          end do
        end do
      end do

    case("XZ_local_smth")

      if(xBoundary /= "periodic") stop "DYNAMIC SMAGORINSKY NOT READY FOR &
          &NON-PERIODIC BOUNDARY CONDITIONS IN X"

      !---------------------------------
      !         Loop over field
      !---------------------------------

      if(nz /= sizeZ) stop " DYNAMIC SMAGORINSKY NOT READY FOR MPI IN Z"

      do k = 1, nz
        ! correct handling of solid and periodic boundaries in z

        if(zBoundary == "solid_wall") then
          kmin = max(1, k - nsmth_DySma)
          kmax = min(nz, k + nsmth_DySma)
        else if(zBoundary == "periodic") then
          kmin = k - nsmth_DySma
          kmax = k + nsmth_DySma
        else
          stop "vertical smoothing: unknown case zBoundary."
        end if

        nsmthv = kmax - kmin + 1

        do j = 1, ny
          do i = 1, nx
            !UAC averaging dyn. Smag. coeff. only over atmosphere cells

            !if(topography) then
            !  nsmthall=0
            !  var3D_DySma(i,j,k)=0.0

            !  do ksmth=kmin,kmax
            !     do ismth=i-nsmth_DySma,i+nsmth_DySma
            !        if(.not.&
            !           topography_mask(i0+ismth,j0+j,ksmth)) then
            !           var3D_DySma(i,j,k)&
            !           =var3D_DySma(i,j,k)&
            !            +var3D_DySma_Extend(ismth,j,ksmth)

            !           nsmthall=nsmthall+1
            !        end if
            !     end do
            !  end do

            !  if(nsmthall > 0) then
            !     var3D_DySma(i,j,k)=var3D_DySma(i,j,k)/nsmthall
            !  end if
            ! else
            !   var3D_DySma(i,j,k)&
            !   =&
            !   sum(&
            !   var3D_DySma_Extend( &
            !   (i-nsmth_DySma):(i+nsmth_DySma), &
            !   j, &
            !   kmin:kmax &
            !   )&
            !   )&
            !   /( (2*nsmth_DySma + 1) * nsmthv)
            !end if
            var3D_DySma(i, j, k) = sum(var3D_DySma_Extend((i - nsmth_DySma):(i &
                &+ nsmth_DySma), j, kmin:kmax)) / ((2 * nsmth_DySma + 1) &
                &* nsmthv)
            !UAE
          end do
        end do
      end do

      ! gagab
    case("YZ_local_smth")

      if(yBoundary /= "periodic") stop "DYNAMIC SMAGORINSKY NOT READY FOR &
          &NON-PERIODIC BOUNDARY CONDITIONS IN Y"

      !---------------------------------
      !         Loop over field
      !---------------------------------

      if(nz /= sizeZ) stop " DYNAMIC SMAGORINSKY NOT READY FOR MPI IN Z"

      do k = 1, nz
        ! correct handling of solid and periodic boundaries in z

        if(zBoundary == "solid_wall") then
          kmin = max(1, k - nsmth_DySma)
          kmax = min(nz, k + nsmth_DySma)
        else if(zBoundary == "periodic") then
          kmin = k - nsmth_DySma
          kmax = k + nsmth_DySma
        else
          stop "vertical smoothing: unknown case zBoundary."
        end if

        nsmthv = kmax - kmin + 1

        do j = 1, ny
          do i = 1, nx
            !UAC averaging dyn. Smag. coeff. only over atmosphere cells

            !if(topography) then
            !  nsmthall=0
            !  var3D_DySma(i,j,k)=0.0

            !  do ksmth=kmin,kmax
            !     do jsmth=j-nsmth_DySma,j+nsmth_DySma
            !        if(.not.&
            !           topography_mask(i0+i,j0+jsmth,ksmth)) then
            !           var3D_DySma(i,j,k)&
            !           =var3D_DySma(i,j,k)&
            !            +var3D_DySma_Extend(i,jsmth,ksmth)

            !           nsmthall=nsmthall+1
            !        end if
            !     end do
            !  end do

            !  if(nsmthall > 0) then
            !     var3D_DySma(i, j, k) = var3D_DySma(i, j, k) / nsmthall
            !  end if
            ! else
            !   var3D_DySma(i,j,k)&
            !   =&
            !   sum(&
            !   var3D_DySma_Extend( i,&
            !   (j-nsmth_DySma):(j+nsmth_DySma), &
            !   kmin:kmax &
            !   )&
            !   )&
            !   /( (2*nsmth_DySma + 1) * nsmthv)
            ! end if
            var3D_DySma(i, j, k) = sum(var3D_DySma_Extend(i, (j &
                &- nsmth_DySma):(j + nsmth_DySma), kmin:kmax)) / ((2 &
                &* nsmth_DySma + 1) * nsmthv)
            !UAE
          end do
        end do
      end do
      ! gagae

    case("X_whole_smth")

      if(xBoundary /= "periodic") stop "DYNAMIC SMAGORINSKY NOT READY FOR &
          &NON-PERIODIC BOUNDARY CONDITIONS IN X"

      !---------------------------------
      !         Loop over field
      !---------------------------------

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            !UAC averaging dyn. Smag. coeff. only over atmosphere cells

            !if(topography) then
            !   var3D_DySma(i,j,k)=0.0
            !   nsmthall=0

            !   do ismth =1,nx
            !      if(.not.topography_mask(i0+ismth,j0+j,k)) then
            !         var3D_DySma(i,j,k)&
            !         = var3D_DySma(i,j,k) &
            !           + var3D_DySma_Extend(ismth,j,k)

            !         nsmthall=nsmthall+1
            !      end if
            !   end do

            !   if(nsmthall > 0) then
            !      var3D_DySma(i,j,k)=var3D_DySma(i,j,k)/nsmthall
            !   end if
            !  else
            !   var3D_DySma(i,j,k)&
            !   =sum(var3D_DySma_Extend( 1:nx, j, k ))/nx
            !end if
            var3D_DySma(i, j, k) = sum(var3D_DySma_Extend(1:nx, j, k)) / nx
            !UAE
          end do
        end do
      end do

    case default
      stop "unknown case homog_dir_DySma."
    end select

    ! deallocate local fields
    deallocate(var3D_DySma_Extend, stat = allocstat); if(allocstat /= 0) stop &
        &"update.f90:dealloc failed"

    return

  end subroutine Var3DSmthDySma

  subroutine setHaloAndBoundary(var3D_HaloBC, nbx_HaloBC, nby_HaloBC, &
      &nbz_HaloBC)
    !--------------------------------------
    ! set Halo and Boundary
    !--------------------------------------

    ! in/out variables
    integer, intent(in) :: nbx_HaloBC, nby_HaloBC, nbz_HaloBC

    real, dimension(- nbx_HaloBC:nx + nbx_HaloBC, - nby_HaloBC:ny &
        &+ nby_HaloBC, - nbz_HaloBC:nz + nbz_HaloBC), intent(inout) :: &
        &var3D_HaloBC

    ! auxiliary fields for "var" with ghost cells (rho)
    real, dimension(nbx_HaloBC, - nby_HaloBC:ny + nby_HaloBC, nz) :: &
        &xRhoSliceLeft_send, xRhoSliceRight_send
    real, dimension(nbx_HaloBC, - nby_HaloBC:ny + nby_HaloBC, nz) :: &
        &xRhoSliceLeft_recv, xRhoSliceRight_recv

    real, dimension(- nbx_HaloBC:nx + nbx_HaloBC, nby_HaloBC, nz) :: &
        &yRhoSliceBack_send, yRhoSliceForw_send
    real, dimension(- nbx_HaloBC:nx + nbx_HaloBC, nby_HaloBC, nz) :: &
        &yRhoSliceBack_recv, yRhoSliceForw_recv

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: i, j, k
    integer :: i0, j0, k0

    !------------------------------
    !          x-direction
    !------------------------------

    if(idim > 1) then
      call mpi_cart_shift(comm, 0, 1, left, right, ierror)

      ! slice size
      sendcount = nbx_HaloBC * (ny + 2 * nby_HaloBC + 1) * nz
      recvcount = sendcount

      ! read slice into contiguous array
      do i = 1, nbx_HaloBC
        xRhoSliceLeft_send(i, - nby_HaloBC:ny + nby_HaloBC, 1:nz) &
            &= var3D_HaloBC(i, - nby_HaloBC:ny + nby_HaloBC, 1:nz)

        xRhoSliceRight_send(i, - nby_HaloBC:ny + nby_HaloBC, 1:nz) &
            &= var3D_HaloBC(nx - nbx_HaloBC + i, - nby_HaloBC:ny + nby_HaloBC, &
            &1:nz)
      end do

      ! left -> right
      source = left
      dest = right
      tag = 100

      i0 = 1; j0 = - nby_HaloBC; k0 = 1

      call mpi_sendrecv(xRhoSliceRight_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, xRhoSliceLeft_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xRhoSliceLeft_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, xRhoSliceRight_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_right, ierror)

      ! write auxiliary slice to var field
      do i = 1, nbx_HaloBC
        ! right halos
        var3D_HaloBC(nx + i, - nby_HaloBC:ny + nby_HaloBC, 1:nz) &
            &= xRhoSliceRight_recv(i, - nby_HaloBC:ny + nby_HaloBC, 1:nz)

        ! left halos
        var3D_HaloBC(- nbx_HaloBC + i, - nby_HaloBC:ny + nby_HaloBC, 1:nz) &
            &= xRhoSliceLeft_recv(i, - nby_HaloBC:ny + nby_HaloBC, 1:nz)
      end do
    else
      do i = 1, nbx_HaloBC
        var3D_HaloBC(nx + i, :, :) = var3D_HaloBC(i, :, :)
        var3D_HaloBC(- i + 1, :, :) = var3D_HaloBC(nx - i + 1, :, :)
      end do
    end if

    !------------------------------
    !          y-direction
    !------------------------------

    if(jdim > 1) then
      call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

      ! slice size
      sendcount = nby_HaloBC * (nx + 2 * nbx_HaloBC + 1) * nz
      recvcount = sendcount

      ! read slice into contiguous array
      do j = 1, nby_HaloBC
        yRhoSliceBack_send(- nbx_HaloBC:nx + nbx_HaloBC, j, 1:nz) &
            &= var3D_HaloBC(- nbx_HaloBC:nx + nbx_HaloBC, j, 1:nz)

        yRhoSliceForw_send(- nbx_HaloBC:nx + nbx_HaloBC, j, 1:nz) &
            &= var3D_HaloBC(- nbx_HaloBC:nx + nbx_HaloBC, ny - nby_HaloBC + j, &
            &1:nz)
      end do

      ! back -> forw
      source = back
      dest = forw
      tag = 100

      i0 = - nbx_HaloBC; j0 = 1; k0 = 1

      call mpi_sendrecv(yRhoSliceForw_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, yRhoSliceBack_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(yRhoSliceBack_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, yRhoSliceForw_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_forw, ierror)

      ! write auxiliary slice to var field
      do j = 1, nby_HaloBC
        ! right halos
        var3D_HaloBC(- nbx_HaloBC:nx + nbx_HaloBC, ny + j, 1:nz) &
            &= yRhoSliceForw_recv(- nbx_HaloBC:nx + nbx_HaloBC, j, 1:nz)

        ! left halos
        var3D_HaloBC(- nbx_HaloBC:nx + nbx_HaloBC, - nby_HaloBC + j, 1:nz) &
            &= yRhoSliceBack_recv(- nbx_HaloBC:nx + nbx_HaloBC, j, 1:nz)
      end do
    else
      do j = 1, nby_HaloBC
        var3D_HaloBC(:, ny + j, :) = var3D_HaloBC(:, j, :)
        var3D_HaloBC(:, - j + 1, :) = var3D_HaloBC(:, ny - j + 1, :)
      end do
    end if

    !------------------------------
    !          z-direction
    !------------------------------

    select case(zBoundary)

    case("periodic")

      do k = 1, nbz_HaloBC
        var3D_HaloBC(:, :, nz + k) = var3D_HaloBC(:, :, k)
        var3D_HaloBC(:, :, - k + 1) = var3D_HaloBC(:, :, nz - k + 1)
      end do

    case("solid_wall")

      do k = 1, nbz_HaloBC
        var3D_HaloBC(:, :, - k + 1) = var3D_HaloBC(:, :, k)
        var3D_HaloBC(:, :, nz + k) = var3D_HaloBC(:, :, nz - k + 1)
      end do

    case default
      stop "setBoundary: unknown case zBoundary"
    end select

    return

  end subroutine setHaloAndBoundary

  !------------------------------------------------------------------------

  subroutine smooth_hor_shapiro(fc_shap, n_shap, flux, var, dt)

    !--------------------------------------------------------------------
    !    horizontal local smoothing of density, winds, pressure,
    !    and density fluctuations
    !    order of shapiro filter given by 2*n_shap
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

    ! in/out variables
    type(var_type), intent(inout) :: var
    real, intent(inout) :: fc_shap
    integer, intent(in) :: n_shap
    type(flux_type), intent(inout) :: flux
    real, intent(in) :: dt

    ! allocatable fields
    real, dimension(:, :, :), allocatable :: field
    type(var_type) :: var_l

    integer :: allocstat
    integer :: i, j, k
    integer :: nsmth
    integer :: iVar, ivmax
    integer :: i_lapl
    integer :: nz_max

    allocate(field(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "smooth_shapiro:alloc failed"

    call allocate_var_type(var_l)
    call reset_var_type(var_l)

    if(timeScheme == "semiimplicit") then
      ! in explicit integration smoothing of density, winds,
      ! and density fluctuations
      ! pressure fluctuations are not smoothened
      ivmax = 6
    else
      ! in explicit integration smoothing of density, and winds
      ivmax = 4
    end if

    ! make sure that boundary conditions are satisfied

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    ! smoothing in x-direction

    if(sizeX > 1) then
      ! 2n-th-order x-derivative of all fields that are to be smoothed

      var_l = var

      ! in case of BLC only filtering of deviations from equilibrium state

      if((TestCase == "baroclinic_LC") .or. (TestCase == "baroclinic_ID")) then
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              if(topography) then
                var_l%rho(i, j, k) = var%rho(i, j, k) - dens_env_pp(i, j, k) &
                    &+ rhoStratTFC(i, j, k)
              else
                var_l%rho(i, j, k) = var%rho(i, j, k) - dens_env_pp(i, j, k) &
                    &+ rhoStrat(k)
              end if
            end do
          end do
        end do

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              var_l%u(i, j, k) = var%u(i, j, k) - u_env_pp(i, j, k)
            end do
          end do
        end do

        if(timeScheme == "semiimplicit") then
          if(topography) then
            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  var%rhop(i, j, k) = var%rhop(i, j, k) - dens_env_pp(i, j, k) &
                      &+ rhoStratTFC(i, j, k)
                end do
              end do
            end do
          else
            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  var_l%rhop(i, j, k) = var%rhop(i, j, k) - dens_env_pp(i, j, &
                      &k) + rhoStrat(k)
                end do
              end do
            end do
          end if
        end if

        ! horizontal boundary conditions so that everything is ready

        !for parallelized directions

        call setHalos(var_l, "var")

        ! non-parallel boundary conditions in x-direction

        select case(xBoundary)
        case("periodic")
          if(idim == 1) call setBoundary_x_periodic(var_l, flux, "var")
        case default
          stop "setBoundary: unknown case xBoundary"
        end select
      end if

      do i_lapl = 1, n_shap
        field = var_l%rho
        do i = 1, nx
          var_l%rho(i, 1:ny, 1:nz) = (field(i - 1, 1:ny, 1:nz) + field(i + 1, &
              &1:ny, 1:nz) - 2.0 * field(i, 1:ny, 1:nz)) / 4.0
        end do

        field = var_l%u
        do i = 1, nx
          var_l%u(i, 1:ny, 1:nz) = (field(i - 1, 1:ny, 1:nz) + field(i + 1, &
              &1:ny, 1:nz) - 2.0 * field(i, 1:ny, 1:nz)) / 4.0
        end do

        field = var_l%v
        do i = 1, nx
          var_l%v(i, 1:ny, 1:nz) = (field(i - 1, 1:ny, 1:nz) + field(i + 1, &
              &1:ny, 1:nz) - 2.0 * field(i, 1:ny, 1:nz)) / 4.0
        end do

        field = var_l%w
        do i = 1, nx
          var_l%w(i, 1:ny, 1:nz) = (field(i - 1, 1:ny, 1:nz) + field(i + 1, &
              &1:ny, 1:nz) - 2.0 * field(i, 1:ny, 1:nz)) / 4.0
        end do

        if(timeScheme == "semiimplicit") then
          field = var_l%rhop
          do i = 1, nx
            var_l%rhop(i, 1:ny, 1:nz) = (field(i - 1, 1:ny, 1:nz) + field(i &
                &+ 1, 1:ny, 1:nz) - 2.0 * field(i, 1:ny, 1:nz)) / 4.0
          end do
        end if

        ! horizontal boundary conditions so that everything is ready for
        ! the next iteration

        !for non-parallelized directions

        call setHalos(var_l, "var")

        ! non-parallel boundary conditions in x-direction

        select case(xBoundary)
        case("periodic")
          if(idim == 1) call setBoundary_x_periodic(var_l, flux, "var")
        case default
          stop "setBoundary: unknown case xBoundary"
        end select
      end do

      ! apply filter

      var%rho(1:nx, 1:ny, 1:nz) = var%rho(1:nx, 1:ny, 1:nz) + fc_shap * (- 1) &
          &** (n_shap + 1) * var_l%rho(1:nx, 1:ny, 1:nz)
      var%u(1:nx, 1:ny, 1:nz) = var%u(1:nx, 1:ny, 1:nz) + fc_shap * (- 1) &
          &** (n_shap + 1) * var_l%u(1:nx, 1:ny, 1:nz)
      var%v(1:nx, 1:ny, 1:nz) = var%v(1:nx, 1:ny, 1:nz) + fc_shap * (- 1) &
          &** (n_shap + 1) * var_l%v(1:nx, 1:ny, 1:nz)
      var%w(1:nx, 1:ny, 1:nz - 1) = var%w(1:nx, 1:ny, 1:nz - 1) + fc_shap * (- &
          &1) ** (n_shap + 1) * var_l%w(1:nx, 1:ny, 1:nz - 1)
      if(timeScheme == "semiimplicit") then
        var%rhop(1:nx, 1:ny, 1:nz) = var%rhop(1:nx, 1:ny, 1:nz) + fc_shap * (- &
            &1) ** (n_shap + 1) * var_l%rhop(1:nx, 1:ny, 1:nz)
      end if

      ! boundary conditions again

      call setHalos(var, "var")
      call setBoundary(var, flux, "var")
    end if

    !testb
    !goto 100
    !teste

    if(sizeY > 1) then
      ! 2n-th-order y-derivative of all fields that are to be smoothed

      var_l = var

      ! in case of BLC only filtering of deviations from equilibrium state

      if((TestCase == "baroclinic_LC") .or. (TestCase == "baroclinic_ID")) then
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              if(topography) then
                var_l%rho(i, j, k) = var%rho(i, j, k) - dens_env_pp(i, j, k) &
                    &+ rhoStratTFC(i, j, k)
              else
                var_l%rho(i, j, k) = var%rho(i, j, k) - dens_env_pp(i, j, k) &
                    &+ rhoStrat(k)
              end if
            end do
          end do
        end do

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              var_l%u(i, j, k) = var%u(i, j, k) - u_env_pp(i, j, k)
            end do
          end do
        end do

        if(timeScheme == "semiimplicit") then
          if(topography) then
            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  var%rhop(i, j, k) = var%rhop(i, j, k) - dens_env_pp(i, j, k) &
                      &+ rhoStratTFC(i, j, k)
                end do
              end do
            end do
          else
            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  var_l%rhop(i, j, k) = var%rhop(i, j, k) - dens_env_pp(i, j, &
                      &k) + rhoStrat(k)
                end do
              end do
            end do
          end if
        end if

        ! horizontal boundary conditions so that everything is ready

        !for parallelized directions

        call setHalos(var_l, "var")

        ! non-parallel boundary conditions in x-direction

        select case(yBoundary)
        case("periodic")
          if(jdim == 1) call setBoundary_y_periodic(var_l, flux, "var")
        case default
          stop "setBoundary: unknown case yBoundary"
        end select
      end if

      do i_lapl = 1, n_shap
        field = var_l%rho
        do j = 1, ny
          var_l%rho(1:nx, j, 1:nz) = (field(1:nx, j - 1, 1:nz) + field(1:nx, j &
              &+ 1, 1:nz) - 2.0 * field(1:nx, j, 1:nz)) / 4.0
        end do

        field = var_l%u
        do j = 1, ny
          var_l%u(1:nx, j, 1:nz) = (field(1:nx, j - 1, 1:nz) + field(1:nx, j &
              &+ 1, 1:nz) - 2.0 * field(1:nx, j, 1:nz)) / 4.0
        end do

        field = var_l%v
        do j = 1, ny
          var_l%v(1:nx, j, 1:nz) = (field(1:nx, j - 1, 1:nz) + field(1:nx, j &
              &+ 1, 1:nz) - 2.0 * field(1:nx, j, 1:nz)) / 4.0
        end do

        field = var_l%w
        do j = 1, ny
          var_l%w(1:nx, j, 1:nz) = (field(1:nx, j - 1, 1:nz) + field(1:nx, j &
              &+ 1, 1:nz) - 2.0 * field(1:nx, j, 1:nz)) / 4.0
        end do

        if(timeScheme == "semiimplicit") then
          field = var_l%rhop
          do j = 1, ny
            var_l%rhop(1:nx, j, 1:nz) = (field(1:nx, j - 1, 1:nz) &
                &+ field(1:nx, j + 1, 1:nz) - 2.0 * field(1:nx, j, 1:nz)) / 4.0
          end do
        end if

        ! horizontal boundary conditions so that everything is ready for
        ! the next iteration

        !for non-parallelized directions

        call setHalos(var_l, "var")

        ! non-parallel boundary conditions in y-direction

        select case(yBoundary)
        case("periodic")
          if(jdim == 1) call setBoundary_y_periodic(var_l, flux, "var")
        case default
          stop "setBoundary: unknown case yBoundary"
        end select
      end do

      ! apply filter

      var%rho(1:nx, 1:ny, 1:nz) = var%rho(1:nx, 1:ny, 1:nz) + fc_shap * (- 1) &
          &** (n_shap + 1) * var_l%rho(1:nx, 1:ny, 1:nz)
      var%u(1:nx, 1:ny, 1:nz) = var%u(1:nx, 1:ny, 1:nz) + fc_shap * (- 1) &
          &** (n_shap + 1) * var_l%u(1:nx, 1:ny, 1:nz)
      var%v(1:nx, 1:ny, 1:nz) = var%v(1:nx, 1:ny, 1:nz) + fc_shap * (- 1) &
          &** (n_shap + 1) * var_l%v(1:nx, 1:ny, 1:nz)
      var%w(1:nx, 1:ny, 1:nz - 1) = var%w(1:nx, 1:ny, 1:nz - 1) + fc_shap * (- &
          &1) ** (n_shap + 1) * var_l%w(1:nx, 1:ny, 1:nz - 1)
      if(timeScheme == "semiimplicit") then
        var%rhop(1:nx, 1:ny, 1:nz) = var%rhop(1:nx, 1:ny, 1:nz) + fc_shap * (- &
            &1) ** (n_shap + 1) * var_l%rhop(1:nx, 1:ny, 1:nz)
      end if

      ! boundary conditions again

      call setHalos(var, "var")
      call setBoundary(var, flux, "var")
    end if

    !testb
    100 continue
    !teste

    ! deallocate local fields

    deallocate(field, stat = allocstat); if(allocstat /= 0) stop &
        &"smooth_shapiro:dealloc failed"

    call deallocate_var_type(var_l)

    return

  end subroutine smooth_hor_shapiro

  !------------------------------------------------------------------------

  subroutine smooth_shapiro(fc_shap, n_shap, flux, var)

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

    ! in/out variables
    type(var_type), intent(inout) :: var
    real, intent(in) :: fc_shap
    integer, intent(in) :: n_shap
    type(flux_type), intent(inout) :: flux

    ! allocatable fields
    real, dimension(:, :, :), allocatable :: field
    type(var_type) :: var_l

    integer :: allocstat
    integer :: i, j, k
    integer :: nsmth
    integer :: iVar, ivmax
    integer :: i_lapl
    integer :: nz_max

    allocate(field(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "smooth_shapiro:alloc failed"

    call allocate_var_type(var_l)
    call reset_var_type(var_l)

    if(timeScheme == "semiimplicit") then
      ! in semiimplicit integration smoothing of density, winds,
      ! and density fluctuations
      ! pressure fluctuations are not smoothened
      ivmax = 6
    else
      ! in explicit integration smoothing of density, and winds
      ivmax = 4
    end if

    ! make sure that boundary conditions are satisfied

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    ! smoothing in x-direction

    if(sizeX > 1) then
      ! 2n-th-order x-derivative of all fields that are to be smoothed

      var_l = var

      ! in case of BLC only filtering of deviations from equilibrium state

      if((TestCase == "baroclinic_LC") .or. (TestCase == "baroclinic_ID")) then
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              if(topography) then
                var_l%rho(i, j, k) = var%rho(i, j, k) - dens_env_pp(i, j, k) &
                    &+ rhoStratTFC(i, j, k)
              else
                var_l%rho(i, j, k) = var%rho(i, j, k) - dens_env_pp(i, j, k) &
                    &+ rhoStrat(k)
              end if
            end do
          end do
        end do

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              var_l%u(i, j, k) = var%u(i, j, k) - u_env_pp(i, j, k)
            end do
          end do
        end do

        if(timeScheme == "semiimplicit") then
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                var_l%rhop(i, j, k) = var%rhop(i, j, k) - dens_env_pp(i, j, k) &
                    &+ rhoStrat(k)
              end do
            end do
          end do
        end if

        ! horizontal boundary conditions so that everything is ready

        !for parallelized directions

        call setHalos(var_l, "var")

        ! non-parallel boundary conditions in x-direction

        select case(xBoundary)
        case("periodic")
          if(idim == 1) call setBoundary_x_periodic(var_l, flux, "var")
        case default
          stop "setBoundary: unknown case xBoundary"
        end select
      end if

      do i_lapl = 1, n_shap
        field = var_l%rho
        do i = 1, nx
          var_l%rho(i, 1:ny, 1:nz) = (field(i - 1, 1:ny, 1:nz) + field(i + 1, &
              &1:ny, 1:nz) - 2.0 * field(i, 1:ny, 1:nz)) / 4.0
        end do

        field = var_l%u
        do i = 1, nx
          var_l%u(i, 1:ny, 1:nz) = (field(i - 1, 1:ny, 1:nz) + field(i + 1, &
              &1:ny, 1:nz) - 2.0 * field(i, 1:ny, 1:nz)) / 4.0
        end do

        field = var_l%v
        do i = 1, nx
          var_l%v(i, 1:ny, 1:nz) = (field(i - 1, 1:ny, 1:nz) + field(i + 1, &
              &1:ny, 1:nz) - 2.0 * field(i, 1:ny, 1:nz)) / 4.0
        end do

        field = var_l%w
        do i = 1, nx
          var_l%w(i, 1:ny, 1:nz) = (field(i - 1, 1:ny, 1:nz) + field(i + 1, &
              &1:ny, 1:nz) - 2.0 * field(i, 1:ny, 1:nz)) / 4.0
        end do

        if(timeScheme == "semiimplicit") then
          field = var_l%rhop
          do i = 1, nx
            var_l%rhop(i, 1:ny, 1:nz) = (field(i - 1, 1:ny, 1:nz) + field(i &
                &+ 1, 1:ny, 1:nz) - 2.0 * field(i, 1:ny, 1:nz)) / 4.0
          end do
        end if

        ! horizontal boundary conditions so that everything is ready for
        ! the next iteration

        !for non-parallelized directions

        call setHalos(var_l, "var")

        ! non-parallel boundary conditions in x-direction

        select case(xBoundary)
        case("periodic")
          if(idim == 1) call setBoundary_x_periodic(var_l, flux, "var")
        case default
          stop "setBoundary: unknown case xBoundary"
        end select
      end do

      ! apply filter

      var%rho(1:nx, 1:ny, 1:nz) = var%rho(1:nx, 1:ny, 1:nz) + fc_shap * (- 1) &
          &** (n_shap + 1) * var_l%rho(1:nx, 1:ny, 1:nz)
      var%u(1:nx, 1:ny, 1:nz) = var%u(1:nx, 1:ny, 1:nz) + fc_shap * (- 1) &
          &** (n_shap + 1) * var_l%u(1:nx, 1:ny, 1:nz)
      var%v(1:nx, 1:ny, 1:nz) = var%v(1:nx, 1:ny, 1:nz) + fc_shap * (- 1) &
          &** (n_shap + 1) * var_l%v(1:nx, 1:ny, 1:nz)
      var%w(1:nx, 1:ny, 1:nz - 1) = var%w(1:nx, 1:ny, 1:nz - 1) + fc_shap * (- &
          &1) ** (n_shap + 1) * var_l%w(1:nx, 1:ny, 1:nz - 1)
      if(timeScheme == "semiimplicit") then
        var%rhop(1:nx, 1:ny, 1:nz) = var%rhop(1:nx, 1:ny, 1:nz) + fc_shap * (- &
            &1) ** (n_shap + 1) * var_l%rhop(1:nx, 1:ny, 1:nz)
      end if

      ! boundary conditions again

      call setHalos(var, "var")
      call setBoundary(var, flux, "var")
    end if

    !testb
    !goto 100
    !teste

    !smoothing in z-direction

    if(sizeY > 1) then
      ! 2n-th-order y-derivative of all fields that are to be smoothed

      var_l = var

      ! in case of BLC only filtering of deviations from equilibrium state

      if((TestCase == "baroclinic_LC") .or. (TestCase == "baroclinic_ID")) then
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              if(topography) then
                var_l%rho(i, j, k) = var%rho(i, j, k) - dens_env_pp(i, j, k) &
                    &+ rhoStratTFC(i, j, k)
              else
                var_l%rho(i, j, k) = var%rho(i, j, k) - dens_env_pp(i, j, k) &
                    &+ rhoStrat(k)
              end if
            end do
          end do
        end do

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              var_l%u(i, j, k) = var%u(i, j, k) - u_env_pp(i, j, k)
            end do
          end do
        end do

        if(timeScheme == "semiimplicit") then
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                var_l%rhop(i, j, k) = var%rhop(i, j, k) - dens_env_pp(i, j, k) &
                    &+ rhoStrat(k)
              end do
            end do
          end do
        end if

        ! horizontal boundary conditions so that everything is ready

        !for parallelized directions

        call setHalos(var_l, "var")

        ! non-parallel boundary conditions in x-direction

        select case(yBoundary)
        case("periodic")
          if(jdim == 1) call setBoundary_y_periodic(var_l, flux, "var")
        case default
          stop "setBoundary: unknown case yBoundary"
        end select
      end if

      do i_lapl = 1, n_shap
        field = var_l%rho
        do j = 1, ny
          var_l%rho(1:nx, j, 1:nz) = (field(1:nx, j - 1, 1:nz) + field(1:nx, j &
              &+ 1, 1:nz) - 2.0 * field(1:nx, j, 1:nz)) / 4.0
        end do

        field = var_l%u
        do j = 1, ny
          var_l%u(1:nx, j, 1:nz) = (field(1:nx, j - 1, 1:nz) + field(1:nx, j &
              &+ 1, 1:nz) - 2.0 * field(1:nx, j, 1:nz)) / 4.0
        end do

        field = var_l%v
        do j = 1, ny
          var_l%v(1:nx, j, 1:nz) = (field(1:nx, j - 1, 1:nz) + field(1:nx, j &
              &+ 1, 1:nz) - 2.0 * field(1:nx, j, 1:nz)) / 4.0
        end do

        field = var_l%w
        do j = 1, ny
          var_l%w(1:nx, j, 1:nz) = (field(1:nx, j - 1, 1:nz) + field(1:nx, j &
              &+ 1, 1:nz) - 2.0 * field(1:nx, j, 1:nz)) / 4.0
        end do

        if(timeScheme == "semiimplicit") then
          field = var_l%rhop
          do j = 1, ny
            var_l%rhop(1:nx, j, 1:nz) = (field(1:nx, j - 1, 1:nz) &
                &+ field(1:nx, j + 1, 1:nz) - 2.0 * field(1:nx, j, 1:nz)) / 4.0
          end do
        end if

        ! horizontal boundary conditions so that everything is ready for
        ! the next iteration

        !for non-parallelized directions

        call setHalos(var_l, "var")

        ! non-parallel boundary conditions in y-direction

        select case(yBoundary)
        case("periodic")
          if(jdim == 1) call setBoundary_y_periodic(var_l, flux, "var")
        case default
          stop "setBoundary: unknown case yBoundary"
        end select
      end do

      ! apply filter

      var%rho(1:nx, 1:ny, 1:nz) = var%rho(1:nx, 1:ny, 1:nz) + fc_shap * (- 1) &
          &** (n_shap + 1) * var_l%rho(1:nx, 1:ny, 1:nz)
      var%u(1:nx, 1:ny, 1:nz) = var%u(1:nx, 1:ny, 1:nz) + fc_shap * (- 1) &
          &** (n_shap + 1) * var_l%u(1:nx, 1:ny, 1:nz)
      var%v(1:nx, 1:ny, 1:nz) = var%v(1:nx, 1:ny, 1:nz) + fc_shap * (- 1) &
          &** (n_shap + 1) * var_l%v(1:nx, 1:ny, 1:nz)
      var%w(1:nx, 1:ny, 1:nz - 1) = var%w(1:nx, 1:ny, 1:nz - 1) + fc_shap * (- &
          &1) ** (n_shap + 1) * var_l%w(1:nx, 1:ny, 1:nz - 1)
      if(timeScheme == "semiimplicit") then
        var%rhop(1:nx, 1:ny, 1:nz) = var%rhop(1:nx, 1:ny, 1:nz) + fc_shap * (- &
            &1) ** (n_shap + 1) * var_l%rhop(1:nx, 1:ny, 1:nz)
      end if

      ! boundary conditions again

      call setHalos(var, "var")
      call setBoundary(var, flux, "var")
    end if

    ! filtering in z-direction

    ! 2nd-order z-derivative of all fields that are to be smoothed

    var_l = var

    ! in case of BLC only filtering of deviations from equilibrium state

    if((TestCase == "baroclinic_LC") .or. (TestCase == "baroclinic_ID")) then
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            if(topography) then
              var_l%rho(i, j, k) = var%rho(i, j, k) - dens_env_pp(i, j, k) &
                  &+ rhoStratTFC(i, j, k)
            else
              var_l%rho(i, j, k) = var%rho(i, j, k) - dens_env_pp(i, j, k) &
                  &+ rhoStrat(k)
            end if
          end do
        end do
      end do

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            var_l%u(i, j, k) = var%u(i, j, k) - u_env_pp(i, j, k)
          end do
        end do
      end do

      if(timeScheme == "semiimplicit") then
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              var_l%rhop(i, j, k) = var%rhop(i, j, k) - dens_env_pp(i, j, k) &
                  &+ rhoStrat(k)
            end do
          end do
        end do
      end if

    end if

    field = var_l%rho
    do k = 1, nz
      var_l%rho(1:nx, 1:ny, k) = (field(1:nx, 1:ny, k - 1) + field(1:nx, 1:ny, &
          &k + 1) - 2.0 * field(1:nx, 1:ny, k)) / 4.0
    end do

    field = var_l%u
    do k = 1, nz
      var_l%u(1:nx, 1:ny, k) = (field(1:nx, 1:ny, k - 1) + field(1:nx, 1:ny, k &
          &+ 1) - 2.0 * field(1:nx, 1:ny, k)) / 4.0
    end do

    field = var_l%v
    do k = 1, nz
      var_l%v(1:nx, 1:ny, k) = (field(1:nx, 1:ny, k - 1) + field(1:nx, 1:ny, k &
          &+ 1) - 2.0 * field(1:nx, 1:ny, k)) / 4.0
    end do

    field = var_l%w
    do k = 1, nz
      var_l%w(1:nx, 1:ny, k) = (field(1:nx, 1:ny, k - 1) + field(1:nx, 1:ny, k &
          &+ 1) - 2.0 * field(1:nx, 1:ny, k)) / 4.0
    end do

    if(timeScheme == "semiimplicit") then
      field = var_l%rhop
      do k = 1, nz
        var_l%rhop(1:nx, 1:ny, k) = (field(1:nx, 1:ny, k - 1) + field(1:nx, &
            &1:ny, k + 1) - 2.0 * field(1:nx, 1:ny, k)) / 4.0
      end do
    end if

    ! apply filter

    var%rho(1:nx, 1:ny, 1:nz) = var%rho(1:nx, 1:ny, 1:nz) + 1.e-2 * fc_shap &
        &* var_l%rho(1:nx, 1:ny, 1:nz)
    var%u(1:nx, 1:ny, 1:nz) = var%u(1:nx, 1:ny, 1:nz) + 1.e-2 * fc_shap &
        &* var_l%u(1:nx, 1:ny, 1:nz)
    var%v(1:nx, 1:ny, 1:nz) = var%v(1:nx, 1:ny, 1:nz) + 1.e-2 * fc_shap &
        &* var_l%v(1:nx, 1:ny, 1:nz)
    var%w(1:nx, 1:ny, 1:nz - 1) = var%w(1:nx, 1:ny, 1:nz - 1) + 1.e-2 &
        &* fc_shap * var_l%w(1:nx, 1:ny, 1:nz - 1)
    if(timeScheme == "semiimplicit") then
      var%rhop(1:nx, 1:ny, 1:nz) = var%rhop(1:nx, 1:ny, 1:nz) + 1.e-2 &
          &* fc_shap * var_l%rhop(1:nx, 1:ny, 1:nz)
    end if

    ! boundary conditions again

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    do k = 1, nz !theta
      do j = 1, ny
        do i = 1, nx
          var%rho(i, j, k) = (var%rho(i, j, k) - rhoStrat(k)) / PStrat(k) !- dens_env_pp(i, j, k) + rhoStrat(k)

        end do
      end do
    end do

    !testb
    100 continue
    !teste

    ! deallocate local fields

    deallocate(field, stat = allocstat); if(allocstat /= 0) stop &
        &"smooth_shapiro:dealloc failed"

    call deallocate_var_type(var_l)

    return

  end subroutine smooth_shapiro
  !UAE

  !---------------------------------------------------------------------

  subroutine BGstate_update(var, flux, dt, m, q_P, q_rho, int_mod, &
      &heating_switch)

    ! (1) update of the reference-atmosphere profile
    ! (2) synchronization of the vertical wind with the assumption that the
    ! w0 according to O'Neill & Klein (2014) is the horizontal-mean
    ! vertical wind

    ! in/out variables
    type(var_type), intent(inout) :: var
    type(flux_type), intent(in) :: flux

    real, intent(in) :: dt
    integer, intent(in) :: m

    real, dimension(- nbz:nz + nbz), intent(inout) :: q_P, q_rho

    character(len = *), intent(in) :: int_mod

    integer, intent(in) :: heating_switch

    character(len = 40) :: w0_mod

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: heat
    real, dimension(- nbz:nz + nbz) :: S_bar, PStratold
    integer :: i, j, k
    real, dimension(1:nz) :: sum_local, sum_global, sum_local2, sum_global2
    real, dimension(- nbz:nz + nbz) :: press0
    real, dimension(- nbz:nz + nbz) :: divPw !, divrhow !UA

    real, dimension(- nbz:nz + nbz) :: rhow_bar, rho_bar

    real, dimension(- nbz:nz + nbz) :: w_0 !UA

    real :: dptopdt
    real :: expo

    real :: sum_d, sum_n

    ! No heating in TFC (FJApr2023)
    if(topography) return

    ! w0_mod = 'Almgrenetal08'
    w0_mod = 'ONK14'

    w_0 = 0.
    S_bar = 0.
    heat = 0.

    divPw = 0.

    ! -S eq(9)  ONeill+Klein2014
    call calculate_heating(var, flux, heat)

    ! calculate horizontal mean of heat(:,:,:)
    do k = 1, nz
      sum_local(k) = sum(heat(1:nx, 1:ny, k))
    end do
    !global sum and average
    call mpi_allreduce(sum_local(1), sum_global(1), nz - 1 + 1, &
        &mpi_double_precision, mpi_sum, comm, ierror)
    sum_global = sum_global / (sizeX * sizeY)

    S_bar(1:nz) = sum_global(1:nz)

    if(w0_mod == 'Almgrenetal08') then
      ! horizontal mean of the vertical density flux

      rhow_bar = 0.
      rho_bar = 0.

      sum_local = 0.
      sum_global = 0.
      sum_local2 = 0.
      sum_global2 = 0.

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx

            sum_local2(k) = sum_local2(k) + var%rho(i, j, k) + rhoStrat(k)

            if(k == 1) then
              sum_local(k) = sum_local(k) + 0.5 * flux%rho(i, j, k, 3)
            else if(k == nz) then
              sum_local(k) = sum_local(k) + 0.5 * flux%rho(i, j, k - 1, 3)
            else
              sum_local(k) = sum_local(k) + 0.5 * (flux%rho(i, j, k - 1, 3) &
                  &+ flux%rho(i, j, k, 3))
            end if
          end do
        end do
      end do
      call mpi_allreduce(sum_local(1), sum_global(1), nz - 1 + 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      sum_global = sum_global / (sizeX * sizeY)

      rhow_bar(1:nz) = sum_global(1:nz)

      call mpi_allreduce(sum_local2(1), sum_global2(1), nz - 1 + 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      sum_global2 = sum_global2 / (sizeX * sizeY)

      rho_bar(1:nz) = sum_global2(1:nz)

    end if

    !non_dim. pressure; eq(12) ONeill+Klein2014

    do k = 1, nz
      press0(k) = PStrat(k) ** gamma
    end do

    !  time derivative of reference pressure at the model top

    sum_d = 0.0
    sum_n = 0.0

    if(w0_mod == 'Almgrenetal08') then
      do k = 1, nz

        expo = exp(- g_ndim * (rho_bar(k) + rhoStrat_s(k)) / (gamma &
            &* press0(k)) * z(k))

        sum_n = sum_n + expo * (- S_bar(k) / PStrat(k) + g_ndim * rhow_bar(k) &
            &/ (gamma * press0(k)))

        sum_d = sum_d + expo / (gamma * press0(k))
      end do
    else if(w0_mod == 'ONK14') then
      do k = 1, nz
        sum_n = sum_n - S_bar(k) / PStrat(k)
        sum_d = sum_d + 1. / (gamma * press0(k))
      end do
    else
      stop 'ERROR: wrong w0_mod'
    end if

    dptopdt = sum_n / sum_d

    ! horizontal-mean vertical wind

    w_0 = 0.

    if(w0_mod == 'Almgrenetal08') then
      do k = 1, nz - 1
        expo = exp(- g_ndim * (rho_bar(k) + rhoStrat_s(k)) / (gamma &
            &* press0(k)) * z(k))

        w_0(k) = w_0(k - 1) + dz * expo * (- S_bar(k) / Pstrat(k) + g_ndim &
            &* rhow_bar(k) / (gamma * press0(k)) - dptopdt / (gamma &
            &* press0(k)))
      end do

      do k = 1, nz - 1
        expo = exp(g_ndim * (rho_bar(k) + rhoStrat_s(k)) / (gamma * press0(k)) &
            &* 0.5 * (z(k) + z(k + 1)))

        w_0(k) = expo * w_0(k)
      end do
    else if(w0_mod == 'ONK14') then
      w_0(1) = dz * (- S_bar(1) / Pstrat(1) - (1. / (gamma * press0(1))) &
          &* dptopdt)

      do k = 2, nz - 1
        w_0(k) = w_0(k - 1) + dz * (- S_bar(k) / Pstrat(k) - (1. / (gamma &
            &* press0(k))) * dptopdt)
      end do
    else
      stop 'ERROR: wrong w0_mod'
    end if

    ! update PStrat

    do k = 1, nz
      divPw(k) = (PstratTilde(k) * w_0(k) - PstratTilde(k - 1) * w_0(k - 1)) &
          &/ dz
    end do

    do k = 1, nz
      if(int_mod == "expl") then
        !init q
        if(m == 1) then
          q_P(k) = 0.
          q_rho(k) = 0.
        end if

        ! update: q(m-1) -> q(m)

        q_P(k) = alphaRK(m) * q_P(k) - dt * divPw(k) - dt * S_bar(k)

        ! update PStrat

        Pstrat(k) = Pstrat(k) + betaRK(m) * q_P(k) !PIold
      else if(int_mod == "impl") then
        if(heating_switch == 0) then
          PStrat(k) = PStrat(k) - dt * divPw(k) - dt * S_bar(k) !PIold
        else
          PStrat(k) = PStrat00(k) - dt * divPw(k) - dt * S_bar(k) !PIold
        end if
      else
        print *, "update.f90/BGstate_update: wrong int_mod"
        stop
      end if

      if(PStrat(k) <= 0.) then
        print *, 'ERROR in BGstate_update: PStrat(', k, ') =', PStrat(k), '<=0.'
        stop
      end if

      piStrat(k) = PStrat(k) ** (kappa / (1.0 - kappa))
    end do

    pStrat(0) = pStrat(1)
    pStrat(- 1) = pStrat(0)
    pStrat(nz + 1) = pStrat(nz)
    pStrat(nz + 2) = pStrat(nz + 1)

    do k = 1, nz
      PstratTilde(k) = 0.5 * (PStrat(k) + PStrat(k + 1))
    end do

    ! the following could most probably be deleted
    ! update of non-dimensional squared Brunt-Vaisala frequency
    ! (this could perhaps be done a bit nicer)

    N2 = 0.

    do k = 1, nz
      if(k == 1) then
        bvsStrat(k) = g_ndim / thetaStrat(k) * (thetaStrat(k + 1) &
            &- thetaStrat(k)) / dz
      else if(k == nz) then
        bvsStrat(k) = g_ndim / thetaStrat(k) * (thetaStrat(k) - thetaStrat(k &
            &- 1)) / dz
      else
        bvsStrat(k) = g_ndim / thetaStrat(k) * (thetaStrat(k + 1) &
            &- thetaStrat(k - 1)) / (2.0 * dz)
      end if

      N2 = max(N2, bvsStrat(k))
    end do

    bvsStrat(- 1) = bvsStrat(1)
    bvsStrat(0) = bvsStrat(1)

    bvsStrat(nz + 1) = bvsStrat(nz)
    bvsStrat(nz + 2) = bvsStrat(nz)

    if(N2 < 0.) then
      stop 'ERROR: N2 < 0'
    else
      NN = sqrt(N2)
    end if

  end subroutine BGstate_update

  !---------------------------------------------------------------------

  subroutine piUpdate(var, pinew, dt, int_mod, flux)
    !-----------------------------------------
    ! in the compressible model pi' is updated
    ! in the explicit and implicit Euler step
    !-----------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: pinew
    real, intent(in) :: dt
    character(len = *), intent(in) :: int_mod
    type(flux_type), intent(in) :: flux

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: heat

    ! local variables
    integer :: i, j, k
    real :: fL, fR, gB, gF, hD, hU
    real :: fluxDiff
    real :: dPdPi

    select case(int_mod)

    case("expl")
      do k = 1, nz + 1
        do j = 1, ny + 1
          do i = 1, nx + 1
            ! using JPu as the carrier flux
            fR = var%u(i, j, k) ! mass flux accros right cell edge
            fL = var%u(i - 1, j, k) ! left
            gF = var%v(i, j, k) ! forward
            gB = var%v(i, j - 1, k) ! backward
            hU = var%w(i, j, k) ! upward
            hD = var%w(i, j, k - 1) ! downward

            ! convective part
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

            fluxDiff = fluxDiff / jac(i, j, k)

            ! (\partial P / \partial \pi) from latest update of P
            dPdPi = 1 / (gamma - 1) * (Rsp / pref) ** (1 - gamma) * var%P(i, &
                &j, k) ** (2 - gamma)

            ! update density
            pinew(i, j, k) = var%pi(i, j, k) - dt * fluxDiff / dPdPi
          end do
        end do
      end do

    case("expl_heating")

      call calculate_heating(var, flux, heat)
      do k = - nbz, nz + nbz
        do j = - nby, ny + nby
          do i = - nbx, nx + nbx
            dPdPi = 1 / (gamma - 1) * (Rsp / pref) ** (1 - gamma) * var%P(i, &
                &j, k) ** (2 - gamma)

            pinew(i, j, k) = var%pi(i, j, k) - dt * heat(i, j, k) / dPdPi
          end do
        end do
      end do

    case default
      stop "piUpdate: Unknown case int_mod."

    end select

  end subroutine piUpdate
  !---------------------------------------------------------------------

  ! SK
  subroutine bvsUpdate(bvs, var)
    !--------------------------------------
    ! updates N^2 in the compressible model
    ! because it is timedependent
    !--------------------------------------
    ! in/out variables
    type(var_type), intent(in) :: var
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, (- 1):(nz + 2)), &
        &intent(inout) :: bvs

    ! local variables
    integer :: i, j, k

    if(testCase == "SkamarockKlemp94") then
      return
    end if

    ! N^2 = gP / ( \rho \Bar{\theta}^2 J) * d J\Bar{\theta} / (dz)
    if(topography) then
      ! bvs = 0.0
      do i = - nbx, nx + nbx
        do j = - nby, ny + nby
          ! Lower boundary.
          bvs(i, j, - 1) = g_ndim * var%P(i, j, 0) / (var%rho(i, j, 0) &
              &+ rhoStratTFC(i, j, 0)) / (thetaStratTFC(i, j, 0) ** 2) &
              &/ jac(i, j, 0) * (thetaStratTFC(i, j, 1) - thetaStratTFC(i, j, &
              &0)) / dz
          bvs(i, j, 0) = bvs(i, j, - 1)
          ! Between boundaries.
          do k = 1, nz
            bvs(i, j, k) = g_ndim * var%P(i, j, k) / (var%rho(i, j, k) &
                &+ rhoStratTFC(i, j, k)) / (thetaStratTFC(i, j, k) ** 2) &
                &/ jac(i, j, k) * 0.5 * (thetaStratTFC(i, j, k + 1) &
                &- thetaStratTFC(i, j, k - 1)) / dz
          end do
          ! Upper boundary.
          bvs(i, j, nz + 1) = g_ndim * var%P(i, j, nz + 1) / (var%rho(i, j, nz &
              &+ 1) + rhoStratTFC(i, j, nz + 1)) / (thetaStratTFC(i, j, nz &
              &+ 1) ** 2) / jac(i, j, nz + 1) * (thetaStratTFC(i, j, nz + 1) &
              &- thetaStratTFC(i, j, nz)) / dz
          bvs(i, j, nz + 2) = bvs(i, j, nz + 1)
        end do
      end do
    end if
  end subroutine bvsUpdate

end module update_module
