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
  public :: timestep
  public :: init_update
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

  subroutine momentumPredictor(var, flux, dt, q, m, mmp_mod, int_mod, &
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

    if(corset == 'constant') then
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

    case("periodic")
      i0 = 0
      i1 = nx
    case default
      stop "momentumPredictor: unknown case xBoundary."
    end select

    if(mmp_mod == "lhs") then
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

  
            

            ! Explicit integration of Coriolis force in TFC.
            if(topography .and. mmp_mod == "lhs") then
              uOldTFC(i, j, k) = var%u(i, j, k)
              vC = 0.5 * (var%v(i, j, k) + var%v(i, j - 1, k))
              vR = 0.5 * (var%v(i + 1, j, k) + var%v(i + 1, j - 1, k))
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
            if(mmp_mod == "lhs") then
              F = - fluxDiff + volForce !200413
            else
              stop 'ERROR: wrong mmp_mod'
            end if

            ! interpolated density
            select case(model)

            case("pseudo_incompressible")

              rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i + 1, j, k))
              rhoM = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k))

              if(topography) then
                rhoStratEdgeR = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i &
                    &+ 1, j, k))
                rhoM_1 = rhoM_1 + rhoStratEdgeR
                rhoM = rhoM + rhoStratEdgeR
              end if
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
              end if

              
              !--- pressure gradient term -> piGrad
              piR = var%pi(i + 1, j, k)
                piL = var%pi(i, j, k)

              if(topography) then
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
                      &/ dx + met13EdgeR * (- piUUEdgeR &
                      &+ 4.0 * piUEdgeR - 3.0 * piEdgeR) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeR = 0.5 * (var%pi(i, j, k - 2) + var%pi(i + 1, j, k &
                      &- 2))
                  piDEdgeR = 0.5 * (var%pi(i, j, k - 1) + var%pi(i + 1, j, k &
                      &- 1))
                  piEdgeR = 0.5 * (var%pi(i, j, k) + var%pi(i + 1, j, k))
                  piGrad = kappaInv * MaInv2 * pEdgeR / rhou * ((piR - piL) &
                      &/ dx + met13EdgeR * (piDDEdgeR - 4.0 &
                      &* piDEdgeR + 3.0 * piEdgeR) * 0.5 / dz)
                else
                  piUEdgeR = 0.5 * (var%pi(i, j, k + 1) + var%pi(i + 1, j, k &
                      &+ 1))
                  piDEdgeR = 0.5 * (var%pi(i, j, k - 1) + var%pi(i + 1, j, k &
                      &- 1))
                  piGrad = kappaInv * MaInv2 * pEdgeR / rhou * ((piR - piL) &
                      &/ dx + met13EdgeR * (piUEdgeR &
                      &- piDEdgeR) * 0.5 / dz)
                end if
              end if

              volfcx = 0.0

              ! ustar
              uhorx = var%u(i, j, k)

              if(topography) then
                ! Coriolis force is integrated on LHS.
                uAst = uhorx + dt * (- piGrad + volfcx / rhou)
              end if

              if(spongeLayer .and. sponge_uv) then
                if(topography) then
                  uAst = uAst - dt * 0.5 * (kr_sp_tfc(i, j, k) + kr_sp_tfc(i &
                      &+ 1, j, k)) * uhorx
                end if
              end if

              usave(i, j, k) = uAst

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
              end if

              
              !--- pressure gradient terms -> piGradx, piGrady
              piR = var%pi(i + 1, j, k)
                piL = var%pi(i, j, k)

              if(topography) then
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
                      &/ dx + met13EdgeR * (- piUUEdgeR &
                      &+ 4.0 * piUEdgeR - 3.0 * piEdgeR) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeR = 0.5 * (var%pi(i, j, k - 2) + var%pi(i + 1, j, k &
                      &- 2))
                  piDEdgeR = 0.5 * (var%pi(i, j, k - 1) + var%pi(i + 1, j, k &
                      &- 1))
                  piEdgeR = 0.5 * (var%pi(i, j, k) + var%pi(i + 1, j, k))
                  piGradX = kappaInv * MaInv2 * pEdgeR / rhou * ((piR - piL) &
                      &/ dx + met13EdgeR * (piDDEdgeR &
                      &- 4.0 * piDEdgeR + 3.0 * piEdgeR) * 0.5 / dz)
                else
                  piUEdgeR = 0.5 * (var%pi(i, j, k + 1) + var%pi(i + 1, j, k &
                      &+ 1))
                  piDEdgeR = 0.5 * (var%pi(i, j, k - 1) + var%pi(i + 1, j, k &
                      &- 1))
                  piGradX = kappaInv * MaInv2 * pEdgeR / rhou * ((piR - piL) &
                      &/ dx + met13EdgeR * (piUEdgeR &
                      &- piDEdgeR) * 0.5 / dz)
                end if
              end if

              volfcx = 0.0
                volfcy = 0.0

              ! ustar
              uhorx = var%u(i, j, k)

              vhory = 0.25 * (var%v(i, j - 1, k) + var%v(i, j, k) + var%v(i &
                  &+ 1, j - 1, k) + var%v(i + 1, j, k))

              facu = 1.0

              if(spongeLayer .and. sponge_uv) then
                if(topography) then
                  facu = facu + dt * 0.5 * (kr_sp_tfc(i, j, k) + kr_sp_tfc(i &
                      &+ 1, j, k))
                end if
              end if

              facv = facu

              if(topography) then
                ! Coriolis force is integrated on LHS.
                uAst = 1.0 / facu * (uhorx + dt * (- piGradX + volfcx / rhou))
              end if

              usave(i, j, k) = uAst

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

    case("periodic")
      j0 = 0
      j1 = ny
    case default
      stop "momentumPredictor: unknown case yBoundary."
    end select

    if(mmp_mod == "lhs") then
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

            ! Explicit integration of Coriolis force in TFC.
            if(topography .and. mmp_mod == "lhs") then
              vOldTFC(i, j, k) = var%v(i, j, k)
              uC = 0.5 * (uOldTFC(i, j, k) + uOldTFC(i - 1, j, k))
              uF = 0.5 * (uOldTFC(i, j + 1, k) + uOldTFC(i - 1, j + 1, k))
              
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
            if(mmp_mod == "lhs") then
              F = - fluxDiff + volForce !UA 200413
            else
              stop 'ERROR: wrong mmp_mod'
            end if

            ! interpolated density
            select case(model)

            case("pseudo_incompressible")

              rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j + 1, k))
              rhoM = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k))

              if(topography) then
                rhoStratEdgeF = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, j &
                    &+ 1, k))
                rhoM_1 = rhoM_1 + rhoStratEdgeF
                rhoM = rhoM + rhoStratEdgeF
              end if

            
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
              end if

              piF = var%pi(i, j + 1, k)
                piB = var%pi(i, j, k)

              if(topography) then
                
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
                      &/ dy + met23EdgeF * (- piUUEdgeF &
                      &+ 4.0 * piUEdgeF - 3.0 * piEdgeF) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeF = 0.5 * (var%pi(i, j, k - 2) + var%pi(i, j + 1, k &
                      &- 2))
                  piDEdgeF = 0.5 * (var%pi(i, j, k - 1) + var%pi(i, j + 1, k &
                      &- 1))
                  piEdgeF = 0.5 * (var%pi(i, j, k) + var%pi(i, j + 1, k))
                  piGrad = kappaInv * MaInv2 * pEdgeF / rhov * ((piF - piB) &
                      &/ dy + met23EdgeF * (piDDEdgeF - 4.0 &
                      &* piDEdgeF + 3.0 * piEdgeF) * 0.5 / dz)
                else
                  piUEdgeF = 0.5 * (var%pi(i, j, k + 1) + var%pi(i, j + 1, k &
                      &+ 1))
                  piDEdgeF = 0.5 * (var%pi(i, j, k - 1) + var%pi(i, j + 1, k &
                      &- 1))
                  piGrad = kappaInv * MaInv2 * pEdgeF / rhov * ((piF - piB) &
                      &/ dy + met23EdgeF * (piUEdgeF &
                      &- piDEdgeF) * 0.5 / dz)
                end if
              else
                piGrad = kappaInv * MaInv2 * Pstrat(k) / rhov * (piF - piB) / dy
              end if

              volfcy = 0.0

              ! vstar
              uhorx = 0.25 * (var%u(i - 1, j, k) + var%u(i - 1, j + 1, k) &
                    &+ var%u(i, j, k) + var%u(i, j + 1, k))

              vhory = var%v(i, j, k)

              f_cor_v = 0.5 * (f_cor_nd(j) + f_cor_nd(j + 1))

              if(topography) then
                ! Coriolis force is integrated on LHS.
                vAst = vhory + dt * (- piGrad + volfcy / rhov)
              end if

              if(spongeLayer .and. sponge_uv) then
                if(topography) then
                  vAst = vAst - dt * 0.5 * (kr_sp_tfc(i, j, k) + kr_sp_tfc(i, &
                      &j + 1, k)) * vhory
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
              end if

              piF = var%pi(i, j + 1, k)
                piB = var%pi(i, j, k)

              if(topography) then
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
                      &/ dy + met23EdgeF * (- piUUEdgeF &
                      &+ 4.0 * piUEdgeF - 3.0 * piEdgeF) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeF = 0.5 * (var%pi(i, j, k - 2) + var%pi(i, j + 1, k &
                      &- 2))
                  piDEdgeF = 0.5 * (var%pi(i, j, k - 1) + var%pi(i, j + 1, k &
                      &- 1))
                  piEdgeF = 0.5 * (var%pi(i, j, k) + var%pi(i, j + 1, k))
                  piGradY = kappaInv * MaInv2 * pEdgeF / rhov * ((piF - piB) &
                      &/ dy + met23EdgeF * (piDDEdgeF &
                      &- 4.0 * piDEdgeF + 3.0 * piEdgeF) * 0.5 / dz)
                else
                  piUEdgeF = 0.5 * (var%pi(i, j, k + 1) + var%pi(i, j + 1, k &
                      &+ 1))
                  piDEdgeF = 0.5 * (var%pi(i, j, k - 1) + var%pi(i, j + 1, k &
                      &- 1))
                  piGradY = kappaInv * MaInv2 * pEdgeF / rhov * ((piF - piB) &
                      &/ dy + met23EdgeF * (piUEdgeF &
                      &- piDEdgeF) * 0.5 / dz)
                end if
              end if

              volfcx = 0.0
                volfcy = 0.0

              ! vstar
              uhorx = 0.25 * (var%u(i - 1, j, k) + var%u(i - 1, j + 1, k) &
                    &+ var%u(i, j, k) + var%u(i, j + 1, k))

              vhory = var%v(i, j, k)

              facv = 1.0

              if(spongeLayer .and. sponge_uv) then
                if(topography) then
                  facv = facv + dt * 0.5 * (kr_sp_tfc(i, j, k) + kr_sp_tfc(i, &
                      &j + 1, k))
                end if
              end if

              facu = facv

              if(topography) then
                ! Coriolis force is integrated on LHS.
                vAst = 1.0 / facv * (vhory + dt * (- piGradY + volfcy / rhov))
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
    case default
      stop "momentumPredictor: unknown case zBoundary."
    end select

    if(mmp_mod == "lhs") then
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

            
            
            ! Explicit integration of Coriolis force in TFC.
            if(topography .and. mmp_mod == "lhs") then
              vC = 0.5 * (vOldTFC(i, j, k) + vOldTFC(i, j - 1, k))
              vU = 0.5 * (vOldTFC(i, j, k + 1) + vOldTFC(i, j - 1, k + 1))
              uC = 0.5 * (uOldTFC(i, j, k) + uOldTFC(i - 1, j, k))
              uU = 0.5 * (uOldTFC(i, j, k + 1) + uOldTFC(i - 1, j, k + 1))
              
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
            if(mmp_mod == "lhs" .and. topography) then
              F = - fluxDiff + volForce
            else if(mmp_mod == "lhs") then
              F = - fluxDiff
            else
              stop 'ERROR: wrong mmp_mod'
            end if

            ! interpolated densities
            select case(model)

            case("pseudo_incompressible")
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
              end if

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
              end if

              
              !--- pressure gradient term -> piGrad
              piU = var%pi(i, j, k + 1)
                piD = var%pi(i, j, k)

              if(topography) then
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
                piREdgeU = (jac(i + 1, j, k + 1) * var%pi(i + 1, j, k) &
                    &+ jac(i + 1, j, k) * var%pi(i + 1, j, k + 1)) / (jac(i &
                    &+ 1, j, k) + jac(i + 1, j, k + 1))
                piLEdgeU = (jac(i - 1, j, k + 1) * var%pi(i - 1, j, k) &
                    &+ jac(i - 1, j, k) * var%pi(i - 1, j, k + 1)) / (jac(i &
                    &- 1, j, k) + jac(i - 1, j, k + 1))
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
              end if

              volfcz = 0.0

              ! wstar
              wvert = var%w(i, j, k)

              if(topography) then
                  buoy = - g_ndim * (jac(i, j, k + 1) * rhopOld(i, j, k) &
                      &/ rho000 / jac(i, j, k) + jac(i, j, k) * rhopOld(i, j, &
                      &k + 1) / rho001 / jac(i, j, k + 1)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))
                end if

             wAst = wvert + dt * (buoy - piGrad + volfcz / rhow)

              if(spongeLayer) then
                if(topography) then
                  wAst = wAst - dt * (jac(i, j, k + 1) * kr_sp_w_tfc(i, j, k) &
                      &+ jac(i, j, k) * kr_sp_w_tfc(i, j, k + 1)) / (jac(i, j, &
                      &k) + jac(i, j, k + 1)) * wvert
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

        

        do k = k0, k1
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
              end if

              

              !--- pressure gradient term -> piGrad
              piU = var%pi(i, j, k + 1)
                piD = var%pi(i, j, k)

              if(topography) then
                
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
                piREdgeU = (jac(i + 1, j, k + 1) * var%pi(i + 1, j, k) &
                    &+ jac(i + 1, j, k) * var%pi(i + 1, j, k + 1)) / (jac(i &
                    &+ 1, j, k) + jac(i + 1, j, k + 1))
                piLEdgeU = (jac(i - 1, j, k + 1) * var%pi(i - 1, j, k) &
                    &+ jac(i - 1, j, k) * var%pi(i - 1, j, k + 1)) / (jac(i &
                    &- 1, j, k) + jac(i - 1, j, k + 1))
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
              end if

              
              ! Metric terms due to topography growth.
              volfcz = 0.0

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
              end if

              facw = 1.0


              if(spongeLayer) then
                if(topography) then
                  facw = facw + dt * (jac(i, j, k + 1) * kr_sp_w_tfc(i, j, k) &
                      &+ jac(i, j, k) * kr_sp_w_tfc(i, j, k + 1)) / (jac(i, j, &
                      &k) + jac(i, j, k + 1))
                end if
              end if

              heat0 = heat(i, j, k)

              heat1 = heat(i, j, k + 1)

              ! Buoyancy is predicted after momentum in implicit steps.
              if(topography) then
                buoy = - g_ndim * (jac(i, j, k + 1) * var%rhop(i, j, k) &
                      &/ rho000 / jac(i, j, k) + jac(i, j, k) * var%rhop(i, j, &
                      &k + 1) / rho001 / jac(i, j, k + 1)) / (jac(i, j, k) &
                      &+ jac(i, j, k + 1))
              end if

              if(topography) then
                uC = 0.5 * (var%u(i, j, k) + var%u(i - 1, j, k))
                  uU = 0.5 * (var%u(i, j, k + 1) + var%u(i - 1, j, k + 1))
                  vC = 0.5 * (var%v(i, j, k) + var%v(i, j - 1, k))
                  vU = 0.5 * (var%v(i, j, k + 1) + var%v(i, j - 1, k + 1))

                wAst = 1.0 / (facw + rhoStratEdgeU / rhow * bvsstw * dt &
                      &** 2.0) * (wvert - dt * piGrad + dt * buoy + dt &
                      &* volfcz / rhow + rhoStratEdgeU / rhow * bvsstw * dt &
                      &** 2.0 * (jac(i, j, k + 1) * (met(i, j, k, 1, 3) * uC &
                      &+ met(i, j, k, 2, 3) * vC) + jac(i, j, k) * (met(i, j, &
                      &k + 1, 1, 3) * uU + met(i, j, k + 1, 2, 3) * vU)) &
                      &/ (jac(i, j, k) + jac(i, j, k + 1)))
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

    

            ! update: q(m-1) -> q(m)
            q(i, j, k) = dt * F + alphaRK(m) * q(i, j, k)

            ! update density
            var%rho(i, j, k) = var%rho(i, j, k) + betaRK(m) * q(i, j, k)
          end do
        end do
      end do
    else if(upd_var == "rhop") then
      if(upd_mod == "lhs") then
        if(int_mod /= "expl") then
          stop 'ERROR: wrong int_mod for upd_mod = lhs'
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
              F = - fluxDiff

              
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
                end if

                if(topography) then
                  ! Momentum is predicted before buoyancy in implicit
                  ! steps.
                  wvrt = 0.5 * (wOldTFC(i, j, k) + wOldTFC(i, j, k - 1))
                end if

                if(topography) then
                  ! Compute P coefficients.
                  pEdgeU = (jac(i, j, k + 1) * pStratTFC(i, j, k) + jac(i, j,&
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
                end if

                facw = 1.0

                if(spongeLayer) then
                  if(topography) then
                    facw = facw + dt * kr_sp_w_tfc(i, j, k)
                  end if
                end if

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
                end if

                if(topography) then
                  wvrt = 0.5 * (vertWindTFC(i, j, k, var) + vertWindTFC(i, &
                        &j, k - 1, var))
                end if

                if(topography) then
                  buoy = - g_ndim * rhop / rho
                  buoy = buoy - dt * rhoStratTFC(i, j, k) / rho &
                        &* bvsStratTFC(i, j, k) * wvrt
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
    else
      stop 'upd_var unknown'
    end if

  end subroutine massUpdate

  !-------------------------------------------------------------------------

  subroutine timestep(var, dt, errFlag)
    !---------------------------------------------
    ! compute time step from stability criteria:
    ! 1) CFL criterion for advection
    ! 2) von Neumann cirterion for dissipation
    ! 3) set maximum time step
    ! 4) buouyancy acceleration ->  1/2 b_max*dt^2 < dz
    ! 5) gravity-wave group condition: c_g_x * dt < dx, ...dy,dz
    !---------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var
    real, intent(out) :: dt
    logical, intent(out) :: errFlag

    ! locals
    real :: uMax, vMax, wMax
    real :: dtConv, dtVisc
    real :: dtConv_loc, dtVisc_loc
    real :: dtMax
    real :: dtWave

    ! Buoyancy time step restriction
    real :: dtBuoy, dtBuoy_loc
    real, dimension(3) :: bMax, bMaxNew, duMax
    real :: buoyMax, buoyMin, buoyMaxNew, the_New, the_max, the_min
    real :: bb, ww

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

      case("pseudo_incompressible")

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
        !  Acceleration condition
        !---------------------------
        bMax = 0.0
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx

              select case(model)

              case("pseudo_incompressible")
                if(topography) then
                  bMaxNew = abs(var%rho(i, j, k)) / (rhoStratTFC(i, j, k) &
                      &+ var%rho(i, j, k)) * vertical
                end if

              case default
                stop "timeStep: unknown case model."
              end select

              if(bMaxNew(1) > bMax(1)) bMax(1) = bMaxNew(1)
              if(bMaxNew(2) > bMax(2)) bMax(2) = bMaxNew(2)
              if(bMaxNew(3) > bMax(3)) bMax(3) = bMaxNew(3)
            end do
          end do
        end do
        bMax = FrInv2 * bMax

        ! check whether acceleration condition is needed
        duMax = bMax * dtConv
        if((duMax(1) > 1.e-2 * uMax .or. duMax(2) > 1.e-2 * vMax .or. duMax(3) &
            &> 1.e-2 * wMax) .and. ((bMax(1) /= 0.) .and. (bMax(2) /= 0.) &
            &.and. (bMax(3) /= 0.))) then

          dtBuoy_loc = max(- uMax / bMax(1) + sqrt((uMax / bMax(1)) ** 2 + 2. &
              &* cfl * dx / bMax(1)), - vMax / bMax(2) + sqrt((vMax / bMax(2)) &
              &** 2 + 2. * cfl * dy / bMax(2)), - wMax / bMax(3) + sqrt((wMax &
              &/ bMax(3)) ** 2 + 2. * cfl * dz / bMax(3)))

          if(dtBuoy_loc * tRef < 1.e-2) then
            print *, "dtBuoy_loc*tRef  = ", dtBuoy_loc * tRef
            print *, "bMax(3) = ", bMax(3) * FrInv2
          end if

          if(topography) then
            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  select case(model)
                  case("pseudo_incompressible")
                    bb = abs(var%rho(i, j, k)) / (rhoStratTFC(i, j, k) &
                        &+ var%rho(i, j, k))
                  case default
                    stop "Error in timestep: Unknown case model!"
                  end select
                  ww = abs(0.5 * (vertWindTFC(i, j, k, var) + vertWindTFC(i, &
                      &j, k - 1, var))) + small
                  dtBuoy_loc = min(dtBuoy_loc, - ww / bb + sqrt((ww / bb) &
                      &** 2.0 + 2.0 * cfl * jac(i, j, k) * dz / bb))
                end do
              end do
            end do
          end if
        else

          dtBuoy_loc = 1.0e20 / tRef ! set to high value if not needed

        end if

        ! find global minimum

        call mpi_reduce(dtBuoy_loc, dtBuoy, 1, mpi_double_precision, mpi_min, &
            &root, comm, ierror)

        call mpi_bcast(dtBuoy, 1, mpi_double_precision, root, comm, ierror)

        !---------------------------
        !   von Neumann condition
        !----------------------------

        dtVisc = 0.5 * min(dx ** 2, dy ** 2, dz ** 2) * Re

        if(topography) then
          dtVisc_loc = dtVisc
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                dtVisc_loc = min(dtVisc_loc, 0.5 * (jac(i, j, k) * dz) ** 2.0 &
                    &* Re)
              end do
            end do
          end do
          call mpi_reduce(dtVisc_loc, dtVisc, 1, mpi_double_precision, &
              &mpi_min, root, comm, ierror)
          call mpi_bcast(dtVisc, 1, mpi_double_precision, root, comm, ierror)
        end if

        !----------------------------
        !    Maximal time step
        !----------------------------
        dtMax = dtMax_dim / tRef

        !------------------------------------
        !    Gravity wave time period
        !------------------------------------

        dtWave = 1. / (NN + small)

        
        !-------------------------------
        !        Make your choice
        !-------------------------------

        dt = min(dtVisc, dtConv, dtMax, dtBuoy)

        !-----------------------------------------
        !     Inform on time step restrictions
        !-----------------------------------------

        if(master) then

          write(*, fmt = "(a25,es15.1,a8)") "dtVisc =", dtVisc * tRef, "seconds"
          write(*, fmt = "(a25,es15.1,a8)") "dtConv =", dtConv * tRef, "seconds"
          write(*, fmt = "(a25,es15.1,a8)") "dtMax =", dtMax * tRef, "seconds"
          write(*, fmt = "(a25,es15.1,a8)") "dtBuoy =", dtBuoy * tRef, "seconds"
          write(*, fmt = "(a25,es15.1,a8)") "dtWave =", dtWave * tRef, "seconds"
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
          else if(dt == dtBuoy) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtBuoy = ", dt * tRef, &
                &"seconds"
          else if(dt == dtWave) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtWave = ", dt * tRef, &
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

end module update_module
