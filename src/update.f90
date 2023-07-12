module update_module

  use type_module
  use timeScheme_module
  use atmosphere_module
  use flux_module
  use algebra_module
  use ice_module
  use poisson_module
  use boundary_module
  use mpi_module
  use output_module

  implicit none

  private
  ! all module objects (variables, functions, subroutines)
  ! are internal to the module by default

  !------------------------
  !   public subroutines
  !------------------------
  public :: momentumPredictor
  public :: massUpdate
  public :: iceUpdate
  public :: thetaUpdate
  public :: tracerUpdate
  public :: timestep
  public :: init_update
  public :: set_spongeLayer
  public :: CoefDySma_update
  public :: Var3DSmthDySma
  public :: ice2Update, ice2Update_source, ice2Update_apb
  public :: setHaloAndBoundary

  public :: smooth_shapiro
  public :: smooth_hor_shapiro

  public :: BGstate_update

  ! TFC FJ
  public :: momentumPredictorTestTFC, massUpdateTestTFC

  !-------------------------------
  !    private module variables
  !------------------------------

  contains

  subroutine set_spongeLayer(var, dt, variable)
    !--------------------------------------
    ! relaxes the predicted solution to
    ! the background state
    !--------------------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var
    real, intent(in) :: dt
    character(len = *), intent(in) :: variable

    ! local variables
    integer :: i, j, k, iVar

    ! relaxation parameters
    real :: alpha, beta
    real :: spongeAlphaZ, spongeDz

    ! variables for rho
    real :: rho_old, rho_bg, rho_new
    real :: uOld, uBG, uNew
    real :: vOld, vBG, vNew
    real :: wOld, wBG, wNew

    ! variables for ice
    real :: nAer_bg, nIce_bg, qIce_bg, qv_bg
    real :: T, p

    real, dimension(1:nz) :: sum_local, sum_global

    real, dimension(1:ny) :: c4_strtd
    real :: yjets, yjetn, dy_hs
    real :: ymin, ymax, yloc, jwdth
    ! integer :: j00

    ! TFC FJ
    integer :: i00, j00
    real :: spongeAlphaX, spongeAlphaY

    ! return if sponge layer with relaxation is switched off
    if(.not. spongeLayer) then
      return
    end if

    ! nondimensionalize relaxation parameter
    spongeAlphaZ = spongeAlphaZ_dim * tRef

    ! thickness of sponge layer
    spongeDz = z(nz) - z(kSponge)

    ! TFC FJ
    ! Meridional dependence is only implemented for semi-implicit procedure!
    c4_strtd = 1.0

    ! TFC FJ
    ! Define parameters needed for TFC sponge layers.
    if(topography .and. spongeTFC .and. lateralSponge) then
      i00 = is + nbx - 1
      j00 = js + nby - 1
      spongeAlphaX = spongeAlphaZ
      spongeAlphaY = spongeAlphaZ
    end if

    select case(variable)

      !UAB 200413
    case("ref")
      ! save total density and subtract the reference-atmosphere density
      ! from this again after the update of the latter

      if(fluctuationMode) then
        do k = 1, nz
          var(:, :, k, 1) = var(:, :, k, 1) + rhoStrat(k)
        end do
      end if

      if((timeScheme == "semiimplicit") .or. auxil_equ) then
        do k = 1, nz
          var(:, :, k, 6) = var(:, :, k, 6) + rhoStrat(k)
        end do
      end if

      do k = kSponge, nz

        alpha = spongeAlphaZ * ((z(k) - zSponge) / spongeDz) ** sponge_order
        beta = 1. / (1. + alpha * 0.5 * dt) ** 2

        rhoStrat(k) = (1. - beta) * rhoStrat_0(k) + beta * rhoStrat(k)
        pStrat(k) = (1. - beta) * pStrat_0(k) + beta * pStrat(k)

        thetaStrat(k) = PStrat(k) / rhoStrat(k)
      end do

      do k = - 1, nz + 1
        PstratTilde(k) = 0.5 * (PStrat(k) + PStrat(k + 1))
        rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k + 1))
        thetaStratTilde(k) = PStratTilde(k) / rhoStratTilde(k)
      end do

      ! adjust density fluctuations to new reference atmosphere
      if(fluctuationMode) then
        do k = 1, nz
          var(:, :, k, 1) = var(:, :, k, 1) - rhoStrat(k)
        end do
      end if

      if((timeScheme == "semiimplicit") .or. auxil_equ) then
        do k = 1, nz
          var(:, :, k, 6) = var(:, :, k, 6) - rhoStrat(k)
        end do
      end if

      ! update of non-dimensional squared Brunt-Vaisala frequency
      ! (this could perhaps be done a bit nicer)

      bvsStrat(- 1) = g_ndim / thetaStrat(0) * (thetaStrat(1) - thetaStrat(0)) &
          / dz

      bvsStrat(0) = g_ndim / thetaStrat(0) * (thetaStrat(1) - thetaStrat(0)) &
          / dz

      N2 = max(bvsStrat(- 1), bvsStrat(0))

      do k = 1, nz
        bvsStrat(k) = g_ndim / thetaStrat(k) * (thetaStrat(k + 1) &
            - thetaStrat(k - 1)) / (2.0 * dz)

        N2 = max(N2, bvsStrat(k))
      end do

      bvsStrat(nz + 1) = g_ndim / thetaStrat(nz + 1) * (thetaStrat(nz + 1) &
          - thetaStrat(nz)) / dz

      N2 = max(N2, bvsStrat(nz + 1))

      if(N2 < 0.) then
        stop 'ERROR: N2 < 0'
      else
        NN = sqrt(N2)
      end if

      !testb
      do k = - 1, nz + 1
        if(master .and. N2 == bvsStrat(k)) print *, 'N2 = max at k =', k
      end do
      !teste
      !UAE 200413

    case("rho")

      ! TFC FJ
      ! No sponge applied to any quantity other than w.
      ! No sponge applied to rho in Boussinesq model.
      if(.not. spongeTFC .and. model /= "Boussinesq") then
        do k = kSponge, nz
          do j = 1, ny
            do i = 1, nx

              if((TestCase == "baroclinic_LC") .or. (TestCase &
                  == "baroclinic_ID")) then
                !UAB 200413
                !rho_bg = dens_env_pp(i, j, k)
                if(fluctuationMode) then
                  if(topography) then
                    ! TFC FJ
                    rho_bg = dens_env_pp(i, j, k) - rhoStratTFC(i, j, k)
                  else
                    rho_bg = dens_env_pp(i, j, k) - rhoStrat(k)
                  end if
                else
                  rho_bg = dens_env_pp(i, j, k)
                end if
                !UAE 200413
              else
                if(fluctuationMode) then
                  rho_bg = 0.0 ! push back to zero perturbation
                else
                  rho_bg = rhoStrat(k)
                end if
              end if

              rho_old = var(i, j, k, 1)
              if(diffusive_sponge) then
                alpha = spongeAlphaZ * exp((z(k) - lz(1)) / zSponge)
              else
                alpha = c4_strtd(j) * spongeAlphaZ * ((z(k) - zSponge) &
                    / spongeDz) ** sponge_order
              end if
              beta = 1. / (1. + alpha * 0.5 * dt) ** 2
              rho_new = (1. - beta) * rho_bg + beta * rho_old

              var(i, j, k, 1) = rho_new

            end do
          end do
        end do
      end if

    case("rhop")

      ! TFC FJ
      ! No sponge applied to any quantity other than w.
      if(.not. spongeTFC) then
        do k = kSponge, nz
          do j = 1, ny
            do i = 1, nx

              if((TestCase == "baroclinic_LC") .or. (TestCase &
                  == "baroclinic_ID")) then
                if(topography) then
                  ! TFC FJ
                  rho_bg = dens_env_pp(i, j, k) - rhoStratTFC(i, j, k)
                else
                  rho_bg = dens_env_pp(i, j, k) - rhoStrat(k)
                end if
              else
                rho_bg = 0.0 ! push back to zero perturbation
              end if

              rho_old = var(i, j, k, 6)
              if(diffusive_sponge) then
                alpha = spongeAlphaZ * exp((z(k) - lz(1)) / zSponge)
              else
                alpha = c4_strtd(j) * spongeAlphaZ * ((z(k) - zSponge) &
                    / spongeDz) ** sponge_order
              end if
              beta = 1. / (1. + alpha * 0.5 * dt) ** 2
              rho_new = (1. - beta) * rho_bg + beta * rho_old

              var(i, j, k, 6) = rho_new

            end do
          end do
        end do
      end if

    case("ice")

      nAer_bg = 0.0 !init_nAer * rhoRef * lRef**3
      nIce_bg = 0.0
      qIce_bg = 0.0
      qv_bg = 0.0

      do k = kSponge, nz
        do j = 1, ny
          do i = 1, nx

            select case(iceTestcase)
            case("homogeneous_qv")
              qv_bg = 0.0 !init_qv
            case("homogeneous_SIce")
              !call find_temperature(T,i,j,k,var)
              !p = press0_dim * ( (PStrat(k)/p0)**gamma_1  +var(i,j,k,5) )**kappaInv
              qv_bg = 0.0 !epsilon0 * init_SIce * p_saturation(T) / p
            end select

            if(diffusive_sponge) then
              alpha = spongeAlphaZ * exp((z(k) - lz(1)) / zSponge)
            else
              alpha = spongeAlphaZ * ((z(k) - zSponge) / spongeDz) &
                  ** sponge_order
            end if
            beta = 1. / (1. + alpha * 0.5 * dt) ** 2
            var(i, j, k, nVar - 3) = (1. - beta) * nAer_bg + beta * var(i, j, &
                k, nVar - 3)
            var(i, j, k, nVar - 2) = (1. - beta) * nIce_bg + beta * var(i, j, &
                k, nVar - 2)
            var(i, j, k, nVar - 1) = (1. - beta) * qIce_bg + beta * var(i, j, &
                k, nVar - 1)
            var(i, j, k, nVar) = (1. - beta) * qv_bg + beta * var(i, j, k, nVar)
            do iVar = 0, 3
              if(var(i, j, k, nVar - iVar) .lt. 0.0) var(i, j, k, nVar - iVar) &
                  = 0.0
            end do

          end do
        end do
      end do

      ! second sponge for ice particles at lower boundary
      ! didn't prove useful
      !
      !if (iceTestcase == "qv_relaxation") then
      !  do k = 0, ceiling(mountainHeight_dim/(dz*lRef))+5
      !    do j = 1,ny
      !      do i = 1,nx
      !        qv_bg = 0.0
      !        alpha = 100*spongeAlphaZ*(1-z(k)/(mountainHeight_dim/lRef+5*dz))
      !        beta = 1./(1.+alpha*0.5*dt)**2
      !        var(i,j,k,nVar-3) = (1.-beta)*nAer_bg + beta*var(i,j,k,nVar-3)
      !        var(i,j,k,nVar-2) = (1.-beta)*nIce_bg + beta*var(i,j,k,nVar-2)
      !        var(i,j,k,nVar-1) = (1.-beta)*qIce_bg + beta*var(i,j,k,nVar-1)
      !        var(i,j,k,nVar)   = (1.-beta) * qv_bg + beta*var(i,j,k,nVar)
      !        do iVar=0,3
      !           if (var(i,j,k,nVar-iVar) .lt. 0.0) var(i,j,k,nVar-iVar) = 0.0
      !        end do
      !      end do
      !    end do
      !  end do
      !end if

    case("uvw")
      ! relax u to:
      !   baroclinic cases (2D or 3D):
      !       0 or
      !       environmental u or
      !       u (no relaxation)
      !   else: horizontal mean

      ! local horizontal sum in the sponge layer

      do k = kSponge, nz
        sum_local(k) = sum(var(1:nx, 1:ny, k, 2))
      end do

      ! global sum and average

      call mpi_allreduce(sum_local(kSponge), sum_global(kSponge), nz - kSponge &
          + 1, mpi_double_precision, mpi_sum, comm, ierror)
      sum_global = sum_global / (sizeX * sizeY)

      ! TFC FJ
      ! No sponge applied to any quantity other than w.
      if(.not. spongeTFC) then
        do k = kSponge, nz
          do j = 1, ny
            !do i = 0,nx
            do i = 1, nx
              if((TestCase == "baroclinic_LC") .or. (TestCase &
                  == "baroclinic_ID")) then
                if(Sponge_Rel_Bal_Type == "hyd") then
                  ! relax to hydrost bal
                  uBG = 0.
                else if(Sponge_Rel_Bal_Type == "env") then
                  ! relax to geostr bal
                  uBG = u_env_pp(i, j, k)
                else
                  ! free development
                  uBG = var(i, j, k, 2)
                end if
              else
                uBG = sum_global(k)
              end if

              uOld = var(i, j, k, 2)
              if(diffusive_sponge) then
                alpha = spongeAlphaZ * exp((z(k) - lz(1)) / zSponge)
              else
                alpha = c4_strtd(j) * spongeAlphaZ * ((z(k) - zSponge) &
                    / spongeDz) ** sponge_order
              end if
              beta = 1. / (1. + alpha * 0.5 * dt) ** 2
              uNew = (1. - beta) * uBG + beta * uOld

              var(i, j, k, 2) = uNew
            end do
          end do
        end do
      end if

      ! relax v to:
      !   baroclinic cases (2D or 3D):
      !       0 or
      !       environmental v or
      !       v (no relaxation)
      !   else: horizontal mean

      ! local horizontal sum in the sponge layer

      do k = kSponge, nz
        sum_local(k) = sum(var(1:nx, 1:ny, k, 3))
      end do

      ! global sum and average

      call mpi_allreduce(sum_local(kSponge), sum_global(kSponge), nz - kSponge &
          + 1, mpi_double_precision, mpi_sum, comm, ierror)
      sum_global = sum_global / (sizeX * sizeY)

      ! TFC FJ
      ! No sponge applied to any quantity other than w.
      if(.not. spongeTFC) then
        do k = kSponge, nz
          !do j = 0,ny !gaga 1->0
          do j = 1, ny
            do i = 1, nx
              if((TestCase == "baroclinic_LC") .or. (TestCase &
                  == "baroclinic_ID")) then
                if(Sponge_Rel_Bal_Type == "hyd") then
                  ! relax to hydrost bal
                  vBG = 0.
                else if(Sponge_Rel_Bal_Type == "env") then
                  ! relax to geostr bal
                  vBG = v_env_pp(i, j, k)
                else
                  ! free development
                  vBG = var(i, j, k, 3)
                end if
              else
                vBG = sum_global(k)
              end if

              vOld = var(i, j, k, 3)
              if(diffusive_sponge) then
                alpha = spongeAlphaZ * exp((z(k) - lz(1)) / zSponge)
              else
                alpha = c4_strtd(j) * spongeAlphaZ * ((z(k) - zSponge) &
                    / spongeDz) ** sponge_order
              end if
              beta = 1. / (1. + alpha * 0.5 * dt) ** 2
              vNew = (1. - beta) * vBG + beta * vOld

              var(i, j, k, 3) = vNew

            end do
          end do
        end do
      end if

      ! achatzb relax w to zero

      if(topography .and. spongeTFC) then
        ! TFC FJ
        ! TFC sponge layers.
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              wOld = var(i, j, k, 4)
              alpha = 0.0
              if(lateralSponge) then
                ! Zonal sponge.
                if(x(i00 + i) <= xSponge0) then
                  alpha = alpha + spongeAlphaX * sin(0.5 * pi * (xSponge0 &
                      - x(i00 + i)) / (xSponge0 - lx(0))) ** 2.0
                else if(x(i00 + i) >= xSponge1) then
                  alpha = alpha + spongeAlphaX * sin(0.5 * pi * (x(i00 + i) &
                      - xSponge1) / (lx(1) - xSponge1)) ** 2.0
                end if
                ! Meridional sponge.
                if(y(j00 + j) <= ySponge0) then
                  alpha = alpha + spongeAlphaY * sin(0.5 * pi * (ySponge0 &
                      - y(j00 + j)) / (ySponge0 - ly(0))) ** 2.0
                else if(y(j00 + j) >= ySponge1) then
                  alpha = alpha + spongeAlphaY * sin(0.5 * pi * (y(j00 + j) &
                      - ySponge1) / (ly(1) - ySponge1)) ** 2.0
                end if
              end if
              ! Vertical sponge.
              if(heightTFC(i, j, k) >= zSponge) then
                alpha = alpha + spongeAlphaZ * sin(0.5 * pi * (heightTFC(i, j, &
                    k) - zSponge) / (lz(1) - zSponge)) ** 2.0
              end if
              ! Adjust for terrain-following velocity.
              alpha = alpha / jac(i, j, k)
              ! Apply total sponge.
              beta = 1.0 / (1.0 + alpha * dt)
              wNew = beta * wOld
              var(i, j, k, 4) = wNew
            end do
          end do
        end do
      else
        do k = kSponge, nz
          wBG = backgroundFlow_dim(3) / uRef !0.0
          do j = 1, ny
            do i = 1, nx
              wOld = var(i, j, k, 4)
              if(diffusive_sponge) then
                alpha = spongeAlphaZ * exp((z(k) - lz(1)) / zSponge)
              else
                alpha = c4_strtd(j) * spongeAlphaZ * ((z(k) - zSponge) &
                    / spongeDz) ** sponge_order
              end if
              beta = 1. / (1. + alpha * 0.5 * dt) ** 2
              wNew = (1. - beta) * wBG + beta * wOld

              var(i, j, k, 4) = wNew
            end do
          end do
        end do
      end if

    case default
      stop "spongeLayer: Unknown variable"
    end select

  end subroutine set_spongeLayer

  !---------------------------------------------------------------------

  subroutine wind_ip(var, xip, yip, zip, wind_recon, u_ip, v_ip, w_ip)

    !------------------------------------------------------------------
    ! winds at interpolation points needed for the implementation of the
    ! immersed-boundary topographic boundary condition
    !------------------------------------------------------------------

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(in) :: var
    real, intent(in) :: xip, yip, zip ! coordinates of interpolation point
    ! wind to be reconstructed using the winds at the interpolation point
    ! (u,v,w)
    character(len = *), intent(in) :: wind_recon
    real, intent(out) :: u_ip, v_ip, w_ip ! interpolated winds

    integer :: i00, j00
    integer :: i, j, k

    i00 = is + nbx - 1
    j00 = js + nby - 1

    if(wind_recon == 'u' .or. wind_recon == 'v') then
      !----------------------------------------
      ! interpolation for reconstructing u or v
      !----------------------------------------

      k = nint((zip - lz(0) + dz / 2.0) / dz)

      ! u at interpolation point

      i = floor((xip - lx(0)) / dx) - i00
      j = floor((yip - ly(0) + dy / 2.) / dy) - j00

      u_ip = ((x(i + i00 + 1) + dx / 2. - xip) * ((y(j + j00 + 1) - yip) &
          * var(i, j, k, 2) + (yip - y(j + j00)) * var(i, j + 1, k, 2)) + (xip &
          - x(i + i00) - dx / 2.) * ((y(j + j00 + 1) - yip) * var(i + 1, j, k, &
          2) + (yip - y(j + j00)) * var(i + 1, j + 1, k, 2))) / (dx * dy)

      ! v at interpolation point

      i = floor((xip - lx(0) + dx / 2.) / dx) - i00
      j = floor((yip - ly(0)) / dy) - j00

      v_ip = ((x(i + i00 + 1) - xip) * ((y(j + j00 + 1) + dy / 2. - yip) &
          * var(i, j, k, 3) + (yip - y(j + j00) - dy / 2.) * var(i, j + 1, k, &
          3)) + (xip - x(i + i00)) * ((y(j + j00 + 1) + dy / 2. - yip) * var(i &
          + 1, j, k, 3) + (yip - y(j + j00) - dy / 2.) * var(i + 1, j + 1, k, &
          3))) / (dx * dy)

      ! w at interpolation point
      ! here and further below the vertical averaging could be coded
      ! more efficiently by averaging var directly (instead of the
      ! horizontally interpolated winds)

      i = floor((xip - lx(0) + dx / 2.) / dx) - i00
      j = floor((yip - ly(0) + dy / 2.) / dy) - j00

      w_ip = 0.5 * (((x(i + i00 + 1) - xip) * ((y(j + j00 + 1) - yip) * var(i, &
          j, k, 4) + (yip - y(j + j00)) * var(i, j + 1, k, 4)) + (xip - x(i &
          + i00)) * ((y(j + j00 + 1) - yip) * var(i + 1, j, k, 4) + (yip - y(j &
          + j00)) * var(i + 1, j + 1, k, 4))) / (dx * dy) + ((x(i + i00 + 1) &
          - xip) * ((y(j + j00 + 1) - yip) * var(i, j, k - 1, 4) + (yip - y(j &
          + j00)) * var(i, j + 1, k - 1, 4)) + (xip - x(i + i00)) * ((y(j &
          + j00 + 1) - yip) * var(i + 1, j, k - 1, 4) + (yip - y(j + j00)) &
          * var(i + 1, j + 1, k - 1, 4))) / (dx * dy))
    else if(wind_recon == 'w') then
      !----------------------------------------
      ! interpolation for reconstructing w
      !----------------------------------------

      k = nint((zip - lz(0)) / dz)

      ! u at interpolation point

      i = floor((xip - lx(0)) / dx) - i00
      j = floor((yip - ly(0) + dy / 2.) / dy) - j00

      u_ip = 0.5 * (((x(i + i00 + 1) + dx / 2. - xip) * ((y(j + j00 + 1) &
          - yip) * var(i, j, k, 2) + (yip - y(j + j00)) * var(i, j + 1, k, 2)) &
          + (xip - x(i + i00) - dx / 2.) * ((y(j + j00 + 1) - yip) * var(i &
          + 1, j, k, 2) + (yip - y(j + j00)) * var(i + 1, j + 1, k, 2))) / (dx &
          * dy) + ((x(i + i00 + 1) + dx / 2. - xip) * ((y(j + j00 + 1) - yip) &
          * var(i, j, k + 1, 2) + (yip - y(j + j00)) * var(i, j + 1, k + 1, &
          2)) + (xip - x(i + i00) - dx / 2.) * ((y(j + j00 + 1) - yip) * var(i &
          + 1, j, k + 1, 2) + (yip - y(j + j00)) * var(i + 1, j + 1, k + 1, &
          2))) / (dx * dy))

      ! v at interpolation point

      i = floor((xip - lx(0) + dx / 2.) / dx) - i00
      j = floor((yip - ly(0)) / dy) - j00

      v_ip = 0.5 * (((x(i + i00 + 1) - xip) * ((y(j + j00 + 1) + dy / 2. &
          - yip) * var(i, j, k, 3) + (yip - y(j + j00) - dy / 2.) * var(i, j &
          + 1, k, 3)) + (xip - x(i + i00)) * ((y(j + j00 + 1) + dy / 2. - yip) &
          * var(i + 1, j, k, 3) + (yip - y(j + j00) - dy / 2.) * var(i + 1, j &
          + 1, k, 3))) / (dx * dy) + ((x(i + i00 + 1) - xip) * ((y(j + j00 &
          + 1) + dy / 2. - yip) * var(i, j, k + 1, 3) + (yip - y(j + j00) - dy &
          / 2.) * var(i, j + 1, k + 1, 3)) + (xip - x(i + i00)) * ((y(j + j00 &
          + 1) + dy / 2. - yip) * var(i + 1, j, k + 1, 3) + (yip - y(j + j00) &
          - dy / 2.) * var(i + 1, j + 1, k + 1, 3))) / (dx * dy))

      ! w at interpolation point

      i = floor((xip - lx(0) + dx / 2.0) / dx) - i00
      j = floor((yip - ly(0) + dy / 2.0) / dy) - j00

      w_ip = ((x(i + i00 + 1) - xip) * ((y(j + j00 + 1) - yip) * var(i, j, k, &
          4) + (yip - y(j + j00)) * var(i, j + 1, k, 4)) + (xip - x(i + i00)) &
          * ((y(j + j00 + 1) - yip) * var(i + 1, j, k, 4) + (yip - y(j + j00)) &
          * var(i + 1, j + 1, k, 4))) / (dx * dy)
    else
      print *, "wrong wind_recon"
      stop
    endif

  end subroutine wind_ip

  !---------------------------------------------------------------------

  subroutine momentumPredictor_nc(var, flux, force, dt, q, m, mmp_mod, &
      int_mod, facray)
    !-------------------------------------
    !  calculates the velocities u^*
    !
    !  Corilolis effect handled explicitly
    !-------------------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var

    real, dimension(- 1:nx, - 1:ny, - 1:nz, 3, nVar), intent(in) :: flux
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

    real, intent(in) :: dt, facray
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 3), &
        intent(inout) :: q
    integer, intent(in) :: m

    logical :: spongeLayer_s, topography_s

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

    ! classical RK3
    real :: rhoM_0, uM_0, vM_0, wM_0, momM_0

    ! non-dimensional Corilois parameter (= inverse Rossby number)
    !FS real :: f_cor_nd
    real, dimension(0:ny + 1) :: f_cor_nd

    !real :: rho, rhop, rhou, rhov, rhow, facu, facv, facw, facr, pstw, buoy
    real :: rho, rhop, rhou, rhov, rhow, facu, facv, facw, facr, pstw, pstw_0, &
        buoy
    real :: rho10, rho01
    real :: rhov0m, rhov00, rhov1m, rhov10
    real :: rhou00, rhoum0, rhou01, rhoum1
    real :: rho000, rho001
    real :: volfcx, volfcy
    real :: bvsstw

    real :: rho_e, rhou_e, rhov_e, rhow_e, pstw_e
    real :: rho10_e, rho01_e
    real :: rhov0m_e, rhov00_e, rhov1m_e, rhov10_e
    real :: rhou00_e, rhoum0_e, rhou01_e, rhoum1_e
    real :: rho000_e, rho001_e
    real :: piR_e, piL_e, piF_e, piB_e, piU_e, piD_e

    real :: f_cor_v

    real :: u_ip, v_ip, w_ip, u_ip_n, v_ip_n, w_ip_n, u_ip_t, v_ip_t, w_ip_t, &
        u_rp, v_rp, w_rp, u_rp_n, v_rp_n, w_rp_n, u_rp_t, v_rp_t, w_rp_t

    integer :: i00, j00

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
            / (ymax - ymin))
      end do
    else if(corset == 'constant') then
      f_cor_nd(0:ny + 1) = f_Coriolis_dim * tRef
    else
      stop 'ERROR: wrong corset'
    end if

    if(topography) then
      i00 = is + nbx - 1
      j00 = js + nby - 1
    end if

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
        topography_s = topography

        spongeLayer = .false.
        topography = .false.
      else if(int_mod == 'impl') then
        kr_sp = kr_sp * facray
        kr_sp_w = kr_sp_w * facray
        alprlx = alprlx * facray
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
            fR = flux(i, j, k, 1, 2)
            fL = flux(i - 1, j, k, 1, 2)
            gF = flux(i, j, k, 2, 2)
            gB = flux(i, j - 1, k, 2, 2)
            hU = flux(i, j, k, 3, 2)
            hD = flux(i, j, k - 1, 3, 2)
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz ! diverg.

            volForce = 0.

            if(mmp_mod == "tot") then
              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piR = var(i + 1, j, k, 5) - var_env(i + 1, j, k, 5)
                piL = var(i, j, k, 5) - var_env(i, j, k, 5)
              else
                piR = var(i + 1, j, k, 5)
                piL = var(i, j, k, 5)
              end if

              piGrad = kappaInv * MaInv2 * Pstrat(k) * (piR - piL) / dx

              if(TestCase == "baroclinic_LC") then !FS
                piGrad = piGrad + kappaInv * MaInv2 * (PStrat(k) &
                    - pStrat_0(k)) * (var_env(i + 1, j, k, 5) - var_env(i, j, &
                    k, 5)) / dx
              end if

              !---- volume forces
              volForce = 0.5 * (force(i, j, k, 1) + force(i + 1, j, k, 1))

              if(TestCase == "baroclinic_LC") then
                if(model == "pseudo_incompressible") then
                  rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i + 1, j, k))

                  if(fluctuationMode) then
                    rhoM_1 = rhoM_1 + rhoStrat(k)
                  end if
                else if(model == "Boussinesq") then
                  rhoM_1 = rho00
                else
                  stop "momentumPredictor: unkown model."
                end if

                volforce = volforce - rhoM_1 * RoInv(j) * 0.25 * (var_env(i, j &
                    - 1, k, 3) + var_env(i + 1, j - 1, k, 3) + var_env(i, j, &
                    k, 3) + var_env(i + 1, j, k, 3))

              end if

              if(topography) then
                ! Rayleigh damping for topography (immersed boundary)

                if(model == "pseudo_incompressible") then
                  rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i + 1, j, k))

                  if(fluctuationMode) then
                    rhoM_1 = rhoM_1 + rhoStrat(k)
                  end if
                else if(model == "Boussinesq") then
                  rhoM_1 = rho00
                else
                  stop "momentumPredictor: unkown model."
                end if

                !UAC if(topography_mask(i00+i,j00+j,k)&
                !   .or.&
                !   topography_mask(i00+i+1,j00+j,k)) then
                !   volForce = volForce - alprlx * rhoM_1*var(i,j,k,2)
                !end if
                if(k < kbl_topo(i, j, 1)) then
                  volForce = volForce - alprlx * rhoM_1 * var(i, j, k, 2)
                else if(k == kbl_topo(i, j, 1)) then
                  call wind_ip(var, x_ip(i, j, 1), y_ip(i, j, 1), z_ip(i, j, &
                      1), 'u', u_ip, v_ip, w_ip)

                  u_ip_n = (u_ip * dhdx(i, j, 1) + v_ip * dhdy(i, j, 1) &
                      - w_ip) * dhdx(i, j, 1) / (1 + dhdx(i, j, 1) ** 2 &
                      + dhdy(i, j, 1) ** 2)

                  u_ip_t = u_ip - u_ip_n

                  u_rp_t = velocity_reconst_t(i, j, 1) * u_ip_t
                  u_rp_n = velocity_reconst_n(i, j, 1) * u_ip_n
                  u_rp = u_rp_t + u_rp_n

                  volForce = volForce - alprlx * rhoM_1 * (var(i, j, k, 2) &
                      - u_rp)
                end if
                !UAE
              end if
            end if

            if(TestCase == "baroclinic_LC") then
              if(background == "HeldSuarez") then
                ! Rayleigh damping

                if(model == "pseudo_incompressible") then
                  rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i + 1, j, k))

                  if(fluctuationMode) then
                    rhoM_1 = rhoM_1 + rhoStrat(k)
                  end if
                else if(model == "Boussinesq") then
                  rhoM_1 = rho00
                else
                  stop "momentumPredictor: unkown model."
                end if

                volForce = volForce - kv_hs(j, k) * rhoM_1 * (var(i, j, k, 2) &
                    - var_env(i, j, k, 2))
              end if
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

            case("pseudo_incompressible")

              rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i + 1, j, k))
              rhoM = 0.5 * (var(i, j, k, 1) + var(i + 1, j, k, 1))

              if(fluctuationMode) then
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
            uM_1 = var(i, j, k, 2)
            momM_1 = rhoM_1 * uM_1

            ! q(m-1) -> q(m)

            q(i, j, k, 1) = dt * F + alpha(m) * q(i, j, k, 1)

            ! rhoU(m-1) -> rhoU(m)
            momM = momM_1 + beta(m) * q(i, j, k, 1)

            ! calc u(m,*)
            uAst = momM / rhoM

            ! uAst -> var
            var(i, j, k, 2) = uAst
          end do
        end do
      end do
    else if(mmp_mod == "rhs") then
      if(int_mod == "expl") then
        do k = 1, nz
          do j = 1, ny
            do i = i0, i1
              rhou = 0.5 * (var(i, j, k, 1) + var(i + 1, j, k, 1))
              if(fluctuationMode) then
                rhou = rhou + rhoStrat(k)
              end if

              if(TestCase == "baroclinic_LC") then
                rhou_e = 0.5 * (var_env(i, j, k, 1) + var_env(i + 1, j, k, 1))
                if(fluctuationMode) then
                  rhou_e = rhou_e + rhoStrat_0(k)
                end if
              end if

              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piR = var(i + 1, j, k, 5) - var_env(i + 1, j, k, 5)
                piL = var(i, j, k, 5) - var_env(i, j, k, 5)

                piR_e = var_env(i + 1, j, k, 5)
                piL_e = var_env(i, j, k, 5)
              else
                piR = var(i + 1, j, k, 5)
                piL = var(i, j, k, 5)
              end if

              piGrad = kappaInv * MaInv2 * Pstrat(k) / rhou * (piR - piL) / dx

              if(TestCase == "baroclinic_LC") then !FS
                piGrad = piGrad + kappaInv * MaInv2 * (Pstrat(k) / rhou &
                    - Pstrat_0(k) / rhou_e) * (piR_e - piL_e) / dx
              end if

              ! gravity-wave forcing
              if(raytracer .or. (testCase == "mountainwave")) then
                volfcx = 0.5 * (force(i, j, k, 1) + force(i + 1, j, k, 1))
              else
                volfcx = 0.0
              end if

              ! ustar
              if(TestCase == "baroclinic_LC") then
                uhorx = var(i, j, k, 2) - var_env(i, j, k, 2)
              else
                uhorx = var(i, j, k, 2)
              end if

              vhory = 0.25 * (var(i, j - 1, k, 3) + var(i, j, k, 3) + var(i &
                  + 1, j - 1, k, 3) + var(i + 1, j, k, 3))

              uAst = uhorx + dt * (f_cor_nd(j) * vhory - piGrad + volfcx / rhou)

              if(topography) then
                ! Rayleigh damping for topography (immersed boundary)

                !UAC if(topography_mask(i00+i,j00+j,k)&
                !   .or.&
                !   topography_mask(i00+i+1,j00+j,k)) then
                !   uAst = uAst - dt* alprlx*uhorx
                !end if
                if(TestCase == "baroclinic_LC") then
                  stop 'combination of topography with baroclinic  LC not &
                      possible yet'
                end if

                if(k < kbl_topo(i, j, 1)) then
                  uAst = uAst - dt * alprlx * uhorx
                else if(k == kbl_topo(i, j, 1)) then
                  call wind_ip(var, x_ip(i, j, 1), y_ip(i, j, 1), z_ip(i, j, &
                      1), 'u', u_ip, v_ip, w_ip)

                  u_ip_n = (u_ip * dhdx(i, j, 1) + v_ip * dhdy(i, j, 1) &
                      - w_ip) * dhdx(i, j, 1) / (1 + dhdx(i, j, 1) ** 2 &
                      + dhdy(i, j, 1) ** 2)

                  u_ip_t = u_ip - u_ip_n

                  u_rp_t = velocity_reconst_t(i, j, 1) * u_ip_t
                  u_rp_n = velocity_reconst_n(i, j, 1) * u_ip_n
                  u_rp = u_rp_t + u_rp_n

                  uAst = uAst - dt * alprlx * (var(i, j, k, 2) - u_rp)
                end if
                !UAE
              end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  uAst = uAst - dt * kv_hs(j, k) * uhorx
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                uAst = uAst - dt * kr_sp(j, k) * uhorx
              end if

              usave(i, j, k) = uAst

              if(TestCase == "baroclinic_LC") then
                usave(i, j, k) = usave(i, j, k) + var_env(i, j, k, 2)
              end if
            end do
          end do
        end do
      else if(int_mod == "impl") then
        do k = 1, nz
          do j = 1, ny
            do i = i0, i1
              rhou = 0.5 * (var(i, j, k, 1) + var(i + 1, j, k, 1))

              rhov0m = 0.5 * (var(i, j, k, 1) + var(i, j - 1, k, 1))
              rhov00 = 0.5 * (var(i, j + 1, k, 1) + var(i, j, k, 1))
              rhov1m = 0.5 * (var(i + 1, j, k, 1) + var(i + 1, j - 1, k, 1))
              rhov10 = 0.5 * (var(i + 1, j + 1, k, 1) + var(i + 1, j, k, 1))

              if(fluctuationMode) then
                rhou = rhou + rhoStrat(k)
                rhov0m = rhov0m + rhoStrat(k)
                rhov00 = rhov00 + rhoStrat(k)
                rhov1m = rhov1m + rhoStrat(k)
                rhov10 = rhov10 + rhoStrat(k)
              end if

              if(TestCase == "baroclinic_LC") then
                rhou_e = 0.5 * (var_env(i, j, k, 1) + var_env(i + 1, j, k, 1))

                rhov0m_e = 0.5 * (var_env(i, j, k, 1) + var_env(i, j - 1, k, 1))
                rhov00_e = 0.5 * (var_env(i, j + 1, k, 1) + var_env(i, j, k, 1))
                rhov1m_e = 0.5 * (var_env(i + 1, j, k, 1) + var_env(i + 1, j &
                    - 1, k, 1))
                rhov10_e = 0.5 * (var_env(i + 1, j + 1, k, 1) + var_env(i + 1, &
                    j, k, 1))

                if(fluctuationMode) then
                  rhou_e = rhou_e + rhoStrat_0(k)
                  rhov0m_e = rhov0m_e + rhoStrat_0(k)
                  rhov00_e = rhov00_e + rhoStrat_0(k)
                  rhov1m_e = rhov1m_e + rhoStrat_0(k)
                  rhov10_e = rhov10_e + rhoStrat_0(k)
                end if
              end if

              !--- pressure gradient terms -> piGradx, piGrady
              if(TestCase == "baroclinic_LC") then
                piR = var(i + 1, j, k, 5) - var_env(i + 1, j, k, 5)
                piL = var(i, j, k, 5) - var_env(i, j, k, 5)

                piR_e = var_env(i + 1, j, k, 5)
                piL_e = var_env(i, j, k, 5)
              else
                piR = var(i + 1, j, k, 5)
                piL = var(i, j, k, 5)
              end if

              piGradx = kappaInv * MaInv2 * Pstrat(k) / rhou * (piR - piL) / dx !FS

              if(TestCase == "baroclinic_LC") then !FS
                piGradx = piGradx + kappaInv * MaInv2 * (PStrat(k) / rhou &
                    - pStrat_0(k) / rhou_e) * (piR_e - piL_e) / dx
              end if

              if(TestCase == "baroclinic_LC") then
                piGrady = kappaInv * MaInv2 * 0.25 * (Pstrat(k) / rhov0m &
                    * (var(i, j, k, 5) - var(i, j - 1, k, 5) - var_env(i, j, &
                    k, 5) + var_env(i, j - 1, k, 5)) / dy + Pstrat(k) / rhov00 &
                    * (var(i, j + 1, k, 5) - var(i, j, k, 5) - var_env(i, j &
                    + 1, k, 5) + var_env(i, j, k, 5)) / dy + Pstrat(k) &
                    / rhov1m * (var(i + 1, j, k, 5) - var(i + 1, j - 1, k, 5) &
                    - var_env(i + 1, j, k, 5) + var_env(i + 1, j - 1, k, 5)) &
                    / dy + Pstrat(k) / rhov10 * (var(i + 1, j + 1, k, 5) &
                    - var(i + 1, j, k, 5) - var_env(i + 1, j + 1, k, 5) &
                    + var_env(i + 1, j, k, 5)) / dy)

                piGrady = piGrady + kappaInv * MaInv2 * 0.25 * ((Pstrat(k) &
                    / rhov0m - Pstrat_0(k) / rhov0m_e) * (var_env(i, j, k, 5) &
                    - var_env(i, j - 1, k, 5)) / dy + (Pstrat(k) / rhov00 &
                    - Pstrat_0(k) / rhov00_e) * (var_env(i, j + 1, k, 5) &
                    - var_env(i, j, k, 5)) / dy + (Pstrat(k) / rhov1m &
                    - Pstrat_0(k) / rhov1m_e) * (var_env(i + 1, j, k, 5) &
                    - var_env(i + 1, j - 1, k, 5)) / dy + (Pstrat(k) / rhov10 &
                    - Pstrat_0(k) / rhov10_e) * (var_env(i + 1, j + 1, k, 5) &
                    - var_env(i + 1, j, k, 5)) / dy)
              else
                piGrady = kappaInv * MaInv2 * 0.25 * (Pstrat(k) / rhov0m &
                    * (var(i, j, k, 5) - var(i, j - 1, k, 5)) / dy + Pstrat(k) &
                    / rhov00 * (var(i, j + 1, k, 5) - var(i, j, k, 5)) / dy &
                    + Pstrat(k) / rhov1m * (var(i + 1, j, k, 5) - var(i + 1, j &
                    - 1, k, 5)) / dy + Pstrat(k) / rhov10 * (var(i + 1, j + 1, &
                    k, 5) - var(i + 1, j, k, 5)) / dy)
              end if

              ! gravity-wave forcing
              if(raytracer .or. (testCase == "mountainwave")) then
                volfcx = 0.5 * (force(i, j, k, 1) + force(i + 1, j, k, 1))
                volfcy = 0.5 * (force(i, j, k, 2) + force(i, j + 1, k, 2))
              else
                volfcx = 0.0
                volfcy = 0.0
              end if

              ! ustar
              if(TestCase == "baroclinic_LC") then
                uhorx = var(i, j, k, 2) - var_env(i, j, k, 2)
              else
                uhorx = var(i, j, k, 2)
              end if

              vhory = 0.25 * (var(i, j - 1, k, 3) + var(i, j, k, 3) + var(i &
                  + 1, j - 1, k, 3) + var(i + 1, j, k, 3))

              facu = 1.0

              if(topography) then
                ! Rayleigh damping for topography (immersed boundary)

                !UAC if(topography_mask(i00+i,j00+j,k)&
                !   .or.&
                !   topography_mask(i00+i+1,j00+j,k)) then
                !   facu = facu + dt*alprlx
                !end if
                if(k < kbl_topo(i, j, 1)) then
                  facu = facu + dt * alprlx
                else if(k == kbl_topo(i, j, 1)) then
                  stop 'implementation topography into semi-implicit time step &
                      still to be done'
                end if
                !UAE
              end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  facu = facu + dt * kv_hs(j, k)
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                facu = facu + dt * kr_sp(j, k)
              end if

              facv = facu

              !UAC
              !uAst &
              !=   1.0/(facu*facv + (f_cor_nd(j)*dt)**2) &
              !    * (  facv * (uhorx + dt*(volfcx/rhou - piGradx)) &
              !  + f_cor_nd(j)*dt &
              !    * (vhory + dt*(volfcy/rhou - piGrady)))
              if(testCase == "SkamarockKlemp94") then
                uAst = 1.0 / facu * (uhorx + dt * (f_cor_nd(j) * vhory &
                    + volfcx / rhou - piGradx) + dt ** 2 * f_cor_nd(j) ** 2 &
                    * backgroundFlow_dim(1) / uRef)
              else
                uAst = 1.0 / facu * (uhorx + dt * (f_cor_nd(j) * vhory &
                    + volfcx / rhou - piGradx))
              end if
              !UAE

              usave(i, j, k) = uAst

              if(TestCase == "baroclinic_LC") then
                usave(i, j, k) = usave(i, j, k) + var_env(i, j, k, 2)
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
            fR = flux(i, j, k, 1, 3)
            fL = flux(i - 1, j, k, 1, 3)
            gF = flux(i, j, k, 2, 3)
            gB = flux(i, j - 1, k, 2, 3)
            hU = flux(i, j, k, 3, 3)
            hD = flux(i, j, k - 1, 3, 3)
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

            volForce = 0.

            if(mmp_mod == "tot") then
              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piF = var(i, j + 1, k, 5) - var_env(i, j + 1, k, 5)
                piB = var(i, j, k, 5) - var_env(i, j, k, 5)
              else
                piF = var(i, j + 1, k, 5)
                piB = var(i, j, k, 5)
              end if

              piGrad = kappaInv * MaInv2 * pStrat(k) * (piF - piB) / dy

              if(TestCase == "baroclinic_LC") then !FS
                piGrad = piGrad + kappaInv * MaInv2 * (PStrat(k) &
                    - pStrat_0(k)) * (var_env(i, j + 1, k, 5) - var_env(i, j, &
                    k, 5)) / dy
              end if

              !---- volume forces
              volForce = 0.5 * (force(i, j, k, 2) + force(i, j + 1, k, 2))

              if(TestCase == "baroclinic_LC") then
                if(model == "pseudo_incompressible") then
                  rhoM = rhoOld(i, j, k)
                  rhoM_1 = rhoOld(i, j + 1, k)

                  if(fluctuationMode) then
                    rhoM = rhoM + rhoStrat(k)
                    rhoM_1 = rhoM_1 + rhoStrat(k)
                  end if
                else if(model == "Boussinesq") then
                  rhoM = rho00
                  rhoM_1 = rho00
                else
                  stop "momentumPredictor: unkown model."
                end if

                volforce = volforce + RoInv(j) * (0.25 * (rhoM * (var_env(i &
                    - 1, j, k, 2) + var_env(i, j, k, 2)) + rhoM_1 * (var_env(i &
                    - 1, j + 1, k, 2) + var_env(i, j + 1, k, 2))))
                !  + 0.5*(RoInv(j)+RoInv(j+1)) &

              end if

              if(topography) then
                ! Rayleigh damping for topography (immersed boundary)

                if(model == "pseudo_incompressible") then
                  rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j + 1, k))

                  if(fluctuationMode) then
                    rhoM_1 = rhoM_1 + rhoStrat(k)
                  end if
                else if(model == "Boussinesq") then
                  rhoM_1 = rho00
                else
                  stop "momentumPredictor: unkown model."
                end if

                !UAC if(topography_mask(i00+i,j00+j,k)&
                !   .or.&
                !   topography_mask(i00+i,j00+j+1,k)) then
                !   volForce = volForce - alprlx * rhoM_1*var(i,j,k,3)
                !end if
                if(k < kbl_topo(i, j, 2)) then
                  volForce = volForce - alprlx * rhoM_1 * var(i, j, k, 3)
                else if(k == kbl_topo(i, j, 2)) then
                  call wind_ip(var, x_ip(i, j, 2), y_ip(i, j, 2), z_ip(i, j, &
                      2), 'v', u_ip, v_ip, w_ip)

                  v_ip_n = (u_ip * dhdx(i, j, 2) + v_ip * dhdy(i, j, 2) &
                      - w_ip) * dhdy(i, j, 2) / (1 + dhdx(i, j, 2) ** 2 &
                      + dhdy(i, j, 2) ** 2)

                  v_ip_t = v_ip - v_ip_n

                  v_rp_t = velocity_reconst_t(i, j, 2) * v_ip_t
                  v_rp_n = velocity_reconst_n(i, j, 2) * v_ip_n
                  v_rp = v_rp_t + v_rp_n

                  volForce = volForce - alprlx * rhoM_1 * (var(i, j, k, 3) &
                      - v_rp)
                end if
                !UAE
              end if
            end if

            if(TestCase == "baroclinic_LC") then
              if(background == "HeldSuarez") then
                ! Rayleigh damping

                if(model == "pseudo_incompressible") then
                  rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j + 1, k))

                  if(fluctuationMode) then
                    rhoM_1 = rhoM_1 + rhoStrat(k)
                  end if
                else if(model == "Boussinesq") then
                  rhoM_1 = rho00
                else
                  stop "momentumPredictor: unkown model."
                end if

                volForce = volForce - 0.5 * (kv_hs(j, k) + kv_hs(j + 1, k)) &
                    * rhoM_1 * (var(i, j, k, 3) - var_env(i, j, k, 3))

              end if
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
              F = - fluxDiff + volForce
            else
              stop 'ERROR: wrong mmp_mod'
            end if

            ! interpolated density
            select case(model)

            case("pseudo_incompressible")

              rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j + 1, k))
              rhoM = 0.5 * (var(i, j, k, 1) + var(i, j + 1, k, 1))

              if(fluctuationMode) then
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
            vM_1 = var(i, j, k, 3)
            momM_1 = rhoM_1 * vM_1

            ! q(m-1) -> q(m)
            q(i, j, k, 2) = dt * F + alpha(m) * q(i, j, k, 2)

            ! rhoV(m-1) -> rhoV(m)
            momM = momM_1 + beta(m) * q(i, j, k, 2)

            ! calc v(m,*)
            vAst = momM / rhoM

            ! vAst -> var
            var(i, j, k, 3) = vAst
          end do
        end do
      end do
    else if(mmp_mod == "rhs") then
      if(int_mod == "expl") then
        do k = 1, nz
          do j = j0, j1
            do i = 1, nx
              rhov = 0.5 * (var(i, j, k, 1) + var(i, j + 1, k, 1))
              if(fluctuationMode) then
                rhov = rhov + rhoStrat(k)
              end if

              if(TestCase == "baroclinic_LC") then
                rhov_e = 0.5 * (var_env(i, j, k, 1) + var_env(i, j + 1, k, 1))

                if(fluctuationMode) then
                  rhov_e = rhov_e + rhoStrat_0(k)
                end if
              end if

              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piF = var(i, j + 1, k, 5) - var_env(i, j + 1, k, 5)
                piB = var(i, j, k, 5) - var_env(i, j, k, 5)

                if(TestCase == "baroclinic_LC") then !FS
                  piF_e = var_env(i, j + 1, k, 5)
                  piB_e = var_env(i, j, k, 5)
                end if
              else
                piF = var(i, j + 1, k, 5)
                piB = var(i, j, k, 5)
              end if

              piGrad = kappaInv * MaInv2 * Pstrat(k) / rhov * (piF - piB) / dy

              if(TestCase == "baroclinic_LC") then !FS
                piGrad = piGrad + kappaInv * MaInv2 * (Pstrat(k) / rhov &
                    - Pstrat_0(k) / rhov_e) * (piF_e - piB_e) / dy
              end if

              ! gravity-wave forcing
              if(raytracer .or. (testCase == "mountainwave")) then
                volfcy = 0.5 * (force(i, j, k, 2) + force(i, j + 1, k, 2))
              else
                volfcy = 0.0
              end if

              ! vstar
              if(TestCase == "baroclinic_LC") then
                uhorx = 0.25 * (var(i - 1, j, k, 2) + var(i - 1, j + 1, k, 2) &
                    - var_env(i - 1, j, k, 2) - var_env(i - 1, j + 1, k, 2) &
                    + var(i, j, k, 2) + var(i, j + 1, k, 2) - var_env(i, j, k, &
                    2) - var_env(i, j + 1, k, 2))
              else
                uhorx = 0.25 * (var(i - 1, j, k, 2) + var(i - 1, j + 1, k, 2) &
                    + var(i, j, k, 2) + var(i, j + 1, k, 2))
              end if

              vhory = var(i, j, k, 3)

              f_cor_v = 0.5 * (f_cor_nd(j) + f_cor_nd(j + 1))

              if(testCase == "SkamarockKlemp94") then
                vAst = vhory + dt * (- f_cor_v * (uhorx &
                    - backgroundFlow_dim(1) / uRef) - piGrad + volfcy / rhov)
              else
                vAst = vhory + dt * (- f_cor_v * uhorx - piGrad + volfcy / rhov)
              end if

              if(topography) then
                ! Rayleigh damping for topography (immersed boundary)

                !UAC if(topography_mask(i00+i,j00+j,k)&
                !   .or.&
                !   topography_mask(i00+i,j00+j+1,k)) then
                !   vAst = vAst - dt* alprlx*vhory
                !end if
                if(k < kbl_topo(i, j, 2)) then
                  vAst = vAst - dt * alprlx * vhory
                else if(k == kbl_topo(i, j, 2)) then
                  call wind_ip(var, x_ip(i, j, 2), y_ip(i, j, 2), z_ip(i, j, &
                      2), 'v', u_ip, v_ip, w_ip)

                  v_ip_n = (u_ip * dhdx(i, j, 2) + v_ip * dhdy(i, j, 2) &
                      - w_ip) * dhdy(i, j, 2) / (1 + dhdx(i, j, 2) ** 2 &
                      + dhdy(i, j, 2) ** 2)

                  v_ip_t = v_ip - v_ip_n

                  v_rp_t = velocity_reconst_t(i, j, 2) * v_ip_t
                  v_rp_n = velocity_reconst_n(i, j, 2) * v_ip_n
                  v_rp = v_rp_t + v_rp_n

                  vAst = vAst - dt * alprlx * (var(i, j, k, 3) - v_rp)
                end if
                !UAE
              end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  vAst = vAst - dt * 0.5 * (kv_hs(j, k) + kv_hs(j + 1, k)) &
                      * vhory
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                vAst = vAst - dt * 0.5 * (kr_sp(j, k) + kr_sp(j + 1, k)) * vhory
              end if

              var(i, j, k, 3) = vAst
            end do
          end do
        end do
      else if(int_mod == "impl") then
        do k = 1, nz
          do j = j0, j1
            do i = 1, nx
              rhoum0 = 0.5 * (var(i, j, k, 1) + var(i - 1, j, k, 1))
              rhou00 = 0.5 * (var(i + 1, j, k, 1) + var(i, j, k, 1))
              rhoum1 = 0.5 * (var(i, j + 1, k, 1) + var(i - 1, j + 1, k, 1))
              rhou01 = 0.5 * (var(i + 1, j + 1, k, 1) + var(i, j + 1, k, 1))

              rhov = 0.5 * (var(i, j, k, 1) + var(i, j + 1, k, 1))

              if(fluctuationMode) then
                rhov = rhov + rhoStrat(k)
                rhoum0 = rhoum0 + rhoStrat(k)
                rhou00 = rhou00 + rhoStrat(k)
                rhoum1 = rhoum1 + rhoStrat(k)
                rhou01 = rhou01 + rhoStrat(k)
              end if

              if(TestCase == "baroclinic_LC") then
                rhoum0_e = 0.5 * (var_env(i, j, k, 1) + var_env(i - 1, j, k, 1))
                rhou00_e = 0.5 * (var_env(i + 1, j, k, 1) + var_env(i, j, k, 1))
                rhoum1_e = 0.5 * (var_env(i, j + 1, k, 1) + var_env(i - 1, j &
                    + 1, k, 1))
                rhou01_e = 0.5 * (var_env(i + 1, j + 1, k, 1) + var_env(i, j &
                    + 1, k, 1))

                rhov_e = 0.5 * (var_env(i, j, k, 1) + var_env(i, j + 1, k, 1))

                if(fluctuationMode) then
                  rhov_e = rhov_e + rhoStrat_0(k)
                  rhoum0_e = rhoum0_e + rhoStrat_0(k)
                  rhou00_e = rhou00_e + rhoStrat_0(k)
                  rhoum1_e = rhoum1_e + rhoStrat_0(k)
                  rhou01_e = rhou01_e + rhoStrat_0(k)
                end if
              end if

              !--- pressure gradient terms -> piGradx, piGrady
              if(TestCase == "baroclinic_LC") then
                piGradx = kappaInv * MaInv2 * 0.25 * (Pstrat(k) / rhou00 &
                    * (var(i + 1, j, k, 5) - var(i, j, k, 5) - var_env(i + 1, &
                    j, k, 5) + var_env(i, j, k, 5)) / dx + Pstrat(k) / rhoum0 &
                    * (var(i, j, k, 5) - var(i - 1, j, k, 5) - var_env(i, j, &
                    k, 5) + var_env(i - 1, j, k, 5)) / dx + Pstrat(k) / rhou01 &
                    * (var(i + 1, j + 1, k, 5) - var(i, j + 1, k, 5) &
                    - var_env(i + 1, j + 1, k, 5) + var_env(i, j + 1, k, 5)) &
                    / dx + Pstrat(k) / rhoum1 * (var(i, j + 1, k, 5) - var(i &
                    - 1, j + 1, k, 5) - var_env(i, j + 1, k, 5) + var_env(i &
                    - 1, j + 1, k, 5)) / dx)

                piGradx = piGradx + kappaInv * MaInv2 * 0.25 * ((Pstrat(k) &
                    / rhou00 - Pstrat_0(k) / rhou00_e) * (var_env(i + 1, j, k, &
                    5) - var_env(i, j, k, 5)) / dx + (Pstrat(k) / rhoum0 &
                    - Pstrat_0(k) / rhoum0_e) * (var_env(i, j, k, 5) &
                    - var_env(i - 1, j, k, 5)) / dx + (Pstrat(k) / rhou01 &
                    - Pstrat_0(k) / rhou01_e) * (var_env(i + 1, j + 1, k, 5) &
                    - var_env(i, j + 1, k, 5)) / dx + (Pstrat(k) / rhoum1 &
                    - Pstrat_0(k) / rhoum1_e) * (var_env(i, j + 1, k, 5) &
                    - var_env(i - 1, j + 1, k, 5)) / dx)
              else
                piGradx = kappaInv * MaInv2 * 0.25 * (Pstrat(k) / rhou00 &
                    * (var(i + 1, j, k, 5) - var(i, j, k, 5)) / dx + Pstrat(k) &
                    / rhoum0 * (var(i, j, k, 5) - var(i - 1, j, k, 5)) / dx &
                    + Pstrat(k) / rhou01 * (var(i + 1, j + 1, k, 5) - var(i, j &
                    + 1, k, 5)) / dx + Pstrat(k) / rhoum1 * (var(i, j + 1, k, &
                    5) - var(i - 1, j + 1, k, 5)) / dx)
              end if

              if(TestCase == "baroclinic_LC") then
                piF = var(i, j + 1, k, 5) - var_env(i, j + 1, k, 5)
                piB = var(i, j, k, 5) - var_env(i, j, k, 5)

                piF_e = var_env(i, j + 1, k, 5)
                piB_e = var_env(i, j, k, 5)
              else
                piF = var(i, j + 1, k, 5)
                piB = var(i, j, k, 5)
              end if

              piGrady = kappaInv * MaInv2 * Pstrat(k) / rhov * (piF - piB) / dy

              if(TestCase == "baroclinic_LC") then !FS
                piGrady = piGrady + kappaInv * MaInv2 * (Pstrat(k) / rhov &
                    - Pstrat_0(k) / rhov_e) * (piF_e - piB_e) / dy
              end if

              ! gravity-wave forcing
              if(raytracer .or. (testCase == "mountainwave")) then
                volfcx = 0.5 * (force(i, j, k, 1) + force(i + 1, j, k, 1))
                volfcy = 0.5 * (force(i, j, k, 2) + force(i, j + 1, k, 2))
              else
                volfcx = 0.0
                volfcy = 0.0
              end if

              ! vstar
              if(TestCase == "baroclinic_LC") then
                uhorx = 0.25 * (var(i - 1, j, k, 2) + var(i - 1, j + 1, k, 2) &
                    - var_env(i - 1, j, k, 2) - var_env(i - 1, j + 1, k, 2) &
                    + var(i, j, k, 2) + var(i, j + 1, k, 2) - var_env(i, j, k, &
                    2) - var_env(i, j + 1, k, 2))
              else
                uhorx = 0.25 * (var(i - 1, j, k, 2) + var(i - 1, j + 1, k, 2) &
                    + var(i, j, k, 2) + var(i, j + 1, k, 2))
              end if

              vhory = var(i, j, k, 3)

              facv = 1.0

              if(topography) then
                ! Rayleigh damping for topography (immersed boundary)

                !UAC if(topography_mask(i00+i,j00+j,k)&
                !   .or.&
                !   topography_mask(i00+i,j00+j+1,k)) then
                !   facv = facv + dt*alprlx
                !end if
                if(k < kbl_topo(i, j, 2)) then
                  facv = facv + dt * alprlx
                else if(k == kbl_topo(i, j, 2)) then
                  stop 'implementation topography into semi-implicit time step &
                      still to be done'
                end if
                !UAE
              end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  facv = facv + dt * 0.5 * (kv_hs(j, k) + kv_hs(j + 1, k))
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                facv = facv + dt * 0.5 * (kr_sp(j, k) + kr_sp(j + 1, k))
              end if

              facu = facv

              !UAC
              !vAst &
              != 1.0 &
              !  /(  facu*facv &
              !    + (0.5*(f_cor_nd(j) + f_cor_nd(j+1))*dt)**2) &
              !  * (- 0.5*(f_cor_nd(j) + f_cor_nd(j+1))*dt &
              !       * (uhorx + dt * (volfcx/rhov - piGradx)) &
              !     + facu * (vhory + dt * (volfcy/rhov - piGrady)))
              f_cor_v = 0.5 * (f_cor_nd(j) + f_cor_nd(j + 1))

              if(testCase == "SkamarockKlemp94") then
                vAst = 1.0 / facv * (vhory + dt * (- f_cor_v * (uhorx &
                    - backgroundFlow_dim(1) / uRef) + volfcy / rhov - piGrady))
              else
                vAst = 1.0 / facv * (vhory + dt * (- f_cor_v * uhorx + volfcy &
                    / rhov - piGrady))
              end if
              !UAE

              var(i, j, k, 3) = vAst
            end do
          end do
        end do
      else
        stop 'ERROR: unknown int_mod'
      end if

      ! now the new u can be put into the proper array
      var(:, :, :, 2) = usave(:, :, :)
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
            fR = flux(i, j, k, 1, 4)
            fL = flux(i - 1, j, k, 1, 4)
            gF = flux(i, j, k, 2, 4)
            gB = flux(i, j - 1, k, 2, 4)
            hU = flux(i, j, k, 3, 4)
            hD = flux(i, j, k - 1, 3, 4)
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

            if(mmp_mod == "tot") then
              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piU = var(i, j, k + 1, 5) - var_env(i, j, k + 1, 5)
                piD = var(i, j, k, 5) - var_env(i, j, k, 5)
              else
                piU = var(i, j, k + 1, 5)
                piD = var(i, j, k, 5)
              end if

              piGrad = 0.5 * kappaInv * MaInv2 * (Pstrat(k) + Pstrat(k + 1)) &
                  * (piU - piD) / dz

              if(TestCase == "baroclinic_LC") then !FS
                piGrad = piGrad + 0.5 * kappaInv * MaInv2 * (Pstrat(k) &
                    + Pstrat(k + 1) - pStrat_0(k) - pStrat_0(k + 1)) &
                    * (var_env(i, j, k + 1, 5) - var_env(i, j, k, 5)) / dz
              end if

              !---- volume forces
              volForce = 0.5 * (force(i, j, k, 3) + force(i, j, k + 1, 3))

              if(TestCase == "baroclinic_LC") then
                if(model == "pseudo_incompressible") then
                  drho_e = 0.5 * (var_env(i, j, k, 1) + var_env(i, j, k + 1, 1))

                  if(.not. fluctuationMode) then
                    drho_e = drho_e - rhoStratTilde(k)
                  end if
                else if(model == "Boussinesq") then
                  stop 'ERROR: baroclinic LC not ready yet for  Boussinesq'
                else
                  stop "momentumPredictor: unkown model."
                end if

                volForce = volForce + FrInv2 * drho_e
              end if

              if(topography) then
                ! Rayleigh damping for topography (immersed boundary)

                if(model == "pseudo_incompressible") then
                  rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j, k + 1))

                  if(fluctuationMode) then
                    rhoM_1 = rhoM_1 + rhoStratTilde(k)
                  end if
                else if(model == "Boussinesq") then
                  rhoM_1 = rho00
                else
                  stop "momentumPredictor: unkown model."
                end if

                !UAC if(topography_mask(i00+i,j00+j,k)&
                !   .or.&
                !   topography_mask(i00+i,j00+j,k+1)) then
                !   volForce = volForce - alprlx * rhoM_1*var(i,j,k,4)
                !end if
                if(k < kbl_topo(i, j, 3)) then
                  volForce = volForce - alprlx * rhoM_1 * var(i, j, k, 4)
                else if(k == kbl_topo(i, j, 3)) then
                  call wind_ip(var, x_ip(i, j, 3), y_ip(i, j, 3), z_ip(i, j, &
                      3), 'w', u_ip, v_ip, w_ip)

                  w_ip_n = (- u_ip * dhdx(i, j, 3) - v_ip * dhdy(i, j, 3) &
                      + w_ip) / (1 + dhdx(i, j, 3) ** 2 + dhdy(i, j, 3) ** 2)

                  w_ip_t = w_ip - w_ip_n

                  w_rp_t = velocity_reconst_t(i, j, 3) * w_ip_t
                  w_rp_n = velocity_reconst_n(i, j, 3) * w_ip_n
                  w_rp = w_rp_t + w_rp_n

                  volForce = volForce - alprlx * rhoM_1 * (var(i, j, k, 4) &
                      - w_rp)
                end if
                !UAE
              end if
            end if

            if(TestCase == "baroclinic_LC") then
              if(background == "HeldSuarez") then
                ! Rayleigh damping
                if(model == "pseudo_incompressible") then
                  rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j, k + 1))

                  if(fluctuationMode) then
                    rhoM_1 = rhoM_1 + rhoStratTilde(k)
                  end if
                else if(model == "Boussinesq") then
                  rhoM_1 = rho00
                else
                  stop "momentumPredictor: unkown model."
                end if

                volForce = volforce - 0.5 * (kw_hs(k) + kw_hs(k + 1)) * rhoM_1 &
                    * var(i, j, k, 4)
              end if
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
              F = - fluxDiff
            else
              stop 'ERROR: wrong mmp_mod'
            end if

            ! interpolated densities
            select case(model)

            case("pseudo_incompressible")

              rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j, k + 1)) !rho(m-1)
              rhoM = 0.5 * (var(i, j, k, 1) + var(i, j, k + 1, 1)) !rho(m)

              if(fluctuationMode) then
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
            wM_1 = var(i, j, k, 4)
            momM_1 = rhoM_1 * wM_1

            ! q(m-1) -> q(m)
            q(i, j, k, 3) = dt * F + alpha(m) * q(i, j, k, 3)

            ! rhoW(m-1) -> rhoW(m)
            momM = momM_1 + beta(m) * q(i, j, k, 3)

            ! calc w(m,*)
            wAst = momM / rhoM

            ! wAst -> var
            var(i, j, k, 4) = wAst
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
              rho000 = var(i, j, k, 1)
              rho001 = var(i, j, k + 1, 1)

              rhow = 0.5 * (var(i, j, k, 1) + var(i, j, k + 1, 1))

              if(fluctuationMode) then
                rho000 = rho000 + rhoStrat(k)
                rho001 = rho001 + rhoStrat(k + 1)

                rhow = rhow + rhoStratTilde(k)
              end if

              if(TestCase == "baroclinic_LC") then
                rho000_e = var_env(i, j, k, 1)
                rho001_e = var_env(i, j, k + 1, 1)

                rhow_e = 0.5 * (var_env(i, j, k, 1) + var_env(i, j, k + 1, 1))

                if(fluctuationMode) then
                  rho000_e = rho000_e + rhoStrat_0(k)
                  rho001_e = rho001_e + rhoStrat_0(k + 1)

                  rhow_e = rhow_e + 0.5 * (rhoStrat_0(k) + rhoStrat_0(k + 1))
                end if
              end if

              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piU = var(i, j, k + 1, 5) - var_env(i, j, k + 1, 5)
                piD = var(i, j, k, 5) - var_env(i, j, k, 5)

                piU_e = var_env(i, j, k + 1, 5)
                piD_e = var_env(i, j, k, 5)
              else
                piU = var(i, j, k + 1, 5)
                piD = var(i, j, k, 5)
              end if

              piGrad = kappaInv * MaInv2 * pstw / rhow * (piU - piD) / dz

              if(TestCase == "baroclinic_LC") then !FS
                piGrad = piGrad + kappaInv * MaInv2 * (pstw / rhow - pstw_e &
                    / rhow_e) * (piU_e - piD_e) / dz
              end if

              ! wstar
              wvert = var(i, j, k, 4)

              if(TestCase == "baroclinic_LC") then
                buoy = - g_ndim * 0.5 * (rhopOld(i, j, k) / rho000 &
                    - var_env(i, j, k, 6) / rho000_e + rhopOld(i, j, k + 1) &
                    / rho001 - var_env(i, j, k + 1, 6) / rho001_e)
              else
                buoy = - g_ndim * 0.5 * (rhopOld(i, j, k) / rho000 &
                    + rhopOld(i, j, k + 1) / rho001)
              end if

              wAst = wvert + dt * (buoy - piGrad)

              if(topography) then
                ! Rayleigh damping for topography (immersed boundary)

                !UAC if(topography_mask(i00+i,j00+j,k)&
                !   .or.&
                !   topography_mask(i00+i,j00+j,k+1)) then
                !   wAst = wAst - dt* alprlx*wvert
                !end if
                if(k < kbl_topo(i, j, 3)) then
                  wAst = wAst - dt * alprlx * wvert
                else if(k == kbl_topo(i, j, 3)) then
                  call wind_ip(var, x_ip(i, j, 3), y_ip(i, j, 3), z_ip(i, j, &
                      3), 'w', u_ip, v_ip, w_ip)

                  w_ip_n = (- u_ip * dhdx(i, j, 3) - v_ip * dhdy(i, j, 3) &
                      + w_ip) / (1 + dhdx(i, j, 3) ** 2 + dhdy(i, j, 3) ** 2)

                  w_ip_t = w_ip - w_ip_n

                  w_rp_t = velocity_reconst_t(i, j, 3) * w_ip_t
                  w_rp_n = velocity_reconst_n(i, j, 3) * w_ip_n
                  w_rp = w_rp_t + w_rp_n

                  wAst = wAst - dt * alprlx * (var(i, j, k, 4) - w_rp)
                end if
                !UAE
              end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping

                  wAst = wAst - dt * 0.5 * (kw_hs(k) + kw_hs(k + 1)) * wvert
                end if
              end if

              if(spongeLayer) then
                wAst = wAst - dt * 0.5 * (kr_sp_w(j, k) + kr_sp_w(j, k + 1)) &
                    * wvert
              end if

              var(i, j, k, 4) = wAst
            end do
          end do
        end do
      else if(int_mod == "impl") then
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

        do k = k0, k1
          pstw = 0.5 * (Pstrat(k) + Pstrat(k + 1))
          pstw_0 = 0.5 * (Pstrat_0(k) + Pstrat_0(k + 1))

          if(TestCase == "baroclinic_LC") then
            pstw_e = 0.5 * (Pstrat_0(k) + Pstrat_0(k + 1))
          end if

          do j = 1, ny
            do i = 1, nx
              rho000 = var(i, j, k, 1)
              rho001 = var(i, j, k + 1, 1)

              rhow = 0.5 * (var(i, j, k, 1) + var(i, j, k + 1, 1))

              if(fluctuationMode) then
                rho000 = rho000 + rhoStrat(k)
                rho001 = rho001 + rhoStrat(k + 1)

                rhow = rhow + rhoStratTilde(k)
              end if

              if(TestCase == "baroclinic_LC") then
                rho000_e = var_env(i, j, k, 1)
                rho001_e = var_env(i, j, k + 1, 1)

                rhow_e = 0.5 * (var_env(i, j, k, 1) + var_env(i, j, k + 1, 1))

                if(fluctuationMode) then
                  rho000_e = rho000_e + rhoStrat_0(k)
                  rho001_e = rho001_e + rhoStrat_0(k + 1)

                  rhow_e = rhow_e + 0.5 * (rhoStrat_0(k) + rhoStrat_0(k + 1))
                end if
              end if

              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piU = var(i, j, k + 1, 5) - var_env(i, j, k + 1, 5)
                piD = var(i, j, k, 5) - var_env(i, j, k, 5)

                piU_e = var_env(i, j, k + 1, 5)
                piD_e = var_env(i, j, k, 5)
              else
                piU = var(i, j, k + 1, 5)
                piD = var(i, j, k, 5)
              end if

              piGrad = kappaInv * MaInv2 * pstw / rhow * (piU - piD) / dz

              if(TestCase == "baroclinic_LC") then !FS
                piGrad = piGrad + kappaInv * MaInv2 * (pstw / rhow - pstw_e &
                    / rhow_e) * (piU_e - piD_e) / dz
              end if

              ! wstar
              wvert = var(i, j, k, 4)

              ! squared Brunt-Vaisala frequency averaged to half
              ! levels
              ! (could be done a bit nicer by determining this without
              ! averaging directly from the reference-atmosphere
              ! density)
              bvsstw = 0.5 * (bvsStrat(k) + bvsStrat(k + 1))

              facw = 1.0

              if(topography) then
                ! Rayleigh damping for topography (immersed boundary)

                !UAC if(topography_mask(i00+i,j00+j,k)&
                !   .or.&
                !   topography_mask(i00+i,j00+j,k+1)) then
                !   facw = facw + alprlx*dt
                !end if
                if(k < kbl_topo(i, j, 3)) then
                  facw = facw + alprlx * dt
                else if(k == kbl_topo(i, j, 3)) then
                  stop 'implementation topography into semi-implicit time step &
                      still to be done'
                end if
              end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping

                  facw = facw + dt * 0.5 * (kw_hs(k) + kw_hs(k + 1))
                end if
              end if

              if(spongeLayer) then
                facw = facw + dt * 0.5 * (kr_sp_w(j, k) + kr_sp_w(j, k + 1))
              end if

              !heat0 &
              != heat(i,j,k) - S_bar(k) &
              !  - Pstrat(k)/g_ndim * bvsStrat(k) &
              !    * 0.5*(w_0(k) + w_0(k-1))
              heat0 = heat(i, j, k)

              !heat1 &
              != heat(i,j,k+1) - S_bar(k+1) &
              !  - Pstrat(k+1)/g_ndim * bvsStrat(k+1) &
              !    * 0.5*(w_0(k+1) + w_0(k))
              heat1 = heat(i, j, k + 1)

              if(TestCase == "baroclinic_LC") then
                wAst = 1.0 / (facw + rhoStratTilde(k) / rhow * pstw / pstw_0 &
                    * bvsstw * dt ** 2) * (wvert - dt * piGrad - dt * g_ndim &
                    * 0.5 * (rhopOld(i, j, k) / rho000 - var_env(i, j, k, 6) &
                    / rho000_e + rhopOld(i, j, k + 1) / rho001 - var_env(i, j, &
                    k + 1, 6) / rho001_e + dt * (rhoStrat(k) / Pstrat_0(k) &
                    * heat0 / rho000 + rhoStrat(k + 1) / Pstrat_0(k + 1) &
                    * heat1 / rho001)))
                !/(  facw &
                !  + rhoStratTilde(k)/rhow * bvsstw * dt**2) &
                !* (  rhoStrat(k)/Pstrat(k) &
                !+ rhoStrat(k+1)/Pstrat(k+1) &
              else
                wAst = 1.0 / (facw + rhoStratTilde(k) / rhow * pstw / pstw_0 &
                    * bvsstw * dt ** 2) * (wvert - dt * piGrad - dt * g_ndim &
                    * 0.5 * (rhopOld(i, j, k) / rho000 + rhopOld(i, j, k + 1) &
                    / rho001 + dt * (rhoStrat(k) / Pstrat_0(k) * heat0 &
                    / rho000 + rhoStrat(k + 1) / Pstrat_0(k + 1) * heat1 &
                    / rho001)))
                !/(  facw &
                !  + rhoStratTilde(k)/rhow * bvsstw * dt**2) &
                !* (  rhoStrat(k)/Pstrat(k) &
                !+ rhoStrat(k+1)/Pstrat(k+1) &
              end if

              var(i, j, k, 4) = wAst
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
        topography = topography_s
      else if(int_mod == 'impl') then
        kr_sp = kr_sp / facray
        kr_sp_w = kr_sp_w / facray
        alprlx = alprlx / facray
      end if
    end if

  end subroutine momentumPredictor_nc
  !end subroutine momentumPredictor

  !---------------------------------------------------------------------

  !subroutine momentumPredictor_wc (var,flux,force,dt,q,m,mmp_mod,int_mod, &
  subroutine momentumPredictor(var, flux, force, dt, q, m, mmp_mod, int_mod, &
      facray)
    !----------------------------------
    !  calculates the velocities u^*
    !----------------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var

    real, dimension(- 1:nx, - 1:ny, - 1:nz, 3, nVar), intent(in) :: flux
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
        intent(inout) :: q
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

    ! TFC FJ
    real :: jacEdgeR, jacEdgeF, jacEdgeU
    real :: pEdgeR, pEdgeF, pEdgeU
    real :: piEdgeR, piUEdgeR, piUUEdgeR, piDEdgeR, piDDEdgeR, piEdgeF, &
        piUEdgeF, piUUEdgeF, piDEdgeF, piDDEdgeF, piREdgeU, piLEdgeU, &
        piFEdgeU, piBEdgeU
    real :: rhoStratEdgeR, rhoStratEdgeF, rhoStratEdgeU
    real :: chris11EdgeU, chris22EdgeU, chris13EdgeU, chris23EdgeU
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
        buoy
    real :: rho10, rho01
    real :: rhov0m, rhov00, rhov1m, rhov10
    real :: rhou00, rhoum0, rhou01, rhoum1
    real :: rho000, rho001
    real :: volfcx, volfcy
    real :: bvsstw

    !UAB
    real :: rho_e, rhou_e, rhov_e, rhow_e, pstw_e
    real :: rho10_e, rho01_e
    real :: rhov0m_e, rhov00_e, rhov1m_e, rhov10_e
    real :: rhou00_e, rhoum0_e, rhou01_e, rhoum1_e
    real :: rho000_e, rho001_e
    real :: piR_e, piL_e, piF_e, piB_e, piU_e, piD_e

    real :: f_cor_v

    real :: u_ip, v_ip, w_ip, u_ip_n, v_ip_n, w_ip_n, u_ip_t, v_ip_t, w_ip_t, &
        u_rp, v_rp, w_rp, u_rp_n, v_rp_n, w_rp_n, u_rp_t, v_rp_t, w_rp_t
    !UAE

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
            / (ymax - ymin))
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
      if(int_mod == 'expl' .and. .not. spongeTFC) then
        ! TFC FJ
        spongeLayer_s = spongeLayer
        ! topography_s = topography

        spongeLayer = .false.
        ! topography = .false.
      else if(int_mod == 'impl') then
        kr_sp = kr_sp * facray
        kr_sp_w = kr_sp_w * facray
        alprlx = alprlx * facray
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
            fR = flux(i, j, k, 1, 2)
            fL = flux(i - 1, j, k, 1, 2)
            gF = flux(i, j, k, 2, 2)
            gB = flux(i, j - 1, k, 2, 2)
            hU = flux(i, j, k, 3, 2)
            hD = flux(i, j, k - 1, 3, 2)
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz ! diverg.

            ! TFC FJ
            ! Adjust zonal momentum flux divergence.
            if(topography) then
              jacEdgeR = 0.5 * (jac(i, j, k) + jac(i + 1, j, k))
              fluxDiff = fluxDiff / jacEdgeR
            end if

            volForce = 0.

            if(mmp_mod == "tot") then
              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piR = var(i + 1, j, k, 5) - var_env(i + 1, j, k, 5)
                piL = var(i, j, k, 5) - var_env(i, j, k, 5)
              else
                piR = var(i + 1, j, k, 5)
                piL = var(i, j, k, 5)
              end if

              if(topography) then
                ! TFC FJ
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) - var_env(:, :, :, 5)
                end if
                ! Compute values at cell edges.
                pEdgeR = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    + pStratTFC(i + 1, j, k) / jac(i + 1, j, k))
                ! Compute pressure gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  piUUEdgeR = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 1, 3) &
                      * var(i, j, k + 2, 5) + jac(i + 1, j, k + 2) * met(i &
                      + 1, j, k + 2, 1, 3) * var(i + 1, j, k + 2, 5))
                  piUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
                      * var(i, j, k + 1, 5) + jac(i + 1, j, k + 1) * met(i &
                      + 1, j, k + 1, 1, 3) * var(i + 1, j, k + 1, 5))
                  piEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * var(i, &
                      j, k, 5) + jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) &
                      * var(i + 1, j, k, 5))
                  piGrad = kappaInv * MaInv2 * pEdgeR * ((jac(i + 1, j, k) &
                      * piR - jac(i, j, k) * piL) / dx + (- piUUEdgeR + 4.0 &
                      * piUEdgeR - 3.0 * piEdgeR) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeR = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 1, 3) &
                      * var(i, j, k - 2, 5) + jac(i + 1, j, k - 2) * met(i &
                      + 1, j, k - 2, 1, 3) * var(i + 1, j, k - 2, 5))
                  piDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, 3) &
                      * var(i, j, k - 1, 5) + jac(i + 1, j, k - 1) * met(i &
                      + 1, j, k - 1, 1, 3) * var(i + 1, j, k - 1, 5))
                  piEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * var(i, &
                      j, k, 5) + jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) &
                      * var(i + 1, j, k, 5))
                  piGrad = kappaInv * MaInv2 * pEdgeR * ((jac(i + 1, j, k) &
                      * piR - jac(i, j, k) * piL) / dx + (piDDEdgeR - 4.0 &
                      * piDEdgeR + 3.0 * piEdgeR) * 0.5 / dz)
                else
                  piUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
                      * var(i, j, k + 1, 5) + jac(i + 1, j, k + 1) * met(i &
                      + 1, j, k + 1, 1, 3) * var(i + 1, j, k + 1, 5))
                  piDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, 3) &
                      * var(i, j, k - 1, 5) + jac(i + 1, j, k - 1) * met(i &
                      + 1, j, k - 1, 1, 3) * var(i + 1, j, k - 1, 5))
                  piGrad = kappaInv * MaInv2 * pEdgeR * ((jac(i + 1, j, k) &
                      * piR - jac(i, j, k) * piL) / dx + (piUEdgeR - piDEdgeR) &
                      * 0.5 / dz)
                end if
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) + var_env(:, :, :, 5)
                end if
              else
                piGrad = kappaInv * MaInv2 * Pstrat(k) * (piR - piL) / dx
              end if

              if(TestCase == "baroclinic_LC") then !FS
                piGrad = piGrad + kappaInv * MaInv2 * (PStrat(k) &
                    - pStrat_0(k)) * (var_env(i + 1, j, k, 5) - var_env(i, j, &
                    k, 5)) / dx
              end if

              !---- volume forces
              volForce = 0.5 * (force(i, j, k, 1) + force(i + 1, j, k, 1))

              if(TestCase == "baroclinic_LC") then
                if(model == "pseudo_incompressible") then
                  rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i + 1, j, k))

                  if(fluctuationMode) then
                    if(topography) then
                      ! TFC FJ
                      rhoM_1 = rhoM_1 + 0.5 * (rhoStratTFC(i, j, k) &
                          + rhoStratTFC(i + 1, j, k))
                    else
                      rhoM_1 = rhoM_1 + rhoStrat(k)
                    end if
                  end if
                else if(model == "Boussinesq") then
                  rhoM_1 = rho00
                else
                  stop "momentumPredictor: unkown model."
                end if

                volforce = volforce - rhoM_1 * RoInv(j) * 0.25 * (var_env(i, j &
                    - 1, k, 3) + var_env(i + 1, j - 1, k, 3) + var_env(i, j, &
                    k, 3) + var_env(i + 1, j, k, 3))

              end if

              ! if (topography) then
              !    ! Rayleigh damping for topography (immersed boundary)
              !
              !    if (model == "pseudo_incompressible") then
              !        rhoM_1 = 0.5 * (rhoOld(i,j,k) + rhoOld(i+1,j,k))
              !
              !        if( fluctuationMode ) then
              !           rhoM_1 = rhoM_1 + rhoStrat(k)
              !        end if
              !       else if (model == "Boussinesq") then
              !        rhoM_1 = rho00
              !       else
              !        stop"momentumPredictor: unkown model."
              !    end if
              !
              !    if (k < kbl_topo(i,j,1)) then
              !       volForce = volForce - alprlx * rhoM_1*var(i,j,k,2)
              !      else if (k == kbl_topo(i,j,1)) then
              !       call wind_ip(var, &
              !                  & x_ip(i,j,1),y_ip(i,j,1),z_ip(i,j,1),&
              !                  & 'u',u_ip,v_ip,w_ip)
              !
              !        u_ip_n &
              !        = (u_ip*dhdx(i,j,1) + v_ip*dhdy(i,j,1) - w_ip)&
              !          *dhdx(i,j,1) &
              !          /(1 + dhdx(i,j,1)**2 + dhdy(i,j,1)**2)
              !
              !        u_ip_t = u_ip - u_ip_n
              !
              !        u_rp_t = velocity_reconst_t(i,j,1)*u_ip_t
              !        u_rp_n = velocity_reconst_n(i,j,1)*u_ip_n
              !        u_rp = u_rp_t + u_rp_n
              !
              !        volForce &
              !        = volForce -  alprlx * rhoM_1*(var(i,j,k,2)-u_rp)
              !    end if
              ! end if

            end if

            if(TestCase == "baroclinic_LC") then
              if(background == "HeldSuarez") then
                ! Rayleigh damping

                if(model == "pseudo_incompressible") then
                  rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i + 1, j, k))

                  if(fluctuationMode) then
                    if(topography) then
                      ! TFC FJ
                      rhoM_1 = rhoM_1 + 0.5 * (rhoStratTFC(i, j, k) &
                          + rhoStratTFC(i + 1, j, k))
                    else
                      rhoM_1 = rhoM_1 + rhoStrat(k)
                    end if
                  end if
                else if(model == "Boussinesq") then
                  rhoM_1 = rho00
                else
                  stop "momentumPredictor: unkown model."
                end if

                volForce = volForce - kv_hs(j, k) * rhoM_1 * (var(i, j, k, 2) &
                    - var_env(i, j, k, 2))

              end if
            end if

            ! TFC FJ
            ! Explicit integration of Coriolis force in TFC.
            if(topography .and. mmp_mod == "lhs") then
              uOldTFC(i, j, k) = var(i, j, k, 2)
              vC = 0.5 * (var(i, j, k, 3) + var(i, j - 1, k, 3))
              vR = 0.5 * (var(i + 1, j, k, 3) + var(i + 1, j - 1, k, 3))
              if(testCase == "baroclinic_LC") then
                vC = vC - 0.5 * (var_env(i, j, k, 3) + var_env(i, j - 1, k, 3))
                vR = vR - 0.5 * (var_env(i + 1, j, k, 3) + var_env(i + 1, j &
                    - 1, k, 3))
              end if
              volForce = volForce + 0.5 * f_cor_nd(j) * ((rhoOld(i, j, k) &
                  + rhoStratTFC(i, j, k)) * vC + (rhoOld(i + 1, j, k) &
                  + rhoStratTFC(i + 1, j, k)) * vR)
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

            case("pseudo_incompressible")

              rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i + 1, j, k))
              rhoM = 0.5 * (var(i, j, k, 1) + var(i + 1, j, k, 1))

              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  ! Adjust for 3D fields.
                  rhoStratEdgeR = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i &
                      + 1, j, k))
                  rhoM_1 = rhoM_1 + rhoStratEdgeR
                  rhoM = rhoM + rhoStratEdgeR
                else
                  rhoM_1 = rhoM_1 + rhoStrat(k)
                  rhoM = rhoM + rhoStrat(k)
                end if
              end if

            case("Boussinesq")
              rhoM_1 = rho00
              rhoM = rho00
            case default
              stop "momentumPredictor: unkown case model."
            end select

            ! velocity and momentum at t(m-1)
            uM_1 = var(i, j, k, 2)
            momM_1 = rhoM_1 * uM_1

            ! q(m-1) -> q(m)

            q(i, j, k, 1) = dt * F + alpha(m) * q(i, j, k, 1)

            ! rhoU(m-1) -> rhoU(m)
            momM = momM_1 + beta(m) * q(i, j, k, 1)

            ! calc u(m,*)
            uAst = momM / rhoM

            ! uAst -> var
            var(i, j, k, 2) = uAst
          end do
        end do
      end do
    else if(mmp_mod == "rhs") then
      if(int_mod == "expl") then
        do k = 1, nz
          do j = 1, ny
            do i = i0, i1
              rhou = 0.5 * (var(i, j, k, 1) + var(i + 1, j, k, 1))
              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  rhoStratEdgeR = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i &
                      + 1, j, k))
                  rhou = rhou + rhoStratEdgeR
                else
                  rhou = rhou + rhoStrat(k)
                end if
              end if

              if(TestCase == "baroclinic_LC") then
                rhou_e = 0.5 * (var_env(i, j, k, 1) + var_env(i + 1, j, k, 1))
                if(fluctuationMode) then
                  if(topography) then
                    ! TFC FJ
                    rhou_e = rhou_e + 0.5 * (rhoStratTFC(i, j, k) &
                        + rhoStratTFC(i + 1, j, k))
                  else
                    rhou_e = rhou_e + rhoStrat_0(k)
                  end if
                end if
              end if

              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piR = var(i + 1, j, k, 5) - var_env(i + 1, j, k, 5)
                piL = var(i, j, k, 5) - var_env(i, j, k, 5)

                piR_e = var_env(i + 1, j, k, 5)
                piL_e = var_env(i, j, k, 5)
              else
                piR = var(i + 1, j, k, 5)
                piL = var(i, j, k, 5)
              end if

              if(topography) then
                ! TFC FJ
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) - var_env(:, :, :, 5)
                end if
                ! Compute values at cell edges.
                pEdgeR = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    + pStratTFC(i + 1, j, k) / jac(i + 1, j, k))
                ! Compute pressure gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  piUUEdgeR = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 1, 3) &
                      * var(i, j, k + 2, 5) + jac(i + 1, j, k + 2) * met(i &
                      + 1, j, k + 2, 1, 3) * var(i + 1, j, k + 2, 5))
                  piUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
                      * var(i, j, k + 1, 5) + jac(i + 1, j, k + 1) * met(i &
                      + 1, j, k + 1, 1, 3) * var(i + 1, j, k + 1, 5))
                  piEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * var(i, &
                      j, k, 5) + jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) &
                      * var(i + 1, j, k, 5))
                  piGrad = kappaInv * MaInv2 * pEdgeR / rhou * ((jac(i + 1, j, &
                      k) * piR - jac(i, j, k) * piL) / dx + (- piUUEdgeR + 4.0 &
                      * piUEdgeR - 3.0 * piEdgeR) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeR = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 1, 3) &
                      * var(i, j, k - 2, 5) + jac(i + 1, j, k - 2) * met(i &
                      + 1, j, k - 2, 1, 3) * var(i + 1, j, k - 2, 5))
                  piDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, 3) &
                      * var(i, j, k - 1, 5) + jac(i + 1, j, k - 1) * met(i &
                      + 1, j, k - 1, 1, 3) * var(i + 1, j, k - 1, 5))
                  piEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * var(i, &
                      j, k, 5) + jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) &
                      * var(i + 1, j, k, 5))
                  piGrad = kappaInv * MaInv2 * pEdgeR / rhou * ((jac(i + 1, j, &
                      k) * piR - jac(i, j, k) * piL) / dx + (piDDEdgeR - 4.0 &
                      * piDEdgeR + 3.0 * piEdgeR) * 0.5 / dz)
                else
                  piUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
                      * var(i, j, k + 1, 5) + jac(i + 1, j, k + 1) * met(i &
                      + 1, j, k + 1, 1, 3) * var(i + 1, j, k + 1, 5))
                  piDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, 3) &
                      * var(i, j, k - 1, 5) + jac(i + 1, j, k - 1) * met(i &
                      + 1, j, k - 1, 1, 3) * var(i + 1, j, k - 1, 5))
                  piGrad = kappaInv * MaInv2 * pEdgeR / rhou * ((jac(i + 1, j, &
                      k) * piR - jac(i, j, k) * piL) / dx + (piUEdgeR &
                      - piDEdgeR) * 0.5 / dz)
                end if
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) + var_env(:, :, :, 5)
                end if
              else
                piGrad = kappaInv * MaInv2 * Pstrat(k) / rhou * (piR - piL) / dx
              end if

              if(TestCase == "baroclinic_LC") then !FS
                if(topography) then
                  ! TFC FJ
                  ! Compute values at cell edges.
                  pEdgeR = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                      + pStratTFC(i + 1, j, k) / jac(i + 1, j, k))
                  ! Compute pressure gradient component.
                  if(k == 1 .and. zBoundary == "solid_wall") then
                    piUUEdgeR = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 1, &
                        3) * var_env(i, j, k + 2, 5) + jac(i + 1, j, k + 2) &
                        * met(i + 1, j, k + 2, 1, 3) * var_env(i + 1, j, k &
                        + 2, 5))
                    piUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, &
                        3) * var_env(i, j, k + 1, 5) + jac(i + 1, j, k + 1) &
                        * met(i + 1, j, k + 1, 1, 3) * var_env(i + 1, j, k &
                        + 1, 5))
                    piEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) &
                        * var_env(i, j, k, 5) + jac(i + 1, j, k) * met(i + 1, &
                        j, k, 1, 3) * var_env(i + 1, j, k, 5))
                    piGrad = piGrad + kappaInv * MaInv2 * (pEdgeR / rhou &
                        - pEdgeR / rhou_e) * ((jac(i + 1, j, k) * piR_e &
                        - jac(i, j, k) * piL_e) / dx + (- piUUEdgeR + 4.0 &
                        * piUEdgeR - 3.0 * piEdgeR) * 0.5 / dz)
                  else if(k == nz .and. zBoundary == "solid_wall") then
                    piDDEdgeR = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 1, &
                        3) * var_env(i, j, k - 2, 5) + jac(i + 1, j, k - 2) &
                        * met(i + 1, j, k - 2, 1, 3) * var_env(i + 1, j, k &
                        - 2, 5))
                    piDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, &
                        3) * var_env(i, j, k - 1, 5) + jac(i + 1, j, k - 1) &
                        * met(i + 1, j, k - 1, 1, 3) * var_env(i + 1, j, k &
                        - 1, 5))
                    piEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) &
                        * var_env(i, j, k, 5) + jac(i + 1, j, k) * met(i + 1, &
                        j, k, 1, 3) * var_env(i + 1, j, k, 5))
                    piGrad = piGrad + kappaInv * MaInv2 * (pEdgeR / rhou &
                        - pEdgeR / rhou_e) * ((jac(i + 1, j, k) * piR_e &
                        - jac(i, j, k) * piL_e) / dx + (piDDEdgeR - 4.0 &
                        * piDEdgeR + 3.0 * piEdgeR) * 0.5 / dz)
                  else
                    piUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, &
                        3) * var_env(i, j, k + 1, 5) + jac(i + 1, j, k + 1) &
                        * met(i + 1, j, k + 1, 1, 3) * var_env(i + 1, j, k &
                        + 1, 5))
                    piDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, &
                        3) * var_env(i, j, k - 1, 5) + jac(i + 1, j, k - 1) &
                        * met(i + 1, j, k - 1, 1, 3) * var_env(i + 1, j, k &
                        - 1, 5))
                    piGrad = piGrad + kappaInv * MaInv2 * (pEdgeR / rhou &
                        - pEdgeR / rhou_e) * ((jac(i + 1, j, k) * piR_e &
                        - jac(i, j, k) * piL_e) / dx + (piUEdgeR - piDEdgeR) &
                        * 0.5 / dz)
                  end if
                else
                  piGrad = piGrad + kappaInv * MaInv2 * (Pstrat(k) / rhou &
                      - Pstrat_0(k) / rhou_e) * (piR_e - piL_e) / dx
                end if
              end if

              ! gravity-wave forcing
              if(raytracer .or. (testCase == "mountainwave")) then
                volfcx = 0.5 * (force(i, j, k, 1) + force(i + 1, j, k, 1))
              else
                volfcx = 0.0
              end if

              ! ustar
              if(TestCase == "baroclinic_LC") then
                uhorx = var(i, j, k, 2) - var_env(i, j, k, 2)
              else
                uhorx = var(i, j, k, 2)
              end if

              if(topography) then
                ! TFC FJ
                ! Coriolis force is integrated on LHS.
                uAst = uhorx + dt * (- piGrad + volfcx / rhou)
              else
                vhory = 0.25 * (var(i, j - 1, k, 3) + var(i, j, k, 3) + var(i &
                    + 1, j - 1, k, 3) + var(i + 1, j, k, 3))

                uAst = uhorx + dt * (f_cor_nd(j) * vhory - piGrad + volfcx &
                    / rhou)
              end if

              ! if (topography) then
              !    ! Rayleigh damping for topography (immersed boundary)
              !
              !    if (TestCase == "baroclinic_LC") then
              !       stop'combination of topography with baroclinic &
              !          & LC not possible yet'
              !    end if
              !
              !    if (k < kbl_topo(i,j,1)) then
              !       uAst = uAst - dt* alprlx*uhorx
              !      else if (k == kbl_topo(i,j,1)) then
              !       call wind_ip(var, &
              !                  & x_ip(i,j,1),y_ip(i,j,1),z_ip(i,j,1),&
              !                  & 'u',u_ip,v_ip,w_ip)
              !
              !        u_ip_n &
              !        = (u_ip*dhdx(i,j,1) + v_ip*dhdy(i,j,1) - w_ip)&
              !          *dhdx(i,j,1) &
              !          /(1 + dhdx(i,j,1)**2 + dhdy(i,j,1)**2)
              !
              !        u_ip_t = u_ip - u_ip_n
              !
              !        u_rp_t = velocity_reconst_t(i,j,1)*u_ip_t
              !        u_rp_n = velocity_reconst_n(i,j,1)*u_ip_n
              !        u_rp = u_rp_t + u_rp_n
              !
              !        uAst &
              !        = uAst - dt*alprlx * (var(i,j,k,2)-u_rp)
              !    end if
              ! end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  uAst = uAst - dt * kv_hs(j, k) * uhorx
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                uAst = uAst - dt * kr_sp(j, k) * uhorx
              end if

              usave(i, j, k) = uAst

              if(TestCase == "baroclinic_LC") then
                usave(i, j, k) = usave(i, j, k) + var_env(i, j, k, 2)
              end if
            end do
          end do
        end do
      else if(int_mod == "impl") then
        do k = 1, nz
          do j = 1, ny
            do i = i0, i1
              rhou = 0.5 * (var(i, j, k, 1) + var(i + 1, j, k, 1))

              rhov0m = 0.5 * (var(i, j, k, 1) + var(i, j - 1, k, 1))
              rhov00 = 0.5 * (var(i, j + 1, k, 1) + var(i, j, k, 1))
              rhov1m = 0.5 * (var(i + 1, j, k, 1) + var(i + 1, j - 1, k, 1))
              rhov10 = 0.5 * (var(i + 1, j + 1, k, 1) + var(i + 1, j, k, 1))

              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  rhoStratEdgeR = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i &
                      + 1, j, k))
                  rhou = rhou + rhoStratEdgeR
                else
                  rhou = rhou + rhoStrat(k)
                  rhov0m = rhov0m + rhoStrat(k)
                  rhov00 = rhov00 + rhoStrat(k)
                  rhov1m = rhov1m + rhoStrat(k)
                  rhov10 = rhov10 + rhoStrat(k)
                end if
              end if

              if(TestCase == "baroclinic_LC") then
                rhou_e = 0.5 * (var_env(i, j, k, 1) + var_env(i + 1, j, k, 1))

                rhov0m_e = 0.5 * (var_env(i, j, k, 1) + var_env(i, j - 1, k, 1))
                rhov00_e = 0.5 * (var_env(i, j + 1, k, 1) + var_env(i, j, k, 1))
                rhov1m_e = 0.5 * (var_env(i + 1, j, k, 1) + var_env(i + 1, j &
                    - 1, k, 1))
                rhov10_e = 0.5 * (var_env(i + 1, j + 1, k, 1) + var_env(i + 1, &
                    j, k, 1))

                if(fluctuationMode) then
                  if(topography) then
                    ! TFC FJ
                    rhou_e = rhou_e + 0.5 * (rhoStratTFC(i, j, k) &
                        + rhoStratTFC(i + 1, j, k))
                  else
                    rhou_e = rhou_e + rhoStrat_0(k)
                    rhov0m_e = rhov0m_e + rhoStrat_0(k)
                    rhov00_e = rhov00_e + rhoStrat_0(k)
                    rhov1m_e = rhov1m_e + rhoStrat_0(k)
                    rhov10_e = rhov10_e + rhoStrat_0(k)
                  end if
                end if
              end if

              !--- pressure gradient terms -> piGradx, piGrady
              if(TestCase == "baroclinic_LC") then
                piR = var(i + 1, j, k, 5) - var_env(i + 1, j, k, 5)
                piL = var(i, j, k, 5) - var_env(i, j, k, 5)

                piR_e = var_env(i + 1, j, k, 5)
                piL_e = var_env(i, j, k, 5)
              else
                piR = var(i + 1, j, k, 5)
                piL = var(i, j, k, 5)
              end if

              if(topography) then
                ! TFC FJ
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) - var_env(:, :, :, 5)
                end if
                ! Compute values at cell edges.
                pEdgeR = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    + pStratTFC(i + 1, j, k) / jac(i + 1, j, k))
                ! Compute pressure gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  piUUEdgeR = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 1, 3) &
                      * var(i, j, k + 2, 5) + jac(i + 1, j, k + 2) * met(i &
                      + 1, j, k + 2, 1, 3) * var(i + 1, j, k + 2, 5))
                  piUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
                      * var(i, j, k + 1, 5) + jac(i + 1, j, k + 1) * met(i &
                      + 1, j, k + 1, 1, 3) * var(i + 1, j, k + 1, 5))
                  piEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * var(i, &
                      j, k, 5) + jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) &
                      * var(i + 1, j, k, 5))
                  piGradX = kappaInv * MaInv2 * pEdgeR / rhou * ((jac(i + 1, &
                      j, k) * piR - jac(i, j, k) * piL) / dx + (- piUUEdgeR &
                      + 4.0 * piUEdgeR - 3.0 * piEdgeR) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeR = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 1, 3) &
                      * var(i, j, k - 2, 5) + jac(i + 1, j, k - 2) * met(i &
                      + 1, j, k - 2, 1, 3) * var(i + 1, j, k - 2, 5))
                  piDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, 3) &
                      * var(i, j, k - 1, 5) + jac(i + 1, j, k - 1) * met(i &
                      + 1, j, k - 1, 1, 3) * var(i + 1, j, k - 1, 5))
                  piEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * var(i, &
                      j, k, 5) + jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) &
                      * var(i + 1, j, k, 5))
                  piGradX = kappaInv * MaInv2 * pEdgeR / rhou * ((jac(i + 1, &
                      j, k) * piR - jac(i, j, k) * piL) / dx + (piDDEdgeR &
                      - 4.0 * piDEdgeR + 3.0 * piEdgeR) * 0.5 / dz)
                else
                  piUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
                      * var(i, j, k + 1, 5) + jac(i + 1, j, k + 1) * met(i &
                      + 1, j, k + 1, 1, 3) * var(i + 1, j, k + 1, 5))
                  piDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, 3) &
                      * var(i, j, k - 1, 5) + jac(i + 1, j, k - 1) * met(i &
                      + 1, j, k - 1, 1, 3) * var(i + 1, j, k - 1, 5))
                  piGradX = kappaInv * MaInv2 * pEdgeR / rhou * ((jac(i + 1, &
                      j, k) * piR - jac(i, j, k) * piL) / dx + (piUEdgeR &
                      - piDEdgeR) * 0.5 / dz)
                end if
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) + var_env(:, :, :, 5)
                end if
              else
                piGradx = kappaInv * MaInv2 * Pstrat(k) / rhou * (piR - piL) &
                    / dx !FS
              end if

              if(TestCase == "baroclinic_LC") then !FS
                if(topography) then
                  ! TFC FJ
                  ! Compute values at cell edges.
                  pEdgeR = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                      + pStratTFC(i + 1, j, k) / jac(i + 1, j, k))
                  ! Compute pressure gradient component.
                  if(k == 1 .and. zBoundary == "solid_wall") then
                    piUUEdgeR = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 1, &
                        3) * var_env(i, j, k + 2, 5) + jac(i + 1, j, k + 2) &
                        * met(i + 1, j, k + 2, 1, 3) * var_env(i + 1, j, k &
                        + 2, 5))
                    piUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, &
                        3) * var_env(i, j, k + 1, 5) + jac(i + 1, j, k + 1) &
                        * met(i + 1, j, k + 1, 1, 3) * var_env(i + 1, j, k &
                        + 1, 5))
                    piEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) &
                        * var_env(i, j, k, 5) + jac(i + 1, j, k) * met(i + 1, &
                        j, k, 1, 3) * var_env(i + 1, j, k, 5))
                    piGradX = piGradX + kappaInv * MaInv2 * (pEdgeR / rhou &
                        - pEdgeR / rhou_e) * ((jac(i + 1, j, k) * piR_e &
                        - jac(i, j, k) * piL_e) / dx + (- piUUEdgeR + 4.0 &
                        * piUEdgeR - 3.0 * piEdgeR) * 0.5 / dz)
                  else if(k == nz .and. zBoundary == "solid_wall") then
                    piDDEdgeR = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 1, &
                        3) * var_env(i, j, k - 2, 5) + jac(i + 1, j, k - 2) &
                        * met(i + 1, j, k - 2, 1, 3) * var_env(i + 1, j, k &
                        - 2, 5))
                    piDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, &
                        3) * var_env(i, j, k - 1, 5) + jac(i + 1, j, k - 1) &
                        * met(i + 1, j, k - 1, 1, 3) * var_env(i + 1, j, k &
                        - 1, 5))
                    piEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) &
                        * var_env(i, j, k, 5) + jac(i + 1, j, k) * met(i + 1, &
                        j, k, 1, 3) * var_env(i + 1, j, k, 5))
                    piGradX = piGradX + kappaInv * MaInv2 * (pEdgeR / rhou &
                        - pEdgeR / rhou_e) * ((jac(i + 1, j, k) * piR_e &
                        - jac(i, j, k) * piL_e) / dx + (piDDEdgeR - 4.0 &
                        * piDEdgeR + 3.0 * piEdgeR) * 0.5 / dz)
                  else
                    piUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, &
                        3) * var_env(i, j, k + 1, 5) + jac(i + 1, j, k + 1) &
                        * met(i + 1, j, k + 1, 1, 3) * var_env(i + 1, j, k &
                        + 1, 5))
                    piDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, &
                        3) * var_env(i, j, k - 1, 5) + jac(i + 1, j, k - 1) &
                        * met(i + 1, j, k - 1, 1, 3) * var_env(i + 1, j, k &
                        - 1, 5))
                    piGradX = piGradX + kappaInv * MaInv2 * (pEdgeR / rhou &
                        - pEdgeR / rhou_e) * ((jac(i + 1, j, k) * piR_e &
                        - jac(i, j, k) * piL_e) / dx + (piUEdgeR - piDEdgeR) &
                        * 0.5 / dz)
                  end if
                else
                  piGradx = piGradx + kappaInv * MaInv2 * (PStrat(k) / rhou &
                      - pStrat_0(k) / rhou_e) * (piR_e - piL_e) / dx
                end if
              end if

              if(TestCase == "baroclinic_LC") then
                piGrady = kappaInv * MaInv2 * 0.25 * (Pstrat(k) / rhov0m &
                    * (var(i, j, k, 5) - var(i, j - 1, k, 5) - var_env(i, j, &
                    k, 5) + var_env(i, j - 1, k, 5)) / dy + Pstrat(k) / rhov00 &
                    * (var(i, j + 1, k, 5) - var(i, j, k, 5) - var_env(i, j &
                    + 1, k, 5) + var_env(i, j, k, 5)) / dy + Pstrat(k) &
                    / rhov1m * (var(i + 1, j, k, 5) - var(i + 1, j - 1, k, 5) &
                    - var_env(i + 1, j, k, 5) + var_env(i + 1, j - 1, k, 5)) &
                    / dy + Pstrat(k) / rhov10 * (var(i + 1, j + 1, k, 5) &
                    - var(i + 1, j, k, 5) - var_env(i + 1, j + 1, k, 5) &
                    + var_env(i + 1, j, k, 5)) / dy)

                piGrady = piGrady + kappaInv * MaInv2 * 0.25 * ((Pstrat(k) &
                    / rhov0m - Pstrat_0(k) / rhov0m_e) * (var_env(i, j, k, 5) &
                    - var_env(i, j - 1, k, 5)) / dy + (Pstrat(k) / rhov00 &
                    - Pstrat_0(k) / rhov00_e) * (var_env(i, j + 1, k, 5) &
                    - var_env(i, j, k, 5)) / dy + (Pstrat(k) / rhov1m &
                    - Pstrat_0(k) / rhov1m_e) * (var_env(i + 1, j, k, 5) &
                    - var_env(i + 1, j - 1, k, 5)) / dy + (Pstrat(k) / rhov10 &
                    - Pstrat_0(k) / rhov10_e) * (var_env(i + 1, j + 1, k, 5) &
                    - var_env(i + 1, j, k, 5)) / dy)
              else
                piGrady = kappaInv * MaInv2 * 0.25 * (Pstrat(k) / rhov0m &
                    * (var(i, j, k, 5) - var(i, j - 1, k, 5)) / dy + Pstrat(k) &
                    / rhov00 * (var(i, j + 1, k, 5) - var(i, j, k, 5)) / dy &
                    + Pstrat(k) / rhov1m * (var(i + 1, j, k, 5) - var(i + 1, j &
                    - 1, k, 5)) / dy + Pstrat(k) / rhov10 * (var(i + 1, j + 1, &
                    k, 5) - var(i + 1, j, k, 5)) / dy)
              end if

              ! gravity-wave forcing
              if(raytracer .or. (testCase == "mountainwave")) then
                volfcx = 0.5 * (force(i, j, k, 1) + force(i + 1, j, k, 1))
                volfcy = 0.5 * (force(i, j, k, 2) + force(i, j + 1, k, 2))
              else
                volfcx = 0.0
                volfcy = 0.0
              end if

              ! ustar
              if(TestCase == "baroclinic_LC") then
                uhorx = var(i, j, k, 2) - var_env(i, j, k, 2)
              else
                uhorx = var(i, j, k, 2)
              end if

              vhory = 0.25 * (var(i, j - 1, k, 3) + var(i, j, k, 3) + var(i &
                  + 1, j - 1, k, 3) + var(i + 1, j, k, 3))

              facu = 1.0

              ! if (topography) then
              !    ! Rayleigh damping for topography (immersed boundary)
              !
              !    if(k < kbl_topo(i,j,1)) then
              !       facu = facu + dt*alprlx
              !      else if(k == kbl_topo(i,j,1)) then
              !        stop'implementation topography into &
              !            &semi-implicit time step still to be done'
              !    end if
              ! end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  facu = facu + dt * kv_hs(j, k)
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                facu = facu + dt * kr_sp(j, k)
              end if

              facv = facu

              if(topography) then
                ! TFC FJ
                ! Coriolis force is integrated on LHS.
                uAst = 1.0 / facu * (uhorx + dt * (- piGradX + volfcx / rhou))
              else
                if(testCase == "SkamarockKlemp94") then
                  uAst = 1.0 / (facu * facv + (f_cor_nd(j) * dt) ** 2) * (facv &
                      * (uhorx + dt * (volfcx / rhou - piGradx)) + f_cor_nd(j) &
                      * dt * (vhory + dt * (volfcy / rhou - piGrady)) &
                      + f_cor_nd(j) ** 2 * dt ** 2 * backgroundFlow_dim(1) &
                      / uRef)
                else
                  uAst = 1.0 / (facu * facv + (f_cor_nd(j) * dt) ** 2) * (facv &
                      * (uhorx + dt * (volfcx / rhou - piGradx)) + f_cor_nd(j) &
                      * dt * (vhory + dt * (volfcy / rhou - piGrady)))
                end if
              end if

              usave(i, j, k) = uAst

              if(TestCase == "baroclinic_LC") then
                usave(i, j, k) = usave(i, j, k) + var_env(i, j, k, 2)
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
            fR = flux(i, j, k, 1, 3)
            fL = flux(i - 1, j, k, 1, 3)
            gF = flux(i, j, k, 2, 3)
            gB = flux(i, j - 1, k, 2, 3)
            hU = flux(i, j, k, 3, 3)
            hD = flux(i, j, k - 1, 3, 3)
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

            ! TFC FJ
            ! Adjust meridional momentum flux divergence.
            if(topography) then
              jacEdgeF = 0.5 * (jac(i, j, k) + jac(i, j + 1, k))
              fluxDiff = fluxDiff / jacEdgeF
            end if

            volForce = 0.

            if(mmp_mod == "tot") then
              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piF = var(i, j + 1, k, 5) - var_env(i, j + 1, k, 5)
                piB = var(i, j, k, 5) - var_env(i, j, k, 5)
              else
                piF = var(i, j + 1, k, 5)
                piB = var(i, j, k, 5)
              end if

              if(topography) then
                ! TFC FJ
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) - var_env(:, :, :, 5)
                end if
                ! Compute values at cell edges.
                pEdgeF = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    + pStratTFC(i, j + 1, k) / jac(i, j + 1, k))
                ! Compute pressure gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  piUUEdgeF = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 2, 3) &
                      * var(i, j, k + 2, 5) + jac(i, j + 1, k + 2) * met(i, j &
                      + 1, k + 2, 2, 3) * var(i, j + 1, k + 2, 5))
                  piUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                      * var(i, j, k + 1, 5) + jac(i, j + 1, k + 1) * met(i, j &
                      + 1, k + 1, 2, 3) * var(i, j + 1, k + 1, 5))
                  piEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) * var(i, &
                      j, k, 5) + jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) &
                      * var(i, j + 1, k, 5))
                  piGrad = kappaInv * MaInv2 * pEdgeF * ((jac(i, j + 1, k) &
                      * piF - jac(i, j, k) * piB) / dy + (- piUUEdgeF + 4.0 &
                      * piUEdgeF - 3.0 * piEdgeF) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeF = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 2, 3) &
                      * var(i, j, k - 2, 5) + jac(i, j + 1, k - 2) * met(i, j &
                      + 1, k - 2, 2, 3) * var(i, j + 1, k - 2, 5))
                  piDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) &
                      * var(i, j, k - 1, 5) + jac(i, j + 1, k - 1) * met(i, j &
                      + 1, k - 1, 2, 3) * var(i, j + 1, k - 1, 5))
                  piEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) * var(i, &
                      j, k, 5) + jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) &
                      * var(i, j + 1, k, 5))
                  piGrad = kappaInv * MaInv2 * pEdgeF * ((jac(i, j + 1, k) &
                      * piF - jac(i, j, k) * piB) / dy + (piDDEdgeF - 4.0 &
                      * piDEdgeF + 3.0 * piEdgeF) * 0.5 / dz)
                else
                  piUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                      * var(i, j, k + 1, 5) + jac(i, j + 1, k + 1) * met(i, j &
                      + 1, k + 1, 2, 3) * var(i, j + 1, k + 1, 5))
                  piDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) &
                      * var(i, j, k - 1, 5) + jac(i, j + 1, k - 1) * met(i, j &
                      + 1, k - 1, 2, 3) * var(i, j + 1, k - 1, 5))
                  piGrad = kappaInv * MaInv2 * pEdgeF * ((jac(i, j + 1, k) &
                      * piF - jac(i, j, k) * piB) / dy + (piUEdgeF - piDEdgeF) &
                      * 0.5 / dz)
                end if
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) + var_env(:, :, :, 5)
                end if
              else
                piGrad = kappaInv * MaInv2 * pStrat(k) * (piF - piB) / dy
              end if

              if(TestCase == "baroclinic_LC") then !FS
                piGrad = piGrad + kappaInv * MaInv2 * (PStrat(k) &
                    - pStrat_0(k)) * (var_env(i, j + 1, k, 5) - var_env(i, j, &
                    k, 5)) / dy
              end if

              !---- volume forces
              volForce = 0.5 * (force(i, j, k, 2) + force(i, j + 1, k, 2))

              if(TestCase == "baroclinic_LC") then
                if(model == "pseudo_incompressible") then
                  rhoM = rhoOld(i, j, k)
                  rhoM_1 = rhoOld(i, j + 1, k)

                  if(fluctuationMode) then
                    if(topography) then
                      ! TFC FJ
                      rhoStratEdgeF = 0.5 * (rhoStratTFC(i, j, k) &
                          + rhoStratTFC(i, j + 1, k))
                      rhoM = rhoM + rhoStratEdgeF
                      rhoM_1 = rhoM_1 + rhoStratEdgeF
                    else
                      rhoM = rhoM + rhoStrat(k)
                      rhoM_1 = rhoM_1 + rhoStrat(k)
                    end if
                  end if
                else if(model == "Boussinesq") then
                  rhoM = rho00
                  rhoM_1 = rho00
                else
                  stop "momentumPredictor: unkown model."
                end if

                volforce = volforce + RoInv(j) * 0.25 * (rhoM * (var_env(i &
                    - 1, j, k, 2) + var_env(i, j, k, 2)) + rhoM_1 * (var_env(i &
                    - 1, j + 1, k, 2) + var_env(i, j + 1, k, 2)))

              end if

              !    if (topography) then
              !       ! Rayleigh damping for topography (immersed boundary)
              !
              !       if (model == "pseudo_incompressible") then
              !           rhoM_1 = 0.5 * (rhoOld(i,j,k) + rhoOld(i,j+1,k))
              !
              !           if( fluctuationMode ) then
              !              rhoM_1 = rhoM_1 + rhoStrat(k)
              !           end if
              !          else if (model == "Boussinesq") then
              !           rhoM_1 = rho00
              !          else
              !           stop"momentumPredictor: unkown model."
              !       end if
              !
              !       if(k < kbl_topo(i,j,2)) then
              !          volForce = volForce - alprlx * rhoM_1*var(i,j,k,3)
              !         else if(k == kbl_topo(i,j,2)) then
              !          call wind_ip(var, &
              !                     & x_ip(i,j,2),y_ip(i,j,2),z_ip(i,j,2),&
              !                     &'v',u_ip,v_ip,w_ip)
              !
              !           v_ip_n &
              !           = (u_ip*dhdx(i,j,2) + v_ip*dhdy(i,j,2) - w_ip)&
              !             *dhdy(i,j,2) &
              !             /(1 + dhdx(i,j,2)**2 + dhdy(i,j,2)**2)
              !
              !           v_ip_t = v_ip-v_ip_n
              !
              !           v_rp_t = velocity_reconst_t(i,j,2)*v_ip_t
              !           v_rp_n = velocity_reconst_n(i,j,2)*v_ip_n
              !           v_rp = v_rp_t + v_rp_n
              !
              !           volForce &
              !           = volForce - alprlx * rhoM_1*(var(i,j,k,3)-v_rp)
              !       end if
              !    end if

            end if

            if(TestCase == "baroclinic_LC") then
              if(background == "HeldSuarez") then
                ! Rayleigh damping

                if(model == "pseudo_incompressible") then
                  rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j + 1, k))

                  if(fluctuationMode) then
                    if(topography) then
                      ! TFC FJ
                      rhoM_1 = rhoM_1 + 0.5 * (rhoStratTFC(i, j, k) &
                          + rhoStratTFC(i, j + 1, k))
                    else
                      rhoM_1 = rhoM_1 + rhoStrat(k)
                    end if
                  end if
                else if(model == "Boussinesq") then
                  rhoM_1 = rho00
                else
                  stop "momentumPredictor: unkown model."
                end if

                volForce = volForce - 0.5 * (kv_hs(j, k) + kv_hs(j + 1, k)) &
                    * rhoM_1 * (var(i, j, k, 3) - var_env(i, j, k, 3))

              end if
            end if

            ! TFC FJ
            ! Explicit integration of Coriolis force in TFC.
            if(topography .and. mmp_mod == "lhs") then
              vOldTFC(i, j, k) = var(i, j, k, 3)
              uC = 0.5 * (uOldTFC(i, j, k) + uOldTFC(i - 1, j, k))
              uF = 0.5 * (uOldTFC(i, j + 1, k) + uOldTFC(i - 1, j + 1, k))
              if(testCase == "SkamarockKlemp94") then
                uC = uC - backgroundFlow_dim(1) / uRef
                uF = uF - backgroundFlow_dim(1) / uRef
              else if(testCase == "baroclinic_LC") then
                uC = uC - 0.5 * (var_env(i, j, k, 2) + var_env(i - 1, j, k, 2))
                uF = uF - 0.5 * (var_env(i, j + 1, k, 2) + var_env(i - 1, j &
                    + 1, k, 2))
              end if
              volForce = volForce - 0.5 * (f_cor_nd(j) * (rhoOld(i, j, k) &
                  + rhoStratTFC(i, j, k)) * uC + f_cor_nd(j + 1) * (rhoOld(i, &
                  j + 1, k) + rhoStratTFC(i, j + 1, k)) * uF)
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

            case("pseudo_incompressible")

              rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j + 1, k))
              rhoM = 0.5 * (var(i, j, k, 1) + var(i, j + 1, k, 1))

              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  ! Adjust for 3D fields.
                  rhoStratEdgeF = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, &
                      j + 1, k))
                  rhoM_1 = rhoM_1 + rhoStratEdgeF
                  rhoM = rhoM + rhoStratEdgeF
                else
                  rhoM_1 = rhoM_1 + rhoStrat(k)
                  rhoM = rhoM + rhoStrat(k)
                end if
              end if

            case("Boussinesq")
              rhoM_1 = rho00
              rhoM = rho00
            case default
              stop "momentumPredictor: unkown case model."
            end select

            ! velocity and momentum at t(m-1)
            vM_1 = var(i, j, k, 3)
            momM_1 = rhoM_1 * vM_1

            ! q(m-1) -> q(m)
            q(i, j, k, 2) = dt * F + alpha(m) * q(i, j, k, 2)

            ! rhoV(m-1) -> rhoV(m)
            momM = momM_1 + beta(m) * q(i, j, k, 2)

            ! calc v(m,*)
            vAst = momM / rhoM

            ! vAst -> var
            var(i, j, k, 3) = vAst
          end do
        end do
      end do
    else if(mmp_mod == "rhs") then
      if(int_mod == "expl") then
        do k = 1, nz
          do j = j0, j1
            do i = 1, nx
              rhov = 0.5 * (var(i, j, k, 1) + var(i, j + 1, k, 1))
              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  rhoStratEdgeF = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, &
                      j + 1, k))
                  rhov = rhov + rhoStratEdgeF
                else
                  rhov = rhov + rhoStrat(k)
                end if
              end if

              if(TestCase == "baroclinic_LC") then
                rhov_e = 0.5 * (var_env(i, j, k, 1) + var_env(i, j + 1, k, 1))

                if(fluctuationMode) then
                  if(topography) then
                    ! TFC FJ
                    rhov_e = rhov_e + 0.5 * (rhoStratTFC(i, j, k) &
                        + rhoStratTFC(i, j + 1, k))
                  else
                    rhov_e = rhov_e + rhoStrat_0(k)
                  end if
                end if
              end if

              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piF = var(i, j + 1, k, 5) - var_env(i, j + 1, k, 5)
                piB = var(i, j, k, 5) - var_env(i, j, k, 5)

                if(TestCase == "baroclinic_LC") then !FS
                  piF_e = var_env(i, j + 1, k, 5)
                  piB_e = var_env(i, j, k, 5)
                end if
              else
                piF = var(i, j + 1, k, 5)
                piB = var(i, j, k, 5)
              end if

              if(topography) then
                ! TFC FJ
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) - var_env(:, :, :, 5)
                end if
                ! Compute values at cell edges.
                pEdgeF = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    + pStratTFC(i, j + 1, k) / jac(i, j + 1, k))
                ! Compute pressure gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  piUUEdgeF = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 2, 3) &
                      * var(i, j, k + 2, 5) + jac(i, j + 1, k + 2) * met(i, j &
                      + 1, k + 2, 2, 3) * var(i, j + 1, k + 2, 5))
                  piUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                      * var(i, j, k + 1, 5) + jac(i, j + 1, k + 1) * met(i, j &
                      + 1, k + 1, 2, 3) * var(i, j + 1, k + 1, 5))
                  piEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) * var(i, &
                      j, k, 5) + jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) &
                      * var(i, j + 1, k, 5))
                  piGrad = kappaInv * MaInv2 * pEdgeF / rhov * ((jac(i, j + 1, &
                      k) * piF - jac(i, j, k) * piB) / dy + (- piUUEdgeF + 4.0 &
                      * piUEdgeF - 3.0 * piEdgeF) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeF = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 2, 3) &
                      * var(i, j, k - 2, 5) + jac(i, j + 1, k - 2) * met(i, j &
                      + 1, k - 2, 2, 3) * var(i, j + 1, k - 2, 5))
                  piDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) &
                      * var(i, j, k - 1, 5) + jac(i, j + 1, k - 1) * met(i, j &
                      + 1, k - 1, 2, 3) * var(i, j + 1, k - 1, 5))
                  piEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) * var(i, &
                      j, k, 5) + jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) &
                      * var(i, j + 1, k, 5))
                  piGrad = kappaInv * MaInv2 * pEdgeF / rhov * ((jac(i, j + 1, &
                      k) * piF - jac(i, j, k) * piB) / dy + (piDDEdgeF - 4.0 &
                      * piDEdgeF + 3.0 * piEdgeF) * 0.5 / dz)
                else
                  piUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                      * var(i, j, k + 1, 5) + jac(i, j + 1, k + 1) * met(i, j &
                      + 1, k + 1, 2, 3) * var(i, j + 1, k + 1, 5))
                  piDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) &
                      * var(i, j, k - 1, 5) + jac(i, j + 1, k - 1) * met(i, j &
                      + 1, k - 1, 2, 3) * var(i, j + 1, k - 1, 5))
                  piGrad = kappaInv * MaInv2 * pEdgeF / rhov * ((jac(i, j + 1, &
                      k) * piF - jac(i, j, k) * piB) / dy + (piUEdgeF &
                      - piDEdgeF) * 0.5 / dz)
                end if
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) + var_env(:, :, :, 5)
                end if
              else
                piGrad = kappaInv * MaInv2 * Pstrat(k) / rhov * (piF - piB) / dy
              end if

              if(TestCase == "baroclinic_LC") then !FS
                if(topography) then
                  ! TFC FJ
                  ! Compute values at cell edges.
                  pEdgeF = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                      + pStratTFC(i, j + 1, k) / jac(i, j + 1, k))
                  ! Compute pressure gradient component.
                  if(k == 1 .and. zBoundary == "solid_wall") then
                    piUUEdgeF = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 2, &
                        3) * var_env(i, j, k + 2, 5) + jac(i, j + 1, k + 2) &
                        * met(i, j + 1, k + 2, 2, 3) * var_env(i, j + 1, k &
                        + 2, 5))
                    piUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, &
                        3) * var_env(i, j, k + 1, 5) + jac(i, j + 1, k + 1) &
                        * met(i, j + 1, k + 1, 2, 3) * var_env(i, j + 1, k &
                        + 1, 5))
                    piEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) &
                        * var_env(i, j, k, 5) + jac(i, j + 1, k) * met(i, j &
                        + 1, k, 2, 3) * var_env(i, j + 1, k, 5))
                    piGrad = piGrad + kappaInv * MaInv2 * (pEdgeF / rhov &
                        - pEdgeF / rhov_e) * ((jac(i, j + 1, k) * piF_e &
                        - jac(i, j, k) * piB_e) / dy + (- piUUEdgeF + 4.0 &
                        * piUEdgeF - 3.0 * piEdgeF) * 0.5 / dz)
                  else if(k == nz .and. zBoundary == "solid_wall") then
                    piDDEdgeF = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 2, &
                        3) * var_env(i, j, k - 2, 5) + jac(i, j + 1, k - 2) &
                        * met(i, j + 1, k - 2, 2, 3) * var_env(i, j + 1, k &
                        - 2, 5))
                    piDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, &
                        3) * var_env(i, j, k - 1, 5) + jac(i, j + 1, k - 1) &
                        * met(i, j + 1, k - 1, 2, 3) * var_env(i, j + 1, k &
                        - 1, 5))
                    piEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) &
                        * var_env(i, j, k, 5) + jac(i, j + 1, k) * met(i, j &
                        + 1, k, 2, 3) * var_env(i, j + 1, k, 5))
                    piGrad = piGrad + kappaInv * MaInv2 * (pEdgeF / rhov &
                        - pEdgeF / rhov_e) * ((jac(i, j + 1, k) * piF_e &
                        - jac(i, j, k) * piB_e) / dy + (piDDEdgeF - 4.0 &
                        * piDEdgeF + 3.0 * piEdgeF) * 0.5 / dz)
                  else
                    piUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, &
                        3) * var_env(i, j, k + 1, 5) + jac(i, j + 1, k + 1) &
                        * met(i, j + 1, k + 1, 2, 3) * var_env(i, j + 1, k &
                        + 1, 5))
                    piDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, &
                        3) * var_env(i, j, k - 1, 5) + jac(i, j + 1, k - 1) &
                        * met(i, j + 1, k - 1, 2, 3) * var_env(i, j + 1, k &
                        - 1, 5))
                    piGrad = piGrad + kappaInv * MaInv2 * (pEdgeF / rhov &
                        - pEdgeF / rhov_e) * ((jac(i, j + 1, k) * piF_e &
                        - jac(i, j, k) * piB_e) / dy + (piUEdgeF - piDEdgeF) &
                        * 0.5 / dz)
                  end if
                else
                  piGrad = piGrad + kappaInv * MaInv2 * (Pstrat(k) / rhov &
                      - Pstrat_0(k) / rhov_e) * (piF_e - piB_e) / dy
                end if
              end if

              ! gravity-wave forcing
              if(raytracer .or. (testCase == "mountainwave")) then
                volfcy = 0.5 * (force(i, j, k, 2) + force(i, j + 1, k, 2))
              else
                volfcy = 0.0
              end if

              ! vstar
              if(TestCase == "baroclinic_LC") then
                uhorx = 0.25 * (var(i - 1, j, k, 2) + var(i - 1, j + 1, k, 2) &
                    - var_env(i - 1, j, k, 2) - var_env(i - 1, j + 1, k, 2) &
                    + var(i, j, k, 2) + var(i, j + 1, k, 2) - var_env(i, j, k, &
                    2) - var_env(i, j + 1, k, 2))
              else
                uhorx = 0.25 * (var(i - 1, j, k, 2) + var(i - 1, j + 1, k, 2) &
                    + var(i, j, k, 2) + var(i, j + 1, k, 2))
              end if

              vhory = var(i, j, k, 3)

              f_cor_v = 0.5 * (f_cor_nd(j) + f_cor_nd(j + 1))

              if(topography) then
                ! TFC FJ
                ! Coriolis force is integrated on LHS.
                vAst = vhory + dt * (- piGrad + volfcy / rhov)
              else
                if(testCase == "SkamarockKlemp94") then
                  vAst = vhory + dt * (- f_cor_v * (uhorx &
                      - backgroundFlow_dim(1) / uRef) - piGrad + volfcy / rhov)
                else
                  vAst = vhory + dt * (- f_cor_v * uhorx - piGrad + volfcy &
                      / rhov)
                end if
              end if

              ! if (topography) then
              !    ! Rayleigh damping for topography (immersed boundary)
              !
              !    if(k < kbl_topo(i,j,2)) then
              !       vAst = vAst - dt* alprlx*vhory
              !      else if(k == kbl_topo(i,j,2)) then
              !       call wind_ip(var, &
              !                  & x_ip(i,j,2),y_ip(i,j,2),z_ip(i,j,2),&
              !                  &'v',u_ip,v_ip,w_ip)
              !
              !        v_ip_n &
              !        = (u_ip*dhdx(i,j,2) + v_ip*dhdy(i,j,2) - w_ip)&
              !          *dhdy(i,j,2) &
              !          /(1 + dhdx(i,j,2)**2 + dhdy(i,j,2)**2)
              !
              !        v_ip_t = v_ip-v_ip_n
              !
              !        v_rp_t = velocity_reconst_t(i,j,2)*v_ip_t
              !        v_rp_n = velocity_reconst_n(i,j,2)*v_ip_n
              !        v_rp = v_rp_t + v_rp_n
              !
              !        vAst &
              !        = vAst - dt*alprlx * (var(i,j,k,3)-v_rp)
              !    end if
              ! end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  vAst = vAst - dt * 0.5 * (kv_hs(j, k) + kv_hs(j + 1, k)) &
                      * vhory
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                vAst = vAst - dt * 0.5 * (kr_sp(j, k) + kr_sp(j + 1, k)) * vhory
              end if

              var(i, j, k, 3) = vAst
            end do
          end do
        end do
      else if(int_mod == "impl") then
        do k = 1, nz
          do j = j0, j1
            do i = 1, nx
              rhoum0 = 0.5 * (var(i, j, k, 1) + var(i - 1, j, k, 1))
              rhou00 = 0.5 * (var(i + 1, j, k, 1) + var(i, j, k, 1))
              rhoum1 = 0.5 * (var(i, j + 1, k, 1) + var(i - 1, j + 1, k, 1))
              rhou01 = 0.5 * (var(i + 1, j + 1, k, 1) + var(i, j + 1, k, 1))

              rhov = 0.5 * (var(i, j, k, 1) + var(i, j + 1, k, 1))

              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  rhoStratEdgeF = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, &
                      j + 1, k))
                  rhov = rhov + rhoStratEdgeF
                else
                  rhov = rhov + rhoStrat(k)
                  rhoum0 = rhoum0 + rhoStrat(k)
                  rhou00 = rhou00 + rhoStrat(k)
                  rhoum1 = rhoum1 + rhoStrat(k)
                  rhou01 = rhou01 + rhoStrat(k)
                end if
              end if

              if(TestCase == "baroclinic_LC") then
                rhoum0_e = 0.5 * (var_env(i, j, k, 1) + var_env(i - 1, j, k, 1))
                rhou00_e = 0.5 * (var_env(i + 1, j, k, 1) + var_env(i, j, k, 1))
                rhoum1_e = 0.5 * (var_env(i, j + 1, k, 1) + var_env(i - 1, j &
                    + 1, k, 1))
                rhou01_e = 0.5 * (var_env(i + 1, j + 1, k, 1) + var_env(i, j &
                    + 1, k, 1))

                rhov_e = 0.5 * (var_env(i, j, k, 1) + var_env(i, j + 1, k, 1))

                if(fluctuationMode) then
                  if(topography) then
                    ! TFC FJ
                    rhov_e = rhov_e + 0.5 * (rhoStratTFC(i, j, k) &
                        + rhoStratTFC(i, j + 1, k))
                  else
                    rhov_e = rhov_e + rhoStrat_0(k)
                    rhoum0_e = rhoum0_e + rhoStrat_0(k)
                    rhou00_e = rhou00_e + rhoStrat_0(k)
                    rhoum1_e = rhoum1_e + rhoStrat_0(k)
                    rhou01_e = rhou01_e + rhoStrat_0(k)
                  end if
                end if
              end if

              !--- pressure gradient terms -> piGradx, piGrady
              if(TestCase == "baroclinic_LC") then
                piGradx = kappaInv * MaInv2 * 0.25 * (Pstrat(k) / rhou00 &
                    * (var(i + 1, j, k, 5) - var(i, j, k, 5) - var_env(i + 1, &
                    j, k, 5) + var_env(i, j, k, 5)) / dx + Pstrat(k) / rhoum0 &
                    * (var(i, j, k, 5) - var(i - 1, j, k, 5) - var_env(i, j, &
                    k, 5) + var_env(i - 1, j, k, 5)) / dx + Pstrat(k) / rhou01 &
                    * (var(i + 1, j + 1, k, 5) - var(i, j + 1, k, 5) &
                    - var_env(i + 1, j + 1, k, 5) + var_env(i, j + 1, k, 5)) &
                    / dx + Pstrat(k) / rhoum1 * (var(i, j + 1, k, 5) - var(i &
                    - 1, j + 1, k, 5) - var_env(i, j + 1, k, 5) + var_env(i &
                    - 1, j + 1, k, 5)) / dx)

                piGradx = piGradx + kappaInv * MaInv2 * 0.25 * ((Pstrat(k) &
                    / rhou00 - Pstrat_0(k) / rhou00_e) * (var_env(i + 1, j, k, &
                    5) - var_env(i, j, k, 5)) / dx + (Pstrat(k) / rhoum0 &
                    - Pstrat_0(k) / rhoum0_e) * (var_env(i, j, k, 5) &
                    - var_env(i - 1, j, k, 5)) / dx + (Pstrat(k) / rhou01 &
                    - Pstrat_0(k) / rhou01_e) * (var_env(i + 1, j + 1, k, 5) &
                    - var_env(i, j + 1, k, 5)) / dx + (Pstrat(k) / rhoum1 &
                    - Pstrat_0(k) / rhoum1_e) * (var_env(i, j + 1, k, 5) &
                    - var_env(i - 1, j + 1, k, 5)) / dx)
              else
                piGradx = kappaInv * MaInv2 * 0.25 * (Pstrat(k) / rhou00 &
                    * (var(i + 1, j, k, 5) - var(i, j, k, 5)) / dx + Pstrat(k) &
                    / rhoum0 * (var(i, j, k, 5) - var(i - 1, j, k, 5)) / dx &
                    + Pstrat(k) / rhou01 * (var(i + 1, j + 1, k, 5) - var(i, j &
                    + 1, k, 5)) / dx + Pstrat(k) / rhoum1 * (var(i, j + 1, k, &
                    5) - var(i - 1, j + 1, k, 5)) / dx)
              end if

              if(TestCase == "baroclinic_LC") then
                piF = var(i, j + 1, k, 5) - var_env(i, j + 1, k, 5)
                piB = var(i, j, k, 5) - var_env(i, j, k, 5)

                piF_e = var_env(i, j + 1, k, 5)
                piB_e = var_env(i, j, k, 5)
              else
                piF = var(i, j + 1, k, 5)
                piB = var(i, j, k, 5)
              end if

              if(topography) then
                ! TFC FJ
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) - var_env(:, :, :, 5)
                end if
                ! Compute values at cell edges.
                pEdgeF = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    + pStratTFC(i, j + 1, k) / jac(i, j + 1, k))
                ! Compute pressure gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  piUUEdgeF = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 2, 3) &
                      * var(i, j, k + 2, 5) + jac(i, j + 1, k + 2) * met(i, j &
                      + 1, k + 2, 2, 3) * var(i, j + 1, k + 2, 5))
                  piUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                      * var(i, j, k + 1, 5) + jac(i, j + 1, k + 1) * met(i, j &
                      + 1, k + 1, 2, 3) * var(i, j + 1, k + 1, 5))
                  piEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) * var(i, &
                      j, k, 5) + jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) &
                      * var(i, j + 1, k, 5))
                  piGradY = kappaInv * MaInv2 * pEdgeF / rhov * ((jac(i, j &
                      + 1, k) * piF - jac(i, j, k) * piB) / dy + (- piUUEdgeF &
                      + 4.0 * piUEdgeF - 3.0 * piEdgeF) * 0.5 / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  piDDEdgeF = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 2, 3) &
                      * var(i, j, k - 2, 5) + jac(i, j + 1, k - 2) * met(i, j &
                      + 1, k - 2, 2, 3) * var(i, j + 1, k - 2, 5))
                  piDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) &
                      * var(i, j, k - 1, 5) + jac(i, j + 1, k - 1) * met(i, j &
                      + 1, k - 1, 2, 3) * var(i, j + 1, k - 1, 5))
                  piEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) * var(i, &
                      j, k, 5) + jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) &
                      * var(i, j + 1, k, 5))
                  piGradY = kappaInv * MaInv2 * pEdgeF / rhov * ((jac(i, j &
                      + 1, k) * piF - jac(i, j, k) * piB) / dy + (piDDEdgeF &
                      - 4.0 * piDEdgeF + 3.0 * piEdgeF) * 0.5 / dz)
                else
                  piUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                      * var(i, j, k + 1, 5) + jac(i, j + 1, k + 1) * met(i, j &
                      + 1, k + 1, 2, 3) * var(i, j + 1, k + 1, 5))
                  piDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) &
                      * var(i, j, k - 1, 5) + jac(i, j + 1, k - 1) * met(i, j &
                      + 1, k - 1, 2, 3) * var(i, j + 1, k - 1, 5))
                  piGradY = kappaInv * MaInv2 * pEdgeF / rhov * ((jac(i, j &
                      + 1, k) * piF - jac(i, j, k) * piB) / dy + (piUEdgeF &
                      - piDEdgeF) * 0.5 / dz)
                end if
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) + var_env(:, :, :, 5)
                end if
              else
                piGrady = kappaInv * MaInv2 * Pstrat(k) / rhov * (piF - piB) &
                    / dy
              end if

              if(TestCase == "baroclinic_LC") then !FS
                if(topography) then
                  ! TFC FJ
                  ! Compute values at cell edges.
                  pEdgeF = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                      + pStratTFC(i, j + 1, k) / jac(i, j + 1, k))
                  ! Compute pressure gradient component.
                  if(k == 1 .and. zBoundary == "solid_wall") then
                    piUUEdgeF = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 2, &
                        3) * var_env(i, j, k + 2, 5) + jac(i, j + 1, k + 2) &
                        * met(i, j + 1, k + 2, 2, 3) * var_env(i, j + 1, k &
                        + 2, 5))
                    piUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, &
                        3) * var_env(i, j, k + 1, 5) + jac(i, j + 1, k + 1) &
                        * met(i, j + 1, k + 1, 2, 3) * var_env(i, j + 1, k &
                        + 1, 5))
                    piEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) &
                        * var_env(i, j, k, 5) + jac(i, j + 1, k) * met(i, j &
                        + 1, k, 2, 3) * var_env(i, j + 1, k, 5))
                    piGradY = piGradY + kappaInv * MaInv2 * (pEdgeF / rhov &
                        - pEdgeF / rhov_e) * ((jac(i, j + 1, k) * piF_e &
                        - jac(i, j, k) * piB_e) / dy + (- piUUEdgeF + 4.0 &
                        * piUEdgeF - 3.0 * piEdgeF) * 0.5 / dz)
                  else if(k == nz .and. zBoundary == "solid_wall") then
                    piDDEdgeF = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 2, &
                        3) * var_env(i, j, k - 2, 5) + jac(i, j + 1, k - 2) &
                        * met(i, j + 1, k - 2, 2, 3) * var_env(i, j + 1, k &
                        - 2, 5))
                    piDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, &
                        3) * var_env(i, j, k - 1, 5) + jac(i, j + 1, k - 1) &
                        * met(i, j + 1, k - 1, 2, 3) * var_env(i, j + 1, k &
                        - 1, 5))
                    piEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) &
                        * var_env(i, j, k, 5) + jac(i, j + 1, k) * met(i, j &
                        + 1, k, 2, 3) * var_env(i, j + 1, k, 5))
                    piGradY = piGradY + kappaInv * MaInv2 * (pEdgeF / rhov &
                        - pEdgeF / rhov_e) * ((jac(i, j + 1, k) * piF_e &
                        - jac(i, j, k) * piB_e) / dy + (piDDEdgeF - 4.0 &
                        * piDEdgeF + 3.0 * piEdgeF) * 0.5 / dz)
                  else
                    piUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, &
                        3) * var_env(i, j, k + 1, 5) + jac(i, j + 1, k + 1) &
                        * met(i, j + 1, k + 1, 2, 3) * var_env(i, j + 1, k &
                        + 1, 5))
                    piDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, &
                        3) * var_env(i, j, k - 1, 5) + jac(i, j + 1, k - 1) &
                        * met(i, j + 1, k - 1, 2, 3) * var_env(i, j + 1, k &
                        - 1, 5))
                    piGradY = piGradY + kappaInv * MaInv2 * (pEdgeF / rhov &
                        - pEdgeF / rhov_e) * ((jac(i, j + 1, k) * piF_e &
                        - jac(i, j, k) * piB_e) / dy + (piUEdgeF - piDEdgeF) &
                        * 0.5 / dz)
                  end if
                else
                  piGrady = piGrady + kappaInv * MaInv2 * (Pstrat(k) / rhov &
                      - Pstrat_0(k) / rhov_e) * (piF_e - piB_e) / dy
                end if
              end if

              ! gravity-wave forcing
              if(raytracer .or. (testCase == "mountainwave")) then
                volfcx = 0.5 * (force(i, j, k, 1) + force(i + 1, j, k, 1))
                volfcy = 0.5 * (force(i, j, k, 2) + force(i, j + 1, k, 2))
              else
                volfcx = 0.0
                volfcy = 0.0
              end if

              ! vstar
              if(TestCase == "baroclinic_LC") then
                uhorx = 0.25 * (var(i - 1, j, k, 2) + var(i - 1, j + 1, k, 2) &
                    - var_env(i - 1, j, k, 2) - var_env(i - 1, j + 1, k, 2) &
                    + var(i, j, k, 2) + var(i, j + 1, k, 2) - var_env(i, j, k, &
                    2) - var_env(i, j + 1, k, 2))
              else
                uhorx = 0.25 * (var(i - 1, j, k, 2) + var(i - 1, j + 1, k, 2) &
                    + var(i, j, k, 2) + var(i, j + 1, k, 2))
              end if

              vhory = var(i, j, k, 3)

              facv = 1.0

              ! if (topography) then
              !    ! Rayleigh damping for topography (immersed boundary)
              !
              !    if(k < kbl_topo(i,j,2)) then
              !       facv = facv + dt*alprlx
              !      else if(k == kbl_topo(i,j,2)) then
              !        stop'implementation topography into &
              !            &semi-implicit time step still to be done'
              !    end if
              ! end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  facv = facv + dt * 0.5 * (kv_hs(j, k) + kv_hs(j + 1, k))
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                facv = facv + dt * 0.5 * (kr_sp(j, k) + kr_sp(j + 1, k))
              end if

              facu = facv

              if(topography) then
                ! TFC FJ
                ! Coriolis force is integrated on LHS.
                vAst = 1.0 / facv * (vhory + dt * (- piGradY + volfcy / rhov))
              else
                if(testCase == "SkamarockKlemp94") then
                  vAst = 1.0 / (facu * facv + (0.5 * (f_cor_nd(j) + f_cor_nd(j &
                      + 1)) * dt) ** 2) * (- 0.5 * (f_cor_nd(j) + f_cor_nd(j &
                      + 1)) * dt * ((uhorx - backgroundFlow_dim(1) / uRef) &
                      + dt * (volfcx / rhov - piGradx)) + facu * (vhory + dt &
                      * (volfcy / rhov - piGrady)))
                else
                  vAst = 1.0 / (facu * facv + (0.5 * (f_cor_nd(j) + f_cor_nd(j &
                      + 1)) * dt) ** 2) * (- 0.5 * (f_cor_nd(j) + f_cor_nd(j &
                      + 1)) * dt * (uhorx + dt * (volfcx / rhov - piGradx)) &
                      + facu * (vhory + dt * (volfcy / rhov - piGrady)))
                end if
              end if

              var(i, j, k, 3) = vAst
            end do
          end do
        end do
      else
        stop 'ERROR: unknown int_mod'
      end if

      ! now the new u can be put into the proper array
      var(:, :, :, 2) = usave(:, :, :)
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
            fR = flux(i, j, k, 1, 4)
            fL = flux(i - 1, j, k, 1, 4)
            gF = flux(i, j, k, 2, 4)
            gB = flux(i, j - 1, k, 2, 4)
            hU = flux(i, j, k, 3, 4)
            hD = flux(i, j, k - 1, 3, 4)
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

            ! TFC FJ
            ! Adjust vertical momentum flux divergence.
            if(topography) then
              ! Adjust Cartesian vertical momentum flux divergence.
              jacEdgeU = 0.5 * (jac(i, j, k) + jac(i, j, k + 1))
              fluxDiff = fluxDiff / jacEdgeU
              ! Compute zonal momentum flux divergences.
              do ll = 0, 1
                do mm = 0, 1
                  fR = flux(i - ll, j, k + mm, 1, 2)
                  fL = flux(i - 1 - ll, j, k + mm, 1, 2)
                  gF = flux(i - ll, j, k + mm, 2, 2)
                  gB = flux(i - ll, j - 1, k + mm, 2, 2)
                  hU = flux(i - ll, j, k + mm, 3, 2)
                  hD = flux(i - ll, j, k - 1 + mm, 3, 2)
                  fluxDiffU(ll, mm) = (fR - fL) / dx + (gF - gB) / dy + (hU &
                      - hD) / dz
                  jacEdgeR = 0.5 * (jac(i - ll, j, k + mm) + jac(i + 1 - ll, &
                      j, k + mm))
                  fluxDiffU(ll, mm) = fluxDiffU(ll, mm) / jacEdgeR
                end do
              end do
              ! Compute meridional momentum flux divergences.
              do ll = 0, 1
                do mm = 0, 1
                  fR = flux(i, j - ll, k + mm, 1, 3)
                  fL = flux(i - 1, j - ll, k + mm, 1, 3)
                  gF = flux(i, j - ll, k + mm, 2, 3)
                  gB = flux(i, j - 1 - ll, k + mm, 2, 3)
                  hU = flux(i, j - ll, k + mm, 3, 3)
                  hD = flux(i, j - ll, k - 1 + mm, 3, 3)
                  fluxDiffV(ll, mm) = (fR - fL) / dx + (gF - gB) / dy + (hU &
                      - hD) / dz
                  jacEdgeF = 0.5 * (jac(i, j - ll, k + mm) + jac(i, j + 1 &
                      - ll, k + mm))
                  fluxDiffV(ll, mm) = fluxDiffV(ll, mm) / jacEdgeF
                end do
              end do
              ! Compute transformed vertical momentum flux divergence.
              ! metEdgeR1 = 0.5 * (met(i, j, k, 1, 3) &
              !             + met(i + 1, j, k, 1, 3))
              ! metEdgeR2 = 0.5 * (met(i, j, k + 1, 1, 3) &
              !             + met(i + 1, j, k + 1, 1, 3))
              ! metEdgeL1 = 0.5 * (met(i, j, k, 1, 3) &
              !             + met(i - 1, j, k, 1, 3))
              ! metEdgeL2 = 0.5 * (met(i, j, k + 1, 1, 3) &
              !             + met(i - 1, j, k + 1, 1, 3))
              ! metEdgeF1 = 0.5 * (met(i, j, k, 2, 3) &
              !             + met(i, j + 1, k, 2, 3))
              ! metEdgeF2 = 0.5 * (met(i, j, k + 1, 2, 3) &
              !             + met(i, j + 1, k + 1, 2, 3))
              ! metEdgeB1 = 0.5 * (met(i, j, k, 2, 3) &
              !             + met(i, j - 1, k, 2, 3))
              ! metEdgeB2 = 0.5 * (met(i, j, k + 1, 2, 3) &
              !             + met(i, j - 1, k + 1, 2, 3))
              ! fluxDiff = fluxDiff / jacEdgeU &
              !            + 0.25 * (metEdgeR1 * fluxDiffU(0, 0) &
              !            + metEdgeR2 * fluxDiffU(0, 1) &
              !            + metEdgeL1 * fluxDiffU(1, 0) &
              !            + metEdgeL2 * fluxDiffU(1, 1)) &
              !            + 0.25 * (metEdgeF1 * fluxDiffV(0, 0) &
              !            + metEdgeF2 * fluxDiffV(0, 1) &
              !            + metEdgeB1 * fluxDiffV(1, 0) &
              !            + metEdgeB2 * fluxDiffV(1, 1))
              ! fluxDiff = fluxDiff / jacEdgeU &
              !            + 0.25 * (met(i, j, k, 1, 3) * (fluxDiffU(0, 0) &
              !            + fluxDiffU(1, 0)) + met(i, j, k + 1, 1, 3) &
              !            * (fluxDiffU(0, 1) + fluxDiffU(1, 1))) &
              !            + 0.25 * (met(i, j, k, 2, 3) * (fluxDiffV(0, 0) &
              !            + fluxDiffV(1, 0)) + met(i, j, k + 1, 2, 3) &
              !            * (fluxDiffV(0, 1) + fluxDiffV(1, 1)))
              fluxDiff = trafoTFC(i, j, k, fluxDiffU(0, 0), fluxDiffU(0, 1), &
                  fluxDiffU(1, 0), fluxDiffU(1, 1), fluxDiffV(0, 0), &
                  fluxDiffV(0, 1), fluxDiffV(1, 0), fluxDiffV(1, 1), fluxDiff, &
                  "tfc")
            end if

            if(mmp_mod == "tot") then
              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piU = var(i, j, k + 1, 5) - var_env(i, j, k + 1, 5)
                piD = var(i, j, k, 5) - var_env(i, j, k, 5)
              else
                piU = var(i, j, k + 1, 5)
                piD = var(i, j, k, 5)
              end if

              if(topography) then
                ! TFC FJ
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) - var_env(:, :, :, 5)
                end if
                ! Compute values at cell edges.
                pEdgeU = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    + pStratTFC(i, j, k + 1) / jac(i, j, k + 1))
                piREdgeU = 0.5 * (jac(i + 1, j, k) * met(i + 1, j, k, 3, 1) &
                    * var(i + 1, j, k, 5) + jac(i + 1, j, k + 1) * met(i + 1, &
                    j, k + 1, 3, 1) * var(i + 1, j, k + 1, 5))
                piLEdgeU = 0.5 * (jac(i - 1, j, k) * met(i - 1, j, k, 3, 1) &
                    * var(i - 1, j, k, 5) + jac(i - 1, j, k + 1) * met(i - 1, &
                    j, k + 1, 3, 1) * var(i - 1, j, k + 1, 5))
                piFEdgeU = 0.5 * (jac(i, j + 1, k) * met(i, j + 1, k, 3, 2) &
                    * var(i, j + 1, k, 5) + jac(i, j + 1, k + 1) * met(i, j &
                    + 1, k + 1, 3, 2) * var(i, j + 1, k + 1, 5))
                piBEdgeU = 0.5 * (jac(i, j - 1, k) * met(i, j - 1, k, 3, 2) &
                    * var(i, j - 1, k, 5) + jac(i, j - 1, k + 1) * met(i, j &
                    - 1, k + 1, 3, 2) * var(i, j - 1, k + 1, 5))
                chris11EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                    1) * var(i, j, k, 5) + pStratTFC(i, j, k + 1) * chris(i, &
                    j, k + 1, 1, 1) * var(i, j, k + 1, 5))
                chris22EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                    2) * var(i, j, k, 5) + pStratTFC(i, j, k + 1) * chris(i, &
                    j, k + 1, 2, 2) * var(i, j, k + 1, 5))
                chris13EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                    3) * met(i, j, k, 1, 3) * var(i, j, k, 5) + pStratTFC(i, &
                    j, k + 1) * chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, 1, &
                    3) * var(i, j, k + 1, 5))
                chris23EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                    3) * met(i, j, k, 2, 3) * var(i, j, k, 5) + pStratTFC(i, &
                    j, k + 1) * chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, 2, &
                    3) * var(i, j, k + 1, 5))
                ! Compute pressure gradient component.
                piGrad = kappaInv * MaInv2 * pEdgeU * ((piREdgeU - piLEdgeU) &
                    * 0.5 / dx + (piFEdgeU - piBEdgeU) * 0.5 / dy + (jac(i, j, &
                    k + 1) * met(i, j, k + 1, 3, 3) * var(i, j, k + 1, 5) &
                    - jac(i, j, k) * met(i, j, k, 3, 3) * var(i, j, k, 5)) &
                    / dz) + kappaInv * MaInv2 * (chris11EdgeU + chris22EdgeU &
                    + 2.0 * (chris13EdgeU + chris23EdgeU))
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) + var_env(:, :, :, 5)
                end if
              else
                piGrad = 0.5 * kappaInv * MaInv2 * (Pstrat(k) + Pstrat(k + 1)) &
                    * (piU - piD) / dz
              end if

              if(TestCase == "baroclinic_LC") then !FS
                piGrad = piGrad + 0.5 * kappaInv * MaInv2 * (Pstrat(k) &
                    + Pstrat(k + 1) - pStrat_0(k) - pStrat_0(k + 1)) &
                    * (var_env(i, j, k + 1, 5) - var_env(i, j, k, 5)) / dz
              end if

              !---- volume forces
              volForce = 0.5 * (force(i, j, k, 3) + force(i, j, k + 1, 3))

              if(TestCase == "baroclinic_LC") then
                if(model == "pseudo_incompressible") then
                  drho_e = 0.5 * (var_env(i, j, k, 1) + var_env(i, j, k + 1, 1))

                  if(.not. fluctuationMode) then
                    drho_e = drho_e - rhoStratTilde(k)
                  end if
                else if(model == "Boussinesq") then
                  stop 'ERROR: baroclinic LC not ready yet for  Boussinesq'
                else
                  stop "momentumPredictor: unkown model."
                end if

                if(topography) then
                  ! TFC FJ
                  volForce = volForce + FrInv2 * drho_e / (0.5 * (jac(i, j, k) &
                      + jac(i, j, k + 1)))
                else
                  volForce = volForce + FrInv2 * drho_e
                end if
              end if

              !    if (topography) then
              !       ! Rayleigh damping for topography (immersed boundary)
              !
              !       if (model == "pseudo_incompressible") then
              !           rhoM_1 = 0.5 * (rhoOld(i,j,k) + rhoOld(i,j,k+1))
              !
              !           if( fluctuationMode ) then
              !              rhoM_1 = rhoM_1 + rhoStratTilde(k)
              !           end if
              !          else if (model == "Boussinesq") then
              !           rhoM_1 = rho00
              !          else
              !           stop"momentumPredictor: unkown model."
              !       end if
              !
              !       if(k < kbl_topo(i,j,3)) then
              !          volForce = volForce - alprlx * rhoM_1*var(i,j,k,4)
              !         else if(k == kbl_topo(i,j,3)) then
              !          call wind_ip(var, &
              !                     & x_ip(i,j,3),y_ip(i,j,3),z_ip(i,j,3),&
              !                     & 'w',u_ip,v_ip,w_ip)
              !
              !           w_ip_n &
              !           = (- u_ip*dhdx(i,j,3) - v_ip*dhdy(i,j,3) + w_ip)&
              !             /(1 + dhdx(i,j,3)**2 + dhdy(i,j,3)**2)
              !
              !           w_ip_t = w_ip-w_ip_n
              !
              !           w_rp_t = velocity_reconst_t(i,j,3)*w_ip_t
              !           w_rp_n = velocity_reconst_n(i,j,3)*w_ip_n
              !           w_rp = w_rp_t+w_rp_n
              !
              !           volForce &
              !           = volForce - alprlx * rhoM_1*(var(i,j,k,4)-w_rp)
              !       end if
              !    end if

            end if

            if(TestCase == "baroclinic_LC") then
              if(background == "HeldSuarez") then
                ! Rayleigh damping
                if(model == "pseudo_incompressible") then
                  rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j, k + 1))

                  if(fluctuationMode) then
                    if(topography) then
                      ! TFC FJ
                      rhoM_1 = rhoM_1 + 0.5 * (rhoStratTFC(i, j, k) &
                          + rhoStratTFC(i, j, k + 1))
                    else
                      rhoM_1 = rhoM_1 + rhoStratTilde(k)
                    end if
                  end if
                else if(model == "Boussinesq") then
                  rhoM_1 = rho00
                else
                  stop "momentumPredictor: unkown model."
                end if

                volForce = volforce - 0.5 * (kw_hs(k) + kw_hs(k + 1)) * rhoM_1 &
                    * var(i, j, k, 4)
              end if
            end if

            ! TFC FJ
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
                vC = vC - 0.5 * (var_env(i, j, k, 3) + var_env(i, j - 1, k, 3))
                vU = vU - 0.5 * (var_env(i, j, k + 1, 3) + var_env(i, j - 1, k &
                    + 1, 3))
                uC = uC - 0.5 * (var_env(i, j, k, 2) + var_env(i - 1, j, k, 2))
                uU = vU - 0.5 * (var_env(i, j, k + 1, 2) + var_env(i - 1, j, k &
                    + 1, 2))
              end if
              volForce = volForce + 0.5 * f_cor_nd(j) * (met(i, j, k, 1, 3) &
                  * (rhoOld(i, j, k) + rhoStratTFC(i, j, k)) * vC + met(i, j, &
                  k + 1, 1, 3) * (rhoOld(i, j, k + 1) + rhoStratTFC(i, j, k &
                  + 1)) * vU) - 0.5 * f_cor_nd(j) * (met(i, j, k, 2, 3) &
                  * (rhoOld(i, j, k) + rhoStratTFC(i, j, k)) * uC + met(i, j, &
                  k + 1, 2, 3) * (rhoOld(i, j, k + 1) + rhoStratTFC(i, j, k &
                  + 1)) * uU)
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

            case("pseudo_incompressible")

              rhoM_1 = 0.5 * (rhoOld(i, j, k) + rhoOld(i, j, k + 1)) !rho(m-1)
              rhoM = 0.5 * (var(i, j, k, 1) + var(i, j, k + 1, 1)) !rho(m)

              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  ! Adjust for 3D fields.
                  rhoStratEdgeU = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, &
                      j, k + 1))
                  rhoM_1 = rhoM_1 + rhoStratEdgeU
                  rhoM = rhoM + rhoStratEdgeU
                else
                  rhoM_1 = rhoM_1 + rhoStratTilde(k)
                  rhoM = rhoM + rhoStratTilde(k)
                end if
              end if

            case("Boussinesq")
              rhoM_1 = rho00
              rhoM = rho00
            case default
              stop "momentumPredictor: unkown case model."
            end select

            ! velocity and momentum at t(m-1)
            wM_1 = var(i, j, k, 4)
            momM_1 = rhoM_1 * wM_1

            ! q(m-1) -> q(m)
            q(i, j, k, 3) = dt * F + alpha(m) * q(i, j, k, 3)

            ! rhoW(m-1) -> rhoW(m)
            momM = momM_1 + beta(m) * q(i, j, k, 3)

            ! calc w(m,*)
            wAst = momM / rhoM

            ! wAst -> var
            var(i, j, k, 4) = wAst
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
              rho000 = var(i, j, k, 1)
              rho001 = var(i, j, k + 1, 1)

              rhow = 0.5 * (var(i, j, k, 1) + var(i, j, k + 1, 1))

              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  rhoStratEdgeU = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, &
                      j, k + 1))
                  rho000 = rho000 + rhoStratTFC(i, j, k)
                  rho001 = rho001 + rhoStratTFC(i, j, k + 1)
                  rhow = rhow + rhoStratEdgeU
                else
                  rho000 = rho000 + rhoStrat(k)
                  rho001 = rho001 + rhoStrat(k + 1)

                  rhow = rhow + rhoStratTilde(k)
                end if
              end if

              if(TestCase == "baroclinic_LC") then
                rho000_e = var_env(i, j, k, 1)
                rho001_e = var_env(i, j, k + 1, 1)

                rhow_e = 0.5 * (var_env(i, j, k, 1) + var_env(i, j, k + 1, 1))

                if(fluctuationMode) then
                  if(topography) then
                    ! TFC FJ
                    rho000_e = rho000_e + rhoStratTFC(i, j, k)
                    rho001_e = rho001_e + rhoStratTFC(i, j, k + 1)
                    rhow_e = rhow_e + 0.5 * (rhoStratTFC(i, j, k) &
                        + rhoStratTFC(i, j, k + 1))
                  else
                    rho000_e = rho000_e + rhoStrat_0(k)
                    rho001_e = rho001_e + rhoStrat_0(k + 1)

                    rhow_e = rhow_e + 0.5 * (rhoStrat_0(k) + rhoStrat_0(k + 1))
                  end if
                end if
              end if

              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piU = var(i, j, k + 1, 5) - var_env(i, j, k + 1, 5)
                piD = var(i, j, k, 5) - var_env(i, j, k, 5)

                piU_e = var_env(i, j, k + 1, 5)
                piD_e = var_env(i, j, k, 5)
              else
                piU = var(i, j, k + 1, 5)
                piD = var(i, j, k, 5)
              end if

              if(topography) then
                ! TFC FJ
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) - var_env(:, :, :, 5)
                end if
                ! Compute values at cell edges.
                pEdgeU = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    + pStratTFC(i, j, k + 1) / jac(i, j, k + 1))
                piREdgeU = 0.5 * (jac(i + 1, j, k) * met(i + 1, j, k, 3, 1) &
                    * var(i + 1, j, k, 5) + jac(i + 1, j, k + 1) * met(i + 1, &
                    j, k + 1, 3, 1) * var(i + 1, j, k + 1, 5))
                piLEdgeU = 0.5 * (jac(i - 1, j, k) * met(i - 1, j, k, 3, 1) &
                    * var(i - 1, j, k, 5) + jac(i - 1, j, k + 1) * met(i - 1, &
                    j, k + 1, 3, 1) * var(i - 1, j, k + 1, 5))
                piFEdgeU = 0.5 * (jac(i, j + 1, k) * met(i, j + 1, k, 3, 2) &
                    * var(i, j + 1, k, 5) + jac(i, j + 1, k + 1) * met(i, j &
                    + 1, k + 1, 3, 2) * var(i, j + 1, k + 1, 5))
                piBEdgeU = 0.5 * (jac(i, j - 1, k) * met(i, j - 1, k, 3, 2) &
                    * var(i, j - 1, k, 5) + jac(i, j - 1, k + 1) * met(i, j &
                    - 1, k + 1, 3, 2) * var(i, j - 1, k + 1, 5))
                chris11EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                    1) * var(i, j, k, 5) + pStratTFC(i, j, k + 1) * chris(i, &
                    j, k + 1, 1, 1) * var(i, j, k + 1, 5))
                chris22EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                    2) * var(i, j, k, 5) + pStratTFC(i, j, k + 1) * chris(i, &
                    j, k + 1, 2, 2) * var(i, j, k + 1, 5))
                chris13EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                    3) * met(i, j, k, 1, 3) * var(i, j, k, 5) + pStratTFC(i, &
                    j, k + 1) * chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, 1, &
                    3) * var(i, j, k + 1, 5))
                chris23EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                    3) * met(i, j, k, 2, 3) * var(i, j, k, 5) + pStratTFC(i, &
                    j, k + 1) * chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, 2, &
                    3) * var(i, j, k + 1, 5))
                ! Compute pressure gradient component.
                piGrad = kappaInv * MaInv2 * pEdgeU / rhow * ((piREdgeU &
                    - piLEdgeU) * 0.5 / dx + (piFEdgeU - piBEdgeU) * 0.5 / dy &
                    + (jac(i, j, k + 1) * met(i, j, k + 1, 3, 3) * var(i, j, k &
                    + 1, 5) - jac(i, j, k) * met(i, j, k, 3, 3) * var(i, j, k, &
                    5)) / dz) + kappaInv * MaInv2 / rhow * (chris11EdgeU &
                    + chris22EdgeU + 2.0 * (chris13EdgeU + chris23EdgeU))
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) + var_env(:, :, :, 5)
                end if
              else
                piGrad = kappaInv * MaInv2 * pstw / rhow * (piU - piD) / dz
              end if

              if(TestCase == "baroclinic_LC") then !FS
                if(topography) then
                  ! TFC FJ
                  ! Compute values at cell edges.
                  pEdgeU = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                      + pStratTFC(i, j, k + 1) / jac(i, j, k + 1))
                  piREdgeU = 0.5 * (jac(i + 1, j, k) * met(i + 1, j, k, 3, 1) &
                      * var_env(i + 1, j, k, 5) + jac(i + 1, j, k + 1) * met(i &
                      + 1, j, k + 1, 3, 1) * var_env(i + 1, j, k + 1, 5))
                  piLEdgeU = 0.5 * (jac(i - 1, j, k) * met(i - 1, j, k, 3, 1) &
                      * var_env(i - 1, j, k, 5) + jac(i - 1, j, k + 1) * met(i &
                      - 1, j, k + 1, 3, 1) * var_env(i - 1, j, k + 1, 5))
                  piFEdgeU = 0.5 * (jac(i, j + 1, k) * met(i, j + 1, k, 3, 2) &
                      * var_env(i, j + 1, k, 5) + jac(i, j + 1, k + 1) &
                      * met(i, j + 1, k + 1, 3, 2) * var_env(i, j + 1, k + 1, &
                      5))
                  piBEdgeU = 0.5 * (jac(i, j - 1, k) * met(i, j - 1, k, 3, 2) &
                      * var_env(i, j - 1, k, 5) + jac(i, j - 1, k + 1) &
                      * met(i, j - 1, k + 1, 3, 2) * var_env(i, j - 1, k + 1, &
                      5))
                  chris11EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                      1) * var_env(i, j, k, 5) + pStratTFC(i, j, k + 1) &
                      * chris(i, j, k + 1, 1, 1) * var_env(i, j, k + 1, 5))
                  chris22EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                      2) * var_env(i, j, k, 5) + pStratTFC(i, j, k + 1) &
                      * chris(i, j, k + 1, 2, 2) * var_env(i, j, k + 1, 5))
                  chris13EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                      3) * met(i, j, k, 1, 3) * var_env(i, j, k, 5) &
                      + pStratTFC(i, j, k + 1) * chris(i, j, k + 1, 1, 3) &
                      * met(i, j, k + 1, 1, 3) * var_env(i, j, k + 1, 5))
                  chris23EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                      3) * met(i, j, k, 2, 3) * var_env(i, j, k, 5) &
                      + pStratTFC(i, j, k + 1) * chris(i, j, k + 1, 2, 3) &
                      * met(i, j, k + 1, 2, 3) * var_env(i, j, k + 1, 5))
                  ! Compute pressure gradient component.
                  piGrad = piGrad + kappaInv * MaInv2 * (pEdgeU / rhow &
                      - pEdgeU / rhow_e) * ((piREdgeU - piLEdgeU) * 0.5 / dx &
                      + (piFEdgeU - piBEdgeU) * 0.5 / dy + (jac(i, j, k + 1) &
                      * met(i, j, k + 1, 3, 3) * var_env(i, j, k + 1, 5) &
                      - jac(i, j, k) * met(i, j, k, 3, 3) * var_env(i, j, k, &
                      5)) / dz) + kappaInv * MaInv2 * (1.0 / rhow - 1.0 &
                      / rhow_e) * (chris11EdgeU + chris22EdgeU + 2.0 &
                      * (chris13EdgeU + chris23EdgeU))
                else
                  piGrad = piGrad + kappaInv * MaInv2 * (pstw / rhow - pstw_e &
                      / rhow_e) * (piU_e - piD_e) / dz
                end if
              end if

              ! wstar
              wvert = var(i, j, k, 4)

              if(TestCase == "baroclinic_LC") then
                if(topography) then
                  ! TFC FJ
                  buoy = - g_ndim * 0.5 * (rhopOld(i, j, k) / rho000 / jac(i, &
                      j, k) - var_env(i, j, k, 6) / rho000_e / jac(i, j, k) &
                      + rhopOld(i, j, k + 1) / rho001 / jac(i, j, k + 1) &
                      - var_env(i, j, k + 1, 6) / rho001_e / jac(i, j, k + 1))
                else
                  buoy = - g_ndim * 0.5 * (rhopOld(i, j, k) / rho000 &
                      - var_env(i, j, k, 6) / rho000_e + rhopOld(i, j, k + 1) &
                      / rho001 - var_env(i, j, k + 1, 6) / rho001_e)
                end if
              else
                if(topography) then
                  ! TFC FJ
                  buoy = - g_ndim * 0.5 * (rhopOld(i, j, k) / rho000 / jac(i, &
                      j, k) + rhopOld(i, j, k + 1) / rho001 / jac(i, j, k + 1))
                else
                  buoy = - g_ndim * 0.5 * (rhopOld(i, j, k) / rho000 &
                      + rhopOld(i, j, k + 1) / rho001)
                end if
              end if

              wAst = wvert + dt * (buoy - piGrad)

              ! if (topography) then
              !    ! Rayleigh damping for topography (immersed boundary)
              !
              !    if(k < kbl_topo(i,j,3)) then
              !       wAst = wAst - dt* alprlx*wvert
              !      else if(k == kbl_topo(i,j,3)) then
              !       call wind_ip(var, &
              !                  & x_ip(i,j,3),y_ip(i,j,3),z_ip(i,j,3),&
              !                  & 'w',u_ip,v_ip,w_ip)
              !
              !       w_ip_n &
              !       = (- u_ip*dhdx(i,j,3) - v_ip*dhdy(i,j,3) + w_ip)&
              !         /(1 + dhdx(i,j,3)**2 + dhdy(i,j,3)**2)
              !
              !       w_ip_t = w_ip - w_ip_n
              !
              !       w_rp_t = velocity_reconst_t(i,j,3)*w_ip_t
              !       w_rp_n = velocity_reconst_n(i,j,3)*w_ip_n
              !       w_rp = w_rp_t+w_rp_n
              !
              !       wAst &
              !       = wAst - dt*alprlx * (var(i,j,k,4)-w_rp)
              !    end if
              ! end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping

                  wAst = wAst - dt * 0.5 * (kw_hs(k) + kw_hs(k + 1)) * wvert
                end if
              end if

              if(spongeLayer) then
                if(topography .and. spongeTFC) then
                  ! TFC FJ
                  wAst = wAst - dt * 0.5 * (alphaTFC(i, j, k) + alphaTFC(i, j, &
                      k + 1)) * wvert
                else
                  wAst = wAst - dt * 0.5 * (kr_sp_w(j, k) + kr_sp_w(j, k + 1)) &
                      * wvert
                end if
              end if

              var(i, j, k, 4) = wAst
            end do
          end do
        end do
      else if(int_mod == "impl") then
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

        do k = k0, k1
          pstw = 0.5 * (Pstrat(k) + Pstrat(k + 1))
          pstw_0 = 0.5 * (Pstrat_0(k) + Pstrat_0(k + 1))

          if(TestCase == "baroclinic_LC") then
            pstw_e = 0.5 * (Pstrat_0(k) + Pstrat_0(k + 1))
          end if

          do j = 1, ny
            do i = 1, nx
              rho000 = var(i, j, k, 1)
              rho001 = var(i, j, k + 1, 1)

              rhow = 0.5 * (var(i, j, k, 1) + var(i, j, k + 1, 1))

              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  rhoStratEdgeU = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, &
                      j, k + 1))
                  rho000 = rho000 + rhoStratTFC(i, j, k)
                  rho001 = rho001 + rhoStratTFC(i, j, k + 1)
                  rhow = rhow + rhoStratEdgeU
                else
                  rho000 = rho000 + rhoStrat(k)
                  rho001 = rho001 + rhoStrat(k + 1)

                  rhow = rhow + rhoStratTilde(k)
                end if
              end if

              if(TestCase == "baroclinic_LC") then
                rho000_e = var_env(i, j, k, 1)
                rho001_e = var_env(i, j, k + 1, 1)

                rhow_e = 0.5 * (var_env(i, j, k, 1) + var_env(i, j, k + 1, 1))

                if(fluctuationMode) then
                  if(topography) then
                    ! TFC FJ
                    rho000_e = rho000_e + rhoStratTFC(i, j, k)
                    rho001_e = rho001_e + rhoStratTFC(i, j, k + 1)
                    rhow_e = rhow_e + 0.5 * (rhoStratTFC(i, j, k) &
                        + rhoStratTFC(i, j, k + 1))
                  else
                    rho000_e = rho000_e + rhoStrat_0(k)
                    rho001_e = rho001_e + rhoStrat_0(k + 1)

                    rhow_e = rhow_e + 0.5 * (rhoStrat_0(k) + rhoStrat_0(k + 1))
                  end if
                end if
              end if

              !--- pressure gradient term -> piGrad
              if(TestCase == "baroclinic_LC") then
                piU = var(i, j, k + 1, 5) - var_env(i, j, k + 1, 5)
                piD = var(i, j, k, 5) - var_env(i, j, k, 5)

                piU_e = var_env(i, j, k + 1, 5)
                piD_e = var_env(i, j, k, 5)
              else
                piU = var(i, j, k + 1, 5)
                piD = var(i, j, k, 5)
              end if

              if(topography) then
                ! TFC FJ
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) - var_env(:, :, :, 5)
                end if
                ! Compute values at cell edges.
                pEdgeU = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    + pStratTFC(i, j, k + 1) / jac(i, j, k + 1))
                piREdgeU = 0.5 * (jac(i + 1, j, k) * met(i + 1, j, k, 3, 1) &
                    * var(i + 1, j, k, 5) + jac(i + 1, j, k + 1) * met(i + 1, &
                    j, k + 1, 3, 1) * var(i + 1, j, k + 1, 5))
                piLEdgeU = 0.5 * (jac(i - 1, j, k) * met(i - 1, j, k, 3, 1) &
                    * var(i - 1, j, k, 5) + jac(i - 1, j, k + 1) * met(i - 1, &
                    j, k + 1, 3, 1) * var(i - 1, j, k + 1, 5))
                piFEdgeU = 0.5 * (jac(i, j + 1, k) * met(i, j + 1, k, 3, 2) &
                    * var(i, j + 1, k, 5) + jac(i, j + 1, k + 1) * met(i, j &
                    + 1, k + 1, 3, 2) * var(i, j + 1, k + 1, 5))
                piBEdgeU = 0.5 * (jac(i, j - 1, k) * met(i, j - 1, k, 3, 2) &
                    * var(i, j - 1, k, 5) + jac(i, j - 1, k + 1) * met(i, j &
                    - 1, k + 1, 3, 2) * var(i, j - 1, k + 1, 5))
                chris11EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                    1) * var(i, j, k, 5) + pStratTFC(i, j, k + 1) * chris(i, &
                    j, k + 1, 1, 1) * var(i, j, k + 1, 5))
                chris22EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                    2) * var(i, j, k, 5) + pStratTFC(i, j, k + 1) * chris(i, &
                    j, k + 1, 2, 2) * var(i, j, k + 1, 5))
                chris13EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                    3) * met(i, j, k, 1, 3) * var(i, j, k, 5) + pStratTFC(i, &
                    j, k + 1) * chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, 1, &
                    3) * var(i, j, k + 1, 5))
                chris23EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                    3) * met(i, j, k, 2, 3) * var(i, j, k, 5) + pStratTFC(i, &
                    j, k + 1) * chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, 2, &
                    3) * var(i, j, k + 1, 5))
                ! Compute pressure gradient component.
                piGrad = kappaInv * MaInv2 * pEdgeU / rhow * ((piREdgeU &
                    - piLEdgeU) * 0.5 / dx + (piFEdgeU - piBEdgeU) * 0.5 / dy &
                    + (jac(i, j, k + 1) * met(i, j, k + 1, 3, 3) * var(i, j, k &
                    + 1, 5) - jac(i, j, k) * met(i, j, k, 3, 3) * var(i, j, k, &
                    5)) / dz) + kappaInv * MaInv2 / rhow * (chris11EdgeU &
                    + chris22EdgeU + 2.0 * (chris13EdgeU + chris23EdgeU))
                if(testCase == "baroclinic_LC") then
                  var(:, :, :, 5) = var(:, :, :, 5) + var_env(:, :, :, 5)
                end if
              else
                piGrad = kappaInv * MaInv2 * pstw / rhow * (piU - piD) / dz
              end if

              if(TestCase == "baroclinic_LC") then !FS
                if(topography) then
                  ! TFC FJ
                  ! Compute values at cell edges.
                  pEdgeU = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                      + pStratTFC(i, j, k + 1) / jac(i, j, k + 1))
                  piREdgeU = 0.5 * (jac(i + 1, j, k) * met(i + 1, j, k, 3, 1) &
                      * var_env(i + 1, j, k, 5) + jac(i + 1, j, k + 1) * met(i &
                      + 1, j, k + 1, 3, 1) * var_env(i + 1, j, k + 1, 5))
                  piLEdgeU = 0.5 * (jac(i - 1, j, k) * met(i - 1, j, k, 3, 1) &
                      * var_env(i - 1, j, k, 5) + jac(i - 1, j, k + 1) * met(i &
                      - 1, j, k + 1, 3, 1) * var_env(i - 1, j, k + 1, 5))
                  piFEdgeU = 0.5 * (jac(i, j + 1, k) * met(i, j + 1, k, 3, 2) &
                      * var_env(i, j + 1, k, 5) + jac(i, j + 1, k + 1) &
                      * met(i, j + 1, k + 1, 3, 2) * var_env(i, j + 1, k + 1, &
                      5))
                  piBEdgeU = 0.5 * (jac(i, j - 1, k) * met(i, j - 1, k, 3, 2) &
                      * var_env(i, j - 1, k, 5) + jac(i, j - 1, k + 1) &
                      * met(i, j - 1, k + 1, 3, 2) * var_env(i, j - 1, k + 1, &
                      5))
                  chris11EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                      1) * var_env(i, j, k, 5) + pStratTFC(i, j, k + 1) &
                      * chris(i, j, k + 1, 1, 1) * var_env(i, j, k + 1, 5))
                  chris22EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                      2) * var_env(i, j, k, 5) + pStratTFC(i, j, k + 1) &
                      * chris(i, j, k + 1, 2, 2) * var_env(i, j, k + 1, 5))
                  chris13EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                      3) * met(i, j, k, 1, 3) * var_env(i, j, k, 5) &
                      + pStratTFC(i, j, k + 1) * chris(i, j, k + 1, 1, 3) &
                      * met(i, j, k + 1, 1, 3) * var_env(i, j, k + 1, 5))
                  chris23EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                      3) * met(i, j, k, 2, 3) * var_env(i, j, k, 5) &
                      + pStratTFC(i, j, k + 1) * chris(i, j, k + 1, 2, 3) &
                      * met(i, j, k + 1, 2, 3) * var_env(i, j, k + 1, 5))
                  ! Compute pressure gradient component.
                  piGrad = piGrad + kappaInv * MaInv2 * (pEdgeU / rhow &
                      - pEdgeU / rhow_e) * ((piREdgeU - piLEdgeU) * 0.5 / dx &
                      + (piFEdgeU - piBEdgeU) * 0.5 / dy + (jac(i, j, k + 1) &
                      * met(i, j, k + 1, 3, 3) * var_env(i, j, k + 1, 5) &
                      - jac(i, j, k) * met(i, j, k, 3, 3) * var_env(i, j, k, &
                      5)) / dz) + kappaInv * MaInv2 * (1.0 / rhow + 1.0 &
                      / rhow_e) * (chris11EdgeU + chris22EdgeU + 2.0 &
                      * (chris13EdgeU + chris23EdgeU))
                else
                  piGrad = piGrad + kappaInv * MaInv2 * (pstw / rhow - pstw_e &
                      / rhow_e) * (piU_e - piD_e) / dz
                end if
              end if

              ! wstar
              wvert = var(i, j, k, 4)

              ! squared Brunt-Vaisala frequency averaged to half
              ! levels
              ! (could be done a bit nicer by determining this without
              ! averaging directly from the reference-atmosphere
              ! density)
              if(topography) then
                ! TFC FJ
                bvsstw = 0.5 * (bvsStratTFC(i, j, k) + bvsStratTFC(i, j, k + 1))
              else
                bvsstw = 0.5 * (bvsStrat(k) + bvsStrat(k + 1))
              end if

              facw = 1.0

              ! if (topography) then
              !    ! Rayleigh damping for topography (immersed boundary)
              !
              !    if(k < kbl_topo(i,j,3)) then
              !       facw = facw + alprlx*dt
              !      else if(k == kbl_topo(i,j,3)) then
              !        stop'implementation topography into &
              !            &semi-implicit time step still to be done'
              !    end if
              ! end if

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping

                  facw = facw + dt * 0.5 * (kw_hs(k) + kw_hs(k + 1))
                end if
              end if

              if(spongeLayer) then
                if(topography .and. spongeTFC) then
                  ! TFC FJ
                  facw = facw + dt * 0.5 * (alphaTFC(i, j, k) + alphaTFC(i, j, &
                      k + 1))
                else
                  facw = facw + dt * 0.5 * (kr_sp_w(j, k) + kr_sp_w(j, k + 1))
                end if
              end if

              !heat0 &
              != heat(i,j,k) - S_bar(k) &
              !  - Pstrat(k)/g_ndim * bvsStrat(k) &
              !    * 0.5*(w_0(k) + w_0(k-1))
              heat0 = heat(i, j, k)

              !heat1 &
              != heat(i,j,k+1) - S_bar(k+1) &
              !  - Pstrat(k+1)/g_ndim * bvsStrat(k+1) &
              !    * 0.5*(w_0(k+1) + w_0(k))
              heat1 = heat(i, j, k + 1)

              ! TFC FJ
              ! Buoyancy is predicted after momentum in implicit steps.
              if(topography) then
                if(testCase == "baroclinic_LC") then
                  buoy = - g_ndim * 0.5 * (var(i, j, k, 6) / rho000 / jac(i, &
                      j, k) - var_env(i, j, k, 6) / rho000_e / jac(i, j, k) &
                      + var(i, j, k + 1, 6) / rho001 / jac(i, j, k + 1) &
                      - var_env(i, j, k + 1, 6) / rho001_e / jac(i, j, k + 1))
                else
                  buoy = - g_ndim * 0.5 * (var(i, j, k, 6) / rho000 / jac(i, &
                      j, k) + var(i, j, k + 1, 6) / rho001 / jac(i, j, k + 1))
                end if
              end if

              if(topography) then
                ! TFC FJ
                uC = 0.5 * (var(i, j, k, 2) + var(i - 1, j, k, 2))
                uU = 0.5 * (var(i, j, k + 1, 2) + var(i - 1, j, k + 1, 2))
                vC = 0.5 * (var(i, j, k, 3) + var(i, j - 1, k, 3))
                vU = 0.5 * (var(i, j, k + 1, 3) + var(i, j - 1, k + 1, 3))
                wAst = 1.0 / (facw + rhoStratEdgeU / rhow * bvsstw * dt &
                    ** 2.0) * (wvert - dt * piGrad + dt * buoy + rhoStratEdgeU &
                    / rhow * bvsstw * dt ** 2.0 * (0.5 * (met(i, j, k, 1, 3) &
                    * uC + met(i, j, k + 1, 1, 3) * uU) + 0.5 * (met(i, j, k, &
                    2, 3) * vC + met(i, j, k + 1, 2, 3) * vU)))
              else
                if(TestCase == "baroclinic_LC") then
                  wAst = 1.0 / (facw + rhoStratTilde(k) / rhow * pstw / pstw_0 &
                      * bvsstw * dt ** 2) * (wvert - dt * piGrad - dt * g_ndim &
                      * 0.5 * (rhopOld(i, j, k) / rho000 - var_env(i, j, k, 6) &
                      / rho000_e + rhopOld(i, j, k + 1) / rho001 - var_env(i, &
                      j, k + 1, 6) / rho001_e + dt * (rhoStrat(k) &
                      / Pstrat_0(k) * heat0 / rho000 + rhoStrat(k + 1) &
                      / Pstrat_0(k + 1) * heat1 / rho001)))
                  !/(  facw &
                  !  + rhoStratTilde(k)/rhow * bvsstw * dt**2) &
                  !* (  rhoStrat(k)/Pstrat(k) &
                  !+ rhoStrat(k+1)/Pstrat(k+1) &
                else
                  wAst = 1.0 / (facw + rhoStratTilde(k) / rhow * pstw / pstw_0 &
                      * bvsstw * dt ** 2) * (wvert - dt * piGrad - dt * g_ndim &
                      * 0.5 * (rhopOld(i, j, k) / rho000 + rhopOld(i, j, k &
                      + 1) / rho001 + dt * (rhoStrat(k) / Pstrat_0(k) * heat0 &
                      / rho000 + rhoStrat(k + 1) / Pstrat_0(k + 1) * heat1 &
                      / rho001)))
                  !/(  facw &
                  !  + rhoStratTilde(k)/rhow * bvsstw * dt**2) &
                  !* (  rhoStrat(k)/Pstrat(k) &
                  !+ rhoStrat(k+1)/Pstrat(k+1) &
                end if
              end if

              var(i, j, k, 4) = wAst
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
      if(int_mod == 'expl' .and. .not. spongeTFC) then
        ! TFC FJ
        spongeLayer = spongeLayer_s
        ! topography = topography_s
      else if(int_mod == 'impl') then
        kr_sp = kr_sp / facray
        kr_sp_w = kr_sp_w / facray
        alprlx = alprlx / facray
      end if
    end if

    !end subroutine momentumPredictor_wc
  end subroutine momentumPredictor

  !-------------------------------------------------------------------------

  subroutine thetaUpdate(var, var0, flux, source, dt, q, m)
    !-----------------------------
    ! adds theta flux to cell theta
    !-----------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(in) :: var0

    real, dimension(- 1:nx, - 1:ny, - 1:nz, 3, nVar), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux)
    ! dir = 1..3 > f-, g- and h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW, fTheta

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(in) :: source

    real, intent(in) :: dt
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        intent(inout) :: q

    integer, intent(in) :: m

    ! local variables
    integer :: i, j, k, l
    real :: fL, fR ! flux Left/Right
    real :: gB, gF ! flux Backward/Forward
    real :: hD, hU ! flux Downward/Upward
    real :: fluxDiff ! convective part
    real :: F ! F(phi)

    ! advection of background
    real :: u, v, w, w_true
    real :: adv

    ! init q
    if(m == 1) q = 0.

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

          fL = flux(i - 1, j, k, 1, 6) ! theta flux accros left cell edge
          fR = flux(i, j, k, 1, 6) ! right
          gB = flux(i, j - 1, k, 2, 6) ! backward
          gF = flux(i, j, k, 2, 6) ! forward
          hD = flux(i, j, k - 1, 3, 6) ! downward
          hU = flux(i, j, k, 3, 6) ! upward

          ! convective part
          fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

          ! advective part of background stratification
          u = 0.5 * (var(i, j, k, 2) + var(i - 1, j, k, 2))
          v = 0.5 * (var(i, j, k, 3) + var(i, j - 1, k, 3))
          w = 0.5 * (var(i, j, k, 4) + var(i, j, k - 1, 4))
          w_true = u * vertical(1) + v * vertical(2) + w * vertical(3)

          adv = w_true * Fr2 * theta00 * N2

          ! diffusive part
          ! diff = ....

          ! F(phi)
          F = - fluxDiff - adv + source(i, j, k, 6)

          select case(timeSchemeType)

          case("lowStorage")

            ! update: q(m-1) -> q(m)
            q(i, j, k) = dt * F + alpha(m) * q(i, j, k)

            ! update potential temperature
            var(i, j, k, 6) = var(i, j, k, 6) + beta(m) * q(i, j, k)

          case("classical")

            var(i, j, k, 6) = rk(1, m) * var0(i, j, k, 6) + rk(2, m) * var(i, &
                j, k, 6) + rk(3, m) * dt * F

          case default
            stop "thetaUpdate: unknown case timeSchemeType"
          end select

        end do
      end do
    end do

    if(verbose .and. master) print *, "update.f90/thetaUpdate: theta(m=", m, &
        ") calculated."

  end subroutine thetaUpdate

  !--------------------------------------------------------------------------

  !UAC subroutine massUpdate (var,flux,dt,q,m,upd_var,upd_mod,int_mod)
  subroutine massUpdate(var, flux, dt, q, m, upd_var, upd_mod, int_mod, facray)
    !-----------------------------
    ! adds mass flux to cell mass
    !-----------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var

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

    real, dimension(- 1:nx, - 1:ny, - 1:nz, 3, nVar), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux)
    ! dir = 1..3 > f-, g- and h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW

    !UAC real, intent(in) :: dt
    real, intent(in) :: dt, facray
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        intent(inout) :: q

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
        pstwm, piU, piD, piGrad

    ! TFC FJ
    real :: pEdgeU, pEdgeD
    real :: piREdgeU, piLEdgeU, piFEdgeU, piBEdgeU, piREdgeD, piLEdgeD, &
        piFEdgeD, piBEdgeD
    real :: chris11EdgeU, chris11EdgeD, chris22EdgeU, chris22EdgeD, &
        chris13EdgeU, chris13EdgeD, chris23EdgeU, chris23EdgeD
    real :: piGradZEdgeU, piGradZEdgeD

    real :: rho_p

    real, dimension(- nbz:nz + nbz) :: w_0
    real, dimension(- nbz:nz + nbz) :: S_bar
    real :: heat_flc

    !UAB
    real :: rho_e, pstw_e, pstwm_e, rhow_e, rhowm_e
    !UAE

    real, dimension(1:nz) :: sum_local, sum_global

    real, dimension(- nbz:nz + nbz) :: rhopw_bar

    real :: ymax, yloc

    ymax = ly_dim(1) / lRef

    if(correctDivError) then
      print *, 'ERROR: correction divergence error not allowed'
      stop
    end if

    ! TFC FJ
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
            fL = flux(i - 1, j, k, 1, 1) ! mass flux accros left cell edge
            fR = flux(i, j, k, 1, 1) ! right
            gB = flux(i, j - 1, k, 2, 1) ! backward
            gF = flux(i, j, k, 2, 1) ! forward
            hD = flux(i, j, k - 1, 3, 1) ! downward
            hU = flux(i, j, k, 3, 1) ! upward

            ! convective part
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

            ! TFC FJ
            ! Adjust mass flux divergence.
            if(topography) then
              fluxDiff = fluxDiff / jac(i, j, k)
            end if

            ! F(phi)
            F = - fluxDiff

            !UAB
            ! density relaxation
            if(dens_relax) then
              if(background /= "HeldSuarez") then
                stop 'ERROR: density relaxation only ready for background &
                    = HeldSuarez'
              end if

              if(fluctuationMode) then
                rho = var(i, j, k, 1) + rhoStrat(k)
              else
                rho = var(i, j, k, 1)
              end if

              rho_e = Pstrat(k) / the_env_pp(i, j, k)

              F = F - kt_hs(j, k) * (rho - rho_e)
            end if
            !UAE

            ! update: q(m-1) -> q(m)
            q(i, j, k) = dt * F + alpha(m) * q(i, j, k)

            ! update density
            var(i, j, k, 1) = var(i, j, k, 1) + beta(m) * q(i, j, k)
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
              fL = flux(i - 1, j, k, 1, 6) ! mass flux accros left cell edge
              fR = flux(i, j, k, 1, 6) ! right
              gB = flux(i, j - 1, k, 2, 6) ! backward
              gF = flux(i, j, k, 2, 6) ! forward
              !hD = flux(i,j,k-1,3,6) - rhopw_bar(k-1) ! downward
              !hU = flux(i,j,k,3,6) - rhopw_bar(k)   ! upward
              hD = flux(i, j, k - 1, 3, 6) ! downward
              hU = flux(i, j, k, 3, 6) ! upward

              ! convective part
              fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

              ! TFC FJ
              if(topography) then
                fluxDiff = fluxDiff / jac(i, j, k)
              end if

              rhop = var(i, j, k, 6)

              rho = var(i, j, k, 1)
              if(fluctuationMode) then
                rho = rho + rhoStrat(k)
              end if

              if(topography) then
                ! TFC FJ
                wvrt = 0.5 * (vertWindTFC(i, j, k, var) + vertWindTFC(i, j, k &
                    - 1, var))
              else
                !wvrt &
                != 0.5 * (var(i,j,k,4)-w_0(k) + var(i,j,k-1,4)-w_0(k-1))
                wvrt = 0.5 * (var(i, j, k, 4) + var(i, j, k - 1, 4))
              end if

              !heat_flc= heat(i,j,k) - S_bar(k)
              heat_flc = heat(i, j, k)

              if(topography) then
                ! TFC FJ
                F = - fluxDiff + rhoStratTFC(i, j, k) / g_ndim &
                    * bvsStratTFC(i, j, k) * wvrt
              else
                ! F(phi)
                !F &
                != - fluxDiff + rhoStrat(k)/g_ndim * bvsStrat(k)*wvrt &
                !  + rhoStrat(k)/Pstrat(k) * heat_flc&
                !  - alprlx * (rhop - rho + rhoStrat(k))
                F = - fluxDiff + PStrat(k) / PStrat_0(k) * rhoStrat(k) &
                    / g_ndim * bvsStrat(k) * wvrt + rhoStrat(k) / Pstrat_0(k) &
                    * heat_flc
              end if

              ! density relaxation
              if(dens_relax) then
                if(background /= "HeldSuarez") then
                  stop 'ERROR: density relaxation only ready for background &
                      = HeldSuarez'
                end if

                if(fluctuationMode) then
                  rho = var(i, j, k, 1) + rhoStrat(k)
                else
                  rho = var(i, j, k, 1)
                end if

                !rho_e = Pstrat(k)/the_env_pp(i,j,k)
                rho_e = Pstrat_0(k) / the_env_pp(i, j, k)

                F = F - kt_hs(j, k) * (rho - rho_e)
              end if

              ! update: q(m-1) -> q(m)
              q(i, j, k) = dt * F + alpha(m) * q(i, j, k)

              ! update density
              var(i, j, k, 6) = var(i, j, k, 6) + beta(m) * q(i, j, k)
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

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              fL = flux(i - 1, j, k, 1, 6) ! mass flux accros left cell edge
              fR = flux(i, j, k, 1, 6) ! right
              gB = flux(i, j - 1, k, 2, 6) ! backward
              gF = flux(i, j, k, 2, 6) ! forward
              !hD = flux(i,j,k-1,3,6) - rhopw_bar(k-1) ! downward
              !hU = flux(i,j,k,3,6) - rhopw_bar(k)   ! upward
              hD = flux(i, j, k - 1, 3, 6) ! downward
              hU = flux(i, j, k, 3, 6) ! upward

              ! convective part
              fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

              ! TFC FJ
              if(topography) then
                fluxDiff = fluxDiff / jac(i, j, k)
              end if

              ! F(phi)
              F = - fluxDiff

              !UAB
              ! density relaxation
              if(dens_relax) then
                if(background /= "HeldSuarez") then
                  stop 'ERROR: density relaxation only ready for background &
                      = HeldSuarez'
                end if

                if(fluctuationMode) then
                  rho = var(i, j, k, 1) + rhoStrat(k)
                else
                  rho = var(i, j, k, 1)
                end if

                !rho_e = Pstrat(k)/the_env_pp(i,j,k)
                rho_e = Pstrat_0(k) / the_env_pp(i, j, k)

                F = F - kt_hs(j, k) * (rho - rho_e)
              end if
              !UAE

              ! update: q(m-1) -> q(m)
              q(i, j, k) = dt * F + alpha(m) * q(i, j, k)

              ! update density
              var(i, j, k, 6) = var(i, j, k, 6) + beta(m) * q(i, j, k)
            end do
          end do
        end do
      else if(upd_mod == "rhs") then
        ! calculate bstar ...

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

        if(int_mod == "impl") then
          ! if(topography) then
          !    i00=is+nbx-1
          !    j00=js+nby-1
          ! end if

          kr_sp = kr_sp * facray
          kr_sp_w = kr_sp_w * facray
          alprlx = alprlx * facray

          do k = 1, nz
            pstw = 0.5 * (Pstrat(k) + Pstrat(k + 1))
            pstwm = 0.5 * (Pstrat(k - 1) + Pstrat(k))

            do j = 1, ny
              do i = 1, nx
                rho = var(i, j, k, 1)
                rhow = 0.5 * (var(i, j, k, 1) + var(i, j, k + 1, 1))
                rhowm = 0.5 * (var(i, j, k - 1, 1) + var(i, j, k, 1))

                if(fluctuationMode) then
                  if(topography) then
                    ! TFC FJ
                    rho = rho + rhoStratTFC(i, j, k)
                    rhow = rhow + 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, &
                        j, k + 1))
                    rhowm = rhowm + 0.5 * (rhoStratTFC(i, j, k) &
                        + rhoStratTFC(i, j, k - 1))
                  else
                    rho = rho + rhoStrat(k)
                    rhow = rhow + rhoStratTilde(k)
                    rhowm = rhowm + rhoStratTilde(k - 1)
                  end if
                end if

                if(topography) then
                  ! TFC FJ
                  ! Momentum is predicted before buoyancy in implicit
                  ! steps.
                  wvrt = 0.5 * (wOldTFC(i, j, k) + wOldTFC(i, j, k - 1))
                else
                  wvrt = 0.5 * (var(i, j, k, 4) + var(i, j, k - 1, 4))
                end if

                !heat_flc &
                != heat(i,j,k) - S_bar(k) &
                !  - Pstrat(k)/g_ndim * bvsStrat(k) &
                !    * 0.5*(w_0(k) + w_0(k-1))
                heat_flc = heat(i, j, k)

                if(topography) then
                  ! TFC FJ
                  if(testCase == "baroclinic_LC") then
                    var(:, :, :, 5) = var(:, :, :, 5) - var_env(:, :, :, 5)
                  end if
                  ! Compute P coefficients.
                  pEdgeU = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                      + pStratTFC(i, j, k + 1) / jac(i, j, k + 1))
                  pEdgeD = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                      + pStratTFC(i, j, k - 1) / jac(i, j, k - 1))
                  ! Interpolate pressure differences.
                  piREdgeU = 0.5 * (jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) &
                      * var(i + 1, j, k, 5) + jac(i + 1, j, k + 1) * met(i &
                      + 1, j, k + 1, 1, 3) * var(i + 1, j, k + 1, 5))
                  piLEdgeU = 0.5 * (jac(i - 1, j, k) * met(i - 1, j, k, 1, 3) &
                      * var(i - 1, j, k, 5) + jac(i - 1, j, k + 1) * met(i &
                      - 1, j, k + 1, 1, 3) * var(i - 1, j, k + 1, 5))
                  piREdgeD = 0.5 * (jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) &
                      * var(i + 1, j, k, 5) + jac(i + 1, j, k - 1) * met(i &
                      + 1, j, k - 1, 1, 3) * var(i + 1, j, k - 1, 5))
                  piLEdgeD = 0.5 * (jac(i - 1, j, k) * met(i - 1, j, k, 1, 3) &
                      * var(i - 1, j, k, 5) + jac(i - 1, j, k - 1) * met(i &
                      - 1, j, k - 1, 1, 3) * var(i - 1, j, k - 1, 5))
                  piFEdgeU = 0.5 * (jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) &
                      * var(i, j + 1, k, 5) + jac(i, j + 1, k + 1) * met(i, j &
                      + 1, k + 1, 2, 3) * var(i, j + 1, k + 1, 5))
                  piBEdgeU = 0.5 * (jac(i, j - 1, k) * met(i, j - 1, k, 2, 3) &
                      * var(i, j - 1, k, 5) + jac(i, j - 1, k + 1) * met(i, j &
                      - 1, k + 1, 2, 3) * var(i, j - 1, k + 1, 5))
                  piFEdgeD = 0.5 * (jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) &
                      * var(i, j + 1, k, 5) + jac(i, j + 1, k - 1) * met(i, j &
                      + 1, k - 1, 2, 3) * var(i, j + 1, k - 1, 5))
                  piBEdgeD = 0.5 * (jac(i, j - 1, k) * met(i, j - 1, k, 2, 3) &
                      * var(i, j - 1, k, 5) + jac(i, j - 1, k - 1) * met(i, j &
                      - 1, k - 1, 2, 3) * var(i, j - 1, k - 1, 5))
                  ! Interpolate Christoffel symbols.
                  chris11EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                      1) * var(i, j, k, 5) + pStratTFC(i, j, k + 1) * chris(i, &
                      j, k + 1, 1, 1) * var(i, j, k + 1, 5))
                  chris11EdgeD = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                      1) * var(i, j, k, 5) + pStratTFC(i, j, k - 1) * chris(i, &
                      j, k - 1, 1, 1) * var(i, j, k - 1, 5))
                  chris22EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                      2) * var(i, j, k, 5) + pStratTFC(i, j, k + 1) * chris(i, &
                      j, k + 1, 2, 2) * var(i, j, k + 1, 5))
                  chris22EdgeD = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                      2) * var(i, j, k, 5) + pStratTFC(i, j, k - 1) * chris(i, &
                      j, k - 1, 2, 2) * var(i, j, k - 1, 5))
                  chris13EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                      3) * met(i, j, k, 1, 3) * var(i, j, k, 5) + pStratTFC(i, &
                      j, k + 1) * chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, &
                      1, 3) * var(i, j, k + 1, 5))
                  chris13EdgeD = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                      3) * met(i, j, k, 1, 3) * var(i, j, k, 5) + pStratTFC(i, &
                      j, k - 1) * chris(i, j, k - 1, 1, 3) * met(i, j, k - 1, &
                      1, 3) * var(i, j, k - 1, 5))
                  chris23EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                      3) * met(i, j, k, 2, 3) * var(i, j, k, 5) + pStratTFC(i, &
                      j, k + 1) * chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, &
                      2, 3) * var(i, j, k + 1, 5))
                  chris23EdgeD = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                      3) * met(i, j, k, 2, 3) * var(i, j, k, 5) + pStratTFC(i, &
                      j, k - 1) * chris(i, j, k - 1, 2, 3) * met(i, j, k - 1, &
                      2, 3) * var(i, j, k - 1, 5))
                  ! Compute pressure gradients.
                  piGradZEdgeU = kappaInv * MaInv2 * pEdgeU / rhow * (0.5 &
                      * (piREdgeU - piLEdgeU) / dx + 0.5 * (piFEdgeU &
                      - piBEdgeU) / dy + (jac(i, j, k + 1) * met(i, j, k + 1, &
                      3, 3) * var(i, j, k + 1, 5) - jac(i, j, k) * met(i, j, &
                      k, 3, 3) * var(i, j, k, 5)) / dz) + kappaInv * MaInv2 &
                      / rhow * (chris11EdgeU + chris22EdgeU + 2.0 &
                      * (chris13EdgeU + chris23EdgeU))
                  piGradZEdgeD = kappaInv * MaInv2 * pEdgeD / rhowm * (0.5 &
                      * (piREdgeD - piLEdgeD) / dx + 0.5 * (piFEdgeD &
                      - piBEdgeD) / dy + (jac(i, j, k) * met(i, j, k, 3, 3) &
                      * var(i, j, k, 5) - jac(i, j, k - 1) * met(i, j, k - 1, &
                      3, 3) * var(i, j, k - 1, 5)) / dz) + kappaInv * MaInv2 &
                      / rhowm * (chris11EdgeD + chris22EdgeD + 2.0 &
                      * (chris13EdgeD + chris23EdgeD))
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
                    var(:, :, :, 5) = var(:, :, :, 5) + var_env(:, :, :, 5)
                    ! Interpolate pressure differences.
                    piREdgeU = 0.5 * (jac(i + 1, j, k) * met(i + 1, j, k, 1, &
                        3) * var_env(i + 1, j, k, 5) + jac(i + 1, j, k + 1) &
                        * met(i + 1, j, k + 1, 1, 3) * var_env(i + 1, j, k &
                        + 1, 5))
                    piLEdgeU = 0.5 * (jac(i - 1, j, k) * met(i - 1, j, k, 1, &
                        3) * var_env(i - 1, j, k, 5) + jac(i - 1, j, k + 1) &
                        * met(i - 1, j, k + 1, 1, 3) * var_env(i - 1, j, k &
                        + 1, 5))
                    piREdgeD = 0.5 * (jac(i + 1, j, k) * met(i + 1, j, k, 1, &
                        3) * var_env(i + 1, j, k, 5) + jac(i + 1, j, k - 1) &
                        * met(i + 1, j, k - 1, 1, 3) * var_env(i + 1, j, k &
                        - 1, 5))
                    piLEdgeD = 0.5 * (jac(i - 1, j, k) * met(i - 1, j, k, 1, &
                        3) * var_env(i - 1, j, k, 5) + jac(i - 1, j, k - 1) &
                        * met(i - 1, j, k - 1, 1, 3) * var_env(i - 1, j, k &
                        - 1, 5))
                    piFEdgeU = 0.5 * (jac(i, j + 1, k) * met(i, j + 1, k, 2, &
                        3) * var_env(i, j + 1, k, 5) + jac(i, j + 1, k + 1) &
                        * met(i, j + 1, k + 1, 2, 3) * var_env(i, j + 1, k &
                        + 1, 5))
                    piBEdgeU = 0.5 * (jac(i, j - 1, k) * met(i, j - 1, k, 2, &
                        3) * var_env(i, j - 1, k, 5) + jac(i, j - 1, k + 1) &
                        * met(i, j - 1, k + 1, 2, 3) * var_env(i, j - 1, k &
                        + 1, 5))
                    piFEdgeD = 0.5 * (jac(i, j + 1, k) * met(i, j + 1, k, 2, &
                        3) * var_env(i, j + 1, k, 5) + jac(i, j + 1, k - 1) &
                        * met(i, j + 1, k - 1, 2, 3) * var_env(i, j + 1, k &
                        - 1, 5))
                    piBEdgeD = 0.5 * (jac(i, j - 1, k) * met(i, j - 1, k, 2, &
                        3) * var_env(i, j - 1, k, 5) + jac(i, j - 1, k - 1) &
                        * met(i, j - 1, k - 1, 2, 3) * var_env(i, j - 1, k &
                        - 1, 5))
                    ! Interpolate Christoffel symbols.
                    chris11EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, &
                        1, 1) * var_env(i, j, k, 5) + pStratTFC(i, j, k + 1) &
                        * chris(i, j, k + 1, 1, 1) * var_env(i, j, k + 1, 5))
                    chris11EdgeD = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, &
                        1, 1) * var_env(i, j, k, 5) + pStratTFC(i, j, k - 1) &
                        * chris(i, j, k - 1, 1, 1) * var_env(i, j, k - 1, 5))
                    chris22EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, &
                        2, 2) * var_env(i, j, k, 5) + pStratTFC(i, j, k + 1) &
                        * chris(i, j, k + 1, 2, 2) * var_env(i, j, k + 1, 5))
                    chris22EdgeD = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, &
                        2, 2) * var_env(i, j, k, 5) + pStratTFC(i, j, k - 1) &
                        * chris(i, j, k - 1, 2, 2) * var_env(i, j, k - 1, 5))
                    chris13EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, &
                        1, 3) * met(i, j, k, 1, 3) * var_env(i, j, k, 5) &
                        + pStratTFC(i, j, k + 1) * chris(i, j, k + 1, 1, 3) &
                        * met(i, j, k + 1, 1, 3) * var_env(i, j, k + 1, 5))
                    chris13EdgeD = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, &
                        1, 3) * met(i, j, k, 1, 3) * var_env(i, j, k, 5) &
                        + pStratTFC(i, j, k - 1) * chris(i, j, k - 1, 1, 3) &
                        * met(i, j, k - 1, 1, 3) * var_env(i, j, k - 1, 5))
                    chris23EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, &
                        2, 3) * met(i, j, k, 2, 3) * var_env(i, j, k, 5) &
                        + pStratTFC(i, j, k + 1) * chris(i, j, k + 1, 2, 3) &
                        * met(i, j, k + 1, 2, 3) * var_env(i, j, k + 1, 5))
                    chris23EdgeD = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, &
                        2, 3) * met(i, j, k, 2, 3) * var_env(i, j, k, 5) &
                        + pStratTFC(i, j, k - 1) * chris(i, j, k - 1, 2, 3) &
                        * met(i, j, k - 1, 2, 3) * var_env(i, j, k - 1, 5))
                    ! Compute pressure gradients.
                    piGradZEdgeU = kappaInv * MaInv2 * (pEdgeU / rhow - pEdgeU &
                        / rhow_e) * (0.5 * (piREdgeU - piLEdgeU) / dx + 0.5 &
                        * (piFEdgeU - piBEdgeU) / dy + (jac(i, j, k + 1) &
                        * met(i, j, k + 1, 3, 3) * var_env(i, j, k + 1, 5) &
                        - jac(i, j, k) * met(i, j, k, 3, 3) * var_env(i, j, k, &
                        5)) / dz) + kappaInv * MaInv2 * (1.0 / rhow - 1.0 &
                        / rhow_e) * (chris11EdgeU + chris22EdgeU + 2.0 &
                        * (chris13EdgeU + chris23EdgeU))
                    piGradZEdgeD = kappaInv * MaInv2 * (pEdgeD / rhowm &
                        - pEdgeD / rhowm_e) * (0.5 * (piREdgeD - piLEdgeD) &
                        / dx + 0.5 * (piFEdgeD - piBEdgeD) / dy + (jac(i, j, &
                        k) * met(i, j, k, 3, 3) * var_env(i, j, k, 5) - jac(i, &
                        j, k - 1) * met(i, j, k - 1, 3, 3) * var_env(i, j, k &
                        - 1, 5)) / dz) + kappaInv * MaInv2 * (1.0 / rhowm &
                        - 1.0 / rhowm_e) * (chris11EdgeD + chris22EdgeD + 2.0 &
                        * (chris13EdgeD + chris23EdgeD))
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
                    piGrad = kappaInv * MaInv2 * 0.5 * (pstw / rhow * (var(i, &
                        j, k + 1, 5) - var(i, j, k, 5) - var_env(i, j, k + 1, &
                        5) + var_env(i, j, k, 5)) / dz + pstwm / rhowm &
                        * (var(i, j, k, 5) - var(i, j, k - 1, 5) - var_env(i, &
                        j, k, 5) + var_env(i, j, k - 1, 5)) / dz)

                    pstw_e = 0.5 * (pStrat_0(k) + pStrat_0(k + 1))
                    pstwm_e = 0.5 * (pStrat_0(k - 1) + pStrat_0(k))

                    rhow_e = 0.5 * (var_env(i, j, k, 1) + var_env(i, j, k + 1, &
                        1))
                    rhowm_e = 0.5 * (var_env(i, j, k - 1, 1) + var_env(i, j, &
                        k, 1))

                    if(fluctuationMode) then
                      rhow_e = rhow_e + 0.5 * (rhoStrat_0(k) + rhoStrat_0(k &
                          + 1))
                      rhowm_e = rhowm_e + 0.5 * (rhoStrat_0(k - 1) &
                          + rhoStrat_0(k))
                    end if

                    piGrad = piGrad + kappaInv * MaInv2 * 0.5 * ((pstw / rhow &
                        - pstw_e / rhow_e) * (var_env(i, j, k + 1, 5) &
                        - var_env(i, j, k, 5)) / dz + (pstwm / rhowm - pstwm_e &
                        / rhowm_e) * (var_env(i, j, k, 5) - var_env(i, j, k &
                        - 1, 5)) / dz)
                  else
                    piGrad = kappaInv * MaInv2 * 0.5 * (pstw / rhow * (var(i, &
                        j, k + 1, 5) - var(i, j, k, 5)) / dz + pstwm / rhowm &
                        * (var(i, j, k, 5) - var(i, j, k - 1, 5)) / dz)
                  end if
                end if

                ! due to damping of wind in land cells (if there is
                ! topography)
                facw = 1.0

                ! if(topography) then
                !    !UAC if(topography_mask(i00+i,j00+j,k)&
                !    !   .or.&
                !    !   topography_mask(i00+i,j00+j,k+1)) then
                !    !   facw = facw + alprlx*dt
                !    !end if
                !    stop'implementation topography into &
                !        &semi-implicit time step still to be done'
                !    !UAE
                ! end if

                if(TestCase == "baroclinic_LC") then
                  if(background == "HeldSuarez") then
                    ! Rayleigh damping

                    facw = facw + dt * kw_hs(k)
                  end if
                end if

                if(spongeLayer) then
                  if(topography .and. spongeTFC) then
                    ! TFC FJ
                    facw = facw + dt * alphaTFC(i, j, k)
                  else
                    facw = facw + dt * kr_sp_w(j, k)
                  end if
                end if

                if(TestCase == "baroclinic_LC") then
                  rho_e = var_env(i, j, k, 1)

                  if(fluctuationMode) then
                    if(topography) then
                      ! TFC FJ
                      rho_e = rho_e + rhoStratTFC(i, j, k)
                    else
                      rho_e = rho_e + rhoStrat_0(k)
                    end if
                  end if

                  if(topography) then
                    ! TFC FJ
                    ! Predict buoyancy.
                    buoy = - g_ndim * (var(i, j, k, 6) / rho - var_env(i, j, &
                        k, 6) / rho_e)
                    buoy = 1.0 / (facw + rhoStratTFC(i, j, k) / rho &
                        * bvsStratTFC(i, j, k) * dt ** 2.0) * (- &
                        rhoStratTFC(i, j, k) / rho * bvsStratTFC(i, j, k) * dt &
                        * jac(i, j, k) * (wvrt - dt * piGrad) + facw * buoy &
                        + rhoStratTFC(i, j, k) / rho * bvsStratTFC(i, j, k) &
                        * dt * jac(i, j, k) * facw * 0.5 * (met(i, j, k, 1, 3) &
                        * (var(i, j, k, 2) + var(i - 1, j, k, 2)) + met(i, j, &
                        k, 2, 3) * (var(i, j, k, 3) + var(i, j - 1, k, 3))))
                  else
                    buoy = 1.0 / (facw + rhoStrat(k) / rho * PStrat(k) &
                        / PStrat_0(k) * bvsStrat(k) * dt ** 2) * (- &
                        rhoStrat(k) / rho * PStrat(k) / PStrat_0(k) &
                        * bvsStrat(k) * dt * (wvrt - dt * piGrad) - facw &
                        * g_ndim * (var(i, j, k, 6) / rho - var_env(i, j, k, &
                        6) / rho_e + dt / rho * rhoStrat(k) / Pstrat_0(k) &
                        * heat_flc))
                    !/( facw &
                    !  + rhoStrat(k)/rho * bvsStrat(k) * dt**2) &
                    !* (-rhoStrat(k)/rho * bvsStrat(k) * dt &
                    !+ dt/rho * rhoStrat(k)/Pstrat(k) &
                  end if

                  buoy = buoy - g_ndim * var_env(i, j, k, 6) / rho_e
                else
                  if(topography) then
                    ! TFC FJ
                    ! Predict buoyancy.
                    buoy = - g_ndim * var(i, j, k, 6) / rho
                    buoy = 1.0 / (facw + rhoStratTFC(i, j, k) / rho &
                        * bvsStratTFC(i, j, k) * dt ** 2.0) * (- &
                        rhoStratTFC(i, j, k) / rho * bvsStratTFC(i, j, k) * dt &
                        * jac(i, j, k) * (wvrt - dt * piGrad) + facw * buoy &
                        + rhoStratTFC(i, j, k) / rho * bvsStratTFC(i, j, k) &
                        * dt * jac(i, j, k) * facw * 0.5 * (met(i, j, k, 1, 3) &
                        * (var(i, j, k, 2) + var(i - 1, j, k, 2)) + met(i, j, &
                        k, 2, 3) * (var(i, j, k, 3) + var(i, j - 1, k, 3))))
                  else
                    buoy = 1.0 / (facw + rhoStrat(k) / rho * PStrat(k) &
                        / PStrat_0(k) * bvsStrat(k) * dt ** 2) * (- &
                        rhoStrat(k) / rho * PStrat(k) / PStrat_0(k) &
                        * bvsStrat(k) * dt * (wvrt - dt * piGrad) - facw &
                        * g_ndim / rho * (var(i, j, k, 6) + dt * rhoStrat(k) &
                        / Pstrat_0(k) * heat_flc))
                    !/( facw &
                    ! + rhoStrat(k)/rho * bvsStrat(k) * dt**2) &
                    !* (-rhoStrat(k)/rho * bvsStrat(k) * dt &
                    !+ dt *  rhoStrat(k)/Pstrat(k) &
                  end if
                end if

                var(i, j, k, 6) = - buoy * rho / g_ndim
              end do
            end do
          end do

          kr_sp = kr_sp / facray
          kr_sp_w = kr_sp_w / facray
          alprlx = alprlx / facray

        else if(int_mod == "expl") then
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                rhop = var(i, j, k, 6)

                rho = var(i, j, k, 1)
                if(fluctuationMode) then
                  if(topography) then
                    ! TFC FJ
                    rho = rho + rhoStratTFC(i, j, k)
                  else
                    rho = rho + rhoStrat(k)
                  end if
                end if

                if(topography) then
                  ! TFC FJ
                  wvrt = 0.5 * (vertWindTFC(i, j, k, var) + vertWindTFC(i, j, &
                      k - 1, var))
                else
                  !wvrt &
                  != 0.5 &
                  !  * (var(i,j,k,4)-w_0(k) + var(i,j,k-1,4)-w_0(k-1))
                  wvrt = 0.5 * (var(i, j, k, 4) + var(i, j, k - 1, 4))
                end if

                !heat_flc= heat(i,j,k) - S_bar(k)
                heat_flc = heat(i, j, k)

                if(topography) then
                  ! TFC FJ
                  buoy = - g_ndim * rhop / rho
                  buoy = buoy - dt * rhoStratTFC(i, j, k) / rho &
                      * bvsStratTFC(i, j, k) * wvrt
                else
                  buoy = - g_ndim * rhop / rho - dt * (rhoStrat(k) / rho &
                      * PStrat(k) / PStrat_0(k) * bvsStrat(k) * wvrt + g_ndim &
                      / rho * rhoStrat(k) / Pstrat_0(k) * heat_flc)
                  !- dt * (  rhoStrat(k)/rho * bvsStrat(k)*wvrt &
                  !* rhoStrat(k)/Pstrat(k) * heat_flc)
                  !UAE
                end if

                var(i, j, k, 6) = - buoy * rho / g_ndim
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

  !-----------------------------------------------------------------------
  subroutine tracerUpdate(var, flux, dt, q, m)

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var

    real, dimension(- 1:nx, - 1:ny, - 1:nz, 3, nVar), intent(in) :: flux

    real, intent(in) :: dt
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        intent(inout) :: q

    integer, intent(in) :: m
    integer :: i00, j00

    ! local variables
    integer :: i, j, k, l
    real :: fL, fR ! flux Left/Right
    real :: gB, gF ! flux Backward/Forward
    real :: hD, hU ! flux Downward/Upward
    real :: fluxDiff ! convective part
    real :: F ! F(phi)
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: rho

    if(correctDivError) then
      print *, 'ERROR: correction divergence error not allowed'
      stop
    end if

    ! init q
    if(m == 1) q = 0.

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          fL = flux(i - 1, j, k, 1, iVart) ! mass flux accros left cell edge
          fR = flux(i, j, k, 1, iVart) ! right
          gB = flux(i, j - 1, k, 2, iVart) ! backward
          gF = flux(i, j, k, 2, iVart) ! forward
          hD = flux(i, j, k - 1, 3, iVart) ! downward
          hU = flux(i, j, k, 3, iVart) ! upward

          ! convective part
          fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

          if(topography) then
            fluxDiff = fluxDiff / jac(i, j, k)
          end if

          ! F(phi)
          F = - fluxDiff

          if(dens_relax) then
            stop "update.f90: dens_relax not implemented in tracerUpdate"
          end if

          ! update: q(m-1) -> q(m)
          q(i, j, k) = dt * F + alpha(m) * q(i, j, k)

          ! update density
          var(i, j, k, iVart) = var(i, j, k, iVart) + beta(m) * q(i, j, k)

        end do
      end do
    end do

  end subroutine tracerUpdate

  !-----------------------------------------------------------------------

  subroutine iceUpdate(var, var0, flux, source, dt, q, m)
    !-----------------------------
    ! adds ice flux to cell ice
    !-----------------------------
    ! mainly analogous to massUpdate

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var0
    real, dimension(- 1:nx, - 1:ny, - 1:nz, 3, nVar), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux)
    ! dir = 1..3 > f-, g- and h-flux in x,y,z-direction
    ! iFlux = 8..11 > Rho_nAer, Rho_nIce, Rho_qIce, Rho_qv

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(in) :: source

    real, intent(in) :: dt
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 4), &
        intent(inout) :: q

    integer, intent(in) :: m

    ! local integer
    integer :: i, j, k, iVar

    ! local variables
    real, dimension(4) :: fL, fR ! flux Left/Right
    real, dimension(4) :: gB, gF ! flux Backward/Forward
    real, dimension(4) :: hD, hU ! flux Downward/Upward
    real, dimension(4) :: fluxDiff ! convective part
    real, dimension(4) :: F ! F(phi)
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: rho
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 4) :: &
        rho_source_term
    real :: T, p, SIce

    var0 = var

    ! init q
    if(m == 1) q = 0.

    if(fluctuationMode) then
      do k = - 1, nz + 1
        rho(:, :, k) = var(:, :, k, 1) + rhoStrat(k)
      end do
    else
      rho = var(:, :, :, 1)
    end if

    if(correctDivError) then
      do k = 0, 3
        rho_source_term(:, :, :, 4 - k) = var(:, :, :, nVar - k) * source(:, &
            :, :, 1)
      end do
    else
      rho_source_term(:, :, :, :) = 0.0
    end if

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

          !UAC if (topography_mask(i+is+nbx-1,j+js+nby-1,k)==.false.) then
          ! topography not used as a condition anymore. This should all
          ! be done by the winds responding to the immersed boundary
          if(k > 0) then
            !UAE

            fL = flux(i - 1, j, k, 1, nVar - 3:nVar) ! mass flux across left cell edge
            fR = flux(i, j, k, 1, nVar - 3:nVar) ! right
            gB = flux(i, j - 1, k, 2, nVar - 3:nVar) ! backward
            gF = flux(i, j, k, 2, nVar - 3:nVar) ! forward
            hD = flux(i, j, k - 1, 3, nVar - 3:nVar) ! downward
            hU = flux(i, j, k, 3, nVar - 3:nVar) ! upward

            ! convective part
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

            ! diffusive part
            ! diff = ....

            ! F(phi)
            F = - fluxDiff

            F(:) = F(:) + rho_source_term(i, j, k, :) + rho(i, j, k) &
                * source(i, j, k, nVar - 3:nVar)

            select case(timeSchemeType)

            case("lowStorage")

              ! update: q(m-1) -> q(m)
              q(i, j, k, :) = dt * F(:) + alpha(m) * q(i, j, k, :)

              ! update variables
              var(i, j, k, nVar - 3:nVar) = var(i, j, k, nVar - 3:nVar) &
                  + beta(m) * q(i, j, k, 1:4) / rho(i, j, k)

            case("classical")

              var(i, j, k, nVar - 3:nVar) = rk(1, m) * var0(i, j, k, nVar &
                  - 3:nVar) + rk(2, m) * var(i, j, k, nVar - 3:nVar) + rk(3, &
                  m) * dt * F(1:4) / rho(i, j, k)

            case default
              stop "iceUpdate: unknown case timeSchemeType"
            end select

            do iVar = nVar - 3, nVar
              ! avoid negative values for all ice variables
              if((var(i, j, k, iVar) .lt. 0.0)) then
                var(i, j, k, iVar) = 0.0
              end if
            end do

          end if
        end do
      end do
    end do

    if(verbose .and. master) print *, "update.f90/iceUpdate: ice(m=", m, ") &
        calculated."

  end subroutine iceUpdate

  !-------------------------------------------------------------------------

  subroutine ice2Update_source(var, flux, source, dt, q, m)
    !-----------------------------
    ! adds ice flux to cell ice field
    !-----------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var
    real, dimension(- 1:nx, - 1:ny, - 1:nz, 3, nVar), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux)
    ! dir = 1..3 > f-, g- and h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(in) :: source

    real, intent(in) :: dt
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVarIce), &
        intent(inout) :: q

    integer, intent(in) :: m

    ! local variables
    integer :: i, j, k, l
    real :: fL, fR ! flux Left/Right
    real :: gB, gF ! flux Backward/Forward
    real :: hD, hU ! flux Downward/Upward
    real :: fluxDiff ! convective part
    real :: F ! F(phi)

    !!$    ! TFC FJ
    !!$    real :: pEdgeU, pEdgeD
    !!$    real :: piREdgeU, piLEdgeU, piFEdgeU, piBEdgeU, &
    !!$         piREdgeD, piLEdgeD, piFEdgeD, piBEdgeD
    !!$    real :: chris11EdgeU, chris11EdgeD, chris22EdgeU, chris22EdgeD, &
    !!$         chris13EdgeU, chris13EdgeD, chris23EdgeU, chris23EdgeD
    !!$    real :: piGradZEdgeU, piGradZEdgeD

    integer :: ii, iVar

    ! init q
    if(m == 1) q = 0.

    do ii = 1, nVarIce
      iVar = iVarIce(ii)

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx

            fL = flux(i - 1, j, k, 1, iVar) ! mass flux accros left cell edge
            fR = flux(i, j, k, 1, iVar) ! right
            gB = flux(i, j - 1, k, 2, iVar) ! backward
            gF = flux(i, j, k, 2, iVar) ! forward
            hD = flux(i, j, k - 1, 3, iVar) ! downward
            hU = flux(i, j, k, 3, iVar) ! upward

            ! convective part
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

            ! TFC FJ
            ! Adjust mass flux divergence.
            if(topography) then
              fluxDiff = fluxDiff / jac(i, j, k)
            end if

            ! F(phi)
            F = - fluxDiff + source(i, j, k, iVar)

            ! update: q(m-1) -> q(m)
            q(i, j, k, ii) = dt * F + alpha(m) * q(i, j, k, ii)

            ! update fields
            var(i, j, k, iVar) = var(i, j, k, iVar) + beta(m) * q(i, j, k, ii)

          end do !i
        end do !j
      end do !k

    end do !ii

  end subroutine ice2Update_source

  !-----------------------------------------------------------------------

  subroutine ice2Update_apb(var, flux, source, dt, q, m, update_type)
    !-----------------------------
    ! adds ice flux to cell ice field
    !-----------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var
    real, dimension(- 1:nx, - 1:ny, - 1:nz, 3, nVar), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux)
    ! dir = 1..3 > f-, g- and h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(in) :: source

    real, intent(in) :: dt
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVarIce), &
        intent(inout) :: q

    integer, intent(in) :: m
    character(len = 3), intent(in) :: update_type

    ! local variables
    integer :: i, j, k, l
    real :: fL, fR ! flux Left/Right
    real :: gB, gF ! flux Backward/Forward
    real :: hD, hU ! flux Downward/Upward
    real :: fluxDiff ! convective part
    real :: F ! F(phi)

    !!$    ! TFC FJ
    !!$    real :: pEdgeU, pEdgeD
    !!$    real :: piREdgeU, piLEdgeU, piFEdgeU, piBEdgeU, &
    !!$         piREdgeD, piLEdgeD, piFEdgeD, piBEdgeD
    !!$    real :: chris11EdgeU, chris11EdgeD, chris22EdgeU, chris22EdgeD, &
    !!$         chris13EdgeU, chris13EdgeD, chris23EdgeU, chris23EdgeD
    !!$    real :: piGradZEdgeU, piGradZEdgeD

    integer :: ii, iVar

    ! init q
    if(m == 1) q = 0.

    do ii = 1, nVarIce
      iVar = iVarIce(ii)

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx

            if(update_type .eq. 'ADV' .or. update_type .eq. 'BOT') then

              fL = flux(i - 1, j, k, 1, iVar) ! mass flux accros left cell edge
              fR = flux(i, j, k, 1, iVar) ! right
              gB = flux(i, j - 1, k, 2, iVar) ! backward
              gF = flux(i, j, k, 2, iVar) ! forward
              hD = flux(i, j, k - 1, 3, iVar) ! downward
              hU = flux(i, j, k, 3, iVar) ! upward

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
              F = - fluxDiff + source(i, j, k, iVar)
            elseif(update_type .eq. 'ADV') then
              F = - fluxDiff
            elseif(update_type .eq. 'PHY') then
              F = source(i, j, k, iVar)
            else
              print *, 'wrong update_type in ice2Update_apb'
              stop

            end if

            ! update: q(m-1) -> q(m)
            q(i, j, k, ii) = dt * F + alpha(m) * q(i, j, k, ii)

            ! update fields
            var(i, j, k, iVar) = var(i, j, k, iVar) + beta(m) * q(i, j, k, ii)

          end do !i
        end do !j
      end do !k

    end do !ii

  end subroutine ice2Update_apb

  !-----------------------------------------------------------------------

  subroutine ice2Update(var, flux, dt, q, m, int_mod, facray)
    !-----------------------------
    ! adds ice flux to cell ice field
    !-----------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var

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
    character(len = *), intent(in) :: int_mod

    real, dimension(- 1:nx, - 1:ny, - 1:nz, 3, nVar), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux)
    ! dir = 1..3 > f-, g- and h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW

    !UAC real, intent(in) :: dt
    real, intent(in) :: dt, facray
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVarIce), &
        intent(inout) :: q

    integer, intent(in) :: m
    integer :: i00, j00

    ! local variables
    integer :: i, j, k, l
    real :: fL, fR ! flux Left/Right
    real :: gB, gF ! flux Backward/Forward
    real :: hD, hU ! flux Downward/Upward
    real :: fluxDiff ! convective part
    real :: F ! F(phi)

    !    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz) :: heat

    !    real :: buoy0, buoy, rho, rhow, rhowm, rhop, wvrt, facw, facr, &
    !         & pstw, pstwm, piU, piD, piGrad

    ! TFC FJ
    real :: pEdgeU, pEdgeD
    real :: piREdgeU, piLEdgeU, piFEdgeU, piBEdgeU, piREdgeD, piLEdgeD, &
        piFEdgeD, piBEdgeD
    real :: chris11EdgeU, chris11EdgeD, chris22EdgeU, chris22EdgeD, &
        chris13EdgeU, chris13EdgeD, chris23EdgeU, chris23EdgeD
    real :: piGradZEdgeU, piGradZEdgeD

    real :: rho_p

    real, dimension(- nbz:nz + nbz) :: w_0
    real, dimension(- nbz:nz + nbz) :: S_bar
    real :: heat_flc

    !UAB
    real :: rho_e, pstw_e, pstwm_e, rhow_e, rhowm_e
    !    !UAE

    !    real, dimension(1:nz) :: sum_local, sum_global

    !    real, dimension(-nbz:nz+nbz) :: rhopw_bar

    real :: ymax, yloc
    integer :: ii, iVar

    ymax = ly_dim(1) / lRef

    ! init q
    if(m == 1) q = 0.

    do ii = 1, nVarIce
      iVar = iVarIce(ii)

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx

            fL = flux(i - 1, j, k, 1, iVar) ! mass flux accros left cell edge
            fR = flux(i, j, k, 1, iVar) ! right
            gB = flux(i, j - 1, k, 2, iVar) ! backward
            gF = flux(i, j, k, 2, iVar) ! forward
            hD = flux(i, j, k - 1, 3, iVar) ! downward
            hU = flux(i, j, k, 3, iVar) ! upward

            ! convective part
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz

            ! TFC FJ
            ! Adjust mass flux divergence.
            if(topography) then
              fluxDiff = fluxDiff / jac(i, j, k)
            end if

            ! F(phi)
            F = - fluxDiff

            ! update: q(m-1) -> q(m)
            q(i, j, k, ii) = dt * F + alpha(m) * q(i, j, k, ii)

            ! update fields
            var(i, j, k, iVar) = var(i, j, k, iVar) + beta(m) * q(i, j, k, ii)

          end do !i
        end do !j
      end do !k

    end do !ii

  end subroutine ice2Update

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
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(in) :: var
    real, intent(out) :: dt
    logical, intent(out) :: errFlag

    ! locals
    real :: uMax, vMax, wMax
    real :: dtConv, dtVisc, dtCond, dtWKB
    real :: dtConv_loc, dtWKB_loc
    real :: dtMax
    real :: dtWave

    ! Buoyancy time step restriction
    real :: dtBuoy, dtBuoy_loc
    real, dimension(3) :: bMax, bMaxNew, duMax
    real :: buoyMax, buoyMin, buoyMaxNew, the_New, the_max, the_min

    ! achatzb test deletion:
    ! ! Gravity wave time stop restriction
    ! real :: lambdaMax    ! max GW length to be time resolved
    ! real :: dtWave, lambdaX, lambdaZ, kk, mm, kMin, cX, cZ
    ! achatze

    ! local integer
    integer :: i, j, k

    ! sponge layer
    real :: ReMin ! min Reynolds number in domain

    ! verbose
    logical, parameter :: giveInfo = .true.

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

      case("Boussinesq", "pseudo_incompressible")

        !----------------------------
        !   Full model time step
        !----------------------------

        !----------------------
        !     CFL condition
        !----------------------
        uMax = maxval(abs(var(1:nx, 1:ny, 1:nz, 2))) + small
        vMax = maxval(abs(var(1:nx, 1:ny, 1:nz, 3))) + small
        wMax = maxval(abs(var(1:nx, 1:ny, 1:nz, 4))) + small

        dtConv_loc = cfl * min(dx / uMax, dy / vMax, dz / wMax)

        ! find global minimum

        call mpi_reduce(dtConv_loc, dtConv, 1, mpi_double_precision, mpi_min, &
            root, comm, ierror)

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
                if(fluctuationMode) then
                  bMaxNew = abs(var(i, j, k, 1)) / (rhoStrat(k) + var(i, j, k, &
                      1)) * vertical
                else
                  bMaxNew = abs(rhoStrat(k) - var(i, j, k, 1)) / var(i, j, k, &
                      1) * vertical
                end if

              case("Boussinesq")
                ! TFC FJ
                ! Boussinesq: density fluctuations are stored in
                ! var(:, :, :, 6)!
                bMaxNew = abs(- var(i, j, k, 6)) / rho00 * vertical

                ! bMaxNew = var(i,j,k,6)/theta00 * vertical

              case default
                stop "timeStep: unknown case model."
              end select

              ! TFC FJ
              if(topography) then
                bMaxNew = bMaxNew / jac(i, j, k)
              end if

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
            > 1.e-2 * wMax) .and. ((bMax(1) /= 0.) .and. (bMax(2) /= 0.) .and. &
            (bMax(3) /= 0.))) then

          dtBuoy_loc = max(- uMax / bMax(1) + sqrt((uMax / bMax(1)) ** 2 + 2. &
              * cfl * dx / bMax(1)), - vMax / bMax(2) + sqrt((vMax / bMax(2)) &
              ** 2 + 2. * cfl * dy / bMax(2)), - wMax / bMax(3) + sqrt((wMax &
              / bMax(3)) ** 2 + 2. * cfl * dz / bMax(3)))

          !xxxx debug
          if(dtBuoy_loc * tRef < 1.e-2) then

            print *, "dtBuoy_loc*tRef  = ", dtBuoy_loc * tRef
            print *, "bMax(3) = ", bMax(3) * FrInv2

          end if
          !xxxx end debug

        else

          dtBuoy_loc = 1.0e20 / tRef ! set to high value if not needed

        end if

        ! find global minimum

        call mpi_reduce(dtBuoy_loc, dtBuoy, 1, mpi_double_precision, mpi_min, &
            root, comm, ierror)

        call mpi_bcast(dtBuoy, 1, mpi_double_precision, root, comm, ierror)

        !---------------------------
        !   von Neumann condition
        !----------------------------

        dtVisc = 0.5 * min(dx ** 2, dy ** 2, dz ** 2) * Re
        dtCond = 0.5 * min(dx ** 2, dy ** 2, dz ** 2) / mu_conduct

        !----------------------------
        !    Maximal time step
        !----------------------------
        dtMax = dtMax_dim / tRef

        !------------------------------------
        !    Gravity wave time period
        !------------------------------------

        !UAB
        !dtWave = pi/(NN+small)
        !FS: to be consistent with Rieper et al. (2013)
        dtWave = 1. / (NN + small) !1.7/(NN+small)
        !UAE

        !------------------------------------
        !     WKB "CFL" criterion
        !------------------------------------

        if(raytracer) then
          dtWKB_loc = dz / (cgz_max + small)

          if(sizeX > 1) dtWKB_loc = min(dtWKB_loc, dx / (cgx_max + small))
          if(sizeY > 1) dtWKB_loc = min(dtWKB_loc, dy / (cgy_max + small))

          dtWKB_loc = cfl_wave * dtWKB_loc

          ! find global minimum

          call mpi_reduce(dtWKB_loc, dtWKB, 1, mpi_double_precision, mpi_min, &
              root, comm, ierror)

          call mpi_bcast(dtWKB, 1, mpi_double_precision, root, comm, ierror)
        end if

        ! if (testCase == "hotBubble_heat")then! .or. testCase == "hotBubble_heatedLayer")then

        !    buoyMax = 0.
        !    the_max = 0.
        !    the_min = 1.e20
        !    do k = 1,nz
        !       do j = 1,ny
        !          do i = 1,nx

        !             buoyMaxNew = -g_ndim*var(i,j,k,1)/(var(i,j,k,1)+rhoStrat(k))

        !             if (buoyMaxNew > buoyMax)then
        !                buoyMax = buoyMaxNew
        !             end if

        !             the_New = PStrat(k)/(var(i,j,k,1)+rhoStrat(k))

        !             if (the_New > the_max)then
        !                the_max = the_New
        !             end if
        !             if (the_New < the_min) then
        !                the_min = the_New
        !             end if

        !          end do
        !       end do
        !    end do

        !    dtBuoy_loc = cfl *sqrt(dx*the_min/(g_ndim*(the_max-the_min)))!* buoyMax
        !    ! find global minimum

        !    call mpi_reduce(dtBuoy_loc, dtBuoy, 1, mpi_double_precision,&
        !         & mpi_min, root, comm, ierror)

        !    call mpi_bcast(dtBuoy, 1, mpi_double_precision, root, comm, &
        !         & ierror)

        !    end if

        !-------------------------------
        !        Make your choice
        !-------------------------------

        if(dtWave_on .and. timeScheme /= 'semiimplicit') then
          dt = min(dtVisc, dtCond, dtConv, dtMax, dtBuoy, dtWave)
        else
          dt = min(dtVisc, dtCond, dtConv, dtMax, dtBuoy)
          !if (timeScheme == 'semiimplicit') then
          !   dt = min(dtVisc,dtCond,dtConv,dtMax,dtBuoy,10.*dtWave)
          !  else
          !   dt = min(dtVisc,dtCond,dtConv,dtMax,dtBuoy)
          !end if
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
          write(*, fmt = "(a25,es15.1,a8)") "dtBuoy =", dtBuoy * tRef, "seconds"
          write(*, fmt = "(a25,es15.1,a8)") "dtWave =", dtWave * tRef, "seconds"
          if(raytracer) then
            write(*, fmt = "(a25,es15.1,a8)") "dtWKB =", dtWKB * tRef, "seconds"
          end if
          print *, ""

          if(dt == dtMax) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtMax = ", dt * tRef, &
                "seconds"
          else if(dt == dtConv) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtConv = ", dt * tRef, &
                "seconds"
          else if(dt == dtVisc) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtVisc = ", dt * tRef, &
                "seconds"
          else if(dt == dtCond) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtCond = ", dt * tRef, &
                "seconds"
          else if(dt == dtBuoy) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtBuoy = ", dt * tRef, &
                "seconds"
          else if(dt == dtWave) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtWave = ", dt * tRef, &
                "seconds"
          else if(dt == dtWKB) then
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = dtWKB =", dt * tRef, &
                "seconds"
          else
            write(*, fmt = "(a25,es15.1,a8)") "--> dt = ????? = ", dt * tRef, &
                "seconds"
          end if
          print *, ""

        end if

      case default
        stop "timestep: unknown case model."
      end select ! WKB / full model

    end if

    ! error handling for too small time steps
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
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var

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
        Sij_smth
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
                grid anisotropies!'
          elseif(dy / dx > 10.) then
            print *, 'WARNING: dy/dx > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                grid anisotropies!'
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

          uL = var(i - 1, j, k, 2)
          uR = var(i, j, k, 2)

          ! UA not sure whether this averaging is necessary.
          ! Compare, e.g. the viscous fluxes. There it is not done.

          ! replied by JW (20160824): For now, there is no revision on
          ! this part, since it is a relatively minor issue at the moment.

          uB = 0.5 * (var(i - 1, j - 1, k, 2) + var(i, j - 1, k, 2))
          uF = 0.5 * (var(i - 1, j + 1, k, 2) + var(i, j + 1, k, 2))
          uD = 0.5 * (var(i - 1, j, k - 1, 2) + var(i, j, k - 1, 2))
          uU = 0.5 * (var(i - 1, j, k + 1, 2) + var(i, j, k + 1, 2))

          vL = 0.5 * (var(i - 1, j, k, 3) + var(i - 1, j - 1, k, 3))
          vR = 0.5 * (var(i + 1, j, k, 3) + var(i + 1, j - 1, k, 3))
          vB = var(i, j - 1, k, 3)
          vF = var(i, j, k, 3)
          vD = 0.5 * (var(i, j, k - 1, 3) + var(i, j - 1, k - 1, 3))
          vU = 0.5 * (var(i, j, k + 1, 3) + var(i, j - 1, k + 1, 3))

          wL = 0.5 * (var(i - 1, j, k - 1, 4) + var(i - 1, j, k, 4))
          wR = 0.5 * (var(i + 1, j, k - 1, 4) + var(i + 1, j, k, 4))
          wB = 0.5 * (var(i, j - 1, k - 1, 4) + var(i, j - 1, k, 4))
          wF = 0.5 * (var(i, j + 1, k - 1, 4) + var(i, j + 1, k, 4))
          wD = var(i, j, k - 1, 4)
          wU = var(i, j, k, 4)

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

          Sij(i, j, k, 1, 1) = 0.5 * (du_dx + du_dx)
          Sij(i, j, k, 1, 2) = 0.5 * (du_dy + dv_dx)
          Sij(i, j, k, 1, 3) = 0.5 * (du_dz + dw_dx)
          Sij(i, j, k, 2, 1) = 0.5 * (dv_dx + du_dy)
          Sij(i, j, k, 2, 2) = 0.5 * (dv_dy + dv_dy)
          Sij(i, j, k, 2, 3) = 0.5 * (dv_dz + dw_dy)
          Sij(i, j, k, 3, 1) = 0.5 * (dw_dx + du_dz)
          Sij(i, j, k, 3, 2) = 0.5 * (dw_dy + dv_dz)
          Sij(i, j, k, 3, 3) = 0.5 * (dw_dz + dw_dz)

          S_norm(i, j, k) = 0.0

          do jw = 1, 3
            do iw = 1, 3
              S_norm(i, j, k) = S_norm(i, j, k) + Sij(i, j, k, iw, jw) &
                  * Sij(i, j, k, iw, jw)
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
          ui_smth(i, j, k, 1) = 0.5 * (var(i, j, k, 2) + var(i - 1, j, k, 2))
          ui_smth(i, j, k, 2) = 0.5 * (var(i, j, k, 3) + var(i, j - 1, k, 3))
          ui_smth(i, j, k, 3) = 0.5 * (var(i, j, k, 4) + var(i, j, k - 1, 4))

          Sn_smth(i, j, k) = S_norm(i, j, k)

          do jw = 1, 3
            do iw = 1, 3
              uiuj_smth(i, j, k, iw, jw) = ui_smth(i, j, k, iw) * ui_smth(i, &
                  j, k, jw)

              Sij_smth(i, j, k, iw, jw) = Sij(i, j, k, iw, jw)

              S_Sij_smth(i, j, k, iw, jw) = S_norm(i, j, k) * Sij(i, j, k, iw, &
                  jw)
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
            "XZ_local_smth")
      end do

      call Var3DSmthDySma(Sn_smth(1:nx, 1:ny, 1:nz), smth_npts1_DySma, &
          "XZ_local_smth")

      do jw = 1, 3
        do iw = 1, 3
          call Var3DSmthDySma(uiuj_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              smth_npts1_DySma, "XZ_local_smth")
          call Var3DSmthDySma(Sij_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              smth_npts1_DySma, "XZ_local_smth")
          call Var3DSmthDySma(S_Sij_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              smth_npts1_DySma, "XZ_local_smth")
        end do
      end do
    elseif(nx .eq. 1) then
      do iw = 1, 3
        call Var3DSmthDySma(ui_smth(1:nx, 1:ny, 1:nz, iw), smth_npts1_DySma, &
            "YZ_local_smth")
      end do

      call Var3DSmthDySma(Sn_smth(1:nx, 1:ny, 1:nz), smth_npts1_DySma, &
          "YZ_local_smth")

      do jw = 1, 3
        do iw = 1, 3
          call Var3DSmthDySma(uiuj_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              smth_npts1_DySma, "YZ_local_smth")
          call Var3DSmthDySma(Sij_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              smth_npts1_DySma, "YZ_local_smth")
          call Var3DSmthDySma(S_Sij_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              smth_npts1_DySma, "YZ_local_smth")
        end do
      end do
    else
      do iw = 1, 3
        call Var3DSmthDySma(ui_smth(1:nx, 1:ny, 1:nz, iw), smth_npts1_DySma, &
            "XYZ_local_smth")
      end do

      call Var3DSmthDySma(Sn_smth(1:nx, 1:ny, 1:nz), smth_npts1_DySma, &
          "XYZ_local_smth")

      do jw = 1, 3
        do iw = 1, 3
          call Var3DSmthDySma(uiuj_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              smth_npts1_DySma, "XYZ_local_smth")
          call Var3DSmthDySma(Sij_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              smth_npts1_DySma, "XYZ_local_smth")
          call Var3DSmthDySma(S_Sij_smth(1:nx, 1:ny, 1:nz, iw, jw), &
              smth_npts1_DySma, "XYZ_local_smth")
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
                  j, k, iw) * ui_smth(i, j, k, jw)

              Mij(i, j, k, iw, jw) = S_Sij_smth(i, j, k, iw, jw) - (2.0 &
                  * smth_npts1_DySma + 1.0) ** 2 * Sn_smth(i, j, k) &
                  * Sij_smth(i, j, k, iw, jw)

              ! allow for grid anisotropy

              if(iw == 3 .or. jw == 3) then
                Mij(i, j, k, iw, jw) = Mij(i, j, k, iw, jw) * delta_vs
              else
                Mij(i, j, k, iw, jw) = Mij(i, j, k, iw, jw) * delta_hs
              end if

              LijMij_smth(i, j, k) = LijMij_smth(i, j, k) + Lij(i, j, k, iw, &
                  jw) * Mij(i, j, k, iw, jw)

              MijMij_smth(i, j, k) = MijMij_smth(i, j, k) + Mij(i, j, k, iw, &
                  jw) * Mij(i, j, k, iw, jw)
            end do
          end do
        end do
      end do
    end do

    if(ny .eq. 1) then
      call Var3DSmthDySma(LijMij_smth(1:nx, 1:ny, 1:nz), smth_npts2_DySma, &
          "XZ_local_smth")
      call Var3DSmthDySma(MijMij_smth(1:nx, 1:ny, 1:nz), smth_npts2_DySma, &
          "XZ_local_smth")
    elseif(nx .eq. 1) then
      call Var3DSmthDySma(LijMij_smth(1:nx, 1:ny, 1:nz), smth_npts2_DySma, &
          "YZ_local_smth")
      call Var3DSmthDySma(MijMij_smth(1:nx, 1:ny, 1:nz), smth_npts2_DySma, &
          "YZ_local_smth")
    else
      call Var3DSmthDySma(LijMij_smth(1:nx, 1:ny, 1:nz), smth_npts2_DySma, &
          "XYZ_local_smth")
      call Var3DSmthDySma(MijMij_smth(1:nx, 1:ny, 1:nz), smth_npts2_DySma, &
          "XYZ_local_smth")
    end if

    !---------------------------------
    !         Get the final results
    !---------------------------------

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

          if(MijMij_smth(i, j, k) /= 0.) then
            CS2_DySma(i, j, k) = 0.5 * LijMij_smth(i, j, k) / MijMij_smth(i, &
                j, k)
          else
            CS2_DySma(i, j, k) = 0.
          end if

          if(CS2_DySma(i, j, k) < 0.0) then
            CS2_DySma(i, j, k) = 0.0
          end if

          var(i, j, k, 7) = CS2_DySma(i, j, k) * S_norm(i, j, k)
        end do
      end do
    end do

    ! *** set values for the ghost celss ***
    call setHaloAndBoundary(var(:, :, :, 7), nbx, nby, nbz)

    ! deallocate local fields
    deallocate(Sij, stat = allocstat); if(allocstat /= 0) stop &
        "update.f90:dealloc failed"
    deallocate(Lij, stat = allocstat); if(allocstat /= 0) stop &
        "update.f90:dealloc failed"
    deallocate(Mij, stat = allocstat); if(allocstat /= 0) stop &
        "update.f90:dealloc failed"
    deallocate(S_norm, stat = allocstat); if(allocstat /= 0) stop &
        "update.f90:dealloc failed"
    deallocate(uiuj_smth, stat = allocstat); if(allocstat /= 0) stop &
        "update.f90:dealloc failed"
    deallocate(S_Sij_smth, stat = allocstat); if(allocstat /= 0) stop &
        "update.f90:dealloc failed"
    deallocate(Sij_smth, stat = allocstat); if(allocstat /= 0) stop &
        "update.f90:dealloc failed"
    deallocate(ui_smth, stat = allocstat); if(allocstat /= 0) stop &
        "update.f90:dealloc failed"
    deallocate(Sn_smth, stat = allocstat); if(allocstat /= 0) stop &
        "update.f90:dealloc failed"
    deallocate(LijMij_smth, stat = allocstat); if(allocstat /= 0) stop &
        "update.f90:dealloc failed"
    deallocate(MijMij_smth, stat = allocstat); if(allocstat /= 0) stop &
        "update.f90:dealloc failed"
    deallocate(CS2_DySma, stat = allocstat); if(allocstat /= 0) stop &
        "update.f90:dealloc failed"

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
        - nsmth_DySma):(ny + nsmth_DySma), (0 - nsmth_DySma):(nz &
        + nsmth_DySma)), stat = allocstat)
    if(allocstat /= 0) stop "Var3DSmthDySma:alloc failed"

    ! set the values for var3D_DySma_Extend

    var3D_DySma_Extend(1:nx, 1:ny, 1:nz) = var3D_DySma(1:nx, 1:ny, 1:nz)

    call setHaloAndBoundary(var3D_DySma_Extend(:, :, :), nsmth_DySma, &
        nsmth_DySma, nsmth_DySma)

    ! start to do the smoothing

    i0 = is + nbx - 1
    j0 = js + nby - 1

    select case(homog_dir_DySma)

    case("XYZ_local_smth")

      if(xBoundary /= "periodic") stop "DYNAMIC SMAGORINSKY NOT READY FOR &
          NON-PERIODIC BOUNDARY CONDITIONS IN X"

      if(yBoundary /= "periodic") stop "DYNAMIC SMAGORINSKY NOT READY FOR &
          NON-PERIODIC BOUNDARY CONDITIONS IN Y"

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
                + nsmth_DySma), (j - nsmth_DySma):(j + nsmth_DySma), &
                kmin:kmax)) / (((2 * nsmth_DySma + 1) ** 2) * nsmthv)
            !UAE
          end do
        end do
      end do

    case("XZ_local_smth")

      if(xBoundary /= "periodic") stop "DYNAMIC SMAGORINSKY NOT READY FOR &
          NON-PERIODIC BOUNDARY CONDITIONS IN X"

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
                + nsmth_DySma), j, kmin:kmax)) / ((2 * nsmth_DySma + 1) &
                * nsmthv)
            !UAE
          end do
        end do
      end do

      ! gagab
    case("YZ_local_smth")

      if(yBoundary /= "periodic") stop "DYNAMIC SMAGORINSKY NOT READY FOR &
          NON-PERIODIC BOUNDARY CONDITIONS IN Y"

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
            var3D_DySma(i, j, k) = var3D_DySma(i, j, k) / nsmthall
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
                - nsmth_DySma):(j + nsmth_DySma), kmin:kmax)) / ((2 &
                * nsmth_DySma + 1) * nsmthv)
            !UAE
          end do
        end do
      end do
      ! gagae

    case("X_whole_smth")

      if(xBoundary /= "periodic") stop "DYNAMIC SMAGORINSKY NOT READY FOR &
          NON-PERIODIC BOUNDARY CONDITIONS IN X"

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
        "update.f90:dealloc failed"

    return

  end subroutine Var3DSmthDySma

  subroutine setHaloAndBoundary(var3D_HaloBC, nbx_HaloBC, nby_HaloBC, &
      nbz_HaloBC)
    !--------------------------------------
    ! set Halo and Boundary
    !--------------------------------------

    ! in/out variables
    integer, intent(in) :: nbx_HaloBC, nby_HaloBC, nbz_HaloBC

    real, dimension(- nbx_HaloBC:nx + nbx_HaloBC, - nby_HaloBC:ny &
        + nby_HaloBC, - nbz_HaloBC:nz + nbz_HaloBC), intent(inout) :: &
        var3D_HaloBC

    ! auxiliary fields for "var" with ghost cells (rho)
    real, dimension(nbx_HaloBC, - nby_HaloBC:ny + nby_HaloBC, nz) :: &
        xRhoSliceLeft_send, xRhoSliceRight_send
    real, dimension(nbx_HaloBC, - nby_HaloBC:ny + nby_HaloBC, nz) :: &
        xRhoSliceLeft_recv, xRhoSliceRight_recv

    real, dimension(- nbx_HaloBC:nx + nbx_HaloBC, nby_HaloBC, nz) :: &
        yRhoSliceBack_send, yRhoSliceForw_send
    real, dimension(- nbx_HaloBC:nx + nbx_HaloBC, nby_HaloBC, nz) :: &
        yRhoSliceBack_recv, yRhoSliceForw_recv

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
            = var3D_HaloBC(i, - nby_HaloBC:ny + nby_HaloBC, 1:nz)

        xRhoSliceRight_send(i, - nby_HaloBC:ny + nby_HaloBC, 1:nz) &
            = var3D_HaloBC(nx - nbx_HaloBC + i, - nby_HaloBC:ny + nby_HaloBC, &
            1:nz)
      end do

      ! left -> right
      source = left
      dest = right
      tag = 100

      i0 = 1; j0 = - nby_HaloBC; k0 = 1

      call mpi_sendrecv(xRhoSliceRight_send(i0, j0, k0), sendcount, &
          mpi_double_precision, dest, tag, xRhoSliceLeft_recv(i0, j0, k0), &
          recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xRhoSliceLeft_send(i0, j0, k0), sendcount, &
          mpi_double_precision, dest, tag, xRhoSliceRight_recv(i0, j0, k0), &
          recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          sts_right, ierror)

      ! write auxiliary slice to var field
      do i = 1, nbx_HaloBC
        ! right halos
        var3D_HaloBC(nx + i, - nby_HaloBC:ny + nby_HaloBC, 1:nz) &
            = xRhoSliceRight_recv(i, - nby_HaloBC:ny + nby_HaloBC, 1:nz)

        ! left halos
        var3D_HaloBC(- nbx_HaloBC + i, - nby_HaloBC:ny + nby_HaloBC, 1:nz) &
            = xRhoSliceLeft_recv(i, - nby_HaloBC:ny + nby_HaloBC, 1:nz)
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
            = var3D_HaloBC(- nbx_HaloBC:nx + nbx_HaloBC, j, 1:nz)

        yRhoSliceForw_send(- nbx_HaloBC:nx + nbx_HaloBC, j, 1:nz) &
            = var3D_HaloBC(- nbx_HaloBC:nx + nbx_HaloBC, ny - nby_HaloBC + j, &
            1:nz)
      end do

      ! back -> forw
      source = back
      dest = forw
      tag = 100

      i0 = - nbx_HaloBC; j0 = 1; k0 = 1

      call mpi_sendrecv(yRhoSliceForw_send(i0, j0, k0), sendcount, &
          mpi_double_precision, dest, tag, yRhoSliceBack_recv(i0, j0, k0), &
          recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(yRhoSliceBack_send(i0, j0, k0), sendcount, &
          mpi_double_precision, dest, tag, yRhoSliceForw_recv(i0, j0, k0), &
          recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          sts_forw, ierror)

      ! write auxiliary slice to var field
      do j = 1, nby_HaloBC
        ! right halos
        var3D_HaloBC(- nbx_HaloBC:nx + nbx_HaloBC, ny + j, 1:nz) &
            = yRhoSliceForw_recv(- nbx_HaloBC:nx + nbx_HaloBC, j, 1:nz)

        ! left halos
        var3D_HaloBC(- nbx_HaloBC:nx + nbx_HaloBC, - nby_HaloBC + j, 1:nz) &
            = yRhoSliceBack_recv(- nbx_HaloBC:nx + nbx_HaloBC, j, 1:nz)
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

  subroutine smooth_shapiro_0(var)

    !--------------------------------------------------------------------
    !    local smoothing of density, winds, pressure,
    !    and density fluctuations
    !    use Shapiro weighting up to nsmth = 4
    !-------------------------------------------------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var

    ! allocatable fields
    real, dimension(:, :, :), allocatable :: field, field_0, field_1

    integer :: allocstat
    integer :: i, j, k
    integer :: nsmth
    integer :: iVar, ivmax

    allocate(field(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        = allocstat)
    if(allocstat /= 0) stop "smooth_shapiro:alloc failed"

    allocate(field_0(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        = allocstat)
    if(allocstat /= 0) stop "smooth_shapiro:alloc failed"

    allocate(field_1(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        = allocstat)
    if(allocstat /= 0) stop "smooth_shapiro:alloc failed"

    if(timeScheme == "semiimplicit") then
      ! in explicit integration smoothing of density, winds,
      ! and density fluctuations
      ! pressure fluctuations are not smoothened
      ivmax = 6
    else
      ! in explicit integration smoothing of density, and winds
      ivmax = 4
    end if

    do iVar = 1, ivmax
      if(iVar == 5) goto 100

      field(:, :, :) = var(:, :, :, iVar)

      ! set the values for field_0

      field_0 = field

      ! start to do the smoothing

      if(sizeX > 1 .and. sizeY > 1 .and. sizeZ > 1) then
        ! 3D smoothing

        if(nbx == 1 .and. nby == 1 .and. nbz == 1) then
          nsmth = 1
        else if(nbx == 2 .and. nby == 2 .and. nbz == 2) then
          nsmth = 2
        else if(nbx == 3 .and. nby == 3 .and. nbz == 3) then
          nsmth = 3
        else if(nbx == 4 .and. nby == 4 .and. nbz == 4) then
          nsmth = 4
        else
          stop 'ERROR: wrong nbx, nby, nbz in smoothing'
        end if

        if(nsmth == 1) then
          ! smooth in x

          do k = - nbz, nz + nbz
            do j = - nby, ny + nby
              do i = 1, nx
                field_0(i, j, k) = (field(i - 1, j, k) + field(i + 1, j, k) &
                    + 2.0 * field(i, j, k)) / 4.0
              end do
            end do
          end do

          ! smooth in y

          do k = - nbz, nz + nbz
            do j = 1, ny
              do i = 1, nx
                field_1(i, j, k) = (field_0(i, j - 1, k) + field_0(i, j + 1, &
                    k) + 2.0 * field_0(i, j, k)) / 4.0
              end do
            end do
          end do

          ! smooth in z

          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                field(i, j, k) = (field_1(i, j, k - 1) + field_1(i, j, k + 1) &
                    + 2.0 * field_1(i, j, k)) / 4.0
              end do
            end do
          end do
        elseif(nsmth == 2) then
          ! smooth in x

          do k = - nbz, nz + nbz
            do j = - nby, ny + nby
              do i = 1, nx
                field_0(i, j, k) = (- field(i - 2, j, k) - field(i + 2, j, k) &
                    + 4.0 * (field(i - 1, j, k) + field(i + 1, j, k)) + 10.0 &
                    * field(i, j, k)) / 16.0
              end do
            end do
          end do

          ! smooth in y

          do k = - nbz, nz + nbz
            do j = 1, ny
              do i = 1, nx
                field_1(i, j, k) = (- field_0(i, j - 2, k) - field_0(i, j + 2, &
                    k) + 4.0 * (field_0(i, j - 1, k) + field_0(i, j + 1, k)) &
                    + 10.0 * field_0(i, j, k)) / 16.0
              end do
            end do
          end do

          ! smooth in z

          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                field(i, j, k) = (- field_1(i, j, k - 2) - field_1(i, j, k &
                    + 2) + 4.0 * (field_1(i, j, k - 1) + field_1(i, j, k + 1)) &
                    + 10.0 * field_1(i, j, k)) / 16.0
              end do
            end do
          end do
        elseif(nsmth == 3) then
          ! smooth in x

          do k = - nbz, nz + nbz
            do j = - nby, ny + nby
              do i = 1, nx
                field_0(i, j, k) = (field(i - 3, j, k) + field(i + 3, j, k) &
                    - 6.0 * (field(i - 2, j, k) + field(i + 2, j, k)) + 15.0 &
                    * (field(i - 1, j, k) + field(i + 1, j, k)) + 44.0 &
                    * field(i, j, k)) / 64.0
              end do
            end do
          end do

          ! smooth in y

          do k = - nbz, nz + nbz
            do j = 1, ny
              do i = 1, nx
                field_1(i, j, k) = (field_0(i, j - 3, k) + field_0(i, j + 3, &
                    k) - 6.0 * (field_0(i, j - 2, k) + field_0(i, j + 2, k)) &
                    + 15.0 * (field_0(i, j - 1, k) + field_0(i, j + 1, k)) &
                    + 44.0 * field_0(i, j, k)) / 64.0
              end do
            end do
          end do

          ! smooth in z

          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                field(i, j, k) = (field_1(i, j, k - 3) + field_1(i, j, k + 3) &
                    - 6.0 * (field_1(i, j, k - 2) + field_1(i, j, k + 2)) &
                    + 15.0 * (field_1(i, j, k - 1) + field_1(i, j, k + 1)) &
                    + 44.0 * field_1(i, j, k)) / 64.0
              end do
            end do
          end do
        elseif(nsmth == 4) then
          ! smooth in x

          do k = - nbz, nz + nbz
            do j = - nby, ny + nby
              do i = 1, nx
                field_0(i, j, k) = (- field(i - 4, j, k) - field(i + 4, j, k) &
                    + 8.0 * (field(i - 3, j, k) + field(i + 3, j, k)) - 28.0 &
                    * (field(i - 2, j, k) + field(i + 2, j, k)) + 56.0 &
                    * (field(i - 1, j, k) + field(i + 1, j, k)) + 186.0 &
                    * field(i, j, k)) / 256.0
              end do
            end do
          end do

          ! smooth in y

          do k = - nbz, nz + nbz
            do j = 1, ny
              do i = 1, nx
                field_1(i, j, k) = (- field_0(i, j - 4, k) - field_0(i, j + 4, &
                    k) + 8.0 * (field_0(i, j - 3, k) + field_0(i, j + 3, k)) &
                    - 28.0 * (field_0(i, j - 2, k) + field_0(i, j + 2, k)) &
                    + 56.0 * (field_0(i, j - 1, k) + field_0(i, j + 1, k)) &
                    + 186.0 * field_0(i, j, k)) / 256.0
              end do
            end do
          end do

          ! smooth in z

          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                field(i, j, k) = (field_1(i, j, k - 4) + field_1(i, j, k + 4) &
                    + 8.0 * (field_1(i, j, k - 3) + field_1(i, j, k + 3)) &
                    - 28.0 * (field_1(i, j, k - 2) + field_1(i, j, k + 2)) &
                    + 56.0 * (field_1(i, j, k - 1) + field_1(i, j, k + 1)) &
                    + 186.0 * field_1(i, j, k)) / 256.0
              end do
            end do
          end do
        end if
      else if(sizeX > 1 .and. sizeY == 1 .and. sizeZ > 1) then
        ! 2D smoothing in x and z

        if(nbx == 1 .and. nbz == 1) then
          nsmth = 1
        else if(nbx == 2 .and. nbz == 2) then
          nsmth = 2
        else if(nbx == 3 .and. nbz == 3) then
          nsmth = 3
        else if(nbx == 4 .and. nbz == 4) then
          nsmth = 4
        else
          stop 'ERROR: wrong nbx, nby, nbz in smoothing'
        end if

        if(nsmth == 1) then
          ! smooth in x

          do k = - nbz, nz + nbz
            do j = 1, ny
              do i = 1, nx
                field_1(i, j, k) = (field_0(i - 1, j, k) + field_0(i + 1, j, &
                    k) + 2.0 * field_0(i, j, k)) / 4.0
              end do
            end do
          end do

          ! smooth in z

          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                field(i, j, k) = (field_1(i, j, k - 1) + field_1(i, j, k + 1) &
                    + 2.0 * field_1(i, j, k)) / 4.0
              end do
            end do
          end do
        elseif(nsmth == 2) then
          ! smooth in x

          do k = - nbz, nz + nbz
            do j = 1, ny
              do i = 1, nx
                field_1(i, j, k) = (- field_0(i - 2, j, k) - field_0(i + 2, j, &
                    k) + 4.0 * (field_0(i - 1, j, k) + field_0(i + 1, j, k)) &
                    + 10.0 * field_0(i, j, k)) / 16.0
              end do
            end do
          end do

          ! smooth in z

          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                field(i, j, k) = (- field_1(i, j, k - 2) - field_1(i, j, k &
                    + 2) + 4.0 * (field_1(i, j, k - 1) + field_1(i, j, k + 1)) &
                    + 10.0 * field_1(i, j, k)) / 16.0
              end do
            end do
          end do
        elseif(nsmth == 3) then
          ! smooth in x

          do k = - nbz, nz + nbz
            do j = 1, ny
              do i = 1, nx
                field_1(i, j, k) = (field_0(i - 3, j, k) + field_0(i + 3, j, &
                    k) - 6.0 * (field_0(i - 2, j, k) + field_0(i + 2, j, k)) &
                    + 15.0 * (field_0(i - 1, j, k) + field_0(i + 1, j, k)) &
                    + 44.0 * field_0(i, j, k)) / 64.0
              end do
            end do
          end do

          ! smooth in z

          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                field(i, j, k) = (field_1(i, j, k - 3) + field_1(i, j, k + 3) &
                    - 6.0 * (field_1(i, j, k - 2) + field_1(i, j, k + 2)) &
                    + 15.0 * (field_1(i, j, k - 1) + field_1(i, j, k + 1)) &
                    + 44.0 * field_1(i, j, k)) / 64.0
              end do
            end do
          end do
        elseif(nsmth == 4) then
          ! smooth in x

          do k = - nbz, nz + nbz
            do j = 1, ny
              do i = 1, nx
                field_1(i, j, k) = (- field_0(i - 4, j, k) - field_0(i + 4, j, &
                    k) + 8.0 * (field_0(i - 3, j, k) + field_0(i + 3, j, k)) &
                    - 28.0 * (field_0(i - 2, j, k) + field_0(i + 2, j, k)) &
                    + 56.0 * (field_0(i - 1, j, k) + field_0(i + 1, j, k)) &
                    + 186.0 * field_0(i, j, k)) / 256.0
              end do
            end do
          end do

          ! smooth in z

          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                field(i, j, k) = (field_1(i, j, k - 4) + field_1(i, j, k + 4) &
                    + 8.0 * (field_1(i, j, k - 3) + field_1(i, j, k + 3)) &
                    - 28.0 * (field_1(i, j, k - 2) + field_1(i, j, k + 2)) &
                    + 56.0 * (field_1(i, j, k - 1) + field_1(i, j, k + 1)) &
                    + 186.0 * field_1(i, j, k)) / 256.0
              end do
            end do
          end do
        end if
      else if(sizeX == 1 .and. sizeY > 1 .and. sizeZ > 1) then
        if(nby == 1 .and. nbz == 1) then
          nsmth = 1
        else if(nby == 2 .and. nbz == 2) then
          nsmth = 2
        else if(nby == 3 .and. nbz == 3) then
          nsmth = 3
        else if(nby == 4 .and. nbz == 4) then
          nsmth = 4
        else
          stop 'ERROR: wrong nbx, nby, nbz in smoothing'
        end if

        if(nsmth == 1) then
          ! smooth in y

          do k = - nbz, nz + nbz
            do j = 1, ny
              do i = 1, nx
                field_1(i, j, k) = (field_0(i, j - 1, k) + field_0(i, j + 1, &
                    k) + 2.0 * field_0(i, j, k)) / 4.0
              end do
            end do
          end do

          ! smooth in z

          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                field(i, j, k) = (field_1(i, j, k - 1) + field_1(i, j, k + 1) &
                    + 2.0 * field_1(i, j, k)) / 4.0
              end do
            end do
          end do
        elseif(nsmth == 2) then
          ! smooth in y

          do k = - nbz, nz + nbz
            do j = 1, ny
              do i = 1, nx
                field_1(i, j, k) = (- field_0(i, j - 2, k) - field_0(i, j + 2, &
                    k) + 4.0 * (field_0(i, j - 1, k) + field_0(i, j + 1, k)) &
                    + 10.0 * field_0(i, j, k)) / 16.0
              end do
            end do
          end do

          ! smooth in z

          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                field(i, j, k) = (- field_1(i, j, k - 2) - field_1(i, j, k &
                    + 2) + 4.0 * (field_1(i, j, k - 1) + field_1(i, j, k + 1)) &
                    + 10.0 * field_1(i, j, k)) / 16.0
              end do
            end do
          end do
        elseif(nsmth == 3) then
          ! smooth in y

          do k = - nbz, nz + nbz
            do j = 1, ny
              do i = 1, nx
                field_1(i, j, k) = (field_0(i, j - 3, k) + field_0(i, j + 3, &
                    k) - 6.0 * (field_0(i, j - 2, k) + field_0(i, j + 2, k)) &
                    + 15.0 * (field_0(i, j - 1, k) + field_0(i, j + 1, k)) &
                    + 44.0 * field_0(i, j, k)) / 64.0
              end do
            end do
          end do

          ! smooth in z

          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                field(i, j, k) = (field_1(i, j, k - 3) + field_1(i, j, k + 3) &
                    - 6.0 * (field_1(i, j, k - 2) + field_1(i, j, k + 2)) &
                    + 15.0 * (field_1(i, j, k - 1) + field_1(i, j, k + 1)) &
                    + 44.0 * field_1(i, j, k)) / 64.0
              end do
            end do
          end do
        elseif(nsmth == 4) then
          ! smooth in y

          do k = - nbz, nz + nbz
            do j = 1, ny
              do i = 1, nx
                field_1(i, j, k) = (- field_0(i, j - 4, k) - field_0(i, j + 4, &
                    k) + 8.0 * (field_0(i, j - 3, k) + field_0(i, j + 3, k)) &
                    - 28.0 * (field_0(i, j - 2, k) + field_0(i, j + 2, k)) &
                    + 56.0 * (field_0(i, j - 1, k) + field_0(i, j + 1, k)) &
                    + 186.0 * field_0(i, j, k)) / 256.0
              end do
            end do
          end do

          ! smooth in z

          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                field(i, j, k) = (field_1(i, j, k - 4) + field_1(i, j, k + 4) &
                    + 8.0 * (field_1(i, j, k - 3) + field_1(i, j, k + 3)) &
                    - 28.0 * (field_1(i, j, k - 2) + field_1(i, j, k + 2)) &
                    + 56.0 * (field_1(i, j, k - 1) + field_1(i, j, k + 1)) &
                    + 186.0 * field_1(i, j, k)) / 256.0
              end do
            end do
          end do
        end if
      else
        stop "ERROR: smoothing not ready for 2D in x and y or 1D"
      end if

      var(:, :, :, iVar) = field(:, :, :)

      100 continue
    end do

    ! deallocate local fields
    deallocate(field, stat = allocstat); if(allocstat /= 0) stop &
        "smooth_shapiro:dealloc failed"
    deallocate(field_0, stat = allocstat); if(allocstat /= 0) stop &
        "smooth_shapiro:dealloc failed"
    deallocate(field_1, stat = allocstat); if(allocstat /= 0) stop &
        "smooth_shapiro:dealloc failed"

    return

  end subroutine smooth_shapiro_0

  !UAB
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
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var
    real, intent(inout) :: fc_shap
    integer, intent(in) :: n_shap
    real, dimension(- 1:nx, - 1:ny, - 1:nz, 3, nVar), intent(inout) :: flux
    real, intent(in) :: dt

    ! allocatable fields
    real, dimension(:, :, :), allocatable :: field
    real, dimension(:, :, :, :), allocatable :: var_l

    integer :: allocstat
    integer :: i, j, k
    integer :: nsmth
    integer :: iVar, ivmax
    integer :: i_lapl
    integer :: nz_max

    allocate(field(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        = allocstat)
    if(allocstat /= 0) stop "smooth_shapiro:alloc failed"

    allocate(var_l(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), stat &
        = allocstat)
    if(allocstat /= 0) stop "smooth_shapiro:alloc failed"

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
              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k) &
                      + rhoStratTFC(i, j, k)
                else
                  var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k) &
                      + rhoStrat(k)
                end if
              else
                var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k)
              end if
            end do
          end do
        end do

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              var_l(i, j, k, 2) = var(i, j, k, 2) - u_env_pp(i, j, k)
            end do
          end do
        end do

        if(timeScheme == "semiimplicit") then
          if(topography) then
            ! TFC FJ
            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  var(i, j, k, 6) = var(i, j, k, 6) - dens_env_pp(i, j, k) &
                      + rhoStratTFC(i, j, k)
                end do
              end do
            end do
          else
            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  var_l(i, j, k, 6) = var(i, j, k, 6) - dens_env_pp(i, j, k) &
                      + rhoStrat(k)
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
        do iVar = 1, ivmax
          if(iVar /= 5) then
            field(:, :, :) = var_l(:, :, :, iVar)

            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  var_l(i, j, k, iVar) = (field(i - 1, j, k) + field(i + 1, j, &
                      k) - 2.0 * field(i, j, k)) / 4.0
                end do
              end do
            end do
          end if
        end do

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

      do iVar = 1, ivmax
        if(iVar /= 5) then
          if(iVar == 4) then
            nz_max = nz - 1
          else
            nz_max = nz
          end if

          do k = 1, nz_max
            do j = 1, ny
              do i = 1, nx

                var(i, j, k, iVar) = var(i, j, k, iVar) + fc_shap * (- 1) &
                    ** (n_shap + 1) * var_l(i, j, k, iVar)

              end do
            end do
          end do
        end if
      end do

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
              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k) &
                      + rhoStratTFC(i, j, k)
                else
                  var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k) &
                      + rhoStrat(k)
                end if
              else
                var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k)
              end if
            end do
          end do
        end do

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              var_l(i, j, k, 2) = var(i, j, k, 2) - u_env_pp(i, j, k)
            end do
          end do
        end do

        if(timeScheme == "semiimplicit") then
          if(topography) then
            ! TFC FJ
            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  var(i, j, k, 6) = var(i, j, k, 6) - dens_env_pp(i, j, k) &
                      + rhoStratTFC(i, j, k)
                end do
              end do
            end do
          else
            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  var_l(i, j, k, 6) = var(i, j, k, 6) - dens_env_pp(i, j, k) &
                      + rhoStrat(k)
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
        do iVar = 1, ivmax
          if(iVar /= 5) then
            field(:, :, :) = var_l(:, :, :, iVar)

            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  var_l(i, j, k, iVar) = (field(i, j - 1, k) + field(i, j + 1, &
                      k) - 2.0 * field(i, j, k)) / 4.0
                end do
              end do
            end do
          end if
        end do

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

      do iVar = 1, ivmax
        if(iVar /= 5) then
          if(iVar == 4) then
            nz_max = nz - 1
          else
            nz_max = nz
          end if

          do k = 1, nz_max
            do j = 1, ny
              do i = 1, nx

                var(i, j, k, iVar) = var(i, j, k, iVar) + fc_shap * (- 1) &
                    ** (n_shap + 1) * var_l(i, j, k, iVar)
              end do
            end do
          end do
        end if
      end do

      ! boundary conditions again

      call setHalos(var, "var")
      call setBoundary(var, flux, "var")
    end if

    !testb
    100 continue
    !teste

    ! deallocate local fields

    deallocate(field, stat = allocstat); if(allocstat /= 0) stop &
        "smooth_shapiro:dealloc failed"
    deallocate(var_l, stat = allocstat); if(allocstat /= 0) stop &
        "smooth_shapiro:dealloc failed"

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
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var
    real, intent(in) :: fc_shap
    integer, intent(in) :: n_shap
    real, dimension(- 1:nx, - 1:ny, - 1:nz, 3, nVar), intent(inout) :: flux

    ! allocatable fields
    real, dimension(:, :, :), allocatable :: field
    real, dimension(:, :, :, :), allocatable :: var_l

    integer :: allocstat
    integer :: i, j, k
    integer :: nsmth
    integer :: iVar, ivmax
    integer :: i_lapl
    integer :: nz_max

    allocate(field(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        = allocstat)
    if(allocstat /= 0) stop "smooth_shapiro:alloc failed"

    allocate(var_l(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), stat &
        = allocstat)
    if(allocstat /= 0) stop "smooth_shapiro:alloc failed"

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
              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k) &
                      + rhoStratTFC(i, j, k)
                else
                  var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k) &
                      + rhoStrat(k)
                end if
              else
                var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k)
              end if
            end do
          end do
        end do

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              var_l(i, j, k, 2) = var(i, j, k, 2) - u_env_pp(i, j, k)
            end do
          end do
        end do

        if(timeScheme == "semiimplicit") then
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                var_l(i, j, k, 6) = var(i, j, k, 6) - dens_env_pp(i, j, k) &
                    + rhoStrat(k)
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
        do iVar = 1, ivmax
          if(iVar /= 5) then

            field(:, :, :) = var_l(:, :, :, iVar)

            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  var_l(i, j, k, iVar) = (field(i - 1, j, k) + field(i + 1, j, &
                      k) - 2.0 * field(i, j, k)) / 4.0
                end do
              end do
            end do
          end if

        end do

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

      do iVar = 1, ivmax
        if(iVar /= 5) then
          if(iVar == 4) then
            nz_max = nz - 1
          else
            nz_max = nz
          end if

          do k = 1, nz_max
            do j = 1, ny
              do i = 1, nx

                var(i, j, k, iVar) = var(i, j, k, iVar) + fc_shap * (- 1) &
                    ** (n_shap + 1) * var_l(i, j, k, iVar)

              end do
            end do
          end do
        end if
      end do

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
              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k) &
                      + rhoStratTFC(i, j, k)
                else
                  var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k) &
                      + rhoStrat(k)
                end if
              else
                var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k)
              end if
            end do
          end do
        end do

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              var_l(i, j, k, 2) = var(i, j, k, 2) - u_env_pp(i, j, k)
            end do
          end do
        end do

        if(timeScheme == "semiimplicit") then
          do k = 1, nz
            do j = 1, ny
              do i = 1, nx
                var_l(i, j, k, 6) = var(i, j, k, 6) - dens_env_pp(i, j, k) &
                    + rhoStrat(k)
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
        do iVar = 1, ivmax
          if(iVar /= 5) then

            field(:, :, :) = var_l(:, :, :, iVar)

            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  var_l(i, j, k, iVar) = (field(i, j - 1, k) + field(i, j + 1, &
                      k) - 2.0 * field(i, j, k)) / 4.0
                end do
              end do
            end do
          end if

        end do

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

      do iVar = 1, ivmax
        if(iVar /= 5) then
          if(iVar == 4) then
            nz_max = nz - 1
          else
            nz_max = nz
          end if

          do k = 1, nz_max
            do j = 1, ny
              do i = 1, nx

                var(i, j, k, iVar) = var(i, j, k, iVar) + fc_shap * (- 1) &
                    ** (n_shap + 1) * var_l(i, j, k, iVar)

              end do
            end do
          end do
        end if
      end do

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
            if(fluctuationMode) then
              if(topography) then
                ! TFC FJ
                var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k) &
                    + rhoStratTFC(i, j, k)
              else
                var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k) &
                    + rhoStrat(k)
              end if
            else
              var_l(i, j, k, 1) = var(i, j, k, 1) - dens_env_pp(i, j, k)
            end if
          end do
        end do
      end do

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            var_l(i, j, k, 2) = var(i, j, k, 2) - u_env_pp(i, j, k)
          end do
        end do
      end do

      if(timeScheme == "semiimplicit") then
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              var_l(i, j, k, 6) = var(i, j, k, 6) - dens_env_pp(i, j, k) &
                  + rhoStrat(k)
            end do
          end do
        end do
      end if

    end if

    do iVar = 1, ivmax
      if(iVar /= 5) then

        field(:, :, :) = var_l(:, :, :, iVar)

        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              var_l(i, j, k, iVar) = (field(i, j, k - 1) + field(i, j, k + 1) &
                  - 2.0 * field(i, j, k)) / 4.0
            end do
          end do
        end do
      end if

    end do

    ! apply filter

    do iVar = 1, ivmax
      if(iVar /= 5) then
        if(iVar == 4) then
          nz_max = nz - 1
        else
          nz_max = nz
        end if

        do k = 1, nz_max
          do j = 1, ny
            do i = 1, nx

              var(i, j, k, iVar) = var(i, j, k, iVar) + 1.e-2 * fc_shap &
                  * var_l(i, j, k, iVar)

            end do
          end do
        end do
      end if
    end do

    ! boundary conditions again

    call setHalos(var, "var")
    call setBoundary(var, flux, "var")

    do k = 1, nz !theta
      do j = 1, ny
        do i = 1, nx
          var(i, j, k, 1) = (var(i, j, k, 1) - rhoStrat(k)) / PStrat(k) !- dens_env_pp(i, j, k) + rhoStrat(k)

        end do
      end do
    end do

    !testb
    100 continue
    !teste

    ! deallocate local fields

    deallocate(field, stat = allocstat); if(allocstat /= 0) stop &
        "smooth_shapiro:dealloc failed"
    deallocate(var_l, stat = allocstat); if(allocstat /= 0) stop &
        "smooth_shapiro:dealloc failed"

    return

  end subroutine smooth_shapiro
  !UAE

  !---------------------------------------------------------------------

  subroutine BGstate_update(var, flux, dt, m, q_P, q_rho, int_mod, &
      heating_switch)

    ! (1) update of the reference-atmosphere profile
    ! (2) synchronization of the vertical wind with the assumption that the
    ! w0 according to O'Neill & Klein (2014) is the horizontal-mean
    ! vertical wind

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, nVar), &
        intent(inout) :: var
    real, dimension(- 1:nx, - 1:ny, - 1:nz, 3, nVar), intent(in) :: flux

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
        mpi_double_precision, mpi_sum, comm, ierror)
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

            if(fluctuationMode) then
              sum_local2(k) = sum_local2(k) + var(i, j, k, 1) + rhoStrat(k)
            else
              sum_local2(k) = sum_local2(k) + var(i, j, k, 1)
            end if

            if(k == 1) then
              sum_local(k) = sum_local(k) + 0.5 * flux(i, j, k, 3, 1)
            else if(k == nz) then
              sum_local(k) = sum_local(k) + 0.5 * flux(i, j, k - 1, 3, 1)
            else
              sum_local(k) = sum_local(k) + 0.5 * (flux(i, j, k - 1, 3, 1) &
                  + flux(i, j, k, 3, 1))
            end if
          end do
        end do
      end do
      call mpi_allreduce(sum_local(1), sum_global(1), nz - 1 + 1, &
          mpi_double_precision, mpi_sum, comm, ierror)
      sum_global = sum_global / (sizeX * sizeY)

      rhow_bar(1:nz) = sum_global(1:nz)

      call mpi_allreduce(sum_local2(1), sum_global2(1), nz - 1 + 1, &
          mpi_double_precision, mpi_sum, comm, ierror)
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
            * press0(k)) * z(k))

        sum_n = sum_n + expo * (- S_bar(k) / PStrat(k) + g_ndim * rhow_bar(k) &
            / (gamma * press0(k)))

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
            * press0(k)) * z(k))

        w_0(k) = w_0(k - 1) + dz * expo * (- S_bar(k) / Pstrat(k) + g_ndim &
            * rhow_bar(k) / (gamma * press0(k)) - dptopdt / (gamma * press0(k)))
      end do

      do k = 1, nz - 1
        expo = exp(g_ndim * (rho_bar(k) + rhoStrat_s(k)) / (gamma * press0(k)) &
            * 0.5 * (z(k) + z(k + 1)))

        w_0(k) = expo * w_0(k)
      end do
    else if(w0_mod == 'ONK14') then
      w_0(1) = dz * (- S_bar(1) / Pstrat(1) - (1. / (gamma * press0(1))) &
          * dptopdt)

      do k = 2, nz - 1
        w_0(k) = w_0(k - 1) + dz * (- S_bar(k) / Pstrat(k) - (1. / (gamma &
            * press0(k))) * dptopdt)
      end do
    else
      stop 'ERROR: wrong w0_mod'
    end if

    ! update PStrat

    do k = 1, nz
      divPw(k) = (PstratTilde(k) * w_0(k) - PstratTilde(k - 1) * w_0(k - 1)) &
          / dz
    end do

    !! save total density and subtract the reference-atmosphere density
    !! from this again after the update of the latter

    !if (fluctuationMode) then
    !   do k = 1,nz
    !      var(:,:,k,1) = var(:,:,k,1) + rhoStrat(k)
    !   end do
    !end if

    !sum_local = 0.
    !sum_global = 0.

    !do k = 1,nz
    !   do j = 1,ny
    !      do i = 1,nx
    !         sum_local(k)  = sum_local(k)  + var(i,j,k,1)
    !      end do
    !   end do
    !end do
    !call mpi_allreduce(sum_local(1),sum_global(1),&
    !     nz,&
    !     mpi_double_precision,mpi_sum,comm,ierror)
    !sum_global = sum_global/(sizeX*sizeY)
    !
    !rhoStrat(1:nz) = rhoStrat_d(1:nz) + sum_global(1:nz)

    do k = 1, nz
      if(int_mod == "expl") then
        !init q
        if(m == 1) then
          q_P(k) = 0.
          q_rho(k) = 0.
        end if

        ! update: q(m-1) -> q(m)

        q_P(k) = alpha(m) * q_P(k) - dt * divPw(k) - dt * S_bar(k)

        ! update PStrat

        Pstrat(k) = Pstrat(k) + beta(m) * q_P(k) !PIold
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

      !update piStrat
      !thetaStrat(k) = PStrat(k)/rhoStrat(k) !FSMar2021
      piStrat(k) = PStrat(k) ** (kappa / (1.0 - kappa))
    end do

    !! adjust stratification so that no unstable layers exist

    !do k = 2,nz
    !   if (rhoStrat(k) > rhoStrat(k-1) * pStrat(k)/pStrat(k-1)) then
    !      rhoStrat(k) = rhoStrat(k-1) * pStrat(k)/pStrat(k-1)
    !      thetaStrat(k) = thetaStrat(k-1)
    !   end if
    !end do

    !sum_local = 0.
    !sum_global = 0.

    !do k = 1,nz
    !   do j = 1,ny
    !      do i = 1,nx
    !         sum_local(k)  = sum_local(k)  + var(i,j,k,1)
    !      end do
    !   end do
    !end do
    !call mpi_allreduce(sum_local(1),sum_global(1),&
    !     nz,&
    !     mpi_double_precision,mpi_sum,comm,ierror)
    !sum_global = sum_global/(sizeX*sizeY)

    !rhoStrat_d = 0.

    !do k=1,nz
    !   rhoStrat_d(k) = rhoStrat(k) - sum_global(k)
    !end do

    pStrat(0) = pStrat(1)
    pStrat(- 1) = pStrat(0)
    pStrat(nz + 1) = pStrat(nz)
    pStrat(nz + 2) = pStrat(nz + 1)

    do k = 1, nz
      PstratTilde(k) = 0.5 * (PStrat(k) + PStrat(k + 1))
      !rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k+1))
      !thetaStratTilde(k) =  PStratTilde(k)/rhoStratTilde(k)
    end do

    !! adjust density fluctuations to new reference atmosphere
    !do k = 1,nz
    !   if (fluctuationMode) then
    !      var(:,:,k,1) = var(:,:,k,1) - rhoStrat(k)
    !   end if

    !   if ((timeScheme == "semiimplicit") .or. auxil_equ) then
    !      var(:,:,k,6) = var(:,:,k,6) - rhoStrat(k)
    !   end if
    !end do

    ! the following could most probably be deleted
    ! update of non-dimensional squared Brunt-Vaisala frequency
    ! (this could perhaps be done a bit nicer)

    N2 = 0.

    do k = 1, nz
      if(k == 1) then
        bvsStrat(k) = g_ndim / thetaStrat(k) * (thetaStrat(k + 1) &
            - thetaStrat(k)) / dz
      else if(k == nz) then
        bvsStrat(k) = g_ndim / thetaStrat(k) * (thetaStrat(k) - thetaStrat(k &
            - 1)) / dz
      else
        bvsStrat(k) = g_ndim / thetaStrat(k) * (thetaStrat(k + 1) &
            - thetaStrat(k - 1)) / (2.0 * dz)
      end if

      N2 = max(N2, bvsStrat(k))
    end do

    bvsStrat(- 1) = bvsStrat(1)
    bvsStrat(0) = bvsStrat(1)

    bvsStrat(nz + 1) = bvsStrat(nz)
    bvsStrat(nz + 2) = bvsStrat(nz)

    ! N2 = max(N2, bvsStrat(nz+1))

    if(N2 < 0.) then
      stop 'ERROR: N2 < 0'
    else
      NN = sqrt(N2)
    end if

    !testb
    ! do k = -1,nz+1
    !    if (master .and. N2 == bvsStrat(k)) print*,'N2 = max at k =',k
    ! end do
    !teste

  end subroutine BGstate_update

  !---------------------------------------------------------------------

  ! TFC FJ
  subroutine momentumPredictorTestTFC(var, flux, force, dMom, int_mod)

    real, dimension((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        + nbz), nVar), intent(inout) :: var
    real, dimension((- 1):nx, (- 1):ny, (- 1):nz, 3, nVar), intent(in) :: flux
    real, dimension(0:(nx + 1), 0:(ny + 1), 0:(nz + 1), 3), intent(in) :: force
    real, dimension((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        + nbz), 3), intent(inout) :: dMom
    character(len = *), intent(in) :: int_mod

    real, dimension((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        + nbz), nVar) :: var_tfc

    var_tfc = var

    if(int_mod == "expl") then
      rhoOld = var(:, :, :, 1)
      call momentumPredictor(var_tfc, flux, force, 1.0, dMom, 1, "tot", &
          "expl", 1.0)
      topography = .false.
      call momentumPredictor(var, flux, force, 1.0, dMom, 1, "tot", "expl", 1.0)
      topography = .true.
      print *, "Momentum predictor difference: ", maxval(abs(var_tfc - var))
    else if(int_mod == "impl") then
      rhoOld = var(:, :, :, 1)
      call momentumPredictor(var_tfc, flux, force, 1.0, dMom, 1, "lhs", &
          "expl", 1.0)
      topography = .false.
      call momentumPredictor(var, flux, force, 1.0, dMom, 1, "lhs", "expl", 1.0)
      topography = .true.
      print *, "Momentum predictor difference (lhs, expl): ", &
          maxval(abs(var_tfc - var))

      var_tfc = var
      rhopOld = var(:, :, :, 6)
      call momentumPredictor(var_tfc, flux, force, 1.0, dMom, 1, "rhs", &
          "impl", 1.0)
      topography = .false.
      call momentumPredictor(var, flux, force, 1.0, dMom, 1, "rhs", "impl", 1.0)
      topography = .true.
      print *, "Momentum predictor difference (rhs, impl): ", &
          maxval(abs(var_tfc - var))

      var_tfc = var
      rhopOld = var(:, :, :, 6)
      call momentumPredictor(var_tfc, flux, force, 1.0, dMom, 1, "rhs", &
          "expl", 1.0)
      topography = .false.
      call momentumPredictor(var, flux, force, 1.0, dMom, 1, "rhs", "expl", 1.0)
      topography = .true.
      print *, "Momentum predictor difference (rhs, expl): ", &
          maxval(abs(var_tfc - var))
    end if

  end subroutine momentumPredictorTestTFC

  !---------------------------------------------------------------------

  ! TFC FJ
  subroutine massUpdateTestTFC(var, flux, dRho, int_mod)

    real, dimension((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        + nbz), nVar), intent(inout) :: var
    real, dimension((- 1):nx, (- 1):ny, (- 1):nz, 3, nVar), intent(in) :: flux
    real, dimension((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        + nbz)), intent(inout) :: dRho
    character(len = *), intent(in) :: int_mod

    real, dimension((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        + nbz), nVar) :: var_tfc

    var_tfc = var

    if(int_mod == "expl") then
      call massUpdate(var_tfc, flux, 1.0, dRho, 1, "rho", "tot", "expl", 1.0)
      call massUpdate(var_tfc, flux, 1.0, dRho, 1, "rhop", "tot", "expl", 1.0)
      topography = .false.
      call massUpdate(var, flux, 1.0, dRho, 1, "rho", "tot", "expl", 1.0)
      call massUpdate(var, flux, 1.0, dRho, 1, "rhop", "tot", "expl", 1.0)
      topography = .true.
      print *, "Mass update difference: ", maxval(abs(var_tfc - var))
    else if(int_mod == "impl") then
      call massUpdate(var_tfc, flux, 1.0, dRho, 1, "rho", "lhs", "expl", 1.0)
      call massUpdate(var_tfc, flux, 1.0, dRho, 1, "rhop", "lhs", "expl", 1.0)
      topography = .false.
      call massUpdate(var, flux, 1.0, dRho, 1, "rho", "lhs", "expl", 1.0)
      call massUpdate(var, flux, 1.0, dRho, 1, "rhop", "lhs", "expl", 1.0)
      topography = .true.
      print *, "Mass update difference (lhs, expl): ", maxval(abs(var_tfc &
          - var))

      var_tfc = var
      wOldTFC = var(:, :, :, 4)
      call massUpdate(var_tfc, flux, 1.0, dRho, 1, "rhop", "rhs", "impl", 1.0)
      topography = .false.
      call massUpdate(var, flux, 1.0, dRho, 1, "rhop", "rhs", "impl", 1.0)
      topography = .true.
      print *, "Mass update difference (rhs, impl): ", maxval(abs(var_tfc &
          - var))

      var_tfc = var
      call massUpdate(var_tfc, flux, 1.0, dRho, 1, "rhop", "rhs", "expl", 1.0)
      topography = .false.
      call massUpdate(var, flux, 1.0, dRho, 1, "rhop", "rhs", "expl", 1.0)
      topography = .true.
      print *, "Mass update difference (rhs, expl): ", maxval(abs(var_tfc &
          - var))
    end if

  end subroutine massUpdateTestTFC

end module update_module
