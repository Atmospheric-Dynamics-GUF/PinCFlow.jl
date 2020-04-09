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
  public :: timestep
  public :: init_update
  public :: set_spongeLayer

  public :: CoefDySma_update         
  public :: Var3DSmthDySma

  public :: setHaloAndBoundary

  public :: smooth_shapiro
  public :: smooth_hor_shapiro

  public :: BGstate_update

  
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
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real, intent(in) :: dt
    character(len=*), intent(in) :: variable
    
    ! local variables
    integer :: i,j,k, iVar

    ! relaxation parameters
    real :: alpha, beta
    real :: spongeAlphaZ, spongeDz
    
    ! variables for rho 
    real    :: rho_old, rho_bg, rho_new
    real    :: uOld, uBG, uNew
    real    :: vOld, vBG, vNew
    real    :: wOld, wBG, wNew
    
    ! variables for ice
    real    :: nAer_bg, nIce_bg, qIce_bg, qv_bg
    real    :: T,p


    real, dimension(1:nz) :: sum_local, sum_global
    
    ! return if sponge layer with relaxation is switched off
    if( .not. spongeLayer ) then 
       return
    end if
    
    ! nondimensionalize relaxation parameter   
    spongeAlphaZ = spongeAlphaZ_dim * tRef

    ! thickness of sponge layer
    spongeDz = z(nz) - z(kSponge) 

    
    select case( variable )

    case( "rho" ) 
       do k = kSponge, nz
          do j = 1,ny
             do i = 1,nx

                if ((TestCase == "baroclinic_LC") &
                  & .or.(TestCase == "baroclinic_ID")) then
                   rho_bg = dens_env_pp(i, j, k)
                  else
                   if( fluctuationMode ) then
                       rho_bg = 0.0   ! push back to zero perturbation
                      else
                       rho_bg = rhoStrat(k)
                   end if
                end if
                
                rho_old = var(i,j,k,1)
                alpha = spongeAlphaZ*(z(k)-zSponge)/spongeDz
                beta = 1./(1.+alpha*0.5*dt)**2
                rho_new = (1.-beta)*rho_bg + beta*rho_old
                
                var(i,j,k,1) = rho_new
                
             end do
          end do
       end do

       
    case( "ice" ) 
    
       nAer_bg = 0.0 !init_nAer * rhoRef * lRef**3 
       nIce_bg = 0.0
       qIce_bg = 0.0
       qv_bg = 0.0
       
       do k = kSponge, nz
          do j = 1,ny
             do i = 1,nx
             
                select case (iceTestcase) 
                   case ("homogeneous_qv") 
                      qv_bg = 0.0 !init_qv         
                   case ("homogeneous_SIce")
                      !call find_temperature(T,i,j,k,var)
                      !p = press0_dim * ( (PStrat(k)/p0)**gamma_1  +var(i,j,k,5) )**kappaInv
                      qv_bg = 0.0 !epsilon0 * init_SIce * p_saturation(T) / p 
                end select
                
                alpha = spongeAlphaZ*(z(k)-zSponge)/spongeDz
                beta = 1./(1.+alpha*0.5*dt)**2
                var(i,j,k,nVar-3) = (1.-beta)*nAer_bg + beta*var(i,j,k,nVar-3)
                var(i,j,k,nVar-2) = (1.-beta)*nIce_bg + beta*var(i,j,k,nVar-2)
                var(i,j,k,nVar-1) = (1.-beta)*qIce_bg + beta*var(i,j,k,nVar-1)
                var(i,j,k,nVar)   = (1.-beta) * qv_bg + beta*var(i,j,k,nVar)
                do iVar=0,3
                  if (var(i,j,k,nVar-iVar) .lt. 0.0) var(i,j,k,nVar-iVar) = 0.0
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
       
       
    case( "uvw" )
       ! relax u to:
       !   baroclinic cases (2D or 3D):
       !       0 or
       !       environmental u or
       !       u (no relaxation)
       !   else: horizontal mean

       ! local horizontal sum in the sponge layer

       do k = kSponge, nz
          sum_local(k) = sum(var(1:nx,1:ny,k,2))
       end do

       ! global sum and average

       call mpi_allreduce(sum_local(kSponge),sum_global(kSponge),&
                          nz-kSponge+1,&
                          mpi_double_precision,mpi_sum,comm,ierror)
       sum_global = sum_global/(sizeX*sizeY)

       do k = kSponge, nz
          do j = 1,ny
             !do i = 0,nx
             do i = 1,nx
                if ((TestCase == "baroclinic_LC") &
                 & .or.(TestCase == "baroclinic_ID")) then
                   if (Sponge_Rel_Bal_Type == "hyd") then
                      ! relax to hydrost bal
                      uBG = 0.
                     else if (Sponge_Rel_Bal_Type == "env") then
                      ! relax to geostr bal
                      uBG = u_env_pp(i, j, k) 
                     else
                      ! free development
                      uBG = var(i,j,k,2)
                   end if
                  else
                   uBG = sum_global(k)
                end if 
                
                uOld = var(i,j,k,2)
                alpha = spongeAlphaZ * (z(k)-zSponge) / spongeDz
                beta = 1./(1.+alpha*0.5*dt)**2
                
                uNew = (1.-beta)*uBG + beta*uOld
                
                var(i,j,k,2) = uNew
             end do
          end do
       end do
       
       ! relax v to:
       !   baroclinic cases (2D or 3D):
       !       0 or
       !       environmental v or
       !       v (no relaxation)
       !   else: horizontal mean
       
       ! local horizontal sum in the sponge layer

       do k = kSponge, nz
          sum_local(k) = sum(var(1:nx,1:ny,k,3))
       end do

       ! global sum and average

       call mpi_allreduce(sum_local(kSponge),sum_global(kSponge),&
                          nz-kSponge+1,&
                          mpi_double_precision,mpi_sum,comm,ierror)
       sum_global = sum_global/(sizeX*sizeY)

       do k = kSponge, nz
          !do j = 0,ny !gaga 1->0
          do j = 1,ny
             do i = 1,nx
                if ((TestCase == "baroclinic_LC") &
                 & .or.(TestCase == "baroclinic_ID")) then
                   if (Sponge_Rel_Bal_Type == "hyd") then
                      ! relax to hydrost bal
                      vBG = 0.
                     else if (Sponge_Rel_Bal_Type == "env") then
                      ! relax to geostr bal
                      vBG = v_env_pp(i, j, k)
                     else
                      ! free development
                      vBG = var(i,j,k,3)
                   end if
                  else  
                   vBG = sum_global(k)
                end if
                
                vOld = var(i,j,k,3)
                alpha = spongeAlphaZ * (z(k)-zSponge) / spongeDz
                beta = 1./(1.+alpha*0.5*dt)**2
                
                vNew = (1.-beta)*vBG + beta*vOld
                
                var(i,j,k,3) = vNew
                
             end do
          end do
       end do
       
       ! achatzb relax w to zero

       do k = kSponge, nz
          wBG = backgroundFlow_dim(3) / uRef !0.0

          do j = 1,ny
             do i = 1,nx
                
                wOld = var(i,j,k,4)
                alpha = spongeAlphaZ * (z(k)-zSponge) / spongeDz
                beta = 1./(1.+alpha*0.5*dt)**2
                
                wNew = (1.-beta)*wBG + beta*wOld

                var(i,j,k,4) = wNew
                
             end do
          end do
       end do
    case default
       stop"spongeLayer: Unknown variable"
    end select
    

  end subroutine set_spongeLayer


!---------------------------------------------------------------------


  subroutine momentumPredictor (var,flux,flux_rhopw,force,dt,q,m,mmp_mod,int_mod)
    !----------------------------------
    !  calculates the velocities u^*
    !----------------------------------
    
    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var

    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux) 
    ! dir = 1..3 > f,g,h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW
    real, dimension(-1:nx,-1:ny,-1:nz), intent(in) :: flux_rhopw

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
    character(len=*), intent(in) :: mmp_mod, int_mod

    ! volume forces 
    ! mmp_mod = tot =>
    ! 1) gravitational / buoyancy force (cell-centered)
    ! 2) Coriolis force (cell-centered)
    ! 3) WKB wave driving (cell-centered)
    ! mmp_mod = rhs =>
    ! 1) WKB wave driving (cell-centered)
    real, dimension(0:nx+1,0:ny+1,0:nz+1,3), intent(in) :: force

    real, intent(in) :: dt 
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,3), &
         & intent(inout) :: q
    integer, intent(in) :: m

    ! local variables
    real :: fL, fR, gB,gF, hD,hU
    ! flux Left/Right, Backward/Forward, Downward/Upward

    ! usave to keep the new u until v has been updated as well
    ! (for mmp_mod = rhs)
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz) :: usave
    
    ! other stuff
    real :: rhoM_1, rhoM               ! rho(m-1), rho(m)
    real :: drho_e
    real :: uM_1, vM_1, wM_1           ! u(m-1), v(m-1) and w(m-1)
    real :: momM_1, momM               ! momentum at t(m-1) and t(m)
    real :: piR,piL, piF,piB, piU,piD  
    real :: fluxDiff                   ! conv. and viscous flux contr.
    real :: piGrad                     ! pressure gradient
    real :: piGradx, piGrady           ! horizontal pressure-gradient
                                       ! components
    real :: F                          ! update part for Runge-Kutta step
    real :: uAst, vAst, wAst           ! predicted velocities u*, v* and w*
    real :: uhorx, vhory, wvert
    
    integer :: i,j,k
    integer :: i0,i1,j0,j1,k0,k1

    ! gravity
    real :: thetaEdge
    real :: volForce

    ! local interpolation values
    real :: pBarEdge   ! stratified background pressure interpolated
    real :: thetaBar   ! stratified pot. temp. interpolated 
    real :: rhoEdge    ! interpolated density to velocity edge

    ! classical RK3
    real :: rhoM_0, uM_0, vM_0, wM_0, momM_0

    ! non-dimensional Corilois parameter (= inverse Rossby number)
    real :: f_cor_nd

    real :: rho, rhop, rhou, rhov, rhow, facu, facv, facw, facr, pstw, buoy
    real :: rho10, rho01
    real :: rhov0m, rhov00, rhov1m, rhov10
    real :: rhou00, rhoum0, rhou01, rhoum1
    real :: rho000, rho001
    real :: volfcx, volfcy
    real :: bvsstw

    integer :: i00,j00

    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz) :: heat

    real :: rhop_0, rhop_1
    real :: rho_p0, rho_p1

    real, dimension(-nbz:nz+nbz) :: w_0 
    real, dimension(-nbz:nz+nbz) :: S_bar 
    real :: heat_flc, heat0, heat1

    ! non-dimensional Corilois parameter (= inverse Rossby number)
    f_cor_nd = f_Coriolis_dim*tRef

    if(topography) then
       i00=is+nbx-1
       j00=js+nby-1
    end if

    if (int_mod == "impl") then
       ! environmental heating
       call calculate_heating(var,flux,heat)

       ! heating by GW entropy-flux convergence
       if (raytracer) heat(:,:,:) = heat(:,:,:) + var(:,:,:,8)
    end if

    if( correctDivError ) then
        print*,'ERROR: correction divergence error not allowed'
        stop
    end if

    ! init q
    if (m == 1) q = 0.

    ! -------------------------------------
    !            predict u -> u*
    ! -------------------------------------

    select case( xBoundary ) 
       
    case( "solid_wall" )
       i0 = 1
       i1 = nx-1
    case( "periodic" ) 
       i0 = 0
       i1 = nx
    case default
       stop"momentumPredictor: unknown case xBoundary."
    end select

    if (mmp_mod == "tot" .or. mmp_mod == "lhs") then
       if (int_mod /= "expl") then
          stop'ERROR: wrong int_mod for mmp_mod = tot or mmp_mod = lhs'
       end if

       do k = 1,nz
          do j = 1,ny
             do i = i0,i1
             
                !--- convective fluxes -> conv
                fR = flux(i  ,j,k,1,2)
                fL = flux(i-1,j,k,1,2)
                gF = flux(i,j  ,k,2,2)
                gB = flux(i,j-1,k,2,2)
                hU = flux(i,j,k  ,3,2)
                hD = flux(i,j,k-1,3,2)
                fluxDiff = (fR-fL)/dx + (gF-gB)/dy + (hU-hD)/dz ! diverg.

                if (mmp_mod == "tot") then
                   !--- pressure gradient term -> piGrad
                   if (TestCase == "baroclinic_LC") then
                      piR = var(i+1,j,k,5) - var_env(i+1,j,k,5)
                      piL = var(i,  j,k,5) - var_env(i,  j,k,5)
                     else
                      piR = var(i+1,j,k,5)
                      piL = var(i,  j,k,5)
                   end if

                   piGrad = kappaInv*MaInv2 * Pstrat(k) * (piR-piL)/dx

                   !---- volume forces
                   volForce = 0.5*( force(i,j,k,1) + force(i+1,j,k,1) )

                   if (TestCase == "baroclinic_LC") then
                      if (model == "pseudo_incompressible") then
                          rhoM_1 = 0.5 * (rhoOld(i,j,k) + rhoOld(i+1,j,k))

                          if( fluctuationMode ) then
                             rhoM_1 = rhoM_1 + rhoStrat(k)
                          end if
                         else if (model == "Boussinesq") then
                          rhoM_1 = rho00
                         else
                          stop"momentumPredictor: unkown model."
                      end if

                      volforce &
                      =   volforce &
                        - rhoM_1*RoInv &
                          * 0.25 &
                          * (  var_env(i,j-1,k,3) + var_env(i+1,j-1,k,3) &
                             + var_env(i,j  ,k,3) + var_env(i+1,j  ,k,3))

                      !UAB
                      if (background == "HeldSuarez") then
                         ! Rayleigh damping in boundary layer

                         if (model == "pseudo_incompressible") then
                             rhoM_1 &
                             = 0.5 * (rhoOld(i,j,k) + rhoOld(i+1,j,k))

                             if( fluctuationMode ) then
                                rhoM_1 = rhoM_1 + rhoStrat(k)
                             end if
                            else if (model == "Boussinesq") then
                             rhoM_1 = rho00
                            else
                             stop"momentumPredictor: unkown model."
                         end if

                         volForce = volForce - kv_hs(k)*rhoM_1*(var(i,j,k,2)-u_env_pp(i,j,k)) !FS
                      end if
                      !UAE
                   end if

                   if (topography) then
                      ! Rayleigh damping in land cells (immersed boundary)

                      if (model == "pseudo_incompressible") then
                          rhoM_1 = 0.5 * (rhoOld(i,j,k) + rhoOld(i+1,j,k))

                          if( fluctuationMode ) then
                             rhoM_1 = rhoM_1 + rhoStrat(k)
                          end if
                         else if (model == "Boussinesq") then
                          rhoM_1 = rho00
                         else
                          stop"momentumPredictor: unkown model."
                      end if

                      if(topography_mask(i00+i,j00+j,k)&
                         .or.&
                         topography_mask(i00+i+1,j00+j,k)) then
                         volForce = volForce - alprlx * rhoM_1*var(i,j,k,2)
                      end if
                   end if
                end if
        
                !--------------------
                !   d/dt ... = F(phi) (RHS of ODE)
                !--------------------
                ! fluxDiff -> convective and viscous fluxes
                ! piGrad   -> pressure gradient along x scaled with 1/Ma^2
                ! volForce -> Gravity, Coriolis
                if (mmp_mod == "tot") then
                   F = -fluxDiff - piGrad + volForce
                  else if (mmp_mod == "lhs") then
                   F = -fluxDiff
                  else
                   stop'ERROR: wrong mmp_mod'
                end if
             
                ! interpolated density
                select case( model ) 

                case( "pseudo_incompressible" )

                   rhoM_1 = 0.5 * ( rhoOld(i,j,k) + rhoOld(i+1,j,k) )
                   rhoM   = 0.5 * ( var(i,j,k,1) + var(i+1,j,k,1) )
                
                   if( fluctuationMode ) then
                      rhoM_1 = rhoM_1 + rhoStrat(k)
                      rhoM   = rhoM   + rhoStrat(k)
                   end if
                
                case( "Boussinesq" )
                   rhoM_1 = rho00
                   rhoM = rho00
                case default
                   stop"momentumPredictor: unkown case model."
                end select
             
                ! velocity and momentum at t(m-1)
                uM_1 = var(i,j,k,2)
                momM_1 = rhoM_1 * uM_1
             
                ! q(m-1) -> q(m)
             
                q(i,j,k,1) = dt * F + alpha(m) * q(i,j,k,1)

                ! rhoU(m-1) -> rhoU(m)
                momM = momM_1 + beta(m) * q(i,j,k,1)

                ! calc u(m,*)
                uAst = momM / rhoM
   
                ! uAst -> var
                var(i,j,k,2) = uAst
             end do
          end do
       end do
      else if (mmp_mod == "rhs") then
       if (int_mod == "expl") then
          do k = 1,nz
             do j = 1,ny
                do i = i0,i1
                   rhou = 0.5 * (var(i,j,k,1) + var(i+1,j,k,1))
                   if( fluctuationMode ) then
                       rhou =  rhou + rhoStrat(k)
                   end if

                   !--- pressure gradient term -> piGrad
                   if (TestCase == "baroclinic_LC") then
                      piR = var(i+1,j,k,5) - var_env(i+1,j,k,5)
                      piL = var(i,  j,k,5) - var_env(i,  j,k,5)
                     else
                      piR = var(i+1,j,k,5)
                      piL = var(i,  j,k,5)
                   end if

                   piGrad = kappaInv*MaInv2 * Pstrat(k)/rhou * (piR-piL)/dx

                   ! gravity-wave forcing
                   if (raytracer .or. (testCase == "mountainwave")) then
                      volfcx = 0.5 * (force(i,j,k,1) + force(i+1,j,k,1))
                     else
                      volfcx = 0.0
                   end if

                   ! ustar
                   if (TestCase == "baroclinic_LC") then
                      uhorx = var(i,j,k,2) - var_env(i,j,k,2)
                     else
                      uhorx = var(i,j,k,2)
                   end if

                   vhory &
                   = 0.25 * (  var(i  ,j-1,k,3) + var(i  ,j,k,3) &
                             + var(i+1,j-1,k,3) + var(i+1,j,k,3))

                   uAst &
                   = uhorx + dt*(f_cor_nd*vhory - piGrad + volfcx/rhou)

                   if (topography) then
                      ! Rayleigh damping in land cells (immersed boundary)
                      if(topography_mask(i00+i,j00+j,k)&
                         .or.&
                         topography_mask(i00+i+1,j00+j,k)) then
                         uAst = uAst - dt* alprlx*uhorx
                      end if
                   end if

                   usave(i,j,k) = uAst

                   if (TestCase == "baroclinic_LC") then
                      usave(i,j,k) = usave(i,j,k) + var_env(i,j,k,2)
                   end if
                end do
             end do
          end do
         else if (int_mod == "impl") then
          do k = 1,nz
             do j = 1,ny
                do i = i0,i1
                   rhou = 0.5 * (var(i,j,k,1) + var(i+1,j,k,1))

                   rhov0m =  0.5 * (var(i  ,j  ,k,1) + var(i  ,j-1,k,1))
                   rhov00 =  0.5 * (var(i  ,j+1,k,1) + var(i  ,j  ,k,1))
                   rhov1m =  0.5 * (var(i+1,j  ,k,1) + var(i+1,j-1,k,1))
                   rhov10 =  0.5 * (var(i+1,j+1,k,1) + var(i+1,j  ,k,1))

                   if( fluctuationMode ) then
                         rhou  =    rhou + rhoStrat(k)
                       rhov0m  =  rhov0m + rhoStrat(k)
                       rhov00  =  rhov00 + rhoStrat(k)
                       rhov1m  =  rhov1m + rhoStrat(k)
                       rhov10  =  rhov10 + rhoStrat(k)
                   end if

                   !--- pressure gradient terms -> piGradx, piGrady
                   if (TestCase == "baroclinic_LC") then
                      piR = var(i+1,j,k,5) - var_env(i+1,j,k,5)
                      piL = var(i,  j,k,5) - var_env(i,  j,k,5)
                     else
                      piR = var(i+1,j,k,5)
                      piL = var(i,  j,k,5)
                   end if

                   piGradx &
                   = kappaInv*MaInv2 * Pstrat(k)/rhou &
                     * (var(i+1,j,k,5) - var(i,j,k,5))/dx

                   if (TestCase == "baroclinic_LC") then
                      piGrady &
                      = kappaInv*MaInv2 &
                        * 0.25 &
                        * (  Pstrat(k)/rhov0m &
                             * (  var(i  ,j  ,k,5) - var(i  ,j-1,k,5) &
                                - var_env(i  ,j  ,k,5) &
                                + var_env(i  ,j-1,k,5))/dy &
                           + Pstrat(k)/rhov00 &
                             * (  var(i  ,j+1,k,5) - var(i  ,j  ,k,5) &
                                - var_env(i  ,j+1,k,5) &
                                + var_env(i  ,j  ,k,5))/dy &
                           + Pstrat(k)/rhov1m &
                             * (  var(i+1,j  ,k,5) - var(i+1,j-1,k,5) &
                                - var_env(i+1,j  ,k,5) &
                                + var_env(i+1,j-1,k,5))/dy &
                           + Pstrat(k)/rhov10 &
                             * (  var(i+1,j+1,k,5) - var(i+1,j  ,k,5) &
                                - var_env(i+1,j+1,k,5) &
                                + var_env(i+1,j  ,k,5))/dy)
                     else
                      piGrady &
                      = kappaInv*MaInv2 &
                        * 0.25 &
                        * (  Pstrat(k)/rhov0m &
                             * (var(i  ,j  ,k,5) - var(i  ,j-1,k,5))/dy &
                           + Pstrat(k)/rhov00 &
                             * (var(i  ,j+1,k,5) - var(i  ,j  ,k,5))/dy &
                           + Pstrat(k)/rhov1m &
                             * (var(i+1,j  ,k,5) - var(i+1,j-1,k,5))/dy &
                           + Pstrat(k)/rhov10 &
                             * (var(i+1,j+1,k,5) - var(i+1,j  ,k,5))/dy)
                   end if

                   ! gravity-wave forcing
                   if (raytracer .or. (testCase == "mountainwave")) then
                      volfcx = 0.5 * (force(i,j,k,1) + force(i+1,j,k,1))
                      volfcy = 0.5 * (force(i,j,k,2) + force(i,j+1,k,2))
                     else
                      volfcx = 0.0
                      volfcy = 0.0
                   end if

                   ! ustar
                   if (TestCase == "baroclinic_LC") then
                      uhorx = var(i,j,k,2) - var_env(i,j,k,2)
                     else
                      uhorx = var(i,j,k,2)
                   end if

                   vhory &
                   = 0.25 * (  var(i  ,j-1,k,3) + var(i  ,j,k,3) &
                             + var(i+1,j-1,k,3) + var(i+1,j,k,3))

                   facu = 1.0

                   if (topography) then
                      ! Rayleigh damping in land cells (immersed boundary)

                      if(topography_mask(i00+i,j00+j,k)&
                         .or.&
                         topography_mask(i00+i+1,j00+j,k)) then
                         facu = facu + dt*alprlx
                      end if
                   end if

                   facv = facu

                   uAst &
                   = 1.0/(facu*facv + (f_cor_nd*dt)**2) &
                     * (  facv * (uhorx + dt*(volfcx/rhou - piGradx)) &
                        + f_cor_nd*dt &
                          * (vhory + dt*(volfcy/rhou - piGrady)))

                   usave(i,j,k) = uAst

                   if (TestCase == "baroclinic_LC") then
                      usave(i,j,k) = usave(i,j,k) + var_env(i,j,k,2)
                   end if
                end do
             end do
          end do
         else
          stop'ERROR: unknown int_mod'
       end if
      else
       stop'ERROR: unknown mmp_mod'
    end if

    ! -------------------------------------
    !            predict v -> v*
    ! -------------------------------------

    select case( yBoundary ) 
       
    case( "solid_wall" )
       j0 = 1
       j1 = ny-1
    case( "periodic" ) 
       j0 = 0
       j1 = ny
    case default
       stop"momentumPredictor: unknown case yBoundary."
    end select
    
    if (mmp_mod == "tot" .or. mmp_mod == "lhs") then
       if (int_mod /= "expl") then
          stop'ERROR: wrong int_mod for mmp_mod = tot or mmp_mod = lhs'
       end if

       do k = 1,nz
          do j = j0,j1
             do i = 1,nx

                !--- convective part -> conv
                fR = flux(i,  j,k,1,3)
                fL = flux(i-1,j,k,1,3)
                gF = flux(i,j,  k,2,3)
                gB = flux(i,j-1,k,2,3)
                hU = flux(i,j,k,  3,3)
                hD = flux(i,j,k-1,3,3)
                fluxDiff = (fR-fL)/dx + (gF-gB)/dy + (hU-hD)/dz

                if (mmp_mod == "tot") then
                   !--- pressure gradient term -> piGrad
                   if (TestCase == "baroclinic_LC") then
                      piF = var(i,j+1,k,5) - var_env(i,j+1,k,5)
                      piB = var(i,j,  k,5) - var_env(i,j,  k,5)
                     else
                      piF = var(i,j+1,k,5)
                      piB = var(i,j,  k,5)
                   end if

                   piGrad = kappaInv*MaInv2 * pStrat(k) * (piF-piB)/dy

                   !---- volume forces
                   volForce = 0.5*( force(i,j,k,2) + force(i,j+1,k,2) )

                   if (TestCase == "baroclinic_LC") then
                      if (model == "pseudo_incompressible") then
                          rhoM   = rhoOld(i,j  ,k)
                          rhoM_1 = rhoOld(i,j+1,k)

                          if( fluctuationMode ) then
                             rhoM   = rhoM   + rhoStrat(k)
                             rhoM_1 = rhoM_1 + rhoStrat(k)
                          end if
                         else if (model == "Boussinesq") then
                          rhoM   = rho00
                          rhoM_1 = rho00
                         else
                          stop"momentumPredictor: unkown model."
                      end if

                      volforce &
                      =   volforce &
                        + RoInv &
                          * 0.25 &
                          * (  rhoM &
                               * (  var_env(i-1,j,k,2) &
                                  + var_env(i  ,j,k,2)) &
                             + rhoM_1 &
                               * (  var_env(i-1,j+1,k,2) &
                                  + var_env(i  ,j+1,k,2)))

                      !UAB
                      if (background == "HeldSuarez") then
                         ! Rayleigh damping in boundary layer

                         if (model == "pseudo_incompressible") then
                             rhoM_1 &
                             = 0.5 * (rhoOld(i,j,k) + rhoOld(i,j+1,k))

                             if( fluctuationMode ) then
                                rhoM_1 = rhoM_1 + rhoStrat(k)
                             end if
                            else if (model == "Boussinesq") then
                             rhoM_1 = rho00
                            else
                             stop"momentumPredictor: unkown model."
                         end if

                         volForce = volForce - kv_hs(k)*rhoM_1*var(i,j,k,3)!FS
                      end if
                      !UAE
                   end if

                   if (topography) then
                      ! Rayleigh damping in land cells (immersed boundary)

                      if (model == "pseudo_incompressible") then
                          rhoM_1 = 0.5 * (rhoOld(i,j,k) + rhoOld(i,j+1,k))

                          if( fluctuationMode ) then
                             rhoM_1 = rhoM_1 + rhoStrat(k)
                          end if
                         else if (model == "Boussinesq") then
                          rhoM_1 = rho00
                         else
                          stop"momentumPredictor: unkown model."
                      end if

                      if(topography_mask(i00+i,j00+j,k)&
                         .or.&
                         topography_mask(i00+i,j00+j+1,k)) then
                         volForce = volForce - alprlx * rhoM_1*var(i,j,k,3)
                      end if
                   end if
                end if
             
                !--------------------
                !   F(phi) = RHS
                !--------------------
                ! fluxDiff -> convective and viscous fluxes
                ! piGrad   -> pressure gradient along x
                ! volForce -> Gravity, Coriolis
                if (mmp_mod == "tot") then
                   F = -fluxDiff - piGrad + volForce
                  else if (mmp_mod == "lhs") then
                   F = -fluxDiff
                  else
                   stop'ERROR: wrong mmp_mod'
                end if

                ! interpolated density
                select case( model ) 
             
                case( "pseudo_incompressible" )

                   rhoM_1 = 0.5*( rhoOld(i,j,k) + rhoOld(i,j+1,k) )
                   rhoM   = 0.5*( var(i,j,k,1) + var(i,j+1,k,1) )

                   if( fluctuationMode ) then
                      rhoM_1 = rhoM_1 + rhoStrat(k)
                      rhoM   = rhoM   + rhoStrat(k)
                   end if
                
                case( "Boussinesq" )
                   rhoM_1 = rho00
                   rhoM = rho00
                case default
                   stop"momentumPredictor: unkown case model."
                end select

                ! velocity and momentum at t(m-1)
                vM_1 = var(i,j,k,3)
                momM_1 = rhoM_1 * vM_1

                ! q(m-1) -> q(m)
                q(i,j,k,2) = dt * F + alpha(m) * q(i,j,k,2)

                ! rhoV(m-1) -> rhoV(m)
                momM = momM_1 + beta(m) * q(i,j,k,2)

                ! calc v(m,*)
                vAst = momM / rhoM

                ! vAst -> var
                var(i,j,k,3) = vAst
             end do
          end do
       end do
      else if (mmp_mod == "rhs") then
       if (int_mod == "expl") then
          do k = 1,nz
             do j = j0,j1
                do i = 1,nx
                   rhov = 0.5 * (var(i,j,k,1) + var(i,j+1,k,1))
                   if( fluctuationMode ) then
                       rhov =  rhov + rhoStrat(k)
                   end if

                   !--- pressure gradient term -> piGrad
                   if (TestCase == "baroclinic_LC") then
                      piF = var(i,j+1,k,5) - var_env(i,j+1,k,5)
                      piB = var(i,  j,k,5) - var_env(i,  j,k,5)
                     else
                      piF = var(i,j+1,k,5)
                      piB = var(i,  j,k,5)
                   end if

                   piGrad = kappaInv*MaInv2 * Pstrat(k)/rhov * (piF-piB)/dy

                   ! gravity-wave forcing
                   if (raytracer .or. (testCase == "mountainwave")) then
                      volfcy = 0.5 * (force(i,j,k,2) + force(i,j+1,k,2))
                     else
                      volfcy = 0.0
                   end if

                   ! vstar
                   if (TestCase == "baroclinic_LC") then
                      uhorx &
                      = 0.25 &
                        * (  var(i-1,j,k,2) + var(i-1,j+1,k,2) &
                           - var_env(i-1,j,k,2) - var_env(i-1,j+1,k,2) &
                           + var(i  ,j,k,2) + var(i  ,j+1,k,2) &
                           - var_env(i  ,j,k,2) - var_env(i  ,j+1,k,2))
                     else
                      uhorx &
                      = 0.25 * (  var(i-1,j,k,2) + var(i-1,j+1,k,2) &
                                + var(i  ,j,k,2) + var(i  ,j+1,k,2))
                   end if

                   vhory = var(i,j,k,3)

                   vAst &
                   = vhory + dt*(-f_cor_nd*uhorx - piGrad + volfcy/rhov)

                   if (topography) then
                      ! Rayleigh damping in land cells (immersed boundary)
                      if(topography_mask(i00+i,j00+j,k)&
                         .or.&
                         topography_mask(i00+i,j00+j+1,k)) then
                         vAst = vAst - dt* alprlx*vhory
                      end if
                   end if

                   var(i,j,k,3) = vAst
                end do
             end do
          end do
         else if (int_mod == "impl") then
          do k = 1,nz
             do j = j0,j1
                do i = 1,nx
                   rhoum0 = 0.5 * (var(i  ,j  ,k,1) + var(i-1,j  ,k,1))
                   rhou00 = 0.5 * (var(i+1,j  ,k,1) + var(i  ,j  ,k,1))
                   rhoum1 = 0.5 * (var(i  ,j+1,k,1) + var(i-1,j+1,k,1))
                   rhou01 = 0.5 * (var(i+1,j+1,k,1) + var(i  ,j+1,k,1))

                    rhov = 0.5 * (var(i,j,k,1) + var(i,j+1,k,1))

                   if( fluctuationMode ) then
                         rhov  =    rhov  + rhoStrat(k)
                       rhoum0  =  rhoum0  + rhoStrat(k)
                       rhou00  =  rhou00  + rhoStrat(k)
                       rhoum1  =  rhoum1  + rhoStrat(k)
                       rhou01  =  rhou01  + rhoStrat(k)
                   end if

                   !--- pressure gradient terms -> piGradx, piGrady
                   if (TestCase == "baroclinic_LC") then
                      piGradx &
                      = kappaInv*MaInv2 &
                        * 0.25 &
                        * (  Pstrat(k)/rhou00 &
                             * (  var(i+1,j  ,k,5) - var(i  ,j  ,k,5) &
                                - var_env(i+1,j  ,k,5) &
                                + var_env(i  ,j  ,k,5))/dx &
                           + Pstrat(k)/rhoum0 &
                             * (  var(i  ,j  ,k,5) - var(i-1,j  ,k,5) &
                                - var_env(i  ,j  ,k,5) &
                                + var_env(i-1,j  ,k,5))/dx &
                           + Pstrat(k)/rhou01 &
                             * (  var(i+1,j+1,k,5) - var(i  ,j+1,k,5) &
                                - var_env(i+1,j+1,k,5) &
                                + var_env(i  ,j+1,k,5))/dx &
                           + Pstrat(k)/rhoum1 &
                             * (  var(i  ,j+1,k,5) - var(i-1,j+1,k,5) &
                                - var_env(i  ,j+1,k,5) &
                                + var_env(i-1,j+1,k,5))/dx)
                     else
                      piGradx &
                      = kappaInv*MaInv2 &
                        * 0.25 &
                        * (  Pstrat(k)/rhou00 &
                             * (var(i+1,j  ,k,5) - var(i  ,j  ,k,5))/dx &
                           + Pstrat(k)/rhoum0 &
                             * (var(i  ,j  ,k,5) - var(i-1,j  ,k,5))/dx &
                           + Pstrat(k)/rhou01 &
                             * (var(i+1,j+1,k,5) - var(i  ,j+1,k,5))/dx &
                           + Pstrat(k)/rhoum1 &
                             * (var(i  ,j+1,k,5) - var(i-1,j+1,k,5))/dx)
                   end if

                   if (TestCase == "baroclinic_LC") then
                      piF = var(i,j+1,k,5) - var_env(i,j+1,k,5)
                      piB = var(i,j  ,k,5) - var_env(i,j  ,k,5)
                     else
                      piF = var(i,j+1,k,5)
                      piB = var(i,j  ,k,5)
                   end if

                   piGrady &
                   = kappaInv*MaInv2 * Pstrat(k)/rhov * (piF-piB)/dy

                   ! gravity-wave forcing
                   if (raytracer .or. (testCase == "mountainwave")) then
                      volfcx = 0.5 * (force(i,j,k,1) + force(i+1,j,k,1))
                      volfcy = 0.5 * (force(i,j,k,2) + force(i,j+1,k,2))
                     else
                      volfcx = 0.0
                      volfcy = 0.0
                   end if

                   ! vstar
                   if (TestCase == "baroclinic_LC") then
                      uhorx &
                      = 0.25 &
                        * (  var(i-1,j,k,2) + var(i-1,j+1,k,2) &
                           - var_env(i-1,j,k,2) - var_env(i-1,j+1,k,2) &
                           + var(i  ,j,k,2) + var(i  ,j+1,k,2) &
                           - var_env(i  ,j,k,2) - var_env(i  ,j+1,k,2))
                     else
                      uhorx &
                      = 0.25 * (  var(i-1,j,k,2) + var(i-1,j+1,k,2) &
                                + var(i  ,j,k,2) + var(i  ,j+1,k,2))
                   end if

                   vhory = var(i,j,k,3)

                   facv = 1.0

                   if (topography) then
                      ! Rayleigh damping in land cells (immersed boundary)

                      if(topography_mask(i00+i,j00+j,k)&
                         .or.&
                         topography_mask(i00+i,j00+j+1,k)) then
                         facv = facv + dt*alprlx
                      end if
                   end if

                   facu = facv

                   vAst &
                   = 1.0/(facu*facv + (f_cor_nd*dt)**2) &
                     * (- f_cor_nd*dt &
                          * (uhorx + dt * (volfcx/rhov - piGradx)) &
                        + facu * (vhory + dt * (volfcy/rhov - piGrady)))

                   var(i,j,k,3) = vAst
                end do
             end do
          end do
         else
          stop'ERROR: unknown int_mod'
       end if

       ! now the new u can be put into the proper array
       var(:,:,:,2) = usave(:,:,:)
      else
       stop'ERROR: unknown mmp_mod'
    end if

    !testb
    !write(42) var
    !stop
    !teste


    ! -------------------------------------
    !            predict w -> w*
    ! -------------------------------------
    
    select case( zBoundary ) 
       
    case( "solid_wall" )
       k0 = 1
       k1 = nz-1
    case( "periodic" ) 
       k0 = 0
       k1 = nz
    case default
       stop"momentumPredictor: unknown case zBoundary."
    end select

    if (mmp_mod == "tot" .or. mmp_mod == "lhs") then
       if (int_mod /= "expl") then
          stop'ERROR: wrong int_mod for mmp_mod = tot or mmp_mod = lhs'
       end if

       do k = k0,k1
          do j = 1,ny
             do i = 1,nx

                !--- convective part -> conv
                fR = flux(i,  j,k,1,4)
                fL = flux(i-1,j,k,1,4)
                gF = flux(i,j,  k,2,4)
                gB = flux(i,j-1,k,2,4)
                hU = flux(i,j,k,  3,4)
                hD = flux(i,j,k-1,3,4)
                fluxDiff = (fR-fL)/dx + (gF-gB)/dy + (hU-hD)/dz

                if (mmp_mod == "tot") then
                   !--- pressure gradient term -> piGrad
                   if (TestCase == "baroclinic_LC") then
                      piU = var(i,j,k+1,5) - var_env(i,j,k+1,5)
                      piD = var(i,j,  k,5) - var_env(i,j,  k,5)
                     else
                      piU = var(i,j,k+1,5)
                      piD = var(i,j,  k,5)
                   end if

                   piGrad = 0.5 * kappaInv*MaInv2* &
                        & ( Pstrat(k)+Pstrat(k+1) ) * (piU-piD)/dz
   
                   !---- volume forces
                   volForce = 0.5*( force(i,j,k,3) + force(i,j,k+1,3) )

                   if (TestCase == "baroclinic_LC") then
                      if (model == "pseudo_incompressible") then
                          drho_e &
                          = 0.5 * (var_env(i,j,k,1) + var_env(i,j,k+1,1))

                          if( .not. fluctuationMode ) then
                             drho_e = drho_e - rhoStratTilde(k)
                          end if
                         else if (model == "Boussinesq") then
                          stop'ERROR: baroclinic LC not ready yet for &
                                  & Boussinesq'
                         else
                          stop"momentumPredictor: unkown model."
                      end if

                      volForce = volForce + FrInv2 * drho_e
                   end if

                   if (topography) then
                      ! Rayleigh damping in land cells (immersed boundary)

                      if (model == "pseudo_incompressible") then
                          rhoM_1 = 0.5 * (rhoOld(i,j,k) + rhoOld(i,j,k+1))

                          if( fluctuationMode ) then
                             rhoM_1 = rhoM_1 + rhoStratTilde(k)
                          end if
                         else if (model == "Boussinesq") then
                          rhoM_1 = rho00
                         else
                          stop"momentumPredictor: unkown model."
                      end if

                      if(topography_mask(i00+i,j00+j,k)&
                         .or.&
                         topography_mask(i00+i,j00+j,k+1)) then
                         volForce = volForce - alprlx * rhoM_1*var(i,j,k,4)
                      end if
                   end if
                end if

                !--------------------
                !   F(phi) = RHS
                !--------------------
                ! fluxDiff -> convective and viscous fluxes
                ! piGrad   -> pressure gradient along x
                ! volForce -> Gravity, Coriolis
                if (mmp_mod == "tot") then
                   F = -fluxDiff - piGrad + volForce
                  else if (mmp_mod == "lhs") then
                   F = -fluxDiff
                  else
                   stop'ERROR: wrong mmp_mod'
                end if

                ! interpolated densities
                select case( model ) 
                
                case( "pseudo_incompressible" )

                   rhoM_1 = 0.5*(rhoOld(i,j,k) + rhoOld(i,j,k+1)) !rho(m-1)
                   rhoM   = 0.5*(var(i,j,k,1) + var(i,j,k+1,1))   !rho(m)

                   if( fluctuationMode ) then
                      rhoM_1 = rhoM_1 + rhoStratTilde(k)
                      rhoM   = rhoM   + rhoStratTilde(k)
                   end if
                
                case( "Boussinesq" )
                   rhoM_1 = rho00
                   rhoM = rho00
                case default
                   stop"momentumPredictor: unkown case model."
                end select
             
                ! velocity and momentum at t(m-1)
                wM_1 = var(i,j,k,4)
                momM_1 = rhoM_1 * wM_1

                ! q(m-1) -> q(m)
                q(i,j,k,3) = dt * F + alpha(m) * q(i,j,k,3)

                ! rhoW(m-1) -> rhoW(m)
                momM = momM_1 + beta(m) * q(i,j,k,3)

                ! calc w(m,*)
                wAst = momM / rhoM

                ! wAst -> var
                var(i,j,k,4) = wAst
             end do
          end do
       end do
      else if (mmp_mod == "rhs") then
       if (int_mod == "expl") then
          do k = k0,k1
             pstw = 0.5*(Pstrat(k) + Pstrat(k+1))

             do j = 1,ny
                do i = 1,nx
                   rho000 = var(i,j,k  ,1)
                   rho001 = var(i,j,k+1,1)

                     rhow = 0.5 * (var(i,j,k,1) + var(i,j,k+1,1))

                   if( fluctuationMode ) then
                      rho000 = rho000 + rhoStrat(k  )
                      rho001 = rho001 + rhoStrat(k+1)

                       rhow =   rhow + rhoStratTilde(k)
                   end if

                   !--- pressure gradient term -> piGrad
                   if (TestCase == "baroclinic_LC") then
                      piU = var(i,j,k+1,5) - var_env(i,j,k+1,5)
                      piD = var(i,  j,k,5) - var_env(i,  j,k,5)
                     else
                      piU = var(i,j,k+1,5)
                      piD = var(i,  j,k,5)
                   end if

                   piGrad = kappaInv*MaInv2 * pstw/rhow * (piU-piD)/dz

                   ! wstar
                   wvert = var(i,j,k,4)

                   if (TestCase == "baroclinic_LC") then
                      buoy  &
                      = - g_ndim &
                          * 0.5 &
                          * (  (rhopOld(i,j,k) - var_env(i,j,k,6))/rho000 &
                             + (  rhopOld(i,j,k+1) &
                                - var_env(i,j,k+1,6))/rho001)
                     else
                      buoy  &
                      = - g_ndim &
                          * 0.5*(  rhopOld(i,j,k  )/rho000 &
                                 + rhopOld(i,j,k+1)/rho001)
                   end if

                   wAst = wvert + dt*(buoy - piGrad)

                   if (topography) then
                      ! Rayleigh damping in land cells (immersed boundary)
                      if(topography_mask(i00+i,j00+j,k)&
                         .or.&
                         topography_mask(i00+i,j00+j,k+1)) then
                         wAst = wAst - dt* alprlx*wvert
                      end if
                   end if

                   var(i,j,k,4) = wAst
                end do
             end do
          end do
         else if (int_mod == "impl") then
          !UAB
          ! heating due to relaxation, entropy diffusion and GWs, its 
          ! horizontal mean and the horizontal-mean vertical wind 
          ! resulting from it

          if (heatingONK14 .or. TurbScheme .or. rayTracer) then
             call heat_w0(var,flux,flux_rhopw,heat,S_bar,w_0)
            else
             heat = 0.
             S_bar = 0.
             w_0 = 0.
          end if
          !UAE

          do k = k0,k1
             pstw = 0.5*(Pstrat(k) + Pstrat(k+1))

             do j = 1,ny
                do i = 1,nx
                   rho000 = var(i,j,k  ,1)
                   rho001 = var(i,j,k+1,1)

                   rhow = 0.5 * (var(i,j,k,1) + var(i,j,k+1,1))

                   if( fluctuationMode ) then
                      rho000 = rho000 + rhoStrat(k  )
                      rho001 = rho001 + rhoStrat(k+1)

                        rhow =  rhow + rhoStratTilde(k)
                   end if

                   !--- pressure gradient term -> piGrad
                   if (TestCase == "baroclinic_LC") then
                      piU = var(i,j,k+1,5) - var_env(i,j,k+1,5)
                      piD = var(i,  j,k,5) - var_env(i,  j,k,5)
                     else
                      piU = var(i,j,k+1,5)
                      piD = var(i,  j,k,5)
                   end if

                   piGrad = kappaInv*MaInv2 * pstw/rhow * (piU-piD)/dz

                   ! wstar
                   wvert = var(i,j,k,4)

                   ! ssquared Brunt-Vaisala frequency averaged to half 
                   ! levels
                   ! (could be done a bit nicer by determining this without
                   ! averaging directly from the reference-atmosphere
                   ! density)
                   bvsstw = 0.5 * (bvsStrat(k) + bvsStrat(k+1))

                   facw = 1.0
                   facr = 1.0 + alprlx*dt

                   if (topography) then
                      ! Rayleigh damping in land cells (immersed boundary)
                      if(topography_mask(i00+i,j00+j,k)&
                         .or.&
                         topography_mask(i00+i,j00+j,k+1)) then
                         facw = facw + alprlx*dt
                      end if
                   end if

                   if (TestCase == "baroclinic_LC") then
                      rhop_0 = rhopOld(i,j,k  ) - var_env(i,j,k  ,6)
                      rhop_1 = rhopOld(i,j,k+1) - var_env(i,j,k+1,6)

                      rho_p0 = rho000 - var_env(i,j,k  ,6) - rhoStrat(k  )
                      rho_p1 = rho001 - var_env(i,j,k+1,6) - rhoStrat(k+1)
                     else
                      rhop_0 = rhopOld(i,j,k  )
                      rhop_1 = rhopOld(i,j,k+1)

                      rho_p0 = rho000 - rhoStrat(k  )
                      rho_p1 = rho001 - rhoStrat(k+1)
                   end if

                   !UAB
                   !if (heatingONK14)then
                   !wAst &
                   != 1.0 &
                   !  /(  facw*facr &
                   !    + rhoStratTilde(k)/rhow * bvsstw * dt**2) &
                   !  * (  facr * (wvert - dt * piGrad) &
                   !     - dt * g_ndim &
                   !       * 0.5*(  rhop_0/rho000 + rhop_1/rho001&
                   !              + dt &
                   !                * (  rhoStrat(k)/Pstrat(k) &
                   !                     * heat(i,j,k)/rho000 &
                   !                   + rhoStrat(k+1)/Pstrat(k+1) &
                   !                     *  heat(i,j,k+1)/rho001 &
                   !                   + 0.5/g_ndim &
                   !                     *(  rhoStrat(k)*bvsStrat(k) &
                   !                         *(w_0(k)+w_0(k-1))/rho000 &
                   !                       + rhoStrat(k+1)*bvsStrat(k+1) &
                   !                         *(w_0(k)+w_0(k+1))/rho001) &
                   !                   + alprlx &
                   !                     * (  rho_p0/rho000 &
                   !                        + rho_p1/rho001)))) 
                   !else
                   !   wAst &
                   !        = 1.0 &
                   !  /(  facw*facr &
                   !    + rhoStratTilde(k)/rhow * bvsstw * dt**2) &
                   !  * (  facr * (wvert - dt * piGrad) &
                   !     - dt * g_ndim &
                   !       * 0.5*(  rhop_0/rho000 + rhop_1/rho001&
                   !              + dt &
                   !                * (  rhoStrat(k)/Pstrat(k) &
                   !                     * heat(i,j,k)/rho000 &
                   !                   + rhoStrat(k+1)/Pstrat(k+1) &
                   !                     * heat(i,j,k+1)/rho001 &
                   !                   + alprlx &
                   !                     * (  rho_p0/rho000 &
                   !                        + rho_p1/rho001))))
                   !end if

                   heat0 &
                   = heat(i,j,k) - S_bar(k) &
                     + g_ndim/Pstrat(k) * bvsStrat(k) &
                       * 0.5*(w_0(k) + w_0(k-1))

                   heat1 &
                   = heat(i,j,k+1) - S_bar(k+1) &
                     + g_ndim/Pstrat(k+1) * bvsStrat(k+1) &
                       * 0.5*(w_0(k+1) + w_0(k))

                   wAst &
                   = 1.0 &
                     /(  facw*facr &
                       + rhoStratTilde(k)/rhow * bvsstw * dt**2) &
                     * (  facr * (wvert - dt * piGrad) &
                        - dt * g_ndim &
                          * 0.5*(  rhop_0/rho000 + rhop_1/rho001&
                                 + dt &
                                   * (  rhoStrat(k)/Pstrat(k) &
                                        * heat0/rho000 &
                                      + rhoStrat(k+1)/Pstrat(k+1) &
                                        * heat1/rho001 &
                                      + alprlx &
                                        * (  rho_p0/rho000 &
                                           + rho_p1/rho001))))
                   !UAE

                   var(i,j,k,4) = wAst
                end do
             end do
          end do
         else
          stop'ERROR: unknown int_mod'
       end if
      else
       stop'ERROR: unknown mmp_mod'
    end if

  end subroutine momentumPredictor


!-------------------------------------------------------------------------


  subroutine thetaUpdate (var,var0,flux,source,dt,q,m)
    !-----------------------------
    ! adds theta flux to cell theta
    !-----------------------------
    
    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var

    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var0

    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux) 
    ! dir = 1..3 > f-, g- and h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW, fTheta

    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: source
    
    real, intent(in) :: dt
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz), &
         & intent(inout) :: q

    integer, intent(in) :: m
    
    ! local variables
    integer :: i,j,k,l
    real    :: fL,fR        ! flux Left/Right
    real    :: gB,gF        ! flux Backward/Forward
    real    :: hD,hU        ! flux Downward/Upward
    real    :: fluxDiff         ! convective part
    real    :: F            ! F(phi)

    ! advection of background
    real :: u,v,w,w_true
    real :: adv
       

    ! init q
    if (m == 1) q = 0.

    do k = 1,nz
       do j = 1,ny
          do i = 1,nx

             fL = flux(i-1,j,k,1,6) ! theta flux accros left cell edge
             fR = flux(i,j,k,1,6)   ! right
             gB = flux(i,j-1,k,2,6) ! backward
             gF = flux(i,j,k,2,6)   ! forward
             hD = flux(i,j,k-1,3,6) ! downward
             hU = flux(i,j,k,3,6)   ! upward

             ! convective part
             fluxDiff = (fR-fL)/dx + (gF-gB)/dy + (hU-hD)/dz

             ! advective part of background stratification
             u = 0.5*( var(i,j,k,2) + var(i-1,j,k,2) )
             v = 0.5*( var(i,j,k,3) + var(i,j-1,k,3) )
             w = 0.5*( var(i,j,k,4) + var(i,j,k-1,4) )
             w_true = u*vertical(1) + v*vertical(2) + w*vertical(3) 

             adv = w_true*Fr2*theta00*N2

             ! diffusive part
             ! diff = ....

             ! F(phi)
             F = -fluxDiff - adv + source(i,j,k,6)

             select case( timeSchemeType ) 

             case( "lowStorage" ) 

                ! update: q(m-1) -> q(m)
                q(i,j,k) = dt*F + alpha(m) * q(i,j,k)

                ! update potential temperature
                var(i,j,k,6) = var(i,j,k,6) + beta(m) * q(i,j,k)

             case( "classical" )

                var(i,j,k,6) = rk(1,m) * var0(i,j,k,6) &
                     &       + rk(2,m) * var (i,j,k,6) &
                     &       + rk(3,m) * dt*F

             case default
                stop"thetaUpdate: unknown case timeSchemeType"
             end select



          end do
       end do
    end do

    if(verbose .and. master) &
    & print*,"update.f90/thetaUpdate: theta(m=",m,") calculated."

  end subroutine thetaUpdate


!--------------------------------------------------------------------------

  subroutine massUpdate (var,flux,flux_rhopw,dt,q,m,upd_var,upd_mod,int_mod)
    !-----------------------------
    ! adds mass flux to cell mass
    !-----------------------------
    
    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var

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
    character(len=*), intent(in) :: upd_var, upd_mod, int_mod

    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux) 
    ! dir = 1..3 > f-, g- and h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW
    real, dimension(-1:nx,-1:ny,-1:nz), intent(in) :: flux_rhopw
    
    real, intent(in) :: dt
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz), &
         & intent(inout) :: q

    integer, intent(in) :: m
    integer :: i00, j00
    
    ! local variables
    integer :: i,j,k,l
    real    :: fL,fR        ! flux Left/Right
    real    :: gB,gF        ! flux Backward/Forward
    real    :: hD,hU        ! flux Downward/Upward
    real    :: fluxDiff     ! convective part
    real    :: F            ! F(phi)

    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz) :: heat

    real :: buoy0, buoy, rho, rhow, rhowm, rhop, wvrt, facw, facr, &
          & pstw, pstwm, piU, piD, piGrad

    real :: rho_p

    real, dimension(-nbz:nz+nbz) :: w_0 
    real, dimension(-nbz:nz+nbz) :: S_bar 
    real :: heat_flc

    !UAB
    real :: rho_e
    !UAE

    real, dimension(1:nz) :: sum_local, sum_global
    real, dimension(-nbz:nz+nbz) :: avgrhopw

    if( correctDivError ) then
        print*,'ERROR: correction divergence error not allowed'
        stop
    end if

    ! init q
    if (m == 1) q = 0.

    if (upd_var == "rho") then
       if (upd_mod /= "tot" .and. upd_mod /= "lhs") then
          print*,'ERROR: wrong upd_mod for upd_var = rho'
          stop
       end if

       if (int_mod /= "expl") stop'ERROR: wrong int_mod for upd_var = rho'

       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                fL = flux(i-1,j,k,1,1) ! mass flux accros left cell edge
                fR = flux(i,j,k,1,1)   ! right
                gB = flux(i,j-1,k,2,1) ! backward
                gF = flux(i,j,k,2,1)   ! forward
                hD = flux(i,j,k-1,3,1) ! downward
                hU = flux(i,j,k,3,1)   ! upward

                ! convective part
                fluxDiff = (fR-fL)/dx + (gF-gB)/dy + (hU-hD)/dz

                ! F(phi)
                F = -fluxDiff

                !UAB
                ! density relaxation
                if (dens_relax) then
                   if (background /= "HeldSuarez") then
                      stop'ERROR: density relaxation only ready for &
                         &background = HeldSuarez'
                   end if

                   if( fluctuationMode )  then
                      rho = var(i,j,k,1) + rhoStrat(k)
                     else   
                      rho = var(i,j,k,1)
                   end if  

                   rho_e = Pstrat(k)/the_env_pp(i,j,k)

                   F = F - kt_hs(j,k) * (rho-rho_e)
                end if
                !UAE

                ! update: q(m-1) -> q(m)
                q(i,j,k) = dt*F + alpha(m) * q(i,j,k)

                ! update density
                var(i,j,k,1) = var(i,j,k,1) + beta(m) * q(i,j,k)
             end do
          end do
       end do
      else if (upd_var == "rhop") then
       if (upd_mod == "tot") then
          if (int_mod /= "expl") then
             stop'ERROR: wrong int_mod for upd_mod = tot'
          end if

          !UAB
          !! environmental heating
          !call calculate_heating(var,flux,heat)
          !! heating by GW entropy-flux convergence
          !if (raytracer) heat(:,:,:) = heat(:,:,:) + var(:,:,:,8)

          ! heating due to relaxation, entropy diffusion and GWs, its 
          ! horizontal mean and the horizontal-mean vertical wind 
          ! resulting from it

          if (heatingONK14 .or. TurbScheme .or. rayTracer) then
             call heat_w0(var,flux,flux_rhopw,heat,S_bar,w_0)
            else
             heat = 0.
             S_bar = 0.
             w_0 = 0.
          end if
          !UAE

          !FS
          sum_local = 0.
          sum_global = 0.
          avgrhopw = 0.
             do k = 1,nz
             sum_local(k) =  sum(flux(1:nx,1:ny,k,3,6))
             end do
          
          !global sum and average
          call mpi_allreduce(sum_local(1),sum_global(1),&
               nz-1+1,&
               mpi_double_precision,mpi_sum,comm,ierror)
          sum_global = sum_global/(sizeX*sizeY)
          
          avgrhopw(1:nz) = sum_global
          
          do k = 1,nz
             do j = 1,ny
                do i = 1,nx
                   fL = flux(i-1,j,k,1,6) ! mass flux accros left cell edge
                   fR = flux(i,j,k,1,6)   ! right
                   gB = flux(i,j-1,k,2,6) ! backward
                   gF = flux(i,j,k,2,6)   ! forward
                   hD = flux(i,j,k-1,3,6) ! downward
                   hU = flux(i,j,k,3,6)   ! upward
   
                   ! convective part
                   fluxDiff = (fR-fL)/dx + (gF-gB)/dy + (hU-hD)/dz &
                        - (avgrhopw(k)-avgrhopw(k-1))/dz !FS
   
                   rhop = var(i,j,k,6)

                   rho = var(i,j,k,1)
                   if( fluctuationMode ) then
                       rho =  rho + rhoStrat(k)
                   end if

                   !UAB
                   !wvrt = 0.5 * (var(i,j,k,4) + var(i,j,k-1,4))
                   wvrt &
                   = 0.5 * (var(i,j,k,4)-w_0(k) + var(i,j,k-1,4)-w_0(k-1)) 

                   heat_flc= heat(i,j,k) - S_bar(k)
                   !UAE

                   ! F(phi)
                   !UAB
                   !F &
                   != - fluxDiff + rhoStrat(k)/g_ndim * bvsStrat(k)*wvrt &
                   !  + rhoStrat(k)/Pstrat(k) * heat(i,j,k) &
                   !  - alprlx * (rhop - rho + rhoStrat(k)) 
                   F &
                   = - fluxDiff + rhoStrat(k)/g_ndim * bvsStrat(k)*wvrt &
                     + rhoStrat(k)/Pstrat(k) * heat_flc&
                     - alprlx * (rhop - rho + rhoStrat(k))
                   !UAE
   
                   !UAB
                   ! density relaxation
                   if (dens_relax) then
                      if (background /= "HeldSuarez") then
                         stop'ERROR: density relaxation only ready for &
                            &background = HeldSuarez'
                      end if

                      if( fluctuationMode )  then
                         rho = var(i,j,k,1) + rhoStrat(k)
                        else   
                         rho = var(i,j,k,1)
                      end if  

                      rho_e = Pstrat(k)/the_env_pp(i,j,k)
   
                      F = F - kt_hs(j,k) * (rho-rho_e)
                   end if
                   !UAE

                   ! update: q(m-1) -> q(m)
                   q(i,j,k) = dt*F + alpha(m) * q(i,j,k)
   
                   ! update density
                   var(i,j,k,6) = var(i,j,k,6) + beta(m) * q(i,j,k)
                end do
             end do
          end do
         else if (upd_mod == "lhs") then
          if (int_mod /= "expl") then
             stop'ERROR: wrong int_mod for upd_mod = lhs'
          end if

          !FS
          sum_local = 0.
          sum_global = 0.
          avgrhopw = 0.
          do k = 1,nz
             sum_local(k) =  sum(flux(1:nx,1:ny,k,3,6))
          end do
          
          !global sum and average
          call mpi_allreduce(sum_local(1),sum_global(1),&
               nz-1+1,&
               mpi_double_precision,mpi_sum,comm,ierror)
          sum_global = sum_global/(sizeX*sizeY)
          
          avgrhopw(1:nz) = sum_global

          do k = 1,nz
             do j = 1,ny
                do i = 1,nx
                   fL = flux(i-1,j,k,1,6) ! mass flux accros left cell edge
                   fR = flux(i,j,k,1,6)   ! right
                   gB = flux(i,j-1,k,2,6) ! backward
                   gF = flux(i,j,k,2,6)   ! forward
                   hD = flux(i,j,k-1,3,6) ! downward
                   hU = flux(i,j,k,3,6)   ! upward
   
                   ! convective part
                   fluxDiff = (fR-fL)/dx + (gF-gB)/dy + (hU-hD)/dz &
                        - (avgrhopw(k)-avgrhopw(k-1))/dz !FS
   
                   ! F(phi)
                   F = -fluxDiff

                   !UAB
                   ! density relaxation
                   if (dens_relax) then
                      if (background /= "HeldSuarez") then
                         stop'ERROR: density relaxation only ready for &
                            &background = HeldSuarez'
                      end if

                      if( fluctuationMode )  then
                         rho = var(i,j,k,1) + rhoStrat(k)
                        else   
                         rho = var(i,j,k,1)
                      end if  

                      rho_e = Pstrat(k)/the_env_pp(i,j,k)
   
                      F = F - kt_hs(j,k) * (rho-rho_e)
                   end if
                   !UAE 

                   ! update: q(m-1) -> q(m)
                   q(i,j,k) = dt*F + alpha(m) * q(i,j,k)
   
                   ! update density
                   var(i,j,k,6) = var(i,j,k,6) + beta(m) * q(i,j,k)
                end do
             end do
          end do
         else if (upd_mod == "rhs") then
          ! calculate bstar ...

          !UAB
          !! environmental heating
          !call calculate_heating(var,flux,heat)
          !! heating by GW entropy-flux convergence
          !if (raytracer) heat(:,:,:) = heat(:,:,:) + var(:,:,:,8)

          ! heating due to relaxation, entropy diffusion and GWs, its 
          ! horizontal mean and the horizontal-mean vertical wind 
          ! resulting from it

          if (heatingONK14 .or. TurbScheme .or. rayTracer) then
             call heat_w0(var,flux,flux_rhopw,heat,S_bar,w_0)
            else
             heat = 0.
             S_bar = 0.
             w_0 = 0.
          end if
          !UAE

          if (int_mod == "impl") then
             if(topography) then
                i00=is+nbx-1
                j00=js+nby-1
             end if

             do k = 1,nz
                pstw  = 0.5*(Pstrat(k  ) + Pstrat(k+1))
                pstwm = 0.5*(Pstrat(k-1) + Pstrat(k  ))

                do j = 1,ny
                   do i = 1,nx
                      if (TestCase == "baroclinic_LC") then
                         rhop = var(i,j,k,6) - var_env(i,j,k,6)
                        else
                         rhop = var(i,j,k,6)
                      end if

                      rho = var(i,j,k,1)
                      rhow  = 0.5*(var(i,j,k  ,1) + var(i,j,k+1,1))
                      rhowm = 0.5*(var(i,j,k-1,1) + var(i,j,k  ,1))

                      if( fluctuationMode ) then
                            rho =   rho +      rhoStrat(k)
                           rhow =  rhow + rhoStratTilde(k)
                          rhowm = rhowm + rhoStratTilde(k-1)
                      end if

                      !UAB
                      !wvrt = 0.5 * (var(i,j,k,4) + var(i,j,k-1,4))
                      wvrt &
                      = 0.5 &
                        * (var(i,j,k,4)-w_0(k) + var(i,j,k-1,4)-w_0(k-1)) 
   
                      heat_flc= heat(i,j,k) - S_bar(k)
                      !UAE

                      if (TestCase == "baroclinic_LC") then
                         piGrad &
                         = kappaInv*MaInv2 &
                           * 0.5*(  pstw/rhow &
                                    * (  var(i,j,k+1,5) &
                                       - var(i,j,  k,5) &
                                       - var_env(i,j,k+1,5) &
                                       + var_env(i,j,  k,5))/dz &
                                  + pstwm/rhowm &
                                    * (  var(i,j,k  ,5) &
                                       - var(i,j,k-1,5) &
                                       - var_env(i,j,k  ,5) &
                                       + var_env(i,j,k-1,5))/dz)
                        else
                         piGrad &
                         = kappaInv*MaInv2 &
                           * 0.5*(  pstw/rhow &
                                    * (var(i,j,k+1,5)-var(i,j,  k,5))/dz &
                                  + pstwm/rhowm &
                                    * (var(i,j,k  ,5)-var(i,j,k-1,5))/dz)
                      end if
   

                      ! due to relaxation of density fluctuation to the
                      ! correspondiung result from density and 
                      ! reference-state density
                      facr = 1.0 + alprlx*dt

                      ! due to damping of wind in land cells (if there is
                      ! topography)
                      facw = 1.0

                      if(topography) then
                         if(topography_mask(i00+i,j00+j,k)&
                            .or.&
                            topography_mask(i00+i,j00+j,k+1)) then
                            facw = facw + alprlx*dt
                         end if
                      end if

                      if (TestCase == "baroclinic_LC") then
                         rho_p = rho - var_env(i,j,k,6) - rhoStrat(k)
                        else
                         rho_p = rho - rhoStrat(k)
                      end if

                      !UAB
                      !buoy &
                      != 1.0 &
                      !  /( facw*facr &
                      !    + rhoStrat(k)/rho * bvsStrat(k) * dt**2) &
                      !  * (-rhoStrat(k)/rho * bvsStrat(k) * dt &
                      !      *(wvrt - dt * piGrad) &
                      !     - facw * g_ndim/rho &
                      !       * (  rhop &
                      !          + dt &
                      !            * (  rhoStrat(k)/Pstrat(k) &
                      !                 * heat(i,j,k) &
                      !               + alprlx * rho_p)))
                      buoy &
                      = 1.0 &
                        /( facw*facr &
                          + rhoStrat(k)/rho * bvsStrat(k) * dt**2) &
                        * (-rhoStrat(k)/rho * bvsStrat(k) * dt &
                            *(wvrt - dt * piGrad) &
                           - facw * g_ndim/rho &
                             * (  rhop &
                                + dt &
                                  * (  rhoStrat(k)/Pstrat(k) * heat_flc &
                                     + alprlx * rho_p)))
                      !UAE

                      var(i,j,k,6) = -buoy * rho/g_ndim

                      if (TestCase == "baroclinic_LC") then
                         var(i,j,k,6) = var(i,j,k,6) + var_env(i,j,k,6)   
                      end if
                   end do
                end do
             end do
            else if (int_mod == "expl") then
             do k = 1,nz
                do j = 1,ny
                   do i = 1,nx
                      if (TestCase == "baroclinic_LC") then
                         rhop = var(i,j,k,6) - var_env(i,j,k,6)
                        else
                         rhop = var(i,j,k,6)
                      end if


                      rho = var(i,j,k,1)
                      if( fluctuationMode ) then
                          rho =  rho + rhoStrat(k)
                      end if

                      !UAB
                      !wvrt = 0.5 * (var(i,j,k,4) + var(i,j,k-1,4))
                      wvrt &
                      = 0.5 &
                        * (var(i,j,k,4)-w_0(k) + var(i,j,k-1,4)-w_0(k-1)) 
   
                      heat_flc= heat(i,j,k) - S_bar(k)
                      !UAE

                      if (TestCase == "baroclinic_LC") then
                         rho_p = rho - var_env(i,j,k,6) - rhoStrat(k)
                        else
                         rho_p = rho - rhoStrat(k)
                      end if

                      !UAB
                      !buoy &
                      != - g_ndim * rhop/rho &
                      !  - dt &
                      !    * (  rhoStrat(k)/rho * bvsStrat(k)*wvrt &
                      !       + g_ndim/rho &
                      !         * (  rhoStrat(k)/Pstrat(k) * heat(i,j,k) &
                      !            - alprlx * (rhop - rho_p)))
                      buoy &
                      = - g_ndim * rhop/rho &
                        - dt &
                          * (  rhoStrat(k)/rho * bvsStrat(k)*wvrt &
                             + g_ndim/rho &
                               * (  rhoStrat(k)/Pstrat(k) * heat_flc &
                                  - alprlx * (rhop - rho_p)))
                      !UAE

                      var(i,j,k,6) = -buoy * rho/g_ndim

                      if (TestCase == "baroclinic_LC") then
                         var(i,j,k,6) = var(i,j,k,6) + var_env(i,j,k,6)   
                      end if
                   end do
                end do
             end do
            else
             stop'int_mod unknown'
          end if
         else
          stop'upd_mod unknown'
       end if
      else
       stop'upd_var unknown'
    end if

  end subroutine massUpdate

!---------------------------------------------------------------------------


  subroutine iceUpdate (var,var0,flux,source,dt,q,m)
    !-----------------------------
    ! adds ice flux to cell ice
    !-----------------------------
    ! mainly analogous to massUpdate 


    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var0
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux) 
    ! dir = 1..3 > f-, g- and h-flux in x,y,z-direction
    ! iFlux = 8..11 > Rho_nAer, Rho_nIce, Rho_qIce, Rho_qv

    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: source
    
    real, intent(in) :: dt
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,4), &
         & intent(inout) :: q

    integer, intent(in) :: m

    ! local integer
    integer :: i,j,k,iVar
    
    ! local variables
    real, dimension(4)    :: fL,fR        ! flux Left/Right
    real, dimension(4)    :: gB,gF        ! flux Backward/Forward
    real, dimension(4)    :: hD,hU        ! flux Downward/Upward
    real, dimension(4)    :: fluxDiff         ! convective part
    real, dimension(4)    :: F            ! F(phi)
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz) :: rho
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,4) :: rho_source_term
    real :: T,p,SIce
    
    var0 = var

    ! init q
    if (m == 1) q = 0.

    if ( fluctuationMode ) then
      do k = -1,nz+1
        rho(:,:,k) = var(:,:,k,1)+rhoStrat(k) 
      end do
    else
      rho = var(:,:,:,1)
    end if 

    if (correctDivError) then
      do k = 0,3
        rho_source_term(:,:,:,4-k) = var(:,:,:,nVar-k) * source(:,:,:,1)
      end do
    else 
      rho_source_term(:,:,:,:) = 0.0
    end if

    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
          
           if (topography_mask(i+is+nbx-1,j+js+nby-1,k)==.false.) then
          
             fL = flux(i-1,j,k,1,nVar-3:nVar) ! mass flux across left cell edge
             fR = flux(i,j,k,1,nVar-3:nVar)   ! right
             gB = flux(i,j-1,k,2,nVar-3:nVar) ! backward
             gF = flux(i,j,k,2,nVar-3:nVar)   ! forward
             hD = flux(i,j,k-1,3,nVar-3:nVar) ! downward
             hU = flux(i,j,k,3,nVar-3:nVar)   ! upward

             ! convective part
             fluxDiff = (fR-fL)/dx + (gF-gB)/dy + (hU-hD)/dz

             ! diffusive part
             ! diff = ....

             ! F(phi)
             F = -fluxDiff

             F(:) = F(:) + rho_source_term(i,j,k,:) + rho(i,j,k)*source(i,j,k,nVar-3:nVar)
             
             select case( timeSchemeType ) 
                
             case( "lowStorage" ) 

                ! update: q(m-1) -> q(m)
                q(i,j,k,:) = dt*F(:) + alpha(m) * q(i,j,k,:)

                ! update variables
                var(i,j,k,nVar-3:nVar) = var(i,j,k,nVar-3:nVar) + beta(m) * q(i,j,k,1:4) / rho(i,j,k)

             case( "classical" )

                var(i,j,k,nVar-3:nVar) = rk(1,m) * var0(i,j,k,nVar-3:nVar) &
                     &       + rk(2,m) * var (i,j,k,nVar-3:nVar) &
                     &       + rk(3,m) * dt*F(1:4) / rho(i,j,k)

             case default
                stop "iceUpdate: unknown case timeSchemeType"
             end select
             
             do iVar = nVar-3,nVar
               ! avoid negative values for all ice variables
               if ((var(i,j,k,iVar).lt. 0.0)) then
                 var(i,j,k,iVar) = 0.0
               end if
             end do
            
           end if
          end do
       end do
    end do

    if(verbose .and. master) print*,"update.f90/iceUpdate: ice(m=",m,") calculated." 

  end subroutine iceUpdate

!-------------------------------------------------------------------------


  subroutine timestep (var,dt,errFlag)
    !---------------------------------------------
    ! compute time step from stability criteria:
    ! 1) CFL criterion for advection
    ! 2) von Neumann cirterion for dissipation
    ! 3) set maximum time step
    ! 4) buouyancy acceleration ->  1/2 b_max*dt^2 < dz
    ! 5) gravity-wave group condition: c_g_x * dt < dx, ...dy,dz
    !---------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),intent(in) &
            & :: var
    real, intent(out) :: dt
    logical, intent(out) :: errFlag

    ! locals
    real :: uMax, vMax, wMax
    real :: dtConv, dtVisc, dtCond, dtWKB
    real :: dtConv_loc, dtWKB_loc
    real :: dtMax
    real :: dtWave

    ! Buoyancy time step restriction
    real              :: dtBuoy, dtBuoy_loc
    real,dimension(3) :: bMax, bMaxNew, duMax

    ! achatzb test deletion:
    ! ! Gravity wave time stop restriction
    ! real :: lambdaMax    ! max GW length to be time resolved
    ! real :: dtWave, lambdaX, lambdaZ, kk, mm, kMin, cX, cZ
    ! achatze
    
    ! local integer
    integer :: i,j,k

    ! sponge layer
    real :: ReMin         ! min Reynolds number in domain
    
    ! verbose
    logical, parameter :: giveInfo = .true.
    
    errFlag = .false.

    !-------------------------------------------
    !              Fixed time step
    !-------------------------------------------

    if( tStepChoice == "fix" ) then  

       dt = dtMax_dim / tRef
       errFlag = .false.

       if( master ) then
       write(*,fmt="(a25,es15.1,a8)") "dt = dtFix = ", dt*tRef, "seconds"
       end if

    else

    !-------------------------------------------
    !           Variable time step
    !-------------------------------------------

       select case( model ) 

       case( "Boussinesq", "pseudo_incompressible")

          !----------------------------
          !   Full model time step
          !----------------------------

          !----------------------
          !     CFL condition
          !----------------------
          uMax = maxval( abs(var(1:nx,1:ny,1:nz,2) )) + small
          vMax = maxval( abs(var(1:nx,1:ny,1:nz,3) )) + small
          wMax = maxval( abs(var(1:nx,1:ny,1:nz,4) )) + small

          dtConv_loc = cfl * min(dx/uMax, dy/vMax, dz/wMax)

          ! find global minimum

          call mpi_reduce(dtConv_loc, dtConv, 1, mpi_double_precision,&
                        & mpi_min, root, comm, ierror)

          call mpi_bcast(dtConv, 1, mpi_double_precision, root, comm, &
                       & ierror)


          !---------------------------
          !  Acceleration condition
          !---------------------------
          bMax = 0.0
          do k = 1,nz
             do j = 1,ny
                do i = 1,nx

                   select case( model ) 
                      
                   case( "pseudo_incompressible" )
                      if( fluctuationMode ) then
                         bMaxNew = abs(var(i,j,k,1)) / &
                              & (rhoStrat(k) + var(i,j,k,1)) * vertical
                      else
                         bMaxNew = abs(rhoStrat(k)-var(i,j,k,1)) &
                              & / var(i,j,k,1) * vertical
                      end if
                      
                   case( "Boussinesq" )
                      bMaxNew = var(i,j,k,6)/theta00 * vertical 

                   case default
                      stop"timeStep: unknown case model."
                   end select
                   
                   if( bMaxNew(1) > bMax(1) ) bMax(1) = bMaxNew(1)
                   if( bMaxNew(2) > bMax(2) ) bMax(2) = bMaxNew(2)
                   if( bMaxNew(3) > bMax(3) ) bMax(3) = bMaxNew(3)
                end do
             end do
          end do
          bMax = FrInv2 * bMax
          
          ! check whether acceleration condition is needed
          duMax = bMax * dtConv
          if( (  duMax(1) > 1.e-2*uMax .or.&
               & duMax(2) > 1.e-2*vMax .or.&
               & duMax(3) > 1.e-2*wMax).and.&
              ((bMax(1) /= 0.).and.(bMax(2) /= 0.).and.(bMax(3) /= 0.)) ) &
              then
             
             dtBuoy_loc = max(&
                & -uMax/bMax(1)+sqrt((uMax/bMax(1))**2+2.*cfl*dx/bMax(1)),&
                & -vMax/bMax(2)+sqrt((vMax/bMax(2))**2+2.*cfl*dy/bMax(2)),&
                & -wMax/bMax(3)+sqrt((wMax/bMax(3))**2+2.*cfl*dz/bMax(3)))
             
             !xxxx debug             
             if( dtBuoy_loc*tRef < 1.e-2 ) then

                print*,"dtBuoy_loc*tRef  = ", dtBuoy_loc*tRef
                print*,"bMax(3) = ", bMax(3)*FrInv2
                
             end if
             !xxxx end debug
             
          else
             
             dtBuoy_loc = 1.0e20/tRef   ! set to high value if not needed

          end if
          
          ! find global minimum

          call mpi_reduce(dtBuoy_loc, dtBuoy, 1, mpi_double_precision,&
                        & mpi_min, root, comm, ierror)

          call mpi_bcast(dtBuoy, 1, mpi_double_precision, root, comm, &
                       & ierror)

          !---------------------------
          !   von Neumann condition
          !----------------------------

          dtVisc = 0.5 * min(dx**2, dy**2, dz**2) * Re
          dtCond = 0.5 * min(dx**2, dy**2, dz**2) / mu_conduct


          !----------------------------
          !    Maximal time step
          !----------------------------
          dtMax = dtMax_dim / tRef



          !------------------------------------
          !    Gravity wave time period
          !------------------------------------
          
          !UAB
          !dtWave = pi/(NN+small)
          dtWave = 1.7/(NN+small)
          !UAE
          
          !------------------------------------
          !     WKB "CFL" criterion
          !------------------------------------

          if (raytracer) then
             dtWKB_loc = dz/(cgz_max + small)

             if(sizeX > 1) dtWKB_loc = min(dtWKB_loc, dx/(cgx_max + small))
             if(sizeY > 1) dtWKB_loc = min(dtWKB_loc, dy/(cgy_max + small))

             dtWKB_loc = cfl_wave * dtWKB_loc

             ! find global minimum

             call mpi_reduce(dtWKB_loc, dtWKB, 1, mpi_double_precision,&
                           & mpi_min, root, comm, ierror)

             call mpi_bcast(dtWKB, 1, mpi_double_precision, root, comm, &
                          & ierror)
          end if

          !-------------------------------
          !        Make your choice
          !-------------------------------

          if ( dtWave_on .and. timeScheme /= 'semiimplicit'  ) then
             dt = min(dtVisc,dtCond,dtConv,dtMax,dtBuoy,dtWave)
            else
             dt = min(dtVisc,dtCond,dtConv,dtMax,dtBuoy)
             !if (timeScheme == 'semiimplicit') then
             !   dt = min(dtVisc,dtCond,dtConv,dtMax,dtBuoy,10.*dtWave)
             !  else
             !   dt = min(dtVisc,dtCond,dtConv,dtMax,dtBuoy)
             !end if
          end if      

          if (raytracer) dt = min(dt, dtWKB)

          !-----------------------------------------
          !     Inform on time step restrictions
          !-----------------------------------------

       if( master ) then

          write(*,fmt="(a25,es15.1,a8)") "dtVisc =", dtVisc*tRef, "seconds"
          write(*,fmt="(a25,es15.1,a8)") "dtCond =", dtCond*tRef, "seconds"
          write(*,fmt="(a25,es15.1,a8)") "dtConv =", dtConv*tRef, "seconds"
          write(*,fmt="(a25,es15.1,a8)") "dtMax =", dtMax*tRef, "seconds"
          write(*,fmt="(a25,es15.1,a8)") "dtBuoy =", dtBuoy*tRef, "seconds"
          write(*,fmt="(a25,es15.1,a8)") "dtWave =", dtWave*tRef, "seconds"
          if (raytracer) then
             write(*,fmt="(a25,es15.1,a8)") "dtWKB =", dtWKB*tRef,"seconds"
          end if
          print*,""

          if(dt == dtMax) then
             write(*,fmt="(a25,es15.1,a8)") "--> dt = dtMax = ", dt*tRef, &
                     & "seconds"
          else if(dt == dtConv) then
             write(*,fmt="(a25,es15.1,a8)") "--> dt = dtConv = ", dt*tRef,&
                    &  "seconds"
          else if(dt == dtVisc) then
             write(*,fmt="(a25,es15.1,a8)") "--> dt = dtVisc = ", dt*tRef,&
                    &  "seconds"
          else if(dt == dtCond) then
             write(*,fmt="(a25,es15.1,a8)") "--> dt = dtCond = ", dt*tRef,&
                    &  "seconds"
          else if(dt == dtBuoy) then
             write(*,fmt="(a25,es15.1,a8)") "--> dt = dtBuoy = ", dt*tRef,&
                    &  "seconds"
          else if(dt == dtWave) then
             write(*,fmt="(a25,es15.1,a8)") "--> dt = dtWave = ", dt*tRef,&
                    &  "seconds"
          else if(dt == dtWKB) then
             write(*,fmt="(a25,es15.1,a8)") "--> dt = dtWKB =", dt*tRef, &
                                          & "seconds"
          else
             write(*,fmt="(a25,es15.1,a8)") "--> dt = ????? = ", dt*tRef, &
                     & "seconds"
          end if
          print*,""

       end if

       case default
          stop"timestep: unknown case model."
       end select   ! WKB / full model

    end if

    ! error handling for too small time steps
    if( dt*tRef < dtMin_dim ) errFlag = .true. 

  end subroutine timestep


!-------------------------------------------------------------------------


  subroutine init_update
    !--------------------------------------
    ! allocate variables for update module
    !--------------------------------------
    
    ! local variables
    integer :: allocstat

    ! number of equations
    nxyz = nx*ny*nz
    
    ! number of non-zeros in the equations system
    nnz = 7*nxyz - 2*nx*ny
    ! xxx: verify this value when changing boundary conditions

          
  end subroutine init_update


!---------------------------------------------------------------------

  subroutine CoefDySma_update(var)
    !--------------------------------------
    ! calculate the Coefficient for Dynamic Smagorinsky Scheme
    !--------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var

    ! more variables
    real :: delta_hs, delta_vs
    real :: uL, uR, uB, uF, uD, uU
    real :: vL, vR, vB, vF, vD, vU
    real :: wL, wR, wB, wF, wD, wU
    real :: du_dx, du_dy, du_dz
    real :: dv_dx, dv_dy, dv_dz
    real :: dw_dx, dw_dy, dw_dz

    ! allocatable fields
    real, dimension(:,:,:,:,:), allocatable :: Sij, Lij, &
                                               Mij
    real, dimension(:,:,:), allocatable :: S_norm

    real, dimension(:,:,:,:,:), allocatable :: uiuj_smth, &
                                               S_Sij_smth, &
                                               Sij_smth
    real, dimension(:,:,:,:), allocatable :: ui_smth
    real, dimension(:,:,:), allocatable :: Sn_smth
    real, dimension(:,:,:), allocatable :: LijMij_smth, &
                                           MijMij_smth
    real, dimension(:,:,:), allocatable :: CS2_DySma

    integer :: allocstat
    integer :: i,j,k
    integer :: iw,jw

    integer :: smth_npts1_DySma, smth_npts2_DySma
    parameter (smth_npts1_DySma = 1)   ! revised by JW (20160824)
    parameter (smth_npts2_DySma = 2)  ! revised by JW (20160824)

    ! Allocate local fields
    allocate(Sij(1:nx,1:ny,1:nz,1:3,1:3), stat=allocstat)
    if(allocstat/=0) stop"CoefDySma_update:alloc failed"

    allocate(Lij(1:nx,1:ny,1:nz,1:3,1:3), stat=allocstat)
    if(allocstat/=0) stop"CoefDySma_update:alloc failed"

    allocate(Mij(1:nx,1:ny,1:nz,1:3,1:3), stat=allocstat)
    if(allocstat/=0) stop"CoefDySma_update:alloc failed"

    allocate(S_norm(1:nx,1:ny,1:nz), stat=allocstat)
    if(allocstat/=0) stop"CoefDySma_update:alloc failed"

    allocate(uiuj_smth(1:nx,1:ny,1:nz,1:3,1:3), stat=allocstat)
    if(allocstat/=0) stop"CoefDySma_update:alloc failed"

    allocate(S_Sij_smth(1:nx,1:ny,1:nz,1:3,1:3), stat=allocstat)
    if(allocstat/=0) stop"CoefDySma_update:alloc failed"

    allocate(Sij_smth(1:nx,1:ny,1:nz,1:3,1:3), stat=allocstat)
    if(allocstat/=0) stop"CoefDySma_update:alloc failed"

    allocate(ui_smth(1:nx,1:ny,1:nz,1:3), stat=allocstat)
    if(allocstat/=0) stop"CoefDySma_update:alloc failed"

    allocate(Sn_smth(1:nx,1:ny,1:nz), stat=allocstat)
    if(allocstat/=0) stop"CoefDySma_update:alloc failed"

    allocate(LijMij_smth(1:nx,1:ny,1:nz), stat=allocstat)
    if(allocstat/=0) stop"CoefDySma_update:alloc failed"

    allocate(MijMij_smth(1:nx,1:ny,1:nz), stat=allocstat)
    if(allocstat/=0) stop"CoefDySma_update:alloc failed"

    allocate(CS2_DySma(1:nx,1:ny,1:nz), stat=allocstat)
    if(allocstat/=0) stop"CoefDySma_update:alloc failed"

    ! calculate delta

    if (TurbScheme) then
       if (ny == 1 .and. nx == 1) then
          stop'ERROR: turbulence assumes either nx > 1 or ny > 1'
         else
          if (nx == 1) then
             delta_hs = dy**2 ! 2D problems in y and z
            else if (ny == 1) then
             delta_hs = dx**2 ! 2D problems in x and z
            else
             delta_hs = dx*dy ! 3D problems

             if (dx/dy > 10.) then
                print*,'WARNING: dx/dy > 10!'
                print*,'The turbulence scheme is not ready for such &
                        & horizontal grid anisotropies!'
               elseif (dy/dx > 10.) then
                print*,'WARNING: dy/dx > 10!'
                print*,'The turbulence scheme is not ready for such &
                        & horizontal grid anisotropies!'
             end if
          end if

          delta_vs = dz**2
       end if
    end if

    ! calculate S_ij

    !---------------------------------
    !         Loop over field
    !---------------------------------

    do k = 1,nz
       do j = 1,ny
          do i = 1,nx

          uL = var(i-1,j,k,2)
          uR = var(i,j,k,2)

          ! UA not sure whether this averaging is necessary. 
          ! Compare, e.g. the viscous fluxes. There it is not done.

          ! replied by JW (20160824): For now, there is no revision on 
          ! this part, since it is a relatively minor issue at the moment. 

          uB = 0.5 * ( var(i-1,j-1,k,2) + var(i,j-1,k,2) )
          uF = 0.5 * ( var(i-1,j+1,k,2) + var(i,j+1,k,2) )
          uD = 0.5 * ( var(i-1,j,k-1,2) + var(i,j,k-1,2) )
          uU = 0.5 * ( var(i-1,j,k+1,2) + var(i,j,k+1,2) )

          vL = 0.5 * ( var(i-1,j,k,3) + var(i-1,j-1,k,3) )
          vR = 0.5 * ( var(i+1,j,k,3) + var(i+1,j-1,k,3) )
          vB = var(i,j-1,k,3)
          vF = var(i,j,k,3)
          vD = 0.5 * ( var(i,j,k-1,3) + var(i,j-1,k-1,3) )
          vU = 0.5 * ( var(i,j,k+1,3) + var(i,j-1,k+1,3) )

          wL = 0.5 * ( var(i-1,j,k-1,4) + var(i-1,j,k,4) )
          wR = 0.5 * ( var(i+1,j,k-1,4) + var(i+1,j,k,4) )
          wB = 0.5 * ( var(i,j-1,k-1,4) + var(i,j-1,k,4) )
          wF = 0.5 * ( var(i,j+1,k-1,4) + var(i,j+1,k,4) )
          wD = var(i,j,k-1,4)
          wU = var(i,j,k,4)

          du_dx = ( uR - uL ) / dx

          !UA in case without averaging, no factor 2!

          ! replied by JW (20160824): For now, there is no revision on 
          ! this part, since it is a relatively minor issue at the moment. 

          du_dy = ( uF - uB ) / ( 2.0 * dy )
          du_dz = ( uU - uD ) / ( 2.0 * dz )

          dv_dx = ( vR - vL ) / ( 2.0 * dx )
          dv_dy = ( vF - vB ) / dy
          dv_dz = ( vU - vD ) / ( 2.0 * dz )

          dw_dx = ( wR - wL ) / ( 2.0 * dx )
          dw_dy = ( wF - wB ) / ( 2.0 * dy )
          dw_dz = ( wU - wD ) / dz

          Sij(i,j,k,1,1) = 0.5 * ( du_dx + du_dx )
          Sij(i,j,k,1,2) = 0.5 * ( du_dy + dv_dx )
          Sij(i,j,k,1,3) = 0.5 * ( du_dz + dw_dx )
          Sij(i,j,k,2,1) = 0.5 * ( dv_dx + du_dy )
          Sij(i,j,k,2,2) = 0.5 * ( dv_dy + dv_dy )
          Sij(i,j,k,2,3) = 0.5 * ( dv_dz + dw_dy )
          Sij(i,j,k,3,1) = 0.5 * ( dw_dx + du_dz )
          Sij(i,j,k,3,2) = 0.5 * ( dw_dy + dv_dz )
          Sij(i,j,k,3,3) = 0.5 * ( dw_dz + dw_dz )

          S_norm(i,j,k) = 0.0

          do jw = 1,3
             do iw = 1,3
                S_norm(i,j,k) &
                = S_norm(i,j,k) + Sij(i,j,k,iw,jw) * Sij(i,j,k,iw,jw)
             end do
          end do

          S_norm(i,j,k) = S_norm(i,j,k) * 2.0
          S_norm(i,j,k) = sqrt( S_norm(i,j,k) )
       end do
    end do
 end do

 !---------------------------------
 !         Loop over field
 !         Prepare the data (before smoothing at the 1st step)
 !---------------------------------

  do k = 1,nz
     do j = 1,ny
        do i = 1,nx
           ui_smth(i,j,k,1) = 0.5 * ( var(i,j,k,2) + var(i-1,j,k,2) )
           ui_smth(i,j,k,2) = 0.5 * ( var(i,j,k,3) + var(i,j-1,k,3) )
           ui_smth(i,j,k,3) = 0.5 * ( var(i,j,k,4) + var(i,j,k-1,4) )

           Sn_smth(i,j,k)  =  S_norm(i,j,k)

           do jw = 1,3
              do iw = 1,3
                 uiuj_smth(i,j,k,iw,jw) &
                 = ui_smth(i,j,k,iw) * ui_smth(i,j,k,jw)

                 Sij_smth(i,j,k,iw,jw) = Sij(i,j,k,iw,jw)


                 S_Sij_smth(i,j,k,iw,jw) &
                 = S_norm(i,j,k) * Sij(i,j,k,iw,jw)
              end do
           end do
        end do
     end do
  end do

  !---------------------------------
  !         Smoothing at the 1st step
  !---------------------------------

  if(ny.eq.1)then
     do iw = 1,3
        call Var3DSmthDySma( ui_smth(1:nx,1:ny,1:nz,iw), &
                             & smth_npts1_DySma, "XZ_local_smth" )
     end do

     call Var3DSmthDySma( Sn_smth(1:nx,1:ny,1:nz), &
                          & smth_npts1_DySma, "XZ_local_smth" )

     do jw = 1,3
        do iw = 1,3
           call Var3DSmthDySma( uiuj_smth(1:nx,1:ny,1:nz,iw,&
                                                & jw), &
                                & smth_npts1_DySma, "XZ_local_smth" )
           call Var3DSmthDySma( Sij_smth(1:nx,1:ny,1:nz,iw,&
                                               & jw),  &
                                & smth_npts1_DySma, "XZ_local_smth" )
           call Var3DSmthDySma( S_Sij_smth(1:nx,1:ny,1:nz,iw,&
                                                   & jw), &
                                & smth_npts1_DySma, "XZ_local_smth" )
        end do
     end do
  elseif (nx.eq.1)then
     do iw = 1,3
        call Var3DSmthDySma( ui_smth(1:nx,1:ny,1:nz,iw), &
                             & smth_npts1_DySma, "YZ_local_smth" )
     end do

     call Var3DSmthDySma( Sn_smth(1:nx,1:ny,1:nz), &
                          & smth_npts1_DySma, "YZ_local_smth" )

     do jw = 1,3
        do iw = 1,3
           call Var3DSmthDySma( uiuj_smth(1:nx,1:ny,1:nz,iw,&
                                                & jw), &
                                & smth_npts1_DySma, "YZ_local_smth" )
           call Var3DSmthDySma( Sij_smth(1:nx,1:ny,1:nz,iw,&
                                               & jw),  &
                                & smth_npts1_DySma, "YZ_local_smth" )
           call Var3DSmthDySma( S_Sij_smth(1:nx,1:ny,1:nz,iw,&
                                                   & jw), &
                                & smth_npts1_DySma, "YZ_local_smth" )
        end do
     end do
  else
     do iw = 1,3
        call Var3DSmthDySma( ui_smth(1:nx,1:ny,1:nz,iw), &
                             & smth_npts1_DySma, "XYZ_local_smth" )
     end do

     call Var3DSmthDySma( Sn_smth(1:nx,1:ny,1:nz), &
                          & smth_npts1_DySma, "XYZ_local_smth" )

     do jw = 1,3
        do iw = 1,3
           call Var3DSmthDySma( uiuj_smth(1:nx,1:ny,1:nz,iw,&
                                                & jw), &
                                & smth_npts1_DySma, "XYZ_local_smth" )
           call Var3DSmthDySma( Sij_smth(1:nx,1:ny,1:nz,iw,&
                                               & jw),  &
                                & smth_npts1_DySma, "XYZ_local_smth" )
           call Var3DSmthDySma( S_Sij_smth(1:nx,1:ny,1:nz,iw,&
                                                   & jw), &
                                & smth_npts1_DySma, "XYZ_local_smth" )
        end do
     end do
  end if

  !---------------------------------
  !         Loop over field
  !         (Smoothing at the 2nd step)
  !---------------------------------

  do k = 1,nz
     do j = 1,ny
        do i = 1,nx
           LijMij_smth(i,j,k) = 0.0
           MijMij_smth(i,j,k) = 0.0

           do jw = 1,3
              do iw = 1,3
                 Lij(i,j,k,iw,jw) &
                 = uiuj_smth(i,j,k,iw,jw) &
                   - ui_smth(i,j,k,iw) * ui_smth(i,j,k,jw)

                 Mij(i,j,k,iw,jw) &
                 = S_Sij_smth(i,j,k,iw,jw) &
                   - (2.0*smth_npts1_DySma+1.0)**2 &
                     * Sn_smth(i,j,k) * Sij_smth(i,j,k,iw,jw)

                 ! allow for grid anisotropy

                 if (iw == 3 .or. jw == 3) then
                    Mij(i,j,k,iw,jw) = Mij(i,j,k,iw,jw) * delta_vs
                   else
                    Mij(i,j,k,iw,jw) = Mij(i,j,k,iw,jw) * delta_hs
                 end if

                 LijMij_smth(i,j,k) &
                 = LijMij_smth(i,j,k) + Lij(i,j,k,iw,jw) * Mij(i,j,k,iw,jw)

                 MijMij_smth(i,j,k) &
                 = MijMij_smth(i,j,k) &
                   + Mij(i,j,k,iw,jw) * Mij(i,j,k,iw,jw)
              end do
           end do
        end do
     end do
  end do

  if(ny.eq.1)then
     call Var3DSmthDySma( LijMij_smth(1:nx,1:ny,1:nz), &
                          & smth_npts2_DySma, "XZ_local_smth" )
     call Var3DSmthDySma( MijMij_smth(1:nx,1:ny,1:nz), &
                          & smth_npts2_DySma, "XZ_local_smth" )
  elseif(nx.eq.1)then
     call Var3DSmthDySma( LijMij_smth(1:nx,1:ny,1:nz), &
                          & smth_npts2_DySma, "YZ_local_smth" )
     call Var3DSmthDySma( MijMij_smth(1:nx,1:ny,1:nz), &
                          & smth_npts2_DySma, "YZ_local_smth" )
  else
     call Var3DSmthDySma( LijMij_smth(1:nx,1:ny,1:nz), &
                          & smth_npts2_DySma, "XYZ_local_smth" )
     call Var3DSmthDySma( MijMij_smth(1:nx,1:ny,1:nz), &
                          & smth_npts2_DySma, "XYZ_local_smth" )
  end if


  !---------------------------------
  !         Get the final results
  !---------------------------------

  do k = 1,nz
     do j = 1,ny
        do i = 1,nx

           if(MijMij_smth(i,j,k) /= 0.) then
              CS2_DySma(i,j,k) &
              = 0.5 * LijMij_smth(i,j,k) / MijMij_smth(i,j,k)
             else
              CS2_DySma(i,j,k)=0.
           end if

           if( CS2_DySma(i,j,k) < 0.0 )then
              CS2_DySma(i,j,k) = 0.0
           end if

           var(i,j,k,7) &
           = CS2_DySma(i,j,k) * S_norm(i,j,k)
        end do
     end do
  end do

  ! *** set values for the ghost celss ***
  call setHaloAndBoundary( var(:,:,:,7), nbx, nby, nbz )

  ! deallocate local fields
  deallocate(Sij, stat=allocstat); if(allocstat/=0) &
       & stop"update.f90:dealloc failed"
  deallocate(Lij, stat=allocstat); if(allocstat/=0) &
       & stop"update.f90:dealloc failed"
  deallocate(Mij, stat=allocstat); if(allocstat/=0) &
       & stop"update.f90:dealloc failed"
  deallocate(S_norm, stat=allocstat); if(allocstat/=0) &
       & stop"update.f90:dealloc failed"
  deallocate(uiuj_smth, stat=allocstat); if(allocstat/=0) &
       & stop"update.f90:dealloc failed"
  deallocate(S_Sij_smth, stat=allocstat); if(allocstat/=0) &
       & stop"update.f90:dealloc failed"
  deallocate(Sij_smth, stat=allocstat); if(allocstat/=0) &
       & stop"update.f90:dealloc failed"
  deallocate(ui_smth, stat=allocstat); if(allocstat/=0) &
       & stop"update.f90:dealloc failed"
  deallocate(Sn_smth, stat=allocstat); if(allocstat/=0) &
       & stop"update.f90:dealloc failed"
  deallocate(LijMij_smth, stat=allocstat); if(allocstat/=0) &
       & stop"update.f90:dealloc failed"
  deallocate(MijMij_smth, stat=allocstat); if(allocstat/=0) &
       & stop"update.f90:dealloc failed"
  deallocate(CS2_DySma, stat=allocstat); if(allocstat/=0) &
       & stop"update.f90:dealloc failed"

  return


  end subroutine CoefDySma_update


  subroutine Var3DSmthDySma(var3D_DySma,nsmth_DySma,homog_dir_DySma)
    !--------------------------------------
    ! calculate the Coefficient for Dynamic Smagorinsky Scheme
    !--------------------------------------

    ! in/out variables
    real, dimension(1:nx,1:ny,1:nz), &
         & intent(inout) :: var3D_DySma
    integer, intent(in) :: nsmth_DySma

    character(len=*), intent(in) :: homog_dir_DySma

    ! allocatable fields
    real, dimension(:,:,:), allocatable :: var3D_DySma_Extend

    integer :: allocstat
    integer :: i,j,k
    integer :: kmin,kmax,nsmthv
    integer :: nsmthall,ismth,jsmth,ksmth
    integer :: i0,j0

    allocate(var3D_DySma_Extend( (0-nsmth_DySma) : (nx+nsmth_DySma), &
                               & (0-nsmth_DySma) : (ny+nsmth_DySma), &
                               & (0-nsmth_DySma) : (nz+nsmth_DySma) ), &
             & stat=allocstat)
    if(allocstat/=0) stop"Var3DSmthDySma:alloc failed"

    ! set the values for var3D_DySma_Extend

    var3D_DySma_Extend(1:nx,1:ny,1:nz) = var3D_DySma(1:nx,1:ny,1:nz)

    call setHaloAndBoundary( var3D_DySma_Extend(:,:,:), nsmth_DySma, &
                             & nsmth_DySma, nsmth_DySma )

    ! start to do the smoothing

    i0=is+nbx-1
    j0=js+nby-1

    select case( homog_dir_DySma )

    case( "XYZ_local_smth" )

       if(xBoundary /= "periodic") &
       stop "DYNAMIC SMAGORINSKY NOT READY FOR NON-PERIODIC BOUNDARY &
             &CONDITIONS IN X"

       if(yBoundary /= "periodic") &
       stop "DYNAMIC SMAGORINSKY NOT READY FOR NON-PERIODIC BOUNDARY &
             &CONDITIONS IN Y"

       !---------------------------------
       !         Loop over field
       !---------------------------------

       if(nz /= sizeZ) stop" DYNAMIC SMAGORINSKY NOT READY FOR MPI IN Z"

       do k = 1,nz
          ! correct handling of solid and periodic boundaries in z

          if(zBoundary == "solid_wall") then
             kmin=max( 1,k-nsmth_DySma)
             kmax=min(nz,k+nsmth_DySma)
            else if (zBoundary == "periodic") then
             kmin=k-nsmth_DySma
             kmax=k+nsmth_DySma
            else
             stop"vertical smoothing: unknown case zBoundary."
          end if

          nsmthv=kmax-kmin+1


          do j = 1,ny
             do i = 1,nx
                ! averaging dyn. Smag. coeff. only over atmosphere cells

                if(topography) then
                   nsmthall=0
                   var3D_DySma(i,j,k)=0.0

                   do ksmth=kmin,kmax
                      do jsmth=j-nsmth_DySma,j+nsmth_DySma
                         do ismth=i-nsmth_DySma,i+nsmth_DySma
                            if(.not.&
                               topography_mask(i0+ismth,j0+jsmth,ksmth)) &
                               then
                               var3D_DySma(i,j,k)&
                               =var3D_DySma(i,j,k)&
                                +var3D_DySma_Extend(ismth,jsmth,ksmth)

                               nsmthall=nsmthall+1
                            end if
                         end do
                       end do
                    end do

                   if(nsmthall > 0) then
                      var3D_DySma(i,j,k)=var3D_DySma(i,j,k)/nsmthall
                   end if
                  else
                   var3D_DySma(i,j,k)&
                   =&
                   sum(&
                   var3D_DySma_Extend( &
                   (i-nsmth_DySma):(i+nsmth_DySma), &
                   (j-nsmth_DySma):(j+nsmth_DySma), &
                   kmin:kmax &
                   )&
                   )&
                   /( ((2*nsmth_DySma + 1)**2) * nsmthv)
                end if
             end do
          end do
       end do

    case( "XZ_local_smth" )

       if(xBoundary /= "periodic") &
       stop "DYNAMIC SMAGORINSKY NOT READY FOR NON-PERIODIC BOUNDARY &
             &CONDITIONS IN X"

       !---------------------------------
       !         Loop over field
       !---------------------------------

       if(nz /= sizeZ) stop" DYNAMIC SMAGORINSKY NOT READY FOR MPI IN Z"

       do k = 1,nz
          ! correct handling of solid and periodic boundaries in z

          if(zBoundary == "solid_wall") then
             kmin=max( 1,k-nsmth_DySma)
             kmax=min(nz,k+nsmth_DySma)
            else if (zBoundary == "periodic") then
             kmin=k-nsmth_DySma
             kmax=k+nsmth_DySma
            else
             stop"vertical smoothing: unknown case zBoundary."
          end if

          nsmthv=kmax-kmin+1
             
          do j = 1,ny
             do i = 1,nx
                ! averaging dyn. Smag. coeff. only over atmosphere cells

                if(topography) then
                  nsmthall=0
                  var3D_DySma(i,j,k)=0.0

                  do ksmth=kmin,kmax
                     do ismth=i-nsmth_DySma,i+nsmth_DySma
                        if(.not.&
                           topography_mask(i0+ismth,j0+j,ksmth)) then
                           var3D_DySma(i,j,k)&
                           =var3D_DySma(i,j,k)&
                            +var3D_DySma_Extend(ismth,j,ksmth)

                           nsmthall=nsmthall+1
                        end if
                     end do
                  end do

                  if(nsmthall > 0) then
                     var3D_DySma(i,j,k)=var3D_DySma(i,j,k)/nsmthall
                  end if
                 else
                   var3D_DySma(i,j,k)&
                   =&
                   sum(&
                   var3D_DySma_Extend( &
                   (i-nsmth_DySma):(i+nsmth_DySma), &
                   j, &
                   kmin:kmax &
                   )&
                   )&
                   /( (2*nsmth_DySma + 1) * nsmthv)
               end if
            end do
         end do
      end do

! gagab
  case( "YZ_local_smth" )

       if(yBoundary /= "periodic") &
       stop "DYNAMIC SMAGORINSKY NOT READY FOR NON-PERIODIC BOUNDARY &
             &CONDITIONS IN Y"

       !---------------------------------
       !         Loop over field
       !---------------------------------

       if(nz /= sizeZ) stop" DYNAMIC SMAGORINSKY NOT READY FOR MPI IN Z"

       do k = 1,nz
          ! correct handling of solid and periodic boundaries in z

          if(zBoundary == "solid_wall") then
             kmin=max( 1,k-nsmth_DySma)
             kmax=min(nz,k+nsmth_DySma)
            else if (zBoundary == "periodic") then
             kmin=k-nsmth_DySma
             kmax=k+nsmth_DySma
            else
             stop"vertical smoothing: unknown case zBoundary."
          end if

          nsmthv=kmax-kmin+1
             
          do j = 1,ny
             do i = 1,nx
                ! averaging dyn. Smag. coeff. only over atmosphere cells

                if(topography) then
                  nsmthall=0
                  var3D_DySma(i,j,k)=0.0

                  do ksmth=kmin,kmax
                     do jsmth=j-nsmth_DySma,j+nsmth_DySma
                        if(.not.&
                           topography_mask(i0+i,j0+jsmth,ksmth)) then
                           var3D_DySma(i,j,k)&
                           =var3D_DySma(i,j,k)&
                            +var3D_DySma_Extend(i,jsmth,ksmth)

                           nsmthall=nsmthall+1
                        end if
                     end do
                  end do

                  if(nsmthall > 0) then
                     var3D_DySma(i,j,k)=var3D_DySma(i,j,k)/nsmthall
                  end if
                 else
                   var3D_DySma(i,j,k)&
                   =&
                   sum(&
                   var3D_DySma_Extend( i,&
                   (j-nsmth_DySma):(j+nsmth_DySma), &
                   kmin:kmax &
                   )&
                   )&
                   /( (2*nsmth_DySma + 1) * nsmthv)
               end if
            end do
         end do
      end do
! gagae

    case( "X_whole_smth" )

       if(xBoundary /= "periodic") &
       stop "DYNAMIC SMAGORINSKY NOT READY FOR NON-PERIODIC BOUNDARY &
             &CONDITIONS IN X"

       !---------------------------------
       !         Loop over field
       !---------------------------------

       do k = 1,nz
          do j = 1,ny
             ! averaging dyn. Smag. coeff. only over atmosphere cells

             do i = 1,nx
                if(topography) then
                   var3D_DySma(i,j,k)=0.0
                   nsmthall=0

                   do ismth =1,nx
                      if(.not.topography_mask(i0+ismth,j0+j,k)) then
                         var3D_DySma(i,j,k)&
                         = var3D_DySma(i,j,k) &
                           + var3D_DySma_Extend(ismth,j,k)

                         nsmthall=nsmthall+1
                      end if
                   end do

                   if(nsmthall > 0) then
                      var3D_DySma(i,j,k)=var3D_DySma(i,j,k)/nsmthall
                   end if
                  else
                   var3D_DySma(i,j,k)&
                   =sum(var3D_DySma_Extend( 1:nx, j, k ))/nx
                end if
             end do
          end do
       end do

    case default
       stop"unknown case homog_dir_DySma."
    end select

    ! deallocate local fields
    deallocate(var3D_DySma_Extend, stat=allocstat); if(allocstat/=0) &
         & stop"update.f90:dealloc failed"

    return

  end subroutine Var3DSmthDySma


  subroutine setHaloAndBoundary(var3D_HaloBC, nbx_HaloBC, nby_HaloBC, &
                                & nbz_HaloBC)
    !--------------------------------------
    ! set Halo and Boundary
    !--------------------------------------

    ! in/out variables
    integer, intent(in) :: nbx_HaloBC,nby_HaloBC,nbz_HaloBC

    real, dimension(-nbx_HaloBC:nx+nbx_HaloBC, -nby_HaloBC:ny+nby_HaloBC, &
                  & -nbz_HaloBC:nz+nbz_HaloBC), &
         & intent(inout) :: var3D_HaloBC

    ! auxiliary fields for "var" with ghost cells (rho)
    real, dimension(nbx_HaloBC,-nby_HaloBC:ny+nby_HaloBC,nz) &
            & :: xRhoSliceLeft_send, xRhoSliceRight_send
    real, dimension(nbx_HaloBC,-nby_HaloBC:ny+nby_HaloBC,nz) &
            & :: xRhoSliceLeft_recv, xRhoSliceRight_recv

    real, dimension(-nbx_HaloBC:nx+nbx_HaloBC,nby_HaloBC,nz) &
            & :: yRhoSliceBack_send, yRhoSliceForw_send
    real, dimension(-nbx_HaloBC:nx+nbx_HaloBC,nby_HaloBC,nz) &
            & :: yRhoSliceBack_recv, yRhoSliceForw_recv

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: i, j, k
    integer :: i0, j0, k0

    !------------------------------
    !          x-direction
    !------------------------------

    if( idim > 1 ) then
       call mpi_cart_shift(comm,0,1,left,right,ierror)

       ! slice size
       sendcount = nbx_HaloBC*(ny+2*nby_HaloBC+1)*nz
       recvcount = sendcount

       ! read slice into contiguous array
       do i = 1,nbx_HaloBC
          xRhoSliceLeft_send (i,-nby_HaloBC:ny+nby_HaloBC,1:nz) &
          = var3D_HaloBC(i,-nby_HaloBC:ny+nby_HaloBC,1:nz)

          xRhoSliceRight_send(i,-nby_HaloBC:ny+nby_HaloBC,1:nz) &
          = var3D_HaloBC(nx-nbx_HaloBC+i,-nby_HaloBC:ny+nby_HaloBC,1:nz)
       end do

       ! left -> right
       source = left
       dest = right
       tag = 100

       i0 = 1; j0 = -nby_HaloBC; k0 = 1

       call mpi_sendrecv(xRhoSliceRight_send(i0,j0,k0), sendcount, &
            & mpi_double_precision, dest, tag, &
            & xRhoSliceLeft_recv(i0,j0,k0), recvcount, &
            & mpi_double_precision, &
            & source, mpi_any_tag, comm, sts_left, ierror)

       ! right -> left
       source = right
       dest = left
       tag = 100

       call mpi_sendrecv(xRhoSliceLeft_send(i0,j0,k0), sendcount, &
            & mpi_double_precision, dest, tag, &
            & xRhoSliceRight_recv(i0,j0,k0), recvcount, &
            & mpi_double_precision, &
            & source, mpi_any_tag, comm, sts_right, ierror)

       ! write auxiliary slice to var field
       do i = 1,nbx_HaloBC
          ! right halos
          var3D_HaloBC(nx+i,-nby_HaloBC:ny+nby_HaloBC,1:nz) &
          = xRhoSliceRight_recv(i,-nby_HaloBC:ny+nby_HaloBC,1:nz)

          ! left halos
          var3D_HaloBC(-nbx_HaloBC+i,-nby_HaloBC:ny+nby_HaloBC,1:nz) &
          = xRhoSliceLeft_recv(i,-nby_HaloBC:ny+nby_HaloBC,1:nz)
       end do
      else
       do i = 1,nbx_HaloBC
         var3D_HaloBC(nx+i,:,:) = var3D_HaloBC(i,:,:)
         var3D_HaloBC(-i+1,:,:) = var3D_HaloBC(nx-i+1,:,:)
       end do
    end if

    !------------------------------
    !          y-direction
    !------------------------------

    if( jdim > 1 ) then
       call mpi_cart_shift(comm,1,1,back,forw,ierror)

       ! slice size
       sendcount = nby_HaloBC*(nx+2*nbx_HaloBC+1)*nz
       recvcount = sendcount

       ! read slice into contiguous array
       do j = 1, nby_HaloBC
          yRhoSliceBack_send(-nbx_HaloBC:nx+nbx_HaloBC,j,1:nz) &
          = var3D_HaloBC(-nbx_HaloBC:nx+nbx_HaloBC,j,1:nz)

          yRhoSliceForw_send(-nbx_HaloBC:nx+nbx_HaloBC,j,1:nz) &
          = var3D_HaloBC(-nbx_HaloBC:nx+nbx_HaloBC,ny-nby_HaloBC+j,1:nz)
       end do

       ! back -> forw
       source = back
       dest = forw
       tag = 100
             
       i0 = -nbx_HaloBC; j0 = 1; k0 = 1

       call mpi_sendrecv(yRhoSliceForw_send(i0,j0,k0), sendcount, &
            & mpi_double_precision, dest, tag, &
            & yRhoSliceBack_recv(i0,j0,k0), recvcount, &
            & mpi_double_precision, &
            & source, mpi_any_tag, comm, sts_back, ierror)

       ! forw -> back
       source = forw
       dest = back
       tag = 100

       call mpi_sendrecv(yRhoSliceBack_send(i0,j0,k0), sendcount, &
            & mpi_double_precision, dest, tag, &
            & yRhoSliceForw_recv(i0,j0,k0), recvcount, &
            & mpi_double_precision, &
            & source, mpi_any_tag, comm, sts_forw, ierror)

       ! write auxiliary slice to var field
       do j = 1, nby_HaloBC
          ! right halos
          var3D_HaloBC(-nbx_HaloBC:nx+nbx_HaloBC,ny+j,1:nz) &
          = yRhoSliceForw_recv(-nbx_HaloBC:nx+nbx_HaloBC,j,1:nz)

          ! left halos
          var3D_HaloBC(-nbx_HaloBC:nx+nbx_HaloBC,-nby_HaloBC+j,1:nz) &
          = yRhoSliceBack_recv(-nbx_HaloBC:nx+nbx_HaloBC,j,1:nz)
       end do
      else
       do j = 1,nby_HaloBC
          var3D_HaloBC(:,ny+j,:) = var3D_HaloBC(:,j,:)
          var3D_HaloBC(:,-j+1,:) = var3D_HaloBC(:,ny-j+1,:)
       end do
    end if

    !------------------------------
    !          z-direction
    !------------------------------


    select case( zBoundary )
       
    case( "periodic" ) 

          do k = 1,nbz_HaloBC
             var3D_HaloBC(:,:,nz+k) = var3D_HaloBC(:,:,k)
             var3D_HaloBC(:,:,-k+1) = var3D_HaloBC(:,:,nz-k+1)
          end do
       
    case( "solid_wall" )

          do k = 1,nbz_HaloBC
             var3D_HaloBC(:,:,-k+1) = var3D_HaloBC(:,:,k)
             var3D_HaloBC(:,:,nz+k) = var3D_HaloBC(:,:,nz-k+1)
          end do
       
    case default
       stop"setBoundary: unknown case zBoundary"
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
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var

    ! allocatable fields
    real, dimension(:,:,:), allocatable :: field, field_0, field_1

    integer :: allocstat
    integer :: i,j,k
    integer :: nsmth
    integer :: iVar, ivmax

    allocate(field(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz), stat=allocstat)
    if(allocstat/=0) stop"smooth_shapiro:alloc failed"

    allocate(field_0(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz), stat=allocstat)
    if(allocstat/=0) stop"smooth_shapiro:alloc failed"

    allocate(field_1(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz), stat=allocstat)
    if(allocstat/=0) stop"smooth_shapiro:alloc failed"

    if (timeScheme == "semiimplicit") then
       ! in explicit integration smoothing of density, winds, 
       ! and density fluctuations
       ! pressure fluctuations are not smoothened
       ivmax = 6
      else
       ! in explicit integration smoothing of density, and winds 
       ivmax = 4
    end if

    do iVar = 1,ivmax
       if (iVar == 5) goto 100

       field(:,:,:) = var(:,:,:,iVar)

       ! set the values for field_0

       field_0 = field

       ! start to do the smoothing

       if (sizeX > 1 .and. sizeY > 1 .and. sizeZ > 1) then
          ! 3D smoothing

          if (nbx == 1 .and. nby == 1 .and. nbz == 1) then
             nsmth = 1
            else if (nbx == 2 .and. nby == 2 .and. nbz == 2) then
             nsmth = 2
            else if (nbx == 3 .and. nby == 3 .and. nbz == 3) then
             nsmth = 3
            else if (nbx == 4 .and. nby == 4 .and. nbz == 4) then
             nsmth = 4
            else
             stop'ERROR: wrong nbx, nby, nbz in smoothing'
          end if

          if (nsmth == 1) then
             ! smooth in x

             do k = -nbz,nz+nbz
                do j = -nby,ny+nby
                   do i = 1,nx
                      field_0(i,j,k) &
                      = (field(i-1,j,k) + field(i+1,j,k) &
                         + 2.0*field(i,j,k))/4.0
                   end do
                end do
             end do

             ! smooth in y

             do k = -nbz,nz+nbz
                do j = 1,ny
                   do i = 1,nx
                      field_1(i,j,k) &
                      = (field_0(i,j-1,k) + field_0(i,j+1,k) &
                         + 2.0*field_0(i,j,k))/4.0
                   end do
                end do
             end do
   
             ! smooth in z

             do k = 1,nz
                do j = 1,ny
                   do i = 1,nx
                      field(i,j,k) &
                      = (field_1(i,j,k-1) + field_1(i,j,k+1) &
                         + 2.0*field_1(i,j,k))/4.0
                   end do
                end do
             end do
            elseif (nsmth == 2) then
             ! smooth in x

             do k = -nbz,nz+nbz
                do j = -nby,ny+nby
                   do i = 1,nx
                      field_0(i,j,k) &
                      = (- field(i-2,j,k) - field(i+2,j,k) &
                         + 4.0*(field(i-1,j,k) + field(i+1,j,k)) &
                         + 10.0*field(i,j,k))/16.0
                   end do
                end do
             end do
   
             ! smooth in y
   
             do k = -nbz,nz+nbz
                do j = 1,ny
                   do i = 1,nx
                      field_1(i,j,k) &
                      = (- field_0(i,j-2,k) - field_0(i,j+2,k) &
                         + 4.0*(field_0(i,j-1,k) + field_0(i,j+1,k)) &
                         + 10.0*field_0(i,j,k))/16.0
                   end do
                end do
             end do

             ! smooth in z

             do k = 1,nz
                do j = 1,ny
                   do i = 1,nx
                      field(i,j,k) &
                      = (- field_1(i,j,k-2) - field_1(i,j,k+2) &
                         + 4.0*(field_1(i,j,k-1) + field_1(i,j,k+1)) &
                         + 10.0*field_1(i,j,k))/16.0
                   end do
                end do
             end do
            elseif (nsmth == 3) then
             ! smooth in x

             do k = -nbz,nz+nbz
                do j = -nby,ny+nby
                   do i = 1,nx
                      field_0(i,j,k) &
                      = (field(i-3,j,k) + field(i+3,j,k) &
                         - 6.0*(field(i-2,j,k) + field(i+2,j,k)) &
                         + 15.0*(field(i-1,j,k) + field(i+1,j,k)) &
                         + 44.0*field(i,j,k))/64.0
                   end do
                end do
             end do

             ! smooth in y

             do k = -nbz,nz+nbz
                do j = 1,ny
                   do i = 1,nx
                      field_1(i,j,k) &
                      = (field_0(i,j-3,k) + field_0(i,j+3,k) &
                         - 6.0*(field_0(i,j-2,k) + field_0(i,j+2,k)) &
                         + 15.0*(field_0(i,j-1,k) + field_0(i,j+1,k)) &
                         + 44.0*field_0(i,j,k))/64.0
                   end do
                end do
             end do

             ! smooth in z

             do k = 1,nz
                do j = 1,ny
                   do i = 1,nx
                      field(i,j,k) &
                      = (field_1(i,j,k-3) + field_1(i,j,k+3) &
                         - 6.0*(field_1(i,j,k-2) + field_1(i,j,k+2)) &
                         + 15.0*(field_1(i,j,k-1) + field_1(i,j,k+1)) &
                         + 44.0*field_1(i,j,k))/64.0
                   end do
                end do
             end do
            elseif (nsmth == 4) then
             ! smooth in x

             do k = -nbz,nz+nbz
                do j = -nby,ny+nby
                   do i = 1,nx
                      field_0(i,j,k) &
                      = (- field(i-4,j,k) - field(i+4,j,k) &
                         + 8.0*(field(i-3,j,k) + field(i+3,j,k)) &
                         - 28.0*(field(i-2,j,k) + field(i+2,j,k)) &
                         + 56.0*(field(i-1,j,k) + field(i+1,j,k)) &
                         + 186.0*field(i,j,k))/256.0
                   end do
                end do
             end do
   
             ! smooth in y
   
             do k = -nbz,nz+nbz
                do j = 1,ny
                   do i = 1,nx
                      field_1(i,j,k) &
                      = (- field_0(i,j-4,k) - field_0(i,j+4,k) &
                         + 8.0*(field_0(i,j-3,k) + field_0(i,j+3,k)) &
                         - 28.0*(field_0(i,j-2,k) + field_0(i,j+2,k)) &
                         + 56.0*(field_0(i,j-1,k) + field_0(i,j+1,k)) &
                         + 186.0*field_0(i,j,k))/256.0
                   end do
                end do
             end do

             ! smooth in z

             do k = 1,nz
                do j = 1,ny
                   do i = 1,nx
                      field(i,j,k) &
                      = (field_1(i,j,k-4) + field_1(i,j,k+4) &
                         + 8.0*(field_1(i,j,k-3) + field_1(i,j,k+3)) &
                         - 28.0*(field_1(i,j,k-2) + field_1(i,j,k+2)) &
                         + 56.0*(field_1(i,j,k-1) + field_1(i,j,k+1)) &
                         + 186.0*field_1(i,j,k))/256.0
                   end do
                end do
             end do
          end if
         else if (sizeX > 1 .and. sizeY == 1 .and. sizeZ > 1) then
          ! 2D smoothing in x and z

          if (nbx == 1 .and. nbz == 1) then
             nsmth = 1
            else if (nbx == 2 .and. nbz == 2) then
             nsmth = 2
            else if (nbx == 3 .and. nbz == 3) then
             nsmth = 3
            else if (nbx == 4 .and. nbz == 4) then
             nsmth = 4
            else
             stop'ERROR: wrong nbx, nby, nbz in smoothing'
          end if

          if (nsmth == 1) then
             ! smooth in x
   
             do k = -nbz,nz+nbz
                do j = 1,ny
                   do i = 1,nx
                      field_1(i,j,k) &
                      = (field_0(i-1,j,k) + field_0(i+1,j,k) &
                         + 2.0*field_0(i,j,k))/4.0
                   end do
                end do
             end do

             ! smooth in z
   
             do k = 1,nz
                do j = 1,ny
                   do i = 1,nx
                      field(i,j,k) &
                      = (field_1(i,j,k-1) + field_1(i,j,k+1) &
                         + 2.0*field_1(i,j,k))/4.0
                   end do
                end do
             end do
            elseif (nsmth == 2) then
             ! smooth in x
   
             do k = -nbz,nz+nbz
                do j = 1,ny
                   do i = 1,nx
                      field_1(i,j,k) &
                      = (- field_0(i-2,j,k) - field_0(i+2,j,k) &
                         + 4.0*(field_0(i-1,j,k) + field_0(i+1,j,k)) &
                         + 10.0*field_0(i,j,k))/16.0
                   end do
                end do
             end do

             ! smooth in z

             do k = 1,nz
                do j = 1,ny
                   do i = 1,nx
                      field(i,j,k) &
                      = (- field_1(i,j,k-2) - field_1(i,j,k+2) &
                         + 4.0*(field_1(i,j,k-1) + field_1(i,j,k+1)) &
                         + 10.0*field_1(i,j,k))/16.0
                   end do
                end do
             end do
            elseif (nsmth == 3) then
             ! smooth in x

             do k = -nbz,nz+nbz
                do j = 1,ny
                   do i = 1,nx
                      field_1(i,j,k) &
                      = (field_0(i-3,j,k) + field_0(i+3,j,k) &
                         - 6.0*(field_0(i-2,j,k) + field_0(i+2,j,k)) &
                         + 15.0*(field_0(i-1,j,k) + field_0(i+1,j,k)) &
                         + 44.0*field_0(i,j,k))/64.0
                   end do
                end do
             end do

             ! smooth in z

             do k = 1,nz
                do j = 1,ny
                   do i = 1,nx
                      field(i,j,k) &
                      = (field_1(i,j,k-3) + field_1(i,j,k+3) &
                         - 6.0*(field_1(i,j,k-2) + field_1(i,j,k+2)) &
                         + 15.0*(field_1(i,j,k-1) + field_1(i,j,k+1)) &
                         + 44.0*field_1(i,j,k))/64.0
                   end do
                end do
             end do
            elseif (nsmth == 4) then
             ! smooth in x

             do k = -nbz,nz+nbz
                do j = 1,ny
                   do i = 1,nx
                      field_1(i,j,k) &
                      = (- field_0(i-4,j,k) - field_0(i+4,j,k) &
                         + 8.0*(field_0(i-3,j,k) + field_0(i+3,j,k)) &
                         - 28.0*(field_0(i-2,j,k) + field_0(i+2,j,k)) &
                         + 56.0*(field_0(i-1,j,k) + field_0(i+1,j,k)) &
                         + 186.0*field_0(i,j,k))/256.0
                   end do
                end do
             end do
   
             ! smooth in z
      
             do k = 1,nz
                do j = 1,ny
                   do i = 1,nx
                      field(i,j,k) &
                      = (field_1(i,j,k-4) + field_1(i,j,k+4) &
                         + 8.0*(field_1(i,j,k-3) + field_1(i,j,k+3)) &
                         - 28.0*(field_1(i,j,k-2) + field_1(i,j,k+2)) &
                         + 56.0*(field_1(i,j,k-1) + field_1(i,j,k+1)) &
                         + 186.0*field_1(i,j,k))/256.0
                   end do
                end do
             end do
          end if
         else if (sizeX == 1 .and. sizeY > 1 .and. sizeZ > 1) then
          if (nby == 1 .and. nbz == 1) then
             nsmth = 1
            else if (nby == 2 .and. nbz == 2) then
             nsmth = 2
            else if (nby == 3 .and. nbz == 3) then
             nsmth = 3
            else if (nby == 4 .and. nbz == 4) then
             nsmth = 4
            else
             stop'ERROR: wrong nbx, nby, nbz in smoothing'
          end if

          if (nsmth == 1) then
             ! smooth in y
   
             do k = -nbz,nz+nbz
                do j = 1,ny
                   do i = 1,nx
                      field_1(i,j,k) &
                      = (field_0(i,j-1,k) + field_0(i,j+1,k) &
                         + 2.0*field_0(i,j,k))/4.0
                   end do
                end do
             end do

             ! smooth in z
   
             do k = 1,nz
                do j = 1,ny
                   do i = 1,nx
                      field(i,j,k) &
                      = (field_1(i,j,k-1) + field_1(i,j,k+1) &
                         + 2.0*field_1(i,j,k))/4.0
                   end do
                end do
             end do
            elseif (nsmth == 2) then
             ! smooth in y

             do k = -nbz,nz+nbz
                do j = 1,ny
                   do i = 1,nx
                      field_1(i,j,k) &
                      = (- field_0(i,j-2,k) - field_0(i,j+2,k) &
                         + 4.0*(field_0(i,j-1,k) + field_0(i,j+1,k)) &
                         + 10.0*field_0(i,j,k))/16.0
                   end do
                end do
             end do

             ! smooth in z

             do k = 1,nz
                do j = 1,ny
                   do i = 1,nx
                      field(i,j,k) &
                      = (- field_1(i,j,k-2) - field_1(i,j,k+2) &
                         + 4.0*(field_1(i,j,k-1) + field_1(i,j,k+1)) &
                         + 10.0*field_1(i,j,k))/16.0
                   end do
                end do
             end do
            elseif (nsmth == 3) then
             ! smooth in y
   
             do k = -nbz,nz+nbz
                do j = 1,ny
                   do i = 1,nx
                      field_1(i,j,k) &
                      = (field_0(i,j-3,k) + field_0(i,j+3,k) &
                         - 6.0*(field_0(i,j-2,k) + field_0(i,j+2,k)) &
                         + 15.0*(field_0(i,j-1,k) + field_0(i,j+1,k)) &
                         + 44.0*field_0(i,j,k))/64.0
                   end do
                end do
             end do

             ! smooth in z

             do k = 1,nz
                do j = 1,ny
                   do i = 1,nx
                      field(i,j,k) &
                      = (field_1(i,j,k-3) + field_1(i,j,k+3) &
                         - 6.0*(field_1(i,j,k-2) + field_1(i,j,k+2)) &
                         + 15.0*(field_1(i,j,k-1) + field_1(i,j,k+1)) &
                         + 44.0*field_1(i,j,k))/64.0
                   end do
                end do
             end do
            elseif (nsmth == 4) then
             ! smooth in y

             do k = -nbz,nz+nbz
                do j = 1,ny
                   do i = 1,nx
                      field_1(i,j,k) &
                      = (- field_0(i,j-4,k) - field_0(i,j+4,k) &
                         + 8.0*(field_0(i,j-3,k) + field_0(i,j+3,k)) &
                         - 28.0*(field_0(i,j-2,k) + field_0(i,j+2,k)) &
                         + 56.0*(field_0(i,j-1,k) + field_0(i,j+1,k)) &
                         + 186.0*field_0(i,j,k))/256.0
                   end do
                end do
             end do

             ! smooth in z
   
             do k = 1,nz
                do j = 1,ny
                   do i = 1,nx
                      field(i,j,k) &
                      = (field_1(i,j,k-4) + field_1(i,j,k+4) &
                         + 8.0*(field_1(i,j,k-3) + field_1(i,j,k+3)) &
                         - 28.0*(field_1(i,j,k-2) + field_1(i,j,k+2)) &
                         + 56.0*(field_1(i,j,k-1) + field_1(i,j,k+1)) &
                         + 186.0*field_1(i,j,k))/256.0
                   end do
                end do
             end do
          end if
         else
          stop"ERROR: smoothing not ready for 2D in x and y or 1D"
       end if

       var(:,:,:,iVar) = field(:,:,:)

100    continue
    end do

    ! deallocate local fields
    deallocate(field, stat=allocstat); if(allocstat/=0) &
         & stop"smooth_shapiro:dealloc failed"
    deallocate(field_0, stat=allocstat); if(allocstat/=0) &
         & stop"smooth_shapiro:dealloc failed"
    deallocate(field_1, stat=allocstat); if(allocstat/=0) &
         & stop"smooth_shapiro:dealloc failed"

    return

  end subroutine smooth_shapiro_0

!UAB
!------------------------------------------------------------------------

  subroutine smooth_hor_shapiro(fc_shap,n_shap,flux,var)
 
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
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real, intent(in) :: fc_shap
    integer, intent(in) :: n_shap
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), &
         & intent(inout) :: flux

    ! allocatable fields
    real, dimension(:,:,:), allocatable :: field
    real, dimension(:,:,:,:), allocatable :: var_l

    integer :: allocstat
    integer :: i,j,k
    integer :: nsmth
    integer :: iVar, ivmax
    integer :: i_lapl
    integer :: nz_max

    allocate(field(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz), stat=allocstat)
    if(allocstat/=0) stop"smooth_shapiro:alloc failed"

    allocate(var_l(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
           & stat=allocstat)
    if(allocstat/=0) stop"smooth_shapiro:alloc failed"

    if (timeScheme == "semiimplicit") then
       ! in explicit integration smoothing of density, winds, 
       ! and density fluctuations
       ! pressure fluctuations are not smoothened
       ivmax = 4!6 !FS
      else
       ! in explicit integration smoothing of density, and winds 
       ivmax = 4
    end if


   
    var(1:nx,1:ny,1:nz,2) = var(1:nx,1:ny,1:nz,2) - u_env_pp(1:nx,1:ny,1:nz) !FS
    
    ! do i = 1,nx
    !    do j = 1, ny
    !       var(i,j,1:nz,1) = var(i,j,1:nz,1) - dens_env_pp(i,j,1:nz) + rhoStrat(1:nz)
    !    end do
    ! end do
 ! make sure that boundary conditions are satisfied

    call setHalos( var, "var" )
    call setBoundary (var,flux,"var")

    ! smoothing in x-direction

    if (sizeX > 1) then
       ! 2n-th-order x-derivative of all fields that are to be smoothed

       var_l = var

       do i_lapl = 1, n_shap
          do iVar = 2,ivmax !FS
             if (iVar /= 5) then
          

                field(:,:,:) = var_l(:,:,:,iVar)

                do k = 1,nz
                   do j = 1,ny
                      do i = 1,nx
                         var_l(i,j,k,iVar) &
                         = (  field(i-1,j,k) + field(i+1,j,k) &
                            - 2.0*field(i,j,k))/4.0
                      end do
                   end do
                end do
             end if
          end do

          ! horizontal boundary conditions so that everything is ready for 
          ! the next iteration

          !for non-parallelized directions

          call setHalos( var_l, "var" )

          ! non-parallel boundary conditions in x-direction

          select case( xBoundary )
          case( "periodic" )
            if ( idim == 1 ) call setBoundary_x_periodic(var_l,flux,"var")
          case default
             stop"setBoundary: unknown case xBoundary"
          end select
       end do

       ! apply filter
    
       do iVar = 2,ivmax!FS
          if (iVar /= 5) then
             if (iVar == 4) then
                nz_max = nz-1
               else
                nz_max = nz
             end if

             do k = 1,nz_max
                do j = 1,ny
                   do i = 1,nx
                     
                      var(i,j,k,iVar) &
                      = var(i,j,k,iVar) &
                        + fc_shap * (-1)**(n_shap + 1) * var_l(i,j,k,iVar)
                     
                   end do
                end do
             end do
          end if
       end do

       ! boundary conditions again

       call setHalos( var, "var" )
       call setBoundary (var,flux,"var")
    end if

    !testb
    !goto 100
    !teste

    if (sizeY > 1) then
       ! 2n-th-order y-derivative of all fields that are to be smoothed

       var_l = var

       do i_lapl = 1, n_shap
          do iVar = 2,ivmax!FS
             if (iVar /= 5) then

                          
                field(:,:,:) = var_l(:,:,:,iVar)

                do k = 1,nz
                   do j = 1,ny
                      do i = 1,nx
                         var_l(i,j,k,iVar) &
                         = (  field(i,j-1,k) + field(i,j+1,k) &
                            - 2.0*field(i,j,k))/4.0
                      end do
                   end do
                end do
             end if

                      
          end do

          ! horizontal boundary conditions so that everything is ready for 
          ! the next iteration

          !for non-parallelized directions

          call setHalos( var_l, "var" )

          ! non-parallel boundary conditions in y-direction
    
          select case( yBoundary )
          case( "periodic" ) 
            if ( jdim == 1 ) call setBoundary_y_periodic(var_l,flux,"var")
          case default
             stop"setBoundary: unknown case yBoundary"
          end select   
       end do

       ! apply filter
    
       do iVar = 2,ivmax!FS
          if (iVar /= 5) then
             if (iVar == 4) then
                nz_max = nz-1
               else
                nz_max = nz
             end if

             do k = 1,nz_max
                do j = 1,ny
                   do i = 1,nx

                      var(i,j,k,iVar) &
                      = var(i,j,k,iVar) &
                        + fc_shap * (-1)**(n_shap + 1) * var_l(i,j,k,iVar)
                    
                   end do
                end do
             end do
          end if
       end do
 
    var(1:nx,1:ny,1:nz,2) = var(1:nx,1:ny,1:nz,2) + u_env_pp(1:nx,1:ny,1:nz) !FS
    
    ! do i = 1,nx
    !    do j = 1, ny
    !       var(i,j,1:nz,1) = var(i,j,1:nz,1) + dens_env_pp(i,j,1:nz) - rhoStrat(1:nz)
    !    end do
    ! end do

       ! boundary conditions again
   
       call setHalos( var, "var" )
       call setBoundary (var,flux,"var")
    end if

    !testb
100 continue
    !teste

    ! deallocate local fields

    deallocate(field, stat=allocstat); if(allocstat/=0) &
         & stop"smooth_shapiro:dealloc failed"
    deallocate(var_l, stat=allocstat); if(allocstat/=0) &
         & stop"smooth_shapiro:dealloc failed"

    return

  end subroutine smooth_hor_shapiro

!------------------------------------------------------------------------

  subroutine smooth_shapiro(fc_shap,n_shap,flux,var)
 
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
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real, intent(in) :: fc_shap
    integer, intent(in) :: n_shap
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), &
         & intent(inout) :: flux

    ! allocatable fields
    real, dimension(:,:,:), allocatable :: field
    real, dimension(:,:,:,:), allocatable :: var_l

    integer :: allocstat
    integer :: i,j,k
    integer :: nsmth
    integer :: iVar, ivmax
    integer :: i_lapl
    integer :: nz_max

    allocate(field(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz), stat=allocstat)
    if(allocstat/=0) stop"smooth_shapiro:alloc failed"

    allocate(var_l(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
           & stat=allocstat)
    if(allocstat/=0) stop"smooth_shapiro:alloc failed"

    if (timeScheme == "semiimplicit") then
       ! in semiimplicit integration smoothing of density, winds, 
       ! and density fluctuations
       ! pressure fluctuations are not smoothened
       ivmax = 4!6!FS 
      else
       ! in explicit integration smoothing of density, and winds 
       ivmax = 4
    end if

    ! make sure that boundary conditions are satisfied

    call setHalos( var, "var" )
    call setBoundary (var,flux,"var")

    ! smoothing in x-direction

    if (sizeX > 1) then
       ! 2n-th-order x-derivative of all fields that are to be smoothed

       var_l = var

       do i_lapl = 1, n_shap 
          do iVar = 2,ivmax!FS
             if (iVar /= 5) then

                        

                field(:,:,:) = var_l(:,:,:,iVar)

                do k = 1,nz
                   do j = 1,ny
                      do i = 1,nx
                         var_l(i,j,k,iVar) &
                         = (  field(i-1,j,k) + field(i+1,j,k) &
                            - 2.0*field(i,j,k))/4.0
                      end do
                   end do
                end do
             end if

                       
          end do

          ! horizontal boundary conditions so that everything is ready for 
          ! the next iteration

          !for non-parallelized directions

          call setHalos( var_l, "var" )

          ! non-parallel boundary conditions in x-direction

          select case( xBoundary )
          case( "periodic" )
            if ( idim == 1 ) call setBoundary_x_periodic(var_l,flux,"var")
          case default
             stop"setBoundary: unknown case xBoundary"
          end select
       end do

       ! apply filter
    
       do iVar = 2,ivmax!FS
          if (iVar /= 5) then
             if (iVar == 4) then
                nz_max = nz-1
               else
                nz_max = nz
             end if

             do k = 1,nz_max
                do j = 1,ny
                   do i = 1,nx
                    
                      var(i,j,k,iVar) &
                      = var(i,j,k,iVar) &
                        + fc_shap * (-1)**(n_shap + 1) * var_l(i,j,k,iVar)
                    
                   end do
                end do
             end do
          end if
       end do

       ! boundary conditions again

       call setHalos( var, "var" )
       call setBoundary (var,flux,"var")
    end if

    !testb
    !goto 100
    !teste

    !smoothing in z-direction

    if (sizeY > 1) then
       ! 2n-th-order y-derivative of all fields that are to be smoothed

       var_l = var
       
       do i_lapl = 1, n_shap
          do iVar = 2,ivmax!FS
             if (iVar /= 5) then

                        

                field(:,:,:) = var_l(:,:,:,iVar)

                do k = 1,nz
                   do j = 1,ny
                      do i = 1,nx
                         var_l(i,j,k,iVar) &
                         = (  field(i,j-1,k) + field(i,j+1,k) &
                            - 2.0*field(i,j,k))/4.0
                      end do
                   end do
                end do
             end if

          end do

          ! horizontal boundary conditions so that everything is ready for 
          ! the next iteration

          !for non-parallelized directions

          call setHalos( var_l, "var" )

          ! non-parallel boundary conditions in y-direction
    
          select case( yBoundary )
          case( "periodic" ) 
            if ( jdim == 1 ) call setBoundary_y_periodic(var_l,flux,"var")
          case default
             stop"setBoundary: unknown case yBoundary"
          end select   
       end do

       ! apply filter
    
       do iVar = 2,ivmax!FS
          if (iVar /= 5) then
             if (iVar == 4) then
                nz_max = nz-1
               else
                nz_max = nz
             end if

             do k = 1,nz_max
                do j = 1,ny
                   do i = 1,nx
                    
                      var(i,j,k,iVar) &
                      = var(i,j,k,iVar) &
                      + fc_shap * (-1)**(n_shap + 1) * var_l(i,j,k,iVar)
                     
                   end do
                end do
             end do
          end if
       end do

       ! boundary conditions again
   
       call setHalos( var, "var" )
       call setBoundary (var,flux,"var")
    end if

    ! filtering in z-direction

    ! 2nd-order z-derivative of all fields that are to be smoothed

    var_l = var

    do iVar = 2,ivmax!FS
       if (iVar /= 5) then

                   

          field(:,:,:) = var_l(:,:,:,iVar)

          do k = 1,nz
             do j = 1,ny
                do i = 1,nx
                   var_l(i,j,k,iVar) &
                   = (  field(i,j,k-1) + field(i,j,k+1) &
                      - 2.0*field(i,j,k))/4.0
                end do
             end do
          end do
       end if

                
    end do

    ! apply filter
    
    do iVar = 2,ivmax!FS
       if (iVar /= 5) then
          if (iVar == 4) then
             nz_max = nz-1
            else
             nz_max = nz
          end if

          do k = 1,nz_max
             do j = 1,ny
                do i = 1,nx
                  
                   var(i,j,k,iVar) &
                   = var(i,j,k,iVar) &
                     + 1.e-2 * fc_shap * var_l(i,j,k,iVar)
                
                end do
             end do
          end do
       end if
    end do

    ! boundary conditions again
   
    call setHalos( var, "var" )
    call setBoundary (var,flux,"var")

    !testb
100 continue
    !teste

    ! deallocate local fields

    deallocate(field, stat=allocstat); if(allocstat/=0) &
         & stop"smooth_shapiro:dealloc failed"
    deallocate(var_l, stat=allocstat); if(allocstat/=0) &
         & stop"smooth_shapiro:dealloc failed"

    return

  end subroutine smooth_shapiro
!UAE

!---------------------------------------------------------------------

  !UAB
  !subroutine BGstate_update(var,flux,dt,m,w_0,q_P,q_rho,int_mod)
  subroutine BGstate_update(var,flux,flux_rhopw,dt,m,q_P,q_rho,int_mod)
  !UAE

  ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), intent(in) :: flux
    real, dimension(-1:nx,-1:ny,-1:nz), intent(in) :: flux_rhopw
    real,intent(in) :: dt
    integer, intent(in) :: m

    !UAB
    !real, dimension(-nbz:nz+nbz),intent(inout) :: w_0 
    !UAE

    real, dimension(-nbz:nz+nbz),intent(inout) :: q_P, q_rho

    character(len=*), intent(in) :: int_mod

    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz) :: heat
    real, dimension(-nbz:nz+nbz) :: S_bar, PStratold
    integer :: i,j,k
    real, dimension(1:nz) :: sum_local, sum_global 
    real, dimension(-nbz:nz+nbz) :: press0 
    real, dimension(-nbz:nz+nbz) :: divPw, divrhow !UA

    real, dimension(-nbz:nz+nbz) :: w_0 !UA

    real, dimension(-nbz:nz+nbz) :: avgrhopw
    real :: dptopdt

    real :: sum_d, sum_n
    

    w_0 = 0.
    S_bar = 0.
    heat = 0.
    sum_local = 0.
    sum_global = 0.

    !UAB
    divPw = 0.
    divrhow = 0.
    !UAE

    ! S eq(9)  ONeill+Klein2014
    !UAB call calculate_heating_ONK14(var,flux,heat) 
    call calculate_heating(var,flux,heat) 
    !UAE


    ! calculate horizontal mean of heat(:,:,:)
    do k = 1,nz
       sum_local(k) = sum(heat(1:nx,1:ny,k))
    end do
    !global sum and average
    call mpi_allreduce(sum_local(1),sum_global(1),&
         nz-1+1,&
         mpi_double_precision,mpi_sum,comm,ierror)
    sum_global = sum_global/(sizeX*sizeY)

    S_bar(1:nz) = (-1)*sum_global(1:nz)


    do k = 1,nz
       press0(k) = PStrat(k)**gamma  
    end do
    !UAE

    sum_d = 0.0
    sum_n = 0.0

    !calculate horizontal mean of vertical rhop flux 
    sum_local(:) = 0.
    sum_global(:) = 0.
    avgrhopw = 0.
 
    if(timeScheme == "semiimplicit")then
       do k = 1,nz
          sum_local(k) =  sum(flux(1:nx,1:ny,k,3,6))
       end do
    else
       do k = 1,nz
          sum_local(k) = sum(flux_rhopw(1:nx,1:ny,k))
       end do
    end if
       
       !global sum and average
    call mpi_allreduce(sum_local(1),sum_global(1),&
         nz-1+1,&
         mpi_double_precision,mpi_sum,comm,ierror)
    sum_global = sum_global/(sizeX*sizeY)
    
    avgrhopw(1:nz) = sum_global

    do k = 1,nz
       sum_n = sum_n + S_bar(k)/PStrat(k) - avgrhopw(k)*g_ndim/(gamma*press0(k)) 
       sum_d = sum_d +  1./(gamma*press0(k))
    end do

    dptopdt = sum_n/sum_d
    !UAE


    w_0(1) = dz*(S_bar(1)/Pstrat(1) - (1./(gamma*press0(1)))*dptopdt- avgrhopw(1)*g_ndim/(gamma*press0(1)) )

    do k = 2,nz-1
       w_0(k) &
       = w_0(k-1) &
         + dz*(S_bar(k)/Pstrat(k) - (1./(gamma*press0(k)))*dptopdt- avgrhopw(k)*g_ndim/(gamma*press0(k))  )
    end do
    !UAE

    ! update PStrat and rhoStrat and thetaStrat

    !UAB
    do k=1,nz
       divPw(k) = (PstratTilde(k)*w_0(k) - PstratTilde(k-1)*w_0(k-1))/dz
       divrhow(k) &
       = (rhoStratTilde(k)*w_0(k) - rhoStratTilde(k-1)*w_0(k-1) + avgrhopw(k) - avgrhopw(k-1))/dz
    end do
    !UAE

    !UAB
    ! save total density and subtract the reference-atmosphere density from 
    ! this again after the update of the latter

    if (fluctuationMode) then
       do k = 1,nz
          var(:,:,k,1) = var(:,:,k,1) + rhoStrat(k)
       end do
    end if

    if (timeScheme == "semiimplicit") then
       do k = 1,nz
          var(:,:,k,6) = var(:,:,k,6) + rhoStrat(k)
       end do
    end if

    ! background only changed where w_0 /= 0
    !do k = -1,nz+1
    do k = 1,nz
    !UAE
       if(int_mod == "expl")then
          !init q
          if (m == 1) then
             q_P(k) = 0.
             q_rho(k) = 0.
          end if
          
          ! update: q(m-1) -> q(m)


          q_P(k) = alpha(m) * q_P(k) - dt*divPw(k) + dt*S_bar(k) 
          !UAE

          ! update PStrat

          Pstrat(k) = Pstrat(k) + beta(m) * q_P(k) 

          ! update: q(m-1) -> q(m)


          q_rho(k) = alpha(m) * q_rho(k) - dt*divrhow(k)
          !UAE

          ! update rhoStrat

          rhostrat(k) = rhostrat(k) + beta(m) * q_rho(k)
         else if (int_mod == "impl")then
          !UAB
      
          PStrat(k) = PStrat(k) - dt*divPw(k) + dt*S_bar(k)
          rhoStrat(k) = rhoStrat(k) - dt*divrhow(k)
          !UAE
         else
          print*, "update.f90/BGstate_update: wrong int_mod"
          stop
       end if
      
       !testb
       !print*,'k,Pstrat(k),rhoStrat(k)'
       !print*,k,Pstrat(k),rhoStrat(k)
       !teste

       !update thetaStrat and piSTrat
       pistrat(k) = PStrat(k)**(kappa/(1.0 - kappa))
       thetaStrat(k) = PStrat(k)/rhoStrat(k)

    end do

    !UAB
    do k = -1,nz+1
       PstratTilde(k) = 0.5*(PStrat(k)+PStrat(k+1))
       rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k+1))
       thetaStratTilde(k) =  PStratTilde(k)/rhoStratTilde(k)
    end do

    ! adjust density fluctuations to new reference atmosphere
    if (fluctuationMode) then
       do k = 1,nz
          var(:,:,k,1) = var(:,:,k,1) - rhoStrat(k)
       end do
    end if

    if (timeScheme == "semiimplicit") then
       do k = 1,nz
          var(:,:,k,6) = var(:,:,k,6) - rhoStrat(k)
       end do
    end if
    !UAE

    !UAB
    ! update of non-dimensional squared Brunt-Vaisala frequency
    ! (this could perhaps be done a bit nicer)

    bvsStrat(-1) &
    = g_ndim/thetaStrat(0) * (thetaStrat(1) - thetaStrat(0))/dz

    bvsStrat(0) &
    = g_ndim/thetaStrat(0) * (thetaStrat(1) - thetaStrat(0))/dz

    N2 = max(bvsStrat(-1),bvsStrat(0))

    do k = 1,nz
       bvsStrat(k) &
       = g_ndim/thetaStrat(k) &
         * (thetaStrat(k+1) - thetaStrat(k-1))/(2.0 * dz)
          
       N2 = max(N2, bvsStrat(k))
    end do

    bvsStrat(nz+1) &
    = g_ndim/thetaStrat(nz+1) * (thetaStrat(nz+1) - thetaStrat(nz))/dz

    N2 = max(N2, bvsStrat(nz+1))

    if(N2 < 0.) then
       stop'ERROR: N2 < 0'
      else
       NN = sqrt(N2)
    end if

    !testb
    do k = -1,nz+1
       if (master .and. N2 == bvsStrat(k)) print*,'N2 = max at k =',k
    end do
    !teste
    !UAE

    
    ! call output_profile(iOut,w_0,'wStrat.dat')
    ! call mpi_barrier(comm,ierror)
    ! call output_profile(iOut,rhoStrat,'rhoStrat.dat')
    ! call mpi_barrier(comm,ierror)
    ! call output_profile(iOut,thetaStrat,'thetaStrat.dat')
    ! call mpi_barrier(comm,ierror)
    ! call output_profile(iOut,PStrat,'PStrat.dat')
    ! call mpi_barrier(comm,ierror)
    ! call output_profile(iOut,piStrat,'piStrat.dat')
    ! call mpi_barrier(comm,ierror)
    ! call output_profile(iOut,bvsStrat,'bvsStrat.dat')
    ! call mpi_barrier(comm,ierror)

    S_bar = (-1)*S_bar


  end subroutine BGstate_update





end module update_module

