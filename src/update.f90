module update_module

  use type_module
  use timeScheme_module
  use atmosphere_module
  use flux_module
  use algebra_module


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

  public :: CoefDySma_update         ! modified by Junhong Wei (20160726)
  public :: Var3DSmthDySma           ! modified by Junhong Wei (20160916)

  public :: setHaloAndBoundary       ! modified by Junhong Wei (20170203)

  
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
    integer :: i,j,k

    ! relaxation parameters
    real :: alpha, beta
    real :: spongeAlphaZ, spongeDz
    
    ! variables for rho 
    real    :: rho_old, rho_bg, rho_new
    real    :: uOld, uBG, uNew
    real    :: vOld, vBG, vNew
    real    :: wOld, wBG, wNew

!   achatzb
    real, dimension(1:nz) :: sum_local, sum_global
!   achatze
    
    ! return if sponge layer with relaxation is switched off
    if( .not. spongeLayer ) then 
       return
!   else
!      print*,"WARNING: Sponge layer with relaxation is on!"
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

                if( fluctuationMode ) then
                   rho_bg = 0.0               ! push back to zero perturbation
                else
                   rho_bg = rhoStrat(k)
                end if
                
                rho_old = var(i,j,k,1)
                alpha = spongeAlphaZ*(z(k)-zSponge)/spongeDz
                beta = 1./(1.+alpha*0.5*dt)**2
                rho_new = (1.-beta)*rho_bg + beta*rho_old
                
                var(i,j,k,1) = rho_new
                
             end do
          end do
       end do

       
       
    case( "uvw" )
       
!      achatzb relax u to horizontal mean
!      ! relax u to background flow
!      do k = kSponge, nz
!         do j = 1,ny
!            do i = 0,nx
!               
!               uOld = var(i,j,k,2)
!               uBG = backgroundFlow(1)
!               alpha = spongeAlphaZ * (z(k)-zSponge) / spongeDz
!               beta = 1./(1.+alpha*0.5*dt)**2
!               
!               uNew = (1.-beta)*uBG + beta*uOld
!               
!               var(i,j,k,2) = uNew
!               
!            end do
!         end do
!      end do

       ! local horizontal sum in the sponge layer

       do k = kSponge, nz
          sum_local(k) = sum(var(1:nx,1:ny,k,2))
       end do

       ! global sum and average

       call mpi_allreduce(sum_local(kSponge),sum_global(kSponge),nz-kSponge+1,&
                          mpi_double_precision,mpi_sum,comm,ierror)
       sum_global = sum_global/(sizeX*sizeY)

       do k = kSponge, nz
          uBG = sum_global(k)

          do j = 1,ny
             do i = 0,nx
                
                uOld = var(i,j,k,2)
                alpha = spongeAlphaZ * (z(k)-zSponge) / spongeDz
                beta = 1./(1.+alpha*0.5*dt)**2
                
                uNew = (1.-beta)*uBG + beta*uOld
                
                var(i,j,k,2) = uNew
                
             end do
          end do
       end do
!      achatze
       
!      achatzb relax v to horizontal mean
!      ! relax v to background flow
!      do k = kSponge, nz
!         do j = 1,ny
!            do i = 1,nx
!               
!               vOld = var(i,j,k,3)
!               vBG = backgroundFlow(2)
!               alpha = spongeAlphaZ * (z(k)-zSponge) / spongeDz
!               beta = 1./(1.+alpha*0.5*dt)**2
!               
!               vNew = (1.-beta)*vBG + beta*vOld
!               
!               var(i,j,k,3) = vNew
!               
!            end do
!         end do
!      end do
       
       ! local horizontal sum in the sponge layer

       do k = kSponge, nz
          sum_local(k) = sum(var(1:nx,1:ny,k,3))
       end do

       ! global sum and average

       call mpi_allreduce(sum_local(kSponge),sum_global(kSponge),nz-kSponge+1,&
                          mpi_double_precision,mpi_sum,comm,ierror)
       sum_global = sum_global/(sizeX*sizeY)

       do k = kSponge, nz
          vBG = sum_global(k)

          do j = 1,ny
             do i = 1,nx
                
                vOld = var(i,j,k,3)
                alpha = spongeAlphaZ * (z(k)-zSponge) / spongeDz
                beta = 1./(1.+alpha*0.5*dt)**2
                
                vNew = (1.-beta)*vBG + beta*vOld
                
                var(i,j,k,3) = vNew
                
             end do
          end do
       end do
!      achatze
       
!      achatzb relax w to zero
!      ! relax w to background flow
!      do k = kSponge, nz
!         do j = 1,ny
!            do i = 1,nx
!               
!               wOld = var(i,j,k,4)
!               wBG = backgroundFlow(3)
!               alpha = spongeAlphaZ * (z(k)-zSponge) / spongeDz
!               beta = 1./(1.+alpha*0.5*dt)**2
!               
!               wNew = (1.-beta)*wBG + beta*wOld

!               var(i,j,k,4) = wNew
!               
!            end do
!         end do
!      end do

       do k = kSponge, nz
          wBG = 0.0

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
!      achatze


    case default
       stop "spongeLayer: Unknown variable"
    end select
    

  end subroutine set_spongeLayer


!---------------------------------------------------------------------


  subroutine momentumPredictor (var,var0,flux,source, force,dt,q,m)
    !----------------------------------
    !  calculates the velocities u^*
    !----------------------------------
    
    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var0

    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux) 
    ! dir = 1..3 > f,g,h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW

    ! source terms
    ! 1) divergence error source terms
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: source

    ! volume forces (cell centered) 
    ! 1) gravitational / buoyancy force
    ! 2) Coriolis force
    real, dimension(0:nx+1,0:ny+1,0:nz+1,3), intent(in) :: force


    real, intent(in) :: dt 
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,3), &
         & intent(inout) :: q
    integer, intent(in) :: m
 

    ! local variables
    real :: fL, fR, gB,gF, hD,hU
    ! flux Left/Right, Backward/Forward, Downward/Upward
    
    ! other stuff
    real :: rhoM_1, rhoM               ! rho(m-1), rho(m)
    real :: uM_1, vM_1, wM_1           ! u(m-1), v(m-1) and w(m-1)
    real :: momM_1, momM               ! momentum at t(m-1) and t(m)
    real :: piR,piL, piF,piB, piU,piD  
    real :: fluxDiff                   ! conv. and viscous flux contribution
    real :: piGrad                     ! pressure gradient
    real :: F                          ! update part for Runge-Kutta step
    real :: uAst, vAst, wAst           ! predicted velocities u*, v* and w*
    
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

    ! test 
    real :: u,v,w,uAbs

!   achatzb
    integer :: i00,j00
!   achatze

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
       stop "momentumPredictor: unknown case xBoundary."
    end select

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
             fluxDiff = (fR-fL)/dx + (gF-gB)/dy + (hU-hD)/dz     ! divergence

             !--- pressure gradient term -> piGrad
             piR = var(i+1,j,k,5)
             piL = var(i,  j,k,5)
             piGrad = kappaInv*MaInv2 * Pstrat(k) * (piR-piL)/dx

             !---- volume forces
             volForce = 0.5*( force(i,j,k,1) + force(i+1,j,k,1) )
        
             !--------------------
             !   d/dt ... = F(phi) (RHS of ODE)
             !--------------------
             ! fluxDiff -> convective and viscous fluxes
             ! piGrad   -> pressure gradient along x scaled with 1/Ma^2
             ! volForce -> Gravity, Coriolis
             F = -fluxDiff - piGrad + volForce
             
!            testb
!            if((i == 0).and.(j == 1).and.(k == 310)) then
!                print*,"fluxDiff, piGrad, volForce at (0,1,310) = 0:", &
!                      & fluxDiff, piGrad, volForce
!                print*,"fL, fR = ",fL, fR
!                print*,"gB, gF = ",gB, gF
!                print*,"hD, hU = ",hD, hU
!                print*,"force(i,j,k,1), force(i+1,j,k,1)", &
!                     & force(i,j,k,1), force(i+1,j,k,1)
!            end if

!            if((i == 1).and.(j == 1).and.(k == 310)) then
!                print*,"fluxDiff, piGrad, volForce at (1,1,310) = 0:", &
!                      & fluxDiff, piGrad, volForce
!                print*,"fL, fR = ",fL, fR
!                print*,"gB, gF = ",gB, gF
!                print*,"hD, hU = ",hD, hU
!                print*,"force(i,j,k,1), force(i+1,j,k,1)", &
!                     & force(i,j,k,1), force(i+1,j,k,1)
!            end if

!            if((i == 32).and.(j == 1).and.(k == 310)) then
!                print*,"fluxDiff, piGrad, volForce at (32,1,310) = 0:", &
!                      & fluxDiff, piGrad, volForce
!                print*,"fL, fR = ",fL, fR
!                print*,"gB, gF = ",gB, gF
!                print*,"hD, hU = ",hD, hU
!                print*,"force(i,j,k,1), force(i+1,j,k,1)", &
!                     & force(i,j,k,1), force(i+1,j,k,1)
!            end if
!            teste

             if( correctDivError ) F = F + source(i,j,k,2)
             
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
                stop "momentumPredictor: unkown case model."
             end select
             
             ! velocity and momentum at t(m-1)
             uM_1 = var(i,j,k,2)
             momM_1 = rhoM_1 * uM_1
             
             ! q(m-1) -> q(m)
             
             select case( timeSchemeType ) 

             case( "lowStorage" ) 

                q(i,j,k,1) = dt * F + alpha(m) * q(i,j,k,1)

                ! rhoU(m-1) -> rhoU(m)
                momM = momM_1 + beta(m) * q(i,j,k,1)

                ! calc u(m,*)
                uAst = momM / rhoM

                ! uAst -> var
                var(i,j,k,2) = uAst

             case( "classical" )

                ! interpolated density at t(0)
                select case( model ) 

                case( "pseudo_incompressible" )
                   rhoM_0   = 0.5 * ( var0(i,j,k,1) + var0(i+1,j,k,1) )
                   if( fluctuationMode ) rhoM_0 = rhoM_0 + rhoStrat(k)
                case( "Boussinesq" )
                   rhoM_0 = rho00
                case default
                   stop "momentumPredictor: unkown case model."
                end select

                ! velocity and momentum at t(0)
                uM_0 = var0(i,j,k,2)
                momM_0 = rhoM_0 * uM_0
                
                ! update momentum
                momM =   rk(1,m) * momM_0 &
                     & + rk(2,m) * momM_1 &
                     & + rk(3,m) * dt*F

                ! calc u(m,*)
                uAst = momM / rhoM

                ! uAst -> var
                var(i,j,k,2) = uAst

             case default
                stop "thetaUpdate: unknown case timeSchemeType"
             end select


          end do
       end do
    end do

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
       stop "momentumPredictor: unknown case yBoundary."
    end select
    
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

             !--- pressure gradient term -> piGrad
             piF = var(i,j+1,k,5)
             piB = var(i,j,  k,5)
             piGrad = kappaInv*MaInv2 * pStrat(k) * (piF-piB)/dy

             !---- volume forces
             volForce = 0.5*( force(i,j,k,2) + force(i,j+1,k,2) )
             
             !--------------------
             !   F(phi) = RHS
             !--------------------
             ! fluxDiff -> convective and viscous fluxes
             ! piGrad   -> pressure gradient along x
             ! volForce -> Gravity, Coriolis
             F = -fluxDiff - piGrad + volForce

             if( correctDivError ) F = F + source(i,j,k,3)

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
                stop "momentumPredictor: unkown case model."
             end select


             ! velocity and momentum at t(m-1)
             vM_1 = var(i,j,k,3)
             momM_1 = rhoM_1 * vM_1


             select case( timeSchemeType ) 

             case( "lowStorage" ) 


                ! q(m-1) -> q(m)
                q(i,j,k,2) = dt * F + alpha(m) * q(i,j,k,2)

                ! rhoV(m-1) -> rhoV(m)
                momM = momM_1 + beta(m) * q(i,j,k,2)

                ! calc v(m,*)
                vAst = momM / rhoM

                ! vAst -> var
                var(i,j,k,3) = vAst

             case( "classical" )

                ! interpolated density at t(0)
                select case( model ) 

                case( "pseudo_incompressible" )
                   rhoM_0   = 0.5 * ( var0(i,j,k,1) + var0(i,j+1,k,1) )
                   if( fluctuationMode ) rhoM_0 = rhoM_0 + rhoStrat(k)
                case( "Boussinesq" )
                   rhoM_0 = rho00
                case default
                   stop "momentumPredictor: unkown case model."
                end select
                
                ! velocity and momentum at t(0)
                vM_0 = var0(i,j,k,3)
                momM_0 = rhoM_0 * vM_0

                ! update momentum
                momM =   rk(1,m) * momM_0 &
                     & + rk(2,m) * momM_1 &
                     & + rk(3,m) * dt*F

                ! calc u(m,*)
                vAst = momM / rhoM

                ! uAst -> var
                var(i,j,k,3) = vAst

!                var(i,j,k,3) = rk(1,m) * var0(i,j,k,3) &
!                     &       + rk(2,m) * var (i,j,k,3) &
!                     &       + rk(3,m) * dt*F


             case default
                stop "Update: unknown case timeSchemeType"
             end select



          end do
       end do
    end do

             
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
       stop "momentumPredictor: unknown case zBoundary."
    end select

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

             !--- pressure gradient term -> piGrad
             piU = var(i,j,k+1,5)
             piD = var(i,j,  k,5)
             piGrad = 0.5 * kappaInv*MaInv2* &
                  & ( Pstrat(k)+Pstrat(k+1) ) * (piU-piD)/dz

             !---- volume forces
             volForce = 0.5*( force(i,j,k,3) + force(i,j,k+1,3) )

             !--------------------
             !   F(phi) = RHS
             !--------------------
             ! fluxDiff -> convective and viscous fluxes
             ! piGrad   -> pressure gradient along x
             ! volForce -> Gravity, Coriolis
             F = -fluxDiff - piGrad + volForce

             if( correctDivError ) F = F + source(i,j,k,4)

             ! interpolated densities
             select case( model ) 
                
             case( "pseudo_incompressible" )

                rhoM_1 = 0.5*( rhoOld(i,j,k) + rhoOld(i,j,k+1) )   ! rho(m-1)
                rhoM   = 0.5*( var(i,j,k,1) + var(i,j,k+1,1) )     ! rho(m)

                if( fluctuationMode ) then
                   rhoM_1 = rhoM_1 + rhoStratTilde(k)
                   rhoM   = rhoM   + rhoStratTilde(k)
                end if
                
             case( "Boussinesq" )
                rhoM_1 = rho00
                rhoM = rho00
             case default
                stop "momentumPredictor: unkown case model."
             end select
             
             
             ! velocity and momentum at t(m-1)
             wM_1 = var(i,j,k,4)
             momM_1 = rhoM_1 * wM_1

             select case( timeSchemeType ) 

             case( "lowStorage" ) 

                ! q(m-1) -> q(m)
                q(i,j,k,3) = dt * F + alpha(m) * q(i,j,k,3)

                ! rhoW(m-1) -> rhoW(m)
                momM = momM_1 + beta(m) * q(i,j,k,3)

                ! calc w(m,*)
                wAst = momM / rhoM

                ! wAst -> var
                var(i,j,k,4) = wAst

             case( "classical" )

                if( model == "pseudo_incompressible" ) then
                   print*,"WARNING: classical RK not implemented for pseud_inc."
                end if

                var(i,j,k,4) = rk(1,m) * var0(i,j,k,4) &
                     &       + rk(2,m) * var (i,j,k,4) &
                     &       + rk(3,m) * dt*F

             case default
                stop "thetaUpdate: unknown case timeSchemeType"
             end select




          end do
       end do
    end do


    if (testCase == "predictorTest") then  
       print*,""
       print*,"update.f90/predict: flux at left and right bound should be equal."
       print*,"fRhoU(-1) = ", flux(-1,1,1,1,2)
       print*,"fRhoU(nx) = ", flux(nx,1,1,1,2)
       print*,""
    end if
             
!   achatzb
!   -------------------------------------
!   in case of topography, 
!   set all velocities normal to the topographic surface to zero,
!   set density in land cells to background density
!   ------------------------------------

    i00=is+nbx-1
    j00=js+nby-1

    if(topography) then
       do k = 0, nz+1
          do j = 0, ny+1
             do i = 0, nx+1
!               u at x interfaces
                if(&
                   topography_mask(i00+i,j00+j,k)&
                   .or.&
                   topography_mask(i00+i+1,j00+j,k)&
                ) then
                   var(i,j,k,2)=0.
                end if

!               v at y interfaces
                if(&
                   topography_mask(i00+i,j00+j,k)&
                   .or.&
                   topography_mask(i00+i,j00+j+1,k)&
                ) then
                   var(i,j,k,3)=0.
                end if

!               w at z interfaces
                if(&
                   topography_mask(i00+i,j00+j,k)&
                   .or.&
                   topography_mask(i00+i,j00+j,k+1)&
                ) then
                   var(i,j,k,4)=0.
                end if

!               density in land cells
                if(topography_mask(i0+i,j0+j,k)) then
                   if( fluctuationMode ) then
                      var(i,j,k,1) = 0.0
                     else
                      var(i,j,k,1) = rhoStrat(k) 
                   end if
                end if

             end do
          end do
       end do
    end if
!   -------------------------------------
!   achatze

  end subroutine momentumPredictor


!---------------------------------------------------------------------------


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
                stop "thetaUpdate: unknown case timeSchemeType"
             end select



          end do
       end do
    end do

!    if(verbose) print*,"update.f90/thetaUpdate: theta(m=",m,") calculated."   ! modified by Junhong Wei (20170216)

    if(verbose .and. master) print*,"update.f90/thetaUpdate: theta(m=",m,") calculated."   ! modified by Junhong Wei (20170216)

  end subroutine thetaUpdate


!---------------------------------------------------------------------------


  subroutine massUpdate (var,var0,flux,source,dt,q,m)
    !-----------------------------
    ! adds mass flux to cell mass
    !-----------------------------
    
    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var0
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux) 
    ! dir = 1..3 > f-, g- and h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW

    ! source terms
    ! 1) divergence error source terms
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
       

    ! init q
    if (m == 1) q = 0.

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

             ! diffusive part
             ! diff = ....

             ! F(phi)
             F = -fluxDiff

             ! subtract divergence error
             if( correctDivError ) then
                
                F = F + source(i,j,k,1)
                
             end if
             
             select case( timeSchemeType ) 
                
             case( "lowStorage" ) 

                ! update: q(m-1) -> q(m)
                q(i,j,k) = dt*F + alpha(m) * q(i,j,k)

                ! update density
                var(i,j,k,1) = var(i,j,k,1) + beta(m) * q(i,j,k)

             case( "classical" )

                var(i,j,k,1) = rk(1,m) * var0(i,j,k,1) &
                     &       + rk(2,m) * var (i,j,k,1) &
                     &       + rk(3,m) * dt*F

             case default
                stop "massUpdate: unknown case timeSchemeType"
             end select

          end do
       end do
    end do

!    if(verbose) print*,"update.f90/massUpdate: rho(m=",m,") calculated."   ! modified by Junhong Wei (20170216)

    if(verbose .and. master) print*,"update.f90/massUpdate: rho(m=",m,") calculated."   ! modified by Junhong Wei (20170216)

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
         & intent(in) :: var0
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), intent(in) :: flux
    ! flux(i,j,k,dir,iFlux) 
    ! dir = 1..3 > f-, g- and h-flux in x,y,z-direction
    ! iFlux = 8..10 > Rho_nIce, Rho_qIce, Rho_qv

    ! source terms
    ! 1) divergence error source terms
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: source
    
    real, intent(in) :: dt
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,3), &
         & intent(inout) :: q

    integer, intent(in) :: m
    
    ! local variables
    integer :: i,j,k,l
    real, dimension(3)    :: fL,fR        ! flux Left/Right
    real, dimension(3)    :: gB,gF        ! flux Backward/Forward
    real, dimension(3)    :: hD,hU        ! flux Downward/Upward
    real, dimension(3)    :: fluxDiff         ! convective part
    real, dimension(3)    :: F            ! F(phi)
    

    ! init q
    if (m == 1) q = 0.

    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
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

             ! subtract divergence error
             if( correctDivError ) then
                
                F(:) = F(:) + source(i,j,k,nVar-3:nVar)
                
             end if
             
             select case( timeSchemeType ) 
                
             case( "lowStorage" ) 

                ! update: q(m-1) -> q(m)
                q(i,j,k,:) = dt*F(:) + alpha(m) * q(i,j,k,:)

                ! update variables
                var(i,j,k,nVar-3:nVar) = var(i,j,k,nVar-3:nVar) + beta(m) * q(i,j,k,1:4)

             case( "classical" )

                var(i,j,k,nVar-3:nVar) = rk(1,m) * var0(i,j,k,nVar-3:nVar) &
                     &       + rk(2,m) * var (i,j,k,nVar-3:nVar) &
                     &       + rk(3,m) * dt*F(1:4)

             case default
                stop "iceUpdate: unknown case timeSchemeType"
             end select

          end do
       end do
    end do

!    if(verbose) print*,"update.f90/massUpdate: rho(m=",m,") calculated."   ! modified by Junhong Wei (20170216)

    if(verbose .and. master) print*,"update.f90/iceUpdate: ice(m=",m,") calculated."   ! modified by Junhong Wei (20170216)

  end subroutine iceUpdate


!---------------------------------------------------------------------------

  subroutine timestep (var,ray,dt,errFlag)
    !---------------------------------------------
    ! compute time step from stability criteria:
    ! 1) CFL criterion for advection
    ! 2) von Neumann cirterion for dissipation
    ! 3) set maximum time step
    ! 4) buouyancy acceleration ->  1/2 b_max*dt^2 < dz
    ! 5) gravity wave phase condition: c_phi_x * dt < dx, ...dy,dz
    !---------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),intent(in) :: var
    type(rayType), dimension(nRay), intent(in)                        :: ray
    real, intent(out) :: dt
    logical, intent(out) :: errFlag

    ! locals
    real :: uMax, vMax, wMax
    real :: dtConv, dtVisc, dtCond
    real :: dtMax

    ! ray tracer velocities
    real :: cgxMax, cgyMax, cgzMax
    real :: dtRay

    ! Buoyancy time step restriction
    real              :: dtBuoy
    real,dimension(3) :: bMax, bMaxNew, duMax

    ! Gravity wave time stop restriction
    real :: lambdaMax                     ! max GW length to be time resolved
    real :: dtWave, lambdaX, lambdaZ, kk, mm, kMin, cX, cZ
    
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

       if( master ) then   ! modified by Junhong Wei (20170216)
       write(*,fmt="(a25,es15.1,a8)") "dt = dtFix = ", dt*tRef, "seconds"
       end if   ! modified by Junhong Wei (20170216)

    else

    !-------------------------------------------
    !           Variable time step
    !-------------------------------------------

       select case( model ) 

       case( "WKB" ) 

          !----------------------------
          !    Ray tracer time step
          !----------------------------

          ! CFL time step restriction for 
          ! Eulerian wave action transport
          
          uMax = maxval( abs(var(1:nx,1:ny,1:nz,2) )) + small
          vMax = maxval( abs(var(1:nx,1:ny,1:nz,3) )) + small
          wMax = maxval( abs(var(1:nx,1:ny,1:nz,4) )) + small
          
          cgxMax = abs(maxval(ray(:)%cgx)) + abs(uMax)
          cgyMax = abs(maxval(ray(:)%cgy)) + abs(vMax)
          cgzMax = abs(maxval(ray(:)%cgz)) + abs(wMax)

          dtRay = cfl * min(dx/cgxMax, dy/cgyMax, dz/cgzMax)


          !----------------------------
          !    Maximal time step
          !----------------------------

          dtMax = dtMax_dim / tRef

          
          !-------------------------------
          !        Make your choice
          !-------------------------------

          dt = min(dtMax,dtRay)
          
          if( dt == dtMax ) then
             write(*,fmt="(a25,es15.1,a8)") "dt = dtMax = ", dt*tRef, "seconds"
          else if( dt == dtRay ) then 
             write(*,fmt="(a25,es15.1,a8)") "dt = dtRay = ", dt*tRef, "seconds"
          end if

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

          dtConv = cfl * min(dx/uMax, dy/vMax, dz/wMax)

!         testb
!         print*,"uMax and cfl = ",uMax*uRef, cfl*dx/uMax*tRef
!         
!         do k=1,nz
!            if(uMax == (maxval( abs(var(1:nx,1:ny,k,2) )) + small)) then
!               print*,"reached at k = ",k
!            end if
!         end do

!         print*,"vMax and cfl = ",vMax*uRef, cfl*dy/vMax*tRef
!         
!         do k=1,nz
!            if(vMax == (maxval( abs(var(1:nx,1:ny,k,3) )) + small)) then
!               print*,"reached at k = ",k
!            end if
!         end do

!         print*,"wMax and cfl = ",wMax*uRef, cfl*dz/wMax*tRef
!         
!         do k=1,nz
!            if(wMax == (maxval( abs(var(1:nx,1:ny,k,4) )) + small)) then
!               print*,"reached at k = ",k
!            end if
!         end do
!         teste


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
                      stop "timeStep: unknown case model."
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
             
             dtBuoy = max(&
                  & -uMax/bMax(1)+sqrt((uMax/bMax(1))**2+2.*cfl*dx/bMax(1)),&
                  & -vMax/bMax(2)+sqrt((vMax/bMax(2))**2+2.*cfl*dy/bMax(2)),&
                  & -wMax/bMax(3)+sqrt((wMax/bMax(3))**2+2.*cfl*dz/bMax(3)))
             
!xxxx debug             
             if( dtBuoy*tRef < 1.e-2 ) then

                print*,"dtBuoy*tRef  = ", dtBuoy*tRef
                print*,"bMax(3) = ", bMax(3)*FrInv2
                
             end if
!xxxx end debug
             
          else
             
             dtBuoy = 1.0e20/tRef   ! set to high value if not needed

          end if
          
!old:          dtBuoy = sqrt(2.0*dz/bMax)


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
          !
          
          dtWave = 1./(NN+small)
          
!!$          ! lambdaMax = lambdaMax_dim/lRef
!!$          
!!$          if( testCase == "wavePacket" ) then
!!$             
!!$             lambdaX = lambdaX_dim / lRef
!!$             lambdaZ = lambdaZ_dim / lRef
!!$             kk = abs(2.*pi/lambdaX) + small
!!$             mm = abs(2.*pi/lambdaZ) + small
!!$             kMin = min(kk,mm)
!!$             
!!$             cX = NN/kMin/sqrt(2.0)
!!$             cZ = cX * kk/mm
!!$
!!$             dtWave = cfl_wave * min(dx/cX, dz/cZ)
!!$             
!!$!             dtWave = cfl_wave * 2.0*pi / (NN+small) * dtWave
!!$             
!!$          end if
          

          !-------------------------------
          !        Make your choice
          !-------------------------------

          dt = min(dtVisc,dtCond,dtConv,dtMax,dtBuoy,dtWave)


          !-----------------------------------------
          !     Inform on time step restrictions
          !-----------------------------------------

       if( master ) then   ! modified by Junhong Wei (20170216)

          write(*,fmt="(a25,es15.1,a8)") "dtVisc = ", dtVisc*tRef, "seconds"
          write(*,fmt="(a25,es15.1,a8)") "dtCond = ", dtCond*tRef, "seconds"
          write(*,fmt="(a25,es15.1,a8)") "dtConv = ", dtConv*tRef, "seconds"
          write(*,fmt="(a25,es15.1,a8)") "dtMax = ", dtMax*tRef, "seconds"
          write(*,fmt="(a25,es15.1,a8)") "dtBuoy = ", dtBuoy*tRef, "seconds"
          write(*,fmt="(a25,es15.1,a8)") "dtWave = ", dtWave*tRef, "seconds"
          print*,""

          if(dt == dtMax) then
             write(*,fmt="(a25,es15.1,a8)") "--> dt = dtMax = ", dt*tRef, "seconds"
          else if(dt == dtConv) then
             write(*,fmt="(a25,es15.1,a8)") "--> dt = dtConv = ", dt*tRef, "seconds"
          else if(dt == dtVisc) then
             write(*,fmt="(a25,es15.1,a8)") "--> dt = dtVisc = ", dt*tRef, "seconds"
          else if(dt == dtCond) then
             write(*,fmt="(a25,es15.1,a8)") "--> dt = dtCond = ", dt*tRef, "seconds"
          else if(dt == dtBuoy) then
             write(*,fmt="(a25,es15.1,a8)") "--> dt = dtBuoy = ", dt*tRef, "seconds"
          else if(dt == dtWave) then
             write(*,fmt="(a25,es15.1,a8)") "--> dt = dtWave = ", dt*tRef, "seconds"
          else
             write(*,fmt="(a25,es15.1,a8)") "--> dt = ????? = ", dt*tRef, "seconds"
          end if
          print*,""

       end if   ! modified by Junhong Wei (20170216)

       case default
          stop "timestep: unknown case model."
       end select   ! WKB / full model

    end if

    ! error handling for too small time steps
    if( dt*tRef < dtMin_dim ) errFlag = .true. 

  end subroutine timestep


!---------------------------------------------------------------------------


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
! modified by Junhong Wei (20160726)  (starting line)

  subroutine CoefDySma_update(var)
    !--------------------------------------
    ! calculate the Coefficient for Dynamic Smagorinsky Scheme
    !--------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var

    ! more variables
    real :: delta_DySma
    real :: uL_DySma, uR_DySma, uB_DySma, uF_DySma, uD_DySma, uU_DySma
    real :: vL_DySma, vR_DySma, vB_DySma, vF_DySma, vD_DySma, vU_DySma
    real :: wL_DySma, wR_DySma, wB_DySma, wF_DySma, wD_DySma, wU_DySma
    real :: du_dx_DySma, du_dy_DySma, du_dz_DySma
    real :: dv_dx_DySma, dv_dy_DySma, dv_dz_DySma
    real :: dw_dx_DySma, dw_dy_DySma, dw_dz_DySma

    ! allocatable fields
    real, dimension(:,:,:,:,:), allocatable :: Sij_DySma, Lij_DySma, Mij_DySma
    real, dimension(:,:,:), allocatable :: SSS_DySma

    real, dimension(:,:,:,:,:), allocatable :: uiuj_smth_DySma, SSS_Sij_smth_DySma, Sij_smth_DySma
    real, dimension(:,:,:,:), allocatable :: ui_smth_DySma
    real, dimension(:,:,:), allocatable :: SSS_smth_DySma
    real, dimension(:,:,:), allocatable :: LijMij_smth_DySma, MijMij_smth_DySma
    real, dimension(:,:,:), allocatable :: CS2_DySma

    integer :: allocstat
    integer :: i,j,k
    integer :: iii_no,jjj_no

    !UA 2 corresponds to averaging over 5 grid cells. Why so many? Not better 1,
    !so that you average over 3?

    ! replied by JW (20160824): It has been revised. 

    integer :: smth_npts1_DySma, smth_npts2_DySma
!    parameter (smth_npts1_DySma = 2)  ! revised by JW (20160824)
    parameter (smth_npts1_DySma = 1)   ! revised by JW (20160824)
    parameter (smth_npts2_DySma = 2)  ! revised by JW (20160824)
!    parameter (smth_npts2_DySma = 5)   ! revised by JW (20160824)

    ! Allocate local fields
    allocate(Sij_DySma(1:nx,1:ny,1:nz,1:3,1:3), stat=allocstat); if(allocstat/=0) &
         & stop "CoefDySma_update:alloc failed"
    allocate(Lij_DySma(1:nx,1:ny,1:nz,1:3,1:3), stat=allocstat); if(allocstat/=0) &
         & stop "CoefDySma_update:alloc failed"
    allocate(Mij_DySma(1:nx,1:ny,1:nz,1:3,1:3), stat=allocstat); if(allocstat/=0) &
         & stop "CoefDySma_update:alloc failed"
    allocate(SSS_DySma(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "CoefDySma_update:alloc failed"
    allocate(uiuj_smth_DySma(1:nx,1:ny,1:nz,1:3,1:3), stat=allocstat); if(allocstat/=0) &
         & stop "CoefDySma_update:alloc failed"
    allocate(SSS_Sij_smth_DySma(1:nx,1:ny,1:nz,1:3,1:3), stat=allocstat); if(allocstat/=0) &
         & stop "CoefDySma_update:alloc failed"
    allocate(Sij_smth_DySma(1:nx,1:ny,1:nz,1:3,1:3), stat=allocstat); if(allocstat/=0) &
         & stop "CoefDySma_update:alloc failed"
    allocate(ui_smth_DySma(1:nx,1:ny,1:nz,1:3), stat=allocstat); if(allocstat/=0) &
         & stop "CoefDySma_update:alloc failed"
    allocate(SSS_smth_DySma(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "CoefDySma_update:alloc failed"
    allocate(LijMij_smth_DySma(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "CoefDySma_update:alloc failed"
    allocate(MijMij_smth_DySma(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "CoefDySma_update:alloc failed"
    allocate(CS2_DySma(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "CoefDySma_update:alloc failed"


    ! calculate delta

  if(ny.eq.1)then
    delta_DySma = (dx*dz)**(1.0/2.0)          ! 2D problems
  else
    delta_DySma = (dx*dy*dz)**(1.0/3.0)       ! 3D problems
  end if


    ! calculate S_ij

       !---------------------------------
       !         Loop over field
       !---------------------------------

       do k = 1,nz
          do j = 1,ny
             do i = 1,nx

             uL_DySma = var(i-1,j,k,2)
             uR_DySma = var(i,j,k,2)
             !UA not sure whether this averaging is necessary. Compare, e.g. the
             !viscous fluxes. There it is not done.

             ! replied by JW (20160824): For now, there is no revision on this part, since it is a relatively minor issue at the moment. 
             uB_DySma = 0.5 * ( var(i-1,j-1,k,2) + var(i,j-1,k,2) )
             uF_DySma = 0.5 * ( var(i-1,j+1,k,2) + var(i,j+1,k,2) )
             uD_DySma = 0.5 * ( var(i-1,j,k-1,2) + var(i,j,k-1,2) )
             uU_DySma = 0.5 * ( var(i-1,j,k+1,2) + var(i,j,k+1,2) )

             vL_DySma = 0.5 * ( var(i-1,j,k,3) + var(i-1,j-1,k,3) )
             vR_DySma = 0.5 * ( var(i+1,j,k,3) + var(i+1,j-1,k,3) )
             vB_DySma = var(i,j-1,k,3)
             vF_DySma = var(i,j,k,3)
             vD_DySma = 0.5 * ( var(i,j,k-1,3) + var(i,j-1,k-1,3) )
             vU_DySma = 0.5 * ( var(i,j,k+1,3) + var(i,j-1,k+1,3) )

             wL_DySma = 0.5 * ( var(i-1,j,k-1,4) + var(i-1,j,k,4) )
             wR_DySma = 0.5 * ( var(i+1,j,k-1,4) + var(i+1,j,k,4) )
             wB_DySma = 0.5 * ( var(i,j-1,k-1,4) + var(i,j-1,k,4) )
             wF_DySma = 0.5 * ( var(i,j+1,k-1,4) + var(i,j+1,k,4) )
             wD_DySma = var(i,j,k-1,4)
             wU_DySma = var(i,j,k,4)


             du_dx_DySma = ( uR_DySma - uL_DySma ) / dx
             !UA in case without averaging, no factor 2!

             ! replied by JW (20160824): For now, there is no revision on this part, since it is a relatively minor issue at the moment. 
             du_dy_DySma = ( uF_DySma - uB_DySma ) / ( 2.0 * dy )
             du_dz_DySma = ( uU_DySma - uD_DySma ) / ( 2.0 * dz )

             dv_dx_DySma = ( vR_DySma - vL_DySma ) / ( 2.0 * dx )
             dv_dy_DySma = ( vF_DySma - vB_DySma ) / dy
             dv_dz_DySma = ( vU_DySma - vD_DySma ) / ( 2.0 * dz )

             dw_dx_DySma = ( wR_DySma - wL_DySma ) / ( 2.0 * dx )
             dw_dy_DySma = ( wF_DySma - wB_DySma ) / ( 2.0 * dy )
             dw_dz_DySma = ( wU_DySma - wD_DySma ) / dz

       !UA replace 0.25 by 0.5 below (typo in my text, sorry)

       ! replied by JW (20160824): It has been revised.

!             Sij_DySma(i,j,k,1,1) = 0.25 * ( du_dx_DySma + du_dx_DySma )  ! revised by JW (20160824)
!             Sij_DySma(i,j,k,1,2) = 0.25 * ( du_dy_DySma + dv_dx_DySma )  ! revised by JW (20160824)
!             Sij_DySma(i,j,k,1,3) = 0.25 * ( du_dz_DySma + dw_dx_DySma )  ! revised by JW (20160824)
!             Sij_DySma(i,j,k,2,1) = 0.25 * ( dv_dx_DySma + du_dy_DySma )  ! revised by JW (20160824)
!             Sij_DySma(i,j,k,2,2) = 0.25 * ( dv_dy_DySma + dv_dy_DySma )  ! revised by JW (20160824)
!             Sij_DySma(i,j,k,2,3) = 0.25 * ( dv_dz_DySma + dw_dy_DySma )  ! revised by JW (20160824)
!             Sij_DySma(i,j,k,3,1) = 0.25 * ( dw_dx_DySma + du_dz_DySma )  ! revised by JW (20160824)
!             Sij_DySma(i,j,k,3,2) = 0.25 * ( dw_dy_DySma + dv_dz_DySma )  ! revised by JW (20160824)
!             Sij_DySma(i,j,k,3,3) = 0.25 * ( dw_dz_DySma + dw_dz_DySma )  ! revised by JW (20160824)


             Sij_DySma(i,j,k,1,1) = 0.5 * ( du_dx_DySma + du_dx_DySma )  ! revised by JW (20160824)
             Sij_DySma(i,j,k,1,2) = 0.5 * ( du_dy_DySma + dv_dx_DySma )  ! revised by JW (20160824)
             Sij_DySma(i,j,k,1,3) = 0.5 * ( du_dz_DySma + dw_dx_DySma )  ! revised by JW (20160824)
             Sij_DySma(i,j,k,2,1) = 0.5 * ( dv_dx_DySma + du_dy_DySma )  ! revised by JW (20160824)
             Sij_DySma(i,j,k,2,2) = 0.5 * ( dv_dy_DySma + dv_dy_DySma )  ! revised by JW (20160824)
             Sij_DySma(i,j,k,2,3) = 0.5 * ( dv_dz_DySma + dw_dy_DySma )  ! revised by JW (20160824)
             Sij_DySma(i,j,k,3,1) = 0.5 * ( dw_dx_DySma + du_dz_DySma )  ! revised by JW (20160824)
             Sij_DySma(i,j,k,3,2) = 0.5 * ( dw_dy_DySma + dv_dz_DySma )  ! revised by JW (20160824)
             Sij_DySma(i,j,k,3,3) = 0.5 * ( dw_dz_DySma + dw_dz_DySma )  ! revised by JW (20160824)


               SSS_DySma(i,j,k) = 0.0
               do jjj_no = 1,3
               do iii_no = 1,3
               SSS_DySma(i,j,k) = SSS_DySma(i,j,k) + ( Sij_DySma(i,j,k,iii_no,jjj_no) * Sij_DySma(i,j,k,iii_no,jjj_no) )
               end do
               end do
               SSS_DySma(i,j,k) = SSS_DySma(i,j,k) * 2.0
               SSS_DySma(i,j,k) = sqrt( SSS_DySma(i,j,k) )


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

        ui_smth_DySma(i,j,k,1) = 0.5 * ( var(i,j,k,2) + var(i-1,j,k,2) )
        ui_smth_DySma(i,j,k,2) = 0.5 * ( var(i,j,k,3) + var(i,j-1,k,3) )
        ui_smth_DySma(i,j,k,3) = 0.5 * ( var(i,j,k,4) + var(i,j,k-1,4) )

        SSS_smth_DySma(i,j,k)  =  SSS_DySma(i,j,k)

        do jjj_no = 1,3
           do iii_no = 1,3

           uiuj_smth_DySma(i,j,k,iii_no,jjj_no) = ui_smth_DySma(i,j,k,iii_no) * ui_smth_DySma(i,j,k,jjj_no)
           Sij_smth_DySma(i,j,k,iii_no,jjj_no)  = Sij_DySma(i,j,k,iii_no,jjj_no)
           SSS_Sij_smth_DySma(i,j,k,iii_no,jjj_no) = SSS_DySma(i,j,k) * Sij_DySma(i,j,k,iii_no,jjj_no)

           end do
        end do

        end do
     end do
  end do

       !---------------------------------
       !         Smoothing at the 1st step
       !---------------------------------

       !UA why only smoothening in x? You should also smoothen in z (and in y if
       !3D)

       ! replied by JW (20160824): The "Var3DSmthDySma_X" has been revised. For now, the smoothing is done over the x direction first, and then it is done again over the z direction. Please check. For the testing stage, I would only focus on the 2.5D simulation. 
       ! Note by Junhong Wei (20160916): The name of "Var3DSmthDySma_X" is changed into "Var3DSmthDySma". The new subroutine includes all kinds of smoothing. I have modified all the related lines for the new subroutine "Var3DSmthDySma". The changes are not marked. 

  if(ny.eq.1)then

  do iii_no = 1,3
    call Var3DSmthDySma( ui_smth_DySma(1:nx,1:ny,1:nz,iii_no), smth_npts1_DySma, "XZ_local_smth" )
  end do

    call Var3DSmthDySma( SSS_smth_DySma(1:nx,1:ny,1:nz),       smth_npts1_DySma, "XZ_local_smth" )

  do jjj_no = 1,3
     do iii_no = 1,3

     call Var3DSmthDySma( uiuj_smth_DySma(1:nx,1:ny,1:nz,iii_no,jjj_no), smth_npts1_DySma, "XZ_local_smth" )
     call Var3DSmthDySma( Sij_smth_DySma(1:nx,1:ny,1:nz,iii_no,jjj_no),  smth_npts1_DySma, "XZ_local_smth" )
     call Var3DSmthDySma( SSS_Sij_smth_DySma(1:nx,1:ny,1:nz,iii_no,jjj_no), smth_npts1_DySma, "XZ_local_smth" )

     end do
  end do

  else

  do iii_no = 1,3
    call Var3DSmthDySma( ui_smth_DySma(1:nx,1:ny,1:nz,iii_no), smth_npts1_DySma, "XYZ_local_smth" )
  end do

    call Var3DSmthDySma( SSS_smth_DySma(1:nx,1:ny,1:nz),       smth_npts1_DySma, "XYZ_local_smth" )

  do jjj_no = 1,3
     do iii_no = 1,3

     call Var3DSmthDySma( uiuj_smth_DySma(1:nx,1:ny,1:nz,iii_no,jjj_no), smth_npts1_DySma, "XYZ_local_smth" )
     call Var3DSmthDySma( Sij_smth_DySma(1:nx,1:ny,1:nz,iii_no,jjj_no),  smth_npts1_DySma, "XYZ_local_smth" )
     call Var3DSmthDySma( SSS_Sij_smth_DySma(1:nx,1:ny,1:nz,iii_no,jjj_no), smth_npts1_DySma, "XYZ_local_smth" )

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

        LijMij_smth_DySma(i,j,k) = 0.0
        MijMij_smth_DySma(i,j,k) = 0.0

        do jjj_no = 1,3
           do iii_no = 1,3

        Lij_DySma(i,j,k,iii_no,jjj_no) = uiuj_smth_DySma(i,j,k,iii_no,jjj_no) &
          &  - ( ui_smth_DySma(i,j,k,iii_no)*ui_smth_DySma(i,j,k,jjj_no) )

        !UA below is not right. Should be delta_DySma**2 * ... - (2.0*smth_npts1_DySma+1.0)**2 * ...

        ! replied by JW (20160824): It has been revised.

!        Mij_DySma(i,j,k,iii_no,jjj_no) = SSS_Sij_smth_DySma(i,j,k,iii_no,jjj_no) - ( SSS_smth_DySma(i,j,k)*Sij_smth_DySma(i,j,k,iii_no,jjj_no) )  ! revised by JW (20160824)
!        Mij_DySma(i,j,k,iii_no,jjj_no) = Mij_DySma(i,j,k,iii_no,jjj_no) * (delta_DySma**2)                                                        ! revised by JW (20160824)

        ! revised by JW (20160824)
        Mij_DySma(i,j,k,iii_no,jjj_no) = SSS_Sij_smth_DySma(i,j,k,iii_no,jjj_no) &
          &  - ( ( (2.0*smth_npts1_DySma+1.0)**2 )*SSS_smth_DySma(i,j,k)*        &
          &      Sij_smth_DySma(i,j,k,iii_no,jjj_no) )
        Mij_DySma(i,j,k,iii_no,jjj_no) = Mij_DySma(i,j,k,iii_no,jjj_no) * (delta_DySma**2)


        LijMij_smth_DySma(i,j,k) = LijMij_smth_DySma(i,j,k) + ( Lij_DySma(i,j,k,iii_no,jjj_no)*Mij_DySma(i,j,k,iii_no,jjj_no) )
        MijMij_smth_DySma(i,j,k) = MijMij_smth_DySma(i,j,k) + ( Mij_DySma(i,j,k,iii_no,jjj_no)*Mij_DySma(i,j,k,iii_no,jjj_no) )

           end do
        end do

        end do
     end do
  end do

    !UA Again, also smoothen in y and z. I would use smth_npts2_DySma= 5, e.g. 

        ! replied by JW (20160824): It has been revised.

!    call Var3DSmthDySma( LijMij_smth_DySma(1:nx,1:ny,1:nz),       smth_npts2_DySma, "X_whole_smth" )
!    call Var3DSmthDySma( MijMij_smth_DySma(1:nx,1:ny,1:nz),       smth_npts2_DySma, "X_whole_smth" )

  if(ny.eq.1)then

    call Var3DSmthDySma( LijMij_smth_DySma(1:nx,1:ny,1:nz),       smth_npts2_DySma, "XZ_local_smth" )
    call Var3DSmthDySma( MijMij_smth_DySma(1:nx,1:ny,1:nz),       smth_npts2_DySma, "XZ_local_smth" )

  else

    call Var3DSmthDySma( LijMij_smth_DySma(1:nx,1:ny,1:nz),       smth_npts2_DySma, "XYZ_local_smth" )
    call Var3DSmthDySma( MijMij_smth_DySma(1:nx,1:ny,1:nz),       smth_npts2_DySma, "XYZ_local_smth" )

  end if

       !---------------------------------
       !         Get the final results
       !---------------------------------

  do k = 1,nz
     do j = 1,ny
        do i = 1,nx

           if(MijMij_smth_DySma(i,j,k) /= 0.) then
              CS2_DySma(i,j,k) &
              = 0.5 * LijMij_smth_DySma(i,j,k) / MijMij_smth_DySma(i,j,k)
             else
             CS2_DySma(i,j,k)=0.
           end if

           ! *** modified by Junhong Wei (20160918) *** (starting line)
           if( CS2_DySma(i,j,k) < 0.0 )then

              CS2_DySma(i,j,k) = 0.0

           end if
           ! *** modified by Junhong Wei (20160918) *** (finishing line)

        var(i,j,k,7) = CS2_DySma(i,j,k) * (delta_DySma**2) * SSS_DySma(i,j,k)

        end do
     end do
  end do


  ! *** modified by Junhong Wei (20160918) *** (starting line)

  ! *** set values for the ghost celss ***

    call setHaloAndBoundary( var(:,:,:,7), nbx, nby, nbz )

  ! *** modified by Junhong Wei (20160918) *** (finishing line)



    ! deallocate local fields
    deallocate(Sij_DySma, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(Lij_DySma, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(Mij_DySma, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(SSS_DySma, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(uiuj_smth_DySma, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(SSS_Sij_smth_DySma, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(Sij_smth_DySma, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(ui_smth_DySma, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(SSS_smth_DySma, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(LijMij_smth_DySma, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(MijMij_smth_DySma, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(CS2_DySma, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"


          return

  end subroutine CoefDySma_update



  subroutine Var3DSmthDySma(var3D_DySma,nsmth_DySma,homog_dir_DySma)
    !--------------------------------------
    ! calculate the Coefficient for Dynamic Smagorinsky Scheme
    ! Note by Junhong Wei: This subroutine has been modified all over again on 20160915. The changes are not marked, since the old version is erased. 
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
!   achatzb
    integer :: nsmthall,ismth,jsmth,ksmth
    integer :: i0,j0
!   achatze


    allocate(var3D_DySma_Extend( (0-nsmth_DySma) : (nx+nsmth_DySma), &
      &                          (0-nsmth_DySma) : (ny+nsmth_DySma), &
      &                          (0-nsmth_DySma) : (nz+nsmth_DySma) ), stat=allocstat)
    if(allocstat/=0) stop "Var3DSmthDySma:alloc failed"

    ! set the values for var3D_DySma_Extend

    var3D_DySma_Extend(1:nx,1:ny,1:nz) = var3D_DySma(1:nx,1:ny,1:nz)

    call setHaloAndBoundary( var3D_DySma_Extend(:,:,:), nsmth_DySma, nsmth_DySma, nsmth_DySma )

! start to do the smoothing

!   achatzb
    i0=is+nbx-1
    j0=js+nby-1
!   achatze

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

     if(nz /= sizeZ) stop " DYNAMIC SMAGORINSKY NOT READY FOR MPI IN Z"

       do k = 1,nz
!         achatzb correct handling of solid and periodic boundaries in z
!         kmin=max( 1,k-nsmth_DySma)
!         kmax=min(nz,k+nsmth_DySma)
!         nsmthv=kmax-kmin+1
          if(zBoundary == "solid_wall") then
             kmin=max( 1,k-nsmth_DySma)
             kmax=min(nz,k+nsmth_DySma)
            else if (zBoundary == "periodic") then
             kmin=k-nsmth_DySma
             kmax=k+nsmth_DySma
            else
             stop "vertical smoothing: unknown case zBoundary."
          end if

          nsmthv=kmax-kmin+1
!         achatze

!         testb
!         if (master) write(*,*)k,kmin,kmax,nsmthv
!         teste
             
          do j = 1,ny
             do i = 1,nx

!               achatzb 
!               averaging dyn. Smag. coeff. only over atmosphere cells

!               var3D_DySma(i,j,k)&
!               =&
!               sum(&
!               var3D_DySma_Extend( &
!               (i-nsmth_DySma):(i+nsmth_DySma), &
!               (j-nsmth_DySma):(j+nsmth_DySma), &
!               kmin:kmax &
!               )&
!               )&
!               /( ((2*nsmth_DySma + 1)**2) * nsmthv)

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
!              achatze

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

     if(nz /= sizeZ) stop " DYNAMIC SMAGORINSKY NOT READY FOR MPI IN Z"

       do k = 1,nz
!         achatzb correct handling of solid and periodic boundaries in z
!         kmin=max( 1,k-nsmth_DySma)
!         kmax=min(nz,k+nsmth_DySma)
!         nsmthv=kmax-kmin+1
          if(zBoundary == "solid_wall") then
             kmin=max( 1,k-nsmth_DySma)
             kmax=min(nz,k+nsmth_DySma)
            else if (zBoundary == "periodic") then
             kmin=k-nsmth_DySma
             kmax=k+nsmth_DySma
            else
             stop "vertical smoothing: unknown case zBoundary."
          end if

          nsmthv=kmax-kmin+1
!         achatze

!         testb
!         if (master) write(*,*)k,kmin,kmax,nsmthv
!         teste
             
          do j = 1,ny
             do i = 1,nx

!               achatzb 
!               averaging dyn. Smag. coeff. only over atmosphere cells

!               var3D_DySma(i,j,k)&
!               =&
!               sum(&
!               var3D_DySma_Extend( &
!               (i-nsmth_DySma):(i+nsmth_DySma), &
!               j, &
!               kmin:kmax &
!               )&
!               )&
!               /( (2*nsmth_DySma + 1) * nsmthv)

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
!              achatze
            end do
         end do
      end do


    case( "X_whole_smth" )

    if(xBoundary /= "periodic") &
    stop "DYNAMIC SMAGORINSKY NOT READY FOR NON-PERIODIC BOUNDARY &
          &CONDITIONS IN X"


       !---------------------------------
       !         Loop over field
       !---------------------------------

       do k = 1,nz
          do j = 1,ny

!            achatzb 
!            averaging dyn. Smag. coeff. only over atmosphere cells

!            var3D_DySma(:,j,k)&
!            =sum(var3D_DySma_Extend( 1:nx, j, k ))/(nx*1.0)

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
!            achatze
          end do
       end do


    case default
       stop "unknown case homog_dir_DySma."
    end select





    ! deallocate local fields
    deallocate(var3D_DySma_Extend, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"

          return

  end subroutine Var3DSmthDySma






  subroutine setHaloAndBoundary(var3D_HaloBC,nbx_HaloBC,nby_HaloBC,nbz_HaloBC)
    !--------------------------------------
    ! set Halo and Boundary
    !--------------------------------------

    ! in/out variables
    integer, intent(in) :: nbx_HaloBC,nby_HaloBC,nbz_HaloBC

    real, dimension(-nbx_HaloBC:nx+nbx_HaloBC,-nby_HaloBC:ny+nby_HaloBC,-nbz_HaloBC:nz+nbz_HaloBC), &
         & intent(inout) :: var3D_HaloBC

    ! auxiliary fields for "var" with ghost cells (rho)
    real, dimension(nbx_HaloBC,-nby_HaloBC:ny+nby_HaloBC,nz) :: xRhoSliceLeft_send, xRhoSliceRight_send
    real, dimension(nbx_HaloBC,-nby_HaloBC:ny+nby_HaloBC,nz) :: xRhoSliceLeft_recv, xRhoSliceRight_recv

    real, dimension(-nbx_HaloBC:nx+nbx_HaloBC,nby_HaloBC,nz) :: yRhoSliceBack_send, yRhoSliceForw_send
    real, dimension(-nbx_HaloBC:nx+nbx_HaloBC,nby_HaloBC,nz) :: yRhoSliceBack_recv, yRhoSliceForw_recv

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
             xRhoSliceLeft_send (i,-nby_HaloBC:ny+nby_HaloBC,1:nz) = var3D_HaloBC(i,-nby_HaloBC:ny+nby_HaloBC,1:nz)
             xRhoSliceRight_send(i,-nby_HaloBC:ny+nby_HaloBC,1:nz) = var3D_HaloBC(nx-nbx_HaloBC+i,-nby_HaloBC:ny+nby_HaloBC,1:nz)
          end do

             ! left -> right
             source = left
             dest = right
             tag = 100

             i0 = 1; j0 = -nby_HaloBC; k0 = 1

             call mpi_sendrecv(xRhoSliceRight_send(i0,j0,k0), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & xRhoSliceLeft_recv(i0,j0,k0), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_left, ierror)

             ! right -> left
             source = right
             dest = left
             tag = 100

             call mpi_sendrecv(xRhoSliceLeft_send(i0,j0,k0), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & xRhoSliceRight_recv(i0,j0,k0), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_right, ierror)

             ! write auxiliary slice to var field
             do i = 1,nbx_HaloBC

                ! right halos
                var3D_HaloBC(nx+i,-nby_HaloBC:ny+nby_HaloBC,1:nz) = xRhoSliceRight_recv(i,-nby_HaloBC:ny+nby_HaloBC,1:nz)

                ! left halos
                var3D_HaloBC(-nbx_HaloBC+i,-nby_HaloBC:ny+nby_HaloBC,1:nz) = xRhoSliceLeft_recv(i,-nby_HaloBC:ny+nby_HaloBC,1:nz)

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
             yRhoSliceBack_send(-nbx_HaloBC:nx+nbx_HaloBC,j,1:nz) = var3D_HaloBC(-nbx_HaloBC:nx+nbx_HaloBC,j,1:nz)
             yRhoSliceForw_send(-nbx_HaloBC:nx+nbx_HaloBC,j,1:nz) = var3D_HaloBC(-nbx_HaloBC:nx+nbx_HaloBC,ny-nby_HaloBC+j,1:nz)
          end do

             ! back -> forw
             source = back
             dest = forw
             tag = 100
             
             i0 = -nbx_HaloBC; j0 = 1; k0 = 1

             call mpi_sendrecv(yRhoSliceForw_send(i0,j0,k0), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & yRhoSliceBack_recv(i0,j0,k0), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_back, ierror)

             ! forw -> back
             source = forw
             dest = back
             tag = 100

             call mpi_sendrecv(yRhoSliceBack_send(i0,j0,k0), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & yRhoSliceForw_recv(i0,j0,k0), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_forw, ierror)

             ! write auxiliary slice to var field
             do j = 1, nby_HaloBC

                ! right halos
                var3D_HaloBC(-nbx_HaloBC:nx+nbx_HaloBC,ny+j,1:nz) = yRhoSliceForw_recv(-nbx_HaloBC:nx+nbx_HaloBC,j,1:nz)

                ! left halos
                var3D_HaloBC(-nbx_HaloBC:nx+nbx_HaloBC,-nby_HaloBC+j,1:nz) = yRhoSliceBack_recv(-nbx_HaloBC:nx+nbx_HaloBC,j,1:nz)

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
       stop "setBoundary: unknown case zBoundary"
    end select

          return

  end subroutine setHaloAndBoundary


! modified by Junhong Wei (20160726)  (finishing line)
!---------------------------------------------------------------------


end module update_module

