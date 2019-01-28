module flux_module

  use type_module
  use xweno_module
  use muscl_module
  use atmosphere_module
  use algebra_module

  implicit none

  private
  ! all module objects (variables, functions, subroutines)
  ! are internal to the module by default

  ! public routines
  public :: reconstruction
  public :: massFlux
  public :: thetaFlux
  public :: momentumFlux
  public :: iceFlux
  public :: volumeForce
! achatzb old topographic scheme removed
! public :: bottomTopography
! achatze
  public :: init_fluxes 
  public :: terminate_fluxes
  public :: thetaSource
  public :: massSource
  public :: momentumSource
  public :: iceSource

  ! private routines
  private :: absDiff
  private :: slopeFunction


  ! internal module variables
  real, dimension(:,:,:),     allocatable :: rhoBar, rhoOld
  real, dimension(:,:,:),     allocatable :: uBar
  real, dimension(:,:,:),     allocatable :: vBar
  real, dimension(:,:,:),     allocatable :: wBar
  real, dimension(:,:,:),     allocatable :: thetaBar
  real, dimension(:,:,:),     allocatable :: nIceBar
  real, dimension(:,:,:),     allocatable :: qIceBar
  real, dimension(:,:,:),     allocatable :: SIceBar

  real, dimension(:,:,:,:,:), allocatable :: rhoTilde
  real, dimension(:,:,:,:,:), allocatable :: uTilde
  real, dimension(:,:,:,:,:), allocatable :: vTilde
  real, dimension(:,:,:,:,:), allocatable :: wTilde
  real, dimension(:,:,:,:,:), allocatable :: thetaTilde
  real, dimension(:,:,:,:,:), allocatable :: nIceTilde
  real, dimension(:,:,:,:,:), allocatable :: qIceTilde
  real, dimension(:,:,:,:,:), allocatable :: SIceTilde

  ! public variables
  ! needed for 
  ! 1) BC correction
  ! 2) explicit boundary setting
  ! 3) update module
  public :: rhoTilde, thetaTilde
  public :: uTilde, vTilde, wTilde
  public :: rhoOld
  public :: nIceTilde, qIceTilde, SIceTilde


  ! phiTilde(i,j,k,dir,Left/Right) with
  ! dir = 1,2,3 for x,y and z reconstruction directions
  ! 0 = Left and 1 = Right for left and right cell edge



contains

! subroutine bottomTopography (var,flux)
!   !------------------------------------------------------
!   ! sets w at k = 1 so that w/u = mountain slope
!   !------------------------------------------------------

!   ! in/out variables
!   real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
!        & intent(inout) :: var

!   real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), &
!        & intent(inout) :: flux
!   ! flux(i,j,k,dir,iFlux) 
!   ! dir = 1..3 > f,g,h-flux in x,y,z-direction
!   ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW

!   ! local vars
!   integer :: i,j,k, kMountain
!   real :: u, w_top, slope 
!   real :: rhoD, rhoU, w, hRho

!   ! mass flux correction
!   real :: wSurf

!   ! debugging
!   real :: checkSum


!   !---------------------------------------------
!   !  topography: set w at k = 1,...,kMountain
!   !---------------------------------------------

!   if( topography ) then

!      print*,"WARNING: topography is on!!!"

!      ! heighest level to be manipulated
!      kMountain = floor(mountainHeight_dim/lRef/dz)

!      do k = 1,kMountain  ! ground level for u
!         do j = 1,ny
!            do i = 1,nx

!               u = 0.5*( var(i,j,k,2)+var(i,j,k+1,2) )
!               slope = slopeFunction( x(i) )
!               w_top = slope * u

!               var(i,j,k,4) = w_top

!            end do
!         end do
!      end do

!   end if ! topography


! end subroutine bottomTopography


  !---------------------------------------------------------------------------


  real function slopeFunction ( xx )
    !
    ! calculates slope of "witch of agnesi" mountain 
    !
    real, intent(in) :: xx           ! x is reserved for grid

    ! local vars
    real :: l,h,s

    ! init vars

    ! scaled quantities
    h = mountainHeight_dim / lRef
    l = mountainWidth_dim / lRef
    s = xx/l        

    if( abs(s) < 5.0 ) then
       slopeFunction = -2.0*h/l * s/(1+s**2)**2
    else
       slopeFunction = 0.0
    end if


  end function slopeFunction


  !---------------------------------------------------------------------------


  subroutine reconstruction (var,variable)
    !--------------------------------------------------
    ! reconstructs "variable" with 
    ! SALD, ALDM, constant, MUSCL
    !-------------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var

    character(len=*), intent(in) :: variable
    integer :: i,j,k

    ! locals
    integer :: dir, lr

    select case (reconstType)



       !----------------------------
       !   Constant reconstruction
       !----------------------------
       
    case( 'constant2' ) 
       
       select case( variable )
          
       case( "rho" )

          do dir = 1,3
             do lr = 0,1
                rhoTilde(:,:,:,dir,lr) = var(:,:,:,1)
             end do
          end do


       case( "uvw" )

          do dir = 1,3
             do lr = 0,1
                uTilde(:,:,:,dir,lr) = var(:,:,:,2)
                vTilde(:,:,:,dir,lr) = var(:,:,:,3)
                wTilde(:,:,:,dir,lr) = var(:,:,:,4)
             end do
          end do


       case( "theta" )

          do dir = 1,3
             do lr = 0,1
                thetaTilde(:,:,:,dir,lr) = var(:,:,:,6)
             end do
          end do

          
       case( "ice" )
         if (include_ice) then
          do dir = 1,3
             do lr = 0,1
                nIceTilde(:,:,:,dir,lr) = var(:,:,:,nVar-2)
                qIceTilde(:,:,:,dir,lr) = var(:,:,:,nVar-1)
                SIceTilde(:,:,:,dir,lr) = var(:,:,:,nVar)
             end do
          end do
         end if

       case default
          stop "reconstruction: unknown case model."
       end select



 

       !---------------------------
       !   MUSCL reconstrcution
       !---------------------------

    case( 'MUSCL' ) 

       select case( variable )

       case( "rho" )

          rhoBar(:,:,:) = var(:,:,:,1)
          call reconstruct_MUSCL(rhoBar,rhoTilde,nxx,nyy,nzz,limiterType1)

       case( "uvw" )

          uBar(:,:,:)   = var(:,:,:,2)
          vBar(:,:,:)   = var(:,:,:,3)
          wBar(:,:,:)   = var(:,:,:,4)

          call reconstruct_MUSCL(uBar,uTilde,nxx,nyy,nzz,limiterType1)
          call reconstruct_MUSCL(vBar,vTilde,nxx,nyy,nzz,limiterType1)
          call reconstruct_MUSCL(wBar,wTilde,nxx,nyy,nzz,limiterType1)

       case( "theta" )

          thetaBar(:,:,:) = var(:,:,:,6)
          call reconstruct_MUSCL(thetaBar,thetaTilde,nxx,nyy,nzz,limiterType1)

       case ( "ice" )

         if (include_ice) then
           
          nIceBar(:,:,:)   = var(:,:,:,nVar-2)
          qIceBar(:,:,:)   = var(:,:,:,nVar-1)
          SIceBar(:,:,:)   = var(:,:,:,nVar)

          call reconstruct_MUSCL(nIceBar,nIceTilde,nxx,nyy,nzz,limiterType1)
          call reconstruct_MUSCL(qIceBar,qIceTilde,nxx,nyy,nzz,limiterType1)
          call reconstruct_MUSCL(SIceBar,SIceTilde,nxx,nyy,nzz,limiterType1)

         end if

       case default
          stop "reconstruction: unknown case model."
       end select




       !---------------------------
       !     no reconstruction
       !---------------------------

    case( 'constant' )

       return


       !---------------------------
       !    SALD reconstruction
       !---------------------------

    case( 'SALD' )  ! simplified ALDM using reconstruct_SALD

       select case( variable )

       case( "rho" )

          rhoBar(:,:,:) = var(:,:,:,1)
          call reconstruct_SALD(rhoBar,rhoTilde)

       case( "uvw" )

          uBar(:,:,:)   = var(:,:,:,2)
          vBar(:,:,:)   = var(:,:,:,3)
          wBar(:,:,:)   = var(:,:,:,4)

          call reconstruct_SALD(uBar,uTilde)
          call reconstruct_SALD(vBar,vTilde)
          call reconstruct_SALD(wBar,wTilde)

       case( "theta" ) 

          thetaBar(:,:,:) = var(:,:,:,6)
          call reconstruct_SALD(thetaBar, thetaTilde)

       case ( "ice" )

         if (include_ice) then
           
          nIceBar(:,:,:)   = var(:,:,:,nVar-2)
          qIceBar(:,:,:)   = var(:,:,:,nVar-1)
          SIceBar(:,:,:)   = var(:,:,:,nVar)

          call reconstruct_SALD(nIceBar,nIceTilde)
          call reconstruct_SALD(qIceBar,qIceTilde)
          call reconstruct_SALD(SIceBar,SIceTilde)

         end if

       case default
          stop "reconstruction: unknown case variable."
       end select


       !---------------------------
       !    ALDM reconstruction
       !---------------------------

    case( 'ALDM' ) ! full 3D reconstruction according to ALDM

       select case( variable ) 

       case( "rho" )

          rhoBar(:,:,:) = var(:,:,:,1)
          call reconstruct_ALDM(rhoBar,rhoTilde)

       case( "uvw" )

          uBar(:,:,:)   = var(:,:,:,2)
          vBar(:,:,:)   = var(:,:,:,3)
          wBar(:,:,:)   = var(:,:,:,4)

          call reconstruct_ALDM(uBar,uTilde)
          call reconstruct_ALDM(vBar,vTilde)
          call reconstruct_ALDM(wBar,wTilde)

       case( "theta" ) 

          thetaBar(:,:,:) = var(:,:,:,6)
          call reconstruct_ALDM(thetaBar,thetaTilde)

       case ( "ice" )

         if (include_ice) then
           
          nIceBar(:,:,:)   = var(:,:,:,nVar-2)
          qIceBar(:,:,:)   = var(:,:,:,nvar-1)
          SIceBar(:,:,:)   = var(:,:,:,nVar)

          call reconstruct_ALDM(nIceBar,nIceTilde)
          call reconstruct_ALDM(qIceBar,qIceTilde)
          call reconstruct_ALDM(SIceBar,SIceTilde)

         end if

       case default
          stop "reconstruction: unknown case variable."
       end select

    case default
       print*, "reconstruction: unknown reconstruction type."
    end select




  end subroutine reconstruction


  !---------------------------------------------------------------------------


  subroutine thetaSource (var,source)
    !---------------------------------------------------------------------
    ! computes theta*div(u) for reconstructed u, which do not satisfy
    ! the divergence constraint -> source term corrects flux difference
    !---------------------------------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var

    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: source



    integer :: i,j,k
    real :: uL,uR       ! L=Left i-1/2, R=Right i+1/2
    real :: vB,vF       ! B=Backward j-1/2, F=Forward j+1/2
    real :: wD,wU       ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurfL, vSurfF, wSurfU     ! velocities at cell surface
    real :: uSurfR, vSurfB, wSurfD     ! velocities at cell surface
    real :: div, theta


    !--------------------------------------------------------
    !                  Divergence correction
    !--------------------------------------------------------
    !   u*grad(theta) = div(theta*u) - theta*div(u)

    do k = 1,nz
       do j = 1,ny
          do i = 1,nx


             select case( fluxType )

             case( "central" )
                return
                
             case( "upwind","ILES" )

                uL = uTilde(i,j,k,1,0)
                uR = uTilde(i,j,k,1,1)
                uSurfR = 0.5*(uL + uR)
                
                uL = uTilde(i-1,j,k,1,0)
                uR = uTilde(i-1,j,k,1,1)
                uSurfL = 0.5*(uL + uR)

                vB = vTilde(i,j,k,2,0)
                vF = vTilde(i,j,k,2,1)
                vSurfF = 0.5*(vB + vF)
                
                vB = vTilde(i,j-1,k,2,0)
                vF = vTilde(i,j-1,k,2,1)
                vSurfB = 0.5*(vB + vF)

                wD = wTilde(i,j,k,3,0)
                wU = wTilde(i,j,k,3,1)
                wSurfU = 0.5*(wD + wU)

                wD = wTilde(i,j,k-1,3,0)
                wU = wTilde(i,j,k-1,3,1)
                wSurfD = 0.5*(wD + wU)

                div = (uSurfR - uSurfL)/dx &
                     & + (vSurfF - vSurfB) / dy &
                     & + (wSurfU - wSurfD) / dz

                theta = var(i,j,k,6)

                source(i,j,k,6) = theta*div
                
                
             case default
                stop "thetaFlux: unknown case fluxType"
             end select

             
          end do
       end do
    end do

    
  end subroutine thetaSource


  !---------------------------------------------------------------------------


  subroutine thetaFlux (var,flux)
    !---------------------------------------------------------------------
    ! computes the theta flux at all cell edges using reconstructed values
    !---------------------------------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var

    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), &
         & intent(out) :: flux
    ! flux(i,j,k,dir,iFlux) 
    ! dir = 1..3 > f,g,h-flux in x,y,z-direction
    ! iFlux = 1..nVar > fRho, fRhoU, rRhoV, fRhoW, fTheta


    integer :: i,j,k,l
    real :: thetaL,thetaR, uL,uR       ! L=Left i-1/2, R=Right i+1/2
    real :: thetaB,thetaF, vB,vF       ! B=Backward j-1/2, F=Forward j+1/2
    real :: thetaD,thetaU, wD,wU       ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurf, vSurf, wSurf    ! velocities at cell surface

    real  :: fTheta, gTheta, hTheta

    ! avoid abs() for linerisation
    real :: delta
    real, parameter :: delta0 = 1.0e-6



    !--------------------------------------------------------
    !                          Advection
    !--------------------------------------------------------



    !-----------------------------------------
    !       Zonal theta fluxes in x: f
    !-----------------------------------------

    do k = 1,nz
       do j = 1,ny
          do i = 0,nx



             select case( fluxType )

             case( "central" )
                thetaL = var(i,j,k,6)
                thetaR = var(i+1,j,k,6)
                uSurf = var(i,j,k,2)
                fTheta = uSurf * 0.5*(thetaL + thetaR)

             case( "upwind" )
                thetaR = thetaTilde(i+1,j,k,1,0)
                thetaL = thetaTilde(i,j,k,1,1)
                uL = uTilde(i,j,k,1,0)
                uR = uTilde(i,j,k,1,1)
                uSurf = 0.5*(uL + uR)
                fTheta = flux_muscl(uSurf,thetaL,thetaR)

             case( "ILES" )
                thetaR = thetaTilde(i+1,j,k,1,0)
                thetaL = thetaTilde(i,j,k,1,1)
                uL = uTilde(i,j,k,1,0)
                uR = uTilde(i,j,k,1,1)
                uSurf = 0.5*(uL + uR)
                fTheta = flux_aldm(thetaL,thetaR,uSurf,&
                     &             thetaL,thetaR,uL,uR,sigmaC)
             case default
                stop "thetaFlux: unknown case fluxType"
             end select

             flux(i,j,k,1,6) = fTheta
          end do
       end do
    end do


    !-----------------------------------------
    !    Meridional theta fluxes in y: g
    !-----------------------------------------

    do k = 1,nz
       do j = 0,ny
          do i = 1,nx



             select case( fluxType )

             case( "central" )
                thetaF = var(i,j+1,k,6)
                thetaB = var(i,j,k,6)
                vSurf = var(i,j,k,3)
                gTheta = vSurf * 0.5*(thetaB + thetaF)
                
             case( "upwind" )
                thetaF = thetaTilde(i,j+1,k,2,0)
                thetaB = thetaTilde(i,j,k,2,1)
                vB = vTilde(i,j,k,2,0)
                vF = vTilde(i,j,k,2,1)
                vSurf = 0.5*(vB + vF)
                gTheta = flux_muscl(vSurf,thetaB,thetaF)
                
             case( "ILES" )
                thetaF = thetaTilde(i,j+1,k,2,0)
                thetaB = thetaTilde(i,j,k,2,1)
                vB = vTilde(i,j,k,2,0)
                vF = vTilde(i,j,k,2,1)
                vSurf = 0.5*(vB + vF)
                gTheta = flux_aldm(thetaB,thetaF,vSurf,&
                     &             thetaB,thetaF,vB,vF,sigmaC)

             case default
                stop "thetaFlux: unknown case fluxType"
             end select

             flux(i,j,k,2,6) = gTheta

          end do
       end do
    end do


    !-----------------------------------------
    !      Vertical theta fluxes in z: h
    !-----------------------------------------

    do k = 0,nz
       do j = 1,ny
          do i = 1,nx

             select case( fluxType )

             case( "central" )
                thetaU = var(i,j,k+1,6)
                thetaD = var(i,j,k,6)
                wSurf = var(i,j,k,4)
                hTheta = wSurf * 0.5*(thetaD + thetaU)                   

             case( "upwind" )
                thetaU = thetaTilde(i,j,k+1,3,0)
                thetaD = thetaTilde(i,j,k,3,1)
                wD = wTilde(i,j,k,3,0)
                wU = wTilde(i,j,k,3,1)
                wSurf = 0.5*(wD + wU)
                hTheta = flux_muscl(wSurf,thetaD,thetaU)

             case( "ILES" )
                thetaU = thetaTilde(i,j,k+1,3,0)
                thetaD = thetaTilde(i,j,k,3,1)
                wD = wTilde(i,j,k,3,0)
                wU = wTilde(i,j,k,3,1)
                wSurf = 0.5*(wD + wU)
                hTheta = flux_aldm(thetaD,thetaU,wSurf,&
                     &             thetaD,thetaU,wU,wD,sigmaC)
                
             case default
                stop "thetaFlux: unknown case fluxType"
             end select

             flux(i,j,k,3,6) = hTheta

          end do
       end do
    end do



    !--------------------------------------------------------------
    !                      Heat conduction
    !--------------------------------------------------------------


    if( mu_conduct > 0.0 ) then


       !-----------------------------------------
       !       Zonal theta fluxes in x: f
       !-----------------------------------------

       do k = 1,nz
          do j = 1,ny
             do i = 0,nx

                thetaL = var(i,j,k,6)
                thetaR = var(i+1,j,k,6)
                fTheta = mu_conduct * (thetaR-thetaL)/dx

                flux(i,j,k,1,6) = flux(i,j,k,1,6) - fTheta
             end do
          end do
       end do


       !-----------------------------------------
       !    Meridional theta fluxes in y: g
       !-----------------------------------------

       do k = 1,nz
          do j = 0,ny
             do i = 1,nx

                thetaF = var(i,j+1,k,6)
                thetaB = var(i,j,k,6)
                gTheta = mu_conduct * (thetaF-thetaB)/dy

                flux(i,j,k,2,6) = flux(i,j,k,2,6) - gTheta
             end do
          end do
       end do


       !-----------------------------------------
       !      Vertical theta fluxes in z: h
       !-----------------------------------------

       do k = 0,nz
          do j = 1,ny
             do i = 1,nx

                thetaU = var(i,j,k+1,6)
                thetaD = var(i,j,k,6)
                hTheta = mu_conduct * (thetaU-thetaD)/dz

                flux(i,j,k,3,6) = flux(i,j,k,3,6) - hTheta
             end do
          end do
       end do

    end if  ! mu_conduct > 0.0

    
    if (verbose) print*,"thetaFlux: &
         &theta fluxes fTheta, gTheta and fTheta calculated"

  end subroutine thetaFlux


  !---------------------------------------------------------------------------


  subroutine massSource (var,source)
    !---------------------------------------------------------------------
    ! computes source term for the conti equation
    ! 1) divergence term rho*div(u)
    !---------------------------------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var

    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: source



    integer :: i,j,k
    real :: uL,uR       ! L=Left i-1/2, R=Right i+1/2
    real :: vB,vF       ! B=Backward j-1/2, F=Forward j+1/2
    real :: wD,wU       ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurfL, vSurfF, wSurfU     ! velocities at cell surface
    real :: uSurfR, vSurfB, wSurfD     ! velocities at cell surface
    real :: PstratU, PstratD
    real :: divPu


    !--------------------------------------------------------
    !                  Divergence correction
    !--------------------------------------------------------
    !   u*grad(rho) = div(rho*u) - rho*div(u)

    do k = 1,nz
       do j = 1,ny
          do i = 1,nx


             select case( fluxType )

             case( "central" )
                return
                
             case( "upwind","ILES" )

                uR = var(i,j,k,2); uL = var(i-1,j,k,2)
                vF = var(i,j,k,3); vB = var(i,j-1,k,3)
                wU = var(i,j,k,4); wD = var(i,j,k-1,4)

                PstratU = PstratTilde(k)
                PstratD = PstratTilde(k-1)

                divPu  = Pstrat(k) * ( (uR-uL)/dx + (vF-vB)/dy ) &
                     & + (PstratU*wU - PstratD*wD)/dz

!!$                uSurfR = var(i,j,k,2)
!!$                uSurfL = var(i-1,j,k,2)
!!$
!!$                vSurfF = var(i,j,k,3)
!!$                vSurfB = var(i,j-1,k,3)
!!$
!!$                wSurfU = var(i,j,k,4)
!!$                wSurfD = var(i,j,k-1,4)
!!$
!!$                div = (uSurfR - uSurfL)/dx &
!!$                     & + (vSurfF - vSurfB) / dy &
!!$                     & + (wSurfU - wSurfD) / dz

                source(i,j,k,1) = divPu/thetaStrat(k)
                
             case default
                stop "rhoFlux: unknown case fluxType"
             end select

             
          end do
       end do
    end do

    
  end subroutine massSource


  !---------------------------------------------------------------------------


  subroutine massFlux (var,flux)
    !---------------------------------------------------------------------
    ! computes the mass flux at all cell edges using reconstructed values
    !---------------------------------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var

    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), &
         & intent(out) :: flux
    ! flux(i,j,k,dir,iFlux) 
    ! dir = 1..3 > f,g,h-flux in x,y,z-direction
    ! iFlux = 1..nVar > fRho, fRhoU, rRhoV, fRhoW, fRho


    integer :: i,j,k,l
    real :: rhoL,rhoR, uL,uR       ! L=Left i-1/2, R=Right i+1/2
    real :: rhoB,rhoF, vB,vF       ! B=Backward j-1/2, F=Forward j+1/2
    real :: rhoD,rhoU, wD,wU       ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurf, vSurf, wSurf    ! velocities at cell surface

    real  :: fRho, gRho, hRho

    ! avoid abs() for linerisation
    real :: delta
    real, parameter :: delta0 = 1.0e-6


!   variables for the Dynamic Smagorinsky scheme      
!   modified by Junhong Wei (20160803)
    real :: coef_t_DySma, Pr_t_DySma, drho_dxi_DySma       
!   achatzb for inclusion density dependent diffusivity and viscosity
    real :: coef_t, drho_dxi
!   achatze

    Pr_t_DySma = 0.5  ! modified by Junhong Wei (20160803)



    !-----------------------------------------
    !       Zonal rho fluxes in x: f
    !-----------------------------------------

    do k = 1,nz
       do j = 1,ny
          do i = 0,nx



             select case( fluxType )

             case( "central" )

                if( fluctuationMode ) then
                   rhoL = var(i,j,k,1) + rhoStrat(k)
                   rhoR = var(i+1,j,k,1) + rhoStrat(k)
                else
                   rhoL = var(i,j,k,1)
                   rhoR = var(i+1,j,k,1)
                end if

                uSurf = var(i,j,k,2)
                fRho = uSurf * 0.5*(rhoL + rhoR)

             case( "upwind" )

                if( fluctuationMode ) then
                   rhoR = rhoTilde(i+1,j,k,1,0) + rhoStrat(k)
                   rhoL = rhoTilde(i,j,k,1,1)   + rhoStrat(k)
                else
                   rhoR = rhoTilde(i+1,j,k,1,0)
                   rhoL = rhoTilde(i,j,k,1,1)
                end if
                
                uSurf = var(i,j,k,2)
                fRho = flux_muscl(uSurf,rhoL,rhoR)

             case( "ILES" )
                rhoR = rhoTilde(i+1,j,k,1,0)
                rhoL = rhoTilde(i,j,k,1,1)
                uL = uTilde(i,j,k,1,0)
                uR = uTilde(i,j,k,1,1)
! variant a)
                uSurf = var(i,j,k,2)

! variant b)
!                uSurf = 0.5*(uL + uR)
                if( fluctuationMode ) then
                   fRho = flux_aldm(rhoL + rhoStrat(k),rhoR + rhoStrat(k),uSurf,&
                        &           rhoL,rhoR,uL,uR,sigmaC)
                else
                   fRho = flux_aldm(rhoL,rhoR,uSurf,&
                        &           rhoL,rhoR,uL,uR,sigmaC)
                end if

             case default
                stop "rhoFlux: unknown case fluxType"
             end select



!               achatzb
!               inclusion density dependent diffusivity

!               ! modified by Junhong Wei (20160803) --- starting line
!               if(DySmaScheme)then

!               coef_t_DySma   = 0.5*( var(i,j,k,7) + var(i+1,j,k,7) )

!               if( fluctuationMode ) then
!                  rhoL = var(i,j,k,1) + rhoStrat(k)
!                  rhoR = var(i+1,j,k,1) + rhoStrat(k)
!               else
!                  rhoL = var(i,j,k,1)
!                  rhoR = var(i+1,j,k,1)
!               end if

!               drho_dxi_DySma = ( rhoR - rhoL ) / dx

!               fRho           &
!               = fRho - ( coef_t_DySma * drho_dxi_DySma / Pr_t_DySma )

!               end if
!               ! modified by Junhong Wei (20160803) --- finishing line

                select case( model )
                  case( "pseudo_incompressible" )
                   coef_t = mu_conduct * rhoStrat(1)/rhoStrat(k)
                  case( "Boussinesq" )
                   coef_t = mu_conduct
                  case default
                   stop "diffusivity: unkown case model."
                end select

                if(DySmaScheme)then
                   coef_t  &
                   = coef_t &
                     + 0.5*( var(i,j,k,7) + var(i+1,j,k,7) )/Pr_t_DySma
                end if

                if( fluctuationMode ) then
                   rhoL = var(i,j,k,1) + rhoStrat(k)
                   rhoR = var(i+1,j,k,1) + rhoStrat(k)
                else
                   rhoL = var(i,j,k,1)
                   rhoR = var(i+1,j,k,1)
                end if

                drho_dxi = ( rhoR - rhoL ) / dx

                fRho = fRho - coef_t * drho_dxi
!               achatze



             flux(i,j,k,1,1) = fRho
          end do
       end do
    end do


    !-----------------------------------------
    !    Meridional rho fluxes in y: g
    !-----------------------------------------

    do k = 1,nz
       do j = 0,ny
          do i = 1,nx

             select case( fluxType )

             case( "central" )
                if( fluctuationMode ) then
                   rhoF = var(i,j+1,k,1) + rhoStrat(k)
                   rhoB = var(i,j,k,1)   + rhoStrat(k)
                else
                   rhoF = var(i,j+1,k,1)
                   rhoB = var(i,j,k,1)
                end if

                vSurf = var(i,j,k,3)
                gRho = vSurf * 0.5*(rhoB + rhoF)

             case( "upwind" )
                if( fluctuationMode ) then
                   rhoF = rhoTilde(i,j+1,k,2,0) + rhoStrat(k)
                   rhoB = rhoTilde(i,j,k,2,1)   + rhoStrat(k)
                else
                   rhoF = rhoTilde(i,j+1,k,2,0)
                   rhoB = rhoTilde(i,j,k,2,1)
                end if

                vSurf = var(i,j,k,3)
                gRho = flux_muscl(vSurf,rhoB,rhoF)

             case( "ILES" )
                rhoF = rhoTilde(i,j+1,k,2,0)
                rhoB = rhoTilde(i,j,k,2,1)
                vB = vTilde(i,j,k,2,0)
                vF = vTilde(i,j,k,2,1)
! variant a)
                vSurf = var(i,j,k,3)
! variant b)
!                vSurf = 0.5*(vF + vB)

                if( fluctuationMode ) then
                   gRho = flux_aldm(rhoB+rhoStrat(k),rhoF+rhoStrat(k),vSurf,&
                        &           rhoB,rhoF,vB,vF,sigmaC)
                else                
                   gRho = flux_aldm(rhoB,rhoF,vSurf,&
                        &           rhoB,rhoF,vB,vF,sigmaC)
                end if

             case default
                stop "rhoFlux: unknown case fluxType"
             end select



!               achatzb
!               inclusion density dependent diffusivity

!               ! modified by Junhong Wei (20160803) --- starting line
!               if(DySmaScheme)then

!               coef_t_DySma   = 0.5*( var(i,j,k,7) + var(i,j+1,k,7) )

!               if( fluctuationMode ) then
!                  rhoF = var(i,j+1,k,1) + rhoStrat(k)
!                  rhoB = var(i,j,k,1)   + rhoStrat(k)
!               else
!                  rhoF = var(i,j+1,k,1)
!                  rhoB = var(i,j,k,1)
!               end if

!               drho_dxi_DySma = ( rhoF - rhoB ) / dy

!               gRho           &
!               = gRho - ( coef_t_DySma * drho_dxi_DySma / Pr_t_DySma )

!               end if
!               ! modified by Junhong Wei (20160803) --- finishing line

                select case( model )
                  case( "pseudo_incompressible" )
                   coef_t = mu_conduct * rhoStrat(1)/rhoStrat(k)
                  case( "Boussinesq" )
                   coef_t = mu_conduct
                  case default
                   stop "diffusivity: unkown case model."
                end select

                if(DySmaScheme)then
                   coef_t  &
                   = coef_t &
                     + 0.5*( var(i,j,k,7) + var(i+1,j,k,7) )/Pr_t_DySma
                end if

                if( fluctuationMode ) then
                   rhoF = var(i,j+1,k,1) + rhoStrat(k)
                   rhoB = var(i,j,k,1)   + rhoStrat(k)
                else
                   rhoF = var(i,j+1,k,1)
                   rhoB = var(i,j,k,1)
                end if

                drho_dxi = ( rhoF - rhoB ) / dy

                gRho = gRho -  coef_t * drho_dxi
!               achatze



             flux(i,j,k,2,1) = gRho

          end do
       end do
    end do


    !-----------------------------------------
    !      Vertical rho fluxes in z: h
    !-----------------------------------------

! xxx fluctuatonMode: check wheter rhoStrat has to be changed at k = nz+1, and k = 0

    do k = 0,nz
       do j = 1,ny
          do i = 1,nx

             select case( fluxType )

             case( "central" )
                if( fluctuationMode ) then
                   rhoU = var(i,j,k+1,1) + rhoStratTilde(k)   
                   ! background rho at half level     
                   rhoD = var(i,j,k,1)   + rhoStratTilde(k)
                else
                   rhoU = var(i,j,k+1,1)
                   rhoD = var(i,j,k,1)
                end if

                wSurf = var(i,j,k,4)
                hRho = wSurf * 0.5*(rhoD + rhoU)                   

             case( "upwind" )

                if( fluctuationMode ) then
                   rhoU = rhoTilde(i,j,k+1,3,0) + rhoStratTilde(k)  ! background at half level
                   rhoD = rhoTilde(i,j,k,3,1)   + rhoStratTilde(k)
                else
                   rhoU = rhoTilde(i,j,k+1,3,0)
                   rhoD = rhoTilde(i,j,k,3,1)
                end if

                wSurf = var(i,j,k,4)
                hRho = flux_muscl(wSurf,rhoD,rhoU)

             case( "ILES" )
                rhoU = rhoTilde(i,j,k+1,3,0)
                rhoD = rhoTilde(i,j,k,3,1)
                wD = wTilde(i,j,k,3,0)
                wU = wTilde(i,j,k,3,1)
! variant a)
                wSurf = var(i,j,k,4)
! variant b)
!                wSurf = 0.5*(wU + wD)    ! crahsed after heavy oscillations

                if( fluctuationMode ) then
                   hRho = flux_aldm(rhoD+rhoStratTilde(k),rhoU+rhoStratTilde(k),wSurf,&
                        &             rhoD,rhoU,wU,wD,sigmaC)
                else
                   hRho = flux_aldm(rhoD,rhoU,wSurf,&
                        &             rhoD,rhoU,wU,wD,sigmaC)
                end if

             case default
                stop "rhoFlux: unknown case fluxType"
             end select




!               achatzb
!               inclusion density dependent diffusivity

!               ! modified by Junhong Wei (20160803) --- starting line
!               if(DySmaScheme)then

!               coef_t_DySma   = 0.5*( var(i,j,k,7) + var(i,j,k+1,7) )

!               !UA here I am not sure whether it would not be better to 
!               !rather take the vertical derivative of the density 
!               !without the contribution from the reference atmosphere. 
!               !Maybe try?

!               ! replied by JW (20160825): I haven't done anything for 
!               ! the above comment. Please check. 

!               if( fluctuationMode ) then
!               !   rhoU = var(i,j,k+1,1) + rhoStrat(k+1)
!               !   ! background rho at half level     
!               !   rhoD = var(i,j,k,1)   + rhoStrat(k)

!                  rhoU = var(i,j,k+1,1)
!                  ! background rho at half level     
!                  rhoD = var(i,j,k,1)
!               else
!               !   rhoU = var(i,j,k+1,1)
!               !   rhoD = var(i,j,k,1)

!                  rhoU = var(i,j,k+1,1) - rhoStrat(k+1)
!                  rhoD = var(i,j,k,1)   - rhoStrat(k)
!               end if

!               drho_dxi_DySma = ( rhoU - rhoD ) / dz

!               hRho           &
!               = hRho - ( coef_t_DySma * drho_dxi_DySma / Pr_t_DySma )

!               end if
!               ! modified by Junhong Wei (20160803) --- finishing line

                select case( model )
                  case( "pseudo_incompressible" )
                   coef_t = mu_conduct * rhoStrat(1)/rhoStrat(k)
                  case( "Boussinesq" )
                   coef_t = mu_conduct
                  case default
                   stop "diffusivity: unkown case model."
                end select

                if(DySmaScheme)then
                   coef_t  &
                   = coef_t &
                     + 0.5*( var(i,j,k,7) + var(i+1,j,k,7) )/Pr_t_DySma
                end if

                if( fluctuationMode ) then
                   rhoU = var(i,j,k+1,1)
                   rhoD = var(i,j,k,1)
                else
                   rhoU = var(i,j,k+1,1) - rhoStrat(k+1)
                   rhoD = var(i,j,k,1)   - rhoStrat(k)
                end if

                drho_dxi = ( rhoU - rhoD ) / dz

                hRho = hRho - coef_t * drho_dxi
!               achatze



             flux(i,j,k,3,1) = hRho

          end do
       end do
    end do

    if (verbose) print*,"rhoFlux: &
         &rho fluxes fRho, gRho and fRho calculated"

  end subroutine massFlux

  
  !---------------------------------------------------------------------------

  subroutine iceFlux (var, flux)
   ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var

    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), &
         & intent(out) :: flux
    ! flux(i,j,k,dir,iFlux) 
    ! dir = 1..3 > f,g,h-flux in x,y,z-direction

    integer :: nqS, dir, k, j, i

    ! All ice particles obey to the general mass flux
    do nqS = 0,2  
      do dir=1,3
        do k = 0,nz
          do j = 1,ny
            do i = 1,nx
              flux(i,j,k,dir,nVar-nqS) = var(i,j,k,nVar-nqS) * flux(i,j,k,dir,1)
            end do
          end do
        end do
      end do
    end do

  end subroutine iceFlux


  subroutine iceSource (var, source)
    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var

    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: source

    integer :: nqS, k, j, i

    do nqS = 0,2  
        do k = 0,nz
          do j = 1,ny
            do i = 1,nx
              source(i,j,k,nVar-nqS) = var(i,j,k,nVar-nqS) * source(i,j,k,1)
            end do
          end do
        end do
    end do
    
  end  subroutine iceSource
  
!-----------------------------------------------------------------------!

  
  subroutine volumeForce (var,time,force)
    !-------------------------------------------------------
    ! computes and adds up all volume forces on a grid cell
    !-------------------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var

    ! volume forces 
    real, dimension(0:nx+1,0:ny+1,0:nz+1,3), intent(out) :: force

    ! local variables
    integer :: i,j,k,l

    ! gravitational forces
    real :: dRho, theta
    real, dimension(3) :: gForce
    
    ! Coriolis force
    real :: n1,n2,n3
    real :: u1,u2,u3
    real :: f1,f2,f3
    real :: rho

!   achatzb
    ! zonal wind relaxation
    real :: ft_relax
    real :: time
    real, dimension(0:nx+1) :: fx_relax
    real :: xextent_relax
    integer :: i0
!   achatze

    
    !--------------------------------------------
    !             Gravitational force 
    !--------------------------------------------
    
    do k = 0,nz+1
       do j = 0,ny+1
          do i = 0,nx+1

             select case( model ) 

             case( "Boussinesq" ) 

                theta = var(i,j,k,6)
                gForce = FrInv2*rho00/theta00 * theta * vertical

             case( "pseudo_incompressible" )

                if( fluctuationMode) then
                   dRho = var(i,j,k,1)
                else
                   dRho = var(i,j,k,1) - rhoStrat(k)
                end if
                
                gForce = -FrInv2 * dRho * vertical
                
             case default
                stop "volumeForce: unknown case model."
             end select

             force(i,j,k,:) = gForce

          end do
       end do
    end do

!   testb
!   print*,"after g:"
!   print*,"f(0,1,310), f(1,1,310) = ", force(0,1,310,1), force(1,1,310,1)
!   print*,"f(1,1,310), f(2,1,310) = ", force(1,1,310,1), force(2,1,310,1)
!   print*,"f(32,1,310), f(33,1,310) = ", &
!        & force(32,1,310,1), force(33,1,310,1)
!   teste

    !--------------------------------------------
    !             Coriolis force 
    !--------------------------------------------
    ! Coriolis force defined in scalar volume cells
    ! -> needs interpolation
    
    if( RoInv > 0.0 ) then

       do k = 0,nz+1
          do j = 0,ny+1
             do i = 0,nx+1

                u1 = 0.5*( var(i,j,k,2) + var(i-1,j,k,2) )
                u2 = 0.5*( var(i,j,k,3) + var(i,j-1,k,3) )
                u3 = 0.5*( var(i,j,k,4) + var(i,j,k-1,4) )

                select case( model ) 

                case( "Boussinesq" ) 
                   rho = rho00

                case( "pseudo_incompressible" )
!                   rho = var(i,j,k,1)

                if( fluctuationMode) then
                rho = var(i,j,k,1) + rhoStrat(k)
                else
                rho = var(i,j,k,1)
                end if

                case default
                   stop "volumeForce: unknown case model."
                end select

                n1 = vertical(1)
                n2 = vertical(2)
                n3 = vertical(3)

                f1 = n2*u3 - n3*u2
                f2 = n3*u1 - n1*u3
                f3 = n1*u2 - n2*u1

                ! Coriolis force normally written in LHS with "+"
                ! gets now a "-" since force is assumed on the RHS

                force(i,j,k,:) = force(i,j,k,:) - ( rho*RoInv*(/f1,f2,f3/) )

             end do
          end do
       end do

    end if ! RoInv > 0.0

!   testb
!   print*,"after Coriolis:"

!   print*,"u(-1,1,310), u(0,1,310) = ", var(-1,1,310,2), var(0,1,310,2)
!   print*,"u(0,1,310), u(1,1,310) = ", var(0,1,310,2), var(1,1,310,2)
!   print*,"u(-1,1,310), u(0,1,310) = ", var(32,1,310,2), var(33,1,310,2)

!   print*,"v(0,1,310), v(0,0,310) = ", var(0,1,310,3), var(0,0,310,3)
!   print*,"v(1,1,310), v(1,0,310) = ", var(1,1,310,3), var(1,0,310,3)
!   print*,"v(32,1,310), v(32,0,310) = ", var(32,1,310,3), var(32,0,310,3)

!   print*,"leads to"

!   print*,"f(0,1,310), f(1,1,310) = ", force(0,1,310,1), force(1,1,310,1)
!   print*,"f(1,1,310), f(2,1,310) = ", force(1,1,310,1), force(2,1,310,1)
!   print*,"f(32,1,310), f(33,1,310) = ", &
!        & force(32,1,310,1), force(33,1,310,1)
!   teste

!   achatzb
    !--------------------------------------------
    !             zonal wind relaxation
    !--------------------------------------------

    i0=is+nbx-1
    
    if( testCase == "mountainwave") then

       if(time < t_ramp) then
          ft_relax = (1.0 - cos(time*pi/(t_ramp*2.0)))/t_relax
         else if (time < t_relax-t_ramp)  then
          ft_relax = 1.0/t_relax
         else if (time < t_relax) then
          ft_relax = (1.0 - cos((t_relax-time)*pi/(t_ramp*2.0)))/t_relax
         else
          ft_relax = 0.0
       end if

       xextent_relax = lx(1) - lx(0) - xextent_norelax

       ! if no-relaxation zone is larger than model domain then no 
       ! relaxation

       if(xextent_relax > 0.0) then
          do i = 0,nx+1
             if(x(i0+i) < lx(0)+0.5*xextent_relax) then
               fx_relax(i) = cos((x(i0+i)-lx(0))*pi/xextent_relax)
             else if(x(i0+i) < lx(1)-0.5*xextent_relax) then
               fx_relax(i) = 0.0
             else
               fx_relax(i) = cos((lx(1)-x(i0+i))*pi/xextent_relax)
             end if
!            testb
!            write(*,*)i,x(i0+i)*lRef,fx_relax(i)
!            teste
          end do
         else
          do i = 0,nx+1
             fx_relax(i) = 0.0
          end do
       end if

       do k = 0,nz+1
          do j = 0,ny+1
             do i = 0,nx+1

                select case( model ) 

                case( "Boussinesq" ) 
                   rho = rho00

                case( "pseudo_incompressible" )

                if( fluctuationMode) then
                rho = 0.5*(var(i,j,k,1)+var(i+1,j,k,1)) + rhoStrat(k)
                else
                rho = 0.5*(var(i,j,k,1)+var(i+1,j,k,1))
                end if

                case default
                   stop "volumeForce: unknown case model."
                end select

                force(i,j,k,1) &
                = force(i,j,k,1) &
                  - rho * ft_relax * fx_relax(i) * (var(i,j,k,2) - u_relax)

             end do
          end do
       end do

    end if ! mountainwave
!   achatze

  end subroutine volumeForce


  !---------------------------------------------------------------------------


  function flux_aldm(uB,uF,vSurf,vL,vR,vBarL,vBarR,sigma)
    !--------------------------------------------
    !   ALDM flux function (cf. Adams, Hickel)
    !--------------------------------------------

    ! in/out arguments
    real, intent(in) :: uB,uF   ! velocity to be transported 
    real, intent(in) :: vSurf   ! averaged cell face transport velocity
    real, intent(in) :: vL, vR  ! reconstructed transport velocity
    real, intent(in) :: vBarL,vBarR ! filtered transport velocities
    real, intent(in) :: sigma   ! ILES parameter
    real             :: flux_aldm


    flux_aldm = 0.5*(uB+uF)*vSurf - sigma*(vR-vL)*abs(vBarR-vBarL)


  end function flux_aldm


  !---------------------------------------------------------------------------


  function flux_muscl(uSurf, phiUp, phiDown)
    !----------------------------
    !   upwind flux function
    !----------------------------

    ! in/out arguments
    real, intent(in) :: uSurf           ! cell face value
    real, intent(in) :: phiUp, phiDown  ! upwind, downwind values
    real             :: flux_muscl

    !    flux_muscl = uSurf*0.5*(phiUp+phiDown) &
    !         &     + abs(uSurf)*0.5*(phiUp-phiDown)


    if( uSurf > 0.0 ) then
       flux_muscl = uSurf * phiUp 
    else
       flux_muscl = uSurf * phiDown
    end if


  end function flux_muscl


  !---------------------------------------------------------------------------


  subroutine momentumFlux (var,flux)
    !---------------------------------------------------------------------------
    ! computes the momentum fluxes at the cell edges using reconstructed values
    !---------------------------------------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var

    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), intent(out) :: flux
    ! flux(i,j,k,dir,iFlux) 
    ! dir = 1..3 > f,g,h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, fRhoV, fRhoW


    ! local variables
    integer :: i,j,k,l

    ! density at cell edge
    real :: rhoEdge

    ! uTilde at cell edges
    real :: uL,uR, vL,vR, wL,wR       ! L=Left at i+1, R=Right at i
    real :: uB,uF, vB,vF, wB,wF       ! B=Backward at j+1, F=Forward at j
    real :: uD,uU, vD,vU, wD,wU       ! D=Downward at k+1, U=Upward at k

    ! cell averaged values at cell centres
    real :: uBarL,uBarR, vBarL,vBarR, wBarL,wBarR   ! at i and i+1
    real :: uBarB,uBarF, vBarB,vBarF, wBarB, wBarF  ! at j and j+1
    real :: uBarD,uBarU, vBarD,vBarU, wBarD, wBarU  ! at k and k+1

    ! upwinding
    real :: uSurf, vSurf, wSurf

    ! local flux variables
    real :: fRhoU, gRhoU, hRhoU      ! rho*U momentum fluxes
    real :: fRhoV, gRhoV, hRhoV      ! rho*V momentum fluxes
    real :: fRhoW, gRhoW, hRhoW      ! rho*W momentum fluxes

    ! viscous fluxes
    real, dimension(0:nx+1,0:ny+1,0:nz+1) :: divU   ! div(u) field of viscosity
    real :: div        
    real :: du_dx, du_dy, du_dz                     ! partial derivatives 
    real :: dv_dx, dv_dy, dv_dz
    real :: dw_dx, dw_dy, dw_dz
    real :: fRhoU_visc, gRhoU_visc, hRhoU_visc      ! viscous momentum fluxes
    real :: fRhoV_visc, gRhoV_visc, hRhoV_visc     
    real :: fRhoW_visc, gRhoW_visc, hRhoW_visc

    ! debugging
    real :: rhoEdge2

    ! avoid abs() for linerisation
    real :: delta

    !for basic-state density
    real :: rhos

    !achatzb for putting together molecular and turbulent vioscosity
    real :: coef_v
    !achatze

    if(verbose) print*,"fluxes.f90/momentumFlux: Entering subroutine..."



    !------------------------------
    !     flux for rho*u
    !------------------------------

    ! flux fRhoU
    do k = 1,nz
       do j = 1,ny
          do i = -1,nx

             ! density at flux point
             select case( model ) 

             case( "Boussinesq" )
                rhoEdge = rho00

             case default

                select case( fluxType ) 

                case( "central" )
                   ! density interpolation consistent with conti eq

                   rhoEdge=0.25*(var(i,j,k,1) &
                        & + var(i+1,j,k,1) &
                        & + var(i+1,j,k,1) &
                        & + var(i+2,j,k,1) )

                case( "upwind", "ILES" )
                   ! density interpolation consistent with conti eq
                   rhoEdge=0.25*(rhoTilde(i,j,k,1,1) &
                        & + rhoTilde(i+1,j,k,1,0) &
                        & + rhoTilde(i+1,j,k,1,1) &
                        & + rhoTilde(i+2,j,k,1,0) )

                case default
                   stop "momentumFlux: unknown fluxType."
                end select
                
                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)
                
             end select ! model

             ! velocity
             select case( fluxType )

             case( "central" )
                uL = var(i,j,k,2)
                uR = var(i+1,j,k,2)

                fRhoU = 0.25*(uL+uR)**2

             case( "upwind", "ILES" )

                uR = uTilde(i+1,j,k,1,0)
                uL = uTilde(i,j,k,1,1) 
                uBarL = uBar(i,j,k)
                uBarR = uBar(i+1,j,k)

                uSurf = 0.5*(uL + uR)

                select case( fluxType ) 

                case( "upwind" )
                   fRhoU = flux_muscl(uSurf,uL,uR)
                case( "ILES" )
                   fRhoU = flux_aldm(uL,uR,uSurf,uL,uR,uBarL,uBarR,sigmaX)
                case default
                   stop "momentumFlux: unknown fluxType."
                end select

             case default
                stop "momentumFlux: unknown fluxType."
             end select

             flux(i,j,k,1,2) = rhoEdge * fRhoU

          end do
       end do
    end do

    !----------------------------------------------------------

    !  flux gRhoU
    do k = 1,nz
       do j = 0,ny
          do i = 0,nx

             ! density at flux point
             select case( model ) 

             case( "Boussinesq" )
                rhoEdge = rho00

             case default

                select case( fluxType ) 

                case( "central" )
                   ! density interpolation consistent with conti eq
                   rhoEdge = 0.25*( var(i+1,j,k,1) &
                        & + var(i+1,j+1,k,1) &
                        & + var(i,j,k,1) &
                        & + var(i,j+1,k,1) )

                case( "upwind", "ILES" )
                   rhoEdge = 0.25*( rhoTilde(i+1,j,k,2,1) &
                        & + rhoTilde(i+1,j+1,k,2,0) &
                        & + rhoTilde(i,j,k,2,1) &
                        & + rhoTilde(i,j+1,k,2,0) )
                case default 
                   stop "momentumFlux: unknown fluxType."
                end select
                
                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)
                
             end select ! model

             ! velocity
             select case( fluxType )

             case( "central" )
                uF = var(i,j+1,k,2)
                uB = var(i,j,k,2)
                vR = var(i+1,j,k,3)
                vL = var(i,j,k,3)

                gRhoU = 0.25*(uB+uF)*(vL+vR)

             case( "upwind", "ILES" )
                uF = uTilde(i,j+1,k,2,0)
                uB = uTilde(i,j,k,2,1)
                vR = vTilde(i+1,j,k,1,0)
                vL = vTilde(i,j,k,1,1)
                uBarF = uBar(i,j+1,k)
                uBarB = uBar(i,j,k)
                vSurf = 0.5*(vR + vL)

                select case( fluxType ) 

                case( "upwind" )
                   gRhoU = flux_muscl(vSurf,uB,uF)
                case( "ILES" )
                   gRhoU = flux_aldm(uB,uF,vSurf,uB,uF,uBarB,uBarF,sigmaX)
                case default
                   stop "momentumFlux: unknown fluxType."
                end select

             case default
                stop "momentumFlux: unknown fluxType."
             end select

             !                call absDiff(vBarR-vBarL,delta)
             !                gRhoU = 0.25*(uB+uF) * (vL+vR) - sigmaX*(vR-vL) * delta

             flux(i,j,k,2,2) = rhoEdge * gRhoU
          end do
       end do
    end do

    !-----------------------------------------------------------------------

    ! flux hRhoU
    do k = 0,nz
       do j = 1,ny
          do i = 0,nx

             ! density at flux point
             select case( model ) 

             case( "Boussinesq" )
                rhoEdge = rho00

             case default

                select case( fluxType ) 

                case( "central" )
                   ! density interpolation consistent with conti eq
                   rhoEdge = 0.25*( var(i,j,k,1) &
                        & + var(i,j,k+1,1) &
                        & + var(i+1,j,k,1) &
                        & + var(i+1,j,k+1,1)  )

                case( "upwind", "ILES" )
                   rhoEdge = 0.25*( rhoTilde(i,j,k,3,1) &
                        & + rhoTilde(i,j,k+1,3,0) &
                        & + rhoTilde(i+1,j,k,3,1) &
                        & + rhoTilde(i+1,j,k+1,3,0)  )
                   
                case default 
                   stop "momentumFlux: unknown fluxType."
                end select
                
                ! comment: for CDS rhoEdge should add rhoStrat
                ! for each var(...1) individually for 100% correctness
                
                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStratTilde(k)
                
             end select ! model

             ! velocity
             select case( fluxType )

             case( "central" )
                uU = var(i,j,k+1,2)
                uD = var(i,j,k,2)
                wR = var(i+1,j,k,4)
                wL = var(i,j,k,4)

                hRhoU = 0.25*(uD + uU)*(wL + wR)

             case( "upwind", "ILES" )
                uU = uTilde(i,j,k+1,3,0)
                uD = uTilde(i,j,k,3,1)
                wR = wTilde(i+1,j,k,1,0)
                wL = wTilde(i,j,k,1,1)
                uBarU = uBar(i,j,k+1)
                uBarD = uBar(i,j,k)
                wSurf = 0.5*(wL + wR)

                select case( fluxType ) 

                case( "upwind" )
                   hRhoU = flux_muscl(wSurf,uD,uU)
                case( "ILES" )
                   hRhoU = flux_aldm(uU,uD,wSurf,uD,uU,uBarD,uBarU,sigmaX)
                case default
                   stop "momentumFlux: unknown fluxType."
                end select

             case default
                stop "momentumFlux: unknown fluxType."
             end select

             flux(i,j,k,3,2) = rhoEdge * hRhoU

          end do
       end do
    end do


    !------------------------------
    !     flux for rho*v
    !------------------------------

    !  flux fRhoV
    do k = 1,nz
       do j = 0,ny
          do i = 0,nx

             ! density at flux point
             select case( model ) 

             case( "Boussinesq" )
                rhoEdge = rho00

             case default

                select case( fluxType ) 

                case( "central" )
                   ! density interpolation consistent with conti eq
                   rhoEdge = 0.25*( var(i,j,k,1) &
                        & + var(i+1,j,k,1) &
                        & + var(i,j+1,k,1) &
                        & + var(i+1,j+1,k,1) )

                case( "upwind", "ILES" )
                   rhoEdge = 0.25*( rhoTilde(i,j,k,1,1) &
                        & + rhoTilde(i+1,j,k,1,0) &
                        & + rhoTilde(i,j+1,k,1,1) &
                        & + rhoTilde(i+1,j+1,k,1,0) )

                case default 
                   stop "momentumFlux: unknown fluxType."
                end select

                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

             end select ! model

             ! velocity
             select case( fluxType )

             case( "central" )
                vR = var(i+1,j,k,3)
                vL = var(i  ,j,k,3)
                uF = var(i,j+1,k,2)
                uB = var(i,j  ,k,2)

                fRhoV = 0.25 * (vL+vR) * (uB+uF)

             case( "upwind", "ILES" )
                vR = vTilde(i+1,j,k,1,0)
                vL = vTilde(i  ,j,k,1,1)
                uF = uTilde(i,j+1,k,2,0)
                uB = uTilde(i,j  ,k,2,1)
                vBarR = vBar(i+1,j,k)
                vBarL = vBar(i,j  ,k)
                uSurf = 0.5*(uB + uF)

                select case( fluxType ) 

                case( "upwind" )
                   fRhoV = flux_muscl(uSurf,vL,vR)
                case( "ILES" )
                   fRhoV = flux_aldm(vL,vR,uSurf,vL,vR,vBarL,vBarR,sigmaY)
                case default
                   stop "momentumFlux: unknown fluxType."
                end select

             case default
                stop "momentumFlux: unknown fluxType."
             end select

             !                call absDiff(uBarF-uBarB,delta)
             !                fRhoV = 0.25*(vL+vR)*(uB+uF) - sigmaY*(uF-uB)*delta

             flux(i,j,k,1,3) = rhoEdge * fRhoV
          end do
       end do
    end do

    !---------------------------------------------------------------------

    ! flux gRhoV
    do k = 1,nz
       do j = -1,ny
          do i = 1,nx

             ! density at flux point
             select case( model ) 

             case( "Boussinesq" )
                rhoEdge = rho00

             case default

                select case( fluxType ) 

                case( "central" )
                   ! density interpolation consistent with conti eq
                   rhoEdge = 0.25*( var(i,j,k,1) &
                        & + var(i,j+1,k,1) &
                        & + var(i,j+1,k,1) &
                        & + var(i,j+2,k,1) )

                case( "upwind", "ILES" )
                   rhoEdge = 0.25*( rhoTilde(i,j,k,2,1) &
                        & + rhoTilde(i,j+1,k,2,0) &
                        & + rhoTilde(i,j+1,k,2,1) &
                        & + rhoTilde(i,j+2,k,2,0) )

                case default 
                   stop "momentumFlux: unknown fluxType."
                end select

                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

             end select ! model

             ! velocity
             select case( fluxType )

             case( "central" )
                vF = var(i,j+1,k,3)
                vB = var(i,j,k,3)

                gRhoV = 0.25*(vB+vF)**2

             case( "upwind", "ILES" )
                vF = vTilde(i,j+1,k,2,0)
                vB = vTilde(i,j,k,2,1)
                vBarF = vBar(i,j+1,k)
                vBarB = vBar(i,j,k)
                vSurf = 0.5*(vB + vF)

                select case( fluxType ) 

                case( "upwind" )
                   gRhoV = flux_muscl(vSurf,vB,vF)
                case( "ILES" )
                   gRhoV = flux_aldm(vB,vF,vSurf,vB,vF,vBarB,vBarF,sigmaY)
                case default
                   stop "momentumFlux: unknown fluxType."
                end select

             case default
                stop "momentumFlux: unknown fluxType."
             end select

             !                call absDiff(vBarF-vBarB,delta)
             !                gRhoV = 0.25*(vB+vF)**2 - sigmaY*(vF-vB)*delta

             flux(i,j,k,2,3) = rhoEdge * gRhoV
          end do
       end do
    end do

    !-------------------------------------------------------------------------

    ! vertical flux hRhoV
    do k = 0,nz
       do j = 0,ny
          do i = 1,nx

             ! density at flux point
             select case( model ) 

             case( "Boussinesq" )
                rhoEdge = rho00

             case default

                select case( fluxType ) 

                case( "central" )
                   ! density interpolation consistent with conti eq
                   rhoEdge = 0.25*(var(i,j,k,1) &
                        & + var(i,j,k+1,1) &
                        & + var(i,j+1,k,1) &
                        & + var(i,j+1,k+1,1) )

                case( "upwind", "ILES" )
                   rhoEdge = 0.25*(rhoTilde(i,j,k,3,1) &
                        & + rhoTilde(i,j,k+1,3,0) &
                        & + rhoTilde(i,j+1,k,3,1) &
                        & + rhoTilde(i,j+1,k+1,3,0) )
                case default 
                   stop "momentumFlux: unknown fluxType."
                end select

                ! comment: for CDS rhoEdge should add rhoStrat for each 
                ! var(...1) individually for 100% consistency with conti eq.
                
                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStratTilde(k)

             end select ! model

             ! velocity
             select case( fluxType )

             case( "central" )
                vU = var(i,j,k+1,3)
                vD = var(i,j,k,3)
                wF = var(i,j+1,k,4)
                wB = var(i,j,k,4)

                hRhoV = 0.25*(vD+vU)*(wB+wF)

             case( "upwind", "ILES" )
                vU = vTilde(i,j,k+1,3,0)
                vD = vTilde(i,j,k,3,1)
                wF = wTilde(i,j+1,k,2,0)
                wB = wTilde(i,j,k,2,1)
                vBarU = vBar(i,j,k+1)
                vBarD = vBar(i,j,k)
                wSurf = 0.5*(wB + wF)

                select case( fluxType ) 

                case( "upwind" )
                   hRhoV = flux_muscl(wSurf,vD,vU)
                case( "ILES" )
                   hRhoV = flux_aldm(vD,vU,wSurf,vD,vU,vBarD,vBarU,sigmaY)
                case default
                   stop "momentumFlux: unknown fluxType."
                end select

             case default
                stop "momentumFlux: unknown fluxType."
             end select

             !                call absDiff(wBarF-wBarB,delta)
             !                hRhoV = 0.25*(vD+vU)*(wB+wF) - sigmaY*(wF-wB)*delta

             flux(i,j,k,3,3) = rhoEdge * hRhoV
          end do
       end do
    end do


    !------------------------------
    !     flux for rho*w
    !------------------------------

    ! flux fRhoW
    do k = 0,nz
       do j = 1,ny
          do i = 0,nx

             ! density at flux point
             select case( model ) 

             case( "Boussinesq" )
                rhoEdge = rho00

             case default

                select case( fluxType ) 

                case( "central" )

                   if( fluctuationMode) then
                      rhoEdge = 0.25*(var(i,j,k,1) &
                           & + var(i+1,j,k,1) &
                           & + var(i,j,k+1,1) &
                           & + var(i+1,j,k+1,1) ) &
                           & + 0.5*(rhoStrat(k) + rhoStrat(k+1))
                   else
                      rhoEdge = 0.25*(var(i,j,k,1) &
                           & + var(i+1,j,k,1) &
                           & + var(i,j,k+1,1) &
                           & + var(i+1,j,k+1,1) )
                   end if

                case( "upwind", "ILES" )

                   if( fluctuationMode ) then
                      rhoEdge = 0.25*(rhoTilde(i,j,k,1,1) &
                           & + rhoTilde(i+1,j,k,1,0) &
                           & + rhoTilde(i,j,k+1,1,1) &
                           & + rhoTilde(i+1,j,k+1,1,0) ) &
                           & + 0.5*(rhoStrat(k) + rhoStrat(k+1))
                   else
                      rhoEdge = 0.25*(rhoTilde(i,j,k,1,1) &
                           & + rhoTilde(i+1,j,k,1,0) &
                           & + rhoTilde(i,j,k+1,1,1) &
                           & + rhoTilde(i+1,j,k+1,1,0) )
                   end if



                case default 
                   stop "momentumFlux: unknown fluxType."
                end select

             end select ! model

             ! velocity
             select case( fluxType )

             case( "central" )
                wR = var(i+1,j,k,4)
                wL = var(i,j,k,4)
                uU = var(i,j,k+1,2)
                uD = var(i,j,k,2)

                fRhoW = 0.25*(wL+wR)*(uD+uU)

             case( "upwind", "ILES" )
                wR = wTilde(i+1,j,k,1,0)
                wL = wTilde(i,j,k,1,1)
                uU = uTilde(i,j,k+1,3,0)
                uD = uTilde(i,j,k,3,1)
                wBarR = wBar(i+1,j,k)
                wBarL = wBar(i,j,k)
                uSurf = 0.5*(uD + uU)

                select case( fluxType ) 

                case( "upwind" )
                   fRhoW = flux_muscl(uSurf,wL,wR)
                case( "ILES" )
                   fRhoW = flux_aldm(wL,wR,uSurf,wL,wR,wBarR,wBarL,sigmaZ)
                case default
                   stop "momentumFlux: unknown fluxType."
                end select

             case default
                stop "momentumFlux: unknown fluxType."
             end select

             !                call absDiff(uBarU-uBarD,delta)
             !                fRhoW = 0.25*(wL+wR)*(uD+uU) - sigmaZ*(uU-uD)*delta

             flux(i,j,k,1,4) = rhoEdge * fRhoW
          end do
       end do
    end do

    !-------------------------------------------------------------------

    ! flux gRhoW
    do k = 0,nz
       do j = 0,ny
          do i = 1,nx

             ! density at flux point
             select case( model ) 

             case( "Boussinesq" )
                rhoEdge = rho00

             case default

                select case( fluxType ) 

                case( "central" )
                   
                   if( fluctuationMode ) then
                      rhoEdge = 0.25*( var(i,j,k,1) &
                           & + var(i,j+1,k,1) &
                           & + var(i,j,k+1,1) &
                           & + var(i,j+1,k+1,1) )&
                           & + 0.5*(rhoStrat(k) + rhoStrat(k+1))
                   else
                      rhoEdge = 0.25*( var(i,j,k,1) &
                           & + var(i,j+1,k,1) &
                           & + var(i,j,k+1,1) &
                           & + var(i,j+1,k+1,1) )
                   end if
                   
                case( "upwind", "ILES" )
                   if( fluctuationMode ) then
                      rhoEdge = 0.25*( rhoTilde(i,j,k,2,1) &
                           & + rhoTilde(i,j+1,k,2,0) &
                           & + rhoTilde(i,j,k+1,2,1) &
                           & + rhoTilde(i,j+1,k+1,2,0) ) &
                           & + 0.5*(rhoStrat(k) + rhoStrat(k+1))
                   else
                      rhoEdge = 0.25*( rhoTilde(i,j,k,2,1) &
                           & + rhoTilde(i,j+1,k,2,0) &
                           & + rhoTilde(i,j,k+1,2,1) &
                           & + rhoTilde(i,j+1,k+1,2,0) )
                   end if
                case default 
                   stop "momentumFlux: unknown fluxType."
                end select

             end select ! model

             ! velocity
             select case( fluxType )

             case( "central" )
                wF = var(i,j+1,k,4)
                wB = var(i,j,k,4)
                vU = var(i,j,k+1,3)
                vD = var(i,j,k,3)

                gRhoW = 0.25*(wB+wF)*(vD+vU)

             case( "upwind", "ILES" )
                wF = wTilde(i,j+1,k,2,0)
                wB = wTilde(i,j,k,2,1)
                vU = vTilde(i,j,k+1,3,0)
                vD = vTilde(i,j,k,3,1)
                wBarF = wBar(i,j+1,k)
                wBarB = wBar(i,j,k)
                vSurf = 0.5*(vU + vD)

                select case( fluxType ) 

                case( "upwind" )
                   gRhoW = flux_muscl(vSurf,wB,wF)
                case( "ILES" )
                   gRhoW = flux_aldm(wB,wF,vSurf,wB,wF,wBarB,wBarF,sigmaZ)
                case default
                   stop "momentumFlux: unknown fluxType."
                end select

             case default
                stop "momentumFlux: unknown fluxType."
             end select

             !                call absDiff(vBarU-vBarD, delta)
             !                gRhoW = 0.25*(wB+wF)*(vD+vU) - sigmaZ*(vU-vD)*delta

             flux(i,j,k,2,4) = rhoEdge * gRhoW
          end do
       end do
    end do

    !-----------------------------------------------------------------------

    ! flux hRhoW
    do k = -1,nz
       do j = 1,ny
          do i = 1,nx

             ! density at flux point
             select case( model ) 

             case( "Boussinesq" )
                rhoEdge = rho00

             case default

                select case( fluxType ) 

                case( "central" )

                   if( fluctuationMode ) then
                      rhoEdge = 0.25*( var(i,j,k,1) + rhoStrat(k) &
                           & + var(i,j,k+1,1) + rhoStrat(k+1) &
                           & + var(i,j,k+1,1) + rhoStrat(k+1) &
                           & + var(i,j,k+2,1) + rhoStrat(k+2) )
                      
                   else
                      rhoEdge = 0.25*( var(i,j,k,1) &
                           & + var(i,j,k+1,1) &
                           & + var(i,j,k+1,1) &
                           & + var(i,j,k+2,1) )
                   end if

                case( "upwind", "ILES" )

                   if( fluctuationMode ) then
                      rhoEdge = 0.25*( rhoTilde(i,j,k,3,1) &
                           & + rhoTilde(i,j,k+1,3,0) &
                           & + rhoTilde(i,j,k+1,3,1) &
                           & + rhoTilde(i,j,k+2,3,0) ) &
                           & + 0.5*(rhoStratTilde(k) + rhoStratTilde(k+1))
                   else
                      rhoEdge = 0.25*( rhoTilde(i,j,k,3,1) &
                           & + rhoTilde(i,j,k+1,3,0) &
                           & + rhoTilde(i,j,k+1,3,1) &
                           & + rhoTilde(i,j,k+2,3,0) )
                   end if
                   
                case default 
                   stop "momentumFlux: unknown fluxType."
                end select

             end select ! model

             ! velocity
             select case( fluxType )

             case( "central" )
                wU = var(i,j,k+1,4)
                wD = var(i,j,k,4)

                hRhoW = 0.25*(wD+wU)**2

             case( "upwind", "ILES" )
                wU = wTilde(i,j,k+1,3,0)
                wD = wTilde(i,j,k,3,1)
                wBarU = wBar(i,j,k+1)
                wBarD = wBar(i,j,k)
                wSurf = 0.5*(wD + wU)

                select case( fluxType ) 

                case( "upwind" )
                   hRhoW = flux_muscl(wSurf,wD,wU)
                case( "ILES" )
                   hRhoW = flux_aldm(wD,wU,wSurf,wD,wU,wBarD,wBarU,sigmaZ)
                case default
                   stop "momentumFlux: unknown fluxType."
                end select

             case default
                stop "momentumFlux: unknown fluxType."
             end select

             flux(i,j,k,3,4) = rhoEdge * hRhoW

          end do
       end do
    end do


    
    
    !-------------------------------------------------------------------
    !                          Viscous Fluxes
    !-------------------------------------------------------------------
    
    
    ! prevent calculation for inviscid flows
!    if( ReInv == 0.0 ) return                ! I choose to still do the calculation for inviscid flows (modified by Junhong Wei, 20160804)


    select case ( model ) 

    case( "Boussinesq" ) 

       !UA I do not quite get why in the Boussinesq case the stress tensor his
       !handled so strangely. See comments below ...

       !------------------------------
       !     flux for rho*u
       !------------------------------

       ! horizontal flux fRhoU
       do k = 1,nz
          do j = 1,ny
             do i = -1,nx

                du_dx = ( var(i+1,j,k,2) - var(i,j,k,2) )/dx     ! du/dx at i+1/2
!                fRhoU_visc = ReInv * du_dx         ! revised by JW (20160825)

                ! modified by Junhong Wei (20160803) --- starting line
                ! This part is revised again by JW (20160825)
                if(DySmaScheme)then

                !UA here not rather replace du_dx by du_dx + du_dx?
                ! Replied by JW (20160825): It is revised again for here and the rest of the code in "Boussinesq" background.

!                fRhoU_visc = ( ReInv + ( rho00*var(i+1,j,k,7) ) ) * du_dx
!                fRhoU_visc = ( ReInv + ( rho00*var(i+1,j,k,7) ) ) * ( du_dx + du_dx )
                fRhoU_visc = ( ReInv + ( rho00*var(i,j,k,7) ) ) * ( du_dx + du_dx )

                else

                fRhoU_visc = ReInv * ( du_dx + du_dx )

                end if
                ! modified by Junhong Wei (20160803) --- finishing line

                flux(i,j,k,1,2) = flux(i,j,k,1,2) - fRhoU_visc
             end do
          end do
       end do

       ! horizontal flux gRhoU
       do k = 1,nz
          do j = 0,ny
             do i = 0,nx
                du_dy = ( var(i,j+1,k,2) - var(i,j,k,2) )/dy     ! du/dy at j+1/2
                dv_dx = ( var(i+1,j,k,3) - var(i,j,k,3) )/dx     ! dv/dx     ! revised by JW (20160825)
!                gRhoU_visc = ReInv * du_dy          ! revised by JW (20160825)

                ! modified by Junhong Wei (20160803) --- starting line
                if(DySmaScheme)then

                !UA here not rather replace du_dy by du_dy + dv_dx (as in case
                !pseudo-incompressible? Same in all further instances below ...
!                gRhoU_visc = ( ReInv + ( rho00*var(i,j+1,k,7) ) ) * ( du_dy + dv_dx ) ! revised by JW (20160825)
                gRhoU_visc = ( ReInv + ( rho00*var(i,j,k,7) ) ) * ( du_dy + dv_dx ) ! revised by JW (20160825)

                else

                gRhoU_visc = ReInv * ( du_dy + dv_dx )          ! revised by JW (20160825)

                end if
                ! modified by Junhong Wei (20160803) --- finishing line

                flux(i,j,k,2,2) = flux(i,j,k,2,2) - gRhoU_visc
             end do
          end do
       end do

       ! vertical flux hRhoU
       do k = 0,nz
          do j = 1,ny
             do i = 0,nx
                du_dz = ( var(i,j,k+1,2) - var(i,j,k,2) )/dz       ! du/dz  at k+1/2
                dw_dx = ( var(i+1,j,k,4) - var(i,j,k,4) )/dx       ! dw/dx        ! revised by JW (20160825)
!                hRhoU_visc = ReInv * du_dz           ! revised by JW (20160825)

                ! modified by Junhong Wei (20160803) --- starting line
                if(DySmaScheme)then

                !UA du_dz + dw_dx?
!                hRhoU_visc = ( ReInv + ( rho00*var(i,j,k+1,7) ) ) * ( du_dz + dw_dx )  ! revised by JW (20160825)
                hRhoU_visc = ( ReInv + ( rho00*var(i,j,k,7) ) ) * ( du_dz + dw_dx )  ! revised by JW (20160825)

                else

                hRhoU_visc = ReInv * ( du_dz + dw_dx )  ! revised by JW (20160825)

                end if
                ! modified by Junhong Wei (20160803) --- finishing line

                flux(i,j,k,3,2) =  flux(i,j,k,3,2) - hRhoU_visc          
             end do
          end do
       end do

       !------------------------------
       !     flux for rho*v
       !------------------------------

       ! horizontal flux fRhoV
       do k = 1,nz
          do j = 0,ny
             do i = 0,nx
                dv_dx = ( var(i+1,j,k,3) - var(i,j,k,3) )/dx    ! dv/dx at i+1/2
                du_dy = ( var(i,j+1,k,2) - var(i,j,k,2) )/dy      ! du/dy       ! revised by JW (20160825)
!                fRhoV_visc = ReInv * dv_dx        ! revised by JW (20160825)

                ! modified by Junhong Wei (20160803) --- starting line
                if(DySmaScheme)then

                !UA dv_dx + du_dy?
!                fRhoV_visc = ( ReInv + ( rho00*var(i+1,j,k,7) ) ) * ( dv_dx + du_dy )     ! revised by JW (20160825)
                fRhoV_visc = ( ReInv + ( rho00*var(i,j,k,7) ) ) * ( dv_dx + du_dy )     ! revised by JW (20160825)

                else

                fRhoV_visc = ReInv * ( dv_dx + du_dy )     ! revised by JW (20160825)

                end if
                ! modified by Junhong Wei (20160803) --- finishing line

                flux(i,j,k,1,3) = flux(i,j,k,1,3) - fRhoV_visc
             end do
          end do
       end do

       ! horizontal flux gRhoV
       do k = 1,nz
          do j = -1,ny
             do i = 1,nx
                dv_dy = ( var(i,j+1,k,3) - var(i,j,k,3) )/dy    ! dv/dy at j+1/2
!                gRhoV_visc = ReInv * dv_dy       ! revised by JW (20160825)

                ! modified by Junhong Wei (20160803) --- starting line
                if(DySmaScheme)then

                !UA dv_dy + dv_dy?
!                gRhoV_visc = ( ReInv + ( rho00*var(i,j+1,k,7) ) ) * ( dv_dy + dv_dy )     ! revised by JW (20160825)
                gRhoV_visc = ( ReInv + ( rho00*var(i,j,k,7) ) ) * ( dv_dy + dv_dy )     ! revised by JW (20160825)

                else

                gRhoV_visc = ReInv * ( dv_dy + dv_dy )      ! revised by JW (20160825)

                end if
                ! modified by Junhong Wei (20160803) --- finishing line

                flux(i,j,k,2,3) = flux(i,j,k,2,3) - gRhoV_visc             
             end do
          end do
       end do

       ! vertical flux hRhoV
       do k = 0,nz
          do j = 0,ny
             do i = 1,nx
                dv_dz = ( var(i,j,k+1,3) - var(i,j,k,3) )/dz    ! dv/dz at k+1/2
                dw_dy = ( var(i,j+1,k,4) - var(i,j,k,4) )/dy    ! dw/dy            ! revised by JW (20160825)
!                hRhoV_visc = ReInv * dv_dz                 ! revised by JW (20160825)

                ! modified by Junhong Wei (20160803) --- starting line
                if(DySmaScheme)then

                !UA dv_dz + dw_dy?
!                hRhoV_visc = ( ReInv + ( rho00*var(i,j,k+1,7) ) ) * ( dv_dz + dw_dy )         ! revised by JW (20160825)
                hRhoV_visc = ( ReInv + ( rho00*var(i,j,k,7) ) ) * ( dv_dz + dw_dy )         ! revised by JW (20160825)

                else

                hRhoV_visc = ReInv * ( dv_dz + dw_dy )         ! revised by JW (20160825)

                end if
                ! modified by Junhong Wei (20160803) --- finishing line

                flux(i,j,k,3,3) = flux(i,j,k,3,3) - hRhoV_visc
             end do
          end do
       end do


       !------------------------------
       !     flux for rho*w
       !------------------------------

       ! horizontal flux fRhoW
       do k = 0,nz
          do j = 1,ny
             do i = 0,nx
                dw_dx = ( var(i+1,j,k,4) - var(i,j,k,4) )/dx   ! dw/dx at i+1/2
                du_dz = ( var(i,j,k+1,2) - var(i,j,k,2) )/dz   ! du/dz             ! revised by JW (20160825)
!                fRhoW_visc = ReInv * dw_dx           ! revised by JW (20160825)

                ! modified by Junhong Wei (20160803) --- starting line
                if(DySmaScheme)then

                !UA dw_dx + du_dz?
!                fRhoW_visc = ( ReInv + ( rho00*var(i+1,j,k,7) ) ) * ( dw_dx + du_dz )       ! revised by JW (20160825)
                fRhoW_visc = ( ReInv + ( rho00*var(i,j,k,7) ) ) * ( dw_dx + du_dz )       ! revised by JW (20160825)

                else

                fRhoW_visc = ReInv * ( dw_dx + du_dz )        ! revised by JW (20160825)

                end if
                ! modified by Junhong Wei (20160803) --- finishing line

                flux(i,j,k,1,4) = flux(i,j,k,1,4) - fRhoW_visc
             end do
          end do
       end do

       ! horizontal flux gRhoW
       do k = 0,nz
          do j = 0,ny
             do i = 1,nx
                dw_dy = ( var(i,j+1,k,4) - var(i,j,k,4) )/dy   ! dw/dy at j+1/2
                dv_dz = ( var(i,j,k+1,2) - var(i,j,k,2) )/dz   ! dv/dz                    ! revised by JW (20160825)
!                gRhoW_visc = ReInv * dw_dy                                               ! revised by JW (20160825)

                ! modified by Junhong Wei (20160803) --- starting line
                if(DySmaScheme)then

                !UA dw_dy + dv_dz?
!                gRhoW_visc = ( ReInv + ( rho00*var(i,j+1,k,7) ) ) * ( dw_dy + dv_dz )     ! revised by JW (20160825)
                gRhoW_visc = ( ReInv + ( rho00*var(i,j,k,7) ) ) * ( dw_dy + dv_dz )     ! revised by JW (20160825)

                else

                gRhoW_visc = ReInv * ( dw_dy + dv_dz )                                    ! revised by JW (20160825)

                end if
                ! modified by Junhong Wei (20160803) --- finishing line

                flux(i,j,k,2,4) =  flux(i,j,k,2,4) - gRhoW_visc
             end do
          end do
       end do

       ! vertical flux hRhoW
       do k = -1,nz
          do j = 1,ny
             do i = 1,nx
                dw_dz = ( var(i,j,k+1,4) - var(i,j,k,4) )/dz   ! dw/dz at k+1/2
!                hRhoW_visc = ReInv * dw_dz                                      ! revised by JW (20160825)

                ! modified by Junhong Wei (20160803) --- starting line
                if(DySmaScheme)then

                !UA dw_dz + dw_dz?
!                hRhoW_visc = ( ReInv + ( rho00*var(i,j,k+1,7) ) ) * ( dw_dz + dw_dz )        ! revised by JW (20160825)
                hRhoW_visc = ( ReInv + ( rho00*var(i,j,k,7) ) ) * ( dw_dz + dw_dz )        ! revised by JW (20160825)

                else

                hRhoW_visc = ReInv * ( dw_dz + dw_dz )                                       ! revised by JW (20160825)

                end if
                ! modified by Junhong Wei (20160803) --- finishing line

                flux(i,j,k,3,4) = flux(i,j,k,3,4) - hRhoW_visc
             end do
          end do
       end do


       !----------------------------------------------------------------------


    case( "pseudo_incompressible") 

       !------------------------------------
       !      calc div(u) for viscosity
       !------------------------------------

       do k = 1,nz
          do j = 0,ny+1                   ! j = 0 and ny+1 are div's at ghost cells
             do i = 0,nx+1                ! same for i
                uR = var(i,j,k,2)
                uL = var(i-1,j,k,2)
                vF = var(i,j,k,3)
                vB = var(i,j-1,k,3)
                wU = var(i,j,k,4)
                wD = var(i,j,k-1,4)

                ! divergence at (i,j,k):
                divU(i,j,k) = (uR-uL)/dx + (vF-vB)/dy + (wU-wD)/dz   
             end do
          end do
       end do
       divU(:,:,0) = divU(:,:,1)       ! set div's in bottom ghost cells
       divU(:,:,nz+1) = divU(:,:,nz)   ! set div's in top ghost cells


       !------------------------------
       !     flux for rho*u
       !------------------------------

       ! horizontal flux fRhoU
       do k = 1,nz
          do j = 1,ny
             do i = -1,nx
                ! du/dx at i+1/2
                du_dx = ( var(i+1,j,k,2) - var(i,j,k,2) )/dx 

                div = divU(i+1,j,k)

!               achatzb slight clean up of density dependent viscosity
!               fRhoU_visc = ReInv * ( du_dx + du_dx - 2./3.*div )

!               if(DySmaScheme)then
!                  if( fluctuationMode ) then
!                     fRhoU_visc &
!                     = (ReInv &
!                        + ((var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7))) &
!                       * ( du_dx + du_dx - 2./3.*div )
!                    else
!                     fRhoU_visc &
!                     = ( ReInv + ( var(i,j,k,1)*var(i,j,k,7) ) ) &
!                       * ( du_dx + du_dx - 2./3.*div )
!                  end if
!               end if

                coef_v = ReInv * rhoStrat(1)

                if(DySmaScheme)then
                   if( fluctuationMode ) then
                       coef_v &
                       = coef_v + (var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7)
                     else
                       coef_v &
                       = coef_v + var(i,j,k,1) * var(i,j,k,7)
                   end if
                end if

                fRhoU_visc = coef_v * ( du_dx + du_dx - 2./3.*div )
!               achatze

                flux(i,j,k,1,2) = flux(i,j,k,1,2) - fRhoU_visc
             end do
          end do
       end do

       ! horizontal flux gRhoU
       do k = 1,nz
          do j = 0,ny
             do i = 0,nx
                ! du/dy at j+1/2
                du_dy = ( var(i,j+1,k,2) - var(i,j,k,2) )/dy
                dv_dx = ( var(i+1,j,k,3) - var(i,j,k,3) )/dx     ! dv/dx

!               achatzb slight clean up of density dependent viscosity
!               gRhoU_visc = ReInv * ( du_dy + dv_dx )

!               if(DySmaScheme)then
!                  if( fluctuationMode ) then
!                     gRhoU_visc &
!                     = (ReInv &
!                        + ((var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7))) &
!                       * ( du_dy + dv_dx )
!                    else
!                     gRhoU_visc &
!                     = ( ReInv + ( var(i,j,k,1)*var(i,j,k,7) ) ) &
!                       * ( du_dy + dv_dx )
!                  end if
!               end if

                coef_v = ReInv * rhoStrat(1)

                if(DySmaScheme)then
                   if( fluctuationMode ) then
                       coef_v &
                       = coef_v + (var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7)
                     else
                       coef_v &
                       = coef_v + var(i,j,k,1) * var(i,j,k,7)
                   end if
                end if

                gRhoU_visc = coef_v * ( du_dy + dv_dx )
!               achatze

                flux(i,j,k,2,2) = flux(i,j,k,2,2) - gRhoU_visc
             end do
          end do
       end do

       ! vertical flux hRhoU
       do k = 0,nz
          do j = 1,ny
             do i = 0,nx
                ! du/dz at k+1/2
                du_dz = ( var(i,j,k+1,2) - var(i,j,k,2) )/dz
                dw_dx = ( var(i+1,j,k,4) - var(i,j,k,4) )/dx       ! dw/dx

!               achatzb slight clean up of density dependent viscosity
!               hRhoU_visc = ReInv * ( du_dz + dw_dx )

!               if(DySmaScheme)then
!                  if( fluctuationMode ) then
!                     hRhoU_visc &
!                     = ( ReInv &
!                        + ((var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7))) &
!                       * ( du_dz + dw_dx )
!                    else
!                     hRhoU_visc &
!                     = ( ReInv + ( var(i,j,k,1)*var(i,j,k,7) ) ) &
!                       * ( du_dz + dw_dx )
!                  end if
!               end if

                coef_v = ReInv * rhoStrat(1)

                if(DySmaScheme)then
                   if( fluctuationMode ) then
                       coef_v &
                       = coef_v + (var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7)
                     else
                       coef_v &
                       = coef_v + var(i,j,k,1) * var(i,j,k,7)
                   end if
                end if

                hRhoU_visc = coef_v * ( du_dz + dw_dx )
!               achatze

                flux(i,j,k,3,2) =  flux(i,j,k,3,2) - hRhoU_visc          
             end do
          end do
       end do

       !------------------------------
       !     flux for rho*v
       !------------------------------

       ! horizontal flux fRhoV
       do k = 1,nz
          do j = 0,ny
             do i = 0,nx
                ! dv/dx at i+1/2
                dv_dx = ( var(i+1,j,k,3) - var(i,j,k,3) )/dx
                du_dy = ( var(i,j+1,k,2) - var(i,j,k,2) )/dy      ! dv/dy

!               achatzb slight clean up of density dependent viscosity
!               fRhoV_visc = ReInv * ( dv_dx + du_dy )

!               if(DySmaScheme)then
!                  if( fluctuationMode ) then
!                     fRhoV_visc &
!                     = ( ReInv &
!                         + ((var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7))) &
!                       * ( dv_dx + du_dy )
!                    else
!                     fRhoV_visc &
!                     = ( ReInv + ( var(i,j,k,1)*var(i,j,k,7) ) ) &
!                       * ( dv_dx + du_dy )
!                  end if
!               end if

                coef_v = ReInv * rhoStrat(1)

                if(DySmaScheme)then
                   if( fluctuationMode ) then
                       coef_v &
                       = coef_v + (var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7)
                     else
                       coef_v &
                       = coef_v + var(i,j,k,1) * var(i,j,k,7)
                   end if
                end if

                fRhoV_visc = coef_v * ( dv_dx + du_dy )
!               achatze

                flux(i,j,k,1,3) = flux(i,j,k,1,3) - fRhoV_visc
             end do
          end do
       end do

       ! horizontal flux gRhoV
       do k = 1,nz
          do j = -1,ny
             do i = 1,nx
                ! dv/dy at j+1/2
                dv_dy = ( var(i,j+1,k,3) - var(i,j,k,3) )/dy

                div = divU(i,j+1,k)

!               achatzb slight clean up of density dependent viscosity
!               gRhoV_visc = ReInv * ( dv_dy + dv_dy - 2./3.*div )

!               if(DySmaScheme)then
!                  if( fluctuationMode ) then
!                     gRhoV_visc &
!                     = ( ReInv &
!                         + ((var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7))) &
!                       * ( dv_dy + dv_dy - 2./3.*div )
!                    else
!                     gRhoV_visc &
!                     = ( ReInv + ( var(i,j,k,1)*var(i,j,k,7) ) ) &
!                       * ( dv_dy + dv_dy - 2./3.*div )
!                  end if
!               end if

                coef_v = ReInv * rhoStrat(1)

                if(DySmaScheme)then
                   if( fluctuationMode ) then
                       coef_v &
                       = coef_v + (var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7)
                     else
                       coef_v &
                       = coef_v + var(i,j,k,1) * var(i,j,k,7)
                   end if
                end if

                gRhoV_visc = coef_v * ( dv_dy + dv_dy - 2./3.*div )
!               achatze

                flux(i,j,k,2,3) = flux(i,j,k,2,3) - gRhoV_visc             
             end do
          end do
       end do

       ! vertical flux hRhoV
       do k = 0,nz
          do j = 0,ny
             do i = 1,nx
                ! dv/dz at k+1/2
                dv_dz = ( var(i,j,k+1,3) - var(i,j,k,3) )/dz
                dw_dy = ( var(i,j+1,k,4) - var(i,j,k,4) )/dy    ! dw/dy

!               achatzb slight clean up of density dependent viscosity
!               hRhoV_visc = ReInv * ( dv_dz + dw_dy )

!               if(DySmaScheme)then
!                  if( fluctuationMode ) then
!                     hRhoV_visc &
!                     = ( ReInv &
!                         + ((var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7))) &
!                       * ( dv_dz + dw_dy )
!                    else
!                     hRhoV_visc &
!                     = ( ReInv + ( var(i,j,k,1)*var(i,j,k,7) ) ) &
!                     * ( dv_dz + dw_dy )
!                  end if
!               end if

                coef_v = ReInv * rhoStrat(1)

                if(DySmaScheme)then
                   if( fluctuationMode ) then
                       coef_v &
                       = coef_v + (var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7)
                     else
                       coef_v &
                       = coef_v + var(i,j,k,1) * var(i,j,k,7)
                   end if
                end if

                hRhoV_visc = coef_v * ( dv_dz + dw_dy )
!               achatze

                flux(i,j,k,3,3) = flux(i,j,k,3,3) - hRhoV_visc
             end do
          end do
       end do


       !------------------------------
       !     flux for rho*w
       !------------------------------

       ! horizontal flux fRhoW
       do k = 0,nz
          do j = 1,ny
             do i = 0,nx
                ! dw/dx at i     +1/2
                dw_dx = ( var(i+1,j,k,4) - var(i,j,k,4) )/dx
                du_dz = ( var(i,j,k+1,2) - var(i,j,k,2) )/dz   ! du/dz

!               achatzb slight clean up of density dependent vicosity
!               fRhoW_visc = ReInv * ( dw_dx + du_dz )

!               if(DySmaScheme)then
!                  if( fluctuationMode ) then
!                     fRhoW_visc &
!                     = ( ReInv &
!                         + ((var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7))) &
!                       * ( dw_dx + du_dz )
!                    else
!                     fRhoW_visc &
!                     = ( ReInv + ( var(i,j,k,1)*var(i,j,k,7) ) ) &
!                       * ( dw_dx + du_dz )
!                  end if
!               end if

                coef_v = ReInv * rhoStrat(1)

                if(DySmaScheme)then
                   if( fluctuationMode ) then
                       coef_v &
                       = coef_v + (var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7)
                     else
                       coef_v &
                       = coef_v + var(i,j,k,1) * var(i,j,k,7)
                   end if
                end if

                fRhoW_visc = coef_v * ( dw_dx + du_dz )
!               achatze

                flux(i,j,k,1,4) = flux(i,j,k,1,4) - fRhoW_visc
             end do
          end do
       end do

       ! horizontal flux gRhoW
       do k = 0,nz
          do j = 0,ny
             do i = 1,nx
                ! dw/dy at j+1/2
                dw_dy = ( var(i,j+1,k,4) - var(i,j,k,4) )/dy
                dv_dz = ( var(i,j,k+1,2) - var(i,j,k,2) )/dz   ! dv/dz

!               achatzb slight clean up of density dependent viscosity
!               gRhoW_visc = ReInv * ( dw_dy + dv_dz)

!               if(DySmaScheme)then
!                  if( fluctuationMode ) then
!                     gRhoW_visc &
!                     = ( ReInv &
!                         + ((var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7))) &
!                       * ( dw_dy + dv_dz)
!                    else
!                     gRhoW_visc &
!                     = ( ReInv + ( var(i,j,k,1)*var(i,j,k,7) ) ) &
!                       * ( dw_dy + dv_dz)
!                  end if
!               end if

                coef_v = ReInv * rhoStrat(1)

                if(DySmaScheme)then
                   if( fluctuationMode ) then
                       coef_v &
                       = coef_v + (var(i,j,k,1)+rhoStrat(k)) * var(i,j,k,7)
                     else
                       coef_v &
                       = coef_v + var(i,j,k,1) * var(i,j,k,7)
                   end if
                end if

                gRhoW_visc = coef_v * ( dw_dy + dv_dz)
!               achatze

                flux(i,j,k,2,4) =  flux(i,j,k,2,4) - gRhoW_visc
             end do
          end do
       end do

       ! vertical flux hRhoW
       do k = -1,nz
          do j = 1,ny
             do i = 1,nx
                ! dw/dz at k+1/2
                dw_dz = ( var(i,j,k+1,4) - var(i,j,k,4) )/dz

                div = divU(i,j,k+1)

!               slight clean up of density dependent viscosity
!               hRhoW_visc = ReInv * ( dw_dz + dw_dz - 2./3.*div )

!               if(DySmaScheme)then
!                  if( fluctuationMode ) then
!                     if(k == -1) then
!                        rhos = rhoStrat(0)
!                       else
!                        rhos = rhoStrat(k)
!                     end if

!                     hRhoW_visc &
!                     = &
!                     ( ReInv + ( (var(i,j,k,1)+rhos) * var(i,j,k,7) ) ) &
!                     * ( dw_dz + dw_dz - 2./3.*div )
!                    else
!                     hRhoW_visc &
!                     = ( ReInv &
!                         + (var(i,j,k,1)*var(i,j,k,7))) &
!                       * ( dw_dz + dw_dz - 2./3.*div )
!                  end if
!               end if

                coef_v = ReInv * rhoStrat(1)

                if(DySmaScheme)then
                   if( fluctuationMode ) then
                      if(k == -1) then
                         rhos = rhoStrat(0)
                        else
                         rhos = rhoStrat(k)
                      end if

                       coef_v &
                       = coef_v + (var(i,j,k,1)+rhos) * var(i,j,k,7)
                     else
                       coef_v &
                       = coef_v + var(i,j,k,1) * var(i,j,k,7)
                   end if
                end if

                hRhoW_visc = coef_v * ( dw_dz + dw_dz - 2./3.*div )
!               achatze

                flux(i,j,k,3,4) = flux(i,j,k,3,4) - hRhoW_visc
             end do
          end do
       end do

    case default
       stop "momentumFlux: unknown case model"
    end select



    if (verbose) print*,"fluxes.f90/momentumFlux: &
         & momentum fluxes fRhoU, fRhoV, fRhoW, &
         & gRhoU, gRhoV, gRhoW&
         & hRhoU, hRhoV, hRhoW calculated."

  end subroutine momentumFlux


  !---------------------------------------------------------------------------


  subroutine momentumSource (var,source)
    !---------------------------------------------------------------------
    ! computes the source terms in the momentum equation
    ! 1) div error correction: rhoU * div(u)
    !---------------------------------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var

    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: source

    integer :: i,j,k
    real :: uL,uR, uB,uF, uD,uU
    real :: vL,vR, vB,vF, vD,vU
    real :: wL,wR, wB,wF, wD,wU
    real :: uSurfL, vSurfF, wSurfU     ! velocities at cell surface
    real :: uSurfR, vSurfB, wSurfD     ! velocities at cell surface
    real :: divPu, u, v, w, theta
    real :: PstratU, PstratD, PstratC

    !--------------------------------------------------------
    !              Divergence correction for rhoU
    !--------------------------------------------------------

    do k = 1,nz
       do j = 1,ny
          do i = 0,nx


             select case( fluxType )

             case( "central" )
                
                ! no correction implemented
                return
                
             case( "upwind","ILES" )
                
                uR = uTilde(i+1,j,k,1,0)
                uL = uTilde(i,j,k,1,1) 
                uSurfR = 0.5*(uL + uR)
                               
                uR = uTilde(i,j,k,1,0)
                uL = uTilde(i-1,j,k,1,1) 
                uSurfL = 0.5*(uL + uR)
                
                vR = vTilde(i+1,j,k,1,0)
                vL = vTilde(i,j,k,1,1)
                vSurfF = 0.5*(vR + vL)

                vR = vTilde(i+1,j-1,k,1,0)
                vL = vTilde(i,j-1,k,1,1)
                vSurfB = 0.5*(vR + vL)
                
                wR = wTilde(i+1,j,k,1,0)
                wL = wTilde(i,j,k,1,1)
                wSurfU = 0.5*(wL + wR)
                
                wR = wTilde(i+1,j,k-1,1,0)
                wL = wTilde(i,j,k-1,1,1)
                wSurfD = 0.5*(wL + wR)

                PstratU = PstratTilde(k)
                PstratD = PstratTilde(k-1)
                PstratC = Pstrat(k)
                
                divPu = PstratC*( (uSurfR - uSurfL)/dx &
                     & + (vSurfF - vSurfB) / dy ) &
                     & + (PstratU*wSurfU - PstratD*wSurfD) / dz
                
                theta = thetaStrat(k)
                
                u   = var(i,j,k,2) 
                
                source(i,j,k,2) = u * divPu / theta
                
             case default
                stop "thetaFlux: unknown case fluxType"
             end select

             
          end do
       end do
    end do

    !--------------------------------------------------------
    !              Divergence correction for rhoV
    !--------------------------------------------------------

    do k = 1,nz
       do j = 0,ny
          do i = 1,nx


             select case( fluxType )

             case( "central" )
                
                ! no correction implemented
                return
                
             case( "upwind","ILES" )

                uF = uTilde(i,j+1,k,2,0)
                uB = uTilde(i,j  ,k,2,1)
                uSurfR = 0.5*(uB + uF)

                uF = uTilde(i-1,j+1,k,2,0)
                uB = uTilde(i-1,j  ,k,2,1)
                uSurfL = 0.5*(uB + uF)

                vF = vTilde(i,j+1,k,2,0)
                vB = vTilde(i,j,k,2,1)
                vSurfF = 0.5*(vB + vF)
                
                vF = vTilde(i,j,k,2,0)
                vB = vTilde(i,j-1,k,2,1)
                vSurfB = 0.5*(vB + vF)
                
                wF = wTilde(i,j+1,k,2,0)
                wB = wTilde(i,j,k,2,1)
                wSurfU = 0.5*(wB + wF)

                wF = wTilde(i,j+1,k-1,2,0)
                wB = wTilde(i,j,k-1,2,1)
                wSurfD = 0.5*(wB + wF)

                PstratU = PstratTilde(k)
                PstratD = PstratTilde(k-1)
                PstratC = Pstrat(k)
                
                divPu = PstratC*( (uSurfR - uSurfL)/dx &
                     & + (vSurfF - vSurfB)/dy ) &
                     & + (PstratU*wSurfU - PstratD*wSurfD)/dz

                theta = thetaStrat(k)
                
                v   = var(i,j,k,3) 
                
                source(i,j,k,3) = v * divPu / theta
                
             case default
                stop "thetaFlux: unknown case fluxType"
             end select

             
          end do
       end do
    end do


    !--------------------------------------------------------
    !               Divergence correction for rhoW
    !--------------------------------------------------------
    
    do k = 1,nz-1           !xxxx assume solid wall here
       do j = 1,ny
          do i = 0,nx


             select case( fluxType )

             case( "central" )
                
                ! no correction implemented
                return
                
             case( "upwind","ILES" )
                
                uU = uTilde(i,j,k+1,3,0)
                uD = uTilde(i,j,k,3,1)
                uSurfR = 0.5*(uD + uU)

                uU = uTilde(i-1,j,k+1,3,0)
                uD = uTilde(i-1,j,k,3,1)
                uSurfL = 0.5*(uD + uU)

                vU = vTilde(i,j,k+1,3,0)
                vD = vTilde(i,j,k,3,1)
                vSurfF = 0.5*(vU + vD)

                vU = vTilde(i,j-1,k+1,3,0)
                vD = vTilde(i,j-1,k,3,1)
                vSurfB = 0.5*(vU + vD)

                wU = wTilde(i,j,k+1,3,0)
                wD = wTilde(i,j,k,3,1)
                wSurfU = 0.5*(wD + wU)
                
                wU = wTilde(i,j,k,3,0)
                wD = wTilde(i,j,k-1,3,1)
                wSurfD = 0.5*(wD + wU)

                PstratU = Pstrat(k+1)
                PstratD = Pstrat(k)
                PstratC = PstratTilde(k)
                                
                divPu = PstratC*( (uSurfR - uSurfL)/dx &
                     & + (vSurfF - vSurfB)/dy ) &
                     & + (PstratU*wSurfU - PstratD*wSurfD) / dz
                
                theta = thetaStratTilde(k)
                
                w   = var(i,j,k,4) 
                
                source(i,j,k,4) = w * divPu / theta
                
             case default
                stop "thetaFlux: unknown case fluxType"
             end select

             
          end do
       end do
    end do


    
  end subroutine momentumSource


  !---------------------------------------------------------------------------


  subroutine init_fluxes
    !---------------------------------------------
    ! 1) set parameter for central or ILES flux
    ! 2) allocate flux module variables
    !---------------------------------------------

    ! local variables
    integer :: allocstat


    ! module variables 

    ! rhoBar
    allocate(rhoBar(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhoBar"

    ! rhoOld
    allocate( rhoOld(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz), stat=allocstat)
    if( allocstat /= 0) stop "init_fluxes: alloc of rhoOld failed"

    ! uBar
    allocate(uBar(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "init_fluxes: could not allocate uBar"

    ! vBar
    allocate(vBar(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate vBar"

    ! wBar
    allocate(wBar(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate wBar"

    ! thetaBar
    allocate(thetaBar(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate thetaBar"

    ! rhoTilde
    allocate(rhoTilde(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,1:3,0:1),&
         & stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhoTilde"

    ! uTilde
    allocate(uTilde(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,1:3,0:1),&
         & stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate uTilde"

    ! vTilde
    allocate(vTilde(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,1:3,0:1),&
         & stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate vTilde"

    ! wTilde
    allocate(wTilde(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,1:3,0:1),&
         & stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate wTilde"

    ! thetaTilde
    allocate(thetaTilde(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,1:3,0:1),&
         & stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate thetaTilde"



    if (verbose) print*,"init_fluxes: &
         &rhoBar, uBar, vBar, wBar, thetaBar, &
         & rhoTilde, uTilde, vTilde, wTilde allocated."
  end subroutine init_fluxes


  !---------------------------------------------------------------------------


  subroutine terminate_fluxes
    !-----------------------------------
    ! deallocates flux module variables
    !-----------------------------------

    ! local variables
    integer :: allocstat

    !---------------- deallocate variables -----------------------

    deallocate(rhoBar,stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate rhoBar"

    deallocate( rhoOld, stat=allocstat)
    if(allocstat /= 0) stop "terminate_fluxes: dealloc of rhoOld failed"

    deallocate(uBar,stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate uBar"

    deallocate(vBar,stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate vBar"

    deallocate(wBar,stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate wBar"

    deallocate(thetaBar,stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate thetaBar"

    deallocate(rhoTilde,stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate rhoTilde"

    deallocate(uTilde,stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate uTilde"

    deallocate(vTilde,stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate vTilde"

    deallocate(wTilde,stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate wTilde"

    deallocate(thetaTilde,stat=allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate thetaTilde"


  end subroutine terminate_fluxes


  ! ----------------------------------------------------

  subroutine absDiff(x, absX)
    !----------------------------------------------------------------------
    ! differentiable approx. of absolute value abs()-function for linFloit
    !----------------------------------------------------------------------

    ! in/out variables
    real :: x
    intent(in) :: x
    real :: absX
    intent(out) :: absX

    ! local vars
    real, parameter :: delta0 = 1.0e-18


    if (abs(x) > delta0) then
       absX = abs(x) 
    else
       absX = (x**2 + delta0**2)/2.0/delta0
    end if


  end subroutine absDiff


end module flux_module
