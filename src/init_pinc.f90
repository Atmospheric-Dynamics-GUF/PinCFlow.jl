module init_module


  use type_module
  use atmosphere_module


  implicit none

  private         ! private module

  public  :: initialise
  public  :: setup
  public  :: init_ParamStudy



contains

  subroutine init_paramStudy
    !--------------------------------
    ! read data for parameter study
    !--------------------------------
    
    ! open input file input.f90
    open (unit=10, file="input.f90", action="read", &
         form="formatted", status="old", position="rewind")
    
    ! read parameter study list
    read (unit=10, nml=parameterList)

    close( unit=10 ) 
  

  end subroutine init_paramStudy

 
! ------------------------------------------------------------------------


  subroutine setup (var,var0,flux,force,dRho,dMom,dTheta)
    !-----------------------------------------
    ! allocate var and flux / read input.f90
    !-----------------------------------------

    ! in/out variables 
    real,dimension(:,:,:,:), allocatable,intent(out) :: var, var0
    real,dimension(:,:,:,:,:), allocatable,intent(out) :: flux
    real,dimension(:,:,:,:), allocatable,intent(out) :: force
    real, dimension(:,:,:), allocatable :: dRho        ! RK-Update for rho
    real, dimension(:,:,:,:), allocatable :: dMom      ! ...rhoU,rhoV,rhoW
    real, dimension(:,:,:), allocatable :: dTheta       ! RK-Update for theta


    integer :: allocstat
    integer :: i, j, k, iVar


    ! constants
    pi = 4*atan(1.0)


    ! open input file input.f90
    open (unit=10, file="input.f90", action="read", &
         form="formatted", status="old", position="rewind")

    ! read grid info
    read (unit=10, nml=grid)
   
    ! total dimension of var fields with ghost cells
    nxx = nx + 2*nbx + 1
    nyy = ny + 2*nby + 1
    nzz = nz + 2*nbz + 1
    

    !---------------------------------------------------
    ! allocate x,y,z - cell centered coordinate fields
    !---------------------------------------------------
    
    allocate(x(-nbx:nx+nbx), stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate x"

    allocate(y(-nby:ny+nby), stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate y"

    allocate(z(-nbz:nz+nbz), stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate z"


    !-------------------------------------
    !      allocate variable fields
    !-------------------------------------


    read (unit=10, nml=variables)

    ! allocate var = (rho,u,v,w,pEx)
    allocate(var(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate var"

    ! allocate var0 = (rho,u,v,w,pEx)
    allocate(var0(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate var"

    ! allocate dRho
    allocate(dRho(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dRho."

    ! allocate dMom
    allocate(dMom(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,3),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dMom."
    
    ! allocate dTheta
    allocate(dTheta(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dTheta."

    ! allocate flux = (f,g,h / fRho, fRhoU, fRhoV, fRhoW, fTheta)
    allocate(flux(-1:nx,-1:ny,-1:nz,3,nVar),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate flux"

    ! allocate force 
    allocate(force(0:nx+1,0:ny+1,0:nz+1,3),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate force"
    

    !-------------------------------------
    !    read name lists from input.f90
    !-------------------------------------

    ! read model equation specifications
    read (unit=10, nml=modelList)

    ! read solver specifications
    read (unit=10, nml=solverList)
   
    ! read Poisson solver specifications
    read (unit=10, nml=poissonSolverList)

    ! read atmosphere specifications
    read (unit=10, nml=atmosphereList)

    ! read topography data
    read (unit=10, nml=topographyList)

    ! read boundary info
    read (unit=10, nml=boundaryList)
    
    ! read boundary info2
    read (unit=10, nml=boundaryList2)

    ! read output specifications
    read (unit=10, nml=outputList)
    
    ! read programme debug parameters
    read (unit=10, nml=debuggingList)
    
    ! read wkb infos
    read (unit=10, nml=wkbList)

    ! read test case name
    read (unit=10, nml=testCaseList)


    ! close input file pinc.f
    close (unit=10)



    !---------------------------------------
    !        Model equation settings
    !---------------------------------------
    
    select case( model ) 
       
    case( "Boussinesq" )
       
       updateMass = .false.
       predictMomentum = .true.
       correctMomentum = .true.
       updateTheta = .true.
       
       ! overwrite unsuitable input settings
       topography = .false. 
       spongeLayer = .false.
       pressureScaling = .false.
       
       ! overwrite atmosphere to uniform
       background = "uniform"

       ! never offset theta 
       thetaOffset = .false.

       
       
    case( "pseudo_incompressible" ) 
       
       updateMass = .true.
       predictMomentum = .true.
       correctMomentum = .true.
       updateTheta = .false.
       
       ! overwrite unsuitable input settings
       if( zBoundary == "periodic" ) then
          print*,"WARNING: zBoundary periodic not possible. Reset to solid_wall!"
          zBoundary = "solid_wall"      
       end if

    case( "WKB" )
      
       updateMass = .false.
       predictMomentum = .false.
       correctMomentum = .false.
       updateTheta = .false.

       if( zBoundary == "periodic" ) stop "WKB not ready for zBoundary = periodic."
       
    case default
       print*,"model = ", model
       stop "initialize: Unknown model" 
    end select

    
    
  end subroutine setup


! ------------------------------------------------------------------------


  subroutine initialise (var)
    !------------------
    ! setup test cases
    !------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
     ! in/out variables
 !   real, dimension(:,:,:,:),allocatable,intent(inout) :: var

    ! local variables
    integer :: i,j,k

    ! greshoVortex
    real :: xu, yu, zu, xv, yv, xw, zw, delX, delY
    real :: x0, y0, r, uVortex
    real :: uPhi, p, rho0
    real :: pInf, u0, r0

    ! densityDisc
    real :: v0, w0, dens, rhoDisc, z0
    
    ! Robert's hot and cold bubble
  real :: dTheta1, a1, sigma1, xCenter1, zCenter1
  real :: dTheta2, a2, sigma2, xCenter2, zCenter2


    
    ! hotBubble
    real :: x_dim, y_dim, z_dim, delZ
    real :: dTheta, dTheta0, theta, rho
!    real :: xCenter_dim, zCenter_dim    ! position of thermal
!    real :: xRadius_dim, zRadius_dim    ! radii of thermal (ellipse)
    real :: dTheta_dim! , dTheta0_dim

    ! wavepacket: all quantities are scaled
    real :: kk, mm, kTot                 ! vertical, zonal, total wave number
    real :: omi, omi2                    ! intrinsic frequency, squared 
    real :: bAmp, uAmp, wAmp, pAmp       ! amplitudes for buoyancy, u, w, Exner function
    real :: delx2, delz2                 ! squared distance from center
    real :: Gauss                        ! Gaussian distribution value
    real :: sigma                        ! width of Gaussian distribution
    real :: xCenter, zCenter             ! center of Gaussian distribution
    real :: phi                          ! local phase 
    real :: u,v,w,b                        ! buoyancy, zonal + vertical velocity
    real :: lambdaX, lambdaZ             ! zonal and vertical wave length

    ! monochromatic wave
    real :: th, f, f2
    real :: lambda
    real :: amp, cot, al
    real :: uRot, vRot, wRot
    
    
    ! open input file input.f90
    open (unit=10, file="input.f90", action="read", &
         form="formatted", status="old", position="rewind")

    
    !---------------------------
    !      General settings
    !---------------------------
    
    select case( model ) 
       
    case( "Boussinesq" )
       var(:,:,:,1) = rho00          ! always constant background density

       ! calc vertical
       vertical(1) = -cos(vert_theta*pi/180.)*cos(vert_alpha*pi/180.)
       vertical(2) = +cos(vert_theta*pi/180.)*sin(vert_alpha*pi/180.)
       vertical(3) = +sin(vert_theta*pi/180.)

    case( "pseudo_incompressible" )
       
       ! set vertical always paralle to z-axis
       vertical = (/0.0, 0.0, 1.0 /)
       
    case default
       stop "initialize: unknown case model."
    end select
    
    
    
    !---------------------------------------
    !        Test case settings
    !---------------------------------------
    
    select case (testCase)

       !-------------------------------------
       !             Boussinesq only
       !-------------------------------------

    case( "uniform_theta" ) 
       
       ! zero atmospheric background flow
       var(:,:,:,2) = backgroundFlow_dim(1) / uRef
       var(:,:,:,3) = backgroundFlow_dim(2) / uRef
       var(:,:,:,4) = backgroundFlow_dim(3) / uRef


       ! constant pressure variable pi' 
       var(:,:,:,5) = 0.0

       ! uniform pot temp perturbation
       var(:,:,:,6) = 0.0 / thetaRef


    case( "monochromeWave" ) 

       

!!$       ! read test case input data
!!$       read (unit=10, nml=monochromeWave)
!!$       
!!$       ! gather quantities
!!$       th = vert_theta * pi/180.0
!!$       al = vert_alpha * pi/180.0
!!$       f = f_Coriolis_dim * tRef
!!$       f2 = f**2
!!$       omi = -sqrt(N2*cos(th)**2 + f2*sin(th)**2)
!!$       
!!$       
!!$       ! wave vector
!!$       lambda = lambda_dim / lRef
!!$       kTot = 2.0*pi / lambda
!!$       
!!$       kk = kTot*cos(th)
!!$       mm = kTot*sin(th)
!!$
!!$       ! amplitudes in non-rotated frame
!!$       amp = amplitudeFactor*(-omi/kk)
!!$       cot = 1.0/tan(th)
       
!!$       do k = 1,nz
!!$          j = 1
!!$          do i = 1,nx
!!$             
!!$             ! phase
!!$             phi = -kk*x(i) - mm*z(k)
!!$             
!!$             ! amplitudes
!!$             u = amp * cos(phi)
!!$             v = amp*f/omi * sin(phi)
!!$             w = amp*cot * cos(phi)
!!$             b = amp*N2/omi*cot * sin(phi)
!!$
!!$            ! project on rotated coordinate system
!!$             
!!$             wRot = dot_product( (/u,v,w/), vertical) 
!!$             vRot = v
!!$             uRot = sqrt(u**2+w**2 - wRot**2)
!!$             theta = Fr2 * theta00 * b
!!$             
!!$             ! assign to var field
!!$             var(i,j,k,2:4) = (/ uRot,vRot,wRot /)
!!$             var(i,j,k,6) = theta
!!$
!!$          end do
!!$       end do
       
       
       ! open input file (double precision)
       open(80,file='bouss_init.datinp',form='unformatted',&
            & access='direct',recl=2*nx*nz, err=1078)
       !
       goto 1079
1078   stop "initialize: error reading data file."
1079   continue

       read(80,rec=1)&
       &           ((var(i,1,k,2),i=1,nx),k=1,nz)
       read(80,rec=2)&
       &           ((var(i,1,k,3),i=1,nx),k=1,nz)
       read(80,rec=3)&
       &           ((var(i,1,k,4),i=1,nx),k=1,nz)
       read(80,rec=4)&
       &           ((var(i,1,k,6),i=1,nx),k=1,nz)

       ! make non-dimensional
       var(:,:,:,2:4) = var(:,:,:,2:4)/uRef


       ! calc theta
       var(:,:,:,6) = var(:,:,:,6)/lRef*tRef**2
       var(:,:,:,6) = Fr2 * theta00 * var(:,:,:,6)

       ! set periodic BC
       var(0,:,:,2) = var(nx,:,:,2)
       var(:,0,:,3) = var(:,ny,:,3)
       var(:,:,0,4) = var(:,:,nz,4)

       close(80)


!!$       !-----------------------------
!!$       !  Interpolate to cell faces
!!$       !-----------------------------
!!$
!!$       ! copy u values to ghost cells
!!$       var(0,:,:,2) = var(nx,:,:,2)
!!$       var(nx+1,:,:,2) = var(1,:,:,2)
!!$       
!!$       ! average zonal velocities to cell face...
!!$       do i = 0,nx
!!$          var(i,:,:,2) = 0.5*( var(i,:,:,2) + var(i+1,:,:,2) )
!!$       end do
!!$
!!$       ! copy v values to ghost cells
!!$       var(:,0,:,3) = var(:,ny,:,3)
!!$       var(:,ny+1,:,3) = var(:,1,:,3)
!!$       
!!$       ! average meridional velocities to cell face...
!!$       do j = 0,ny
!!$          var(:,j,:,3) = 0.5*( var(:,j,:,3) + var(:,j+1,:,3) )
!!$       end do
!!$
!!$       ! copy w values to ghost cells
!!$       var(:,:,0,4) = var(:,:,nz,4)
!!$       var(:,:,nz+1,4) = var(:,:,1,4)
!!$       
!!$       ! average vertical velocities to cell face...
!!$       do k = 0,nz
!!$          var(:,:,k,4) = 0.5*( var(:,:,k,4) + var(:,:,k+1,4) )
!!$       end do
!!$

       
       
       !----------------------------------------
       !          WKB theory 
       !----------------------------------------

    case( 'wavePacket' )              ! 1D/2D wave packet

       ! read test case input data
       read (unit=10, nml=wavePacket)

       ! scale input data
       lambdaX = lambdaX_dim/lRef     ! non-dim zonal wave length
       lambdaZ = lambdaZ_dim/lRef     !         vert. wave length
       
       xCenter = xCenter_dim/lRef     ! scaled position of wave packtet
       zCenter = zCenter_dim/lRef 

       sigma = sigma_dim/lRef         ! sigma width of Gaussian distribution
       
       ! wave numbers
       mm = -2.0*pi/lambdaZ
       kk = 2.0*pi/lambdaX 
       kTot = sqrt(mm**2 + kk**2)
       
       ! intrinsic frequency
       omi = sqrt(N2)*kk/kTot
       omi2 = omi**2
       
       ! amplitude coefficients
       bAmp = amplitudeFactor * N2/mm                  ! buoyancy
       uAmp = mm/kk * omi/N2 * bAmp
       wAmp = omi/N2 * bAmp
       pAmp = kappa*Ma2 * mm/kk**2 * omi2/N2 * bAmp    ! Exner pressure


       do k = 1,nz
          j = 1
          do i = 1,nx
             
             ! Gaussian profile: 1D and 2D
             if( wavePacketDim == 1 ) then
                delx2 = 0.0
             else
                delx2 = (x(i)-xCenter)**2
             end if
             delz2 = (z(k)-zCenter)**2
             Gauss = exp(-(delx2 + delz2)/2./sigma**2) 

             ! common phase term: 1D and 2D
             phi = kk*x(i) + mm*z(k)
             
             ! buoyancy
             b = bAmp*Gauss*cos(phi)
             
             ! cell centered values
             rho = 1./(1.+Fr2*b) * rhoStrat(k) 
             u = uAmp*Gauss*cos(phi-pi/2.)
             w = wAmp*Gauss*cos(phi+pi/2.)
             p = pAmp*Gauss*cos(phi-pi/2.) / thetaStrat(k)
             theta = Fr2 * theta00 * b
                          

             ! write to field
             select case( model ) 
             case( "pseudo_incompressible" )
                var(i,j,k,1) = rho
                
             case( "Boussinesq" ) 
                var(i,j,k,6) = theta

             case default
                stop "initialize: unknown case model"
             end select
             
             var(i,j,k,2) = u
             var(i,j,k,4) = w
             var(i,j,k,5) = p
             
          end do
       end do
       
       ! copy u values to ghost cells
       var(0,:,:,2) = var(nx,:,:,2)
       var(nx+1,:,:,2) = var(1,:,:,2)
       
       ! average zonal velocities to cell face...
       do i = 0,nx
          var(i,:,:,2) = 0.5*( var(i,:,:,2) + var(i+1,:,:,2) )
       end do
       
       ! average vertical velocities to cell faces
       do k = 1,nz-1
          var(:,:,k,4) = 0.5*( var(:,:,k,4) + var(:,:,k+1,4) )
       end do
       var(:,:,nz,4) = 0.0        ! reset velocity at wall to zero
       

       !-------------------------------------------------------------------
       

    case( 'wavePacket_raytracer' )
       
       ! read test case input data
       read (unit=10, nml=wavePacket)
       
       
       ! initial background flow
       var(:,:,:,2) = meanFlowX_dim / uRef
       var(:,:,:,3) = 0.0
       var(:,:,1:nz-1,4) = meanFlowZ_dim / uRef
       var(:,:,0,4 ) = 0.0
       var(:,:,nz,4) = 0.0
       
       ! constant pressure variable pi' 
       var(:,:,:,5) = 0.0
       
       ! stratified density
       do k = 1,nz
          var(:,:,k,1) = rhoStrat(k)
       end do
       


       !---------------------------------------------------
       !                Hot and cold bubbles
       !---------------------------------------------------

       case( 'Robert_Bubble' )
          ! read test case input data
          read (unit=10, nml=robert_bubble)
          
          ! atmospheric background flow
          var(:,:,:,2) = backgroundFlow_dim(1)/uRef
          var(:,:,:,3) = backgroundFlow_dim(2)/uRef
          var(:,:,:,4) = backgroundFlow_dim(3)/uRef
          
          ! constant pressure variable pi' 
          var(:,:,:,5) = 0.0

          ! non-dimensionalize input parameters
          dTheta1 = dTheta1_dim/thetaRef
          a1 = a1_dim/lRef
          sigma1 = sigma1_dim/lRef
          xCenter1 = xCenter1_dim/lRef
          zCenter1 = zCenter1_dim/lRef
          
          dTheta2 = dTheta2_dim/thetaRef
          a2 = a2_dim/lRef
          sigma2 = sigma2_dim/lRef
          xCenter2 = xCenter2_dim/lRef
          zCenter2 = zCenter2_dim/lRef
          
          
          ! potential temperature and density
          do k = 1,nz
             do j = 1,ny
                do i = 1,nx
                
                   ! init theta off set
                   dTheta = 0.0
                   
                   ! hot bubble
                   delX = (x(i) - xCenter1) 
                   delZ = (z(k) - zCenter1)
                   r = sqrt(delX**2 + delZ**2)  ! scaled radius
                     
                   if( r<= a1 ) then
                      dTheta = dTheta1
                   else
                      Gauss = exp(-(r-a1)**2/sigma1**2)
                      dTheta = dTheta1*Gauss
                   end if
                   
                   ! cold bubble
                   delX = (x(i) - xCenter2) 
                   delZ = (z(k) - zCenter2)
                   r = sqrt(delX**2 + delZ**2)  ! scaled radius
                     
                   if( r<= a2 ) then
                      dTheta = dTheta + dTheta2
                   else
                      Gauss = exp(-(r-a2)**2/sigma2**2)
                      dTheta = dTheta + dTheta2*Gauss
                   end if

                   
                   ! total potential temperature
                   theta = thetaStrat(k) + dTheta
                   
                   select case( model ) 
                      
                   case( "pseudo_incompressible" ) 
                      
                      ! calc pseudo-incompressible density rho*
                      if ( referenceQuantities == "SI" ) then
                         rho = p0**kappa/Rsp * Pstrat(k) / theta
                      else
                         rho = Pstrat(k) / theta
                      end if
                      var(i,j,k,1) = rho
                      
                   case( "Boussinesq" )
                      
                      ! set pot Temp deviation
                      var(i,j,k,6) = dTheta
                      
                   case default
                      stop "initialize: unknown model."
                   end select
                   
                
                end do
             end do
          end do


!-----------------------------------------------------------------------


    case( 'coldBubble' )

       ! read test case input data
       read (unit=10, nml=bubble)

       if (referenceQuantities == "SI" ) stop "initialize: SI units not allowed"

       ! zero start velocity 
       var(:,:,:,2) = 0.0
       var(:,:,:,3) = 0.0
       var(:,:,:,4) = 0.0
       
       ! constant pressure variable pi' 
       var(:,:,:,5) = 0.0


       ! potential temperature and density
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                x_dim = x(i) * lRef       ! dimensional lenghts
                z_dim = z(k) * lRef
                
                delX = (x_dim - xCenter_dim) / xRadius_dim
                delZ = (z_dim - zCenter_dim) / zRadius_dim
                
                r = sqrt(delX**2 + delZ**2)  ! scaled radius
                
                if( r<=1.0 ) then  ! inside bubble

                   dTheta_dim = 0.5*dTheta0_dim * (1.0 + (cos(pi*r/2.0))**2)
                   theta = thetaStrat(k) - dTheta_dim / thetaRef
                   
                   ! calc pseudo-incompressible density rho*
                   if ( referenceQuantities == "SI" ) then
                      rho = p0**kappa/Rsp * Pstrat(k) / theta
                   else
                      rho = Pstrat(k) / theta
                   end if
                   var(i,j,k,1) = rho

                else  ! outside bubble
                   ! keep background density
                   var(i,j,k,1) = rhoStrat(k)
                end if
                
             end do
          end do
       end do
   

!-----------------------------------------------------------------------


    case( 'hotBubble' )

       ! read test case input data
       read (unit=10, nml=bubble)
       
       ! zero start velocity 
       var(:,:,:,2) = backgroundFlow(1)
       var(:,:,:,3) = 0.0
       var(:,:,:,4) = 0.0
       
       ! constant pressure variable pi' 
       var(:,:,:,5) = 0.0


       ! potential temperature and density
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                x_dim = x(i) * lRef       ! dimensional lenghts
                z_dim = z(k) * lRef
                
                delX = (x_dim - xCenter_dim) / xRadius_dim
                delZ = (z_dim - zCenter_dim) / zRadius_dim
                delZ = zExcentricity * delZ
                
                r = sqrt(delX**2 + delZ**2)  ! scaled radius
                
                if( r<=1.0 ) then
                   !--------------------
                   !    inside bubble
                   !--------------------

                   dTheta_dim = dTheta0_dim * (cos(pi*r/2.0))**2

                   theta = thetaStrat(k) + dTheta_dim / thetaRef
                   
                   ! calc pseudo-incompressible density rho*
                   if ( referenceQuantities == "SI" ) then
                      rho = p0**kappa/Rsp * Pstrat(k) / theta
                   else
                      rho = Pstrat(k) / theta
                   end if
                   var(i,j,k,1) = rho
                   

                else
                   !------------------------------------------
                   !  outside bubble keep background density
                   !------------------------------------------
                   var(i,j,k,1) = rhoStrat(k) 
                end if
                
             end do
          end do
       end do


!-----------------------------------------------------------------------


 case( 'hotBubble3D' )

       ! read test case input data
       read (unit=10, nml=bubble)
       
       if (referenceQuantities == "SI" ) stop "initialize: SI units not allowed"

       ! zero start velocity 
       var(:,:,:,2) = 0.0
       var(:,:,:,3) = 0.0
       var(:,:,:,4) = 0.0
       
       ! constant pressure variable pi' 
       var(:,:,:,5) = 0.0


       ! potential temperature and density
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                x_dim = x(i) * lRef       ! dimensional lenghts
                y_dim = y(j) * lRef
                z_dim = z(k) * lRef

                delX = (x_dim - xCenter_dim) / xRadius_dim
                delY = (y_dim - xCenter_dim) / xRadius_dim
                delZ = (z_dim - zCenter_dim) / zRadius_dim
                
                r = sqrt(delX**2 + delY**2 + delZ**2)  ! scaled radius
                
                if( r<=1.0 ) then          ! inside bubble

                   dTheta_dim = dTheta0_dim * (cos(pi*r/2.0))**2

                   theta = thetaStrat(k) + dTheta_dim / thetaRef
                   
                   ! calc pseudo-incompressible density rho*
                   if ( referenceQuantities == "SI" ) then
                      rho = p0**kappa/Rsp * Pstrat(k) / theta
                   else
                      rho = Pstrat(k) / theta
                   end if
                   var(i,j,k,1) = rho

                else  ! outside bubble
                   ! keep background density
                   var(i,j,k,1) = rhoStrat(k)
                end if
                
             end do
          end do
       end do

       !------------------------------------------------
       !               Bottom topography
       !------------------------------------------------

    case( 'agnesiMountain' )
       
       if (referenceQuantities == "SI" ) stop "initialize: SI units not allowed"

       ! density
       do j = 1,ny
          do i = 1,nx
             var(i,j,1:nz,1) = rhoStrat(1:nz)
          end do
       end do

       ! uniform horizontal initial velocity
       var(:,:,:,2) = backgroundFlow(1)
       var(:,:,:,3) = 0.0
       var(:,:,:,4) = 0.0

       ! constant pressure variable pi' 
       var(:,:,:,5) = 0.0

    case( 'steadyFlow' )

       ! uniform horizontal initial velocity
       var(:,:,:,2) = backgroundFlow(1)
       var(:,:,:,3) = 0.0
       var(:,:,:,4) = 0.0

       ! constant pressure variable pi' 
       var(:,:,:,5) = 0.0


       ! potential temperature and density
       do j = 1,ny
          do i = 1,nx

             var(i,j,1:nz,1) = rhoStrat(1:nz)

          end do
       end do


! ------------- without stratification / gravity ---------------------

       
    case ("densitySineTransport")
       
       rho0 = 1.0;  u0 = 1.0;  v0 = 0.0;  w0 = 0.0

       ! density sine curve in x
       do i = -nbx, nx+nbx
          var(i,:,:,1) = rho0 + 0.1*sin( 2*pi * x(i) )
       end do
       ! constant velocity vector in x-direction
       var(:,:,:,2) = u0
       var(:,:,:,3) = v0
       var(:,:,:,4) = w0
       
       ! constant Exner pressure
       var(:,:,:,5) = p0


! -------------------------------------------------------------------------


    case( 'densityDiscXYZ' )
       ! center and radius of disc
       x0 = 0.5;  y0 = 0.5;  z0 = 0.5;  r0 = 1./3.   
       u0 = 1.0;  v0 = 1.0;  w0 = 1.0       ! constant backround velocity
       rhoDisc = 0.1                        ! density of disc

       ! density 
       var(:,:,:,1) = 1.0
       do i = 1,nx
          do j = 1,ny
             do k = 1,nz
                delX = x(i) - x0
                delY = y(j) - y0
                delZ = z(k) - z0
                r = sqrt(delX**2 + delY**2 + delZ**2)
                if( r >= r0 ) var(i,j,k,1) = var(i,j,k,1) + rhoDisc
             end do
          end do
       end do

       ! horizontal velocity u
       var(:,:,:,2) = u0
       
       ! horizontal velocity v
       var(:,:,:,3) = v0
       
       ! vertical velocity w
       var(:,:,:,4) = w0
       var(:,:,-nbz:2,4) = 0.0            ! zero near walls
       var(:,:,nz-2:nz+nbz,4) = 0.0

       ! pressure variable pi'
       var(:,:,:,5) = 0.0


! -------------------------------------------------------------------------


    case( 'densityDiscXZ' )
       ! transport of a density disc with constant speed in x
       x0 = 0.5;  z0 = 0.5;  r0 = 1./3.   ! center and radius of disc
       u0 = 1.0;  v0 = 0.0;  w0 = 0.0     ! constant backround velocity
       rhoDisc = 0.1                         ! density of disc

       ! density 
       var(:,:,:,1) = 1.0
       do i = 1,nx
          do j = 1,ny
             do k = 1,nz
                delX = x(i) - x0
                delZ = z(k) - z0
                r = sqrt(delX**2 + delZ**2)
                if( r >= r0 ) var(i,j,k,1) = var(i,j,k,1) + rhoDisc
             end do
          end do
       end do

       ! horizontal velocity u
       var(:,:,:,2) = u0
       
       ! horizontal velocity v
       var(:,:,:,3) = v0
       
       ! vertical velocity w
       var(:,:,:,4) = w0
       if( w0 /= 0.0 ) then
          var(:,:,-nbz:2,4) = 0.0            ! zero near walls
          var(:,:,nz-2:nz+nbz,4) = 0.0
       end if

       ! pressure variable pi'
       var(:,:,:,5) = 0.0


! -------------------------------------------------------------------------


    case( 'densityDiscXY' )
       x0 = 0.5;  y0 = 0.5;  r0 = 1./3.   ! center and radius of disc
       u0 = 1.0;  v0 = 1.0;  w0 = 0.0     ! constant backround velocity
       rhoDisc = 0.1                         ! density of disc

       ! density 
       var(:,:,:,1) = 1.0
       do i = 1,nx
          do j = 1,ny
             do k = 1,nz
                delX = x(i) - x0
                delY = y(j) - y0
                r = sqrt(delX**2 + delY**2)
                if( r >= r0 ) var(i,j,k,1) = var(i,j,k,1) + rhoDisc
             end do
          end do
       end do


       ! horizontal velocity u
       var(:,:,:,2) = u0
       
       ! horizontal velocity v
       var(:,:,:,3) = v0
       
       ! vertical velocity w
       var(:,:,:,4) = w0

       ! pressure variable pi'
       var(:,:,:,5) = 0.0


! --------------------------------------------------------------------------

       
    case ("greshoVortexXY")
       updateMass = .true.

       ! setup
       x0 = 0.5; y0 = 0.5;  z0 = 0.5        ! center vortex
       uVortex = 1.0;
       r0 = 1./3.
       ! background
       pInf = 10.0
       rho0 = 1.0;  u0 = 0.0;  v0 = 0.0;  w0 = 0.0

       ! density 
       var(:,:,:,1) = rho0

       ! velocities: background
       var(:,:,:,2) = u0;  var(:,:,:,3) = v0;  var(:,:,:,4) = w0

       ! vortex velocity u
       do i = 0,nx
          do j = 1,ny
             do k = 1,nz
                xu = x(i) + 0.5*dx
                yu = y(j)
                delX = xu-x0
                delY = yu-y0
                r = sqrt(delX**2 + delY**2)
                uPhi = greshoVelocity(uVortex, r/r0)
                var(i,j,k,2) = var(i,j,k,2) - uPhi * delY/r
             end do
          end do
       end do


       ! vortex velocity v
       do i = 1,nx
          do j = 0,ny
             do k = 1,nz
                xv = x(i) 
                yv = y(j) + 0.5*dy
                delX = xv-x0
                delY = yv-y0
                r = sqrt(delX**2 + delY**2)
                uPhi = greshoVelocity(uVortex, r/r0)
                var(i,j,k,3) = var(i,j,k,3) + uPhi * delX/r
             end do
          end do
       end do

       ! pressure variable
       ! set pi'(inf) including ghost cells
       var(:,:,:,5) = kappaInv * pInf**kappa        
       do i = 1,nx
          do j = 1,ny
             do k = 1,nz
                delX = x(i) - x0
                delY = y(j) - y0
                r = sqrt(delX**2 + delY**2)
                p = greshoPressure(pInf, uVortex, r/r0)

                ! pressure variable pi' (ref. Klein)
                var(i,j,k,5) = kappaInv * ( p**kappa - pInf**kappa ) 
             end do
          end do
       end do


! --------------------------------------------------------------------------
      
       
    case ("greshoVortexXZ")
       updateMass = .true.

       ! setup
       x0 = 0.0;  z0 = 0.5        ! center vortex
       uVortex = 1.0;
       r0 = 1./3.
       ! background
       pInf = 10.0
       rho0 = 1.0;  u0 = 1.0;  v0 = 0.0;  w0 = 0.0

       var(:,:,:,1) = rho0

       ! velocities: background
       var(:,:,:,2) = u0;  var(:,:,:,3) = v0;  var(:,:,:,4) = w0

       ! vortex velocity u
       do i = 0,nx
          do j = 1,ny
             do k = 1,nz
                xu = x(i) + 0.5*dx
                zu = z(k)
                delX = xu-x0
                delZ = zu-z0
                r = sqrt(delX**2 + delZ**2)
                uPhi = greshoVelocity(uVortex, r/r0)
                var(i,j,k,2) = var(i,j,k,2) - uPhi * delZ/r
             end do
          end do
       end do

       ! vortex velocity w
       do i = 1,nx
          do j = 1,ny
             do k = 0,nz
                xw = x(i) 
                zw = z(k) + 0.5*dz
                delX = xw-x0
                delZ = zw-z0
                r = sqrt(delX**2 + delZ**2)
                uPhi = greshoVelocity(uVortex, r/r0)
                var(i,j,k,4) = var(i,j,k,4) + uPhi * delX/r
             end do
          end do
       end do

       ! pressure variable
       ! set pi'(inf) including ghost cells
       var(:,:,:,5) = kappaInv * pInf**kappa      
       do i = 1,nx
          do j = 1,ny
             do k = 1,nz
                delX = x(i) - x0
                delZ = z(k) - z0
                r = sqrt(delX**2 + delZ**2)
                p = greshoPressure(pInf, uVortex, r/r0)

                ! pressure variable pi' (ref. Klein)
                var(i,j,k,5) = kappaInv * ( p**kappa - pInf**kappa ) 
             end do
          end do
       end do


! --------------------------------------------------------------------------


    case ("projectionTest")
       ! test the projection step
       
       var(:,:,:,1) = 1.0       ! constant density

       var(:,:,:,2) = 1.0       ! constant velocity in x
       
       ! velocity sine curve along x
       !do i = 0,nx
       !   xu = x(i) + dx/2.
       !   var(i,:,:,2) = 1 + 0.0001*sin( 2*pi*xu)
       !end do
       
       var(:,:,:,3:4) = 0.0     ! zero velocity in y,z

       !var(:,:,:,5) = 1.0       ! constant pressure variable
       ! pressure sine curve 
       do i = 0,nx
          var(i,:,:,5) = 1.0 + 0.01*sin( 2*pi*x(i) )
       end do


       ! pressure fluctuation
!       var(nx/2,ny/2,1,5) = 1.5


! --------------------------------------------------------------------------


    case ("matrixStructureTest")
       ! test the momentum transport at constant density and 
       ! pressure along x. 

       ! constant density
       var(:,:,:,1) = 2.0

!!$       ! velocity sine curve along x
!!$       do i = 0,nx
!!$          xu = x(i) + dx/2.
!!$          var(i,:,:,2) = sin( 2*pi*xu)
!!$       end do
!!$       var(:,:,:,3) = 0.0
!!$       var(:,:,:,4) = 0.0

!!$       ! velocity sine curve along y 
!!$       do j = 0,ny
!!$          yv = y(j) + dy/2.0
!!$          var(:,j,:,3) = sin( 2*pi*yv)
!!$       end do
!!$       var(:,:,:,2) = 0.0
!!$       var(:,:,:,4) = 0.0


       ! velocity sine curve along z
       do k = 0,nz
          zw = z(k) + dz/2.0
          var(:,:,k,4) = sin( 2*pi*zw)
       end do
       var(:,:,:,2) = 0.0
       var(:,:,:,3) = 0.0


       ! constant Exner pressure
       var(:,:,:,5) = 1.0


! --------------------------------------------------------------------------
       
       
    case ("momentumTransportX")
       ! test the momentum transport at constant density and 
       ! pressure along x. 
       updateMass = .false.
       correctMomentum = .false.

       ! constant density
       var(:,:,:,1) = 2.0

       ! velocity sine curve along x
       do i = 0,nx
          xu = x(i) + dx/2.
          var(i,:,:,2) = 1.0 + 1.0e-5 * sin( 2*pi*xu)
       end do
     
       ! zero velocity along y and z
       var(:,:,:,3) = 0.0
       var(:,:,:,4) = 0.0

       ! constant Exner pressure
       var(:,:,:,5) = 1.0


! --------------------------------------------------------------------------


    case ("momentumTransportY")
       ! test the momentum transport at constant density and 
       ! pressure along y. 
       updateMass = .false.
       correctMomentum = .false.
       
       ! constant density
       var(:,:,:,1) = 2.0

       ! velocity sine curve along y 
       do j = 0,ny
          yv = y(j) + dy/2.0
          var(:,j,:,3) = 50.0 + sin( 2*pi*yv)
       end do
       
       ! zero velocity along x and z
       var(:,:,:,2) = 0.0
       var(:,:,:,4) = 0.0

       ! constant Exner pressure
       var(:,:,:,5) = 1.0


! --------------------------------------------------------------------------


    case ("momentumFluxTest")
       ! constant density
       var(:,:,:,1) = 2.0

       ! constant velocity fields
       var(:,:,:,2) = 1.0
       var(:,:,:,3) = 2.0
       var(:,:,:,4) = 3.0

       ! constant Exner pressure
       var(:,:,:,5) = 1.0


! --------------------------------------------------------------------------


    case ("densityTransportX")
       ! density sine curve in x
!       predictMomentum = .false.
!       correctMomentum = .false.
       
       u0 = 1.0;  v0 = 0.0;  w0 = 0.0

       do i = 1,nx
          var(i,:,:,1) = sin( 2*pi * x(i) )
       end do
       ! constant velocity vector in x-direction
       var(:,:,:,2) = u0
       var(:,:,:,3) = v0
       var(:,:,:,4) = w0
       var(:,:,:,5) = p0


! --------------------------------------------------------------------------


    case ("densityTransportY")
!       predictMomentum = .false.
!       correctMomentum = .false.
       ! density sine curve in y

       u0 = 0.0;  v0 = 1.0;  w0 = 0.0
       do j = 1,ny
          var(:,j,:,1) = sin( 2*pi * (y(j)-0.25) )
       end do
       ! constant velocity vector in y-direction
       var(:,:,:,2) = u0
       var(:,:,:,3) = v0
       var(:,:,:,4) = w0

       ! constant Exner pressure
       var(:,:,:,5) = p0
       

! --------------------------------------------------------------------------


    case ("massFluxTest")
       ! constant density
       var(:,:,:,1) = 1.5

       ! constant velocity field:
       var(:,:,:,2) = 1.0
       var(:,:,:,3) = 2.0
       var(:,:,:,4) = 3.0

       ! constant Exner pressure
       var(:,:,:,5) = 1.0   

    case ("standard")

       do k = 1,nz
          var(:,:,k,:) = z(k)   ! Schichtung
       end do

    case("xparabel")
       ! rho parabolic in x
       do i = -nbx, nx+nbx
          var(i,:,:,1) = (x(i)+dx/2)**3 - (x(i)-dx/2)**3
          var(i,:,:,1) = 1.0/3.0/dx * var(i,:,:,1) 
       end do

    case("yparabel")
       ! u-velocity parabolic in y
       do j = -nby, ny+nby
          var(:,j,:,2) = (y(j)+dy/2)**3 - (y(j)-dy/2)**3
          var(:,j,:,2) = 1.0/3.0/dy * var(:,j,:,2) 
       end do

    case("zparabel")
       ! u-velocity parabolic in z
       do k = -nbz, nz+nbz
          var(:,:,k,2) = (z(k)+dz/2)**3 - (z(k)-dz/2)**3
          var(:,:,k,2) = 1.0/3.0/dz * var(:,:,k,2) 
       end do

    case default

       print *,"init.f90/initialise: testCase = ", testCase
       stop "init.f90/initialise: This testCase is not valid. Stop."

    end select

    
    ! close input file pinc.f
    close (unit=10)


    
    !----------------------------------
    !     Output system settings 
    !----------------------------------

    print*,""
    print*,"Initializing System: "
    print*,"  1) Reference quantities: "

    write(*,fmt="(a25,f7.3,a)") "rhoRef = ", rhoRef, " kg/m^3"
    write(*,fmt="(a25,f7.3,a)") "pRef = ", pRef*1.e-3, " kPa"
    write(*,fmt="(a25,f7.3,a)") "aRef,uRef = ", uRef, " m/s"    
    write(*,fmt="(a25,f7.3,a)") "lRef = ", lRef*1.e-3, " km"    
    write(*,fmt="(a25,f7.3,a)") "tRef = ", tRef, " s"
    write(*,fmt="(a25,f7.3,a)") "thetaRef = ", thetaRef, " K"
    write(*,fmt="(a25,f7.3,a)") "FRef = ", FRef, " N/m^3"
    print*,""

    print*,"  2) Non-dimensional numbers: "

    write(*,fmt="(a25,es10.3)") "Ma = ", Ma
    write(*,fmt="(a25,es10.3)") "Fr = ", Fr 
    write(*,fmt="(a25,es10.3)") "Re = ", Re
    print*,""

    
    print*,"  3) Extreme values: "
    write(*,fmt="(a25,es10.1,a,f5.1,a)") "PStrat = ", PStrat(nz)*pRef/1000.0,&
         & " kPa at z = ", z(nz) * lRef/1000.0, " km"
    write(*,fmt="(a25,es10.1,a,f5.1,a)") "rhoStrat = ", rhoStrat(nz)*rhoRef,&
         & " kg/m3 at z = ", z(nz) * lRef/1000.0, " km"
    write(*,fmt="(a25,es10.1,a,f5.1,a)") "thetaStrat = ", thetaStrat(nz)*thetaRef,&
         & " K at z = ", z(nz) * lRef/1000.0, "km"
    print*,""

    print*,"  4) Constants: "
    write(*,fmt="(a25,f7.3,a)") "gamma = ", gamma, " "
    write(*,fmt="(a25,f7.3,a)") "g = ", g, " m/s^2"
    write(*,fmt="(a25,f7.3,a)") "R_sp = ", Rsp, " J/kg/K"
    write(*,fmt="(a25,f7.3,a)") "f_Coriolis = ", f_Coriolis_dim, " 1/s"
    write(*,fmt="(a25,f7.3,a)") "mu_viscous = ", mu_viscous_dim, " m^2/s"
    write(*,fmt="(a25,f7.3,a)") "mu_conduct = ", mu_conduct_dim, " m^2/s"
    print*,""


    print*,"  5) Background: "
    select case (background)
    case ("isothermal") 
       write(*,fmt="(a25,a15)") "background = ", background
       write(*,fmt="(a25,f7.3,a)") "T0 = ", Temp0_dim, " K"
       write(*,fmt="(a25,f7.4,a)") "N = ", sqrt(N2)/tRef, " 1/s"
    case ("isentropic")
       write(*,fmt="(a25,a7)") "background = ", background
       write(*,fmt="(a25,f7.3,a)") "theta_0 = ", theta0_dim, " K"
    case ("const-N")
       write(*,fmt="(a25,a7)") "background = ", background
       write(*,fmt="(a25,f7.3,a)") "N = ", N_BruntVaisala_dim, " 1/s"
    case ( "uniform" )
       write(*,fmt="(a25,a15)") "background = ", background
       write(*,fmt="(a25,f7.3,a)") "rho0 = ", rhoStrat(1)*rhoRef, " kg/m^3"
       write(*,fmt="(a25,f7.1,a)") "theta0 = ", thetaStrat(1)*thetaRef, " K"
       write(*,fmt="(a25,f7.2,a)") "N_Brunt-Vaisala = ", NN/tRef, " 1/s"
    end select
    write(*,fmt="(a25,f7.1,a)") "u0 = ", &
         & backgroundFlow(1)*uRef," m/s"
    write(*,fmt="(a25,f7.1,a)") "v0 = ", &
         & backgroundFlow(2)*uRef," m/s"
    write(*,fmt="(a25,f7.1,a)") "w0 = ", &
         & backgroundFlow(3)*uRef," m/s"
    print*,""
    
    
    print*,"  6) Model equations: "
    write(*,fmt="(a25,a15)") "model = ", model
    write(*,fmt="(a25)") "inclination angle: "
    write(*,fmt="(a25,f5.1,a)") "theta = ", vert_theta, " deg"
    write(*,fmt="(a25,f5.1,a)") "alpha = ", vert_alpha, " deg"
    print*,""


   
    
    print*,"  7) Domain: "
    write(*,fmt="(a25,a,f6.1,a,f6.1,a)") "[xMin, xMax] = ", &
         & "[",lx(0)*lRef/1000.0,"km, ",lx(1)*lRef/1000.0,"km ]" 
    write(*,fmt="(a25,a,f6.1,a,f6.1,a)") "[yMin, yMax] = ", &
         & "[",ly(0)*lRef/1000.0,"km, ",ly(1)*lRef/1000.0,"km ]" 
    write(*,fmt="(a25,a,f6.1,a,f6.1,a)") "[zMin, zMax] = ", &
         & "[",lz(0)*lRef/1000.0,"km, ",lz(1)*lRef/1000.0,"km ]"
    write(*,fmt="(a25,i4,a,i4,a,i4)") "nx x ny x nz = ", &
         & nx," x ", ny, " x ", nz
    print*,""


    print*,"  8) Boundary conditions: "
    write(*,fmt="(a25,a)") "xBoundary = ", trim(xBoundary)
    write(*,fmt="(a25,a)") "yBoundary = ", trim(yBoundary)
    write(*,fmt="(a25,a)") "zBoundary = ", trim(zBoundary)
    if( spongeLayer ) then
       write(*,fmt="(a25,a)") "sponge layer = ", "on"
       write(*,fmt="(a25,f5.1,a)") "height = ", &
            & (lz(1)-lz(0))*spongeHeight*lRef/1000.0," km"
       write(*,fmt="(a25,es7.1,a)") "relaxation  = ", spongeAlphaZ_dim, " 1/s"
    else
       write(*,fmt="(a25,a)") "sponge layer = ", "off"
    end if
    print*,""

    
    
    print*,"  9) Poisson Solver: "
    write(*,fmt="(a25,a)") "solver = ", poissonSolverType
    write(*,fmt="(a25,es7.1)") "tolPoisson = ", tolPoisson
    print*,""
    
    
    print*," 10) Topography: "
    if( topography) then 
       write(*,fmt="(a25,a)") "topography = ", "on"
       write(*,fmt="(a25,f6.1,a)") "mountain height = ", mountainHeight_dim, " m"
       write(*,fmt="(a25,f6.1,a)") "mountain width = ", mountainWidth_dim, " m"
    else
       write(*,fmt="(a25,a)") "topography = ", "off"
    end if
    print*,""
    
    


!---------------------------------------------------------------------------


  contains

    function discDensity(r)
      ! in/out
      real :: discDensity
      real, intent(in) :: r
      
      ! local variables
      real :: dens
      
      
      if (r >= 1.0) then   ! outer region 
         dens = 0.0
      else if (r < 0.5) then ! 0 <= r < r0/2  region
         dens = 4.0 * r**2
      else ! 1/2 <= r < 1
         dens = 4.0 * (1.0 - r)**2
      end if

      discDensity = dens
      

    end function discDensity


!---------------------------------------------------------------------------


    function greshoVelocity (u0, r) ! gresho vortex with normed radius r0=1
      ! in/out
      real :: greshoVelocity
      real, intent(in) :: u0, r

      ! local variables
      real :: u 
      
      if (r >= 1.0) then   ! outer region 
         u = 0.0
      else if (r < 0.5) then ! 0 <= r < r0/2  region
         u = 2.0 * r
      else ! 1/2 <= r < 1
         u = 2.0*(1.0-r)
      end if

      greshoVelocity = u*u0

    end function greshoVelocity


!---------------------------------------------------------------------------

    
    function greshoPressure(pInf, u0, r)
      ! in/out
      real :: greshoPressure
      real, intent(in) :: pInf, u0, r
      
      ! local variables
      real :: p
      
      
      if (r >= 1.0) then   ! outer region 
         p = 0.0
      else if (r < 0.5) then ! 0 <= r < r0/2  region
         p = 4.0 * r**2
      else ! 1/2 <= r < 1
         p = 4.0 * (1.0 - r)**2
      end if

      greshoPressure = pInf - 0.5 * u0**2 * p

    end function greshoPressure

      
  end subroutine initialise


end module init_module
