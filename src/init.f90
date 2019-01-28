module init_module


  use type_module
  use atmosphere_module

  implicit none

  private         ! private module

  public  :: initialise
  public  :: setup
  public  :: init_ParamStudy

  private :: cphase


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


  subroutine setup (var,var0,flux,force,source,dRho,dMom,dTheta,dIce)
    !-----------------------------------------
    ! allocate var and flux / read input.f90
    !-----------------------------------------

    ! in/out variables 
    real,dimension(:,:,:,:), allocatable,intent(out) :: var, var0, source
    real,dimension(:,:,:,:,:), allocatable,intent(out) :: flux
    real,dimension(:,:,:,:), allocatable,intent(out) :: force
    real,dimension(:,:,:), allocatable :: dRho        ! RK-Update for rho
    real,dimension(:,:,:,:), allocatable :: dMom      ! ...rhoU,rhoV,rhoW
    real,dimension(:,:,:), allocatable :: dTheta       ! RK-Update for theta
    real,dimension(:,:,:,:), allocatable :: dIce       ! RK-Update for nIce,qIce,SIce


    integer :: allocstat
    integer :: i, j, k, iVar


    ! constants
    pi = 4*atan(1.0)


    ! open input file input.f90
    open (unit=10, file="input.f90", action="read", &
         form="formatted", status="old", position="rewind")

    ! read grid info
!    read (unit=10, nml=grid)   ! modified by Junhong Wei (20161110)

    ! total dimension of var fields with ghost cells
!    nxx = nx + 2*nbx + 1   ! modified by Junhong Wei (20161110)
!    nyy = ny + 2*nby + 1   ! modified by Junhong Wei (20161110)
!    nzz = nz + 2*nbz + 1   ! modified by Junhong Wei (20161110)


    !---------------------------------------------------
    ! allocate x,y,z - cell centered coordinate fields
    !---------------------------------------------------

!    allocate(x(-nbx:nx+nbx), stat=allocstat)   ! modified by Junhong Wei (20161104)
    allocate(x(-nbx:sizeX+nbx), stat=allocstat)   ! modified by Junhong Wei (20161104)
    if(allocstat /= 0) stop "init.f90: could not allocate x"

!    allocate(y(-nby:ny+nby), stat=allocstat)   ! modified by Junhong Wei (20161104)
    allocate(y(-nby:sizeY+nby), stat=allocstat)   ! modified by Junhong Wei (20161104)
    if(allocstat /= 0) stop "init.f90: could not allocate y"

!    allocate(z(-nbz:nz+nbz), stat=allocstat)   ! modified by Junhong Wei (20161104)
    allocate(z(-nbz:sizeZ+nbz), stat=allocstat)   ! modified by Junhong Wei (20161104)
    if(allocstat /= 0) stop "init.f90: could not allocate z"

!   achatzb
    !---------------------------------------------------
    ! allocate topography mask and surface
    !---------------------------------------------------

    allocate(&
    topography_mask(-nbx+1:sizeX+nbx,-nby+1:sizeY+nby,-nbz+1:sizeZ+nbz), &
    stat=allocstat&
    )
    if(allocstat /= 0) stop "init.f90: could not allocate topgoraphy_mask"
    allocate(&
    topography_surface(-nbx+1:sizeX+nbx,-nby+1:sizeY+nby), stat=allocstat&
    )
    if(allocstat /= 0) &
    stop "init.f90: could not allocate topography_surface"
!   achatze


    !-------------------------------------
    !      allocate variable fields
    !-------------------------------------


    read (unit=10, nml=variables)
    if (include_ice) nVar = nVar+3

    ! allocate var = (rho,u,v,w,pEx)
    allocate(var(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate var"

    ! allocate var0 = (rho,u,v,w,pEx)
    allocate(var0(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate var0"

    ! allocate source for rho,u,v,w,pEx,theta
    allocate(source(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate source"

    ! allocate dRho
    allocate(dRho(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dRho."

    ! allocate dMom
    allocate(dMom(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,3),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dMom."

    ! allocate dTheta
    allocate(dTheta(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dTheta."

    ! allocate dIce
    if (include_ice) then
      allocate(dIce(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,3),stat=allocstat)
      if(allocstat /= 0) stop "init.f90: Could not allocate dIce."
    end if

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

    ! allocate varIn, varOut and offset
    allocate(varIn(nVar),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate varIn"

    allocate(varOut(nVar),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate varOut"

    allocate(offset(nVar-2),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate offset"

    ! read output specifications
    read (unit=10, nml=outputList)
    if (include_ice) then
      forall (i=0:2)
        varOut(nVar-i) = 1
        varIn(nVar-i) = 1
      end forall
    end if

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

    ! open info file input.f90
    open (unit=90, file="info.txt", action="write", &
         form="formatted", status="replace")

    ! write model info file
    write(90,fmt="(a)") ""
    write(90,fmt="(a)") "Model: "//trim(model)
    write(90,fmt="(a)") ""


    select case( model ) 

    case( "anelastic" )

       pressureScaling = .false.


    case( "Boussinesq" )

       updateMass = .false.
       predictMomentum = .true.
       correctMomentum = .true.
       updateTheta = .true.
       updateIce = .false.

       ! overwrite unsuitable input settings
       topography = .false. 
       spongeLayer = .false.
       pressureScaling = .false.

       ! do not work with fluctuations
       fluctuationMode = .false.

       ! never offset theta 
       thetaOffset = .false.



    case( "pseudo_incompressible" ) 

       updateMass = .true.
       predictMomentum = .true.
       correctMomentum = .true.
       updateTheta = .false.
       updateIce = .true.

       if( master ) then   ! modified by Junhong Wei (20170216)
       write(90,"(a25)",advance = "no") "updateMass != "
       write(90,*) updateMass
       write(90,"(a25)",advance = "no") "predictMomentum != "
       write(90,*) predictMomentum
       write(90,"(a25)",advance = "no") "correctMomentum != "
       write(90,*) correctMomentum
       write(90,"(a25)",advance = "no") "updateTheta != "
       write(90,*) updateTheta
       write(90,*) ""
       write(90,"(a25)",advance = "no") "updateIce != "
       write(90,*) updateIce
       write(90,*) ""
       end if   ! modified by Junhong Wei (20170216)

       ! overwrite unsuitable input settings
       if( zBoundary == "periodic" ) then
          print*,"WARNING: zBoundary periodic not possible. &
               & Reset to solid_wall!"
          zBoundary = "solid_wall"
          write(90,"(a)") "WARNING: zBoundary periodic not possible. &
               & Reset to solid_wall!"
       end if

    case( "WKB" )

       raytracer = .true.
       updateMass = .false.
       predictMomentum = .false.
       correctMomentum = .false.
       updateTheta = .false.
       updateIce = .false.

       fluctuationMode = .false.       ! Work with full density
       background = "isothermal"       ! Theory only for isothermal bg.

       !-----------------
       ! Write info file
       !-----------------
       write(90,"(a25)",advance = "no") "raytracer != "
       write(90,*) raytracer
       write(90,"(a25)",advance = "no") "updateMass != "
       write(90,*) updateMass
       write(90,"(a25)",advance = "no") "predictMomentum != "
       write(90,*) predictMomentum
       write(90,"(a25)",advance = "no") "correctMomentum != "
       write(90,*) correctMomentum
       write(90,"(a25)",advance = "no") "updateTheta != "
       write(90,*) updateTheta
       write(90,"(a25)",advance = "no") "fluctuationMode != "
       write(90,*) fluctuationMode
       write(90,"(a25)",advance = "no") "background != "
       write(90,"(a)") background
       write(90,"(a25)",advance = "no") "updateIce != "
       write(90,"(a)") updateIce
       write(90,*) ""


    case default
       print*,"model = ", model
       stop "initialize: Unknown model" 
    end select

    close(90)  ! info file

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

    integer :: i0, j0   ! modified by Junhong Wei (20161121)

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
    real :: ll_3DWP                      ! meridional wave number ! modified by Junhong Wei for 3DWP (20170828)
    real :: omi, omi2                    ! intrinsic frequency, squared 
    real :: bAmp, uAmp, wAmp, pAmp       ! amplitudes for buoyancy, u, w, Exner function
    real :: delx2, delz2                 ! squared distance from center
    real :: envel                        ! envelope of wave packet
    real :: Gauss                        ! Gaussian distribution value
    real :: sigma                        ! width of Gaussian distribution
    real :: sigma_hor                    ! width of Gaussian distribution (horizontal direction)   ! modified by Junhong Wei (20170214)
    real :: sigma_hor_yyy                ! width of Gaussian distribution   ! modified by Junhong Wei for 3DWP (20170828)
    real :: L_cos                        ! width of Cosine distribution
    real :: xCenter, zCenter             ! center of distribution
    real :: phi                          ! local phase 
    real :: u,v,w,b                        ! buoyancy, zonal + vertical velocity
    real :: lambdaX, lambdaZ             ! zonal and vertical wave length

    real :: lambdaY, yCenter ! variables for 3DWP   ! modified by Junhong Wei for 3DWP (20170828)


    ! wave 2 
!    complex, dimension(0:nx+1,0:ny+1,0:nz+1,4,0:2)  :: Psi ! modified by Junhong Wei
    complex, dimension(0:nx+1,0:ny+1,0:nz+1,5,0:2)  :: Psi ! modified by Junhong Wei
    real :: u1,w1,b1,p1
    real :: u2,w2,b2,p2
!    logical, parameter :: initWave2 = .true.      ! modified by Junhong Wei
    logical, parameter :: initWave2 = .false.      ! modified by Junhong Wei


    ! monochromatic wave
    real :: th, f, f2
    real :: lambda
    real :: amp, cot
    real :: uRot, vRot, wRot

    ! testing with Stefan
    real :: xx

    ! random noise on background
    real, dimension(0:nx+1,0:ny+1,0:nz+1) :: randNoise         ! noise on density
    real, parameter           :: randAmp = 0.0     ! amplitude of random noise

    ! jet stream 
    real :: u_jet                        ! jet stream velocity
    real :: L_jet
    real :: z0_jet

    ! debug   ! modified by Junhong Wei (20161121)
    integer :: i00, j00   ! modified by Junhong Wei (20161121)

    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar) :: flux ! modified by Junhong Wei (20161205)


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

       if( master ) then ! modified by Junhong Wei (20161121)
       print*,"WARNING: only wave packet test case is &
            & changed to density fluctuation mode!"
       end if ! modified by Junhong Wei (20161121)
       
       ! set vertical always paralle to z-axis
       vertical = (/0.0, 0.0, 1.0 /)

    case( "WKB" ) 
       ! 

    case default
       stop "initialize: unknown case model."
    end select

! modified by Junhong Wei (20161121) *** starting line ***
    !-----------------------
    !      MPI stuff
    !-----------------------
    i0 = is + nbx - 1   ! 0 index, replace i -> i + i0 in x and y fields
    j0 = js + nby - 1
! modified by Junhong Wei (20161121) *** finishing line ***
    

! on default there is no initial ice or humidity in the atmosphere
       if (include_ice) var(:,:,:,nVar-2:nVar) = 0.0
!---------------------------------------------------------------

    select case (testCase)

       !-------------------------------------
       !             Boussinesq only
       !-------------------------------------


    case( "sinus" ) 
       var(:,:,:,1) = 1.0
       var(:,:,:,2:4) = 0.0


       do i = -2,nx+3

          xx = lx(0) + (real(i-1)+0.5)*dx
          theta = sin(pi*xx*lRef)
          var(i,1,1,6) = theta

       end do


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

       ! read test case input data
       read (unit=10, nml=monochromeWave)

       ! gather quantities
       th = vert_theta * pi/180.0
       f = f_Coriolis_dim * tRef
       f2 = f**2
       omi = -sqrt(N2*cos(th)**2 + f2*sin(th)**2)


       ! wave vector
       lambda = lambda_dim / lRef
       kTot = 2.0*pi / lambda

       kk = kTot*cos(th)
       mm = kTot*sin(th)

       ! amplitudes in non-rotated frame
       amp = amplitudeFactor*(-omi/kk)
       cot = 1.0/tan(th)

       do k = 1,nz
          j = 1
          do i = 1,nx

             ! phase
             phi = -kk*x(i) - mm*z(k)

             ! amplitudes
             u = amp * cos(phi)
             v = amp*f/omi * sin(phi)
             w = amp*cot * cos(phi)
             b = amp*N2/omi*cot * sin(phi)

             ! project on rotated coordinate system

             wRot = dot_product( (/u,v,w/), vertical) 
             vRot = v
             uRot = sqrt(u**2+w**2 - wRot**2)
             theta = Fr2 * theta00 * b

             ! assign to var field
             var(i,j,k,2:4) = (/uRot,vRot,wRot/)
             var(i,j,k,6) = theta

          end do
       end do


       !-----------------------------
       !  Interpolate to cell faces
       !-----------------------------

       ! copy u values to ghost cells
       var(0,:,:,2) = var(nx,:,:,2)
       var(nx+1,:,:,2) = var(1,:,:,2)

       ! average zonal velocities to cell face...
       do i = 0,nx
          var(i,:,:,2) = 0.5*( var(i,:,:,2) + var(i+1,:,:,2) )
       end do

       ! copy v values to ghost cells
       var(:,0,:,3) = var(:,ny,:,3)
       var(:,ny+1,:,3) = var(:,1,:,3)

       ! average meridional velocities to cell face...
       do j = 0,ny
          var(:,j,:,3) = 0.5*( var(:,j,:,3) + var(:,j+1,:,3) )
       end do

       ! copy w values to ghost cells
       var(:,:,0,4) = var(:,:,nz,4)
       var(:,:,nz+1,4) = var(:,:,1,4)

       ! average vertical velocities to cell face...
       do k = 0,nz
          var(:,:,k,4) = 0.5*( var(:,:,k,4) + var(:,:,k+1,4) )
       end do

       
       !----------------------------------------
       !          WKB theory 
       !----------------------------------------

    case( 'wavePacket' )              ! 1D/2D wave packet

       !---------------------
       ! set up random noise
       !---------------------

       call random_number(randNoise)
       do k = 1,nz
          randNoise(:,:,k) = randNoise(:,:,k)*randAmp
       end do

       !--------------------
       ! set up jet stream
       !--------------------

       ! read test case input data
       read (unit=10, nml=wavePacket)

       u0 = u0_jet_dim / uRef             ! amplitude of jet
       L_jet = L_jet_dim / lRef           ! half width of cos profile
       z0_jet = z0_jet_dim / lRef         ! center of jet

       do k = 1,nz
          delz = (z(k)-z0_jet)

          ! Cosine
          if( abs(delz) .le. L_jet ) then
             u_jet = 0.5*u0*(1.0 + cos(pi*delz/L_jet))
          else
             u_jet = 0.0
          end if

          var(:,:,k,2) = u_jet

       end do


       !--------------------
       !     set up GWP
       !--------------------
       close(10)
       call init_GWP(Psi,kk,mm, ll_3DWP )   ! modified by J. Wei for 3DWP
       open (unit=10, file="input.f90", action="read", &
         form="formatted", status="old", position="rewind")
       ! added close and open statements to avoid unclear position specifiers in init_GWP

       do k = 0,(nz+1)   ! modified by Junhong Wei for 3DWP (20171204)
         do j = 0,(ny+1)   ! modified by Junhong Wei for 3DWP (20170922)
          do i = 0,(nx+1) ! modified by Junhong Wei (20161207)

              phi = kk*x(i+i0) + mm*z(k) + ll_3DWP*y(j+j0)
             
             ! wave 1
             u1 = real( Psi(i,j,k,1,1) * exp(phi*imag) )
             w1 = real( Psi(i,j,k,2,1) * exp(phi*imag) ) 
             b1 = real( Psi(i,j,k,3,1) * exp(phi*imag) ) 
             p1 = real( Psi(i,j,k,4,1) * exp(phi*imag) ) 


             ! wave 2
             if( initWave2 ) then
                u2 = real( Psi(i,j,k,1,2) * exp(2.*phi*imag) )
                w2 = real( Psi(i,j,k,2,2) * exp(2.*phi*imag) ) 
                b2 = real( Psi(i,j,k,3,2) * exp(2.*phi*imag) ) 
                p2 = real( Psi(i,j,k,4,2) * exp(2.*phi*imag) ) 
             end if


             ! sum of wave 1 and 2
             if( initWave2 ) then
                b = b1 + b2
                u = u1 + u2
                w = w1 + w2
                p = p1 + p2
             else
                b = b1
                u = u1
                w = w1
                p = p1                
             end if


             ! additional vars
             rho = 1./(1.+Fr2*b) * rhoStrat(k) 
             theta = Fr2 * theta00 * b


             ! write to field
             select case( model ) 
             case( "pseudo_incompressible" )

                ! add random noise
                rho = rho + randNoise(i,j,k)

                ! subtract background for fluctuation mode
                if( fluctuationMode ) rho = rho - rhoStrat(k)

                ! write to field
                var(i,j,k,1) = rho 

             case( "Boussinesq" ) 
                var(i,j,k,6) = theta

             case default
                stop "initialize: unknown case model"
             end select

             var(i,j,k,2) = var(i,j,k,2) + u
             var(i,j,k,4) = w
             var(i,j,k,5) = p

             var(i,j,k,3) = real( Psi(i,j,k,5,1) * exp(phi*imag) )

          end do
         end do   ! modified by Junhong Wei for 3DWP (20170922)
       end do

       ! average zonal velocities to cell face...
       do i = 0,nx
          var(i,:,:,2) = 0.5*( var(i,:,:,2) + var(i+1,:,:,2) )
       end do

       ! average vertical velocities to cell faces
       do k = 0,nz   ! modified by Junhong Wei for 3DWP (20171204)
          var(:,:,k,4) = 0.5*( var(:,:,k,4) + var(:,:,k+1,4) )
       end do

             select case( model ) 
             case( "pseudo_incompressible" )

          var(:,:,0,4)  = 0.0        ! reset velocity at wall to zero
          var(:,:,nz,4) = 0.0        ! reset velocity at wall to zero

             case( "Boussinesq" ) 


             case default
                stop "initialize: unknown case model"
             end select


       ! average meridional velocities to cell face...
       do j = 0,ny
          var(:,j,:,3) = 0.5*( var(:,j,:,3) + var(:,j+1,:,3) )
       end do

       ! no initial ice in the atmosphere, but supersaturation SIce=1.5
       if (include_ice) then 
         var(:,:,:,nVar-2:nVar-1) = 0.0
         var(:,:,:,nVar) = 1.5
       end if

!---------------------------------------------------------------

!   achatzb
    !   -----------------------------------------------------------------

    !   read parameters for temporary wind relaxation
    !   zero-wind initial state for montain-wave simulations 

    case( 'mountainwave' )  
       ! read parameters for temporary wind relaxation

       read (unit=10, nml=mountainwavelist)

       ! nondimensionalization

       u_relax = u_relax/uRef

       t_relax = t_relax/tRef
       t_ramp = t_ramp/tRef

       xextent_norelax = xextent_norelax/lRef

       ! increase relaxation wind u_relax so that u = u_relax after the
       ! relaxation period (in zonally symmetric case without topography)

       u_relax = u_relax/(1.0 - exp(4.0*t_ramp/(pi*t_relax) - 1.0))

       ! zero wind
  
       var(:,:,:,2) = 0.0 ! u
       var(:,:,:,3) = 0.0 ! v
       var(:,:,:,4) = 0.0 ! w


       ! density, potential temperature, and pressure

       do k = 0,(nz+1)
          do j = 0,(ny+1)
             do i = 0,(nx+1)
                select case( model ) 
                   case( "pseudo_incompressible" )
                      ! initialization density = background density
                      ! subtract background for fluctuation mode

                      if( fluctuationMode ) then
                         rho = 0.0
                        else
                         rho = rhoStrat(k) 
                      end if

                      ! write to field
                      var(i,j,k,1) = rho 

                   case( "Boussinesq" ) 
                      ! initialization zero buoyancy fluctuations

                      var(i,j,k,6) = 0.0

                   case default
                      stop "initialize: unknown case model"
                end select

                ! initialization zero pressure fluctuations

                var(i,j,k,5) = 0.0
             end do
          end do
       end do
!   achatze

    !-------------------------------------------------------------------


    case( 'wavePacket_raytracer' )

       ! for the full set up see setup_wkb routine in wkb.f90

       ! read test case input data
       read (unit=10, nml=wavePacket)


       ! stratified density
       do k = 1,nz
          var(:,:,k,1) = rhoStrat(k)
       end do

       ! initial background flow
       var(:,:,:,2) = meanFlowX_dim / uRef
       var(:,:,:,3) = 0.0
       var(:,:,1:nz-1,4) = meanFlowZ_dim / uRef
       var(:,:,0,4 ) = 0.0
       var(:,:,nz,4) = 0.0

       ! constant pressure variable pi' 
       var(:,:,:,5) = 0.0


       !--------------------
       ! set up jet stream
       !--------------------
       u0 = u0_jet_dim / uRef                          ! amplitude of jet
       L_jet = L_jet_dim / lRef                        ! half width of cos profile
       z0_jet = z0_jet_dim / lRef                      ! center of jet

       do k = 1,nz
          delz = (z(k)-z0_jet)

          ! Cosine profile
          if( abs(delz) .le. L_jet ) then
             u_jet = 0.5*u0*(1.0 + cos(pi*delz/L_jet))
          else
             u_jet = 0.0
          end if

          var(:,:,k,2) = var(:,:,k,2) + u_jet

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

                   ! subtract background for fluctuation mode
                   if( fluctuationMode ) rho = rho - rhoStrat(k)

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

       ! start velocity 
       var(:,:,:,2) = backgroundFlow(1)
       var(:,:,:,3) = backgroundFlow(2)
       var(:,:,:,4) = backgroundFlow(3)

       ! constant pressure variable pi' 
       var(:,:,:,5) = 0.0

       ! zero potential temperature
       var(:,:,:,6) = 0.0

       ! potential temperature and density

       dTheta0 = dTheta0_dim / thetaRef

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

                   dTheta = dTheta0 * (cos(pi*r/2.0))**2

                   theta = thetaStrat(k) + dTheta

                   select case( model ) 

                   case( "pseudo_incompressible" ) 

                      if( fluctuationMode )  then

                         ! calc pseudo-incompressible density rho*
                         if ( referenceQuantities == "SI" ) then
                            rho = p0**kappa/Rsp * Pstrat(k) / theta - rhoStrat(k)
                         else
                            rho = Pstrat(k) / theta - rhoStrat(k)
                         end if
                      else

                         ! calc pseudo-incompressible density rho*
                         if ( referenceQuantities == "SI" ) then
                            rho = p0**kappa/Rsp * Pstrat(k) / theta
                         else
                            rho = Pstrat(k) / theta
                         end if
                      end if

                      var(i,j,k,1) = rho

                   case( "Boussinesq" )

                      ! set pot Temp deviation
                      var(i,j,k,6) = dTheta

                   case default
                      stop "initialize: unknown model."
                   end select

                else
                   !------------------------------------------
                   !  outside bubble keep background density
                   !------------------------------------------

                   if( fluctuationMode ) then
                      var(i,j,k,1) = 0.0
                   else
                      var(i,j,k,1) = rhoStrat(k) 
                   end if

                end if

             end do
          end do
       end do

       ! modified by Junhong Wei (20161122) *** starting line ***

       !-----------------------------------------------------------------------

       
    case( 'hotBubble2D' )

!       ! read test case input data
!       if( master ) read (unit=10, nml=bubble)
!
!       ! broadcase input data
!       call mpi_bcast(dTheta0_dim,1,mpi_double_precision,0,comm, ierror)
!       call mpi_bcast(xRadius_dim,1,mpi_double_precision,0,comm, ierror)
!       call mpi_bcast(yRadius_dim,1,mpi_double_precision,0,comm, ierror)
!       call mpi_bcast(zRadius_dim,1,mpi_double_precision,0,comm, ierror)
!       call mpi_bcast(xCenter_dim,1,mpi_double_precision,0,comm, ierror)
!       call mpi_bcast(yCenter_dim,1,mpi_double_precision,0,comm, ierror)
!       call mpi_bcast(zCenter_dim,1,mpi_double_precision,0,comm, ierror)
       !       call mpi_bcast(zExcentricity,1,mpi_double_precision,0,comm, ierror)

       ! read test case input data   ! modified by Junhong Wei (20161122)
       read (unit=10, nml=wavePacket)   ! modified by Junhong Wei (20161122)

       read (unit=10, nml=bubble)
       

       if (referenceQuantities == "SI" ) stop "initialize: SI units not allowed"

       ! start velocity
!       var(:,:,:,2) = backgroundFlow_dim(1) / uRef
!       var(:,:,:,3) = backgroundFlow_dim(2) / uRef
!       var(:,:,:,4) = backgroundFlow_dim(3) / uRef

       var(:,:,:,2) = 0.0
       var(:,:,:,3) = 0.0
       var(:,:,:,4) = 0.0
       
       ! constant pressure variable pi' 
       var(:,:,:,5) = 0.0

       ! set local index
       i00 = is + nbx - 1

       ! potential temperature and density
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                x_dim = x(i+i00) * lRef       ! dimensional lenghts
                z_dim = z(k) * lRef

                delX = (x_dim - xCenter_dim) / xRadius_dim
                delZ = (z_dim - zCenter_dim) / zRadius_dim

                r = sqrt(delX**2 + delZ**2)  ! scaled radius

                if( r<=1.0 ) then          ! inside bubble

                   dTheta_dim = dTheta0_dim * (cos(pi*r/2.0))**2
                   theta = thetaStrat(k) + dTheta_dim / thetaRef

                   if( fluctuationMode )  then

                      ! calc pseudo-incompressible density rho*
                      if ( referenceQuantities == "SI" ) then
                         rho = p0**kappa/Rsp * Pstrat(k) / theta - rhoStrat(k)
                      else
                         rho = Pstrat(k) / theta - rhoStrat(k)
                      end if

                   else

                      ! calc pseudo-incompressible density rho*
                      if ( referenceQuantities == "SI" ) then
                         rho = p0**kappa/Rsp * Pstrat(k) / theta 
                      else
                         rho = Pstrat(k) / theta 
                      end if

                   end if  ! fluctuation mode

                   select case( model ) 

                   case( "pseudo_incompressible" ) 

                       var(i,j,k,1) = rho

                   case( "Boussinesq" ) 

                       var(i,j,k,1) = rhoStrat(k)
                       var(i,j,k,6) = dTheta_dim / thetaRef

                   case default
                      stop "initialize: unknown model."
                   end select
                    
                   
                else  ! outside bubble
                   
                   select case( model ) 

                   case( "pseudo_incompressible" ) 

                     if( fluctuationMode ) then
                        var(i,j,k,1) = 0.0
                     else
                        var(i,j,k,1) = rhoStrat(k)
                     end if

                   case( "Boussinesq" ) 

                     var(i,j,k,1) = rhoStrat(k)
                     var(i,j,k,6) = 0.0

                   case default
                      stop "initialize: unknown model."
                   end select
                   
                end if
                
             end do
          end do
       end do

       
       ! modified by Junhong Wei (20161122) *** finishing line ***
       
       
       
       !-----------------------------------------------------------------------


    case( 'hotBubble3D' )

       ! read test case input data   ! modified by Junhong Wei (20161122)
       read (unit=10, nml=wavePacket)   ! modified by Junhong Wei (20161122)

       ! read test case input data
       read (unit=10, nml=bubble)

       if (referenceQuantities == "SI" ) stop "initialize: SI units not allowed"

       ! zero start velocity 
       var(:,:,:,2) = 0.0
       var(:,:,:,3) = 0.0
       var(:,:,:,4) = 0.0

       ! constant pressure variable pi' 
       var(:,:,:,5) = 0.0


       ! modified by Junhong Wei (20161122) *** starting line ***
       
       
!       ! potential temperature and density
!       do k = 1,nz
!          do j = 1,ny
!             do i = 1,nx
!                x_dim = x(i) * lRef       ! dimensional lenghts
!                y_dim = y(j) * lRef
!                z_dim = z(k) * lRef
!
!                delX = (x_dim - xCenter_dim) / xRadius_dim
!                delY = (y_dim - xCenter_dim) / xRadius_dim
!                delZ = (z_dim - zCenter_dim) / zRadius_dim
!
!                r = sqrt(delX**2 + delY**2 + delZ**2)  ! scaled radius
!
!                if( r<=1.0 ) then          ! inside bubble
!
!                   dTheta_dim = dTheta0_dim * (cos(pi*r/2.0))**2
!
!                   theta = thetaStrat(k) + dTheta_dim / thetaRef
!
!                   ! calc pseudo-incompressible density rho*
!                   if ( referenceQuantities == "SI" ) then
!                      rho = p0**kappa/Rsp * Pstrat(k) / theta
!                   else
!                      rho = Pstrat(k) / theta
!                   end if
!                   var(i,j,k,1) = rho
!
!                else  ! outside bubble
!                   ! keep background density
!                   var(i,j,k,1) = rhoStrat(k)
!                end if
!
!             end do
!          end do
!       end do


       ! set local index
       i00 = is + nbx - 1
       j00 = js + nby - 1

       ! potential temperature and density
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                x_dim = x(i+i00) * lRef       ! dimensional lengh
                y_dim = y(j+j00) * lRef
                z_dim = z(k) * lRef

                delX = (x_dim - xCenter_dim) / xRadius_dim
                delY = (y_dim - xCenter_dim) / xRadius_dim
                delZ = (z_dim - zCenter_dim) / zRadius_dim

                r = sqrt(delX**2 + delY**2 + delZ**2)  ! scaled radius

                if( r<=1.0 ) then          ! inside bubble

                   dTheta_dim = dTheta0_dim * (cos(pi*r/2.0))**2
                   theta = thetaStrat(k) + dTheta_dim / thetaRef

                   if( fluctuationMode )  then

                      ! calc pseudo-incompressible density rho*
                      if ( referenceQuantities == "SI" ) then
                         rho = p0**kappa/Rsp * Pstrat(k) / theta - rhoStrat(k)
                      else
                         rho = Pstrat(k) / theta - rhoStrat(k)
                      end if

                   else

                      ! calc pseudo-incompressible density rho*
                      if ( referenceQuantities == "SI" ) then
                         rho = p0**kappa/Rsp * Pstrat(k) / theta
                      else
                         rho = Pstrat(k) / theta
                      end if

                   end if  ! fluctuation mode

                   select case( model )

                   case( "pseudo_incompressible" )

                       var(i,j,k,1) = rho

                   case( "Boussinesq" )

                       var(i,j,k,1) = rhoStrat(k)
                       var(i,j,k,6) = dTheta_dim / thetaRef

                   case default
                      stop "initialize: unknown model."
                   end select

                else  ! outside bubble

                   select case( model )

                   case( "pseudo_incompressible" )

                     if( fluctuationMode ) then
                        var(i,j,k,1) = 0.0
                     else
                        var(i,j,k,1) = rhoStrat(k)
                     end if

                   case( "Boussinesq" )

                     var(i,j,k,1) = rhoStrat(k)
                     var(i,j,k,6) = 0.0

                   case default
                      stop "initialize: unknown model."
                   end select

                end if

             end do
          end do
       end do       
       
       ! modified by Junhong Wei (20161122) *** finishing line ***       
       
       
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
       updateIce = .true.

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
       updateIce = .true.

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
       updateIce = .false.

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
       updateIce = .false.

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


!   achatzb
!   -------------------------------------
!   in case of topography, 
!   set all velocities normal to the topographic surface to zero,
!   set density in land cells to background density
!   ------------------------------------

    i0=is+nbx-1
    j0=js+nby-1

    if(topography) then
       do k = 0, nz+1
          do j = 0, ny+1
             do i = 0, nx+1
!               u at x interfaces
                if(&
                   topography_mask(i0+i,j0+j,k)&
                   .or.&
                   topography_mask(i0+i+1,j0+j,k)&
                ) then
                   var(i,j,k,2)=0.
                end if

!               v at y interfaces
                if(&
                   topography_mask(i0+i,j0+j,k)&
                   .or.&
                   topography_mask(i0+i,j0+j+1,k)&
                ) then
                   var(i,j,k,3)=0.
                end if

!               w at z interfaces
                if(&
                   topography_mask(i0+i,j0+j,k)&
                   .or.&
                   topography_mask(i0+i,j0+j,k+1)&
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


    ! close input file pinc.f
    close (unit=10)



    !----------------------------------
    !     Output system settings 
    !----------------------------------

    if( master ) then   ! modified by Junhong Wei (20170216)

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
    write(*,fmt="(a25,es10.3,a,f5.1,a)") "PStrat = ", PStrat(nz)*pRef/1000.0,&
         & " kPa at z = ", z(nz) * lRef/1000.0, " km"
    write(*,fmt="(a25,es10.3,a,f5.1,a)") "rhoStrat = ", rhoStrat(nz)*rhoRef,&
         & " kg/m3 at z = ", z(nz) * lRef/1000.0, " km"
    write(*,fmt="(a25,es10.3,a,f5.1,a)") "thetaStrat = ", thetaStrat(nz)*thetaRef,&
         & " K at z = ", z(nz) * lRef/1000.0, "km"
    print*,""

    print*,"  4) Constants: "
    write(*,fmt="(a25,f7.3,a)") "gamma = ", gamma, " "
    write(*,fmt="(a25,f7.3,a)") "g = ", g, " m/s^2"
    write(*,fmt="(a25,f7.3,a)") "R_sp = ", Rsp, " J/kg/K"
    write(*,fmt="(a25,es10.3,a)") "f_Coriolis = ", f_Coriolis_dim, " 1/s"
    write(*,fmt="(a25,es10.3,a)") "mu_viscous = ", mu_viscous_dim, " m^2/s"
    write(*,fmt="(a25,es10.3,a)") "mu_conduct = ", mu_conduct_dim, " m^2/s"
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
       write(*,fmt="(a25,es8.1,a)") "relaxation  = ", spongeAlphaZ_dim, " 1/s"
    else
       write(*,fmt="(a25,a)") "sponge layer = ", "off"
    end if
    print*,""



    print*,"  9) Poisson Solver: "
    write(*,fmt="(a25,a)") "solver = ", poissonSolverType
    write(*,fmt="(a25,es8.1)") "tolPoisson = ", tolPoisson
!   achatzb
    write(*,fmt="(a25,es8.1)") "tolCond = ", tolCond
!   achatze
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

    end if   ! modified by Junhong Wei (20170216)



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


    !------------------------------------------------------------------------------


!    subroutine init_GWP(Psi,kk,mm)
    subroutine init_GWP(Psi,kk,mm, ll_3DWP )   ! modified by Junhong Wei for 3DWP (20170828)

      !------------------------------------------------
      !  calculate complex amplitudes for
      !    1) first harmonics,  leading order: Psi(:,:,:,0)
      !    2) second harmonics, leading order: Psi(:,:,:,1)
      !------------------------------------------------

      ! in/out variables
      !real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
      !     & intent(in) :: var                         ! mean flow velocities    

      !    type(rayType), dimension(nRay)        :: ray
      !    real, dimension(0:nx+1,0:ny+1,0:nz+1) :: waveAct 

      ! wave amplitude
!      complex, dimension(0:nx+1,0:ny+1,0:nz+1,4,0:2), intent(out) :: Psi ! modified by Junhong Wei
      complex, dimension(0:nx+1,0:ny+1,0:nz+1,5,0:2), intent(out) :: Psi ! modified by Junhong Wei
      real, intent(out) :: kk,mm
      real, intent(out) :: ll_3DWP   ! modified by Junhong Wei for 3DWP (20170828)

      ! local variables
      real :: A, rho
      real :: omi
      real :: kk2,mm2,omi2, kTot2
      complex :: u10,w10,b11,pi12
      integer :: iRay
      integer :: i,j,k
      integer :: ii,jj

      ! local amplitude values and derivatives
      complex :: u10_r, u10_c, u10_l, u10_t, u10_b
      complex :: w10_r, w10_c, w10_l, w10_t, w10_b
      complex :: theta11_r, theta11_c, theta11_l, theta11_t, theta11_b
      complex :: pi12_c 
      complex :: pi0_t, pi0_c, pi0_b
      complex :: theta0_c
      complex :: du10_dx, du10_dz, dw10_dx, dw10_dz 
      complex :: dpi0_dz, dtheta11_dx, dtheta11_dz
      complex :: Div
      complex :: Press, PressU, PressW
      complex :: d1u10, d1w10, d1theta11
      complex :: coeff, aux1
      complex :: M11,M12,M13,M14
      complex :: M21,M22,M23,M24
      complex :: M31,M32,M33,M34
      complex :: M41,M42,M43,M44
      complex, dimension(4) :: RHS, Sol
      complex, dimension(4,4) :: M2, M2inv
      complex :: u21, w21,b22, pi23

      ! mean value calculation
      real :: rho0, rho0_t, rho0_b, d_dz, ypsi


      ! debugging stuff
      complex :: summe
      complex, dimension(4) :: term

      ! local fields
      !    real, dimension(nx,ny,nz) :: omegaMean   ! cell averaged intrinsic freq.
      !    real, dimension(nx,ny,nz) :: kMean, mMean ! cell mean wave numbers

      ! more debugging stuff
      real :: B11_pinc
      real :: D1TH11

    integer :: i0, j0   ! modified by Junhong Wei (20161201)

      complex :: tmp_var_3DWP   ! modified by Junhong Wei for 3DWP (20170921)

! modified by Junhong Wei (20161201) *** starting line ***
    !-----------------------
    !      MPI stuff
    !-----------------------
    i0 = is + nbx - 1   ! 0 index, replace i -> i + i0 in x and y fields
    j0 = js + nby - 1
! modified by Junhong Wei (20161201) *** finishing line ***

      !------------------------
      !    Init data
      !-----------------------

      ! open input file input.f90
      open (unit=10, file="input.f90", action="read", &
           form="formatted", status="old", position="rewind")

      ! read test case input data
      read (unit=10, nml=wavePacket)

      ! scale input data
      lambdaX = lambdaX_dim/lRef     ! non-dim zonal wave length
      lambdaZ = lambdaZ_dim/lRef     !         vert. wave length

      xCenter = xCenter_dim/lRef     ! scaled position of wave packtet
      zCenter = zCenter_dim/lRef 

      sigma = sigma_dim/lRef         ! sigma width of Gaussian distribution
      sigma_hor = sigma_hor_dim/lRef ! sigma width of Gaussian distribution  ! modified by Junhong Wei (20170214)
      L_cos = L_cos_dim/lRef         ! half length of cosine profile

      lambdaY = lambdaY_dim/lRef     !         meridional wave length   ! modified by Junhong Wei for 3DWP (20170828)
      yCenter = yCenter_dim/lRef     ! scaled position of wave packtet  ! modified by Junhong Wei for 3DWP (20170828)
      sigma_hor_yyy = sigma_hor_yyy_dim/lRef ! sigma width of Gaussian distribution  ! modified by Junhong Wei for 3DWP (20170828)

            if( (ABS(lambdaY_dim)) .GT. 0.1 ) then   ! modified by Junhong Wei for 3DWP (20170921)
               ll_3DWP = 2.0*pi/lambdaY     ! modified by Junhong Wei for 3DWP (20170828)
            else                            ! modified by Junhong Wei for 3DWP (20170828)
               ll_3DWP = 0.0                ! modified by Junhong Wei for 3DWP (20170828)
            end if                          ! modified by Junhong Wei for 3DWP (20170828)

      ! wave numbers 
!      kk = 2.0*pi/lambdaX    !xxx  new signs by Ulrich, 14.9.2012   ! modified by Junhong Wei for 3DWP (20171128)

            if( (ABS(lambdaX_dim)) .GT. 0.1 ) then   ! modified by Junhong Wei for 3DWP (20171128)
               kk = 2.0*pi/lambdaX     ! modified by Junhong Wei for 3DWP (20171128)
            else                            ! modified by Junhong Wei for 3DWP (20171128)
               kk = 0.0                ! modified by Junhong Wei for 3DWP (20171128)
            end if                          ! modified by Junhong Wei for 3DWP (20171128)

      mm = 2.0*pi/lambdaZ   
      kk2 = kk**2
      mm2 = mm**2
!      kTot2 = kk2 + mm2   ! modified by Junhong Wei for 3DWP (20170828)
      kTot2 = kk2 + mm2 + ( ll_3DWP * ll_3DWP )   ! modified by Junhong Wei for 3DWP (20170828)
      kTot = sqrt(kTot2)

      ! intrinsic frequency
!      omi = omiSign * sqrt(N2)*kk/kTot  ! modified by Junhong Wei
!      omi = omiSign * sqrt( ( (N2*kk*kk) + (RoInv*RoInv*mm*mm) ) )/kTot  ! modified by Junhong Wei   ! modified by Junhong Wei for 3DWP (20170828)
       omi = omiSign * sqrt( ( ( N2 * ((kk*kk)+(ll_3DWP*ll_3DWP)) ) + (RoInv*RoInv*mm*mm) ) )/kTot  ! modified by Junhong Wei for 3DWP (20170828)
      omi2 = omi**2

      ! amplitude coefficients for wave 1
      bAmp = amplitudeFactor * N2/mm                  ! buoyancy
      uAmp = mm/kk * omi/N2 * bAmp
      wAmp = omi/N2 * bAmp
      pAmp = kappa*Ma2 * mm/kk**2 * omi2/N2 * bAmp    ! Exner pressure


      close(10)

       !----------------------
       !  output of init data
       !----------------------

    if( master ) then   ! modified by Junhong Wei (20170216)

!xxx

print*,"omi = ", omi/tRef
print*,"mm = ", mm/lRef

print*,"RoInv = ", RoInv/tRef   ! modified by Junhong Wei

       print*,""
       print*,"  0) Test case: "
       write(*,fmt="(a25,a35)") "Test case  = ", "wave packet (full model)"
       write(*,fmt="(a25,f10.1,a)") "lambda_x = ", lambdaX_dim, " m"
       write(*,fmt="(a25,f10.1,a)") "lambda_z = ", lambdaZ_dim, " m"
       write(*,fmt="(a25,f10.1a7)") "c_x  = ", omi/kk*uRef, " m/s"
       write(*,fmt="(a25,f10.1,a7)") "c_z  = ", omi/mm*uRef, " m/s"
       ! new sign by Ulrich Achatz, 14.9.2012
       write(*,fmt="(a25,f10.1,a7)") "cg_x  = ", -NN*mm**2/kTot**3 * uRef, " m/s"
       write(*,fmt="(a25,f10.1,a7)") "cg_z  = ", NN*mm*kk/kTot**3 * uRef, " m/s"
       write(*,fmt="(a25,f10.1,a7)") "u_jet  = ", u0_jet_dim, " m/s"
       print*,""

    end if   ! modified by Junhong Wei (20170216)



      !---------------------------------------
      !        calc amplitude Psi_1^0 
      !     (first harmonic, leading order)
      !---------------------------------------

!      do k = 1,nz
      do k = 0,(nz+1)   ! modified by Junhong Wei for 3DWP (20171204)
!         j = 1   ! modified by Junhong Wei for 3DWP (20170921)
          do j = 0,(ny+1)   ! modified by Junhong Wei for 3DWP (20170921)
!         do i = 1,nx ! modified by Junhong Wei (20161207)
         do i = 0,(nx+1) ! modified by Junhong Wei (20161207)


            ! profile: 1D and 2D
            if( wavePacketDim == 1 ) then
               delx = 0.0
            else
!               delx = (x(i)-xCenter)   ! modified by Junhong Wei (20161201)
               delx = ( x(i+i0) -xCenter)   ! modified by Junhong Wei (20161201)
            end if

            if( wavePacketDim == 3 ) then
               dely = ( y(j+j0) -yCenter)
            else
               dely = 0.0
            end if

            delz = (z(k)-zCenter)

            select case(wavePacketType) 

            case(1)

               ! Gaussian
!               envel = exp(-(delx**2 + delz**2)/2./sigma**2)   ! modified by Junhong Wei (20170214)
!               envel = ( exp(-(delz**2)/2./sigma**2) ) * ( exp(-(delx**2)/2./sigma_hor**2) )   ! modified by Junhong Wei for 3DWP (20170214)

!           achatzb cosine profile horizontally so that fields are zero 
!           at the horizontal boundaries
!           in case of zero sigma in x or y direction use infinity
!           envel &
!           = ( exp(-(delz**2)/2./sigma**2) ) &
!           * ( exp(-(delx**2)/2./sigma_hor**2) ) &
!           * ( exp(-(dely**2)/2./sigma_hor_yyy**2) )   &
!           ! modified by Junhong Wei for 3DWP (20170921)

            if(sigma_hor == 0.0) then
               envel = 1.0
              else if(abs(delx) < sigma_hor) then
               envel &
               = 1.0 - amp_mod_x + amp_mod_x *cos(delx*pi/(sigma_hor*2.0))
              else
               envel = 1.0 - amp_mod_x
            end if

            if(sigma_hor_yyy == 0.0) then
               envel = 1.0 * envel
              else if(abs(dely) < sigma_hor_yyy) then
               envel &
               = (1.0 - amp_mod_y &
                  + amp_mod_y * cos(dely*pi/(sigma_hor_yyy*2.0))) &
                 * envel
              else
               envel = envel * (1.0 - amp_mod_y)
            end if

            envel = envel * exp(-(delz**2)/2./sigma**2)
!           achatze

            case(2) 

               ! Cosine
               if( abs(delz) .le. L_cos ) then
                  envel = 0.5*(1.0 + cos(pi*delz/L_cos))
               else
                  envel = 0.0
               end if

            case default
               stop "init.f90: unknown wavePacketType. Stop."
            end select



            b11 = cmplx(envel*bAmp, 0.0 )

!           Modified by Junhong Wei

            theta0 = thetaStrat(k)

            tmp_var_3DWP = cmplx( 0.0,  ((omi*omi)-N2) / ( mm*N2*( (omi*omi)-(RoInv*RoInv) ) )    )   ! modified by Junhong Wei for 3DWP (20170921)

!            u10 = cmplx(0.0, -mm/kk * omi/N2) * b11
!            u10 = cmplx(0.0, ((omi*omi)-N2)*kk*omi / ( mm*N2*( (omi*omi)-(RoInv*RoInv) ) )  ) * b11   ! modified by Junhong Wei for 3DWP (20170921)

            u10 = tmp_var_3DWP * cmplx( kk*omi, ll_3DWP*RoInv ) * b11   ! modified by Junhong Wei for 3DWP (20170921)

            w10 = cmplx(0.0, omi/N2) * b11

!            pi12 = cmplx(0.0, -kappa*Ma2* mm/kk**2 * omi**2/N2 / theta0) * b11
            pi12 = cmplx(0.0, kappa*Ma2*( (omi*omi)-N2 ) / N2 / mm / theta0) * b11

!            Psi(i,j,k,:,1) = (/u10, w10, b11, pi12/) ! modified by Junhong Wei

!            Psi(i,j,k,:,1) = (/u10, w10, b11, pi12, ( cmplx( 0.0, 0.0 ) * b11 ) /) ! modified by Junhong Wei
!            Psi(i,j,k,:,1) = (/u10, w10, b11, pi12, ( cmplx( ((omi*omi)-N2)*kk*RoInv / ( mm*N2*( (omi*omi)-(RoInv*RoInv) ) )  , 0.0 ) * b11 ) /) ! modified by Junhong Wei   ! modified by Junhong Wei for 3DWP (20170921)


            Psi(i,j,k,:,1) = (/u10, w10, b11, pi12, ( tmp_var_3DWP * cmplx( ll_3DWP*omi, kk*RoInv*(0.0-1.0) ) * b11 ) /) ! modified by Junhong Wei for 3DWP (20170921)

!           Modified by Junhong Wei

         end do
       end do   ! modified by Junhong Wei for 3DWP (20170921)
      end do


! modified by Junhong Wei (20161207) *** starting line ***
!
!      !---------------------------------------
!      !    set ghost cell values for Psi_1^0 
!      !    x/y: Periodic, z: solid wall 
!      !---------------------------------------
!
!      ! periodic in x
!      Psi(0,:,:,:,1) = Psi(nx,:,:,:,1)
!      Psi(nx+1,:,:,:,1) = Psi(1,:,:,:,1)
!
!      ! periodic in y
!      ! implement for 3D
!
!      ! solid wall -> reflect u10 and w10 with change of sign
!      Psi(:,:,nz+1,2,1) = -Psi(:,:,nz,2,1)
!      Psi(:,:,0,2,1) = -Psi(:,:,1,2,1)
!
!      ! solid wall: reflect b11 and pi12 without change of sign 
!      Psi(:,:,nz+1,3:4,1) = Psi(:,:,nz,3:4,1)
!      Psi(:,:,0,3:4,1) = Psi(:,:,1,3:4,1)
!
! modified by Junhong Wei (20161207) *** finishing line ***


      !---------------------------------------
      !        calc amplitude Psi_2^1 
      !     (second harmonic, first order)
      !---------------------------------------

      do k = 1,nz
         j = 1
         do i = 1,nx

            !-------------------------
            !       Set up RHS
            !-------------------------

            ! zonal velocities right, center, left, top, bottom
            u10_r = Psi(i+1,j, k ,1,1)
            u10_c = Psi( i ,j, k ,1,1)
            u10_l = Psi(i-1,j, k ,1,1)
            u10_t = Psi( i ,j,k+1,1,1)
            u10_b = Psi( i ,j,k-1,1,1)

            ! vertical velocities top, center, bottom
            w10_r = Psi(i+1,j, k ,2,1)
            w10_l = Psi(i-1,j, k ,2,1)          
            w10_t = Psi( i ,j,k+1,2,1)
            w10_c = Psi( i ,j, k ,2,1)
            w10_b = Psi( i ,j,k-1,2,1)


            ! Buoyancy and pot. temp.
            theta11_r = Fr2 * thetaStrat(k) * Psi(i+1,j,k,3,1)
            theta11_c = Fr2 * thetaStrat(k) * Psi(i,j,k,3,1)
            theta11_l = Fr2 * thetaStrat(k) * Psi(i-1,j,k,3,1)
            theta11_t = Fr2 * thetaStrat(k+1) * Psi(i,j,k+1,3,1)
            theta11_b = Fr2 * thetaStrat(k-1) * Psi(i,j,k-1,3,1)

            ! Second order Exner pressure
            pi12_c = Psi(i,j,k,4,1)

            ! Background Exner pressure and pot. temp.
            pi0_t = (Pstrat(k+1)/p0)**gamma_1
            pi0_c = (Pstrat(k)/p0)**gamma_1
            pi0_b = (Pstrat(k-1)/p0)**gamma_1
            theta0_c = thetaStrat(k)


            ! derivatives          
            du10_dx = 0.5*(u10_r - u10_l)/dx
            du10_dz = 0.5*(u10_t - u10_b)/dz
            dw10_dx = 0.5*(w10_r - w10_l)/dx
            dw10_dz = 0.5*(w10_t - w10_b)/dz
            dpi0_dz = 0.5*(pi0_t - pi0_b)/dz
            dtheta11_dx = 0.5*(theta11_r - theta11_l)/dx
            dtheta11_dz = 0.5*(theta11_t - theta11_b)/dz

            ! divergence term -> Eq. (7.21)
            Div = -du10_dx - dw10_dz - (1.-kappa)/kappa*w10_c/pi0_c*dpi0_dz

            ! Pressure terms
            Press = 0.5*kappaInv*MaInv2*imag*theta11_c*pi12_c
            PressU = kk*Press
            PressW = mm*Press


            ! intermediate terms
            d1u10 = 0.5*(u10_c*du10_dx + w10_c*du10_dz + Div*u10_c)
            d1w10 = 0.5*(u10_c*dw10_dx + w10_c*dw10_dz + Div*w10_c)
            d1theta11 = 0.5*(u10_c*dtheta11_dx + w10_c*dtheta11_dz + Div*theta11_c)

            RHS(1) = -d1u10 - PressU
            RHS(2) = -d1w10 - PressW
            RHS(3) = -FrInv2/NN/theta0_c * d1theta11
            RHS(4) = (0.0, 0.0)


            !----------------------------------------------------
            !       Set up inverted system matrix M(2om,2k,2m)
            !----------------------------------------------------

            coeff = 1./(4.*omi2*kTot2 - kk2*N2)
            aux1 = 4.*omi2-N2

            M11 = 2.*imag*mm2*omi; M12 = -2.*imag*kk*mm*omi; M13 = kk*mm*NN; M14 = -0.5*imag*kk*aux1
            M21 = M12;           M22 = 2.*imag*kk2*omi;    M23 = -kk2*NN;  M24 = -2.*imag*omi2*mm
            M31 = -M13;          M32 = -M23;             M33 = 2.*imag*omi*kTot2; M34 = -omi*mm*NN
            M41 = M14;           M42 = M24;              M43 = -M34;     M44 = -0.5*imag*omi*aux1


            ! inverted matrix: check ok
            M2inv(1,:) = (/M11,M12,M13,M14/)
            M2inv(2,:) = (/M21,M22,M23,M24/)
            M2inv(3,:) = (/M31,M32,M33,M34/)
            M2inv(4,:) = (/M41,M42,M43,M44/)
            M2inv = coeff*M2inv


            !---------------------------------------
            !   Solve linear System -> save in Psi
            !---------------------------------------

            sol = matmul(M2inv,RHS)

            u21 = sol(1)
            w21 = sol(2)
            b22 = sol(3) * NN
            pi23 = sol(4) * kappa*Ma2 / theta0_c

!            Psi(i,j,k,:,2) = (/u21,w21,b22,pi23/) ! modified by Junhong Wei

!            Psi(i,j,k,:,2) = (/u21,w21,b22,pi23,u21/) ! modified by Junhong Wei

            Psi(i,j,k,:,2) = (/u21,w21,b22,pi23, ( cmplx( 0.0, 0.0 ) * b11 ) /) ! modified by Junhong Wei

         end do
      end do



    end subroutine init_GWP


  end subroutine initialise


  !--------------------------------------------------------------------------    


  function cphase(c) result(phi)

    !------------------------------------
    !  phase of complex number
    !------------------------------------

    ! in/out
    complex, intent(in) :: c
    real                :: phi

    ! local vars
    real :: a,b

    a = real(c)
    b = aimag(c)

    ! first quadrant
    if( a>=0 .and. b>=0) then
       phi = atan(b/a)

       ! second quadrant
    else if( a<0. .and. b>=0. ) then
       phi = pi - atan(-b/a)

       ! third quadrant
    else if( a<0. .and. b<0.) then
       phi = -pi + atan(b/a)

       ! fourth quadrant
    else if( a>=0. .and. b<0 ) then
       phi = -atan(-b/a)
    else
       stop "wkb.f90/cphase: case not included. Stop."
    end if


    phi = phi*180./pi


  end function cphase




end module init_module
