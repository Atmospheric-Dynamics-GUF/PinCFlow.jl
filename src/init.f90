module init_module


  use type_module
  use atmosphere_module
  use ice_module
  use mpi_module
  use boundary_module

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

    ! open the namelist file
    open (unit=10, file=file_namelist, action="read", &
         form="formatted", status="old", position="rewind")

    ! read parameter study list
    read (unit=10, nml=parameterList)

    close( unit=10 ) 


  end subroutine init_paramStudy



  ! ----------------------------------------------------------------------


  !UAB
  !subrouine setup (var,var0,var1,flux,force,source,dRho,dRhop,dMom,dTheta)
  subroutine setup (var,var0,var1,flux,flux_rhopw,force,source,dRho,dRhop,dMom, &
                  & dTheta,dPStrat,drhoStrat,w_0,dIce)
  
  !UAE
    !-----------------------------------------
    ! allocate var and flux / read the namelist
    !-----------------------------------------

    ! in/out variables 
    real,dimension(:,:,:,:), allocatable,intent(out) :: var, var0, var1, &
                                                      & source
    real,dimension(:,:,:,:,:), allocatable,intent(out) :: flux
    real, dimension(:,:,:), allocatable, intent(out) :: flux_rhopw
    real,dimension(:,:,:,:), allocatable,intent(out) :: force
    real,dimension(:,:,:), allocatable :: dRho,dRhop   ! RK-Update for rho
    real,dimension(:,:,:,:), allocatable :: dMom ! ...rhoU,rhoV,rhoW
    real,dimension(:,:,:), allocatable :: dTheta ! RK-Update for theta
    real,dimension(:,:,:,:), allocatable :: dIce       ! RK-Update for nIce,qIce,qAer,qv

    !UAB
    real, dimension(:), allocatable :: dPStrat, drhoStrat ! RK-Update for P
    real, dimension(:), allocatable :: w_0 !! w_0 from ONK14
    !UAE

    integer :: allocstat
    integer :: i, j, k, iVar

    ! constants
    pi = 4*atan(1.0)

    ! set default values of all namelist parameters
    call default_values

    ! open input file input.f90
    open (unit=10, file=file_namelist, action="read", &
         form="formatted", status="old", position="rewind")

    !---------------------------------------------------
    ! allocate x,y,z - cell centered coordinate fields
    !---------------------------------------------------
  
    allocate(x(-nbx:sizeX+nbx), stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate x"

    allocate(y(-nby:sizeY+nby), stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate y"

    allocate(z(-nbz:sizeZ+nbz), stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate z"

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


    !-------------------------------------
    !      allocate variable fields
    !-------------------------------------

    read (unit=10, nml=variables)
    if (include_ice) nVar = nVar+4

    ! allocate var = (rho,u,v,w,pEx)
    allocate(var(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate var"

    ! allocate var0 = (rho,u,v,w,pEx)
    allocate(var0(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate var0"

    ! allocate var1 = (rho,u,v,w,pEx)
    allocate(var1(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate var1"

    ! allocate source for rho,u,v,w,pEx,theta
    allocate(source(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate source"

    ! allocate dRho
    allocate(dRho(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dRho."

    ! allocate dRhop
    allocate(dRhop(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dRhop."

    ! allocate dMom
    allocate(dMom(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,3),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dMom."

    ! allocate dTheta
    allocate(dTheta(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dTheta."

    ! allocate dIce
    if (include_ice) then
      allocate(dIce(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,4),stat=allocstat)
      if(allocstat /= 0) stop "init.f90: Could not allocate dIce."
    end if

    !UAB
    ! allocate dPStrat
    allocate(dPStrat(-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate dPStrat."

    ! allocate drhoStrat
    allocate(drhoStrat(-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: Could not allocate drhoStrat."

    ! allocate w_0
    allocate(w_0(-nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate w_0"
    !UAE

    ! allocate flux = (f,g,h / fRho, fRhoU, fRhoV, fRhoW, fTheta)
    allocate(flux(-1:nx,-1:ny,-1:nz,3,nVar),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate flux"

    ! allocate flux_rhopw 
    allocate(flux_rhopw(-1:nx,-1:ny,-1:nz),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate flux_rhopw"

    ! allocate force 
    allocate(force(0:nx+1,0:ny+1,0:nz+1,3),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate force"

    ! allocate environm. pot. temp. perturbations
    !allocate(the_env(-nbx+1:sizeX+nbx,-nby+1:sizeY+nby,-nbz+1:sizeZ+nbz), stat=allocstat)
    !if(allocstat /= 0) stop "init.f90: could not allocate the_env"

    !allocate(dens_env(-nbx+1:sizeX+nbx,-nby+1:sizeY+nby,-nbz+1:sizeZ+nbz), stat=allocstat)
    !if(allocstat /= 0) stop "init.f90: could not allocate dens_env"
    
    !allocate(u_env(-nbx+1:sizeX+nbx,-nby+1:sizeY+nby,-nbz+1:sizeZ+nbz), stat=allocstat)
    !if(allocstat /= 0) stop "init.f90: could not allocate u_env"
    
    
    allocate(p_env_pp(1:nx,1:ny,0:nz+1),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate p_env_pp"
    
    allocate(the_env_pp(1:nx,1:ny,0:nz+1),stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate the_env_pp"
    
    allocate(dens_env_pp(1:nx,1:ny,0:nz+1), stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate dens_env_pp"
    
    allocate(u_env_pp(0:nx,1:ny,0:nz+1), stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate u_env"
    
    allocate(v_env_pp(1:nx,0:ny,0:nz+1), stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate v_env"
    
    allocate(u_const(0:nx,1:ny), stat=allocstat)
    if(allocstat /= 0) stop "init.f90: could not allocate u_const"
    
    !UAb
    !if (TestCase == "baroclinic_LC") then
    !   allocate(var_env(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
    !           & stat=allocstat)
    !   if(allocstat /= 0) stop "init.f90: could not allocate var_env"
    !end if
    !UAE

    !-------------------------------------
    !    read name lists
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
      forall (i=0:3) 
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

    if (include_ice) then
      ! read ice physics parametrization
      read(unit=10, nml=iceLIst)
    end if 

    ! close input file pinc.f
    close (unit=10)

    ! decide if heating of the reference atmosphere is applied or not
    heatingRefAtmo = (heatingONK14 .or. TurbScheme .or. rayTracer)

    !UAB
    if (TestCase == "baroclinic_LC") then
       allocate(var_env(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
               & stat=allocstat)
       if(allocstat /= 0) stop "init.f90: could not allocate var_env"

       if (background == "HeldSuarez") then
          allocate(kt_hs(1:ny,0:nz+1), stat=allocstat)
          if(allocstat /= 0) stop "init.f90: could not allocate kt_hs"

          allocate(kv_hs(0:nz+1), stat=allocstat)
          if(allocstat /= 0) stop "init.f90: could not allocate kv_hs"
       end if
    end if
    !UAE


    !---------------------------------------
    !        Model equation settings
    !---------------------------------------

    ! open info file
    open (unit=90, file="info.txt", action="write", &
         form="formatted", status="replace")

    ! write model info file
    write(90,fmt="(a)"), ""
    write(90,fmt="(a)"), "Model: "//trim(model)
    write(90,fmt="(a)"), ""


    select case( model ) 

    case( "anelastic" )

       pressureScaling = .false.


    case( "Boussinesq" )

       updateMass = .false.
       predictMomentum = .true.
       correctMomentum = .true.
       updateTheta = .true.
       updateIce = .true.

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
       write(90,*), updateMass
       write(90,"(a25)",advance = "no") "predictMomentum != "
       write(90,*), predictMomentum
       write(90,"(a25)",advance = "no") "correctMomentum != "
       write(90,*), correctMomentum
       write(90,"(a25)",advance = "no") "updateTheta != "
       write(90,*), updateTheta
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

    case default
       print*,"model = ", model
       stop "initialize: Unknown model" 
    end select

    close(90)  ! info file

  end subroutine setup


  ! --------------------------------------------------------------------


  subroutine initialise (var)
    !------------------
    ! setup test cases
    !------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var

    ! local variables
    integer :: i,j,k

    integer :: i0, j0, k_test, k_1, k_2    ! modified by Junhong Wei (20161121)

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
    real :: dTheta_dim! , dTheta0_dim

    ! wavepacket: all quantities are scaled
    real :: kk, mm, kTot            ! vertical, zonal, total wave number
    real :: ll                      ! meridional wave number
    real :: omi, omi2               ! intrinsic frequency, squared 
    real :: bAmp, uAmp, wAmp, pAmp  ! ampl. for buoyancy, u, w, Exner pr.
    real :: delx2, delz2            ! squared distance from center
    real :: envel                   ! envelope of wave packet
    real :: Gauss                   ! Gaussian distribution value
    real :: sigma_z                 ! vert. width of Gaussian distribution
    real :: sigma_x                 ! hor. width Gaussian distr. (x dir.)
    real :: sigma_y                 ! hor. width Gaussian distr. (y dir.)
    real :: L_cos                   ! width of Cosine distribution
    real :: xCenter, yCenter        ! center of distribution (hor.)
    real :: zCenter                 ! center of distribution (vert.)
    real :: phi                     ! local phase 
    real :: u,v,w,b                 ! buoyancy, velocities
    real :: lambdaX, lambdaY        ! hor. wave lengths
    real :: lambdaZ                 ! vert. wave lengths


    complex, dimension(0:nx+1,0:ny+1,0:nz+1,5,0:2)  :: Psi
    real :: u1,w1,b1,p1
    real :: u2,w2,b2,p2

    logical, parameter :: initWave2 = .false.

    ! monochromatic wave
    real :: th, f, f2
    real :: lambda
    real :: amp, cot
    real :: uRot, vRot, wRot

    ! testing with Stefan
    real :: xx

    ! random noise on background
    real, dimension(0:nx+1,0:ny+1,0:nz+1) :: randNoise
    real, parameter           :: randAmp = 0.0

    ! jet stream 
    real :: u_jet                        ! jet stream velocity
    real :: L_jet
    real :: z0_jet

    ! baroclinic life cycle
    real  :: term_B1, term_B2, term_A1, term_A2, term_P1, term_P2
    real  :: alpha_t, alpha_s, z_diff 
    real  :: F_1, F_2, F_a, dFadz, z_1, z_2, smooth_diff, L_z
    real  :: T_c, term_a, term_b, DgDy, dbdy, dT0dy, dTcdy, dZdy, dady
    real  :: dpdy, dpsidy, T0_s, DgDy_tr, Ly_end, Lx_end, bar_sigma_x
    real  :: tanhyP1, tanhyP2, ztdif, signmtanh_P1, signmtanh_P2
    real  :: p_dim, tmp_dim, pressure, hor_wind, the_dim, rhe_dim, the, &
             ampl_rho
    real, dimension(1:nx,1:ny) :: P_ref, dPrefdy, Z_trop, T0_t
    real, dimension(1:nx,1:ny) :: gamma_tr, gamma_st
    real, dimension(1:nx,1:ny) :: dpdy_trop, P_trop
    real, dimension(0:nz+1) :: F_0, dF0dz, rho_n
    real:: dthdz, streamfunc, cp, dexndy, dexndz, ex_pr_pert, dstdz
    real  :: theta_bar_0, noise_mag
    real, dimension(1:nx,1:ny,1:nz) :: noise
    real, dimension(0:ny+1,0:nz+1) :: pi_pr_yz
    real, dimension(0:nx+1,0:nz+1) :: pi_pr_xz
    real, dimension(1:nx,1:ny,1:nz-1) :: br_vais_sq, balance, balance1, &
             balance2,balance3, balance4
    character(len=20) :: noise_var, defth

    integer :: i00, j00   ! modified by Junhong Wei (20161121)

    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar) :: flux

    real :: rho_int_m0, rho_int_00, rho_int_mp, rho_int_0p, &
            rho_int_0m, rho_int_pm, rho_int_p0

    real :: ymin, ymax, z_trpp0, deltht, thet0, ntrp, nstr, jwdth
    real :: z_trpp, fstrpp, fstht, xloc, yloc, yjet0, yjet, zloc, z_baro
    real :: zeta_baro, buoy0, buoy1, rho1, rhop0, rhop1
    real :: rptb, thtptb, rhotot

    real :: ptptb_x, ptptb_y, ptptb_z, ptptb_dh, ptptb_dz, ptptb_amp

    real, dimension(1:ny) :: s2_strtd, c2_strtd, c4_strtd
    real :: yjets, yjetn, dy_hs, tempev, pistar, thetastar, sig_pr, facsig
    real :: ka_hs, ks_hs, kf_hs

    ! open the namelist file
    open (unit=10, file=file_namelist, action="read", &
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

    !-----------------------
    !      MPI stuff
    !-----------------------
    i0 = is + nbx - 1   ! 0 index, replace i -> i + i0 in x and y fields
    j0 = js + nby - 1
    

! on default there is no initial ice, humidity or aerosols in the atmosphere
       if (include_ice) var(:,:,:,nVar-3:nVar) = 0.0
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


       ! initialize the ice variables according to iceTestcase
       if (include_ice) call setup_ice(var)


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

       !----------------------------------------------------------
       !          Gravity Waves: wave resolving simulations or WKB
       !----------------------------------------------------------

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
       call init_GWP(Psi,kk,mm,ll)

       do k = 0,(nz+1)
          do j = 0,(ny+1)
             do i = 0,(nx+1)
                phi = kk*x(i+i0) + mm*z(k) + ll*y(j+j0)
             
                ! wave 1
                u1 = real( Psi(i,j,k,1,1) * exp(phi*imag) )
                w1 = real( Psi(i,j,k,2,1) * exp(phi*imag) ) 
                b1 = real( Psi(i,j,k,3,1) * exp(phi*imag) ) 
                p1 = real( Psi(i,j,k,4,1) * exp(phi*imag) ) 


                ! wave 2
                if( initWave2 ) then
                   stop 'ERROR: 2ndary wave not ready for 2D or 3D wave p.'
                   u2 = real( Psi(i,j,k,1,2) * exp(2.*phi*imag) )
                   w2 = real( Psi(i,j,k,2,2) * exp(2.*phi*imag) ) 
                   b2 = real( Psi(i,j,k,3,2) * exp(2.*phi*imag) ) 
                   p2 = real( Psi(i,j,k,4,2) * exp(2.*phi*imag) ) 
                end if


                ! sum of wave 1 and 2
                if( initWave2 ) then
                   stop 'ERROR: 2ndary wave not ready for 2D or 3D wave p.'
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

       ! initialize the ice variables according to iceTestcase
       if (include_ice) call setup_ice(var)

!---------------------------------------------------------------


    case( 'mountainwave' )  
       ! for wave resolving simulation of mountain waves:
       ! read parameters for temporary wind relaxation
       ! zero-wind initial state for montain-wave simulations 

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
       
       
       ! initialize the ice variables according to iceTestcase
       if (include_ice) call setup_ice(var)
       


       !-----------------------------------------------------------------

    case( 'raytracer' )
       ! WKB simulations: Wave packet or mountain waves
       ! for the full set up see routine setup_wkb

       if (.not. raytracer) stop 'raytracer not set correctly'

       ! read namelist for wkb ray tracer
       read (unit=10, nml=LagrangeRayTracing)

       if (sizeX == 1) then
          print*,'sizeX = 1, hence fac_dk_init = 0'
          fac_dk_init = 0.0
       end if

       if (sizeY == 1) then
          print*,'sizeY = 1, hence fac_dl_init = 0'
          fac_dl_init = 0.0
       end if

       if (fac_dk_init == 0.0) then
          dk_init = 0.0
         else if (wlrx_init /= 0.0) then
          dk_init = fac_dk_init * 2.0*pi/wlrx_init
         else if (wlry_init /= 0.0) then
          dk_init = fac_dk_init * 2.0*pi/wlry_init
         else
          stop 'ERROR: BOTH WLRX_INIT and WLRY_INIT = 0.0'
       end if

       if (fac_dl_init == 0.0) then
          dl_init = 0.0
         else if (wlry_init /= 0.0) then
          dl_init = fac_dl_init * 2.0*pi/wlry_init
         else if (wlrx_init /= 0.0) then
          dl_init = fac_dl_init * 2.0*pi/wlrx_init
         else
          stop 'ERROR: BOTH WLRX_INIT and WLRY_INIT = 0.0'
       end if

       zmin_wkb = zmin_wkb_dim/lRef

       ! in WKB mountain-wave case read parameters for wind relaxation

       if (case_wkb == 3) then
          read (unit=10, nml=mountainwavelist)

          ! nondimensionalization

          u_relax = u_relax/uRef

          t_relax = t_relax/tRef
          t_ramp = t_ramp/tRef

          xextent_norelax = xextent_norelax/lRef

          ! increase relaxation wind u_relax so that u = u_relax after the
          ! relaxation period (in x-independent case without topography)

          u_relax = u_relax/(1.0 - exp(4.0*t_ramp/(pi*t_relax) - 1.0))
       end if

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


       !------------------------------------------------------------------


    case( 'coldBubble' )

       ! read test case input data
       read (unit=10, nml=bubble)

       if (referenceQuantities == "SI" ) then
          stop "initialize: SI units not allowed"
       end if

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
                   dTheta_dim &
                   = 0.5*dTheta0_dim * (1.0 + (cos(pi*r/2.0))**2)
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


       !-----------------------------------------------------------------


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
                            rho &
                            = p0**kappa/Rsp * Pstrat(k) / theta &
                              - rhoStrat(k)
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

       !------------------------------------------------------------------



    case( 'hotBubble2D' )

       ! read test case input data 
       read (unit=10, nml=wavePacket)

       read (unit=10, nml=bubble)
       

       if (referenceQuantities == "SI" ) then
          stop "initialize: SI units not allowed"
       end if

       ! start velocity
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
                         rho &
                         = p0**kappa/Rsp * Pstrat(k) / theta - rhoStrat(k)
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

       
       !-------------------------------------------------------------------

    case( 'hotBubble3D' )

       ! read test case input data 
       read (unit=10, nml=wavePacket)

       ! read test case input data
       read (unit=10, nml=bubble)

       if (referenceQuantities == "SI" ) then
          stop "initialize: SI units not allowed"
       end if

       ! zero start velocity 
       var(:,:,:,2) = 0.0
       var(:,:,:,3) = 0.0
       var(:,:,:,4) = 0.0

       ! constant pressure variable pi' 
       var(:,:,:,5) = 0.0

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
                         rho &
                         = p0**kappa/Rsp * Pstrat(k) / theta - rhoStrat(k)
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
       
       !----------------------------------------------------------------
       !               Baroclinic life cycle: realistic
       !
       ! either setup from Kuehnlein et al (2012): background = 'const-N'
       ! or setup of Held & Suarez (1994): background = 'HeldSuarez'
       !----------------------------------------------------------------

   case( 'baroclinic_LC' )

       !UAB
       if ((background /= "const-N") .and. (background /= "HeldSuarez")) &
       & then
          stop 'ERROR: baroclinic_LC needs for background either const-N &
             & or HeldSuarez'
       end if
       !UAE

       ! read test case input data   
       read (unit=10, nml=baroclinic_LC) 
!*****************************************************************
! initialize fields
!*****************************************************************
       var(:,:,:,1) = 0.0
       var(:,:,:,2) = 0.0
       var(:,:,:,3) = 0.0
       var(:,:,:,4) = 0.0
       var(:,:,:,5) = 0.0
       var(:,:,:,6) = 0.0

        the_env_pp(:,:,:) = 0.0
       dens_env_pp(:,:,:) = 0.0
          p_env_pp(:,:,:) = 0.0

       u_env_pp(:,:,:) = 0.0
       v_env_pp(:,:,:) = 0.0

       ! non-dimensional quantities

       ymin = ly_dim(0)/lRef
       ymax = ly_dim(1)/lRef

       jwdth = jwdth_dim/lRef       ! jet width

       !UAB
       if (background == "const-N") then
       !UAE
          z_trpp0 = z_trpp0_dim/lRef ! mean tropopause height

          z_baro = z_baro_dim/lRef   ! alt. above which the atmosph. is 
                                     ! barotropic

          deltht = 3.e1/thetaRef     ! merid. potential-temp. contrast

          thet0 = thet0_dim/thetaRef ! characteristic potential temperature

          ntrp = ntrp_dim*tRef       ! Brunt-Vaisala frequency troposphere
          nstr = nstr_dim*tRef       ! Brunt-Vaisala frequency stratosphere
       !UAB
         else if (background == "HeldSuarez") then
          yjets = 0.75*ymin + 0.25*ymax
          yjetn = 0.25*ymin + 0.75*ymax

          dy_hs = 2.0*jwdth
         else
          stop 'ERROR: baroclinic_LC needs for background either const-N &
             & or HeldSuarez'
       end if
       !UAE

       if (master .and. jwdth > 0.5*(ymax-ymin)) then
          stop 'ERROR: jet width too large'
       end if

       if (sizeY <= 1) stop 'ERROR: Barocl LC expects sizeY > 1'

       ! set local index

       i00 = is + nbx - 1   ! 0 index, 
                            ! replace i -> i + i0 in x and y fields
       j00 = js + nby - 1

       !UAB
       if (background == "const-N") then
       !UAE
          ! potential-temperature field

          do j = 1, ny
             yloc = y(j+j00)

             ! tropopause height

             if (yloc < 0.5*(ymax+ymin)) then
                yjet0 = 0.75*ymin + 0.25*ymax

                if (yloc - yjet0 < -jwdth) then
                   fstrpp = -1.e0
                  else if (      yloc - yjet0 >= -jwdth &
                          &.and. yloc - yjet0 < jwdth) then
                   fstrpp = sin(0.5*pi*(yloc - yjet0)/jwdth)
                  else if (yloc - yjet0 >= jwdth) then
                   fstrpp = 1.e0
                end if
               else
                yjet0 = 0.25*ymin + 0.75*ymax

                if (yloc - yjet0 < -jwdth) then
                   fstrpp = 1.e0
                  else if (      yloc - yjet0 >= -jwdth &
                          &.and. yloc - yjet0 < jwdth) then
                   fstrpp = -sin(0.5*pi*(yloc - yjet0)/jwdth)
                  else if (yloc - yjet0 >= jwdth) then
                   fstrpp = -1.e0
                end if
             end if

             z_trpp &
             =   z_trpp0 &
               + g_ndim*deltht/(2.0*thet0 * (nstr**2 - ntrp**2)) * fstrpp

             do k = 0, nz+1
                zloc = z(k)

                if (zloc < z_trpp) then
                   ! below tropopause

                   if (yloc < 0.5*(ymax+ymin)) then
                      yjet = yjet0 - kaptpp*zloc

                      if (yloc - yjet < -jwdth) then
                         fstht = -1.e0
                        else if (      yloc - yjet >= -jwdth &
                                &.and. yloc - yjet < jwdth) then
                         fstht = sin(0.5*pi*(yloc - yjet)/jwdth)
                        else if (yloc - yjet >= jwdth) then
                         fstht = 1.e0
                      end if
                     else
                      yjet = yjet0 + kaptpp*zloc

                      if (yloc - yjet < -jwdth) then
                         fstht = 1.e0
                        else if (      yloc - yjet >= -jwdth &
                                &.and. yloc - yjet < jwdth) then
                         fstht = -sin(0.5*pi*(yloc - yjet)/jwdth)
                        else if (yloc - yjet >= jwdth) then
                         fstht = -1.e0
                      end if
                   end if

                   the_env_pp(:,j,k) &
                   =   thet0 * (1.e0 + ntrp**2/g_ndim *zloc) &
                     + 0.5*deltht * fstht
                  else
                   ! above tropopause

                   if (yloc < 0.5*(ymax+ymin)) then
                      yjet = yjet0 - kaptpp*z_trpp

                      if (yloc - yjet < -jwdth) then
                         fstht = -1.e0
                        else if (      yloc - yjet >= -jwdth &
                                &.and. yloc - yjet < jwdth) then
                         fstht = sin(0.5*pi*(yloc - yjet)/jwdth)
                        else if (yloc - yjet >= jwdth) then
                         fstht = 1.e0
                      end if
                     else
                      yjet = yjet0 + kaptpp*z_trpp

                      if (yloc - yjet < -jwdth) then
                         fstht = 1.e0
                        else if (      yloc - yjet >= -jwdth &
                                &.and. yloc - yjet < jwdth) then
                         fstht = -sin(0.5*pi*(yloc - yjet)/jwdth)
                        else if (yloc - yjet >= jwdth) then
                         fstht = -1.e0
                      end if
                   end if

                   if (zloc < z_baro) then
                      zeta_baro &
                      = sin(0.5*pi * (z_baro - zloc)/(z_baro - z_trpp))
                     else
                      zeta_baro = 0.0
                   end if

                   the_env_pp(:,j,k) &
                   =   thet0 * (1.e0 + nstr**2/g_ndim *zloc) &
                     + zeta_baro &
                       * (- thet0 * (nstr**2 - ntrp**2)/g_ndim * z_trpp &
                          + 0.5*deltht * fstht &
                            * (1.0 + 5.0*(z_trpp - zloc)/z_trpp))
                end if

                dens_env_pp(:,j,k) = Pstrat(k) /  the_env_pp(:,j,k)

                if (fluctuationMode) then
                   var(1:nx,j,k,1) = dens_env_pp(1:nx,j,k) - rhoStrat(k)
                  else
                   var(1:nx,j,k,1) = dens_env_pp(1:nx,j,k)
                end if
             end do
          end do

          ! Exner pressure just below the bottom
          ! (so that the bottom pressure vanishes)

          do j = 1, ny
             do i = 1, nx
                if (fluctuationMode) then
                   rhop0 = var(i,j,0,1)
                   rhop1 = var(i,j,1,1)

                   rho0 = var(i,j,0,1) + rhoStrat(0)
                   rho1 = var(i,j,1,1) + rhoStrat(1)
                  else
                   rhop0 = var(i,j,0,1) - rhoStrat(0)
                   rhop1 = var(i,j,1,1) - rhoStrat(1)

                   rho0 = var(i,j,0,1)
                   rho1 = var(i,j,1,1)
                end if

                buoy0 = - g_ndim * rhop0/rho0
                buoy1 = - g_ndim * rhop1/rho1

                rho = 0.5*(rho0 + rho1)

                var(i,j,0,5) &
                = - 0.25*dz * Ma2*kappa * rho/PstratTilde(0) &
                  * 0.5*(buoy0 + buoy1)
             end do
          end do

          ! Exner pressure up to just above the lid

          do k = 0, nz
             do j = 1,ny
                do i = 1,nx
                   if (fluctuationMode) then
                      rhop0 = var(i,j,k,1)
                      rhop1 = var(i,j,k+1,1)
   
                      rho0 = var(i,j,k,1) + rhoStrat(k)
                      rho1 = var(i,j,k+1,1) + rhoStrat(k+1)
                     else
                      rhop0 = var(i,j,k,1) - rhoStrat(k)
                      rhop1 = var(i,j,k+1,1) - rhoStrat(k+1)
   
                      rho0 = var(i,j,k,1)
                      rho1 = var(i,j,k+1,1)
                   end if

                   buoy0 = - g_ndim * rhop0/rho0
                   buoy1 = - g_ndim * rhop1/rho1

                   rho = 0.5*(rho0 + rho1)

                   var(i,j,k+1,5) &
                   = var(i,j,k,5) &
                     + 0.5*dz * Ma2*kappa * rho/PstratTilde(k) &
                     * 0.5*(buoy0 + buoy1)
                end do
             end do
          end do
         else if (background == "HeldSuarez") then
          ! define stretched squared sine and cosine field
          ! define stretched quadrupled cosine field

          do j = 1, ny
             yloc = y(j+j00)

             if (yloc < ymin) then
                stop 'ERROR: y < ymin'
               else if ((yloc >= ymin) .and. (yloc < yjets - 0.5*dy_hs)) &
                & then
                s2_strtd(j) = 1.0
                c2_strtd(j) = 0.0
                c4_strtd(j) = 0.0
               else if ((yloc >= yjets - 0.5*dy_hs) &
                &       .and. (yloc < yjets + 0.5*dy_hs)) then
                s2_strtd(j) &
                = sin((yloc - (yjets + 0.5*dy_hs))/dy_hs * 0.5*pi)**2

                c2_strtd(j) &
                = cos((yloc - (yjets + 0.5*dy_hs))/dy_hs * 0.5*pi)**2

                c4_strtd(j) &
                = cos((yloc - (yjets + 0.5*dy_hs))/dy_hs * 0.5*pi)**4
               else if ((yloc >= yjets + 0.5*dy_hs) &
                &       .and. (yloc < yjetn - 0.5*dy_hs)) then
                s2_strtd(j) = 0.0
                c2_strtd(j) = 1.0
                c4_strtd(j) = 1.0
               else if ((yloc >= yjetn - 0.5*dy_hs) &
                &       .and. (yloc < yjetn + 0.5*dy_hs)) then
                s2_strtd(j) &
                = sin((yloc - (yjetn - 0.5*dy_hs))/dy_hs * 0.5*pi)**2

                c2_strtd(j) &
                = cos((yloc - (yjetn - 0.5*dy_hs))/dy_hs * 0.5*pi)**2

                c4_strtd(j) &
                = cos((yloc - (yjetn - 0.5*dy_hs))/dy_hs * 0.5*pi)**4
               else if ((yloc >= yjetn + 0.5*dy_hs) .and. (yloc <= ymax)) &
                & then
                s2_strtd(j) = 1.0
                c2_strtd(j) = 0.0
                c4_strtd(j) = 0.0
               else if (yloc > ymax) then
                stop 'ERROR: y > ymax'
             end if
          end do

          ! y- and z-dependent thermal relaxation rate
          ! z-dependent Rayleigh-damping-rate

          if (ta_hs_dim == 0.0) then
             ka_hs = 0.0
            else
             ka_hs = tRef/ta_hs_dim
          end if

          if (ts_hs_dim == 0.0) then
             ks_hs = 0.0
            else
             ks_hs = tRef/ts_hs_dim
          end if

          if (tf_hs_dim == 0.0) then
             kf_hs = 0.0
            else
             kf_hs = tRef/tf_hs_dim
          end if

          do k=0,nz+1
             sig_pr = pistrat(k)**(1.0/kappa)

             facsig = max(0.0, (sig_pr - sigb_hs)/(1.0 - sigb_hs))

             !testb
             !print*,'k,pistrat(k),sig_pr,facsig,kf_hs'
             !print*,k,pistrat(k),sig_pr,facsig,kf_hs
             !kv_hs(k) = 0.0
             !teste

             kv_hs(k) = kf_hs*facsig 
             
             do j=1,ny
                kt_hs(j,k) = ka_hs + (ks_hs - ka_hs)*facsig*c4_strtd(j) 
             end do
          end do

          do j=1,ny
             ! (total) Exner pressure just above and below the surface so 
             ! that it = 1 at the surface, by an Euler integration of 
             ! hydrostatic equilibrium

             tempev &
             = max( tp_strato, tp_srf_trp - tpdiffhor_tropo*s2_strtd(j))

             var(1:nx,j,0,5) = 1.0 + 0.5*dz * kappa/tempev
             var(1:nx,j,1,5) = 1.0 - 0.5*dz * kappa/tempev

             ! potential temperature just below and above the surface

             do k=0,1
                do i=1,nx
                   tempev &
                   = max( tp_strato, &
                          var(i,j,k,5) &
                          * (tp_srf_trp - tpdiffhor_tropo * s2_strtd(j) &
                             - ptdiffvert_tropo/kappa * log(var(i,j,k,5)) &
                               * c2_strtd(j)))
             
                   the_env_pp(i,j,k) = tempev/var(i,j,k,5)
                end do
             end do

             ! for k > 1:
             ! Exner pressure and potential temperature by upward 
             ! integration of hydrostatic equilibrium, using a trapezoidal 
             ! leapfrog

             do k = 2, nz+1
                do i=1,nx
                   pistar &
                   = var(i,j,k-2,5) - 2.0*dz * kappa/the_env_pp(i,j,k-1)

                   tempev &
                   = max( tp_strato, &
                          pistar &
                          * (tp_srf_trp - tpdiffhor_tropo * s2_strtd(j) &
                             - ptdiffvert_tropo/kappa * log(pistar) &
                               * c2_strtd(j)))

                   thetastar = tempev/pistar

                   var(i,j,k,5) &
                   =   var(i,j,k-1,5) &
                     - 0.5*dz &
                       * (kappa/thetastar + kappa/the_env_pp(i,j,k-1))

                   tempev &
                   = max( tp_strato, &
                          var(i,j,k,5) &
                          * (tp_srf_trp - tpdiffhor_tropo * s2_strtd(j) &
                             - ptdiffvert_tropo/kappa * log(var(i,j,k,5)) &
                               * c2_strtd(j)))

                   the_env_pp(i,j,k) = tempev/var(i,j,k,5)
                end do
             end do

             ! density

             do k = 0,nz+1
                dens_env_pp(:,j,k) = Pstrat(k) /  the_env_pp(:,j,k)

                if (fluctuationMode) then
                   var(1:nx,j,k,1) = dens_env_pp(1:nx,j,k) - rhoStrat(k)
                  else
                   var(1:nx,j,k,1) = dens_env_pp(1:nx,j,k)
                end if
             end do
          end do

          ! subtract reference-atmosphere Exner pressure from the total

          do k = 0, nz+1
             var(1:nx,1:ny,k,5) = var(1:nx,1:ny,k,5) - pistrat(k) !FS
          end do
         else
          stop 'ERROR: wrong background for baroclinic_LC'
       end if

       p_env_pp(1:nx,1:ny,0:nz+1) = var(1:nx,1:ny,0:nz+1,5)

       call setHalos( var, "var" )
       call setBoundary (var, flux, "var")

       ! determine horizontal wind from density and Exner-pressure 
       ! fluctuations

       do k=0,nz+1
          do j=1,ny
             do i = 1,nx
                ! u-wind

                rho_int_0m = 0.5*(var(i  ,j-1,k,1) + var(i  ,j  ,k,1))
                rho_int_00 = 0.5*(var(i  ,j  ,k,1) + var(i  ,j+1,k,1))
                rho_int_pm = 0.5*(var(i+1,j-1,k,1) + var(i+1,j  ,k,1))
                rho_int_p0 = 0.5*(var(i+1,j  ,k,1) + var(i+1,j+1,k,1))

                rho = var(i,j,k,1)

                if (fluctuationMode) then
                   rho_int_0m = rho_int_0m + rhoStrat(k)
                   rho_int_00 = rho_int_00 + rhoStrat(k)
                   rho_int_pm = rho_int_pm + rhoStrat(k)
                   rho_int_p0 = rho_int_p0 + rhoStrat(k)

                   rho = rho + rhoStrat(k)
                end if

                var(i,j,k,2) &
                = - Ro/(Ma2*kappa*dy) &
                    * Pstrat(k) &
                    * 0.25 &
                    * (  (var(i  ,j  ,k,5) - var(i  ,j-1,k,5)) &
                         / rho_int_0m &
                       + (var(i  ,j+1,k,5) - var(i  ,j  ,k,5)) &
                         / rho_int_00 &
                       + (var(i+1,j  ,k,5) - var(i+1,j-1,k,5)) &
                         / rho_int_pm &
                       + (var(i+1,j+1,k,5) - var(i+1,j  ,k,5)) &
                         / rho_int_p0)

                u_env_pp(i,j,k) = var(i,j,k,2)
                v_env_pp(i,j,k) = 0.
             end do
          end do
       end do

  
    ! density fluctuations

       if (timeScheme == "semiimplicit" .or. auxil_equ) then
          if( fluctuationMode ) then
              do k = -nbz,nz+nbz
                 var(:,:,k,6) = var(:,:,k,1)
              end do
             else
              do k = -nbz,nz+nbz
                 var(:,:,k,6) = var(:,:,k,1) - rhoStrat(k)
              end do
           end if
        end if

       ! store environmental state also in var_env

       call setHalos( var, "var" )
       call setBoundary (var, flux, "var")

       var_env = var

       !-----------------------------------------------------------
       ! add local potential-temperature perturbation
       !-----------------------------------------------------------
       
       if (add_ptptb) then
          ptptb_x = ptptb_x_dim/lRef
          ptptb_y = ptptb_y_dim/lRef
          ptptb_z = ptptb_z_dim/lRef

          ptptb_dh = ptptb_dh_dim/lRef
          ptptb_dz = ptptb_dz_dim/lRef

          ptptb_amp = ptptb_amp_dim/thetaRef

          do k = 1,nz
             zloc = z(k)
   
             do j = 1,ny
                yloc = y(j00+j)
   
                do i = 1, nx
                   xloc = x(i00+i)
  

                  if (ptptb_dh <= 0.) then
                      rptb &
                      = sqrt(((zloc - ptptb_z)/ptptb_dz)**2)
                   else 

                      rptb &
                           = sqrt(  ((xloc - ptptb_x)/ptptb_dh)**2 &
                           + ((yloc - ptptb_y)/ptptb_dh)**2 &
                           + ((zloc - ptptb_z)/ptptb_dz)**2)
                  end if
   
                  if (rptb <= 1.0) then
                     thtptb = ptptb_amp * cos(0.5*pi*rptb)**2
                  else
                     thtptb = 0.0
                  end if
                    
                                      
                 if (fluctuationMode) then
                     rho = var(i,j,k,1) + rhoStrat(k)!Pstrat(k)
                  else
                     rho = var(i,j,k,1)
                  end if
   
                  theta = Pstrat(k)/rho + thtptb
      
                  if (fluctuationMode) then
                     var(i,j,k,1) =   Pstrat(k)/theta - rhostrat(k)!PStrat(k)
                    else
                     var(i,j,k,1) =   Pstrat(k)/theta
                  end if
                end do
             end do
          end do

          ! add local PT perturbation on SH !FS
          ptptb_y = (-1.)*ptptb_y
          do k = 1,nz
             zloc = z(k)
   
             do j = 1,ny
                yloc = y(j00+j)
   
                do i = 1, nx
                   xloc = x(i00+i)
   
                   rptb &
                   = sqrt(  ((xloc - ptptb_x)/ptptb_dh)**2 &
                          + ((yloc - ptptb_y)/ptptb_dh)**2 &
                          + ((zloc - ptptb_z)/ptptb_dz)**2)
   
                  if (rptb <= 1.0) then
                     thtptb = ptptb_amp * cos(0.5*pi*rptb)**2
                    else
                      thtptb = 0.0
                  end if

   
                  if (fluctuationMode) then
                  rho = var(i,j,k,1) + rhoStrat(k)!Pstrat(k)
                    else
                     rho = var(i,j,k,1)
                  end if
   
                  theta = Pstrat(k)/rho - thtptb
      
                  if (fluctuationMode) then
                     var(i,j,k,1) =  Pstrat(k)/theta - rhoStrat(k)!Pstrat(k)
                    else
                     var(i,j,k,1) =   Pstrat(k)/theta
                  end if
                end do
             end do
          end do
        end if


       !------------------------------------------
       !        Adding the noise
       !------------------------------------------

       ! noise added to density fluctuations in a relative manner

       noise_var = "none"
       noise_mag = 1.0

       if (add_noise) then
          call Random_Seed()
          call noise_array(noise_mag, noise_var, noise)

          do k=1,nz
             do j=1,ny
                do i=1,nx
                   if (fluctuationMode) then
                      var(i,j,k,1) = var(i,j,k,1) * (1.0 + noise(i,j,k))
                     else
                      var(i,j,k,1) &
                      = rhoStrat(k) &
                        + (var(i,j,k,1) - rhoStrat(k)) &
                          * (1.0 + noise(i,j,k))
                   end if
                end do
             end do
          end do
       end if    


       !------------------------------------------------
       !               Baroclinic life cycle: idealistic
       !------------------------------------------------       
   case( 'baroclinic_ID' )

       ! read test case input data   
       read (unit=10, nml=baroclinic_ID) 
!*****************************************************************
! initialize fields
!*****************************************************************
       var(:,:,:,1) = 0.0
       var(:,:,:,2) = 0.0
       var(:,:,:,3) = 0.0
       var(:,:,:,4) = 0.0
       var(:,:,:,5) = 0.0
       var(:,:,:,6) = 0.0

       the_env_pp(:,:,:) = 0.0

       dens_env_pp(:,:,:) = 0.0

       ! achatzc: would be much nicer if everything was coded in 
       ! non-dimensional units

       ! achatzc: it is unfortunate that the below is so far most of the 
       ! time applied with isothermal atmosphere
       ! would be better to have a constant lapse rate of 6K/km, e.g., 
       ! for the troposphere and isothermal above the tropopause

       Ly_end = ly_dim(1)

       ! to also enable test of 2D dynamics in x-z plane
       if (sizeX > 1 .and. sizeY == 1) then
          Lx_end = lx_dim(1)
          bar_sigma_x = bar_sigma_y
       end if

       L_z = z(nz) * lRef

       k_2 = int((z_tr_dim - 0.5*dz* lRef)/(dz* lRef)) + 1
       z_2 = z(k_2) * lRef

       theta_bar_0 = thetaStrat(k_2)*thetaRef  ! at the tropopause 

       if (master) print *, 'theta_bar_0 = ', theta_bar_0

       ! to also enable test of 2D dynamics in x-z plane
       if (sizeX > 1 .and. sizeY == 1) then
          alpha_t &
          = u_strength*f_Coriolis_dim*theta_bar_0*Lx_end*bar_sigma_x &
            /(dTh_atm*g*z_2*pi)
         else
          alpha_t &
          = u_strength*f_Coriolis_dim*theta_bar_0*Ly_end*bar_sigma_y &
            /(dTh_atm*g*z_2*pi)
       end if

       if (master) then
          print*,"Barocl. test. instability measure alpha*dTh_atm/dz = ", &
                & alpha_t*dTh_atm/(dz*lRef)
          print*,"Barocl. test. magnitude of \theta alpha*dTh_atm = ", &
                & alpha_t*dTh_atm
       end if  
           
       if (init_2Dto3D) then
          call init_data2D3D(lastrecordnum, var)

          ! achatzc: here is an isue, since the semi-implicit time-stepping
          ! scheme stores the density fluctuations in var(... , 6)

          ! store env pot temp: local
          the_env_pp(1:nx,1:ny,1:nz) = var(1:nx,1:ny,1:nz,6) 

          ! store env dens: local
          dens_env_pp(1:nx,1:ny,1:nz) = var(1:nx,1:ny,1:nz,1)         

          ! actually not used. left as backup for sponge.
          u_env_pp(0:nx,1:ny,1:nz) = var(0:nx,1:ny,1:nz,2)            
         else 
          !***************************************************************
          ! tropopause height: Z_trop, 
          ! temperature at the ground: T0_t, 
          ! lapse rates in the troposphere and stratosphere: gamma_tr, 
          !                                                  gamma_st, 
          ! pressure at a reference height z_ref: P_ref
          ! exact pressure at tropopause: P_trop 
          !***************************************************************

          ! set local index
          i00 = is + nbx - 1   ! 0 index, 
                               ! replace i -> i + i0 in x and y fields
          j00 = js + nby - 1
      
          if (master) then
             print*, 'k_2 = ', k_2 , 'Numerical tropopause is at z = ', &
                    & z_2, ', while H_t = ', z_tr_dim
          end if

          cp = Rsp/kappa

          if (init_bal == "geostr_id") then
             ! determine the Exner-pressure fluctuations

             do k = 0,nz+1
                z_dim = z(k) * lRef

                if (z_dim.lt.z_2) then   ! Troposphere
                   F_0(k) &
                   = alpha_t &
                     * (2.*z_dim * sin(pi*z_dim/(2.*z_2)) - z_dim**2/z_2)
                  else
                   F_0(k) &
                   = alpha_t &
                     * (  2.*(z_dim-L_z) &
                          * sin(pi*(z_dim-L_z)/(2.*(z_2-L_z))) &
                        - (z_dim-L_z)**2/(z_2-L_z)) &
                     * z_2/(z_2-L_z)
                end if

                F_a = F_0(k)

                ! below also enables test of 2D dynamics in x-z plane

                if (sizeX > 1 .and. sizeY == 1) then
                   do i = 1,nx
                      x_dim = x(i+i00) * lRef
   
                      if (x_dim < Lx_end * (0.25 - 0.5*bar_sigma_x)) then
                         term_a = 1.0
                        else if(x_dim >= Lx_end*(0.25 - 0.5*bar_sigma_x) &
                         & .and. x_dim < Lx_end*(0.25 + 0.5*bar_sigma_x)) &
                         & then
                         term_a &
                         = 0.5 &
                           * (  1.0 &
                              - sin(pi &
                                    * ((x_dim - 0.25*Lx_end) &
                                       /(Lx_end*bar_sigma_x))))
                        else if(x_dim >= Lx_end*(0.25 + 0.5*bar_sigma_x) &
                         .and. x_dim < Lx_end*(0.75 - 0.5*bar_sigma_x)) &
                         & then
                         term_a = 0.0
                        else if(x_dim >= Lx_end*(0.75 - 0.5*bar_sigma_x) &
                         .and. x_dim < Lx_end*(0.75 + 0.5*bar_sigma_x)) &
                         & then
                         term_a &
                         = 0.5 &
                           * (  1.0 &
                              + sin(pi &
                                    * ((x_dim - 0.75*Lx_end) &
                                       /(Lx_end*bar_sigma_x))))
                        else
                         term_a = 1.0
                      end if

                      term_b = (dTh_atm - 2.*dTh_atm*term_a)/theta_bar_0

                      if (balance_eq == 'QG') then
                         stop 'ERROR: balance_eq == QG not provided'
                        else
                         streamfunc = g*F_a*term_b/(f_Coriolis_dim)
 
                         pi_pr_xz(i,k) &
                         = f_Coriolis_dim*streamfunc &
                           /(cp*thetaRef*thetaStrat(k))

                         if (balance_eq == 'QG') then
                            stop 'ERROR: balance_eq == QG not provided'
                           else
                            var(i,:,k,5) = pi_pr_xz(i,k)
                         end if
                      end if
                   end do
                  else
                   do j = 1,ny
                      y_dim = y(j+j00) * lRef
   
                      if (y_dim < Ly_end * (0.25 - 0.5*bar_sigma_y)) then
                         term_a = 1.0
                        else if(y_dim >= Ly_end*(0.25 - 0.5*bar_sigma_y) &
                         & .and. y_dim < Ly_end*(0.25 + 0.5*bar_sigma_y)) &
                         & then
                         term_a &
                         = 0.5 &
                           * (  1.0 &
                              - sin(pi &
                                    * ((y_dim - 0.25*Ly_end) &
                                       /(Ly_end*bar_sigma_y))))
                        else if(y_dim >= Ly_end*(0.25 + 0.5*bar_sigma_y) &
                         .and. y_dim < Ly_end*(0.75 - 0.5*bar_sigma_y)) &
                         & then
                         term_a = 0.0
                        else if(y_dim >= Ly_end*(0.75 - 0.5*bar_sigma_y) &
                         .and. y_dim < Ly_end*(0.75 + 0.5*bar_sigma_y)) &
                         & then
                         term_a &
                         = 0.5 &
                           * (  1.0 &
                              + sin(pi &
                                    * ((y_dim - 0.75*Ly_end) &
                                       /(Ly_end*bar_sigma_y))))
                        else
                         term_a = 1.0
                      end if

                      term_b = (dTh_atm - 2.*dTh_atm*term_a)/theta_bar_0

                      if (balance_eq == 'QG') then
                         stop 'ERROR: balance_eq == QG not provided'
                        else
                         streamfunc = g*F_a*term_b/(f_Coriolis_dim)
 
                         pi_pr_yz(j,k) &
                         = f_Coriolis_dim*streamfunc &
                           /(cp*thetaRef*thetaStrat(k))

                         if (balance_eq == 'QG') then
                            stop 'ERROR: balance_eq == QG not provided'
                           else
                            var(:,j,k,5) = pi_pr_yz(j,k)
                         end if
                      end if
                   end do
                end if
             end do

             ! determine density from the Exner-pressure fluctuations

             do k=1,nz
                do j=1,ny
                   do i = 1,nx
                      ! density

                      var(i,j,k,1) &
                      = - cp * thetaRef /(g * lRef * dz) &
                          * 0.5 &
                          * (  PstratTilde(k-1) &
                               * (var(i,j,k,5) - var(i,j,k-1,5)) &
                             + PstratTilde(k) &
                               * (var(i,j,k+1,5) - var(i,j,k,5))) 
     
                      if (.not. fluctuationMode) then
                         var(i,j,k,1) = var(i,j,k,1) + rhoStrat(k)
                      end if
                   end do
                end do
             end do

             call setHalos( var, "var" )
             call setBoundary (var, flux, "var")

             ! determine horizontal wind from density and Exner-pressure 
             ! fluctuations

             ! below also enables test of 2D dynamics in x-z plane

             if (sizeX > 1 .and. sizeY == 1) then
                do k=1,nz
                   do j=1,ny
                      do i = 1,nx
                         ! v-wind

                         rho_int_m0 &
                         = 0.5*(var(i-1,j  ,k,1) + var(i  ,j  ,k,1)) &
                           * rhoRef
                         rho_int_00 &
                         = 0.5*(var(i  ,j  ,k,1) + var(i+1,j  ,k,1)) &
                           * rhoRef
                         rho_int_mp &
                         = 0.5*(var(i-1,j+1,k,1) + var(i  ,j+1,k,1)) &
                           * rhoRef
                         rho_int_0p &
                         = 0.5*(var(i  ,j+1,k,1) + var(i+1,j+1,k,1)) &
                           * rhoRef

                         rho =  var(i,j,k,1)

                         if (fluctuationMode) then
                            rho_int_m0 = rho_int_m0 + rhoStrat(k)* rhoRef
                            rho_int_00 = rho_int_00 + rhoStrat(k)* rhoRef
                            rho_int_mp = rho_int_mp + rhoStrat(k)* rhoRef
                            rho_int_0p = rho_int_0p + rhoStrat(k)* rhoRef

                            rho = rho + rhoStrat(k)* rhoRef
                         end if

                         var(i,j,k,3) &
                         =  1.0/f_Coriolis_dim &
                            * cp * thetaRef*rhoRef/(uRef*lRef*dx) &
                            * Pstrat(k) &
                            * 0.25 &
                            * (  (var(i  ,j  ,k,5) - var(i-1,j  ,k,5)) &
                                 / rho_int_m0 &
                               + (var(i+1,j  ,k,5) - var(i  ,j  ,k,5)) &
                                 / rho_int_00 &
                               + (var(i  ,j+1,k,5) - var(i-1,j+1,k,5)) &
                                 / rho_int_mp &
                               + (var(i+1,j+1,k,5) - var(i  ,j+1,k,5)) &
                                 / rho_int_0p)

                         ! store env pot temp: local
                         the_env_pp(i,j,k) = Pstrat(k)/rho * rhoRef 

                         ! store env dens: local
                         dens_env_pp(i,j,k) = rho / rhoRef   

                         ! for sponge.
                         u_env_pp(i,j,k) = 0.
                         v_env_pp(i,j,k) = var(i,j,k,3)

                         if (      timeScheme /= "semiimplicit" &
                          & .and. .not.auxil_equ) &
                         & then
                            ! Here it is pot temp (re-calculated below)
                            var(i,j,k,6) = the_env_pp(i,j,k)
                         end if
                      end do
                   end do
                end do

                !testb
                !var(:,:,:,2) = - var(:,:,:,3)
                !var(:,:,:,3) = 0.

                !u_env_pp = - v_env_pp
                !v_env_pp = 0.
                !teste
               else
                do k=1,nz
                   do j=1,ny
                      do i = 1,nx
                         ! u-wind

                         rho_int_0m &
                         = 0.5*(var(i  ,j-1,k,1) + var(i  ,j  ,k,1)) &
                           * rhoRef
                         rho_int_00 &
                         = 0.5*(var(i  ,j  ,k,1) + var(i  ,j+1,k,1)) &
                           * rhoRef
                         rho_int_pm &
                         = 0.5*(var(i+1,j-1,k,1) + var(i+1,j  ,k,1)) &
                           * rhoRef
                         rho_int_p0 &
                         = 0.5*(var(i+1,j  ,k,1) + var(i+1,j+1,k,1)) &
                           * rhoRef

                         rho = var(i,j,k,1) * rhoRef

                         if (fluctuationMode) then
                            rho_int_0m = rho_int_0m + rhoStrat(k)* rhoRef
                            rho_int_00 = rho_int_00 + rhoStrat(k)* rhoRef
                            rho_int_pm = rho_int_pm + rhoStrat(k)* rhoRef
                            rho_int_p0 = rho_int_p0 + rhoStrat(k)* rhoRef

                            rho = rho + rhoStrat(k)* rhoRef
                         end if

                         var(i,j,k,2) &
                         = - 1.0/f_Coriolis_dim &
                             * cp * thetaRef*rhoRef/(uRef*lRef*dy) &
                             * Pstrat(k) &
                             * 0.25 &
                             * (  (var(i  ,j  ,k,5) - var(i  ,j-1,k,5)) &
                                  / rho_int_0m &
                                + (var(i  ,j+1,k,5) - var(i  ,j  ,k,5)) &
                                  / rho_int_00 &
                                + (var(i+1,j  ,k,5) - var(i+1,j-1,k,5)) &
                                  / rho_int_pm &
                                + (var(i+1,j+1,k,5) - var(i+1,j  ,k,5)) &
                                  / rho_int_p0)

                         ! store env pot temp: local
                         the_env_pp(i,j,k) = Pstrat(k)/rho * rhoRef 

                         ! store env dens: local
                         dens_env_pp(i,j,k) = rho / rhoRef   

                         ! for sponge.
                         u_env_pp(i,j,k) = var(i,j,k,2)
                         v_env_pp(i,j,k) = 0.

                         if (      timeScheme /= "semiimplicit" &
                           & .and. .not.auxil_equ) &
                         & then
                            ! Here it is pot temp (re-calculated below)
                            var(i,j,k,6) = the_env_pp(i,j,k)
                         end if
                      end do
                   end do
                end do

                !testb
                !var(:,:,:,3) = var(:,:,:,2)
                !var(:,:,:,2) = 0.

                !v_env_pp = u_env_pp
                !u_env_pp = 0.
                !teste
             end if

             ! achatzc: this loop probably not needed?
             do j=1,ny
                u_const(0:nx,j) = var(0:nx,j,nz,2)
             end do
            else if (init_bal == "hydrost") then
             do k = 1,nz
                do j = 1,ny
                   do i = 1,nx
                      if( fluctuationMode )  then
                          var(i,j,k,1) = 0.  ! density
                        else
                         var(i,j,k,1) = rhoStrat(k)  ! density
                      end if                        

                      var(i,j,k,2) = 0.0  ! u
                      var(i,j,k,3) = 0.0  ! v
                      var(i,j,k,4) = 0.0  ! w
                      var(i,j,k,5) = 0.0  ! deviation of Exner pressure

                      ! achatzc: this is potentially an issue since 
                      ! var (... , 6) is used in semi-implicit time 
                      ! stepping for the density fluctuations

                      if (      timeScheme /= "semiimplicit" &
                        & .and. .not.auxil_equ) then
                         var(i,j,k,6) = thetaStrat(k)  ! pot temp
                      end if
                   end do
                end do
             end do
            else
             stop "initialize: init_bal not def. for this model."
          end if
       end if

       do k = 1, nz-1
          do j = 1, ny
             do i = 1, nx
                br_vais_sq(i,j,k) &
                = g &
                  * (the_env_pp(i,j,k+1) - the_env_pp(i,j,k)) &
                  / (dz*the_env_pp(i,j,k)* lRef)
             end do
          end do
       end do
    
       !------------------------------------------
       !        Output of background \theta
       !------------------------------------------

       if (output_theta_bgr) then
          call output_background(the_env_pp, nz, 'theta_bgr.dat', &
                               & thetaRef)
       end if

       if (output_rho_bgr) then
          call output_background(dens_env_pp, nz, 'rho_bgr.dat', rhoRef)
       end if

       !------------------------------------------
       !        Output of background \theta
       !------------------------------------------

       if (output_br_vais_sq) then
          call output_background(br_vais_sq, nz-1, 'Br_Va_bgr.dat', 1.)
       end if

       do k = 1, nz-1
          do j = 1, ny
             do i = 1, nx
                balance(i,j,k) &
                =  uRef * (var(i,j,k+1,2) - var(i,j,k,2)) / (dz*lRef) &
                 + g * (var(i,j+1,k,6) - var(i,j,k,6)) &
                   / (f_Coriolis_dim * dy*lRef * thetaStrat(k))
             end do
          end do
       end do

       do k = 1, nz-1
          do j = 1, ny
             do i = 1, nx
                balance1(i,j,k) &
                =   f_Coriolis_dim*uRef*var(i,j,k,2) &
                  + thetaRef*Rsp*var(i,j,k,6) &
                    * (var(i,j+1,k,5) - var(i,j,k,5))/(kappa*dy*lRef)

                balance3(i,j,k) &
                = g*(var(i,j,k,6) - thetaStrat(k))/thetaStrat(k)

                balance4(i,j,k) &
                = thetaRef*Rsp*var(i,j,k,6) &
                  * (var(i,j,k+1,5) - var(i,j,k,5))/(kappa*dz*lRef)

                balance2(i,j,k) &
                =   thetaRef*Rsp*var(i,j,k,6) &
                    * (var(i,j,k+1,5) - var(i,j,k,5))/(kappa*dz*lRef) &
                  - g*(var(i,j,k,6) - thetaStrat(k))/thetaStrat(k)
             end do
          end do
       end do

       call output_background(balance, nz-1, 'balance.dat', 1.)
       call output_background(balance1, nz-1, 'balance1.dat', 1.)
       call output_background(balance2, nz-1, 'balance2.dat', 1.)
       call output_background(balance3, nz-1, 'balance3.dat', 1.)
       call output_background(balance4, nz-1, 'balance4.dat', 1.)

       !------------------------------------------
       !        Adding the noise
       !------------------------------------------

       ! noise added to density fluctuations in a relative mainer

       noise_var = "none"
       noise_mag = 1.0

       if (add_noise) then
          call Random_Seed()
          call noise_array(noise_mag, noise_var, noise)

          do k=1,nz
             do j=1,ny
                do i=1,nx
                   if (fluctuationMode) then
                      var(i,j,k,1) = var(i,j,k,1) * (1.0 + noise(i,j,k))
                     else
                      var(i,j,k,1) &
                      = rhoStrat(k) &
                        + (var(i,j,k,1) - rhoStrat(k)) &
                          * (1.0 + noise(i,j,k))
                   end if
                end do
             end do
          end do
       end if    

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


       ! ----------------------------------------------------------------

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

       ! ----------------------------------------------------------------

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


       ! ----------------------------------------------------------------


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


       ! -----------------------------------------------------------------


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

       ! -----------------------------------------------------------------

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

    case( "nIce_w_test" )  ! for pseudo_incompressible case
       ! This testcase emulates Peter Spichtinger's box model.  
         
       updateMass = .false.
       predictMomentum = .false.
       correctMomentum = .false.
       
       init_SIce=SIce_crit(T_nuc)
       
       ! initial atmospheric background flow
       var(:,:,:,2) = backgroundFlow_dim(1) / uRef
       var(:,:,:,3) = backgroundFlow_dim(2) / uRef
       var(:,:,:,4) = backgroundFlow_dim(3) / uRef


       ! constant pressure variable pi' 
       var(:,:,:,5) = 0.0

       ! background pot temp
       do k=0,nz
         var(:,:,k,6) = thetaStrat(k)
         if( fluctuationMode ) then
            var(:,:,k,1) = 0.0
         else
            var(:,:,k,1) = rhoStrat(k) 
         end if
       end do


       ! initialize the ice variables according to iceTestcase
       if (include_ice) then
         call setup_ice(var)
       else 
         stop "Error: ice variables must be included for nIce_w_test"
       end if
       
       !---------------------------------------------------------------------------
       
    case ("projectionTest")
       ! test the projection step

       var(:,:,:,1) = 1.0       ! constant density

       var(:,:,:,2) = 1.0       ! constant velocity in x

       var(:,:,:,3:4) = 0.0     ! zero velocity in y,z

       ! pressure sine curve 
       do i = 0,nx
          var(i,:,:,5) = 1.0 + 0.01*sin( 2*pi*x(i) )
       end do


       ! -----------------------------------------------------------------

    case ("matrixStructureTest")
       ! test the momentum transport at constant density and 
       ! pressure along x. 

       ! constant density
       var(:,:,:,1) = 2.0

       ! velocity sine curve along z
       do k = 0,nz
          zw = z(k) + dz/2.0
          var(:,:,k,4) = sin( 2*pi*zw)
       end do
       var(:,:,:,2) = 0.0
       var(:,:,:,3) = 0.0


       ! constant Exner pressure
       var(:,:,:,5) = 1.0


       ! ----------------------------------------------------------------


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


       ! ----------------------------------------------------------------


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


       ! ----------------------------------------------------------------


    case ("momentumFluxTest")
       ! constant density
       var(:,:,:,1) = 2.0

       ! constant velocity fields
       var(:,:,:,2) = 1.0
       var(:,:,:,3) = 2.0
       var(:,:,:,4) = 3.0

       ! constant Exner pressure
       var(:,:,:,5) = 1.0


       ! ----------------------------------------------------------------


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


       ! -----------------------------------------------------------------


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


       ! ------------------------------------------------------------


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
    write(*,fmt="(a25,es10.1,a,f5.1,a)") "PStrat = ", &
            & PStrat(nz)*pRef/1000.0,&
            & " kPa at z = ", z(nz) * lRef/1000.0, " km"
    write(*,fmt="(a25,es10.1,a,f5.1,a)") "rhoStrat = ", &
            & rhoStrat(nz)*rhoRef,&
            & " kg/m3 at z = ", z(nz) * lRef/1000.0, " km"
    write(*,fmt="(a25,es10.1,a,f5.1,a)") "thetaStrat = ", &
            & thetaStrat(nz)*thetaRef,&
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
       write(*,fmt="(a25,f7.1,a)") "theta0 = ", &
               & thetaStrat(1)*thetaRef, " K"
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
    !write(*,fmt="(a25,a,f6.1,a,f6.1,a)") "[xMin, xMax] = ", &
    !     & "[",lx(0)*lRef/1000.0,"km, ",lx(1)*lRef/1000.0,"km ]" 
    print*,"[xMin, xMax] = ", &
         & "[",lx(0)*lRef/1000.0,"km, ",lx(1)*lRef/1000.0,"km ]" 
    !write(*,fmt="(a25,a,f6.1,a,f6.1,a)") "[yMin, yMax] = ", &
    !     & "[",ly(0)*lRef/1000.0,"km, ",ly(1)*lRef/1000.0,"km ]" 
    print*,"[yMin, yMax] = ", &
         & "[",ly(0)*lRef/1000.0,"km, ",ly(1)*lRef/1000.0,"km ]" 
    !write(*,fmt="(a25,a,f6.1,a,f6.1,a)") "[zMin, zMax] = ", &
    !     & "[",lz(0)*lRef/1000.0,"km, ",lz(1)*lRef/1000.0,"km ]"
    print*,"[zMin, zMax] = ", &
         & "[",lz(0)*lRef/1000.0,"km, ",lz(1)*lRef/1000.0,"km ]"
    !write(*,fmt="(a25,i4,a,i4,a,i4)") "nx x ny x nz = ", &
    !     & nx," x ", ny, " x ", nz
    print*,"nx x ny x nz = ", nx," x ", ny, " x ", nz
    print*,""


    print*,"  8) Boundary conditions: "
    write(*,fmt="(a25,a)") "xBoundary = ", trim(xBoundary)
    write(*,fmt="(a25,a)") "yBoundary = ", trim(yBoundary)
    write(*,fmt="(a25,a)") "zBoundary = ", trim(zBoundary)
    if( spongeLayer ) then
       write(*,fmt="(a25,a)") "sponge layer = ", "on"
       write(*,fmt="(a25,f5.1,a)") "height = ", &
            & (lz(1)-lz(0))*spongeHeight*lRef/1000.0," km"
       write(*,fmt="(a25,es8.1,a)") "relaxation  = ", &
               & spongeAlphaZ_dim, " 1/s"
    else
       write(*,fmt="(a25,a)") "sponge layer = ", "off"
    end if
    print*,""



    print*,"  9) Poisson Solver: "
    write(*,fmt="(a25,a)") "solver = ", poissonSolverType
    write(*,fmt="(a25,es8.1)") "tolPoisson = ", tolPoisson
    write(*,fmt="(a25,es8.1)") "tolCond = ", tolCond
    print*,""


    print*," 10) Topography: "
    if( topography) then 
       write(*,fmt="(a25,a)") "topography = ", "on"
       write(*,fmt="(a25,f6.1,a)") "mountain height = ", mountainHeight_dim, " m"
       write(*,fmt="(a25,es8.1,a)") "mountain width = ", mountainWidth_dim, " m"
    else
       write(*,fmt="(a25,a)") "topography = ", "off"
    end if
    print*,""

    end if   ! modified by Junhong Wei (20170216)



    !-------------------------------------------------------------------


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


    !--------------------------------------------------------------------


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


    !--------------------------------------------------------------------


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


    !---------------------------------------------------------------------


    subroutine init_GWP(Psi,kk,mm,ll)

      !------------------------------------------------
      !  calculate complex amplitudes for
      !    1) first harmonics,  leading order: Psi(:,:,:,0)
      !    2) second harmonics, leading order: Psi(:,:,:,1)
      !------------------------------------------------

      ! WARNING: 
      ! 2nd harmonics probably not ready for 2D and 3D wave packets

      ! in/out variables
      ! wave amplitude
      complex, dimension(0:nx+1,0:ny+1,0:nz+1,5,0:2), intent(out) :: Psi
      real, intent(out) :: kk,mm
      real, intent(out) :: ll

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

      ! more debugging stuff
      real :: B11_pinc
      real :: D1TH11

      integer :: i0, j0   ! modified by Junhong Wei (20161201)

      complex :: tmp_var_3DWP

      !-----------------------
      !      MPI stuff
      !-----------------------
      i0 = is + nbx - 1   ! 0 index, replace i -> i + i0 in x and y fields
      j0 = js + nby - 1

      !------------------------
      !    Init data
      !-----------------------

      ! scale input data
      lambdaX = lambdaX_dim/lRef     ! non-dim wave length in x dir.
      lambdaY = lambdaY_dim/lRef     ! non-dim wave length in y dir.
      lambdaZ = lambdaZ_dim/lRef     ! non-dim vert. wave length

      xCenter = xCenter_dim/lRef     ! scaled position wave packtet x dir.
      yCenter = yCenter_dim/lRef     ! scaled position wave packtet y dir.
      zCenter = zCenter_dim/lRef     ! scaled position wave packtet z dir.

      sigma_x = sigma_hor_dim/lRef   ! x width of Gaussian distribution
      sigma_y &
      = sigma_hor_yyy_dim/lRef       ! y width of Gaussian distribution
      sigma_z = sigma_dim/lRef       ! vert. width of Gaussian distribution



      L_cos = L_cos_dim/lRef         ! half length of cosine profile

      if(ABS(lambdaY_dim) /= 0.0) then
         ll = 2.0*pi/lambdaY 
        else
         ll = 0.0
      end if

      if(ABS(lambdaX_dim) /= 0.0) then
         kk = 2.0*pi/lambdaX
        else
         kk = 0.0
      end if


      mm = 2.0*pi/lambdaZ   
      kk2 = kk**2
      mm2 = mm**2
      kTot2 = kk2 + mm2 +  ll * ll
      kTot = sqrt(kTot2)

      ! intrinsic frequency
      omi = omiSign * sqrt(N2 * (kk*kk + ll*ll) + RoInv*RoInv*mm*mm)/kTot
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


      if( master ) then
         print*,"omi = ", omi/tRef
         print*,"mm = ", mm/lRef

         print*,"RoInv = ", RoInv/tRef   ! modified by Junhong Wei

        print*,""
        print*,"  0) Test case: "
        write(*,fmt="(a25,a35)") "Test case  = ", &
                & "wave packet (full model)"
        write(*,fmt="(a25,f10.1,a)") "lambda_x = ", lambdaX_dim, " m"
        write(*,fmt="(a25,f10.1,a)") "lambda_z = ", lambdaZ_dim, " m"
        write(*,fmt="(a25,f10.1a7)") "c_x  = ", omi/kk*uRef, " m/s"
        write(*,fmt="(a25,f10.1,a7)") "c_z  = ", omi/mm*uRef, " m/s"
        write(*,fmt="(a25,f10.1,a7)") "cg_x  = ", &
                & -NN*mm**2/kTot**3 * uRef, " m/s"
        write(*,fmt="(a25,f10.1,a7)") "cg_z  = ", &
                & NN*mm*kk/kTot**3 * uRef, " m/s"
        write(*,fmt="(a25,f10.1,a7)") "u_jet  = ", u0_jet_dim, " m/s"
        print*,""
     end if   ! modified by Junhong Wei (20170216)

     !---------------------------------------
     !        calc amplitude Psi_1^0 
     !     (first harmonic, leading order)
     !---------------------------------------

      do k = 0,(nz+1)
         do j = 0,(ny+1)
            do i = 0,(nx+1)
               ! profile: 1D and 2D

               if( wavePacketDim == 1 ) then
                  delx = 0.0
                 else
                  delx = ( x(i+i0) -xCenter)
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
               ! cosine profile horizontally so that fields are zero 
               ! at the horizontal boundaries
               ! in case of zero sigma in x or y direction use infinity


               if(sigma_x == 0.0) then
                  envel = 1.0
                 else if(abs(delx) < sigma_x) then
                  envel &
                  = 1.0 - amp_mod_x + amp_mod_x *cos(delx*pi/(sigma_x*2.0))
                 else
                  envel = 1.0 - amp_mod_x
               end if

               if(sigma_y == 0.0) then
                  envel = 1.0 * envel
                 else if(abs(dely) < sigma_y) then
                  envel &
                  = (1.0 - amp_mod_y &
                     + amp_mod_y * cos(dely*pi/(sigma_y*2.0))) &
                    * envel
                 else
                  envel = envel * (1.0 - amp_mod_y)
               end if

               envel = envel * exp(-(delz**2)/2./sigma_z**2)

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

               theta0 = thetaStrat(k)

               tmp_var_3DWP &
               = cmplx( 0.0,  &
                        (omi*omi - N2) / (mm*N2*(omi*omi - RoInv*RoInv)))

               u10 = tmp_var_3DWP * cmplx( kk*omi, ll*RoInv ) * b11

               w10 = cmplx(0.0, omi/N2) * b11

               pi12 &
               = cmplx(0.0, &
                       kappa*Ma2*(omi*omi - N2) / N2 / mm / theta0) * b11

               Psi(i,j,k,:,1) &
               = (/u10, w10, b11, pi12, &
                   ( tmp_var_3DWP * cmplx( ll*omi, -kk*RoInv) * b11 ) /)

            end do
         end do
      end do

      !---------------------------------------
      !        calc amplitude Psi_2^1 
      !     (second harmonic, first order)
      !---------------------------------------

      ! WARNING:
      ! this part would probably still tp be adjusted fpor 2D or 3D wave p.

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
            d1theta11 &
            = 0.5*(u10_c*dtheta11_dx + w10_c*dtheta11_dz + Div*theta11_c)

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

            Psi(i,j,k,:,2) &
            = (/u21,w21,b22,pi23, ( cmplx( 0.0, 0.0 ) * b11 ) /)
         end do
      end do

    end subroutine init_GWP

  end subroutine initialise


  !-------------------------------------------------------------------    


  function cphase(c) result(phi)   ! currently not being used

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


!------------------------------------------------------------------
 subroutine init_data2D3D( &
       & iIn, &
       & var )
    !-------------------------------
    !  reads data from file pf_all_in.dat
    !-------------------------------
    
    ! input counter
    integer, intent(in) :: iIn
    
    ! argument fields
    real,dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),intent(out)  &
    & :: var
    
    ! local and global output field and record nbs.
    real*4,dimension(nx,ny) :: field_prc
    integer irc_prc,irc_out
    real*4,dimension(1,sizeY) :: field_in
    
    ! local variables
    integer :: i,j,k, iVar

    ! needed for output to screen
    character (len = 20)  :: fmt, form


    integer :: i_prc,i_mst,i_out,j_prc,j_mst,j_out

    if( master ) then
       print*,""
       print*," Input from File ", fileinitstate2D
       print*,""
       write(*,fmt="(a25,i15)") " reading record no. ", iIn
    end if

    
    ! open input file
   if(master) then
       open(40,file=fileinitstate2D,form="unformatted",access='direct',&
            & recl=1*SizeY)
!      print*,"pf_all_in.dat opened"
    end if

    var = 0.0

    !---------------------------------------
    !       layerwise input and non-dimensionalizing
    !---------------------------------------

    irc_prc = 0
    do iVar = 1, nVar
       if (varIn(iVar)==1) irc_prc = irc_prc + 1
    end do
    irc_prc = irc_prc * iIn * nz

    do iVar = 1, nVar
       if (varIn(iVar)==1) then
          do k = 1, nz
             ! read data layerwise

             irc_prc = irc_prc + 1

             if(master) then
                read(40,rec=irc_prc) field_in


                do j=1,ny
                   j_mst=j

                   do j_prc= 1,nprocy
                      j_out=ny*(j_prc-1)+j

                      do i_prc=1,nprocx
                         do i=1,nx
                            i_out = 1 !nx*(i_prc-1)+i
   
                            i_mst=nprocy*nx*(i_prc-1)+(j_prc-1)*nx+i


                            field_mst(i_mst,j_mst)=field_in(i_out,j_out)
                         end do
                      end do
                   end do
                end do
             end if


             call mpi_barrier(comm,ierror)

             do j = 1, ny
                ! data distributed over all processors


                call mpi_scatter(field_mst(1,j),nx,mpi_real,&
                                 field_prc(1,j),nx,mpi_real,0,comm,&
                                 ierror)
                do i = 1, nx
                   ! non-dimensionalization

                   select case (iVar) 

                      case(1) ! density
                        if(fluctuationMode) then
                           if(rhoOffset) then
                              var(i,j,k,iVar) = field_prc(i,j) / rhoRef
                             else
                              var(i,j,k,iVar) = field_prc(i,j) / rhoRef &
                                                - rhoStrat(k)
                           end if
                          else
                           if(rhoOffset) then
                              var(i,j,k,iVar) = field_prc(i,j) / rhoRef &
                                                + rhoStrat(k)
                             else
                              var(i,j,k,iVar) = field_prc(i,j) / rhoRef
                           end if
                        end if
                      
                      ! interpolate velocities to cell faces
                      
                      case(2) ! u velocity
                        var(i,j,k,iVar) = field_prc(i,j) / uRef &
                                          + offset(iVar)


                      case(3) ! v velocity
                        var(i,j,k,iVar) = field_prc(i,j) / uRef &
                                          + offset(iVar)


                      case(4) ! w velocity
                        var(i,j,k,iVar) = field_prc(i,j) / uRef &
                                          + offset(iVar)

                      case(5) ! Exner function pi' 
                              !(deviation from background)
                        var(i,j,k,iVar) = field_prc(i,j)

                      case(6) ! potential temperature theta' 
                              ! (deviation from background, Boussinesq)
                              ! potential temperature not used 
                              ! (only density  needed)
                        var(i,j,k,iVar) = field_prc(i,j)/thetaRef


                      case(7) ! dynamic Smagorinsky coefficient
                              !(deviation from background)

                        var(i,j,k,iVar) = field_prc(i,j) / (uRef*lRef)

                      case default
                        stop "tec360: unkown iVar"
                   end select ! iVar
                end do ! i
             end do ! j
          end do ! k
       end if
    end do ! iVar

    !------------------------------------
    !              close file
    !------------------------------------
    if(master) close(unit=40)
    
  end subroutine init_data2D3D
!------------------------------------------------------------------
 subroutine therm_rel_param
    ! local variables
    integer :: k
    integer :: allocstat
    real :: tau_sc, spongeAlphaZ_inv
    real, dimension(1:nz)  :: tau_z
    !-------------------------------
    !  creates height-dependent thermal relaxation function
    !-------------------------------
    ! allocate relaxation parameter: hight-dependent
    !allocate(tau_z(1:nz),stat=allocstat)
    !if(allocstat /= 0) stop "init.f90: could not allocate tau_z"
    


      ! nondimensionalize heating relaxation parameter  
      tau_sc = tau_relax/tref
      ! nondimensionalize sponge relaxation parameter   
      spongeAlphaZ_inv = 1./(spongeAlphaZ_dim * tRef)    

      do k = 1,nz
        if (spongeLayer) then
            if (k.ge.kSponge) then
!            print *, 'k=', k, ' z(k)=', z(k), ' z_dim(k)=', lref*z(k)
                  select case(Sponge_Rel_Type)
                    case("constant")
                        tau_z(k) = tau_sc
                    case("linear")
                        tau_z(k) = ((z(k) - z(kSponge))*spongeAlphaZ_inv + (z(nz) - z(k))*tau_sc)/(z(nz) - z(kSponge))                        
                    case default
                        stop "init: relaxation is not defined."
                end select
            else
                tau_z(k) = tau_sc
            end if
        else 
            tau_z(k) = tau_sc
        end if
        if (tau_z(k).le.1.e-10) then
            stop "init: small thermal relaxation parameter."
        end if    
      end do

 
    
  end subroutine therm_rel_param
  !-------------------------------------------------------------------------
    subroutine noise_array(amplitude, ntovar, noise)

    ! local variables
    integer :: i,j,k
    integer :: allocstat, root
    real :: valRef, amplitude, perturbVal, sum_loc, sum_glob
    character(len=20) :: ntovar
    real :: rand   ! random number
    real :: Lx, Ly, Lz
    real, dimension(1:nx,1:ny,1:nz)  :: noise

    !-------------------------------
    !  add noise
    !-------------------------------

    select case( ntovar ) 
        case( "rho" )
            valRef = rhoref
        case( "vel" )
            valRef = uref    
        case( "th" )
            valRef = thetaRef
        case( "none" )
            valRef = 1.0
        case default
            stop "noise_array: Unknown variable" 
    end select        

    do k=1,nz
       do j = 1,ny
          do i = 1,nx                    
             ! calculate a random number between 0 and 1
             call Random_Number(rand) 

             ! calc. perturbation scaled by amplitude
             perturbVal = 2.0*(0.5-rand)*proc_noise*amplitude

             ! non-dimensionalization ...

             if ( referenceQuantities == "SI" ) then 
                noise(i,j,k) = perturbVal
               else
                noise(i,j,k) = perturbVal/valRef
             end if
          enddo
       enddo
    enddo   

  end subroutine noise_array
  
   !-------------------------------------------------------------------------
   subroutine clean_noise(noise)
    ! local variables
    integer :: i,j,k
    integer :: allocstat, root
    real :: valRef, amplitude, perturbVal, sum_loc, sum_glob
    character(len=20) :: ntovar
    real :: rand   ! random number
    real :: Lx, Ly, Lz
    real, dimension(1:nx,1:ny,1:nz)  :: noise


sum_loc = 0.
sum_glob = 0.
    
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx                    
            sum_loc = sum_loc + noise(i,j,k)
        enddo
      enddo
    enddo   

            
    call mpi_allreduce(sum_loc, sum_glob, 1, &
         & mpi_double_precision, mpi_sum, comm, ierror)

    !if(master) then
    !    print*,"sum_glob = ", sum_glob
    !end if         
    
    sum_loc = 0.

    do k=1,nz
      do j = 1,ny
        do i = 1,nx                    
                noise(i,j,k) = noise(i,j,k) - sum_glob/(sizeX*sizeY*sizeZ) ! add noise to 3D-pot. temperature field                                    
                sum_loc = sum_loc + noise(i,j,k)
        enddo
      enddo
    enddo   

    sum_glob = 0.            
    call mpi_allreduce(sum_loc, sum_glob, 1, &
         & mpi_double_precision, mpi_sum, comm, ierror)

    !if(master) then
    !    print*,"sum_glob after norm = ", sum_glob
    !end if        
       
       
  end subroutine clean_noise
  

     !-------------------------------------------------------------------------
   subroutine find_magnitude_y(slice, calc_magn, ntovar)
    ! local variables
    integer :: xslice, zslice
    real :: calc_magn, max_loc, min_loc, max_glob, min_glob, valRef
    !real, dimension(1:nx, 1:ny, 1:nz)  :: field
    real, dimension(1:ny)  :: slice
    character(len=20) :: ntovar

    !slice(1:ny) = field(xslice, 1:ny, zslice)

    max_loc = maxval(slice)
    min_loc = minval(slice)

    call mpi_allreduce(max_loc, max_glob, 1, &
         & mpi_double_precision, mpi_max, comm, ierror)
    call mpi_allreduce(min_loc, min_glob, 1, &
         & mpi_double_precision, mpi_min, comm, ierror)


    calc_magn = max_glob - min_glob

    select case( ntovar ) 
        case( "rho" )
            valRef = rhoref
        case( "vel" )
            valRef = uref    
        case( "th" )
            valRef = thetaRef
        case default
            stop "noise_array: Unknown variable" 
    end select

    if (master) then
     print *, "Noise is added to variable ", ntovar
     if ( referenceQuantities == "SI" ) then 
          print *, "Noise magnitude = ", calc_magn
     else
          print *,  "Noise magnitude = ", calc_magn*valRef
     end if
    end if
  end subroutine find_magnitude_y
   !-------------------------------------------------------------------------
   subroutine averagevel(var)
    ! local variables
    integer :: i,j,k
    real,dimension(:,:,:,:), allocatable,intent(out) :: var

       !-----------------------------
       !  Interpolate to cell faces
       !-----------------------------
print *, "nx = ", nx
print *, "var0=",  var(0,1,1,2)
       ! average zonal velocities to cell face...
       do i = 0,nx
          var(i,:,:,2) = 0.5*( var(i,:,:,2) + var(i+1,:,:,2) )
       end do

       ! average meridional velocities to cell face...
       do j = 0,ny
          var(:,j,:,3) = 0.5*( var(:,j,:,3) + var(:,j+1,:,3) )
       end do

       ! average vertical velocities to cell face...
       do k = 0,nz
          var(:,:,k,4) = 0.5*( var(:,:,k,4) + var(:,:,k+1,4) )
       end do
 
    
    u_env_pp(:,:,:) = var(1:nx,1:ny,1:nz,2)     
  end subroutine averagevel
 !------------------------------------------------------------------------- 
 
  subroutine output_background(th_bgr, max_nz, filename_bgr, ref)
    !-------------------------------
    !  writes background theta
    !-------------------------------
    
    integer max_nz
    ! argument fields
    real,dimension(1:nx,1:ny,1:max_nz),intent(in)  :: th_bgr
    character(len=*) :: filename_bgr

    ! local and global output field and record nbs.
    real*4,dimension(nx,ny) :: field_prc
!   real*4,dimension(SizeX,SizeY) :: field_out
    integer irc_prc,irc_out
    
    ! local variables
    integer :: i,j,k

    integer :: i0,i1, j0,j1, k0,k1

    ! hotBubble output
    real :: rho, theta_dim, ref


    integer :: i_prc,i_mst,i_out,j_prc,j_mst,j_out


    if( master ) then 
        if ((TestCase == "baroclinic_LC").or.(TestCase == "baroclinic_ID")) then
            print*,""
            print*," Output of background into file ", filename_bgr
            print*,""
        else
            print*,""
            print*," There is no environmental state for this testcase, EXIT "
            print*,""        
            stop
        end if 
    end if    
    
    !------------------------------
    !   prepare output file
    !------------------------------
 
    ! open output file

!   open(40,file=dataFile,form="unformatted",access='direct',recl=nx*ny)
    if(master) then
       open(51,file=filename_bgr,form="unformatted",access='direct',&
recl=SizeX*SizeY)
    end if

    
    !---------------------------------------
    !       dimensionalising and layerwise output
    !---------------------------------------


    irc_prc =  0
          do k = 1, max_nz
             ! dimensionalization
             do j = 1, ny
                do i = 1, nx
                        select case(model) 

                           case("pseudo_incompressible")
                       
                             if (referenceQuantities == "SI" ) then
                                theta_dim = th_bgr(i,j,k)
                               else
                                theta_dim = th_bgr(i,j,k) * ref
                             end if
                        
                             field_prc(i,j) = real(theta_dim, kind=4)
                         
                           case( "Boussinesq" )
                              stop "output_background: background undefined"

                           case( "WKB" )
                              stop "output_background: background undefined"
                           case default
                              stop "output_background: unknown model"
                        end select ! model
                      
                end do ! i
                call mpi_gather(field_prc(1,j),nx,mpi_real,&
                                field_mst(1,j),nx,mpi_real,0,comm,ierror)
             end do ! j

             ! layerwise output
             irc_prc=irc_prc+1
             call mpi_barrier(comm,ierror)
             if(master) then
                do j=1,ny
                   j_mst=j

                   do j_prc= 1,nprocy
                      j_out=ny*(j_prc-1)+j

                      do i_prc=1,nprocx
                         do i=1,nx
                            i_out=nx*(i_prc-1)+i

                            i_mst=nprocy*nx*(i_prc-1)+(j_prc-1)*nx+i

                            field_out(i_out,j_out)=field_mst(i_mst,j_mst)
                         end do
                      end do
                   end do
                end do

                write(51,rec=irc_prc) field_out
             end if
          end do ! k



    !------------------------------------
    !              close file
    !------------------------------------
    if(master) close(unit=51)
    

  end subroutine output_background
  !--------------------------------------------------------------------------
end module init_module
