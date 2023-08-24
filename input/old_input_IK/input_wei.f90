!------------------------------------------------
!   This is input file for pincFloit
!   The ending .f90 is chosen to facilitate
!   the f90 mode with Emacs
!------------------------------------------------


!------------------------------------------------
!                  Grid and Domain
!------------------------------------------------

&domain

  sizeX = 512,                   ! nb of global grid cells
  sizeY = 1,
  sizeZ = 1000,
  nbx = 3,             ! nb. of ghost cells
  nby = 3,
  nbz = 3,
  lx_dim = 0.0, 9000000.0,       ! domain lenths in m
  ly_dim = 0.0, 40000.0, 
  lz_dim = 0.0, 100000.0,
  nprocx = {nprocx},
  nprocy = {nprocy},

&end



!------------------------------------------------
!                Number of variables 
!------------------------------------------------

&variables

  nVar = 7,         ! number of dependent variables
  nOptVar = 4,

&end


!-----------------------------------------------------------------  
!                          Model equations
!-----------------------------------------------------------------  
  
&modelList

  model = "pseudo_incompressible"    ! pseudo_incompressible / Boussinesq / WKB
  vert_theta = 90.0 !deg    angle of rotation about y
  vert_alpha = 0.0 ! det    angle of rotation about z'

&end


!------------------------------------------------
!                   Pinc Solver
!------------------------------------------------

&solverList

  cfl = 0.5
  cfl_wave = 0.5                 ! passage rate of phase throuh a cell
  dtMax_dim = 1.0 !s              ! max time step in s
  tStepChoice = "cfl"             ! "fix" -> time step dtMax_dim is taken
                                  ! "cfl" -> stability criteria used
  timeScheme = "LS_Will_RK3"      ! LS_Will_RK3 -> Williamson / Euler /
                                  ! LS_TVD_RK3 / CL_TVD_RK3
  fluxType   = "central"           ! ILES / central / upwind
  reconstType = "constant"           ! MUSCL / constant / SALD / ALDM
  limiterType1 = "MCVariant"      ! minmod / Cada / MCVariant
  fluctuationMode = .true.        ! use rho' as primary variable
  DySmaScheme = .true.            ! Dynamic Smagorinsky Scheme

&end

!-------------------------------------------------
!                  Poisson solver
!-------------------------------------------------

&poissonSolverList

  tolPoisson = 1.0e-2          ! abbort 
  maxIterPoisson = 500
  poissonSolverType = "bicgstab"         ! "bicgstab" / "gcr" / "adi" / "hypre"
  storageType = "opr"          ! "csr" (compressed sparse row) 
                                     !  "opr" (lin operator)

  preconditioner = "no"        ! for operator-Solver: "no" / "adi" / "adi_z" 
                                      ! for csr-storage: "ilu", "diag", "no"
  dtau = 4.0e-4                 ! time parameter for ADI (imperical value)
  maxIterADI = 2                ! nb of iterations for ADI preconditioner

  initialCleaning = .false.     ! makes initial projection
  pressureScaling = .true.      ! .true. / .false. Scaling with PStrat
  useNAG = .false.              ! use NAG routine for TDMA algorithm 
  correctMomentum = .true.      ! turn velocity projection on/off 
  correctDivError = .false.      ! true -> subtract rho*div(u)

&end	


!------------------------------------------------
!                    Atmosphere 
!------------------------------------------------

&atmosphereList

  referenceQuantities = "Klein"  ! Klein / WKB / SI / general
  
  specifyReynolds = .false.    ! false -> give mu_viscous, true-> give ReInv
  ReInv = 0.0      !           inverse Reynolds number, ReInv = 0 -> inviscid flow
  mu_viscous_dim = 0.0 ! m^2/s     kinematic viscosity
  mu_conduct_dim = 0.0 ! m^2/s     heat conductivity
  
  background = "isothermal"
  
  !                     const-N    -> set N_BruntVaisala_dim
  !                     isothermal -> set Temp0_dim in K
  !                     isentropic -> set theta0_dim in K 
  !                     uniform    -> constant density (Boussinesq)
  !                     realistic  -> isentr. troposphere / isoth. stratosphere
  !
  !                        for realistic background...
  z_tr_dim = 12000.0 !m    height of tropopause
  theta_tr_dim = 240 !K    const pot. temp. of troposphere

  theta0_dim = 240 ! K      
  !                     isentropic -> background pot. temp.
  !                     const-N    -> ground pot. temp. 
  !                     uniform    -> background pot temp for Boussinesq
  
  Temp0_dim = 240  ! K      
  !                     isothermal -> background temperature
  
  press0_dim =  68880.0
  !                     ground pressure (at z=0) in Pa
  
  N_BruntVaisala_dim = 0.020
  !                     Brunt-Vaisala frequency for 
  !                     1) "const-N" atmosphere in 1/s
  !                     2) "unifrom" Boussinesq 
  
  backgroundFlow_dim =  0.0, 0.0, 0.0 !m/s
  !                     zonal background flow velocity u
  f_Coriolis_dim = 0.00010 ! 1/s       Coriolis parameter
  
&end


!-----------------------------------------------------------
!          Bottom topography (not correctly implemented)
!-----------------------------------------------------------

&topographyList

  topography = .false.   ! set w explicitly at k=1
  mountainHeight_dim = 630       ! in m
  mountainWidth_dim = 1000       ! m

&end
  
  
!-----------------------------------------------------------
!                          Boundary
!-----------------------------------------------------------

&boundaryList

  ! correction of solid wall boundary
  rhoFluxCorr = .false.   ! replace vertical mass flux by CDS at k=1, k=nz-1
  uFluxCorr = .false.
  vFluxCorr = .false.
  wFluxCorr = .false.
  thetaFluxCorr = .false.   ! replace vertical theta flux by CDS at k=1, k=nz-1
  nbCellCorr = 1
  
  ! sponge layer at upper boundary
  spongeLayer = .false.  ! sponge with relaxation to background
  spongeHeight = 0.5     ! relative height of sponge layer
  spongeAlphaZ_dim = 1.67e-3 ! relaxation rate coeff in 1/s 
  
&end

&boundaryList2

  ! boundary types
  xBoundary = "periodic"   ! periodic
  yBoundary = "periodic"   ! periodic
  zBoundary = "solid_wall" ! periodic / solid_wall
  
&end




!------------------------------------------------
!                   Input / Output 
!------------------------------------------------

&outputList

  outputType = "time"    ! timeStep / time
  
  nOutput = 1            ! output every nOutput's time step 
  !                        for outputType = "timeStep"

  maxIter = 2            ! stop after maxIter time steps

  outputTimeDiff =  86400.0  !s    ! output every ... seconds
  maxTime = 86400.0 !s          ! stop after maxTime seconds

  dataFileName = ""      ! empty string "" -> dataFileName = testCase
  restartFile = "restart.ref"   ! restart file in TEC360 format
  restart = .false.      ! true / false

  dimOut = .true.,.false.,.true.      ! 2D(x,z)-plot dimOut = 1,0,1, 3D with 1,1,1

  varOut = 1,1,1,1,1,0,0   ! 1 = output, 0 = no output 
  !                       primary variables: rho,u,v,w,pi',theta'

  offset = 0.0, 0.0, 0.0, 0.0, 0.0 ! offset for primary variables
  rhoOffset = .false.               ! subtract background

  ! optional variables
  optVarOut = 0,0,0, 0,0,0           ! 1 = output, 0 = no output for 
  !                         1) p-pBar in kPa 
  !                         2) buoyancy
  !                         3) background pressure pBar in kPa
  !                         4) background density rhoBar in kg/m^3
  !                         5) div(Pu)
  !                         6) stratification perturbation db/dz
  thetaOffset = .true.               ! subtract background

  ! WKB variables
  wkbVarOut = 0,0,0,0, 0,0,0,0, 0,0,0,0
  !                         1) Wave action amplitude
  !                         2) u00 in m/s
  !                         3) th02
  !                         4) pi02

  !                         5) u10 in m/s
  !                         6) w10 in m/s
  !                         7) b11 in m/s^2
  !                         8) pi12
 
  !                         9) u21 in m/s
  !                         10) w21 in m/s
  !                         11) b22 in m/s^2
  !                         12) pi23


  solutionTime = .true.     ! TECPLOT's "solution time" out yes/no
  solutionTimeUnit = "min"    ! "s", "min" or "h"
  showGhostCellsX = .false. ! shows ghost cells in x
  showGhostCellsY = .false. ! shows ghost cells in y 
  showGhostCellsZ = .false. ! shows ghost cells in z

&end


!------------------------------------------------
!                 Parameter study 
!------------------------------------------------

&parameterList

  parameterStudy = .false.    ! .true. / .false. 
  startParam = 1
  endParam = 10        stepParam = 1
  paramName = ""

&end


!------------------------------------------------
!           Debugging & Error handling 
!------------------------------------------------

&debuggingList

  verbose = .false.
  dtMin_dim = 1.0e-5       ! stop program if dt < dtMin

&end


!------------------------------------------------
!                 WKB & Ray tracing
!------------------------------------------------

&wkbList

  rayTracer = .false.           ! set up ray tracer
!  nRayRatio = 5                ! nb of rays per finite-volume-cell
  ! waveFluxType = "upwind"      ! upwind / central (both 2nd order)
  ! limiterType = "MCVariant"    ! limter upwind transport operator
                                      ! minmod / thirdOrder / Cada / MCVariant
&end


!------------------------------------------------
!  Testcases: initial and boundary conditions
!------------------------------------------------

! general
&testCaseList

  testCase = "wavePacket"
  ! Boussinesq: uniform_theta, wavePacket
  ! agnesiMountain -> see topography

&end

&monochromeWave

  lambdaZ_dim = 6000.0 !m       vertical wave length
  
&end

!----------------
! wave packets
!----------------

!  test cases: wavePacket, wavePacket_raytracer
&wavePacket
  
  wavePacketType = 1            ! 1 = Gaussian, 2 = Cosine
  wavePacketDim = 2             ! 1 = 1D, 2 = 2D, 3 = 3D ! modified by Junhong Wei for 3DWP (20170828)
  ! If working on the 2.5D Wave Packet, please use wavePacketDim = 2   ! modified by Junhong Wei for 3DWP (20170921)
  lambdaX_dim = 300000.0          ! zonal wave length in m ! modified by Junhong Wei for 3DWP (20170828)
  !                               if lambdaX = 0.0 --> kk0 = 0 ! modified by Junhong Wei for 3DWP (20170828)
  lambdaY_dim = 0.0          ! meridional wave length in m ! modified by Junhong Wei for 3DWP (20170828)
  ! if the absolute value of lambdaY_dim is less than or equal to 0.1 --> ll0 = 0 ! modified by Junhong Wei for 3DWP (20170921)
  lambdaZ_dim = -1000.0         ! vertical wave length in m
  amplitudeFactor = 0.5         ! normalilized buoyancy amplitude
  xCenter_dim = 4500000.0         ! zonal center of wave packet in m
  yCenter_dim = 20000.0          ! meridional center of wave packet in m ! modified by Junhong Wei for 3DWP (20170828)
  zCenter_dim = 30000.0         ! vertical...
  sigma_dim = 5000.0            ! Gaussian distribution width (in the vertical direction)
  sigma_hor_dim = 1500000.0        ! Gaussian distribution width (in the horizontal direction)
  sigma_hor_yyy_dim = 1500000.0    ! Gaussian distribution width (in the horizontal direction) ! modified by Junhong Wei for 3DWP (20170828)
  L_cos_dim = 10000.0           ! half width of cosine profile of GWP
  meanFlowX_dim = 0.0           ! mean flow in m/s / jet flow amplitude
  meanFlowZ_dim = 0.0           ! mean vertical flow m/s
  u0_jet_dim = 0.0              ! amplitude max of jet velocity
  L_jet_dim = 5000.0            ! half width of cosine profile of jet in 1m
  z0_jet_dim = 50000.0          ! center of jet stream
  omiSign = 1                 ! frequency branch

&end

!----------------
!    bubbles
!----------------

! test cases:  ->  Robter_Bubble
&Robert_Bubble

  dTheta1_dim = 0.5 !K          pot temp offset
  a1_dim = 150.0 !m             radius of plateau
  sigma1_dim = 50.0 !m          Gaussian edge profile
  xCenter1_dim = 0.0 !m      
  zCenter1_dim = 300.0 !m

  dTheta2_dim = -0.15 !K        pot temp offset
  a2_dim = 0.0 !m               radius of plateau
  sigma2_dim = 50.0 !m          Gaussian edge profile
  xCenter2_dim = 60.0 !      
  zCenter2_dim = 640.0 !m

&end

! test cases: hotBubble, coldBubble, hotBubble3D
&bubble

  dTheta0_dim = 6.0 !K
  xRadius_dim = 1000.0 !m
  zRadius_dim = 1000.0 !m
  xCenter_dim = 0.0 !m
  zCenter_dim = 0.0 !m
  zExcentricity = 1.0

&end

! ---------- available test cases: -----------
        
        ! ------- test nonlinear WKB theory ---------
        ! wavePacket1D and wavePacket1D_raytracer
        ! wavePacket2D and wavePacket1D_raytracer
        
        
        ! ------- with bottom topography
        ! steadyFlow: test steady background flow
        ! 
        ! agnesiMountain (rf. Smolarkiewiecz) - does not work!!!
        ! gamma = 1.069... N = 0.01 1/s  120km x 120km x 60km

	! ------- with gravity & stratification --------
	! coldBubble: same input variables as hotBubble
	! hotBubble: dTheta0_dim, xRadius_dim, zRadius_dim, 
        ! ...xCenter_dim, zCenter_dim
	! hotBubble3D: like hotBubble with y data taken from x dimension

	! --------------- without gravity --------------
 	! densitySineTransport
	! densityDiscXYZ
	! densityDiscXY
	! densityDiscXZ
	! greshoVortexXY
	! greshoVortexXZ
	! projectionTest       
	! matrixStructureTest, set maxIter = 1, nx = ny = nz = 3

	! ------------------------------------------------------------	

	! momentumTransportX, set maxIter = 100, nx = 50, ny = nz = 1
	! momentumTransportY, set maxIter = 100, ny = 50, nx = nz = 1, 
	! momentumFluxTest, nx = ny = nz = 2
        ! cfl <= 0.60 for RK3	
	
	! ---------------------------------------------

	! densityTransportX, nx = 50, ny = nz = 1 
	! densityTransportY, ny = 50, nx = nz = 1
	! massFluxTest,     set nx = ny = nz = 2

   	! ---------------------------------------------	  

	! xparabel 

