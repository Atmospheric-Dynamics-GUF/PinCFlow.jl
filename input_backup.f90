!------------------------------------------------
!   This is input file for pincFloit
!   The ending .f90 is chosen to facilitate
!   the f90 mode with Emacs
!------------------------------------------------


!------------------------------------------------
!                  Grid and Domain
!------------------------------------------------

&domain

  sizeX = 100,          ! nb of internal cells in global domain
  sizeY = 100,          
  sizeZ = 100,
  nbx = 3,              ! nb. of ghost cells
  nby = 3,
  nbz = 3,

&end

&grid

  nx = 100,             ! nb of grid internal cells 
  ny = 1,               ! max 40^3 or 170^2 on acerix with csr storage
  nz = 100,              ! at least nz = 3!
  lx_dim = -500.0, 500.0,       ! domain lenths in m
  ly_dim = 0.0, 1000.0, 
  lz_dim = 0.0, 1000.0,
 
&end



!------------------------------------------------
!                Number of variables 
!------------------------------------------------

&variables

  nVar = 6,         ! number of dependent variables
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

  cfl = 0.9
  cfl_wave = 1.0                  ! 1 -> GW phases passes a cell / time step
  dtMax_dim = 1.0 !s              ! max time step in s
  tStepChoice = "cfl"             ! "fix" -> time step dtMax_dim is taken
                                  ! "cfl" -> stability criteria used
  timeScheme = "LS_Will_RK3"      ! LS_Will_RK3 -> Williamson / Euler /
                                  ! LS_TVD_RK3 / CL_TVD_RK3
  fluxType   = "upwind"           ! ILES / central / upwind
  reconstType = "MUSCL"           ! MUSCL / constant / SALD / ALDM
  limiterType1 = "MCVariant"      ! minmod / Cada / MCVariant
  fluctuationMode = .true.        ! use rho' as primary variable
                                  ! not rho=rhoStrat + rho'

&end

!-------------------------------------------------
!                  Poisson solver
!-------------------------------------------------

&poissonSolverList

  tolPoisson = 1.0e-9          ! abbort 
  maxIterPoisson = 5000
  poissonSolverType = "bicgstab"         ! "bicgstab" / "gcr" / "adi"
  storageType = "opr"          ! "csr" (compressed sparse row) 
                                     !  "opr" (lin operator)

  preconditioner = "no"        ! for operator-Solver: "no" / "adi" / "adi_z" 
                                      ! for csr-storage: "ilu", "diag", "no"
  dtau = 4.0e-4                 ! time parameter for ADI (imperical value)
  maxIterADI = 2                ! nb of iterations for ADI preconditioner

  initialCleaning = .false.     ! makes initial projection
  pressureScaling = .true.      ! .true. / .false. Scaling with PStrat
  useNAG = .false.               ! use NAG routine for TDMA algorithm 
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
  
  background = "isentropic"
  
  !                     const-N    -> set N_BruntVaisala_dim
  !                     isothermal -> set Temp0_dim in K
  !                     isentropic -> set theta0_dim in K 
  !                     uniform    -> constant density (Boussinesq)
  
  theta0_dim = 300 ! K      
  !                     isentropic -> background pot. temp.
  !                     const-N    -> ground pot. temp. 
  !                     uniform    -> background pot temp for Boussinesq
  
  Temp0_dim = 300  ! K      
  !                     isothermal -> background temperature
  
  press0_dim =  101325.0    
  !                     ground pressure (at z=0) in Pa
  
  N_BruntVaisala_dim = 0.01 
  !                     Brunt-Vaisala frequency for 
  !                     1) "const-N" atmosphere in 1/s
  !                     2) "unifrom" Boussinesq 
  
  backgroundFlow_dim =  0.0, 0.0, 0.0 !m/s
  !                     zonal background flow velocity u
  f_Coriolis_dim = 0.0e-2 ! 1/s       Coriolis parameter
  
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

  outputTimeDiff =  300.0  !s    ! output every ... seconds
  maxTime = 1500.0 !s          ! stop after maxTime seconds

  dataFileName = ""      ! empty string "" -> dataFileName = testCase
  restartFile = "restart.ref"   ! restart file in TEC360 format
  restart = .false.      ! true / false

  dimOut = 1,0,1        ! 2D(x,z)-plot dimOut = 1,0,1, 3D with 1,1,1

  varOut = 1,1,0, 1,1,1   ! 1 = output, 0 = no output 
  !                       primary variables: rho,u,v,w,pi',theta'

  offset = 0.0, 0.0, 0.0, 0.0, 0.0 ! offset for primary variables
  rhoOffset = .true.               ! subtract background

  ! optional variables
  optVarOut = 0,1,0, 0,0,0           ! 1 = output, 0 = no output for 
  !                         1) p-pBar in kPa 
  !                         2) buoyancy
  !                         3) background pressure pBar in kPa
  !                         4) background density rhoBar in kg/m^3
  !                         5) div(Pu)
  !                         6) stratification perturbation db/dz
  thetaOffset = .true.               ! subtract background

  ! WKB variables
  wkbVarOut = 0,0,0, 0,0,0, 0  
  !                         1) Wave action amplitude
  !                         2) u01 in m/s
  !                         4) ... 
  !                         5) 
  !                         6) 
  !                         7) 


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

  rayTracer = .false.          ! set up ray tracer
  outComplex = .true.          ! output amplitudes complex / absolute value
  nRayRatio = 2                ! nb of rays per finite-volume-cell
  waveFluxType = "upwind"      ! upwind / central (both 2nd order)
  limiterType = "MCVariant"    ! limter upwind transport operator
                                      ! minmod / thirdOrder / Cada / MCVariant
&end


!------------------------------------------------
!  Testcases: initial and boundary conditions
!------------------------------------------------

! general
&testCaseList

  testCase = "Robert_Bubble"
  ! Boussinesq: uniform_theta
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

  wavePacketDim = 1             ! 1 = 1D, 2 = 2D
  lambdaX_dim = 3000.0          ! horizontal wave length in m
  !                               if lambdaX = 0.0 --> kk0 = 0
  lambdaZ_dim = 1000.0          ! vertical wave length in m
  amplitudeFactor = 0.1        ! normalilized buoyancy amplitude
  xCenter_dim = 0.0          ! horizontal center of wave packet in m
  zCenter_dim = 30000.0          ! vertical...
  sigma_dim = 5000.0             ! Gaussian distribution width
  meanFlowX_dim = 0.0          ! m/s
  meanFlowZ_dim = 0.0          ! m/s

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

  dTheta2_dim = -0.15 !K          pot temp offset
  a2_dim = 0.0 !m             radius of plateau
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
  zExcentricity = .25

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

