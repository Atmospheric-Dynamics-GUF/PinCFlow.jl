!------------------------------------------------
!   This is input file for pincFloit
!   The ending .f90 is chosen to facilitate
!   the f90 mode with Emacs
!------------------------------------------------


!------------------------------------------------
!                  Grid and Domain
!------------------------------------------------

&grid
	nx = 20,             ! nb of grid internal cells 
	ny = 1,               ! max 40^3 or 170^2 on acerix with csr storage
	nz = 20,              ! at least nz = 3!
        nbx = 3,              ! nb. of ghost cells
	nby = 3,
	nbz = 3,
	lx_dim = -20000, 20000.0,       ! domain lenths in m
	ly_dim = 0.0, 40000.0, 
	lz_dim = 0.0, 40000.0,
&end



!------------------------------------------------
!                Number of variables 
!------------------------------------------------

&variables
	nVar = 5,         ! number of dependent variables
	nOptVar = 4,
&end


!------------------------------------------------
!                   Pinc Solver
!------------------------------------------------

&solverList
                                          ! Problems with SI units: Exner too large. 
        pFactor = 1.0 !                   ! testFactor for p-Gradient: optimal for 1.126...1.127
        cfl = 0.9
	dtMax_dim = 30                    ! max time step in s
        tStepChoice = "cfl"               ! "fix" -> time step dtMax_dim is taken
                                          ! "cfl" -> stability criteria used
	timeScheme = "LS_Will_RK3"        ! LS_Will_RK3 / Euler / LS_TVD_RK3
	fluxType = "ILES"              ! ILES / central
	reconstType = "SALD"          ! xweno / constant / SALD
&end

!-------------------------------------------------
!                  Poisson solver
!-------------------------------------------------

&poissonSolverList
	tolPoisson = 1.0e-8          ! abbort 
	maxIterPoisson = 1000
        poissonSolverType = "bicgstab"         ! "bicgstab" / "gcr" / "adi"
	storageType = "opr"          ! "csr" (compressed sparse row) 
                                     !  "opr" (lin operator)

	preconditioner = "no"        ! for operator-Solver: "no" / "adi" / "adi_z" 
                                      ! for csr-storage: "ilu", "diag", "no"
        dtau = 4.0e-4                 ! time parameter for ADI (imperical value)
        maxIterADI = 2                ! nb of iterations for ADI preconditioner

	initialCleaning = .false.     ! makes initial projection
        pressureScaling = .true.      ! .true. / .false. Scaling with PStrat
        pressurePrediction = .false.   ! advective prediction of pressure field
        useNAG = .false.               ! use NAG routine for TDMA algorithm 
        correctMomentum = .true.      ! turn velocity projection on/off 
&end	


!------------------------------------------------
!                    Atmosphere 
!------------------------------------------------

&atmosphereList

        referenceQuantities = "Klein"   ! Klein / WKB / SI / general
        
	ReInv = 0.0                     ! inverse Reynolds number, 
                                        ! ReInv = 0 -> inviscid flow

	background = "isothermal"
	
	! stable     -> set N_BruntVaisala_dim
	! isothermal -> set Temp0_dim in K
	! isentropic -> set theta0_dim in K 
	! constant   -> no stratification, g is set to 0.0

	theta0_dim = 300 ! K      ! background pot. temp. in K for 'isentropic'
                                  ! ground pot. temp. in K for 'stable' atmosphere
        Temp0_dim = 300  ! K      ! background temperature in K for 'isothermal' 

        press0_dim =  101325.0    ! ground pressure (at z=0) in Pa

	N_BruntVaisala_dim = 0.01 ! Brunt-Vaisala frequency for 
                                  ! 'stable' atmosphere in 1/s

        backgroundFlow = 0.0      ! zonal background flow velocity u in m/s
        
&end


!-----------------------------------------------------------
!          Bottom topography (not correctly implemented)
!-----------------------------------------------------------

&topographyList

        topography = .false.       ! via boundary treatment 
        topography2 = .false.      ! via conservative force
        gForceMonitor = .false.    ! prints gForce on screen at i = nx/2

        mountainHeight = 630       ! in m

        mountainWidth = 1000       ! m

        topForceFactor = 1.0       ! topographic force proportional factor
&end


!------------------------------------------------
!                   Input / Output 
!------------------------------------------------

&outputList
	outputType = "timeStep"    ! timeStep / time

	nOutput = 1            ! output every nOutput's time step 
                               ! for outputType = "timeStep"

	maxIter = 1           ! stop after maxIter time steps
	
	outputTimeDiff =  360  !s    ! output every ... seconds
	maxTime = 10800 !s          ! stop after maxTime seconds
	
	dataFileName = ""      ! empty string "" -> dataFileName = testCase
        restartFile = "restart.ref"   ! restart file in TEC360 format
        restart = .false.      ! true / false
        
	dimOut = 1,0,1        ! 2D(x,z)-plot dimOut = 1,0,1, 3D with 1,1,1
        
        ! primary variables
	varOut = 0,1,0,0,0    ! 1 = output, 0 = no output for (rho,u,v,w,pi)
                              ! rho, u,v,w, ExnerFunction
        offset = 0.0, 0.0, 0.0, 0.0, 0.0 ! offset for primary variables
        rhoOffset = .true.               ! subtract background

        ! optional variables
	optVarOut = 0,0,0, 0,0,1           ! 1 = output, 0 = no output for 
        !                                  1) p-pBar in kPa 
        !                                  2) pot. temperature theta in K 
        !                                  3) background pressure pBar in kPa
        !                                  4) background density rhoBar in kg/m^3
        !                                  5) div(Pu)
        !                                  6) stratification perturbation db/dz
        thetaOffset = .true.               ! subtract background
        
        ! WKB variables
        wkbVarOut = 1,1,0, 0,0,0, 0  
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
        paramName = "gForceFactor"
&end


!------------------------------------------------
!           Debugging & Error handling 
!------------------------------------------------

&debuggingList
	verbose = .false.
        dtMin = 1.0e-5       ! stop program if dt < dtMin
&end


!------------------------------------------------
!                 WKB & Ray tracing
!------------------------------------------------

&wkbList
         rayTracer = .true.           ! set up ray tracer
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
        testCase = "wavePacket_raytracer"
&end

!----------------
! wave packets
!----------------
&wavePacket
!        wavePacket, wavePacket_raytracer
         wavePacketDim = 1             ! 1 = 1D, 2 = 2D
         lambdaX_dim = 3000.0          ! horizontal wave length in m
                                      ! if lambdaX = 0.0 --> kk0 = 0
         lambdaZ_dim = 1000.0          ! vertical wave length in m
         amplitudeFactor = 0.1        ! normalilized buoyancy amplitude
         xCenter_dim = 0.0          ! horizontal center of wave packet in m
         zCenter_dim = 20000.0          ! vertical...
         sigma_dim = 5000.0             ! Gaussian distribution width
         meanFlowX_dim = 0.0          ! m/s
         meanFlowZ_dim = 0.0          ! m/s
&end

!----------------
!    bubbles
!----------------
&bubble
	dTheta0_dim = 6.5 !K
	xRadius_dim = 2000.0 !m
	zRadius_dim = 2000.0 !m
	xCenter_dim = 0.0 !m
	zCenter_dim = 2750.0 !m
&end

        ! agnesi Mountain variables
        

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

