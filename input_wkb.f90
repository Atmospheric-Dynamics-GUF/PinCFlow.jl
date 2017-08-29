
! ------------ grid ------------------------------

&grid
	nx = 2,             ! nb of grid internal cells 
	ny = 1,               ! max 40^3 or 170^2 on acerix with csr storage
	nz = 2,              ! at least nz = 3!
        nbx = 3,              ! nb. of ghost cells
	nby = 3,
	nbz = 3,
	lx_dim = -10.0, 10.0,       ! domain lenths in m
	ly_dim = 0.0, 20.0, 
	lz_dim = 0.0, 20.0,
&end



! ---------------- nb. of variables ---------------

&variables
	nVar = 5,         ! number of dependent variables
	nOptVar = 4,
&end



! -------------------- solver -------------------
&solverList
	cfl = 0.9
	dtMax_dim = 30.0                   ! max time step in s
	timeScheme = "LS_Will_RK3"        ! LS_Will_RK3 / Euler / LS_TVD_RK3
	fluxType = "central"              ! ILES / central
	reconstType = "constant"          ! ALDM / constant / SALD
&end



! ------------------ Poisson solver ---------------
&poissonSolverList
	tolPoisson = 1.0e-6          ! abbort: norm(res)/norm(b) < tolPoisson
	maxIterPoisson = 1000
        poissonSolverType = "bicgstab"         ! "bicgstab" / "gcr" / "adi"
	storageType = "opr"           ! "csr" (compressed sparse row) / "opr" (lin operator)

	preconditioner = "no"        ! for operator-Solver: "no" / "adi" / "adi_z" 
                                      ! for csr-storage: "ilu", "diag", "no"
        dtau = 4.0e-4                 ! time parameter for ADI (imperical value)
        maxIterADI = 2                ! nb of iterations for preconditioner ADI

	initialCleaning = .false.     ! makes initial projection
        pressureScaling = .true.     ! .true. / .false. Scaling with PStrat
        useNAG = .false.               ! use NAG routine for tridiagonal system in ADI
        correctMomentum = .true.      ! turn velocity projection on/off 
&end	



! ------------------- atmosphere ------------------
&atmosphereList

	ReInv = 0.0                ! inverse Reynolds number, ReInv = 0 -> inviscid flow

	background = "isothermal"	
	! const-N -> set N_BruntVaisala_dim
	! uniform   -> no stratification
	! isothermal -> T = const, rho = p
	! isentropic -> set theta0_dim in 1K 

	theta0_dim = 300 ! K                ! background pot. temp. in K for 'isentropic'
                                            ! ground pot. temp. in K for 'const-N' atmosphere
        Temp0_dim = 273  ! K                ! background temperature in K for 'isothermal' 

        press0_dim =  101325.0  ! Pa        ! ground pressure (at z=0)

	N_BruntVaisala_dim = 0.01 ! 1/s     ! Brunt-Vaisala frequency for 'const-N' atmosphere
        
&end

! ---------------- bottom topography ------------------------
&topographyList

        topography = .false.       ! via boundary treatment 
        mountainHeight = 630        ! in m
        mountainWidth = 1000        ! m

&end


! -------------------- input / output ----------------------

&outputList
	outputType = "timeStep"    ! timeStep / time

	nOutput = 1           ! output every nOutput's time step
	maxIter = 10           ! stop after maxIter time steps
	
	outputTimeDiff =  300 !s       ! output every ... seconds
	maxTime = 9000 !s         ! stop after maxTime seconds
	
	dataFileName = ""    ! empty string "" -> dataFileName = testCase
        restartFile = "restart.ref"   ! restart file in TEC360 format
        restart = .false.      ! true / false
        
	dimOut = 1,0,1        ! 2D(x,z)-plot dimOut = 1,0,1, 3D with 1,1,1
	varOut = 1,1,0,1,1    ! 1 = output, 0 = no output for (rho,u,v,w,pi)
                              ! rho, u,v,w, pExner

        offset = 0.0, 0.0, 0.0, 0.0, 0.0    ! offset for primary variables
	optVarOut = 0,0,0,0,0,0,0           ! 1 = output, 0 = no output for 
	                                    ! 1) rho, 2) theta, 3) div(Pu), 4) PStrat, 5) rhoStrat
                                            ! 6) Force_x, 7) Force_z

	solutionTime = .true.  ! TECPLOT's "solution time" out yes/no
        solutionTimeUnit = "min"  ! "s", "min" or "h"
	showGhostCellsX = .false. ! shows ghost cells in x
	showGhostCellsY = .false. ! shows ghost cells in y 
	showGhostCellsZ = .false. ! shows ghost cells in z
&end

! ---------------- parameter study ------------------
&parameterList
        parameterStudy = .false.    ! .true. / .false. 
        startParam = 1
        endParam = 10        stepParam = 1
        paramName = "gForceFactor"
&end

! -------------------- debugging & error handling --------

&debuggingList
	verbose = .false.
        dtMin_dim = 1.0e-5       ! stop program if dt < dtMin
&end


!------------------- wkb / ray tracer -------------------
&wkbList
         rayTracer = .true.           ! set up ray tracer
         nRayRatio = 2                ! nb of rays per finite-volume-cell
         lambdaX_dim = 1.0            ! horizontal wave length in m
         lambdaZ_dim = 1.0            ! vertical wave length in m
         meanFlowX_dim = 10.0         ! m/s
         meanFlowZ_dim = 0.0          ! m/s
&end


! ----------- testcase: initial and boundary conditions ---

&testCaseList
        testCase = "wavePacket"

        ! hotBubble variables
	dTheta0_dim = 6.5 !K
	xRadius_dim = 2000.0 !m
	zRadius_dim = 2000.0 !m
	xCenter_dim = 0.0 !m
	zCenter_dim = 2750.0 !m

        ! agnesi Mountain variables
        backgroundFlow_dim = 10.0, 0.0, 0.0 ! m/2     ! horizontal background flow velocity (along x)

	! ---------- available test cases: -----------
        
        ! ------- test nonlinear WKB theory ---------
        ! wavePacket
        
        
        ! ------- with bottom topography
        ! steadyFlow: test steady background flow
        ! 
        ! agnesiMountain (rf. Smolarkiewiecz) - does not work!!!
        ! gamma = 1.069... N = 0.01 1/s  120km x 120km x 60km

	! ------- with gravity & stratification --------
	! coldBubble: same input variables as hotBubble
	! hotBubble: dTheta0_dim, xRadius_dim, zRadius_dim, xCenter_dim, zCenter_dim
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
&end
