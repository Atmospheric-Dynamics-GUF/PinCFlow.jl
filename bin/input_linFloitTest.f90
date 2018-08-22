
! ------------ grid ------------------------------

&grid
	nx = 80,             ! nb of grid internal cells 
	ny = 1,               ! max 40^3 or 170^2 on acerix with csr storage
	nz = 40,             ! at least nz = 3!
	nbx = 3,              ! nb. of ghost cells
	nby = 3,
	nbz = 3,
	lx_dim = -10000.0, 10000.0,       ! domain lenths in m
	ly_dim = 0.0, 10000.0, 
	lz_dim = 0.0, 10000.0,
&end



! ---------------- nb. of variables ---------------

&variables
	nVar = 5,         ! number of dependent variables
	nOptVar = 4,
&end



! -------------------- solver -------------------
&solverList
	cfl = 1.0
	dtMax = 1.0                    ! max time step
	timeScheme = "LS_Will_RK3"     ! LS_Will_RK3 / Euler / LS_TVD_RK3
	fluxType = "ILES"           ! ILES / central
	reconstType = "SALD"       ! ALDM / constant / SALD
&end



! ------------------ Poisson solver ---------------
&poissonSolverList

	tolPoisson = 1.0e-3          ! abbort: norm(res)/norm(b) < tolPoisson
	maxIterPoisson = 1000
	storageType = "csr"           ! full / csr (compressed sparse row)
	preconditioner = "diag"         ! no / ilu / diag
	initialCleaning = .false.     ! makes initial projection
&end	



! ------------------- atmosphere ------------------
&atmosphereList

	Re = 1.0e10                ! Reynolds number

	background = "isentropic"	
	theta0_dim = 300 !K                ! background pot. temp. in K for 'isentropic'

	N_BruntVaisala_dim = 0.01 !1/s     ! Brunt-Vaisala frequency for 'const-N' atmosphere

	
	! const-N -> set N_BruntVaisala_dim
	! uniform   -> no stratification
	! isothermal -> T = const, rho = p
	! isentropic -> set theta0_dim in 1K 
&end



! -------------------- output ----------------------

&outputList
	outputType = "timeStep"    ! timeStep / time

	nOutput = 1           ! output every nOutput's time step
	maxIter = 1           ! stop after maxIter time steps
	
	outputTime = 60 !s       ! output every ... seconds
	maxTime = 420 !s         ! stop after maxTime seconds
	
	dataFileName = ""    ! empty string "" -> dataFileName = testCase
        restartFile = "restart.ref"   ! restart file in TEC360 format
        restart = .true.      ! true / false

	dimOut = 1,0,1        ! 2D(x,z)-plot dimOut = 1,0,1, 3D with 1,1,1
	varOut = 1,1,1,1,1    ! 1 = output, 0 = no output for (rho,u,v,w,pi)
                              ! rho, u,v,w, pExner
        offset = 0.0, 0.0, 0.0, 0.0, 0.0  ! offset for primary variables
	optVarOut = 1,1,1     ! 1 = output, 0 = no output for 
	                      ! div(Pu), rho, theta
	solutionTime = .true.  ! TECPLOT's "solution time" out yes/no
	showGhostCellsX = .false. ! shows ghost cells in x
	showGhostCellsY = .false. ! shows ghost cells in y 
	showGhostCellsZ = .false. ! shows ghost cells in z
&end



! -------------------- debugging --------------------

&debuggingList
	verbose = .f.
&end



! ----------- testcase: initial and boundary conditions ---

&testCaseList
	testCase = "hotBubble",    ! see init.f90 for test cases
	dTheta0_dim = 6.5 !K
	xRadius_dim = 2000.0 !m
	zRadius_dim = 2000.0 !m
	xCenter_dim = 0.0 !m
	zCenter_dim = 2750.0 !m

	! ---------- available test cases: -----------

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

