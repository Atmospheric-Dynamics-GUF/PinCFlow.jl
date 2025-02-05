module type_module

    use mpi
    !-----------------------------------------------------------------
    !    Definition of data types and variables accessible
    !    throughout the code! Use care when working on these
    !    data fields. When adding a new namelist parameter,
    !    set its default value in the subroutine default_values
    !    and add an appropriate line in the subroutine
    !    write_namelists.
    !-----------------------------------------------------------------
  
    implicit none
  
    public
  
    !-----------------------------------------------------------------
    !                     MPI & Domain composition
    !-----------------------------------------------------------------
  
    integer :: sizeX, sizeY, sizeZ
    integer :: nbx, nby, nbz
    integer :: nprocx, nprocy
    real, dimension(0:1) :: lx_dim, ly_dim, lz_dim ! dimensional domain
  
    namelist / domain / sizeX, sizeY, sizeZ, nbx, nby, nbz, lx_dim, ly_dim, &
        &lz_dim, nprocx, nprocy
  
    integer :: sizeXX, sizeYY, sizeZZ
    integer :: nx1, in1, ny1, jn1
    integer :: iStart, jStart
    integer :: is, ie, js, je ! local start and end indices
    logical :: verboseMPI
  
    ! MPI include (parameters needed below)
    !include 'mpif.h'
  
    ! MPI variables
    integer :: ierror
    integer, dimension(2) :: dims, coords
    logical, dimension(2) :: periods
    integer :: back, forw, right, left, rank
    integer :: idim, jdim, icoord, jcoord
    integer :: nbProc, comm
    logical :: master
    integer, dimension(mpi_status_size) :: sts_left, sts_right, sts_back, sts_forw
  
    integer, parameter :: root = 0
  
    !-----------------------------------------------------------------
    !                      (local) Grid & Domain
    !-----------------------------------------------------------------
  
    integer :: nx, ny, nz
    real, dimension(0:1) :: lx, ly, lz ! scaled domain
  
    real, dimension(:), allocatable :: x, y, z
    real :: dx, dy, dz
    integer :: nxx, nyy, nzz ! grid size inclusive ghost cells
    integer :: nxyz ! = nx*ny*nz
  
    !--------------------------------------------------------------
    !              for output global fields to single file
    !--------------------------------------------------------------
  
    real * 4, dimension(:, :), allocatable :: field_out, field_mst
  
    !-----------------------------------------------------------------
    ! for sponge: maximum damping rate in 1/dt
    !-----------------------------------------------------------------
  
    real :: alpspg
  
    !-----------------------------------------------------------------
    ! for
    ! (1) Rayleigh damping in land cells and
    ! (2) density-fluctuation relaxation in semi-implicit time stepping
    !-----------------------------------------------------------------
  
    real :: alprlx
  
    type var_type
      real, dimension(:, :, :), allocatable :: rho ! Density
      real, dimension(:, :, :), allocatable :: u ! Zonal wind
      real, dimension(:, :, :), allocatable :: v ! Meridional wind
      real, dimension(:, :, :), allocatable :: w ! Vertical wind
      real, dimension(:, :, :), allocatable :: pi ! Exner pressure
      real, dimension(:, :, :), allocatable :: rhop ! Density fluctuations
    end type var_type
  
    type flux_type
      real, dimension(:, :, :, :), allocatable :: rho ! Density fluxes
      real, dimension(:, :, :, :), allocatable :: u ! Zonal-wind fluxes
      real, dimension(:, :, :, :), allocatable :: v ! Meridional-wind fluxes
      real, dimension(:, :, :, :), allocatable :: w ! Vertical-wind fluxes
      real, dimension(:, :, :, :), allocatable :: rhop ! Dens-fluct. fluxes
    end type flux_type
  
    !-----------------------------------------------------------------
    !                    Input / Output variables
    !-----------------------------------------------------------------
    integer :: iOut ! output counter ; gagarina: moved from pinc
  
    character(len = 256) :: file_namelist
    integer, parameter :: nVar = 20 ! Maximum number of I/O variables
  
    character(len = 10), dimension(1:20) :: atmvarOut
  
    logical :: prepare_restart
    logical :: restart
  
    integer :: iIn
  
    character(len = 40) :: runName
  
    character(len = 20) :: outputType ! "time" or "timeStep"
    integer :: nOutput ! output every nOutput's time step
    integer :: maxIter ! max nb. of time steps
    real :: outputTimeDiff ! output every ... seconds
    real :: maxTime ! max time in seconds
  
    logical :: detailedinfo
    logical :: RHS_diagnostics
  
    logical :: fancy_namelists
  
    namelist / outputList / atmvarOut, &
        &prepare_restart, restart, iIn, runName, outputType, nOutput, maxIter, &
        &outputTimeDiff, maxTime, detailedinfo, RHS_diagnostics, fancy_namelists
    !achatzb
    !achatze
  
    !-----------------------------------------------------------------
    !                    Debugging & Error handling
    !-----------------------------------------------------------------
  
    logical :: verbose
    real :: dtMin_dim
    namelist / debuggingList / verbose, dtMin_dim
  
    !-----------------------------------------------------------------
    !                           Test cases
    !-----------------------------------------------------------------
  
    ! general
    character(len = 50) :: testCase
    namelist / testCaseList / testCase
  
  
    ! vertical direction
  
    character(len = 25) :: model
  
    namelist / modelList / model
  
    !-----------------------------------------------------------------
    !                               Solver
    !-----------------------------------------------------------------
  
    real :: cfl
    real :: dtMax_dim
  
    character(len = 20) :: tStepChoice ! "cfl", "fix"
    character(len = 20) :: timeScheme ! LS_Will_RK3 / Euler
    character(len = 20) :: timeSchemeType ! lowStorage / classical
    character(len = 20) :: limiterType1 ! minmod / ...
    logical :: auxil_equ ! auxiliary equation for the
    ! density fluctuations to be
    ! used in the explicit
    ! integration of the
    ! pseudo-incompressible system
    logical :: dtWave_on ! gagarina:
    ! time step for gravity wave
    ! resolution
    namelist / solverList / cfl, dtMax_dim, tStepChoice, timeScheme, &
        &auxil_equ, limiterType1, dtWave_on
    !UAC & dens_relax, shap_dts_dim, n_shap
  
    integer :: nStages
    logical :: updateMass ! transport of mass=var(1)  on/off
    logical :: predictMomentum ! transport of momentum=var(2-4) on/off
  
    !-----------------------------------------------------------------
    !                          Poisson solver
    !-----------------------------------------------------------------
  
    real :: tolPoisson
    real :: tolCond, tolref, abs_tol, scaled_atol, alpha_tol, b_norm
    integer :: maxIterPoisson
    character(len = 10) :: preconditioner
    character(len = 10) :: tolcrit
    real :: dtau
    integer :: maxIterADI
    logical :: initialCleaning
    logical :: correctMomentum ! false -> momentumCorrector off
    logical :: correctDivError ! true -> subtract rho*div(u)
    namelist / poissonSolverList / tolPoisson, abs_tol, tolCond, maxIterPoisson, &
        &preconditioner, dtau, maxIterADI, initialCleaning, correctMomentum, &
        &correctDivError, tolcrit
    integer :: nnz ! number of nonzeros
  
    ! hypre and bicgstab objects
    ! integer * 8 grid_hypre, stencil_e, stencil_i, A_hp_e, A_hp_i, b_hp_e, &
    !     b_hp_i, x_hp_e, x_hp_i, solver_hp_e, solver_hp_i
    ! real, dimension (:), allocatable :: values_e
    ! real, dimension (:), allocatable :: values_i
    !UAC real, dimension(:,:,:), allocatable :: ac_b, al_b,ar_b, ab_b,af_b, &
    real, dimension(:, :, :), allocatable :: ac_b, acv_b, ach_b, al_b, ar_b, &
        &ab_b, af_b, ad_b, au_b, alb_b, alf_b, arb_b, arf_b
    !UAE
  
    ! TFC FJ
    real, dimension(:, :, :), allocatable :: aru_b, ard_b, alu_b, ald_b, afu_b, &
        &afd_b, abu_b, abd_b, auu_b, add_b, aruu_b, ardd_b, aluu_b, aldd_b, &
        &afuu_b, afdd_b, abuu_b, abdd_b
  
    !-----------------------------------------------------------------
    !                           Constants
    !-----------------------------------------------------------------
  
    real :: pi ! you know...
    real, parameter :: small = 1.0e-20 ! to devision by zero
    complex, parameter :: imag = (0.0, 1.0) ! imaginary unit
  
    !-----------------------------------------------------------------
    !                         Flux specification
    !-----------------------------------------------------------------
  
    ! ILES (ALDM) parameter
    real :: sigmaC
    real :: sigma0
    real :: sigmaX, sigmaY, sigmaZ
  
    !-----------------------------------------------------------------
    !                             Atmosphere
    !-----------------------------------------------------------------
  
    character(len = 10) :: referenceQuantities ! set of reference quantities
    ! for the hydrostatically
    ! balanced background state
    logical :: specifyReynolds ! choose Reynolds or mu_viscous
    real :: ReInv ! reciprocal Reynolds number
    real :: mu_viscous_dim ! kinematic viscosity
    character(len = 30) :: background ! isentropic / isothermal /
    ! const-N / diflapse / HeldSuarez
    real :: Temp0_dim ! isoth. backgr. temp. in K
    real :: press0_dim ! pressure at z=0 in Pa
    real, dimension(3) :: backgroundFlow_dim
    real :: f_Coriolis_dim ! Coriolis parameter
    character(len = 30) :: corset ! constant / periodic
  
    namelist / atmosphereList / referenceQuantities, specifyReynolds, ReInv, &
        &mu_viscous_dim, background, &
        &Temp0_dim, press0_dim, backgroundFlow_dim, f_Coriolis_dim, &
        &corset
  
    real, dimension(3) :: backgroundFlow
    real :: theta00, rho00, P00 ! background values for Boussinesq
  
    !-----------------------------------------------------------------
    !                         Bottom topography
    !-----------------------------------------------------------------
  
    logical :: topography ! via k = 1
  
    ! Resolved topography
    real, dimension(:, :), allocatable :: topography_surface, &
        &final_topography_surface
  
    ! Layers
    real, dimension(:, :, :), allocatable :: zTFC, zTildeTFC
  
    ! Stretched vertical grid
    real, dimension(:), allocatable :: zS, zTildeS
  
    integer :: ipolTFC
  
    real :: topographyTime
  
    real :: mountainHeight_dim
    real :: mountainWidth_dim
    integer :: mountain_case
    real :: range_factor
    integer :: spectral_modes
    real :: envelope_reduction
    real :: stretch_exponent
  
    namelist / topographyList / topography, ipolTFC, topographyTime, &
        &mountainHeight_dim, mountainWidth_dim, mountain_case, &
        &range_factor, spectral_modes, envelope_reduction, stretch_exponent
    !UAB
    !UAE
  
    !-----------------------------------------------------------------
    !                         Boundary
    !-----------------------------------------------------------------
  
    ! sponge layer
    logical :: spongeLayer, sponge_uv
    real :: spongeHeight
    integer :: kSponge
    real :: zSponge
    real :: spongeAlphaZ_dim, spongeAlphaZ_fac
  
    ! Unified and lateral sponge layers (FJJun2023)
    logical :: unifiedSponge, lateralSponge
    real, dimension(:, :, :), allocatable :: alphaUnifiedSponge
    real :: xSponge0, ySponge0, xSponge1, ySponge1
    real :: dxSponge, dySponge, dzSponge
  
    ! Vertical sponge layer
    character(len = 50) :: spongeType
  
    ! Order of polynomial sponge
    integer :: spongeOrder
  
    ! Damping time of COSMO sponge (in time steps)
    integer :: cosmoSteps
  
    logical :: relax_to_mean
  
    real :: relaxation_period
    real :: relaxation_amplitude
  
    ! gaga: backup, delete later
    real, dimension(:, :), allocatable :: u_const ! constant wind for baroclinic life cycle experiments
  
    namelist / boundaryList / spongeLayer, sponge_uv, &
        &spongeHeight, spongeAlphaZ_dim, spongeAlphaZ_fac, unifiedSponge, &
        &lateralSponge, spongeType, spongeOrder, cosmoSteps, relax_to_mean, &
        &relaxation_period, relaxation_amplitude
  
    ! boundary types
    character(len = 15) :: xBoundary
    character(len = 15) :: yBoundary
    character(len = 15) :: zBoundary
  
    namelist / boundaryList2 / xBoundary, yBoundary, zBoundary
  
    !-----------------------------------------------------------------
    !                     WKB arrays and variables
    !-----------------------------------------------------------------
  
  
    contains
  
    subroutine default_values
  
      !-----------------------------------------------------------------
      !                     Set default values
      !-----------------------------------------------------------------
  
      ! Input/Output
      atmvarOut = ""
      prepare_restart = .false.
      restart = .false.
      iIn = - 1
      runName = "runName"
      outputType = "time"
      nOutput = 1
      maxIter = 1
      outputTimeDiff = 3600.0
      maxTime = 3600.0
      detailedinfo = .false.
      RHS_diagnostics = .false.
      fancy_namelists = .true.
  
      ! Debugging
      dtMin_dim = 1.0e-6
  
      ! Test cases
      testCase = "mountainwave"
  
      ! Model equations
      model = "pseudo_incompressible"
  
      ! Solver
      cfl = 0.5
      dtMax_dim = 1.0e3
      tStepChoice = "cfl"
      timeScheme = "LS_Will_RK3"
      auxil_equ = .false.
      limiterType1 = "MCVariant"
      dtWave_on = .true.
  
      ! Poisson solver
      tolPoisson = 1.0e-8
      abs_tol = 0.0
      tolCond = 1.0e-23
      maxIterPoisson = 1000
      preconditioner = "yes"
      dtau = 4.0e-4
      maxIterADI = 2
      initialCleaning = .true.
      correctMomentum = .true.
      correctDivError = .false.
      tolcrit = "abs"
  
      ! Atmosphere
      referenceQuantities = "Klein"
      specifyReynolds = .false.
      ReInv = 0.0
      mu_viscous_dim = 0.0
      background = "isothermal"
      Temp0_dim = 300.0
      press0_dim = 100000.0
      backgroundFlow_dim = 0.0
      f_Coriolis_dim = 0.0
      corset = "constant"
  
      ! Topography
      topography = .false.
      ipolTFC = 2
      topographyTime = 0.0
      mountainHeight_dim = 0.1 * (lz_dim(1) - lz_dim(0))
      mountainWidth_dim = 0.1 * (lx_dim(1) - lx_dim(0))
      mountain_case = 1
      range_factor = 1.0
      spectral_modes = 1
      envelope_reduction = 0.0
      stretch_exponent = 1.0
  
      ! Boundaries
      spongeLayer = .false.
      sponge_uv = .false.
      spongeHeight = 0.33
      spongeAlphaZ_dim = 0.01
      spongeAlphaZ_fac = 0.01
      unifiedSponge = .false.
      lateralSponge = .false.
      spongeType = "polynomial"
      spongeOrder = 1
      cosmoSteps = 1
      relax_to_mean = .true.
      relaxation_period = 0.0
      relaxation_amplitude = 0.0
      xBoundary = "periodic"
      yBoundary = "periodic"
      zBoundary = "solid_wall"
  
    end subroutine default_values
  
    subroutine write_namelists
  
      ! Write all namelists.
  
      implicit none
  
      integer, parameter :: parameter_space = 24, value_space = 23, &
          &comment_space = 72 - parameter_space - value_space
      character(len = 50) :: character_format, integer_format, logical_format, &
          &comment_format
      character(len = 2) :: counter
      integer :: iVar, jVar
  
      ! Adjust atmVarOut.
      jVar = 0
      do iVar = 1, size(atmVarOut)
        if(atmVarOut(iVar) == "" .or. (iVar > 1 .and. any(atmVarOut(:iVar - 1) &
            &== atmVarOut(iVar)))) cycle
        jVar = jVar + 1
        atmVarOut(jVar) = atmVarOut(iVar)
      end do
      atmVarOut(jVar + 1:) = ""
  
      ! Write namelists in standard format.
      if(master .and. .not. fancy_namelists) then
        open(unit = 90, file = "namelists.txt", action = "write", form &
            &= "formatted", status = "replace")
        write(unit = 90, nml = domain)
        write(90, "(a)") ""
        write(unit = 90, nml = outputList)
        write(90, "(a)") ""
        write(unit = 90, nml = debuggingList)
        write(90, "(a)") ""
        write(unit = 90, nml = testCaseList)
        write(90, "(a)") ""
        write(unit = 90, nml = modelList)
        write(90, "(a)") ""
        write(unit = 90, nml = solverList)
        write(90, "(a)") ""
        write(unit = 90, nml = poissonSolverList)
        write(90, "(a)") ""
        write(unit = 90, nml = atmosphereList)
        write(90, "(a)") ""
        write(unit = 90, nml = topographyList)
        write(90, "(a)") ""
        write(unit = 90, nml = boundaryList)
        write(90, "(a)") ""
        write(unit = 90, nml = boundaryList2)
        write(90, "(a)") ""
        close(90)
        return
      end if
  
      ! Check if comment space is too small.
      if(comment_space <= 0) stop "Comment space too small! Decrease parameter &
          &space and/or value space!"
  
      ! Define character, integer and logical formats.
      write(counter, "(i2)") parameter_space
      character_format = "(2x, a" // trim(adjustl(counter)) // ", a3, a"
      integer_format = "(2x, a" // trim(adjustl(counter)) // ", a3, i"
      logical_format = "(2x, a" // trim(adjustl(counter)) // ", a3, l"
      write(counter, "(i2)") value_space
      character_format = trim(adjustl(character_format)) &
          &// trim(adjustl(counter)) // ", a3, a)"
      integer_format = trim(adjustl(integer_format)) // trim(adjustl(counter)) &
          &// ", a3, a)"
      logical_format = trim(adjustl(logical_format)) // trim(adjustl(counter)) &
          &// ", a3, a)"
  
      ! Define comment format.
      write(counter, "(i2)") parameter_space + value_space + 5
      comment_format = "(" // trim(adjustl(counter)) // "x, a3, a)"
  
      ! Open info file.
      if(master) then
        open(unit = 90, file = "namelists.txt", action = "write", form &
            &= "formatted", status = "replace")
  
        ! Write the header.
        write(90, "(a)") "!" // repeat("-", 78) // "!"
        write(90, "(a)") "!" // repeat(" ", 30) // "PincFlow Namelists" &
            &// repeat(" ", 30) // "!"
        write(90, "(a)") "!" // repeat("-", 78) // "!"
        write(90, "(a)") ""
  
        ! Write domain namelist.
        write(90, "(a)") "&domain"
        call write_integer("sizeX", sizeX, "Cells in x")
        call write_integer("sizeY", sizeY, "Cells in y")
        call write_integer("sizeZ", sizeZ, "Cells in z")
        call write_integer("nbx", nbx, "Halo/ghost cells in x")
        call write_integer("nby", nby, "Halo/ghost cells in y")
        call write_integer("nbz", nbz, "Ghost cells in z")
        call write_float("lx_dim(0)", lx_dim(0), "Minimum of x")
        call write_float("lx_dim(1)", lx_dim(1), "Maximum of x")
        call write_float("ly_dim(0)", ly_dim(0), "Minimum of y")
        call write_float("ly_dim(1)", ly_dim(1), "Maximum of y")
        call write_float("lz_dim(0)", lz_dim(0), "Minimum of z")
        call write_float("lz_dim(1)", lz_dim(1), "Maximum of z")
        call write_integer("nprocx", nprocx, "Processors in x")
        call write_integer("nprocy", nprocy, "Processors in y")
        write(90, "(a)") "&end"
        write(90, "(a)") ""
  
        ! Write output namelist.
        write(90, "(a)") "&outputList"
        do iVar = 1, size(atmvarOut)
          if(iVar == 1 .or. atmvarOut(iVar) /= "") then
            write(counter, "(i2)") iVar
            call write_character("atmvarOut(" // trim(adjustl(counter)) // ")", &
                &atmvarOut(iVar), "Atmospheric output variable")
          end if
        end do
        call write_logical("prepare_restart", prepare_restart, "Save everything &
            &needed for a restart")
        call write_logical("restart", restart, "Restart the model from state in &
            &previous simulation")
        call write_integer("iIn", iIn, "Restart at time step iIn")
        call write_character("runName", runName, "Run name for netCDF file")
        call write_character("outputType", outputType, "'timeStep' or 'time'")
        call write_integer("nOutput", nOutput, "Output every nOutput time steps &
            &for outputType = 'timeStep'")
        call write_integer("maxIter", maxIter, "Stop after maxIter time steps &
            &for outputType = 'timeStep'")
        call write_float("outputTimeDiff", outputTimeDiff, "Output every &
            &outputTimeDiff seconds for outputType = 'time'")
        call write_float("maxTime", maxTime, "Stop after maxTime seconds for &
            &outputType = 'time'")
        call write_logical("detailedinfo", detailedinfo, "Provide info on the &
            &final state of Poisson solver")
        call write_logical("RHS_diagnostics", RHS_diagnostics, "Provide info &
            &about the RHS of Poisson equation")
        call write_logical("fancy_namelists", fancy_namelists, "Write all &
            &namelists with comments")
        write(90, "(a)") "&end"
        write(90, "(a)") ""
  
        ! Write debugging namelist.
        write(90, "(a)") "&debuggingList"
        call write_float("dtMin_dim", dtMin_dim, "Stop if dt < dtMin")
        write(90, "(a)") "&end"
        write(90, "(a)") ""
  
        ! Write testcase namelist.
        write(90, "(a)") "&testCaseList"
        call write_character("testCase", testCase, "Predefined test case")
        write(90, "(a)") "&end"
        write(90, "(a)") ""
  
        ! Write model namelist.
        write(90, "(a)") "&modelList"
        call write_character("model", model, "Dynamic equations")
        write(90, "(a)") "&end"
        write(90, "(a)") ""
  
        ! Write solver namelist.
        write(90, "(a)") "&solverList"
        call write_float("cfl", cfl, "CFL number")
        all write_float("dtMax_dim", dtMax_dim, "Maximum time step")
        call write_character("tStepChoice", tStepChoice, "'fix' or 'cfl'")
        call write_character("timeScheme", timeScheme, "'LS_Will_RK3' or &
            &'semiimplicit'")
        call write_logical("auxil_equ", auxil_equ, "Buoyancy equation")
        call write_character("limiterType1", limiterType1, "'minmod', &
            &'MCVariant' or 'Cada'")
        call write_logical("dtWave_on", dtWave_on, "Limit time step by inverse &
            &buoyancy frequency")
        write(90, "(a)") "&end"
        write(90, "(a)") ""
  
        ! Write Poisson solver namelist.
        write(90, "(a)") "&poissonSolverList"
        call write_float("tolPoisson", tolPoisson, "Abort criterion")
        call write_float("abs_tol", abs_tol, "Lower bound for tolerance")
        call write_float("tolCond", tolCond, "Preconditioner tolerance")
        call write_integer("maxIterPoisson", maxIterPoisson, "Maximum iterations")
        call write_character("preconditioner", preconditioner, "'no' or 'yes'")
        call write_float("dtau", dtau, "Time parameter for preconditioner")
        call write_integer("maxIterADI", maxIterADI, "Preconditioner iterations")
        call write_logical("initialCleaning", initialCleaning, "Enforce initial &
            &non-divergence")
        call write_logical("correctMomentum", correctMomentum, "Correct momentum &
            &so that divergence constraint is fulfilled")
        call write_logical("correctDivError", correctDivError, "Subtract &
            &divergence")
        call write_character("tolcrit", tolcrit, "'abs' or 'rel'")
        write(90, "(a)") "&end"
        write(90, "(a)") ""
  
        ! Write atmosphere namelist.
        write(90, "(a)") "&atmosphereList"
        call write_character("referenceQuantities", referenceQuantities, &
            &"'Klein'")
        call write_logical("specifyReynolds", specifyReynolds, "Use inverse &
            &Reynolds number")
        call write_float("ReInv", ReInv, "Inverse Reynolds number")
        call write_float("mu_viscous_dim", mu_viscous_dim, "Kinematic viscosity")
        call write_character("background", background, "'isothermal'")
        call write_float("Temp0_dim", Temp0_dim, "Background temperature for &
            &'isothermal'")
        call write_float("press0_dim", press0_dim, "Ground pressure")
        do iVar = 1, 3
          write(counter, "(i2)") iVar
          call write_float("backgroundFlow_dim(" // trim(adjustl(counter)) &
              &// ")", backgroundFlow_dim(iVar), "Initial wind")
        end do
        call write_float("f_Coriolis_dim", f_Coriolis_dim, "Coriolis frequency")
        call write_character("corset", corset, "'constant'")
        write(90, "(a)") "&end"
        write(90, "(a)") ""
  
        ! Write topography namelist.
        write(90, "(a)") "&topographyList"
        call write_logical("topography", topography, "Terrain-following &
            &coordinates")
        call write_integer("ipolTFC", ipolTFC, "Interpolation in the &
            &transformation of w")
        call write_float("topographyTime", topographyTime, "Topography growth &
            &time")
        call write_float("mountainHeight_dim", mountainHeight_dim, "Maximum &
            &height")
        call write_float("mountainWidth_dim", mountainWidth_dim, "Half width")
        call write_integer("mountain_case", mountain_case, "Predefined &
            &topography")
        call write_float("range_factor", range_factor, "Ratio between large and &
            &small scales")
        call write_integer("spectral_modes", spectral_modes, "Number of spectral &
            &modes")
        call write_float("envelope_reduction", envelope_reduction, "Relative &
            &reduction of the envelope (between 0 and 1)")
        call write_float("stretch_exponent", stretch_exponent, "Exponent of &
            &vertical grid stretching (1 for no stretching)")
        write(90, "(a)") "&end"
        write(90, "(a)") ""
  
        ! Write boundary namelist.
        write(90, "(a)") "&boundaryList"
        call write_logical("spongeLayer", spongeLayer, "General sponge layer &
            &switch")
        call write_logical("sponge_uv", sponge_uv, "Sponge layer for horizontal &
            &wind if unifiedSponge = .false.")
        call write_float("spongeHeight", spongeHeight, "Relative height of lower &
            &sponge layer edge (scale height for unifiedSponge = .true. and &
            &spongeType = 'exponential')")
        call write_float("spongeAlphaZ_dim", spongeAlphaZ_dim, "Maximum &
            &relaxation rate for unifiedSponge = .true.")
        call write_float("spongeAlphaZ_fac", spongeAlphaZ_fac, "Sponge layer &
            &factor for unifiedSponge = .false.")
        call write_logical("unifiedSponge", unifiedSponge, "Unified sponge for &
            &both time schemes, applied to wind and density")
        call write_logical("lateralSponge", lateralSponge, "Lateral sponge for &
            &unifiedSponge = .true.")
        call write_character("spongeType", spongeType, "Sponge layer profile for &
            &unifiedSponge = .true.")
        call write_integer("spongeOrder", spongeOrder, "Order of polynomial &
            &sponge")
        call write_integer("cosmoSteps", cosmoSteps, "Relative strength of COSMO &
            &sponge")
        call write_logical("relax_to_mean", relax_to_mean, "Relax the wind to &
            &its (terrain-following) horizontal mean (otherwise, relax to the &
            &initial state) if unifiedSponge == .true.")
        call write_float("relaxation_period", relaxation_period, "Period of &
            &an oscillation superposed on the background wind if unifiedSponge &
            &== .true. and relax_to_mean == .false. (0 for no oscillation)")
        call write_float("relaxation_amplitude", relaxation_amplitude, &
            &"Relative amplitude of an oscillation superposed on the background &
            &wind if unifiedSponge == .true. and relax_to_mean == .false.")
        write(90, "(a)") "&end"
        write(90, "(a)") ""
  
        ! Write second boundary namelist.
        write(90, "(a)") "&boundaryList2"
        call write_character("xBoundary", xBoundary, "Boundary conditions in x &
            &('periodic' only)")
        call write_character("yBoundary", yBoundary, "Boundary conditions in y &
            &('periodic' only)")
        call write_character("zBoundary", zBoundary, "Boundary conditions in z &
            &('periodic' or 'solid_wall')")
        write(90, "(a)") "&end"
        write(90, "(a)") ""
  
        ! Close info file.
        close(90)
      end if
  
      contains
  
      subroutine write_character(parameter, value, comment)
  
        character(len = *), intent(in) :: parameter, value, comment
        character(len = len_trim(adjustl(parameter))) :: parameter_trimmed
        character(len = len_trim(adjustl(value))) :: value_trimmed
        character(len = len_trim(adjustl(comment))) :: comment_trimmed
        character(len = parameter_space) :: parameter_output
        character(len = value_space) :: value_output
        character(len = comment_space) :: comment_output
        integer :: position
  
        ! Trim input.
        parameter_trimmed = trim(adjustl(parameter))
        value_trimmed = trim(adjustl(value))
        comment_trimmed = trim(adjustl(comment))
  
        ! Check if parameter space and value space are large enough.
        if(len(parameter_trimmed) > parameter_space) stop "Parameter space too &
            &small!"
        if(len(value_trimmed) + 2 > value_space) stop "Value space too small!"
  
        ! Set parameter output and value output.
        parameter_output = adjustl(parameter_trimmed)
        value_output = repeat(" ", value_space - len(value_trimmed) - 2) // "'" &
            &// value_trimmed // "'"
  
        ! Set comment output.
        if(len(comment_trimmed) > comment_space) then
          position = max(index(comment_trimmed(:comment_space), " ", back &
              &= .true.), index(comment_trimmed(:comment_space), "-", back &
              &= .true.), index(comment_trimmed(:comment_space), "/", back &
              &= .true.))
          if(position == 0) stop "Hyphenation required!"
          comment_output = adjustl(comment_trimmed(:position))
        else
          comment_output = adjustl(comment_trimmed)
        end if
  
        ! Write namelist line.
        write(90, character_format) parameter_output, " = ", value_output, " ! &
            &", trim(comment_output)
  
        ! Continue comment.
        if(len(comment_trimmed) > comment_space) then
          call write_comment(comment_trimmed, position)
        end if
  
      end subroutine write_character
  
      subroutine write_integer(parameter, value, comment)
  
        character(len = *), intent(in) :: parameter, comment
        integer, intent(in) :: value
        character(len = len_trim(adjustl(parameter))) :: parameter_trimmed
        character(len = len_trim(adjustl(comment))) :: comment_trimmed
        character(len = parameter_space) :: parameter_output
        character(len = comment_space) :: comment_output
        integer :: position
  
        ! Trim input.
        parameter_trimmed = trim(adjustl(parameter))
        comment_trimmed = trim(adjustl(comment))
  
        ! Check if parameter space and value space are large enough.
        if(len(parameter_trimmed) > parameter_space) stop "Parameter space too &
            &small!"
        if(integer_size(value) > value_space) stop "Value space too small!"
  
        ! Set parameter output.
        parameter_output = adjustl(parameter_trimmed)
  
        ! Set comment output.
        if(len(comment_trimmed) > comment_space) then
          position = max(index(comment_trimmed(:comment_space), " ", back &
              &= .true.), index(comment_trimmed(:comment_space), "-", back &
              &= .true.), index(comment_trimmed(:comment_space), "/", back &
              &= .true.))
          if(position == 0) stop "Hyphenation required!"
          comment_output = adjustl(comment_trimmed(:position))
        else
          comment_output = adjustl(comment_trimmed)
        end if
  
        ! Write namelist line.
        write(90, integer_format) parameter_output, " = ", value, " ! ", &
            &trim(comment_output)
  
        ! Continue comment.
        if(len(comment_trimmed) > comment_space) then
          call write_comment(comment_trimmed, position)
        end if
  
      end subroutine write_integer
  
      subroutine write_float(parameter, value, comment)
  
        character(len = *), intent(in) :: parameter, comment
        real, intent(in) :: value
        character(len = len_trim(adjustl(parameter))) :: parameter_trimmed
        character(len = len_trim(adjustl(comment))) :: comment_trimmed
        character(len = parameter_space) :: parameter_output
        character(len = 50) :: float_format
        character(len = comment_space) :: comment_output
        character(len = 2) counter
        integer, dimension(1:2) :: digit_count
        integer :: position
  
        ! Trim input.
        parameter_trimmed = trim(adjustl(parameter))
        comment_trimmed = trim(adjustl(comment))
  
        ! Compute significant digits.
        digit_count = significant_digits(value)
  
        ! Check if parameter space and value space are large enough.
        if(len(parameter_trimmed) > parameter_space) stop "Parameter space too &
            &small!"
        if(sum(digit_count) + 5 > value_space) stop "Value space too small!"
  
        ! Set parameter output.
        parameter_output = adjustl(parameter_trimmed)
  
        ! Define float format.
        write(counter, "(i2)") parameter_space
        float_format = "(2x, a" // trim(adjustl(counter)) // ", a3, es"
        write(counter, "(i2)") value_space
        float_format = trim(adjustl(float_format)) // trim(adjustl(counter)) &
            &// "."
        write(counter, "(i2)") digit_count(1)
        float_format = trim(adjustl(float_format)) // trim(adjustl(counter)) &
            &// "e"
        write(counter, "(i2)") digit_count(2)
        float_format = trim(adjustl(float_format)) // trim(adjustl(counter)) &
            &// ", a3, a)"
  
        ! Set comment output.
        if(len(comment_trimmed) > comment_space) then
          position = max(index(comment_trimmed(:comment_space), " ", back &
              &= .true.), index(comment_trimmed(:comment_space), "-", back &
              &= .true.), index(comment_trimmed(:comment_space), "/", back &
              &= .true.))
          if(position == 0) stop "Hyphenation required!"
          comment_output = adjustl(comment_trimmed(:position))
        else
          comment_output = adjustl(comment_trimmed)
        end if
  
        ! Write namelist line.
        write(90, float_format) parameter_output, " = ", value, " ! ", &
            &trim(comment_output)
  
        ! Continue comment.
        if(len(comment_trimmed) > comment_space) then
          call write_comment(comment_trimmed, position)
        end if
  
      end subroutine write_float
  
      subroutine write_logical(parameter, value, comment)
  
        character(len = *), intent(in) :: parameter, comment
        logical, intent(in) :: value
        character(len = len_trim(adjustl(parameter))) :: parameter_trimmed
        character(len = len_trim(adjustl(comment))) :: comment_trimmed
        character(len = parameter_space) :: parameter_output
        character(len = comment_space) :: comment_output
        integer :: position
  
        ! Trim input.
        parameter_trimmed = trim(adjustl(parameter))
        comment_trimmed = trim(adjustl(comment))
  
        ! Check if parameter space is large enough.
        if(len(parameter_trimmed) > parameter_space) stop "Parameter space too &
            &small!"
  
        ! Set parameter output.
        parameter_output = adjustl(parameter_trimmed)
  
        ! Set comment output.
        if(len(comment_trimmed) > comment_space) then
          position = max(index(comment_trimmed(:comment_space), " ", back &
              &= .true.), index(comment_trimmed(:comment_space), "-", back &
              &= .true.), index(comment_trimmed(:comment_space), "/", back &
              &= .true.))
          if(position == 0) stop "Hyphenation required!"
          comment_output = adjustl(comment_trimmed(:position))
        else
          comment_output = adjustl(comment_trimmed)
        end if
  
        ! Write namelist line.
        write(90, logical_format) parameter_output, " = ", value, " ! ", &
            &trim(comment_output)
  
        ! Continue comment.
        if(len(comment_trimmed) > comment_space) then
          call write_comment(comment_trimmed, position)
        end if
  
      end subroutine write_logical
  
      subroutine write_comment(comment, position)
  
        character(len = *), intent(in) :: comment
        integer, intent(in) :: position
        character(len = comment_space) :: comment_output
        integer :: left, right, shift
  
        left = position + 1
        do while(len_trim(comment(left:)) > 0)
          comment_output = adjustl(comment(left:))
          shift = max(index(comment_output, " ", back = .true.), &
              &index(comment_output, "-", back = .true.), index(comment_output, &
              &"/", back = .true.))
          if(shift == 0) stop "Hyphenation required!"
          comment_output = adjustl(comment_output(:shift))
          write(90, comment_format) " ! ", trim(comment_output)
          right = left + shift - 1
          left = right + 1
        end do
  
      end subroutine write_comment
  
    end subroutine write_namelists
  
    function trim_float(input) result(output)
  
      real, intent(in) :: input
      integer, dimension(1:2) :: digit_count
      character(len = 50) :: float_format
      character(len = 2) :: counter
      character(len = max(0, int(sign(1.0, - input))) &
          &+ sum(significant_digits(input)) + 4) :: output
  
      digit_count = significant_digits(input)
  
      write(counter, "(i2)") len(output)
      float_format = "(es" // trim(adjustl(counter)) // "."
      write(counter, "(i2)") digit_count(1)
      float_format = trim(adjustl(float_format)) // trim(adjustl(counter)) // "e"
      write(counter, "(i2)") digit_count(2)
      float_format = trim(adjustl(float_format)) // counter // ")"
  
      write(output, float_format) input
  
    end function trim_float
  
    function trim_integer(input) result(output)
  
      integer, intent(in) :: input
      character(len = 20) :: number
      character(len = integer_size(input)) :: output
  
      write(number, "(i20)") input
      output = trim(adjustl(number))
  
    end function trim_integer
  
    function trim_logical(input) result(output)
  
      logical, intent(in) :: input
      character(len = 1) :: output
  
      write(output, "(l1)") input
  
    end function trim_logical
  
    pure function convert_case(input, choice) result(output)
  
      ! Convert either lower to upper or upper to lower case.
  
      implicit none
  
      character(len = *), intent(in) :: input, choice
      character(len = len(input)) :: output
      character(len = *), parameter :: lower = "abcdefghijklmnopqrstuvwxyz"
      character(len = *), parameter :: upper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      integer :: i, j
  
      output = input
      if(choice == "upper") then
        do i = 1, len(output)
          j = index(lower, output(i:i))
          if(j > 0) output(i:i) = upper(j:j)
        end do
      else if(choice == "lower") then
        do i = 1, len(output)
          j = index(upper, output(i:i))
          if(j > 0) output(i:i) = lower(j:j)
        end do
      end if
  
    end function
  
    pure function significant_digits(input) result(digit_count)
  
      ! Determine the significant digits of a 64-bit float (the 16th digit is
      ! rounded).
  
      implicit none
  
      real, intent(in) :: input
      character(len = 23) :: number
      integer, dimension(1:2) :: digit_count
      integer :: position
  
      write(number, "(es23.15e3)") input
      digit_count = 0
      do position = 21, 23
        if(number(position:position) /= "0") exit
        digit_count(2) = digit_count(2) + 1
      end do
      digit_count(2) = max(1, 3 - digit_count(2))
      do position = 18, 4, - 1
        if(number(position:position) /= "0") exit
        digit_count(1) = digit_count(1) + 1
      end do
      digit_count(1) = max(1, 15 - digit_count(1))
  
    end function significant_digits
  
    pure function integer_size(input) result(size)
  
      ! Determine the size of a 64-bit integer.
  
      implicit none
  
      integer, intent(in) :: input
      character(len = 20) :: number
      integer :: size
  
      write(number, "(i20)") input
      size = len_trim(adjustl(number))
  
    end function integer_size
  
    subroutine allocate_var_type(var)
  
      implicit none
  
      type(var_type), intent(inout) :: var
      integer :: allocstat
  
      ! Allocate basic variables.
      allocate(var%rho(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "Allocation of var_type component rho failed!"
      allocate(var%u(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "Allocation of var_type component u failed!"
      allocate(var%v(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "Allocation of var_type component v failed!"
      allocate(var%w(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "Allocation of var_type component w failed!"
      allocate(var%pi(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "Allocation of var_type component pi failed!"
      allocate(var%rhop(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "Allocation of var_type component rhop failed!"
  
  
    end subroutine allocate_var_type
  
    subroutine allocate_flux_type(flux)
  
      implicit none
  
      type(flux_type), intent(inout) :: flux
      integer :: allocstat
  
      ! Allocate basic fluxes.
      allocate(flux%rho(- 1:nx, - 1:ny, - 1:nz, 1:3), stat = allocstat)
      if(allocstat /= 0) stop "Allocation of flux_type component rho failed!"
      allocate(flux%u(- 1:nx, - 1:ny, - 1:nz, 1:3), stat = allocstat)
      if(allocstat /= 0) stop "Allocation of flux_type component u failed!"
      allocate(flux%v(- 1:nx, - 1:ny, - 1:nz, 1:3), stat = allocstat)
      if(allocstat /= 0) stop "Allocation of flux_type component v failed!"
      allocate(flux%w(- 1:nx, - 1:ny, - 1:nz, 1:3), stat = allocstat)
      if(allocstat /= 0) stop "Allocation of flux_type component w failed!"
      allocate(flux%rhop(- 1:nx, - 1:ny, - 1:nz, 1:3), stat = allocstat)
      if(allocstat /= 0) stop "Allocation of flux_type component rhop failed!"
  
  
    end subroutine allocate_flux_type
  
    subroutine reset_var_type(var)
  
      implicit none
  
      type(var_type), intent(inout) :: var
  
      ! Reset basic variables.
      var%rho = 0.0
      var%u = 0.0
      var%v = 0.0
      var%w = 0.0
      var%pi = 0.0
      var%rhop = 0.0
  
    end subroutine reset_var_type
  
    subroutine reset_flux_type(flux)
  
      implicit none
  
      type(flux_type), intent(inout) :: flux
  
      ! Reset basic fluxes.
      flux%rho = 0.0
      flux%u = 0.0
      flux%v = 0.0
      flux%w = 0.0
      flux%rhop = 0.0
  
    end subroutine reset_flux_type
  
    subroutine deallocate_var_type(var)
  
      implicit none
  
      type(var_type), intent(inout) :: var
      integer :: allocstat
  
      ! Deallocate basic variables.
      deallocate(var%rho, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of var_type component rho failed!"
      deallocate(var%u, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of var_type component u failed!"
      deallocate(var%v, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of var_type component v failed!"
      deallocate(var%w, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of var_type component w failed!"
      deallocate(var%pi, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of var_type component pi failed!"
      deallocate(var%rhop, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of var_type component rhop failed!"
  
    end subroutine deallocate_var_type
  
    subroutine deallocate_flux_type(flux)
  
      implicit none
  
      type(flux_type), intent(inout) :: flux
      integer :: allocstat
  
      ! Deallocate basic fluxes.
      deallocate(flux%rho, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of flux_type component rho failed!"
      deallocate(flux%u, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of flux_type component u failed!"
      deallocate(flux%v, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of flux_type component v failed!"
      deallocate(flux%w, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of flux_type component w failed!"
      deallocate(flux%rhop, stat = allocstat)
      if(allocstat /= 0) stop "Deallocation of flux_type component rhop failed!"
  
    end subroutine deallocate_flux_type
  
  end module type_module
  