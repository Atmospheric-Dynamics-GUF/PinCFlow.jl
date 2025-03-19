!------------------------------------------------------------------------------!
!                              PincFlow Namelists                              !
!------------------------------------------------------------------------------!

&domain
  sizeX                    =                     300 ! Cells in x
  sizeY                    =                       1 ! Cells in y
  sizeZ                    =                     100 ! Cells in z
  nbx                      =                       3 ! Halo/ghost cells in x
  nby                      =                       3 ! Halo/ghost cells in y
  nbz                      =                       3 ! Ghost cells in z
  lx_dim(0)                =                  0.0E+0 ! Minimum of x
  lx_dim(1)                =                  6.0E+4 ! Maximum of x
  ly_dim(0)                =                  0.0E+0 ! Minimum of y
  ly_dim(1)                =                  4.0E+4 ! Maximum of y
  lz_dim(0)                =                  0.0E+0 ! Minimum of z
  lz_dim(1)                =                  2.0E+4 ! Maximum of z
  nprocx                   =                {nprocx} ! Processors in x
  nprocy                   =                {nprocy} ! Processors in y
&end

&outputList
  atmvarOut(1)             =                     'w' ! Atmospheric output
                                                     ! variable
  prepare_restart          =                       F ! Save everything needed
                                                     ! for a restart
  restart                  =                       F ! Restart the model from
                                                     ! state in previous
                                                     ! simulation
  iIn                      =                      -1 ! Restart at time step iIn
  runName                  =          'mountainwave' ! Run name for netCDF file
  outputType               =                  'time' ! 'timeStep' or 'time'
  nOutput                  =                       1 ! Output every nOutput
                                                     ! time steps for
                                                     ! outputType = 'timeStep'
  maxIter                  =                       1 ! Stop after maxIter time
                                                     ! steps for outputType =
                                                     ! 'timeStep'
  outputTimeDiff           =                  3.6E+3 ! Output every
                                                     ! outputTimeDiff seconds
                                                     ! for outputType = 'time'
  maxTime                  =                  3.6E+3 ! Stop after maxTime
                                                     ! seconds for outputType =
                                                     ! 'time'
  fancy_namelists          =                       T ! Write all namelists with
                                                     ! comments
&end

&debuggingList
  dtMin_dim                =                  1.0E-5 ! Stop if dt < dtMin
&end

&testCaseList
  testCase                 =          'mountainwave' ! Predefined test case
&end

&modelList
  model                    = 'pseudo_incompressible' ! Dynamic equations
&end

&solverList
  cfl                      =                  5.0E-1 ! CFL number
  dtMax_dim                =                  6.0E+1 ! Maximum time step
  tStepChoice              =                   'cfl' ! 'fix' or 'cfl'
  limiterType1             =             'MCVariant' ! Flux limiter
                                                     ! ('MCVariant')
&end

&poissonSolverList
  tolPoisson               =                  1.0E-8 ! Abort criterion
  maxIterPoisson           =                    5000 ! Maximum iterations
  preconditioner           =                   'yes' ! 'no' or 'yes'
  dtau                     =                  4.0E+0 ! Time parameter for
                                                     ! preconditioner
  maxIterADI               =                       2 ! Preconditioner iterations
  initialCleaning          =                       T ! Enforce initial non-
                                                     ! divergence
  correctMomentum          =                       T ! Correct momentum so that
                                                     ! divergence constraint is
                                                     ! fulfilled
  tolcrit                  =                   'abs' ! 'abs' or 'rel'
&end

&atmosphereList
  specifyReynolds          =                       F ! Use inverse Reynolds
                                                     ! number
  ReInv                    =                  0.0E+0 ! Inverse Reynolds number
  mu_viscous_dim           =                  0.0E+0 ! Kinematic viscosity
  background               =            'isothermal' ! 'isothermal'
  Temp0_dim                =                  3.0E+2 ! Background temperature
                                                     ! for 'isothermal'
  press0_dim               =                  1.0E+5 ! Ground pressure
  backgroundFlow_dim(1)    =                  1.0E+1 ! Initial wind
  backgroundFlow_dim(2)    =                  0.0E+0 ! Initial wind
  backgroundFlow_dim(3)    =                  0.0E+0 ! Initial wind
  f_Coriolis_dim           =                  0.0E+0 ! Coriolis frequency
  corset                   =              'constant' ! 'constant'
&end

&topographyList
  mountainHeight_dim       =                  4.0E+2 ! Maximum height
  mountainWidth_dim        =                  1.0E+3 ! Half width
  mountain_case            =                       3 ! Predefined topography
  range_factor             =                  1.0E+1 ! Ratio between large and
                                                     ! small scales
  spectral_modes           =                       1 ! Number of spectral modes
  envelope_reduction       =                  0.0E+0 ! Relative reduction of
                                                     ! the envelope (between 0
                                                     ! and 1)
  stretch_exponent         =                  1.0E+0 ! Exponent of vertical
                                                     ! grid stretching (1 for
                                                     ! no stretching)
&end

&boundaryList
  spongeLayer              =                       T ! General sponge layer
                                                     ! switch
  sponge_uv                =                       F ! Sponge layer for
                                                     ! horizontal wind if
                                                     ! unifiedSponge = .false.
  spongeHeight             =                  5.0E-1 ! Relative height of lower
                                                     ! sponge layer edge (scale
                                                     ! height for unifiedSponge
                                                     ! = .true. and spongeType
                                                     ! = 'exponential')
  spongeAlphaZ_dim         =                  1.0E-2 ! Maximum relaxation rate
                                                     ! for unifiedSponge =
                                                     ! .true.
  spongeAlphaZ_fac         =                  1.0E+0 ! Sponge layer factor for
                                                     ! unifiedSponge = .false.
  unifiedSponge            =                       T ! Unified sponge for both
                                                     ! time schemes, applied to
                                                     ! wind and density
  lateralSponge            =                       T ! Lateral sponge for
                                                     ! unifiedSponge = .true.
  spongeType               =            'sinusoidal' ! Sponge layer profile for
                                                     ! unifiedSponge = .true.
  spongeOrder              =                       1 ! Order of polynomial
                                                     ! sponge
  cosmoSteps               =                       1 ! Relative strength of
                                                     ! COSMO sponge
  relax_to_mean            =                       T ! Relax the wind to its
                                                     ! (terrain-following)
                                                     ! horizontal mean
                                                     ! (otherwise, relax to the
                                                     ! initial state) if
                                                     ! unifiedSponge == .true.
  relaxation_period        =                  0.0E+0 ! Period of an oscillation
                                                     ! superposed on the
                                                     ! background wind if
                                                     ! unifiedSponge == .true.
                                                     ! and relax_to_mean ==
                                                     ! .false. (0 for no
                                                     ! oscillation)
  relaxation_amplitude     =                  0.0E+0 ! Relative amplitude of an
                                                     ! oscillation superposed
                                                     ! on the background wind
                                                     ! if unifiedSponge ==
                                                     ! .true. and relax_to_mean
                                                     ! == .false.
  xBoundary                =              'periodic' ! Boundary conditions in x
                                                     ! ('periodic' only)
  yBoundary                =              'periodic' ! Boundary conditions in y
                                                     ! ('periodic' only)
  zBoundary                =            'solid_wall' ! Boundary conditions in z
                                                     ! ('solid_wall')
&end