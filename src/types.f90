module type_module

  !-----------------------------------------------------------------  
  !    Definition of data types and variables accessible
  !    throughout the code! Use care when working on these
  !    data fields.
  !-----------------------------------------------------------------  
 

  !-----------------------------------------------------------
  implicit none
  !-----------------------------------------------------------

  public 
  ! all variables declared herein are public by default


  ! modified by Junhong Wei (20161106) *** starting line ***

  !-----------------------------------------------------------------
  !                     MPI & Domain composition
  !-----------------------------------------------------------------
  
  integer :: sizeX, sizeY, sizeZ
  integer :: nbx, nby, nbz
  integer :: nprocx, nprocy
  real, dimension(0:1) :: lx_dim, ly_dim, lz_dim  ! dimensional domain
  
  namelist / domain / sizeX, sizeY, sizeZ, nbx, nby, nbz, &
       & lx_dim, ly_dim, lz_dim, nprocx, nprocy
  
  integer :: sizeXX, sizeYY, sizeZZ
  integer :: nx1, in1, ny1, jn1
  integer :: iStart, jStart
  integer :: is, ie, js, je      ! local start and end indices
  logical :: verboseMPI


  ! MPI include (parameters needed below)
  include 'mpif.h' 
  
  ! MPI variables
  integer :: ierror, count
  integer, dimension(2) :: dims, coords
  logical, dimension(2) :: period
  integer :: back, forw, right, left, rank
  integer :: idim, jdim, icoord, jcoord
  integer :: nbProc, comm, root             
  logical :: master
  integer, dimension(mpi_status_size) :: sts_left, sts_right, sts_back, sts_forw
  




  
  !-----------------------------------------------------------------  
  !                      (local) Grid & Domain
  !-----------------------------------------------------------------  
  
  integer  :: nx, ny, nz
  real, dimension(0:1) :: lx, ly, lz              ! scaled domain
  
  real, dimension(:), allocatable  :: x,y,z
  real :: dx, dy, dz
  integer :: nxx, nyy, nzz  ! grid size inclusive ghost cells
  integer :: nxyz           ! = nx*ny*nz
  
  !--------------------------------------------------------------
  !              for output global fields to single file
  !--------------------------------------------------------------
  
  ! real*4, dimension(:,:), allocatable :: field_out
  real*4, dimension(:,:), allocatable :: field_out, field_mst
  
  
!  !-----------------------------------------------------------------!  
!  !                         Grid & Domain
!  !-----------------------------------------------------------------!  
!  
!  integer  :: nx, ny, nz
!  integer :: nbx, nby, nbz
!  real, dimension(0:1) :: lx, ly, lz              ! scaled domain
!  real, dimension(0:1) :: lx_dim, ly_dim, lz_dim  ! dimensional domain
!
!  
!  namelist / grid / nx, ny, nz, nbx, nby, nbz, &
!       & lx_dim, ly_dim, lz_dim
!  
!  real, dimension(:), allocatable  :: x,y,z
!  real :: dx, dy, dz
!  integer :: nxx, nyy, nzz  ! grid size inclusive ghost cells
!  integer :: nxyz           ! = nx*ny*nz

  ! modified by Junhong Wei (20161106) *** fnishing line ***


  !-----------------------------------------------------------------  
  !                          Variables 
  !-----------------------------------------------------------------  
  integer :: nVar, nOptVar
  logical :: include_ice ! controls use of additional ice variables nIce,qIce and SIce
  namelist / variables / nVar, nOptVar, include_ice

 

  !-----------------------------------------------------------------  
  !                    Input / Output variables 
  !-----------------------------------------------------------------  
  
  character (len=40)    :: dataFileName
  character (len=40)    :: restartFile   ! Tecplot file with restart data
  logical, dimension(3) :: dimOut     ! (/1,0,1/) = 2D (x,z), (/1,1,1) = 3D
  integer, dimension(10) :: varOut     ! 1 = output, 0 = no output           

! achatzb
! data available in restart file
  integer, dimension(10) :: varIn      ! 1 = available, 0 = not available
! record to be read from restart file (starting from 0!)
  integer :: iIn
! achatze

  real, dimension(8)    :: offset     ! offset for rho, u,v,w, p, nIce, qIce, SIce
  integer, dimension(6) :: optVarOut  ! 1 = output, 0 = no output
  integer, dimension(13):: wkbVarOut  !           --"--
  ! increase dimension of optVarOut for new optional variables

  ! subtract background yes/no
  logical               :: thetaOffset
  logical               :: rhoOffset

  character(len=20)     :: outputType ! "time" or "timeStep"
  integer               :: nOutput    ! output every nOutput's time step 
  integer               :: maxIter    ! max nb. of time steps
  real                  :: outputTimeDiff ! output every ... seconds
  real                  :: maxTime    ! max time in seconds
  logical               :: solutionTime ! TECPLOT's "solution time"
  character(len=3)      :: solutionTimeUnit ! 
  logical :: restart             ! reads restartFile in TECPLOT format
  logical :: showGhostCellsX     ! plots include ghost cells in x-direction
  logical :: showGhostCellsY     ! y-direction
  logical :: showGhostCellsZ     ! z-direction

  namelist / outputList / &
       & dataFileName, restartFile, dimOut, varOut, &
       !achatzb
       & varIn, iIn, &
       !achatze
       & offset, optVarOut, wkbVarOut, outputType, nOutput, maxIter, &
       & outputTimeDiff, maxTime, solutionTime, solutionTimeUnit, &
       & restart, showGhostCellsX, showGhostCellsY, showGhostCellsZ, &
       & thetaOffset, rhoOffset
       

  !-----------------------------------------------------------------  
  !                         Parameter study 
  !-----------------------------------------------------------------  
  
   logical :: parameterStudy                   ! .true. / .false. 
  integer :: startParam, endParam, stepParam
  character(len=50) :: paramName
  namelist / parameterList / parameterStudy, startParam, endParam, &
       & stepParam, paramName


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
  character (len=50) :: testCase
  namelist / testCaseList / testCase

  !------------------------
  !  monochromatic wave
  !------------------------

  real :: lambda_dim       ! wave length in (inclined) z direction
  
  namelist / monochromeWave / lambda_dim, amplitudeFactor


! modified by Junhong Wei for 3DWP (20170828) *** starting line ***

  !---------------------------------------
  !         Gravity wave packet
  !---------------------------------------
  ! test cases: 
  ! 1) wavePacket1D and wavePacket1D_raytracer
  ! 2) wavePacket2D and wavePacket2D_raytracer
  integer :: wavePacketDim
  integer :: wavePacketType
  real :: lambdaX_dim, lambdaY_dim, lambdaZ_dim
  real :: meanFlowX_dim, meanFlowZ_dim
  real :: amplitudeFactor
  real :: sigma_dim
  real :: sigma_hor_dim, sigma_hor_yyy_dim   ! modified by Junhong Wei (20170214)
  real :: L_cos_dim
  integer :: omiSign
  real :: u0_jet_dim
  real :: z0_jet_dim
  real :: L_jet_dim
  real :: xCenter_dim, yCenter_dim, zCenter_dim
! achatzb
  real :: amp_mod_x, amp_mod_y
! achatze
  
  namelist / wavePacket / &
       &       wavePacketType, &
       &       wavePacketDim, &
       &       lambdaX_dim, lambdaY_dim, lambdaZ_dim, &
       &       meanFlowX_dim, meanFlowZ_dim, &
       &       amplitudeFactor, &
       &       xCenter_dim, yCenter_dim, zCenter_dim, sigma_dim, &
       &       sigma_hor_dim, &
! achatzb
       &       amp_mod_x, &
! achatze
       &       sigma_hor_yyy_dim, &
! achatzb
       &       amp_mod_y, &
! achatze
       &       L_cos_dim, omiSign, &
       &       u0_jet_dim, z0_jet_dim, L_jet_dim

! modified by Junhong Wei for 3DWP (20170828) *** finishing line ***
  
  ! hotBubble, coldBubble, hotBubble3D
!  real :: dTheta0_dim, xRadius_dim, zRadius_dim, xCenter_dim, zCenter_dim   ! modified by Junhong Wei for 3DWP (20170922)
  real :: dTheta0_dim, xRadius_dim, zRadius_dim   ! modified by Junhong Wei for 3DWP (20170922)
  real :: zExcentricity
  namelist / bubble / &
       &       dTheta0_dim, xRadius_dim, zRadius_dim, &
       &       xCenter_dim, zCenter_dim, zExcentricity
    

  ! hot and cold bubble by Robert
  real :: dTheta1_dim, a1_dim, sigma1_dim, xCenter1_dim, zCenter1_dim
  real :: dTheta2_dim, a2_dim, sigma2_dim, xCenter2_dim, zCenter2_dim
  namelist / robert_bubble / &
       &    dTheta1_dim, a1_dim, sigma1_dim, xCenter1_dim, zCenter1_dim, &
       &    dTheta2_dim, a2_dim, sigma2_dim, xCenter2_dim, zCenter2_dim


! achatzb
  !-----------------------------------------------------------------
  !                          Mountain Wave
  !-----------------------------------------------------------------

  ! zonal wind to be attained by temporary wind relexation
  real :: u_relax
  ! total relaxation time
  real :: t_relax
  ! duration of ramping up/down the relaxation
  real :: t_ramp
  ! zonal extent of region without wind relaxation
  real :: xextent_norelax
  namelist / mountainwavelist / u_relax,t_relax,t_ramp,xextent_norelax
! achatze
  !-----------------------------------------------------------------  
  !                          Model equations
  !-----------------------------------------------------------------  
  
  ! vertical direction
  real, dimension(3) :: vertical
  real :: vert_theta, vert_alpha
  
  character(len=25) :: model

  namelist / modelList / model, vert_theta, vert_alpha
  

  !-----------------------------------------------------------------  
  !                               Solver 
  !-----------------------------------------------------------------  

  real :: cfl
  real :: cfl_wave
  real :: dtMax_dim 
  character(len=20) :: tStepChoice         ! "cfl", "fix"
  character(len=20) :: timeScheme          ! LS_Will_RK3 / Euler
  character(len=20) :: timeSchemeType      ! lowStorage / classical
  character(len=20) :: fluxType            ! ILES / central / upwind
  character(len=20) :: reconstType         ! ALDM / constant / SALD / MUSCL
  character(len=20) :: limiterType1        ! minmod / ...
  logical           :: fluctuationMode     ! split rho = rhoStrat + rho'
  logical           :: DySmaScheme         ! Dynamic Smagorinsky Scheme ! modified by Junhong Wei (20160722)
  namelist / solverList / cfl, cfl_wave, dtMax_dim, &
       & tStepChoice, timeScheme, &
       & fluxType, reconstType, limiterType1, &
!       & fluctuationMode                         ! modified by Junhong Wei (20160722)
       & fluctuationMode, DySmaScheme             ! modified by Junhong Wei (20160722)

  integer :: nStages
  logical :: updateMass         ! transport of mass=var(1)  on/off
  logical :: predictMomentum    ! transport of momentum=var(2-4) on/off
  logical :: updateTheta        ! transport of theta=var(6) on/off
  logical :: updateIce          ! transport of ice=var(8-10) on/off

 
  !-----------------------------------------------------------------  
  !                          Poisson solver 
  !-----------------------------------------------------------------

  real :: tolPoisson
! achatzb
  real :: tolCond,tolref
! achatze
  integer :: maxIterPoisson
  character(len=20) :: poissonSolverType
  character(len=20) :: storageType
  character(len=10) :: preconditioner
  real :: dtau                             
  integer :: maxIterADI
  logical :: initialCleaning
  logical :: pressureScaling
  logical :: useNAG
  logical :: correctMomentum       ! false -> momentumCorrector off
  logical :: correctDivError       ! true -> subtract rho*div(u)
  namelist / poissonSolverList / tolPoisson, tolCond, maxIterPoisson, &
       & poissonSolverType, storageType, preconditioner, dtau, maxIterADI, &
       & initialCleaning, pressureScaling, useNAG, &
       & correctMomentum, correctDivError
  integer :: nnz      ! number of nonzeros 


  !-----------------------------------------------------------------
  !                           Constants 
  !-----------------------------------------------------------------
 
  real               :: pi                 ! you know...
  real, parameter    :: small = 1.0e-20    ! to devision by zero
  complex, parameter :: imag = (0.0,1.0)   ! imaginary unit



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
  
  character(len=10) :: referenceQuantities     ! set of reference quantities
  
  logical :: specifyReynolds           ! choose Reynolds or mu_viscous
  real :: ReInv                        ! reciprocal Reynolds number
  real :: mu_viscous_dim               ! kinematic viscosity
  real :: mu_conduct_dim               ! heat conductivity
  character(len=30) :: background      ! isentropic / isothermal / const-N
  real :: N_BruntVaisala_dim           ! Brunt-Vaisala frequency in 1/s
  real :: theta0_dim                   ! isentr. background pot. temp. in K
  real :: Temp0_dim                    ! isothermal background temperature in K
  real :: press0_dim                   ! pressure at z=0 in Pa
  real,dimension(3) :: backgroundFlow_dim
  real :: f_Coriolis_dim               ! Coriolis parameter
  real :: z_tr_dim                     ! height of topopause
  real :: theta_tr_dim                 ! const pot. temp. of troposphere
  namelist / atmosphereList / referenceQuantities, &
       & specifyReynolds, ReInv, mu_viscous_dim, mu_conduct_dim,&
       & background, &
       & N_BruntVaisala_dim, theta0_dim, &
       & Temp0_dim, press0_dim, &
       & backgroundFlow_dim,&
       & f_Coriolis_dim,&
       & z_tr_dim, theta_tr_dim

  real,dimension(3) :: backgroundFlow
  real :: theta00, rho00, P00          ! background values for Boussinesq

  real :: mu_conduct


  !-----------------------------------------------------------------
  !                         Bottom topography 
  !-----------------------------------------------------------------
  
  logical :: topography    ! via k = 1
! achatzb
! topography_mask = .true.  if cell is below topographic surface
! topography_mask = .false. if cell is above topographic surface
! topography_surface x-y-dependent mountain surface
  logical, dimension(:,:,:), allocatable :: topography_mask
  real, dimension(:,:), allocatable :: topography_surface
! achatze
  real :: mountainHeight_dim
  real :: mountainWidth_dim
  namelist / topographyList / topography, &
       & mountainHeight_dim, mountainWidth_dim
  
  
  !-----------------------------------------------------------------
  !                         Boundary
  !-----------------------------------------------------------------

  logical :: rhoFluxCorr, uFluxCorr, vFluxCorr, wFluxCorr, thetaFluxCorr
  integer :: nbCellCorr
  
  ! sponge layer
  logical :: spongeLayer
  real    :: spongeHeight
  integer :: kSponge
  real    :: zSponge
  real    :: spongeAlphaZ_dim

  ! boundary types
!  real :: xBoundary
!  character(len=15) :: xxBoundary
!  character(len=15) :: yBoundary
!  character(len=15) :: zBoundary

  namelist / boundaryList / rhoFluxCorr, uFluxCorr, &
       & vFluxCorr, wFluxCorr, thetaFluxCorr, nbCellCorr, &
       & spongeLayer, spongeHeight, &
       & zSponge, &
       & spongeAlphaZ_dim

  ! boundary types
  character(len=15) :: xBoundary
  character(len=15) :: yBoundary
  character(len=15) :: zBoundary

  namelist / boundaryList2 / xBoundary, yBoundary, zBoundary

  

  !-----------------------------------------------------------------
  !                             WKB  
  !-----------------------------------------------------------------
  
  type rayType
     real :: x             ! ray position
     real :: y 
     real :: z
     real :: k             ! ray wave vector
     real :: l
     real :: m
     real :: cgx           ! intrinsic group velocity
     real :: cgy
     real :: cgz
     real :: omega         ! intrinsic frequency
     integer :: iCell      ! position within finite-volume grid
     integer :: jCell
     integer :: kCell
     logical :: right      ! position within grid cell: true = right half
     logical :: up         !                            true = upper half
  end type rayType
   

  ! total number of rays
  integer :: nRay
  
  
  ! namelist "rayTracer"
  logical :: rayTracer ! run ray tracer
  logical :: outComplex     ! output complex amplitudes
  integer :: nRayRatioX, nRayRatioY, nRayRatioZ
  character(len=20) :: waveFluxType
  character(len=20) :: limiterType
  
  namelist / wkbList / rayTracer, outComplex, nRayRatioX, nRayRatioY, nRayRatioZ, &
       &     waveFluxType, limiterType
  

end module type_module
