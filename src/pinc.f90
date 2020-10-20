program pinc_prog

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !% calculates pseudo-incompressible model     %
  !% by Durran with ILES by Adams and Hickel    %
  !% cylflow (c) 2009 Felix Rieper              %
  !% latest revision: January 2010              %
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
  use type_module
  use mpi_module
  use timeScheme_module
  use init_module
  use debug_module
  use wkb_module    
  use output_module
  use xweno_module
  use atmosphere_module
  use boundary_module
  use flux_module
  use update_module
  use poisson_module
  use finish_module
  use ice_module
  use hypretools_module
 
  ! test
  use algebra_module

  !----------------------------------------------------
  implicit none
  !---------------------------------------------------

  integer                     :: iTime, Ice_RKstage
  integer                     :: ice_time_steps
  real                        :: time, dt

  ! CPU Time
  integer                     :: rate, startTimeCount, timeCount  
  real                        :: cpuTime 

  ! MPI stuff            
  logical :: error_flag  
  real :: dt_local     

  ! fields
  real, dimension(:,:,:,:), allocatable :: var, var0, var1
  real, dimension(:,:,:,:), allocatable :: source
  ! var(i,j,k,iVar) iVar = 1..5 > rho,u,v,w,pExner

  real, dimension(:,:,:), allocatable :: dRho, dRhop  ! RK-Update for rho
  real, dimension(:,:,:,:), allocatable :: dMom    ! RK for rhoU,rhoV,rhoW
  real, dimension(:,:,:), allocatable :: dTheta     ! RK-Update for theta
  real, dimension(:,:,:,:), allocatable :: dIce     ! RK-Update for nAer,nIce,qIce,qv

  !UAB
  real, dimension(:), allocatable :: dPStrat, drhoStrat !RK update for P  
  real, dimension(:), allocatable :: w_0
  !UAE


  real, dimension(:,:,:,:,:), allocatable :: flux
  ! flux(i,j,k,dir,iFlux) 
  ! dir = 1..3 > f,g,h-flux in x,y,z-direction
  ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW

 
  !--------------------
  !   WKB variables
  !--------------------
  type(rayType), dimension(:,:,:,:), allocatable     :: ray
  real, dimension(:,:,:,:), allocatable :: ray_var3D


  ! topography via force field
  real, dimension(:,:,:,:), allocatable :: force ! volume forces
  ! force(i,j,k,forceVector)

  ! output per timeStep
  logical :: output
  real    :: nextOutputTime   ! scaled time for next output

  ! general
  integer :: i,j,k,l
  integer :: ix, jy, kz

  ! restart
  logical :: scale

  ! parameter study
  integer :: iParam

  ! classical RK
  real    :: dt_Poisson

  ! error handling and controlling
  logical :: errFlagBiCG, errFlagTStep
  integer :: nIterBicg, nTotalBicg
  real    :: nAverageBicg
  character (len=8)  :: datum
  character (len=10) :: niceDatum
  character (len=10) :: zeit
  character (len=9)  :: niceZeit
  real :: a2

  real :: tolpoisson_s
  real :: fc_shap, shap_dts  
  
  integer :: iVar

  integer :: heating_switch

  !-------------------------------------------------
  !                    Set up 
  !-------------------------------------------------

  file_namelist = 'input.f90'   ! default
  if ( command_argument_count() /= 0 )  &
    &  call get_command_argument(1,file_namelist)  ! take the given one
    
  ! init counter and time
  iOut = 0
  iTime = 0; time = 0.0; cpuTime=0.0

  call init_mpi(error_flag)

  if( error_flag )then
     print*,"error in the init_mpi"
     goto 666
  end if

  if( master ) then
     call date_and_time(date=datum, time=zeit)

     niceDatum = datum(7:8)//"."//datum(5:6)//"."//datum(1:4)
     niceZeit = zeit(1:2)//":"//zeit(3:4)//" Uhr"

     print*,""
     print*,""
     print*,""
     print*,"=============================================================="
     print*," PincFloit (c) 2010 F. Rieper  "
     print*," modified by many others (2018)  "
     print*,""
     write(*,fmt="(a25,a10)") " Today : ", niceDatum
     write(*,fmt="(a25,a9)") " Time  : ", niceZeit

     ! init clock
     call system_clock (count_rate=rate)
     call system_clock (count=startTimeCount)
  end if

  !-------------------------------------------------
  !                    Set up 
  !-------------------------------------------------

  ! init counter and time
  nTotalBicg = 0

  ! 1) allocate variables 
  ! 2) read input.f90
  !UAB
  !call setup (var,var0,var1,flux,force,source,dRho,dRhop,dMom,dTheta)
    call setup(var,var0,var1,flux,force,source,dRho,dRhop,dMom,dTheta, &
            & dPStrat,drhoStrat,w_0,dIce)
  !UAE

  call init_atmosphere      ! set atmospheric background state 
  call init_output

  call initialise (var)     ! set initial conditions

  var(:,:,:,8) = 0.0   ! Heating due to GWs in the rotating atmosphere

  if (poissonSolverType == 'hypre') then
     call SetUpHypre       ! Set Up Hypre objects
    else if (poissonSolverType == 'bicgstab') then
     call SetUpBiCGStab       ! Set BiCGStab arrays
    else
     stop 'ERROR: only HYPRE and BiCGStab ready to be used'
  end if

  call init_xweno       ! set ILES parameters 
  call init_fluxes      ! allocate tilde variables
  call init_update
  call init_timeScheme  ! define Runge-Kutta parameters
  call init_poisson     ! allocate dp
  
  ice_time_steps = 1 
     

  !-------------------------------------------------
  !              Read initial data
  !-------------------------------------------------
  if (restart) then
     !restart from previous output

     call read_data ( iIn,var)

     call setHalos( var, "var" )
     call setBoundary (var, flux, "var")
  end if

  !testb
  !call output_data(iOut, var, iTime, time, cpuTime)
  !teste

  !---------------------------------------------
  !        Initial divergence cleaning
  !---------------------------------------------

  if ( initialCleaning ) then
     call setHalos( var, "var" )
     call setBoundary (var,flux,"var")

     ! 1) allocate variables 
     ! 2) read the namelist
     !call setup(var,var0,var1,flux,force,source,dRho,dRhop,dMom,dTheta, &
     !       & dPStrat,drhoStrat,w_0,dIce)
     
     flux = 0.
     dMom = 0.

     tolpoisson_s = tolPoisson
     if (testCase == 'baroclinic_LC') then
        tolPoisson = 1.e-6
     else
        tolPoisson = 1.e-8
     end if
     call Corrector ( var, flux, dMom, 1.0 , errFlagBicg, nIterBicg, &
                    & 1, "expl")
     tolPoisson = tolpoisson_s
  end if

  !---------------------------------------------
  !               Init ray tracer
  !---------------------------------------------

  if (rayTracer) then
     ! allocate ray fields

     call setup_wkb(ray, ray_var3D, var) 
  end if

  !UAB 200413
  !-------------------------------------------------------------------
  ! in diabatic case store initial reference atmosphere for the sponge
  !-------------------------------------------------------------------

  if (heatingONK14 .or. TurbScheme .or. rayTracer) then
     rhoStrat_0 = rhoStrat
     pStrat_0 = pStrat
  end if
  !UAE 200413

  !------------------------------------------
  !              Initial output
  !------------------------------------------

  call output_data(iOut, var, iTime, time, cpuTime)
  if (testCase == 'baroclinic_LC') then
     call output_profile(iOut, PStrat,'Pstrat.dat')
     call output_profile(iOut, rhoStrat,'rhostrat.dat')
     call output_profile(iOut, thetaStrat,'thetastrat.dat')
     call output_profile(iOut, bvsstrat,'bvsstrat.dat')
  end if
     

  if (rayTracer) then
      call output_wkb(iOut, ray, ray_var3D)
  end if

  output = .false.
  ! set time for first output
  nextOutputTime = time*tRef + outputTimeDiff
  ! and consider restart time


  !-----------------------------------------------------
  !                        Time loop
  !-----------------------------------------------------

  if( master ) then
     print *, "...starting time loop"
  end if 

  if( outputType == "time" ) maxIter = 2**30

 
  time_loop: do iTime = 1,maxIter

     if( master ) then
        print*,""
        print*,""
        print*,""
        print*,"--------------------------------------------------------"
        write(*,fmt="(a25,i15)") " Time step = ", iTime
        write(*,fmt="(a25,f15.1,a8)") " Time = ", time*tRef, "seconds"
        print*,"--------------------------------------------------------"
     end if 


     !----------------------------------
     !         Calc time step
     !----------------------------------

     call timestep (var,dt, errFlagTstep)

     !! in order to obtain an equilibrated pressure do first an 
     !! extremely short time step

     !! testb
     !!goto 20
     !! teste

     !if (iTime == 1) dt = max(1.e-6*dt,dtMin_dim/tRef)

  20 continue

     !testb
     !if (master) print*,'dt * tRef =',dt * tRef
     !teste

     ! abort if time step too small
     if( dt*tRef < dtMin_dim ) then
        if( master ) then
           print*," TimeStep routine: "
           write(*,fmt="(a25,es15.4,a8)") &
                &"dt < dtMin!!! dtMin = ", dtMin_dim, "seconds"
           write(*,fmt="(a25,es15.4,a8)") "dt * tRef = ", dt*tRef, &
                   & "seconds"
           write(*,fmt="(a25,i2.2)") "cycle paramLoop. iPara = ", iParam

           call system_clock(count=timeCount)
           cpuTime = (timeCount - startTimeCount)/real(rate)
        end if

        ! final output
        call output_data(iOut, var, iTime, time, cpuTime)
        if (testCase == 'baroclinic_LC') then
           call output_profile(iOut, PStrat,'Pstrat.dat')
           call output_profile(iOut, rhoStrat,'rhostrat.dat')
           call output_profile(iOut, thetaStrat,'thetastrat.dat')
           call output_profile(iOut, bvsstrat,'bvsstrat.dat')
        end if
        
        if (rayTracer) then
            call output_wkb(iOut, ray, ray_var3D)
        end if

        call mpi_barrier(comm, ierror)
        call mpi_finalize(ierror)
        print*,"pinc.f90: time step too small. dt = ", dt*tRef
        stop
     end if

     if( timeSchemeType == "classical" .and. iTime == 1 ) dt = 0.5*dt

     ! correct dt to hit desired output time for outputType 'time'
     if ( outputType == 'time' ) then
        !testb
        !if (master) then
        !   print*,'(time+dt) * tRef =',(time+dt) * tRef
        !   print*,'(time+dt) * tRef + dtMin_dim =', &
        !         & (time+dt) * tRef + dtMin_dim
        !   print*,'nextOutputTime =',nextOutputTime
        !end if
        !teste
        if( (time+dt) * tRef + dtMin_dim > nextOutputTime ) then
           dt = nextOutputTime/tRef - time      
           output = .true.
           if(master) then
              write(*,fmt="(a25,es15.4,a8)") "dt for output = ", &
                   & dt*tRef, "seconds"
           end if
        end if
     end if

     time = time + dt

     !-----------------------------------------------------------------
     ! relaxation rate for 
     ! (1) Rayleigh damping in land cells and 
     ! (2) density-fluctuation relaxation in semi-implicit time stepping
     !-----------------------------------------------------------------

     if (testCase == 'baroclinic_LC') then
        alprlx = 0.
     else
        alprlx = 5.e-1/dt
     end if

     !testb
     !alprlx = 0.0
     !teste

     !------------------------------------
     !    Call wave saturation scheme
     !------------------------------------
        
     if (rayTracer) then
        if( lsaturation ) then
            call saturation_3D(ray,dt)
        end if
     end if

     !---------------------------------------------------------------
     !          Runge-Kutta stages or semi-implicit time step
     !---------------------------------------------------------------

     if (timeScheme == "semiimplicit") then
        ! testb
        !call test_hypre
        !stop
        ! teste

        ! just for safety ...

        if( correctDivError ) then
            print*,'ERROR: correction divergence error not &
                  & implemented properly'
            stop
        end if

        if( updateTheta ) then
           print*,'ERROR: semiimplicit time stepping does not allow &
                  & updateTheta = .true.'
           stop
        end if

        ! initialize zero volume force 
        ! (to be filled by ray tracer and wind relaxation at the 
        ! horizontal boundaries)

        force = 0.0

        ! Lagrangian WKB model (position-wavenumber space method)

        if (rayTracer) then
           do RKstage = 1, nStages
              call transport_rayvol(var, ray, dt, RKstage)  

              if (RKstage == nStages) then
                 call boundary_rayvol(ray)
                 call split_rayvol(ray)
                 call shift_rayvol(ray)
                 call merge_rayvol(ray)

                 ! GW effects are put into force(...,1/2) and var(...,8)
                 call calc_meanFlow_effect(ray,var,force,ray_var3D)
              end if
           end do
        end if

        ! wind relaxation at horizontal boundaries

        if(       (testCase == "mountainwave") &
         & .or. (raytracer .and. case_wkb == 3)) then
           call volumeForce (var,time,force)        
        end if

        ! set density fluctuation (synchronization step)

        ! testb
        ! if (iTime > 1) goto 50
        ! teste

        if( fluctuationMode ) then
           do kz = -nbz,nz+nbz
              var(:,:,kz,6) = var(:,:,kz,1)
           end do
          else
           do kz = -nbz,nz+nbz
              var(:,:,kz,6) = var(:,:,kz,1) - rhoStrat(kz)
           end do
        end if

        !testb
        !call output_data(iOut, var, iTime, time, cpuTime)
        !teste
           
        ! turbulence scheme:
        ! either prescribed damping time scale for smallest spatial scales
        ! or dynamic Smagorinsky scheme
        ! diffusion coefficient (normalized by squared grid length scale)
        ! stored in var (...,7)

        if (TurbScheme) then
           if(DySmaScheme) then
              call CoefDySma_update(var)

              ! limit Smagorinsky coefficient so that the damping time 
              ! scale for the 2dx-wave is shorter than a time step
              var(:,:,:,7) = min(var(:,:,:,7), 1.e0/(dt * pi**2))
              !var(:,:,:,7) = min(var(:,:,:,7), 5.e0/(dt * pi**2))
             else
              var(:,:,:,7) = tRef/turb_dts
              !var(:,:,:,7) = 0.1* tRef/dt
           end if
        end if

        ! put initial state into var0 in order to save the advecting 
        ! velocities

        var0 = var
        !FS for update of P and rhoStrat in (5)
        PStrat00 = PStrat 
        rhoStrat00 = rhoStrat
        heating_switch = 0

        !FS18082020
        thetaStrat00 = thetaStrat
        bvsStrat00 = bvsStrat
        thetaStratTilde00 = thetaStratTilde
        rhoStratTilde00 = rhoStratTilde
        PStratTilde00 = PStratTilde
        
        call setHalos( var0, "var" )
        call setBoundary (var0,flux,"var")

        ! (1) explicit integration of convective and 
        !     viscous-diffusive/turbulent fluxes over half a time step,
        !     with the advetion velocity kept constant
        !     \psi^# = \psi^n + A^{dt/2} (\psi^n, v^n)

        !testb
        if (master) print*,'beginning a semi-implicit time step'
        if (master) print*,'(1) explicit integration lhs over dt/2'
        !teste


        do RKstage = 1, nStages
           ! Reconstruction

           call setHalos( var, "var" )
           call setBoundary (var, flux, "var")

           call reconstruction (var, "rho") 
           call reconstruction (var, "rhop")         
           call reconstruction (var, "uvw")
           
           call setHalos( var, "varTilde" )
           call setBoundary (var, flux, "varTilde" ) 

           ! Fluxes and Forces
           
           call massFlux (var0,var,flux,"lin",PStrat00,PStratTilde00)
           call momentumFlux (var0,var,flux,"lin",PStrat00,PStratTilde00) 

           call setBoundary (var, flux, "flux") 

           ! RK step for density and density fluctuations
           
           rhoOld = var(:,:,:,1)  ! rhoOld for momentum predictor

           
           call massUpdate(var, flux, 0.5*dt, dRho, RKstage, &
                         & "rho", "tot", "expl")
         
           call massUpdate(var,flux, 0.5*dt, dRhop, RKstage, &
                         & "rhop", "lhs", "expl")
          
              
           ! RK step for momentum
        
           call setHalos( var, "var" )
           call setBoundary (var, flux, "var")

           
           call momentumPredictor(var, flux, force, 0.5*dt, dMom, &
                                & RKstage, "lhs", "expl")
        
           
        end do



        ! (2) implicit integration of the linear right-hand sides of the
        !     equations for density fluctuations and momentum over half a 
        !     time step, under consideration of the divergence constraint 
        !     \psi^{n+1/2} = \psi^# + dt/2 Q(\psi^{n+1/2})

        !testb
        if (master) print*,'(2) implicit integration rhs over dt/2'
        !teste

        call setHalos( var, "var" )
        call setBoundary (var, flux, "var")


        rhopOld = var(:,:,:,6)  ! rhopOld for momentum predictor

        ! update density fluctuations (rhopStar)

        call massUpdate(var,flux,0.5*dt,dRhop,RKstage,"rhop","rhs","impl")

        ! update winds (uStar, vStar, wStar)

        call momentumPredictor(var,flux,force,0.5*dt,dMom,RKstage,"rhs",&
                             & "impl")
        
              
        ! corrector: rhopStar, uStar, vStar, wStar 
        !            -> new rhop, u, v, w
           
        call setHalos( var, "var" )
        call setBoundary (var,flux,"var")



        !UAB
        if (shap_dts_dim > 0.) then
           ! smoothing of the fields in order to limit grid-point noise


           shap_dts = shap_dts_dim/tRef

           fc_shap = min( 1.0, 0.5*dt/shap_dts)
   
           call smooth_hor_shapiro(fc_shap,n_shap,flux,var,0.5*dt)
           !call smooth_shapiro(fc_shap,n_shap,flux,var)
        end if
        !UAE
       


        ! GBcorr -> FS
        if (heatingONK14 .or. TurbScheme .or. rayTracer) then
        !if (heating) then
           if (model == 'Boussinesq') then
              print*, "main:ONeill+Klein2014 heating only for &
                     & pseudo-incompressible dyn."
              stop
           end if

           RKstage = 1
           dPStrat = 0.
           drhoStrat = 0.           
           call BGstate_update(var,flux,0.5*dt,RKstage,dPStrat,drhoStrat, &
                             & "impl",heating_switch)
        end if

        call Corrector ( var, flux, dMom, 0.5*dt, errFlagBicg, nIterBicg, &
                       & RKstage, "impl")
        

        nTotalBicg = nTotalBicg + nIterBicg

        call setHalos( var, "var" ) 
        call setBoundary (var,flux,"var")

        ! put new state into var1 in order to save the advection velocities

        var1 = var     
        !FS18082020
        PStrat01 = PStrat
        rhoStrat01 = rhoStrat
        thetaStrat01 = thetaStrat
        bvsStrat01 = bvsStrat
        thetaStratTilde01 = thetaStratTilde
        rhoStratTilde01 = rhoStratTilde
        PStratTilde01 = PStratTilde



        ! (3) explicit integration of the linear right-hand sides of the
        !     equations for density fluctuations and momentum over half a 
        !     time step, under consideration of the divergence constraint 
        !     \psi^\ast = \psi^n + dt/2 Q(\psi^n)
        !     
        !     could also be replaced by an implicit time step

        !testb
        if (master) print*,'(3) explicit integration rhs over dt/2'
        !teste

        var = var0
        !FS18082020
        PStrat = PStrat00
        rhoStrat = rhoStrat00
        thetaStrat = thetaStrat00
        bvsStrat = bvsStrat00
        thetaStratTilde = thetaStratTilde00
        rhoStratTilde = rhoStratTilde00
        PStratTilde = PStratTilde00
      
        call setHalos( var, "var" )
        call setBoundary (var, flux, "var")
      

        rhopOld = var(:,:,:,6)  ! rhopOld for momentum predictor

        ! update density fluctuations (rhopStar)
        call massUpdate(var,flux,0.5*dt,dRhop,RKstage,"rhop","rhs","expl")
        

        ! update winds (uStar, vStar, wStar)

        call momentumPredictor(var,flux,force,0.5*dt,dMom,RKstage,"rhs",&
                             & "expl")

        ! corrector: uStar, vStar, wStar 
        !            -> new u, v, w

           
        call setHalos( var, "var" )
        call setBoundary (var,flux,"var")

        if (shap_dts_dim > 0.) then
           ! smoothing of the fields in order to limit grid-point noise


           shap_dts = shap_dts_dim/tRef

           fc_shap = min( 1.0, 0.5*dt/shap_dts)

           call smooth_hor_shapiro(fc_shap,n_shap,flux,var,0.5*dt)
           !call smooth_shapiro(fc_shap,n_shap,flux,var)
        end if
       

        ! GBcorr -> FS
        !FS auskommentiert Anfang
        ! if (heatingONK14 .or. TurbScheme .or. rayTracer) then
        ! !if (heating) then
        !    if (model == 'Boussinesq') then
        !       print*, "main:ONeill+Klein2014 heating only for &
        !              & pseudo-incompressible dyn."
        !       stop
        !    end if

        !    !'impl' chosen below because this is not in an RK sub-step

        !    RKstage = 1
        !    dPStrat = 0.
        !    drhoStrat = 0.           
        !    call BGstate_update(var,flux,0.5*dt,RKstage,dPStrat,drhoStrat, &
        !                      & "impl")
        ! end if
         
        ! call Corrector ( var, flux, dMom, 0.5*dt, errFlagBicg, nIterBicg, &
        !      & RKstage, "expl")
        
        ! ! !FS auskommentiert Ende
        !  nTotalBicg = nTotalBicg + nIterBicg



        ! (4) explicit integration of convective and 
        !     viscous-diffusive/turbulent fluxes over a full time step,
        !     with the advection velocity kept constant
        !     \psi^{\ast\ast} = \psi^\ast + A^dt (\psi^\ast, v^{n+1/2})

        !testb
        if (master) print*,'(4) explicit integration lhs over dt'
        !teste

        var0 = var1
        !FS
        PStrat = PStrat01
        rhoStrat = rhoStrat01
        thetaStrat = thetaStrat01
        bvsStrat = bvsStrat01
        thetaStratTilde = thetaStratTilde01
        rhoStratTilde = rhoStratTilde01
        PStratTilde = PStratTilde01

        call setHalos( var0, "var" )
        call setBoundary (var0,flux,"var")

        do RKstage = 1, nStages
           ! Reconstruction

           call setHalos( var, "var" )
           call setBoundary (var, flux, "var")

           call reconstruction (var, "rho")  
           call reconstruction (var, "rhop")         
           call reconstruction (var, "uvw")
           
           call setHalos( var, "varTilde" )
           call setBoundary (var, flux, "varTilde" ) 

           ! Fluxes and Forces
           
           call massFlux (var0,var,flux,"lin",PStrat01,PStratTilde01)
           call momentumFlux (var0,var,flux,"lin",PStrat01,PStratTilde01) 

           call setBoundary (var, flux, "flux") 


           ! RK step for density and density fluctuations
           
           rhoOld = var(:,:,:,1)  ! rhoOld for momentum predictor

         
           call massUpdate(var, flux, dt, dRho, RKstage, &
                        & "rho", "tot", "expl")

           call massUpdate(var, flux, dt, dRhop, RKstage, &
                        & "rhop", "lhs", "expl")
           

           ! RK step for momentum
           call setHalos( var, "var" )
           call setBoundary (var, flux, "var")

           
           call momentumPredictor(var, flux, &
                &                 force, dt, dMom, RKstage, "lhs", "expl")
  
        end do


        ! (5) implicit integration of the linear right-hand sides of the
        !     equations for density fluctuations and momentum over half a 
        !     time step, under consideration of the divergence constraint 
        !     \psi^{n+1} = \psi^{\ast\ast} + dt/2 Q(\psi^{n+1})

        !testb
        if (master) print*,'(5) implicit integration rhs over dt/2'
        !teste

        call setHalos( var, "var" )
        call setBoundary (var, flux, "var")

      
        rhopOld = var(:,:,:,6)  ! rhopOld for momentum predictor

        ! update density fluctuations (rhopStar)

       
        call massUpdate(var,flux,0.5*dt,dRhop,RKstage,"rhop","rhs","impl")
        

        ! update winds (uStar, vStar, wStar)
        call setHalos( var, "var" )
        call setBoundary (var,flux,"var")

       
        call momentumPredictor(var,flux,force,0.5*dt,dMom,RKstage,"rhs",&
                             & "impl")
              
        ! corrector: rhopStar, uStar, vStar, wStar 
        !            -> new rhop, u, v, w
           
        call setHalos( var, "var" )
        call setBoundary (var,flux,"var")

        !UAB
        if (shap_dts_dim > 0.) then
           ! smoothing of the fields in order to limit grid-point noise


           shap_dts = shap_dts_dim/tRef

           fc_shap = min( 1.0, 0.5*dt/shap_dts) 

           call smooth_hor_shapiro(fc_shap,n_shap,flux,var,0.5*dt)
           !call smooth_shapiro(fc_shap,n_shap,flux,var)
        end if
        !UAE

        ! GBcorr -> FS
        if (heatingONK14 .or. TurbScheme .or. rayTracer) then
        !if (heating) then
           if (model == 'Boussinesq') then
              print*, "main:ONeill+Klein2014 heating only for &
                     & pseudo-incompressible dyn."
              stop
           end if

           RKstage = 1
           dPStrat = 0.
           drhoStrat = 0. 
           heating_switch = 1
           call BGstate_update(var,flux,dt,RKstage,dPStrat,drhoStrat, &
                             & "impl",heating_switch) !FS 0.5*dt -> dt
        end if
!FSE

        call Corrector ( var, flux, dMom, 0.5*dt, errFlagBicg, nIterBicg, &
                       & RKstage, "impl")
        


        nTotalBicg = nTotalBicg + nIterBicg


        !testb
        if (master) print*,'semi-implicit time step done'
        !teste

     else ! (timeScheme /= "semiimplicit") explicit time stepping

        Runge_Kutta_Loop: do RKstage = 1, nStages

            ! turbulence scheme:
           ! either prescribed damping time scale for smallest spatial 
           ! scales
           ! or dynamic Smagorinsky scheme
           ! diffusion coefficient (normali. by squared grid length scale)
           ! stored in var (...,7)

           if (TurbScheme) then
              if(DySmaScheme) then
                 call CoefDySma_update(var)
                 !call CoefDySma_update(var,dt)
                else
                 var(:,:,:,7) = tRef/turb_dts
              end if
           end if

           if( timeSchemeType == "classical" ) then
              if( RKstage == 1 ) then
                 var0 = var
              end if
           end if

           ! initialize density fluctuations for the integration

           !if (auxil_equ .and. iTime == 1) then
           !if (auxil_equ) then
           if (     auxil_equ &
              &.or. heatingONK14 .or. TurbScheme .or. rayTracer) then
           !UAE 200413
              if (model /= "pseudo_incompressible") then
                 print*,'auxiliary equation only ready for &
                       & pseudo-incompressible'
                 stop
              end if

              alprlx = 0.

              !testb
              !if (iTime > 1) goto 450
              !teste

              if( fluctuationMode ) then
                 do kz = -nbz,nz+nbz
                    var(:,:,kz,6) = var(:,:,kz,1)
                 end do
                else
                 do kz = -nbz,nz+nbz
                    var(:,:,kz,6) = var(:,:,kz,1) - rhoStrat(kz)
                 end do
              end if

              var0 = var !UA 200413

       450    continue
           end if

           ! initialize zero volume force

           force = 0.0

           ! Lag Ray tracer (position-wavenumber space method)

           if (rayTracer) then
              call calc_meanFlow_effect(ray,var,force,ray_var3D)
              call transport_rayvol(var, ray, dt, RKstage)  
              if (RKstage == nStages) then
                 call boundary_rayvol(ray)
                 call split_rayvol(ray)
                 call shift_rayvol(ray)
                 call merge_rayvol(ray)
              end if
           end if
           
           ! Reconstruction

           call setHalos( var, "var" )
           call setBoundary (var, flux, "var")

           if( updateMass .or. (testcase=="nIce_w_test") ) call reconstruction (var, "rho")         
           if( updateMass ) call reconstruction (var, "rho")  
           if( updateMass .and. auxil_equ ) then
               call reconstruction (var, "rhop")         
           end if

           if( updateTheta ) call reconstruction (var, "theta") 
           if( predictMomentum .or. (testcase=="nIce_w_test") ) call reconstruction (var, "uvw")
           if(( include_ice ) .and. ( updateIce )) call reconstruction(var, "ice")
           
           call setHalos( var, "varTilde" )
           call setBoundary (var, flux, "varTilde" ) 

           ! Fluxes and Forces
           
           if( updateMass ) then 
              call massFlux (var,var,flux,"nln",PStrat,PStratTilde)
              if( correctDivError ) then
                  print*,'ERROR: correction divergence error not &
                        & implemented properly'
                  stop
              end if
           end if

           if( updateTheta ) then 
              call thetaFlux (var,flux)
              call thetaSource (var,source)
           end if

           if( predictMomentum ) then
              call momentumFlux (var,var,flux,"nln",PStrat,PStratTilde)        
              call volumeForce (var,time,force)        
           end if

           call setBoundary (var, flux, "flux") 
           
           ! Evolve in time
           
        
           if (include_ice .and. updateIce) then
              !--------------------------------------
              !               ice_new
              !--------------------------------------
    
            ! find number of ice time steps per dynamic time step
               ice_time_steps = int(dt * tRef / dt_ice)
               if (ice_time_steps .lt. 1) then
                 ice_time_steps=1
               end if
               dt_ice = dt*tRef / ice_time_steps         

               if (RKstage == 1) then 
                do j = 1,ice_time_steps ! microphysical time steps
                  ! check if time-dependent ice physics is turned on
                  if (iceTestcase_specifics(time+(j-1)*dt_ice/tRef,var)) then
                    dIce = 0.0 ! init q
                    do Ice_RKstage = 1, 3 ! Runge-Kutta loop                    
                      call setBoundary(var,flux,"ice") 
                      call setHalos(var,"ice")
                      call reconstruction(var, "ice")
                      !call setHalos(var,"iceTilde")
                      call setBoundary(var,flux,"iceTilde")
                      call iceSource (var, source)
                      call iceFlux (var, flux)
                      call setBoundary(var,flux,"iceFlux")
                      call iceUpdate(var, var0, flux, source, dt_ice/tRef, dIce, Ice_RKstage)
                      call set_spongeLayer(var, stepFrac(RKstage)*dt_ice/tRef, "ice")
                    end do
                  end if
                end do
              end if
              
           else if(iTime==1 .and. RKstage==1 .and. master) then
              print *,"main: IceUpdate off!"
           end if


           if( updateMass ) then
              ! rho_new

              if ( RKstage==1 ) dRho = 0.                    ! init q
           
              rhoOld = var(:,:,:,1)  ! rhoOld for momentum predictor
              call massUpdate(var, flux, dt, dRho, RKstage, &
                            & "rho", "tot", "expl")
              if (testCase /= 'baroclinic_LC') then
                 call set_spongeLayer(var, stepFrac(RKstage)*dt, "rho")
              end if
              !UAE 200413

              if (auxil_equ) then
                 if ( RKstage==1 ) dRhop = 0.                    ! init q

                 rhopOld = var(:,:,:,6)  ! rhopOld for momentum predictor
                 call massUpdate(var, flux, dt, dRhop, RKstage, &
                               & "rhop", "tot", "expl")
                 if (testCase /= 'baroclinic_LC') then
                    call set_spongeLayer(var, stepFrac(RKstage)*dt, "rhop")
                 end if
              !UAE 200413
              end if

             else 
              if(iTime==1 .and. RKstage==1 .and. master) &
              & print *,"main: MassUpdate off!"
           end if

              
           if( updateTheta ) then
              ! theta_new
           
              call setHalos( var, "var" )
              call setBoundary (var, flux, "var")

              if (RKstage == 1) dTheta = 0.                    ! init q

              call thetaUpdate(var, var0, flux, source, dt, dTheta, &
                      & RKstage)
             else 
              if(iTime==1 .and. RKstage==1 .and. master ) &
              & print *,"main: ThetaUpdate off!"
           end if


           if ( predictMomentum ) then
              ! predictor: uStar
           
              call setHalos( var, "var" )
              call setBoundary (var, flux, "var")

              if ( RKstage==1 ) dMom = 0.                    ! init q

              call momentumPredictor(var, flux, force, dt, dMom, RKstage, &
                                   & "tot", "expl")
              if (testCase /= 'baroclinic_LC') then
                 call set_spongeLayer(var, stepFrac(RKstage)*dt, "uvw")
              end if
             else 
              if(iTime==1 .and. RKstage==1 .and. master) &
              & print *,"main: MomentumUpdate off!"
           end if



           !UAB
           if (shap_dts_dim > 0.) then
              ! smoothing of the fields in order to limit grid-point noise

              select case( timeSchemeType ) 
                 case( "lowStorage" ) 
                    dt_Poisson = beta(RKstage)*dt
                 case( "classical" )
                    dt_Poisson = rk(3,RKstage)*dt
                 case default
                    stop"thetaUpdate: unknown case timeSchemeType"
              end select

              shap_dts = shap_dts_dim/tRef

              fc_shap = min( 1.0, dt_Poisson/shap_dts)

              call smooth_hor_shapiro(fc_shap,n_shap,flux,var,dt_Poisson)
              !call smooth_shapiro(fc_shap,n_shap,flux,var)
           end if
           !UAE

           ! implementation of heating ONeill and Klein 2014

           !UAB
           !if (heatingONK14) then
           ! GBcorr -> FS
           if (heatingONK14 .or. TurbScheme .or. rayTracer) then
           !if (heating) then
           !UAE

              if (model == 'Boussinesq') then
                 print*, "main:ONeill+Klein2014 heating only for &
                        & pseudo-incompressible dyn."
                 stop
              end if

              if ( RKstage==1 ) dPStrat = 0.                    ! init q
              if ( RKstage==1 ) drhoStrat = 0.           

              call BGstate_update(var,flux,dt,RKstage,dPStrat,drhoStrat, &
                                & "expl",heating_switch)
              
           end if

           if( correctMomentum ) then
              ! corrector: dp, du -> u_new, p_new
           
              call setHalos( var, "var" )
              call setBoundary (var,flux,"var")

              select case( timeSchemeType ) 
                 case( "lowStorage" ) 
                    dt_Poisson = beta(RKstage)*dt
                 case( "classical" )
                    dt_Poisson = rk(3,RKstage)*dt
                 case default
                    stop"thetaUpdate: unknown case timeSchemeType"
              end select
              
              call Corrector ( var, flux, dMom, dt_Poisson, errFlagBicg, &
                             & nIterBicg, RKstage, "expl" )
              

              nTotalBicg = nTotalBicg + nIterBicg
           
              ! error handling
              if( errFlagBiCG ) then
                 print*," Momentum Corrector error"
                 write(*,fmt="(a25,i2.2)") "maxIterPoisson! iPara=", &
                                           & iParam
                 nAverageBicg = real(nTotalBicg) / real(iTime) / 3.0

                 call system_clock(count=timeCount)
                 cpuTime = (timeCount - startTimeCount)/real(rate)

                 call output_data(iOut, var, iTime, time, cpuTime)
                 call output_profile(iOut, PStrat,'Pstrat')
                 call output_profile(iOut, rhoStrat,'rhostrat')
                 call output_profile(iOut, thetaStrat,'thetastrat')
                 call output_profile(iOut, bvsstrat,'bvsstrat')
                 
                 if (rayTracer) then
                     call output_wkb(iOut, ray, ray_var3D)
                 end if

                 go to 10   ! dealloc fields
              end if
             else
              if(iTime==1 .and. RKstage==1) &
                   & print *,"main: MomentumCorrector off!"
           end if


        end do Runge_Kutta_Loop
     end if ! timeScheme

     !--------------------------------------------------------------
     !                           Output
     !--------------------------------------------------------------
111 continue
     select case( outputType )
        case( 'time' )
           if (output) then
              if( master ) then
                 call system_clock(count=timeCount)
                 cpuTime = (timeCount - startTimeCount)/real(rate)
              end if

              call output_data(&
                   & iOut, &
                   & var,&
                   & iTime, time, cpuTime)

              if (testCase == 'baroclinic_LC') then
                 call output_profile(iOut, PStrat,'Pstrat.dat')
                 call output_profile(iOut, rhoStrat,'rhostrat.dat')
                 call output_profile(iOut, thetaStrat,'thetastrat.dat')
                 call output_profile(iOut, bvsstrat,'bvsstrat.dat')
              end if
              
              if (rayTracer) then
                 call output_wkb(iOut, ray, ray_var3D)
              end if

              output = .false.
              nextOutputTime = nextOutputTime + outputTimeDiff
              if (nextOutputTime >= maxTime) nextOutputTime = maxTime
           end if
        case( 'timeStep' )
           if (modulo(iTime,nOutput) == 0) then
              if( master ) then   ! modified by Junhong Wei for MPI
                 call system_clock(count=timeCount)
                 cpuTime = (timeCount - startTimeCount)/real(rate)
              end if              ! modified by Junhong Wei for MPI

              call output_data(&
                   & iOut, &
                   & var,&
                   & iTime, time, cpuTime)

              if (testCase == 'baroclinic_LC') then
                 call output_profile(iOut, PStrat,'Pstrat.dat')
                 call output_profile(iOut, rhoStrat,'rhostrat.dat')
                 call output_profile(iOut, thetaStrat,'thetastrat.dat')
                 call output_profile(iOut, bvsstrat,'bvsstrat.dat')
              end if
              
              if (rayTracer) then
                 call output_wkb(iOut, ray, ray_var3D)
              end if
           end if
        case default
           stop "main: unknown outputType"
     end select


     !-------------------------------------------
     !              Abort criteria
     !-------------------------------------------

     if( outputType == "time" ) then
        if( time*tRef >= maxTime ) then 
           if( master ) then
              if(timeScheme == 'semiimplicit')then
                 nAverageBicg = real(nTotalBicg) / real(iTime) / 2.0
              else
                 nAverageBicg = real(nTotalBicg) / real(iTime) / 3.0
              end if

              call system_clock(count=timeCount)
              cpuTime = (timeCount - startTimeCount)/real(rate)

              print*,""
              print*,"=================================================="
              print*," Pinc: Resume"
              print*,""
              write(*,fmt="(a31,es15.1)") &
                   & "average Poisson iterations = ", nAverageBicg
              print*,"=================================================="
              print*,""
           end if

           exit      ! leave time_loop
        end if
     end if

  end do time_loop

     
  !-------------------------------------------
  !      Final output for timeStep
  !-------------------------------------------
     
10 if( master ) then   ! modified by Junhong Wei for MPI (20161103)
     if( outputType == "timeStep" ) then
        nAverageBicg = real(nTotalBicg) / real(iTime-1) / nStages

        call system_clock(count=timeCount)

        cpuTime = (timeCount - startTimeCount)/real(rate)

        print*,""
        print*,"========================================================"
        print*," Pinc: Resume"
        print*,""
        write(*,fmt="(a31,f15.1)") "average Poisson iterations = ", &
                                   & nAverageBicg
        if (preconditioner /= "no") &
        & write(*,fmt="(a25,es15.1)") "ADI: dtau = ", dtau

        print*,"========================================================"
        print*,""
     end if

  end if

  !-------------------------------------
  !       Deallocate variables
  !-------------------------------------

  call terminate_fluxes
  call terminate_poisson
  call terminate (var,var0,var1,flux,force,source,dRho,dRhop,dMom,dTheta,dIce)                         
  call terminate_atmosphere
  call terminate_output

  if (poissonSolverType == 'hypre') then
     call CleanUpHypre       ! Clean Up Hypre objects
    else if (poissonSolverType == 'bicgstab') then
     call CleanUpBiCGSTab       ! Clean Up BiCGSTAB arrays
    else 
     stop 'ERROR: HYPRE or BiCGSTab expected as Poisson solvers'
  end if

666  if( master ) then
        print*,""
        print*,"---------------------------------------"
        print*, "          pincFloit finished          "
        print*,"---------------------------------------"
     end if

  call mpi_finalize(ierror)

end program pinc_prog
