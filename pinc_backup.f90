program pinc_prog

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !% calculates pseudo-incompressible model     %
  !% by Durran with ILES by Adams and Hickel    %
  !% cylflow (c) 2009 Felix Rieper              %
  !% latest revision: January 2010              %
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  use type_module
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
 
  ! test
  use algebra_module

  !----------------------------------------------------
  implicit none
  !---------------------------------------------------

  integer                     :: iTime
  real                        :: time, dt

  ! CPU Time
  integer                     :: rate, startTimeCount, timeCount  
  real                        :: cpuTime 

  ! fields
  real, dimension(:,:,:,:), allocatable :: var, var0
  ! var(i,j,k,iVar) iVar = 1..5 > rho,u,v,w,pExner

  real, dimension(:,:,:), allocatable :: dRho      ! RK-Update for rho
  real, dimension(:,:,:,:), allocatable :: dMom    ! RK for rhoU,rhoV,rhoW
  real, dimension(:,:,:), allocatable :: dTheta     ! RK-Update for theta


  real, dimension(:,:,:,:,:), allocatable :: flux
  ! flux(i,j,k,dir,iFlux) 
  ! dir = 1..3 > f,g,h-flux in x,y,z-direction
  ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW

 
  !--------------------
  !   WKB variables
  !--------------------
  real,          dimension(:,:,:), allocatable :: waveAct, waveActOld
  type(rayType), dimension(:), allocatable     :: ray
  ! complex wave amplitudes
  complex, dimension(:,:,:,:,:), allocatable :: Psi
  ! Psi(i,j,k,variable b|u|w|pi,harmonic 0|1|2)

  

  ! topography via force field
  real, dimension(:,:,:,:), allocatable :: force ! volume forces
  ! force(i,j,k,forceVector)

  ! output per timeStep
  logical :: output
  real :: nextOutputTime      ! scaled time for next output
  integer :: iOut             ! output counter

  ! general
  integer :: i,j,k,l

  ! restart
  logical :: scale
  ! test
  !real, dimension(-1:1) :: a
  !real, dimension(3) :: b

  ! parameter study
  integer :: iParam

  ! error handling and controlling
  logical :: errFlagBiCG, errFlagTStep
  integer :: nIterBicg, nTotalBicg
  real :: nAverageBicg
  character (len=8) :: datum
  character (len=10) :: niceDatum
  character (len=10) :: zeit
  character (len=9) :: niceZeit
  real :: a2



  call date_and_time(date=datum, time=zeit)
  niceDatum = datum(7:8)//"."//datum(5:6)//"."//datum(1:4)
  niceZeit = zeit(1:2)//":"//zeit(3:4)//" Uhr"
  print*,""
  print*,""
  print*,""
  print*,"=============================================================="
  print*," PincFloit (c) 2010 F. Rieper  "
  print*,""
  write(*,fmt="(a25,a10)") " Today : ", niceDatum
  write(*,fmt="(a25,a9)") " Time  : ", niceZeit




  !----------------------------------------------
  !              Parameter study
  !----------------------------------------------

  call init_paramStudy             ! read parameter data from input.f90

  if( .not. parameterStudy ) then
     startParam = 1
     endParam = 1
     stepParam = 1
  end if


  paramLoop: do iParam = startParam, endParam, stepParam



     !-------------------------------------------------
     !                    Set up 
     !-------------------------------------------------

     ! init counter and time
     iOut = 0
     iTime = 0; time = 0.0; cpuTime=0.0
     nTotalBicg = 0

     ! init clock
     call system_clock (count_rate=rate)
     call system_clock (count=startTimeCount)

     ! 1) allocate variables 
     ! 2) read input.f90
     call setup (var,var0,flux,force,dRho,dMom,dTheta)


     if( iParam == startParam ) then
        write(*,fmt="(a25,a)") "Test Case : ", testCase
        print*,"=============================================================="
     end if


     call init_atmosphere      ! set atmospheric background state 
     call initialise (var)     ! set initial conditions

     call init_xweno       ! set ILES parameters 
     call init_fluxes      ! allocate tilde variables
     call init_update
     call init_timeScheme  ! define Runge-Kutta parameters
     call init_poisson     ! allocate dp
     call init_output
     
     !-------------------------------------------------
     !              read initial data
     !-------------------------------------------------
     if (restart) then
        scale = .true.                       ! scale with reference quantities
        call readtec360( var,restartFile,time, scale)
        if( maxTime < time*tRef ) stop"restart error: maxTime < current time"
     end if

     !-------------------------------------------------
     !              set up parameter study 
     !-------------------------------------------------
     if( parameterStudy ) then
        print*,""
        print*,""
        print*,""
        print*,"=============================================================="
        print*," Parameter study "
        print*,""
        write(*,fmt="(a25,i15)") " iParam = ", iParam
        print*,"=============================================================="
     end if


     ! ---------------------------------------------
     !        initial divergence cleaning
     ! ---------------------------------------------

     if( initialCleaning ) then

        ! correct u,v,w and pi
        call momentumCorrector (var, 1.0, errFlagBicg, nIterBicg,1)

        call horizontalBoundary (var, "uvw")     ! periodic BC for u,v,w
        call verticalBoundary(var,flux,"uvw")    ! solid wall BC for u,v,w
        call reconstruction (var, "uvw")         ! reconstruct u,v,w

        call horizontalBoundary(var, "p")     ! periodic BC for p
        call verticalBoundary(var,flux,"p")   ! set pGrad at vert. bounds to 0
     end if


     !---------------------------------------------
     !               Init ray tracer
     !---------------------------------------------

     call setup_wkb(ray, waveAct, waveActOld, Psi)   ! allocate ray fields


     !------------------------------------------
     !            initial output
     !------------------------------------------

     call tec360 (&          
          & iOut,&
          & var,&
          & ray, waveAct, Psi,&
          & iTime, time, cpuTime, dt, &
          & iParam)


     
     output = .false.
     nextOutputTime = time*tRef + outputTimeDiff     ! set time for first output
     ! and consider restart time


     !-----------------------------------------------------
     !                        Time loop
     !-----------------------------------------------------

  
     print *, "...starting time loop"
     if( outputType == "time" ) maxIter = 2**30
 
     time_loop: do iTime = 1,maxIter

        print*,""
        print*,""
        print*,""
        print*,"--------------------------------------------------------------"
        write(*,fmt="(a25,i15)") " Time step = ", iTime
        write(*,fmt="(a25,f15.1,a8)") " Time = ", time*tRef, "seconds"
        print*,"--------------------------------------------------------------"


        ! calc time step
        call timestep (var,ray,dt,errFlagTstep)

        ! error handling for parameter study
        if( errFlagTstep ) then
           print*," TimeStep routine: "
           write(*,fmt="(a25,es15.4,a8)") &
                &"dt < dtMin!!! dtMin = ", dtMin_dim, "seconds"
           write(*,fmt="(a25,es15.4,a8)") "dt * tRef = ", dt*tRef, "seconds"
           write(*,fmt="(a25,i2.2)") "cycle paramLoop. iPara = ", iParam
           call system_clock(count=timeCount)
           cpuTime = (timeCount - startTimeCount)/real(rate)

           ! final output
           call tec360(&
                & iOut, &
                & var,&
                & ray, waveAct, Psi, &
                & iTime, time, cpuTime, dt,&
                & iParam)   

           go to 10    ! deallocate fields
        end if


        ! correct dt to hit desired output time for outputType 'time'
        if ( outputType == 'time' ) then
           if( (time+dt) * tRef >= nextOutputTime ) then
              dt = nextOutputTime/tRef - time      
              output = .true.
              write(*,fmt="(a25,es15.4,a8)") "dt for output = ", &
                   & dt*tRef, "seconds"
           end if
        end if



        ! update will be for new time level
        time = time + dt


        !---------------------------------------------------------------
        !                     Runge-Kutta stages
        !---------------------------------------------------------------

        Runge_Kutta_Loop: do RKstage = 1, nStages

           if( verbose ) then
              print*,""
              print*,"-----------------------------------------"
              write(*,fmt="(a25,i1)") "Runge-Kutta stage ", RKstage
              print*,"-----------------------------------------"
           end if


           !-----------------------------
           !         ray tracer
           !-----------------------------
           if (raytracer) then
              
              ! save current wave action for ray transport
              waveActOld = waveAct              
              
              call transport_waveAction(var,ray,waveAct, RKstage, dt)
              
              ! use waveActOld if needed:
              call transport_ray(var, ray, dt, RKstage)  
              call calc_cellIndex(ray)
              call calc_waveAmplitude(ray, waveAct, Psi)
              call calc_meanFlow(Psi,dt,RKstage)

           end if ! ray tracer



           select case( model ) 
              
           case( "pseudo_incompressible" ) 
              !-------------------------------------------
              !         reconstruction of density
              !-------------------------------------------

              call horizontalBoundary (var, "rho")     ! periodic BC for rho
              call verticalBoundary (var, flux, "rho" )! solid wall BC for rho
              call reconstruction (var, "rho")         

              ! periodic BC for rhoTilde
              call horizontalBoundary (var,"rhoTilde")


           case( "Boussinesq" )
              !-------------------------------------------
              !  reconstruction of potential temperature
              !-------------------------------------------

              call horizontalBoundary (var, "theta")     ! periodic BC for rho
              call verticalBoundary (var, flux, "theta" )! solid wall BC for rho
              ! reconstruct rhoTilde(RKstage+1) and save rhoTilde(RKstage):
              call reconstruction (var, "theta")         

              ! periodic BC for rhoTilde
              call horizontalBoundary (var,"thetaTilde")
              
           case default
              stop"pinc: unknown model"
           end select

              



           !-------------------------------------------
           !         reconstruction of velocity
           !-------------------------------------------

           ! set BC for u,v,w
           call horizontalBoundary (var, "uvw")    ! perdiodic BC for u,v,w
           call verticalBoundary (var,flux,"uvw")  ! solid wall BC for u,v,w

           call reconstruction (var, "uvw")         ! reconstruct u,v,w




           !-------------------------------------------------
           !                Transport of scalars
           !-------------------------------------------------

           select case( model )
              
           case( "pseudo_incompressible" )

              !---------------------------------------
              !               rho(RK_stage+1)
              !---------------------------------------

              if( updateMass ) then

                 call massFlux (var,flux)                 ! calc fRho,...

                 ! set hRho = 0 at vert. bounds
                 call verticalBoundary (var,flux,"rhoFlux")

                 ! use CDS at boundary
                 call verticalBoundary(var,flux,"rhoFluxCorr")

                 ! set hRho for topographic boundary (topograpy1)
                 call bottomTopography(var,flux,"rhoFlux")

                 if (RKstage == 1) dRho = 0.                    ! init q

                 ! save old density -> rhoOld and make RK-update
                 rhoOld = var(:,:,:,1)


                 call massUpdate(var, flux, dt,&
                      & dRho, &          ! corresponds to q(RKstage) or k(RKstage)
                      & RKstage)         ! Runge-Kutta stage 

                 ! apply sponge layer to density (relaxation to background)
                 call set_spongeLayer(var, stepFrac(RKstage)*dt, "rho")


              else 
                 if(iTime==1 .and. RKstage==1) print *,"pinc.f90: MassUpdate off!"
              end if


           case( "Boussinesq" )
              
              !---------------------------------------
              !               theta(RK_stage+1)
              !---------------------------------------
              
              if( updateTheta ) then

                 call thetaFlux (var,flux)                 ! calc fTheta,...

                 ! set hRho = 0 at vert. bounds
                 call verticalBoundary (var,flux,"thetaFlux")
                 
                 ! use CDS at boundary
                 call verticalBoundary(var,flux,"thetaFluxCorr")

                 if (RKstage == 1) dTheta = 0.                    ! init q
                 ! save old theta -> thetaOld and make RK-update
                 thetaOld = var(:,:,:,6)

                 call thetaUpdate(var, flux, dt,&
                      & dTheta, &        ! corresponds to q(RKstage) or k(RKstage)
                      & RKstage)         ! Runge-Kutta stage 

                 ! apply sponge layer to density (relaxation to background)
!                 call set_spongeLayer(var, stepFrac(RKstage)*dt, "theta")


              else 
                 if(iTime==1 .and. RKstage==1) print *,"pinc.f90: MassUpdate off!"
              end if

           case default
              stop"main: unknown case model"
           end select


           !-----------------------------------------
           !             predictor: uStar(RKstage+1)
           !-----------------------------------------

           if ( predictMomentum ) then

              call momentumFlux (var,flux)          ! calc fRhoU,...
              call volumeForce (var,force)          ! calc volume forces
              
              ! set w at k=1 to imitate bottom topography
              call bottomTopography(var,flux,"w")   
              
              ! set hRhoU,V,W = 0 at vert. bounds
              call verticalBoundary(var,flux,"uvwFlux")

              ! use CDS at boundary
              call verticalBoundary(var,flux,"uvwFluxCorr")

              if ( RKstage==1 ) dMom = 0.                    ! init q

              call momentumPredictor(var, flux, force, dt, &
                   & dMom, &          ! corresponds to q(RKstage) or k(RKstage)
                   & RKstage)         ! Runge-Kutta stage

              ! apply sponge layer to velocity
              call set_spongeLayer(var, stepFrac(RKstage)*dt, "uvw")
              
              call horizontalBoundary (var,"uvw")    ! periodic BC for u,v,w
              call verticalBoundary(var,flux,"uvw")  ! solid wall BC for u,v,w

           else 
              if(iTime==1 .and. RKstage==1) &
                   & print *,"pinc.f90: MomentumUpdate off!"
           end if


           ! ---------------------------------------------
           !           corrector: dp, du -> u(RKstage), p(RKstage-1)
           ! ---------------------------------------------

           if( correctMomentum ) then

              ! periodic BC for rho needed for PoissonSolver
              call horizontalBoundary (var, "rho")     
              call momentumCorrector ( var, beta(RKstage)*dt, &
                   & errFlagBicg, nIterBicg, RKstage)

              nTotalBicg = nTotalBicg + nIterBicg


              ! error handling
              if( errFlagBiCG ) then
                 print*," Momentum Corrector error"
                 write(*,fmt="(a25,i2.2)") "maxIterPoisson! iPara=", iParam
                 nAverageBicg = real(nTotalBicg) / real(iTime) / 3.0
                 call system_clock(count=timeCount)
                 cpuTime = (timeCount - startTimeCount)/real(rate)
                 call tec360(&
                      & iOut,&
                      & var,&
                      & ray, waveAct, Psi, &
                      & iTime, time, cpuTime, dt,&
                      & iParam)

                 go to 10   ! dealloc fields
              end if

              ! set BC for p
              call horizontalBoundary (var, "p")   ! periodic BC for p
              call verticalBoundary (var,flux,"p") !set pGrad=0 at vert. bounds

           else
              if(iTime==1 .and. RKstage==1) &
                   & print *,"pinc.f90: MomentumCorrector off!"
           end if


        end do Runge_Kutta_Loop



        !------------------------------------------
        !                 output
        !------------------------------------------

        select case( outputType )

        case( 'time' )

           if (output) then
              call system_clock(count=timeCount)
              cpuTime = (timeCount - startTimeCount)/real(rate)
              call tec360(&
                   & iOut, &
                   & var,&
                   & ray, waveAct, Psi, &
                   & iTime, time, cpuTime, dt,&
                   & iParam)

              output = .false.
              nextOutputTime = nextOutputTime + outputTimeDiff
              if (nextOutputTime >= maxTime) nextOutputTime = maxTime
           end if

        case( 'timeStep' )

           if (modulo(iTime,nOutput) == 0) then
              call system_clock(count=timeCount)
              cpuTime = (timeCount - startTimeCount)/real(rate)
              call tec360(&
                   & iOut,&
                   & var,&
                   & ray, waveAct, Psi, &
                   & iTime, time, cpuTime, dt,&
                   & iParam)
           end if

        case default
           stop"pinc.f90: unknown outputType"
        end select


        !-------------------------------------------
        !              Abort criteria
        !-------------------------------------------

        if( outputType == "time" ) then
           if( time*tRef >= maxTime ) then 
              nAverageBicg = real(nTotalBicg) / real(iTime) / 3.0
              call system_clock(count=timeCount)
              cpuTime = (timeCount - startTimeCount)/real(rate)
              print*,""
              print*,"=========================================================="
              print*," Pinc: Resume"
              print*,""
              write(*,fmt="(a31,es15.1)") "average Poisson iterations = ", nAverageBicg
              print*,"=========================================================="
              print*,""
              exit      ! leave time_loop
           end if
        end if


     end do time_loop

     
     !-------------------------------------------
     ! Final output for timeStep
     !-------------------------------------------
     
     if( outputType == "timeStep" ) then
        nAverageBicg = real(nTotalBicg) / real(iTime-1) / 3.0
        call system_clock(count=timeCount)
        cpuTime = (timeCount - startTimeCount)/real(rate)
        print*,""
        print*,"=========================================================="
        print*," Pinc: Resume"
        print*,""
        write(*,fmt="(a31,f15.1)") "average Poisson iterations = ", nAverageBicg
        if (preconditioner /= "no") write(*,fmt="(a25,es15.1)") "ADI: dtau = ", dtau
        print*,"=========================================================="
        print*,""
     end if




     !-------------------------------------
     !       deallocate variables
     !-------------------------------------

10   call terminate_fluxes
     call terminate_poisson
     call terminate (var,var0,dRho,dMom,dTheta)                         
     call terminate_atmosphere
     call terminate_output
     call finish_wkb(ray, waveAct, waveActOld, Psi)

  end do paramLoop





  print*,""
  print*,"---------------------------------------"
  print*, "          pincFloit finished          "
  print*,"---------------------------------------"

end program pinc_prog
