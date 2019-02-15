program pinc_prog

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !% calculates pseudo-incompressible model     %
  !% by Durran with ILES by Adams and Hickel    %
  !% cylflow (c) 2009 Felix Rieper              %
  !% latest revision: January 2010              %
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
  use type_module
  use mpi_module        ! modified by Junhong Wei for MPI (20161102)
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

  ! MPI stuff             ! modified by Junhong Wei for MPI (20161102)
  logical :: error_flag   ! modified by Junhong Wei for MPI (20161102)
  real :: dt_local        ! modified by Junhong Wei for MPI (20161102)

  ! fields
  real, dimension(:,:,:,:), allocatable :: var, var0
  real, dimension(:,:,:,:), allocatable :: source
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
  real    :: nextOutputTime   ! scaled time for next output
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


  file_namelist = 'input.f90'   ! default
  if ( command_argument_count() /= 0 )  &
    &  call get_command_argument(1,file_namelist)  ! take the given one

  !-------------------------------------------------   ! modified by Junhong Wei for MPI (20161103)
  !                    Set up   ! modified by Junhong Wei for MPI (20161103)
  !-------------------------------------------------   ! modified by Junhong Wei for MPI (20161103)

  ! init counter and time   ! modified by Junhong Wei for MPI (20161103)
  iOut = 0   ! modified by Junhong Wei for MPI (20161103)
  iTime = 0; time = 0.0; cpuTime=0.0   ! modified by Junhong Wei for MPI (20161103)



  call init_mpi(error_flag)      ! modified by Junhong Wei for MPI (20161102)

  ! modified by Junhong Wei for MPI (20161102) *** starting line ***
!  if( error_flag ) goto 10       ! modified by Junhong Wei for MPI (20161102)
  if( error_flag )then
  print*,"error in the init_mpi"
  goto 666
  end if


  if( master ) then   ! modified by Junhong Wei for MPI (20161102)
  call date_and_time(date=datum, time=zeit)
  niceDatum = datum(7:8)//"."//datum(5:6)//"."//datum(1:4)
  niceZeit = zeit(1:2)//":"//zeit(3:4)//" Uhr"
  print*,""
  print*,""
  print*,""
  print*,"=============================================================="
  print*," PincFloit (c) 2010 F. Rieper  "
  print*," modified by Junhong Wei (2018)  "   ! modified by Junhong Wei
  print*,""
  write(*,fmt="(a25,a10)") " Today : ", niceDatum
  write(*,fmt="(a25,a9)") " Time  : ", niceZeit

     ! init clock   ! modified by Junhong Wei for MPI (20161102)
     call system_clock (count_rate=rate)        ! modified by Junhong Wei for MPI (20161102)
     call system_clock (count=startTimeCount)   ! modified by Junhong Wei for MPI (20161102)

  end if   ! modified by Junhong Wei for MPI (20161102)




  !----------------------------------------------
  !              Parameter study
  !----------------------------------------------

  call init_paramStudy             ! read parameter data from the namelist

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
!     iOut = 0   ! modified by Junhong Wei for MPI (20161103)
!     iTime = 0; time = 0.0; cpuTime=0.0   ! modified by Junhong Wei for MPI (20161103)
     nTotalBicg = 0

     ! init clock
!     call system_clock (count_rate=rate)   ! modified by Junhong Wei for MPI (20161103)
!     call system_clock (count=startTimeCount)   ! modified by Junhong Wei for MPI (20161103)

     ! 1) allocate variables 
     ! 2) read the namelist
     call setup (var,var0,flux,force,source,dRho,dMom,dTheta)


  if( master ) then   ! modified by Junhong Wei for MPI (20161103)
     if( iParam == startParam ) then
        write(*,fmt="(a25,a)") "Test Case : ", testCase
        print*,"=============================================================="
     end if
  end if              ! modified by Junhong Wei for MPI (20161103)


     call init_atmosphere      ! set atmospheric background state 
     call init_output
     call initialise (var)     ! set initial conditions

     call init_xweno       ! set ILES parameters 
     call init_fluxes      ! allocate tilde variables
     call init_update
     call init_timeScheme  ! define Runge-Kutta parameters
     call init_poisson     ! allocate dp
     
     !-------------------------------------------------
     !              Read initial data
     !-------------------------------------------------
     if (restart) then
        !achatzb
        !restart from previous output
        ! scale = .true.  ! scale with reference quantities
        ! call readtec360( var,restartFile,time, scale)
        ! if( maxTime < time*tRef ) &
        ! & stop "restart error: maxTime < current time"
 
        call read_data ( iIn,var)

        call setHalos( var, "var" )
        call setBoundary (var, flux, "var")

!       testb
!       print*,"r(0,1,310),r(nx,1,310) =",var(0,1,310,1),var(nx,1,310,1)
!       print*,"r(1,1,310),r(nx+1,1,310) =",var(1,1,310,1),var(nx+1,1,310,1)

!       print*,"u(0,1,310),u(nx,1,310) =",var(0,1,310,2),var(nx,1,310,2)
!       print*,"u(1,1,310),u(nx+1,1,310) =",var(1,1,310,2),var(nx+1,1,310,2)

!       print*,"v(0,1,310),v(nx,1,310) =",var(0,1,310,3),var(nx,1,310,3)
!       print*,"v(1,1,310),v(nx+1,1,310) =",var(1,1,310,3),var(nx+1,1,310,3)

!       print*,"w(0,1,310),w(nx,1,310) =",var(0,1,310,4),var(nx,1,310,4)
!       print*,"w(1,1,310),w(nx+1,1,310) =",var(1,1,310,4),var(nx+1,1,310,4)

!       print*,"p(0,1,310),p(nx,1,310) =",var(0,1,310,5),var(nx,1,310,5)
!       print*,"p(1,1,310),p(nx+1,1,310) =",var(1,1,310,5),var(nx+1,1,310,5)
!       teste
        ! achatze
     end if

     !-------------------------------------------------
     !              Set up parameter study 
     !-------------------------------------------------
  if( master ) then   ! modified by Junhong Wei for MPI (20161103)
     if( parameterStudy ) then
        print*,""
        print*,""
        print*,""
        print*,"==========================================================="
        print*," Parameter study "
        print*,""
        write(*,fmt="(a25,i15)") " iParam = ", iParam
        print*,"==========================================================="
     end if
  end if              ! modified by Junhong Wei for MPI (20161103)


     !---------------------------------------------
     !        Initial divergence cleaning
     !---------------------------------------------

     if( initialCleaning ) then
        ! correct u,v,w and pi
        call momentumCorrector (var, dMom, 1.0, errFlagBicg, nIterBicg,1,'initial')
        call setHalos( var, "var" )
        call setBoundary (var,flux,"var")
     end if


     !---------------------------------------------
     !               Init ray tracer
     !---------------------------------------------

!    achatzb
!    initialization WKB would have to be redone!
!    call setup_wkb(ray, waveAct, waveActOld, Psi)   ! allocate ray fields
!    achatze


     !------------------------------------------
     !              Initial output
     !------------------------------------------

!    achatzb
!    reduced argument list for output
!    changed name of output routine
!    call tec360(&
!         & iOut, &
!         & var,&
!         & ray, waveAct, Psi, &
!         & iTime, time, cpuTime, dt,&
!         & iParam)   
     call output_data(&
          & iOut, &
          & var,&
          & iTime, time, cpuTime)
!    achatze

     output = .false.
     nextOutputTime = time*tRef + outputTimeDiff     ! set time for first output
     ! and consider restart time


     !-----------------------------------------------------
     !                        Time loop
     !-----------------------------------------------------

  if( master ) then   ! modified by Junhong Wei for MPI (20161103)
     print *, "...starting time loop"
  end if              ! modified by Junhong Wei for MPI (20161103)

     if( outputType == "time" ) maxIter = 2**30
 
     time_loop: do iTime = 1,maxIter

  if( master ) then   ! modified by Junhong Wei for MPI (20161103)
        print*,""
        print*,""
        print*,""
        print*,"--------------------------------------------------------------"
        write(*,fmt="(a25,i15)") " Time step = ", iTime
        write(*,fmt="(a25,f15.1,a8)") " Time = ", time*tRef, "seconds"
        print*,"--------------------------------------------------------------"
  end if              ! modified by Junhong Wei for MPI (20161103)


        !----------------------------------
        !         Calc time step
        !----------------------------------
!        call timestep (var,ray,dt,errFlagTstep)   ! modified by Junhong Wei for MPI (20161103)

       call timestep (var,ray,dt_local,errFlagTstep)   ! modified by Junhong Wei for MPI (20161103)


! modified by Junhong Wei for MPI (20161103)   *** starting line ***
     ! find global maximum
     call mpi_reduce(dt_local, dt, 1, mpi_double_precision,&
          & mpi_min, root, comm, ierror)

     call mpi_bcast(dt, 1, mpi_double_precision, root, comm, ierror)

     ! abort if time step too small
     if( dt*tRef < dtMin_dim ) then


           if( master ) then

           print*," TimeStep routine: "
           write(*,fmt="(a25,es15.4,a8)") &
                &"dt < dtMin!!! dtMin = ", dtMin_dim, "seconds"
           write(*,fmt="(a25,es15.4,a8)") "dt * tRef = ", dt*tRef, "seconds"
           write(*,fmt="(a25,i2.2)") "cycle paramLoop. iPara = ", iParam
           call system_clock(count=timeCount)
           cpuTime = (timeCount - startTimeCount)/real(rate)

           end if


           ! final output

!          achatzb
!          reduced argument list for output
!          changed name of output routine
!          call tec360(&
!               & iOut, &
!               & var,&
!               & ray, waveAct, Psi, &
!               & iTime, time, cpuTime, dt,&
!               & iParam)   
           call output_data(&
                & iOut, &
                & var,&
                & iTime, time, cpuTime)
!          achatze

        call mpi_barrier(comm, ierror)
        call mpi_finalize(ierror)
        print*,"pinc.f90: time step too small. dt = ", dt*tRef
        stop
     end if
! modified by Junhong Wei for MPI (20161103)   *** finishing line ***


        if( timeSchemeType == "classical" .and. iTime == 1 ) dt = 0.5*dt

        ! correct dt to hit desired output time for outputType 'time'
        if ( outputType == 'time' ) then
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


        !---------------------------------------------------------------
        !                     Runge-Kutta stages
        !---------------------------------------------------------------

        Runge_Kutta_Loop: do RKstage = 1, nStages

        ! modified by Junhong Wei (20161007) --- starting line
        if(DySmaScheme)then
        call CoefDySma_update(var)
        end if
        ! modified by Junhong Wei (20161007) --- finishing line

           if( timeSchemeType == "classical" ) then
              if( RKstage == 1 ) then
                 var0 = var
              end if
           end if

  if( master ) then   ! modified by Junhong Wei for MPI (20161103)
           if( verbose ) then
              print*,""
              print*,"-----------------------------------------"
              write(*,fmt="(a25,i1)") "Runge-Kutta stage ", RKstage
              print*,"-----------------------------------------"
           end if
  end if              ! modified by Junhong Wei for MPI (20161103)


           !------------------------------------------------------
           !                      Ray tracer (WKB part)
           !------------------------------------------------------
           if (raytracer) then
              
              ! save current wave action for ray transport
              waveActOld = waveAct         
              call transport_waveAction(var,ray,waveAct, RKstage, dt)
              call transport_ray(var, ray, dt, RKstage)  
              call calc_cellIndex(ray)
              call calc_waveAmplitude(ray, waveAct, Psi)
              call calc_meanFlow(Psi,var,dt,RKstage)

              ! reset rays to original position
              ! poor results -> do not uncomment
              !if( RKstage == nStages ) then
              !   call calc_waveNumber(ray)
              !   call reset_ray(ray)
              !end if

           end if


           !---------------------------------------------------------
           !                      Reconstruction
           !---------------------------------------------------------

           call setHalos( var, "var" )
           call setBoundary (var, flux, "var")

           if( updateMass ) call reconstruction (var, "rho")         
           if( updateTheta ) call reconstruction (var, "theta") 
           if( predictMomentum ) call reconstruction (var, "uvw")
           
           
           !------------------------------------------------------------
           !                     Fluxes and Forces
           !------------------------------------------------------------
           
           call setHalos( var, "varTilde" )
           call setBoundary (var, flux, "varTilde" ) 

           if( updateMass ) then 
              call massFlux (var,flux)
              if( correctDivError ) call massSource (var,source)
           end if
           if( updateTheta ) then 
              call thetaFlux (var,flux)
              call thetaSource (var,source)
           end if
           if( predictMomentum ) then
              call momentumFlux (var,flux)        
              call volumeForce (var,time,force)        
!             achatzb deactivated old implementation of topography
!             call bottomTopography(var,flux)
!             achatze
!xxxx wrong implemented / call might not be necessary
              if( correctDivError) call momentumSource (var,source)
           end if
           call setBoundary (var, flux, "flux") 

           
           !------------------------------------------------------------
           !                        Evolve in time
           !------------------------------------------------------------
           
           if( updateMass ) then
              !---------------------------------------
              !               rho_new
              !---------------------------------------
              
              if (RKstage == 1) dRho = 0.0                ! init q
              rhoOld = var(:,:,:,1)  ! rhoOld for momentum predictor
              call massUpdate(var, var0, flux, source, dt, dRho, RKstage)
              call set_spongeLayer(var, stepFrac(RKstage)*dt, "rho")
           else 
!              if(iTime==1 .and. RKstage==1) print *,"main: MassUpdate off!"   ! modified by Junhong Wei for MPI (20161103)
              if(iTime==1 .and. RKstage==1 .and. master) print *,"main: MassUpdate off!"   ! modified by Junhong Wei for MPI (20161103)
           end if
           
              
           if( updateTheta ) then
              !---------------------------------------
              !               theta_new
              !---------------------------------------
              
              call setHalos( var, "var" )
              call setBoundary (var, flux, "var")

              if (RKstage == 1) dTheta = 0.                    ! init q
              call thetaUpdate(var, var0, flux, source, dt, dTheta, RKstage) 

           else 
!              if(iTime==1 .and. RKstage==1) print *,"main: ThetaUpdate off!"
              if(iTime==1 .and. RKstage==1 .and. master ) print *,"main: ThetaUpdate off!"
           end if
           

           if ( predictMomentum ) then
              !-----------------------------------------
              !             predictor: uStar
              !-----------------------------------------
              
              call setHalos( var, "var" )
              call setBoundary (var, flux, "var")

              if ( RKstage==1 ) dMom = 0.                    ! init q
              call momentumPredictor(var, var0, flux, source, &
                   &                 force, dt, dMom, RKstage)
              call set_spongeLayer(var, stepFrac(RKstage)*dt, "uvw")
           else 
!              if(iTime==1 .and. RKstage==1) print *,"main: MomentumUpdate off!"   ! modified by Junhong Wei for MPI (20161103)
              if(iTime==1 .and. RKstage==1 .and. master) print *,"main: MomentumUpdate off!"   ! modified by Junhong Wei for MPI (20161103)
           end if
           
           
           if( correctMomentum ) then
              !---------------------------------------------
              !      corrector: dp, du -> u_new, p_new
              !---------------------------------------------
              
              call setHalos( var, "var" )
              call setBoundary (var,flux,"var")

!xxxx new time scheme
              select case( timeSchemeType ) 
                 
              case( "lowStorage" ) 
                 
                 dt_Poisson = beta(RKstage)*dt
                 
              case( "classical" )
                 
                 dt_Poisson = rk(3,RKstage)*dt

              case default
                 stop "thetaUpdate: unknown case timeSchemeType"
              end select
!xxxx end
              
              call momentumCorrector ( var, dMom, dt_Poisson, &
                   & errFlagBicg, nIterBicg, RKstage,'')

              nTotalBicg = nTotalBicg + nIterBicg
              
              ! error handling
              if( errFlagBiCG ) then
                 print*," Momentum Corrector error"
                 write(*,fmt="(a25,i2.2)") "maxIterPoisson! iPara=", iParam
                 nAverageBicg = real(nTotalBicg) / real(iTime) / 3.0

                 call system_clock(count=timeCount)
                 cpuTime = (timeCount - startTimeCount)/real(rate)

!                achatzb
!                reduced argument list for output
!                changed name of output routine
!                call tec360(&
!                     & iOut, &
!                     & var,&
!                     & ray, waveAct, Psi, &
!                     & iTime, time, cpuTime, dt,&
!                     & iParam)   
                 call output_data(&
                      & iOut, &
                      & var,&
                      & iTime, time, cpuTime)
!                achatze

                 go to 10   ! dealloc fields
              end if
           else
              if(iTime==1 .and. RKstage==1) &
                   & print *,"main: MomentumCorrector off!"
           end if


        end do Runge_Kutta_Loop




        !--------------------------------------------------------------
        !                           Output
        !--------------------------------------------------------------

        select case( outputType )

        case( 'time' )

           if (output) then
              if( master ) then   ! modified by Junhong Wei for MPI
                 call system_clock(count=timeCount)
                 cpuTime = (timeCount - startTimeCount)/real(rate)
              end if              ! modified by Junhong Wei for MPI

!             achatzb
!             reduced argument list for output
!             changed name of output routine
!             call tec360(&
!                  & iOut, &
!                  & var,&
!                  & ray, waveAct, Psi, &
!                  & iTime, time, cpuTime, dt,&
!                  & iParam)   
              call output_data(&
                   & iOut, &
                   & var,&
                   & iTime, time, cpuTime)
!             achatze

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

!             achatzb
!             reduced argument list for output
!             changed name of output routine
!             call tec360(&
!                  & iOut, &
!                  & var,&
!                  & ray, waveAct, Psi, &
!                  & iTime, time, cpuTime, dt,&
!                  & iParam)   
              call output_data(&
                   & iOut, &
                   & var,&
                   & iTime, time, cpuTime)
!             achatze
           end if

        case default
           stop "main: unknown outputType"
        end select


        !-------------------------------------------
        !              Abort criteria
        !-------------------------------------------

        if( outputType == "time" ) then
           if( time*tRef >= maxTime ) then 
              if( master ) then   ! modified by Junhong Wei for MPI
                 nAverageBicg = real(nTotalBicg) / real(iTime) / 3.0

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
              end if              ! modified by Junhong Wei for MPI

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
        print*,"=========================================================="
        print*," Pinc: Resume"
        print*,""
        write(*,fmt="(a31,f15.1)") "average Poisson iterations = ", nAverageBicg
        if (preconditioner /= "no") write(*,fmt="(a25,es15.1)") "ADI: dtau = ", dtau
        print*,"=========================================================="
        print*,""
     end if
  end if   ! modified by Junhong Wei for MPI (20161103)



     
     !-------------------------------------
     !       Deallocate variables
     !-------------------------------------

!10   call terminate_fluxes   ! modified by Junhong Wei for MPI (20161102)
     call terminate_fluxes    ! modified by Junhong Wei for MPI (20161102)
     call terminate_poisson
     call terminate(var,var0,flux,force,source,dRho,dMom,dTheta)
     call terminate_atmosphere
     call terminate_output
!    achatzb
!    WKB switched off
!    call finish_wkb(ray, waveAct, waveActOld, Psi)
!    achatze

  end do paramLoop




666  if( master ) then   ! modified by Junhong Wei for MPI (20161102)
  print*,""
  print*,"---------------------------------------"
  print*, "          pincFloit finished          "
  print*,"---------------------------------------"
  end if   ! modified by Junhong Wei for MPI (20161102)

  call mpi_finalize(ierror)   ! modified by Junhong Wei for MPI (20161102)

end program pinc_prog
