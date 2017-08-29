module output_module

  use type_module
  use atmosphere_module
  use wkb_module, only: cabs


  implicit none

  private

  ! public subroutines
  public :: tec360
  public :: readtec360

  ! internal subroutines (listed for completeness)
  public :: init_output
  public :: terminate_output

  ! internal module variables (output variables)
  real, dimension(:,:,:), allocatable :: optVar


contains
  
  subroutine tec360( &
       & iOut, &
       & var,&
       & ray,waveAct,Psi,&
       & iTime,time,cpuTime,dt, &
       & iParam )
    !-------------------------------
    !  creates tecplot .dat-files
    !-------------------------------
    
    ! output counter
    integer, intent(inout) :: iOut
    
    ! argument fields
    real,dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),intent(in)  :: var
    
    ! wkb arguments
    type(rayType), dimension(nRay), intent(in) :: ray
    real,dimension(0:nx+1,0:ny+1,0:nz+1),intent(in)  :: waveAct
    complex, dimension(0:nx+1,0:ny+1,0:nz+1,4,0:2) :: Psi

    ! argument parameters
    integer, intent(in)            :: iTime
    real,intent(in)                :: time, cpuTime
    real, intent(in) :: dt
    integer, intent(in) :: iParam

    ! local variables
    integer :: i,j,k, iVar
    integer :: nOutVar    ! nb. of total output variables 
    integer :: nDimOut    ! nb. of output dimensions
    real :: time_dim
    real :: pressure_dim  ! dimensional pressure

    ! mean wave length over cells
    real, dimension(nx,ny,nz,2) :: lambda
    integer, dimension(nx,ny,nz) :: nRayPerCell
    integer :: iRay
    real :: kk, mm

    ! file name
    integer :: length
    character (len = 20)  :: fmt, form
    character (len = 40)  :: cpuTimeChar, gridName, auxChar
!    character (len = 200) :: title, dataFile   ! modified by Junhong Wei (20161110)
    character (len = 200) :: title, dataFile, zoneChar   ! modified by Junhong Wei (20161110)
    character (len = 1000) :: variablesChar

    ! CPU Time
    integer    :: days, hours, mins, secs
    real       :: timeVar

    ! optional variables
    real :: uR,uL, vF,vB, wU,wD
    real :: pStratU, pStratD, divPu, divU
    real :: p, pEx, pBar

    ! ghostCell output
    integer :: zoneSizeX, zoneSizeY, zoneSizeZ
    integer :: i0,i1, j0,j1, k0,k1

    ! greshoVortex output
    real :: delX, delY, r

    ! hotBubble output
    real :: rho, theta_dim

    ! buoyancy
    real :: b, b_dim, theta

    ! stratification perturbation
    real :: rho_u, rho_d
    real :: theta0_u, theta0_d
    real :: theta_u, theta_d
    real :: dtheta_u, dtheta_d
    real :: Pstrat_u, Pstrat_d
    real :: db_dz

    ! ray scatter plot
    real :: maxWaveAct

    ! WKB amplitudes
    real :: A,u,w,th
    integer :: harmonic

    ! complex output
!    logical, parameter :: outComplex = .true.   -> now in input.f90 file
    real :: u_real, u_imag
    real :: w_real, w_imag
    real :: b_real, b_imag
    real :: p_real, p_imag

    ! MPI stuff   ! modified by Junhong Wei (20161110)
    integer :: i00, j00   ! modified by Junhong Wei (20161110)

    time_dim = time * tRef

    if( master ) then ! modified by Junhong Wei (20161110)
    print*,""
    print*," Output into File "
    print*,""
    write(*,fmt="(a25,i15)") " at time step = ", iTime
    write(*,fmt="(a25,f15.1,a8)") " at physical time = ", time_dim, " seconds"
    end if ! modified by Junhong Wei (20161110)
    
    !------------------------------
    !   prepare output file
    !------------------------------
    
    
    !-------------------------
    !   generate file name
    !-------------------------
    dataFile = ""                                         ! create empty string

    write(unit=auxChar, fmt="(i2.2)") iParam              ! param study tag
    if( parameterStudy ) then 
       dataFile = trim(dataFile)//trim(auxChar)//"-"//trim(paramName)
    else
       dataFile = trim(dataFile)//trim(testCase)          ! test case 
    end if

    ! grid 
    write(unit=form, fmt="(a,i1,a,i1,a,i1,a)") &
         & "(i",int(log10(real(nx)))+1, &
         & ",a,i",int(log10(real(ny)))+1, &
         & ",a,i",int(log10(real(nz)))+1, &
         & ")" ! format for nx x ny x nz

    ! generate text string for grid size
    write(unit=gridName, fmt=form) nx,"x",ny,"x",nz
    dataFile = trim(dataFile)//"-"//trim(gridName)      ! grid size

    if( restart ) dataFile = trim(dataFile)//"-rst"     ! restart tag

    write(unit=auxChar, fmt="(i6.6)") iOut              ! output counter
    dataFile = trim(dataFile)//"-"//trim(auxChar)


    ! modified by Junhong Wei (20161111) *** starting line ***
    ! create MPI-proc related file name
    write(unit=auxChar, fmt="(i3.3)") rank              ! output counter
    dataFile = trim(dataFile)//"_pid"//trim(auxChar)
    ! modified by Junhong Wei (20161111) *** finishing line ***


    dataFile = trim(dataFile)//".dat"                    ! file ending

    !    write(unit=dataFile, fmt="(a)") trim(dataFile)     ! dataFileName -> dataFile

    open( unit=40, file=dataFile, action="write", &
         form="formatted", status="replace")


    !------------------------------
    !   generate tecplot header
    !------------------------------

    ! calc cpu-time in days/hours/minutes/seconds
    timeVar = cpuTime   
    days = floor(timeVar / 86400.0)
    timeVar = timeVar - 86400.0 * days
    hours = floor(timeVar / 3600.0)
    timeVar = timeVar - 3600.0 * hours
    mins = floor(timeVar / 60.0)
    timeVar = timeVar - 60.0 * mins
    secs = int(timeVar)

    write(unit=cpuTimeChar, fmt="(i2,a,2(i2.2,a),i2.2)") &
       & days, " ", hours, ":",mins,":",secs
    
!    write(*,fmt="(a25,a25)")  "CPU time = ", cpuTimeChar
    if( master ) write(*,fmt="(a25,a25)")  "CPU time = ", cpuTimeChar
    
    write(unit=title, fmt="(a,i10.10,a,es10.3,3a)") "TITLE = ""time-step: ", &
         & iTime, " | time:", time_dim, "s | cpu-time:", trim(cpuTimeChar), """"



    ! obtain line for variables
    if( dimOut(1) .and. dimOut(2) .and. dimOut(3) ) then    ! xyz plot
       variablesChar = "Variables =  ""x [m]"" ""y [m]"" ""z [m]"" "
       
    else if( dimOut(1) .and. .not.dimOut(2) .and. dimOut(3) ) then ! xz plot
       variablesChar = "Variables =  ""x [m]"" ""z [m]"" "

    else if( .not.dimOut(1) .and. .not.dimOut(2) .and. dimOut(3) ) then ! z-plot
       variablesChar = "Variables =  ""x [m]"" ""z [m]"" "

    else
       stop"tec360: specify case for line of variables"
    end if

    ! NEW: always write all space variables (necessary for ray tracer)
    if( dimOut(1) .and. dimOut(2) .and. dimOut(3) ) then    ! xyz plot
       nDimOut = 3
       
    else if( dimOut(1) .and. .not.dimOut(2) .and. dimOut(3) ) then ! xz plot
       nDimOut = 2

    else if( .not.dimOut(1) .and. .not.dimOut(2) .and. dimOut(3) ) then ! z-plot
       nDimOut = 2

    else
       stop"tec360: specify case for mDimOut"
    end if

    ! nb of variables initialised with nb of spatial vars
    nOutVar = nDimOut

    !-------------------
    ! primary variables
    !-------------------
    do iVar = 1,nVar
       if (varOut(iVar)==1) then
          nOutVar = nOutVar + 1
          if (iVar==1) then 
             if( rhoOffset ) then
                variablesChar = trim(variablesChar)//&
                     &" ""<greek>r</greek>' [kg/m<sup>3</sup>]""  "
             else
                variablesChar = trim(variablesChar)//&
                     &" ""<greek>r</greek> [kg/m<sup>3</sup>]""  "
             end if
          end if
          if (iVar==2) variablesChar = trim(variablesChar)//" ""u [m/s]""  "
          if (iVar==3) variablesChar = trim(variablesChar)//" ""v [m/s]""  "
          if (iVar==4) variablesChar = trim(variablesChar)//" ""w [m/s]""  "
          if (iVar==5) variablesChar = trim(variablesChar)//" ""<greek>p</greek>' ""  "
          if (iVar==6) variablesChar = trim(variablesChar)//" ""<greek>q</greek>' [K]""  " 
       end if
    end do ! iVar

    !--------------------
    ! optional variables
    !--------------------
    do iVar = 1,size(optVarOut)
       if (optVarOut(iVar)==1) then
          nOutVar = nOutVar + 1
          if (iVar==1) variablesChar = trim(variablesChar)//&
               & " ""p' [kPa]""  "
          if (iVar==2) then
             variablesChar = trim(variablesChar)//&
                  & " ""b [m/s<sup>2</sup>]""  "
          end if
          if (iVar==3) variablesChar = trim(variablesChar)//" ""p<sub>bg</sub> [kPa] ""  "
          if (iVar==4) variablesChar = trim(variablesChar)//&
               & " ""<greek>r</greek><sub>bg</sub> [kg/m<sup>3</sup>] ""  "
          if (iVar==5) variablesChar = trim(variablesChar)//" ""div(Pu) [Pa/s]""  "
          if (iVar==6) variablesChar = trim(variablesChar)//&
               & " ""b<sub>z</sub>/N<sup>2</sup> ""  "
       end if
    end do ! iVar

    !------------------------
    ! WKB optional variables
    !------------------------
    do iVar = 1,size(wkbVarOut)
       if (wkbVarOut(iVar)==1) then
          nOutVar = nOutVar + 1

          ! wave action and mean flow
          if (iVar==1) variablesChar = trim(variablesChar)//&
               & " ""A [Js/m<sup>3</sup>]""  "
          if (iVar==2) variablesChar = trim(variablesChar)//&
               & " ""u<sub>0</sub><sup>(0)</sup>' [m/s]""  "
          if (iVar==3) variablesChar = trim(variablesChar)//&
               & " ""b<sub>0</sub><sup>(2)</sup> [m/s<sup>2</sup>]""  "
          if (iVar==4) variablesChar = trim(variablesChar)//&
               & " ""<greek>q</greek><sub>0</sub>&
               & <greek>p</greek><sub>0</sub><sup>(2)</sup> ""  "

          ! first harmonic
          if (iVar==5) variablesChar = trim(variablesChar)//&
               & " ""u<sub>1</sub><sup>(0)</sup>' [m/s]""  "
          if (iVar==6) variablesChar = trim(variablesChar)//&
               & " ""w<sub>1</sub><sup>(0)</sup>' [m/s]""  "
          if (iVar==7) variablesChar = trim(variablesChar)//&
               & " ""b<sub>1</sub><sup>(1)</sup> [m/s<sup>2</sup>]""  "
          if (iVar==8) variablesChar = trim(variablesChar)//&
               & " ""<greek>q</greek><sub>0</sub>&
               & <greek>p</greek><sub>1</sub><sup>(2)</sup> ""  "

          ! second harmonic
          if (iVar==9) variablesChar = trim(variablesChar)//&
               & " ""u<sub>2</sub><sup>(1)</sup>' [m/s]""  "
          if (iVar==10) variablesChar = trim(variablesChar)//&
               & " ""w<sub>2</sub><sup>(1)</sup>' [m/s]""  "
          if (iVar==11) variablesChar = trim(variablesChar)//&
               & " ""b<sub>2</sub><sup>(2)</sup> [m/s<sup>2</sup>]""  "
          if (iVar==12) variablesChar = trim(variablesChar)//&
               & " ""<greek>q</greek><sub>0</sub>&
               & <greek>p</greek><sub>2</sub><sup>(3)</sup> ""  "

          ! ray information: vertical wave length
          if (iVar==13) variablesChar = trim(variablesChar)//&
               & " ""<greek>l</greek><sub>z</sub> [m]""  "

       end if

    end do ! iVar


    write(unit=40,fmt="(a)") title, trim(variablesChar)

    if( showGhostCellsX ) then
       zoneSizeX = nx+2*nbx+2
    else
       zoneSizeX = nx+1
    end if

    if( showGhostCellsY ) then
       zoneSizeY = ny+2*nby+2
    else
       zoneSizeY = ny+1
    end if

    if( showGhostCellsZ) then
       zoneSizeZ = nz+2*nbz+2
    else
       zoneSizeZ = nz+1
    end if

    !-----------------------------
    !      Zone description
    !-----------------------------
    write(unit=40,fmt="(a,i5)") " Zone "

    ! modified by Junhong Wei (20161110) *** starting line ***
    
!    ! name of zone
!    write(unit=40,fmt="(a)") " T = ""fields""  "

    ! name of zone
    write(unit=auxChar, fmt="(i3.3)") rank              ! output counter
    zoneChar = " T = ""zone_" // trim(auxChar) // """"
    
    write(unit=40,fmt="(a)") trim(zoneChar)

    ! modified by Junhong Wei (20161110) *** finishing line ***
    
    ! size of zone
    if( dimOut(1) .and. dimOut(2) .and. dimOut(3) ) then    ! xyz-plot
       write(unit=40,fmt="(3(a,i5))") " i =",zoneSizeX, &
            &" j =", zoneSizeY, " k =", zoneSizeZ
       
    else if( dimOut(1) .and. .not.dimOut(2) .and. dimOut(3) ) then ! xz-plot
       write(unit=40,fmt="(2(a,i5))") " i =",zoneSizeX, " k =", zoneSizeZ
       
    else if( .not.dimOut(1) .and. .not.dimOut(2) .and. dimOut(3) ) then ! z-plot
       write(unit=40,fmt="(a,i5)") " i =", zoneSizeZ                   ! i must be present

    else
       stop"tec360: this type of dimension specification is not allowed."
    end if

    write(unit=40,fmt="(a)") " Datapacking = Block "

    ! variable location
    write(unit=40,fmt="(a)",advance="no") "VarLocation = ( "
    do i = 1,nDimOut
       write(unit=40,fmt="(a)",advance="no") " Nodal,"
    end do
    
    do i = 1,nOutVar-nDimOut-1
       write(unit=40,fmt="(a)",advance="no") " CellCentered,"
    end do
    write(unit=40,fmt="(a)",advance="no") " CellCentered"   ! last entry without comma
    write(unit=40,fmt="(a)",advance="yes") ")"

    ! double precision for tecplot
    write(unit=40,fmt="(a)",advance="no") "DT = ( "
    do i = 1,nOutVar
       write(unit=40,fmt="(a)",advance="no") "Single "
    end do
    write(unit=40,fmt="(a)",advance="yes") ")"


    ! solution time in seconds, minutes or hours
    if (solutionTime) then
       if ( solutionTimeUnit == "s" ) then
          write(unit=40,fmt="(a,/,es19.10)") "solutiontime = ", time_dim
       else if ( solutionTimeUnit == "min" ) then
          write(unit=40,fmt="(a,/,es19.10)") "solutiontime = ", time_dim/60.0
       else if ( solutionTimeUnit == "h" ) then
          write(unit=40,fmt="(a,/,es19.10)") "solutiontime = ", time_dim/3600.0
       end if
    end if


    !---------------------------------------
    !   write xy-coordinates of nodes       
    !---------------------------------------

! modified by Junhong Wei (20161110) *** starting line ***
    
    ! get MPI local start index 
    i00 = is + nbx - 1   ! 0 index, replace i -> i + i0 in x and y fields
    j00 = js + nby - 1    

! modified by Junhong Wei (20161110) *** finishing line ***    
    
    if (dimOut(1)) then
       if( showGhostCellsX ) then
          i0 = -nbx-1;  i1 = nx+nbx
       else
          i0 = 0;  i1 = nx
       end if
    else
       i0 = 1; i1 = 1  
    end if
    if (dimOut(2)) then
       if( showGhostCellsY ) then
          j0 = -nby-1;  j1 = ny+nby
       else
          j0 = 0;  j1 = ny
       end if
    else
       j0 = 1; j1 = 1 
    end if
    if (dimOut(3)) then
       if( showGhostCellsZ ) then
          k0 = -nbz-1;  k1 = nz+nbz
       else
          k0 = 0;  k1 = nz
       end if
    else
       k0 = 1; k1 = 1
    end if
    
    
    ! always write x-coordinate
    do k = k0, k1
       do j = j0, j1
!          do i = i0, i1   ! modified by Junhong Wei (20161110)
          do i = i0 + i00, i1 + i00  ! MPI shift for local x coords   ! modified by Junhong Wei (20161110)
             if (i == -nbx-1) then
!                write(unit=40,fmt="(4es17.6)") (x(-nbx)-dx/2) * lRef
                write(unit=40,fmt="(4es17.6)") (x(-nbx+i00)-dx/2) * lRef   ! modified by Junhong Wei (20161110)
             else
                write(unit=40,fmt="(4es17.6)") (x(i)+dx/2) * lRef
             end if
          end do
       end do
    end do
    write(unit=40,fmt="(a)" ) ""


    ! write y-coordinate in case
    if (dimOut(2)) then
       do k = k0, k1
!          do j = j0, j1   ! modified by Junhong Wei (20161110)
          do j = j0 + j00, j1 + j00     ! MPI shift for local y coords   ! modified by Junhong Wei (20161110)
          do i = i0, i1
                if (j == -nby-1) then
!                   write(unit=40,fmt="(4es17.6)") (y(-nby)-dy/2) * lRef   ! modified by Junhong Wei (20161110)
                   write(unit=40,fmt="(4es17.6)") (y(-nby+j00)-dy/2) * lRef   ! modified by Junhong Wei (20161110)
                else
                   write(unit=40,fmt="(4es17.6)") (y(j)+dy/2) * lRef
                end if
             end do
          end do
       end do
       write(unit=40,fmt="(a)" ) ""
    end if

    ! always write z-coordinate
    do k = k0, k1
       do j = j0, j1
          do i = i0, i1
             if (k == -nbz-1) then
                write(unit=40,fmt="(4es17.6)") (z(-nbz)-dz/2) * lRef
             else
                write(unit=40,fmt="(4es17.6)") (z(k)+dz/2) * lRef
             end if
          end do
       end do
    end do
    write(unit=40,fmt="(a)" ) ""


    !-------------------------------------------
    !         write primary variables
    !-------------------------------------------

    if (dimOut(1)) then
       i0 = 1; i1 = nx
    else
       i0 = 1; i1 = 1 
    end if
    if (dimOut(2)) then
          j0 = 1; j1 = ny
    else
       j0 = 1; j1 = 1 
    end if
    if (dimOut(3)) then
          k0 = 1; k1 = nz
    else
       k0 = 1; k1 = 1  
    end if

    do iVar = 1, nVar
       if (varOut(iVar)==1) then
          do k = k0, k1
             do j = j0, j1
                do i = i0, i1

                   !---------------------------------------
                   !       dimensionalising and output
                   !---------------------------------------
                   select case (iVar) 

                   case(1) ! density

                      if( fluctuationMode) then
                         
                         if( rhoOffset ) then
                            write(unit=40,fmt="(4es17.6)" ) var(i,j,k,iVar) * rhoRef
                         else
                            write(unit=40,fmt="(4es17.6)" ) &
                                 &(var(i,j,k,iVar) + rhoStrat(k)) * rhoRef
                         end if
                         
                      else
                         
                         if( rhoOffset ) then
                            write(unit=40,fmt="(4es17.6)" ) &
                                 &(var(i,j,k,iVar) - rhoStrat(k)) * rhoRef
                         else
                            write(unit=40,fmt="(4es17.6)" ) var(i,j,k,iVar) * rhoRef
                         end if
                         
                      end if
                      
                      ! average velocities to cell center
                      
                   case(2) ! u velocity
                      write(unit=40,fmt="(4es17.6)" ) &
                           & (0.5*(var(i,j,k,iVar)+var(i-1,j,k,iVar)) &
                           & - offset(iVar)) * uRef

                   case(3) ! v velocity
                      write(unit=40,fmt="(4es17.6)" ) &
                           & (0.5*(var(i,j,k,iVar)+var(i,j-1,k,iVar)) &
                           & - offset(iVar)) * uRef
                      
                   case(4) ! w velocity
                      write(unit=40,fmt="(4es17.6)" ) &
                           & (0.5*(var(i,j,k,iVar)+var(i,j,k-1,iVar)) &
                           & - offset(iVar)) * uRef
                      
                   case(5) ! Exner function pi' (deviation from background)
                      write(unit=40,fmt="(4es17.6)" ) var(i,j,k,iVar)

                   case(6) ! potential temperature theta' 
                           ! (deviation from background, Boussinesq)

                      select case( model ) 

                      case( "pseudo_incompressible")

                         if( fluctuationMode ) then
                            rho = var(i,j,k,1) + rhoStrat(k)
                         else
                            rho = var(i,j,k,1)
                         end if
                         
                         if ( referenceQuantities == "SI" ) then
                            theta_dim = 1.0/Rsp * Pstrat(k) / rho * thetaRef
                         else
                            theta_dim = Pstrat(k) / rho * thetaRef
                         end if
                         
                         if( thetaOffset ) theta_dim = theta_dim - thetaStrat(k)*thetaRef
                         write(unit=40,fmt="(4es17.6)" ) theta_dim
                         
                         
                      case( "Boussinesq" )
                         
                         write(unit=40,fmt="(4es17.6)" ) var(i,j,k,iVar)*thetaRef

                      case( "WKB" )
                         
                         
                      case default
                         stop"tec360: unknown model"                         
                      end select ! model
                      
                   case default
                      stop"tec360: unkown iVar"
                   end select ! iVar

                end do
             end do
          end do
       end if
       write(unit=40,fmt="(a)" ) ""
    end do


    !------------------------------------------------------------
    !                   Calc optional variables
    !------------------------------------------------------------

    if (dimOut(1)) then
       i0 = 1; i1 = nx
    else
       i0 = 1; i1 = 1    
    end if
    if (dimOut(2)) then
       j0 = 1; j1 = ny
    else
       j0 = 1; j1 = 1    
    end if
    if (dimOut(3)) then
       k0 = 1; k1 = nz
    else
       k0 = 1; k1 = 1   
    end if


    !----------------------------------
    !  1) pressure fluctuation in kPa
    !----------------------------------
    
    if (optVarOut(1)==1) then
       
       do k = k0, k1
          do j  = j0, j1
             do i = i0, i1
                
                if( fluctuationMode ) then
                   rho = var(i,j,k,1) + rhoStrat(k)
                else
                   rho = var(i,j,k,1)
                end if

                pBar = p0**(1.0-gamma) * Pstrat(k)**gamma
                pEx = var(i,j,k,5)      ! Exner function (deviation from background)

                p = p0*(pEx + (pBar/p0)**kappa )**kappaInv

                optVar(i,j,k) = (p-pBar) * pRef / 1.0e3
             end do
          end do
       end do
       
       ! write to file
       do k = k0, k1
          do j = j0, j1
             do i = i0, i1
                write(unit=40,fmt="(4es17.6)" ) optVar(i,j,k)
             end do
          end do
       end do
       write(unit=40,fmt="(a)") ""
    end if
    

    !-------------------------------
    !    2) buoyancy
    !-------------------------------
    
    if( optVarOut(2) == 1 ) then
       
       do k = k0, k1
          do j  = j0, j1
             do i = i0, i1

                select case( model ) 
                   
                case( "pseudo_incompressible" )

                   if( fluctuationMode) then
                      rho = var(i,j,k,1) + rhoStrat(k)
                   else 
                      rho = var(i,j,k,1)
                   end if
                   theta = Pstrat(k) / rho
                   b = FrInv2 * (theta-thetaStrat(k))/thetaStrat(k)
                   b_dim = b*lRef/tRef**2
                   
                case( "Boussinesq" )
                   theta = var(i,j,k,6)
                   b = FrInv2 * theta / theta00
                   b_dim = b*lRef/tRef**2
                case default
                   stop"tec360: unknown case model"
                end select
                
                optVar(i,j,k) = b_dim
             end do
          end do
       end do
       
       ! write to file
       do k = k0, k1
          do j = j0, j1
             do i = i0, i1
                write(unit=40,fmt="(4es17.6)" ) optVar(i,j,k)
             end do
          end do
       end do
       write(unit=40,fmt="(a)") ""
    end if


    !--------------------------------
    !  3) background pressure in kPa
    !--------------------------------

    if( optVarOut(3) == 1 ) then
       
       do k = k0, k1
          do j  = j0, j1
             do i = i0, i1
                optVar(i,j,k) = p0**(1.0-gamma)*PStrat(k)**gamma * pRef / 1.0e3
             end do
          end do
       end do
       
       ! write to file
       do k = k0, k1
          do j = j0, j1
             do i = i0, i1
                write(unit=40,fmt="(4es17.6)" ) optVar(i,j,k)
             end do
          end do
       end do
       write(unit=40,fmt="(a)") ""
    end if


    !----------------------------------
    !  4) background density in kg/m^3
    !----------------------------------

    if( optVarOut(4) == 1 ) then
       
       do k = k0, k1
          do j  = j0, j1
             do i = i0, i1
                optVar(i,j,k) = rhoStrat(k) * rhoRef
             end do
          end do
       end do
       
       ! write to file
       do k = k0, k1
          do j = j0, j1
             do i = i0, i1
                write(unit=40,fmt="(4es17.6)" ) optVar(i,j,k)
             end do
          end do
       end do
       write(unit=40,fmt="(a)") ""
    end if

    
    !-------------------------------
    !     5) divergence of Pu
    !-------------------------------
    
    
    if( optVarOut(5) == 1 ) then

       do k = k0, k1
          do j  = j0, j1
             do i = i0, i1
                uR = var(i,j,k,2); uL = var(i-1,j,k,2)
                vF = var(i,j,k,3); vB = var(i,j-1,k,3)
                wU = var(i,j,k,4); wD = var(i,j,k-1,4)

                PStratU = 0.5*(PStrat(k+1) + PStrat(k  ))
                PStratD = 0.5*(PStrat(k  ) + PStrat(k-1))

                divPu = PStrat(k) * (uR-uL)/dx + &
                     & PStrat(k) * (vF-vB)/dy + &
                     & (PStratU*wU - PStratD*wD)/dz

                optVar(i,j,k) = divPu * pRef * uRef / lRef 
             end do
          end do
       end do

       ! write to file
       do k = k0, k1
          do j = j0, j1
             do i = i0, i1
                write(unit=40,fmt="(4es17.6)" ) optVar(i,j,k)
             end do
          end do
       end do
       write(unit=40,fmt="(a)") ""
    end if

    !------------------------------------------
    ! 6)  Stratifiaction perturbation db/dz/N2
    !------------------------------------------
    
    
    if( optVarOut(6) == 1 ) then

       do k = k0, k1
          do j  = j0, j1
             do i = i0, i1
                
                rho_u = var(i,j,k+1,1)
                rho_d = var(i,j,k-1,1)
                theta0_u = thetaStrat(k+1)
                theta0_d = thetaStrat(k-1)
                PStrat_u = PStrat(k+1)
                Pstrat_d = Pstrat(k-1)

                
                ! fluctuation mode
                if( fluctuationMode) then
                   rho_u = rho_u + rhoStrat(k+1)
                   rho_d = rho_d + rhoStrat(k-1)
                end if
                
                if ( referenceQuantities == "SI" ) then
                   theta_u = 1.0/Rsp * Pstrat_u / rho_u
                   theta_d = 1.0/Rsp * Pstrat_d / rho_d
                else
                   theta_u = Pstrat_u / rho_u
                   theta_d = Pstrat_d / rho_d
                end if
                
                dtheta_u = theta_u - theta0_u
                dtheta_d = theta_d - theta0_d
                
                db_dz = FrInv2*0.5/dz * (dtheta_u/theta0_u - dtheta_d/theta0_d)
                optVar(i,j,k) = db_dz / N2

                ! boundary treatment: prob
                if( k==1 .or. k==2 .or. k==nz .or. k==nz-1) then
                   optVar(i,j,k) = 0.0
                end if
                
             end do
          end do
       end do

       ! write to file
       do k = k0, k1
          do j = j0, j1
             do i = i0, i1
                write(unit=40,fmt="(4es17.6)" ) optVar(i,j,k)
             end do
          end do
       end do
       write(unit=40,fmt="(a)") ""
    end if
    

    
    
    !-----------------------------------------------------------
    !                 Calc WKB  variables
    !-----------------------------------------------------------
    
    if( raytracer ) then



       !-----------------------------------------------
       !             Mean flow & Wave action
       !-----------------------------------------------
       harmonic = 0       

       !-------------------------------
       ! 1) Wave action 
       !-------------------------------
       iVar = 1
       if( wkbVarOut(iVar) == 1 ) then

          do k = k0, k1
             do j  = j0, j1
                do i = i0, i1

                   A = waveAct(i,j,k)*rhoRef * uRef**2
                   write(unit=40,fmt="(4es17.6)" ) A

                end do
             end do
          end do
          write(unit=40,fmt="(a)") ""
       end if

       !-------------------------------
       ! 2) Zonal velocity u00
       !-------------------------------
       iVar = 2

       if( wkbVarOut(iVar) == 1 ) then

          do k = k0, k1
             do j  = j0, j1
                do i = i0, i1
                   u = real(Psi(i,j,k,1,harmonic)) * uRef
                   write(unit=40,fmt="(4es17.6)" ) u
                end do
             end do
          end do
          write(unit=40,fmt="(a)") ""
       end if


       !--------------------------------
       ! 3) Potential temperature th01
       !    currently: term containing th01 and pi02
       !--------------------------------

       iVar = 3
       if( wkbVarOut(iVar) == 1 ) then

          do k = k0, k1
             do j  = j0, j1
                do i = i0, i1
                   b = real(Psi(i,j,k,3,harmonic)) * lRef/tRef**2
                   write(unit=40,fmt="(4es17.6)" ) b
                end do
             end do
          end do
          write(unit=40,fmt="(a)") ""
       end if


       !-------------------------------
       ! 4)  Exner pressure pi02
       !-------------------------------

       iVar = 4
       if( wkbVarOut(iVar) == 1 ) then

          do k = k0, k1
             do j  = j0, j1
                do i = i0, i1
                   p = real(Psi(i,j,k,4,harmonic))
                   write(unit=40,fmt="(4es17.6)" ) p
                end do
             end do
          end do
          write(unit=40,fmt="(a)") ""
       end if




       !----------------------------------------------
       !                First harmonic
       !-----------------------------------------------
       harmonic = 1


       !-------------------------------
       ! 5) Zonal velocity u10
       !-------------------------------
       iVar = 5

       if( wkbVarOut(iVar) == 1 ) then

          if( outComplex ) then   ! complex data output for matlab analysis
             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      u_real = real(Psi(i,j,k,1,harmonic)) * uRef
                      write(unit=40,fmt="(4es17.6)" ) u_real
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""

             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      u_imag = aimag(Psi(i,j,k,1,harmonic)) * uRef
                      write(unit=40,fmt="(4es17.6)" ) u_imag
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""


          else  ! output of absolute value for tecplot 

             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      u = cabs(Psi(i,j,k,1,harmonic)) * uRef
                      write(unit=40,fmt="(4es17.6)" ) u
                   end do
                end do
             end do
          end if

          write(unit=40,fmt="(a)") ""
       end if


       !-------------------------------
       ! 6) Vertical velocity w10
       !-------------------------------

       iVar = 6
       if( wkbVarOut(iVar) == 1 ) then

          
          if( outComplex ) then
             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      w_real = real(Psi(i,j,k,2,harmonic)) * uRef
                      write(unit=40,fmt="(4es17.6)" ) w_real
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""

             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      w_imag = aimag(Psi(i,j,k,2,harmonic)) * uRef
                      write(unit=40,fmt="(4es17.6)" ) w_imag
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""


          else  ! output of absolute value

             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      w = cabs(Psi(i,j,k,2,harmonic)) * uRef
                      write(unit=40,fmt="(4es17.6)" ) w
                   end do
                end do
             end do

          end if
          
          write(unit=40,fmt="(a)") ""

       end if


       !--------------------------------
       ! 7)  Buoyancy b11
       !--------------------------------

       iVar = 7
       if( wkbVarOut(iVar) == 1 ) then


          if( outComplex ) then
             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      b_real = real(Psi(i,j,k,3,harmonic)) * lRef/tRef**2
                      write(unit=40,fmt="(4es17.6)" ) b_real
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""

             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      b_imag = aimag(Psi(i,j,k,3,harmonic)) * lRef/tRef**2
                      write(unit=40,fmt="(4es17.6)" ) b_imag
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""


          else  ! output of absolute value

             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      b = cabs(Psi(i,j,k,3,harmonic)) * lRef/tRef**2
                      write(unit=40,fmt="(4es17.6)" ) b
                   end do
                end do
             end do

          end if

          write(unit=40,fmt="(a)") ""
       end if


       !-------------------------------
       ! 8)  Exner pressure
       !-------------------------------
       
       iVar = 8
       if( wkbVarOut(iVar) == 1 ) then

          if( outComplex ) then
             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      p_real = real(Psi(i,j,k,4,harmonic)) 
                      write(unit=40,fmt="(4es17.6)" ) p_real
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""

             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      p_imag = aimag(Psi(i,j,k,4,harmonic))
                      write(unit=40,fmt="(4es17.6)" ) p_imag
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""


          else  ! output of absolute value

             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      p = cabs(Psi(i,j,k,4,harmonic))
                      write(unit=40,fmt="(4es17.6)" ) p
                   end do
                end do
             end do

          end if

          write(unit=40,fmt="(a)") ""
       end if



       !----------------------------------------------
       !                Second harmonic
       !-----------------------------------------------
       harmonic = 2


       !-------------------------------
       ! 9) Zonal velocity u21
       !-------------------------------
       iVar = 9

       if( wkbVarOut(iVar) == 1 ) then

          if( outComplex ) then
             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      u_real = real(Psi(i,j,k,1,harmonic)) * uRef
                      write(unit=40,fmt="(4es17.6)" ) u_real
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""

             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      u_imag = aimag(Psi(i,j,k,1,harmonic)) * uRef
                      write(unit=40,fmt="(4es17.6)" ) u_imag
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""


          else  ! output of absolute value
             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      u = cabs(Psi(i,j,k,1,harmonic)) * uRef
                      write(unit=40,fmt="(4es17.6)" ) u
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""
          end if ! complex

       end if



       !-------------------------------
       ! 10) Vertical velocity w21
       !-------------------------------
       
       iVar = 10
       if( wkbVarOut(iVar) == 1 ) then

          if( outComplex ) then
             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      w_real = real(Psi(i,j,k,2,harmonic)) * uRef
                      write(unit=40,fmt="(4es17.6)" ) w_real
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""

             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      w_imag = aimag(Psi(i,j,k,2,harmonic)) * uRef
                      write(unit=40,fmt="(4es17.6)" ) w_imag
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""


          else  ! output of absolute value
             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      w = cabs(Psi(i,j,k,2,harmonic)) * uRef
                      write(unit=40,fmt="(4es17.6)" ) w
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""
          end if ! complex

       end if


       !--------------------------------
       ! 11)      Buoyancy b22
       !--------------------------------
       
       iVar = 11
       if( wkbVarOut(iVar) == 1 ) then

          if( outComplex ) then
             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      b_real = real(Psi(i,j,k,3,harmonic)) * lRef/tRef**2
                      write(unit=40,fmt="(4es17.6)" ) b_real
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""

             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      b_imag = aimag(Psi(i,j,k,3,harmonic)) * lRef/tRef**2
                      write(unit=40,fmt="(4es17.6)" ) b_imag
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""


          else  ! output of absolute value
             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      b = cabs(Psi(i,j,k,3,harmonic)) * lRef/tRef**2
                      write(unit=40,fmt="(4es17.6)" ) b
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""
          end if ! complex

       end if


       !-------------------------------
       ! 12)  Exner pressure
       !-------------------------------
       
       iVar = 12
       if( wkbVarOut(iVar) == 1 ) then

          if( outComplex ) then
             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      p_real = real(Psi(i,j,k,4,harmonic)) 
                      write(unit=40,fmt="(4es17.6)" ) p_real
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""

             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      p_imag = aimag(Psi(i,j,k,4,harmonic))
                      write(unit=40,fmt="(4es17.6)" ) p_imag
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""


          else  ! output of absolute value
             do k = k0, k1
                do j  = j0, j1
                   do i = i0, i1
                      p = cabs(Psi(i,j,k,4,harmonic))
                      write(unit=40,fmt="(4es17.6)" ) p
                   end do
                end do
             end do
             write(unit=40,fmt="(a)") ""
          end if ! complex

       end if



       !----------------------------------------------
       !              Ray quantities
       !-----------------------------------------------



       !-------------------------------
       ! 13)  lambda_z
       !-------------------------------
       iVar = 13
       
       if( wkbVarOut(iVar) == 1 ) then

          ! Calculation of wave lengths
          
          lambda(:,:,:,:) = 0.0
          nRayPerCell(:,:,:) = 0

          ! ray loop
          do iRay = 1, nRay
             i = ray(iRay)%iCell
             j = 1
             k = ray(iRay)%kCell
             kk = ray(iRay)%k
             mm = ray(iRay)%m
             
             lambda(i,j,k,2) = lambda(i,j,k,2) + 2.0*pi/mm
             nRayPerCell(i,j,k) = nRayPerCell(i,j,k) + 1
          end do
          
          ! cell loop
          do k = 1,nz
             j = 1
             do i = 1,nx
                lambda(i,j,k,2) = lambda(i,j,k,2) / nRayPerCell(i,j,k)
             end do
          end do

          ! output lambda_z into file
          do k = k0, k1
             do j  = j0, j1
                do i = i0, i1
                   write(unit=40,fmt="(4es17.6)" ) lambda(i,j,k,2)*lRef
                end do
             end do
          end do
          write(unit=40,fmt="(a)") ""

       end if



    end if ! raytracer
    


    !------------------------------------------
    !          plot rays as scatter plot
    !------------------------------------------

    
    if (rayTracer) then
       !-----------------------------
       !      Zone description
       !-----------------------------
       write(unit=40,fmt="(a,i5)") " Zone "

       ! name of zone
       write(unit=40,fmt="(a)") " T = ""rays""  "

       ! size of zone: nb. of rays
       write(unit=40,fmt="(a,i10)") " i =", nRay 

       ! data packing
       write(unit=40,fmt="(a)") " Datapacking = Block "

       ! variable location
!       write(unit=40,fmt="(a)",advance="no") "VarLocation = ( Nodal, Nodal, Nodal, "

       write(unit=40,fmt="(a)",advance="no") "VarLocation = ( "
       do i = 1,nDimOut
          write(unit=40,fmt="(a)",advance="no") " Nodal,"
       end do

       !write(unit=40,fmt="(a)",advance="no") "VarLocation = ( "
       !if (dimOut(1)) write(unit=40,fmt="(a)",advance="no") " Nodal,"
       !if (dimOut(2)) write(unit=40,fmt="(a)",advance="no") " Nodal,"       
       !if (dimOut(3)) write(unit=40,fmt="(a)",advance="no") " Nodal,"

       do i = 1,nOutVar-nDimOut-1
          write(unit=40,fmt="(a)",advance="no") " CellCentered,"
       end do
       write(unit=40,fmt="(a)",advance="no") " CellCentered"   ! last entry without comma
       write(unit=40,fmt="(a)",advance="yes") ")"

       ! double precision declaration for tecplot
       write(unit=40,fmt="(a)",advance="no") "DT = ( "
       do i = 1,nOutVar
          write(unit=40,fmt="(a)",advance="no") "Single "
       end do
       write(unit=40,fmt="(a)",advance="yes") ")"

       ! list of passive variables (not given)
       if (nOutVar < 10) then
          write(unit=40,fmt="(a,i1,a1,i1,a)") " passivevarlist = [",nDimOut+1,"-", nOutVar, "]"
       else
          write(unit=40,fmt="(a,i1,a1,i2,a)") " passivevarlist = [",nDimOut+1,"-", nOutVar, "]"
       end if

       ! solution time in seconds, minutes or hours
       if (solutionTime) then
          if ( solutionTimeUnit == "s" ) then
             write(unit=40,fmt="(a,/,es19.10)") " solutiontime = ", time_dim
          else if ( solutionTimeUnit == "min" ) then
             write(unit=40,fmt="(a,/,es19.10)") " solutiontime = ", time_dim/60.0
          else if ( solutionTimeUnit == "h" ) then
             write(unit=40,fmt="(a,/,es19.10)") " solutiontime = ", time_dim/3600.0
          end if
       end if

       !----------------------------------
       !      x, y and z coordinates of rays
       !----------------------------------

       maxWaveAct = maxval(abs(waveAct))

       ! x
       if( dimOut(1) ) then
          do i = 1,nRay
             write(unit=40,fmt="(4es17.6)") ray(i)%x * lRef
          end do
       else                                       ! height of ray symbols in 
          do i = 1,nRay
             write(unit=40,fmt="(4es17.6)") 0.05
          end do
       end if
       
       write(unit=40,fmt="(a)") ""

       ! y
       if( dimOut(2) ) then
          do i = 1,nRay
             write(unit=40,fmt="(4es17.6)") ray(i)%y * lRef
          end do
          write(unit=40,fmt="(a)") ""
       end if
       
       ! z
       do i = 1,nRay
          write(unit=40,fmt="(4es17.6)") ray(i)%z * lRef
       end do

!!$ ! check coordinates: ok
!!$print*,"check in tec360:"
!!$     do i = 1,nRay
!!$        write(*,fmt="(i3,2f10.2)") i, ray(i)%x*lRef, ray(i)%z*lRef
!!$     end do
!!$stop"test stop"

    end if ! raytracer



    !------------------------------------
    !              close file
    !------------------------------------
    close(unit=40)
    
    ! set counter
    iOut = iOut + 1


  end subroutine tec360


  !--------------------------------------------------------------------------


  subroutine readtec360( var, inputFile, time, scale )
    !-----------------------------------
    ! read tec360 data file for restart
    !-----------------------------------

    ! in/out variables
    real,dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),intent(out)  :: var
    character(len=*), intent(in)    :: inputFile

    real,intent(out)      :: time   ! non-dimensional time
    logical, intent(in)   :: scale  ! true: scale with reference quantities 
    !                                 in atmosphere_module
    !                                false: read in data in SI units -> without scaling


    ! local variables
    integer :: i,j,k, iVar
    integer :: nOutVar    ! nb. of total output variables 
    integer :: nDimOut    ! nb. of output dimensions
    real :: time_dim

    ! ghostCell output
    integer :: zoneSizeX, zoneSizeY, zoneSizeZ
    integer :: i0,i1, j0,j1, k0,k1

    ! vars specific to read routine
    real :: realTrash
    character(len=400) :: charTrash
    integer :: openStatus, readStatus, nTextLines

    ! notify user of new data-file format 
    print*,"output.f90/readtec360: Data file header was changed for ray tracing!"


    ! open input file
    open( unit=40, iostat=openStatus, file=inputFile, action="read", &
         form="formatted", status="old", position="rewind")
    if (openStatus /= 0) then
       print*, "iostat = ", openStatus
       print*, "output.f90/readtec360: problem opening file ", trim(inputFile)
       stop
    else
       print*,"successfully opened file ", trim(inputFile)
       print*,""
    end if

    ! get nb. of text lines of TECPLOT header  
    nTextLines = 0
    do i = 1, 100
       read( unit=40, fmt="(4es17.6)", iostat=readStatus) realTrash
       if (readStatus /=0) then ! error -> no real number
          nTextLines = nTextLines + 1
       else ! first real number was read 
          exit
       end if
    end do
    close( unit=40 )
    print*,"nb. of text lines = ", nTextLines


    ! start full read
    open( unit=40, iostat=openStatus, file=inputFile, action="read", &
         form="formatted", status="old", position="rewind")
    do i = 1, nTextLines
       read( unit=40, fmt="(a)") charTrash
       print*, trim(charTrash)
    end do
    print*,""

    ! read solution time and scale it
    if( solutionTime ) then
       read( unit=40, fmt="(4es17.6)") time_dim

       if ( solutionTimeUnit == "s" ) then
          print*,"solution time = ", time_dim," seconds"
          time = time_dim / tRef

       else if ( solutionTimeUnit == "min" ) then
          print*,"solution time = ", time_dim," minutes"
          time = time_dim*60.0 / tRef          
       else if ( solutionTimeUnit == "h" ) then
          print*,"solution time = ", time_dim," hours"
          time = time_dim*3600.0 / tRef          
       end if
    else
       time = 0.0
    end if
    print*,""; print*,""


    !---------------------------------------
    !   read xy-coordinates to trash
    !---------------------------------------

    if (dimOut(1)) then
       if( showGhostCellsX ) then
          i0 = -nbx-1;  i1 = nx+nbx
       else
          i0 = 0;  i1 = nx
       end if
    else
       i0 = 1; i1 = 1  
    end if
    if (dimOut(2)) then
       if( showGhostCellsY ) then
          j0 = -nby-1;  j1 = ny+nby
       else
          j0 = 0;  j1 = ny
       end if
    else
       j0 = 1; j1 = 1 
    end if
    if (dimOut(3)) then
       if( showGhostCellsZ ) then
          k0 = -nbz-1;  k1 = nz+nbz
       else
          k0 = 0;  k1 = nz
       end if
    else
       k0 = 1; k1 = 1
    end if

    ! x-coordinate
    if (dimOut(1)) then
       do k = k0, k1
          do j = j0, j1
             do i = i0, i1
                if (i == -nbx-1) then
                   read(unit=40,fmt="(4es17.6)") realTrash
                else
                   read(unit=40,fmt="(4es17.6)") realTrash
                end if
             end do
          end do
       end do
       read(unit=40,fmt="(a)" ) charTrash
    end if


    ! y-coordinate
    if (dimOut(2)) then
       do k = k0, k1
          do j = j0, j1
             do i = i0, i1
                if (j == -nby-1) then
                   read(unit=40,fmt="(4es17.6)") realTrash
                else
                   read(unit=40,fmt="(4es17.6)") realTrash
                end if
             end do
          end do
       end do
       read(unit=40,fmt="(a)" ) charTrash
    end if


    ! z-coordinate
    if (dimOut(3)) then
       do k = k0, k1
          do j = j0, j1
             do i = i0, i1
                if (k == -nbz-1) then
                   read(unit=40,fmt="(4es17.6)") realTrash
                else
                   read(unit=40,fmt="(4es17.6)") realTrash
                end if
             end do
          end do
       end do
       read(unit=40,fmt="(a)" ) charTrash
    end if



    !-------------------------------------------
    !         read primary variables
    !-------------------------------------------

    ! set field dimensions depending on ghost cell and output dimension
    if (dimOut(1)) then
       if( showGhostCellsX ) then
          i0 = -nbx;  i1 = nx+nbx
       else
          i0 = 1; i1 = nx
       end if
    else
       i0 = 1; i1 = 1 
    end if
    if (dimOut(2)) then
       if( showGhostCellsY ) then
          j0 = -nby;  j1 = ny+nby
       else
          j0 = 1; j1 = ny
       end if
    else
       j0 = 1; j1 = 1 
    end if
    if (dimOut(3)) then
       if( showGhostCellsZ ) then
          k0 = -nbz;  k1 = nz+nbz
       else
          k0 = 1; k1 = nz
       end if
    else
       k0 = 1; k1 = 1 
    end if


    !---------------------------------------
    !       read and make non-dimensional
    !---------------------------------------
    do iVar = 1, nVar
       if (varOut(iVar)==1) then
          do k = k0, k1
             do j = j0, j1
                do i = i0, i1

                   select case (iVar) 

                   case(1) ! density

                      if(showGhostCellsZ) then   ! rhoStrat not defined in ghost cells!!!
                         read(unit=40,fmt="(4es17.6)" ) var(i,j,k,iVar) 
                         if (scale) var(i,j,k,iVar) = var(i,j,k,iVar) / rhoRef  
                      else
                         read(unit=40,fmt="(4es17.6)" ) var(i,j,k,iVar)
                         if (scale) then
                            var(i,j,k,iVar) = var(i,j,k,iVar) / rhoRef
                            var(i,j,k,iVar) = var(i,j,k,iVar) + rhoStrat(k)
                         end if
                      end if

                   case(2:4) ! velocity
                      read(unit=40,fmt="(4es17.6)" ) var(i,j,k,iVar)
                      var(i,j,k,iVar) = var(i,j,k,iVar) + offset(iVar)
                      if (scale) var(i,j,k,iVar) = var(i,j,k,iVar) / uRef


                   case(5) ! Exner function (no scaling needed)
                      read(unit=40,fmt="(4es17.6)" ) var(i,j,k,iVar) 
                      
                   case(6) ! potential temperature deviation (Boussinesq)
                      read(unit=40,fmt="(4es17.6)" ) var(i,j,k,iVar)
                      if (scale) var(i,j,k,iVar) = var(i,j,k,iVar) / thetaRef
                      
                   end select

                end do
             end do
          end do
       end if
       read(unit=40,fmt="(a)" ) charTrash
    end do



    !------------------------------------
    !             close file
    !------------------------------------
    close(unit=40)


  end subroutine readtec360


  !--------------------------------------------------------------------------


  subroutine init_output
    !-----------------------
    ! allocate optVar field
    !-----------------------

    ! local variables
    integer :: allocstat

    ! divPu
    allocate(optVar(-nbx:nx+nbx, -nby:ny+nby, -nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "output.f90/init_output: could not allocate optVar. Stop."

  end subroutine init_output


  !--------------------------------------------------------------------------


  subroutine terminate_output
    !-------------------------
    ! deallocate optVar field
    !-------------------------

    ! local variables
    integer :: allocstat


    !---------------- deallocate variables -----------------------

    deallocate(optVar,stat=allocstat)
    if(allocstat /= 0) stop "output.f90/init_output: could not deallocate optVar. Stop."


  end subroutine terminate_output
  
  
  !--------------------------------------------------------------------------



end module output_module
