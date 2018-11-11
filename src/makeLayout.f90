program teclayout_program
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !% creates a tec360 layout for animating
  !% 1D-plots 
  !% 
  !% 
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  

  
  use type_module
  

  !----------------------------------------
  implicit none
  !----------------------------------------
  
  ! get grid and other info from input file
  ! open input file input.f90
  open (unit=10, file="input.f90", action="read", &
       form="formatted", status="old", position="rewind")
  
  ! read grid info
  read (unit=10, nml=grid)
  
  ! read output info
  read (unit=10, nml=outputList)
  read (unit=10, nml=testCaseList)

  
  close(unit=10) 
  

  call tecLayout(1,10,3)
  
  
  
  
  
contains  
  
  subroutine tecLayout(firstFile, lastFile, skipFile)

    !-----------------------------------
    !  create layout file for tecplot
    !------------------------------------

    ! in/out variables
    integer, intent(in) :: firstFile, lastFile,skipFile
    

    ! local variables
    character(len=*), parameter :: layoutFile = "anime_made.lay"
    character(len=*), parameter :: layoutTemplate = "layout.template"

    integer :: i,j, iFile
    integer :: nDimOut, nOutVar
    character(len=1000) :: line
    character(len=500)  :: variablesChar
    character(len=50)   :: dataFile, auxChar
    character(len=20)   :: form, gridName
    integer :: ios                        ! error status of read
    integer, parameter :: nbHeaderLines=4   ! number of header lines
    integer :: linePos                    ! position of Linemap
    integer :: iLineMap


    !---------------------
    !     open files
    !---------------------
    open( unit=50, file=layoutTemplate, action="read", &
         form="formatted", status="old")
    
    open( unit=60, file=layoutFile, action="write", &
         form="formatted", status="replace")

    !------------------------------
    !     advance input file 
    !------------------------------
    do i = 1,nbHeaderLines
       read( unit=50,fmt="(a)") line 
    end do

    !------------------------------
    !         write header
    !------------------------------

    !------------
    ! first line
    !------------

    ! name of data files
    dataFile = ""                                      ! create empty string
    dataFile = trim(dataFile)//trim(testCase)          ! test case 
    

    ! grid 
    write(unit=form, fmt="(a,i1,a,i1,a,i1,a)") &
         & "(i",int(log10(real(nx)))+1, &
         & ",a,i",int(log10(real(ny)))+1, &
         & ",a,i",int(log10(real(nz)))+1, &
         & ")" ! format for nx x ny x nz

    ! generate text string for grid size
    write(unit=gridName, fmt=form) nx,"x",ny,"x",nz
    dataFile = trim(dataFile)//"-"//trim(gridName)      ! grid size
    
    line = "$!VarSet |LFDSFN1| = """
    do i = firstFile, lastFile, skipFile
       write(unit=auxChar, fmt="(i6.6,a4)") i,".plt"
       line = trim(line)//" """//trim(dataFile)//"-"//trim(auxChar)//""" "
    end do
    line = trim(line)//" "" "
      
    write( unit=60,fmt="(a)") trim(line)


    !--------------
    ! second line
    !--------------

    if( dimOut(1) .and. dimOut(2) .and. dimOut(3) ) then    ! xyz plot
       nDimOut = 3
       
    else if( dimOut(1) .and. .not.dimOut(2) .and. dimOut(3) ) then ! xz plot
       nDimOut = 2

    else if( .not.dimOut(1) .and. .not.dimOut(2) .and. dimOut(3) ) then ! z-plot
       nDimOut = 2

    else
       stop "tec360: specify case for mDimOut"
    end if

    line = "$!VarSet |LFDSVL1| = '"

    ! obtain line for variables
    if( dimOut(1) .and. dimOut(2) .and. dimOut(3) ) then    ! xyz plot
       variablesChar = """x [m]"" ""y [m]"" ""z [m]"" "

    else if( dimOut(1) .and. .not.dimOut(2) .and. dimOut(3) ) then ! xz plot
       variablesChar = """x [m]"" ""z [m]"" "

    else if( .not.dimOut(1) .and. .not.dimOut(2) .and. dimOut(3) ) then ! z-plot
       variablesChar = """x [m]"" ""z [m]"" "

    else
       stop "tec360: specify case for line of variables"
    end if


 ! nb of variables initialised with nb of spatial vars
    nOutVar = nDimOut

    !-------------------
    ! primary variables
    !-------------------
    do i = 1,5
       if (varOut(i)==1) then
          nOutVar = nOutVar + 1
          if (i==1) then 
             if( rhoOffset ) then
                variablesChar = trim(variablesChar)//&
                     &" ""<greek>r</greek>' [kg/m<sup>3</sup>]""  "
             else
                variablesChar = trim(variablesChar)//&
                     &" ""<greek>r</greek> [kg/m<sup>3</sup>]""  "
             end if
          end if
          if (i==2) variablesChar = trim(variablesChar)//" ""u [m/s]""  "
          if (i==3) variablesChar = trim(variablesChar)//" ""v [m/s]""  "
          if (i==4) variablesChar = trim(variablesChar)//" ""w [m/s]""  "
          if (i==5) variablesChar = trim(variablesChar)//" ""<greek>p</greek>' ""  "
       end if
    end do

    !--------------------
    ! optional variables
    !--------------------
    do i = 1,size(optVarOut)
       if (optVarOut(i)==1) then
          nOutVar = nOutVar + 1
          if (i==1) variablesChar = trim(variablesChar)//&
               & " ""p' [kPa]""  "
          if (i==2) then
             if( thetaOffset ) then
                variablesChar = trim(variablesChar)//&
                     & " ""<greek>q</greek>' [K]""  "
             else
                variablesChar = trim(variablesChar)//&
                     & " ""<greek>q</greek> [K]""  "
             end if
          end if
          if (i==3) variablesChar = trim(variablesChar)//" ""p<sub>bg</sub> [kPa] ""  "
          if (i==4) variablesChar = trim(variablesChar)//&
               & " ""<greek>r</greek><sub>bg</sub> [kg/m<sup>3</sup>] ""  "
          if (i==5) variablesChar = trim(variablesChar)//" ""div(Pu) [Pa/s]""  "
          if (i==6) variablesChar = trim(variablesChar)//&
               & " ""b<sub>z</sub>/N<sup>2</sup> ""  "
       end if
    end do

    !------------------------
    ! WKB optional variables
    !------------------------
    do i = 1,size(wkbVarOut)
       if (wkbVarOut(i)==1) then
          nOutVar = nOutVar + 1
          if (i==1) variablesChar = trim(variablesChar)//&
               & " ""A [Js/m<sup>3</sup>]""  "

          if (i==2) variablesChar = trim(variablesChar)//&
               & " ""u<sub>0</sub><sup>(1)</sup>' [m/s]""  "
       end if
    end do
    !----------------------
    ! write variable line
    !----------------------
    line = trim(line)//trim(variablesChar)//"'"
    write( unit=60,fmt="(a)") trim(line)
!    print*,"variable line = ", trim(line)
!    stop "test stop"



    !-------------------------
    ! read and write template 
    !-------------------------
    
    do i = 1,5000
       read( unit=50,fmt="(a)", iostat=ios) line 
       if( ios/=0 ) stop "Stop: problem reading layout.template"
       
       ! identify first line for map description
       if( line(1:9) == "$!LINEMAP" ) then
          iLineMap = i + nbHeaderLines
          print*,"1. Linemap at line = ", iLineMap
          exit
       else
          write( unit=60,fmt="(a)") trim(line)
       end if
    end do

    !---------------------------
    !       create maps
    !---------------------------

    file_loop: do iFile = firstFile, lastFile, skipFile

       close(unit=50)
       
       ! open template
       open( unit=50, file=layoutTemplate, action="read", &
            form="formatted", status="old", position="rewind")

       ! advance to linemap line
       do j = 1, iLineMap
          read(unit=50,fmt="(a)", iostat=ios) line
          if( ios/=0 ) stop "stop pos 1"
       end do
       
       write( unit=60,fmt="(a9,a3,i2,a1)") line(1:9),"  [",iFile,"]"

       read(unit=50,fmt="(a)",iostat=ios) line
       if( ios/=0 ) stop "stop pos 2"
       write( unit=60,fmt="(a13,i2,a1)") line(1:13),iFile,"'"

       ! write map lines
       j_loop: do j = 1,1000
          read(unit=50,fmt="(a)", iostat=ios) line
          if( ios/=0 ) stop "stop pos 1"
          
          if( line(1:2) == "$!" ) then
             exit j_loop
          else          
             write( unit=60,fmt="(a)") trim(line)
          end if
       end do j_loop
       
    end do file_loop
    
    
    !------------------------
    !    write final part
    !------------------------

    ! advance to last part
    do i = 1,10000
       read(unit=50,fmt="(a)", iostat=ios) line
       if( ios/=0 ) stop "stop pos 1"
       if( line(1:12) == "$!XYLINEAXIS" ) exit
    end do
    
    do i = 1,10000
       write( unit=60,fmt="(a)") trim(line)
       read(unit=50,fmt="(a)", iostat=ios, end=100) line
       if( ios/=0 ) stop "stop pos 1"
    end do

100    close( unit=50 )
       close( unit=60 )

  end subroutine tecLayout


!--------------------------------------------------------------------------


end program teclayout_program
