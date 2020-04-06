module output_module

  use type_module
  use atmosphere_module

  implicit none

  private

  ! public subroutines
  public :: output_data
  public :: read_data
  public :: output_wkb
  public :: output_field
  public :: output_profile

  ! internal subroutines (listed for completeness)
  public :: init_output
  public :: terminate_output

  ! internal module variables (output variables)
  real, dimension(:,:,:), allocatable :: optVar


contains
  
  subroutine output_data( &
       & iOut, &
       & var,&
       & iTime, time, cpuTime )

    !-------------------------------
    !  writes data to file pf_all.dat
    !-------------------------------

    ! output counter
    integer, intent(inout) :: iOut
    
    ! argument fields
    real,dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),intent(in)  &
    & :: var
    
    ! local and global output field and record nbs.
    real*4,dimension(nx,ny) :: field_prc
    integer irc_prc,irc_out
    
    ! argument parameters
    integer, intent(in)            :: iTime
    real,intent(in)                :: time, cpuTime

    ! local variables
    integer :: i,j,k, iVar
    real :: time_dim

    ! needed for output to screen
    character (len = 20)  :: fmt, form
    character (len = 40)  :: cpuTimeChar

    ! CPU Time
    integer    :: days, hours, mins, secs
    real       :: timeVar

    ! hotBubble output
    real :: rho, theta_dim

    ! buoyancy
    real :: b, b_dim, theta

    integer :: i_prc,i_mst,i_out,j_prc,j_mst,j_out

    time_dim = time * tRef

    if( master ) then ! modified by Junhong Wei (20161110)
       print*,""
       print*," Output into File "
       print*,""
       write(*,fmt="(a25,i15)") " at time step = ", iTime
       write(*,fmt="(a25,f15.1,a8)") " at physical time = ", time_dim, &
                                      &" seconds"
    end if ! modified by Junhong Wei (20161110)
    
    !------------------------------
    !   prepare output file
    !------------------------------
    
    ! open output file

    if(master) then
       open(41,file='pf_all.dat',form="unformatted",access='direct',&
            & recl=SizeX*SizeY)
    end if

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
    
    if( master ) write(*,fmt="(a25,a25)")  "CPU time = ", cpuTimeChar
    
    !---------------------------------------
    !       dimensionalising and layerwise output
    !---------------------------------------

    irc_prc = 0
    do iVar = 1, nVar
       if (varOut(iVar)==1) irc_prc = irc_prc + 1
    end do
    irc_prc = irc_prc * iOut * nz

    do iVar = 1, nVar
       if (varOut(iVar)==1) then
          do k = 1, nz
             ! dimensionalization

             do j = 1, ny
                do i = 1, nx
                   select case (iVar) 

                      case(1) ! density
                        if(fluctuationMode) then
                           if(rhoOffset) then
                              field_prc(i,j) = var(i,j,k,iVar) * rhoRef 
                             else
                              field_prc(i,j) &
                              = (var(i,j,k,iVar) + rhoStrat(k)) * rhoRef 
                           end if
                          else
                           if(rhoOffset) then
                              field_prc(i,j) &
                              = (var(i,j,k,iVar) - rhoStrat(k)) * rhoRef 
                             else
                              field_prc(i,j) = var(i,j,k,iVar) * rhoRef 
                           end if
                        end if
                      
                      ! average velocities to cell center
                      
                      case(2) ! u velocity
                        field_prc(i,j) &
                        = (var(i,j,k,iVar) - offset(iVar)) &
                          * uRef

                      case(3) ! v velocity
                        field_prc(i,j) &
                        = (var(i,j,k,iVar) - offset(iVar)) &
                         * uRef 
                      
                      case(4) ! w velocity
                        field_prc(i,j) &
                        = (var(i,j,k,iVar) - offset(iVar)) &
                          * uRef 
                   
                      case(5) ! Exner function pi' 
                              !(deviation from background)
                        field_prc(i,j) = (var(i,j,k,iVar)) !+ piStrat(k))**(1./kappa)!FS

                      case(6) ! potential temperature theta' 
                              ! (deviation from background, Boussinesq)

                        select case(model) 

                           case("pseudo_incompressible")
                             if (fluctuationMode ) then
                                rho = var(i,j,k,1) + rhoStrat(k)
                               else
                                rho = var(i,j,k,1)
                             end if
                        
                             theta_dim = Pstrat(k) / rho * thetaRef 
                        
                             if( thetaOffset ) then
                                theta_dim = theta_dim - thetaStrat(k)*thetaRef
                             end if

                             field_prc(i,j) = theta_dim

                           case( "Boussinesq" )
                             field_prc(i,j) = var(i,j,k,iVar)*thetaRef 

                           case( "WKB" )
                        
                           case default
                              stop"tec360: unknown model"
                        end select ! model
                      
                      case(7) ! dynamic Smagorinsky coefficient
                              !(deviation from background)
                        field_prc(i,j) = var(i,j,k,iVar) * uRef * lRef 

                        case default
                           stop"tec360: unkown iVar"
                   end select ! iVar
                end do ! i
                call mpi_gather(field_prc(1,j),nx,mpi_real,&
                                field_mst(1,j),nx,mpi_real,0,comm,ierror)
             end do ! j

             ! layerwise output

             irc_prc=irc_prc+1
!            write(40,rec=irc_prc) field_prc
             call mpi_barrier(comm,ierror)
             if(master) then
                do j=1,ny
                   j_mst=j

                   do j_prc= 1,nprocy
                      j_out=ny*(j_prc-1)+j

                      do i_prc=1,nprocx
                         do i=1,nx
                            i_out=nx*(i_prc-1)+i

                            i_mst=nprocy*nx*(i_prc-1)+(j_prc-1)*nx+i

                            field_out(i_out,j_out)=field_mst(i_mst,j_mst)
                         end do
                      end do
                   end do
                end do

                write(41,rec=irc_prc) field_out
             end if
          end do ! k
       end if
    end do ! iVar

    !------------------------------------
    !              close file
    !------------------------------------

    if(master) close(unit=41)
    
    ! set counter
    iOut = iOut + 1

  end subroutine output_data

  !-------------------------------------------------------------------------

  subroutine read_data( &
       & iIn, &
       & var )

    !-------------------------------
    !  reads data from file pf_all_in.dat
    !-------------------------------
    
    ! input counter
    integer, intent(in) :: iIn
    
    ! argument fields
    real,dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar),intent(out)  &
    & :: var
    
    ! local and global output field and record nbs.
    real*4,dimension(nx,ny) :: field_prc
    integer irc_prc,irc_out
    
    ! local variables
    integer :: i,j,k, iVar

    ! needed for output to screen
    character (len = 20)  :: fmt, form

    ! hotBubble input
    real :: rho, theta_dim

    ! buoyancy
    real :: b, b_dim, theta

    integer :: i_prc,i_mst,i_out,j_prc,j_mst,j_out

    if( master ) then
       print*,""
       print*," Input from File "
       print*,""
       write(*,fmt="(a25,i15)") " reading record no. ", iIn
    end if
    
    ! open input file

    if(master) then
       open(40,file='pf_all_in.dat',form="unformatted",access='direct',&
            & recl=SizeX*SizeY)
!      print*,"pf_all_in.dat opened"
    end if

    var = 0.0

    !---------------------------------------
    !       layerwise input and non-dimensionalizing
    !---------------------------------------

    irc_prc = 0
    do iVar = 1, nVar
       if (varIn(iVar)==1) irc_prc = irc_prc + 1
    end do
    irc_prc = irc_prc * iIn * nz

    do iVar = 1, nVar
       if (varIn(iVar)==1) then
          do k = 1, nz
             ! read data layerwise

             irc_prc = irc_prc + 1

             if(master) then
                read(40,rec=irc_prc) field_out

                do j=1,ny
                   j_mst=j

                   do j_prc= 1,nprocy
                      j_out=ny*(j_prc-1)+j

                      do i_prc=1,nprocx
                         do i=1,nx
                            i_out=nx*(i_prc-1)+i
   
                            i_mst=nprocy*nx*(i_prc-1)+(j_prc-1)*nx+i

                            field_mst(i_mst,j_mst)=field_out(i_out,j_out)
                         end do
                      end do
                   end do
                end do
             end if

             call mpi_barrier(comm,ierror)

             do j = 1, ny
                ! data distributed over all processors

                call mpi_scatter(field_mst(1,j),nx,mpi_real,&
                                 field_prc(1,j),nx,mpi_real,0,comm,&
                                 ierror)
                do i = 1, nx
                   ! non-dimensionalization

                   select case (iVar) 

                      case(1) ! density
                        if(fluctuationMode) then
                           if(rhoOffset) then
                              var(i,j,k,iVar) = field_prc(i,j) / rhoRef
                             else
                              var(i,j,k,iVar) = field_prc(i,j) / rhoRef &
                                                - rhoStrat(k)
                           end if
                          else
                           if(rhoOffset) then
                              var(i,j,k,iVar) = field_prc(i,j) / rhoRef &
                                                + rhoStrat(k)
                             else
                              var(i,j,k,iVar) = field_prc(i,j) / rhoRef
                           end if
                        end if
                      
                      ! interpolate velocities to cell faces
                      
                      case(2) ! u velocity
                        var(i,j,k,iVar) = field_prc(i,j) / uRef &
                                          + offset(iVar)

                      case(3) ! v velocity
                        var(i,j,k,iVar) = field_prc(i,j) / uRef &
                                          + offset(iVar)

                      case(4) ! w velocity
                        var(i,j,k,iVar) = field_prc(i,j) / uRef &
                                          + offset(iVar)

                      case(5) ! Exner function pi' 
                              !(deviation from background)
                        var(i,j,k,iVar) = field_prc(i,j)

                      case(6) ! potential temperature theta' 
                              ! (deviation from background, Boussinesq)
                              ! potential temperature not used 
                              ! (only density  needed)

                      case(7) ! dynamic Smagorinsky coefficient
                              !(deviation from background)

                              var(i,j,k,iVar) = field_prc(i,j) / (uRef*lRef)

                        case default
                           stop"tec360: unkown iVar"
                   end select ! iVar
                end do ! i
             end do ! j
          end do ! k
       end if
    end do ! iVar

    !------------------------------------
    !              close file
    !------------------------------------

    if(master) close(unit=40)
    
  end subroutine read_data

  !-------------------------------------------------------------------------

  subroutine output_wkb( &
       & iOut, &
       & ray,&
       & ray_var3D )

    !--------------------------------------------------
    !  WKB output to pf_wkb_mean.dat and pf_wkb_ray.dat
    !--------------------------------------------------
    
    ! output counter
    integer, intent(inout) :: iOut
    
    ! wkb arguments
    type(rayType), dimension(nray_wrk,0:nx+1,0:ny+1,-1:nz+2), &
    & intent(in) :: ray

    real, dimension(0:nx+1,0:ny+1,0:nz+1,1:6), intent(in) :: ray_var3D

    ! local variables
    integer :: i,j,k, iVar

    ! local and global output field and record nbs.
    real*4,dimension(nx,ny) :: field_prc
    integer irc_prc,irc_out
    
    integer :: i_prc,i_mst,i_out,j_prc,j_mst,j_out

    ! set counter
    iOut = iOut - 1    
    
    !------------------------------
    !   prepare output file
    !------------------------------
    
    ! open output file

    if(master) then
       open(40,file='pf_wkb_mean.dat',form="unformatted",access='direct',&
            & recl=SizeX*SizeY)
    end if

    !---------------------------------------
    !       dimensionalising and layerwise output
    !---------------------------------------

    irc_prc = 6 * iOut * nz

    do iVar = 1, 6
       do k = 1, nz
          ! dimensionalization

          do j = 1, ny
             do i = 1, nx
                select case (iVar) 

                   case(1) ! (du/dt)_GW
                     field_prc(i,j) = ray_var3D(i,j,k,1) * uRef/tRef
                      
                   case(2) ! (dv/dt)_GW
                     field_prc(i,j) = ray_var3D(i,j,k,2) * uRef/tRef

                   case(3) ! (dtheta/dt)_GW
                     field_prc(i,j) = ray_var3D(i,j,k,3) * thetaRef/tRef

                   case(4) ! elastic term in x direction
                     field_prc(i,j) = ray_var3D(i,j,k,4) * uRef/tRef
                      
                   !case(5) ! elastic term in y direction
                   !  field_prc(i,j) = ray_var3D(i,j,k,5) * uRef/tRef
                   case(5) ! u'w' momentum flux , output chage by FDK
                     field_prc(i,j) = ray_var3D(i,j,k,5) * uRef**2

                   case(6) ! GW energy
                     field_prc(i,j) &
                     = ray_var3D(i,j,k,6) * rhoRef*uRef**2! /tRef deleted by FDK
                   
                   case default
                     stop"output_wkb: unkown iVar"
                end select ! iVar
             end do ! i
             call mpi_gather(field_prc(1,j),nx,mpi_real,&
                             field_mst(1,j),nx,mpi_real,0,comm,ierror)
          end do ! j

          ! layerwise output

          irc_prc=irc_prc+1
          call mpi_barrier(comm,ierror)
          if(master) then
             do j=1,ny
                j_mst=j

                do j_prc= 1,nprocy
                   j_out=ny*(j_prc-1)+j

                   do i_prc=1,nprocx
                      do i=1,nx
                         i_out=nx*(i_prc-1)+i

                         i_mst=nprocy*nx*(i_prc-1)+(j_prc-1)*nx+i

                         field_out(i_out,j_out)=field_mst(i_mst,j_mst)
                      end do
                   end do
                end do
             end do

             write(40,rec=irc_prc) field_out
          end if
       end do ! k
    end do ! iVar

    !------------------------------------
    !              close file
    !------------------------------------

    if(master) close(unit=40)

    !--------------------------------------------------------------
    ! here, after the parallelization of MS-GWaM, should also be 
    ! optional output for all ray-volume quantities (x,y,z,k,l,m,N)
    !--------------------------------------------------------------

    ! set counter
    iOut = iOut + 1

  end subroutine output_wkb


  !---------------------------------------------------------------------


  subroutine init_output
    !-----------------------
    ! allocate optVar field
    !-----------------------

    ! local variables
    integer :: allocstat

    ! divPu
    allocate(optVar(-nbx:nx+nbx, -nby:ny+nby, -nbz:nz+nbz),stat=allocstat)
    if(allocstat /= 0) stop "output.f90/init_output: could not allocate optVar. Stop."

    ! for output global fields
    allocate(field_mst(sizeX*nprocy,ny),stat=allocstat)
    if(allocstat /= 0) stop "output.f90/init_output: could not allocate field_mst. Stop."
    allocate(field_out(sizeX,sizeY),stat=allocstat)
    if(allocstat /= 0) stop "output.f90/init_output: could not allocate field_out. Stop."
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

  !------------------------------------------------------------------------
  ! output some singe field in time
    subroutine output_field( &
       & iOut, &
       & field, filename, scaling_coef )

    !-------------------------------
    !  writes data to file filename.dat
    !-------------------------------
    
    ! output counter
    integer, intent(in) :: iOut
    
    ! argument fields
    real,dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz),intent(in)  &
    & :: field
    
    ! local and global output field and record nbs.
    real*4,dimension(nx,ny) :: field_prc
    integer irc_prc,irc_out
    

    ! local variables
    integer :: i,j,k, iVar
    real :: time_dim

    ! needed for output to screen
    character (len = 20)  :: fmt, form
    character (len = 40)  :: cpuTimeChar

    ! CPU Time
    integer    :: days, hours, mins, secs
    real       :: timeVar

    ! hotBubble output
    real :: rho, theta_dim

    ! scale the variable
    real,intent(in)                :: scaling_coef
    character(len=*) :: filename

    ! buoyancy
    real :: b, b_dim, theta

    integer :: i_prc,i_mst,i_out,j_prc,j_mst,j_out


    if( master ) then ! modified by Junhong Wei (20161110)
       print*,""
       print*," Output into File "
       print*, filename

    end if ! modified by Junhong Wei (20161110)
    
    !------------------------------
    !   prepare output file
    !------------------------------
    
    ! open output file

    if(master) then
       open(21,file=filename,form="unformatted",access='direct',&
            & recl=SizeX*SizeY)
    end if


    
    !---------------------------------------
    !       dimensionalising and layerwise output
    !---------------------------------------

    irc_prc = 1
    irc_prc = irc_prc * iOut * nz

          do k = 1, nz
             ! dimensionalization

             do j = 1, ny
                do i = 1, nx
                    field_prc(i,j) = real(field(i,j,k) * scaling_coef &
                                & , kind=4)
                end do ! i
                call mpi_gather(field_prc(1,j),nx,mpi_real,&
                                field_mst(1,j),nx,mpi_real,0,comm,ierror)
             end do ! j

             ! layerwise output

             irc_prc=irc_prc+1
!            write(40,rec=irc_prc) field_prc
             call mpi_barrier(comm,ierror)
             if(master) then
                do j=1,ny
                   j_mst=j

                   do j_prc= 1,nprocy
                      j_out=ny*(j_prc-1)+j

                      do i_prc=1,nprocx
                         do i=1,nx
                            i_out=nx*(i_prc-1)+i

                            i_mst=nprocy*nx*(i_prc-1)+(j_prc-1)*nx+i

                            field_out(i_out,j_out)=field_mst(i_mst,j_mst)
                         end do
                      end do
                   end do
                end do

                write(21,rec=irc_prc) field_out
             end if
          end do ! k

    !------------------------------------
    !              close file
    !------------------------------------

    if(master) close(unit=21)
    

  end subroutine output_field

  !--------------------------------------------------------------------------
  ! output some singe field in time
    subroutine output_profile( &
       & iOut, &
       & field, filename )

    !-------------------------------
    !  writes data to file filename.dat
    !-------------------------------
    
    ! output counter
    integer, intent(inout) :: iOut
    
    ! argument fields
    real,dimension(-nbz:nz+nbz),intent(in)  &
    & :: field
    
    ! local and global output field and record nbs.
    real*4,dimension(nx,ny) :: field_prc
    integer irc_prc,irc_out
    

    ! local variables
    integer :: i,j,k, iVar
    real :: time_dim


    ! needed for output to screen
    character (len = 20)  :: fmt, form
    character (len = 40)  :: cpuTimeChar

    ! CPU Time
    integer    :: days, hours, mins, secs
    real       :: timeVar

    ! hotBubble output
    real :: rho, theta_dim

    character(len=*) :: filename

    ! buoyancy
    real :: b, b_dim, theta

    integer :: i_prc,i_mst,i_out,j_prc,j_mst,j_out

  !  iOut = iOut - 1 


    if( master ) then ! modified by Junhong Wei (20161110)
       print*,""
       print*," Output into File "
       print*, filename

    end if ! modified by Junhong Wei (20161110)
    
    !------------------------------
    !   prepare output file
    !------------------------------
    
    ! open output file

    if(master) then
       open(21,file=filename,form="unformatted",access='direct',&
            & recl=SizeX*SizeY)
    end if


    
    !---------------------------------------
    !       dimensionalising and layerwise output
    !---------------------------------------

    irc_prc = 1
    irc_prc = irc_prc * iOut * nz

          do k = 1, nz
             ! dimensionalization

             do j = 1, ny
                do i = 1, nx

                    field_prc(i,j) = real(field(k)  &
                                & , kind=4)
                end do ! i
                call mpi_gather(field_prc(1,j),nx,mpi_real,&
                                field_mst(1,j),nx,mpi_real,0,comm,ierror)
             end do ! j

             ! layerwise output

             irc_prc=irc_prc+1
!            write(40,rec=irc_prc) field_prc
             call mpi_barrier(comm,ierror)
             if(master) then
                do j=1,ny
                   j_mst=j

                   do j_prc= 1,nprocy
                      j_out=ny*(j_prc-1)+j

                      do i_prc=1,nprocx
                         do i=1,nx
                            i_out=nx*(i_prc-1)+i

                            i_mst=nprocy*nx*(i_prc-1)+(j_prc-1)*nx+i

                            field_out(i_out,j_out)=field_mst(i_mst,j_mst)

                         end do
                      end do
                   end do
                end do

                write(21,rec=irc_prc) field_out
             end if
          end do ! k

    !------------------------------------
    !              close file
    !------------------------------------

    if(master) close(unit=21)


   ! iOut = iOut + 1 

  end subroutine output_profile

end module output_module
