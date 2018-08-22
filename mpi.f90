module mpi_module

! Original: 
! F. Rieper 2011: MPI for square domains in x-y
!
! Modifications:
! G. Bölöni 27.10.2016
!  - MPI for any rectangle domain in x-y 
!    (see also boundary.f90)
!  - add Theta to setHalos
!  - cleaning of setHalos (case "pressure")


  use type_module
  use flux_module

  implicit none

  private         ! private module

  public :: init_mpi
  public :: setHalos
  public :: dot_product3D_glob
  public :: abort


contains

  subroutine abort( message )
    character(len=*), intent(in) :: message
    
    if( master ) then
       print*, message
    end if
    
    call mpi_finalize(ierror)
    stop
    
  end subroutine abort
  
  
!----------------------------------------------------------------------------------


  subroutine setHalos( var, option )
    !-------------------------------
    !  set values in halo cells 
    !-------------------------------

    ! in/out variables      
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    character(len=*), intent(in) :: option

    ! auxiliary fields for "var" with ghost cells (rho)
    real, dimension(nbx,-nby:ny+nby,nz) :: xRhoSliceLeft_send, xRhoSliceRight_send
    real, dimension(nbx,-nby:ny+nby,nz) :: xRhoSliceLeft_recv, xRhoSliceRight_recv

    real, dimension(-nbx:nx+nbx,nby,nz) :: yRhoSliceBack_send, yRhoSliceForw_send
    real, dimension(-nbx:nx+nbx,nby,nz) :: yRhoSliceBack_recv, yRhoSliceForw_recv

    ! auxiliary fields for "pressure", u,v,w 
    real, dimension(nbx,ny,nz) :: xSliceLeft_send, xSliceRight_send
    real, dimension(nbx,ny,nz) :: xSliceLeft_recv, xSliceRight_recv

    real, dimension(nx,nby,nz) :: ySliceBack_send, ySliceForw_send
    real, dimension(nx,nby,nz) :: ySliceBack_recv, ySliceForw_recv

    ! field for staggered quants
    real, dimension(0:nx,nby,nz) :: ySliceBackU_send, ySliceForwU_send
    real, dimension(0:nx,nby,nz) :: ySliceBackU_recv, ySliceForwU_recv

    real, dimension(nbx,0:ny,nz) :: xSliceLeftV_send, xSliceRightV_send
    real, dimension(nbx,0:ny,nz) :: xSliceLeftV_recv, xSliceRightV_recv

    ! aux fields for "varTilde"
    real, dimension(ny,nz) :: x1SliceLeft_send, x1SliceRight_send
    real, dimension(ny,nz) :: x1SliceLeft_recv, x1SliceRight_recv

    real, dimension(nx,nz) :: y1SliceBack_send, y1SliceForw_send
    real, dimension(nx,nz) :: y1SliceBack_recv, y1SliceForw_recv



    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: i, j, k, iVar
    integer :: i0, j0, k0


    !-----------------------------
    !     Find neighbour procs
    !-----------------------------

    if (idim > 1) call mpi_cart_shift(comm,0,1,left,right,ierror)
    if (jdim > 1) call mpi_cart_shift(comm,1,1,back,forw,ierror)

    select case( option )

    case( "var" ) 

       !------------------------------
       !          x-direction
       !------------------------------

       !------------------------------
       !   rho transport: iVar = 1
       !------------------------------
       !if( updateMass ) then
          iVar = 1

          ! slice size
          sendcount = nbx*(ny+2*nby+1)*nz
          recvcount = sendcount


          ! read slice into contiguous array
          do i = 1,nbx
             xRhoSliceLeft_send (i,-nby:ny+nby,1:nz) = var(i,-nby:ny+nby,1:nz,iVar)
             xRhoSliceRight_send(i,-nby:ny+nby,1:nz) = var(nx-nbx+i,-nby:ny+nby,1:nz,iVar)
          end do

          if ( idim > 1 ) then

             ! left -> right
             source = left
             dest = right
             tag = 100

             i0 = 1; j0 = -nby; k0 = 1

             call mpi_sendrecv(xRhoSliceRight_send(i0,j0,k0), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & xRhoSliceLeft_recv(i0,j0,k0), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_left, ierror)

             ! right -> left
             source = right
             dest = left
             tag = 100

             call mpi_sendrecv(xRhoSliceLeft_send(i0,j0,k0), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & xRhoSliceRight_recv(i0,j0,k0), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_right, ierror)

             ! write auxiliary slice to var field
             do i = 1,nbx

                ! right halos
                var(nx+i,-nby:ny+nby,1:nz,iVar) = xRhoSliceRight_recv(i,-nby:ny+nby,1:nz)

                ! left halos
                var(-nbx+i,-nby:ny+nby,1:nz,iVar) = xRhoSliceLeft_recv(i,-nby:ny+nby,1:nz)

             end do

          end if

          if(verbose .and. master) print*,"horizontalHalos: &
               & x-horizontal halos copied."

       !end if ! updateMass


       !------------------------------
       !   u transport: iVar = 2
       !------------------------------
       if( predictMomentum .or. correctMomentum ) then

          iVar = 2

          ! slice size
          sendcount = nbx*ny*nz
          recvcount = sendcount

          ! read slice into contiguous array
          !------------------------------------------------------------------
          ! Note: value at the local domain interface is not interchanged
          ! If this causes a problem: interchange values and take their mean
          !------------------------------------------------------------------
          do i = 1,nbx
             xSliceLeft_send (i,1:ny,1:nz) = var(i,1:ny,1:nz,iVar)
             xSliceRight_send(i,1:ny,1:nz) = var(nx-nbx-1+i,1:ny,1:nz,iVar)
          end do

          if( idim > 1 ) then

             ! left -> right
             source = left
             dest = right
             tag = 100

             call mpi_sendrecv(xSliceRight_send(1,1,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & xSliceLeft_recv(1,1,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_left, ierror)

             ! right -> left
             source = right
             dest = left
             tag = 100

             call mpi_sendrecv(xSliceLeft_send(1,1,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & xSliceRight_recv(1,1,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_right, ierror)

             ! write auxiliary slice to var field
             do i = 1,nbx

                ! right halos
                var(nx+i,1:ny,1:nz,iVar) = xSliceRight_recv(i,1:ny,1:nz)

                ! left halos
                var(-nbx-1+i,1:ny,1:nz,iVar) = xSliceLeft_recv(i,1:ny,1:nz)

             end do

          end if

          if(verbose .and. master) print*,"horizontalHalos: &
               & x-horizontal halos copied."

       end if  ! predictMomentum


       !------------------------------
       !    v transport: iVar = 3
       !------------------------------
       if( predictMomentum  .or. correctMomentum ) then
          iVar = 3

          ! slice size
          sendcount = nbx*(1+ny)*nz
          recvcount = sendcount


          ! read slice into contiguous array
          do i = 1,nbx
             xSliceLeftV_send (i,0:ny,1:nz) = var(i,0:ny,1:nz,iVar)
             xSliceRightV_send(i,0:ny,1:nz) = var(nx-nbx+i,0:ny,1:nz,iVar)
          end do

          if( idim > 1 ) then

             ! left -> right
             source = left
             dest = right
             tag = 100

             call mpi_sendrecv(xSliceRightV_send(1,0,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & xSliceLeftV_recv(1,0,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_left, ierror)

             ! right -> left
             source = right
             dest = left
             tag = 100

             call mpi_sendrecv(xSliceLeftV_send(1,0,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & xSliceRightV_recv(1,0,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_right, ierror)

             ! write auxiliary slice to var field
             do i = 1,nbx

                ! right halos
                var(nx+i,0:ny,1:nz,iVar) = xSliceRightV_recv(i,0:ny,1:nz)

                ! left halos
                var(-nbx+i,0:ny,1:nz,iVar) = xSliceLeftV_recv(i,0:ny,1:nz)

             end do

          end if

          if(verbose .and. master) print*,"horizontalHalos: &
               & x-horizontal halos copied."

       end if ! predictMomentum


       !------------------------------
       !    w transport: iVar = 4
       !------------------------------
       if( predictMomentum .or. correctMomentum ) then
          iVar = 4

          ! slice size
          sendcount = nbx*ny*nz
          recvcount = sendcount


          ! read slice into contiguous array
          do i = 1,nbx
             xSliceLeft_send (i,1:ny,1:nz) = var(i,1:ny,1:nz,iVar)
             xSliceRight_send(i,1:ny,1:nz) = var(nx-nbx+i,1:ny,1:nz,iVar)
          end do

          if( idim > 1 ) then

             ! left -> right
             source = left
             dest = right
             tag = 100

             call mpi_sendrecv(xSliceRight_send(1,1,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & xSliceLeft_recv(1,1,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_left, ierror)

             ! right -> left
             source = right
             dest = left
             tag = 100

             call mpi_sendrecv(xSliceLeft_send(1,1,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & xSliceRight_recv(1,1,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_right, ierror)

             ! write auxiliary slice to var field
             do i = 1,nbx

                ! right halos
                var(nx+i,1:ny,1:nz,iVar) = xSliceRight_recv(i,1:ny,1:nz)

                ! left halos
                var(-nbx+i,1:ny,1:nz,iVar) = xSliceLeft_recv(i,1:ny,1:nz)

             end do

          end if

          if(verbose .and. master) print*,"horizontalHalos: &
               & x-horizontal halos copied."

       end if ! predictMomentum

       !------------------------------
       !   pressure: iVar = 5
       !------------------------------
       iVar = 5

       ! slice size
       sendcount = nbx*ny*nz
       recvcount = sendcount


       ! read slice into contiguous array
       do i = 1,nbx
          xSliceLeft_send (i,1:ny,1:nz) = var(i,1:ny,1:nz,iVar)
          xSliceRight_send(i,1:ny,1:nz) = var(nx-nbx+i,1:ny,1:nz,iVar)
       end do

       if( idim > 1 ) then

          ! left -> right
          source = left
          dest = right
          tag = 100

          call mpi_sendrecv(xSliceRight_send(1,1,1), sendcount, &
               & mpi_double_precision, dest, tag, &
               & xSliceLeft_recv(1,1,1), recvcount, mpi_double_precision, &
               & source, mpi_any_tag, comm, sts_left, ierror)

          ! right -> left
          source = right
          dest = left
          tag = 100

          call mpi_sendrecv(xSliceLeft_send(1,1,1), sendcount, &
               & mpi_double_precision, dest, tag, &
               & xSliceRight_recv(1,1,1), recvcount, mpi_double_precision, &
               & source, mpi_any_tag, comm, sts_right, ierror)

          ! write auxiliary slice to var field
          do i = 1,nbx

             ! right halos
             var(nx+i,1:ny,1:nz,iVar) = xSliceRight_recv(i,1:ny,1:nz)

             ! left halos
             var(-nbx+i,1:ny,1:nz,iVar) = xSliceLeft_recv(i,1:ny,1:nz)

          end do

       end if

       if(verbose .and. master) print*,"horizontalHalos: &
            & x-horizontal halos copied."

       !------------------------------
       !   Theta transport: iVar = 6
       !------------------------------
       !if( updateTheta ) then
          iVar = 6

          ! slice size
          sendcount = nbx*(ny+2*nby+1)*nz
          recvcount = sendcount


          ! read slice into contiguous array
          do i = 1,nbx
             xRhoSliceLeft_send (i,-nby:ny+nby,1:nz) = var(i,-nby:ny+nby,1:nz,iVar)
             xRhoSliceRight_send(i,-nby:ny+nby,1:nz) = var(nx-nbx+i,-nby:ny+nby,1:nz,iVar)
          end do

          if ( idim > 1 ) then

             ! left -> right
             source = left
             dest = right
             tag = 100

             i0 = 1; j0 = -nby; k0 = 1

             call mpi_sendrecv(xRhoSliceRight_send(i0,j0,k0), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & xRhoSliceLeft_recv(i0,j0,k0), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_left, ierror)

             ! right -> left
             source = right
             dest = left
             tag = 100

             call mpi_sendrecv(xRhoSliceLeft_send(i0,j0,k0), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & xRhoSliceRight_recv(i0,j0,k0), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_right, ierror)

             ! write auxiliary slice to var field
             do i = 1,nbx

                ! right halos
                var(nx+i,-nby:ny+nby,1:nz,iVar) = xRhoSliceRight_recv(i,-nby:ny+nby,1:nz)

                ! left halos
                var(-nbx+i,-nby:ny+nby,1:nz,iVar) = xRhoSliceLeft_recv(i,-nby:ny+nby,1:nz)

             end do

          end if

          if(verbose .and. master) print*,"horizontalHalos: &
               & x-horizontal halos copied."

       !end if ! updateTheta


       !------------------------------
       !          y-direction
       !------------------------------

       !------------------------------
       !   rho transport: iVar = 1
       !------------------------------
       !if( updateMass ) then
          iVar = 1

          ! slice size
          sendcount = nby*(nx+2*nbx+1)*nz
          recvcount = sendcount


          ! read slice into contiguous array
          do j = 1, nby
             yRhoSliceBack_send(-nbx:nx+nbx,j,1:nz) = var(-nbx:nx+nbx,j,1:nz,iVar)
             yRhoSliceForw_send(-nbx:nx+nbx,j,1:nz) = var(-nbx:nx+nbx,ny-nby+j,1:nz,iVar)
          end do

          if( jdim > 1 ) then

             ! back -> forw
             source = back
             dest = forw
             tag = 100
             
             i0 = -nbx; j0 = 1; k0 = 1

             call mpi_sendrecv(yRhoSliceForw_send(i0,j0,k0), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & yRhoSliceBack_recv(i0,j0,k0), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_back, ierror)

             ! forw -> back
             source = forw
             dest = back
             tag = 100

             call mpi_sendrecv(yRhoSliceBack_send(i0,j0,k0), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & yRhoSliceForw_recv(i0,j0,k0), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_forw, ierror)

             ! write auxiliary slice to var field
             do j = 1, nby

                ! right halos
                var(-nbx:nx+nbx,ny+j,1:nz,iVar) = yRhoSliceForw_recv(-nbx:nx+nbx,j,1:nz)

                ! left halos
                var(-nbx:nx+nbx,-nby+j,1:nz,iVar) = yRhoSliceBack_recv(-nbx:nx+nbx,j,1:nz)

             end do

          end if

          if(verbose .and. master) print*,"horizontalHalos: &
               & x-horizontal halos copied."


       !end if ! updateMass


       !------------------------------
       !   u transport: iVar = 2
       !------------------------------
       if( predictMomentum  .or. correctMomentum ) then
          iVar = 2

          ! slice size
          sendcount = nby*(1+nx)*nz
          recvcount = sendcount


          ! read slice into contiguous array
          do j = 1, nby
             ySliceBackU_send(0:nx,j,1:nz) = var(0:nx,j,1:nz,iVar)
             ySliceForwU_send(0:nx,j,1:nz) = var(0:nx,ny-nby+j,1:nz,iVar)
          end do

          if( jdim > 1 ) then

             ! back -> forw
             source = back
             dest = forw
             tag = 100

             call mpi_sendrecv(ySliceForwU_send(0,1,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & ySliceBackU_recv(0,1,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_back, ierror)

             ! forw -> back
             source = forw
             dest = back
             tag = 100

             call mpi_sendrecv(ySliceBackU_send(0,1,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & ySliceForwU_recv(0,1,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_forw, ierror)

             ! write auxiliary slice to var field
             do j = 1, nby

                ! right halos
                var(0:nx,ny+j,1:nz,iVar) = ySliceForwU_recv(0:nx,j,1:nz)

                ! left halos
                var(0:nx,-nby+j,1:nz,iVar) = ySliceBackU_recv(0:nx,j,1:nz)

             end do

          end if

          if(verbose .and. master) print*,"horizontalHalos: &
               & x-horizontal halos copied."

       end if ! predictMomentum


       !------------------------------
       !   v transport: iVar = 3
       !------------------------------
       
       if( predictMomentum .or. correctMomentum ) then
          iVar = 3

          ! slice size
          sendcount = nby*nx*nz
          recvcount = sendcount

          ! read slice into contiguous array
          !xxx: I guess nby-1 should be enough: 2 ghost cells are needed
          do j = 1, nby
             ySliceBack_send(1:nx,j,1:nz) = var(1:nx,j,1:nz,iVar)
             ySliceForw_send(1:nx,j,1:nz) = var(1:nx,ny-nby-1+j,1:nz,iVar)
          end do

          if( jdim > 1 ) then

             ! back -> forw
             source = back
             dest = forw
             tag = 100

             call mpi_sendrecv(ySliceForw_send(1,1,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & ySliceBack_recv(1,1,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_back, ierror)

             ! forw -> back
             source = forw
             dest = back
             tag = 100

             call mpi_sendrecv(ySliceBack_send(1,1,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & ySliceForw_recv(1,1,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_forw, ierror)

             ! write auxiliary slice to var field -> into halo region
             do j = 1, nby

                ! forward halos
                var(1:nx,ny+j,1:nz,iVar) = ySliceForw_recv(1:nx,j,1:nz)

                ! backward halos
                var(1:nx,-nby-1+j,1:nz,iVar) = ySliceBack_recv(1:nx,j,1:nz)

             end do

          end if

          if(verbose .and. master) print*,"horizontalHalos: &
               & x-horizontal halos copied."

       end if ! predictMomentum


       !------------------------------
       !   w transport: iVar = 4
       !------------------------------
       if( predictMomentum .or. correctMomentum ) then
          iVar = 4

          ! slice size
          sendcount = nby*nx*nz
          recvcount = sendcount


          ! read slice into contiguous array
          do j = 1, nby
             ySliceBack_send(1:nx,j,1:nz) = var(1:nx,j,1:nz,iVar)
             ySliceForw_send(1:nx,j,1:nz) = var(1:nx,ny-nby+j,1:nz,iVar)
          end do

          if( jdim > 1 ) then

             ! back -> forw
             source = back
             dest = forw
             tag = 100

             call mpi_sendrecv(ySliceForw_send(1,1,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & ySliceBack_recv(1,1,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_back, ierror)

             ! forw -> back
             source = forw
             dest = back
             tag = 100

             call mpi_sendrecv(ySliceBack_send(1,1,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & ySliceForw_recv(1,1,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_forw, ierror)

             ! write auxiliary slice to var field
             do j = 1, nby

                ! right halos
                var(1:nx,ny+j,1:nz,iVar) = ySliceForw_recv(1:nx,j,1:nz)

                ! left halos
                var(1:nx,-nby+j,1:nz,iVar) = ySliceBack_recv(1:nx,j,1:nz)

             end do

          end if

          if(verbose .and. master) print*,"horizontalHalos: &
               & x-horizontal halos copied."


       end if ! predictMomentum

       !------------------------------
       !   pressure: iVar = 5
       !------------------------------
       iVar = 5

       ! slice size
       sendcount = nby*nx*nz
       recvcount = sendcount


       ! read slice into contiguous array
       do j = 1, nby
          ySliceBack_send(1:nx,j,1:nz) = var(1:nx,j,1:nz,iVar)
          ySliceForw_send(1:nx,j,1:nz) = var(1:nx,ny-nby+j,1:nz,iVar)
       end do

       if( jdim > 1 ) then

          ! back -> forw
          source = back
          dest = forw
          tag = 100

          call mpi_sendrecv(ySliceForw_send(1,1,1), sendcount, &
               & mpi_double_precision, dest, tag, &
               & ySliceBack_recv(1,1,1), recvcount, mpi_double_precision, &
               & source, mpi_any_tag, comm, sts_back, ierror)

          ! forw -> back
          source = forw
          dest = back
          tag = 100

          call mpi_sendrecv(ySliceBack_send(1,1,1), sendcount, &
               & mpi_double_precision, dest, tag, &
               & ySliceForw_recv(1,1,1), recvcount, mpi_double_precision, &
               & source, mpi_any_tag, comm, sts_right, ierror)

          ! write auxiliary slice to var field
          do j = 1, nby

             ! right halos
             var(1:nx,ny+j,1:nz,iVar) = ySliceForw_recv(1:nx,j,1:nz)

             ! left halos
             var(1:nx,-nby+j,1:nz,iVar) = ySliceBack_recv(1:nx,j,1:nz)

          end do

       end if

       if(verbose .and. master) print*,"horizontalHalos: &
            & x-horizontal halos copied."

       !------------------------------
       !   Theta transport: iVar = 6
       !------------------------------
       !if( updateTheta ) then
          iVar = 6

          ! slice size
          sendcount = nby*(nx+2*nbx+1)*nz
          recvcount = sendcount


          ! read slice into contiguous array
          do j = 1, nby
             yRhoSliceBack_send(-nbx:nx+nbx,j,1:nz) = var(-nbx:nx+nbx,j,1:nz,iVar)
             yRhoSliceForw_send(-nbx:nx+nbx,j,1:nz) = var(-nbx:nx+nbx,ny-nby+j,1:nz,iVar)
          end do

          if( jdim > 1 ) then

             ! back -> forw
             source = back
             dest = forw
             tag = 100

             i0 = -nbx; j0 = 1; k0 = 1

             call mpi_sendrecv(yRhoSliceForw_send(i0,j0,k0), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & yRhoSliceBack_recv(i0,j0,k0), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_back, ierror)

             ! forw -> back
             source = forw
             dest = back
             tag = 100

             call mpi_sendrecv(yRhoSliceBack_send(i0,j0,k0), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & yRhoSliceForw_recv(i0,j0,k0), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_forw, ierror)

             ! write auxiliary slice to var field
             do j = 1, nby

                ! right halos
                var(-nbx:nx+nbx,ny+j,1:nz,iVar) = yRhoSliceForw_recv(-nbx:nx+nbx,j,1:nz)

                ! left halos
                var(-nbx:nx+nbx,-nby+j,1:nz,iVar) = yRhoSliceBack_recv(-nbx:nx+nbx,j,1:nz)

             end do

          end if

          if(verbose .and. master) print*,"horizontalHalos: &
               & x-horizontal halos copied."


       !end if ! updateTheta


    case( "varTilde" )


       !-----------------------------------
       !             varTilde
       !-----------------------------------


       !------------------------------
       !          x-direction
       !------------------------------
       ! reconstructed density needed in ghost cell i = nx+2
       ! ...in ghost cell i = -1


       !if( updateMass ) then

          ! slice size
          sendcount = ny*nz
          recvcount = sendcount

          ! read slice into contiguous array
          x1SliceLeft_send (1:ny,1:nz) = rhoTilde(2,1:ny,1:nz,1,0)
          x1SliceRight_send(1:ny,1:nz) = rhoTilde(nx-1,1:ny,1:nz,1,1)


          if( idim > 1 ) then

             ! left -> right
             source = left
             dest = right
             tag = 100

             call mpi_sendrecv(x1SliceRight_send(1,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & x1SliceLeft_recv(1,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_left, ierror)

             ! right -> left
             source = right
             dest = left
             tag = 100
             call mpi_sendrecv(x1SliceLeft_send(1,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & x1SliceRight_recv(1,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_right, ierror)


             ! write auxiliary slice to var field

             ! right halos
             rhoTilde(nx+2,1:ny,1:nz,1,0) = x1SliceRight_recv(1:ny,1:nz)

             ! left halos
             rhoTilde(-1,1:ny,1:nz,1,1) = x1SliceLeft_recv(1:ny,1:nz)

          end if

          if(verbose .and. master) print*,"horizontalHalos: &
               & x-horizontal halos copied."

       !end if

       !------------------------------
       !          y-direction
       !------------------------------
       ! reconstructed density needed in ghost cell j = ny+2
       ! ...in ghost cell j = -1


       !if( updateMass ) then

          ! slice size
          sendcount = nx*nz
          recvcount = sendcount

          ! read slice into contiguous array
          y1SliceBack_send(1:nx,1:nz) = rhoTilde(1:nx,2,1:nz,2,0)
          y1SliceForw_send(1:nx,1:nz) = rhoTilde(1:nx,ny-1,1:nz,2,1)

          if( jdim > 1 ) then

             ! back -> forw
             source = back
             dest = forw
             tag = 100

             call mpi_sendrecv(y1SliceForw_send(1,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & y1SliceBack_recv(1,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_back, ierror)

             ! forw -> back
             source = forw
             dest = back
             tag = 100
             call mpi_sendrecv(y1SliceBack_send(1,1), sendcount, &
                  & mpi_double_precision, dest, tag, &
                  & y1SliceForw_recv(1,1), recvcount, mpi_double_precision, &
                  & source, mpi_any_tag, comm, sts_right, ierror)


             ! write auxiliary slice to var field

             ! right halos
             rhoTilde(1:nx,ny+2,1:nz,2,0) = y1SliceForw_recv(1:nx,1:nz)

             ! left halos
             rhoTilde(1:nx,-1,1:nz,2,1) = y1SliceBack_recv(1:nx,1:nz)

          end if

          if(verbose .and. master) print*,"horizontalHalos: &
               & x-horizontal halos copied."


       !end if ! updateMass

    case default

       if( master ) then
          print*,"setHalo: Unknown case. Stop."
          call mpi_barrier(comm,ierror)
          call mpi_finalize(ierror)
          stop
       end if

    end select  ! option "var", "varTilde", ...


  end subroutine setHalos


  !------------------------------------------------------------------


  subroutine init_mpi(error_flag)

    ! in/out vars
    logical, intent(out) :: error_flag 

    ! local vars
    integer :: root

    include 'mpif.h' 

    error_flag = .false. 
    verboseMPI = .false.


    !-----------------------
    !    basic MPI info
    !-----------------------

    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, rank, ierror)
    if( rank .eq. 0 ) then
       master = .true. 
    else
       master = .false.
    end if
    call mpi_comm_size(mpi_comm_world, nbProc, ierror)


    !-----------------------------
    !    Domain decomposition
    !-----------------------------

    ! open input file input.f90 by master
    if( master ) then
       open (unit=10, file="input.f90", action="read", &
            form="formatted", status="old", position="rewind")

       ! read grid info
       read (unit=10, nml=domain)
       close(10)

    end if

    ! broadcast input file data
    count = 1
    root = 0
    call mpi_bcast(sizeX,count,mpi_integer,root, mpi_comm_world, ierror)
    call mpi_bcast(sizeY,count,mpi_integer,root, mpi_comm_world, ierror)
    call mpi_bcast(sizeZ,count,mpi_integer,root, mpi_comm_world, ierror)
    call mpi_bcast(nbx,count,mpi_integer,root, mpi_comm_world, ierror)
    call mpi_bcast(nby,count,mpi_integer,root, mpi_comm_world, ierror)
    call mpi_bcast(nbz,count,mpi_integer,root, mpi_comm_world, ierror)
    call mpi_bcast(lx_dim,2*count,mpi_double_precision,root, mpi_comm_world, ierror)
    call mpi_bcast(ly_dim,2*count,mpi_double_precision,root, mpi_comm_world, ierror)
    call mpi_bcast(lz_dim,2*count,mpi_double_precision,root, mpi_comm_world, ierror)
    call mpi_bcast(nprocx,count,mpi_integer,root, mpi_comm_world, ierror)
    call mpi_bcast(nprocy,count,mpi_integer,root, mpi_comm_world, ierror)

    ! domain size with ghost cells
    sizeXX = sizeX + 2*nbx + 1
    sizeYY = sizeY + 2*nby + 1 

    iStart = -nbx + 1        ! left start index of global domain
    jStart = -nby + 1        ! back start index of global domain 


    ! set up no. of cpus in x-y direction
    idim = nprocx
    jdim = nprocy
    if (idim*jdim/=nbProc) then
       if (master) print *, 'Virtual topology wrong idim*jdim should be equal to nbProc!'
       error_flag = .true.
       return
    end if
    if( master ) then
!       print"(a,i3,a,i3,a)","Virtual topology: (idim, jdim) = (",idim, " ,",jdim,")"
       print"(a,i4,a,i4,a)","Virtual topology: (idim, jdim) = (",idim, " ,",jdim,")"
    end if

    dims(1) = idim
    dims(2) = jdim
    period(1)=.true.
    period(2)=.true.

    call mpi_cart_create(mpi_comm_world,2,dims,&
         &               period,.true.,comm,ierror)
    call mpi_comm_rank(comm, rank, ierror)
    call mpi_cart_coords(comm, rank, 2, coords, ierror) 
    icoord = coords(1) 
    jcoord = coords(2) 

    !--------------------------------------
    !     even-sized domain composition
    !--------------------------------------  

    ! whole y indecees   |------------- isize = 1+imax-istart ------------|
    !    start/end index  ^-istart                                  imax-^
    ! 
    ! 1. interval        |--- iouter1---|   
    !                    |--|--------|--|
    !                     b1  iinner1 b1 
    !    start/end index  ^-is      ie-^ 
    ! 2. interval                 |--|--------|--|   
    ! 3. interval                          |--|--------|--| 
    ! 4. interval                                   |--|-------|--| 
    !                                                   iinner0 
    ! 5. interval = idim's interval                         |--|-------|--|

    ! In each iteration on each interval, the inner area is computed
    ! by using the values of the last iteration in the whole outer area. 

    ! icoord = number of the interval - 1
    ! 
    ! To fit exactly into isize, we use in1 intervals of with iinner1
    ! and (idim-in1) intervals of with (iinner1 - 1) 

    ! isize <= 2*b1 + idim * iinner1  

    ! ==>   iinner1 >= (isize - 2*b1) / idim
    ! ==>   smallest valid integer value:

    ! Replacements: HLRS code -> pincFloit
    ! iouter = nxx
    ! iinner = nx
    ! b1 = nbx

    ! composition along x coordinate
    if (idim==1) then
      nx = sizeX
      is = istart
    else
      nx1 = (sizeX - 1) / idim + 1         ! local domain size
      in1 = sizeX - idim * (nx1-1)         ! 

      if( icoord .lt. in1 ) then
         nx = nx1
         is = istart + icoord * nx
      else
         nx = nx1 - 1
         is = istart + in1*nx1 + (icoord-in1)*nx
      end if
    end if
    nxx = nx + 2*nbx + 1
    ie  = is + nxx - 1

    ! same for y coordinate:
    if (jdim==1) then
      ny = sizeY
      js = jstart
    else
      ny1 = (sizeY - 1) / jdim + 1         ! local domain size
      jn1 = sizeY - jdim * (ny1-1)         ! 

      if( jcoord .lt. jn1 ) then
         ny = ny1
         js = jstart + jcoord * ny
      else
         ny = ny1 - 1
         js = jstart + jn1*ny1 + (jcoord-jn1)*ny
      end if
    end if
    nyy = ny + 2*nby + 1
    je  = js + nyy - 1

    ! serial z coordinate
    nz = sizeZ
    nzz = nz + 2*nbz + 1


    ! Some consistency check between domain size and cpus
    if ( sizeX < idim) then
       if( master ) then
          print *, 'Stopping! Do not use more than',&
               &     sizeX, &
               &    'cpus to compute this application' 
       end if
       error_flag = .true.
       return
    end if
    if ( sizeY < jdim) then
       if( master ) then
          print *, 'Stopping! Do not use more than',&
               &     sizeY, &
               &    'cpus to compute this application'      
       end if
       error_flag = .true.
       return
    end if


    ! display the decomposition 
    if( verboseMPI ) then
       print"(a,a,i3,a,i3,a,i3,a,i3,a,i3,a)","(rank,icoord,jcoord) -> (nx,ny)", &
            & "(", rank, ",", icoord, &
            &  ", ", jcoord, ") -> (", nx, ", ", ny, ")"
    end if

  end subroutine init_mpi


  !--------------------------------------------------------------------------


  function dot_product3D_glob(a,b)
    !---------------------------
    ! dot product for 3D arrays
    !---------------------------

    ! in/out variables
    real :: dot_product3D_glob
    real, dimension(:,:,:), intent(in) :: a,b

    ! locals
    integer, dimension(3) :: aSize, bSize
    integer :: i,j,k

    ! MPI stuff
    real :: dot_product3D_loc
    integer :: root

    aSize = shape(a)
    bSize = shape(b)

    do i = 1,3
       if( aSize(i) .ne. bSize(i) ) stop"dot_product3D failure."
    end do

    dot_product3D_loc = 0.0
    do k = 1, aSize(3)
       do j = 1, aSize(2)
          dot_product3D_loc = dot_product3D_loc + dot_product( a(:,j,k),b(:,j,k) )
       end do
    end do

    !MPI sum over all procs
    root = 0
    call mpi_reduce(dot_product3D_loc, dot_product3D_glob, 1, mpi_double_precision,&
         & mpi_sum, root, comm, ierror)

    call mpi_bcast(dot_product3D_glob, 1, mpi_double_precision, root, comm, ierror)


  end function dot_product3D_glob




end module mpi_module
