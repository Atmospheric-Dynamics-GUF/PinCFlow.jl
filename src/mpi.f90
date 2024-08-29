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

  private ! private module

  public :: init_mpi
  public :: setHalos
  public :: dot_product3D_glob
  public :: abort_message

  contains

  subroutine abort_message(message)
    character(len = *), intent(in) :: message

    if(master) then
      print *, message
    end if

    call mpi_finalize(ierror)
    stop

  end subroutine abort_message

  !----------------------------------------------------------------------------------

  subroutine setHalos(var, option)
    !-------------------------------
    !  set values in halo cells
    !-------------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var
    character(len = *), intent(in) :: option

    ! aux fields for "varTilde"
    real, dimension(ny, nz) :: x1SliceLeft_send, x1SliceRight_send
    real, dimension(ny, nz) :: x1SliceLeft_recv, x1SliceRight_recv

    real, dimension(nx, nz) :: y1SliceBack_send, y1SliceForw_send
    real, dimension(nx, nz) :: y1SliceBack_recv, y1SliceForw_recv

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: i, j, k, iVar, iTilde, nTilde
    integer :: i0, j0, k0

    select case(option)

    case("var")

      ! Set halos of basic variables.
      call setHalosOfField(var%rho)
      call setHalosOfField(var%u)
      call setHalosOfField(var%v)
      call setHalosOfField(var%w)
      call setHalosOfField(var%pi)
      call setHalosOfField(var%rhop)

      ! Set halos of mass-weighted potential temperature.
      if(model == "compressible") then
        call setHalosOfField(var%P)
      end if

      ! Set halos of ice variables.
      if(include_ice .and. updateIce) then
        do iVar = 1, 4
          call setHalosOfField(var%ICE(:, :, :, iVar))
        end do
      end if

    case("ice")

      do iVar = 1, 4
        call setHalosOfField(var%ICE(:, :, :, iVar))
      end do

    case("ice2")

      do iVar = 1, nVarIce
        call setHalosOfField(var%ICE2(:, :, :, iVar))
      end do

    case("tracer")

      call setHalosOfField(var%chi)

    case("varTilde")

      !-----------------------------------
      !             varTilde
      !-----------------------------------

      if(idim > 1) call mpi_cart_shift(comm, 0, 1, left, right, ierror)
      if(jdim > 1) call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

      !      achatz: this part seems to prohibit using less than two ghost cells.
      !      Can it be removed?

      ! if only for rhoTilde
      nTilde = 1
      ! if for another VarTilde set nTilde = 2
      !if ( model == "compressible" ) then
      !  nTilde = 2
      !end if

      do iTilde = 1, nTilde
        !------------------------------
        !          x-direction
        !------------------------------
        ! reconstructed density needed in ghost cell i = nx+2
        ! ...in ghost cell i = -1

        !if( updateMass ) then

        ! slice size
        sendcount = ny * nz
        recvcount = sendcount

        ! read slice into contiguous array
        if(iTilde == 1) then
          x1SliceLeft_send(1:ny, 1:nz) = rhoTilde(2, 1:ny, 1:nz, 1, 0)
          x1SliceRight_send(1:ny, 1:nz) = rhoTilde(nx - 1, 1:ny, 1:nz, 1, 1)
        else if(iTilde == 2) then
          x1SliceLeft_send(1:ny, 1:nz) = PTilde(2, 1:ny, 1:nz, 1, 0)
          x1SliceRight_send(1:ny, 1:nz) = PTilde(nx - 1, 1:ny, 1:nz, 1, 1)
        end if

        if(idim > 1) then

          ! left -> right
          source = left
          dest = right
          tag = 100

          call mpi_sendrecv(x1SliceRight_send(1, 1), sendcount, &
              &mpi_double_precision, dest, tag, x1SliceLeft_recv(1, 1), &
              &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              &sts_left, ierror)

          ! right -> left
          source = right
          dest = left
          tag = 100
          call mpi_sendrecv(x1SliceLeft_send(1, 1), sendcount, &
              &mpi_double_precision, dest, tag, x1SliceRight_recv(1, 1), &
              &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              &sts_right, ierror)

          ! write auxiliary slice to var field
          if(iTilde == 1) then ! rhoTilde
            ! right halos
            rhoTilde(nx + 2, 1:ny, 1:nz, 1, 0) = x1SliceRight_recv(1:ny, 1:nz)

            ! left halos
            rhoTilde(- 1, 1:ny, 1:nz, 1, 1) = x1SliceLeft_recv(1:ny, 1:nz)
          else if(iTilde == 2) then ! PTilde
            ! right halos
            PTilde(nx + 2, 1:ny, 1:nz, 1, 0) = x1SliceRight_recv(1:ny, 1:nz)

            ! left halos
            PTilde(- 1, 1:ny, 1:nz, 1, 1) = x1SliceLeft_recv(1:ny, 1:nz)
          end if

        end if

        if(verbose .and. master) print *, "horizontalHalos:  x-horizontal &
            &halos copied."

        !end if

        !------------------------------
        !          y-direction
        !------------------------------
        ! reconstructed density needed in ghost cell j = ny+2
        ! ...in ghost cell j = -1

        !if( updateMass ) then

        ! slice size
        sendcount = nx * nz
        recvcount = sendcount

        ! read slice into contiguous array
        if(iTilde == 1) then
          y1SliceBack_send(1:nx, 1:nz) = rhoTilde(1:nx, 2, 1:nz, 2, 0)
          y1SliceForw_send(1:nx, 1:nz) = rhoTilde(1:nx, ny - 1, 1:nz, 2, 1)
        else if(iTilde == 2) then
          y1SliceBack_send(1:nx, 1:nz) = PTilde(1:nx, 2, 1:nz, 2, 0)
          y1SliceForw_send(1:nx, 1:nz) = PTilde(1:nx, ny - 1, 1:nz, 2, 1)
        end if

        if(jdim > 1) then

          ! back -> forw
          source = back
          dest = forw
          tag = 100

          call mpi_sendrecv(y1SliceForw_send(1, 1), sendcount, &
              &mpi_double_precision, dest, tag, y1SliceBack_recv(1, 1), &
              &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              &sts_back, ierror)

          ! forw -> back
          source = forw
          dest = back
          tag = 100
          call mpi_sendrecv(y1SliceBack_send(1, 1), sendcount, &
              &mpi_double_precision, dest, tag, y1SliceForw_recv(1, 1), &
              &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              &sts_right, ierror)

          ! write auxiliary slice to var field

          if(iTilde == 1) then ! rhoTilde
            ! right halos
            rhoTilde(1:nx, ny + 2, 1:nz, 2, 0) = y1SliceForw_recv(1:nx, 1:nz)

            ! left halos
            rhoTilde(1:nx, - 1, 1:nz, 2, 1) = y1SliceBack_recv(1:nx, 1:nz)
          else if(iTilde == 2) then ! PTilde
            ! right halos
            PTilde(1:nx, ny + 2, 1:nz, 2, 0) = y1SliceForw_recv(1:nx, 1:nz)

            ! left halos
            PTilde(1:nx, - 1, 1:nz, 2, 1) = y1SliceBack_recv(1:nx, 1:nz)
          end if

        end if

        if(verbose .and. master) print *, "horizontalHalos:  y-horizontal &
            &halos copied."

        !end if ! updateMass
      end do

    case default

      if(master) then
        print *, "setHalo: Unknown case. Stop."
        call mpi_barrier(comm, ierror)
        call mpi_finalize(ierror)
        stop
      end if

    end select ! option "var", "varTilde", ...

  end subroutine setHalos

  !------------------------------------------------------------------

  subroutine init_mpi(error_flag)

    ! in/out vars
    logical, intent(out) :: error_flag

    include 'mpif.h'

    error_flag = .false.
    verboseMPI = .false.

    !-----------------------
    !    basic MPI info
    !-----------------------

    call mpi_init(ierror)
    call mpi_comm_rank(mpi_comm_world, rank, ierror)
    if(rank .eq. 0) then
      master = .true.
    else
      master = .false.
    end if
    call mpi_comm_size(mpi_comm_world, nbProc, ierror)

    !-----------------------------
    !    Domain decomposition
    !-----------------------------

    ! open the namelist by master
    if(master) then
      open(unit = 10, file = file_namelist, action = "read", form &
          &= "formatted", status = "old")

      ! Read domain namelist.
      rewind(unit = 10)
      read(unit = 10, nml = domain)
      close(10)

    end if

    ! broadcast the namelist values
    call mpi_bcast(sizeX, 1, mpi_integer, root, mpi_comm_world, ierror)
    call mpi_bcast(sizeY, 1, mpi_integer, root, mpi_comm_world, ierror)
    call mpi_bcast(sizeZ, 1, mpi_integer, root, mpi_comm_world, ierror)
    call mpi_bcast(nbx, 1, mpi_integer, root, mpi_comm_world, ierror)
    call mpi_bcast(nby, 1, mpi_integer, root, mpi_comm_world, ierror)
    call mpi_bcast(nbz, 1, mpi_integer, root, mpi_comm_world, ierror)
    call mpi_bcast(lx_dim, 2, mpi_double_precision, root, mpi_comm_world, &
        &ierror)
    call mpi_bcast(ly_dim, 2, mpi_double_precision, root, mpi_comm_world, &
        &ierror)
    call mpi_bcast(lz_dim, 2, mpi_double_precision, root, mpi_comm_world, &
        &ierror)
    call mpi_bcast(nprocx, 1, mpi_integer, root, mpi_comm_world, ierror)
    call mpi_bcast(nprocy, 1, mpi_integer, root, mpi_comm_world, ierror)

    ! domain size with ghost cells
    sizeXX = sizeX + 2 * nbx + 1
    sizeYY = sizeY + 2 * nby + 1

    iStart = - nbx + 1 ! left start index of global domain
    jStart = - nby + 1 ! back start index of global domain

    ! set up no. of cpus in x-y direction
    idim = nprocx
    jdim = nprocy
    if(idim * jdim /= nbProc) then
      if(master) print *, 'Virtual topology wrong idim*jdim should be equal to &
          &nbProc!'
      error_flag = .true.
      return
    end if

    dims(1) = idim
    dims(2) = jdim
    periods(1) = .true.
    periods(2) = .true.

    call mpi_cart_create(mpi_comm_world, 2, dims, periods, .true., comm, ierror)
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
    if(idim == 1) then
      nx = sizeX
      is = istart
    else
      nx1 = (sizeX - 1) / idim + 1 ! local domain size
      in1 = sizeX - idim * (nx1 - 1) !

      if(icoord .lt. in1) then
        nx = nx1
        is = istart + icoord * nx
      else
        nx = nx1 - 1
        is = istart + in1 * nx1 + (icoord - in1) * nx
      end if
    end if
    nxx = nx + 2 * nbx + 1
    ie = is + nxx - 1

    ! same for y coordinate:
    if(jdim == 1) then
      ny = sizeY
      js = jstart
    else
      ny1 = (sizeY - 1) / jdim + 1 ! local domain size
      jn1 = sizeY - jdim * (ny1 - 1) !

      if(jcoord .lt. jn1) then
        ny = ny1
        js = jstart + jcoord * ny
      else
        ny = ny1 - 1
        js = jstart + jn1 * ny1 + (jcoord - jn1) * ny
      end if
    end if
    nyy = ny + 2 * nby + 1
    je = js + nyy - 1

    ! serial z coordinate
    nz = sizeZ
    nzz = nz + 2 * nbz + 1

    ! Some consistency check between domain size and cpus
    if(sizeX < idim) then
      if(master) then
        print *, 'Stopping! Do not use more than', sizeX, 'cpus to compute &
            &this application'
      end if
      error_flag = .true.
      return
    end if
    if(sizeY < jdim) then
      if(master) then
        print *, 'Stopping! Do not use more than', sizeY, 'cpus to compute &
            &this application'
      end if
      error_flag = .true.
      return
    end if

    ! Number of halo cells should not be larger than number of physical cells.
    if(idim > 1 .and. nbx > nx) then
      if(master) stop "Error in init_mpi: idim > 1 and nbx > nx!"
    end if
    if(jdim > 1 .and. nby > ny) then
      if(master) stop "Error in init_mpi: jdim > 1 and nby > ny!"
    end if

    ! display the decomposition
    if(verboseMPI) then
      print "(a,a,i3,a,i3,a,i3,a,i3,a,i3,a)", "(rank,icoord,jcoord) &
          &-> (nx,ny)", "(", rank, ",", icoord, ", ", jcoord, ") -> (", nx, ", &
          &", ny, ")"
    end if

  end subroutine init_mpi

  !------------------------------------------------------------------

  function dot_product3D_glob(a, b)
    !---------------------------
    ! dot product for 3D arrays
    !---------------------------

    ! in/out variables
    real :: dot_product3D_glob
    real, dimension(:, :, :), intent(in) :: a, b

    ! locals
    integer, dimension(3) :: aSize, bSize
    integer :: i, j, k

    ! MPI stuff
    real :: dot_product3D_loc

    aSize = shape(a)
    bSize = shape(b)

    do i = 1, 3
      if(aSize(i) .ne. bSize(i)) stop "dot_product3D failure."
    end do

    dot_product3D_loc = 0.0
    do k = 1, aSize(3)
      do j = 1, aSize(2)
        dot_product3D_loc = dot_product3D_loc + dot_product(a(:, j, k), b(:, &
            &j, k))
      end do
    end do

    !MPI sum over all procs
    call mpi_reduce(dot_product3D_loc, dot_product3D_glob, 1, &
        &mpi_double_precision, mpi_sum, root, comm, ierror)

    call mpi_bcast(dot_product3D_glob, 1, mpi_double_precision, root, comm, &
        &ierror)

  end function dot_product3D_glob

end module mpi_module
