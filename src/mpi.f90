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
  use mpi

  implicit none

  private ! private module

  public :: init_mpi
  public :: setHalos
  public :: setHalosOfField
  public :: setHalosOfField2D
  public :: dot_product3D_glob

  contains

  subroutine setHalos(var, option)
    !-------------------------------
    !  set values in halo cells
    !-------------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var
    character(len = *), intent(in) :: option

    integer :: iVar

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

    case("ice")

      do iVar = 1, nVarIce
        call setHalosOfField(var%ICE(:, :, :, iVar))
      end do

    case("tracer")

      call setHalosOfField(var%chi)

    case default

      if(master) then
        print *, "setHalo: Unknown case. Stop."
        call mpi_barrier(comm, ierror)
        call mpi_finalize(ierror)
        stop
      end if

    end select

  end subroutine setHalos

  !------------------------------------------------------------------

  subroutine setHalosOfField(field)

    !-------------------------------
    !  set values in halo cells
    !-------------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: field

    ! auxiliary fields
    real, dimension(nbx, - nby:ny + nby, - nbz:nz + nbz) :: xSliceLeft_send, &
        &xSliceRight_send
    real, dimension(nbx, - nby:ny + nby, - nbz:nz + nbz) :: xSliceLeft_recv, &
        &xSliceRight_recv

    real, dimension(- nbx:nx + nbx, nby, - nbz:nz + nbz) :: ySliceBack_send, &
        &ySliceForw_send
    real, dimension(- nbx:nx + nbx, nby, - nbz:nz + nbz) :: ySliceBack_recv, &
        &ySliceForw_recv

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: i, j, k
    integer :: i0, j0, k0

    !-----------------------------
    !     find neighbour procs
    !-----------------------------

    if(idim > 1) call mpi_cart_shift(comm, 0, 1, left, right, ierror)
    if(jdim > 1) call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

    !------------------------------
    !          x-direction
    !------------------------------

    if(idim > 1) then

      ! slice size
      sendcount = nbx * (ny + 2 * nby + 1) * (nz + 2 * nbz + 1)
      recvcount = sendcount

      ! read slice into contiguous array
      do i = 1, nbx
        xSliceLeft_send(i, :, :) = field(i, :, :)
        xSliceRight_send(i, :, :) = field(nx - nbx + i, :, :)
      end do

      ! left -> right
      source = left
      dest = right
      tag = 100

      i0 = 1; j0 = - nby; k0 = - nbz

      call mpi_sendrecv(xSliceRight_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, xSliceLeft_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xSliceLeft_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, xSliceRight_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_right, ierror)

      ! write auxiliary slice to var field
      do i = 1, nbx
        ! right halos
        field(nx + i, :, :) = xSliceRight_recv(i, :, :)
        ! left halos
        field(- nbx + i, :, :) = xSliceLeft_recv(i, :, :)
      end do

    else

      do i = 1, nbx
        field(nx + i, :, :) = field(i, :, :)
        field(- i + 1, :, :) = field(nx - i + 1, :, :)
      end do

    end if

    !------------------------------
    !          y-direction
    !------------------------------

    if(jdim > 1) then

      ! slice size
      sendcount = nby * (nx + 2 * nbx + 1) * (nz + 2 * nbz + 1)
      recvcount = sendcount

      ! read slice into contiguous array
      do j = 1, nby
        ySliceBack_send(:, j, :) = field(:, j, :)
        ySliceForw_send(:, j, :) = field(:, ny - nby + j, :)
      end do

      ! back -> forw
      source = back
      dest = forw
      tag = 100

      i0 = - nbx; j0 = 1; k0 = - nbz

      call mpi_sendrecv(ySliceForw_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, ySliceBack_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(ySliceBack_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, ySliceForw_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_forw, ierror)

      ! write auxiliary slice to var field
      do j = 1, nby
        ! right halos
        field(:, ny + j, :) = ySliceForw_recv(:, j, :)
        ! left halos
        field(:, - nby + j, :) = ySliceBack_recv(:, j, :)
      end do

    else

      do j = 1, nby
        field(:, ny + j, :) = field(:, j, :)
        field(:, - j + 1, :) = field(:, ny - j + 1, :)
      end do

    end if

  end subroutine setHalosOfField

  !------------------------------------------------------------------

  subroutine setHalosOfField2D(field)

    ! Subroutine needed for halos of topography.

    !-------------------------------
    !  set values in halo cells
    !-------------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby), intent(inout) :: field

    ! auxiliary fields
    real, dimension(nbx, - nby:ny + nby) :: xSliceLeft_send, xSliceRight_send
    real, dimension(nbx, - nby:ny + nby) :: xSliceLeft_recv, xSliceRight_recv
    real, dimension(- nbx:nx + nbx, nby) :: ySliceBack_send, ySliceForw_send
    real, dimension(- nbx:nx + nbx, nby) :: ySliceBack_recv, ySliceForw_recv

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: i, j, k
    integer :: i0, j0, k0

    !-----------------------------
    !     find neighbour procs
    !-----------------------------

    if(idim > 1) call mpi_cart_shift(comm, 0, 1, left, right, ierror)
    if(jdim > 1) call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

    !------------------------------
    !          x-direction
    !------------------------------

    if(idim > 1) then

      ! slice size
      sendcount = nbx * (ny + 2 * nby + 1)
      recvcount = sendcount

      ! read slice into contiguous array
      do i = 1, nbx
        xSliceLeft_send(i, :) = field(i, :)
        xSliceRight_send(i, :) = field(nx - nbx + i, :)
      end do

      ! left -> right
      source = left
      dest = right
      tag = 100

      i0 = 1; j0 = - nby

      call mpi_sendrecv(xSliceRight_send(i0, j0), sendcount, &
          &mpi_double_precision, dest, tag, xSliceLeft_recv(i0, j0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xSliceLeft_send(i0, j0), sendcount, &
          &mpi_double_precision, dest, tag, xSliceRight_recv(i0, j0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_right, ierror)

      ! write auxiliary slice to var field
      do i = 1, nbx
        ! right halos
        field(nx + i, :) = xSliceRight_recv(i, :)
        ! left halos
        field(- nbx + i, :) = xSliceLeft_recv(i, :)
      end do

    else

      do i = 1, nbx
        field(nx + i, :) = field(i, :)
        field(- i + 1, :) = field(nx - i + 1, :)
      end do

    end if

    !------------------------------
    !          y-direction
    !------------------------------

    if(jdim > 1) then

      ! slice size
      sendcount = nby * (nx + 2 * nbx + 1)
      recvcount = sendcount

      ! read slice into contiguous array
      do j = 1, nby
        ySliceBack_send(:, j) = field(:, j)
        ySliceForw_send(:, j) = field(:, ny - nby + j)
      end do

      ! back -> forw
      source = back
      dest = forw
      tag = 100

      i0 = - nbx; j0 = 1

      call mpi_sendrecv(ySliceForw_send(i0, j0), sendcount, &
          &mpi_double_precision, dest, tag, ySliceBack_recv(i0, j0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(ySliceBack_send(i0, j0), sendcount, &
          &mpi_double_precision, dest, tag, ySliceForw_recv(i0, j0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_forw, ierror)

      ! write auxiliary slice to var field
      do j = 1, nby
        ! right halos
        field(:, ny + j) = ySliceForw_recv(:, j)
        ! left halos
        field(:, - nby + j) = ySliceBack_recv(:, j)
      end do

    else

      do j = 1, nby
        field(:, ny + j) = field(:, j)
        field(:, - j + 1) = field(:, ny - j + 1)
      end do

    end if

  end subroutine setHalosOfField2D

  !------------------------------------------------------------------

  subroutine init_mpi(error_flag)

    ! in/out vars
    logical, intent(out) :: error_flag

    !include 'mpif.h'

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
      if(master) then
        print "(a)", "Virtual topology wrong, idim * jdim should be equal to &
            &nbProc!"
        print "(a)", "[idim, jdim] = [" // trim_integer(idim) // ", " &
            &// trim_integer(jdim), "]"
        print "(a)", "nbProc = " // trim_integer(nbProc)
      end if
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
