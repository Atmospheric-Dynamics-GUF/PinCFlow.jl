module bicgstab_tools_module

  use type_module
  use mpi_module ! modified by Junhong Wei (20161106)
  use timeScheme_module
  use atmosphere_module
  use algebra_module

  implicit none

  private
  ! all module objects (variables, functions, subroutines)
  ! are internal to the module by default

  !------------------------
  !   public subroutines
  !------------------------
  ! public :: SetUpHypre
  public :: SetUpBiCGStab
  public :: CleanUpBiCGStab
  ! public :: CleanUpHypre
  ! public :: test_hypre

  !------------------------
  !   private subroutines
  !------------------------
  !private ::

  contains

  !----------------------------------------------------------------------

  subroutine SetUpBICGStab
    ! --------------------------------------
    !    generate arrays for BiCGStab
    !---------------------------------------

    integer :: allocstat

    allocate(ac_b(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:alloc failed"

    !UAB
    allocate(acv_b(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:alloc failed"

    allocate(ach_b(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:alloc failed"
    !UAE

    allocate(al_b(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:alloc failed"

    allocate(ar_b(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:alloc failed"

    allocate(ab_b(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:alloc failed"

    allocate(af_b(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:alloc failed"

    allocate(ad_b(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:alloc failed"

    allocate(au_b(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:alloc failed"

    if(timeScheme == "semiimplicit") then
      allocate(alb_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(alf_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(arb_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(arf_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"
    end if

    ! TFC FJ
    ! Allocate additional matrix elements for TFC.
    if(topography) then
      allocate(aru_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(ard_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(alu_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(ald_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(afu_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(afd_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(abu_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(abd_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(auu_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(add_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(aruu_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(ardd_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(aluu_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(aldd_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(afuu_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(afdd_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(abuu_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"

      allocate(abdd_b(1:nx, 1:ny, 1:nz), stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:alloc failed"
    end if

    return

  end subroutine SetUpBICGStab

  !----------------------------------------------------------------------

  subroutine CleanUpBICGStab
    ! --------------------------------------
    !    destroy arrays for BiCGStab
    !---------------------------------------

    integer :: allocstat

    deallocate(ac_b, stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:dealloc failed"

    !UAB
    deallocate(acv_b, stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:dealloc failed"

    deallocate(ach_b, stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:dealloc failed"
    !UAE

    deallocate(al_b, stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:dealloc failed"

    deallocate(ar_b, stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:dealloc failed"

    deallocate(ab_b, stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:dealloc failed"

    deallocate(af_b, stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:dealloc failed"

    deallocate(ad_b, stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:dealloc failed"

    deallocate(au_b, stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:dealloc failed"

    if(timeScheme == "semiimplicit") then
      deallocate(alb_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(alf_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(arb_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(arf_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"
    end if

    ! TFC FJ
    ! Deallocate additional matrix elements for TFC.
    if(topography) then
      deallocate(aru_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(ard_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(alu_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(ald_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(afu_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(afd_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(abu_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(abd_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(auu_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(add_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(aruu_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(ardd_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(aluu_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(aldd_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(afuu_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(afdd_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(abuu_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"

      deallocate(abdd_b, stat = allocstat)
      if(allocstat /= 0) stop "BiCGStab:dealloc failed"
    end if

    return

  end subroutine CleanUpBICGStab

  !----------------------------------------------------------------------

  	 !   subroutine SetUpHypre
  	 !     ! --------------------------------------
  	 !     !    HYPRE
  	 !     !---------------------------------------

  	 !     ! variables due to the use of HYPRE

  	 !     integer :: ndim_hypre
  	 !     parameter (ndim_hypre = 3)

  	 !     integer ierr_hypre
  	 !     integer npts_x_periodic_hypre, npts_y_periodic_hypre, npts_z_periodic_hypre
  	 !     integer ihp
  	 !     integer offsets_e(7, 3)
  	 !     integer offsets_i(11, 3)
  	 !     integer :: nvalues_e
  	 !     integer :: nvalues_i
  	 !     integer :: allocstat
  	 !     integer :: index_count_hypre
  	 !     integer :: k_sing
  	 !     integer :: i0, j0

  	 !     ! local variables
  	 !     integer :: i, j, k

  	 !     real :: pStratU, pStratD, rhoEdge
  	 !     real :: AL, AR, AB, AF, AD, AU, AC
  	 !     real :: dx2, dy2, dz2
  	 !     real :: atol

  	 !     ! Local parameters
  	 !     integer :: maxIt

  	 !     ! verbose
  	 !     logical, parameter :: giveInfo = .true.

  	 !     if (timeScheme == "semiimplicit") then
  	 !       nvalues_i = nx * ny * nz * ne_hypre_i
  	 !     else
  	 !       nvalues_e = nx * ny * nz * ne_hypre_e
  	 !     end if

  	 !     allocate (bvalue_vector_hypre(1:(nx * ny * nz)), stat = allocstat)
  	 !     if (allocstat /= 0) stop "hypre:alloc failed"

  	 !     allocate (xvalue_vector_hypre(1:(nx * ny * nz)), stat = allocstat)
  	 !     if (allocstat /= 0) stop "hypre:alloc failed"

  	 !     if (timeScheme == "semiimplicit") then
  	 !       allocate (values_i(1:nvalues_i), stat = allocstat)
  	 !       if (allocstat /= 0) stop "hypre:alloc failed"

  	 !       allocate (stencil_indices_i(1:ne_hypre_i), stat = allocstat)
  	 !       if (allocstat /= 0) stop "hypre:alloc failed"

  	 !       stencil_indices_i(:) = (/0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10/)
  	 !     else
  	 !       allocate (values_e(1:nvalues_e), stat = allocstat)
  	 !       if (allocstat /= 0) stop "hypre:alloc failed"

  	 !       allocate (stencil_indices_e(1:ne_hypre_e), stat = allocstat)
  	 !       if (allocstat /= 0) stop "hypre:alloc failed"

  	 !       stencil_indices_e(:) = (/0, 1, 2, 3, 4, 5, 6/)
  	 !     end if

  	 !     ! Set parameters
  	 !     maxIt = maxIterPoisson

  	 !     !---------
  	 !     ! step one
  	 !     ! Set up a grid. Each processor describes the piece of the grid that
  	 !     ! it owns.

  	 !     ! Create an empty 2D grid object
  	 !     call HYPRE_StructGridCreate(mpi_comm_world, ndim_hypre, grid_hypre, &
  	 !         ierr_hypre)

  	 !     !----------------
  	 !     !  x-Boundary
  	 !     !----------------
  	 !     select case (xBoundary)

  	 !     case ("periodic")
  	 !       npts_x_periodic_hypre = sizeX
  	 !     case default
  	 !       stop "pressureBoundaryCondition: unknown case xBoundary."
  	 !     end select

  	 !     !----------------
  	 !     !   y-Boundary
  	 !     !----------------
  	 !     select case (yBoundary)

  	 !     case ("periodic")
  	 !       npts_y_periodic_hypre = sizeY
  	 !     case default
  	 !       stop "pressureBoundaryCondition: unknown case yBoundary."
  	 !     end select

  	 !     !----------------
  	 !     !   z-Boundary
  	 !     !----------------
  	 !     select case (zBoundary)

  	 !     case ("periodic")
  	 !       npts_z_periodic_hypre = sizeZ
  	 !     case ("solid_wall")
  	 !       npts_z_periodic_hypre = 0
  	 !     case default
  	 !       stop "pressureBoundaryCondition: unknown case zBoundary."
  	 !     end select

  	 !     ! Add boxes to the grid
  	 !     call HYPRE_StructGridSetExtents(grid_hypre, (/(icoord) * nx + 1, (jcoord) &
  	 !         * ny + 1, 1/), (/(icoord + 1) * nx, (jcoord + 1) * ny, nz/), ierr_hypre)

  	 !     call HYPRE_StructGridSetPeriodic(grid_hypre, (/npts_x_periodic_hypre, &
  	 !         npts_y_periodic_hypre, npts_z_periodic_hypre/), ierr_hypre)

  	 !     ! This is a collective call finalizing the grid assembly.
  	 !     ! The grid is now ``ready to be used''
  	 !     call HYPRE_StructGridAssemble(grid_hypre, ierr_hypre)

  	 !     !---------
  	 !     ! step two
  	 !     ! Define the discretization stencil

  	 !     if (timeScheme == "semiimplicit") then
  	 !       ! Create an empty 3D, 11-pt stencil object
  	 !       call HYPRE_StructStencilCreate(3, 11, stencil_i, ierr_hypre)
  	 !     else
  	 !       ! Create an empty 3D, 7-pt stencil object
  	 !       call HYPRE_StructStencilCreate(3, 7, stencil_e, ierr_hypre)
  	 !     end if

  	 !     ! Define the geometry of the stencils. Each represents a relative
  	 !     ! offset (in the index space).

  	 !     if (timeScheme == "semiimplicit") then
  	 !       offsets_i(1, :) = (/0, 0, 0/)
  	 !       offsets_i(2, :) = (/- 1, 0, 0/)
  	 !       offsets_i(3, :) = (/1, 0, 0/)
  	 !       offsets_i(4, :) = (/0, - 1, 0/)
  	 !       offsets_i(5, :) = (/0, 1, 0/)
  	 !       offsets_i(6, :) = (/0, 0, - 1/)
  	 !       offsets_i(7, :) = (/0, 0, 1/)
  	 !       offsets_i(8, :) = (/- 1, - 1, 0/)
  	 !       offsets_i(9, :) = (/- 1, 1, 0/)
  	 !       offsets_i(10, :) = (/1, - 1, 0/)
  	 !       offsets_i(11, :) = (/1, 1, 0/)

  	 !       do ihp = 0, 10
  	 !         call HYPRE_StructStencilSetElement(stencil_i, ihp, offsets_i((ihp &
  	 !             + 1), :), ierr_hypre)
  	 !       end do
  	 !     else
  	 !       offsets_e(1, :) = (/0, 0, 0/)
  	 !       offsets_e(2, :) = (/- 1, 0, 0/)
  	 !       offsets_e(3, :) = (/1, 0, 0/)
  	 !       offsets_e(4, :) = (/0, - 1, 0/)
  	 !       offsets_e(5, :) = (/0, 1, 0/)
  	 !       offsets_e(6, :) = (/0, 0, - 1/)
  	 !       offsets_e(7, :) = (/0, 0, 1/)

  	 !       do ihp = 0, 6
  	 !         call HYPRE_StructStencilSetElement(stencil_e, ihp, offsets_e((ihp &
  	 !             + 1), :), ierr_hypre)
  	 !       end do
  	 !     end if

  	 !     !---------
  	 !     ! step three
  	 !     ! Set up a Struct Matrix

  	 !     ! Create an empty matrix object

  	 !     if (timeScheme == "semiimplicit") then
  	 !       call HYPRE_StructMatrixCreate(mpi_comm_world, grid_hypre, stencil_i, &
  	 !           A_hp_i, ierr_hypre)
  	 !     else
  	 !       call HYPRE_StructMatrixCreate(mpi_comm_world, grid_hypre, stencil_e, &
  	 !           A_hp_e, ierr_hypre)
  	 !     end if

  	 !     ! Indicate that the matrix coefficients are ready to be set

  	 !     if (timeScheme == "semiimplicit") then
  	 !       call HYPRE_StructMatrixInitialize(A_hp_i, ierr_hypre)
  	 !     else
  	 !       call HYPRE_StructMatrixInitialize(A_hp_e, ierr_hypre)
  	 !     end if

  	 !     !---------
  	 !     ! step four
  	 !     ! Set up Struct Vectors for b and x.  Each processor sets the vectors
  	 !     ! corresponding to its boxes.

  	 !     ! Create an empty vector object

  	 !     if (timeScheme == "semiimplicit") then
  	 !       call HYPRE_StructVectorCreate(mpi_comm_world, grid_hypre, b_hp_i, &
  	 !           ierr_hypre)
  	 !       call HYPRE_StructVectorCreate(mpi_comm_world, grid_hypre, x_hp_i, &
  	 !           ierr_hypre)
  	 !     else
  	 !       call HYPRE_StructVectorCreate(mpi_comm_world, grid_hypre, b_hp_e, &
  	 !           ierr_hypre)
  	 !       call HYPRE_StructVectorCreate(mpi_comm_world, grid_hypre, x_hp_e, &
  	 !           ierr_hypre)
  	 !     end if

  	 !     ! Indicate that the vector coefficients are ready to be set

  	 !     if (timeScheme == "semiimplicit") then
  	 !       call HYPRE_StructVectorInitialize(b_hp_i, ierr_hypre)
  	 !       call HYPRE_StructVectorInitialize(x_hp_i, ierr_hypre)
  	 !     else
  	 !       call HYPRE_StructVectorInitialize(b_hp_e, ierr_hypre)
  	 !       call HYPRE_StructVectorInitialize(x_hp_e, ierr_hypre)
  	 !     end if

  	 !     !---------
  	 !     ! the commentary "!!*" stand for what I also plan to include here,
  	 !     ! but first other check
  	 !     ! step five
  	 !     ! Set up and use a solver (See the Reference Manual for descriptions
  	 !     ! of all of the options.)

  	 !     ! Create an empty Hybrid Struct solver

  	 !     if (timeScheme == "semiimplicit") then
  	 !       call HYPRE_StructHybridCreate(mpi_comm_world, solver_hp_i, ierr_hypre)
  	 !     else
  	 !       call HYPRE_StructHybridCreate(mpi_comm_world, solver_hp_e, ierr_hypre)
  	 !     end if

  	 !     ! Set print level (set 0 for complete silence or 2 for info on
  	 !     ! iterations)

  	 !     if (timeScheme == "semiimplicit") then
  	 !       call HYPRE_StructHybridSetPrintLevel(solver_hp_i, 0, ierr_hypre)
  	 !     else
  	 !       call HYPRE_StructHybridSetPrintLevel(solver_hp_e, 0, ierr_hypre)
  	 !     end if

  	 !     ! Set the type of Krylov solver to use.
  	 !     ! Current krylov methods set by solver type are:
  	 !     ! 1 DSCG (default)
  	 !     ! 2 GMRES
  	 !     ! 3 BiCGSTAB
  	 !     ! Junhong used 2

  	 !     if (timeScheme == "semiimplicit") then
  	 !       call HYPRE_StructHybridSetSolverType(solver_hp_i, 3, ierr_hypre)
  	 !     else
  	 !       call HYPRE_StructHybridSetSolverType(solver_hp_e, 3, ierr_hypre)
  	 !     end if

  	 !     ! Set tolerance value for relative! residual (unlike in Rieper's BiGSTAB
  	 !     ! where the absolute! residual is used).
  	 !     ! tol = (tolPoisson in the namelist)/tolref
  	 !     !!*  call HYPRE_StructHybridSetTol(solver_hp_e, tol, ierr_hypre)

  	 !     ! Uncomment the line below to have the tolerance for absolute! residual
  	 !     ! (it works only for BiCGSTAB, i.e. solver_type = 3)
  	 !     !call HYPRE_StructHybridSetStopCrit(solver_hp_e, 0, ierr_hypre)

  	 !     ! Set convergence criterion for activating a preconditioner.
  	 !     ! A number between 0 and 1 is accepted, the smaller the number, the
  	 !     ! sooner a preconditioner is activated.

  	 !     if (timeScheme == "semiimplicit") then
  	 !       call HYPRE_StructHybridSetConvergenc(solver_hp_i, tolCond, ierr_hypre)
  	 !     else
  	 !       call HYPRE_StructHybridSetConvergenc(solver_hp_e, tolCond, ierr_hypre)
  	 !     end if

  	 !     ! Set maximum number of iterations. maxIt = maxIterPoisson in the
  	 !     ! namelist

  	 !     if (timeScheme == "semiimplicit") then
  	 !       call HYPRE_StructHybridSetPCGMaxIter(solver_hp_i, maxIt, ierr_hypre)
  	 !     else
  	 !       call HYPRE_StructHybridSetPCGMaxIter(solver_hp_e, maxIt, ierr_hypre)
  	 !     end if

  	 !     ! Setup and solve

  	 !     if (timeScheme == "semiimplicit") then
  	 !       call HYPRE_StructHybridSetup(solver_hp_i, A_hp_i, b_hp_i, x_hp_i, &
  	 !           ierr_hypre)
  	 !     else
  	 !       call HYPRE_StructHybridSetup(solver_hp_e, A_hp_e, b_hp_e, x_hp_e, &
  	 !           ierr_hypre)
  	 !     end if

  	 !     return

  	 !   end subroutine SetUpHypre

  !-----------------------------------------------------------------------

  	 !   subroutine CleanUpHypre
  	 !     ! --------------------------------------
  	 !     !    HYPRE
  	 !     !---------------------------------------

  	 !     ! in/out variables

  	 !     ! variables due to the use of HYPRE

  	 !     integer :: ndim_hypre
  	 !     parameter (ndim_hypre = 3)

  	 !     integer ierr_hypre
  	 !     integer npts_x_periodic_hypre, npts_y_periodic_hypre, npts_z_periodic_hypre

  	 !     integer ihp
  	 !     integer offsets_e(7, 3)
  	 !     integer offsets_i(11, 3)

  	 !     integer :: nvalues_e
  	 !     integer :: nvalues_i

  	 !     integer :: allocstat

  	 !     integer :: index_count_hypre

  	 !     integer :: k_sing
  	 !     integer :: i0, j0

  	 !     ! local variables
  	 !     integer :: i, j, k
  	 !     real :: pStratU, pStratD, rhoEdge
  	 !     real :: AL, AR, AB, AF, AD, AU, AC

  	 !     real :: dx2, dy2, dz2
  	 !     real :: atol

  	 !     ! Local parameters
  	 !     integer :: maxIt

  	 !     ! verbose
  	 !     logical, parameter :: giveInfo = .true.

  	 !     real :: sol_mean_hypre

  	 !     ! Free memory

  	 !     call HYPRE_StructGridDestroy(grid_hypre, ierr_hypre)

  	 !     if (timeScheme == "semiimplicit") then
  	 !       call HYPRE_StructStencilDestroy(stencil_i, ierr_hypre)
  	 !       call HYPRE_StructMatrixDestroy(A_hp_i, ierr_hypre)
  	 !       call HYPRE_StructVectorDestroy(b_hp_i, ierr_hypre)
  	 !       call HYPRE_StructVectorDestroy(x_hp_i, ierr_hypre)
  	 !       call HYPRE_StructHybridDestroy(solver_hp_i, ierr_hypre)
  	 !     else
  	 !       call HYPRE_StructStencilDestroy(stencil_e, ierr_hypre)
  	 !       call HYPRE_StructMatrixDestroy(A_hp_e, ierr_hypre)
  	 !       call HYPRE_StructVectorDestroy(b_hp_e, ierr_hypre)
  	 !       call HYPRE_StructVectorDestroy(x_hp_e, ierr_hypre)
  	 !       call HYPRE_StructHybridDestroy(solver_hp_e, ierr_hypre)
  	 !     end if

  	 !     ! deallocate local fields

  	 !     deallocate (bvalue_vector_hypre, stat = allocstat); if (allocstat /= 0) &
  	 !         stop "hypre:dealloc failed"
  	 !     deallocate (xvalue_vector_hypre, stat = allocstat); if (allocstat /= 0) &
  	 !         stop "hypre:dealloc failed"

  	 !     if (timeScheme == "semiimplicit") then
  	 !       deallocate (values_i, stat = allocstat)
  	 !       if (allocstat /= 0) stop "hypre:dealloc failed"
  	 !     else
  	 !       deallocate (values_e, stat = allocstat)
  	 !       if (allocstat /= 0) stop "hypre:dealloc failed"
  	 !     end if

  	 !     return

  	 !   end subroutine CleanUpHypre

  !-------------------------------------------------------------------------

  	 !   subroutine test_hypre

  	 !     ! --------------------------------------
  	 !     !    test HYPRE in simple setting
  	 !     !---------------------------------------

  	 !     ! variables due to the use of HYPRE

  	 !     integer :: ndim_t
  	 !     parameter (ndim_t = 3)

  	 !     integer, parameter :: ne_hypre_t = 11
  	 !     !integer, parameter :: ne_hypre_t = 7

  	 !     integer * 8 grid_t, stencil_t, A_hp_t, b_hp_t, x_hp_t, solver_hp_t

  	 !     integer ierr_t
  	 !     integer npts_x_periodic_t, npts_y_periodic_t, npts_z_periodic_t
  	 !     integer ihp
  	 !     integer offsets_t(11, 3)
  	 !     !integer offsets_t(7,3)
  	 !     integer :: nvalues_t
  	 !     integer :: allocstat
  	 !     integer :: index_count_t

  	 !     ! local variables
  	 !     integer :: i, j, k

  	 !     real, dimension (:), allocatable :: values_t
  	 !     real, dimension (:), allocatable :: bvalue_vector_t, xvalue_vector_t

  	 !     integer, dimension (:), allocatable :: stencil_indices_t

  	 !     real, dimension (1:nx, 1:ny, 1:nz) :: b ! RHS
  	 !     real :: bmean
  	 !     real, dimension (1:nx, 1:ny, 1:nz) :: sol
  	 !     real :: res ! residual
  	 !     real :: tol_t

  	 !     real :: AC, AL, AR, AB, AF, AD, AU, ALB, ALF, ARB, ARF

  	 !     real :: solmean, xmean

  	 !     integer :: nIter
  	 !     logical :: errFlag

  	 !     ! define test solution

  	 !     solmean = 0.0

  	 !     do k = 1, nz
  	 !       do j = 1, ny
  	 !         do i = 1, nx
  	 !           sol(i, j, k) = sin(real (i * j * k))
  	 !           solmean = solmean + sol(i, j, k)

  	 !           print *, i, j, k, 'sol =', sol(i, j, k)
  	 !         end do ! i_loop
  	 !       end do ! j_loop
  	 !     end do ! k_loop

  	 !     solmean = solmean / (nx * ny * nz)

  	 !     sol(:, :, :) = sol(:, :, :) - solmean

  	 !     nvalues_t = nx * ny * nz * ne_hypre_t

  	 !     allocate (bvalue_vector_t(1:(nx * ny * nz)), stat = allocstat)
  	 !     if (allocstat /= 0) stop "hypre:alloc failed"

  	 !     allocate (xvalue_vector_t(1:(nx * ny * nz)), stat = allocstat)
  	 !     if (allocstat /= 0) stop "hypre:alloc failed"

  	 !     allocate (values_t(1:nvalues_t), stat = allocstat)
  	 !     if (allocstat /= 0) stop "hypre:alloc failed"

  	 !     allocate (stencil_indices_t(1:ne_hypre_t), stat = allocstat)
  	 !     if (allocstat /= 0) stop "hypre:alloc failed"

  	 !     stencil_indices_t(:) = (/0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10/)
  	 !     !stencil_indices_t(:) = (/0,1,2,3,4,5,6/)

  	 !     !---------
  	 !     ! step one
  	 !     ! Set up a grid. Each processor describes the piece of the grid that
  	 !     ! it owns.

  	 !     ! Create an empty 2D grid object
  	 !     call HYPRE_StructGridCreate(mpi_comm_world, ndim_t, grid_t, ierr_t)

  	 !     ! implement periodic BC in x and y

  	 !     !----------------
  	 !     !  x-Boundary
  	 !     !----------------

  	 !     npts_x_periodic_t = sizeX

  	 !     !----------------
  	 !     !   y-Boundary
  	 !     !----------------

  	 !     npts_y_periodic_t = sizeY

  	 !     !----------------
  	 !     !   z-Boundary
  	 !     !----------------

  	 !     npts_z_periodic_t = 0

  	 !     call HYPRE_StructGridSetPeriodic(grid_t, (/npts_x_periodic_t, &
  	 !         npts_y_periodic_t, npts_z_periodic_t/), ierr_t)

  	 !     ! Add boxes to the grid
  	 !     call HYPRE_StructGridSetExtents(grid_t, (/(icoord) * nx + 1, (jcoord) * ny &
  	 !         + 1, 1/), (/(icoord + 1) * nx, (jcoord + 1) * ny, nz/), ierr_t)

  	 !     !! Add boxes to the grid
  	 !     !call HYPRE_StructGridSetExtents(grid_t, &
  	 !     !        & (/ (icoord-1)*nx, (jcoord-1)*ny , 0 /), &
  	 !     !        & (/ icoord*nx-1 , jcoord*ny-1 , nz-1 /), ierr_t)

  	 !     ! This is a collective call finalizing the grid assembly.
  	 !     ! The grid is now ``ready to be used''
  	 !     call HYPRE_StructGridAssemble(grid_t, ierr_t)

  	 !     !---------
  	 !     ! step two
  	 !     ! Define the discretization stencil

  	 !     ! Create an empty 3D, 11-pt stencil object
  	 !     call HYPRE_StructStencilCreate(3, 11, stencil_t, ierr_t)
  	 !     !call HYPRE_StructStencilCreate(3, 7, stencil_t, ierr_t)

  	 !     ! Define the geometry of the stencils. Each represents a relative
  	 !     ! offset (in the index space).

  	 !     offsets_t(1, :) = (/0, 0, 0/)
  	 !     offsets_t(2, :) = (/- 1, 0, 0/)
  	 !     offsets_t(3, :) = (/1, 0, 0/)
  	 !     offsets_t(4, :) = (/0, - 1, 0/)
  	 !     offsets_t(5, :) = (/0, 1, 0/)
  	 !     offsets_t(6, :) = (/0, 0, - 1/)
  	 !     offsets_t(7, :) = (/0, 0, 1/)

  	 !     offsets_t(8, :) = (/- 1, - 1, 0/)
  	 !     offsets_t(9, :) = (/- 1, 1, 0/)
  	 !     offsets_t(10, :) = (/1, - 1, 0/)
  	 !     offsets_t(11, :) = (/1, 1, 0/)

  	 !     do ihp = 0, 10
  	 !       !do ihp=0,6
  	 !       call HYPRE_StructStencilSetElement(stencil_t, ihp, offsets_t((ihp + 1), &
  	 !           :), ierr_t)
  	 !     end do

  	 !     !---------
  	 !     ! step three
  	 !     ! Set up a Struct Matrix

  	 !     ! Create an empty matrix object

  	 !     call HYPRE_StructMatrixCreate(mpi_comm_world, grid_t, stencil_t, A_hp_t, &
  	 !         ierr_t)

  	 !     ! Indicate that the matrix coefficients are ready to be set

  	 !     call HYPRE_StructMatrixInitialize(A_hp_t, ierr_t)

  	 !     ! set matrix values
  	 !     ! determine RHS belonging to solution

  	 !     do k = 1, nz
  	 !       do j = 1, ny
  	 !         do i = 1, nx
  	 !           index_count_t = i
  	 !           index_count_t = index_count_t + ((j - 1) * nx)
  	 !           index_count_t = index_count_t + ((k - 1) * nx * ny)

  	 !           index_count_t = (index_count_t * ne_hypre_t) - ne_hypre_t + 1

  	 !           if (i > 1) then
  	 !             AL = sin(real ((i - 1) * j * k + 2))
  	 !           else
  	 !             AL = sin(real (nx * j * k + 2))
  	 !           end if
  	 !           !AL = 1.0
  	 !           values_t(index_count_t + 1) = AL

  	 !           AR = sin(real (i * j * k + 2))
  	 !           !AR = 1.0
  	 !           values_t(index_count_t + 2) = AR

  	 !           if (j > 1) then
  	 !             AB = sin(real (i * (j - 1) * k + 4))
  	 !           else
  	 !             AB = sin(real (i * ny * k + 4))
  	 !           end if
  	 !           !AB = 1.0
  	 !           values_t(index_count_t + 3) = AB

  	 !           AF = sin(real (i * j * k + 4))
  	 !           !AF = 1.0
  	 !           values_t(index_count_t + 4) = AF

  	 !           if (k > 1) then
  	 !             AD = sin(real (i * j * (k - 1) + 6))
  	 !             !AD= sin(real(i*j*k+ 5))
  	 !             !AD= 1.0
  	 !           else
  	 !             AD = 0.0
  	 !           end if
  	 !           values_t(index_count_t + 5) = AD

  	 !           if (k < nz) then
  	 !             AU = sin(real (i * j * k + 6))
  	 !             !AU = 1.0
  	 !           else
  	 !             AU = 0.0
  	 !           end if
  	 !           values_t(index_count_t + 6) = AU

  	 !           !values_t(index_count_t   ) = sin(real(i*j*k   ))
  	 !           AC = - AR - AL - AF - AB - AU - AD
  	 !           values_t(index_count_t) = AC

  	 !           ALB = 0.0
  	 !           ALF = 0.0
  	 !           ARB = 0.0
  	 !           ARF = 0.0
  	 !           values_t(index_count_t + 7) = ALB
  	 !           values_t(index_count_t + 8) = ALF
  	 !           values_t(index_count_t + 9) = ARB
  	 !           values_t(index_count_t + 10) = ARF

  	 !           !testb
  	 !           !print*,i,j,k,'AC  =',values_t(index_count_t   )
  	 !           !print*,i,j,k,'AL  =',values_t(index_count_t+ 1)
  	 !           !print*,i,j,k,'AR  =',values_t(index_count_t+ 2)
  	 !           !print*,i,j,k,'AB  =',values_t(index_count_t+ 3)
  	 !           !print*,i,j,k,'AF  =',values_t(index_count_t+ 4)
  	 !           !print*,i,j,k,'AD  =',values_t(index_count_t+ 5)
  	 !           !print*,i,j,k,'AU  =',values_t(index_count_t+ 6)

  	 !           !print*,i,j,k,'ALB =',values_t(index_count_t+ 7)
  	 !           !print*,i,j,k,'ALF =',values_t(index_count_t+ 8)
  	 !           !print*,i,j,k,'ARB =',values_t(index_count_t+ 9)
  	 !           !print*,i,j,k,'ARF =',values_t(index_count_t+10)
  	 !           !teste

  	 !           b(i, j, k) = AC * sol(i, j, k)

  	 !           if (i > 1) then
  	 !             b(i, j, k) = b(i, j, k) + AL * sol(i - 1, j, k)
  	 !           else
  	 !             b(i, j, k) = b(i, j, k) + AL * sol(nx, j, k)
  	 !           end if

  	 !           if (i < nx) then
  	 !             b(i, j, k) = b(i, j, k) + AR * sol(i + 1, j, k)
  	 !           else
  	 !             b(i, j, k) = b(i, j, k) + AR * sol(1, j, k)
  	 !           end if

  	 !           if (j > 1) then
  	 !             b(i, j, k) = b(i, j, k) + AB * sol(i, j - 1, k)
  	 !           else
  	 !             b(i, j, k) = b(i, j, k) + AB * sol(i, ny, k)
  	 !           end if

  	 !           if (j < ny) then
  	 !             b(i, j, k) = b(i, j, k) + AF * sol(i, j + 1, k)
  	 !           else
  	 !             b(i, j, k) = b(i, j, k) + AF * sol(i, 1, k)
  	 !           end if

  	 !           if (k > 1) then
  	 !             b(i, j, k) = b(i, j, k) + AD * sol(i, j, k - 1)
  	 !           end if

  	 !           if (k < nz) then
  	 !             b(i, j, k) = b(i, j, k) + AU * sol(i, j, k + 1)
  	 !           end if
  	 !         end do ! i_loop
  	 !       end do ! j_loop
  	 !     end do ! k_loop

  	 !     ! This is a collective call finalizing the matrix assembly.
  	 !     ! The matrix is now ``ready to be used'

  	 !     call HYPRE_StructMatrixSetBoxValues(A_hp_t, (/(icoord) * nx + 1, (jcoord) &
  	 !         * ny + 1, 1/), (/(icoord + 1) * nx, (jcoord + 1) * ny, nz/), &
  	 !         ne_hypre_t, stencil_indices_t, values_t, ierr_t)

  	 !     !call HYPRE_StructMatrixSetBoxValues(A_hp_t, &
  	 !     !     & (/ (icoord-1)*nx, (jcoord-1)*ny , 0 /), &
  	 !     !     & (/ icoord*nx-1 , jcoord*ny-1 , nz-1 /), &
  	 !     !     & ne_hypre_t, stencil_indices_t, values_t, ierr_t)

  	 !     if (master) print *, 'HYPRE_StructMatrixSetBoxValues done'

  	 !     call HYPRE_StructMatrixAssemble(A_hp_t, ierr_t)

  	 !     if (master) print *, 'HYPRE_StructMatrixAssemble done'

  	 !     !---------
  	 !     ! step four
  	 !     ! Set up Struct Vectors for b and x.  Each processor sets the vectors
  	 !     ! corresponding to its boxes.

  	 !     ! Create an empty vector object

  	 !     call HYPRE_StructVectorCreate(mpi_comm_world, grid_t, b_hp_t, ierr_t)
  	 !     call HYPRE_StructVectorCreate(mpi_comm_world, grid_t, x_hp_t, ierr_t)

  	 !     ! Indicate that the vector coefficients are ready to be set

  	 !     call HYPRE_StructVectorInitialize(b_hp_t, ierr_t)
  	 !     call HYPRE_StructVectorInitialize(x_hp_t, ierr_t)

  	 !     ! Set the vector coefficients

  	 !     bmean = 0.0

  	 !     do k = 1, nz
  	 !       do j = 1, ny
  	 !         do i = 1, nx
  	 !           bmean = bmean + b(i, j, k)
  	 !         end do
  	 !       end do
  	 !     end do

  	 !     bmean = bmean / real (nx * ny * nz)

  	 !     print *, 'bmean =', bmean

  	 !     index_count_t = 1
  	 !     do k = 1, nz
  	 !       do j = 1, ny
  	 !         do i = 1, nx
  	 !           b(i, j, k) = b(i, j, k) - bmean

  	 !           bvalue_vector_t(index_count_t) = b(i, j, k)
  	 !           index_count_t = index_count_t + 1

  	 !           print *, i, j, k, 'b =', b(i, j, k)
  	 !         end do
  	 !       end do
  	 !     end do

  	 !     index_count_t = 1
  	 !     do k = 1, nz
  	 !       do j = 1, ny
  	 !         do i = 1, nx
  	 !           sol(i, j, k) = 0.0

  	 !           xvalue_vector_t(index_count_t) = sol(i, j, k)
  	 !           index_count_t = index_count_t + 1
  	 !         end do
  	 !       end do
  	 !     end do

  	 !     call HYPRE_StructVectorSetBoxValues(b_hp_t, (/(icoord) * nx + 1, (jcoord) &
  	 !         * ny + 1, 1/), (/(icoord + 1) * nx, (jcoord + 1) * ny, nz/), &
  	 !         bvalue_vector_t, ierr_t)

  	 !     !call HYPRE_StructVectorSetBoxValues(b_hp_t, &
  	 !     !   & (/ (icoord-1)*nx, (jcoord-1)*ny , 0 /), &
  	 !     !   & (/ icoord*nx-1 , jcoord*ny-1 , nz-1 /), &
  	 !     !   & bvalue_vector_t, ierr_t)

  	 !     if (master) print *, 'HYPRE_StructVectorSetBoxValues done'

  	 !     call HYPRE_StructVectorSetBoxValues(x_hp_t, (/(icoord) * nx + 1, (jcoord) &
  	 !         * ny + 1, 1/), (/(icoord + 1) * nx, (jcoord + 1) * ny, nz/), &
  	 !         xvalue_vector_t, ierr_t)

  	 !     !call HYPRE_StructVectorSetBoxValues(x_hp_t, &
  	 !     !   & (/ (icoord-1)*nx, (jcoord-1)*ny , 0 /), &
  	 !     !   & (/ icoord*nx-1 , jcoord*ny-1 , nz-1 /), &
  	 !     !   & xvalue_vector_t, ierr_t)

  	 !     if (master) print *, 'HYPRE_StructVectorSetBoxValues done'

  	 !     ! This is a collective call finalizing the vector assembly.
  	 !     ! The vectors are now ``ready to be used''

  	 !     call HYPRE_StructVectorAssemble(b_hp_t, ierr_t)

  	 !     if (master) print *, 'HYPRE_StructVectorAssemble done'

  	 !     call HYPRE_StructVectorAssemble(x_hp_t, ierr_t)

  	 !     if (master) print *, 'HYPRE_StructVectorAssemble done'

  	 !     !---------
  	 !     ! step five
  	 !     ! Set up and use a solver (See the Reference Manual for descriptions
  	 !     ! of all of the options.)

  	 !     ! Create an empty Hybrid Struct solver

  	 !     call HYPRE_StructHybridCreate(mpi_comm_world, solver_hp_t, ierr_t)

  	 !     ! Set print level (set 0 for complete silence or 2 for info on
  	 !     ! iterations)

  	 !     call HYPRE_StructHybridSetPrintLevel(solver_hp_t, 0, ierr_t)

  	 !     ! Set the type of Krylov solver to use.
  	 !     ! Current krylov methods set by solver type are:
  	 !     ! 1 DSCG (default)
  	 !     ! 2 GMRES
  	 !     ! 3 BiCGSTAB

  	 !     call HYPRE_StructHybridSetSolverType(solver_hp_t, 2, ierr_t)

  	 !     ! Set convergence criterion for activating a preconditioner.
  	 !     ! A number between 0 and 1 is accepted, the smaller the number, the
  	 !     ! sooner a preconditioner is activated.

  	 !     call HYPRE_StructHybridSetConvergenc(solver_hp_t, 1.e-23, ierr_t)

  	 !     ! Set maximum number of iterations.

  	 !     call HYPRE_StructHybridSetPCGMaxIter(solver_hp_t, 500, ierr_t)

  	 !     ! Setup

  	 !     call HYPRE_StructHybridSetup(solver_hp_t, A_hp_t, b_hp_t, x_hp_t, ierr_t)

  	 !     ! Finalize set up and use a solver
  	 !     ! (See the Reference Manual for descriptions of all of the options.)

  	 !     tol_t = 1.e-10
  	 !     call HYPRE_StructHybridSetTol(solver_hp_t, tol_t, ierr_t)

  	 !     if (master) print *, 'HYPRE_StructHybridSetTol done'

  	 !     ! show matrix elements

  	 !     call HYPRE_StructMatrixPrint(A_hp_t, 0, ierr_t)

  	 !     ! use solver

  	 !     call HYPRE_StructHybridSolve(solver_hp_t, A_hp_t, b_hp_t, x_hp_t, ierr_t)

  	 !     if (master) print *, 'HYPRE_StructHybridSolve done'

  	 !     ! Get the results

  	 !     call HYPRE_StructVectorGetBoxValues(x_hp_t, (/(icoord) * nx + 1, (jcoord) &
  	 !         * ny + 1, 1/), (/(icoord + 1) * nx, (jcoord + 1) * ny, nz/), &
  	 !         xvalue_vector_t, ierr_t)

  	 !     !call HYPRE_StructVectorGetBoxValues(x_hp_t, &
  	 !     !   & (/ (icoord-1)*nx, (jcoord-1)*ny , 0 /), &
  	 !     !   & (/ icoord*nx-1 , jcoord*ny-1 , nz-1 /), &
  	 !     !   & xvalue_vector_t, ierr_t)

  	 !     if (master) print *, 'HYPRE_StructVectorGetBoxValues done'

  	 !     index_count_t = 1

  	 !     solmean = 0.0

  	 !     do k = 1, nz
  	 !       do j = 1, ny
  	 !         do i = 1, nx
  	 !           sol(i, j, k) = xvalue_vector_t(index_count_t)
  	 !           index_count_t = index_count_t + 1

  	 !           solmean = solmean + sol(i, j, k)
  	 !         end do
  	 !       end do
  	 !     end do

  	 !     solmean = solmean / (nx * ny * nz)
  	 !     sol(:, :, :) = sol(:, :, :) - solmean

  	 !     do k = 1, nz
  	 !       do j = 1, ny
  	 !         do i = 1, nx
  	 !           print *, i, j, k, 'sol =', sol(i, j, k)
  	 !         end do
  	 !       end do
  	 !     end do

  	 !     call HYPRE_StructHybridGetNumIterati(solver_hp_t, nIter, ierr_t)

  	 !     if (master) print *, 'HYPRE_StructHybridGetNumIterati done'
  	 !     print *, 'nIter =', nIter

  	 !     index_count_t = 1
  	 !     call HYPRE_StructHybridGetFinalRelat(solver_hp_t, res, ierr_t)

  	 !     if (master) print *, 'HYPRE_StructHybridGetFinalRelat done'
  	 !     print *, 'res =', res

  	 !     return

  	 !   end subroutine test_hypre

end module bicgstab_tools_module
