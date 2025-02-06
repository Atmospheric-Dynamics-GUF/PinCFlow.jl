module bicgstab_tools_module

  use type_module
  use mpi_module
  use timeScheme_module
  use atmosphere_module
  use mpi

  implicit none

  private

  !------------------------
  !   public subroutines
  !------------------------

  public :: SetUpBiCGStab
  public :: CleanUpBiCGStab

  !------------------------
  !   private subroutines
  !------------------------

  contains

  subroutine SetUpBICGStab
    ! --------------------------------------
    !    generate arrays for BiCGStab
    !---------------------------------------

    integer :: allocstat

    allocate(ac_b(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:alloc failed"

    allocate(acv_b(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:alloc failed"

    allocate(ach_b(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:alloc failed"

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

    deallocate(acv_b, stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:dealloc failed"

    deallocate(ach_b, stat = allocstat)
    if(allocstat /= 0) stop "BiCGStab:dealloc failed"

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

    return

  end subroutine CleanUpBICGStab

end module bicgstab_tools_module
