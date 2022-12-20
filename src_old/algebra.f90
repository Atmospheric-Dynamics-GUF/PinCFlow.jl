module algebra_module


  implicit none

  public :: bicgstab_full          ! ref. Meister, pp. 172
  public :: bicgstab_csr      ! BiCGSTab using CSR storage format
  public :: precond_bicgstab  ! ref. Meister, p. 210
  public :: pre_bicgstab_csr  ! ilu(0)-preconditioned BiCGStab, CSR storage

  public :: lu            ! LU with L_ii = 1, ref. my notes
  public :: ilu           ! ILU(0) with no fill in, U_ii = 1, Meister pp.199
  public :: ilu0          ! ILU(0) with no fill in, L_ii = 1, ref. my notes
  public :: ilu_csr           ! ILU(0) for csr matrices A_csr, LU_csr

  public :: mldivide      ! forward / backward subst. for L,U in full format
  public :: mldivide_csr  ! forward / backward substs for L, U in CSR format


  ! basic algebra routines
  public :: norm3D        ! Frobenius norm / l2-norm for arrays
  public :: dot_product3D ! dot product for 3D arrays
  public :: norm          ! vector-2-norm
  public :: matrix_norm   ! Frobenius matrix norm
  public :: matmul_csr    ! matrix-vector product for CSR stored matrices
  public :: full2csr      ! converts A to CSR storage format
  public :: csr2full      ! converts CSR-stored A_csr to A in full format
  public :: printMatrix   ! print full matrix to screen in nice format

  interface mldivide
     module procedure mldivide
  end interface


  logical, parameter :: giveInfo = .true.


contains

  subroutine bicgstab_csr(A_csr,colInd,rowPtr,diagPtr,b,x,&
       & tol,res,max_iter_bicgstab)
    ! -----------------------------------------
    !  BiCGStab, cf. Meister pp. 172
    !  using compressed sparse matrix (CSR)
    ! -----------------------------------------

    real, dimension(:), intent(in) :: A_csr    ! matrix in CSR storage
    integer, dimension(:), intent(in) :: colInd, rowPtr,diagPtr
    real, dimension(:), intent(in) :: b          ! b=RHS,
    real, dimension(:), intent(inout) :: x       ! x = x0
    real, intent(in) :: tol                      ! target residual
    real, intent(out) :: res                     ! residual
    integer, intent(in) :: max_iter_bicgstab


    ! local variables
    integer :: j, n, allocstat
    real, dimension(:), allocatable :: p,r0,rOld,r,s,t,v
    real :: alpha, beta, omega
    real :: normB

    ! verbose
    logical, parameter :: giveInfo = .true.


    n = size(b)

    ! allocate local fields
    allocate(p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(r0(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(rOld(n), stat=allocstat);if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(r(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(s(n), stat=allocstat);  if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(t(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(v(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"

    ! Init
    r0 = b! - matmul_csr(A_csr,colInd,rowPtr,diagPtr,x,'low+up')
    p = r0
    r = r0
    normB = norm(b)


    ! check initial residual
    res = norm(r)
    if (giveInfo) write(*,fmt="(a25,es25.14)") " Initial residual: res0 = ", res
    if (giveInfo) write(*,fmt="(a25,es25.14)") " tol = ", tol
    if (res <= tol) then
       if(giveInfo) print*," ==> no iteration needed."
       return
    end if

    ! Loop
    iteration: do j = 1,max_iter_bicgstab
       v = matmul_csr(A_csr,colInd,rowPtr,diagPtr,p,'low+up')
       alpha = dot_product(r,r0) / dot_product(v,r0)
       s = r - alpha*v
       t = matmul_csr(A_csr,colInd,rowPtr,diagPtr, s,'low+up')

       omega = dot_product(t,s) / dot_product(t,t)
       x = x + alpha*p + omega*s

       rOld = r
       r = s - omega*t

       res = norm(r)
       if (res <= tol) then
          if (giveInfo)  then
             write(*,fmt="(a25,i25)") " Nb.of iterations: j = ", j
             write(*,fmt="(a25,es25.14)") " Final residual: res = ", res
             print*,"--------------------------------------------------"
             print*,""
          end if
          return
       end if

       beta = alpha/omega * dot_product(r,r0) / dot_product(rOld,r0)
       p = r + beta*(p-omega*v)
    end do iteration

    write(*,fmt="(a25,i25)") " BiCGStab_csr: max iterations!!!", &
         & max_iter_bicgstab
    print*,"--------------------------------------------------"; print*,""



    ! deallocate local fields
    deallocate(p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(r0, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(rOld,stat=allocstat);if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(r, stat=allocstat);if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(s, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(t, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(v, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"

  end subroutine bicgstab_csr


!--------------------------------------------------------------------------


  subroutine pre_bicgstab_csr(A_csr, LU_csr, colInd, rowPtr, diagPtr, &
       & b, x, tol, res, max_iter_bicgstab)
    !----------------------------------
    !  ILU-preconditioned BICGSTAB
    !  cf. PincFloit-Doku
    !  cf. Meister 1997
    !----------------------------------

    ! in / out
    real, dimension(:), intent(in) :: A_csr  ! A in CSR format
    real, dimension(:), intent(in) :: LU_csr ! incomplete L,U matrices, CSR
    ! CSR relevant info:
    integer, dimension(:), intent(in) :: colInd, rowPtr, diagPtr

    real, dimension(:), intent(in) :: b          ! b=RHS,
    real, dimension(:), intent(inout) :: x       ! x = x0 on input
    real, intent(in) :: tol                      ! target residual
    real, intent(out) :: res                     ! residual
    integer, intent(in) :: max_iter_bicgstab     ! max. nb. of iterations

    ! local variables
    integer :: j, n, allocstat
    real, dimension(:), allocatable :: p_p, r,r0_p,r_p, &
         & s,s_p, t,t_p, x_p, v,v_p, aux
    real :: alpha_p, beta_p, omega_p, rho_p, rho_p_old
    real :: normB

    logical, parameter :: giveInfo = .true.

    n = size(diagPtr)     ! or size(b), size(rowPtr)-1,

    ! allocate local fields
    allocate(p_p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(r(n), stat=allocstat); if(allocstat/=0) &
         &stop"algebra.f90/bicgstab:alloc failed"
    allocate(r0_p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(r_p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(s(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(s_p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(t(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(t_p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(x_p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(v(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(v_p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(aux(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"

    ! Init
    ! r -> r(0) = r0:
    r = b - matmul_csr(A_csr, colInd, rowPtr, diagPtr, x, 'low+up')
    r0_p = mldivide_csr(LU_csr, colInd, rowPtr, diagPtr, r, 'low')! P_l*r(0)
    p_p = r0_p                        ! p_p(0) = r0_p
    rho_p = dot_product(r0_p,r0_p)    ! rho_p(0)
    x_p = matmul_csr(LU_csr, colInd, rowPtr, diagPtr, x, 'up') ! x_p = U*x0
    normB = norm(b)                   ! norm of b for rel. residual

    ! prevent divsion by 0 = normB
    if (normB <= tol) then
       normB = 1.0
       if (giveInfo) print*,"algebra.f90/pre_bicgstab_csr: &
            & norm(b) < tol --> normB := 1.0"
    end if

    ! check initial residual
    if (norm(r)/normB <= tol) then
       res = norm(r)
       if(giveInfo) print*,"res = ", res
       if(giveInfo) print*,"algebra.f90/pre_bicgstab_csr: &
            & No iteration needed."
       return
    end if

    ! Loop
    iteration: do j = 1, max_iter_bicgstab
       ! mldivide(U,b,'up') = P_r * b:
       aux = mldivide_csr(LU_csr,colInd,rowPtr,diagPtr,p_p,'up')
       v = matmul_csr(A_csr,colInd,rowPtr,diagPtr, aux, 'low+up')
       ! mldivide(L,b,'low') = P_l * b:
       v_p = mldivide_csr(LU_csr,colInd,rowPtr,diagPtr,v,'low')

       alpha_p = rho_p / dot_product(v_p,r0_p)
       s = r - alpha_p * v
       s_p = mldivide_csr(LU_csr,colInd,rowPtr,diagPtr, s, 'low')

       aux =  mldivide_csr(LU_csr,colInd,rowPtr,diagPtr, s_p, 'up')
       t = matmul_csr(A_csr,colInd,rowPtr,diagPtr, aux, 'low+up')
       t_p = mldivide_csr(LU_csr,colInd,rowPtr,diagPtr, t, 'low')

       omega_p = dot_product(t_p,s_p) / dot_product(t_p,t_p)

       x_p = x_p + alpha_p * p_p + omega_p*s_p

       r = s - omega_p*t
       r_p = s_p - omega_p * t_p

       rho_p_old = rho_p
       rho_p = dot_product(r_p,r0_p)
       beta_p = alpha_p / omega_p * rho_p/rho_p_old

       p_p = r_p + beta_p * (p_p - omega_p * v_p)

       res = norm(r)
       if (res/normB <= tol) then
          if(giveInfo) then
             print*,"algebra.f90/precond_bicgstab: &
                  & Number of iterations needed: ", j
             print*,"relResPrecond = ", norm(r_p)/normB
             print*,"relRes = ", res/normB
             print*,""
          end if
          ! solve U*x = x_p for x:
          x = mldivide_csr(LU_csr,colInd,rowPtr,diagPtr, x_p, 'up')
          return
       end if
       !       if(giveInfo) print*,"res = ", res

    end do iteration


    ! max iterations
    print*,"algebra.f90/precond_bicgstab: max iterations reached: ", &
         & max_iter_bicgstab
    print*,"relResPrecond = ", norm(r_p)/normB
    print*,"relRes = ", res/normB
    print*,""
    x = mldivide_csr(LU_csr,colInd,rowPtr,diagPtr, x_p, 'up')               ! solve U*x = x_p for x


    ! deallocate local fields
    deallocate(p_p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(r, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(r0_p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(r_p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(s, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(s_p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(t, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(t_p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(x_p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(v, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(v_p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(aux, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"

  end subroutine pre_bicgstab_csr


!--------------------------------------------------------------------------


  subroutine ilu_csr(A_csr, colInd, rowPtr, diagPtr, LU_csr)
    !-------------------------------------------------
    ! calculates incomplete LU decomposition
    ! in compressed sparse row (CSR) storage format
    !-------------------------------------------------

    ! in/out
    real, dimension(:), intent(in) :: A_csr
    integer, dimension(:), intent(in) :: colInd
    integer, dimension(:), intent(in) :: rowPtr
    integer, dimension(:), intent(in) :: diagPtr
    real, dimension(:), allocatable, intent(out) :: LU_csr

    ! locals
    integer :: i,k,m,n
    integer :: nnz                          ! number of non-zeros
    integer :: allocstat
    integer :: csrInd, lInd, uInd
    real :: U_ik, L_im, U_mk, U_ii, A_ik
    real :: L_ki, L_km, U_mi, A_ki, sum


    ! get array size and allocate
    nnz = size(A_csr)
    n = size(diagPtr)

    allocate( LU_csr(nnz), stat=allocstat)
    if(allocstat /= 0) stop"algebra.f90/ilu_csr: alloc failed"
    LU_csr = 0.0

    ! ILU(0) for csr stored matrices

    iLoop: do i = 1, n

       ! -----------------------------------
       !         U_ik for k = i to n                     ! i-th row of U
       ! -----------------------------------

       kLoop1: do csrInd = diagPtr(i), rowPtr(i+1)-1    ! k = i to n, U_ik
          k = colInd(csrInd)                            ! column index of U

          ! determine sum             ! Sum( L_im * U_mk, m = 1 to i-1)
          sum = 0.0
          mLoop1: do lInd = rowPtr(i), diagPtr(i)-1
             m = colInd(lInd)
             L_im = LU_csr(lInd)

             ! search in row m for U_mk
             do uInd = diagPtr(m), rowPtr(m+1)-1
                if( colInd(uInd) == k) then
                   U_mk = LU_csr(uInd)
                   sum = sum + L_im*U_mk
                   exit
                end if
             end do
          end do mLoop1

          ! calc U_ik
          A_ik = A_csr(csrInd)
          U_ik = A_ik - sum
          LU_csr(csrInd) = U_ik

       end do kLoop1

       ! -------------------------------------
       !         L_ki for k = i+1 to n                ! i-th column of L
       ! -------------------------------------

       kLoop2: do k = i+1, n                 ! k = i+1 to i+q, q band width
          do csrInd = rowPtr(k), diagPtr(k)-1
             if( colInd(csrInd) == i) then   ! A_ki /= 0 -> calculate L_ki

                ! determine sum              ! Sum( L_km*U_mi, m = 1 to i-1)
                sum = 0.0
                mLoop2: do lInd = rowPtr(k), diagPtr(k)-1
                   m = colInd(lInd)
                   L_km = LU_csr(lInd)
                   ! search row m for U_mi:
                   do uInd = diagPtr(m), rowPtr(m+1)-1
                      if( colInd(uInd) == i) then
                         U_mi = LU_csr(uInd)
                         sum = sum + L_km * U_mi
                         exit
                      end if
                   end do
                end do mLoop2

                ! calc L_ki
                U_ii = LU_csr(diagPtr(i))
                A_ki = A_csr(csrInd)
                L_ki = (A_ki - sum) / U_ii
                LU_csr(csrInd) = L_ki

             end if
          end do
       end do kLoop2

    end do iLoop

  end subroutine ilu_csr


!---------------------------------------------------------------------------


  subroutine full2csr(A,A_csr,colInd,rowPtr,diagPtr)
    !--------------------------------
    !  transforms a full matrix to
    !  Matrix in CSR storage format
    !--------------------------------

    ! in/out
    real, dimension(:,:), intent(in) :: A
    real, dimension(:), allocatable, intent(out) :: A_csr
    integer, dimension(:), allocatable, intent(out) :: colInd
    integer, dimension(:), allocatable, intent(out) :: rowPtr
    integer, dimension(:), allocatable, intent(out) :: diagPtr

    ! local variables
    integer :: n        ! A = (n,n)
    integer :: i, j, csrInd
    integer :: nnz      ! number of non-zeros in A
    integer :: allocstat


    ! Find size of fields and allocate
    n = size(A,1)
    nnz = 0
    do j = 1,n
       do i = 1,n
          if ( A(i,j) /= 0) nnz = nnz + 1
       end do
    end do

    allocate( A_csr(nnz), stat=allocstat)
    if(allocstat /= 0) stop"algebra.f90/full2csr: alloc failed"
    allocate( colInd(nnz), stat=allocstat)
    if(allocstat /= 0) stop"algebra.f90/full2csr: alloc failed"
    allocate( rowPtr(n+1),   stat=allocstat)
    if(allocstat /= 0) stop"algebra.f90/full2csr: alloc failed"
    allocate( diagPtr(n),   stat=allocstat)
    if(allocstat /= 0) stop"algebra.f90/full2csr: alloc failed"

    ! write csr-format
    csrInd = 1
    do i = 1, n
       rowPtr(i) = csrInd
       do j = 1, n
          if( A(i,j) /= 0) then
             A_csr(csrInd) = A(i,j)
             colInd(csrInd) = j
             if (i==j) diagPtr(i) = csrInd
             csrInd = csrInd + 1
          end if
       end do
    end do

    rowPtr(n+1) = csrInd

  end subroutine full2csr


!---------------------------------------------------------------------------


  subroutine csr2full(A_csr, colInd, rowPtr, A)
    !------------------------------------------
    ! transforms a full matrix into a matrix
    ! in compact sparse row (CSR) format
    !------------------------------------------

    ! in/out
    real, dimension(:), intent(in) :: A_csr
    integer, dimension(:), intent(in) :: colInd
    integer, dimension(:), intent(in) :: rowPtr
    real, dimension(:,:), intent(out), allocatable :: A

    ! local variables
    integer :: n        ! A = (n,n)
    integer :: i, j, csrInd
    integer :: allocstat

    ! Find size of fields and allocate
    n = size(rowPtr)-1

    allocate( A(n,n), stat=allocstat)
    if(allocstat /= 0) stop"algebra.f90/csr2full: alloc failed"

    ! write zeros
    A = 0.0

    ! write non-zeros
    do i = 1, n
       do csrInd = rowPtr(i), rowPtr(i+1)-1
          j = colInd(csrInd)
          A(i,j) = A_csr(csrInd)
       end do
    end do

  end subroutine csr2full


!---------------------------------------------------------------------------


  function matmul_csr(M_csr, colInd, rowPtr, diagPtr, x, form)
    !------------------------------------------------------
    ! matrix vector multiplication for CSR stored matrices
    !------------------------------------------------------

    ! in / out
    real, dimension(:), allocatable :: matmul_csr
    real, dimension(:), intent(in) :: M_csr
    integer, dimension(:), intent(in) :: colInd, rowPtr, diagPtr
    real, dimension(:), intent(in) :: x
    character(len=*), intent(in) :: form   ! low, up or low+up


    ! local variables
    real, dimension(:), allocatable :: y
    integer :: i,k,n, allocstat

    ! allocate fields
    n = size(x)
    allocate( matmul_csr(n), stat=allocstat)
    if(allocstat/=0) stop"algebra.f90/matmul_csr: alloc failed"
    allocate( y(n), stat=allocstat)
    if(allocstat/=0) stop"algebra.f90/matmul_csr: alloc failed"


    select case(form)

    case('low+up')              ! M_csr is treated as full matrix
       y = 0.0
       do i = 1,n
          do k = rowPtr(i), rowPtr(i+1)-1
             y(i) = y(i) + M_csr(k) * x(colInd(k))
          end do
       end do


       case('low')   ! M_csr is treated as normed lower matrix L, L_ii = 1.0
          y = x      ! accounting for L_ii = 1.0
          do i = 1, n
             do k = rowPtr(i), diagPtr(i)-1
                y(i) = y(i) + M_csr(k) * x(colInd(k))
             end do
          end do

       case('up')     ! M_csr is treated as upper matrix U, U_ii <> 1.0
          y = 0.0
          do i = 1, n
             do k = diagPtr(i), rowPtr(i+1)-1
                y(i) = y(i) + M_csr(k) * x(colInd(k))
             end do
          end do

       case default
          stop"algebra.f90/matmul_csr: &
               & form should be 'up', 'low', or 'low+up'. Stop"

       end select
       matmul_csr = y

  end function matmul_csr


!---------------------------------------------------------------------------


  function mldivide_csr(LU_csr, colInd, rowPtr, diagPtr, b, form)
    ! in / out
    real, dimension(:), allocatable :: mldivide_csr
    real, dimension(:), intent(in) :: LU_csr
    integer, dimension(:), intent(in) :: colInd, rowPtr, diagPtr
    real, dimension(:), intent(in) :: b
    character(len=*), intent(in) :: form

    ! local variables
    real, dimension(:), allocatable :: x
    integer :: i, csrInd, m, n
    integer :: allocstat
    real :: sum, L_im, x_m
    real :: U_nn, U_ii, U_im


    ! get system size and allocate
    n = size(diagPtr)
    allocate( x(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/mldivide_csr:alloc failed"
    allocate( mldivide_csr(n),stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/mldivide_csr:alloc failed"

    select case(form)

    case('low')   ! use lower part of LU_csr = L, assuming L_ii = 1.0

       x(1) = b(1)
       do i = 2, n         ! loop over x_i
          sum = 0.0
          do csrInd = rowPtr(i), diagPtr(i)-1
             m = colInd(csrInd)
             L_im = LU_csr(csrInd)
             sum = sum + L_im*x(m)
          end do
          x(i) = b(i) - sum
       end do
       mldivide_csr = x

       case('up') ! use upper part of LU_csr = U, U_ii <> 1.0

          U_nn = LU_csr(diagPtr(n))
          x(n) = b(n) / U_nn

          do i = n-1, 1, -1
             sum = 0.0
             do csrInd = diagPtr(i)+1, rowPtr(i+1)-1
                m = colInd(csrInd)
                U_im = LU_csr(csrInd)
                sum = sum + U_im*x(m)
             end do
             U_ii = LU_csr(diagPtr(i))
             x(i) = (b(i) - sum) / U_ii
          end do
          mldivide_csr = x

       case default
          stop"algebra.f90/mldivide_csr: use form = 'low' or 'up'. Stop."

       end select

     end function mldivide_csr


!---------------------------------------------------------------------------


  function mldivide(A,b,form)
    real, dimension(:), allocatable :: mldivide
    character(len=*), intent(in) :: form
    real, dimension(:,:), intent(in) :: A
    real, dimension(:), intent(in) :: b

    ! local variables
    real, dimension(:), allocatable :: x
    integer :: i, n, allocstat

    n = size(A,1)

    allocate( x(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/mldivide:alloc failed"
    allocate( mldivide(n),stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/mldivide:alloc failed"

    select case(form)

    case('up')   ! A = upper triangular matrix

       do i = n,1,-1
          x(i) = (b(i) - dot_product(A(i,i+1:n),x(i+1:n)) ) / A(i,i)
       end do

    case('low')    ! A = lower triangular matrix
       !
       do i = 1,n
          x(i) = ( b(i) - dot_product(A(i,1:i-1),x(1:i-1)) ) / A(i,i)
       end do

    case default
       stop"algebra.f90/mldivide: matrix form incorrect"

    end select

    mldivide = x

    deallocate(x,stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/mldivide:dealloc failed"

  end function mldivide


!---------------------------------------------------------------------------


  subroutine lu(A,L,U) ! LU with L_ii = 1, see docu
    !-----------------------------------
    !     standard LU decomposition
    !-----------------------------------

    ! in/out variables
    integer :: n       ! size of A,L and U
    real, dimension(:,:), intent(in) :: A
    real, dimension(:,:), intent(out) :: L,U

    ! local variables
    integer :: i,k,m
    real :: sum

    ! init
    L = 0.0
    U = 0.0

    n = size(A,1)

    outerLoop: do i = 1,n

       ! ----------------------------
       !       U_ik, k = i to n
       ! ----------------------------

       do k = i, n
          ! Sum( L_im * U_mk, m = 1 to i-1
          sum = 0.0
          do m = 1, i-1
             sum = sum + L(i,m) * U(m,k)
          end do
          U(i,k) = A(i,k) - sum
       end do

       ! ----------------------------
       !       L_ki, k = i+1 to n
       ! ----------------------------

       do k = i+1, n
          ! Sum( L_km * U_mi, m = 1 to i-1
          sum = 0.0
          do m = 1, i-1
             sum = sum + L(k,m) * U(m,i)
          end do
          L(k,i) = (A(k,i) - sum) / U(i,i)
       end do
       L(i,i) = 1.0

    end do outerLoop

  end subroutine lu


!---------------------------------------------------------------------------


  subroutine ilu0(A,L,U)
    !-----------------------------
    ! incomplete LU decomposition
    ! for full matrices
    !-----------------------------

    ! in/out variables
    integer :: n           ! size of A,L and U
    real, dimension(:,:), intent(in) :: A
    real, dimension(:,:), intent(out) :: L,U

    ! local variables
    integer :: i,k,m
    real :: sum

    ! init
    L = 0.0
    U = 0.0

    n = size(A,1)

    outerLoop: do i = 1,n

       ! ----------------------------
       !       U_ik, k = i to n
       ! ----------------------------

       do k = i, n
          if ( A(i,k) /= 0.) then

             ! Sum( L_im * U_mk, m = 1 to i-1
             sum = 0.0
             do m = 1, i-1
                if ( A(i,m) /= 0. .and. A(m,k) /=0. ) then
                   sum = sum + L(i,m) * U(m,k)
                end if
             end do
             U(i,k) = A(i,k) - sum
          end if
       end do

       ! ----------------------------
       !       L_ki, k = i+1 to n
       ! ----------------------------

       do k = i+1, n
          if ( A(k,i) /= 0.0 ) then

             ! Sum( L_km * U_mi, m = 1 to i-1
             sum = 0.0
             do m = 1, i-1
                if ( A(k,m) /=0. .and. A(m,i) /=0. ) then
                   sum = sum + L(k,m) * U(m,i)
                end if
             end do
             L(k,i) = (A(k,i) - sum) / U(i,i)
          end if

       end do
       L(i,i) = 1.0

    end do outerLoop

  end subroutine ilu0


!--------------------------------------------------------------------------


  subroutine ilu(A,L,U)
    !--------------------------------------
    !   LU decomposition cf. Meister 1997
    !--------------------------------------

    ! in/out variables
    integer :: n            ! size of A,L and U
    real, dimension(:,:), intent(in) :: A
    real, dimension(:,:), intent(out) :: L,U

    ! local variables
    integer :: i,k,m
    real :: sum

    ! init
    L = 0.0
    U = 0.0

    n = size(A,1)

    outerLoop: do i = 1,n

       do k = i,n
          if ( A(k,i) /= 0.) then
             sum = 0.0
             do m = 1,i-1
                if ( A(k,m) /= 0. .and. A(m,i) /=0. ) then
                   sum = sum + L(k,m)*U(m,i)
                end if
             end do
             L(k,i) = A(k,i) - sum
          end if
       end do

       do k = i,n
          if ( A(i,k) /=0. ) then
             sum = 0.0
             do m = 1,i-1
                if ( A(i,m) /=0. .and. A(m,k) /=0. ) then
                   sum = sum + L(i,m)*U(m,k)
                end if
             end do
             U(i,k) = 1./L(i,i) * ( A(i,k) - sum)
          end if
       end do

    end do outerLoop

  end subroutine ilu


!---------------------------------------------------------------------------


  subroutine precond_bicgstab(A,L,U,b,x,tol,res,max_iter_bicgstab)
    !-----------------------------------------
    ! precond. bicgstab, Meister 1997, p. 210
    !-----------------------------------------

    ! in/out variables
    real, dimension(:,:), intent(in) :: A        ! A
    real, dimension(:,:), intent(in) :: L,U      ! incomplete LU matrices
    real, dimension(:), intent(in) :: b          ! b=RHS,
    real, dimension(:), intent(inout) :: x       ! x = x0
    real, intent(in) :: tol                      ! target residual
    real, intent(out) :: res                     ! residual
    integer, intent(in) :: max_iter_bicgstab     ! max. nb. of iterations

    ! local variables
    integer :: j, n, allocstat
    real, dimension(:), allocatable :: p_p, r,r0_p,r_p, s,s_p, &
         & t,t_p, x_p, v,v_p
    real :: alpha_p, beta_p, omega_p, rho_p, rho_p_old
    real :: normB

    n = size(A,1)

    ! allocate local fields
    allocate(p_p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(r(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(r0_p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(r_p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(s(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(s_p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(t(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(t_p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(x_p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(v(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(v_p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"

    ! Init
    r = b - matmul(A,x)               ! r -> r(0) = r0
    r0_p = mldivide(L,r,'low')        ! P_l*r(0)
    p_p = r0_p                        ! p_p(0) = r0_p
    rho_p = dot_product(r0_p,r0_p)    ! rho_p(0)
    x_p = matmul(U,x)                 ! x_p = U*x0
    normB = norm(b)                   ! norm of b for rel. residual

    ! check initial residual
    if (norm(r)/normB <= tol) then
       res = norm(r)
       if(giveInfo) print*,"res = ", res
       if(giveInfo) print*,"algebra.f90/bicgstab: No iteration needed."
       return
    end if

    ! Loop
    iteration: do j = 1, max_iter_bicgstab
       v = matmul(A, mldivide(U,p_p,'up') )  ! mldivide(U,b,'up') = P_r * b
       v_p = mldivide(L,v,'low')             ! mldivide(L,b,'low') = P_l * b

       alpha_p = rho_p / dot_product(v_p,r0_p)
       s = r - alpha_p * v
       s_p = mldivide(L,s,'low')

       t = matmul(A, mldivide(U,s_p,'up'))
       t_p = mldivide(L,t,'low')

       omega_p = dot_product(t_p,s_p) / dot_product(t_p,t_p)

       x_p = x_p + alpha_p * p_p + omega_p*s_p

       r = s - omega_p*t
       r_p = s_p - omega_p * t_p

       rho_p_old = rho_p
       rho_p = dot_product(r_p,r0_p)
       beta_p = alpha_p / omega_p * rho_p/rho_p_old

       p_p = r_p + beta_p * (p_p - omega_p * v_p)

       res = norm(r)
       if (res/normB <= tol) then
          if(giveInfo) then
             print*,"algebra.f90/precond_bicgstab: &
                  & Number of iterations needed: ", j
             print*,"relResPrecond = ", norm(r_p)/normB
             print*,"relRes = ", res/normB
             print*,""
          end if
          x = mldivide(U,x_p,'up')
          return
       end if
       !       if(giveInfo) print*,"res = ", res


    end do iteration

    ! max iterations
    print*,"algebra.f90/precond_bicgstab: max iterations reached: ", &
         & max_iter_bicgstab
    x = mldivide(U,x_p,'up')

    ! deallocate local fields
    deallocate(p_p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(r, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(r0_p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(r_p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(s, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(s_p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(t, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(t_p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(x_p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(v, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(v_p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"

  end subroutine precond_bicgstab


!---------------------------------------------------------------------------


  subroutine bicgstab_full(A,b,x,tol,res,max_iter_bicgstab)
    !-----------------------------------------
    ! BICGSTAB, Meister 1997, p. 210
    !-----------------------------------------

    ! in/out variables
    real, dimension(:,:), intent(in) :: A
    real, dimension(:), intent(in) :: b          ! b=RHS,
    real, dimension(:), intent(inout) :: x       ! x = x0
    real, intent(in) :: tol                      ! target residual
    real, intent(out) :: res                     ! residual
    integer, intent(in) :: max_iter_bicgstab

    ! local variables
    integer :: j, n, allocstat
    real, dimension(:), allocatable :: p,r0,rOld,r,s,t,v
    real :: alpha, beta, omega
    real :: normB

    n = size(A,1)

    ! allocate local fields
    allocate(p(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(r0(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(rOld(n), stat=allocstat);if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(r(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(s(n), stat=allocstat);  if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(t(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"
    allocate(v(n), stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:alloc failed"

    ! Init
    r0 = b - matmul(A,x)
    p = r0
    r = r0
    normB = norm(b)

    ! check initial residual
    res = norm(r)
    if (res/normB <= tol) then
       if(giveInfo) print*,"relRes = ", res/normB
       if(giveInfo) print*,"algebra.f90/bicgstab: No iteration needed."
       return
    end if

    ! Loop
    iteration: do j = 1,max_iter_bicgstab
       v = matmul(A,p)
       alpha = dot_product(r,r0) / dot_product(v,r0)
       s = r - alpha*v
       t = matmul(A,s)

       omega = dot_product(t,s) / dot_product(t,t)
       x = x + alpha*p + omega*s

       rOld = r
       r = s - omega*t

       res = norm(r)
       if (res/normB <= tol) then
          if(giveInfo) print*,"algebra.f90/bicgstab: &
               & Number of iterations needed: ", j
          if(giveInfo) print*,"relResidual = ", res/normB
          return
       end if

       beta = alpha/omega * dot_product(r,r0) / dot_product(rOld,r0)
       p = r + beta*(p-omega*v)
    end do iteration

    print*,"algebra.f90/bicgstab: max iterations reached", max_iter_bicgstab

    ! deallocate local fields
    deallocate(p, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(r0, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(rOld,stat=allocstat);if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(r, stat=allocstat);if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(s, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(t, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"
    deallocate(v, stat=allocstat); if(allocstat/=0) &
         & stop"algebra.f90/bicgstab:dealloc failed"

  end subroutine bicgstab_full


!--------------------------------------------------------------------------


  function dot_product3D(a,b)
    !---------------------------
    ! dot product for 3D arrays
    !---------------------------

    ! in/out variables
    real :: dot_product3D
    real, dimension(:,:,:), intent(in) :: a,b

    ! locals
    integer, dimension(3) :: aSize, bSize
    integer :: i,j,k

    aSize = shape(a)
    bSize = shape(b)

    do i = 1,3
       if( aSize(i) .ne. bSize(i) ) stop"dot_product3D failure."
    end do

    dot_product3D = 0.0
    do k = 1, aSize(3)
       do j = 1, aSize(2)
          dot_product3D = dot_product3D + dot_product( a(:,j,k),b(:,j,k) )
       end do
    end do

  end function dot_product3D


!--------------------------------------------------------------------------


  function norm3D(a)
    !------------------------------
    ! Frobinius-norm for 3D arrays
    !------------------------------

    ! in/out variables
    real :: norm3D
    real, dimension(:,:,:), intent(in) :: a

    ! locals
    integer, dimension(3) :: aSize
    integer :: i,j,k

    aSize = shape(a)

    norm3D = 0.0
    do k = 1, aSize(3)
       do j = 1, aSize(2)
          norm3D = norm3D + dot_product( a(:,j,k),a(:,j,k) )
       end do
    end do
    norm3D = sqrt(norm3D)

  end function norm3D


!--------------------------------------------------------------------------


  function norm(a)
    ! 2-Norm for vectors
    real :: norm
    real, dimension(:), intent(in) :: a

    norm = sqrt( dot_product(a,a) )

  end function norm


!--------------------------------------------------------------------------


  function matrix_norm(A)
    ! Frobenius norm for matrices
    real :: matrix_norm
    real, dimension(:,:), intent(in) :: A
    integer :: i,n

    matrix_norm = 0.
    n = size(A,1)

    do i = 1,n
       matrix_norm = matrix_norm + dot_product( A(i,:), A(i,:) )
    end do
    matrix_norm = sqrt(matrix_norm)

  end function matrix_norm


!--------------------------------------------------------------------------


  subroutine printMatrix(A,name,inputForm)
    !----------------------------------------------
    ! print full matrix in nice form: line by line
    !----------------------------------------------

    ! in / out variables
    real, dimension(:,:), intent(in) :: A
    character(len=*), intent(in),optional :: name
    character(len=*), intent(in), optional :: inputForm

    ! locals
    real, dimension(:), allocatable :: Aline
    integer :: m, n, i
    integer :: allocstat
    character(len=20) :: form

    ! get size and allocate


    n = size(A,1)
    m = size(A,2)
    if ( n /= m ) then
       print*,"n <> m. returning"
       return
    end if

    allocate(Aline(n), stat=allocstat)
    if (allocstat /= 0) stop"testAlgebra.f90/printMatrix: alloc failed."

    ! define output format using n...es25.14
    if ( present(inputForm) ) then
       form = inputForm
    else
       write(unit=form,fmt="(a,i3,a)") "(", n, "es25.14)"
    end if

    ! print matrix
    if ( present(name) ) then
       write(*,fmt="(2a)") name, " = "
    else
       print*,"Matrix = "
    end if

    do i = 1, n
       Aline = A(i,:)
       write(*,fmt="(4es25.14)") Aline
    end do

    deallocate(Aline, stat=allocstat)
    if (allocstat /= 0) stop"testAlgebra.f90/printMatrix: dealloc failed."

  end subroutine printMatrix


end module algebra_module
