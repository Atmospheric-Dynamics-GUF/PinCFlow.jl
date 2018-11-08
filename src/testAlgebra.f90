program test_program
! programme for testing algbra_module subroutines

use algebra_module

  implicit none


  ! common variables and parameter
  integer :: max_iter_bicgstab
  integer :: i,j

  ! bicgstab test
  logical :: bicgstab_test = .false.
  real, dimension(3,3) :: A
  real, dimension(3) :: b, x, xSol, line
  real :: tol, res
  
  ! operator_bicgstab
  logical :: opbicgstab_test = .true.

  ! ilu  test0
  logical :: ilu_test = .false.
  real, dimension(3,3) :: L1, U1
  real, dimension(3) :: v

  ! pre_bicgstab test1
  integer,parameter :: n = 15
  integer :: ndiag = 6
  logical :: precond_test = .false.
  real, dimension(n,n) :: H, L, U
  real, dimension(n) :: c, y, ySol
  
  ! ilu_csr test2
  logical, parameter :: ilucsr_test = .false.
  real, dimension(4,4) :: A2, L2, U2
  real, dimension(:,:), allocatable :: LU2
  real, dimension(:), allocatable :: A2_csr
  integer, dimension(:), allocatable :: colInd2
  integer, dimension(:), allocatable :: rowPtr2
  integer, dimension(:), allocatable :: diagPtr2

  real, dimension(:), allocatable :: LU2_csr

  ! mldivide_csr test2
  logical, parameter :: mldivide_csr_test = .false.
  real, dimension(4) :: x2, y2, b2

  ! matmul_csr test
  logical, parameter :: matmul_csr_test = .false.
  real, dimension(4) :: Lx2, Ux2, Ax2

  ! pre_bicgstab_csr test4
  logical, parameter :: pre_bicgstab_csr_test = .false.
  integer, parameter :: n4 = 100
  real, dimension(n4,n4) :: A4
  real, dimension(:), allocatable :: A4_csr, LU4_csr
  integer, dimension(:), allocatable :: colInd4, rowPtr4, diagPtr4
  real, dimension(n4) :: x4, b4, xSol4
  integer :: off1, off2, maxIterPoisson4
  real :: tol4, res4

  ! test reshape
  logical, parameter :: reshapeTest = .false.
  real, dimension(9) :: line1
  real, dimension(3,3) :: field1








  !--------------------------
  !     operator_bicgstab test
  !--------------------------


  if (opbicgstab_test) then
     print*,""; print*,"       ------------- operator bicgstab test -----------------"; print*,""
     A(1,:) = (/ 5., 2., 3. /)
     A(2,:) = (/ 4., 15., 6. /)
     A(3,:) = (/ 7., 8., 29. /)
     xSol = (/ 1.0, 2.0, 3.0 /)

     b = matmul(A,xSol)
     

     print*,"A = "
     do i = 1,3
        line = A(i,:)
        print*,line
     end do
     print*,"b = "
     print*,b  
  
     x = 0.0        ! start vector
     tol = 1.0e-10                  ! final residual
     max_iter_bicgstab = 150

     call operator_bicgstab(A,b,x,tol,res,max_iter_bicgstab)
     print*, "x = "
     print*, x
     print*, "|x-xSol| = "
     print*, norm(xSol-x)
     print*, "res = "
     print*, res
  end if








  ! ---------------------------------
  !             reshape test
  ! --------------------------------
  
  if (reshapeTest) then
     
     do i = 1,9; line1(i) = real(i)
     end do
     print*,"line1 = ", line1

     field1 = reshape(line1,(/ 3,3 /) )
     do i = 1,3
        do j = 1,3
           print*,"i,j,field1(i,j)", i, j, field1(i,j)
        end do
     end do

  end if

  !----------------------------------
  !   pre_bicgstab_csr test4
  !----------------------------------
  ! successful test: comparison with matlab code precond_bicgstab.m

  if (pre_bicgstab_csr_test) then
    print*,""; print*,"       ------------- pre_bicgstab_csr test -----------------"; print*,""

    ! setting up banded matrix
    off1 = 1
    off2 = 3 ! offset for super and sub diagonals

    A4 = 0.0
    do i = 1, n4
      A4(i,i) = -4.0            ! diagonal
      j = i + off1              ! first super diagonal
      if (j<=n4) A4(i,j) = 1.0
      j = i + off2              ! second super diagonal
      if (j<=n4) A4(i,j) = 1.0
      j = i-off1                ! first subdiagonal
      if (j>=1) A4(i,j) = 1.0
      j = i-off2                ! second subdiagonal
      if (j>=1) A4(i,j) = 1.0
   end do
   
   ! set up problem
   xSol4 = 0.0              ! exact solution
   b4 = matmul(A4,xSol4)    ! according right hand side
   x4 = 10.0                 ! start solution
      
   ! CSR stored matrix
   call full2csr(A4,A4_csr, colInd4, rowPtr4, diagPtr4)

   ! envoke bicgstab
   tol4 = 1.0e-5
   maxIterPoisson4 = 150

   call ilu_csr(A4_csr, colInd4, rowPtr4, diagPtr4, LU4_csr)
   call pre_bicgstab_csr(A4_csr, LU4_csr, colInd4, rowPtr4, diagPtr4, b4, x4, tol4, res4, maxIterPoisson4)
   
   ! assess results
   print*,"|x-xSol|  = ", norm(x4-xSol4)
   print*,"res = ", res4
   
   
   



  end if


  !-------------------------------
  !       matmul_csr test2
  !-------------------------------

 if (matmul_csr_test) then
    print*,""; print*,"       ------------- matmul_csr test -----------------"; print*,""
     A2(1,:) = (/ -4., 2., 0., 0. /)
     A2(2,:) = (/ 2., -4., 2., 0. /)
     A2(3,:) = (/ 0., 2., -4., 2. /)
     A2(4,:) = (/ 0., 0., 2., -4. /)

     x2 = 1.0
     
     call full2csr(A2,A2_csr, colInd2, rowPtr2, diagPtr2)
     call ilu_csr(A2_csr, colInd2, rowPtr2, diagPtr2, LU2_csr)
     
     ! results for CSR stored matrices
     Lx2 = matmul_csr(LU2_csr,colInd2,rowPtr2,diagPtr2, x2, 'low')
     Ux2 = matmul_csr(LU2_csr,colInd2,rowPtr2,diagPtr2, x2, 'up')
     Ax2 = matmul_csr(A2_csr,colInd2,rowPtr2,diagPtr2, x2, 'low+up')
     print*,""
     print*,"matmul_csr / L*x2 = "; print*, Lx2
     print*,"matmul_csr / U*x2 = "; print*, Ux2
     print*,"matmul_csr / A*x2 = "; print*, Ax2
     print*,"------------------------------------"; print*,""


     ! compare with matmul for full matrices
     call ilu0(A2,L2,U2)
     Lx2 = matmul(L2,x2)
     Ux2 = matmul(U2,x2)
     Ax2 = matmul(A2,x2)
     print*,"matmul / L*x2 = "; print*, Lx2
     print*,"matmul / U*x2 = "; print*, Ux2
     print*,"matmul / A*x2 = "; print*, Ax2
     print*,"------------------------------------"; print*,""
  
  end if



  !--------------------------------
  !       mldivide_csr test2
  !--------------------------------
  
  if (mldivide_csr_test) then
    print*,""; print*,"       ------------- mldivide_csr test -----------------"; print*,""
     A2(1,:) = (/ -4., 2., 0., 0. /)
     A2(2,:) = (/ 2., -4., 2., 0. /)
     A2(3,:) = (/ 0., 2., -4., 2. /)
     A2(4,:) = (/ 0., 0., 2., -4. /)

     b2 = 1.0
     
     call full2csr(A2,A2_csr, colInd2, rowPtr2, diagPtr2)
     call ilu_csr(A2_csr, colInd2, rowPtr2, diagPtr2, LU2_csr)
     
     ! results for CSR stored matrices
     y2 = mldivide_csr(LU2_csr, colInd2, rowPtr2, diagPtr2, b2, 'low')
     x2 = mldivide_csr(LU2_csr, colInd2, rowPtr2, diagPtr2, y2, 'up')
     print*,""
     print*,"mldivide_csr/y2 = "; print*, y2
     print*,"mldivide_csr/x2 = "; print*, x2
     print*,"------------------------------------"; print*,""


     ! compare with mldivide for full matrices
     call ilu0(A2,L2,U2)
     y2 = mldivide(L2,b2,'low')
     x2 = mldivide(U2,y2,'up')
     print*,"mldivide/y2 = "; print*, y2
     print*,"mldivide/x2 = "; print*, x2
     print*,"------------------------------------"; print*,""
  
  end if




  !--------------------------------
  !      ilu_csr test (succesful)
  !--------------------------------
  
  if (ilucsr_test) then
    print*,""; print*,"       ------------- ilucsr_csr test -----------------"; print*,""
     A2(1,:) = (/ -4., 2., 0., 0. /)
     A2(2,:) = (/ 2., -4., 2., 0. /)
     A2(3,:) = (/ 0., 2., -4., 2. /)
     A2(4,:) = (/ 0., 0., 2., -4. /)

     call full2csr(A2,A2_csr, colInd2, rowPtr2, diagPtr2)
     call ilu_csr(A2_csr, colInd2, rowPtr2, diagPtr2, LU2_csr)
     call csr2full(LU2_csr, colInd2, rowPtr2, LU2)
     print*,"ilu_csr: "
     call printMatrix(LU2, name = "LU2")
     print*,"-----------------------------------"

     print*,"my ilu: "
     call ilu0(A2,L2,U2)
     call printMatrix(A2, name="A2")
     call printMatrix(L2, name="L2")
     call printMatrix(U2, name="U2")
     print*, "|A2-L2*U2| = "
     print*, matrix_norm(A2-matmul(L2,U2))
     print*,"------------------------------------"
     
  
  end if


  !--------------------------
  !     bicgstab test
  !--------------------------


  if (bicgstab_test) then
    print*,""; print*,"       ------------- bicgstab test -----------------"; print*,""
     A(1,:) = (/ 5, 2, 3 /)
     A(2,:) = (/ 4, 15, 6 /)
     A(3,:) = (/ 7, 8, 29 /)
     xSol = (/ 0.3, 0.7, 0.9 /)
     b = matmul(A,xSol)
     

     print*,"A = "
     do i = 1,3
        line = A(i,:)
        print*,line
     end do
     print*,"b = "
     print*,b  
  
     x = 0.0        ! start vector
     tol = 1.0e-10                  ! final residual
     max_iter_bicgstab = 150

     call bicgstab(A,b,x,tol,res,max_iter_bicgstab)
     print*, "x = "
     print*, x
     print*, "|x-xSol| = "
     print*, norm(xSol-x)
     print*, "res = "
     print*, res
  end if
  
  

  !--------------------------
  !     ilu test
  !--------------------------
  
  if (ilu_test) then
     print*,"testing incomplete LU decomposition..."
     
     A(1,:) = (/ 5, 2, 3 /)
     A(2,:) = (/ 4, 15, 6 /)
     A(3,:) = (/ 7, 8, 29 /)
     xSol = (/ 0.3, 0.7, 0.9 /)
     b = matmul(A,xSol)
     L1 = 0.
     U1 = 0.

     call ilu(A,L,U)
     print*,"A = "
     do i = 1,3
        line = A(i,:)
        print*,line
     end do
     print*,"L = "
     do i = 1,3
        line = L1(i,:)
        print*,line
     end do

     print*,"U = "
     do i = 1,3
        line = U1(i,:)
        print*,line
     end do

     print*, "|A-L1*U1| = "
     print*, matrix_norm(A-matmul(L1,U1))
     
!     call mldivide(L1,b,v,'lower')
     v = mldivide(L1,b,'low')
!     call mldivide(U1,v,x,'upper')
     x = mldivide(U1,v,'up')
     
     
    
     print*, "|x-xSol| = "
     print*, norm(xSol-x)

     
  end if

  !----------------------------------
  !   preconditioned bicgstab test
  !----------------------------------

  if (precond_test) then
     print*,"testing preconditioned bicgstab..."

     if (1==2) then
        ! set up Hilbert matrix
        do i = 1,n
           do j = 1,n
              H(i,j) = 1. / real(i+j-1)
           end do
        end do
     else
        ! set up simple matrix
        do i = 1,n
           do j = 1,n
              H(i,j) = sqrt(real(i)* real(j)) / (real(i) + real(j) )
           end do
        end do
     end if
  


     ! create band matrix
     do i = 1,n
        do j = 1,n
           if (j > i+ndiag) H(i,j) = 0.0
           if (j < i-ndiag) H(i,j) = 0.0
        end do
     end do

     ySol = 1.0
     c = matmul(H,ySol)
     y = 0.0
     tol = 1.0e-5
     max_iter_bicgstab = 150

     call ilu(H,L,U)
     print*, "|H| = ", matrix_norm(H)
     print*, "|L| = ", matrix_norm(L)
     print*, "|U| = ", matrix_norm(U)
     print*, "|H-L*U| = ", matrix_norm(H - matmul(L,U))

     call precond_bicgstab(H,L,U,c,y,tol,res,max_iter_bicgstab)
!     print*, "y = "
!     print*, y
     print*, "|y-ySol| = "
     print*, norm(ySol-y)
     print*, "res = "
     print*, res

  end if



      
      
      
   


end program test_program
