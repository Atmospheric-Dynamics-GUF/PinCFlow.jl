  subroutine pre_bicgstab_csr(A_csr, LU_csr, colInd, rowPtr, diagPtr, &
       & b, x, tol, res, max_iter_bicgstab) 
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
         & stop "algebra.f90/bicgstab:alloc failed"
    allocate(r(n), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:alloc failed"
    allocate(r0_p(n), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:alloc failed"
    allocate(r_p(n), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:alloc failed"
    allocate(s(n), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:alloc failed"
    allocate(s_p(n), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:alloc failed"
    allocate(t(n), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:alloc failed"
    allocate(t_p(n), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:alloc failed"
    allocate(x_p(n), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:alloc failed"
    allocate(v(n), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:alloc failed"    
    allocate(v_p(n), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:alloc failed"
    allocate(aux(n), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:alloc failed"

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
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(r, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(r0_p, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(r_p, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(s, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(s_p, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(t, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(t_p, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(x_p, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(v, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"    
    deallocate(v_p, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(aux, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"

  end subroutine pre_bicgstab_csr
