program test_prog

  ! ------------------------------------
  ! programme to improve efficiency of 
  ! xweno subroutines
  ! ------------------------------------

  

  implicit none
  
  integer, parameter :: n = 5
  real, dimension(n) :: diag, rhs
  real, dimension(n-1) :: lower_diag, upper_diag
  integer :: info

  ! test modulo

  print*,"  modulo(1,10)= ", modulo(1,10)
  print*,"  modulo(11,10)= ", modulo(11,10)
  print*,"  modulo(0,10)= ", modulo(0,10)
  print*,"  modulo(-1,10)= ", modulo(-1,10)
stop


!!$  ! test NAG  
!!$
!!$  upper_diag = (/ 2.1,  -1.0,   1.9,   8.0 /)
!!$  diag =   (/ 3.0,   2.3,  -5.0,  -0.9,   7.1 /)
!!$  lower_diag = (/ 3.4,   3.6,   7.0,  -6.0 /)
!!$  rhs = (/ 2.7,  -0.5,   2.6,   0.6,   2.7 /)
!!$
!!$
!!$  CALL DGTSV(n,1,lower_diag,diag,upper_diag,rhs,n,info)
!!$
!!$
!!$  IF (INFO.EQ.0) THEN
!!$
!!$     print*,'Solution'
!!$     print*, rhs
!!$
!!$  ELSE
!!$     print*,'The (', info, ',', info, ')', &
!!$          & ' element of the factor U is zero'
!!$  END IF


end program test_prog
