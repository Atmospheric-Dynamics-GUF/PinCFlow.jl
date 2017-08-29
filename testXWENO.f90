program testxweno_prog

  ! ------------------------------------
  ! programme to improve efficiency of 
  ! xweno subroutines
  ! ------------------------------------

  use xweno_module
  use type_module
  use debug_module

  implicit none

  integer, parameter            :: phiSize = 10
  integer, parameter            :: nSamples = 1
  real, dimension(phiSize)      :: phiBar
  real, dimension(phiSize,-1:1) :: phiTilde, phiTilde2


  integer :: i,j,pos


  ! CPU Time
  integer                     :: rate, startTimeCount, timeCount  
  real                        :: cpuTime 

  ! test cases
  logical, parameter :: compare_old_new = .false.    ! reconstruct1D_old vs. reconstruct1D_new

  logical, parameter :: symmetry = .true.            ! symmetry check for reccontruct1D
  logical, parameter :: centered = .false.




  call init_xweno       ! set ILES parameters 


  


if( symmetry ) then

   call initSymm1D(phiBar, phiSize)
   call reconstruct1D  (phiBar, phiSize, phiTilde, centered)
   call checkSymmTilde1D(phiTilde, phiSize)

end if












  if( compare_old_new ) then

     ! init phiBar

     do i = 1,phiSize/2
        phiBar(i) = real(i)**2
     end do
     do i = phiSize/2+1,phiSize
        phiBar(i) = real(i)**4
     end do



     phiTilde = 0.0


     ! run warm
     !  do j = 1,50
     !     call reconstruct1D  (phiBar,phiSize,phiTilde,centered)
     !  end do


     ! ----------------------------------
     !         test new routine
     ! ----------------------------------


     ! init clock
     call system_clock (count_rate=rate)
     call system_clock (count=startTimeCount)

     do j = 1,nSamples
        call reconstruct1D  (phiBar,phiSize,phiTilde,centered)
     end do

     ! stop and print time
     call system_clock(count=timeCount)
     cpuTime = (timeCount - startTimeCount)/real(rate) / nSamples

     print*,"CPU time of reconstruct1D: ", cpuTime

     open( unit=50, file="rec1D_new.txt", action="write", &
          form="formatted", status="replace")

     ! output
     do pos = -1,1
        do i = 1,phiSize
           write(unit=50,fmt="(es25.14)") phiTilde(i,pos)
        end do
     end do

     close(unit=50)

     ! ----------------------------------
     !         test old routine
     ! ----------------------------------


     ! init
     phiTilde = 0.0

     ! use old subroutine

     ! init clock
     call system_clock (count_rate=rate)
     call system_clock (count=startTimeCount)

     do j = 1,nSamples
        call reconstruct1D_old  (phiBar,phiSize,phiTilde,centered)
     end do

     ! stop and print time
     call system_clock(count=timeCount)
     cpuTime = (timeCount - startTimeCount)/real(rate) / nSamples
     print*,"CPU time of reconstruct1D_old: ", cpuTime



     open( unit=60, file="rec1D_old.txt", action="write", &
          form="formatted", status="replace")

     ! output
     do pos = -1,1
        do i = 1,phiSize
           write(unit=60,fmt="(es25.14)") phiTilde(i,pos)
        end do
     end do

     close( unit=60) 


  end if















end program testxweno_prog
