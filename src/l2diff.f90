program diff_prog

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !%   calculates the l2-norm of date stored      %
  !%   in tec360 format in file1 and file2        %
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  use type_module
  use timeScheme_module
  use init_module
  use debug_module
  use output_module
  use xweno_module
  use atmosphere_module
  use boundary_module
  use flux_module
  use muscl_module
  use update_module
  use finish_module
  use algebra_module

   
  !----------------------------------------------------
  implicit none
  !---------------------------------------------------
 
  integer                     :: iTime
  real                        :: time1, time2, dt

  ! CPU Time
  integer                     :: rate, startTimeCount, timeCount  
  real                        :: cpuTime 

  ! fields
  real, dimension(:,:,:,:), allocatable :: var, var0
  ! var(i,j,k,iVar) iVar = 1..5 > rho,u,v,w,pExner

  real, dimension(:,:,:), allocatable :: dRho      ! RK-Update for rho
  real, dimension(:,:,:,:), allocatable :: dMom    ! RK for rhoU,rhoV,rhoW
  
  
  real, dimension(:,:,:,:,:), allocatable :: flux
  ! flux(i,j,k,dir,iFlux) 
  ! dir = 1..3 > f,g,h-flux in x,y,z-direction
  ! iFlux = 1..4 > fRho, fRhoU, rRhoV, fRhoW

  ! topography via force field
  real, dimension(:,:,:,:), allocatable :: force ! volume forces
  ! force(i,j,k,forceVector)
 
  ! output per timeStep
  logical :: output
  real :: nextOutputTime      ! scaled time for next output
  
  ! general
  integer :: i,j,k,l,m

  ! do not scale variables -> SI units
  logical :: scale = .false. 
  
  ! file to be compared
  character(len=40) :: file1 = "hotBubble-2x1x2-000002.dat"
  character(len=40) :: file2 = "hotBubble-2x1x2-000001.dat"

  ! L2 differences
  real :: diffRho, diffU, diffV, diffW, diffP

  ! test field
  real, dimension(3,3,3) :: b

  !-------------------------------------------------
  !                init
  !-------------------------------------------------
  
  ! allocate variables
  call setup (var,var0,flux,force,dRho,dMom)    ! read input.f90 -> allocate var,flux

  call init_atmosphere      ! set atmospheric background state  
!  call initialise (var)     ! set initial conditions
  

!  call init_xweno       ! set ILES parameters 
!  call init_fluxes      ! allocate tilde variables
!  call init_update      ! allocate dp
!  call init_output

  var = 0.0
  !------------------------
  !   read data of file1
  !------------------------

  ! reset reference quantities

  call readtec360(var,file1,time1,scale) 

  print*,"rho = ", var(:,:,:,1)
stop

  var0 = 0.0
  !------------------------
  !   read data of file2
  !------------------------

  call readtec360(var0,file2,time2,scale)


  !---------------------------
  !   output rel. deviations
  !---------------------------

  diffRho = norm3D( var(:,:,:,1) - var0(:,:,:,1) ) / norm3D( var(:,:,:,1))
  write(*,fmt="(a35,es10.3)") "relDiffRho = ", diffRho
  
  diffU = norm3D( var(1:nx,1:ny,1:nz,2) - var0(1:nx,1:ny,1:nz,2) ) / norm3D( var(1:nx,1:ny,1:nz,2))
  write(*,fmt="(a35,es10.3)") "relDiffU = ", diffU
  !  write(*,fmt="(a25,es10.3)") "refDiffU at nx/2, nz/2 = ", ( var(nx/2,1,nz/2,2) - var0(nx/2,1,nz/2,2)) / var(nx/2,1,nz/2,2)

  
!  diffW = norm3D( var(:,:,:,4) - var0(:,:,:,4) ) / norm3D( var(:,:,:,4))
  diffW = norm3D( var(:,:,:,4) )
  write(*,fmt="(a35,es10.3)") "relDiffW = ", diffW
  !  write(*,fmt="(a25,es10.3)") "refDiffW at nx/2, nz/2 = ", ( var(nx/2,1,nz/2,4) - var0(nx/2,1,nz/2,4)) / var(nx/2,1,nz/2,4)


  diffP = norm3D( var(:,:,:,5) - var0(:,:,:,5) ) / norm3D( var(:,:,:,5))
  write(*,fmt="(a35,es10.3)") "relDiffP = ", diffP
  !  write(*,fmt="(a25,es10.3)") "refDiffP at nx/2, nz/2 = ", ( var(nx/2,1,nz/2,5) - var0(nx/2,1,nz/2,5)) / var(nx/2,1,nz/2,5)



  
  !-------------------------------------
  !       deallocate variables
  !-------------------------------------
  
!  call terminate_fluxes
!  call terminate (var,var0,dRho,dMom)                         
!  call terminate_atmosphere
!  call terminate_output
  

end program diff_prog
