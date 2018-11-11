module muscl_module
  !-----------------------------------------------------
  ! MUSCL reconstruction routines 
  ! This module is independent of the rest of the code
  ! 1) it always assumes 1 additional cells on both ends
  ! of the array for interpolation
  !-----------------------------------------------------

  implicit none

  private
  ! all module variables are internal to the module
  ! if not delared otherwise with the public attribute
  
  public :: reconstruct_MUSCL
  public :: muscl_reconstruct1D
  public :: muscl_reconstruct2D
  public :: muscl_reconstruct3D
  public :: limit
  public :: minmod
  
 
 
  

contains
  
  subroutine reconstruct_MUSCL(&
       & uBar, uTilde, &
       & sizeX, sizeY, sizeZ, &
       & limiter)
    
    !-------------------------------------------------------------
    !  MUSCL reconstruction of a single 3D data field u 
    !  of size uSizeX * uSizeY * sizeZ into data field uTilde with
    !  limiter of limterType - in ALL spatial directions 
    !  Note: all fields have a single ghost cell at the boundary
    !-------------------------------------------------------------
    
    ! in/out
    integer, intent(in) :: sizeX, sizeY, sizeZ
    real, dimension(sizeX,sizeY,sizeZ), intent(in) :: uBar
    real, dimension(sizeX,sizeY,sizeZ,1:3,0:1), intent(out) :: uTilde
    ! 1:3 -> spation dimension to be reconstrueced, 1=x, 2=y, 3=z
    ! 0:1 -> 0 = left/right or backward/forward or down/up
    character(len=*), intent(in) :: limiter
    
    ! local fields
    
    ! local vars

    ! ------------- reconstruction to x +/- 1/2 --------------
    
    call muscl_reconstruct3D(&
         & uBar, &
         & sizeX, sizeY, sizeZ, &
         & uTilde, &
         & limiter, "x")


    ! ------------- reconstruction to y +/- 1/2 --------------
    
    call muscl_reconstruct3D(&
         & uBar, &
         & sizeX, sizeY, sizeZ, &
         & uTilde, &
         & limiter, "y")

    ! ------------- reconstruction to z +/- 1/2 --------------
    
    call muscl_reconstruct3D(&
         & uBar, &
         & sizeX, sizeY, sizeZ, &
         & uTilde, &
         & limiter, "z")

  end subroutine reconstruct_MUSCL


!------------------------------------------------------------------
  

  subroutine muscl_reconstruct3D(u, sizeX, sizeY, sizeZ, uTilde, &
       &                         limiterType, direction)

    !-------------------------------------------------------------
    !  MUSCL reconstruction of a single 3D data field u 
    !  of size uSizeX * uSizeY * sizeZ into data field uTilde with
    !  limiter of limterType - in a SINGLE specified direction
    !  Note: all fields have a single ghost cell at the boundary
    !-------------------------------------------------------------

    ! in/out
    integer, intent(in) :: sizeX, sizeY, sizeZ
    real, dimension(sizeX,sizeY,sizeZ), intent(in) :: u
    real, dimension(sizeX,sizeY,sizeZ,1:3,0:1), intent(out) :: uTilde
    ! 1:3 -> spation dimension to be reconstrueced, 1=x, 2=y, 3=z
    ! 0:1 -> 0 = left/right or backward/forward or down/up
    character(len=*), intent(in) :: limiterType
    character(len=1), intent(in) :: direction

    ! local fields
    real, dimension(sizeX)     :: phiX
    real, dimension(sizeX,0:1) :: phiTildeX
    real, dimension(sizeY)     :: phiY
    real, dimension(sizeY,0:1) :: phiTildeY
    real, dimension(sizeZ)     :: phiZ
    real, dimension(sizeZ,0:1) :: phiTildeZ
    
    ! local vars
    integer :: i,j,k
    
    ! debuggung
    logical, parameter :: debugging = .true.

    select case( direction )

    case( 'x' )

       ! reconstruction in x-direction
       do k = 2, sizeZ-1
          do j = 2, sizeY-1
             phiX = u(:,j,k)
             call muscl_reconstruct1D( phiX, sizeX, phiTildeX, limiterType)
             uTilde(:,j,k,1,:) = phiTildeX
          end do
       end do
       
    case( 'y' )

       ! reconstruction in y-direction
       do k = 2, sizeZ-1
          do i = 2, sizeX-1
             phiY = u(i,:,k)
             call muscl_reconstruct1D(phiY, sizeY, phiTildeY, limiterType)
             uTilde(i,:,k,2,:) = phiTildeY
          end do
       end do


    case( 'z' )

       ! reconstruction in z-direction
       do j = 2, sizeY-1
          do i = 2, sizeX-1
             phiZ = u(i,j,:)
             call muscl_reconstruct1D(phiZ, sizeZ, phiTildeZ, limiterType)
             uTilde(i,j,:,3,:) = phiTildeZ
          end do
       end do

       
    case default
       stop "musc_reconstruct3D: unknown direction.Stop."
    end select


  end subroutine muscl_reconstruct3D


  ! --------------------------------------------------------------------------------------


  subroutine muscl_reconstruct2D( u, uSizeX, uSizeY, uTilde, limiterType, direction )

    !-------------------------------------------------------------
    !  MUSCL reconstruction of a single 2D data field phi 
    !  of size uSizeX times uSizeY into data field phiTilde with
    !  limiter of limterType - only in a specified direction
    !-------------------------------------------------------------

    ! in/out
    integer, intent(in) :: uSizeX, uSizeY
    real, dimension(uSizeX,uSizeY), intent(in) :: u
    real, dimension(uSizeX,uSizeY,1:2,0:1), intent(out) :: uTilde
    ! 1:2 -> direction of reconstruction, 0:1 -> left/right or backward/forward 
    character(len=*), intent(in) :: limiterType
    character(len=1), intent(in) :: direction

    ! local fields
    real, dimension(uSizeX)     :: phiX
    real, dimension(uSizeX,0:1) :: phiTildeX
    real, dimension(uSizeY)     :: phiY
    real, dimension(uSizeY,0:1) :: phiTildeY
          
    ! local vars
    integer :: i,j,k
    
    ! init uTilde
    uTilde = 1000.0

    select case( direction )

    case( 'x' )

       ! reconstruction in x-direction
       do j = 2, uSizeY-1
          phiX = u(:,j)
          call muscl_reconstruct1D( phiX, uSizeX, phiTildeX, limiterType)
          uTilde(:,j,1,:) = phiTildeX
       end do

    case( 'y' )

       ! reconstruction in y-direction
       do i = 2, uSizeX-1
          phiY = u(i,:)
          call muscl_reconstruct1D(phiY, uSizeY, phiTildeY, limiterType)
          uTilde(i,:,2,:) = phiTildeY
       end do

    case default
       stop "muscl.f90/reconstruct2D: unknown direction.Stop."
    end select


  end subroutine muscl_reconstruct2D


  ! --------------------------------------------------------------------------------------


  subroutine muscl_reconstruct1D( phi, phiSize, phiTilde, limiterType )
    
    !-------------------------------------
    !  MUSCL reconstruction of a single
    !  1D data field phi of size phiSize
    !  into data field phiTilde with
    !  limiter of limterType
    !-------------------------------------


    ! in/out
    integer, intent(in) :: phiSize
    real, dimension(phiSize), intent(in)      :: phi
    real, dimension(phiSize,0:1), intent(out) :: phiTilde
    character(len=*), intent(in) :: limiterType

    ! local varibales
    real :: sigmaR, sigmaL, sigma
    real :: deltaR, deltaL, theta
      
    ! local vars
    integer :: i,j,k
    
    ! init phiTilde with 
    phiTilde = 1000.0
    
    do i = 2, phiSize-1

       
       deltaL = phi(i) - phi(i-1)
       deltaR = phi(i+1) - phi(i)
       
       if (deltaR == 0.) then
           phiTilde(i,1) = phi(i)
           phiTilde(i,0) = phi(i)
          else if(deltaL == 0.) then
           theta = deltaL / deltaR
       
           phiTilde(i,1) = phi(i) + 0.5*limit(theta,limiterType) * deltaR
           phiTilde(i,0) = phi(i)
          else
           theta = deltaL / deltaR
       
           phiTilde(i,1) = phi(i) + 0.5*limit(theta,limiterType) * deltaR
           phiTilde(i,0) = phi(i) &
                           - 0.5*limit(1./theta, limiterType) * deltaL
       end if
       
       ! there is an alternative formulation without calculation of 1/theta
       ! ONLY valid for symmetric limiter functions...
       ! ref. ALGORITMY 2009, Cada paper

    end do

  end subroutine muscl_reconstruct1D


  ! ---------------------------------------------------------------------------------------


  function limit(theta, limiterType )  result(sigma)

    !--------------------------------------------
    !  definition of various limiter functions
    !  in normal form
    !--------------------------------------------

    ! in/out
    real, intent(in)              :: theta
    character(len=*), intent(in)  :: limiterType
    real                          :: sigma

    ! local variables 
    real             :: s


    select case( limiterType )
       
    case( 'minmod')
       sigma = minmod(theta, 1.0)

    case( 'MCVariant' ) 
       s = (2.0 + theta)/3.0
       sigma = max(0.0, min(2*theta,s, 2.0))
       
    case( 'Cada' )
       s = (2.0 + theta)/3.0
       sigma = max(0.0,min(s,max(-theta/2.0,min(2*theta,s,1.5))))
       
    case default
       stop "muscl.f90/limit: unknown limiter type.Stop."        
    end select


  end function limit


  ! --------------------------------------------------------------------------------------------


  function minmod(a,b) result(c)
    
    !------------------------------
    !  classical minmod limiter
    !------------------------------
    
    ! in/out
    real, intent(in) :: a,b
    real :: c


    if( a*b > 0.0 ) then
       if( abs(a) < abs(b) ) then
          c = a
       else
          c = b
       end if
    else
       c = 0.0
    end if

  end function minmod


  ! --------------------------------------------------------------------------------------------


end module muscl_module
