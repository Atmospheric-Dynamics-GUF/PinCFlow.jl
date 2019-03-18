module boundary_module
  
  use type_module
  use flux_module        ! BC for rhoTilde, uTilde,...
  use atmosphere_module

  implicit none

  private

  !-------------------
  !  public routines
  !-------------------
  public :: setBoundary


  !--------------------
  !  private routines
  !--------------------
  private:: setBoundary_x_periodic
  private ::setBoundary_y_periodic
  private ::setBoundary_z_periodic
  private ::setBoundary_z_solidWall 


contains

  subroutine setBoundary (var,flux,option)
    !-------------------------------------------
    ! 1) set values in ghost cells according to 
    !    boundary conditions
    ! 2) set fluxes explicitly to zero at walls
    !-------------------------------------------

    ! in/out variables      
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), &
         & intent(inout) :: flux
    character(len=*), intent(in) :: option
    
    !------------------------------
    !          x-direction
    !------------------------------
    select case( xBoundary )
       
    case( "periodic" )

       ! modified by Junhong Wei (20161109) *** starting line ***
       
!       call setBoundary_x_periodic(var,flux,option)

      if ( idim > 1 ) then
!        if (master) print*, "periodic x-boundaries implied by MPI halo exchange"   ! modified by Junhong Wei (20170216)
      else
        call setBoundary_x_periodic(var,flux,option)
      endif

      ! modified by Junhong Wei (20161109) *** finishing line ***
      
    case default
       stop "setBoundary: unknown case xBoundary"
    end select

    
    !------------------------------
    !          y-direction
    !------------------------------
    select case( yBoundary )
       
    case( "periodic" ) 

       ! modified by Junhong Wei (20161109) *** starting line ***

!       call setBoundary_y_periodic(var,flux,option)

      if ( jdim > 1 ) then
!        if (master) print*, "periodic y-boundaries implied by MPI halo exchange"   ! modified by Junhong Wei (20170216)
      else
        call setBoundary_y_periodic(var,flux,option)
      endif

      ! modified by Junhong Wei (20161109) *** finishing line ***
       
    case default
       stop "setBoundary: unknown case yBoundary"
    end select

    
    !------------------------------
    !          z-direction
    !------------------------------
    select case( zBoundary )
       
    case( "periodic" ) 
       call setBoundary_z_periodic(var,flux,option)
       
    case( "solid_wall" )
       call setBoundary_z_solidWall(var,flux,option)
       
    case default
       stop "setBoundary: unknown case zBoundary"
    end select



  end subroutine setBoundary


  !----------------------------------------------------------------------------


  subroutine setBoundary_x_periodic(var,flux,option)
    !-------------------------------------------
    ! 1) set values in ghost cells according to 
    !    boundary conditions
    ! 2) set fluxes explicitly to zero at walls
    !-------------------------------------------

    ! in/out variables      
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), &
         & intent(inout) :: flux
    character(len=*), intent(in) :: option

    ! local variables
    integer :: i,j,k,iVar

    real, dimension(-nby:ny+nby,-nbz:nz+nbz) :: uBound

    select case( option )

    case( "var" )
       !-----------------------------------
       !                var
       !-----------------------------------

       if( updateMass ) then
          ! density -> iVar = 1
          do i = 1,nbx
             var(nx+i,:,:,1) = var(i,:,:,1)
             var(-i+1,:,:,1) = var(nx-i+1,:,:,1)
          end do

! modified by Junhong Wei (20161109) *** starting line ***
          
          if(verbose .and. master) print*,"horizontalBoundary: &
               & x-horizontal BC for rho set."

! modified by Junhong Wei (20161109) *** finishing line ***          
          
       end if


      if( updateIce ) then
          ! ice variables -> iVar = nVar-3, nVar
          do iVar = nVar-3,nVar
            do i = 1,nbx
               var(nx+i,:,:,iVar) = var(i,:,:,iVar)
               var(-i+1,:,:,iVar) = var(nx-i+1,:,:,iVar)
            end do
          end do

          if(verbose .and. master) print*,"horizontalBoundary: &
               & x-horizontal BC for ice variables set."        
          
       end if


       if( updateTheta ) then
          ! pot. temp.  -> iVar = 6
          do i = 1,nbx
             var(nx+i,:,:,6) = var(i,:,:,6)
             var(-i+1,:,:,6) = var(nx-i+1,:,:,6)
          end do

! modified by Junhong Wei (20161109) *** starting line ***
          
          if(verbose .and. master ) print*,"horizontalBoundary: &
               & x-horizontal BC for theta set."

! modified by Junhong Wei (20161109) *** finishing line ***          
          
       end if


       if( predictMomentum ) then
          ! velocity u (staggered along x) -> iVar = 2

!         achatzb
!         uBound causes problems in restarts, hence removed
!         uBound = 0.5*( var(0,:,:,2) + var(nx,:,:,2) )  ! velocity at bound
!         var(0,:,:,2) = uBound
!         var(nx,:,:,2) = uBound
          var(0,:,:,2) = var(nx,:,:,2)
!         achatze

          do i = 1,nbx
             ! velocity u (staggered along x) -> iVar = 2
             var(nx+i,:,:,2) = var(i,:,:,2)
             var(-i,:,:,2) = var(nx-i,:,:,2)
             ! velocity v -> iVar = 3                    ! ghost cells
             var(nx+i,:,:,3) = var(i,:,:,3)              ! right
             var(-i+1,:,:,3) = var(nx-i+1,:,:,3)         ! left 
             ! velocity w -> iVar = 4
             var(nx+i,:,:,4) = var(i,:,:,4)
             var(-i+1,:,:,4) = var(nx-i+1,:,:,4)
          end do

! modified by Junhong Wei (20161109) *** starting line ***
          
          if(verbose .and. master) print *,"boundary.f90/horizontalBoundary: &
               & x-horizontal BC for u, v, w set."

! modified by Junhong Wei (20161109) *** finishing line ***

       end if


       if( correctMomentum ) then
          ! pressure Variable -> iVar = 5
          var(nx+1,:,:,5) = var(1,:,:,5)            ! right ghost cell
          var(0,:,:,5) = var(nx,:,:,5)             ! left ghost cells

! modified by Junhong Wei (20161109) *** starting line ***

          if(verbose .and. master) print*,"boundary.f90/horizontalBoundary: &
               & x-horizontal BC for p set."

! modified by Junhong Wei (20161109) *** finishing line ***

       end if



    case( "varTilde" ) 
       !-----------------------------------
       !             varTilde
       !-----------------------------------
       
       if( updateMass ) then
          ! reconstructed density needed in ghost cell i = nx+2
          rhoTilde(nx+2,:,:,1,0) = rhoTilde(2,:,:,1,0)

          ! ...in ghost cell i = -1
          rhoTilde(-1,:,:,1,1) = rhoTilde(nx-1,:,:,1,1)

! modified by Junhong Wei (20161109) *** starting line ***
          
          if(verbose .and. master) print*,"horizontalBoundary: &
               & x-horizontal BC for rhoTilde set."

! modified by Junhong Wei (20161109) *** finishing line ***          
          
       end if


       if ( updateIce ) then
         ! reconstructed density needed in ghost cell i = nx+2
         nAerTilde(nx+2,:,:,1,0) = nAerTilde(2,:,:,1,0)
         nIceTilde(nx+2,:,:,1,0) = nIceTilde(2,:,:,1,0)
         qIceTilde(nx+2,:,:,1,0) = qIceTilde(2,:,:,1,0)
         qvTilde(nx+2,:,:,1,0) = qvTilde(2,:,:,1,0)

         ! ...in ghost cell i = -1
         nAerTilde(-1,:,:,1,1) = nAerTilde(nx-1,:,:,1,1)
         nIceTilde(-1,:,:,1,1) = nIceTilde(nx-1,:,:,1,1)
         qIceTilde(-1,:,:,1,1) = qIceTilde(nx-1,:,:,1,1)
         qvTilde(-1,:,:,1,1) = qvTilde(nx-1,:,:,1,1)

         if(verbose .and. master) print*,"horizontalBoundary: &
               & x-horizontal BC for iceTilde variables set."

       end if


       if( updateTheta ) then
          ! reconstructed density needed in ghost cell i = nx+2
          thetaTilde(nx+2,:,:,1,0) = thetaTilde(2,:,:,1,0)

          ! ...in ghost cell i = -1
          thetaTilde(-1,:,:,1,1) = thetaTilde(nx-1,:,:,1,1)

! modified by Junhong Wei (20161109) *** starting line ***          
          
          if(verbose .and. master) print*,"horizontalBoundary: &
               & x-horizontal BC for thetaTilde set."

! modified by Junhong Wei (20161109) *** finishing line ***          
          
       end if


    case( "flux" )
       !-----------------------------------
       !              flux
       !-----------------------------------
       
       return
       


    case default
       stop "setBoundary_x: unknown option."
    end select





  end subroutine setBoundary_x_periodic


  !--------------------------------------------------------------------


  subroutine setBoundary_y_periodic(var,flux,option)
    !-------------------------------------------
    ! 1) set values in ghost cells according to 
    !    boundary conditions
    ! 2) set fluxes explicitly to zero at walls
    !-------------------------------------------

    ! in/out variables      
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), &
         & intent(inout) :: flux
    character(len=*), intent(in) :: option

    ! local variables
    integer :: i,j,k,iVar

    real, dimension(-nbx:nx+nbx,-nbz:nz+nbz) :: vBound


    select case( option )


    case( "var" )
       !-----------------------------------
       !                var
       !-----------------------------------

       if( updateMass ) then
          ! density -> iVar = 1
          do j = 1,nby
             var(:,ny+j,:,1) = var(:,j,:,1)
             var(:,-j+1,:,1) = var(:,ny-j+1,:,1)
          end do

! modified by Junhong Wei (20161109) *** starting line ***                    
          
          if(verbose .and. master) print *,"horizontalBoundary: &
               & y-horizontal BC for rho set."

! modified by Junhong Wei (20161109) *** finishing line ***                    
          
       end if


       if( updateIce ) then
          ! ice variables -> iVar = nVar-3, nVar
          do iVar = nVar-3,nVar
            do j = 1,nby
               var(:,ny+j,:,iVar) = var(:,j,:,iVar)
               var(:,-j+1,:,iVar) = var(:,ny-j+1,:,iVar)
            end do
          end do                
          
          if(verbose .and. master) print *,"horizontalBoundary: &
               & y-horizontal BC for ice variables set."
      
       end if


       if( updateTheta ) then
          ! density -> iVar = 6
          do j = 1,nby
             var(:,ny+j,:,6) = var(:,j,:,6)
             var(:,-j+1,:,6) = var(:,ny-j+1,:,6)
          end do

! modified by Junhong Wei (20161109) *** starting line ***          
          
          if(verbose .and. master) print *,"horizontalBoundary: &
               & y-horizontal BC for theta set."

! modified by Junhong Wei (20161109) *** finishing line ***          
          
       end if


       if( predictMomentum ) then
          ! velocity v (staggared along y) -> iVar = 3

!         achatzb
!         use of vBound causes problems in restarts, hence removed
!         vBound = 0.5*( var(:,0,:,3) + var(:,ny,:,3) )
!         var(:,0,:,3) = vBound
!         var(:,ny,:,3) = vBound
          var(:,0,:,3) = var(:,ny,:,3)
!         achatze

          do j = 1,nby
             ! velocity u -> iVar = 2
             var(:,ny+j,:,2) = var(:,j,:,2)
             var(:,-j+1,:,2) = var(:,ny-j+1,:,2)
             ! velocity v (staggared along y) -> iVar = 3
             var(:,ny+j,:,3) = var(:,j,:,3)
             var(:,-j,:,3) = var(:,ny-j,:,3)
             ! velocity w -> iVar = 4
             var(:,ny+j,:,4) = var(:,j,:,4)
             var(:,-j+1,:,4) = var(:,ny-j+1,:,4)
          end do

! modified by Junhong Wei (20161109) *** starting line ***          
          
          if(verbose .and. master) print *,"boundary.f90/horizontalBoundary: &
               & y-horizontal BC for u, v, w set."

! modified by Junhong Wei (20161109) *** finishing line ***          
          
       end if


       if( correctMomentum ) then
          ! pressure variable -> iVar = 5
          var(:,ny+1,:,5) = var(:,1,:,5)    ! forward ghost cell
          var(:,0,:,5) = var(:,ny,:,5)      ! backward

! modified by Junhong Wei (20161109) *** starting line ***
          
          if(verbose .and. master) print *,"boundary.f90/horizontalBoundary: &
               & y-horizontal BC for p set."

! modified by Junhong Wei (20161109) *** finishing line ***          
          
       end if


    case( "varTilde" ) 
       !-----------------------------------
       !              varTilde
       !-----------------------------------
       
       if( updateMass ) then
          ! reconstructed density needed in ghost cell j = ny+2
!          rhoTilde(:,ny+2,:,1,0) = rhoTilde(:,2,:,1,0)   ! modified by Junhong Wei (suggestion from Mark Schlutow)
          rhoTilde(:,ny+2,:,2,0) = rhoTilde(:,2,:,2,0)   ! modified by Junhong Wei (suggestion from Mark Schlutow)

          ! ...in ghost cell j = -1
!          rhoTilde(:,-1,:,1,1) = rhoTilde(:,ny-1,:,1,1)   ! modified by Junhong Wei (suggestion from Mark Schlutow)
          rhoTilde(:,-1,:,2,1) = rhoTilde(:,ny-1,:,2,1)   ! modified by Junhong Wei (suggestion from Mark Schlutow)

! modified by Junhong Wei (20161109) *** starting line ***
          
          if(verbose .and. master) print*,"horizontalBoundary: &
               & y-horizontal BC for rhoTilde set."

! modified by Junhong Wei (20161109) *** finishing line ***
          
       end if


       if ( updateIce ) then
         ! see above, similar to rho
         nAerTilde(:,ny+2,:,2,0) = nAerTilde(:,2,:,2,0)
         nIceTilde(:,ny+2,:,2,0) = nIceTilde(:,2,:,2,0)
         qIceTilde(:,ny+2,:,2,0) = qIceTilde(:,2,:,2,0)
         qvTilde(:,ny+2,:,2,0) = qvTilde(:,2,:,2,0)

         ! see above, similar to rho
         nAerTilde(:,-1,:,2,1) = nAerTilde(:,ny-1,:,2,1)   
         nIceTilde(:,-1,:,2,1) = nIceTilde(:,ny-1,:,2,1)   
         qIceTilde(:,-1,:,2,1) = qIceTilde(:,ny-1,:,2,1)   
         qvTilde(:,-1,:,2,1) = qvTilde(:,ny-1,:,2,1)   

         if(verbose .and. master) print*,"horizontalBoundary: &
               & y-horizontal BC for iceTilde variables set."

       end if


       if( updateTheta ) then
          ! reconstructed density needed in ghost cell j = ny+2
!          thetaTilde(:,ny+2,:,1,0) = thetaTilde(:,2,:,1,0)   ! modified by Junhong Wei (suggestion from Mark Schlutow)
          thetaTilde(:,ny+2,:,2,0) = thetaTilde(:,2,:,2,0)   ! modified by Junhong Wei (suggestion from Mark Schlutow)

          ! ...in ghost cell j = -1
!          thetaTilde(:,-1,:,1,1) = thetaTilde(:,ny-1,:,1,1)   ! modified by Junhong Wei (suggestion from Mark Schlutow)
          thetaTilde(:,-1,:,2,1) = thetaTilde(:,ny-1,:,2,1)   ! modified by Junhong Wei (suggestion from Mark Schlutow)

! modified by Junhong Wei (20161109) *** starting line ***
          
          if(verbose .and. master) print*,"horizontalBoundary: &
               & y-horizontal BC for thetaTilde set."

! modified by Junhong Wei (20161109) *** finishing line ***
          
       end if


    case( "flux" )
       !-----------------------------------
       !              flux
       !-----------------------------------
       return
       

    case default
       stop "setBoundary_y: unknown option."
    end select


  end subroutine setBoundary_y_periodic


!-------------------------------------------------------------------------------


  subroutine setBoundary_z_periodic(var,flux,option)
    !-------------------------------------------
    ! 1) set values in ghost cells according to 
    !    boundary conditions
    ! 2) set fluxes explicitly to zero at walls
    !-------------------------------------------

    ! in/out variables      
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), &
         & intent(inout) :: flux
    character(len=*), intent(in) :: option

    ! local variables
    integer :: i,j,k,iVar

    real, dimension(-nbx:nx+nbx,-nby:ny+nby) :: wBound

    select case( option )
       
    case( "var" )
       !-----------------------------------
       !                var
       !-----------------------------------

       if( updateMass ) then
          ! density -> iVar = 1
          do k = 1,nbz
             var(:,:,nz+k,1) = var(:,:,k,1)
             var(:,:,-k+1,1) = var(:,:,nz-k+1,1)
          end do

! modified by Junhong Wei (20161109) *** starting line ***          
          
          if(verbose .and. master) print *,"setBoundary_z_periodic: &
               & z-periodic BC for rho set."

! modified by Junhong Wei (20161109) *** finishing line ***
          
       end if


       if( updateIce ) then
          ! ice variables -> iVar = nVar-3, nVar
          do iVar = nVar-3,nVar
            do k = 1,nbz
               var(:,:,nz+k,iVar) = var(:,:,k,iVar)
               var(:,:,-k+1,iVar) = var(:,:,nz-k+1,iVar)
            end do
          end do

          if(verbose .and. master) print *,"setBoundary_z_periodic: &
               & z-periodic BC for ice variables set."
        
       end if


       if( updateTheta ) then
          ! density -> iVar = 6
          do k = 1,nbz
             var(:,:,nz+k,6) = var(:,:,k,6)
             var(:,:,-k+1,6) = var(:,:,nz-k+1,6)
          end do

! modified by Junhong Wei (20161109) *** starting line ***
          
          if(verbose .and. master) print *,"setBoundary_z_periodic: &
               & z-perdiodic BC for theta set."

! modified by Junhong Wei (20161109) *** finishing line ***
          
       end if


       if( predictMomentum ) then
          ! velocity w (staggared along z) -> iVar = 4

!         achatzb
!         use of wBound causes problems in restarts, hence removed
!         wBound = 0.5*( var(:,:,0,4) + var(:,:,nz,4) )
!         var(:,:,0,4) = wBound
!         var(:,:,nz,4) = wBound
          var(:,:,0,4) = var(:,:,nz,4)
!         achatze

          do k = 1,nbz
             ! velocity u -> iVar = 2
             var(:,:,nz+k,2) = var(:,:,k,2)
             var(:,:,-k+1,2) = var(:,:,nz-k+1,2)
             ! velocity v -> iVar = 3
             var(:,:,nz+k,3) = var(:,:,k,3)
             var(:,:,-k+1,3) = var(:,:,nz-k+1,3)
             ! velocity w (staggared along z) -> iVar = 4
             var(:,:,nz+k,4) = var(:,:,k,4)
             var(:,:,-k,4) = var(:,:,nz-k,4)
          end do

! modified by Junhong Wei (20161109) *** starting line ***
          
          if(verbose .and. master) print *,"setBoundary_z_periodic: &
               & z-periodic BC for u, v, w set."

! modified by Junhong Wei (20161109) *** finishing line ***
          
       end if
       

       if( correctMomentum ) then
          ! pressure variable -> iVar = 5
          var(:,:,nz+1,5) = var(:,:,1,5)    ! forward ghost cell
          var(:,:,0,5) = var(:,:,nz,5)      ! backward

! modified by Junhong Wei (20161109) *** starting line ***
          
          if(verbose .and. master) print *,"setBoundary_z_periodic: &
               & z-periodic BC for p set."

! modified by Junhong Wei (20161109) *** finishing line ***
          
       end if


    case( "varTilde" ) 
       !-----------------------------------
       !              varTilde
       !-----------------------------------
       
       if( updateMass ) then
          ! reconstructed density needed in ghost cell k = nz+2
!          rhoTilde(:,:,nz+2,1,0) = rhoTilde(:,:,2,1,0)   ! modified by Junhong Wei (suggestion from Mark Schlutow)
          rhoTilde(:,:,nz+2,3,0) = rhoTilde(:,:,2,3,0)   ! modified by Junhong Wei (suggestion from Mark Schlutow)

          ! ...in ghost cell j = -1
!          rhoTilde(:,:,-1,1,1) = rhoTilde(:,:,nz-1,1,1)   ! modified by Junhong Wei (suggestion from Mark Schlutow)
          rhoTilde(:,:,-1,3,1) = rhoTilde(:,:,nz-1,3,1)   ! modified by Junhong Wei (suggestion from Mark Schlutow)

! modified by Junhong Wei (20161109) *** starting line ***
          
          if(verbose .and. master) print*,"setBoundary_z_periodic: &
               & z-periodic BC for rhoTilde set."

! modified by Junhong Wei (20161109) *** finishing line ***
          
       end if


       if ( updateIce ) then
         ! see above, similar to rho
         nAerTilde(:,:,nz+2,3,0) = nAerTilde(:,:,2,3,0) 
         nIceTilde(:,:,nz+2,3,0) = nIceTilde(:,:,2,3,0) 
         qIceTilde(:,:,nz+2,3,0) = qIceTilde(:,:,2,3,0) 
         qvTilde(:,:,nz+2,3,0) = qvTilde(:,:,2,3,0) 

         ! see above, similar to rho
         nAerTilde(:,:,-1,3,1) = nAerTilde(:,:,nz-1,3,1)  
         nIceTilde(:,:,-1,3,1) = nIceTilde(:,:,nz-1,3,1) 
         qIceTilde(:,:,-1,3,1) = qIceTilde(:,:,nz-1,3,1)  
         qvTilde(:,:,-1,3,1) = qvTilde(:,:,nz-1,3,1)  

         if(verbose .and. master) print*,"setBoundary_z_periodic:  &
               & z-periodic BC for iceTilde variables set."

       end if


       if( updateTheta ) then
          ! reconstructed density needed in ghost cell k = nz+2
!          thetaTilde(:,:,nz+2,1,0) = thetaTilde(:,:,2,1,0)   ! modified by Junhong Wei (suggestion from Mark Schlutow)
          thetaTilde(:,:,nz+2,3,0) = thetaTilde(:,:,2,3,0)   ! modified by Junhong Wei (suggestion from Mark Schlutow)

          ! ...in ghost cell k = -1
!          thetaTilde(:,:,-1,1,1) = thetaTilde(:,:,nz-1,1,1)   ! modified by Junhong Wei (suggestion from Mark Schlutow)
          thetaTilde(:,:,-1,3,1) = thetaTilde(:,:,nz-1,3,1)   ! modified by Junhong Wei (suggestion from Mark Schlutow)

! modified by Junhong Wei (20161109) *** starting line ***
          
          if(verbose .and. master) print*,"setBoundary_z_periodic: &
               & z-periodic BC for thetaTilde set."

! modified by Junhong Wei (20161109) *** finishing line ***
          
       end if


    case( "flux" )
       !-----------------------------------
       !              flux
       !-----------------------------------
       return
       

    case default
       stop "setBoundary_y: unknown option."
    end select


  end subroutine setBoundary_z_periodic


!-------------------------------------------------------------------------------


  subroutine setBoundary_z_solidWall(var,flux,option)
    !-------------------------------------------
    ! 1) set values in ghost cells according to 
    !    boundary conditions
    ! 2) set fluxes explicitly to zero at walls
    !-------------------------------------------

    ! in/out variables      
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), &
         & intent(inout) :: flux
    character(len=*), intent(in) :: option

    ! local variables
    integer :: i,j,k,iVar
    real, dimension(-nby:ny+nby,-nbz:nz+nbz) :: uBound
    real, dimension(-nbx:nx+nbx,-nbz:nz+nbz) :: vBound

    ! rho flux correction
    real :: rhoU, rhoD, hRho, wSurf

    ! uvw flux correction
    real :: rhoEdge
    real :: hRhoU, hRhoV, hRhoW
    real :: uD, uU, wL, wR
    real :: vD, vU, wB, wF, wD, wU

    ! theta flux correction
    real :: thetaU, thetaD, hTheta

    ! uvw flux correction
    real :: thetaEdge
    real :: hThetaU, hThetaV, hThetaW

    
    select case( option )

    case( "var" )
       !-----------------------------------
       !                var
       !-----------------------------------


       if( updateMass ) then
          ! reflect at boundary without change of sign
          ! rho -> iVar = 1
          do k = 1,nbz
             var(:,:,-k+1,1) = var(:,:,k,1)
             var(:,:,nz+k,1) = var(:,:,nz-k+1,1)
          end do
       end if

       if( updateIce ) then
          ! reflect at boundary without change of sign
          ! ice variables -> iVar = nVar-3, nVar
          do iVar = nVar-3, nVar
            do k = 1,nbz
              var(:,:,-k+1,iVar) = var(:,:,k,iVar)
              var(:,:,nz+k,iVar) = var(:,:,nz-k+1,iVar)
            end do
          end do
       end if

       if( updateTheta ) then
          ! reflect at boundary withOUT change of sign
          ! theta -> iVar = 6
          do k = 1,nbz
             var(:,:,-k+1,6) = var(:,:,k,6)
             var(:,:,nz+k,6) = var(:,:,nz-k+1,6)
          end do
       end if


       if( predictMomentum ) then
          ! transverse velocities u,v: 
          ! -> reflect at bound. and change sign (no slip)

          ! w -> set to zero at bound, 
          !      reflect at bound with change of sign

!         achatzb
          var(:,:,0,4) = 0.0
          var(:,:,nz,4) = 0.0
!         achatze

          do k = 1,nbz
             ! u
             var(:,:,-k+1,2) = -var(:,:,k,2)
             var(:,:,nz+k,2) = -var(:,:,nz-k+1,2)
             ! v
             var(:,:,-k+1,3) = -var(:,:,k,3)
             var(:,:,nz+k,3) = -var(:,:,nz-k+1,3)
             ! w 
             var(:,:,-k,4) = -var(:,:,k,4)
             var(:,:,nz+k,4) = -var(:,:,nz-k,4)
          end do
       end if


       if( correctMomentum ) then
          ! set gradient at vertical boundary to 0
          var(:,:,0,5) = var(:,:,1,5)       ! z = 0
          var(:,:,nz+1,5) = var(:,:,nz,5)   ! z = zMax
       end if


       !---------------------------------------------------------------         
       
    case( "varTilde " )
       !-----------------------------------
       !             varTilde
       !-----------------------------------
       
       return



    case( "flux" )
       !-----------------------------------
       !                flux
       !-----------------------------------
       
       if( updateMass ) then
          ! set vertical fluxes at wall to 0

          ! density
          flux(:,:,0,3,1) = 0.0
          flux(:,:,nz,3,1) = 0.0

          ! modified by Junhong Wei (20161109) *** starting line ***
          
          if (verbose .and. master) print*,"boundary.f90/verticalBoundary: &
               &vertical BC for rho set"

          ! modified by Junhong Wei (20161109) *** finishing line ***
          
          ! replace flux by CDS fluxes at upper / lower region
          if( rhoFluxCorr ) then

             ! vertical mass flux at the bottom region
             do k = 1, nbCellCorr
                do j = 1,ny
                   do i = 1,nx

                      if( fluctuationMode ) then
                         rhoU = var(i,j,k+1,1) + rhoStrat(k+1)
                         rhoD = var(i,j,k,1)   + rhoStrat(k)
                      else
                         rhoU = var(i,j,k+1,1)
                         rhoD = var(i,j,k,1)
                      end if
                      
                      wSurf = var(i,j,k,4)

                      hRho = wSurf * 0.5*(rhoD + rhoU)

                      flux(i,j,k,3,1) = hRho
                   end do
                end do
             end do

             ! vertical mass flux at the top region
             do k = nz-1,nz-nbCellCorr,-1
                do j = 1,ny
                   do i = 1,nx

                      if( fluctuationMode ) then
                         rhoU = var(i,j,k+1,1) + rhoStrat(k+1)
                         rhoD = var(i,j,k,1)   + rhoStrat(k)
                      else
                         rhoU = var(i,j,k+1,1)
                         rhoD = var(i,j,k,1)
                      end if
                      
                      wSurf = var(i,j,k,4)

                      hRho = wSurf * 0.5*(rhoD + rhoU)

                      flux(i,j,k,3,1) = hRho
                   end do
                end do
             end do
          end if ! rhoFluxCorr
       end if ! updateMass


       if( updateTheta ) then

          ! set fluxes accros solid wall boundary -> 0
          flux(:,:,0,3,6) = 0.0
          flux(:,:,nz,3,6) = 0.0

          ! modified by Junhong Wei (20161109) *** starting line ***
          
          if (verbose .and. master) print*,"verticalBoundary: &
               &vertical BC for theta set"

          ! modified by Junhong Wei (20161109) *** finishing line ***          
          
          if( thetaFluxCorr ) then

             ! replace flux by CDS fluxes at upper / lower region

             ! vertical mass flux at the bottom region
             do k = 1, nbCellCorr
                do j = 1,ny
                   do i = 1,nx

                      thetaU = var(i,j,k+1,6)
                      thetaD = var(i,j,k,6)
                      wSurf = var(i,j,k,4)

                      hTheta = wSurf * 0.5*(thetaD + thetaU)

                      flux(i,j,k,3,6) = hTheta
                   end do
                end do
             end do


             ! vertical mass flux at the top region
             do k = nz-1,nz-nbCellCorr,-1
                do j = 1,ny
                   do i = 1,nx

                      thetaU = var(i,j,k+1,6)
                      thetaD = var(i,j,k,6)
                      wSurf = var(i,j,k,4)

                      hTheta = wSurf * 0.5*(thetaD + thetaU)

                      flux(i,j,k,3,6) = hTheta
                   end do
                end do
             end do
          end if ! thetaFluxCorr
       end if ! updateTheta


       if( predictMomentum ) then

          ! momentum rho*u
          flux(:,:,0,3,2) = 0.0
          flux(:,:,nz,3,2) = 0.0

          ! momentum rho*v
          flux(:,:,0,3,3) = 0.0
          flux(:,:,nz,3,3) = 0.0

          ! momentum rho*w
          flux(:,:,-1,3,4) = 0.0
          flux(:,:,nz,3,4) = 0.0

          ! modified by Junhong Wei (20161109) *** starting line ***
          
          if (verbose .and. master) print*,"boundary.f90/verticalBoundary: &
               &vertical flux-BC for u,v,w set"

          ! modified by Junhong Wei (20161109) *** finishing line ***
          
          ! replace flux by CDS fluxes at upper / lower region

          if( uFluxCorr ) then

             ! vertical flux hRhoU at bottom boundary
             do k = 1,nbCellCorr
                do j = 1,ny
                   do i = 0,nx

                      ! consistent with density interpolation in conti eq
                      rhoEdge = 0.25*( var(i,j,k,1) &
                           & + var(i,j,k+1,1) &
                           & + var(i+1,j,k,1) &
                           & + var(i+1,j,k+1,1)  )
                      
                      if( fluctuationMode ) rhoEdge = rhoEdge + rhoStratTilde(k)
                      
                      uD = var(i,j,k+1,2)
                      uU = var(i,j,k,2)
                      wL = var(i+1,j,k,4)
                      wR = var(i,j,k,4)

                      hRhoU = 0.25*(uD + uU)*(wL + wR)

                      flux(i,j,k,3,2) = rhoEdge * hRhoU
                   end do
                end do
             end do

             ! vertical flux hRhoU at top boundary
             do k = nz-1,nz-nbCellCorr
                do j = 1,ny
                   do i = 0,nx

                      ! consistent with density interpolation in conti eq
                      rhoEdge = 0.25*( var(i,j,k,1) &
                           & + var(i,j,k+1,1) &
                           & + var(i+1,j,k,1) &
                           & + var(i+1,j,k+1,1)  )
                      
                      if( fluctuationMode ) rhoEdge = rhoEdge + rhoStratTilde(k)
                      
                      uD = var(i,j,k+1,2)
                      uU = var(i,j,k,2)
                      wL = var(i+1,j,k,4)
                      wR = var(i,j,k,4)

                      hRhoU = 0.25*(uD + uU)*(wL + wR)

                      flux(i,j,k,3,2) = rhoEdge * hRhoU
                   end do
                end do
             end do

          end if

          !------------------------------------------------------------

          if( vFluxCorr ) then

             ! vertical flux correction hRhoV at bottom
             !            do k = 0,nz
             do k = 1,nbCellCorr
                do j = 0,ny
                   do i = 1,nx

                      ! consistent with density interpolation in conti eq
                      rhoEdge = 0.25*(var(i,j,k,1) &
                           & + var(i,j,k+1,1) &
                           & + var(i,j+1,k,1) &
                           & + var(i,j+1,k+1,1) )
                      
                      if( fluctuationMode ) rhoEdge = rhoEdge + rhoStratTilde(k)
                      
                      vD = var(i,j,k+1,3)
                      vU = var(i,j,k,3)
                      wB = var(i,j+1,k,4)
                      wF = var(i,j,k,4)

                      hRhoV = 0.25*(vD+vU)*(wB+wF)

                      flux(i,j,k,3,3) = rhoEdge * hRhoV
                   end do
                end do
             end do

             ! vertical flux hRhoV at top
             !            do k = 0,nz
             do k = nz-1,nz-nbCellCorr
                do j = 0,ny
                   do i = 1,nx

                      ! consistent with density interpolation in conti eq
                      rhoEdge = 0.25*(var(i,j,k,1) &
                           & + var(i,j,k+1,1) &
                           & + var(i,j+1,k,1) &
                           & + var(i,j+1,k+1,1) )
                      
                      if( fluctuationMode ) rhoEdge = rhoEdge + rhoStratTilde(k)
                      
                      vD = var(i,j,k+1,3)
                      vU = var(i,j,k,3)
                      wB = var(i,j+1,k,4)
                      wF = var(i,j,k,4)

                      hRhoV = 0.25*(vD+vU)*(wB+wF)

                      flux(i,j,k,3,3) = rhoEdge * hRhoV
                   end do
                end do
             end do

          end if


          !----------------------------------------------------------


          if( wFluxCorr ) then

             ! vertical flux hRhoW at bottom
             do k = 0, 1-nbCellCorr
                do j = 1,ny
                   do i = 1,nx

                      ! consistent with density interpolation in conti eq
                      if( fluctuationMode ) then
                         rhoEdge = 0.25*( var(i,j,k,1) + rhoStrat(k) &
                              & + var(i,j,k+1,1) + rhoStrat(k+1) &
                              & + var(i,j,k+1,1) + rhoStrat(k+1) &
                              & + var(i,j,k+2,1) + rhoStrat(k+2) )

                      else
                         rhoEdge = 0.25*( var(i,j,k,1) &
                              & + var(i,j,k+1,1) &
                              & + var(i,j,k+1,1) &
                              & + var(i,j,k+2,1) )
                      end if

                      wD = var(i,j,k+1,4)
                      wU = var(i,j,k,4)

                      hRhoW = 0.25*(wD+wU)**2

                      flux(i,j,k,3,4) = rhoEdge * hRhoW
                   end do
                end do
             end do

             ! vertical flux hRhoW at top
             do k = nz-1,nz-nbCellCorr
                do j = 1,ny
                   do i = 1,nx

                      ! consistent with density interpolation in conti eq
                      if( fluctuationMode ) then
                         rhoEdge = 0.25*( var(i,j,k,1) + rhoStrat(k) &
                              & + var(i,j,k+1,1) + rhoStrat(k+1) &
                              & + var(i,j,k+1,1) + rhoStrat(k+1) &
                              & + var(i,j,k+2,1) + rhoStrat(k+2) )
                         
                      else
                         rhoEdge = 0.25*( var(i,j,k,1) &
                              & + var(i,j,k+1,1) &
                              & + var(i,j,k+1,1) &
                              & + var(i,j,k+2,1) )
                      end if
                      
                      wD = var(i,j,k+1,4)
                      wU = var(i,j,k,4)

                      hRhoW = 0.25*(wD+wU)**2
                      
                      flux(i,j,k,3,4) = rhoEdge * hRhoW
                   end do
                end do
             end do
          end if ! wFluxCorr
       end if ! predictMomentum
       
    case default
       stop "setBoundary_z: unknown option."
    end select


  end subroutine setBoundary_z_solidWall



end module boundary_module
