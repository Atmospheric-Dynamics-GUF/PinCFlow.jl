module boundary_module
  
  use type_module
  use flux_module        ! BC for rhoTilde, uTilde,...
  use atmosphere_module
  use mpi_module

  implicit none

  private

  !-------------------
  !  public routines
  !-------------------
  public :: setBoundary
  public :: setBoundary_x_periodic
  public :: setBoundary_y_periodic


  !--------------------
  !  private routines
  !--------------------
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

      if ( idim > 1 ) then
         ! boundary conditions taken care of by setHalos
        else
         call setBoundary_x_periodic(var,flux,option)
      endif

    case default
       stop "setBoundary: unknown case xBoundary"
    end select

    
    !------------------------------
    !          y-direction
    !------------------------------
    select case( yBoundary )
       
    case( "periodic" ) 

      if ( jdim > 1 ) then
         ! boundary conditions taken care of by setHalos
        else
         call setBoundary_y_periodic(var,flux,option)
      endif

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


  !-----------------------------------------------------------------------


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

          if(verbose .and. master) then
             print*,"horizontalBoundary: x-horizontal BC for rho set."
          end if
       end if

       if (timeScheme == "semiimplicit" .or. auxil_equ) then
          ! density fluctuations -> iVar = 6
          do i = 1,nbx
             var(nx+i,:,:,6) = var(i,:,:,6)
             var(-i+1,:,:,6) = var(nx-i+1,:,:,6)
          end do

          if(verbose .and. master) then
             print*,"horizontalBoundary: x-horizontal BC for rhop set."
          end if
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

          if (timeScheme == "semiimplicit") then
             print*,'ERROR: updateTheta = .true. not allowed for &
                   & timeScheme == "semiimplicit"'
             stop
          end if

          do i = 1,nbx
             var(nx+i,:,:,6) = var(i,:,:,6)
             var(-i+1,:,:,6) = var(nx-i+1,:,:,6)
          end do

          if(verbose .and. master ) print*,"horizontalBoundary: &
               & x-horizontal BC for theta set."
       end if


       if( predictMomentum ) then
          ! velocity u (staggered along x) -> iVar = 2

          var(0,:,:,2) = var(nx,:,:,2)

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

          if (verbose .and. master) then
             print *,"boundary.f90/horizontalBoundary: &
                    & x-horizontal BC for u, v, w set."
          end if
       end if


       if( correctMomentum ) then
          ! pressure Variable -> iVar = 5
          var(nx+1,:,:,5) = var(1,:,:,5)            ! right ghost cell
          var(0,:,:,5) = var(nx,:,:,5)             ! left ghost cells

          if (verbose .and. master) then
              print*,"boundary.f90/horizontalBoundary: &
                    & x-horizontal BC for p set."
          end if
       end if
       
    case( "ice")
     ! ice variables -> iVar = nVar-3, nVar
        do iVar = nVar-3,nVar
            do i = 1,nbx
               var(nx+i,:,:,iVar) = var(i,:,:,iVar)
               var(-i+1,:,:,iVar) = var(nx-i+1,:,:,iVar)
            end do
        end do

        if(verbose .and. master) print*,"horizontalBoundary: &
               & x-horizontal BC for ice variables set."   
               
    case( "iceTilde" )
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

    case( "iceFlux")
    ! nothing
               
    case( "varTilde" ) 
       !-----------------------------------
       !             varTilde
       !-----------------------------------
       
       ! the following three boundary-condition calls can probably be 
       ! removed
       ! probably only necessary for ALDM (that is not used any more 
       ! anyway)

       if( updateMass ) then
          ! reconstructed density needed in ghost cell i = nx+2
          rhoTilde(nx+2,:,:,1,0) = rhoTilde(2,:,:,1,0)

          ! ...in ghost cell i = -1
          rhoTilde(-1,:,:,1,1) = rhoTilde(nx-1,:,:,1,1)

          if(verbose .and. master) then
             print*,"horizontalBoundary: x-horizontal BC for rhoTilde set."
          end if
       end if

       if( timeScheme == "semiimplicit" .or. auxil_equ ) then
          ! reconstructed density fluctuation needed in ghost cell i = nx+2
          rhopTilde(nx+2,:,:,1,0) = rhopTilde(2,:,:,1,0)

          ! ...in ghost cell i = -1
          rhopTilde(-1,:,:,1,1) = rhopTilde(nx-1,:,:,1,1)

          if(verbose .and. master) then
             print*,"horizontalBoundary: x-horizontal BC for rhoTilde set."
          end if
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

          if (verbose .and. master) then
             print*,"horizontalBoundary: &
                   & x-horizontal BC for thetaTilde set."
          end if
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

          if (verbose .and. master) then
             print *,"horizontalBoundary: y-horizontal BC for rho set."
          end if
       end if

       if (timeScheme == "semiimplicit" .or. auxil_equ) then
          ! density fluctuations -> iVar = 6
          do j = 1,nby
             var(:,ny+j,:,6) = var(:,j,:,6)
             var(:,-j+1,:,6) = var(:,ny-j+1,:,6)
          end do

          if (verbose .and. master) then
             print *,"horizontalBoundary: y-horizontal BC for rho set."
          end if
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
          ! potential temperature -> iVar = 6

          if (timeScheme == "semiimplicit") then
             print*,'ERROR: updateTheta = .true. not allowed for &
                   & timeScheme == "semiimplicit"'
             stop
          end if

          do j = 1,nby
             var(:,ny+j,:,6) = var(:,j,:,6)
             var(:,-j+1,:,6) = var(:,ny-j+1,:,6)
          end do

          if (verbose .and. master) then
             print *,"horizontalBoundary: y-horizontal BC for theta set."
          end if
       end if

       if( predictMomentum ) then
          ! velocity v (staggared along y) -> iVar = 3

          var(:,0,:,3) = var(:,ny,:,3)

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

          if (verbose .and. master) then
             print *,"boundary.f90/horizontalBoundary: &
                   & y-horizontal BC for u, v, w set."
          end if
       end if

       if( correctMomentum ) then
          ! pressure variable -> iVar = 5
          var(:,ny+1,:,5) = var(:,1,:,5)    ! forward ghost cell
          var(:,0,:,5) = var(:,ny,:,5)      ! backward

          if (verbose .and. master) then
             print *,"boundary.f90/horizontalBoundary: &
                   & y-horizontal BC for p set."
          end if
       end if

       
    case( "ice")
    ! ice variables -> iVar = nVar-3, nVar
        do iVar = nVar-3,nVar
            do j = 1,nby
               var(:,ny+j,:,iVar) = var(:,j,:,iVar)
               var(:,-j+1,:,iVar) = var(:,ny-j+1,:,iVar)
            end do
        end do                
        
        if(verbose .and. master) print *,"horizontalBoundary: &
               & y-horizontal BC for ice variables set."
               
    case( "iceTilde" )
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
               
    case( "iceFlux")
    ! nothing
               

    case( "varTilde" ) 
       !-----------------------------------
       !              varTilde
       !-----------------------------------
       
       ! the following three boundary-condition calls can probably be 
       ! removed
       ! probably only necessary for ALDM (that is not used any more 
       ! anyway)

       if( updateMass ) then
          ! reconstructed density needed in ghost cell j = ny+2
          rhoTilde(:,ny+2,:,2,0) = rhoTilde(:,2,:,2,0)

          ! ...in ghost cell j = -1
          rhoTilde(:,-1,:,2,1) = rhoTilde(:,ny-1,:,2,1)

          if (verbose .and. master) then
             print*,"horizontalBoundary: y-horizontal BC for rhoTilde set."
          end if
       end if

       if( timeScheme == "semiimplicit" .or. auxil_equ ) then
          ! reconstructed density fluctuations needed in ghost cell 
          ! j = ny+2
          rhopTilde(:,ny+2,:,2,0) = rhopTilde(:,2,:,2,0)

          ! ...in ghost cell j = -1
          rhopTilde(:,-1,:,2,1) = rhopTilde(:,ny-1,:,2,1)

          if (verbose .and. master) then
             print*,"horizontalBoundary: y-horizontal BC for rhoTilde set."
          end if
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
          thetaTilde(:,ny+2,:,2,0) = thetaTilde(:,2,:,2,0)

          ! ...in ghost cell j = -1
          thetaTilde(:,-1,:,2,1) = thetaTilde(:,ny-1,:,2,1)

          if (verbose .and. master) then
             print*,"horizontalBoundary: &
                  & y-horizontal BC for thetaTilde set."
          end if
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


!-------------------------------------------------------------------------


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

          if (verbose .and. master) then
             print *,"setBoundary_z_periodic: z-periodic BC for rho set."
          end if
       end if

       if( timeScheme == "semiimplicit" .or. auxil_equ ) then
          ! density fluctuations -> iVar = 6
          do k = 1,nbz
             var(:,:,nz+k,6) = var(:,:,k,6)
             var(:,:,-k+1,6) = var(:,:,nz-k+1,6)
          end do

          if (verbose .and. master) then
             print *,"setBoundary_z_periodic: z-periodic BC for rhop set."
          end if
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

          if (timeScheme == "semiimplicit") then
             print*,'ERROR: updateTheta = .true. not allowed for &
                   & timeScheme == "semiimplicit"'
             stop
          end if

          do k = 1,nbz
             var(:,:,nz+k,6) = var(:,:,k,6)
             var(:,:,-k+1,6) = var(:,:,nz-k+1,6)
          end do

          if (verbose .and. master) then
             print *,"setBoundary_z_periodic: &
                    & z-perdiodic BC for theta set."
          end if
       end if


       if( predictMomentum ) then
          ! velocity w (staggared along z) -> iVar = 4

          var(:,:,0,4) = var(:,:,nz,4)

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

          if (verbose .and. master) then
             print *,"setBoundary_z_periodic: &
                    & z-periodic BC for u, v, w set."
          end if
       end if

       if( correctMomentum ) then
          ! pressure variable -> iVar = 5
          var(:,:,nz+1,5) = var(:,:,1,5)    ! forward ghost cell
          var(:,:,0,5) = var(:,:,nz,5)      ! backward

          if (verbose .and. master) then
             print *,"setBoundary_z_periodic: &
                    & z-periodic BC for p set."
          end if
       end if


    case( "ice" )
    ! ice variables -> iVar = nVar-3, nVar
        do iVar = nVar-3,nVar
            do k = 1,nbz
               var(:,:,nz+k,iVar) = var(:,:,k,iVar)
               var(:,:,-k+1,iVar) = var(:,:,nz-k+1,iVar)
            end do
        end do

        if(verbose .and. master) print *,"setBoundary_z_periodic: &
               & z-periodic BC for ice variables set."
               
    case( "iceTilde" )
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
               
    case( "iceFlux")
    ! nothing
       
       
    case( "varTilde" ) 
       !-----------------------------------
       !              varTilde
       !-----------------------------------
       
       ! the following three boundary-condition calls can probably be 
       ! removed
       ! probably only necessary for ALDM (that is not used any more 
       ! anyway)

       if( updateMass ) then
          ! reconstructed density needed in ghost cell k = nz+2
          rhoTilde(:,:,nz+2,3,0) = rhoTilde(:,:,2,3,0)

          ! ...in ghost cell j = -1
          rhoTilde(:,:,-1,3,1) = rhoTilde(:,:,nz-1,3,1)

          if (verbose .and. master) then
             print*,"setBoundary_z_periodic: &
                   & z-periodic BC for rhoTilde set."
          end if
       end if

       if( timeScheme == "semiimplicit" .or. auxil_equ ) then
          ! reconstructed density needed in ghost cell k = nz+2
          rhopTilde(:,:,nz+2,3,0) = rhopTilde(:,:,2,3,0)

          ! ...in ghost cell j = -1
          rhopTilde(:,:,-1,3,1) = rhopTilde(:,:,nz-1,3,1)

          if (verbose .and. master) then
             print*,"setBoundary_z_periodic: &
                   & z-periodic BC for rhoTilde set."
          end if
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
          thetaTilde(:,:,nz+2,3,0) = thetaTilde(:,:,2,3,0)

          ! ...in ghost cell k = -1
          thetaTilde(:,:,-1,3,1) = thetaTilde(:,:,nz-1,3,1)

          if (verbose .and. master) then
             print*,"setBoundary_z_periodic: &
                   & z-periodic BC for thetaTilde set."
          end if
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

!-------------------------------------------------------------------------

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
          ! reflect at boundary with change of sign
          ! rho -> iVar = 1

          ! in case of baroclinic life-cycle simulation only the 
          ! deviations from the environmental state are reflected

          if (testCase == "baroclinic_LC") then 
             if (fluctuationMode) then
                var(1:nx,1:ny,0,1) &
                = - var(1:nx,1:ny,1,1) + dens_env_pp(:,:,1) - rhoStrat(1) &
                  + dens_env_pp(:,:,0) - rhoStrat(0)

                var(1:nx,1:ny,nz+1,1) &
                = - var(1:nx,1:ny,nz,1) &
                  + dens_env_pp(:,:,nz) - rhoStrat(nz) &
                  + dens_env_pp(:,:,nz+1) - rhoStrat(nz+1)
               else
                var(1:nx,1:ny,0,1) &
                = - var(1:nx,1:ny,1,1) + dens_env_pp(:,:,1) &
                  + dens_env_pp(:,:,0)

                var(1:nx,1:ny,nz+1,1) &
                = - var(1:nx,1:ny,nz,1) + dens_env_pp(:,:,nz) &
                  + dens_env_pp(:,:,nz+1)
             end if
            else
             do k = 1,nbz
                var(:,:,-k+1,1) = -var(:,:,k,1)
                var(:,:,nz+k,1) = -var(:,:,nz-k+1,1)
             end do
          end if

          if (timeScheme == "semiimplicit" .or. auxil_equ) then
             ! vertical boundary condition for density fluctuations

             if (testCase == "baroclinic_LC") then 
                var(1:nx,1:ny,0,6) &
                = - var(1:nx,1:ny,1,6) + dens_env_pp(:,:,1) - rhoStrat(1) &
                  + dens_env_pp(:,:,0) - rhoStrat(0)

                var(1:nx,1:ny,nz+1,6) &
                = - var(1:nx,1:ny,nz,6) &
                  + dens_env_pp(:,:,nz) - rhoStrat(nz) &
                  + dens_env_pp(:,:,nz+1) - rhoStrat(nz+1)
               else
                do k = 1,nbz
                   var(:,:,-k+1,6) = -var(:,:,k,6)
                   var(:,:,nz+k,6) = -var(:,:,nz-k+1,6)
                end do
             end if
          end if
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

          if (timeScheme == "semiimplicit") then
             print*,'ERROR: updateTheta = .true. not allowed for &
                   & timeScheme == "semiimplicit"'
             stop
          end if

          if (testCase == 'baroclinic_LC') then 
             print*,'ERROR: updateTheta = .true. not allowed for &
                   & testCase == "baroclinic_LC"'
             stop
          end if

          do k = 1,nbz
             var(:,:,-k+1,6) = var(:,:,k,6)
             var(:,:,nz+k,6) = var(:,:,nz-k+1,6)
          end do
       end if

       if( predictMomentum ) then
          ! w -> set to zero at bound, 
          !      reflect at bound with change of sign

          var(:,:,0,4) = 0.0
          var(:,:,nz,4) = 0.0

          do k = 1,nbz
             var(:,:,-k,4) = -var(:,:,k,4)
             var(:,:,nz+k,4) = -var(:,:,nz-k,4)
          end do

          ! transverse velocities u,v: 
          ! -> general case: reflect at bound. and change sign (no slip)
          !    life-cycle simulation: only deviations from environmental 
          !                           state reflected,
          !                           sign not changed (free slip)

          if (testCase == "baroclinic_LC") then 
             ! u
             var(1:nx,1:ny,0,2) &
             = var(1:nx,1:ny,1,2) - u_env_pp (:,:,1) + u_env_pp (:,:,0) 
             var(1:nx,1:ny,nz+1,2) &
             = var(1:nx,1:ny,nz,2) - u_env_pp (:,:,nz) + u_env_pp (:,:,nz+1)

             ! v
             var(1:nx,1:ny,0,3) &
             = var(1:nx,1:ny,1,3) - v_env_pp (:,:,1) + v_env_pp (:,:,0) 
             var(1:nx,1:ny,nz+1,3) &
             = var(1:nx,1:ny,nz,3) - v_env_pp (:,:,nz) + v_env_pp (:,:,nz+1)
            else
             do k = 1,nbz
                ! u
                var(:,:,-k+1,2) = -var(:,:,k,2)
                var(:,:,nz+k,2) = -var(:,:,nz-k+1,2)

                ! v
                var(:,:,-k+1,3) = -var(:,:,k,3)
                var(:,:,nz+k,3) = -var(:,:,nz-k+1,3)
             end do
          end if
       end if

       if( correctMomentum ) then
          ! set gradient at vertical boundary to 0
          ! life-cycle simulation: only gradient of deviations from 
          !                        environmental pressure are set to 0

          if (testCase == "baroclinic_LC") then 
             ! z = 0
             var(1:nx,1:ny,0,5) &
             = var(1:nx,1:ny,1,5) - p_env_pp(:,:,1) + p_env_pp(:,:,0)

             ! z = zMax
             var(1:nx,1:ny,nz+1,5) &
             = var(1:nx,1:ny,nz,5) - p_env_pp(:,:,nz) + p_env_pp(:,:,nz+1)  
            else
             ! z = 0
             var(:,:,0,5) = var(:,:,1,5)       

             ! z = zMax
             var(:,:,nz+1,5) = var(:,:,nz,5)   
          end if
       end if

       ! additional call to setHalos in case baroclinic_LC in order to set 
       ! the overlap betwen vertical and horizontal halos right

       if (testCase == "baroclinic_LC") call setHalos( var, "var" )


    case( "ice" )
     ! reflect at boundary without change of sign
     ! ice variables -> iVar = nVar-3, nVar
        do iVar = nVar-3, nVar
            do k = 1,nbz
              var(:,:,-k+1,iVar) = var(:,:,k,iVar)
              var(:,:,nz+k,iVar) = var(:,:,nz-k+1,iVar)
            end do
        end do
        
    case( "iceTilde" )
    !nothing
        

    case( "iceFlux" )
    ! set vertical fluxes at wall to 0
          ! analog to mass flux modification above

          ! ice variables iVar=nVar-3,nVar
   !       do iVar = nVar-3, nVar
   !         flux(:,:,0,3,iVar) = 0.0
   !         flux(:,:,nz,3,iVar) = 0.0
   !       end do
          
   !       if (verbose .and. master) print*,"boundary.f90/verticalBoundary: &
   !            &vertical BC for ice variables set"
               
          ! replace flux by CDS fluxes at upper / lower region
   !       if( iceFluxCorr ) then


             ! vertical ice fluxes in the bottom region
    !         do k = 1, nbCellCorr
    !           do j = 1,ny
    !             do i = 1,nx
    !               do iVar = nVar-3, nVar
    !                 if( fluctuationMode ) then
    !                   rhoU = var(i,j,k+1,1) + rhoStrat(k+1)
    !                   rhoD = var(i,j,k,1)   + rhoStrat(k)
    !                 else
    !                   rhoU = var(i,j,k+1,1)
    !                   rhoD = var(i,j,k,1)
    !                 end if
    !                 rhoU = rhoU*var(i,j,k+1,iVar) ! reuse of above variables
    !                 ! names do not necessarily refer to rho
    !                 rhoD = rhoD*var(i,j,k,iVar)   
    !                 wSurf = var(i,j,k,4)
    !                 hRho = wSurf * 0.5*(rhoD + rhoU)
    !                 flux(i,j,k,3,iVar) = hRho
    !               end do
    !             end do
    !           end do
    !         end do

             ! vertical ice fluxes at the top region
    !         do k = nz-1,nz-nbCellCorr,-1
    !           do j = 1,ny
    !             do i = 1,nx
    !               do iVar = nVar-3, nVar 
    !                 if( fluctuationMode ) then
    !                   rhoU = var(i,j,k+1,1) + rhoStrat(k+1)
    !                   rhoD = var(i,j,k,1)   + rhoStrat(k)
    !                 else
    !                   rhoU = var(i,j,k+1,1)
    !                   rhoD = var(i,j,k,1)
    !                 end if
    !                 rhoU = rhoU*var(i,j,k+1,iVar)
    !                 rhoD = rhoD*var(i,j,k,iVar)
    !                 wSurf = var(i,j,k,4)
    !                 hRho = wSurf * 0.5*(rhoD + rhoU)
    !                 flux(i,j,k,3,iVar) = hRho
    !               end do
    !             end do
    !           end do
    !         end do
    !      end if ! iceFluxCorr

    
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

          ! these settings correspond to solid-wall condition (w = 0) at 
          ! the bottom and top, whence advective and turbulent mass 
          ! fluxes must vanish there 
          ! a different b.c. seems appropriate if molecular diffusion is 
          ! taken into account!

          ! density
          flux(:,:,0,3,1) = 0.0
          flux(:,:,nz,3,1) = 0.0

          ! diffusive potential-temperature fluxes
          flux(:,:,0,3,5) = 0.0
          flux(:,:,nz,3,5) = 0.0

          if (verbose .and. master) then
             print*,"boundary.f90/verticalBoundary: &
                    &vertical BC for rho set"
          end if

          ! replace flux by CDS fluxes at upper / lower region
          if( rhoFluxCorr ) stop 'ERROR: rhoFluxCorr = .false. expected'
       end if ! updateMass

       
       !if( updateIce ) then
          ! set vertical fluxes at wall to 0
          ! analog to mass flux modification above
          
       if (timeScheme == "semiimplicit" .or. auxil_equ) then
          ! density fluctuations
          flux(:,:,0,3,6) = 0.0
          flux(:,:,nz,3,6) = 0.0
       end if

          ! ice variables iVar=nVar-3,nVar
       !   do iVar = nVar-3, nVar
       !     flux(:,:,0,3,iVar) = 0.0
       !     flux(:,:,nz,3,iVar) = 0.0
       !   end do
          
       !   if (verbose .and. master) print*,"boundary.f90/verticalBoundary: &
       !        &vertical BC for ice variables set"
       !        
          ! replace flux by CDS fluxes at upper / lower region
       !   if( iceFluxCorr ) then

             ! vertical ice fluxes in the bottom region
       !      do k = 1, nbCellCorr
       !        do j = 1,ny
       !          do i = 1,nx
       !            do iVar = nVar-3, nVar
       !              if( fluctuationMode ) then
       !                rhoU = var(i,j,k+1,1) + rhoStrat(k+1)
       !                rhoD = var(i,j,k,1)   + rhoStrat(k)
       !              else
       !                rhoU = var(i,j,k+1,1)
       !                rhoD = var(i,j,k,1)
       !              end if
       !              rhoU = rhoU*var(i,j,k+1,iVar) ! reuse of above variables
       !              ! names do not necessarily refer to rho
       !              rhoD = rhoD*var(i,j,k,iVar)   
       !              wSurf = var(i,j,k,4)
       !              hRho = wSurf * 0.5*(rhoD + rhoU)
       !              flux(i,j,k,3,iVar) = hRho
       !            end do
       !          end do
       !        end do
       !      end do

             ! vertical ice fluxes at the top region
       !      do k = nz-1,nz-nbCellCorr,-1
       !        do j = 1,ny
       !          do i = 1,nx
       !            do iVar = nVar-3, nVar 
       !              if( fluctuationMode ) then
       !                rhoU = var(i,j,k+1,1) + rhoStrat(k+1)
       !                rhoD = var(i,j,k,1)   + rhoStrat(k)
       !              else
       !                rhoU = var(i,j,k+1,1)
       !                rhoD = var(i,j,k,1)
       !              end if
       !              rhoU = rhoU*var(i,j,k+1,iVar)
       !              rhoD = rhoD*var(i,j,k,iVar)
       !             wSurf = var(i,j,k,4)
       !              hRho = wSurf * 0.5*(rhoD + rhoU)
       !              flux(i,j,k,3,iVar) = hRho
       !            end do
       !          end do
       !        end do
       !      end do
       !   end if ! iceFluxCorr
       !end if ! updateIce
       
       
       if( updateTheta ) then
          if (timeScheme == "semiimplicit") then
             print*,'ERROR: updateTheta = .true. not allowed for &
                   & timeScheme == "semiimplicit"'
             stop
          end if

          ! set fluxes accros solid wall boundary -> 0

          flux(:,:,0,3,6) = 0.0
          flux(:,:,nz,3,6) = 0.0

          if (verbose .and. master) then
             print*,"verticalBoundary: vertical BC for theta set"
          end if

          if( thetaFluxCorr ) stop 'ERROR: thetaFluxCorr = .false. expected'
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

          if (verbose .and. master) then
             print*,"boundary.f90/verticalBoundary: &
                    &vertical flux-BC for u,v,w set"
          end if

          ! replace flux by CDS fluxes at upper / lower region

          if( uFluxCorr ) stop 'ERROR: uFluxCorr = .false. expected'

          if( vFluxCorr ) stop 'ERROR: vFluxCorr = .false. expected'

          if( wFluxCorr ) stop 'ERROR: wFluxCorr = .false. expected'
       end if ! predictMomentum
       
    case default
       stop "setBoundary_z: unknown option."
    end select

  end subroutine setBoundary_z_solidWall

end module boundary_module
