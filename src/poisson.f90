module poisson_module

  use type_module
  use mpi_module   ! modified by Junhong Wei (20161106)
  use timeScheme_module
  use atmosphere_module
  use algebra_module


  implicit none

  private
  ! all module objects (variables, functions, subroutines)
  ! are internal to the module by default


  !------------------------
  !   public subroutines
  !------------------------
  public :: momentumCorrector
  public :: init_poisson
  public :: terminate_poisson

  !------------------------
  !   private subroutines       
  !------------------------
  private :: getIndex
  private :: poissonSolver_csr
  private :: pressureBoundaryCondition
  private :: correctorStep
  private :: linOpr
  private :: linOprXYZ
  private :: bicgstab
  private :: hypre                     ! modified by Junhong Wei
  private :: gcr
  private :: poissonSolver
  private :: adi
  private :: adi_z
  private :: thomas   ! tridiagonal matrix algoritm (TDMA) / Thomas algorithm




  !-------------------------------
  !    private module variables
  !------------------------------

  ! pressure correction
  real, dimension(:,:,:), allocatable :: dp

  ! solution to Poisson problem
  real, dimension(:,:,:), allocatable :: sol_old1, sol_old2

  ! predicted pressure
  real, dimension(:,:,:), allocatable :: p_pred

  ! tolerance for initial divergence cleaning
  real, parameter :: tolInitial = 1.0e-9

  real :: tol
    


contains

  subroutine momentumCorrector( var,dMom,dt,errFlagBicg, nIter, m,opt)
    ! -------------------------------------------------
    !              correct uStar and p
    ! -------------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,3), &
         & intent(inout) :: dMom

    real, intent(in)                :: dt
    logical, intent(out)            :: errFlagBicg
    integer, intent(out)            :: nIter
    integer, intent(in)             :: m
    character(len=*), intent(in)    :: opt


    ! local variables
    real, dimension(1:nx,1:ny,1:nz) :: rhs        ! RHS
    logical :: onlyinfo 
   


    ! Note: dp is a module variable
    ! calc dp 
    select case( storageType )

    case( "opr" )

       ! Calc RHS of Poisson problem
       onlyinfo = .false.
       call calc_RHS( rhs,var,dt,onlyinfo )
       
       call poissonSolver( rhs,var,dt,errFlagBicg, nIter, m, opt )

    case( "csr" )
       call poissonSolver_csr( var, dt, opt )

    case default
       stop "update.f90/momentumCorrector: unknown storageType"
    end select

    ! set horizontal and vertical BC for dp
    call pressureBoundaryCondition

    ! correct p and uStar with dp
    call correctorStep (var, dMom, dt, m)
    
    ! give info on final divergence
    onlyinfo = .true.
    call calc_RHS( rhs,var,dt,onlyinfo )
       


  end subroutine momentumCorrector


  !---------------------------------------------------------------------------


  subroutine linOpr( var, sIn, Ls )
    ! --------------------------------------
    !   Linear Operator in Poisson problem
    !   Functions as A*x
    ! --------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         &intent(in) :: var
    real, dimension(1:nx,1:ny,1:nz), intent(out) :: Ls
    real, dimension(1:nx,1:ny,1:nz), intent(in)  :: sIn

    ! local field (extended by ghost cells)
    real, dimension(0:nx+1,0:ny+1,0:nz+1) :: s

    ! modified by Junhong Wei (20161106) *** starting line ***

    ! auxiliary fields for "dp"
    real, dimension(ny,nz) :: xSliceLeft_send, xSliceRight_send
    real, dimension(ny,nz) :: xSliceLeft_recv, xSliceRight_recv

    real, dimension(nx,nz) :: ySliceBack_send, ySliceForw_send
    real, dimension(nx,nz) :: ySliceBack_recv, ySliceForw_recv

    ! modified by Junhong Wei (20161106) *** finishing line ***

    real, dimension(-nby:ny+nby,-nbz:nz+nbz) :: dvar


    ! local variables
    integer :: i,j,k
    real :: pStratU, pStratD, rhoEdge
    real :: AL,AR, AB,AF, AD,AU, AC
    real :: sL,sR, sB,sF, sD,sU, sC

    real :: dx2, dy2, dz2

    ! modified by Junhong Wei (20161106) *** starting line ***

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

!   achatzb
    if(topography) stop 'linOpr not ready for topography!'
!   achatze

    ! modified by Junhong Wei (20161106) *** finishing line ***

    ! auxiliary variables
    dx2 = 1.0/dx**2
    dy2 = 1.0/dy**2
    dz2 = 1.0/dz**2

    ! work with auxiliary field s
    s(1:nx,1:ny,1:nz) = sIn

    ! modified by Junhong Wei (20161106) *** starting line ***

    ! Find neighbour procs
    if (idim > 1) call mpi_cart_shift(comm,0,1,left,right,ierror)
    if (jdim > 1) call mpi_cart_shift(comm,1,1,back,forw,ierror)

    ! modified by Junhong Wei (20161106) *** finishing line ***


    select case( model ) 

       !----------------------------------------
       !       Pseudo-incompressible model
       !----------------------------------------
    case( "pseudo_incompressible" )

       ! modified by Junhong Wei (20161106) *** starting line ***

!       !----------------------------
!       !     Boundary conditions
!       !----------------------------
!
!       ! perdiodic in x
!       s(0,:,:) = s(nx,:,:)
!       s(nx+1,:,:) = s(1,:,:)
!
!       ! periodic in y
!       s(:,0,:) = s(:,ny,:)
!       s(:,ny+1,:) = s(:,1,:)

       !----------------------------
       !   set Halo cells: xSlice
       !----------------------------

       if( xBoundary == "periodic" ) then

         ! slice size
         sendcount = ny*nz
         recvcount = sendcount

         ! read slice into contiguous array
         xSliceLeft_send (:,:) = s(1, 1:ny,1:nz)
         xSliceRight_send(:,:) = s(nx,1:ny,1:nz)

         if ( idim > 1 ) then

            ! left -> right
            source = left
            dest = right
            tag = 100

            call mpi_sendrecv(xSliceRight_send(1,1), sendcount, &
                 & mpi_double_precision, dest, tag, &
                 & xSliceLeft_recv(1,1), recvcount, mpi_double_precision, &
                 & source, mpi_any_tag, comm, sts_left, ierror)

            ! right -> left
            source = right
            dest = left
            tag = 100

            call mpi_sendrecv(xSliceLeft_send(1,1), sendcount, &
                 & mpi_double_precision, dest, tag, &
                 & xSliceRight_recv(1,1), recvcount, mpi_double_precision, &
                 & source, mpi_any_tag, comm, sts_right, ierror)

            ! right halos
            s(nx+1,1:ny,1:nz) = xSliceRight_recv(1:ny,1:nz)
            !          var(nx+i,1:ny,1:nz,iVar) = xSliceRight_recv(i,1:ny,1:nz)

            ! left halos
            s(0,1:ny,1:nz) = xSliceLeft_recv(1:ny,1:nz)
            !             var(-nbx+i,1:ny,1:nz,iVar) = xSliceLeft_recv(i,1:ny,1:nz)

         else

           s(0,:,:) = s(nx,:,:)
           s(nx+1,:,:) = s(1,:,:)

         end if

       else
         stop "Poisson: unknown case xBoundary"
       endif

       if(verbose .and. master) print*,"horizontalHalos: &
            & x-horizontal halos copied."


       !------------------------------
       !   set Halo cells: ySlice 
       !------------------------------

       if( yBoundary == "periodic" ) then

         ! slice size
         sendcount = nx*nz
         recvcount = sendcount

         ! read slice into contiguous array
         ySliceBack_send(:,:) = s(1:nx, 1,1:nz)
         ySliceForw_send(:,:) = s(1:nx,ny,1:nz)

         if( jdim > 1 ) then

            ! back -> forw
            source = back
            dest = forw
            tag = 100

            call mpi_sendrecv(ySliceForw_send(1,1), sendcount, &
                 & mpi_double_precision, dest, tag, &
                 & ySliceBack_recv(1,1), recvcount, mpi_double_precision, &
                 & source, mpi_any_tag, comm, sts_back, ierror)

            ! forw -> back
            source = forw
            dest = back
            tag = 100

            call mpi_sendrecv(ySliceBack_send(1,1), sendcount, &
                 & mpi_double_precision, dest, tag, &
                 & ySliceForw_recv(1,1), recvcount, mpi_double_precision, &
                 & source, mpi_any_tag, comm, sts_right, ierror)

            ! forward halos
            s(1:nx,ny+1,1:nz) = ySliceForw_recv(1:nx,1:nz)

            ! backward halos
            s(1:nx,0,1:nz) = ySliceBack_recv(1:nx,1:nz)

         else

           s(:,0,:) = s(:,ny,:)
           s(:,ny+1,:) = s(:,1,:)

         end if
       else
         stop "Poisson: unknown case xBoundary"
       end if

       if(verbose .and. master) print*,"horizontalHalos: &
            & x-horizontal halos copied."

       ! modified by Junhong Wei (20161106) *** finishing line ***


       !---------------------------------
       !         Loop over field
       !---------------------------------

       k_loop: do k = 1,nz
          j_loop: do j = 1,ny
             i_loop: do i = 1,nx

                ! ------------------ A(i+1,j,k) ------------------
                rhoEdge = 0.5*( var(i+1,j,k,1) + var(i,j,k,1) )
                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

                AR = dx2 * pStrat(k)**2/rhoEdge
                sR = s(i+1,j,k)

                ! ------------------- A(i-1,j,k) --------------------
                rhoEdge = 0.5*( var(i,j,k,1) + var(i-1,j,k,1) )
                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)                

                AL = dx2 * pStrat(k)**2 / rhoEdge
                sL = s(i-1,j,k)

                ! -------------------- A(i,j+1,k) ----------------------
                rhoEdge = 0.5*( var(i,j+1,k,1) + var(i,j,k,1) )
                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

                AF = dy2 * pStrat(k)**2 / rhoEdge
                sF = s(i,j+1,k)

                ! --------------------- A(i,j-1,k) -----------------------
                rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j-1,k,1) )
                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

                AB = dy2 * pStrat(k)**2 / rhoEdge
                sB = s(i,j-1,k)

                ! ---------------------- A(i,j,k+1) ------------------------
                if (k<nz) then
                   rhoEdge = 0.5* ( var(i,j,k+1,1) + var(i,j,k,1) )
                   if( fluctuationMode ) rhoEdge = rhoEdge + rhoStratTilde(k)

                   pStratU = 0.5* ( pStrat(k+1) + pStrat(k) )
                   AU = dz2 * pStratU**2 / rhoEdge
                   sU = s(i,j,k+1)
                else ! k = nz -> upwad boundary (solid wall)
                   ! A(i,j,nz+1) = 0
                   AU = 0.0
                   sU = 0.0
                end if

                ! ----------------------- A(i,j,k-1) ------------------------
                if (k>1) then 
                   rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j,k-1,1) )
                   if( fluctuationMode ) rhoEdge = rhoEdge + rhoStratTilde(k-1)
                   
                   pStratD = 0.5* ( pStrat(k) + pStrat(k-1) )
                   AD = dz2 * pStratD**2 / rhoEdge
                   sD = s(i,j,k-1)
                else ! k = 1 -> downward boundary (solid wall)
                   ! A(i,j,0) = 0
                   AD = 0.0
                   sD = 0.0
                end if

                ! ----------------------- A(i,j,k) --------------------------
                AC = - AR - AL - AF - AB - AU - AD
                sC = s(i,j,k)


                ! -------------------- apply Operator ---------------------
                Ls(i,j,k) = AL*sL + AR*sR + AF*sF + AB*sB + AU*sU + AD*sD + AC*sC


                ! ----------------- scale with thetaStrat --------------------
                if( pressureScaling ) then
                   Ls(i,j,k) = Ls(i,j,k) / PStrat(k)
                end if

             end do i_loop
          end do j_loop
       end do k_loop



    case( "Boussinesq" )
       !----------------------------------------
       !             Boussinesq model
       !----------------------------------------


       ! modified by Junhong Wei (20161106) *** starting line ***
!       
!       !----------------------------
!       !     Boundary conditions
!       !----------------------------
!
!       if( xBoundary == "periodic" ) then
!          ! perdiodic in x
!          s(0,:,:) = s(nx,:,:)
!          s(nx+1,:,:) = s(1,:,:)
!       end if
!
!       if( yBoundary == "periodic" ) then
!          ! periodic in y
!          s(:,0,:) = s(:,ny,:)
!          s(:,ny+1,:) = s(:,1,:)
!       end if
!
!       if( zBoundary == "periodic" ) then
!          ! periodic in z
!          s(:,:,0) = s(:,:,nz)
!          s(:,:,nz+1) = s(:,:,1)
!       end if

       !------------------------------
       !   set Halo cells: xSlice 
       !------------------------------

       if( xBoundary == "periodic" ) then

         ! slice size
         sendcount = ny*nz
         recvcount = sendcount

         ! read slice into contiguous array
         xSliceLeft_send (:,:) = s(1, 1:ny,1:nz)
         xSliceRight_send(:,:) = s(nx,1:ny,1:nz)

         if ( idim > 1 ) then

            ! left -> right
            source = left
            dest = right
            tag = 100

            call mpi_sendrecv(xSliceRight_send(1,1), sendcount, &
                 & mpi_double_precision, dest, tag, &
                 & xSliceLeft_recv(1,1), recvcount, mpi_double_precision, &
                 & source, mpi_any_tag, comm, sts_left, ierror)

            ! right -> left
            source = right
            dest = left
            tag = 100

            call mpi_sendrecv(xSliceLeft_send(1,1), sendcount, &
                 & mpi_double_precision, dest, tag, &
                 & xSliceRight_recv(1,1), recvcount, mpi_double_precision, &
                 & source, mpi_any_tag, comm, sts_right, ierror)

            ! right halos
            s(nx+1,1:ny,1:nz) = xSliceRight_recv(1:ny,1:nz)
            !          var(nx+i,1:ny,1:nz,iVar) = xSliceRight_recv(i,1:ny,1:nz)

            ! left halos
            s(0,1:ny,1:nz) = xSliceLeft_recv(1:ny,1:nz)
            !             var(-nbx+i,1:ny,1:nz,iVar) =
            !             xSliceLeft_recv(i,1:ny,1:nz)

         else

            s(0,:,:) = s(nx,:,:)
            s(nx+1,:,:) = s(1,:,:)

         end if
       else
         stop "Poisson: unknown case xBoundary"
       end if

       !------------------------------
       !   set Halo cells: ySlice 
       !------------------------------

       if( yBoundary == "periodic" ) then

         ! slice size
         sendcount = nx*nz
         recvcount = sendcount

         ! read slice into contiguous array
         ySliceBack_send(:,:) = s(1:nx, 1,1:nz)
         ySliceForw_send(:,:) = s(1:nx,ny,1:nz)

         if( jdim > 1 ) then

            ! back -> forw
            source = back
            dest = forw
            tag = 100

            call mpi_sendrecv(ySliceForw_send(1,1), sendcount, &
                 & mpi_double_precision, dest, tag, &
                 & ySliceBack_recv(1,1), recvcount, mpi_double_precision, &
                 & source, mpi_any_tag, comm, sts_back, ierror)

            ! forw -> back
            source = forw
            dest = back
            tag = 100

            call mpi_sendrecv(ySliceBack_send(1,1), sendcount, &
                 & mpi_double_precision, dest, tag, &
                 & ySliceForw_recv(1,1), recvcount, mpi_double_precision, &
                 & source, mpi_any_tag, comm, sts_right, ierror)

            ! forward halos
            s(1:nx,ny+1,1:nz) = ySliceForw_recv(1:nx,1:nz)

            ! backward halos
            s(1:nx,0,1:nz) = ySliceBack_recv(1:nx,1:nz)

         else

            s(:,0,:) = s(:,ny,:)
            s(:,ny+1,:) = s(:,1,:)

         endif
       else
         stop "Poisson: unknown case yBoundary"
       end if

       if( zBoundary == "periodic" ) then
          ! periodic in z
          s(:,:,0) = s(:,:,nz)
          s(:,:,nz+1) = s(:,:,1)
       end if

      ! modified by Junhong Wei (20161106) *** finishing line ***       



       !---------------------------------
       !         Loop over field
       !---------------------------------

       do k = 1,nz
          do j = 1,ny
             do i = 1,nx

                ! ------------------ A(i+1,j,k) ------------------
                AR = dx2
                sR = s(i+1,j,k)

                ! ------------------- A(i-1,j,k) --------------------
                AL = dx2
                sL = s(i-1,j,k)

                ! -------------------- A(i,j+1,k) ----------------------
                AF = dy2
                sF = s(i,j+1,k)

                ! --------------------- A(i,j-1,k) -----------------------
                AB = dy2
                sB = s(i,j-1,k)

                ! ---------------------- A(i,j,k+1) ------------------------
                AU = dz2 
                sU = s(i,j,k+1)

                ! ----------------------- A(i,j,k-1) ------------------------
                AD = dz2
                sD = s(i,j,k-1)

                ! ----------------------- A(i,j,k) --------------------------
                AC = - AR - AL - AF - AB - AU - AD
                sC = s(i,j,k)

                ! -------------------- apply Operator ---------------------
                Ls(i,j,k) = AL*sL + AR*sR + AF*sF + AB*sB + AU*sU + AD*sD + AC*sC


             end do
          end do
       end do

       !--------------------------------------
       !  Correction for solid wall boundary
       !--------------------------------------

       !---------------------------------
       !          zBoundary
       !---------------------------------

       if( zBoundary == "solid_wall" ) then

          do k = 1,nz,nz-1      ! k = 1 and nz
             do j = 1,ny
                do i = 1,nx

                   ! ------------------ A(i+1,j,k) ------------------
                   AR = dx2
                   sR = s(i+1,j,k)

                   ! ------------------- A(i-1,j,k) --------------------
                   AL = dx2
                   sL = s(i-1,j,k)

                   ! -------------------- A(i,j+1,k) ----------------------
                   AF = dy2
                   sF = s(i,j+1,k)

                   ! --------------------- A(i,j-1,k) -----------------------
                   AB = dy2
                   sB = s(i,j-1,k)

                   ! ---------------------- A(i,j,k+1) ------------------------
                   ! k = nz -> upwad boundary (solid wall)
                   ! A(i,j,nz+1) = 0
                   AU = 0.0
                   sU = 0.0

                   ! ----------------------- A(i,j,k-1) ------------------------
                   ! k = 1 -> downward boundary (solid wall)
                   ! A(i,j,0) = 0
                   AD = 0.0
                   sD = 0.0

                   ! ----------------------- A(i,j,k) --------------------------

                   AC = - AR - AL - AF - AB - AU - AD
                   sC = s(i,j,k)


                   ! -------------------- apply Operator ---------------------
                   Ls(i,j,k) = AL*sL + AR*sR + AF*sF + AB*sB + AU*sU + AD*sD + AC*sC

                end do
             end do
          end do

       end if ! zBoundary


    case default
       stop "linOpr: unknown case model"
    end select



  end subroutine linOpr


  !---------------------------------------------------------------------------


  subroutine calc_RHS( b,var,dt,onlyinfo )

    !----------------------------------------
    !   calculates the RHS of the 
    !   Poisson problem
    !----------------------------------------


    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         &intent(inout) :: var
    real, intent(in) :: dt
    real, dimension(1:nx,1:ny,1:nz), intent(out) :: b       ! RHS
    logical, intent(in) :: onlyinfo                         ! give info in div
    
    ! local vars
    real :: uR,uL, vF,vB, wU,wD
    real :: pStratU, pStratD

    integer :: i,j,k
    real :: div, divSum, divSumScaled

!   achatzb
    real :: bu,bv,bw,bl2loc,divL2_norm,divL2_norm_local
!   achatze

!   testb
!   integer :: i_maxb,j_maxb,k_maxb
!   teste

    ! extrapolation parameter for sol
    
    ! check L2-norm of divergence
    real :: divL2, divMax

    ! verbose
    logical, parameter :: giveInfo = .true.

    ! modified by Junhong Wei (20161106) *** starting line ***

    ! MPI stuff
    real :: divL2_local, divSum_local
    integer :: root

!   achatzb
    integer :: i0,j0
!   achatze

!xxxx
    if( giveInfo .and. master ) then
       print*,""
       print*,"----------------------------------------------"
       print*,"   calc_RHS: computing RHS... "
       print*,"----------------------------------------------"
       print*,""
    end if

    ! modified by Junhong Wei (20161106) *** finishing line ***

    ! ----------------------------------
    !          Poisson Problem
    ! ----------------------------------

    !    if(verbose) print*,"update.f90/poissonSolver: &   ! modified by Junhong Wei (20161106)
    if( master .and. verbose ) print*,"update.f90/poissonSolver: &
         & Setting up Poisson problem."   ! modified by Junhong Wei (20161106)

!   achatzb
    i0=is+nbx-1
    j0=js+nby-1
!   achatze
    
    
    !--------------------------------------------------
    !    setup b = Ma^2 * P * u^*  (right hand side)
    !--------------------------------------------------

    divSum = 0.0
    divSumScaled = 0.0
    divL2  = 0.0
    divMax = 0.0

    divSum_local = 0.0   ! modified by Junhong Wei (20161107)
    divL2_local = 0.0   ! modified by Junhong Wei (20161107)

!   achatzb
    divL2_norm = 0.0
    divL2_norm_local = 0.0
!   achatze

    select case( model ) 

    case( "pseudo_incompressible" ) 

       do k = 1,nz
          do j  = 1,ny
             do i = 1,nx

                uR = var(i,j,k,2); uL = var(i-1,j,k,2)
                vF = var(i,j,k,3); vB = var(i,j-1,k,3)
                wU = var(i,j,k,4); wD = var(i,j,k-1,4)

!xxxx 10.5.2011
!                pStratU = 0.5*(pStrat(k+1) + pStrat(k  ))
!                pStratD = 0.5*(pStrat(k  ) + pStrat(k-1))
                
                PstratU = PstratTilde(k)
                PstratD = PstratTilde(k-1)
                
!xxxx end
                
!               achatzb
!               introduce a reference norm for the RHS
!               b(i,j,k) = Pstrat(k) * ( (uR-uL)/dx + (vF-vB)/dy ) &
!                    & + (PstratU*wU - PstratD*wD)/dz
                bu = Pstrat(k) * (uR-uL)/dx
                bv = Pstrat(k) * (vF-vB)/dy 
                bw = (PstratU*wU - PstratD*wD)/dz
                b(i,j,k) = bu + bv + bw
                bl2loc = bu**2 + bv**2 + bw**2
!               achatze

!               testb
!               if((i == 1).and.(j == 1).and.(k == 310)) then
!                  print*,"at (i,j,k) = (1,1,310):"
!                  print*,"bu = ",bu
!                  print*,"bv = ",bv
!                  print*,"bw = ",bw
!                  print*,"b = ",b(i,j,k)

!                  print*,"uL,uR = ",uL,uR

!                  print*,"u(0,1,310),u(8,1,310) = ", &
!                        & var(0,1,310,2),var(nx,1,310,2)

!                  print*,"u(1,1,310),u(9,1,310) = ", &
!                        & var(1,1,310,2),var(nx+1,1,310,2)
!               end if
!               teste

!               achatzb
                if(topography) then
                   if(topography_mask(i0+i,j0+j,k)) then
                      b(i,j,k) = 0.0
                      bl2loc = 0.0
                   end if
                end if
!               achatze

                ! L2-norm of the divergence div(Pu)
                ! divL2 = divL2 + b(i,j,k)**2

                divL2_local = divL2_local + b(i,j,k)**2

!               achatzb
                divL2_norm_local = divL2_norm_local + bl2loc
!               achatze

                ! max norm of div(Pu)
                if( abs(b(i,j,k)) > divMax ) then
                   divMax = abs(b(i,j,k))
!                  testb
!                  i_maxb = i
!                  j_maxb = j
!                  k_maxb = k
!                  teste
                end if
                
                ! scale RHS with Ma^2 * kappa
                b(i,j,k) = b(i,j,k) * Ma**2 * kappa
                
                ! check sum for solvability criterion
                ! divSum = divSum + b(i,j,k)
                divSum_local = divSum_local + b(i,j,k)
                
                ! Skalierung mit thetaStrat
                if ( pressureScaling ) then
                   b(i,j,k) = b(i,j,k) / PStrat(k)
                   divSumScaled = divSumScaled + b(i,j,k)
                end if

             end do
          end do
       end do

!      testb
!      print*,"divMax = ",divMax," reached at ",i_maxb,j_maxb,k_maxb
!      teste
       
!       ! scale by number of cells   ! modified by Junhong Wei (20161107)
!       divL2 = sqrt(divL2/nx/ny/nz)   ! modified by Junhong Wei (20161107)

       ! modified by Junhong Wei (20161107)   *** starting line ***

       !MPI: sum divSum_local over all procs
       root = 0
       call mpi_reduce(divSum_local, divSum, 1, &
            & mpi_double_precision, mpi_sum, root, comm, ierror)
       
       call mpi_bcast(divSum, 1, &
            & mpi_double_precision, root, comm, ierror)

       !MPI: sum divL2_local over all procs
       root = 0
       call mpi_reduce(divL2_local, divL2, 1, &
            & mpi_double_precision, mpi_sum, root, comm, ierror)
       
       call mpi_bcast(divL2, 1, &
            & mpi_double_precision, root, comm, ierror)

!      achatzb
       !MPI: sum divL2_norm_local over all procs
       root = 0
       call mpi_reduce(divL2_norm_local, divL2_norm, 1, &
            & mpi_double_precision, mpi_sum, root, comm, ierror)
       
       call mpi_bcast(divL2_norm, 1, &
            & mpi_double_precision, root, comm, ierror)
!      achatze

       ! scale div
       divL2_local = sqrt(divL2_local/nx/ny/nz)
       divL2 = sqrt(divL2/sizeX/sizeY/sizeZ)

!      achatzb
       divL2_norm_local = sqrt(divL2_norm_local/nx/ny/nz)
       divL2_norm = sqrt(divL2_norm/sizeX/sizeY/sizeZ)

       tolref=divL2/divL2_norm
       if(master) print*,"tolref = ",tolref

!      root=0
!      call mpi_bcast(tolref, 1, mpi_double_precision, root, comm, ierror)
!      call mpi_barrier(comm,ierror)
!      achatze

       ! modified by Junhong Wei (20161107)   *** finishing line ***


    case( "Boussinesq" )

       do k = 1,nz
          do j  = 1,ny
             do i = 1,nx

                uR = var(i,j,k,2); uL = var(i-1,j,k,2)
                vF = var(i,j,k,3); vB = var(i,j-1,k,3)
                wU = var(i,j,k,4); wD = var(i,j,k-1,4)

!               achatzb
!               introduce a reference norm for the RHS
!               div = (uR-uL)/dx + (vF-vB)/dy + (wU-wD)/dz
                bu = (uR-uL)/dx
                bv = (vF-vB)/dy
                bw = (wU-wD)/dz
                div = bu + bv + bw
                bl2loc = bu**2 + bv**2 + bw**2
!               achatze

                b(i,j,k) =  Ma**2 * kappa / theta00 * div

!               achatzb
                if(topography) then
                   if(topography_mask(i0+i,j0+j,k)) then
                      b(i,j,k) = 0.0
                      bl2loc = 0.0
                   end if
                end if
!               achatze

                ! check sum for solvability criterion (shoud be zero)
                ! divSum = divSum + b(i,j,k)
                divSum_local = divSum_local + b(i,j,k)
                
                ! divL2 = divL2 + div**2
                divL2_local = divL2_local + div**2

!               achatzb
                divL2_norm_local = divL2_norm_local + bl2loc
!               achatze

                if( abs(div) > divMax ) divMax = abs(div)

             end do
          end do
       end do


!       ! scale by number of cells   ! modified by Junhong Wei (20161107)
!       divL2 = sqrt(divL2/nx/ny/nz)   ! modified by Junhong Wei (20161107)

       ! modified by Junhong Wei (20161107)   *** starting line ***

              !MPI: sum divSum_local over all procs
       root = 0
       call mpi_reduce(divSum_local, divSum, 1, &
            & mpi_double_precision, mpi_sum, root, comm, ierror)

       call mpi_bcast(divSum, 1, &
            & mpi_double_precision, root, comm, ierror)

       !MPI: sum divL2_local over all procs
       root = 0
       call mpi_reduce(divL2_local, divL2, 1, &
            & mpi_double_precision, mpi_sum, root, comm, ierror)

       call mpi_bcast(divL2, 1, &
            & mpi_double_precision, root, comm, ierror)

!      achatzb
       !MPI: sum divL2_norm_local over all procs
       root = 0
       call mpi_reduce(divL2_norm_local, divL2_norm, 1, &
            & mpi_double_precision, mpi_sum, root, comm, ierror)
       
       call mpi_bcast(divL2_norm, 1, &
            & mpi_double_precision, root, comm, ierror)
!      achatze

       ! scale div
       divL2_local = sqrt(divL2_local/nx/ny/nz)
       divL2 = sqrt(divL2/sizeX/sizeY/sizeZ)

       ! scale by number of cells 
       divL2_local = sqrt(divL2_local/nx/ny/nz)
       divL2 = sqrt(divL2/sizeX/sizeY/sizeZ)

!      achatzb
       divL2_norm_local = sqrt(divL2_norm_local/nx/ny/nz)
       divL2_norm = sqrt(divL2_norm/sizeX/sizeY/sizeZ)

       tolref=divL2/divL2_norm
       if(master) print*,"tolref = ",tolref

!      root=0
!      call mpi_bcast(tolref, 1, mpi_double_precision, root, comm, ierror)
!      call mpi_barrier(comm,ierror)
!      achatze

       ! modified by Junhong Wei (20161107)   *** finishing line ***
     

    case default
       stop "poissonSolver: unknown case model."
    end select


    !-------------------------------
    !     Display info on screen
    !-------------------------------

    !    if( onlyinfo ) then   ! modified by Junhong Wei (20161107)
        if( master .and. onlyinfo ) then   ! modified by Junhong Wei (20161107)
       
       ! Information on divergence

       print*,""
       print*," Poisson Solver ", trim(poissonSolverType), ": Final state"
       print*,""

       select case ( model ) 
       case( "Boussinesq" ) 
          
          write(*,fmt="(a25,es17.6)") "L2(div(u)) [1/s] = ", &
               & divL2*uRef/lRef

          write(*,fmt="(a25,es17.6)") "max(div(u)) [1/s] = ", &
               & divMax*uRef/lRef

!         achatzb
          write(*,fmt="(a25,es17.6)") "rms terms (div(u)) [1/s] = ", &
               & divL2_norm*uRef/lRef

          write(*,fmt="(a25,es17.6)") "normalized L2(div(u)) = ", &
               & divL2/divL2_norm
!         achatze
          
       case( "pseudo_incompressible" )

          write(*,fmt="(a25,es17.6)") "L2(div(Pu)) [Pa/s] = ", &
               & divL2*rhoRef*thetaRef*uRef/lRef

          write(*,fmt="(a25,es17.6)") "max(div(Pu)) [Pa/s] = ", &
               & divMax*rhoRef*thetaRef*uRef/lRef
          
!         achatzb
          write(*,fmt="(a25,es17.6)") "rms terms (div(Pu)) [Pa/s] = ", &
               & divL2_norm*rhoRef*thetaRef*uRef/lRef

          write(*,fmt="(a25,es17.6)") "normalized L2(div(Pu)) = ", &
               & divL2/divL2_norm
!         achatze
          
       case default
          stop "calc_RHS: unkown case model"
       end select

       print*,""
       print*,"-------------------------------------------------"
       print*,""


          
    else

       !       if( giveInfo ) then   ! modified by Junhong Wei (20161107)
              if( master .and. giveInfo ) then   ! modified by Junhong Wei (20161107)
          print*,""
          print*," Poisson Solver ", trim(poissonSolverType), ": Initial state"
          print*,""
          write(*,fmt="(a25,es17.6)") "Sum over RHS = ", divSum
          write(*,fmt="(a25,es17.6)") "Sum over scaled RHS  = ", divSumScaled

          ! Information on divergence
          select case ( model ) 
          case( "Boussinesq" ) 

             write(*,fmt="(a25,es17.6)") "L2(div(u)) [1/s] = ", &
                  & divL2*uRef/lRef

             write(*,fmt="(a25,es17.6)") "max(div(u)) [1/s] = ", &
                  & divMax*uRef/lRef

!            achatzb
             write(*,fmt="(a25,es17.6)") "rms terms (div(u)) [1/s] = ", &
                  & divL2_norm*uRef/lRef

             write(*,fmt="(a25,es17.6)") "normalized L2(div(u)) = ", &
                  & divL2/divL2_norm
!            achatze
          
          case( "pseudo_incompressible" )

             write(*,fmt="(a25,es17.6)") "L2(div(Pu)) [Pa/s] = ", &
                  & divL2*rhoRef*thetaRef*uRef/lRef

             write(*,fmt="(a25,es17.6)") "max(div(Pu)) [Pa/s] = ", &
                  & divMax*rhoRef*thetaRef*uRef/lRef

!            achatzb
             write(*,fmt="(a25,es17.6)") "rms terms (div(Pu)) [Pa/s] = ", &
                  & divL2_norm*rhoRef*thetaRef*uRef/lRef
   
             write(*,fmt="(a25,es17.6)") "normalized L2(div(Pu)) = ", &
                  & divL2/divL2_norm
!            achatze
          
          case default
             stop "calc_RHS: unkown case model"
          end select

       end if

    end if

  end subroutine calc_RHS


  !---------------------------------------------------------------------------


  subroutine poissonSolver( b,var,dt,errFlagBicg,nIter,m,opt )
    ! -------------------------------------------------
    ! solves the Poisson problem with 
    ! application of linear operator L
    ! -------------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         &intent(inout) :: var
    real, intent(in) :: dt
    logical, intent(out) :: errFlagBicg
    integer, intent(out) :: nIter
    integer, intent(in) :: m
    real, dimension(1:nx,1:ny,1:nz),intent(in) :: b        ! RHS
    character(len=*), intent(in)    :: opt



    ! local vars
    real, dimension(1:nx,1:ny,1:nz) :: sol      ! solution of Poisson problem
    real :: res
    
    real :: dtInv

    ! verbose
    logical, parameter :: giveInfo = .true.



    ! Init
    if (dt == 0.0) stop "poissonSolver: dt = 0.0. Stopping."
    dtInv = 1.0/dt


    !--------------------------------
    !     Linear equation solver
    !--------------------------------

    sol = 0.0

    select case( poissonSolverType ) 

    case( "bicgstab" )

       call bicgstab(var, b, dt, sol, res, nIter, errFlagBicg, opt)

    case( "gcr" ) 

       call gcr(var, b, dt, sol, res, nIter, errFlagBicg)

    case( "adi" ) 

       call adi(var, b, dt, sol, res, nIter, errFlagBicg)

    case( "hypre" )                                           ! modified by Junhong Wei (2016/07/07)

       call hypre(var, b, dt, sol, res, nIter, errFlagBicg, opt)   ! modified by Junhong Wei (2016/07/07)

    case default
       stop "Unknown PoissonSolver. Stop"
    end select


    dp(1:nx,1:ny,1:nz) = dtInv * sol      ! pass solution to pressure corrector


  end subroutine poissonSolver


  !---------------------------------------------------------------------------


  subroutine linOprXYZ( var, q, Lq, direction )
    ! --------------------------------------
    !   Linear Operator in Poisson problem
    !   Functions as A*x
    ! --------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         &intent(in) :: var
    real, dimension(1:nx,1:ny,1:nz), intent(out) :: Lq
    real, dimension(0:nx+1,0:ny+1,0:nz+1), intent(inout)  :: q ! with ghost cells
    character(len=1), intent(in) :: direction 


    ! local variables
    integer :: i,j,k
    real :: pStratU, pStratD, rhoEdge
    real :: AL,AR, AB,AF, AD,AU, AC
    real :: qL,qR, qB,qF, qD,qU, qC

    real :: dx2, dy2, dz2

!   achatzb
    if(topography) stop 'linOprXYZ not ready for topography!'
!   achatze

    ! auxiliary variables
    dx2 = 1.0/dx**2
    dy2 = 1.0/dy**2
    dz2 = 1.0/dz**2

    ! ----------------------------
    !     boundary conditions
    ! ----------------------------

    ! perdiodic in x
    q(0,:,:) = q(nx,:,:)
    q(nx+1,:,:) = q(1,:,:)

    ! periodic in y
    q(:,0,:) = q(:,ny,:)
    q(:,ny+1,:) = q(:,1,:)


    select case( direction ) 

       !--------------------------
       !       operator Lx
       !--------------------------

    case( "x" )

       k_loop_x: do k = 1,nz
          j_loop_x: do j = 1,ny
             i_loop_x: do i = 1,nx

                ! ------------------ A(i+1,j,k) ------------------
                rhoEdge = 0.5*( var(i+1,j,k,1) + var(i,j,k,1) )
                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

                AR = dx2 * pStrat(k)**2/rhoEdge
                qR = q(i+1,j,k)

                ! ------------------- A(i-1,j,k) --------------------
                rhoEdge = 0.5*( var(i,j,k,1) + var(i-1,j,k,1) )
                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)
                
                AL = dx2 * pStrat(k)**2 / rhoEdge
                qL = q(i-1,j,k)

                ! ----------------------- A(i,j,k) --------------------------
                AC = - AL - AR 
                qC = q(i,j,k)

                ! -------------------- apply Operator ---------------------
                Lq(i,j,k) = AL*qL + AC*qC + AR*qR 

                ! ---------------------- scale with PStrat --------------------
                if( pressureScaling ) then
                   Lq(i,j,k) = Lq(i,j,k) / Pstrat(k)
                end if

             end do i_loop_x
          end do j_loop_x
       end do k_loop_x



       !--------------------------
       !       operator Ly
       !--------------------------

    case( "y" )

       k_loop_y: do k = 1,nz
          j_loop_y: do j = 1,ny
             i_loop_y: do i = 1,nx

                ! -------------------- A(i,j+1,k) ----------------------
                rhoEdge = 0.5*( var(i,j+1,k,1) + var(i,j,k,1) )
                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

                AF = dy2 * pStrat(k)**2 / rhoEdge
                qF = q(i,j+1,k)

                ! --------------------- A(i,j-1,k) -----------------------
                rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j-1,k,1) )
                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)
                
                AB = dy2 * pStrat(k)**2 / rhoEdge
                qB = q(i,j-1,k)

                ! ----------------------- A(i,j,k) --------------------------
                AC = - AF - AB 
                qC = q(i,j,k)

                ! -------------------- apply Operator ---------------------
                Lq(i,j,k) = AF*qF + AB*qB + AC*qC

                ! ---------------------- scale with PStrat --------------------
                if( pressureScaling ) then
                   Lq(i,j,k) = Lq(i,j,k) / Pstrat(k)
                end if

             end do i_loop_y
          end do j_loop_y
       end do k_loop_y


       !--------------------------
       !       operator Lz
       !--------------------------

    case( "z" )

       k_loop_z: do k = 1,nz
          j_loop_z: do j = 1,ny
             i_loop_z: do i = 1,nx

                ! ---------------------- A(i,j,k+1) ------------------------
                if (k<nz) then
                   rhoEdge = 0.5* ( var(i,j,k+1,1) + var(i,j,k,1) )
                   if( fluctuationMode ) rhoEdge = rhoEdge + rhoStratTilde(k)

                   pStratU = 0.5* ( pStrat(k+1) + pStrat(k) )
                   AU = dz2 * pStratU**2 / rhoEdge
                   qU = q(i,j,k+1)
                else ! k = nz -> upwad boundary (solid wall)
                   ! A(i,j,nz+1) = 0
                   AU = 0.0
                   qU = 0.0
                end if

                ! ----------------------- A(i,j,k-1) ------------------------
                if (k>1) then 
                   rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j,k-1,1) )
                   if( fluctuationMode ) rhoEdge = rhoEdge + rhoStratTilde(k)

                   pStratD = 0.5* ( pStrat(k) + pStrat(k-1) )
                   AD = dz2 * pStratD**2 / rhoEdge
                   qD = q(i,j,k-1)
                else ! k = 1 -> downward boundary (solid wall)
                   ! A(i,j,0) = 0
                   AD = 0.0
                   qD = 0.0
                end if

                ! ----------------------- A(i,j,k) --------------------------
                AC = - AU - AD
                qC = q(i,j,k)

                ! -------------------- apply Operator ---------------------
                Lq(i,j,k) = AU*qU + AD*qD + AC*qC

                ! ---------------------- scale with PStrat --------------------
                if( pressureScaling ) then
                   Lq(i,j,k) = Lq(i,j,k) / Pstrat(k)
                end if

             end do i_loop_z
          end do j_loop_z
       end do k_loop_z


    case default
       stop "linOprXYZ: unknown direction"
    end select


  end subroutine linOprXYZ


  !---------------------------------------------------------------------------


  subroutine adi(var, b, dt, sol, res, nIter, errFlag)
    ! --------------------------------------
    !   ADI scheme for 3D
    !   cf. Roache, Computational Fluid Dynamics
    !---------------------------------------

    ! In/Out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var
    real, dimension(1:nx,1:ny,1:nz), intent(in) :: b        ! RHS 
    real, intent(in) :: dt
    real, dimension(1:nx,1:ny,1:nz), intent(inout) :: sol
    real, intent(out) :: res                     ! residual
    integer, intent(out) :: nIter
    logical, intent(out) :: errFlag

    ! Local variables
    integer :: i,j,k, allocstat, iter, maxIt
    real :: pStratU, pStratD

    ! Local field
    real, dimension(1:nx,1:ny,1:nz) :: sol0


    ! allocatable fields
    real, dimension(:,:,:), allocatable :: r,Lx3D
    real, dimension(:,:,:), allocatable :: Lx,Ly,Lz
    real, dimension(:,:,:), allocatable :: q,qn,q01,q02
    real, dimension(:),     allocatable :: left_a,center,right_a,rhs

    real :: rhoEdge, AL,AR,AB,AF,AD,AU,AC
    real :: dx2, dy2, dz2

    ! verbose
    logical, parameter :: giveInfo = .true.

    ! test variables
    real :: res0, res01, res02, res03

!   achatzb
    if(topography) stop 'adi not ready for topography!'
!   achatze

    ! auxiliary variables
    dx2 = 1.0/dx**2
    dy2 = 1.0/dy**2
    dz2 = 1.0/dz**2

    if( pressureScaling .and. poissonSolverType=="adi" ) &
         & stop "ADI as poissonSolver: no pressure scaling allowed. Stop."

    ! Set parameter for ADI as Poisson Solver
    if( poissonSolverType == "adi" ) then
       maxIt = maxIterPoisson
       tol = tolPoisson
    else
       maxIt = maxIterADI       ! as preconditioner only 1 or 2 iterations
       tol = 0.0
    end if

    ! Save input solution
    sol0 = sol


    !---- Allocate fields without ghost cells: r, Lx3D
    allocate(r(1:nx,1:ny,1:nz), stat=allocstat)
    if(allocstat/=0) stop "adi:alloc failed"
    allocate(Lx3D(1:nx,1:ny,1:nz), stat=allocstat)
    if(allocstat/=0) stop "adi:alloc failed"

    !--- Allocate Lx, Ly, Lz
    allocate(Lx(1:nx,1:ny,1:nz), stat=allocstat)
    if(allocstat/=0) stop "adi:alloc failed"
    allocate(Ly(1:nx,1:ny,1:nz), stat=allocstat)
    if(allocstat/=0) stop "adi:alloc failed"
    allocate(Lz(1:nx,1:ny,1:nz), stat=allocstat)
    if(allocstat/=0) stop "adi:alloc failed"


    !---- Allocate fields with ghost cells: q, qn, q01, q02
    allocate(q(0:nx+1,0:ny+1,0:nz+1), stat=allocstat)
    if(allocstat/=0) stop "adi:alloc failed"
    allocate(qn(0:nx+1,0:ny+1,0:nz+1), stat=allocstat)
    if(allocstat/=0) stop "adi:alloc failed"
    allocate(q01(0:nx+1,0:ny+1,0:nz+1), stat=allocstat)
    if(allocstat/=0) stop "adi:alloc failed"
    allocate(q02(0:nx+1,0:ny+1,0:nz+1), stat=allocstat)
    if(allocstat/=0) stop "adi:alloc failed"



    ! Error flag
    errFlag = .false.    

    ! Init
    call linOpr( var, sol, Lx3D ) ! -> r = L(x) - b         NOT b - L(x) !!!
    r = Lx3D - b


    !---- Check initial residual
    res = dt*maxval(abs(r))             ! abort criterium by Smolarkiewicz
    ! res = maxval(abs(r))  
    res0 = res
    print*,"  ADI "
    if (giveInfo) write(*,fmt="(a25,es17.6)") " Init residual: res0 = ", res
    if (giveInfo) write(*,fmt="(a25,es17.6)") " tol = ", tol
    ! comment by Smolarkiewicz: run at least one iteration
!!$    if (res <= tol) then
!!$       if(giveInfo) print*," ==> no iteration needed."
!!$       nIter = 0
!!$       return
!!$    end if


    ! Loop
    iteration: do iter = 1,maxIt

       ! Init qn with sol
       qn(1:nx,1:ny,1:nz) = sol

       !
       !---------------------------------
       !    sweep in x direction -> q01
       !---------------------------------
       !

       !---- prepare explicit terms Lx, Ly, Lz for rhs

       call linOprXYZ( var, qn, Lx, "x" )
       call linOprXYZ( var, qn, Ly, "y" )
       call linOprXYZ( var, qn, Lz, "z" )


       !---- set up matrix diagonals for Thomas algorithm

       ! Allocate variables for x sweeps
       allocate(left_a(1:nx),stat=allocstat)
       if(allocstat/=0) stop "adi:alloc failed"
       allocate(center(1:nx),stat=allocstat) 
       if(allocstat/=0) stop "adi:alloc failed"
       allocate(right_a(1:nx), stat=allocstat)
       if(allocstat/=0) stop "adi:alloc failed"
       allocate(rhs(1:nx), stat=allocstat)
       if(allocstat/=0) stop "adi:alloc failed"

       k_loop_x: do k = 1,nz
          j_loop_x: do j = 1,ny

             !--- matrix diagonals
             i_loop_x: do i = 1,nx

                ! A(i+1,j,k)
                rhoEdge = 0.5*( var(i+1,j,k,1) + var(i,j,k,1) )
                AR = dx2 * pStrat(k)**2 / rhoEdge
                if( pressureScaling ) then
                   stop "ADI: Please check correct implementation of &
                        & pressureScaling for ADI."
                   AR = AR / Pstrat(k) 
                end if
                
                ! A(i-1,j,k)
                rhoEdge = 0.5*( var(i,j,k,1) + var(i-1,j,k,1) )
                AL = dx2 * pStrat(k)**2 / rhoEdge
                if( pressureScaling ) AL = AL / Pstrat(k)

                ! A(i,j,k)
                AC = - AL - AR

                left_a(i) = -0.5 * AL
                center(i) = 1.0/dtau - 0.5*AC
                right_a(i) = -0.5 * AR

             end do i_loop_x

             !---- Right hand side rhs

             rhs =  qn(1:nx,j,k) / dtau + &
                  & 0.5 * Lx(1:nx,j,k) + &
                  & Ly(1:nx,j,k) + &
                  & Lz(1:nx,j,k) - &
                  & b(1:nx,j,k)

             ! test: switch off BC
             ! left_a / right_a boundary terms 
             rhs(1) = rhs(1) - left_a(1)*qn(nx,j,k)
             rhs(nx) = rhs(nx) - right_a(nx)*qn(1,j,k)
             ! end test

             !---- solve tridiagonal system
             call thomas( left_a, center, right_a, q01(1:nx,j,k), rhs )

          end do j_loop_x
       end do k_loop_x


       ! test quality of q01
!!$
!!$    call linOpr( var, q01(1:nx,1:ny,1:nz), Lx3D ) ! -> r = L(x) - b
!!$        NOT b - L(x) !!!
!!$    r = Lx3D - b
!!$    res01 = dt*maxval(abs(r))             ! abort criterium by Smolarkiewicz
!!$    if (giveInfo) write(*,fmt="(a25,es17.6)") " res(q01) = ", res01
!!$
       ! end test       

       !---- Deallocate variables for x sweep
       deallocate(left_a,stat=allocstat)
       if(allocstat/=0) stop "adi:dealloc failed"
       deallocate(center,stat=allocstat)
       if(allocstat/=0) stop "adi:dealloc failed"
       deallocate(right_a,stat=allocstat)
       if(allocstat/=0) stop "adi:dealloc failed"
       deallocate(rhs,stat=allocstat)
       if(allocstat/=0) stop "adi:dealloc failed"

       !
       !---------------------------------
       !    sweep in y direction -> q02 
       !---------------------------------
       !

       !---- prepare explicit terms Lx, Ly, Lz for rhs
       q = 0.5*(q01+qn)
       call linOprXYZ( var, q, Lx, "x" )
       ! Ly, Lz can be re-used

       !---- set up matrix diagonals for Thomas algorithm

       ! Allocate variables for x sweeps
       allocate(left_a(1:ny), stat=allocstat)
       if(allocstat/=0) stop "adi:alloc failed"
       allocate(center(1:ny), stat=allocstat)
       if(allocstat/=0) stop "adi:alloc failed"
       allocate(right_a(1:ny), stat=allocstat)
       if(allocstat/=0) stop "adi:alloc failed"
       allocate(rhs(1:ny), stat=allocstat)
       if(allocstat/=0) stop "adi:alloc failed"

       k_loop_y: do k = 1,nz
          i_loop_y: do i = 1,nx

             !--- matrix diagonals
             j_loop_y: do j = 1,ny

                ! A(i,j+1,k)
                rhoEdge = 0.5*( var(i,j+1,k,1) + var(i,j,k,1) )
                AF = dy2 * pStrat(k)**2 / rhoEdge
                if( pressureScaling ) AF = AF / Pstrat(k)


                ! A(i,j-1,k)
                rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j-1,k,1) )
                AB = dy2 * pStrat(k)**2 / rhoEdge
                if( pressureScaling ) AB = AB / Pstrat(k)

                ! A(i,j,k)
                AC = -AB - AF

                left_a(j) = -0.5*AB
                center(j) = 1.0/dtau - 0.5*AC
                right_a(j) = -0.5*AF

             end do j_loop_y

             !---- Right hand side rhs

             rhs =  qn(i,1:ny,k) / dtau + &
                  & Lx(i,1:ny,k) + &
                  & 0.5*Ly(i,1:ny,k) + &
                  & Lz(i,1:ny,k) - &
                  & b(i,1:ny,k)

             ! forward / backward boundary terms 
             rhs(1) = rhs(1) - left_a(1)*q01(i,ny,k)
             rhs(ny) = rhs(ny) - right_a(ny)*q01(i,1,k)

             !---- solve tridiagonal system
             call thomas( left_a, center, right_a, q02(i,1:ny,k), rhs )


          end do i_loop_y
       end do k_loop_y

!!$! test quality of q02
!!$
!!$    call linOpr( var, q02(1:nx,1:ny,1:nz), Lx3D ) ! -> r = L(x) - b                    NOT b - L(x) !!!
!!$    r = Lx3D - b
!!$    res02 = dt*maxval(abs(r))             ! abort criterium by Smolarkiewicz
!!$    if (giveInfo) write(*,fmt="(a25,es17.6)") " res(q02) = ", res02
!!$
!!$! end test       


       !---- Deallocate variables for x sweep
       deallocate(left_a,stat=allocstat)
       if(allocstat/=0) stop "adi:dealloc failed"
       deallocate(center,stat=allocstat)
       if(allocstat/=0) stop "adi:dealloc failed"
       deallocate(right_a,stat=allocstat)
       if(allocstat/=0) stop "adi:dealloc failed"
       deallocate(rhs,stat=allocstat)
       if(allocstat/=0) stop "adi:dealloc failed"



       !
       !---------------------------------
       !    sweep in z direction -> sol
       !---------------------------------
       !

       !---- prepare explicit terms Lx, Ly, Lz for rhs
       q =  0.5*(q02+qn)
       call linOprXYZ( var, q, Ly, "y" )
       ! Lx, Lz can be re-used

       !---- set up matrix diagonals for Thomas algorithm

       ! Allocate variables for x sweeps
       allocate(left_a(1:nz), stat=allocstat)
       if(allocstat/=0) stop "adi:alloc failed"
       allocate(center(1:nz), stat=allocstat)
       if(allocstat/=0) stop "adi:alloc failed"
       allocate(right_a(1:nz), stat=allocstat)
       if(allocstat/=0) stop "adi:alloc failed"
       allocate(rhs(1:nz), stat=allocstat)
       if(allocstat/=0) stop "adi:alloc failed"

       j_loop_z: do j = 1,ny
          i_loop_z: do i = 1,nx

             !--- matrix diagonals
             k_loop_z: do k = 1,nz

                !  A(i,j,k+1) 
                if (k<nz) then

                   rhoEdge = 0.5* ( var(i,j,k+1,1) + var(i,j,k,1) )
                   AU = dz2 * pStratU**2 / rhoEdge
                   if( pressureScaling ) AU = AU / PstratTilde(k)

                else ! k = nz -> upwad boundary (solid wall) -> A(i,j,nz+1) = 0
                   AU = 0.0
                end if

                !  A(i,j,k-1) 
                if (k>1) then 

                   rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j,k-1,1) )
                   AD = dz2 * pStratD**2 / rhoEdge
                   if( pressureScaling ) AD = AD / PstratTilde(k-1)

                else ! k = 1 -> downward boundary (solid wall) -> A(i,j,0) = 0
                   AD = 0.0
                end if

                !  A(i,j,k) 
                AC = - AU - AD

                left_a(k) = -0.5 * AD
                center(k) = 1.0/dtau - 0.5*AC
                right_a(k) = -0.5 * AU

             end do k_loop_z

             !---- Right hand side rhs


             rhs =  qn(i,j,1:nz) / dtau + &
                  & Lx(i,j,1:nz) + &
                  & Ly(i,j,1:nz) + &
                  & 0.5*Lz(i,j,1:nz) - &
                  & b(i,j,1:nz)

             !---- solve tridiagonal system
             call thomas( left_a, center, right_a, sol(i,j,1:nz), rhs )

          end do i_loop_z
       end do j_loop_z

       ! test quality of sol
!!$
!!$    call linOpr( var, sol(1:nx,1:ny,1:nz), Lx3D ) ! -> r = L(x) - b
!!$  NOT b - L(x) !!!
!!$    r = Lx3D - b
!!$    res03 = dt*maxval(abs(r))             ! abort criterium by Smolarkiewicz
!!$    if (giveInfo) write(*,fmt="(a25,es17.6)") " res(sol) = ", res03
!!$    print*,"---------------------------------------------------------"
!!$
       ! end test       


       !---- Deallocate variables for x sweep
       deallocate(left_a,stat=allocstat)
       if(allocstat/=0) stop "adi:dealloc failed"
       deallocate(center,stat=allocstat)
       if(allocstat/=0) stop "adi:dealloc failed"
       deallocate(right_a,stat=allocstat)
       if(allocstat/=0) stop "adi:dealloc failed"
       deallocate(rhs,stat=allocstat)
       if(allocstat/=0) stop "adi:dealloc failed"


       !
       !---------------------------------
       !    Evaluate result
       !---------------------------------
       !

       call linOpr( var, sol, Lx3D ) ! -> r = L(x) - b, NOT b - L(x) !!!
       r = Lx3D - b
       res = dt*maxval(abs(r))    ! abort criterion by Smolarkiewiecz
       !  res = maxval(abs(r))            

       !---- Abort criterium 
       if( poissonSolverType=="adi" ) then ! Abort only for ADI as Poisson solver 
          if (res <= tol) then
             if (giveInfo)  then 
                write(*,fmt="(a25,i25)") " Nb.of iterations: iter = ", iter
                write(*,fmt="(a25,es17.6)") " Final residual: res = ", res
                print*,"--------------------------------------------------"
                print*,""
             end if
             nIter = iter
             return
          end if
       end if


       ! Treatment for ADI as preconditioner
       if( preconditioner == "adi") then
          write(*,fmt="(a25,es17.6)") " ADI residual = ", res
          write(*,fmt="(a25,f25.18)") " -> residual reduction = ", res/res0
          ! Return old solution if no residual reduction
          if( res/res0 >= 1.0 ) then
             print*,"==> no residual reduction!!!"
             !             dtau = 0.9*dtau
             !             write(*,fmt="(a25,f25.18)") " reduce dtau by 10%. dtau = ", dtau
             errFlag = .true.
             exit iteration
          end if
       end if



       !  emergency exit 
       if( poissonSolverType=="adi" ) then 
          if (res > 1.0) then
             if (giveInfo)  then 
                print*,"Emergency exit from ADI: no convergence "
                write(*,fmt="(a25,i25)") " Nb.of iterations: iter = ", iter
                write(*,fmt="(a25,es17.6)") " Final residual: res = ", res
                print*,"----------------- stopping -------------------------"
                print*,""
                stop
             end if
             nIter = iter
             return
          end if
       end if
       !  emergency exit


    end do iteration

    ! Max iterations
    if( poissonSolverType=="adi" ) then  
       write(*,fmt="(a25,i25)") " ADI: max iterations!!!", &
            & maxIt
       call linOpr( var, sol, Lx3D ) ! -> r = L(x) - b     NOT b - L(x) !!!
       r = Lx3D - b
       res = dt*maxval(abs(r))     ! abort crit by Smolarkiewicz, unstable!?
       !  res = maxval(abs(r))            
       write(*,fmt="(a25,es17.6)") " Final residual: res = ", res
       print*,"--------------------------------------------------"; print*,""
       errFlag = .true.
       nIter = iter
    end if

    !---- Deallocate fields without ghost cells
    deallocate(r, stat=allocstat); if(allocstat/=0) stop "adi:dealloc failed"
    deallocate(Lx3D, stat=allocstat); if(allocstat/=0) stop "adi:dealloc failed"

    !--- Deallocate Lx, Ly, Lz
    deallocate(Lx, stat=allocstat); if(allocstat/=0) stop "adi:dealloc failed"
    deallocate(Ly, stat=allocstat); if(allocstat/=0) stop "adi:dealloc failed"
    deallocate(Lz, stat=allocstat); if(allocstat/=0) stop "adi:dealloc failed"

    !---- Deallocate fields with ghost cells: q, q01, q02
    deallocate(q, stat=allocstat); if(allocstat/=0) stop "adi:dealloc failed"
    deallocate(qn, stat=allocstat); if(allocstat/=0) stop "adi:dealloc failed"
    deallocate(q01, stat=allocstat); if(allocstat/=0) stop "adi:dealloc failed"
    deallocate(q02, stat=allocstat); if(allocstat/=0) stop "adi:dealloc failed"



  end subroutine adi


  !---------------------------------------------------------------------------


  subroutine adi_z(var, b, dt, sol, res, nIter, errFlag)
    ! --------------------------------------
    !   ADI sweep in z-direction
    !   cf. Skamarock et al., MWR 1997 
    !   can only be used as preconditioner
    !---------------------------------------

    ! In/Out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var
    real, dimension(1:nx,1:ny,1:nz), intent(in) :: b        ! RHS 
    real, intent(in) :: dt
    real, dimension(1:nx,1:ny,1:nz), intent(inout) :: sol
    real, intent(out) :: res                     ! residual
    integer, intent(out) :: nIter
    logical, intent(out) :: errFlag

    ! Local variables
    integer :: i,j,k, allocstat, iter, maxIt
    real :: pStratU, pStratD

    ! Local field
    real, dimension(1:nx,1:ny,1:nz) :: sol0, solNAG



    ! allocatable fields
    real, dimension(1:nx,1:ny,1:nz) :: r, Lx3D
    real, dimension(1:nx,1:ny,1:nz)  :: Lx, Ly, Lz
    real, dimension(0:nx+1,0:ny+1,0:nz+1)  :: q, qn
    real, dimension(1:nz)      :: left_a, center, right_a, rhs

    real :: rhoEdge, AD, AU, AC
    real :: dx2, dy2, dz2

    ! verbose
    logical, parameter :: giveInfo = .true.

    ! test variables
    real :: res0

    ! NAG routine variables
    integer :: info
    real, dimension(1:nz) :: diag, rhs0 
    real, dimension(1:nz-1) :: lower_diag, upper_diag

!   achatzb
    if(topography) stop 'adi_z not ready for topography!'
!   achatze

    ! define auxiliary variables
    dx2 = 1.0/dx**2
    dy2 = 1.0/dy**2
    dz2 = 1.0/dz**2

    maxIt = maxIterADI       ! as preconditioner only 1 or 2 iterations
    tol = 0.0

    ! Save input solution
    sol0 = sol


    ! Error flag
    errFlag = .false.    

    ! Init
    call linOpr( var, sol, Lx3D ) ! -> r = L(x) - b    NOT b - L(x) !!!
    r = Lx3D - b


    !---- Check initial residual
    res = dt*maxval(abs(r))             ! abort criterium by Smolarkiewicz
    !   res = maxval(abs(r))
    res0 = res
    print*,"  ADI "
    if (giveInfo) write(*,fmt="(a25,es17.6)") " Init residual: res0 = ", res
    if (giveInfo) write(*,fmt="(a25,es17.6)") " tol = ", tol
!!$    if (res <= tol) then
!!$       if(giveInfo) print*," ==> no iteration needed."
!!$       nIter = 0
!!$       return
!!$    end if


    ! Loop
    iteration: do iter = 1,maxIt

       ! Init qn with sol
       qn(1:nx,1:ny,1:nz) = sol

       !---- prepare explicit terms Lx, Ly for rhs

       call linOprXYZ( var, qn, Lx, "x" )
       call linOprXYZ( var, qn, Ly, "y" )

       !
       !---------------------------------
       !    sweep in z direction -> sol
       !---------------------------------
       !

       !---- set up matrix diagonals for Thomas algorithm
       q = qn

       j_loop_z: do j = 1,ny
          i_loop_z: do i = 1,nx

             !--- matrix diagonals
             k_loop_z: do k = 1,nz

                !  A(i,j,k+1) 
                if (k<nz) then

                   rhoEdge = 0.5* ( var(i,j,k+1,1) + var(i,j,k,1) )
                   AU = dz2 * pStratU**2 / rhoEdge
                   if( pressureScaling ) AU = AU / PstratTilde(k)

                else ! k = nz -> upwad boundary (solid wall) -> A(i,j,nz+1) = 0
                   AU = 0.0
                end if

                !  A(i,j,k-1) 
                if (k>1) then 

                   rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j,k-1,1) )
                   AD = dz2 * pStratD**2 / rhoEdge
                   if( pressureScaling ) AD = AD / PstratTilde(k-1)

                else ! k = 1 -> downward boundary (solid wall) -> A(i,j,0) = 0
                   AD = 0.0
                end if

                !  A(i,j,k) 
                AC = - AU - AD

                left_a(k) = - AD
                center(k) = 1.0/dtau - AC
                right_a(k) = - AU

             end do k_loop_z


             !---- Right hand side rhs

             rhs =  qn(i,j,1:nz) / dtau + &
                  & Lx(i,j,1:nz) + &
                  & Ly(i,j,1:nz) - &
                  & b(i,j,1:nz)



             !---- solve tridiagonal system
             if ( useNAG ) then

!                call DGTSV(nz,1,left_a(2:nz),center,right_a(1:nz-1),rhs,nz,info)
!                sol(i,j,1:nz) = rhs

                 print*,"no NAG option!!! ==> exit!"
                 STOP

             else
                call thomas( left_a, center, right_a, sol(i,j,1:nz), rhs )
             end if



          end do i_loop_z
       end do j_loop_z


       !
       !---------------------------------
       !    Evaluate result
       !---------------------------------
       !

       call linOpr( var, sol, Lx3D ) ! -> r = L(x) - b       NOT b - L(x) !!!
       r = Lx3D - b
       res = dt*maxval(abs(r)) ! abort crit by Smolarkiewicz, unstable!?
       !  res = maxval(abs(r))            
    end do iteration
    nIter = iter - 1

    write(*,fmt="(a25,es17.6)") " ADI residual = ", res
    write(*,fmt="(a25,f25.18)") " -> residual reduction = ", res/res0
    ! Return input / start solution if no residual reduction
    if( res/res0 >= 1.0 ) then
       print*,"==> no residual reduction!!!"
       dtau = 0.9*dtau
       write(*,fmt="(a25,f25.18)") " reduce dtau by 10%. dtau = ", dtau
       errFlag = .true.
    end if


  end subroutine adi_z


  !---------------------------------------------------------------------------


  subroutine thomas( l, c, r, sol, b ) 
    ! -------------------------------------
    !      Thomas algorithmus
    !      using diagonals l,c,r
    ! -------------------------------------

    ! In/Out variables
    real, dimension(:), intent(inout) :: l,c,r   ! matrix diagonals: 
    ! left, center, right
    real, dimension(:), intent(inout) :: b       ! rhs
    real, dimension(:), intent(out)   :: sol     ! solution

    ! Local vars
    integer :: i,n


    n = size(b)

    ! Forward sweep
    do i = 2,n
       c(i) = c(i) - l(i)/c(i-1)*r(i-1)
       b(i) = b(i) - l(i)/c(i-1)*b(i-1)       
    end do

    ! Backward sweep
    do i = n-1,1,-1
       b(i) = b(i) - r(i)/c(i+1)*b(i+1)
    end do

    ! Solve 
    do i = 1,n
       sol(i) = b(i)/c(i)
    end do


  end subroutine thomas


  !---------------------------------------------------------------------------


  subroutine gcr(var, b, dt, sol, res, nIter, errFlag)
    ! --------------------------------------
    !   GCR(k) scheme using linear operator
    !   cf. Skamarock et al., AMS 1997
    !---------------------------------------

    ! In/Out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var
    real, dimension(1:nx,1:ny,1:nz), intent(in) :: b        ! RHS 
    real, intent(in) :: dt
    real, dimension(1:nx,1:ny,1:nz), intent(inout) :: sol
    real, intent(out) :: res                     ! residual
    integer, intent(out) :: nIter
    logical, intent(out) :: errFlag

    ! Local parameters
    integer  :: maxIt

    ! Local variables
    integer :: i,j,k, allocstat
    real, dimension(:,:,:), allocatable :: r,Lx,p,Lp,q,Lq
    real :: alpha, beta

    ! Preconditioner variables
    logical :: errFlagPrec
    integer :: nIterPrec
    real :: dtPrec, resPrec


    ! verbose
    logical, parameter :: giveInfo = .true.


    !---- Set parameters
    maxIt = maxIterPoisson
    tol = tolPoisson


    !---- Allocate local fields
    allocate(r(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/gcr:alloc failed"
    allocate(Lx(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/gcr:alloc failed"
    allocate(p(1:nx,1:ny,1:nz), stat=allocstat);if(allocstat/=0) &
         & stop "algebra.f90/gcr:alloc failed"
    allocate(Lp(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/gcr:alloc failed"
    allocate(q(1:nx,1:ny,1:nz), stat=allocstat);  if(allocstat/=0) &
         & stop "algebra.f90/gcr:alloc failed"
    allocate(Lq(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/gcr:alloc failed"

    ! Error flag
    errFlag = .false.    

    ! Init
    call linOpr( var, sol, Lx ) ! -> r = b - L(x)
    r = Lx - b

    select case( preconditioner )

    case( "adi" )
       dtPrec = 1.0
       p = 0.0
       call adi(var, r, dtPrec, p, resPrec, nIterPrec, errFlagPrec)
       if (errFlagPrec) p = r   ! ignore adi in case of error
    case( "adi_z" )
       dtPrec = 1.0
       p = 0.0
       call adi(var, r, dtPrec, p, resPrec, nIterPrec, errFlagPrec)
       if (errFlagPrec) p = r    ! ignore adi in case of error
    case default
       p = r       ! no preconditioning
    end select

    call linOpr( var, p, Lp )   ! -> Lp


    !---- Check initial residual
    res = dt*maxval(abs(r))             ! abort criterium by Smolarkiewicz
    !    res = maxval(abs(r)) 
    print*," GCR "
    if (giveInfo) write(*,fmt="(a25,es17.6)") " Init GCR residual = ", res
    if (giveInfo) write(*,fmt="(a25,es17.6)") " tol = ", tol
!!$    if (res <= tol) then
!!$       if(giveInfo) print*," ==> no iteration needed."
!!$       nIter = 0
!!$       return
!!$    end if


    ! Loop
    iteration: do j = 1,maxIt

       beta = - dot_product3D(r,Lp) / dot_product3D(Lp,Lp)

       sol = sol + beta*p

       r = r + beta*Lp       

       !---- Abort criterium 
       res = dt*maxval(abs(r))    ! by Smolarkiewicz, unstable?!
       !  res = maxval(abs(r))            
       if (res <= tol) then
          if (giveInfo)  then 
             write(*,fmt="(a25,i25)") " Nb.of GCR iterations = ", j
             write(*,fmt="(a25,es17.6)") " Final GCR residual = ", res
             print*,"--------------------------------------------------"
             print*,""
          end if
          nIter = j
          return
       end if

       !---- Precondition residual Lp*q = r
       ! -> b = r
       ! -> sol = q
       ! -> nIter = 
       ! -> dt = 1.0
       select case( preconditioner )

       case( "adi" )
          dtPrec = 1.0
          q = 0.0       ! initial guess 
          call adi(var, r, dtPrec, q, resPrec, nIterPrec, errFlagPrec)
          if (errFlagPrec) q = r    ! ignore adi in case of error
       case( "adi_z" )
          dtPrec = 1.0
          q = 0.0       ! initial guess
          call adi_z(var, r, dtPrec, q, resPrec, nIterPrec, errFlagPrec)
          if (errFlagPrec) q = r    ! ignore adi in case of error
       case default
          q = r       ! no preconditioning
       end select


       call linOpr( var, q, Lq )   ! -> Lq

       alpha = - dot_product3D(Lq,Lp) / dot_product3D(Lp,Lp)

       p = q + alpha*p

       Lp = Lq + alpha*Lp

    end do iteration

    ! Max iterations
    write(*,fmt="(a25,i25)") " GCR: max iterations!!!", &
         & maxIt
    if (giveInfo)  then 
       ! res = dt*maxval(abs(r)) ! unstable
       res = maxval(abs(r))
       write(*,fmt="(a25,es17.6)") " Final GCR residual = ", res
    end if
    print*,"--------------------------------------------------"; print*,""
    errFlag = .true.
    nIter = j

    !---- Deallocate local fields: r,Lx,p,Lp,q,Lq
    deallocate(r, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/gcr:dealloc failed"
    deallocate(Lx, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/gcr:dealloc failed"
    deallocate(p,stat=allocstat);if(allocstat/=0) &
         & stop "algebra.f90/gcr:dealloc failed"
    deallocate(Lp, stat=allocstat);if(allocstat/=0) &
         & stop "algebra.f90/gcr:dealloc failed"
    deallocate(q, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/gcr:dealloc failed"
    deallocate(Lq, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/gcr:dealloc failed"


  end subroutine gcr


  !---------------------------------------------------------------------------


  subroutine bicgstab(var, b, dt, sol, res, nIter, errFlag, opt)
    ! --------------------------------------
    !    BiCGStab using linear operator
    !---------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var
    real, dimension(1:nx,1:ny,1:nz), intent(in) :: b        ! RHS 
    real, intent(in) :: dt
    real, dimension(1:nx,1:ny,1:nz), intent(inout) :: sol
    real, intent(out) :: res                     ! residual
    integer, intent(out) :: nIter
    logical, intent(out) :: errFlag
    character(len=*), intent(in)    :: opt
    

    ! Local parameters
    integer  :: maxIt

    ! local variables
    integer :: i,j,k, allocstat
    real, dimension(:,:,:), allocatable :: p,r0,rOld,r,s,t,v, matVec
    real :: alpha, beta, omega

    ! verbose
    logical, parameter :: giveInfo = .true.

    ! modified by Junhong Wei (20161107) *** starting line ***

    ! MPI stuff
    integer :: root
    real :: res_local

!xxxx
    if( giveInfo .and. master ) then
       print*,""
       print*,"----------------------------------------------"
       print*,"BICGSTAB: solving linear system... "
       print*,"----------------------------------------------"
       print*,""
    end if

    ! modified by Junhong Wei (20161107) *** finishing line ***

    ! Set parameters
    maxIt = maxIterPoisson
    if( opt == 'initial' ) then
       print*,"bicgstab: use tolInitial = ", tolInitial
       tol = tolInitial
    else
       tol = tolPoisson
    end if
    

    ! Allocate local fields
    allocate(p(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:alloc failed"
    allocate(r0(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:alloc failed"
    allocate(rOld(1:nx,1:ny,1:nz), stat=allocstat);if(allocstat/=0) &
         & stop "bicgstab:alloc failed"
    allocate(r(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:alloc failed"
    allocate(s(1:nx,1:ny,1:nz), stat=allocstat);  if(allocstat/=0) &
         & stop "bicgstab:alloc failed"
    allocate(t(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:alloc failed"
    allocate(v(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:alloc failed"
    allocate(matVec(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:alloc failed"

    ! error flag
    errFlag = .false.    

    ! Init
    ! r0 = b - Ax
    call linOpr( var, sol, matVec )
    r0 = b - matVec    
    p = r0
    r = r0

!    ! check initial residual   ! modified by Junhong Wei (20161107)
!    ! old    res = norm3D(r)   ! modified by Junhong Wei (20161107)
!    res = dt*maxval(abs(r))             ! abort criterion by Smolarkiewicz   ! modified by Junhong Wei (20161107)
!    !   res = maxval(abs(r))    ! modified by Junhong Wei (20161107)
!    print*,""   ! modified by Junhong Wei (20161107)
!    print*," BiCGStab solver: "   ! modified by Junhong Wei (20161107)
!    if (giveInfo) write(*,fmt="(a25,es17.6)") &   ! modified by Junhong Wei (20161107)
!         & " Initial residual: res = ", res   ! modified by Junhong Wei (20161107)
!    if (giveInfo) write(*,fmt="(a25,es17.6)") " tol = ", tol   ! modified by Junhong Wei (20161107)
! !!$    if (res <= tol) then   ! modified by Junhong Wei (20161107)
! !!$       if(giveInfo) print*," ==> no iteration needed."   ! modified by Junhong Wei (20161107)
! !!$       nIter = 0   ! modified by Junhong Wei (20161107)
! !!$       return   ! modified by Junhong Wei (20161107)
    ! !!$    end if   ! modified by Junhong Wei (20161107)

    

    ! modified by Junhong Wei (20161107) *** starting line ***

    !--------------------------------------
    !   Abort criterion by Smolarkiewicz
    !--------------------------------------
    ! check initial residual
    ! old    res = norm3D(r)
    
    res_local = dt*maxval(abs(r)) 
    
    !MPI find global residual
    root = 0
    call mpi_reduce(res_local, res, 1, mpi_double_precision,&
         & mpi_max, root, comm, ierror)
    
    call mpi_bcast(res, 1, mpi_double_precision, root, comm, ierror)

    !xxxx
    print"(a,i3,2es10.1)","bicgstab: rank, res_local, res", rank, res_local, res
    
    if( master ) then
       print*,""
       print*," BiCGStab solver: "
       if (giveInfo) write(*,fmt="(a25,es17.6)") &
            & " Initial residual: res = ", res
       if (giveInfo) write(*,fmt="(a25,es17.6)") " tol = ", tol
    end if
    if (res <= tol) then
       if(master .and. giveInfo) print*," ==> no iteration needed."
       nIter = 0
!xxx   wait for all processes and return
!       call mpi_barrier(comm,ierror)
       return
    end if

    ! modified by Junhong Wei (20161107) *** finishing line ***

    
    ! Loop
    iteration: do j = 1,maxIt

       ! v = A*p
       call linOpr( var, p, matVec )
       v = matVec

       !       alpha = dot_product3D(r,r0) / dot_product3D(v,r0) ! modified by Junhong Wei (20161107)
              alpha = dot_product3D_glob(r,r0) / dot_product3D_glob(v,r0) ! modified by Junhong Wei (20161107)
       s = r - alpha*v

       ! t = A*s
       call linOpr( var, s, matVec )
       t = matVec

       !       omega = dot_product3D(t,s) / dot_product3D(t,t) ! modified by Junhong Wei (20161107)
              omega = dot_product3D_glob(t,s) / dot_product3D_glob(t,t) ! modified by Junhong Wei (20161107)
       sol = sol + alpha*p + omega*s

       rOld = r
       r = s - omega*t

       ! modified by Junhong Wei (20161107) *** starting line ***

!       !---- Abort criterium 
!       ! old    res = norm3D(r)
!       res = dt*maxval(abs(r))     ! by Smolarkiewicz
!       !  res = maxval(abs(r))
!       if (res <= tol) then
!          if (giveInfo)  then 
!             write(*,fmt="(a25,i25)") " Nb.of iterations: j = ", j
!             write(*,fmt="(a25,es17.6)") " Final residual: res = ", res
!             print*,""
!          end if
!          nIter = j
!          return
!       end if
!
!       beta = alpha/omega * dot_product3D(r,r0) / dot_product3D(rOld,r0)
       !       p = r + beta*(p-omega*v)

       !-----------------------
       !   Abort criterium 
       !-----------------------
       res_local = dt*maxval(abs(r))     ! by Smolarkiewicz

       !MPI find global residual
       root = 0
       call mpi_reduce(res_local, res, 1, mpi_double_precision,&
            & mpi_max, root, comm, ierror)
       
       call mpi_bcast(res, 1, mpi_double_precision, root, comm, ierror)

      
       if (res <= tol) then
          if (master .and. giveInfo)  then 
             write(*,fmt="(a25,i25)") " Nb.of iterations: j = ", j
             write(*,fmt="(a25,es17.6)") " Final residual: res = ", res
             print*,""
          end if
          nIter = j
!xxx      wait for all processes and return
!          call mpi_barrier(comm, ierror)
          return
       end if

       beta = alpha/omega * dot_product3D_glob(r,r0) / dot_product3D_glob(rOld,r0)
       p = r + beta*(p-omega*v)

       ! modified by Junhong Wei (20161107) *** finishing line ***

    end do iteration

    ! max iterations
    if( master ) then ! modified by Junhong Wei (20161107)
    write(*,fmt="(a25,i25)") " Bicgstab: max iterations!!!", &
         & maxIt
    if (giveInfo)  then 
       res = dt*maxval(abs(r))
       write(*,fmt="(a25,es17.6)") " Final BICGSTAB residual = ", res
    end if
    print*,"--------------------------------------------------"; print*,""
    end if ! modified by Junhong Wei (20161107)
    errFlag = .true.
    nIter = j


    ! deallocate local fields
    deallocate(p, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(r0, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(rOld,stat=allocstat);if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(r, stat=allocstat);if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(s, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(t, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(v, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"



  end subroutine bicgstab

  ! modified by Junhong Wei (2016/07/08) (starting line)

  !---------------------------------------------------------------------------


  subroutine hypre(var, b, dt, sol, res, nIter, errFlag, opt)
    ! --------------------------------------
    !    HYPRE
    !---------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var
    real, dimension(1:nx,1:ny,1:nz), intent(in) :: b        ! RHS 
    real, intent(in) :: dt
    real, dimension(1:nx,1:ny,1:nz), intent(inout) :: sol
    real, intent(out) :: res                     ! residual
    integer, intent(out) :: nIter
    logical, intent(out) :: errFlag
    character(len=*), intent(in)    :: opt

    ! variables due to the use of HYPRE

    integer :: ndim_hypre
    parameter (ndim_hypre = 3)

    integer*8 grid_hypre, stencil_hypre, A_hypre, b_hypre, x_hypre, solver_hypre
    integer ierr_hypre
    integer npts_x_periodic_hypre, npts_y_periodic_hypre, npts_z_periodic_hypre

    integer ii_entry_hypre
    integer offsets_hypre(7,3)

    integer :: nentries_hypre
    parameter (nentries_hypre = 7)
    integer :: nvalues_hypre

    integer :: allocstat

    real, dimension(:), allocatable :: values_hypre

    integer :: stencil_indices_hypre(nentries_hypre)

    integer :: index_count_hypre

!   achatzb
    integer :: k_sing
    integer :: i0,j0
!   achatze

    ! local variables
    integer :: i,j,k
    real :: pStratU, pStratD, rhoEdge
    real :: AL,AR, AB,AF, AD,AU, AC

    real :: dx2, dy2, dz2

    real, dimension(:), allocatable :: bvalue_vector_hypre, xvalue_vector_hypre

    ! Local parameters
    integer  :: maxIt

    ! verbose
    logical, parameter :: giveInfo = .true.

    real :: sol_mean_hypre

    ! auxiliary variables
    dx2 = 1.0/dx**2
    dy2 = 1.0/dy**2
    dz2 = 1.0/dz**2

    ! Set parameters
    maxIt = maxIterPoisson
    ! Set parameters
    if( opt == 'initial' ) then
       print*,"hypre: use tolInitial = ", tolInitial
       tol = tolInitial
    else
!      achatzb
!      modified convergence criterion so that iterations stop when the
!      norm of the normalized RHS falls below tolPoisson
!      tol = tolPoisson
       if (tolref /= 0.0 ) then
           tol = tolPoisson/tolref
          else
           tol = tolPoisson
       end if
       if(master) print*,"hypre uses tolerance tol = ", tol
!      achatze
    end if

    ! Allocate local fields
    nvalues_hypre = nx * ny * nz * nentries_hypre

    allocate(values_hypre(1:nvalues_hypre), stat=allocstat); if(allocstat/=0) &
         & stop "hypre:alloc failed"
    allocate(bvalue_vector_hypre(1: (nx * ny * nz) ), stat=allocstat); if(allocstat/=0) &
         & stop "hypre:alloc failed"
    allocate(xvalue_vector_hypre(1: (nx * ny * nz) ), stat=allocstat); if(allocstat/=0) &
         & stop "hypre:alloc failed"

    ! error flag
    errFlag = .false.

!---------
! step one
! Set up a grid. Each processor describes the piece of the grid that it owns.

! Create an empty 2D grid object
   call HYPRE_StructGridCreate(mpi_comm_world, ndim_hypre, grid_hypre, ierr_hypre)

    !----------------
    !  x-Boundary
    !----------------
    select case( xBoundary )

    case( "periodic" )
!       npts_x_periodic_hypre = nx
       npts_x_periodic_hypre = sizeX
    case default
       stop "pressureBoundaryCondition: unknown case xBoundary."
    end select

    !----------------
    !   y-Boundary
    !----------------
    select case( yBoundary )

    case( "periodic" )
!       npts_y_periodic_hypre = ny
       npts_y_periodic_hypre = sizeY
    case default
       stop "pressureBoundaryCondition: unknown case yBoundary."
    end select

    !----------------
    !   z-Boundary
    !----------------
    select case( zBoundary )

    case( "periodic" )
!       npts_z_periodic_hypre = nz
       npts_z_periodic_hypre = sizeZ
    case( "solid_wall" )
       npts_z_periodic_hypre = 0
    case default
       stop "pressureBoundaryCondition: unknown case zBoundary."
    end select


   call HYPRE_StructGridSetPeriodic(grid_hypre,(/npts_x_periodic_hypre, npts_y_periodic_hypre, npts_z_periodic_hypre/),ierr_hypre)

! Add boxes to the grid
   call HYPRE_StructGridSetExtents(grid_hypre, (/ (icoord-1)*nx , (jcoord-1)*ny , 0 /), &
     & (/ (icoord*nx)-1 , (jcoord*ny)-1 , nz-1 /), ierr_hypre)

! This is a collective call finalizing the grid assembly. The grid is now ``ready to be used''
   call HYPRE_StructGridAssemble(grid_hypre, ierr_hypre)

!---------
! step two
! Define the discretization stencil

! Create an empty 3D, 7-pt stencil object
   call HYPRE_StructStencilCreate(3, 7, stencil_hypre, ierr_hypre)

! Define the geometry of the stencil. Each represents a relative offset (in the index space).
   offsets_hypre(1,:) = (/ 0, 0, 0/)
   offsets_hypre(2,:) = (/-1, 0, 0/)
   offsets_hypre(3,:) = (/ 1, 0, 0/)
   offsets_hypre(4,:) = (/ 0,-1, 0/)
   offsets_hypre(5,:) = (/ 0, 1, 0/)
   offsets_hypre(6,:) = (/ 0, 0,-1/)
   offsets_hypre(7,:) = (/ 0, 0, 1/)
   
   do ii_entry_hypre=0,6
     call HYPRE_StructStencilSetElement( stencil_hypre, ii_entry_hypre, offsets_hypre( (ii_entry_hypre+1), : ), ierr_hypre)
   end do

!---------
! step three
! Set up a Struct Matrix

! Create an empty matrix object
   call HYPRE_StructMatrixCreate(mpi_comm_world, grid_hypre, stencil_hypre, A_hypre, ierr_hypre)

! Indicate that the matrix coefficients are ready to be set
   call HYPRE_StructMatrixInitialize(A_hypre, ierr_hypre)

! Set the matrix coefficients.  Each processor assigns coefficients
! for the boxes in the grid that it owns. Note that the coefficients
! associated with each stencil entry may vary from grid point to grid
! point if desired.  Here, we first set the same stencil entries for
! each grid point.  Then we make modifications to grid points near
! the boundary.

! labels for the stencil entries - these correspond to the offsets defined above

    stencil_indices_hypre(:) = (/0,1,2,3,4,5,6/)

! We have *** grid points, each with 7 stencil entrie

!   achatzb
    i0=is+nbx-1
    j0=js+nby-1
!   achatze

    select case( model ) 

       !----------------------------------------
       !       Pseudo-incompressible model
       !----------------------------------------
    case( "pseudo_incompressible" )


       !---------------------------------
       !         Loop over field
       !---------------------------------

       k_loop: do k = 1,nz
          j_loop: do j = 1,ny
             i_loop: do i = 1,nx

!               achatzb
!               ! ------------------ A(i+1,j,k) ------------------
!               rhoEdge = 0.5*( var(i+1,j,k,1) + var(i,j,k,1) )
!               if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

!               AR = dx2 * pStrat(k)**2/rhoEdge

!               ! ------------------- A(i-1,j,k) --------------------
!               rhoEdge = 0.5*( var(i,j,k,1) + var(i-1,j,k,1) )
!               if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

!               AL = dx2 * pStrat(k)**2 / rhoEdge

!               ! -------------------- A(i,j+1,k) ----------------------
!               rhoEdge = 0.5*( var(i,j+1,k,1) + var(i,j,k,1) )
!               if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

!               AF = dy2 * pStrat(k)**2 / rhoEdge

!               ! --------------------- A(i,j-1,k) ----------------------
!               rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j-1,k,1) )
!               if( fluctuationMode ) then
!                  rhoEdge = rhoEdge + rhoStrat(k)
!               end if

!               AB = dy2 * pStrat(k)**2 / rhoEdge

!               ! ---------------------- A(i,j,k+1) ---------------------
!               if (k<nz) then
!                  rhoEdge = 0.5* ( var(i,j,k+1,1) + var(i,j,k,1) )
!                  if( fluctuationMode ) then
!                     rhoEdge = rhoEdge + rhoStratTilde(k)
!                  end if

!                  pStratU = 0.5* ( pStrat(k+1) + pStrat(k) )
!                  AU = dz2 * pStratU**2 / rhoEdge
!               else ! k = nz -> upwad boundary (solid wall)
!                  ! A(i,j,nz+1) = 0
!                  AU = 0.0
!               end if

!               ! ----------------------- A(i,j,k-1) ---------------------
!               if (k>1) then 
!                  rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j,k-1,1) )
!                  if( fluctuationMode ) then
!                      rhoEdge = rhoEdge + rhoStratTilde(k-1)
!                  end if
!                  
!                  pStratD = 0.5* ( pStrat(k) + pStrat(k-1) )
!                  AD = dz2 * pStratD**2 / rhoEdge
!               else ! k = 1 -> downward boundary (solid wall)
!                  ! A(i,j,0) = 0
!                  AD = 0.0
!               end if

!               ! ----------------- scale with thetaStrat ----------------
!               if( pressureScaling ) then
!                  AL = AL / PStrat(k)
!                  AR = AR / PStrat(k)
!                  AB = AB / PStrat(k)
!                  AF = AF / PStrat(k)
!                  AD = AD / PStrat(k)
!                  AU = AU / PStrat(k)
!               end if

!               ! ----------------------- A(i,j,k) -----------------------
!               AC = - AR - AL - AF - AB - AU - AD

                if(topography) then
                   ! stencil with topography
                   ! ------------------ A(i+1,j,k) ------------------
                   if((topography_mask(i0+i,j0+j,k)&
                       .and.(.not.topography_mask(i0+i+1,j0+j,k)))&
                      .or.&
                      ((.not.topography_mask(i0+i,j0+j,k))&
                       .and.(topography_mask(i0+i+1,j0+j,k)))) then
                      AR=0.0
                     else
                      rhoEdge = 0.5*( var(i+1,j,k,1) + var(i,j,k,1) )
                      if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

                      AR = dx2 * pStrat(k)**2/rhoEdge
                   end if
   
                   ! ------------------- A(i-1,j,k) --------------------
                   if((topography_mask(i0+i-1,j0+j,k)&
                       .and.(.not.topography_mask(i0+i,j0+j,k)))&
                      .or.&
                      ((.not.topography_mask(i0+i-1,j0+j,k))&
                       .and.(topography_mask(i0+i,j0+j,k)))) then
                      AL=0.0
                     else
                      rhoEdge = 0.5*( var(i,j,k,1) + var(i-1,j,k,1) )
                      if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

                      AL = dx2 * pStrat(k)**2 / rhoEdge
                   end if

                   ! -------------------- A(i,j+1,k) ----------------------
                   if((topography_mask(i0+i,j0+j,k)&
                       .and.(.not.topography_mask(i0+i,j0+j+1,k)))&
                      .or.&
                      ((.not.topography_mask(i0+i,j0+j,k))&
                       .and.(topography_mask(i0+i,j0+j+1,k)))) then
                      AF=0.0
                     else
                      rhoEdge = 0.5*( var(i,j+1,k,1) + var(i,j,k,1) )
                      if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

                      AF = dy2 * pStrat(k)**2 / rhoEdge
                   end if

                   ! --------------------- A(i,j-1,k) -------------------
                   if((topography_mask(i0+i,j0+j-1,k)&
                       .and.(.not.topography_mask(i0+i,j0+j,k)))&
                      .or.&
                      ((.not.topography_mask(i0+i,j0+j-1,k))&
                       .and.(topography_mask(i0+i,j0+j,k)))) then
                      AB=0.0
                     else
                      rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j-1,k,1) )
                      if( fluctuationMode ) then
                         rhoEdge = rhoEdge + rhoStrat(k)
                      end if

                      AB = dy2 * pStrat(k)**2 / rhoEdge
                   end if

                   ! ---------------------- A(i,j,k+1) ------------------
                   if((topography_mask(i0+i,j0+j,k)&
                       .and.(.not.topography_mask(i0+i,j0+j,k+1)))&
                      .or.&
                      ((.not.topography_mask(i0+i,j0+j,k))&
                       .and.(topography_mask(i0+i,j0+j,k+1)))&
                      .or.&
                      (k == nz)) then
                      AU=0.0
                     else
                      rhoEdge = 0.5* ( var(i,j,k+1,1) + var(i,j,k,1) )
                      if( fluctuationMode ) then
                         rhoEdge = rhoEdge + rhoStratTilde(k)
                      end if

                      pStratU = 0.5* ( pStrat(k+1) + pStrat(k) )
                      AU = dz2 * pStratU**2 / rhoEdge
                   end if

                   ! ----------------------- A(i,j,k-1) -----------------
                   if((topography_mask(i0+i,j0+j,k-1)&
                       .and.(.not.topography_mask(i0+i,j0+j,k)))&
                      .or.&
                      ((.not.topography_mask(i0+i,j0+j,k-1))&
                       .and.(topography_mask(i0+i,j0+j,k)))&
                      .or.&
                      (k == 1)) then
                      AD=0.0
                     else
                      rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j,k-1,1) )
                      if( fluctuationMode ) then
                          rhoEdge = rhoEdge + rhoStratTilde(k-1)
                      end if
                   
                      pStratD = 0.5* ( pStrat(k) + pStrat(k-1) )
                      AD = dz2 * pStratD**2 / rhoEdge
                   end if

                   ! ----------------- scale with thetaStrat ------------
                   if( pressureScaling ) then
                      AL = AL / PStrat(k)
                      AR = AR / PStrat(k)
                      AB = AB / PStrat(k)
                      AF = AF / PStrat(k)
                      AD = AD / PStrat(k)
                      AU = AU / PStrat(k)
                   end if

                   ! ----------------------- A(i,j,k) -------------------
                   ! avoid singularities in case a point is surrounded by
                   ! topographic borders
                   if((AR == 0.0).and.&
                      (AL == 0.0).and.&
                      (AF == 0.0).and.&
                      (AB == 0.0).and.&
                      (AU == 0.0).and.&
                      (AD == 0.0)) then
                      ! take density from lowermost non-topographic cell
                      if(fluctuationMode ) then
                         rhoEdge = var(i,j,k+1,1) + rhoStratTilde(k)
                        else
                         rhoEdge = var(i,j,k+1,1)
                      end if
                      
                      AC = dz2 * pStrat(k)**2 / rhoEdge
                     else
                      AC = - AR - AL - AF - AB - AU - AD
                   end if
                  else
                   ! stencil without topography
                   ! ------------------ A(i+1,j,k) ------------------
                   rhoEdge = 0.5*( var(i+1,j,k,1) + var(i,j,k,1) )
                   if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

                   AR = dx2 * pStrat(k)**2/rhoEdge
   
                   ! ------------------- A(i-1,j,k) --------------------
                   rhoEdge = 0.5*( var(i,j,k,1) + var(i-1,j,k,1) )
                   if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

                   AL = dx2 * pStrat(k)**2 / rhoEdge

                   ! -------------------- A(i,j+1,k) ----------------------
                   rhoEdge = 0.5*( var(i,j+1,k,1) + var(i,j,k,1) )
                   if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)

                   AF = dy2 * pStrat(k)**2 / rhoEdge

                   ! --------------------- A(i,j-1,k) -------------------
                   rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j-1,k,1) )
                   if( fluctuationMode ) then
                      rhoEdge = rhoEdge + rhoStrat(k)
                   end if

                   AB = dy2 * pStrat(k)**2 / rhoEdge

                   ! ---------------------- A(i,j,k+1) ------------------
                   if (k<nz) then
                      rhoEdge = 0.5* ( var(i,j,k+1,1) + var(i,j,k,1) )
                      if( fluctuationMode ) then
                         rhoEdge = rhoEdge + rhoStratTilde(k)
                      end if

                      pStratU = 0.5* ( pStrat(k+1) + pStrat(k) )
                      AU = dz2 * pStratU**2 / rhoEdge
                   else ! k = nz -> upwad boundary (solid wall)
                      ! A(i,j,nz+1) = 0
                      AU = 0.0
                   end if

                   ! ----------------------- A(i,j,k-1) -----------------
                   if (k>1) then 
                      rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j,k-1,1) )
                      if( fluctuationMode ) then
                          rhoEdge = rhoEdge + rhoStratTilde(k-1)
                      end if
                   
                      pStratD = 0.5* ( pStrat(k) + pStrat(k-1) )
                      AD = dz2 * pStratD**2 / rhoEdge
                   else ! k = 1 -> downward boundary (solid wall)
                      ! A(i,j,0) = 0
                      AD = 0.0
                   end if

                   ! ----------------- scale with thetaStrat ------------
                   if( pressureScaling ) then
                      AL = AL / PStrat(k)
                      AR = AR / PStrat(k)
                      AB = AB / PStrat(k)
                      AF = AF / PStrat(k)
                      AD = AD / PStrat(k)
                      AU = AU / PStrat(k)
                   end if

                   ! ----------------------- A(i,j,k) -------------------
                   AC = - AR - AL - AF - AB - AU - AD
                end if
!               achatze


                ! ------------------- define matrix A -------------------

!                  index_count_hypre 
!                  = ( i * j * k * nentries_hypre ) - nentries_hypre + 1

                  index_count_hypre = i
                  index_count_hypre = index_count_hypre + ( (j-1)*nx )
                  index_count_hypre = index_count_hypre + ( (k-1)*nx*ny )

                  index_count_hypre &
                  = ( index_count_hypre * nentries_hypre ) &
                    - nentries_hypre + 1

                  values_hypre(index_count_hypre)   = AC
                  values_hypre(index_count_hypre+1) = AL
                  values_hypre(index_count_hypre+2) = AR
                  values_hypre(index_count_hypre+3) = AB
                  values_hypre(index_count_hypre+4) = AF
                  values_hypre(index_count_hypre+5) = AD
                  values_hypre(index_count_hypre+6) = AU


             end do i_loop
          end do j_loop
       end do k_loop



    case( "Boussinesq" )
       !----------------------------------------
       !             Boussinesq model
       !----------------------------------------


       !---------------------------------
       !         Loop over field
       !---------------------------------

!      achatzb
!      do k = 1,nz
!         do j = 1,ny
!            do i = 1,nx

!               ! ------------------ A(i+1,j,k) ------------------
!               AR = dx2

!               ! ------------------- A(i-1,j,k) --------------------
!               AL = dx2

!               ! -------------------- A(i,j+1,k) ----------------------
!               AF = dy2

!               ! --------------------- A(i,j-1,k) -----------------------
!               AB = dy2

!               ! ---------------------- A(i,j,k+1) -----------------------
!               AU = dz2

!               ! ----------------------- A(i,j,k-1) ----------------------
!               AD = dz2

!               ! ----------------------- A(i,j,k) ------------------------
!               AC = - AR - AL - AF - AB - AU - AD

!               ! ------------------- define matrix A --------------------

!                  index_count_hypre &
!                  = ( i * j * k * nentries_hypre ) - nentries_hypre + 1

!                 index_count_hypre = i
!                 index_count_hypre = index_count_hypre + ( (j-1)*nx )
!                 index_count_hypre = index_count_hypre + ( (k-1)*nx*ny )

!                 index_count_hypre &
!                 = ( index_count_hypre * nentries_hypre ) &
!                   - nentries_hypre + 1

!                 values_hypre(index_count_hypre)   = AC
!                 values_hypre(index_count_hypre+1) = AL
!                 values_hypre(index_count_hypre+2) = AR
!                 values_hypre(index_count_hypre+3) = AB
!                 values_hypre(index_count_hypre+4) = AF
!                 values_hypre(index_count_hypre+5) = AD
!                 values_hypre(index_count_hypre+6) = AU

!            end do
!         end do
!      end do

!      !--------------------------------------
!      !  Correction for solid wall boundary
!      !--------------------------------------

!      !---------------------------------
!      !          zBoundary
!      !---------------------------------

!      if( zBoundary == "solid_wall" ) then

!         do k = 1,nz,nz-1      ! k = 1 and nz
!            do j = 1,ny
!               do i = 1,nx

!                  ! ------------------ A(i+1,j,k) ------------------
!                  AR = dx2

!                  ! ------------------- A(i-1,j,k) --------------------
!                  AL = dx2

!                  ! -------------------- A(i,j+1,k) ----------------------
!                  AF = dy2

!                  ! --------------------- A(i,j-1,k) ---------------------
!                  AB = dy2

!                  ! ---------------------- A(i,j,k+1) --------------------
!                  ! k = nz -> upwad boundary (solid wall)
!                  ! A(i,j,nz+1) = 0
!                  AU = 0.0

!                  ! ----------------------- A(i,j,k-1) ------------------
!                  ! k = 1 -> downward boundary (solid wall)
!                  ! A(i,j,0) = 0
!                  AD = 0.0

!                  ! ----------------------- A(i,j,k) --------------------

!                  AC = - AR - AL - AF - AB - AU - AD

!               ! ------------------- define matrix A --------------------

!                  index_count_hypre &
!                  = ( i * j * k * nentries_hypre ) - nentries_hypre + 1

!                 index_count_hypre = i
!                 index_count_hypre = index_count_hypre + ( (j-1)*nx )
!                 index_count_hypre = index_count_hypre + ( (k-1)*nx*ny )

!                 index_count_hypre &
!                 = ( index_count_hypre * nentries_hypre ) &
!                   - nentries_hypre + 1

!                 values_hypre(index_count_hypre)   = AC
!                 values_hypre(index_count_hypre+1) = AL
!                 values_hypre(index_count_hypre+2) = AR
!                 values_hypre(index_count_hypre+3) = AB
!                 values_hypre(index_count_hypre+4) = AF
!                 values_hypre(index_count_hypre+5) = AD
!                 values_hypre(index_count_hypre+6) = AU

!               end do
!            end do
!         end do

!      end if ! zBoundary

       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                if(topography) then
                   ! stencil with topography

                   ! ------------------ A(i+1,j,k) ------------------
                   if((topography_mask(i0+i,j0+j,k)&
                       .and.(.not.topography_mask(i0+i+1,j0+j,k)))&
                      .or.&
                      ((.not.topography_mask(i0+i,j0+j,k))&
                       .and.(topography_mask(i0+i+1,j0+j,k)))) then
                      AR=0.0
                     else
                      AR = dx2
                   end if
   
                   ! ------------------- A(i-1,j,k) --------------------
                   if((topography_mask(i0+i-1,j0+j,k)&
                       .and.(.not.topography_mask(i0+i,j0+j,k)))&
                      .or.&
                      ((.not.topography_mask(i0+i-1,j0+j,k))&
                       .and.(topography_mask(i0+i,j0+j,k)))) then
                      AL=0.0
                     else
                      AL = dx2
                   end if
   
                   ! -------------------- A(i,j+1,k) ----------------------
                   if((topography_mask(i0+i,j0+j,k)&
                       .and.(.not.topography_mask(i0+i,j0+j+1,k)))&
                      .or.&
                      ((.not.topography_mask(i0+i,j0+j,k))&
                       .and.(topography_mask(i0+i,j0+j+1,k)))) then
                      AF=0.0
                     else
                      AF = dy2
                   end if
   
                   ! --------------------- A(i,j-1,k) ---------------------
                   if((topography_mask(i0+i,j0+j-1,k)&
                       .and.(.not.topography_mask(i0+i,j0+j,k)))&
                      .or.&
                      ((.not.topography_mask(i0+i,j0+j-1,k))&
                       .and.(topography_mask(i0+i,j0+j,k)))) then
                      AB=0.0
                     else
                      AB = dy2
                   end if

                   ! ---------------------- A(i,j,k+1) --------------------
                   if((topography_mask(i0+i,j0+j,k)&
                       .and.(.not.topography_mask(i0+i,j0+j,k+1)))&
                      .or.&
                      ((.not.topography_mask(i0+i,j0+j,k))&
                       .and.(topography_mask(i0+i,j0+j,k+1)))) then
                      AU=0.0
                     else if((zBoundary == "solid_wall").and.(k == nz)) &
                     then
                      AU=0.0
                     else
                      AU = dz2
                   end if
   
                   ! ----------------------- A(i,j,k-1) -------------------
                   if((topography_mask(i0+i,j0+j,k-1)&
                       .and.(.not.topography_mask(i0+i,j0+j,k)))&
                      .or.&
                      ((.not.topography_mask(i0+i,j0+j,k-1))&
                       .and.(topography_mask(i0+i,j0+j,k)))) then
                      AD=0.0
                     else if((zBoundary == "solid_wall").and.(k == 1)) &
                     then
                      AD=0.0
                     else
                      AD = dz2
                   end if
   
                   ! ----------------------- A(i,j,k) ---------------------
                   ! avoid singularities in case a point is surrounded by
                   ! topographic borders
                   if((AR == 0.0).and.&
                      (AL == 0.0).and.&
                      (AF == 0.0).and.&
                      (AB == 0.0).and.&
                      (AU == 0.0).and.&
                      (AD == 0.0)) then
                      AC = dz2
                     else
                      AC = - AR - AL - AF - AB - AU - AD
                   end if
                  else
                   ! stencil without topography

                   ! ------------------ A(i+1,j,k) ------------------
                   AR = dx2

                   ! ------------------- A(i-1,j,k) --------------------
                   AL = dx2

                   ! -------------------- A(i,j+1,k) ----------------------
                   AF = dy2

                   ! --------------------- A(i,j-1,k) --------------------
                   AB = dy2

                   ! ---------------------- A(i,j,k+1) --------------------
                   if((zBoundary == "solid_wall").and.(k == nz)) then
                      AU = 0.0
                     else
                      AU = dz2
                   end if

                   ! ----------------------- A(i,j,k-1) -------------------
                   if((zBoundary == "solid_wall").and.(k == 1)) then
                      AD = 0.0
                     else
                      AD = dz2
                   end if

                   ! ----------------------- A(i,j,k) ---------------------
                   AC = - AR - AL - AF - AB - AU - AD
                end if

                ! ------------------- define matrix A --------------------

!                  index_count_hypre &
!                  = ( i * j * k * nentries_hypre ) - nentries_hypre + 1

                index_count_hypre = i
                index_count_hypre = index_count_hypre + ( (j-1)*nx )
                index_count_hypre = index_count_hypre + ( (k-1)*nx*ny )

                index_count_hypre &
                = ( index_count_hypre * nentries_hypre ) &
                  - nentries_hypre + 1

                values_hypre(index_count_hypre)   = AC
                values_hypre(index_count_hypre+1) = AL
                values_hypre(index_count_hypre+2) = AR
                values_hypre(index_count_hypre+3) = AB
                values_hypre(index_count_hypre+4) = AF
                values_hypre(index_count_hypre+5) = AD
                values_hypre(index_count_hypre+6) = AU
             end do
          end do
       end do
!      achatze


    case default
       stop "linOpr: unknown case model"
    end select


    call HYPRE_StructMatrixSetBoxValues(A_hypre, (/ (icoord-1)*nx , (jcoord-1)*ny , 0 /),  &
      & (/ (icoord*nx)-1 , (jcoord*ny)-1 , nz-1 /), nentries_hypre, stencil_indices_hypre, &
      & values_hypre, ierr_hypre)

! This is a collective call finalizing the matrix assembly. The matrix is now ``ready to be used'

call HYPRE_StructMatrixAssemble(A_hypre,ierr_hypre)


!---------
! step four
! Set up Struct Vectors for b and x.  Each processor sets the vectors corresponding to its boxes.

! Create an empty vector object
  call HYPRE_StructVectorCreate(mpi_comm_world, grid_hypre, b_hypre, ierr_hypre)
  call HYPRE_StructVectorCreate(mpi_comm_world, grid_hypre, x_hypre, ierr_hypre)

! Indicate that the vector coefficients are ready to be set
  call HYPRE_StructVectorInitialize(b_hypre, ierr_hypre)
  call HYPRE_StructVectorInitialize(x_hypre, ierr_hypre)

! Set the vector coefficients

       index_count_hypre = 1
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx

             bvalue_vector_hypre(index_count_hypre) = b(i,j,k)
             index_count_hypre = index_count_hypre + 1

             end do
          end do
       end do


       index_count_hypre = 1
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx

             xvalue_vector_hypre(index_count_hypre) = sol(i,j,k)
             index_count_hypre = index_count_hypre + 1

             end do
          end do
       end do


   call HYPRE_StructVectorSetBoxValues(b_hypre, (/ (icoord-1)*nx , (jcoord-1)*ny , 0 /), &
     & (/ (icoord*nx)-1 , (jcoord*ny)-1 , nz-1 /), bvalue_vector_hypre, ierr_hypre)
   call HYPRE_StructVectorSetBoxValues(x_hypre, (/ (icoord-1)*nx , (jcoord-1)*ny , 0 /), &
     & (/ (icoord*nx)-1 , (jcoord*ny)-1 , nz-1 /), xvalue_vector_hypre, ierr_hypre)

! This is a collective call finalizing the vector assembly. The vectors are now ``ready to be used''
   call HYPRE_StructVectorAssemble(b_hypre, ierr_hypre)
   call HYPRE_StructVectorAssemble(x_hypre, ierr_hypre)


!---------
! step five
! Set up and use a solver (See the Reference Manual for descriptions of all of the options.)

! Create an empty Hybrid Struct solver
  call HYPRE_StructHybridCreate(mpi_comm_world, solver_hypre, ierr_hypre)

! Set print level (set 0 for complete silence or 2 for info on iterations)
  call HYPRE_StructHybridSetPrintLevel(solver_hypre, 0, ierr_hypre)

! Set the type of Krylov solver to use. 
! Current krylov methods set by solver type are: 1 – DSCG (default) 2 – GMRES 3 – BiCGSTAB
! Junhong used 2
  call HYPRE_StructHybridSetSolverType(solver_hypre, 2, ierr_hypre)
! call HYPRE_StructHybridSetSolverType(solver_hypre, 2, ierr_hypre)

! Set tolerance value for relative! residual (unlike in Rieper's BiGSTAB
! where the absolute! residual is used). 
! tol = (tolPoisson in the namelist)/tolref
  call HYPRE_StructHybridSetTol(solver_hypre, tol, ierr_hypre)

! Uncomment the line below to have the tolerance for absolute! residual 
! (it works only for BiCGSTAB, i.e. solver_type = 3)
  !call HYPRE_StructHybridSetStopCrit(solver_hypre, 0, ierr_hypre)

! Set convergence criterion for activating a preconditioner. 
! A number between 0 and 1 is accepted, the smaller the number, the 
! sooner a preconditioner is activated.
! achatzb
! call HYPRE_StructHybridSetConvergenc(solver_hypre, 1.e-23, ierr_hypre)
  call HYPRE_StructHybridSetConvergenc(solver_hypre, tolCond, ierr_hypre)
! achatze

! Set maximum number of iterations. maxIt = maxIterPoisson in the namelist
  call HYPRE_StructHybridSetPCGMaxIter(solver_hypre, maxIt, ierr_hypre)

! Setup and solve
  call HYPRE_StructHybridSetup(solver_hypre, A_hypre, b_hypre, x_hypre, ierr_hypre)
  call HYPRE_StructHybridSolve(solver_hypre, A_hypre, b_hypre, x_hypre, ierr_hypre)

! Get the results
  call HYPRE_StructVectorGetBoxValues(x_hypre, (/ (icoord-1)*nx , (jcoord-1)*ny , 0 /), &
    & (/ (icoord*nx)-1 , (jcoord*ny)-1 , nz-1 /), xvalue_vector_hypre, ierr_hypre)

       index_count_hypre = 1
!       sol_mean_hypre    = 0.0
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
               sol(i,j,k) = xvalue_vector_hypre(index_count_hypre)
               index_count_hypre = index_count_hypre + 1
             end do
          end do
       end do

      call HYPRE_StructHybridGetNumIterati( solver_hypre, nIter, ierr_hypre )
      call HYPRE_StructHybridGetFinalRelat( solver_hypre, res, ierr_hypre )

!          if (giveInfo)  then 
          if ( giveInfo .and. master )  then 
             print*,""
             print*," HYPRE solver: "
             write(*,fmt="(a25,i25)") " Nb.of iterations: j = ", nIter
             write(*,fmt="(a25,es17.6)") " Final residual: res = ", res
             print*,""
          end if

! Free memory

  call HYPRE_StructGridDestroy(grid_hypre, ierr_hypre)
  call HYPRE_StructStencilDestroy(stencil_hypre, ierr_hypre)
  call HYPRE_StructMatrixDestroy(A_hypre, ierr_hypre)
  call HYPRE_StructVectorDestroy(b_hypre, ierr_hypre)
  call HYPRE_StructVectorDestroy(x_hypre, ierr_hypre)
  call HYPRE_StructHybridDestroy(solver_hypre, ierr_hypre)


    ! deallocate local fields
    deallocate(values_hypre, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(bvalue_vector_hypre, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"
    deallocate(xvalue_vector_hypre, stat=allocstat); if(allocstat/=0) &
         & stop "algebra.f90/bicgstab:dealloc failed"


          return

  end subroutine hypre

  ! modified by Junhong Wei (2016/07/08) (finishing line)


  !---------------------------------------------------------------------------


  subroutine poissonSolver_csr( var, dt, opt )
    ! -------------------------------------------------
    ! solves the Poisson problem with A stored in the
    ! compressed sparse row (csr) format
    ! -------------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         &intent(inout) :: var
    real, intent(in) :: dt
    character(len=*), intent(in)    :: opt
    

    ! local fields
    real, dimension(nnz) :: A_csr       ! compressed sparse row matrix (CSR)
    ! L and U combined in one matrix in CSR format:
    real, dimension(:), allocatable :: LU_csr  
    character(len=2), dimension(nnz) :: TM_csr ! test matrix with CSR
    integer, dimension(nnz) :: colInd          ! column index 
    integer, dimension(nxyz+1) :: rowPtr       ! points at new rows
    integer, dimension(nxyz) :: diagPtr        ! points at diagonal element
    real, dimension(nxyz) :: b

    ! local auxiliary fields
    real, dimension(-5:5) :: diag      ! aux field for A: collects elements 
    ! of a row out of CSR order -> needs reordering
    ! aux field for TM -> needs reordering:
    character(len=2), dimension(-5:5) :: TMdiag  
    integer, dimension(-5:5) :: auxColInd  ! aux colInd -> needs reordering

    ! local variables
    integer :: i,j,k, l,m, ii, csrInd
    real :: uR,uL, vF,vB, wU,wD
    real :: pStratU, pStratD, rhoEdge
    real :: AL,AR, AB,AF, AD,AU

    ! diag preconditioner
    real :: A_ii

    ! test bicgstab
    real, dimension(nxyz) :: sol
    real :: res
    real :: divSum

    ! consistent divergence
    real :: fL, fR, gF, gB, hU, hD
    real :: thetaL, thetaR, thetaF, thetaB, thetaU, thetaD
    real :: uSurf, vSurf, wSurf
    real :: rhoL, rhoR, rhoF, rhoB, rhoU, rhoD

    ! verbose
    logical, parameter :: giveInfo = .false.

!   achatzb
    if(topography) stop 'poissonSolver_csr not ready for topography!'
!   achatze


    ! ----------------------------------
    !          Poisson Problem
    ! ----------------------------------

    if(verbose) print*,"update.f90/poissonSolver: &
         & Setting up Poisson problem."


    ! ------------------------------------
    !    setup b = Pu*  (right hand side)
    ! ------------------------------------

    b = 10.0
    divSum = 0.0

    do k = 1,nz
       do j  = 1,ny
          do i = 1,nx

             l = getIndex(i,j,k)     ! 1D index for b

             uR = var(i,j,k,2); uL = var(i-1,j,k,2)
             vF = var(i,j,k,3); vB = var(i,j-1,k,3)
             wU = var(i,j,k,4); wD = var(i,j,k-1,4)

             pStratU = 0.5*(pStrat(k+1) + pStrat(k  ))
             pStratD = 0.5*(pStrat(k  ) + pStrat(k-1))
             b(l) = pStrat(k) * ( (uR-uL)/dx + (vF-vB)/dy ) &
                  & + (pStratU*wU - pStratD*wD)/dz

             ! scale RHS with Ma^2
             b(l) = b(l)*Ma**2*kappa

             divSum = divSum + b(l)

          end do
       end do
    end do

    !--- give info
    if( giveInfo ) then
       print*,""
       print*," Poisson Solver_csr: Initial state"
       print*,""
       write(*,fmt="(a25,es17.6)") "Sum over RHS = ", divSum
       write(*,fmt="(a25,es17.6)") "max(divPu*)[Pa/s] = ", &
            & maxval(abs(b))*pref*uref/lref
    end if


    ! ---------------------------------
    !  setup A_csr: coefficient matrix 
    ! ---------------------------------

    ! init
    A_csr = 1000.0        ! all values should be overwritten 
    TM_csr = "--"         !  - " - 
    rowPtr = 0            
    rowPtr(1) = 1         ! first element starts first row (no zero row)
    diagPtr = 1000
    diagPtr(1) = 1        ! first element is diagonal element
    csrInd = 0            ! counter for elements in A_csr

    ! loop over equation index l -> over lines of A
    k_loop: do k = 1,nz
       j_loop: do j = 1,ny
          i_loop: do i = 1,nx
             l = getIndex(i,j,k)  ! row index of A

             ! init aux fields
             auxColInd = 0

             ! ------------------ A(i+1,j,k) ------------------
             if (i<nx) then
                m = getIndex(i+1,j,k)
                rhoEdge = 0.5*( var(i+1,j,k,1) + var(i,j,k,1) )
                AR = dt/dx**2 * pStrat(k)**2/rhoEdge
                diag(1) = AR
                if(testCase == "matrixStructureTest") TMdiag(1) = "AR"
                auxColInd(1) = m
             else ! i = nx  right boundary (periodic)
                ! nx+1 -> 1
                m = getIndex(1,j,k)
                !                rhoEdge = 0.5*( var(i+1,j,k,1) + var(i,j,k,1) )
                rhoEdge = 0.5*( var(1,j,k,1) + var(nx,j,k,1) )
                AR = dt/dx**2 * pStrat(k)**2/rhoEdge
                diag(-2) = AR
                if(testCase == "matrixStructureTest") TMdiag(-2) = "AR"
                auxColInd(-2) = m
             end if

             ! ------------------- A(i-1,j,k) --------------------
             if (i>1) then
                m = getIndex(i-1,j,k)
                rhoEdge = 0.5*( var(i,j,k,1) + var(i-1,j,k,1) )
                AL = dt/dx**2 * pStrat(k)**2 / rhoEdge
                diag(-1) = AL
                if(testCase == "matrixStructureTest") TMdiag(-1) = "AL"
                auxColInd(-1) = m
             else ! i = 1 left boundary (perdiodic)
                ! i-1 -> nx 
                m = getIndex(nx,j,k)
                !                rhoEdge = 0.5*( var(i,j,k,1) + var(i-1,j,k,1) )
                rhoEdge = 0.5*( var(1,j,k,1) + var(nx,j,k,1) )
                AL = dt/dx**2 * pStrat(k)**2 / rhoEdge
                diag(2) = AL
                if(testCase == "matrixStructureTest") TMdiag(2) = "AL"
                auxColInd(2) = m
             end if

             ! -------------------- A(i,j+1,k) ----------------------
             if (j<ny) then
                m = getIndex(i,j+1,k)
                rhoEdge = 0.5*( var(i,j+1,k,1) + var(i,j,k,1) )
                AF = dt/dy**2 * pStrat(k)**2 / rhoEdge
                diag(3) = AF
                if(testCase == "matrixStructureTest") TMdiag(3) = "AF"
                auxColInd(3) = m
             else ! j = ny -> forward boundary (periodic)
                ! j+1 -> 1
                m = getIndex(i,1,k)
                !                rhoEdge = 0.5*( var(i,j+1,k,1) + var(i,j,k,1) )
                rhoEdge = 0.5*( var(i,1,k,1) + var(i,ny,k,1) )
                AF = dt/dy**2 * pStrat(k)**2 / rhoEdge
                diag(-4) = AF
                if(testCase == "matrixStructureTest") TMdiag(-4) = "AF"
                auxColInd(-4) = m
             end if

             ! --------------------- A(i,j-1,k) -----------------------
             if (j>1) then
                m = getIndex(i,j-1,k)
                rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j-1,k,1) )
                AB = dt/dy**2 * pStrat(k)**2 / rhoEdge
                diag(-3) = AB
                if(testCase == "matrixStructureTest") TMdiag(-3) = "AB"
                auxColInd(-3) = m
             else ! j = 1 -> backward boundary (periodic)
                ! j-1 -> ny
                m = getIndex(i,ny,k)
                !                rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j-1,k,1) )
                rhoEdge = 0.5* ( var(i,1,k,1) + var(i,ny,k,1) )
                AB = dt/dy**2 * pStrat(k)**2 / rhoEdge
                diag(4) = AB
                if(testCase == "matrixStructureTest") TMdiag(4) = "AB"
                auxColInd(4) = m
             end if


             ! ------------------ vertical elements AU, AD ---------------


             ! ---------------------- A(i,j,k+1) ------------------------
             if (k<nz) then
                m = getIndex(i,j,k+1)
                rhoEdge = 0.5* ( var(i,j,k+1,1) + var(i,j,k,1) )
                pStratU = 0.5* ( pStrat(k+1) + pStrat(k) )
                AU = dt/dz**2 * pStratU**2 / rhoEdge
                diag(5) = AU
                if(testCase == "matrixStructureTest") TMdiag(5) = "AU"
                auxColInd(5) = m
             else ! k = nz -> upwad boundary (solid wall)
                ! A(i,j,nz+1) = 0
                AU = 0.0
             end if

             ! ----------------------- A(i,j,k-1) ------------------------
             if (k>1) then 
                m = getIndex(i,j,k-1)
                rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j,k-1,1) )
                pStratD = 0.5* ( pStrat(k) + pStrat(k-1) )
                AD = dt/dz**2 * pStratD**2 / rhoEdge
                diag(-5) = AD
                if(testCase == "matrixStructureTest") TMdiag(-5) = "AD"
                auxColInd(-5) = m
             else ! k = 1 -> downward boundary (solid wall)
                ! A(i,j,0) = 0
                AD = 0.0
             end if

             ! ----------------------- A(i,j,k) --------------------------
             m = getIndex(i,j,k)
             diag(0) = - AR - AL - AF - AB - AU - AD
             if(testCase == "matrixStructureTest") TMdiag(0) = "AP"
             auxColInd(0) = m



             ! ---------------------------------------
             !        reorder entries for A_CSR
             ! ---------------------------------------
             ! sequence in A_csr is re-established as in A

             do ii = -5,5
                if (auxColInd(ii) /= 0) then
                   csrInd = csrInd + 1
                   A_csr(csrInd) = diag(ii)
                   TM_csr(csrInd) = TMdiag(ii)
                   colInd(csrInd) = auxColInd(ii)
                   if( l==colInd(csrInd) ) diagPtr(l) = csrInd  
                   ! pointer for diagonal entry

                end if
             end do
             rowPtr(l+1) = csrInd + 1        ! new row starts at csrInd+1
          end do i_loop
       end do j_loop
    end do k_loop


    !--------------------------------
    !     Linear equation solver
    !--------------------------------
    
    if( opt == "initial" ) then
       tol = tolInitial
       print*,"poisson_csr: use tolInitial = ", tolInitial
    else
       tol = tolPoisson
    end if
    
    select case (preconditioner)

    case ('no')

       sol = 0.0
       print*,""
       print*," BiCGStab without Preconditioner "
       print*,""
       
       
       
       call bicgstab_csr(A_csr,colInd,rowPtr,diagPtr,b,sol,&
            & tol,res,maxIterPoisson)


    case( 'diag' )

       sol = 0.0
       print*,""
       print*," BiCGStab with Diagonal Preconditioner "
       print*,""

       do l = 1, nxyz
          A_ii = A_csr(diagPtr(l))
          do csrInd = rowPtr(l), rowPtr(l+1)-1
             A_csr(csrInd) = A_csr(csrInd) / A_ii
          end do
          b(l) = b(l) / A_ii
       end do
       call bicgstab_csr(A_csr,colInd,rowPtr,diagPtr,b,sol,&
            & tol,res,maxIterPoisson)


    case ('ilu')

       ! warning for precond. BiCGSTAB
       if( ny == 1 ) then
          print*,"WARNING: ILU preconditioner does not work for ny = 1"
       end if

       sol = 0.0
       print*,""; print*,"-------------- ilu-precond. bicgstab_csr ------------"

       call ilu_csr(A_csr, colInd, rowPtr, diagPtr, LU_csr)

       call pre_bicgstab_csr(A_csr,LU_csr,colInd,rowPtr,diagPtr,b,sol,&
            & tol,res,maxIterPoisson)

    case default

       stop "update.f90/poissonSolver_csr: preconditioner not definded."

    end select



    ! ------------------------------------------
    !   Reshape pressure correction: sol -> dp
    ! ------------------------------------------

    dp = 0.0
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx

             l = getIndex(i,j,k)
             dp(i,j,k) = sol(l)

          end do
       end do
    end do

    ! --------------------------
    !       fix p at a point
    ! --------------------------
    !    dp = dp - dp(1,1,1)


    ! --------------------------
    !       give info
    ! --------------------------
    if( giveInfo ) then
       print*," Momentum Corrector"
       print*,""
       write(*,fmt="(a25,es17.6)") " dpMax = ", MAXVAL(dp)
       write(*,fmt="(a25,es17.6)") " dpMin = ", MINVAL(dp)
    end if


    !---------------------------------
    !       matrix setup test
    !---------------------------------

    ! print Testmatrix on Screen
    if (testCase == "matrixStructureTest") then

       print*,"update.f90/poissonSolver_full: TestMatrix TM_csr = "
       write(*,fmt="(27a4)") TM_csr
       print*,"colInd = ", colInd
       print*,"rowPtr = ", rowPtr
       print*,"diagPtr = ", diagPtr
    end if

    ! print A on screen in matrix format
    if (testCase == "matrixStructureTest")  then
       print*,"update.f90/poissonSolver_csr: Matrix A_csr = "
       write(*,fmt="(10es10.1)") A_csr
    end if

  end subroutine poissonSolver_csr


  !-------------------------------------------------------------------------


  subroutine pressureBoundaryCondition
    !--------------------------------------------------
    ! set pressure correction dp in ghost cells for BC
    !--------------------------------------------------

    ! modified by Junhong Wei (20161107) *** starting line ***

    ! auxiliary fields for "dp"
    real, dimension(ny,nz) :: xSliceLeft_send, xSliceRight_send
    real, dimension(ny,nz) :: xSliceLeft_recv, xSliceRight_recv
    
    real, dimension(nx,nz) :: ySliceBack_send, ySliceForw_send
    real, dimension(nx,nz) :: ySliceBack_recv, ySliceForw_recv
    
    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount


    ! Find neighbour procs
    if ( idim > 1 ) call mpi_cart_shift(comm,0,1,left,right,ierror)
    if ( jdim > 1 ) call mpi_cart_shift(comm,1,1,back,forw,ierror)

!xxxx
    if( giveInfo .and. master ) then
       print*,""
       print*,"----------------------------------------------"
       print*,"pressureBoundaryCondition: setting dp in Halos..."
       print*,"----------------------------------------------"
       print*,""
    end if

    !----------------------------
    !   set Halo cells: xSlice
    !----------------------------

    ! slice size
    sendcount = ny*nz
    recvcount = sendcount

    ! read slice into contiguous array
    xSliceLeft_send (:,:) = dp(1, 1:ny,1:nz)
    xSliceRight_send(:,:) = dp(nx,1:ny,1:nz)

    if( idim > 1 ) then

       ! left -> right
       source = left
       dest = right
       tag = 100

       call mpi_sendrecv(xSliceRight_send(1,1), sendcount, &
            & mpi_double_precision, dest, tag, &
            & xSliceLeft_recv(1,1), recvcount, mpi_double_precision, &
            & source, mpi_any_tag, comm, sts_left, ierror)

       ! right -> left
       source = right
       dest = left
       tag = 100

       call mpi_sendrecv(xSliceLeft_send(1,1), sendcount, &
            & mpi_double_precision, dest, tag, &
            & xSliceRight_recv(1,1), recvcount, mpi_double_precision, &
            & source, mpi_any_tag, comm, sts_right, ierror)

       ! right halos
       dp(nx+1,1:ny,1:nz) = xSliceRight_recv(1:ny,1:nz)

       ! left halos
       dp(0,1:ny,1:nz) = xSliceLeft_recv(1:ny,1:nz)

    else

       !if( master ) then
       !   print*,"idim = 1!  Routine setHelos not checked for this case. Stop."
       !   call mpi_barrier(comm,ierror)
       !   call mpi_finalize(ierror)
       !   stop
       !end if

       dp(0,:,:) = dp(nx,:,:)
       dp(nx+1,:,:) = dp(1,:,:)

    end if
    
    if(verbose .and. master) print*,"horizontalHalos: &
         & x-horizontal halos copied."

    !------------------------------
    !   set Halo cells: ySlice 
    !------------------------------

    ! slice size
    sendcount = nx*nz
    recvcount = sendcount

    ! read slice into contiguous array
    ySliceBack_send (:,:) = dp(1:nx, 1,1:nz)
    ySliceForw_send(:,:) = dp(1:nx,ny,1:nz)

    if( jdim > 1 ) then

       ! back -> forw
       source = back
       dest = forw
       tag = 100

       call mpi_sendrecv(ySliceForw_send(1,1), sendcount, &
            & mpi_double_precision, dest, tag, &
            & ySliceBack_recv(1,1), recvcount, mpi_double_precision, &
            & source, mpi_any_tag, comm, sts_back, ierror)

       ! forw -> back
       source = forw
       dest = back
       tag = 100

       call mpi_sendrecv(ySliceBack_send(1,1), sendcount, &
            & mpi_double_precision, dest, tag, &
            & ySliceForw_recv(1,1), recvcount, mpi_double_precision, &
            & source, mpi_any_tag, comm, sts_right, ierror)

       ! forward halos
       dp(1:nx,ny+1,1:nz) = ySliceForw_recv(1:nx,1:nz)

       ! backward halos
       dp(1:nx,0,1:nz) = ySliceBack_recv(1:nx,1:nz)

    else

       !if( master ) then
       !   print*,"jdim = 1!  Routine setHelos not checked for this case. Stop."
       !   call mpi_barrier(comm,ierror)
       !   call mpi_finalize(ierror)
       !   stop
       !end if

       dp(:,0,:) = dp(:,ny,:)
       dp(:,ny+1,:) = dp(:,1,:)

    end if

    if(verbose .and. master) print*,"horizontalHalos: &
         & x-horizontal halos copied."


    !------------------------------------
    ! SERIAL version:
    !
    !     
    !----------------
    !  x-Boundary
    !----------------
    !    dp(0,:,:) = dp(nx,:,:)
    !    dp(nx+1,:,:) = dp(1,:,:)

    !----------------    
    !   y-Boundary
    !----------------
    !    dp(:,0,:) = dp(:,ny,:)
    !    dp(:,ny+1,:) = dp(:,1,:)

!    !----------------
!    !  x-Boundary
!    !----------------
!    select case( xBoundary ) 
!
!    case( "periodic" ) 
!       dp(0,:,:) = dp(nx,:,:)
!       dp(nx+1,:,:) = dp(1,:,:)
!
!    case default
!       stop "pressureBoundaryCondition: unknown case xBoundary."
!    end select
!
!
!    !----------------    
!    !   y-Boundary
!    !----------------
!    select case( yBoundary ) 
!
!    case( "periodic" ) 
!       dp(:,0,:) = dp(:,ny,:)
!       dp(:,ny+1,:) = dp(:,1,:)
!    case default
!       stop "pressureBoundaryCondition: unknown case yBoundary."
!    end select

        ! modified by Junhong Wei (20161107) *** finishing line ***

    !----------------
    !   z-Boundary
    !----------------
    select case( zBoundary )

    case( "periodic" )
       dp(:,:,0) = dp(:,:,nz)
       dp(:,:,nz+1) = dp(:,:,1)

    case( "solid_wall" )
       ! zero gradient   
       dp(:,:,0) = dp(:,:,1)
       dp(:,:,nz+1) = dp(:,:,nz)
    case default
       stop "pressureBoundaryCondition: unknown case yBoundary."
    end select


  end subroutine pressureBoundaryCondition


  !-------------------------------------------------------------------------


  subroutine correctorStep( var, dMom, dt, RKstage)
    !-------------------------------------------------
    !         correct pressure & velocity
    !-------------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         &intent(inout) :: var
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,3), &
         & intent(inout) :: dMom
    real, intent(in) :: dt
    integer, intent(in) :: RKstage

    ! local variables 
    integer :: i,j,k
    real :: pEdge, rhoEdge
    real :: pGradX, pGradY, pGradZ
    real :: du, dv, dw

    ! divergence check
    real :: maxDivPu, divPu
    real :: uR, uL, vF, vB, wU, wD
    real :: pStratU, pStratD 

    ! verbose
    logical, parameter :: giveInfo = .true.

!   achatzb
    integer :: i0,j0
!   achatze

    ! modified by Junhong Wei (20161108) *** starting line ***

        if( giveInfo .and. master ) then
       print*,""
       print*,"----------------------------------------------"
       print*,"correctorStep: correcting p, u, v, and w..."
       print*,"----------------------------------------------"
       print*,""
    end if
       

!   achatzb
    i0=is+nbx-1
    j0=js+nby-1
!   achatze

    ! --------------------------------------
    !             calc p + dp
    ! --------------------------------------

!   achatzb
!   var(0:nx+1,0:ny+1,1:nz,5) &
!   = var(0:nx+1,0:ny+1,1:nz,5) + dp(0:nx+1,0:ny+1,1:nz)
    ! pressure correction only for atmosphere cells
    if(topography) then
       do k = 0,nz+1
          do j = 0,ny+1
             do i = 0,nx+1
                if(.not.topography_mask(i0+i,j0+j,k)) then
                   var(i,j,k,5) = var(i,j,k,5) + dp(i,j,k)
                end if
             end do
          end do
       end do
      else
       var(0:nx+1,0:ny+1,1:nz,5) &
       = var(0:nx+1,0:ny+1,1:nz,5) + dp(0:nx+1,0:ny+1,1:nz)
    end if
!   achatze

!    ! --------------------------------------
!    !             calc p + dp
!    ! --------------------------------------
!
!    var(1:nx,1:ny,1:nz,5) = var(1:nx,1:ny,1:nz,5) + dp(1:nx,1:ny,1:nz)

    ! modified by Junhong Wei (20161108) *** finishing line ***

    ! --------------------------------------
    !           calc du and u + du
    ! --------------------------------------

    do k = 1,nz
       do j = 1,ny
          do i = 0,nx

             rhoEdge = 0.5 * ( var(i,j,k,1) + var(i+1,j,k,1) )
             if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)


             pGradX = ( dp(i+1,j,k) - dp(i,j,k) ) / dx
             du = -dt * kappaInv * MaInv2 * pStrat(k) / rhoEdge * pGradX     

             var(i,j,k,2) = var(i,j,k,2) + du

             ! correct x-momentum tendency
             dMom(i,j,k,1) = dMom(i,j,k,1) + rhoEdge*du/beta(RKstage)

          end do
       end do
    end do

    !--------------------------------------
    !         calc dv and v + dv
    !--------------------------------------

    do k = 1,nz
       do j = 0,ny
          do i = 1,nx
             rhoEdge = 0.5 * ( var(i,j,k,1) + var(i,j+1,k,1) )
             if( fluctuationMode ) rhoEdge = rhoEdge + rhoStrat(k)
             
             pGradY = ( dp(i,j+1,k) - dp(i,j,k) ) / dy
             dv = -dt * kappaInv * MaInv2 * pStrat(k) / rhoEdge * pGradY   

             var(i,j,k,3) = var(i,j,k,3) + dv

             ! correct y-momentum tendency
             dMom(i,j,k,2) = dMom(i,j,k,2) + rhoEdge*dv/beta(RKstage)

          end do
       end do
    end do


    !--------------------------------------
    !         calc w and  w + dw
    !--------------------------------------

    ! solid wall 
    ! do k = 0, nz                 !  boundary values remain unchanged 
    do k = 1, nz-1        
       do j = 1,ny
          do i = 1,nx

             
             if( fluctuationMode ) then
                pEdge = pStratTilde(k)
             else
                pEdge = 0.5 * ( pStrat(k) + pStrat(k+1) )
             end if

             
             rhoEdge = 0.5 * ( var(i,j,k,1) + var(i,j,k+1,1) )
             if( fluctuationMode ) rhoEdge = rhoEdge + rhoStratTilde(k)
             
             pGradZ = ( dp(i,j,k+1) - dp(i,j,k) ) / dz

             dw = -dt * kappaInv * MaInv2 * pEdge / rhoEdge * pGradz   

             var(i,j,k,4) = var(i,j,k,4) + dw

             ! correct z-momentum tendency
             dMom(i,j,k,3) = dMom(i,j,k,3) + rhoEdge*dw/beta(RKstage)


          end do
       end do
    end do


    ! periodic in z
    if( zBoundary == "periodic" ) then

       do k = 0,nz,nz
          do j = 1,ny
             do i = 1,nx


                if( fluctuationMode ) then
                   pEdge = pStratTilde(k)
                else
                   pEdge = 0.5 * ( pStrat(k) + pStrat(k+1) )
                end if

                rhoEdge = 0.5 * ( var(i,j,k,1) + var(i,j,k+1,1) )
                if( fluctuationMode ) rhoEdge = rhoEdge + rhoStratTilde(k)
                
                pGradZ = ( dp(i,j,k+1) - dp(i,j,k) ) / dz

                dw = -dt * kappaInv * MaInv2 * pEdge / rhoEdge * pGradz   

                var(i,j,k,4) = var(i,j,k,4) + dw

                ! correct z-momentum tendency
                dMom(i,j,k,3) = dMom(i,j,k,3) + rhoEdge*dw/beta(RKstage)
                
             end do
          end do
       end do

    end if

!   achatzb
!   -------------------------------------
!   in case of topography, 
!   set all velocities normal to the topographic surface to zero
!   set densities in land cesll to background density
!   ------------------------------------

    if(topography) then
       do k = 0, nz+1
          do j = 0, ny+1
             do i = 0, nx+1
!               u at x interfaces
                if(&
                   topography_mask(i0+i,j0+j,k)&
                   .or.&
                   topography_mask(i0+i+1,j0+j,k)&
                ) then
                   var(i,j,k,2)=0.
                end if

!               v at y interfaces
                if(&
                   topography_mask(i0+i,j0+j,k)&
                   .or.&
                   topography_mask(i0+i,j0+j+1,k)&
                ) then
                   var(i,j,k,3)=0.
                end if

!               w at z interfaces
                if(&
                   topography_mask(i0+i,j0+j,k)&
                   .or.&
                   topography_mask(i0+i,j0+j,k+1)&
                ) then
                   var(i,j,k,4)=0.
                end if

!               density in land cells
                if(topography_mask(i0+i,j0+j,k)) then
                   if( fluctuationMode ) then
                      var(i,j,k,1) = 0.0
                     else
                      var(i,j,k,1) = rhoStrat(k) 
                   end if
                end if

             end do
          end do
       end do
    end if
!   -------------------------------------
!   achatze



  end subroutine correctorStep


  !-------------------------------------------------------------------------


  function getIndex(i,j,k)
    !----------------------------------------------------
    ! compute matrix index l from spatial indicies i,j,k
    !----------------------------------------------------

    ! in/out variables
    integer :: getIndex 
    integer, intent(in) :: i,j,k 

    getIndex = (k-1)*nx*ny + (j-1)*nx + i
  end function getIndex


  !-------------------------------------------------------------------------


  subroutine init_poisson
    !-----------------------------------
    ! allocate poisson module variables
    !-----------------------------------

    ! local variables
    integer :: allocstat

    ! allocate fields
    allocate( dp(0:nx+1,0:ny+1,0:nz+1), stat=allocstat)
    if( allocstat /= 0) stop "init_poisson: alloc failed"

    allocate( sol_old1(1:nx,1:ny,1:nz), stat=allocstat)
    if( allocstat /= 0) stop "init_poisson: alloc failed"

    allocate( sol_old2(1:nx,1:ny,1:nz), stat=allocstat)
    if( allocstat /= 0) stop "init_poisson: alloc failed"

    allocate( p_pred(1:nx,1:ny,1:nz), stat=allocstat)
    if( allocstat /= 0) stop "init_poisson: alloc failed"


  end subroutine init_poisson


  !---------------------------------------------------------------------------


  subroutine terminate_poisson
    !-----------------------------------
    ! deallocate poisson module variables
    !-----------------------------------

    ! local variables
    integer :: allocstat

    ! deallocate fields
    deallocate( dp, stat=allocstat)
    if( allocstat /= 0) stop "init_poisson: dealloc failed"

    deallocate( sol_old1, stat=allocstat)
    if( allocstat /= 0) stop "init_poisson: dealloc failed"

    deallocate( sol_old2, stat=allocstat)
    if( allocstat /= 0) stop "init_poisson: dealloc failed"

    deallocate( p_pred, stat=allocstat)
    if( allocstat /= 0) stop "init_poisson: dealloc failed"

  end subroutine terminate_poisson

end module poisson_module

