module poisson_module

  use type_module
  use mpi_module   ! modified by Junhong Wei (20161106)
  use timeScheme_module
  use atmosphere_module
  use algebra_module
  use hypretools_module
  use output_module

  implicit none

  private
  ! all module objects (variables, functions, subroutines)
  ! are internal to the module by default


  !------------------------
  !   public subroutines
  !------------------------
  public :: Corrector
  public :: init_poisson
  public :: terminate_poisson
  public :: calculate_heating
  public :: heat_w0

  !------------------------
  !   private subroutines       
  !------------------------
  private :: getIndex
!  private :: poissonSolver_csr
  private :: pressureBoundaryCondition
  private :: correctorStep
  private :: linOpr
  private :: linOprXYZ
  private :: bicgstab
  private :: hypre                     ! modified by Junhong Wei
  private :: poissonSolver
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

  !UAB
  !subroutine Corrector( var,flux,dMom,dt,errFlagBicg,nIter,m,opt,w_0)
  subroutine Corrector( var,flux,dMom,dt,errFlagBicg,nIter,m,opt)
  !UAE
    ! -------------------------------------------------
    !              correct uStar, bStar, and p
    ! -------------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), intent(in) :: flux
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,3), &
         & intent(inout) :: dMom

    real, intent(in)                :: dt
    logical, intent(out)            :: errFlagBicg
    integer, intent(out)            :: nIter
    integer, intent(in)             :: m

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction 
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction 
    ! of the winds and density fluctuations
    character(len=*), intent(in)    :: opt

    !UAB
    !! w_0 = horizontal-mean vertical wind induced by heating
    !real, dimension(-nbz:nz+nbz), intent(in) :: w_0
    !UAE

    ! local variables
    real, dimension(1:nx,1:ny,1:nz) :: rhs        ! RHS
    logical :: onlyinfo 

    ! Note: dp is a module variable
    ! calc dp 

    ! Calc RHS of Poisson problem
    onlyinfo = .false.

    !UAB
    !call calc_RHS(rhs, var, flux, dt, onlyinfo, w_0)
    call calc_RHS(rhs, var, flux, dt, onlyinfo)
    !UAE

    call poissonSolver( rhs, var, dt, errFlagBicg, nIter, m, opt )

    ! set horizontal and vertical BC for dp
    call pressureBoundaryCondition

    ! correct p, rhopStar, and uStar with dp
    call correctorStep (var, dMom, dt, m, opt)
    
    !UAB
    !if (detailedinfo) call calc_RHS( rhs,var,flux,dt,detailedinfo,w_0 )
    if (detailedinfo) call calc_RHS( rhs,var,flux,dt,detailedinfo)
    !UAE

  end subroutine Corrector

  !----------------------------------------------------------------------

  subroutine preCond( sIn, sOut)
    ! --------------------------------------
    !   preconditioner for BiCGStab
    !   solves vertical problem exploiting its tri-diagonal character
    !   (Isaacson & Keller 1966, see also Durran's book appendix A.2)
    ! --------------------------------------

    ! in/out variables
    real, dimension(1:nx,1:ny,1:nz), intent(out) :: sOut
    real, dimension(1:nx,1:ny,1:nz), intent(in)  :: sIn

    ! local field
    real, dimension(1:nx,1:ny,1:nz) :: s_pc, q_pc
    real, dimension(1:nx,1:ny) :: p_pc

    ! local variables
    integer :: k
    integer :: i,j

    ! work with auxiliary field s_pc

    s_pc = sIn

    ! upward sweep

    do j=1,ny
       do i=1,nx
          au_b(i,j,nz) = 0.0
       end do
    end do

    do j=1,ny
       do i=1,nx
          q_pc(i,j,1) = - au_b(i,j,1)/ac_b(i,j,1)
          s_pc(i,j,1) =   s_pc(i,j,1)/ac_b(i,j,1)
       end do
    end do

    do k=2,nz
       do j=1,ny
          do i=1,nx
             p_pc(i,j) = 1.0/(ac_b(i,j,k) + ad_b(i,j,k)*q_pc(i,j,k-1))

             q_pc(i,j,k) = - au_b(i,j,k) * p_pc(i,j)

             s_pc(i,j,k) &
             = (s_pc(i,j,k) - ad_b(i,j,k)*s_pc(i,j,k-1)) * p_pc(i,j)
          end do
       end do
    end do

    ! backward pass

    do k=nz-1,1,-1
       do j=1,ny
          do i=1,nx
             s_pc(i,j,k) = s_pc(i,j,k) + q_pc(i,j,k)*s_pc(i,j,k+1)
          end do
       end do
    end do

    ! final result

    sOut = s_pc
    
    return

  end subroutine preCond


  !----------------------------------------------------------------------

  subroutine linOpr( sIn, Ls, opt )
    ! --------------------------------------
    !   Linear Operator in Poisson problem
    !   Functions as A*x
    ! --------------------------------------

    ! in/out variables
    real, dimension(1:nx,1:ny,1:nz), intent(out) :: Ls
    real, dimension(1:nx,1:ny,1:nz), intent(in)  :: sIn

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction 
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction 
    ! of the winds and density fluctuations
    character(len=*), intent(in)    :: opt

    ! local field (extended by ghost cells)
    real, dimension(0:nx+1,0:ny+1,0:nz+1) :: s

    ! auxiliary fields for "dp"
    real, dimension(0:ny+1,0:nz+1) :: xSliceLeft_send, xSliceRight_send
    real, dimension(0:ny+1,0:nz+1) :: xSliceLeft_recv, xSliceRight_recv

    real, dimension(0:nx+1,0:nz+1) :: ySliceBack_send, ySliceForw_send
    real, dimension(0:nx+1,0:nz+1) :: ySliceBack_recv, ySliceForw_recv

    ! local variables
    integer :: i,j,k
    real :: AL,AR, AB,AF, AD,AU, AC, ALB,ALF, ARB,ARF
    real :: sL,sR, sB,sF, sD,sU, sC, sLB,sLF, sRB,sRF

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! work with auxiliary field s
    s(1:nx,1:ny,1:nz) = sIn

    ! Find neighbour procs
    if (idim > 1) call mpi_cart_shift(comm,0,1,left,right,ierror)
    if (jdim > 1) call mpi_cart_shift(comm,1,1,back,forw,ierror)

    select case( model ) 

       !----------------------------------------
       !       Pseudo-incompressible model
       !----------------------------------------
    case( "pseudo_incompressible" )

       !----------------------------
       !   set Halo cells: xSlice
       !----------------------------

       if( xBoundary == "periodic" ) then
         if ( idim > 1 ) then
            ! slice size
            sendcount = (ny+2)*(nz+2)
            recvcount = sendcount

            ! read slice into contiguous array
            xSliceLeft_send (:,:) = s(1, :,:)
            xSliceRight_send(:,:) = s(nx,:,:)


            ! left -> right
            source = left
            dest = right
            tag = 100

            call mpi_sendrecv(xSliceRight_send(0,0), sendcount, &
                 & mpi_double_precision, dest, tag, &
                 & xSliceLeft_recv(0,0), recvcount, mpi_double_precision, &
                 & source, mpi_any_tag, comm, sts_left, ierror)

            ! right -> left
            source = right
            dest = left
            tag = 100

            call mpi_sendrecv(xSliceLeft_send(0,0), sendcount, &
                 & mpi_double_precision, dest, tag, &
                 & xSliceRight_recv(0,0), recvcount, &
                 & mpi_double_precision, &
                 & source, mpi_any_tag, comm, sts_right, ierror)

            ! right halos
            s(nx+1,:,:) = xSliceRight_recv(:,:)

            ! left halos
            s(0   ,:,:) =  xSliceLeft_recv(:,:)
         else
           s(0   ,:,:) = s(nx,:,:)
           s(nx+1,:,:) = s(1 ,:,:)
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
         if( jdim > 1 ) then
            ! slice size
            sendcount = (nx+2)*(nz+2)
            recvcount = sendcount

            ! read slice into contiguous array
            ySliceBack_send(:,:) = s(:, 1,:)
            ySliceForw_send(:,:) = s(:,ny,:)


            ! back -> forw
            source = back
            dest = forw
            tag = 100

            call mpi_sendrecv(ySliceForw_send(0,0), sendcount, &
                 & mpi_double_precision, dest, tag, &
                 & ySliceBack_recv(0,0), recvcount, mpi_double_precision, &
                 & source, mpi_any_tag, comm, sts_back, ierror)

            ! forw -> back
            source = forw
            dest = back
            tag = 100

            call mpi_sendrecv(ySliceBack_send(0,0), sendcount, &
                 & mpi_double_precision, dest, tag, &
                 & ySliceForw_recv(0,0), recvcount, mpi_double_precision, &
                 & source, mpi_any_tag, comm, sts_right, ierror)

            ! forward halos
            s(:,ny+1,:) = ySliceForw_recv(:,:)

            ! backward halos
            s(:,0   ,:) = ySliceBack_recv(:,:)
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

                AR = ar_b(i,j,k)
                sR = s(i+1,j,k)

                ! ------------------- A(i-1,j,k) --------------------

                AL = al_b(i,j,k)
                sL = s(i-1,j,k)

                ! -------------------- A(i,j+1,k) ----------------------

                AF = af_b(i,j,k)
                sF = s(i,j+1,k)

                ! --------------------- A(i,j-1,k) -----------------------

                AB = ab_b(i,j,k)
                sB = s(i,j-1,k)

                ! --------------------- A(i,j,k+1) ------------------------

                if (k<nz) then
                   AU = au_b(i,j,k)
                   sU = s(i,j,k+1)
                else ! k = nz -> upwad boundary (solid wall)
                   ! A(i,j,nz+1) = 0
                   AU = 0.0
                   sU = 0.0
                end if

                ! --------------------- A(i,j,k-1) ------------------------

                if (k>1) then 
                   AD = ad_b(i,j,k)
                   sD = s(i,j,k-1)
                else ! k = 1 -> downward boundary (solid wall)
                   ! A(i,j,0) = 0
                   AD = 0.0
                   sD = 0.0
                end if

                ! -------------------- A(i,j,k) --------------------------

                AC = ac_b(i,j,k)
                sC = s(i,j,k)


                ! -------------------- apply Operator ---------------------

                Ls(i,j,k) &
                = AL*sL + AR*sR + AF*sF + AB*sB + AU*sU + AD*sD + AC*sC

                if (timeScheme == "semiimplicit") then
                   if (opt == 'impl') then
                      ! -------------------- A(i,j,k) ---------------------

                      ALB = alb_b(i,j,k)
                      sLB = s(i-1,j-1,k)

                      ! -------------------- A(i,j,k) ---------------------

                      ALF = alf_b(i,j,k)
                      sLF = s(i-1,j+1,k)

                      ! -------------------- A(i,j,k) ---------------------

                      ARB = arb_b(i,j,k)
                      sRB = s(i+1,j-1,k)

                      ! -------------------- A(i,j,k) ---------------------

                      ARF = arf_b(i,j,k)
                      sRF = s(i+1,j+1,k)

                      Ls(i,j,k) &
                      = Ls(i,j,k) + ALB*sLB + ALF*sLF + ARB*sRB + ARF*sRF
                     else if (opt /= 'expl') then
                      stop'ERROR: linOpr expects opt = expl or opt = impl'
                   end if
                end if

                ! ---------------- scale with thetaStrat ------------------
                if( pressureScaling ) then
                   stop'ERROR: pressure scaling disabled'
                end if
             end do i_loop
          end do j_loop
       end do k_loop

    case( "Boussinesq" )
       !----------------------------------------
       !             Boussinesq model
       !----------------------------------------

       stop'ERROR: linOpr not ready for Boussinesq model'

    case default
       stop "linOpr: unknown case model"
    end select

  end subroutine linOpr


  !-----------------------------------------------------------------------


  !UAB
  !subroutine calc_RHS( b,var,flux,dt,onlyinfo,w_0 )
  subroutine calc_RHS( b,var,flux,dt,onlyinfo)
  !UAE

    !----------------------------------------
    !   calculates the RHS of the 
    !   Poisson problem
    !----------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         &intent(in) :: var 
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), intent(in) :: flux
    real, intent(in) :: dt
    real, dimension(1:nx,1:ny,1:nz), intent(out) :: b  ! RHS
    logical, intent(in) :: onlyinfo                    ! give info in div
    !UAB
    !real, dimension(-nbz:nz+nbz),intent(in) :: w_0 !w_0 due to heating
    !UAE
    
    ! local vars
    real :: uR,uL, vF,vB, wU,wD
    real :: pStratU, pStratD
    real, dimension(1:nz) :: sum_local, sum_global !UA

    integer :: i,j,k
    real :: div, divSum, divSumScaled

    real :: bu,bv,bw,bl2loc,divL2_norm,divL2_norm_local

    ! for some diagnostics ...
    real :: dPudx_norm_local, dPudx_norm, dPvdy_norm_local, dPvdy_norm
    real :: dPwdz_norm_local, dPwdz_norm, Q_norm_local, Q_norm
    real :: Q1_norm_local, Q1_norm, Q2_norm_local, Q2_norm

    ! check L2-norm of divergence
    real :: divL2, divMax

    ! verbose
    logical, parameter :: giveInfo = .true.

    ! MPI stuff
    real :: divL2_local, divSum_local
    integer :: root

    integer :: i0,j0

    real :: rho, the
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz) :: heat
    real, dimension(-nbz:nz+nbz) :: w_0 
    real, dimension(-nbz:nz+nbz) :: S_bar 

    if( giveInfo .and. master ) then
       print*,""
       print*,"----------------------------------------------"
       print*,"   calc_RHS: computing RHS... "
       print*,"----------------------------------------------"
       print*,""
    end if

    ! ----------------------------------
    !          Poisson Problem
    ! ----------------------------------

    if( master .and. verbose ) print*,"update.f90/poissonSolver: &
                                     & Setting up Poisson problem." 

    i0=is+nbx-1
    j0=js+nby-1

    ! heating due to thermal relaxation, molecular or turbulent diffusion,
    ! or gravity waves, and the horizontal-mean vertical wind due to it

    !UAB
    !heat(:, :, :) = 0.
    !call calculate_heating(var,flux,heat)

    !if (raytracer) heat(:,:,:) = heat(:,:,:) + var(:,:,:,8)

    ! GBcorr -> FS
    if (heatingONK14 .or. TurbScheme .or. rayTracer) then
    !if (heating) then
       call heat_w0(var,flux,heat,S_bar,w_0)
      else
       heat = 0.
       S_bar = 0.
       w_0 = 0.
    end if

    ! subtreact horizontal mean of heat(:,:,:)

    do i = 1,nx
       do j = 1,ny
          heat(i,j,1:nz) = heat(i,j,1:nz) - S_bar(1:nz)
       end do
    end do
    !UAE

    ! scale RHS with Ma^2 * kappa, hence ...
    heat(:,:,:) = heat(:,:,:) * Ma**2 * kappa
    
    !--------------------------------------------------
    !    setup b = Ma^2 * P * u^*  (right hand side)
    !--------------------------------------------------

    divSum = 0.0
    divSumScaled = 0.0
    divL2  = 0.0
    divMax = 0.0

    divSum_local = 0.0  
    divL2_local = 0.0

    divL2_norm = 0.0
    divL2_norm_local = 0.0

    if (RHS_diagnostics) then
       dPudx_norm_local = 0.0
       dPudx_norm = 0.0
    
       dPvdy_norm_local = 0.0
       dPvdy_norm = 0.0
    
       dPwdz_norm_local = 0.0
       dPwdz_norm = 0.0
    
       Q_norm_local = 0.0
       Q_norm = 0.0
    
       Q1_norm_local = 0.0
       Q1_norm = 0.0
    
       Q2_norm_local = 0.0
       Q2_norm = 0.0
    end if

    select case( model ) 

    case( "pseudo_incompressible" ) 

       do k = 1,nz
          do j  = 1,ny
             do i = 1,nx

                uR = var(i,j,k,2); uL = var(i-1,j,k,2)
                vF = var(i,j,k,3); vB = var(i,j-1,k,3)
                !UAB
                !wU = var(i,j,k,4); wD = var(i,j,k-1,4)
                wU = var(i,j,k,4) - w_0(k); wD = var(i,j,k-1,4) - w_0(k-1)
                !UAE

                PstratU = PstratTilde(k)
                PstratD = PstratTilde(k-1)
                
                ! scale RHS with Ma^2 * kappa, hence ...
                
                bu = Pstrat(k) * (uR-uL)/dx * Ma**2 * kappa
                bv = Pstrat(k) * (vF-vB)/dy * Ma**2 * kappa
                !UAB
                !if(heatingONK14)then
                !   bw &
                !   = (PstratU*(wU-w_0(k)) - PstratD*(wD-w_0(k-1)))/dz &
                !     * Ma**2 * kappa
                !  else
                !   (PstratU*wU - PstratD*wD)/dz * Ma**2 * kappa
                !end if
                bw = (PstratU*wU - PstratD*wD)/dz * Ma**2 * kappa
                !UAE

                b(i,j,k) = bu + bv + bw + heat(i,j,k)

                ! L2-norm of the divergence div(Pu)
                ! divL2 = divL2 + b(i,j,k)**2

                divL2_local = divL2_local + b(i,j,k)**2

                ! introduce a reference norm for the RHS
                ! b(i,j,k) =   Pstrat(k) * ( (uR-uL)/dx + (vF-vB)/dy ) &
                !          & + (PstratU*wU - PstratD*wD)/dz &
                !          & + heat

                bl2loc = bu**2 + bv**2 + bw**2 + (heat(i,j,k))**2
                divL2_norm_local = divL2_norm_local + bl2loc

                if (RHS_diagnostics) then
                   dPudx_norm_local = dPudx_norm_local + bu**2
                   dPvdy_norm_local = dPvdy_norm_local + bv**2
                   dPwdz_norm_local = dPwdz_norm_local + bw**2
                   Q_norm_local = Q_norm_local + (heat(i,j,k))**2
                end if

                ! max norm of div(Pu)
                if( abs(b(i,j,k)) > divMax ) then
                   divMax = abs(b(i,j,k))
                end if
                
                ! check sum for solvability criterion
                ! divSum = divSum + b(i,j,k)
                divSum_local = divSum_local + b(i,j,k)
                
                ! Skalierung mit thetaStrat
                if ( pressureScaling ) then
                   stop'ERROR: pressure scaling disabled'
                end if

             end do
          end do
       end do

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

       !MPI: sum divL2_norm_local over all procs
       root = 0
       call mpi_reduce(divL2_norm_local, divL2_norm, 1, &
            & mpi_double_precision, mpi_sum, root, comm, ierror)
       
       call mpi_bcast(divL2_norm, 1, &
            & mpi_double_precision, root, comm, ierror)

       if (RHS_diagnostics) then
          call mpi_allreduce(dPudx_norm_local, dPudx_norm, 1, &
               & mpi_double_precision, mpi_sum, comm, ierror)
          call mpi_allreduce(dPvdy_norm_local, dPvdy_norm, 1, &
               & mpi_double_precision, mpi_sum, comm, ierror)
          call mpi_allreduce(dPwdz_norm_local, dPwdz_norm, 1, &
               & mpi_double_precision, mpi_sum, comm, ierror)
          call mpi_allreduce(Q_norm_local, Q_norm, 1, &
               & mpi_double_precision, mpi_sum, comm, ierror)
            
          dPudx_norm = sqrt(dPudx_norm/sizeX/sizeY/sizeZ)            
          dPvdy_norm = sqrt(dPvdy_norm/sizeX/sizeY/sizeZ)
          dPwdz_norm = sqrt(dPwdz_norm/sizeX/sizeY/sizeZ)
          Q_norm = sqrt(Q_norm/sizeX/sizeY/sizeZ)
       end if

       ! scale div
       divL2_local = sqrt(divL2_local/nx/ny/nz)
       divL2 = sqrt(divL2/sizeX/sizeY/sizeZ)

       divL2_norm_local = sqrt(divL2_norm_local/nx/ny/nz)
       divL2_norm = sqrt(divL2_norm/sizeX/sizeY/sizeZ)
   
       alpha_tol = divL2_norm
       b_norm = divL2

       if (divL2_norm /= 0.0) then
          tolref = divL2/divL2_norm
         else
          if (divL2 == 0.0) then
             tolref = 1.0
            else
             stop'ERROR: divL2_norm = 0 while divL2 /= 0'
          end if
       end if

       if(master) then
           print*,"tolref = ",tolref

           if(RHS_diagnostics) then
              print*,"RHS diagnostics:"
              print*,"dPudx_norm = ",dPudx_norm
              print*,"dPvdy_norm = ",dPvdy_norm
              print*,"dPwdz_norm = ",dPwdz_norm
              print*,"Q_norm = ",Q_norm
              print*,"alpha = ",divL2_norm
              print*,"unscaled |b|=", divL2
           end if
       end if 

!      root=0
!      call mpi_bcast(tolref, 1, mpi_double_precision, root, comm, ierror)
!      call mpi_barrier(comm,ierror)

    case( "Boussinesq" )

       do k = 1,nz
          do j  = 1,ny
             do i = 1,nx

                uR = var(i,j,k,2); uL = var(i-1,j,k,2)
                vF = var(i,j,k,3); vB = var(i,j-1,k,3)
                wU = var(i,j,k,4); wD = var(i,j,k-1,4)

                bu = Ma**2 * kappa / theta00 * (uR-uL)/dx
                bv = Ma**2 * kappa / theta00 * (vF-vB)/dy
                bw = Ma**2 * kappa / theta00 * (wU-wD)/dz

                div = bu + bv + bw

!               introduce a reference norm for the RHS
!               div = (uR-uL)/dx + (vF-vB)/dy + (wU-wD)/dz
                bl2loc = bu**2 + bv**2 + bw**2

                if(topography) then
                   stop'ERROR: topography needs implicit time stepping &
                      & that is not ready yet for Boussinesq'
                end if

                ! check sum for solvability criterion (shoud be zero)
                ! divSum = divSum + b(i,j,k)
                divSum_local = divSum_local + b(i,j,k)
                
                ! divL2 = divL2 + div**2
                divL2_local = divL2_local + div**2

                divL2_norm_local = divL2_norm_local + bl2loc

                if( abs(div) > divMax ) divMax = abs(div)

             end do
          end do
       end do

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

       !MPI: sum divL2_norm_local over all procs
       root = 0
       call mpi_reduce(divL2_norm_local, divL2_norm, 1, &
            & mpi_double_precision, mpi_sum, root, comm, ierror)
       
       call mpi_bcast(divL2_norm, 1, &
            & mpi_double_precision, root, comm, ierror)

       ! scale div
       divL2_local = sqrt(divL2_local/nx/ny/nz)
       divL2 = sqrt(divL2/sizeX/sizeY/sizeZ)

       ! scale by number of cells 
       divL2_local = sqrt(divL2_local/nx/ny/nz)
       divL2 = sqrt(divL2/sizeX/sizeY/sizeZ)

       divL2_norm_local = sqrt(divL2_norm_local/nx/ny/nz)
       divL2_norm = sqrt(divL2_norm/sizeX/sizeY/sizeZ)

       alpha_tol = divL2_norm
       b_norm = divL2

       if (divL2_norm /= 0.0) then
          tolref = divL2/divL2_norm
         else
          if (divL2 == 0.0) then
             tolref = 1.0
            else
             stop'ERROR: divL2_norm = 0 while divL2 /= 0'
          end if
       end if


       if(master) print*,"tolref = ",tolref

!      root=0
!      call mpi_bcast(tolref, 1, mpi_double_precision, root, comm, ierror)
!      call mpi_barrier(comm,ierror)
!      achatze

    case default
       stop "poissonSolver: unknown case model."
    end select

    !-------------------------------
    !     Display info on screen
    !-------------------------------

    if( master .and. onlyinfo ) then
       ! Information on divergence

       print*,""
       print*," Poisson Solver ", trim(poissonSolverType), &
            & ": Final state"
       print*,""

       select case ( model ) 
       case( "Boussinesq" ) 
          
          write(*,fmt="(a25,es17.6)") "L2(div(u)) [1/s] = ", &
               & divL2*uRef/lRef

          write(*,fmt="(a25,es17.6)") "max(div(u)) [1/s] = ", &
               & divMax*uRef/lRef

          write(*,fmt="(a25,es17.6)") "rms terms (div(u)) [1/s] = ", &
               & divL2_norm*uRef/lRef

          write(*,fmt="(a25,es17.6)") "normalized L2(div(u)) = ", &
               & divL2/divL2_norm
         
       case( "pseudo_incompressible" )

          write(*,fmt="(a25,es17.6)") "L2(div(Pu)) [Pa/s] = ", &
               & divL2*rhoRef*thetaRef*uRef/lRef

          write(*,fmt="(a25,es17.6)") "max(div(Pu)) [Pa/s] = ", &
               & divMax*rhoRef*thetaRef*uRef/lRef
         
          write(*,fmt="(a25,es17.6)") "rms terms (div(Pu)) [Pa/s] = ", &
                   & divL2_norm*rhoRef*thetaRef*uRef/lRef

          write(*,fmt="(a25,es17.6)") "normalized L2(div(Pu)) = ", &
                   & tolref
       case default
          stop "calc_RHS: unkown case model"
       end select

       print*,""
       print*,"-------------------------------------------------"
       print*,""
      else
       if( master .and. giveInfo ) then
          print*,""
          print*," Poisson Solver ", trim(poissonSolverType), &
               & ": Initial state"
          print*,""
          write(*,fmt="(a25,es17.6)") "Sum over RHS = ", divSum
          write(*,fmt="(a25,es17.6)") "Sum over scaled RHS  = ", &
                  & divSumScaled

          ! Information on divergence
          select case ( model ) 
          case( "Boussinesq" ) 

             write(*,fmt="(a25,es17.6)") "L2(div(u)) [1/s] = ", &
                  & divL2*uRef/lRef

             write(*,fmt="(a25,es17.6)") "max(div(u)) [1/s] = ", &
                  & divMax*uRef/lRef

             write(*,fmt="(a25,es17.6)") "rms terms (div(u)) [1/s] = ", &
                  & divL2_norm*uRef/lRef

             write(*,fmt="(a25,es17.6)") "normalized L2(div(u)) = ", &
                  & divL2/divL2_norm
        
          case( "pseudo_incompressible" )

             write(*,fmt="(a25,es17.6)") "L2(div(Pu)) [Pa/s] = ", &
                  & divL2*rhoRef*thetaRef*uRef/lRef

             write(*,fmt="(a25,es17.6)") "max(div(Pu)) [Pa/s] = ", &
                  & divMax*rhoRef*thetaRef*uRef/lRef

             write(*,fmt="(a25,es17.6)") "rms terms (div(Pu)) &
                                       & [Pa/s] = ", &
                                       & divL2_norm*rhoRef*thetaRef &
                                       & *uRef/lRef
 
             write(*,fmt="(a25,es17.6)") "normalized L2(div(Pu)) = ", &
                  & tolref
        
          case default
             stop "calc_RHS: unkown case model"
          end select
       end if
    end if

  end subroutine calc_RHS

  !----------------------------------------------------------------------


  subroutine poissonSolver( b,var,dt,errFlagBicg,nIter,m,opt )
    ! -------------------------------------------------
    ! solves the Poisson problem with 
    ! application of linear operator L
    ! -------------------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         &intent(in) :: var
    real, intent(in) :: dt

    logical, intent(out) :: errFlagBicg
    integer, intent(out) :: nIter

    integer, intent(in) :: m
    real, dimension(1:nx,1:ny,1:nz),intent(in) :: b        ! RHS

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction 
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction 
    ! of the winds and density fluctuations
    character(len=*), intent(in)    :: opt

    ! local vars
    real, dimension(1:nx,1:ny,1:nz) :: sol ! solution of Poisson problem
    real :: res
    
    real :: dtInv

    ! verbose
    logical, parameter :: giveInfo = .true.

    ! Init
    if (dt == 0.0) stop "poissonSolver: dt = 0.0. Stopping."
    dtInv = 1.0/dt

    !--------------------------------
    !     Linear equation solver
    !     solve for dt * dp ...
    !--------------------------------

    sol = 0.0

    select case( poissonSolverType ) 

    ! bicgstab solver
    case( "bicgstab" )

       select case( model ) 

       case( "pseudo_incompressible" )

        call val_PsIn(var, dt, opt)

       case( "Boussinesq" )

        stop'ERROR: BiCGStab still to be made ready for Boussinesq'

       case default
          stop "linOpr: unknown case model"
       end select

       call bicgstab(b, dt, sol, res, nIter, errFlagBicg, opt)

    case( "gcr" ) 

       stop'ERROR: no gcr provided anymore'

    case( "adi" ) 

       stop'ERROR: no adi provided anymore'

    ! hypre solver
    case( "hypre" )    

       select case( model ) 

       case( "pseudo_incompressible" )

        call val_PsIn(var, dt, opt)

       case( "Boussinesq" )

        call val_hypre_Bous

       case default
          stop "linOpr: unknown case model"
       end select

       call hypre(b, dt, sol, res, nIter, errFlagBicg, opt)
  
       case default
          stop "Unknown PoissonSolver. Stop"
    end select

    ! now get dp from dt * dp ...

    dp(1:nx,1:ny,1:nz) = dtInv * sol ! pass solution to pressure corrector

  end subroutine poissonSolver

  !----------------------------------------------------------------------

  subroutine linOprXYZ( var, q, Lq, direction )
    ! --------------------------------------
    !   Linear Operator in Poisson problem
    !   Functions as A*x
    ! --------------------------------------

    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz, nVar), &
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


  !----------------------------------------------------------------------


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

  
  !----------------------------------------------------------------------

  subroutine bicgstab(b, dt, sol, res, nIter, errFlag, opt)
    ! --------------------------------------
    !    BiCGStab using linear operator
    !    preconditioner applied via A M^-1 M x = b
    !---------------------------------------

    ! in/out variables
    real, dimension(1:nx,1:ny,1:nz), intent(in) :: b        ! RHS 
    real, intent(in) :: dt
    real, dimension(1:nx,1:ny,1:nz), intent(inout) :: sol
    real, intent(out) :: res                     ! residual
    integer, intent(out) :: nIter
    logical, intent(out) :: errFlag

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction 
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction 
    ! of the winds and density fluctuations
    character(len=*), intent(in)    :: opt

    ! Local parameters
    integer  :: maxIt

    ! local variables
    integer :: i,j,k, allocstat
    integer :: j_b
    real, dimension(:,:,:), allocatable :: p,r0,rOld,r,s,t,v, matVec, v_pc
    real :: alpha, beta, omega

    ! verbose
    logical, parameter :: giveInfo = .true.

    ! MPI stuff
    integer :: root
    real :: res_local

    if( giveInfo .and. master ) then
       print*,""
       print*,"----------------------------------------------"
       print*,"BICGSTAB: solving linear system... "
       print*,"----------------------------------------------"
       print*,""
    end if

    sol = 0.0

    ! Set parameters
    maxIt = maxIterPoisson

    !if( opt == 'initial' ) then
    !   print*,"bicgstab: use tolInitial = ", tolInitial
    !   tol = tolInitial
    !else
    !   tol = tolPoisson
    !end if

    ! modified convergence criterion so that iterations stop when either
    ! (a) tolcrit = abs  =>  |Ax - b| < eps b_*
    !     with b_* a suitable norm deciding whether b (= the divergence 
    !     criterion for the winds from the predictor) is small or not
    ! (b) tolcrit = rel  =>  |Ax - b| < eps |b|
    ! here eps = tolPoisson is the user-set convergence criterion
    ! hypre has the criterion |Ax - b| < tol * |b|, hence, with
    ! tolref = divL2/divL2_norm = |b|/b_*

    if (tolcrit == "abs") then
       tol = tolPoisson/tolref
      else if (tolcrit == "rel") then
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
    allocate(v_pc(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:alloc failed"

    ! error flag
    errFlag = .false.    

    ! Init
    ! r0 = b - Ax
    call linOpr( sol, matVec, opt )
    r0 = b - matVec    
    p = r0
    r = r0

    res_local = 0.0
    do k=1,nz
       do j=1,ny
          do i=1,nx
             res_local = res_local + r(i,j,k)**2
          end do
       end do
    end do

    !MPI find global residual
    root = 0
    call mpi_reduce(res_local, res, 1, mpi_double_precision,&
         & mpi_sum, root, comm, ierror)
    
    call mpi_bcast(res, 1, mpi_double_precision, root, comm, ierror)

    res = sqrt(res/sizeX/sizeY/sizeZ)

    if( master ) then
       print*,""
       print*," BiCGStab solver: "
       if (res == 0.0) then
          if (giveInfo) write(*,fmt="(a25,es17.6)") &
            & " Initial residual: res = ", res
         else
          if (giveInfo) write(*,fmt="(a25,es17.6)") &
            & " Initial residual: res = ", res/b_norm
       end if
       if (giveInfo) write(*,fmt="(a25,es17.6)") " tol = ", tol
    end if

    if (res == 0.0 .or. res/b_norm <= tol) then
       if(master .and. giveInfo) print*," ==> no iteration needed."
       nIter = 0
       return
    end if

    ! Loop

    iteration: do j_b = 1,maxIt

       ! v = A*p
       if (preconditioner == 'yes') then
          call preCond( p, v_pc)
         else
          v_pc = p
       end if
       call linOpr( v_pc, matVec, opt )
       v = matVec

       alpha = dot_product3D_glob(r,r0) / dot_product3D_glob(v,r0)
       s = r - alpha*v

       ! t = A*s
       if (preconditioner == 'yes') then
          call preCond( s, v_pc)
         else
          v_pc = s
       end if
       call linOpr( v_pc, matVec, opt )
       t = matVec

       omega = dot_product3D_glob(t,s) / dot_product3D_glob(t,t)
       sol = sol + alpha*p + omega*s

       rOld = r
       r = s - omega*t

       !-----------------------
       !   Abort criterion 
       !-----------------------

       res_local = 0.0
       do k=1,nz
          do j=1,ny
             do i=1,nx
                res_local = res_local + r(i,j,k)**2
             end do
          end do
       end do

       !MPI find global residual
       root = 0
       call mpi_reduce(res_local, res, 1, mpi_double_precision,&
            & mpi_sum, root, comm, ierror)
    
       call mpi_bcast(res, 1, mpi_double_precision, root, comm, ierror)
      
       res = sqrt(res/sizeX/sizeY/sizeZ)

       if (res/b_norm <= tol) then
          if (master .and. giveInfo)  then 
             write(*,fmt="(a25,i25)") " Nb.of iterations: j = ", j_b
             write(*,fmt="(a25,es17.6)") " Final residual: res = ", &
                                         & res/b_norm
             print*,""
          end if

          nIter = j_b

          if (preconditioner == 'yes') then
             s = sol
             call preCond( s, sol)
          end if

          return
       end if

       beta &
       = alpha/omega &
         * dot_product3D_glob(r,r0) / dot_product3D_glob(rOld,r0)
       p = r + beta*(p-omega*v)

    end do iteration

    ! max iteration

    if( master ) then ! modified by Junhong Wei (20161107)
       write(*,fmt="(a25,i25)") " Bicgstab: max iterations!!!", maxIt
       write(*,fmt="(a25,es17.6)") " Final BICGSTAB residual = ", &
                                 & res/b_norm
       print*,"--------------------------------------------------" 
       print*,""
    end if ! modified by Junhong Wei (20161107)

    errFlag = .true.
    nIter = j_b

    ! deallocate local fields
    deallocate(p, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc p failed"
    deallocate(r0, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc r0 failed"
    deallocate(rOld,stat=allocstat);if(allocstat/=0) &
         & stop "bicgstab:dealloc rOld failed"
    deallocate(r, stat=allocstat);if(allocstat/=0) &
         & stop "bicgstab:dealloc r failed"
    deallocate(s, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc s failed"
    deallocate(t, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc t failed"
    deallocate(v, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc v failed"
    deallocate(v_pc, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc v_pcfailed"
    deallocate(matVec, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc matvec failed"

  end subroutine bicgstab

  !----------------------------------------------------------------------

  subroutine bicgstab_2(b, dt, sol, res, nIter, errFlag, opt)
    ! --------------------------------------
    !    BiCGStab using linear operator
    !    preconditioner applied via M^-1 A x = M^-1 b
    !---------------------------------------

    ! in/out variables
    real, dimension(1:nx,1:ny,1:nz), intent(in) :: b        ! RHS 
    real, intent(in) :: dt
    real, dimension(1:nx,1:ny,1:nz), intent(inout) :: sol
    real, intent(out) :: res                     ! residual
    integer, intent(out) :: nIter
    logical, intent(out) :: errFlag

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction 
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction 
    ! of the winds and density fluctuations
    character(len=*), intent(in)    :: opt

    ! Local parameters
    integer  :: maxIt

    ! local variables
    integer :: i,j,k, allocstat
    integer :: j_b
    real, dimension(:,:,:), allocatable :: p,r0,rOld,r,s,t,v, matVec, &
                                           v_pc, b_int
    real :: alpha, beta, omega

    ! verbose
    logical, parameter :: giveInfo = .true.

    ! MPI stuff
    integer :: root
    real :: res_local
    real :: bi_norm_local, bi_norm

    if( giveInfo .and. master ) then
       print*,""
       print*,"----------------------------------------------"
       print*,"BICGSTAB: solving linear system... "
       print*,"----------------------------------------------"
       print*,""
    end if

    sol = 0.0

    ! Set parameters
    maxIt = maxIterPoisson

    !if( opt == 'initial' ) then
    !   print*,"bicgstab: use tolInitial = ", tolInitial
    !   tol = tolInitial
    !else
    !   tol = tolPoisson
    !end if

    ! modified convergence criterion so that iterations stop when either
    ! (a) tolcrit = abs  =>  |Ax - b| < eps b_*
    !     with b_* a suitable norm deciding whether b (= the divergence 
    !     criterion for the winds from the predictor) is small or not
    ! (b) tolcrit = rel  =>  |Ax - b| < eps |b|
    ! here eps = tolPoisson is the user-set convergence criterion
    ! hypre has the criterion |Ax - b| < tol * |b|, hence, with
    ! tolref = divL2/divL2_norm = |b|/b_*

    if (tolcrit == "abs") then
       tol = tolPoisson/tolref
      else if (tolcrit == "rel") then
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
    allocate(v_pc(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:alloc failed"
    allocate(b_int(1:nx,1:ny,1:nz), stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:alloc failed"

    ! error flag
    errFlag = .false.    

    ! redefine RHS and its norm for preconditioner
    if (preconditioner == 'yes') then
       call preCond( b, b_int)

       bi_norm_local = 0.0
       do k=1,nz
          do j=1,ny
             do i=1,nx
                bi_norm_local = bi_norm_local + b_int(i,j,k)**2
             end do
          end do
       end do

       !MPI find global residual
       root = 0
       call mpi_reduce(bi_norm_local, bi_norm, 1, mpi_double_precision,&
            & mpi_sum, root, comm, ierror)
    
       call mpi_bcast(bi_norm, 1, mpi_double_precision, root, comm, ierror)

       bi_norm = sqrt(bi_norm/sizeX/sizeY/sizeZ)
      else
       b_int = b
       bi_norm = b_norm
    end if

    ! Init
    ! r0 = b - Ax
    call linOpr( sol, matVec, opt )
    if (preconditioner == 'yes') then
       v = matVec
       call preCond( v, matVec)
    end if
    r0 = b_int - matVec    
    p = r0
    r = r0

    res_local = 0.0
    do k=1,nz
       do j=1,ny
          do i=1,nx
             res_local = res_local + r(i,j,k)**2
          end do
       end do
    end do

    !MPI find global residual
    root = 0
    call mpi_reduce(res_local, res, 1, mpi_double_precision,&
         & mpi_sum, root, comm, ierror)
    
    call mpi_bcast(res, 1, mpi_double_precision, root, comm, ierror)

    res = sqrt(res/sizeX/sizeY/sizeZ)

    if( master ) then
       print*,""
       print*," BiCGStab solver: "
       if (res == 0.0) then
          if (giveInfo) write(*,fmt="(a25,es17.6)") &
            & " Initial residual: res = ", res
         else
          if (giveInfo) write(*,fmt="(a25,es17.6)") &
            & " Initial residual: res = ", res/bi_norm
       end if
       if (giveInfo) write(*,fmt="(a25,es17.6)") " tol = ", tol
    end if

    if (res == 0.0 .or. res/bi_norm <= tol) then
       if(master .and. giveInfo) print*," ==> no iteration needed."
       nIter = 0
       return
    end if

    ! Loop

    iteration: do j_b = 1,maxIt

       ! v = A*p
       call linOpr( p, matVec, opt )
       if (preconditioner == 'yes') then
          v_pc = matVec
          call preCond( v_pc, matVec)
       end if
       v = matVec

       alpha = dot_product3D_glob(r,r0) / dot_product3D_glob(v,r0)
       s = r - alpha*v

       ! t = A*s
       call linOpr( s, matVec, opt )
       if (preconditioner == 'yes') then
          v_pc = matVec
          call preCond( v_pc, matVec)
       end if
       t = matVec

       omega = dot_product3D_glob(t,s) / dot_product3D_glob(t,t)
       sol = sol + alpha*p + omega*s

       rOld = r
       r = s - omega*t

       !-----------------------
       !   Abort criterion 
       !-----------------------

       res_local = 0.0
       do k=1,nz
          do j=1,ny
             do i=1,nx
                res_local = res_local + r(i,j,k)**2
             end do
          end do
       end do

       !MPI find global residual
       root = 0
       call mpi_reduce(res_local, res, 1, mpi_double_precision,&
            & mpi_sum, root, comm, ierror)
    
       call mpi_bcast(res, 1, mpi_double_precision, root, comm, ierror)
      
       res = sqrt(res/sizeX/sizeY/sizeZ)

       if (res/bi_norm <= tol) then
          if (master .and. giveInfo)  then 
             write(*,fmt="(a25,i25)") " Nb.of iterations: j = ", j_b
             write(*,fmt="(a25,es17.6)") " Final int. residual = ", &
                                         & res/bi_norm
             print*,""
          end if

          call linOpr( sol, matVec, opt )
          r = b - matVec    

          res_local = 0.0
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   res_local = res_local + r(i,j,k)**2
                end do
             end do
          end do

          !MPI find global residual
          root = 0
          call mpi_reduce(res_local, res, 1, mpi_double_precision,&
               & mpi_sum, root, comm, ierror)
    
          call mpi_bcast(res, 1, mpi_double_precision, root, comm, ierror)

          res = sqrt(res/sizeX/sizeY/sizeZ)

          if (master .and. giveInfo)  then 
             write(*,fmt="(a25,es17.6)") " Final residual: res = ", &
                                         & res/b_norm
             print*,""
          end if

          nIter = j_b

          return
       end if

       beta &
       = alpha/omega &
         * dot_product3D_glob(r,r0) / dot_product3D_glob(rOld,r0)
       p = r + beta*(p-omega*v)

    end do iteration

    ! max iteration

    if( master ) then ! modified by Junhong Wei (20161107)
       write(*,fmt="(a25,i25)") " Bicgstab: max iterations!!!", maxIt
       write(*,fmt="(a25,es17.6)") " Final BICGSTAB residual = ", &
                                 & res/bi_norm
       print*,"--------------------------------------------------" 
       print*,""
    end if ! modified by Junhong Wei (20161107)

    errFlag = .true.
    nIter = j_b

    ! deallocate local fields
    deallocate(p, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc p failed"
    deallocate(r0, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc r0 failed"
    deallocate(rOld,stat=allocstat);if(allocstat/=0) &
         & stop "bicgstab:dealloc rOld failed"
    deallocate(r, stat=allocstat);if(allocstat/=0) &
         & stop "bicgstab:dealloc r failed"
    deallocate(s, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc s failed"
    deallocate(t, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc t failed"
    deallocate(v, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc v failed"
    deallocate(v_pc, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc v_pcfailed"
    deallocate(b_int, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc v_pcfailed"
    deallocate(matVec, stat=allocstat); if(allocstat/=0) &
         & stop "bicgstab:dealloc matvec failed"

  end subroutine bicgstab_2

  ! modified by Junhong Wei (2016/07/08) (starting line)

  !----------------------------------------------------------------------

  subroutine hypre(b, dt, sol, res, nIter, errFlag, opt)

    ! --------------------------------------
    !    HYPRE
    !---------------------------------------

    ! in/out variables

    real, dimension(1:nx,1:ny,1:nz), intent(in) :: b        ! RHS 
    real, intent(in) :: dt
    real, dimension(1:nx,1:ny,1:nz), intent(out) :: sol
    real, intent(out) :: res                                ! residual
    integer, intent(out) :: nIter
    logical, intent(out) :: errFlag

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction 
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction 
    ! of the winds and density fluctuations
    character(len=*), intent(in)    :: opt

    ! variables due to the use of HYPRE

    integer :: ndim_hypre
    parameter (ndim_hypre = 3)

    integer ierr_hypre
    integer npts_x_periodic_hypre, npts_y_periodic_hypre, &
          & npts_z_periodic_hypre

    integer ii_entry_hypre

    integer :: index_count_hypre

    integer :: i0,j0, i, j, k

    real :: atol

    ! Local parameters
    integer  :: maxIt

    ! verbose
    logical, parameter :: giveInfo = .true.

    real :: sol_mean_hypre
    
    ! Set parameters

    ! modified convergence criterion so that iterations stop when either
    ! (a) tolcrit = abs  =>  |Ax - b| < eps b_*
    !     with b_* a suitable norm deciding whether b (= the divergence 
    !     criterion for the winds from the predictor) is small or not
    ! (b) tolcrit = rel  =>  |Ax - b| < eps |b|
    ! here eps = tolPoisson is the user-set convergence criterion
    ! hypre has the criterion |Ax - b| < tol * |b|, hence, with
    ! tolref = divL2/divL2_norm = |b|/b_*

    if (tolcrit == "abs") then
       tol = tolPoisson/tolref
      else if (tolcrit == "rel") then
       tol = tolPoisson
    end if
    
    ! atol = abs_tol/b_norm  !not scaled as it is a lower bound

    ! if(master) then
    !    print*,"hypre uses max from ", " tol = ", tol, "atol = ", atol
    ! end if

    ! tol = max(tol, atol)  

    if(master)  then 
        print*,"hypre uses tolerance = ", tol, &
             & "meaning |Ax - b|/|b| <= tol "
        print*,"|Ax - b| <= ", tol*b_norm
    end if

    ! error flag
    errFlag = .false.

    i0=is+nbx-1
    j0=js+nby-1

    if (timeScheme == "semiimplicit") then
       ! This is a collective call finalizing the matrix assembly. 
       ! The matrix is now ``ready to be used'

       call HYPRE_StructMatrixSetBoxValues(A_hp_i, &
            & (/ (icoord)*nx+1, (jcoord)*ny+1 , 1 /), &
            & (/ (icoord+1)*nx , (jcoord+1)*ny , nz /), &
            & ne_hypre_i, stencil_indices_i, values_i, ierr_hypre)

       ! testb
       if (master) print*,'HYPRE_StructMatrixSetBoxValues done'
       ! teste

       call HYPRE_StructMatrixAssemble(A_hp_i,ierr_hypre)

       ! testb
       if (master) print*,'HYPRE_StructMatrixAssemble done'
       ! teste

       ! Set up Struct Vectors for b and x.  Each processor sets the 
       ! vectors corresponding to its boxes.

       ! Set the vector coefficients

       index_count_hypre = 1
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                bvalue_vector_hypre(index_count_hypre) = b(i,j,k)
                index_count_hypre = index_count_hypre + 1
                ! testb
                ! print*,i,j,k,'b =',b(i,j,k)
                ! teste
             end do
          end do
       end do

       index_count_hypre = 1
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                xvalue_vector_hypre(index_count_hypre) = sol(i,j,k)
                index_count_hypre = index_count_hypre + 1
                ! testb
                ! print*,i,j,k,'sol =',sol(i,j,k)
                ! teste
             end do
          end do
       end do

       call HYPRE_StructVectorSetBoxValues(b_hp_i, &
          & (/ (icoord)*nx+1, (jcoord)*ny+1 , 1 /), &
          & (/ (icoord+1)*nx , (jcoord+1)*ny , nz /), &
          & bvalue_vector_hypre, ierr_hypre)

       ! testb
       if (master) print*,'HYPRE_StructVectorSetBoxValues done'
       ! teste

       call HYPRE_StructVectorSetBoxValues(x_hp_i, &
          & (/ (icoord)*nx+1, (jcoord)*ny+1 , 1 /), &
          & (/ (icoord+1)*nx , (jcoord+1)*ny , nz /), &
          & xvalue_vector_hypre, ierr_hypre)

       ! testb
       if (master) print*,'HYPRE_StructVectorSetBoxValues done'
       ! teste

       ! This is a collective call finalizing the vector assembly. 
       ! The vectors are now ``ready to be used''

       call HYPRE_StructVectorAssemble(b_hp_i, ierr_hypre)

       ! testb
       if (master) print*,'HYPRE_StructVectorAssemble done'
       ! teste

       call HYPRE_StructVectorAssemble(x_hp_i, ierr_hypre)

       ! testb
       if (master) print*,'HYPRE_StructVectorAssemble done'
       ! teste

       ! Finalize set up and use a solver 
       ! (See the Reference Manual for descriptions of all of the options.)

       call HYPRE_StructHybridSetTol(solver_hp_i, tol, ierr_hypre)

       ! testb
       if (master) print*,'HYPRE_StructHybridSetTol done'
       ! teste

       call HYPRE_StructHybridSolve(solver_hp_i, A_hp_i, b_hp_i, x_hp_i, &
                                  & ierr_hypre)

       ! testb
       if (master) print*,'HYPRE_StructHybridSolve done'
       ! teste

       ! Get the results

       call HYPRE_StructVectorGetBoxValues(x_hp_i, &
          & (/ (icoord)*nx+1, (jcoord)*ny+1 , 1 /), &
          & (/ (icoord+1)*nx , (jcoord+1)*ny , nz /), &
          & xvalue_vector_hypre, ierr_hypre)

       ! testb
       if (master) print*,'HYPRE_StructVectorGetBoxValues done'
       ! teste

       index_count_hypre = 1

       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                sol(i,j,k) = xvalue_vector_hypre(index_count_hypre)
                index_count_hypre = index_count_hypre + 1
             end do
          end do
       end do

       call HYPRE_StructHybridGetNumIterati( solver_hp_i, nIter, &
                                           & ierr_hypre )

       ! testb
       if (master) print*,'HYPRE_StructHybridGetNumIterati done'
       ! teste

       index_count_hypre = 1
       call HYPRE_StructHybridGetFinalRelat( solver_hp_i, res, ierr_hypre )

       ! testb
       if (master) print*,'HYPRE_StructHybridGetFinalRelat done'
       ! teste
      else
       ! This is a collective call finalizing the matrix assembly. 
       ! The matrix is now ``ready to be used'

       call HYPRE_StructMatrixSetBoxValues(A_hp_e, &
            & (/ (icoord)*nx+1, (jcoord)*ny+1 , 1 /), &
            & (/ (icoord+1)*nx , (jcoord+1)*ny , nz /), &
            & ne_hypre_e, stencil_indices_e, values_e, ierr_hypre)

       call HYPRE_StructMatrixAssemble(A_hp_e,ierr_hypre)

       ! Set up Struct Vectors for b and x.  Each processor sets the 
       ! vectors corresponding to its boxes.

       ! Set the vector coefficients

       index_count_hypre = 1
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                bvalue_vector_hypre(index_count_hypre) = b(i,j,k)
                index_count_hypre = index_count_hypre + 1
                ! testb
                ! print*,i,j,k,'b =',b(i,j,k)
                ! teste
             end do
          end do
       end do

       index_count_hypre = 1
       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                xvalue_vector_hypre(index_count_hypre) = sol(i,j,k)
                index_count_hypre = index_count_hypre + 1
                ! testb
                ! print*,i,j,k,'sol =',sol(i,j,k)
                ! teste
             end do
          end do
       end do

       call HYPRE_StructVectorSetBoxValues(b_hp_e, &
          & (/ (icoord)*nx+1, (jcoord)*ny+1 , 1 /), &
          & (/ (icoord+1)*nx , (jcoord+1)*ny , nz /), &
          & bvalue_vector_hypre, ierr_hypre)

       call HYPRE_StructVectorSetBoxValues(x_hp_e, &
          & (/ (icoord)*nx+1, (jcoord)*ny+1 , 1 /), &
          & (/ (icoord+1)*nx , (jcoord+1)*ny , nz /), &
          & xvalue_vector_hypre, ierr_hypre)

       ! This is a collective call finalizing the vector assembly. 
       ! The vectors are now ``ready to be used''

       call HYPRE_StructVectorAssemble(b_hp_e, ierr_hypre)
       call HYPRE_StructVectorAssemble(x_hp_e, ierr_hypre)

       ! Finalize set up and use a solver 
       ! (See the Reference Manual for descriptions of all of the options.)

       call HYPRE_StructHybridSetTol(solver_hp_e, tol, ierr_hypre)
       call HYPRE_StructHybridSolve(solver_hp_e, A_hp_e, b_hp_e, x_hp_e, &
                                  & ierr_hypre)

       ! Get the results

       call HYPRE_StructVectorGetBoxValues(x_hp_e, &
          & (/ (icoord)*nx+1, (jcoord)*ny+1 , 1 /), &
          & (/ (icoord+1)*nx , (jcoord+1)*ny , nz /), &
          & xvalue_vector_hypre, ierr_hypre)

       index_count_hypre = 1

       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                sol(i,j,k) = xvalue_vector_hypre(index_count_hypre)
                index_count_hypre = index_count_hypre + 1
             end do
          end do
       end do

       call HYPRE_StructHybridGetNumIterati( solver_hp_e, nIter, &
                                           & ierr_hypre )
       call HYPRE_StructHybridGetFinalRelat( solver_hp_e, res, ierr_hypre )
    end if

    if ( giveInfo .and. master )  then 
       print*,""
       print*," HYPRE solver: "
       write(*,fmt="(a25,i25)") " Nb.of iterations: j = ", nIter
       write(*,fmt="(a25,es17.6)") " Final residual: res = ", res
       print*,""
    end if

    return

  end subroutine hypre

  !-----------------------------------------------------------------------

  subroutine pressureBoundaryCondition
    !--------------------------------------------------
    ! set pressure correction dp in ghost cells for BC
    !--------------------------------------------------

    ! modified by Junhong Wei (20161107) *** starting line ***

    ! auxiliary fields for "dp"
    real, dimension(0:ny+1,0:nz+1) :: xSliceLeft_send, xSliceRight_send
    real, dimension(0:ny+1,0:nz+1) :: xSliceLeft_recv, xSliceRight_recv
    
    real, dimension(0:nx+1,0:nz+1) :: ySliceBack_send, ySliceForw_send
    real, dimension(0:nx+1,0:nz+1) :: ySliceBack_recv, ySliceForw_recv
    
    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount


    ! Find neighbour procs
    if ( idim > 1 ) call mpi_cart_shift(comm,0,1,left,right,ierror)
    if ( jdim > 1 ) call mpi_cart_shift(comm,1,1,back,forw,ierror)

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

    if( idim > 1 ) then
       ! slice size
       sendcount = (ny+2)*(nz+2)
       recvcount = sendcount

       ! read slice into contiguous array
       xSliceLeft_send (:,:) = dp(1, 0:ny+1,:)
       xSliceRight_send(:,:) = dp(nx,0:ny+1,:)



       ! left -> right
       source = left
       dest = right
       tag = 100

       call mpi_sendrecv(xSliceRight_send(0,0), sendcount, &
            & mpi_double_precision, dest, tag, &
            & xSliceLeft_recv(0,0), recvcount, mpi_double_precision, &
            & source, mpi_any_tag, comm, sts_left, ierror)

       ! right -> left
       source = right
       dest = left
       tag = 100

       call mpi_sendrecv(xSliceLeft_send(0,0), sendcount, &
            & mpi_double_precision, dest, tag, &
            & xSliceRight_recv(0,0), recvcount, mpi_double_precision, &
            & source, mpi_any_tag, comm, sts_right, ierror)

       ! right halos
       dp(nx+1,0:ny+1,:) = xSliceRight_recv(0:ny+1,:)

       ! left halos
       dp(0   ,0:ny+1,:) =  xSliceLeft_recv(0:ny+1,:)

    else

       dp(0,:,:) = dp(nx,:,:)
       dp(nx+1,:,:) = dp(1,:,:)

    end if
    
    if(verbose .and. master) print*,"horizontalHalos: &
         & x-horizontal halos copied."

    !------------------------------
    !   set Halo cells: ySlice 
    !------------------------------

    if( jdim > 1 ) then
       ! slice size
       sendcount = (nx+2)*(nz+2)

       recvcount = sendcount

       ! read slice into contiguous array
       ySliceBack_send(:,:) = dp(0:nx+1, 1,:)

       ySliceForw_send(:,:) = dp(0:nx+1,ny,:)

       ! back -> forw
       source = back
       dest = forw
       tag = 100

       call mpi_sendrecv(ySliceForw_send(0,0), sendcount, &
            & mpi_double_precision, dest, tag, &
            & ySliceBack_recv(0,0), recvcount, mpi_double_precision, &
            & source, mpi_any_tag, comm, sts_back, ierror)

       ! forw -> back
       source = forw
       dest = back
       tag = 100

       call mpi_sendrecv(ySliceBack_send(0,0), sendcount, &
            & mpi_double_precision, dest, tag, &
            & ySliceForw_recv(0,0), recvcount, mpi_double_precision, &
            & source, mpi_any_tag, comm, sts_right, ierror)

       ! forward halos
       dp(0:nx+1,ny+1,:) = ySliceForw_recv(0:nx+1,:)

       ! backward halos
       dp(0:nx+1,0   ,:) = ySliceBack_recv(0:nx+1,:)

    else

       dp(:,0,:) = dp(:,ny,:)
       dp(:,ny+1,:) = dp(:,1,:)

    end if

    if(verbose .and. master) print*,"horizontalHalos: &
         & x-horizontal halos copied."

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
       stop "pressureBoundaryCondition: unknown case zBoundary."
    end select


  end subroutine pressureBoundaryCondition

  !-----------------------------------------------------------------------

  subroutine correctorStep( var, dMom, dt, RKstage, opt)

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

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction 
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction 
    ! of the winds and density fluctuations
    character(len=*), intent(in)    :: opt

    ! local variables 
    integer :: i,j,k
    real :: pEdge, rhoEdge, rhou, rhov, rho, rho10, rho01
    real :: pGradX, pGradY, pGradZ
    real :: du, dv, dw, db
    real :: facu, facv, facw, facr
    !real :: f_cor_nd
    real, dimension(0:ny+1) :: f_cor_nd
    real :: bvsstw

    real :: rhov0m, rhov00, rhov1m, rhov10
    real :: rhoum0, rhou00, rhoum1, rhou01
    real :: rhow0, rhowm

    ! divergence check
    real :: maxDivPu, divPu
    real :: uR, uL, vF, vB, wU, wD
    real :: pStratU, pStratD 

    ! verbose
    logical, parameter :: giveInfo = .true.

    integer :: i0,j0
    real :: yloc, ymax

    if (model == "Boussinesq") then
       print*,'ERROR: correctorStep not ready for Boussinesq mode'
       stop
    end if

    i0=is+nbx-1
    j0=js+nby-1

    ! non-dimensional Corilois parameter (= inverse Rossby number) 
    if (TestCase == "baroclinic_LC")then !FS
       ymax = ly_dim(1)/lRef  
       do j = 1,ny
          yloc = y(j+j0)
          f_cor_nd(j) = f_Coriolis_dim*tRef*sin(pi*yloc/ymax)
       end do
    else
       f_cor_nd(:) = f_Coriolis_dim*tRef
    end if
    

    ! --------------------------------------
    !             calc p + dp
    ! --------------------------------------

    var(0:nx+1,0:ny+1,0:nz+1,5) &
    = var(0:nx+1,0:ny+1,0:nz+1,5) + dp(0:nx+1,0:ny+1,0:nz+1)

    ! --------------------------------------
    !           calc du and u + du
    ! --------------------------------------

    if (timeScheme == "semiimplicit") then
       if (opt == "impl") then
          do k = 1,nz
             do j = 1,ny
                do i = 0,nx
                   facu = 1.0

                   if (topography) then
                      if (   topography_mask(i0+i,j0+j,k) &
                        & .or. topography_mask(i0+i+1,j0+j,k)) then
                         facu = facu + dt*alprlx
                      end if
                   end if

                   facv = facu

                   rhov0m = 0.5 * (var(i,j-1,k,1) + var(i,j,k,1))
                   if( fluctuationMode ) rhov0m = rhov0m + rhoStrat(k)

                   rhov00 = 0.5 * (var(i,j,k,1) + var(i,j+1,k,1))
                   if( fluctuationMode ) rhov00 = rhov00 + rhoStrat(k)

                   rhov1m = 0.5 * (var(i+1,j-1,k,1) + var(i+1,j,k,1))
                   if( fluctuationMode ) rhov1m = rhov1m + rhoStrat(k)

                   rhov10 = 0.5 * (var(i+1,j,k,1) + var(i+1,j+1,k,1))
                   if( fluctuationMode ) rhov10 = rhov10 + rhoStrat(k)

                   rhou = 0.5*( var(i+1,j,k,1) + var(i,j,k,1) )
                   if( fluctuationMode ) rhou = rhou + rhoStrat(k)
                   pGradX &
                   = kappaInv*MaInv2 * pStrat(k)/rhou &
                     * ( dp(i+1,j,k) - dp(i,j,k) ) / dx

                   pGradY &
                   = kappaInv*MaInv2 &
                     * 0.25 &
                     * (  pStrat(k)/rhov0m * (dp(i,j,k) - dp(i,j-1,k))/dy &
                        + pStrat(k)/rhov00 * (dp(i,j+1,k) - dp(i,j,k))/dy &
                        + pStrat(k)/rhov1m &
                          * (dp(i+1,j,k) - dp(i+1,j-1,k))/dy &
                        + pStrat(k)/rhov10 &
                          * (dp(i+1,j+1,k) - dp(i+1,j,k))/dy)

                          du = - dt/(facu*facv + (f_cor_nd(j)*dt)**2) &
                          * (facv * pGradX + f_cor_nd(j)*dt * pGradY)

                   var(i,j,k,2) = var(i,j,k,2) + du
                end do
             end do
          end do
         else if (opt == "expl") then
          do k = 1,nz
             do j = 1,ny
                do i = 0,nx
                   rhou = 0.5 * ( var(i,j,k,1) + var(i+1,j,k,1) )
                   if( fluctuationMode ) rhou = rhou + rhoStrat(k)

                   pGradX &
                   = kappaInv*MaInv2 * pStrat(k)/rhou &
                     * (dp(i+1,j,k) - dp(i,j,k))/dx

                   du = - dt * pGradX

                   var(i,j,k,2) = var(i,j,k,2) + du
                end do
             end do
          end do
         else
          stop'ERROR: wrong opt in correctorStep'
       end if
      else
       do k = 1,nz
          do j = 1,ny
             do i = 0,nx
                rhou = 0.5 * ( var(i,j,k,1) + var(i+1,j,k,1) )
                if( fluctuationMode ) rhou = rhou + rhoStrat(k)

                pGradX &
                = kappaInv*MaInv2 * pStrat(k)/rhou &
                  * (dp(i+1,j,k) - dp(i,j,k))/dx

                du = - dt *  pGradX

                var(i,j,k,2) = var(i,j,k,2) + du

                ! correct x-momentum tendency
                dMom(i,j,k,1) = dMom(i,j,k,1) + rhou*du/beta(RKstage)
             end do
          end do
       end do
    end if

    !--------------------------------------
    !         calc dv and v + dv
    !--------------------------------------

    if (timeScheme == "semiimplicit") then
       if (opt == "impl") then
          do k = 1,nz
             do j = 0,ny
                do i = 1,nx
                   facv = 1.0

                   if (topography) then
                      if (   topography_mask(i0+i,j0+j,k) &
                        & .or. topography_mask(i0+i,j0+j+1,k)) then
                         facv = facv + dt*alprlx
                      end if
                   end if

                   facu = facv

                   rhou00 = 0.5 * (var(i,j,k,1) + var(i+1,j,k,1))
                   if( fluctuationMode ) rhou00 = rhou00 + rhoStrat(k)

                   rhoum0 = 0.5 * (var(i-1,j,k,1) + var(i,j,k,1))
                   if( fluctuationMode ) rhoum0 = rhoum0 + rhoStrat(k)

                   rhou01 = 0.5 * (var(i,j+1,k,1) + var(i+1,j+1,k,1))
                   if( fluctuationMode ) rhou01 = rhou01 + rhoStrat(k)

                   rhoum1 = 0.5 * (var(i-1,j+1,k,1) + var(i,j+1,k,1))
                   if( fluctuationMode ) rhoum1 = rhoum1 + rhoStrat(k)

                   rhov = 0.5*( var(i,j+1,k,1) + var(i,j,k,1) )
                   if( fluctuationMode ) rhov = rhov + rhoStrat(k)

                   pGradX &
                   = kappaInv*MaInv2 &
                     * 0.25 &
                     * (  pStrat(k)/rhou00 * (dp(i+1,j,k) - dp(i,j,k))/dx &
                        + pStrat(k)/rhoum0 * (dp(i,j,k) - dp(i-1,j,k))/dx &
                        + pStrat(k)/rhou01 &
                          * (dp(i+1,j+1,k) - dp(i,j+1,k))/dx &
                        + pStrat(k)/rhoum1 &
                          * (dp(i,j+1,k) - dp(i-1,j+1,k))/dx)

                   pGradY &
                   = kappaInv*MaInv2 * pStrat(k)/rhov &
                     * (dp(i,j+1,k) - dp(i,j,k))/dy

                   dv = - dt/(facu*facv + (f_cor_nd(j)*dt)**2) &
                          * (- f_cor_nd(j)*dt * pGradX + facu * pGradY)

                   var(i,j,k,3) = var(i,j,k,3) + dv
                end do
             end do
          end do
         else if (opt == "expl") then
          do k = 1,nz
             do j = 0,ny
                do i = 1,nx
                   rhov = 0.5 * (var(i,j,k,1) + var(i,j+1,k,1))
                   if( fluctuationMode ) rhov = rhov + rhoStrat(k)
             
                   pGradY &
                   = kappaInv*MaInv2 * pStrat(k)/rhov &
                     * (dp(i,j+1,k) - dp(i,j,k))/dy

                   dv = - dt * pGradY

                   var(i,j,k,3) = var(i,j,k,3) + dv
                end do
             end do
          end do
         else
          stop'ERROR: wrong opt in correctorStep'
       end if
      else
       do k = 1,nz
          do j = 0,ny
             do i = 1,nx
                rhov = 0.5 * (var(i,j,k,1) + var(i,j+1,k,1))
                if( fluctuationMode ) rhov = rhov + rhoStrat(k)
             
                pGradY &
                = kappaInv*MaInv2 * pStrat(k)/rhov &
                  * (dp(i,j+1,k) - dp(i,j,k))/dy

                dv = - dt * pGradY

                var(i,j,k,3) = var(i,j,k,3) + dv

                ! correct y-momentum tendency
                dMom(i,j,k,2) = dMom(i,j,k,2) + rhov*dv/beta(RKstage)
             end do
          end do
       end do
    end if


    !--------------------------------------
    !         calc w and  w + dw
    !--------------------------------------

    if (timeScheme == "semiimplicit") then
       if (opt == "impl") then
          ! solid wall implies zero change of w at the bottom and top

          do k = 1, nz-1        
             do j = 1,ny
                do i = 1,nx
                   facw = 1.0
                   facr = 1.0 + alprlx*dt

                   if (topography) then
                      if (   topography_mask(i0+i,j0+j,k) &
                        & .or. topography_mask(i0+i,j0+j,k+1)) then
                         facw = facw + dt*alprlx
                      end if
                   end if

                   pEdge = 0.5 * ( pStrat(k+1) + pStrat(k) )

                   rhoEdge = 0.5 * ( var(i,j,k,1) + var(i,j,k+1,1) )
                   if( fluctuationMode ) then
                       rhoEdge = rhoEdge + rhoStratTilde(k)
                   end if
             
                   bvsstw = 0.5 * (bvsStrat(k) + bvsStrat(k+1))

                   pGradZ = ( dp(i,j,k+1) - dp(i,j,k) ) / dz

                   dw = - dt * kappaInv*MaInv2 * facr &
                          /(  facw*facr &
                            + rhoStratTilde(k)/rhoEdge * bvsstw * dt**2) &
                          * pEdge/rhoEdge * pGradz   

                   var(i,j,k,4) = var(i,j,k,4) + dw
                end do
             end do
          end do

          ! periodic in z
          if( zBoundary == "periodic" ) then
              stop'ERROR: period. vert. bound. not ready in correctorStep'
          end if
         else if (opt == "expl") then
          ! solid wall implies zero change of w at the bottom and top

          do k = 1, nz-1        
             do j = 1,ny
                do i = 1,nx
                   if( fluctuationMode ) then
                      pEdge = pStratTilde(k)
                   else
                      pEdge = 0.5 * ( pStrat(k) + pStrat(k+1) )
                   end if
             
                   rhoEdge = 0.5 * ( var(i,j,k,1) + var(i,j,k+1,1) )
                   if( fluctuationMode ) then
                       rhoEdge = rhoEdge + rhoStratTilde(k)
                   end if
             
                   pGradZ = ( dp(i,j,k+1) - dp(i,j,k) ) / dz

                   dw = -dt * kappaInv * MaInv2 * pEdge/rhoEdge * pGradz   

                   var(i,j,k,4) = var(i,j,k,4) + dw
                end do
             end do
          end do

          ! if periodic in z
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
                      if( fluctuationMode ) then
                          rhoEdge = rhoEdge + rhoStratTilde(k)
                      end if
                
                      pGradZ = ( dp(i,j,k+1) - dp(i,j,k) ) / dz

                      dw = -dt * kappaInv*MaInv2 * pEdge/rhoEdge * pGradz

                      var(i,j,k,4) = var(i,j,k,4) + dw
                   end do
                end do
             end do
          end if
         else
          stop'ERROR: wrong opt in correctorStep'
       end if
      else
       ! solid wall implies zero change of w at the bottom and top

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

       ! if periodic in z
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
                   if( fluctuationMode ) then
                       rhoEdge = rhoEdge + rhoStratTilde(k)
                   end if
                
                   pGradZ = ( dp(i,j,k+1) - dp(i,j,k) ) / dz

                   dw = -dt * kappaInv*MaInv2 * pEdge/rhoEdge * pGradz   

                   var(i,j,k,4) = var(i,j,k,4) + dw

                   ! correct z-momentum tendency
                   dMom(i,j,k,3) = dMom(i,j,k,3) + rhoEdge*dw/beta(RKstage)
                end do
             end do
          end do
       end if
    end if

    !------------------------------------------------------------------
    !         calc rhop and rhop + drhop (only for implicit time step)
    !------------------------------------------------------------------

    if (timeScheme == "semiimplicit") then
       if (opt == "impl") then
          do k = 1, nz        
             do j = 1,ny
                do i = 1,nx
                   facw = 1.0
                   facr = 1.0 + alprlx*dt

                   if (topography) then
                      if (   topography_mask(i0+i,j0+j,k) &
                        & .or. topography_mask(i0+i,j0+j,k+1)) then
                         facw = facw + dt*alprlx
                      end if
                   end if

                   rho = var(i,j,k,1)
                   if( fluctuationMode ) rho = rho + rhoStrat(k)

                   rhowm = 0.5 * (var(i,j,k-1,1) + var(i,j,k,1))
                   if( fluctuationMode ) rhowm = rhowm + rhoStratTilde(k-1)

                   rhow0 = 0.5 * (var(i,j,k,1) + var(i,j,k+1,1))
                   if( fluctuationMode ) rhow0 = rhow0 + rhoStratTilde(k)
             
                   if (k == 1) then
                      pGradZ &
                      = kappaInv*MaInv2 &
                        * 0.5*(pStratTilde(k)/rhow0 &
                               * (dp(i,j,k+1) - dp(i,j,k))/dz)
                     else if (k == nz) then
                      pGradZ &
                      = kappaInv*MaInv2 &
                        * 0.5*(pStratTilde(k-1)/rhowm &
                               * (dp(i,j,k) - dp(i,j,k-1))/dz)
                     else
                      pGradZ &
                      = kappaInv*MaInv2 &
                        * 0.5*(  pStratTilde(k)/rhow0 &
                                 * (dp(i,j,k+1) - dp(i,j,k))/dz &
                               + pStratTilde(k-1)/rhowm &
                                 * (dp(i,j,k) - dp(i,j,k-1))/dz)
                   end if

                   db = rhoStrat(k)/rho * bvsStrat(k) * dt**2 &
                        /(  facw*facr &
                          + rhoStrat(k)/rho * bvsStrat(k) * dt**2) &
                        * pGradz   

                   var(i,j,k,6) = var(i,j,k,6) - rho/g_ndim * db
                end do
             end do
          end do

          ! periodic in z
          if( zBoundary == "periodic" ) then
              stop'ERROR: period. vert. bound. not ready in correctorStep'
          end if
       end if
    end if

  end subroutine correctorStep

  !----------------------------------------------------------------------

  function getIndex(i,j,k)
    !----------------------------------------------------
    ! compute matrix index l from spatial indicies i,j,k
    !----------------------------------------------------

    ! in/out variables
    integer :: getIndex 
    integer, intent(in) :: i,j,k 

    getIndex = (k-1)*nx*ny + (j-1)*nx + i

  end function getIndex

  !--------------------------------------------------------------------

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

  !----------------------------------------------------------------------

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

  !!------------------------------------------------------------------------
!
 !  subroutine calculate_heating_0(var,flux,heat)

 !  !-----------------------------------------------
 !  ! calculate heating in the divergence constraint
 !  !-----------------------------------------------

 !   real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz, nVar), &
 !        &intent(in) :: var     
 !   real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), intent(in) :: flux

 !   real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz), &
 !        &intent(out) :: heat  !, term1, term2

 !   ! local variables
 !   integer :: i, j, k
 !   real    :: rho, the, theta_bar_0, rho_e
 !   real, dimension(1:nz) :: tau_relax_i  ! inverse scaled relaxation 
 !                                         ! time for barocl. l. c.
 !   real    :: tau_jet_sc, tau_relax_sc  ! Klein scaling for relax param.

 !   heat = 0.0

 !   !-------------------------------------------------------
 !   ! calculate the environment-induced negative (!) heating
 !   !-------------------------------------------------------

 !   if (     (TestCase == "baroclinic_LC") &
 !       .or. (TestCase == "baroclinic_ID")) then
 !      if (RelaxHeating /= 0)  then
 !         !UAB
 !         if (background /= "HeldSuarez") then
 !         !UAE
 !            do k = 1,nz
 !               if ( referenceQuantities == "SI" ) then
 !                  tau_relax_sc = tau_relax
 !                  tau_jet_sc = tau_jet
 !                 else
 !                  tau_relax_sc = tau_relax/tref
 !                  tau_jet_sc = tau_jet/tref
 !               end if
           
 !               if (RelaxHeating == 1) then
 !                  tau_relax_i(k) = 1./(tau_relax_sc)
 !                  else
 !                  !UAB
 !                  !tau_relax_i = 0.
 !                  tau_relax_i(k) = 0.
 !                  !UAE
 !               end if
 !            end do
 !         !UAB
 !         end if
 !         !UAE

 !         if( master ) then 
 !            print*,""
 !            print*," Poisson Solver, Thermal Relaxation is on: "
 !            print*," Relaxation factor: Div = - rho(Th - Th_e)/tau: "
 !            print*, "tau = ", tref/tau_relax_i(1), " s"
 !         end if

 !         theta_bar_0 = (thetaStrat(1) + thetaStrat(nz))/2.
   
 !         do k = 1,nz
 !            do j = 1,ny
 !               do i = 1,nx
 !                  if( fluctuationMode )  then
 !                     rho = var(i,j,k,1) + rhoStrat(k)
 !                    else   
 !                     rho = var(i,j,k,1)
 !                  end if  

 !                  the = (Pstrat(k))/rho 
 !                  rho_e = (Pstrat(k))/the_env_pp(i,j,k)
 !                  !UAB
 !                  if (background == "HeldSuarez") then
 !                     heat(i,j,k) &
 !                     = - theta_bar_0*(rho - rho_e)*kt_hs(j,k)
 !                    else
 !                  !UAE
 !                     heat(i,j,k) &
 !                     = - theta_bar_0*(rho - rho_e)*tau_relax_i(k)
 !                  !UAB
 !                  end if
 !                  !UAE
 !               end do
 !            end do
 !         end do
 !      end if
 !   end if

 !  ! if (output_heat) then
 !      call output_field( &
 !      & iOut, &
 !      & heat,&
 !      & 'heat_prof.dat', thetaRef*rhoRef/tref )
 !  ! end if

 !   !testb
 !   return
 !   !teste

 !   !------------------------------------------------------------------
 !   ! supplement by negative (!) heating due to molecular and turbulent 
 !   ! diffusivity
 !   !------------------------------------------------------------------

 !   do k=1,nz
 !      do j=1,ny
 !         do i=1,nx
 !            heat(i,j,k) &
 !            = heat(i,j,k) &
 !              + rhoStrat(k) &
 !                * (  (flux(i,j,k,1,5) - flux(i-1,j  ,k  ,1,5))/dx &
 !                   + (flux(i,j,k,2,5) - flux(i  ,j-1,k  ,2,5))/dy &
 !                   + (flux(i,j,k,3,5) - flux(i  ,j  ,k-1,3,5))/dz)
 !         end do
 !      end do
 !   end do

 ! end subroutine calculate_heating_0

 !----------------------------------------------------------------------

 subroutine calculate_heating(var,flux,heat)

   !-----------------------------------------------------------------
   ! calculate heating in the divergence constraint
   ! supplemented by 'heating' due to turbulent and GW entropy fluxes
   !-----------------------------------------------------------------

    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz, nVar), &
         &intent(in) :: var     
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), intent(in) :: flux

    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz), &
         &intent(out) :: heat  !, term1, term2

    ! local variables
    integer :: i, j, k
    real    :: rho, the, theta_bar_0, rho_e
    real, dimension(1:ny,0:nz) :: tau_relax_i  ! inverse scaled relaxation 
                         !FS                 ! time for barocl. l. c.
    real    :: tau_jet_sc, tau_relax_sc  ! Klein scaling for relax param.

    real :: ymax, ymin, yloc
    integer :: j00

    real, dimension(1:nz) :: sum_local, sum_global

    heat = 0.0

    ymin = ly_dim(0)/lRef
    ymax = ly_dim(1)/lRef

    j00 = js+nby-1

    !-------------------------------------------------------
    ! calculate the environment-induced negative (!) heating
    !-------------------------------------------------------

    if (     (TestCase == "baroclinic_LC") &
        .or. (TestCase == "baroclinic_ID")) then

       !UAB
       if (background == "HeldSuarez") then
          do k=1,nz 
             do j=1,ny
                tau_relax_i(j,k) = kt_hs(j,k)
             end do
          end do
         else
       !UAE
          do k = 1,nz
             if ( referenceQuantities == "SI" ) then
                tau_relax_sc = tau_relax  !tau_z(k)  !tau_relax
                tau_jet_sc = tau_jet
               else
                tau_relax_sc = tau_relax/tref !tau_z(k)  !tau_relax/tref
                tau_jet_sc = tau_jet/tref
             end if
         

             do j = 1,ny
                yloc = y(j+j00)

                if (yloc > 0.5*(ymax+ymin)) then
                   ! meridionally dependent tau_sc

                   tau_relax_sc &
                   =   tau_relax/tref &
                     + (tau_relax_low/tref - tau_relax/tref) &
                       *( 1.0 &
                         -0.5&
                          *( tanh((yloc/ymax-0.1)/sigma_tau) &
                            -tanh((yloc/ymax-0.9)/sigma_tau)))  
                  else
                   tau_relax_sc &
                    =   tau_relax/tref &
                     + (tau_relax_low/tref - tau_relax/tref) &
                       *( 1.0 &
                         -0.5 &
                          *( tanh(((-1)*yloc/ymax-0.1)/sigma_tau) &
                            -tanh(((-1)*yloc/ymax-0.9)/sigma_tau))) 
                end if

                tau_relax_i(j,k) = 1./(tau_relax_sc)
             end do
          end do
       !UAB
       end if

       if (dens_relax) tau_relax_i = 0.0
       !UAE

       if( master ) then 
          print*,""
          print*," Poisson Solver, Thermal Relaxation is on: "
          print*," Relaxation factor: Div = - rho(Th - Th_e)/tau: "
          print*, "tau = ", tref/tau_relax_i(1,1), " s"
       end if

       !UA theta_bar_0 = (thetaStrat(1) + thetaStrat(nz))/2.
 
       !UAB 
       ! do k = 1,nz
       !   do j = 1,ny
       !      do i = 1,nx
       !         if( fluctuationMode )  then
       !            rho = var(i,j,k,1) + rhoStrat(k)
       !           else   
       !            rho = var(i,j,k,1)
       !         end if  

       !         the = (Pstrat(k))/rho 
       !         rho_e = (Pstrat(k))/the_env_pp(i,j,k)

       !         !UAB
       !         heat(i,j,k) &
       !         = rho*(Pstrat(k)/rho - Pstrat(k)/rho_e)*tau_relax_i(j,k)
       !         !UAE
       !      end do
       !   end do
       ! end do

       ! only the deviation of the potential-temperaure difference from
       ! the horizontal mean is taken ...

       ! potential-temperature difference

       do k = 1,nz 
          do j = 1,ny
             do i = 1,nx
                if( fluctuationMode )  then
                   rho = var(i,j,k,1) + rhoStrat(k)
                  else   
                   rho = var(i,j,k,1)
                end if  

                the = Pstrat(k)/rho 

                heat(i,j,k) = the - the_env_pp(i,j,k)
             end do
          end do
       end do

       ! subtract horizontal mean of the potential-temperature difference

      ! do k = 1,nz
      !     sum_local(k) = sum(heat(1:nx,1:ny,k))
      !  end do

      !  call mpi_allreduce(sum_local(1),sum_global(1),&
      !       nz, mpi_double_precision,mpi_sum,comm,ierror)
      !  sum_global = sum_global/(sizeX*sizeY)

      !  do k=1,nz
      !     heat(1:nx,1:ny,k) = heat(1:nx,1:ny,k) - sum_global(k)
      !  end do

       ! heating

       do k = 1,nz 
          do j = 1,ny
             do i = 1,nx
                if( fluctuationMode )  then
                   rho = var(i,j,k,1) + rhoStrat(k)
                  else   
                   rho = var(i,j,k,1)
                end if  

                heat(i,j,k) = rho*heat(i,j,k)*tau_relax_i(j,k)
             end do
          end do
       end do

    end if

    !UAB
    !------------------------------------------------------------------
    ! supplement by negative (!) heating due to molecular and turbulent 
    ! diffusivity
    !------------------------------------------------------------------

    if (TurbScheme) then
       do k=1,nz 
          do j=1,ny
             do i=1,nx
                heat(i,j,k) &
                = heat(i,j,k) &
                  + rhoStrat(k) &
                    * (  (flux(i,j,k,1,5) - flux(i-1,j  ,k  ,1,5))/dx &
                       + (flux(i,j,k,2,5) - flux(i  ,j-1,k  ,2,5))/dy &
                       + (flux(i,j,k,3,5) - flux(i  ,j  ,k-1,3,5))/dz)
             end do
          end do
       end do
    end if

    !-------------------------------------------------------------------
    ! supplement by negative (!) heating due GW entropy-flux convergence
    !-------------------------------------------------------------------

    if (raytracer) heat(:,:,:) = heat(:,:,:) + var(:,:,:,8)
    !UAE

    if (output_heat) then
       call output_field( &
       & iOut, &
       & heat,&
       & 'heat_prof.dat', thetaRef*rhoRef/tref )


    end if

  end subroutine calculate_heating
  
 !============================================================== 

  subroutine val_PsIn(var, dt, opt)

    ! calculates the matrix values for the pressure solver
    ! the solver solves for dt * dp, hence no dt in the matrix elements

    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz, nVar), &
         & intent(in) :: var
    real, intent(in) :: dt

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction 
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction 
    ! of the winds and density fluctuations
    character(len=*), intent(in)    :: opt

    ! local variables
    real :: dx2, dy2, dz2, dxy 
    real :: pStratU, pStratD, rhoEdge
    real :: AL,AR, AB,AF, AD,AU, ALB,ALF, ARB,ARF, AC
    real :: facu, facv, facw, facr
    real :: acontr
    real :: bvsstw

    ! non-dimensional Corilois parameter (= inverse Rossby number)
    real, dimension(0:ny+1) :: f_cor_nd

    integer :: i0,j0, i, j, k
    integer :: index_count_hypre
    real :: yloc, ymax

    i0=is+nbx-1
    j0=js+nby-1
    
    ! non-dimensional Corilois parameter (= inverse Rossby number)
    !f_cor_nd = f_Coriolis_dim*tRef
    if (TestCase == "baroclinic_LC")then !FS
       ymax = ly_dim(1)/lRef  
       do j = 1,ny
          yloc = y(j+j0)
          f_cor_nd(j) = f_Coriolis_dim*tRef*sin(pi*yloc/ymax)
       end do
    else
       f_cor_nd(:) = f_Coriolis_dim*tRef
    end if

    ! auxiliary variables
    dx2 = 1.0/dx**2
    dy2 = 1.0/dy**2
    dz2 = 1.0/dz**2    
    dxy = 1.0/(dx*dy)
    

    if (.not. fluctuationMode) stop'ERROR: must use fluctuationMode'

    !---------------------------------
    !         Loop over field
    !---------------------------------

    if (opt == "expl") then
       if (timeScheme == "semiimplicit") then
          do k = 1,nz
             do j = 1,ny
                do i = 1,nx
                   ! ------------------ A(i+1,j,k) ------------------

                   rhoEdge = 0.5*( var(i+1,j,k,1) + var(i,j,k,1) )
                   rhoEdge = rhoEdge + rhoStrat(k)

                   AR = dx2 * pStrat(k)**2/rhoEdge
   
                   ! ------------------- A(i-1,j,k) --------------------

                   rhoEdge = 0.5*( var(i,j,k,1) + var(i-1,j,k,1) )
                   rhoEdge = rhoEdge + rhoStrat(k)

                   AL = dx2 * pStrat(k)**2/rhoEdge

                   ! -------------------- A(i,j+1,k) ----------------------

                   rhoEdge = 0.5*( var(i,j+1,k,1) + var(i,j,k,1) )
                   rhoEdge = rhoEdge + rhoStrat(k)

                   AF = dy2 * pStrat(k)**2 / rhoEdge

                   ! --------------------- A(i,j-1,k) -------------------

                   rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j-1,k,1) )
                   rhoEdge = rhoEdge + rhoStrat(k)

                   AB = dy2 * pStrat(k)**2 / rhoEdge

                   ! ---------------------- A(i,j,k+1) ------------------

                   if (k == nz) then
                      AU=0.0
                     else
                      rhoEdge = 0.5* ( var(i,j,k+1,1) + var(i,j,k,1) )
                      rhoEdge = rhoEdge + rhoStratTilde(k)
   
                      pStratU = 0.5* ( pStrat(k+1) + pStrat(k) )
                      AU = dz2 * pStratU**2 / rhoEdge
                   end if

                   ! ----------------------- A(i,j,k-1) -----------------

                   if (k == 1) then
                      AD=0.0
                     else
                      rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j,k-1,1) )
                      rhoEdge = rhoEdge + rhoStratTilde(k-1)
                   
                      pStratD = 0.5* ( pStrat(k) + pStrat(k-1) )
                      AD = dz2 * pStratD**2 / rhoEdge
                   end if

                   if( pressureScaling ) then
                      stop'ERROR: no pressure scaling allowed'
                   end if

                   ! ----------------------- A(i,j,k) -------------------
                   !testb
                   !AL = 0.0
                   !AR = 0.0
                   !AB = 0.0
                   !AF = 0.0
                   !teste

                   AC = - AR - AL - AF - AB - AU - AD
 

                   ! ------------------ define matrix A -------------------

                   if (poissonSolverType == 'hypre') then
                      ! index_count_hypre 
                      ! = ( i * j * k * ne_hypre_i ) - ne_hypre_i + 1

                      index_count_hypre = i
                      index_count_hypre = index_count_hypre + ( (j-1)*nx )
                      index_count_hypre = index_count_hypre &
                                          + ( (k-1)*nx*ny )

                      index_count_hypre &
                      = ( index_count_hypre * ne_hypre_i ) &
                        - ne_hypre_i + 1

                      values_i(index_count_hypre   ) = AC
                      values_i(index_count_hypre+ 1) = AL
                      values_i(index_count_hypre+ 2) = AR
                      values_i(index_count_hypre+ 3) = AB
                      values_i(index_count_hypre+ 4) = AF
                      values_i(index_count_hypre+ 5) = AD
                      values_i(index_count_hypre+ 6) = AU
                      values_i(index_count_hypre+ 7) = 0.0
                      values_i(index_count_hypre+ 8) = 0.0
                      values_i(index_count_hypre+ 9) = 0.0
                      values_i(index_count_hypre+10) = 0.0


                      ! testb
                      !values_i(index_count_hypre+ 7) &
                      != sin(real(index_count_hypre+ 7))
                      !values_i(index_count_hypre+ 8) &
                      != sin(real(index_count_hypre+ 8))
                      !values_i(index_count_hypre+ 9) &
                      != sin(real(index_count_hypre+ 9))
                      !values_i(index_count_hypre+10) &
                      != sin(real(index_count_hypre+ 10))
                      ! teste

                      !testb
                      !print*,i,j,k,'AC =',values_i(index_count_hypre   )
                      !print*,i,j,k,'AL =',values_i(index_count_hypre+ 1)
                      !print*,i,j,k,'AR =',values_i(index_count_hypre+ 2)
                      !print*,i,j,k,'AB =',values_i(index_count_hypre+ 3)
                      !print*,i,j,k,'AF =',values_i(index_count_hypre+ 4)
                      !print*,i,j,k,'AD =',values_i(index_count_hypre+ 5)
                      !print*,i,j,k,'AU =',values_i(index_count_hypre+ 6)
                      !print*,i,j,k,'00 =',values_i(index_count_hypre+ 7)
                      !!print*,i,j,k,'00 =',values_i(index_count_hypre+ 8)
                      !print*,i,j,k,'00 =',values_i(index_count_hypre+ 9)
                      !print*,i,j,k,'00 =',values_i(index_count_hypre+10)
                      !teste
                     else if (poissonSolverType == 'bicgstab') then
                      ac_b(i,j,k) = AC

                      al_b(i,j,k) = AL
                      ar_b(i,j,k) = AR

                      ab_b(i,j,k) = AB
                      af_b(i,j,k) = AF

                      ad_b(i,j,k) = AD
                      au_b(i,j,k) = AU

                      alb_b(i,j,k) = 0.0
                      alf_b(i,j,k) = 0.0
                      arb_b(i,j,k) = 0.0
                      arf_b(i,j,k) = 0.0
                     else
                      stop'ERROR: val_PsIn expects hypre or bicgstab'
                   end if
                end do ! i_loop
             end do ! j_loop
          end do ! k_loop
         else
          do k = 1,nz
             do j = 1,ny
                do i = 1,nx
                   ! ------------------ A(i+1,j,k) ------------------

                   rhoEdge = 0.5*( var(i+1,j,k,1) + var(i,j,k,1) )
                   rhoEdge = rhoEdge + rhoStrat(k)

                   AR = dx2 * pStrat(k)**2/rhoEdge
   
                   ! ------------------- A(i-1,j,k) --------------------

                   rhoEdge = 0.5*( var(i,j,k,1) + var(i-1,j,k,1) )
                   rhoEdge = rhoEdge + rhoStrat(k)

                   AL = dx2 * pStrat(k)**2/rhoEdge

                   ! -------------------- A(i,j+1,k) ----------------------

                   rhoEdge = 0.5*( var(i,j+1,k,1) + var(i,j,k,1) )
                   rhoEdge = rhoEdge + rhoStrat(k)

                   AF = dy2 * pStrat(k)**2 / rhoEdge

                   ! --------------------- A(i,j-1,k) -------------------

                   rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j-1,k,1) )
                   rhoEdge = rhoEdge + rhoStrat(k)

                   AB = dy2 * pStrat(k)**2 / rhoEdge

                   ! ---------------------- A(i,j,k+1) ------------------

                   if (k == nz) then
                      AU=0.0
                     else
                      rhoEdge = 0.5* ( var(i,j,k+1,1) + var(i,j,k,1) )
                      rhoEdge = rhoEdge + rhoStratTilde(k)
   
                      pStratU = 0.5* ( pStrat(k+1) + pStrat(k) )
                      AU = dz2 * pStratU**2 / rhoEdge
                   end if

                   ! ----------------------- A(i,j,k-1) -----------------

                   if (k == 1) then
                      AD=0.0
                     else
                      rhoEdge = 0.5* ( var(i,j,k,1) + var(i,j,k-1,1) )
                      rhoEdge = rhoEdge + rhoStratTilde(k-1)
                   
                      pStratD = 0.5* ( pStrat(k) + pStrat(k-1) )
                      AD = dz2 * pStratD**2 / rhoEdge
                   end if

                   if( pressureScaling ) then
                      stop'ERROR: no pressure scaling allowed'
                   end if

                   ! ----------------------- A(i,j,k) -------------------

                   AC = - AR - AL - AF - AB - AU - AD
 

                   ! ------------------ define matrix A -------------------

                   if (poissonSolverType == 'hypre') then
                      ! index_count_hypre 
                      ! = ( i * j * k * ne_hypre_e ) - ne_hypre_e + 1

                      index_count_hypre = i
                      index_count_hypre = index_count_hypre + ( (j-1)*nx )
                      index_count_hypre = index_count_hypre &
                                          + ( (k-1)*nx*ny )

                      index_count_hypre &
                      = ( index_count_hypre * ne_hypre_e ) &
                        - ne_hypre_e + 1

                      values_e(index_count_hypre)   = AC
                      values_e(index_count_hypre+1) = AL
                      values_e(index_count_hypre+2) = AR
                      values_e(index_count_hypre+3) = AB
                      values_e(index_count_hypre+4) = AF
                      values_e(index_count_hypre+5) = AD
                      values_e(index_count_hypre+6) = AU

                      !testb
                      !print*,i,j,k,'AC =',values_e(index_count_hypre   )
                      !print*,i,j,k,'AL =',values_e(index_count_hypre+ 1)
                      !print*,i,j,k,'AR =',values_e(index_count_hypre+ 2)
                      !print*,i,j,k,'AB =',values_e(index_count_hypre+ 3)
                      !print*,i,j,k,'AF =',values_e(index_count_hypre+ 4)
                      !print*,i,j,k,'AD =',values_e(index_count_hypre+ 5)
                      !print*,i,j,k,'AU =',values_e(index_count_hypre+ 6)
                      !teste
                     else if (poissonSolverType == 'bicgstab') then
                      ac_b(i,j,k) = AC

                      al_b(i,j,k) = AL
                      ar_b(i,j,k) = AR

                      ab_b(i,j,k) = AB
                      af_b(i,j,k) = AF

                      ad_b(i,j,k) = AD
                      au_b(i,j,k) = AU
                     else
                      stop'ERROR: val_PsIn expects hypre or bicgstab'
                   end if
                end do ! i_loop
             end do ! j_loop
          end do ! k_loop
       end if
      else if (opt == "impl") then
       if (timeScheme /= "semiimplicit") then
          stop'ERROR: for opt = impl must have timeScheme = semiimplicit'
       end if

       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                AL = 0.0
                AR = 0.0
                AB = 0.0
                AF = 0.0
                AD = 0.0
                AU = 0.0
                AC = 0.0
                ALB = 0.0
                ALF = 0.0
                ARB = 0.0
                ARF = 0.0

                ! ------------------- from P UR/dx ------------------------

                facu = 1.0

                if (topography) then
                   if (   topography_mask(i0+i,j0+j,k) &
                     & .or. topography_mask(i0+i+1,j0+j,k)) then
                      facu = facu + dt*alprlx
                   end if
                end if

                facv = facu

                ! A(i+1,j,k) and A(i,j,k)

                rhoEdge = 0.5*( var(i+1,j,k,1) + var(i,j,k,1) )
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = dx2 * pStrat(k)**2/rhoEdge &
                         * facv/(facu*facv + (f_cor_nd(j)*dt)**2)

                AR = AR + acontr
                AC = AC - acontr

                ! A(i,j,k) and A(i,j-1,k)

                rhoEdge = 0.5*( var(i,j-1,k,1) + var(i,j,k,1) )
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = 0.25 * dxy * pStrat(k)**2/rhoEdge &
                         * f_cor_nd(j)*dt/(facu*facv + (f_cor_nd(j)*dt)**2)

                AC = AC + acontr
                AB = AB - acontr

                ! A(i,j,k) and A(i,j+1,k)

                rhoEdge = 0.5*( var(i,j,k,1) + var(i,j+1,k,1) )
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = 0.25 * dxy * pStrat(k)**2/rhoEdge &
                         * f_cor_nd(j)*dt/(facu*facv + (f_cor_nd(j)*dt)**2)

                AF = AF + acontr
                AC = AC - acontr

                ! A(i+1,j,k) and A(i+1,j-1,k)

                rhoEdge = 0.5*( var(i+1,j-1,k,1) + var(i+1,j,k,1) )
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = 0.25 * dxy * pStrat(k)**2/rhoEdge &
                         * f_cor_nd(j)*dt/(facu*facv + (f_cor_nd(j)*dt)**2)

                AR = AR + acontr
                ARB = ARB - acontr

                ! A(i+1,j,k) and A(i+1,j+1,k)

                rhoEdge = 0.5*( var(i+1,j,k,1) + var(i+1,j+1,k,1) )
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = 0.25 * dxy * pStrat(k)**2/rhoEdge &
                         * f_cor_nd(j)*dt/(facu*facv + (f_cor_nd(j)*dt)**2)

                ARF = ARF + acontr
                AR = AR - acontr

                ! ------------------- from - P UL/dx --------------------

                facu = 1.0

                if (topography) then
                   if (   topography_mask(i0+i-1,j0+j,k) &
                     & .or. topography_mask(i0+i,j0+j,k)) then
                      facu = facu + dt*alprlx
                   end if
                end if

                facv = facu

                ! A(i,j,k) and A(i-1,j,k)

                rhoEdge = 0.5*( var(i,j,k,1) + var(i-1,j,k,1) )
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = - dx2 * pStrat(k)**2/rhoEdge &
                           * facv/(facu*facv + (f_cor_nd(j)*dt)**2)

                AC = AC + acontr
                AL = AL - acontr

                ! A(i-1,j,k) and A(i-1,j-1,k)

                rhoEdge = 0.5*( var(i-1,j-1,k,1) + var(i-1,j,k,1) )
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = -0.25 * dxy * pStrat(k)**2/rhoEdge &
                          * f_cor_nd(j)*dt/(facu*facv + (f_cor_nd(j)*dt)**2)

                AL = AL + acontr
                ALB = ALB - acontr

                ! A(i-1,j,k) and A(i-1,j+1,k)

                rhoEdge = 0.5*( var(i-1,j,k,1) + var(i-1,j+1,k,1) )
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = - 0.25 * dxy * pStrat(k)**2/rhoEdge &
                           * f_cor_nd(j)*dt/(facu*facv + (f_cor_nd(j)*dt)**2)

                ALF = ALF + acontr
                AL = AL - acontr

                ! A(i,j,k) and A(i,j-1,k)

                rhoEdge = 0.5*( var(i,j-1,k,1) + var(i,j,k,1) )
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = - 0.25 * dxy * pStrat(k)**2/rhoEdge &
                           * f_cor_nd(j)*dt/(facu*facv + (f_cor_nd(j)*dt)**2)

                AC = AC + acontr
                AB = AB - acontr

                ! A(i,j,k) and A(i,j+1,k)

                rhoEdge = 0.5*( var(i,j,k,1) + var(i,j+1,k,1) )
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = - 0.25 * dxy * pStrat(k)**2/rhoEdge &
                           * f_cor_nd(j)*dt/(facu*facv + (f_cor_nd(j)*dt)**2)

                AF = AF + acontr
                AC = AC - acontr

                ! ------------------- from P VF/dy ------------------------

                facv = 1.0

                if (topography) then
                   if (   topography_mask(i0+i,j0+j,k) &
                     & .or. topography_mask(i0+i,j0+j+1,k)) then
                      facv = facv + dt*alprlx
                   end if
                end if

                facu = facv

                ! A(i+1,j,k) and A(i,j,k)

                rhoEdge = 0.5 * (var(i,j,k,1) + var(i+1,j,k,1))
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = 0.25 * dxy * pStrat(k)**2/rhoEdge &
                         * (-f_cor_nd(j)*dt)/(facu*facv + (f_cor_nd(j)*dt)**2)

                AR = AR + acontr
                AC = AC - acontr

                ! A(i,j,k) and A(i-1,j,k)

                rhoEdge = 0.5 * (var(i-1,j,k,1) + var(i,j,k,1))
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = 0.25 * dxy * pStrat(k)**2/rhoEdge &
                         * (-f_cor_nd(j)*dt)/(facu*facv + (f_cor_nd(j)*dt)**2)

                AC = AC + acontr
                AL = AL - acontr

                ! A(i+1,j+1,k) and A(i,j+1,k)

                rhoEdge = 0.5 * (var(i,j+1,k,1) + var(i+1,j+1,k,1))
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = 0.25 * dxy * pStrat(k)**2/rhoEdge &
                         * (-f_cor_nd(j)*dt)/(facu*facv + (f_cor_nd(j)*dt)**2)

                ARF = ARF + acontr
                AF = AF - acontr

                ! A(i,j+1,k) and A(i-1,j+1,k)

                rhoEdge = 0.5 * (var(i-1,j+1,k,1) + var(i,j+1,k,1))
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = 0.25 * dxy * pStrat(k)**2/rhoEdge &
                         * (-f_cor_nd(j)*dt)/(facu*facv + (f_cor_nd(j)*dt)**2)

                AF = AF + acontr
                ALF = ALF - acontr

                ! A(i,j+1,k) and A(i,j,k)

                rhoEdge = 0.5 * (var(i,j+1,k,1) + var(i,j,k,1))
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = dy2 * pStrat(k)**2/rhoEdge &
                         * facu/(facu*facv + (f_cor_nd(j)*dt)**2)

                AF = AF + acontr
                AC = AC - acontr

                ! ------------------- from - P VB/dy ---------------------

                facv = 1.0

                if (topography) then
                   if (   topography_mask(i0+i,j0+j-1,k) &
                     & .or. topography_mask(i0+i,j0+j,k)) then
                      facv = facv + dt*alprlx
                   end if
                end if

                facu = facv

                ! A(i+1,j-1,k) and A(i,j-1,k)

                rhoEdge = 0.5 * (var(i,j-1,k,1) + var(i+1,j-1,k,1))
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = - 0.25 * dxy * pStrat(k)**2/rhoEdge &
                          * (-f_cor_nd(j)*dt)/(facu*facv + (f_cor_nd(j)*dt)**2)

                ARB = ARB + acontr
                AB = AB - acontr

                ! A(i,j-1,k) and A(i-1,j-1,k)

                rhoEdge = 0.5 * (var(i-1,j-1,k,1) + var(i,j-1,k,1))
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = - 0.25 * dxy * pStrat(k)**2/rhoEdge &
                           * (-f_cor_nd(j)*dt)/(facu*facv + (f_cor_nd(j)*dt)**2)

                AB = AB + acontr
                ALB = ALB - acontr

                ! A(i+1,j,k) and A(i,j,k)

                rhoEdge = 0.5 * (var(i,j,k,1) + var(i+1,j,k,1))
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = - 0.25 * dxy * pStrat(k)**2/rhoEdge &
                           * (-f_cor_nd(j)*dt)/(facu*facv + (f_cor_nd(j)*dt)**2)

                AR = AR + acontr
                AC = AC - acontr

                ! A(i,j,k) and A(i-1,j,k)

                rhoEdge = 0.5 * (var(i-1,j,k,1) + var(i,j,k,1))
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = - 0.25 * dxy * pStrat(k)**2/rhoEdge &
                           * (-f_cor_nd(j)*dt)/(facu*facv + (f_cor_nd(j)*dt)**2)

                AC = AC + acontr
                AL = AL - acontr

                ! A(i,j,k) and A(i,j-1,k)

                rhoEdge = 0.5 * (var(i,j,k,1) + var(i,j-1,k,1))
                rhoEdge = rhoEdge + rhoStrat(k)

                acontr = - dy2 * pStrat(k)**2/rhoEdge &
                           * facu/(facu*facv + (f_cor_nd(j)*dt)**2)

                AC = AC + acontr
                AB = AB - acontr

                ! ------------------- from PU WU/dz ---------------------

                if (k == nz) then
                   AU = 0.0
                  else
                   facw = 1.0
                   facr = 1.0 + alprlx*dt

                   if (topography) then
                      if (   topography_mask(i0+i,j0+j,k) &
                        & .or. topography_mask(i0+i,j0+j,k+1)) then
                         facw = facw + dt*alprlx
                      end if
                   end if

                   bvsstw = 0.5 * (bvsStrat(k) + bvsStrat(k+1))

                   ! A(i,j,k+1) and A(i,j,k)

                   rhoEdge = 0.5*( var(i,j,k+1,1) + var(i,j,k,1) )
                   rhoEdge = rhoEdge + rhoStratTilde(k)
   
                   pStratU = 0.5 * ( pStrat(k+1) + pStrat(k) )

                   AU = dz2 * pStratU**2 / rhoEdge &
                        * facr &
                          /(  facr*facw &
                            + rhoStratTilde(k)/rhoEdge * bvsstw * dt**2)
   
                   AC = AC - AU
                end if

                ! ------------------- from - PD WD/dz ---------------------

                if (k == 1) then
                   AD = 0.0
                  else
                   facw = 1.0
                   facr = 1.0 + alprlx*dt

                   if (topography) then
                      if (   topography_mask(i0+i,j0+j,k-1) &
                        & .or. topography_mask(i0+i,j0+j,k)) then
                         facw = facw + dt*alprlx
                      end if
                   end if

                   bvsstw = 0.5 * (bvsStrat(k-1) + bvsStrat(k))

                   ! A(i,j,k) and A(i,j,k-1)

                   rhoEdge = 0.5*( var(i,j,k,1) + var(i,j,k-1,1) )
                   rhoEdge = rhoEdge + rhoStratTilde(k-1)
   
                   pStratD = 0.5* ( pStrat(k) + pStrat(k-1) )

                   AD = dz2 * pStratD**2 / rhoEdge &
                        * facr &
                          /(  facr*facw &
                            + rhoStratTilde(k-1)/rhoEdge * bvsstw * dt**2)

                   AC = AC - AD
                end if

                if( pressureScaling ) then
                   stop'ERROR: no pressure scaling allowed'
                end if


                ! ------------------- define matrix A -------------------

                if (poissonSolverType == 'hypre') then
                   ! index_count_hypre 
                   ! = ( i * j * k * ne_hypre_i ) - ne_hypre_i + 1

                   index_count_hypre = i
                   index_count_hypre = index_count_hypre + ( (j-1)*nx )
                   index_count_hypre = index_count_hypre + ( (k-1)*nx*ny )

                   index_count_hypre &
                   = ( index_count_hypre * ne_hypre_i ) &
                     - ne_hypre_i + 1

                   values_i(index_count_hypre   ) = AC
                   values_i(index_count_hypre+ 1) = AL
                   values_i(index_count_hypre+ 2) = AR
                   values_i(index_count_hypre+ 3) = AB
                   values_i(index_count_hypre+ 4) = AF
                   values_i(index_count_hypre+ 5) = AD
                   values_i(index_count_hypre+ 6) = AU
                   values_i(index_count_hypre+ 7) = ALB
                   values_i(index_count_hypre+ 8) = ALF
                   values_i(index_count_hypre+ 9) = ARB
                   values_i(index_count_hypre+10) = ARF
                  else if (poissonSolverType == 'bicgstab') then
                   ac_b(i,j,k) = AC

                   al_b(i,j,k) = AL
                   ar_b(i,j,k) = AR

                   ab_b(i,j,k) = AB
                   af_b(i,j,k) = AF

                   ad_b(i,j,k) = AD
                   au_b(i,j,k) = AU

                   alb_b(i,j,k) = ALB
                   alf_b(i,j,k) = ALF
                   arb_b(i,j,k) = ARB
                   arf_b(i,j,k) = ARF
                  else
                   stop'ERROR: val_PsIn expects hypre or bicgstab'
                end if
             end do ! i_loop
          end do ! j_loop
       end do ! k_loop
      else
       stop'ERROR: wrong opt'
    end if

  return

  end subroutine val_PsIn
  
  !============================================================== 

  subroutine  val_hypre_Bous
         
    ! local variables
    integer :: i,j,k
    real :: pStratU, pStratD, rhoEdge
    real :: AL,AR, AB,AF, AD,AU, AC
    real :: dx2, dy2, dz2 
    integer :: index_count_hypre

    if (topography) then
       print*,'ERROR: no topography allowed in Boussinesq mode'
       print*,'(would require semi-implicit time stepping)'
       print*,'(could be implemented easily by simplifying the &
             & pseudo-incompressible case accordingly)'
       stop
    end if

    ! auxiliary variables
    dx2 = 1.0/dx**2
    dy2 = 1.0/dy**2
    dz2 = 1.0/dz**2    

    !---------------------------------
    !         Loop over field
    !---------------------------------

    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
 
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

             ! ------------------- define matrix A --------------------

             ! index_count_hypre &
             ! = ( i * j * k * ne_hypre_e ) - ne_hypre_e + 1

             index_count_hypre = i
             index_count_hypre = index_count_hypre + ( (j-1)*nx )
             index_count_hypre = index_count_hypre + ( (k-1)*nx*ny )

             index_count_hypre &
             = ( index_count_hypre * ne_hypre_e ) &
               - ne_hypre_e + 1

             values_e(index_count_hypre)   = AC
             values_e(index_count_hypre+1) = AL
             values_e(index_count_hypre+2) = AR
             values_e(index_count_hypre+3) = AB
             values_e(index_count_hypre+4) = AF
             values_e(index_count_hypre+5) = AD
             values_e(index_count_hypre+6) = AU
          end do
       end do
    end do

  return

  end subroutine val_hypre_Bous

!------------------------------------------------------------------

  subroutine heat_w0(var,flux,heat,S_bar,w_0)


  ! negative (!) heating, its horizontal mean, 
  ! and the thereby induced vertical wind

  ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var
    !UAB 200413
    real, dimension(-1:nx,-1:ny,-1:nz,3,nVar), intent(in) :: flux
    !UAE 200413
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz), &
         & intent(out) :: heat
    real, dimension(-nbz:nz+nbz),intent(out) :: w_0 
    real, dimension(-nbz:nz+nbz), intent(out) :: S_bar

    integer :: i,j,k
    real, dimension(1:nz) :: sum_local, sum_global 
    real, dimension(-nbz:nz+nbz) :: press0 

    real :: dptopdt

    real :: sum_d, sum_n

    w_0 = 0.
    S_bar = 0.
    heat = 0.

    ! negative (!) heating, i.e. -S eq(9)  ONeill+Klein2014
    call calculate_heating(var,flux,heat) 

    ! calculate horizontal mean of heat(:,:,:)
    do k = 1,nz
       sum_local(k) = sum(heat(1:nx,1:ny,k))
    end do
    !global sum and average
    call mpi_allreduce(sum_local(1),sum_global(1),&
         nz-1+1,&
         mpi_double_precision,mpi_sum,comm,ierror)
    sum_global = sum_global/(sizeX*sizeY)

    S_bar(1:nz) = sum_global(1:nz)

    !non_dim. pressure; eq(12) ONeill+Klein2014

    do k = 1,nz
       press0(k) = PStrat(k)**gamma  
    end do
    
    !  eq(B.14) ONeill+Klein2014

    sum_d = 0.0
    sum_n = 0.0

    do k = 1,nz
      sum_n = sum_n - S_bar(k)/PStrat(k)
      sum_d = sum_d +  1./(gamma*press0(k))
    end do

    dptopdt = sum_n/sum_d

    w_0(1) = dz*(- S_bar(1)/Pstrat(1) - (1./(gamma*press0(1)))*dptopdt)

    do k = 2,nz-1
      w_0(k) &
      = w_0(k-1) &
        + dz*(- S_bar(k)/Pstrat(k) - (1./(gamma*press0(k)))*dptopdt)
    end do


  end subroutine heat_w0

end module poisson_module
