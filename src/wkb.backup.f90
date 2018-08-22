module wkb_module
  !-----------------------------------------------------------
  ! This module is rather independent of the pinc floit code
  ! it integrates the leading order terms of the 
  ! non-linear WKB Theory by Achatz et al., i.e.
  ! 1) the induced mean flow
  ! 2) the first and second harmonics
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  !
  !                    Known problems
  !
  ! 1) no rays in staggered cell for cgx or cgz calculation
  ! -> increase nb. of rays / decrease u0-jet
  ! 2) rays end up with z < 0
  ! -> reason: time scheme allows negative movement
  ! -> use Euler / do not allow 
  !
  !------------------------------------------------------------
  
  use type_module
  use timeScheme_module
  use atmosphere_module
  use muscl_module


  implicit none

  private                        ! all module variables are internal to the module


  !----------------------
  !   public routines
  !----------------------
  public :: setup_wkb
  public :: setup_wkb_dummy
  public :: transport_ray
  public :: calc_cellIndex
  public :: setBoundary_wkb
  public :: transport_waveAction
  public :: finish_wkb
  public :: finish_wkb_dummy
  public :: calc_waveAmplitude
  public :: calc_meanFlow
  
  public :: cabs                     ! calc absolute value of complex numbers
  


  !----------------------
  !   private routines
  !----------------------
  private :: reconstruct_waveAction
  private :: calc_waveFlux
  private :: add_waveFlux
  private :: cphase                   ! calc phase of complex numbers in radians

  
  !------------------------------
  !   private module variables 
  !------------------------------
  real :: lambdaX, lambdaZ        ! nondimensional wave length

  ! displacement vector of rays
  real, dimension(:,:), allocatable :: dxRay, dkRay

  ! reconstructed wave amplitude
  real, dimension(:,:,:,:,:), allocatable :: waveActTilde ! reconstructed wave action
  !              (i,j,k,direction: 1=x, 2=y, 3=z, position: 0=left, 1=right)

  ! wave action flux 
  real, dimension(:,:,:,:), allocatable  :: waveFlux
  
  ! counter of rays in grid cells
  integer, dimension(:,:,:), allocatable :: nRayPerCell

  ! Runge-Kutta tendency of wave action
  real, dimension(:,:,:), allocatable, public :: delWaveAct

  ! Runge-Kutta tendency of mean flow U00
  complex, dimension(:,:,:), allocatable, private :: delU00
  
  
  integer :: i,j,k
  integer :: iRay                  ! index of Ray
  integer :: iCell, jCell, kCell   ! index of FV grid cells
  real    :: xRay, zRay            ! non-dim position of ray
  real    :: kk0, mm0              ! non-dim initial wave numbers
  integer :: nxRay, nyRay, nzRay   ! nb. of rays along x,y,z,

  real :: i_dx, del_i              ! real position within grid / cell
  real :: k_dz, del_k
  character(len=10) :: left_right_pos, up_down_pos
  logical, parameter :: debugging = .false.  

  
  
contains

  subroutine calc_meanFlow(Psi,var,dt,RKstage)
    ! in/out variables
    complex, dimension(0:nx+1,0:ny+1,0:nz+1,4,0:2), intent(inout) :: Psi
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    
    real, intent(in) :: dt
    integer, intent(in) :: RKstage
    
    ! local variables
    complex :: u10_t, u10_b, w10_t, w10_b
    real    :: u00, u, delU
    real    :: rho0_t, rho0_b, rho0_c
    real    :: F, d_dz

!!$    ! mean value calculation
!!$    real :: rho0, rho0_t, rho0_b, d_dz, ypsi
!!$
!!$    ! mean value calculation for 2D GWP
!!$    complex :: w10_l_conj, w10_r_conj
!!$    real    :: d_dx


    ! init tendency
    if( RKStage == 1) delU00 = 0.0
    
    
    do k = 1,nz
       j = 1
       do i = 1,nx
          
          ! calc RHS F of ODE dphi/dt = F
          
          ! values needed for d/dz
          u10_t = Psi(i,j,k+1,1,1)
          u10_b = Psi(i,j,k-1,1,1)
          w10_t = Psi(i,j,k+1,2,1)
          w10_b = Psi(i,j,k-1,2,1)

          u00 = Psi(i,j,k,1,0)
          
          rho0_t = rhoStrat(k+1)
          rho0_b = rhoStrat(k-1)
          rho0_c = rhoStrat(k)

          d_dz = 0.5/dz*( rho0_t*real(u10_t*conjg(w10_t)) &
               &        - rho0_b*real(u10_b*conjg(w10_b)) )
                    
          F = -0.5/rho0_c * d_dz

          ! update tendency and stage value

          delU00(i,j,k) = dt*F + alpha(RKstage)*delU00(i,j,k)
          u00 = u00 + beta(RKstage)*delU00(i,j,k)

          Psi(i,j,k,1,0) = u00

         

       end do
    end do

    !------------------------------
    !xxx update of velocity field
    !------------------------------
    
    ! set boundary values
    delU00(0,:,:) = delU00(nx,:,:)
    delU00(nx+1,:,:) = delU00(1,:,:)
    
    
    do k = 1,nz
       j = 1
       do i = 0, nx
          delU = 0.5*( delU00(i+1,j,k) + delU00(i,j,k) )
          u    = var(i,j,k,2)
          u    = u + beta(RKstage)*delU

          ! save value
          var(i,j,k,2) = u
       end do
    end do


!!$    !---------------------------------------
!!$    !             Calc mean values
!!$    !     (zeroth harmonic, leading order)
!!$    !---------------------------------------
!!$    
!!$    do k = 1,nz
!!$       j = 1
!!$       do i = 1,nx
!!$
!!$          b11 = Psi(i,j,k,3,1)
!!$          rho0 = rhoStrat(k)
!!$
!!$          rho0_t = rhoStrat(k+1)
!!$          rho0_b = rhoStrat(k-1)
!!$          w10_t = Psi(i,j,k+1,2,1)
!!$          w10_b = Psi(i,j,k-1,2,1)
!!$          
!!$          d_dz = 0.5*(rho0_t*cabs(w10_t)**2 - rho0_b*cabs(w10_b)**2) / dz
!!$
!!$          ! 2D wave packet -> include d_dx
!!$          u10_l = Psi(i-1,j,k,1,1)
!!$          u10_r = Psi(i+1,j,k,1,1)
!!$          w10_l_conj = conjg(Psi(i-1,j,k,2,1))
!!$          w10_r_conj = conjg(Psi(i+1,j,k,2,1))
!!$          d_dx = rho0 * 0.5*(real(u10_r*w10_r_conj) - real(u10_l*w10_l_conj)) / dz
!!$          
!!$          ypsi = -0.5*( Fr2*cabs(b11)**2 + 1./rho0 * (d_dx + d_dz) )
!!$          
!!$          ! save sum of pi02-term + b02-term in place of b02
!!$          Psi(i,j,k,3,0) = ypsi
!!$          
!!$          
!!$       end do
!!$    end do




  end subroutine calc_meanFlow


!------------------------------------------------------------------------------

  
  subroutine calc_waveAmplitude(ray,waveAct,Psi)
    
    !------------------------------------------------
    !  calculate complex amplitudes for
    !    1) first harmonics, leading order: Psi(:,:,:,0)
    !    2) second harmonics, leading order: Psi(:,:,:,1)
    !------------------------------------------------
    
    ! in/out variables
    !real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
    !     & intent(in) :: var                         ! mean flow velocities    
    
    type(rayType), dimension(nRay), intent(in)           :: ray
    real, dimension(0:nx+1,0:ny+1,0:nz+1), intent(inout) :: waveAct 
    ! wave amplitude
    complex, dimension(0:nx+1,0:ny+1,0:nz+1,4,0:2), intent(out) :: Psi 
    
    ! local variables
    real :: A, rho
    real :: kk, mm, omi
    real :: kk2,mm2,omi2, kTot2
    complex :: u10,w10,b11,pi12
    integer :: iRay
    integer :: i,j,k
    integer :: ii,jj

    ! local amplitude values and derivatives
    complex :: u10_r, u10_c, u10_l, u10_t, u10_b
    complex :: w10_r, w10_c, w10_l, w10_t, w10_b
    complex :: theta11_r, theta11_c, theta11_l, theta11_t, theta11_b
    complex :: pi12_c 
    complex :: pi0_t, pi0_c, pi0_b
    complex :: theta0_c
    complex :: du10_dx, du10_dz, dw10_dx, dw10_dz 
    complex :: dpi0_dz, dtheta11_dx, dtheta11_dz
    complex :: Div
    complex :: Press, PressU, PressW
    complex :: d1u10, d1w10, d1theta11
    complex :: coeff, aux1
    complex :: M11,M12,M13,M14
    complex :: M21,M22,M23,M24
    complex :: M31,M32,M33,M34
    complex :: M41,M42,M43,M44
    complex, dimension(4) :: RHS, Sol
    complex, dimension(4,4) :: M2, M2inv
    complex :: u21, w21,b22, pi23


    ! debugging stuff
    complex :: summe
    complex, dimension(4) :: term
    
    ! local fields
    real, dimension(nx,ny,nz) :: omegaMean   ! cell averaged intrinsic freq.
    real, dimension(nx,ny,nz) :: kMean, mMean ! cell mean wave numbers

    ! more debugging stuff
    real :: B11_pinc
    real :: D1TH11
    

    !---------------------------------------
    !   calc intrinsic frequency cell mean
    !---------------------------------------

    ! init mean fields
    omegaMean = 0.0
    kMean = 0.0
    mMean = 0.0
    nRayPerCell = 0
    
    ! loop over all rays
    rayloop: do iRay = 1, nRay
       i = ray(iRay)%iCell
       j = ray(iRay)%jCell
       k = ray(iRay)%kCell
     
       omegaMean(i,j,k) = omegaMean(i,j,k) + ray(iRay)%omega
       kMean(i,j,k) = kMean(i,j,k) + ray(iRay)%k
       mMean(i,j,k) = mMean(i,j,k) + ray(iRay)%m

       nRayPerCell(i,j,k) = nRayPerCell(i,j,k) + 1       
    
    end do rayloop

    ! average over cells
    do k = 1,nz
       j = 1
       do i = 1,nx
          
          omegaMean(i,j,k) = omegaMean(i,j,k) / nRayPerCell(i,j,k)
          kMean(i,j,k) = kMean(i,j,k) / nRayPerCell(i,j,k)
          mMean(i,j,k) = mMean(i,j,k) / nRayPerCell(i,j,k)

       end do
    end do
    
    
    !---------------------------------------
    !        calc amplitude Psi_1^0 
    !     (first harmonic, leading order)
    !---------------------------------------

    do k = 1,nz
       j = 1
       do i = 1,nx
          
          ! buoyancy from wave action
          A = waveAct(i,j,k)
          omi = omegaMean(i,j,k)
          rho = rhoStrat(k)
          b11 = cmplx( sqrt(2.0*N2*omi/rho*A), 0.0 )
          
          ! velocities, pressure
          mm = mMean(i,j,k)
          kk = kMean(i,j,k)
          theta0 = thetaStrat(k)

          u10 = cmplx(0.0, -mm/kk * omi/N2) * b11
          
          w10 = cmplx(0.0, omi/N2) * b11

          pi12 = cmplx(0.0, -kappa*Ma2* mm/kk**2 * omi**2/N2 / theta0) * b11

          Psi(i,j,k,:,1) = (/u10, w10, b11, pi12/) 


       end do
    end do



    !---------------------------------------
    !    set ghost cell values for Psi_1^0 
    !    x/y: Periodic, z: solid wall 
    !---------------------------------------
    
    ! periodic in x
    Psi(0,:,:,:,1) = Psi(nx,:,:,:,1)
    Psi(nx+1,:,:,:,1) = Psi(1,:,:,:,1)
    
    ! periodic in y
    ! implement for 3D

    ! solid wall -> reflect u10 and w10 with change of sign
    Psi(:,:,nz+1,2,1) = -Psi(:,:,nz,2,1)
    Psi(:,:,0,2,1) = -Psi(:,:,1,2,1)
    
    ! solid wall: reflect b11 and pi12 without change of sign 
    Psi(:,:,nz+1,3:4,1) = Psi(:,:,nz,3:4,1)
    Psi(:,:,0,3:4,1) = Psi(:,:,1,3:4,1)
    
    

    !---------------------------------------
    !        calc amplitude Psi_2^1 
    !     (second harmonic, first order)
    !---------------------------------------
    
    do k = 1,nz
       j = 1
       do i = 1,nx

          !-------------------------
          !       Set up RHS
          !-------------------------
          
          ! averaged wavenumbers and intrinsic freqencies
          kk = kMean(i,j,k)
          kk2 = kk**2
          mm = mMean(i,j,k)
          mm2 = mm**2
          omi = omegaMean(i,j,k)
          omi2 = omi**2
          kTot2 = kk2 + mm2


          ! zonal velocities right, center, left, top, bottom
          u10_r = Psi(i+1,j, k ,1,1)
          u10_c = Psi( i ,j, k ,1,1)
          u10_l = Psi(i-1,j, k ,1,1)
          u10_t = Psi( i ,j,k+1,1,1)
          u10_b = Psi( i ,j,k-1,1,1)

          ! vertical velocities top, center, bottom
          w10_r = Psi(i+1,j, k ,2,1)
          w10_l = Psi(i-1,j, k ,2,1)          
          w10_t = Psi( i ,j,k+1,2,1)
          w10_c = Psi( i ,j, k ,2,1)
          w10_b = Psi( i ,j,k-1,2,1)
          

          ! Buoyancy and pot. temp.
          theta11_r = Fr2 * thetaStrat(k) * Psi(i+1,j,k,3,1)
          theta11_c = Fr2 * thetaStrat(k) * Psi(i,j,k,3,1)
          theta11_l = Fr2 * thetaStrat(k) * Psi(i-1,j,k,3,1)
          theta11_t = Fr2 * thetaStrat(k+1) * Psi(i,j,k+1,3,1)
          theta11_b = Fr2 * thetaStrat(k-1) * Psi(i,j,k-1,3,1)

          ! Second order Exner pressure
          pi12_c = Psi(i,j,k,4,1)

          ! Background Exner pressure and pot. temp.
          pi0_t = (Pstrat(k+1)/p0)**gamma_1
          pi0_c = (Pstrat(k)/p0)**gamma_1
          pi0_b = (Pstrat(k-1)/p0)**gamma_1
          theta0_c = thetaStrat(k)

         
          ! derivatives          
          du10_dx = 0.5*(u10_r - u10_l)/dx
          du10_dz = 0.5*(u10_t - u10_b)/dz
          dw10_dx = 0.5*(w10_r - w10_l)/dx
          dw10_dz = 0.5*(w10_t - w10_b)/dz
          dpi0_dz = 0.5*(pi0_t - pi0_b)/dz
          dtheta11_dx = 0.5*(theta11_r - theta11_l)/dx
          dtheta11_dz = 0.5*(theta11_t - theta11_b)/dz

          ! divergence term -> Eq. (7.21)
          Div = -du10_dx - dw10_dz - (1.-kappa)/kappa*w10_c/pi0_c*dpi0_dz
          
          ! Pressure terms
          Press = 0.5*kappaInv*MaInv2*imag*theta11_c*pi12_c
          PressU = kk*Press
          PressW = mm*Press
          

          ! intermediate terms
          d1u10 = 0.5*(u10_c*du10_dx + w10_c*du10_dz + Div*u10_c)
          d1w10 = 0.5*(u10_c*dw10_dx + w10_c*dw10_dz + Div*w10_c)
          d1theta11 = 0.5*(u10_c*dtheta11_dx + w10_c*dtheta11_dz + Div*theta11_c)

          RHS(1) = -d1u10 - PressU
          RHS(2) = -d1w10 - PressW
          RHS(3) = -FrInv2/NN/theta0_c * d1theta11
          RHS(4) = (0.0, 0.0)
          

!-------- debugging wave 2 amplitude
!!$if( k == nz/2 ) then
!!$   !B11_pinc = cabs(Psi(i,j,k,3,1))*lRef/tRef**2
!!$   !print"(a,es9.1)","B11 = ",B11_pinc
!!$   D1TH11 = cabs(-FrInv2/NN/theta0_c * d1theta11 ) * lRef/tRef**2
!!$   print"(a,es10.3)","Div = ", cabs(Div/tRef)
!!$   print"(a,es10.3)","u10_c*dtheta11_dx = ", cabs(u10_c*dtheta11_dx * thetaRef/tRef)
!!$   print"(a,es10.3)","w10_c*dtheta11_dz = ", cabs(w10_c*dtheta11_dz * thetaRef/tRef)
!!$   print"(a,es10.3)","d1theta11 = ", cabs(d1theta11) * thetaRef/tRef
!!$   print"(a,es10.3)","-FrInv2/NN/theta0_c * d1theta11 = ", &
!!$        & cabs(-FrInv2/NN/theta0_c * d1theta11) * lRef/tRef**2
!!$     
!!$end if
          
          !----------------------------------------------------
          !       Set up inverted system matrix M(2om,2k,2m)
          !----------------------------------------------------

          coeff = 1./(4.*omi2*kTot2 - kk2*N2)
          aux1 = 4.*omi2-N2
          
          M11 = 2.*imag*mm2*omi; M12 = -2.*imag*kk*mm*omi; M13 = kk*mm*NN; M14 = -0.5*imag*kk*aux1
          M21 = M12;           M22 = 2.*imag*kk2*omi;    M23 = -kk2*NN;  M24 = -2.*imag*omi2*mm
          M31 = -M13;          M32 = -M23;             M33 = 2.*imag*omi*kTot2; M34 = -omi*mm*NN
          M41 = M14;           M42 = M24;              M43 = -M34;     M44 = -0.5*imag*omi*aux1


          ! inverted matrix: check ok
          M2inv(1,:) = (/M11,M12,M13,M14/)
          M2inv(2,:) = (/M21,M22,M23,M24/)
          M2inv(3,:) = (/M31,M32,M33,M34/)
          M2inv(4,:) = (/M41,M42,M43,M44/)
          M2inv = coeff*M2inv


          !---------------------------------------
          !   Solve linear System -> save in Psi
          !---------------------------------------
          
          sol = matmul(M2inv,RHS)
          
          u21 = sol(1)
          w21 = sol(2)
          b22 = sol(3) * NN
          pi23 = sol(4) * kappa*Ma2 / theta0_c
          
          Psi(i,j,k,:,2) = (/u21,w21,b22,pi23/)


          
          
          ! analyse order of magnitude for terms leading to b22
          if( i == nx/2 .and. k == 362) then
             do ii = 1,4
                term(ii) = M2inv(3,ii) * RHS(ii) * NN
             end do
             summe = sum(term)
             print*,"k = ", k
             print "(6a9)","term1", "term2", "term3", "term4", "summe", "b22"
             print "(6es9.1)",(cabs(term(jj))*lRef/tRef**2, jj=1,4), &
                  & cabs(summe*lRef/tRef**2), cabs(b22*lRef/tRef**2)

             print "(3a9)", "M31", "M32", "M33"
             print "(3es9.1)", (cabs(M2inv(3,ii))*tRef, ii = 1,3)

             print "(3a9)", "RHS11", "RHS2", "RHS3"
             print "(3es9.1)", (cabs(RHS(ii))*lRef/tRef**2, ii = 1,3)
             
          end if
          

       end do
    end do

    
  end subroutine calc_waveAmplitude


  
  !---------------------------------------------------------------------------
  

  subroutine transport_waveAction(var,ray,waveAct,RKstage,dt)
    
    !-----------------------------------------------
    !  compound subroutine calling subroutines for:
    !  1) set boundary conditions
    !  2) reconstruct wave action in cells
    !  3) calc wave action fluxes
    !  4) update wave action
    !-----------------------------------------------
    
    ! in/out variables
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var                           ! mean flow velocities    
    
    type(rayType), dimension(nRay), intent(in)           :: ray
    real, dimension(0:nx+1,0:ny+1,0:nz+1), intent(inout) :: waveAct ! wave action
    
    integer, intent(in) :: RKstage
    real, intent(in)    :: dt
    
    
    !------------------------------------------------
    !   set boundary conditions for reconstruction
    !------------------------------------------------

    call setBoundary_wkb(waveAct, 'meanValue')


    !-----------------------------
    !   reconstruct wave action
    !-----------------------------

    call reconstruct_waveAction(waveAct)


    
    !------------------------------------------------
    !   set boundary conditions for flux calculation
    !------------------------------------------------

    call setBoundary_wkb(waveAct, 'reconstValue' )

   
    
    !-----------------------------
    !   calc wave action fluxes
    !-----------------------------
    
    call calc_waveFlux(var,ray,waveAct)



    !---------------------------
    !    add wave fluxes
    !---------------------------
    
    call add_waveFlux(waveAct,RKstage,dt)
    

    
  end subroutine transport_waveAction
  
  
  !-------------------------------------------------------------------------------------
  
  
  subroutine add_waveFlux(waveAct, RKstage, dt)
    
    !-----------------------------------------------
    !        Add wave fluxes to wave action
    !-----------------------------------------------
    
    ! in/out fields
    real, dimension(0:nx+1,0:ny+1,0:nz+1), intent(inout) :: waveAct ! wave action
    integer, intent(in) :: RKstage
    real, intent(in)    :: dt
    
    ! local variables
    integer :: i,j,k
    real    :: fR, fL, hU, hD         ! local wave action fluxes
    real    :: fluxDiff               ! flux difference
    real    :: F                      ! Runge-Kutta RHS 
    
    
    ! init RK tendency
    if( RKstage == 1)  delWaveAct = 0.0
    

    !------------------------------------
    !    Make Runge-Kutta flux update
    !------------------------------------

    do k = 1,nz
       j = 1     ! change for 3D
       do i = 1,nx
          
          ! fluxes
          fR = waveFlux(i,j,k,1)
          fL = waveFlux(i-1,j,k,1)
          hU = waveFlux(i,j,k,3)
          hD = waveFlux(i,j,k-1,3)
          
          ! divergence
          fluxDiff = (fR-fL)/dx + (hU-hD)/dz

          
          ! Runge-Kutta 
          F = -fluxDiff
          
          delWaveAct(i,j,k) = dt*F + alpha(RKstage) * delWaveAct(i,j,k)
          
          waveAct(i,j,k) = waveAct(i,j,k) + beta(RKstage) * delWaveAct(i,j,k)
          
       end do
    end do

  end subroutine add_waveFlux
    
  
  !-------------------------------------------------------------------------------------
  
  
  subroutine calc_waveFlux(var,ray, waveAct)
    
    !-----------------------------------------------
    !        Compute wave action fluxes 
    !-----------------------------------------------
    
    ! in/out fields
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var
    type(rayType), dimension(nRay), intent(in)        :: ray
    real, dimension(0:nx+1,0:ny+1,0:nz+1), intent(in) :: waveAct ! wave action

    
    
    ! local variables
    integer :: i,j,k
    real    :: fFlux, gFlux          ! local wave action fluxes
    real    :: cgx, cgz              ! group velocity components
    
    
    !----------------------------------------------------------------------
    !   calc average intrinsic group velocities at cell edges -> waveFlux
    !----------------------------------------------------------------------
    
    ! reset wave flux and ray counter
    waveFlux = 0.0
    nRayPerCell = 0
    
    !----------------------------
    !       mean cgx
    !----------------------------

    ! 1) add up all group velocities left and right of cell faces
    rayloop1: do iRay = 1, nRay
       
       if( ray(iRay)%right ) then  ! ray in right half of grid cell
          i = ray(iRay)%iCell
       else
          i = ray(iRay)%iCell-1    ! ray in left half of grid cell
       end if
       
       j = ray(iRay)%jCell
       k = ray(iRay)%kCell
       
       waveFlux(i,j,k,1) = waveFlux(i,j,k,1) + ray(iRay)%cgx
       nRayPerCell(i,j,k) = nRayPerCell(i,j,k) + 1           ! count nb of rays 

    end do rayloop1

    ! check that all (x-staggered) cells have rays
    do k = 1,nz
       j = 1
       do i = 0,nx
          if( nRayPerCell(i,j,k) == 0 ) then
             print*,"STOP. In cgx calculation: i,k,nb of rays = ", &
                  & i,k,nRayPerCell(i,j,k)
             stop
!             print*,"WARNING. Increase number of rays?"
          end if
       end do
    end do
    
    ! 2a) devide sum of cgx by nb of rays
    do k = 1,nz           ! loop without left/right boundary
       j = 1              ! change for 3D
       do i = 1,nx-1      ! only internal cell faces
          
          waveFlux(i,j,k,1) = waveFlux(i,j,k,1) / nRayPerCell(i,j,k)
       end do
    end do
    
    ! 2b) devide sum of cgx by nb of rays
    do k = 1,nz            ! loop over left/right boundary
       j = 1               ! change for 3D
       
       ! check: ok
       
       waveFlux(0,j,k,1) = (waveFlux(0,j,k,1) + waveFlux(nx,j,k,1)) / (nRayPerCell(0,j,k)+nRayPerCell(nx,j,k))
       waveFlux(nx,j,k,1) = waveFlux(0,j,k,1)
    end do
    
    ! 3) add mean flow
    do k = 1,nz
       j = 1    ! change for 3D
       do i = 0,nx
          waveFlux(i,j,k,1) = waveFlux(i,j,k,1) + var(i,j,k,2)
       end do
    end do
    

    !xxx check mean cgx: ok
!!$    do k = nz/2,nz/2
!!$       j = 1
!!$       do i = nx/2,nx/2
!!$          write(*,fmt="(a,2i5,f7.3,a)") "i, k, mean cgx = ", &
!!$               & i,k,waveFlux(i,j,k,1)*uRef, " m/s"
!!$       end do
!!$    end do
    

    !----------------------------
    !       mean cgz
    !----------------------------
    
    ! init ray counter
    nRayPerCell = 0

    ! add up all group velocities above and below cell face
    ! for a staggered cell of rays -> group velocity at cell face
    rayLoop2: do iRay = 1, nRay

       if( ray(iRay)%up ) then    ! ray in upper half of grid cell
          k = ray(iRay)%kCell
       else
          k = ray(iRay)%kCell - 1 ! ray in lower half of grid cell
       end if

       i = ray(iRay)%iCell
       j = ray(iRay)%jCell

       waveFlux(i,j,k,3) = waveFlux(i,j,k,3) + ray(iRay)%cgz
       nRayPerCell(i,j,k) = nRayPerCell(i,j,k) + 1           ! count nb of rays 

    end do rayLoop2

    ! check that all (z-staggered) cells have rays
    do k = 0,nz
       j = 1
       do i = 1,nx
          if( nRayPerCell(i,j,k) == 0) then
             print*,"STOP. In cgz calculation: i,k,nb of rays = ", &
                  & i,k,nRayPerCell(i,j,k)
             stop
!             stop"stopping. Increase number of rays."
          end if
       end do
    end do
    
    ! devide sum of cgz by nb of rays
    do k = 1,nz-1         ! only internal cell faces
       j = 1              ! change for 3D
       do i = 1,nx 
          waveFlux(i,j,k,3) = waveFlux(i,j,k,3) / nRayPerCell(i,j,k)
       end do
    end do

    ! devide sum of cgz by nb of rays
    do i = 1,nx             ! loop over lower/upper boundary
       j = 1                ! change for 3D

       ! k = 0
       waveFlux(i,j,0,3) = 0.0  ! no influx over lower boundary
       waveFlux(i,j,nz,3) = waveFlux(i,j,nz,3) / nRayPerCell(i,j,nz)
    end do
    
    ! add mean flow
    do k = 0,nz
       j = 1    ! change for 3D
       do i = 1,nx
          waveFlux(i,j,k,3) = waveFlux(i,j,k,3) + var(i,j,k,4)
       end do
    end do
    
    
    !xxx check mean cgz: ok
!!$    do k = nz/2,nz/2
!!$       j = 1
!!$       do i = nx/2,nx/2
!!$          write(*,fmt="(a,2i5,f7.3,a)") "i, k, mean cgz = ", &
!!$               & i,k,waveFlux(i,j,k,3)*uRef, " m/s"
!!$       end do
!!$    end do
    
    !-----------------------------------
    !   Calc zonal fluxes: cgx*A
    !-----------------------------------
    
    select case( waveFluxType ) 

    case( "upwind" ) 

       do k = 1, nz
          j = 1   ! change for 3D
          do i = 0,nx

             cgx = waveFlux(i,j,k,1)            ! cgx transitionally stored in waveFlux

             ! upwind decision
             if( cgx > 0.0 ) then                 
                ! take left cell value, reconstructed to the right 
                fFlux = cgx * waveActTilde(i,j,k,1,1)
!                fFlux = cgx * waveAct(i,j,k)
             else                                 
                ! take right cell value, reconstructed to the left
                fFlux = cgx * waveActTilde(i+1,j,k,1,0)
!                fFlux = cgx * waveAct(i+1,j,k)
             end if

             waveFlux(i,j,k,1) = fFlux

          end do
       end do

    case( "central" )

       do k = 1, nz
          j = 1   ! change for 3D
          do i = 0,nx

             cgx = waveFlux(i,j,k,1)            ! cgx transitionally stored in waveFlux
             waveFlux(i,j,k,1) = cgx * 0.5*( waveAct(i,j,k) + waveAct(i+1,j,k) )
             
          end do
       end do
       
    case default
       stop"calc_waveFlux: inknown waveFluxType. Stop."
    end select


    !-----------------------------------
    !   Calc meridional fluxes: cgy*A       
    !-----------------------------------

    !  (only for 3D version)


    
    !-----------------------------------
    !   Calc vertical fluxes: cgz*A
    !-----------------------------------
    
    select case( waveFluxType ) 

    case( "upwind" ) 

       do k = 0,nz
          j = 1      ! change for 3D
          do i = 1,nx

             cgz = waveFlux(i,j,k,3)        ! cgz transitionally stored in waveFlux

             ! upwind decision
             if( cgz > 0.0 ) then   ! take lower cell value, reconstructed to the top
                gFlux = cgz * waveActTilde(i,j,k,3,1)
!                gFlux = cgz * waveAct(i,j,k)
             else               ! take upper cell value, reconstructed to the bottom
                gFlux = cgz * waveActTilde(i,j,k+1,3,0)
!                gFlux = cgz * waveAct(i,j,k+1)
             end if

             waveFlux(i,j,k,3) = gFlux

          end do
       end do


    case( "central" )

       do k = 0,nz
          j = 1      ! change for 3D
          do i = 1,nx

             cgz = waveFlux(i,j,k,3)  ! cgz transitionally stored in waveFlux
             waveFlux(i,j,k,3) = cgz * 0.5*(  waveAct(i,j,k) + waveAct(i,j,k+1) )
             
          end do
       end do
       
    case default
       stop"calc_waveFlux: inknown waveFluxType. Stop."
    end select
    



  end subroutine calc_waveFlux
  
  
  !-------------------------------------------------------------------------------------
  
  
  subroutine setBoundary_wkb(waveAct, reconstCase)
    
    !-----------------------------------------------
    !  1) set periodic boundary by copying values
    !     of wave action into ghost cells
    !  2) set solid wall boundary by setting values
    !     into ghost cells (extrapolation, etc.)
    !-----------------------------------------------
    
    ! in/out variables
    real, dimension(0:nx+1,0:ny+1,0:nz+1), intent(inout) :: waveAct
    character(len=*), intent(in) :: reconstCase

    select case( reconstCase ) 


    case( 'meanValue' )

       !-----------------------------
       !   periodic boundaries
       !-----------------------------

       ! left/right
       waveAct(0,:,:) = waveAct(nx,:,:)
       waveAct(nx+1,:,:) = waveAct(1,:,:)

       ! forward/backward
       waveAct(:,0,:) = waveAct(:,ny,:)
       waveAct(:,ny+1,:) = waveAct(:,1,:)

       !-----------------------------------
       !   solid  wall / open boundaries
       !-----------------------------------

!xxx use linear extrapolation for non-reflective boundaries???

       ! lower boundary
       waveAct(:,:,0) = waveAct(:,:,1)

       ! upper boundary
       waveAct(:,:,nz+1) = waveAct(:,:,nz)


    case( 'reconstValue' ) 

       !-----------------------------
       !   periodic boundaries
       !-----------------------------
       
       ! left/right
       waveActTilde(0,:,:,1,:) = waveActTilde(nx,:,:,1,:)
       waveActTilde(nx+1,:,:,1,:) = waveActTilde(1,:,:,1,:)

       ! forward/backward
       waveActTilde(:,0,:,2,:) = waveActTilde(:,ny,:,2,:)
       waveActTilde(:,ny+1,:,2,:) = waveActTilde(:,1,:,2,:)
       
       !-----------------------------------
       !   solid  wall / open boundaries
       !-----------------------------------

       ! upper boundary
       waveActTilde(:,:,nz+1,3,0) = 0.0          ! no wave action outside domain
                                                 ! in case of donward winds

       ! lower boundary
       ! flux is set to zero in calc_waveFlux


    case default
       stop"setBoundary_wkb/wkb: unknown reconstCase"
    end select


  end subroutine setBoundary_wkb


  !-------------------------------------------------------------------------------------
  

  subroutine reconstruct_waveAction(waveAct)
    !----------------------------------------
    !  reconstruct wave action
    !  using MUSCL upwind reconstructions
    !  to avoid oscillations
    !----------------------------------------
    
    ! in/out variables
    real, dimension(0:nx+1,0:ny+1,0:nz+1), intent(inout) :: waveAct
    
    ! used module variables
    ! waveActTilde
    
    !--------------------------
    !   reconstruct along x
    !--------------------------
    
    call muscl_reconstruct3D(waveAct,nx+2,ny+2,nz+2,waveActTilde, &
         &                   limiterType, "x")
    
    !--------------------------
    !   reconstruct along y
    !--------------------------

    !call muscl_reconstruct3D(waveAct,nx+2,ny+2,nz+2,waveActTilde, &
    !     &                   limiterType, "y")
    

    !--------------------------
    !   reconstruct along z
    !--------------------------

    call muscl_reconstruct3D(waveAct,nx+2,ny+2,nz+2,waveActTilde, &
         &                   limiterType, "z")

    

    



    
  end subroutine reconstruct_waveAction
 
  
  !-------------------------------------------------------------------------------------
  
   
  subroutine calc_cellIndex(ray)
    !--------------------------------------------
    ! calculates the indices of the rays
    ! with respect to the FV-grid
    ! - periodic bc in x and z
    ! - re-init of wave number at lower boundary
    !--------------------------------------------

    ! argument list
    type(rayType), dimension(nRay), intent(inout) :: ray
    
    ! local variables


    !-----------------------------
    !       zonal indices
    !-----------------------------
    
    do iRay = 1,nRay

       xRay = ray(iRay)%x
       i_dx = (xRay -lx(0)) / dx
       iCell = ceiling( i_dx )
       
       ! rays leaving right boundary
       if (iCell >nx) then
          iCell = iCell - nx
          xRay = xRay - (lx(1)-lx(0))
       end if
       
       
       ! rays leaving left boundary
       if (iCell<1) then
          iCell = nx - iCell
          xRay = xRay + (lx(1)-lx(0))
       end if

       
       !---------------------------------------
       !  left/right position within grid cell
       !---------------------------------------
       
       i_dx = (xRay -lx(0)) / dx
       del_i = i_dx - iCell + 1
       if( del_i > 0.5 ) then
          ray(iRay)%right = .true.
       else
          ray(iRay)%right = .false.
       end if

       ! store in ray 
       ray(iRay)%iCell = iCell
       ray(iRay)%x = xRay
       
    end do
    
    
    !-----------------------------
    !       vertical indices
    !-----------------------------
    
    do iRay = 1,nRay

       zRay = ray(iRay)%z
       k_dz = (zRay-lz(0)) / dz
       kCell = ceiling( k_dz )
       
       ! rays leaving upper boundary
       if (kCell > nz) then
          kCell = kCell - nz
          ! re-insert ray in the middle of cell 1: 
          zRay = zRay - (lz(1) - lz(0)) + dz/4.0
          
          
          ray(iRay)%k = kk0               ! reset wave numbers to initial value
          ray(iRay)%m = mm0
       end if

      
       ! rays leaving lower boundary: stop
       if (kCell < 1) then 
          print*,"WARNING (wkb.f90): rays leaving lower boundary!"
          print*,"iRay = ", iRay
          print*,"zRay = ", zRay
          print*,"kCell = ", kCell
          print*,"iCell = ", ray(iRay)%iCell
          print*,"cgz   = ", ray(iRay)%cgz
!          stop"Stopping."
       end if
       
       !------------------------------------
       !  up/down position within grid cell
       !------------------------------------
       
       k_dz = (zRay-lz(0)) / dz
       del_k = k_dz - kCell + 1
       if( del_k > 0.5 ) then
          ray(iRay)%up = .true.
       else
          ray(iRay)%up = .false.
       end if
       
       ! save ray index and position
       ray(iRay)%kCell = kCell
       ray(iRay)%z = zRay
             
    end do
    
    !---------------------------------
    !      Check ray positions
    !---------------------------------
    
    if( debugging ) then
       print*,"RK Stage: ", rkStage
       write(*,fmt="(3a4,2a7)") "iRay", "i_c", "k_c", &
            & "lr_pos", "ud_pos"
       do iRay = 1,nRay
          if( ray(iRay)%right) then 
             left_right_pos = "right"
          else
             left_right_pos = "left"
          end if
          if( ray(iRay)%up) then 
             up_down_pos = "up"
          else
             up_down_pos = "down"
          end if
          write(*,fmt="(3i4,2a7)") iRay, ray(iRay)%iCell, ray(iRay)%kCell, &
               & trim(left_right_pos), trim(up_down_pos)
       end do
       print*,""
    end if

    !---------------------------------
    !        Check rays values
    !---------------------------------
    ! 1) grid cell position: ok
    ! 2) left/right position: ok
    if( debugging ) then
       print*,"RK Stage: ", rkStage
       write(*,fmt="(3a4,2a7,4a7)") "iRay", "i_c", "k_c", &
            & "lr_pos", "lr_pos", "k_wave", "m_wave", "cgx", "cgz"
       do iRay = 1,nRay
          if( ray(iRay)%right) then 
             left_right_pos = "right"
          else
             left_right_pos = "left"
          end if
          if( ray(iRay)%up) then 
             up_down_pos = "up"
          else
             up_down_pos = "down"
          end if
          write(*,fmt="(3i4,2a7,4f7.3)") iRay, ray(iRay)%iCell, ray(iRay)%kCell, &
               & trim(left_right_pos), trim(up_down_pos), &
               & ray(iRay)%k, ray(iRay)%m, ray(iRay)%cgx, ray(iRay)%cgz
       end do
       print*,"x-position of ray (4,4) = ", ray(16)%x*lRef
       print*,"z-position of ray (4,4) = ", ray(16)%z*lRef
       print*,""
    end if

  end subroutine calc_cellIndex


  !-------------------------------------------------------------------------------------


  subroutine transport_ray(var,ray,dt,rkStage)
    !------------------------------------------
    ! transports ray along group velocity
    ! trajectory
    !------------------------------------------

    ! argument list
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var
    type(rayType), dimension(nRay), intent(inout) :: ray

    real, intent(in) :: dt
    integer, intent(in) :: rkStage


    ! transport of ray
    real :: uR, uL, u    ! left/right velocities for interpolation velocity u
    real :: wT, wB, w    ! bottom/top velocities for interpolation velocity w
    real :: delta        ! linear interpolation weight
    real :: kk,mm,kTot   ! horizontal, vertical and total wave number
    real :: F            ! RK slope
    real :: cgx, cgy, cgz ! group velocity components

    ! transport of wave vector
    real :: dm
    real :: uT, uB       ! top, bottom velocity for vertical gradient of u
    real :: du_dz        ! vertical velocity gradient
    logical :: right, up ! local position of ray within cell
    
    ! group velocity
    real :: cg           ! common group velocity factor
    

    ! initialize RK-tendencies at first RK stage
    if( RKstage==1 ) then 
       dxRay = 0.0
       dkRay = 0.0
    end if
    

    ray_loop: do iRay = 1, nRay

       ! find surrounding cell
       i = ray(iRay)%iCell
       j = ray(iRay)%jCell
       k = ray(iRay)%kCell

      
       ! get position within cell
       right = ray(iRay)%right
       up = ray(iRay)%up
       
       !-----------------------------
       !     zonal displacement
       !-----------------------------

       ! interpolate background velocities to ray
       uR = var(i,j,k,2)
       uL = var(i-1,j,k,2)
       delta = (ray(iRay)%x - x(i))/dx + 0.5
       u = delta*uR + (1.-delta)*uL
       
       ! zonal group velocity component
       cgx = ray(iRay)%cgx

       ! check velocity interpolation: ok
       ! write(*,fmt="(i3,3f10.2)") iRay, ray(iRay)%x*lRef, ray(iRay)%z*lRef, u*lRef


       ! RK update
       F = cgx + u
       dxRay(1,iRay) = dt*F + alpha(rkStage)*dxRay(1,iRay)
       ! print*,"dxRay(1,iRay) = ", dxRay(1,iRay)
       ray(iRay)%x = ray(iRay)%x + beta(RKstage)*dxRay(1,iRay)


       !-----------------------------
       !     vertical displacement
       !-----------------------------

       ! interpolate background velocities to ray
       wT = var(i,j,k,4)
       wB = var(i,j,k-1,4)
       delta = (ray(iRay)%z - z(k))/dz + 0.5
       w = delta*wT + (1.-delta)*wB

       
       ! vertical group velocity component
       cgz = ray(iRay)%cgz
       
       ! RK update
       F = cgz + w
       dxRay(3,iRay) = dt*F + alpha(rkStage)*dxRay(3,iRay)
       !print*,"dxRay(3,iRay) = ", dxRay(3,iRay)
       ray(iRay)%z = ray(iRay)%z + beta(RKstage)*dxRay(3,iRay)


       !-------------------------------
       !    transport of wavenumber m
       !-------------------------------
       
       ! vertical velocity gradient u
       if( .not.right .and. up) then          ! position a)
          uT = var(i-1,j,k+1,2)
          uB = var(i-1,j,k+1,2)
       else if( right .and. up) then          ! position b)
          uT = var(i,j,k+1,2)
          uB = var(i,j,k,2)
       else if( right .and. .not. up) then         ! position d)
          uT = var(i,j,k,2)
          uB = var(i,j,k-1,2)
       else
          uT = var(i-1,j,k,2)                 ! position c) 
          uB = var(i-1,j,k-1,2)
       end if
       
       du_dz = (uT-uB)/dz
       
       
       ! set velocity gradient zero at lower boundary
       if (ray(iRay)%kCell == 1) then
          du_dz = 0.0
       end if

       if( debugging .and. iRay == 16 ) then 
          write(*,fmt = "(a10,es25.14)") "du_dz = ", du_dz
          write(*,fmt = "(a10,es25.14)") "m = ", ray(iRay)%m
       end if
          
       ! zonal wave number
       kk = ray(iRay)%k            

       ! RK procedure
       F = -kk*du_dz
       
       dkRay(3,iRay) = dt*F + alpha(rkStage)*dkRay(3,iRay)

       ray(iRay)%m = ray(iRay)%m + beta(rkStage)*dkRay(3,iRay)


       !----------------------------------------
       !     Update intrinsic group velocities
       !----------------------------------------
              
       ! get wave numbers
       mm = ray(iRay)%m             ! vertical
       kk = ray(iRay)%k             ! zonal wave number
       kTot = sqrt(mm**2 + kk**2)   ! total (in 2D)
       cg = NN*mm/kTot**3       ! common group velocity factor

       !xxx calculate cg according to Ulrich
       cgx = -cg*mm                 ! zonal group velocity
       cgz =  cg*kk                 ! vertical group velocity
       
       ray(iRay)%cgx = cgx          ! save in ray variable
       ray(iRay)%cgz = cgz

      
    end do ray_loop

    ! test ray displacement
    !    do iRay = 1,nRay
    !       write(*,fmt="(i3,2f10.2)") iRay, dxRay(1,iRay)*lRef, dxRay(2,iRay)*lRef
    !    end do
    
   

  end subroutine transport_ray


  !-------------------------------------------------------------------------------------


  subroutine setup_wkb(ray, waveAct, waveActOld, Psi)

    !------------------------------------------------
    ! allocate ray field
    ! initialize position, wave vector and frequency
    ! for rays
    !------------------------------------------------

    ! argument list
    type(rayType), dimension(:), allocatable, intent(out) :: ray
    real, dimension(:,:,:), allocatable, intent(out) :: waveAct ! wave action
    real, dimension(:,:,:), allocatable, intent(out) :: waveActOld ! wave action
    complex, dimension(:,:,:,:,:), allocatable, intent(out) :: Psi ! wave amplitudes

    ! local variables
    integer ::  allocstat
    real ::     dxRay0, dyRay0, dzRay0    ! non-dim distance between rays
    real ::     kHor, kTot                ! horizontal and total wave number
    real ::     omegaInt                  ! intrinsic frequency
    real ::     cg                        ! common group velocity coefficient

    real :: i_dx, del_i                   ! real position within grid / cell
    real :: k_dz, del_k

    ! wave packet
    real :: lambdaX, lambdaZ
    real :: delx, delz, sigma, L_cos, xCenter, zCenter, envel
    real :: bAmp, bAmp0

    

    !-------------------------
    !    Allocate fields
    !-------------------------

    if(.not.raytracer) then 

       nRayRatioX = 1 ! reduce amount of storage if raytracer is off
       nRayRatioY = 1 
       nRayRatioZ = 1 
       
       ! set some values to avoid problems during initialization 
       ! in mode raytracer = .false.

!xxxx what was that good for???
!       lambdaZ_dim = 1.0
!       lambdaX_dim = 1.0
       nxRay = 1
       nyRay = 1
       nzRay = 1
       nRay  = 1
       
    else

       nxRay = nRayRatioX * nx
       ! nyRay = nRayRatio * ny
       nyRay = 1
       nzRay = nRayRatioZ * nz
       nRay = nxRay * nyRay * nzRay

    end if
    
    ! field of rays
    allocate( ray(nRay), stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate ray"

    ! ray displacement vector
    allocate( dxRay(3,nRay), stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate dxRay"

    ! wave vector increment
    allocate( dkRay(3,nRay), stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate dkRay"
    
    ! wave action field
    allocate( waveAct(0:nx+1,0:ny+1,0:nz+1), stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate waveAct"

    ! wave action field - save old Runge-Kutta stage
    allocate( waveActOld(0:nx+1,0:ny+1,0:nz+1), stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate waveActOld"

    ! Runge-Kutta tendency of wave action
    allocate( delWaveAct(nx,ny,nz), stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate delWaveAct"
    
    ! reconstructed wave action
    allocate( waveActTilde(0:nx+1,0:ny+1,0:nz+1,1:3,0:1), stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate waveActTilde"

    ! wave action flux
    allocate( waveFlux(0:nx,0:ny,0:nz,3), stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate waveFlux"

    ! number of rays per finite volume cell 
    allocate( nRayPerCell(0:nx,0:ny,0:nz), stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate nRayPerCell"

    ! complex wave amplitudes
    allocate( Psi(0:nx+1,0:ny+1,0:nz+1,4,0:2), stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate Psi"

    ! Runge-Kutta tendency of mean Flow delU00
    allocate( delU00(0:nx+1,ny,nz), stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not allocate delU00"
    
    

    !--------------------------
    !   Init ray coordinates
    !--------------------------

    dxRay0 = dx/real(nRayRatioX)   ! initial distance between rays
    ! dyRay0 = dy/real(nRayRatioY) ! not needed unless 3D
    dzRay0 = dz/real(nRayRatioZ)


    ! x coordinates 
    do i = 1,nxRay
       xRay = lx(0) + real(i-1)*dxRay0 + dxRay0/2.0

       do k = 1,nzRay   ! all rays with this x coordinate
          iRay = i + real(k-1)*nxRay
          ray(iRay)%x = xRay
       end do
    end do

    ! y coordinate to half cell width
    do i = 1,nRay
       ray(i)%y = 0.5*dy
    end do

    ! z coordinate
    do k = 1,nzRay
       zRay = lz(0) + real(k-1)*dzRay0 + dzRay0/2.0

       do i = 1,nxRay
          iRay = i + real(k-1)*nxRay
          ray(iRay)%z = zRay
       end do
    end do


    ! check coordinates: ok
!!$    do k = 1,nzRay
!!$       do i = 1,nxRay
!!$          iRay = i + real(k-1)*nxRay
!!$          write(*,fmt="(2i3,2f10.2)") i,k, ray(iRay)%x*lRef, ray(iRay)%z*lRef
!!$       end do
!!$    end do


    !---------------------
    !     Calc ijk-position 
    !---------------------

    call calc_cellIndex(ray)
          
    ! only jCell, all the rest done by calc_cellIndex
    ray(1:nRay)%jCell = 1


    !--------------------------
    !     Calc wave vector
    !--------------------------
    lambdaX = lambdaX_dim/lRef     ! non-dim zonal wave length
    lambdaZ = lambdaZ_dim/lRef     !         vert. wave length
    
    ! zonal wave number
    if( lambdaX == 0.0) then       ! zonally uniform, lambdaX = infty
       kk0 = 0.0                   ! --> cgz = 0
    else
       kk0 = 2.0*pi/lambdaX           
    end if
    
    ! vertical wave number
    if( lambdaZ == 0.0 ) then       ! vertically uniform, lambdaZ = infty
       print*,"setup_wkb: mm0 = 0 is impossible -> cg = 0!"
    else
       mm0 = 2.0*pi/lambdaZ        !xxx m > 0 according to Ulrich for jet case
    end if
    
    ray(:)%k = kk0
    ray(:)%l = 0.0
    ray(:)%m = mm0
    
    ! non-uniform case not needed
!!$    do i = 1,nRay
!!$       ray(i)%k = 2.0*pi/lambdaX
!!$       ray(i)%l = 0.0              ! wave number zero
!!$       ray(i)%m = -2.0*pi/lambdaZ
!!$       ! check wave length
!!$       ! write(*,fmt="(i3,2f10.2)") i, 2.*pi*lRef/ray(i)%k, 2.*pi*lRef/ray(i)%m
!!$    end do


    !-------------------------
    !   Intrinsic frequency
    !-------------------------
   
    do i = 1,nRay
       kHor = ray(i)%k
       kTot = sqrt(kHor**2 + ray(i)%m**2)
!xxx
!old       omegaInt = NN*kHor / kTot
       !new        
       omegaInt = omiSign * NN*kHor / kTot   !xxx omega < 0 acc. Ulrich Achatz for jet 

       ray(i)%omega = omegaInt
       ! check frequency: 
       ! write(*,fmt="(i3,f10.5)") i, ray(i)%omega/tRef
    end do


    !---------------------------------------
    !    Intrinsic group velocity
    !---------------------------------------
    
    kTot = sqrt(mm0**2 + kk0**2)     ! total wave number
    cg = NN*mm0/kTot**3          ! common group velocity coefficient 
    
    !xxx velocity calculation according to Ulrich
    ray(:)%cgx = -cg*mm0
    ray(:)%cgy =  0.0
    ray(:)%cgz =  cg*kk0
    
    
    !-----------------------------------
    !        Initialize wave action
    !-----------------------------------
   
    ! scale envelope parameters
    xCenter = xCenter_dim/lRef     ! scaled position of wave packtet
    zCenter = zCenter_dim/lRef 
    sigma = sigma_dim/lRef         ! sigma width of Gaussian distribution
    L_cos = L_cos_dim/lRef         ! half length of cosine profile

    
    ! common amplitude coefficient
    bAmp0 = amplitudeFactor * N2/mm0
    
    do k = 1,nz
       j = 1
       do i = 1,nx
          
          ! Gaussian profile
          if( wavePacketDim == 1 ) then
             delx = 0.0
          else
             delx = (x(i)-xCenter)
          end if
          delz = (z(k)-zCenter)


          select case(wavePacketType) 

          case(1) ! Gaussian

             envel = exp(-(delx**2 + delz**2)/2./sigma**2)

          case(2)  ! Cosine

             if( abs(delz) .le. L_cos ) then
                envel = 0.5*(1.0 + cos(pi*delz/L_cos))
             else
                envel = 0.0
             end if

          case default
             stop"init.f90: unknown wavePacketType. Stop."
          end select
          
          ! buoyancy
          bAmp = envel * bAmp0

          ! wave action       !xxx new: enforce positive wave action
!          waveAct(i,j,k) = dabs(0.5*rhoStrat(k)*bAmp**2 / N2 / omegaInt)
          waveAct(i,j,k) = 0.5*rhoStrat(k)*bAmp**2 / N2 / omegaInt

       end do
    end do


    

    
    
    !-------------------------------
    !       Feedback on Screen
    !-------------------------------
    ! Feedback to user
    print*,""
    print*," 11) Ray tracer: "
    
    if( raytracer ) then
       write(*,fmt="(a25,a)") "ray tracer = ", "on"
       write(*,fmt="(a25,f7.4, a)") "horizontal wave number = ", kk0/lRef, " 1/m"
       write(*,fmt="(a25,f7.4, a)") "vertical wave number = ", mm0/lRef, " 1/m"
       write(*,fmt="(a25,f7.2, a)") "zonal mean flow = ", meanFlowX_dim, " m/s"
       
       ! group velocity
       write(*,fmt="(a25,f7.2, a)") "group vel: cg_x = ",-cg*mm0*uRef, " m/s"
       write(*,fmt="(a25,f7.2, a)") "group vel: cg_z = ", cg*kk0*uRef, " m/s"

       ! phase velocity 
       write(*,fmt="(a25,f7.2, a)") "phase vel: c_x = ",omegaInt/kk0*uRef, " m/s"
       write(*,fmt="(a25,f7.2, a)") "phase vel: c_z = ", omegaInt/mm0*uRef, " m/s"
       
       write(*,fmt="(a25,i7)") "nb of rays = ", nRay
       write(*,fmt="(a25,f7.1)") "nb rays / grid cell  = ", real(nRay)/nx/ny/nz
    else
       write(*,fmt="(a25,a)") "ray tracer = ", "off"
    end if


  end subroutine setup_wkb


  !-------------------------------------------------------------------------------


  subroutine setup_wkb_dummy(ray, waveAct,waveActOld)
    !---------------------------------------------------
    ! allocate ray field with one ray
    ! initialize position, wave vector and frequency
    ! for a single dummy ray 
    ! (to allow tec360 to work if no ray tracer is used)
    !---------------------------------------------------

    ! argument list
    type(rayType), dimension(:), allocatable, intent(out) :: ray
    real, dimension(:,:,:), allocatable, intent(out) :: waveAct
    real, dimension(:,:,:), allocatable, intent(out) :: waveActOld
    
    ! local variables
    integer ::  allocstat


    !---------------------------------------
    !    allocate fields with dummy size 1
    !---------------------------------------

    ! field of rays
    allocate( ray(1), stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb_dummy: could not allocate ray"
    
    ! field of wave action
    allocate( waveAct(1,1,1), stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb_dummy: could not allocate waveAct"
   
    ! field of wave action
    allocate( waveActOld(1,1,1), stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb_dummy: could not allocate waveActOld"
    


    !--------------------------
    !   init ray coordinates
    !--------------------------

    ray(1)%x = 0.0
    ray(1)%y = 0.0
    ray(1)%z = 0.0

    !---------------------
    !     ijk-position
    !---------------------
    
       ray(1)%iCell = 1
       ray(1)%jCell = 1
       ray(1)%kCell = 1

    !----------------------
    !     wave vector
    !----------------------
    
    ray(1)%k = 1.0
    ray(1)%l = 1.0
    ray(1)%m = 1.0

    !-------------------------
    !  intrinsic frequency
    !-------------------------

    ray(1)%omega = 1.0

  end subroutine setup_wkb_dummy

  
  !--------------------------------------------------------------------------    


  subroutine finish_wkb(ray, waveAct,waveActOld, Psi)

    !----------------------
    !   deallocate fields
    !----------------------

    ! in/out variables
    type(rayType), dimension(:), allocatable   :: ray
    real, dimension(:,:,:), allocatable        :: waveAct
    real, dimension(:,:,:), allocatable        :: waveActOld
    complex, dimension(:,:,:,:,:), allocatable :: Psi

    ! local variables
    integer :: allocstat

    ! field of rays
    deallocate( ray, stat=allocstat)
    if(allocstat /= 0) stop "finish_wkb: could not deallocate ray"

    ! ray displacement vector
    deallocate( dxRay, stat=allocstat)
    if(allocstat /= 0) stop "finish_wkb: could not deallocate dxRay"

    ! ray wave number tendencies
    deallocate( dkRay, stat=allocstat)
    if(allocstat /= 0) stop "finish_wkb: could not deallocate dkRay"
    
    ! wave amplitude field
    deallocate( waveAct, stat=allocstat)
    if(allocstat /= 0) stop "finish_wkb: could not deallocate waveAct"

    ! wave amplitude field - old Runge-Kutta stage
    deallocate( waveActOld, stat=allocstat)
    if(allocstat /= 0) stop "finish_wkb: could not deallocate waveActOld"

    ! reconstructed wave action
    deallocate( waveActTilde, stat=allocstat)
    if(allocstat /= 0) stop "finish_wkb: could not deallocate waveActTilde"

    deallocate( delWaveAct, stat=allocstat)
    if(allocstat /= 0) stop "finish_wkb: could not deallocate delWaveAct"

    ! wave action flux
    deallocate( waveFlux, stat=allocstat)
    if(allocstat /= 0) stop "finish_wkb: could not deallocate waveFlux"

    ! number of rays per cell
    deallocate( nRayPerCell, stat=allocstat)
    if(allocstat /= 0) stop "finish_wkb: could not deallocate nRayPerCell"

    ! complex wave amplitudes
    deallocate( Psi, stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not deallocate Psi"

    ! Runge-Kutta tendency of mean Flow delU00
    deallocate( delU00, stat=allocstat)
    if(allocstat /= 0) stop "setup_wkb: could not deallocate delU00"


        
  end subroutine finish_wkb


  !--------------------------------------------------------------------------    


  subroutine finish_wkb_dummy(ray, waveAct,waveActOld)

    !----------------------
    !   deallocate fields
    !----------------------

    ! in/out variables
    type(rayType), dimension(:), allocatable   :: ray
    real, dimension(:,:,:), allocatable        :: waveAct
    real, dimension(:,:,:), allocatable        :: waveActOld


    ! local variables
    integer :: allocstat

    deallocate( ray, stat=allocstat)
    if(allocstat /= 0) stop "finish_wkb: could not deallocate ray"

    ! wave amplitude field
    deallocate( waveAct, stat=allocstat)
    if(allocstat /= 0) stop "finish_wkb: could not deallocate waveAct"

    ! wave amplitude field - old Runge-Kutta stage
    deallocate( waveActOld, stat=allocstat)
    if(allocstat /= 0) stop "finish_wkb: could not deallocate waveActOld"

        
  end subroutine finish_wkb_dummy


  !--------------------------------------------------------------------------    


  function cabs(c) result(c_abs)
    
    !------------------------------------
    !  absulute value of complex number
    !------------------------------------
    
    ! in/out
    complex, intent(in) :: c
    real                :: c_abs

    ! local vars
    real :: a,b

    a = real(c)
    b = aimag(c)
    
    c_abs = sqrt(a**2 + b**2 )

  end function cabs
    

  !--------------------------------------------------------------------------    


  function cphase(c) result(phi)
    
    !------------------------------------
    !  phase of complex number
    !------------------------------------
    
    ! in/out
    complex, intent(in) :: c
    real                :: phi

    ! local vars
    real :: a,b

    a = real(c)
    b = aimag(c)
    
    ! first quadrant
    if( a>=0 .and. b>=0) then
       phi = atan(b/a)
       
       ! second quadrant
    else if( a<0. .and. b>=0. ) then
       phi = pi - atan(-b/a)

       ! third quadrant
    else if( a<0. .and. b<0.) then
       phi = -pi + atan(b/a)

       ! fourth quadrant
    else if( a>=0. .and. b<0 ) then
       phi = -atan(-b/a)
    else
       stop"wkb.f90/cphase: case not included. Stop."
    end if
       

    
  end function cphase






end module wkb_module
