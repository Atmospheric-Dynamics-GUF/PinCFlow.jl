module ice_module

  use type_module
  use atmosphere_module
  
  implicit none

  private 

  public :: setup_ice, iceTestcase_specifics
  public :: NUCn, DEPq
  public :: SIce_crit, p_saturation, SIce_threshold, awi
  public :: find_temperature, pIce
  public :: latent_heat_ice
  public :: terminal_v_nIce, terminal_v_qIce

  real, parameter :: Mole_mass_water = 18.01528e-3, Mole_mass_dryAir = 28.9644e-3
  real, parameter :: epsilon0 = Mole_mass_water/Mole_mass_dryAir
  real, parameter :: Rv = Rsp/epsilon0  ! specific gas constant for water vapor
  real, parameter :: cpv = 3.5*Rv ! c_p for water vapor
  real, parameter :: rhob = 0.81e3 ! mean snowflake density in kg / m**3
  real :: ice_crystal_volume
  integer :: ISSR_center
  integer :: first_nuc = 0
  real :: t_ramp_qv, t_relax_qv, t_start ! for qv_relaxation iceTestcase

!coefficients for linear fit of nucleation rate
  real, parameter :: afit0=-62.192670609121556
  real, parameter :: afit1=254.77490427507394  
  real, parameter :: j0 = 16.0

!correction term from Koop & Murray 2016
  real, parameter :: delta=1.522

  public :: epsilon0,first_nuc, Rv, t_ramp_qv, t_relax_qv, t_start, ISSR_center

contains

  subroutine setup_ice(var)
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    real :: SIce,T,p,m_ice,rho
    integer :: i,j,k, ISSR_width, ISSR_radius

    ! ice crystal volume and correction factor from log-normal distribution
    ice_crystal_volume = exp(9. * log(sigma_r)**2 / 2.) &
               & * 4./3*pi*radius_solution**3. 

    ! iceTestcases
    select case (iceTestcase)

    case ("homogeneous_qv") 
      var(:,:,:,nVar-3) = init_nAer * rhoRef * lRef**3
      do i=0,nx
        do j=0,ny
          do k=1,nz
             if ( fluctuationMode ) then
               rho = var(i,j,k,1)+rhoStrat(k) 
             else
               rho = var(i,j,k,1)
             end if 
            var(i,j,k,nVar-3) = init_nAer / rho * lRef**3
          end do
        end do
      end do
      var(:,:,:,nVar) = init_qv         
      var(:,:,:,nVar-2:nVar-1) = 0.0

    case ("homogeneous_SIce")        
      var(:,:,:,nVar-2:nVar-1) = 0.0
      do i=0,nx
        do j=0,ny
          do k=1,nz
             if ( fluctuationMode ) then
               rho = var(i,j,k,1)+rhoStrat(k) 
             else
               rho = var(i,j,k,1)
             end if 
            var(i,j,k,nVar-3) = init_nAer / rho * lRef**3
            call find_temperature(T,i,j,k,var)
            p = press0_dim * ( (PStrat(k)/p0)**gamma_1  +var(i,j,k,5) )**kappaInv
            var(i,j,k,nVar) = epsilon0 * init_SIce * SIce_crit(T) * pIce(T) / p 
          end do
        end do
      end do
            
    case ("qv_relaxation")
      init_qv = 0.0
      ISSR_center = int(ISSR_top*mountainHeight_dim/(dz*lRef))
      t_start = t_relax
      t_ramp_qv = t_ramp/10.0 
      t_relax_qv = t_relax/10.0
      do i=0,nx
        do j=0,ny
          do k=0,nz
            if ( fluctuationMode ) then
               rho = var(i,j,k,1)+rhoStrat(k) 
             else
               rho = var(i,j,k,1)
            end if
            var(i,j,k,nVar-3) = init_nAer / rho * lRef**3
          end do
        end do
      end do
            
    case ("stratification")
      init_SIce = SIce_crit(Temp0_dim)
      ISSR_center = ceiling(kSponge/2.)
      var(:,:,:,nVar-2:nVar-1) = 0
      p = press0_dim * ( (PStrat(ISSR_center)/p0)**gamma_1  +var(1,1,ISSR_center,5) )**kappaInv
      init_qv = epsilon0 * init_SIce * pIce(Temp0_dim) / p
      
      do i=0,nx
        do j=0,ny
          do k=0,nz
            if ( fluctuationMode ) then
               rho = var(i,j,k,1)+rhoStrat(k) 
             else
               rho = var(i,j,k,1)
            end if
            var(i,j,k,nVar-3) = init_nAer / rho * lRef**3
            var(i,j,k,nVar) = init_qv*(1+0.5*tanh(2.5*(k-ISSR_center)/ISSR_center))
          end do
        end do
      end do
     
     
    case ("ice_cloud")
      if (init_SIce > 1.3) init_SIce = 1.0
      var(:,:,:,nVar-3) = 0.0
      ISSR_center = ceiling(ISSR_top*kSponge)
      do i=0,nx
        do j=0,ny
          do k=0,nz
            if ( fluctuationMode ) then
               rho = var(i,j,k,1)+rhoStrat(k) 
             else
               rho = var(i,j,k,1)
            end if
            var(i,j,k,nVar-2) = 1.e6 !1.e-20/(rhoRef*lRef**3)*exp(-(k-ISSR_center)**2/200.0)
            var(i,j,k,nVar-1) = init_m_ice*var(i,j,k,nVar-2)
            !call find_temperature(T,i,j,k,var)
            T = Temp0_dim
            p = press0_dim * ( (PStrat(k)/p0)**gamma_1 )**kappaInv !+var(i,j,k,5) )**kappaInv
            var(i,j,k,nVar) = epsilon0 * init_SIce * pIce(T) / p
          end do
        end do
      end do
      
            
      
    case ("1D_ISSR")
      init_SIce = SIce_crit(T_nuc)
      ISSR_center = ceiling(kSponge/4.)
      ISSR_width = ceiling(5 / (dz*lRef) ) !ceiling(kSponge*0.02)
      ISSR_radius = ceiling(30 / (dz*lRef))
      if ((background=='const-N').and.(testcase=="nIce_w_test")) then
        ISSR_center = ceiling( 1./(dz*lRef) * g/N_BruntVaisala_dim**2 * &
             & log((T_nuc-kappa*g**2/(Rsp*N_BruntVaisala_dim**2))/(theta0_dim-kappa*g**2/(Rsp*N_BruntVaisala_dim**2)))) 
        print*, "!!!!!!!!!!! ISSR_center = ", ISSR_center
      end if
      if ((background=='isentropic').and.(testcase=="nIce_w_test")) then
        ISSR_center = ceiling( 1./(dz*lRef) * Rsp * theta0_dim / (kappa*g) * (1.-T_nuc/theta0_dim))
        print*, "!!!!!!!!!!! ISSR_center = ", ISSR_center
      end if
              
      var(:,:,:,nVar-2:nVar-1) = 0.
      p = press0_dim * ( (PStrat(ISSR_center)/p0)**gamma_1  +var(1,1,ISSR_center,5) )**kappaInv
      init_qv = epsilon0 * init_SIce * pIce(T_nuc) / p
      ISSR_center = ISSR_center - (ISSR_radius+2*ISSR_width) !ceiling(backgroundFlow(3) / dz)
      do i=0,nx
        do j=0,ny
        
          do k=1,ISSR_center-ISSR_radius
             if ( fluctuationMode ) then
               rho = var(i,j,k,1)+rhoStrat(k) 
             else
               rho = var(i,j,k,1)
             end if
            var(i,j,k,nVar-3) = init_nAer / rho * lRef**3
            ! Gaussian profile:
            var(i,j,k,nVar) = init_qv * exp( -1.*(k-ISSR_center+1)**2. / (2. * ISSR_width**2))
          end do
          
          do k=ISSR_center-ISSR_radius+1,ISSR_center+ISSR_radius-1
             if ( fluctuationMode ) then
               rho = var(i,j,k,1)+rhoStrat(k) 
             else
               rho = var(i,j,k,1)
             end if 
            var(i,j,k,nVar-3) = init_nAer / rho * lRef**3 
            var(i,j,k,nVar) = init_qv
          end do
          
          do k=ISSR_center+ISSR_radius,nz
             if ( fluctuationMode ) then
               rho = var(i,j,k,1)+rhoStrat(k) 
             else
               rho = var(i,j,k,1)
             end if 
            var(i,j,k,nVar-3) = init_nAer / rho * lRef**3 
            ! Gaussian profile:
            var(i,j,k,nVar) = init_qv * exp( -1.*(k-ISSR_center-1)**2. / (2. * ISSR_width**2))
          end do
          
        end do
      end do

    case default
      stop "iceTestcase: unknown case name"

    end select

  end subroutine setup_ice

!------------------------------------------------------------------------------------

  logical function iceTestcase_specifics(time, var)
  ! in/out variables
    real, intent(in) :: time
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(inout) :: var
    
  ! local integers
    integer :: i,k,i0
    
  ! additional variables
    real :: T,p,frac
  
    i0=is+nbx-1
  
    select case (iceTestcase)
          
    ! iceTestcase: qvrelaxation
      case("qv_relaxation")
      
        if (time > t_start + t_relax_qv) then
          iceTestcase_specifics = .true.
        else 
            if (time > t_start + t_relax_qv - t_ramp_qv) then
                do k=0,nz
                          if (k<ISSR_center-8) then
                            frac = 0.4
                          else if (k>ISSR_center) then
                            frac = 0.05
                          else 
                            frac = 0.98
                          end if
                          do i=0,nx
                            if (x(i0+i)<lx(0)+0.1*(lx(1)-lx(0))) then  
                              p = press0_dim * ( (PStrat(k)/p0)**gamma_1  +var(i,1,k,5) )**kappaInv
                              call find_temperature(T,i,1,k,var)
                              init_SIce = SIce_crit(T)
                              init_qv = epsilon0 * frac*init_SIce * pIce(T) / p
                              var(i,1,k,nVar) = init_qv * cos(pi*(time-t_start)/(2*t_ramp_qv))
                              !var(i,1,k,nVar-3) = init_nAer * cos(pi*(time-t_start)/(2*t_ramp_qv))
                            end if
                            if (x(i0+i)>lx(1)-0.1*(lx(1)-lx(0))) then 
                              p = press0_dim * ( (PStrat(k)/p0)**gamma_1  +var(i,1,k,5) )**kappaInv
                              call find_temperature(T,i,1,k,var)
                              init_SIce = SIce_crit(T)
                              init_qv = epsilon0 * frac*init_SIce * pIce(T) / p
                              var(i,1,k,nVar) = init_qv * cos(pi*(time-t_start)/(2*t_ramp_qv))
                              !var(i,1,k,nVar-3) = init_nAer * cos(pi*(time-t_start)/(2*t_ramp_qv))
                            end if
                          end do
                        end do
                iceTestcase_specifics = .true.
            else
                if (time > t_start + t_ramp_qv) then
                    do k=0,nz
                          if (k<ISSR_center-8) then
                            frac = 0.4
                          else if (k>ISSR_center) then
                            frac = 0.05
                          else 
                            frac = 0.98
                          end if
                          do i=0,nx
                            if (x(i0+i)<lx(0)+0.1*(lx(1)-lx(0))) then 
                              p = press0_dim * ( (PStrat(k)/p0)**gamma_1  +var(i,1,k,5) )**kappaInv
                              call find_temperature(T,i,1,k,var)
                              init_SIce = SIce_crit(T)
                              init_qv = epsilon0 * frac*init_SIce * pIce(T) / p
                              var(i,1,k,nVar) = init_qv
                              !var(i,1,k,nVar-3) = init_nAer
                            end if
                            if (x(i0+i)>lx(1)-0.1*(lx(1)-lx(0))) then 
                              p = press0_dim * ( (PStrat(k)/p0)**gamma_1  +var(i,1,k,5) )**kappaInv
                              call find_temperature(T,i,1,k,var)
                              init_SIce = SIce_crit(T)
                              init_qv = epsilon0 * frac*init_SIce * pIce(T) / p
                              var(i,1,k,nVar) = init_qv
                              !var(i,1,k,nVar-3) = init_nAer
                            end if
                          end do
                        end do
                    iceTestcase_specifics = .true.
                else 
                    if (time > t_start) then
                        if (init_qv==0.0) then
                            print*,"#### water vapor injection start ####"
                        end if
                        do k=0,nz
                          if (k<ISSR_center-8) then
                            frac = 0.4
                          else if (k>ISSR_center) then
                            frac = 0.05
                          else 
                            frac = 0.98
                          end if
                          do i=0,nx
                            if (x(i0+i)<lx(0)+0.1*(lx(1)-lx(0))) then 
                              p = press0_dim * ( (PStrat(k)/p0)**gamma_1  +var(i,1,k,5) )**kappaInv
                              call find_temperature(T,i,1,k,var)
                              init_SIce = SIce_crit(T)
                              init_qv = epsilon0 * frac*init_SIce * pIce(T) / p
                              var(i,1,k,nVar) = init_qv * sin(pi*(time-t_start)/(2*t_ramp_qv))
                              !var(i,1,k,nVar-3) = init_nAer * sin(pi*(time-t_start)/(2*t_ramp_qv))
                            end if
                            if (x(i0+i)>lx(1)-0.1*(lx(1)-lx(0))) then 
                              p = press0_dim * ( (PStrat(k)/p0)**gamma_1  +var(i,1,k,5) )**kappaInv
                              call find_temperature(T,i,1,k,var)
                              init_SIce = SIce_crit(T)
                              init_qv = epsilon0 * frac*init_SIce * pIce(T) / p
                              var(i,1,k,nVar) = init_qv * sin(pi*(time-t_start)/(2*t_ramp_qv))
                              !var(i,1,k,nVar-3) = init_nAer * sin(pi*(time-t_start)/(2*t_ramp_qv))
                            end if
                          end do
                        end do
                        iceTestcase_specifics = .true.
                    else 
                        iceTestcase_specifics = .false.
                    end if
                end if
            end if
        end if
        
      case default
        iceTestcase_specifics = .true.
        
    end select
  
  end function iceTestcase_specifics
  
!------------------------------------------------------------------------------------
! returns terminal terminal_velocity_correction  
  real function terminal_velocity_correction(T,p)
    ! in/out variables
    real, intent(in) :: T, p
    
    terminal_velocity_correction = (p/ 30000.)**(-0.178) * (T/233.)**(-0.394) 
    
  end function terminal_velocity_correction

!------------------------------------------------------------------------------------
! returns terminal sedimentation velocity of nIce in m/s
  real function terminal_v_nIce(m_ice, T, p)
    ! in/out variables
    real, intent(in) :: m_ice, T, p

    ! parameters from fit
    real, parameter :: r0 = 3.
    real, parameter :: a = 63292.37, b = 0.5727273
    real, parameter :: c = 1.1
    real, parameter :: m0 = 2.35e-8
    real, parameter :: an = a * r0**( 0.5 * b * (b-1.) )
    real, parameter :: ex1 = b*c
    real, parameter :: ex2 = 1./c
    real :: corr
    
    corr = terminal_velocity_correction(T,p) 
      
    if (sedimentation_on) then
      terminal_v_nIce = corr * an * m_ice**b * ( m0**ex1 / (m_ice**ex1 + m0**ex1) )**ex2
    else 
      terminal_v_nIce = 0.0
    end if
    
  end function terminal_v_nIce

!----------------------------------------------

! returns terminal sedimentation velocity of qIce in m/s
  real function terminal_v_qIce(m_ice, T, p)
    ! in/out variables
    real, intent(in) :: m_ice, T, p

    ! parameters from fit
    real, parameter :: r0 = 3.
    real, parameter :: a = 63292.37, b = 0.5727273
    real, parameter :: c = 1.1
    real, parameter :: m0 = 8.e-9
    real, parameter :: aq= a * r0**( 0.5 * b * (b+1.) )
    real, parameter :: ex1 = b*c
    real, parameter :: ex2 = 1./c
    real :: corr
      
    corr = terminal_velocity_correction(T,p) 
        
    if (sedimentation_on) then
      terminal_v_qIce = corr * aq * m_ice**b * ( m0**ex1 / (m_ice**ex1 + m0**ex1) )**ex2
    else 
      terminal_v_qIce = 0.0
    end if
  
  end function terminal_v_qIce

!----------------------------------------------

  real function latent_heat_ice(T)
    ! in/out variables
    real, intent(in) :: T
    
    latent_heat_ice = ( 46782.5+35.8925*T-0.07414*T**2 &
                & +541.5*exp(-1.*(T/123.75)**2) ) / Mole_mass_water 

  end function latent_heat_ice

!----------------------------------------------

  real function approximation_model(SIce,T)
    ! in/out variables
    real, intent(in) :: T, SIce 

!    coefficients for Koop-Polynom P3
!    correction +6. is due to the conversion to SI units
      real, parameter :: pk0=-906.7+6.
      real, parameter :: pk1=8502.0
      real, parameter :: pk2=-26924.0
      real, parameter :: pk3=29180.0
    
      real :: delta_aw, awi0

   awi0 = awi(T)
 
   select case (NUC_approx_type)

     case ( "linFit" )
       delta_aw = (SIce - 1.0) * awi0
       approximation_model = afit0 - delta + afit1 * delta_aw

     case ( "Koop" )
       delta_aw = (SIce - 1.0) * awi0
       approximation_model = pk0 + pk1*delta_aw + pk2*delta_aw**2 &
        &  + pk3*delta_aw**3 - delta

     case ( "threshold" )
       approximation_model = afit1 * &
                           & (SIce - SIce_threshold(T,awi0)) * &
                           & awi0 + j0

     case default
       stop "NUC_approx_type: unknown case model."

   end select
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (approximation_model .ge. 100.0) then
  
    print*,"WARNING: approximation_model = ",approximation_model,", the value has been adjusted to 100"
    approximation_model = 100.0
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end function approximation_model

!----------------------------------------------

  real function NUCn(i,j,k,var,SIce,T,p,m_ice)
    ! in/out variables
    real, intent(in) :: SIce,T,p,m_ice
    integer, intent(in) :: i,j,k
    integer :: iVar
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var
    real :: rho
 
    if (nucleation_on) then
     if (SIce .ge. SIce_crit(T)) then 
      if (first_nuc == 0) first_nuc = k
      NUCn = tRef * var(i,j,k,nVar-3) * ice_crystal_volume * 10.0**approximation_model(SIce,T)
      if (.true.) then !((first_nuc == k) .and. (first_nuc .ne. 0)) then
        print*," "
        print*, "#########################"
        print*, "Nucleation event! k = ", k,", i = ",i+is+nbx-1
        print*, "nIce = ", var(i,j,k,nVar-2) / (rhoRef*lRef**3)
        print*, "qice = ", var(i,j,k,nVar-1)
        if ( fluctuationMode ) then
               rho = var(i,j,k,1)+rhoStrat(k) 
            else
               rho = var(i,j,k,1)
        end if 
        print*, "rho = ", rho * rhoRef
        print*, "nAer = ", var(i,j,k,nVar-3) *rho / lRef**3
        print*, "T = ", T
        print*, "p = ", p
        print*, "qv = ",var(i,j,k,nVar)
        print*, "SIce = ", SIce
        print*, "SIce_crit = ",SIce_crit(T)
        print*, "NUC = ", NUCn / (rhoRef*lRef**3 * tRef)
      end if
     else 
      NUCn = 0.0
     end if
    else
      NUCn = 0.0
    end if
    
  end function NUCn

!----------------------------------------------

  real function DEPq(i,j,k,var,SIce,T,p,m_ice)
    ! in/out variables
    real, intent(in) :: SIce,T,p,m_ice
    integer, intent(in) :: i,j,k
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var

   ! parameter for calculation
   real :: lambda, kT, dv, cm, mu, dvstar, how, Schmidt_term
   real :: corr, lat_heat, fkin, a, b, r, schmidt23, gv, bd, bk
   real, parameter :: cunn = 0.7, alpham = 0.5, alphat = 1.
   real, parameter :: c1 = 0.8198373822 ! = 1.42*r0**(-0.5)
   real, parameter :: c2 = 1.5
   real, parameter :: fac_1 = 0.014038499
   real, parameter :: fac_2 = 0.293369825
   ! where fac_i = a_i * r0 ** ( 0.5 * b_i * (b_i - 1.) )
   ! here: a_1 = 0.015755, a_2 = 0.33565, r0 = 3.0 and
   real, parameter :: b_1 = 0.3, b_2 = 0.43 
   ! TODO: where do these numbers come from?
   real, parameter :: av=1.1e6
   real, parameter :: bv=0.51
   real, parameter :: cv=0.57
   real, parameter :: m0=4.4e-9
   real, parameter :: gv0=0.148564433505542
   real :: rho
   
    !if ( fluctuationMode ) then
    !     rho = var(i,j,k,1)+rhoStrat(k) 
    !else
    !     rho = var(i,j,k,1)
    !end if 
    ! rho = rho*rhoRef

    rho = p / (Rsp*T)

    ! #### find lambda #### !
    lambda = 6.6e-8 * T/293.15 * 101325./p

    ! #### find kT #### !
    if (kT_linFit) then
    kT = 0.00122990325719493 + 8.43749062552794e-05*T   
    else    
    kT = 0.002646*T**1.5 / (T + 245.*10**(-12./T) )
    end  if

    ! #### find dv #### !
    if (dv_exp2) then
    dv = 2.142e-05*(T/273.15)**2. * 101325./p
    else
    dv = 2.11e-5 * (T/273.15)**1.94 * 101325./p
    end if

    ! #### find cm #### !
    if (cm_dryAir) then
        cm = sqrt(8.*RSp*T/pi)
    else
        cm = sqrt(8.*Rv*T/pi)
    end if
    
    ! #### find mu #### !
    if (mu_linFit) then
        mu = 2.14079e-6 + 5.57139e-8 * T
    else
        mu = ( 1.458e-6 * T**1.5 ) / ( T + 110.4 )
    end if

    ! #### find dvstar #### !
    a = lambda * cunn
    b = 4. * dv / (alpham * cm)
    r = (3.* c1 * m_ice / (4.*pi*rhob) )**(1./3.)
    fkin = ( r**2 + a*r ) / ( r**2 + b*r + a*b )
    dvstar = dv * fkin
    
    ! #### find correction for small crystals #### !
    r = (3.* m_ice / (4.*pi*rhob) )**(1./3.)
    bd = 4.*dv / (alpham*sqrt(8.*Rv*T/pi))
    bk = 4.*kt / (alphat*sqrt(8.*Rsp*T/pi)*rho*3.5*Rsp)
    corr = (r*r+bk*r+a*bk)/(r*r+bd*r+a*bd)
    
    ! #### find how #### !
    lat_heat = latent_heat_ice(T) / (Rv * T)
    how = 1. / ( (lat_heat-1.) * lat_heat * corr * dv / kT + Rv * T / pIce(T) )
    
    ! #### find Schmidt_term #### !   
    corr = terminal_velocity_correction(T,p) 
    schmidt23 = (mu / ( rho * dv )) **(2./3.)
    gv = gv0 * schmidt23 * corr * rho / mu
    Schmidt_term = 1. + gv * av * (m_ice*c2)**(bv+cv) * m0**cv / &
                    & ( m_ice**cv + m0**cv )
    
    DEPq = var(i,j,k,nVar-2) * 4*pi * how * &
            & (fac_1 * m_ice**b_1 + fac_2 * m_ice**b_2) * &
            & dvstar * Schmidt_term * (SIce-1.0)  &
            & / ( rhoRef * lRef**3 )* tRef
            
        
    if ((DEPq .lt. 0.0) .and. (var(i,j,k,nVar).le.1.e-10)) then
      DEPq = 0.0
    end if
    
  end function DEPq

!----------------------------------------------

  ! from Murphy and Koop, 2005
  real function pIce(T)  ! ice pressure
    ! in/out variables
    real,intent(in) :: T
      
    pIce=exp(9.550426-5723.265/T+3.53068*log(T)-0.00728332*T)

  end function pIce

!----------------------------------------------

  ! from Murphy and Koop, 2005
  real function p_saturation(T)  ! water vapor saturation pressure
    ! in/out variables
    real,intent(in) :: T

    p_saturation=exp(54.842763-6763.22/T-4.210*log(T)+0.000367*T &
        &   +tanh(0.0415*(T-218.8))*(53.878-1331.22/T-9.44523*log(T) &
       &    +0.014025*T))

  end function p_saturation

!----------------------------------------------

  real function SIce_crit(T)
    ! in/out variables
    real, intent(in) :: T
    real :: x0, awi0
    ! parameter
  !  real :: s0,s1

  ! s1 = -3.4057422903191969E-3 
  ! s0 = 2.2525409204521596
  ! SIce_crit = s0 + s1*T

    x0 = 0
    awi0 = awi(T)
    SIce_crit = (x0-j0) / (afit1 * awi0) + SIce_threshold(T,awi0)
  end function SIce_crit

!----------------------------------------------

  real function SIce_threshold(T, awi0)
   ! in/out variables
    real, intent(in) :: T, awi0

    !coefficients for treshold fit
      real, parameter :: s10=2.27697
      real, parameter :: s11=-0.00347231
      real, parameter :: s20=1.67469
      real, parameter :: s21=0.00228125
      real, parameter :: s22=-1.36989e-05
  
      select case ( SIce_threshold_type )
        case ( "linFit" )
          SIce_threshold = s10 + s11 * T

        case ( "quadFit" )
          SIce_threshold = s20 + s21 * T + s22 * T**2

        case ( "exact" )
          SIce_threshold = ( j0 - afit0 + delta ) / ( afit1 * awi0 ) +1.
  
        case default
      stop "SIce_threshold_type: unknown case model."

  end select

  end function SIce_threshold

!----------------------------------------------

  real function awi(T)
    ! in/out variables
    real, intent(in) :: T
      
    ! fit coefficients
    real, parameter :: awi00=0.574312
    real, parameter :: awi10=-0.197855
    real, parameter :: awi11=0.00367699
    real, parameter :: awi20=1.4962
    real, parameter :: awi21=-0.0125061
    real, parameter :: awi22=3.85311e-05

    select case ( awi_type )
   
      case ( "const" )
        awi = awi00

      case ( "linFit" )
        awi = awi10 + awi11 * T

      case ( "quadFit" )
        awi =awi20 + awi21 * T + awi22 * T**2 

      case ( "exact" )
        awi = pIce(T) / p_saturation (T)

      case default
        stop "awi_type: unknown case model."

    end select


  end function awi

!----------------------------------------------

  ! calculate the current absolute temperature in Kelvin
  subroutine find_temperature(T,i,j,k,var)
    ! in/out variables
    real, intent(inout) :: T
    integer, intent(in) :: i,j,k
    real, dimension(-nbx:nx+nbx,-nby:ny+nby,-nbz:nz+nbz,nVar), &
         & intent(in) :: var
    real :: rho

    if (fluctuationMode ) then
      rho = var(i,j,k,1) + rhoStrat(k)
    else
      rho = var(i,j,k,1)
    end if

    select case ( model )
      case("pseudo_incompressible")
        T = ( Pstrat(k) / rho ) * thetaRef &
          & * ( (PStrat(k)/p0)**gamma_1  +var(i,j,k,5) )

      case("Boussinesq")
        T = ( thetaStrat(k) + var(i,j,k,6) ) * thetaRef &
           & * ( (PStrat(k)/p0)**gamma_1  +var(i,j,k,5) )
  
      case default
         stop "find_temperature: undefined model."

    end select

end subroutine find_temperature


end module ice_module
