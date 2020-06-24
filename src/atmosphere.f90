module atmosphere_module

  use type_module

  implicit none

  public   ! all objects are known to programmes using this module

  !--------------
  !    public
  !--------------
  public :: init_atmosphere
  public :: terminate_atmosphere

  real, dimension(:), allocatable :: PStrat, rhoStrat, thetaStrat, &
                                     bvsStrat
  !UAB
  real, dimension(:), allocatable :: PStrat_0, rhoStrat_0
  real, dimension(:), allocatable :: pistrat
  !UAE
  real, dimension(:), allocatable :: PStratTilde, rhoStratTilde, &
                                     thetaStratTilde

  real, dimension(:), allocatable :: Ro, RoInv

  real, dimension(:), allocatable :: PStrat00, PStrat01, rhoStrat00, rhoStrat01, thetaStrat00, thetaStrat01, bvsStrat00, &
        bvsStrat01, PStratTilde00, PStratTilde01, rhoStratTilde00, rhoStratTilde01, thetaStratTilde00, thetaStratTilde01 



  ! reference quantites
  real :: rhoRef     ! reference density
  real :: pRef       ! reference pressure
  real :: uRef, aRef ! reference velocity / speed of sound
  real :: lRef       ! reference length / pressure scale height h_sc
  real :: tRef       ! ref time (time acoustic signal needs to pass domain)
  real :: thetaRef   ! ref (potential) temperature   (Tref  = thetaRef)
  real :: Fref       ! reference force


  ! natural constants
  real :: gamma                     ! adiabatic index 
  real :: gammaInv                  ! 1/gamma
  real :: gamma_1                   ! gamma-1
  real :: kappa                     ! (gamma-1)/gamma
  real :: kappaInv                  ! 1/kappa

  real, parameter  :: g = 9.81      ! gravitational constant in m/s^2
  real, parameter  :: Rsp = 287.0   ! spec. gas const. for dry air in J/kg/K
  real :: g_ndim                    ! nondimensional gravitational constant

  ! flow parameters
  real :: Re                     ! Reynolds number (calc. from input 1/Re)
  real :: Ma, MaInv2,Ma2         ! Mach number and 1/Ma^2, Ma^2
  real :: Fr, FrInv2,Fr2         ! Froude number Fr and 1/Fr^2, Fr^2
  real :: sig                    ! Ma^2/Fr^2
  ! real :: Ro, RoInv         ! Rossby number and its inverse

  ! pressure scale height
  real :: hp

  ! pot. temperature scale height
  real :: hTheta

  ! isentropic atmosphere
  real :: theta0
  real :: term

  ! stable atmosphere
  real :: N2                     ! scaled square of Brunt-Vaisala frequency
  real :: NN                     ! scaled of Brunt-Vaisala frequency
  real :: coeff                  ! long coefficient

  !UAB
  ! Held-Suarez atmosphere
  real :: tp_strato               ! stratosphere temperature
  real :: tp_srf_trp              ! tropical surface temperature
  real :: tpdiffhor_tropo         ! tropospheric temperature difference 
                                  ! between poles and tropics
  real :: ptdiffvert_tropo        ! vertical potential-temperature 
                                  ! difference in troposphere
  !UAE

  ! isothermal
  real :: T0                     ! scaled background temperature

  ! general 
  real :: p0                     ! scaled reference pressure at z=0

  real :: zk                     ! zk = z(k) 

! achatzb
  real :: mountainHeight,mountainWidth,k_mountainw
  real :: x_center,y_center
! achatze



contains


  subroutine init_atmosphere
    !----------------------------------------------------------------------
    ! allocate and initialize: coordinate system and background atmosphere
    ! set reference and thermodynamic quantities
    !----------------------------------------------------------------------

    ! local variables
    integer :: allocstat
    integer :: i,j,k
    real :: pBar, thetaBar
    real :: p1, p2, power   ! exponents
    real :: gammaInv
    real :: eps
    real :: zmax
    real :: zk_half

    ! realistic atmosphere: isentropic troposphere / isothermal stratosphere
    real :: z_tr                   ! scaled height of troposphere
    real :: theta_tr               ! scaled const. pot. temp. of troposphere
    real :: press_tr               ! pressure at tropopause
    real :: P_tr                   ! pressure variable at tropopause
    real :: T_tr                   ! temperature at tropopause
    real :: delZ                   ! distance to tropopause

!gagarinab
    ! for baroclinic case
    real :: T_c_b1, pow_t, pow_s, p_t_b  ! tropopause quantities
    real :: T_bar, T_c_b, p_bar          ! calculated quantities
    real :: tp_sponge !FS
!gagarinae

    !UAB
    real :: pistar, thetastar
    !UAC

    ! debugging
    integer,parameter :: errorlevel = 10 ! 0 -> no output

    integer :: j00
    real :: yloc, ymax
    real, dimension(0:ny+1) :: f_Coriolis_y

    ! allocate PStrat
    if( .not. allocated(pStrat) ) then
       allocate( Pstrat(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate pStrat"
    end if

     !UAB 200413
    ! allocate pStrat_0
    if( .not. allocated(pStrat_0) ) then
       allocate( pStrat_0(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) then
          stop "atmosphere.f90: could not allocate pStrat_0"
       end if
    end if
    !UAE 200413

        ! allocate PStrat
    if( .not. allocated(pStrat00) ) then
       allocate( Pstrat00(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate pStrat"
    end if

    ! allocate PStrat
    if( .not. allocated(pStrat01) ) then
       allocate( Pstrat01(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate pStrat"
    end if

    !UAB
    ! allocate pistrat
    if( .not. allocated(pistrat) ) then
       allocate( pistrat(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate pistrat"
    end if
    !UAE

    ! allocate pStratTilde -> P at half levels
    if( .not. allocated(pStratTilde) ) then
       allocate( PstratTilde(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not all. pStratTilde"
    end if

        ! allocate pStratTilde -> P at half levels
    if( .not. allocated(pStratTilde00) ) then
       allocate( PstratTilde00(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not all. pStratTilde"
    end if

    ! allocate pStratTilde -> P at half levels
    if( .not. allocated(pStratTilde01) ) then
       allocate( PstratTilde01(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not all. pStratTilde"
    end if


    ! allocate rhoStrat
    if( .not. allocated(rhoStrat) ) then
       allocate( rhoStrat(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate rhoStrat"
    end if

    !UAB 200413
    ! allocate rhoStrat_0
    if( .not. allocated(rhoStrat_0) ) then
       allocate( rhoStrat_0(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) then
          stop "atmosphere.f90: could not allocate rhoStrat_0"
       end if
    end if
    !UAE 200413

    ! allocate rhoStrat
    if( .not. allocated(rhoStrat00) ) then
       allocate( rhoStrat00(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate rhoStrat"
    end if

    ! allocate rhoStratTilde 
    if( .not. allocated(rhoStratTilde00) ) then
       allocate( rhoStratTilde00(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could n. all. rhoStratTilde"
    end if



    ! allocate rhoStrat
    if( .not. allocated(rhoStrat01) ) then
       allocate( rhoStrat01(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate rhoStrat"
    end if

    ! allocate rhoStratTilde 
    if( .not. allocated(rhoStratTilde01) ) then
       allocate( rhoStratTilde01(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could n. all. rhoStratTilde"
    end if

    ! allocate rhoStratTilde 
    if( .not. allocated(rhoStratTilde) ) then
       allocate( rhoStratTilde(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could n. all. rhoStratTilde"
    end if


    ! allocate thetaStrat
    if( .not. allocated(thetaStrat) ) then
       allocate( thetaStrat(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not all. thetaStrat"
    end if

    if( .not. allocated(thetaStrat00) ) then
       allocate( thetaStrat00(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not all. thetaStrat"
    end if

    ! allocate squared Brunt-Vaisala frequency bvsStrat
    if( .not. allocated(bvsStrat) ) then
       allocate( bvsStrat(-1:nz+1),stat=allocstat)
       if(allocstat /= 0) &
         & stop "atmosphere.f90: could not allocate bvsStrat"
    end if

    ! allocate thetaStratTilde
    if( .not. allocated(thetaStratTilde) ) then
       allocate( thetaStratTilde(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: c. n. all. thetaStratTilde"
    end if

    if( .not. allocated(thetaStratTilde00) ) then
       allocate( thetaStratTilde00(-1:nz+2),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: c. n. all. thetaStratTilde"
    end if

    ! allocate Ro
    if( .not. allocated(Ro) ) then
       allocate( Ro(0:ny+1),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate Ro"
    end if

    ! allocate RoInv
    if( .not. allocated(RoInv) ) then
       allocate( RoInv(0:ny+1),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate RoInv"
    end if


    !----------------------------------
    !       auxiliary quantities
    !----------------------------------

    if ( testCase == "agnesiMountain" ) then
       gamma = 1.0696864111498257            ! cf. Klein 2008
    else
       gamma = 1.4
    end if
    gamma_1 = gamma-1.0
    kappa = (gamma-1.0)/gamma                ! = R/c_p = Gamma at Klein 2008
    kappaInv = 1.0/kappa             
    gammaInv = 1.0/gamma


    !----------------------------------
    !      reference quantities & 
    !      nondimensional numbers
    !----------------------------------

    select case (referenceQuantities)

    case( "general" )
       ! free references 
       rhoRef = 1.184             ! in kg/m^3
       pRef = 101325.0            ! in Pa = kg/m/s^2
       lRef = pRef/rhoRef/g       ! in m
       Ma = 0.1

       ! dependent references 
       aRef = sqrt(pRef/rhoRef)   ! in m/s
       uRef = Ma*aRef             ! in m/s
       tRef = lRef/uRef           ! in s
       thetaRef = aRef**2/Rsp     ! in K
       Fr = uRef / sqrt(g*lRef)   ! 

    case( "Klein" )
       rhoRef = 1.184             ! in kg/m^3
       pRef = 101325.0            ! in Pa = kg/m/s^2
       aRef = sqrt(pRef/rhoRef)   ! in m/s
       uRef = aRef                ! - "" -
       lRef = pRef/rhoRef/g       ! in m
       tRef = lRef/aRef           ! in s
       thetaRef = aRef**2/Rsp     ! in K
       FRef = rhoRef*uRef**2/lRef ! in N/m^3 = reference force

       ! nondim numbers
       Ma = uRef / aRef           ! Ma = 1 
       Fr = uRef / sqrt(g*lRef)   ! Fr = 1
       !                          ! sig = Ma^2/Fr^2 = 1

    case( "WKB" )
       eps = 0.1         ! asymptotic parameter ( Ma = Fr for eps = kappa! )
       rhoRef = 1.184             ! in kg/m^3
       pRef = 101325.0            ! in Pa = kg/m/s^2

       aRef = sqrt(pRef/rhoRef)   ! in m/s
       thetaRef = aRef**2/Rsp     ! in K
       Ma = eps/sqrt(kappa)
       Fr = sqrt(eps)
       uRef = Ma*aRef             ! in m/2
       lRef = eps*pRef / (kappa*rhoRef*g)   ! in m
       tRef = lRef / uRef

    case( "SI" )
       stop "init_atmosphere: Problems w. Exner pr.. Use Klein's scaling."
       ! with this scaling all quantities are 
       ! in SI units
       ! Note that in this case the thermodynamic
       ! equations have a different form
       ! while the Euler equations are unchanged
       rhoRef = 1.0               ! in kg/m^3
       pRef = 1.0                 ! in Pa = kg/m/s^2
       lRef = 1.0                 ! in m
       tRef = 1.0                 ! in s
       uRef = 1.0                 ! in m/s
       aRef = 1.0                 ! in m/s
       thetaRef = 1.0             ! in K

       ! nondim numbers
       Ma = 1.0                   ! normally 1/sqrt(kappa) 
       !                          ! but we use Exner-pr./kappa like R. Klein
       Fr = 1.0/sqrt(g)           ! Fr = uRef/sqrt(g*lRef)

    case default
       print*,"referenceQuantities = ", referenceQuantities
       stop "init_atmosphere: unknown referenceQuantities. Stopping."
    end select

    ! auxiliary nondimensionals
    Fr2 = Fr**2
    Ma2 = Ma**2
    MaInv2 = 1.0/Ma**2
    FrInv2 = 1.0/Fr**2
    sig = Ma**2/Fr**2   

    ! Rossby number
!   achatzb correction for zero Coriolis
!   Ro = uRef/f_Coriolis_dim/lRef  
!   RoInv = 1.0/Ro   
    j00=js+nby-1
    if (TestCase == "baroclinic_LC") then
       ymax = ly_dim(1)/lRef  
       do j = 0,ny+1
          yloc = y(j+j00)
          f_Coriolis_y(j) = f_Coriolis_dim*sin(pi*yloc/ymax)
          if(f_Coriolis_y(j) /= 0.0) then
             Ro(j) = uRef/f_Coriolis_y(j)/lRef  
             RoInv(j) = 1.0/Ro(j)       
          else
             Ro(j) = 1.d40
             RoInv(j)=0.0
          end if
       end do
    else
       if(f_Coriolis_dim /= 0.0) then
          Ro(:) = uRef/f_Coriolis_dim/lRef  
          RoInv(:) = 1.0/Ro(:) 
       else
          Ro(:) = 1.d40
          RoInv(:) = 0.0
       end if
    end if

    ! Reynolds number
    if( .not. specifyReynolds ) ReInv = mu_viscous_dim/uRef/lRef

    if ( ReInv < 1.0e-20 ) then
       Re = 1.0e20             ! only used in timestep calculation routine
    else
       Re = 1.0/ReInv
    end if

    ! Heat conduction
    mu_conduct = mu_conduct_dim/uRef/lRef + 1.0e-20

    ! scaled reference pressure at z = 0
    p0 = press0_dim / pRef

    ! scaled background flow
    backgroundFlow = backgroundFlow_dim / uRef

    ! nondimensional gravitational constant
    g_ndim = g/(uRef**2 / lRef)


    !----------------------------------
    !            setup domain
    !----------------------------------

    ! scale the domain by reference length lRef
    lx = lx_dim / lRef
    ly = ly_dim / lRef
    lz = lz_dim / lRef

    ! init cell size

    dx = ( lx(1) - lx(0) ) / real(sizeX)
    dy = ( ly(1) - ly(0) ) / real(sizeY)
    dz = ( lz(1) - lz(0) ) / real(sizeZ)

    ! init cell coordinates
     do i = -nbx, sizeX+nbx   ! modified by Junhong Wei (20161104)
       x(i) = lx(0) + real(i-1)*dx + dx/2.0
    end do

    do j = -nby, sizeY+nby   ! modified by Junhong Wei (20161104)
       y(j) = ly(0) + real(j-1)*dy + dy/2.0
    end do

    do k = -nbz, sizeZ+nbz   ! modified by Junhong Wei (20161104)
       z(k) = lz(0) + real(k-1)*dz + dz/2.0
    end do

!   achatzb
    !----------------------------------
    !            setup topography
    !----------------------------------

    mountainHeight = mountainHeight_dim/lRef
    mountainWidth  = mountainWidth_dim/lRef

    x_center = 0.5*(lx(1) + lx(0)) 
    y_center = 0.5*(ly(1) + ly(0)) 

    k_mountainw = pi/mountainWidth

    do j = -nby+1, sizeY+nby   ! modified by Junhong Wei (20161104)
       do i = -nbx+1, sizeX+nbx
          if(abs(x(i)-x_center) <= mountainWidth) then
             if(mountain_case == 1) then
                topography_surface(i,j) &
                = 0.5*mountainHeight &
                    * (1. + cos(k_mountainw*(x(i)-x_center)))
              else if(mountain_case == 2) then ! modififed by FDK (20191216)
                topography_surface(i,j) &
                = 0.5*mountainHeight &
                  * (1. + cos(k_mountainw*(x(i)-x_center))) &
                  * (cos(range_fac*k_mountainw*(x(i)-x_center)))
                  if (topography_surface(i,j) .lt. lz(0) .and. topography ) then
                    print*,"Topography has negative values"
                    print*, "lz_dim(0) should be set accordingly "
                    stop
                  end if
            end if
           else
             topography_surface(i,j) = 0.
          end if
       end do
    end do

    do k = -nbz+1, sizeZ+nbz
       do j = -nby+1, sizeY+nby
          do i = -nbx+1, sizeX+nbx
             if(z(k) < topography_surface(i,j)) then        
                topography_mask(i,j,k) = .true.
               else
                topography_mask(i,j,k) = .false.
             end if
          end do
       end do
    end do
!   achatze



    !---------------------------------------------
    !   Set up Sponge layer
    !---------------------------------------------

    if( spongeLayer ) then
       kSponge = nz - ceiling(spongeHeight*real(nz))
       zSponge = z(kSponge)
    end if

    select case( model ) 
       !--------------------------------------------------------------------
       !             Atmospheres for pseudo-incompressible
       !--------------------------------------------------------------------

    case( "pseudo_incompressible", "WKB" )


       select case (background)


          !-----------------------------------------------------------
          !   Isentropic troposphere / isothermal stratosphere
          !-----------------------------------------------------------

       case( 'realistic' )
          
          ! not implemented versions
          if ( referenceQuantities == "SI" ) &
               & stop "atmosphere.f90: referenceQuantities = SI not impl.."
          if( .not. fluctuationMode ) &
               &  stop "atmosphere.f90: only fluctuationMode = TRUE impl.."

          ! quantities at tropopause
          z_tr = z_tr_dim / lRef
          theta_tr = theta_tr_dim / thetaRef
          ! presure
          press_tr = p0*(1.0-kappa*sig/theta_tr * z_tr)**(1/kappa)
          ! temperature
          T_tr = theta_tr * (press_tr/p0)**kappa
          
          ! quantities at full levels
          do k = -1,nz+2
             zk = z(k)
             delZ = zk - z_tr
             
             if( zk <= z_tr ) then ! isentropic troposphere

                thetaStrat(k) = theta_tr

                power = 1.0/(gamma-1.0)
                term = kappa*sig/theta_tr * zk
                if (term > 1.0)stop "init_atmosphere: negative term. Stop."
                
                Pstrat(k) = p0*( 1.0 - term)**power
                rhoStrat(k) = PStrat(k) / thetaStrat(k)
                
             else                  ! isothermal stratosphere
                
                thetaStrat(k) = theta_tr * exp(kappa*sig/T_tr*delZ)
                Pstrat(k) &
                = p0**kappa * press_tr**(1/gamma) &
                  * exp(-sig/gamma/T_tr*delZ)
                rhoStrat(k) = PStrat(k) / thetaStrat(k)

             end if

          end do ! k loop
          
          ! quantities at half levels
          do k = -1,nz+2
             if (k == nz+2) then
                zk_half = z(k) + dz/2.
                delZ = zk_half - z_tr

                if( zk_half <= z_tr ) then ! isentropic troposphere
                   thetaStratTilde(k) = theta_tr

                   power = 1.0/(gamma-1.0)
                   term = kappa*sig/theta_tr * zk_half
                   if (term > 1.0) stop "init_atmosphere: negative term. Stop."

                   PstratTilde(k) = p0 * (1.0 - term)**power
                   rhoStratTilde(k) = PStratTilde(k) / thetaStratTilde(k)
                  else                  ! isothermal stratosphere
                   thetaStratTilde(k) = theta_tr * exp(kappa*sig/T_tr*delZ)
                   PstratTilde(k) &
                   = p0**kappa * press_tr**(1/gamma) &
                     * exp(-sig/gamma/T_tr*delZ)
                   rhoStratTilde(k) = PStratTilde(k) / thetaStratTilde(k)
                end if
               else
                thetaStratTilde(k) = 0.5 * (thetaStrat(k) + thetaStrat(k+1))
                PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k+1))
                rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k+1))
             end if
          end do ! k loop 
          
          !-----------------------------------------------------------
          !                    Isothermal atmosphere
          !-----------------------------------------------------------

       case( 'isothermal' )     

          if ( referenceQuantities == "SI" ) then  
             !------------------------------------
             !   original equations in SI units
             !------------------------------------

             if( fluctuationMode ) &
                  &stop "init_atmosphere: fluctuationMode not impl. for SI!" 
             T0 = Temp0_dim / thetaRef   ! T0 in K
             N2 = kappa*g**2/Rsp/T0      ! isothermal Brunt-Vaisala fr.^2
             NN = sqrt(N2)               ! 

             hp = Rsp*T0/g               ! pressure scale height
             hTheta = hp/kappa           ! pot. temperature scale height

             do k = -1,nz+2
                PStrat(k) = p0*exp( -z(k)/gamma/hp )
                thetaStrat(k) = T0*exp( z(k) / hTheta )
                rhoStrat(k) = 1.0/Rsp * PStrat(k) / thetaStrat(k)
             end do

          else
             !-----------------------------------------
             !    with reference quantities 
             !-----------------------------------------

             T0 = Temp0_dim / thetaRef      
             N2 = Ma**2/Fr**4*kappa/T0   ! isothermal Brunt-Vaisala fr.^2
             NN = sqrt(N2)                 

             do k = -1,nz+2
                PStrat(k) = p0*exp( -sig*z(k) / gamma / T0 )       
                thetaStrat(k) = T0 * exp( kappa*sig/T0 * z(k) )
                rhoStrat(k) = PStrat(k) / thetaStrat(k)
             end do

             ! rhoStrat and PStrat at half levels
             do k = -1,nz+2
                if (k == nz+2) then
                   zk_half = z(k) + 0.5*dz  ! half level height
                   PStratTilde(k) = p0*exp( -sig*zk_half / gamma / T0 )
                   thetaStratTilde(k) = T0 * exp( kappa*sig/T0 * zk_half )
                   rhoStratTilde(k) = PStratTilde(k) / thetaStratTilde(k)
                  else
                   thetaStratTilde(k) &
                   = 0.5 * (thetaStrat(k) + thetaStrat(k+1))
                   PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k+1))
                   rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k+1))
                end if
             end do

             !xxx need rhoStratTilde(-1) in fluxes.f90, line 1927 \pm
             rhoStratTilde(-1) = rhoStratTilde(0)
             rhoStratTilde(nz+1) = rhoStratTilde(nz)

          end if

          bvsStrat = N2

          !-------------------------------------------------------------------
          !                        Isentropic atmosphere
          !-------------------------------------------------------------------

       case( 'isentropic' )
       
          if (include_ice .and. (iceTestcase=="1D_ISSR")) then
             theta0_dim = T_nuc + 5.0
             press0_dim = p_nuc * (theta0_dim/T_nuc)**kappaInv
             ! scaled reference pressure at z = 0
             p0 = press0_dim / pRef
          end if
             
          if ( referenceQuantities == "SI" ) then  
             !------------------------------------
             !   original equations in SI units
             !------------------------------------

             if( fluctuationMode ) then
             stop "init_atmosphere: fluctuationMode not implmented for SI!" 
             end if

             NN = 0.0
             N2 = 0.0
             theta0 = theta0_dim / thetaRef       ! theta0 in K
             power = 1.0/(gamma-1.0)
             do k = -1,nz+2
                term = kappa*g/Rsp/theta0*z(k)
                if (term > 1.0) then
                   stop "init_atmosphere: negative term with power.Stop."
                end if
                PStrat(k) = p0 * (1.0 - term)**power
                thetaStrat(k) = theta0
                rhoStrat(k) = 1.0/Rsp * PStrat(k) / thetaStrat(k)
             end do

          else
             !-----------------------------------------
             !    with reference quantities 
             !-----------------------------------------
             thetaStrat = theta0_dim / thetaRef
             theta0 = theta0_dim / thetaRef
             NN = 0.0
             N2 = 0.0
             ! hydrostatic background pressure


             do k = -1,nz+2
                power = 1.0/(gamma-1.0)
                term = kappa*sig/theta0 * z(k)
                if (term > 1.0) then
                   print*, "init_atmosphere: neg. term with power.Stopping."
                   print*,"term = ", term
                   print*,"kappa = ", kappa
                   print*,"sig = ", sig
                   print*,"theta0 = ", theta0
                   print*,"z(nz) = ", z(k)
                   print*,"z(nz)*l = ", z(k)
                   print*,"lRef = ", lRef
                   stop "stopping."                
                end if
                Pstrat(k) = p0*( 1.0 - term)**power

                rhoStrat(k) = PStrat(k) / theta0

             end do

             ! quantities at half levels
             do k = -1,nz+2
                if (k == nz+2) then
                   zk_half = z(k) + dz/2.
                   power = 1.0/(gamma-1.0)
                   term = kappa*sig/theta0 * zk_half
                   if (term > 1.0) then
                      print*, "init_atmosphere: neg. term with power."
                      print*,"term = ", term
                      print*,"kappa = ", kappa
                      print*,"sig = ", sig
                      print*,"theta0 = ", theta0
                      print*,"z(nz) = ", z(k)
                      print*,"z(nz)*l = ", z(k)
                      print*,"lRef = ", lRef
                      stop "stopping."                
                   end if
                   PstratTilde(k) = p0*( 1.0 - term)**power
                   rhoStratTilde(k) = PStratTilde(k) / theta0
                  else
                   PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k+1))
                   rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k+1))
                end if
             end do


          end if

          bvsStrat = N2

          !----------------------------------------------------------------
          !                   Stable atmosphere with constant N
          !----------------------------------------------------------------

       case( 'const-N' )
       
          if (include_ice .and. (iceTestcase=="1D_ISSR")) then
             theta0_dim = T_nuc + 5.0
             term = kappa*g**2 / (Rsp*N_BruntVaisala_dim**2)
             press0_dim = p_nuc / (1.+term/theta0_dim*((theta0_dim-term)/(T_nuc-term)-1.))**kappaInv
             ! scaled reference pressure at z = 0
             p0 = press0_dim / pRef
          end if

          if ( referenceQuantities == "SI" ) then  
             !------------------------------------
             !   original equations in SI units
             !------------------------------------
             if( fluctuationMode ) then
               stop "init_atmosphere: fluctuationMode not implmented for SI!"
             end if

             theta0 = theta0_dim / thetaRef       ! theta0 at z=0 in K
             NN = N_BruntVaisala_dim * tRef
             N2 = NN**2  ! Brunt-Vaisala N0^2 in 1/s^2

             coeff = kappa * g**2 / (Rsp * N2 * theta0)
             power = 1./(gamma-1.)

             do k = -1,nz+2
                term =  exp( -N2/g*z(k) )
                
                !if( term > 1.0 ) stop "init_atmosphere: root of a negative number." 
                if (1.0 + coeff*( term - 1.0) < 0.0) stop "init_atmosphere: root of a negative number."

                PStrat(k) = p0 * (1.0 + coeff*( term - 1.0) )**power
                thetaStrat(k) = theta0 * exp(N2/g*z(k))
                rhoStrat(k) = 1.0/Rsp * PStrat(k) / thetaStrat(k)
             end do

          else
             !-----------------------------------------
             !    with reference quantities 
             !-----------------------------------------

             theta0 = theta0_dim / thetaRef         ! theta0 at z=0
             N2 = (N_BruntVaisala_dim * tRef)**2    ! Brunt-Vaisala fr.^2
             NN = sqrt(N2)

             power = 1.0 / (gamma-1.0)

             ! potential temperature and pressure Pstrat
             do k = -1,nz+2
                thetaStrat(k) = theta0 * exp(Fr**2*N2*z(k))

                term = exp( -Fr2 * N2 * z(k) )

                if( 1.0 + FrInv2*kappa*sig/N2/theta0*(term-1.0) < 0.0) then
                    print*,"init_atmosphere: power of a neg. number." 
                    stop 'top of atmosphere too high for const-N'
                end if

                Pstrat(k) &
                = p0 &
                  * (1.0 + FrInv2*kappa*sig/N2/theta0 * (term-1.0) )**power

                rhoStrat(k) = pStrat(k) / thetaStrat(k)             
             end do

             ! rhoStrat at half levels
             do k = -1,nz+2
                if (k == nz+2) then
                   zk_half = z(k) + dz/2.
                   thetaStratTilde(k) = theta0 * exp(Fr**2*N2*zk_half)

                   term = exp( -Fr2 * N2 * zk_half )

                   if( term > 1.0 ) then
                       stop "init_atmosphere: root of a negative number." 
                   end if

                   if( 1.0 + FrInv2*kappa*sig/N2/theta0*(term-1.0) < 0.0) &
                   & then
                      print*,"init_atmosphere: power of a neg. number." 
                      stop 'top of atmosphere too high for const-N'
                   end if

                   PstratTilde(k) &
                   = p0 &
                     * (  1.0 &
                        + FrInv2*kappa*sig/N2/theta0 * (term-1.0) )**power
                   rhoStratTilde(k) = PstratTilde(k) / thetaStratTilde(k)
                  else
                   thetaStratTilde(k) &
                   = 0.5 * (thetaStrat(k) + thetaStrat(k+1))
                   PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k+1))
                   rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k+1))
                end if
             end do

          end if

          bvsStrat = N2

         !-----------------------------------------------------------
         ! Setting for a baroclinic life cycle. Different lapse rates 
         ! in troposphere and stratosphere that are also y-dependent.
         ! Change: Elena Gagarina, 15.03.2018
         !-----------------------------------------------------------

       case( 'diflapse' )

          if ( referenceQuantities == "SI" ) then
             stop "atmosphere.f90: referenceQuantities = SI not impl."
          end if

          ! quantities at tropopause

          if (gamma_t /= 0.) pow_t = g/(Rsp*gamma_t)
          if (gamma_s /= 0.) pow_s = g/(Rsp*gamma_s)

          if (gamma_t /= 0.) then
             p_t_b = press0_dim * (1. - gamma_t*z_tr_dim/Temp0_dim)**pow_t
            else
             p_t_b = press0_dim * exp(- z_tr_dim/(Rsp*Temp0_dim/g))
          end if

          T_c_b = Temp0_dim - gamma_t*z_tr_dim
          
          ! quantities at full levels:

          do k = -1, nz+2
             zk = z(k) * lRef

             if (zk < z_tr_dim) then
                T_bar = Temp0_dim - gamma_t*zk

                if (gamma_t /= 0.) then
                   p_bar = press0_dim * (1. - gamma_t*zk/Temp0_dim)**pow_t 
                  else
                   p_bar = press0_dim * exp(- zk/(Rsp*Temp0_dim/g))
                end if
               else
                T_bar &
                = Temp0_dim - gamma_t*z_tr_dim - gamma_s*(zk - z_tr_dim)

                if (gamma_s /= 0.) then
                   p_bar &
                   = p_t_b * (1. - gamma_s*(zk - z_tr_dim)/T_c_b)**pow_s
                  else
                   p_bar &
                   = p_t_b * exp(- (zk - z_tr_dim)/(Rsp*T_c_b/g))
                end if
             endif

             thetaStrat(k) = T_bar * (press0_dim/p_bar)**kappa / thetaRef

             rhoStrat(k) = p_bar/(Rsp*T_bar) / rhoRef

             Pstrat(k) = rhoStrat(k)*thetaStrat(k)
          enddo
          
          ! quantities at half levels

          ! with the exception of the uppermost ghost layer the
          ! values there are obtained by linear interpolation, in order
          ! to be consistent with the handling of the half levels in the
          ! semi-implicit time stepping and the corresponding pressure
          ! solver
          
          do k = -1, nz+2
             if (k == nz+2) then
                zk_half = (z(k) + dz/2.) * lRef

                if (zk_half < z_tr_dim) then
                   T_bar = Temp0_dim - gamma_t*zk_half

                   if (gamma_t /= 0.) then
                      p_bar &
                      = press0_dim &
                        * (1. - gamma_t*zk_half/Temp0_dim)**pow_t 
                     else
                      p_bar = press0_dim * exp(- zk_half/(Rsp*Temp0_dim/g))
                   end if
                  else
                   T_bar &
                   =   Temp0_dim - gamma_t*z_tr_dim &
                     - gamma_s*(zk_half - z_tr_dim)

                   if (gamma_s /= 0.) then
                      p_bar &
                      = p_t_b &
                        * (1. - gamma_s*(zk_half - z_tr_dim)/T_c_b)**pow_s
                     else
                      p_bar &
                      = p_t_b * exp(- (zk_half - z_tr_dim)/(Rsp*T_c_b/g))
                   end if
                endif 

                thetaStratTilde(k) &
                = T_bar * (press0_dim/p_bar)**(kappa) / thetaRef

                rhoStratTilde(k) = p_bar/(Rsp*T_bar) / rhoRef

                PstratTilde(k) = rhoStratTilde(k)*thetaStratTilde(k)
               else
                PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k+1))

                rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k+1))

                thetaStratTilde(k) = PStratTilde(k)/rhoStratTilde(k)
             end if
          enddo 

       !UAB
         !-----------------------------------------------------------
         ! Setting for an atmmosphere according to Held & Suarez (1994)
         !-----------------------------------------------------------

       case( 'HeldSuarez' )

          if ( referenceQuantities /= "Klein" ) then
             stop "atmosphere.f90: referenceQuantities = Klein expected!"
          end if

          ! non-dimensional parameters of the Held-Suarez reference 
          ! atmosphere:

          ! stratosphere temperature
          tp_strato = tp_strato_dim/thetaRef
          ! tropical surface temperature
          tp_srf_trp = tp_srf_trp_dim/thetaRef 
          ! tropospheric temperature difference between poles and tropics
          tpdiffhor_tropo = tpdiffhor_tropo_dim/thetaRef
          ! vertical potential-temperature difference in troposphere
          ptdiffvert_tropo = ptdiffvert_tropo_dim/thetaRef

          ! Exner pressure just above and below the surface so that it 
          ! = 1 at the surface, by an Euler integration of hydrostatic 
          ! equilibrium

          T_bar = max( tp_strato, tp_srf_trp - 0.5*tpdiffhor_tropo)

          !testb
          !print*,tp_strato, tp_srf_trp - 0.5*tpdiffhor_tropo, T_bar
          !teste

          pistrat(0) = 1.0 + 0.5*dz * kappa/T_bar
          pistrat(1) = 1.0 - 0.5*dz * kappa/T_bar

          ! potential temperature just below and above the surface
          ! P and density determined as well

          do k=0,1
             T_bar &
             = max( tp_strato, &
                    pistrat(k) &
                    * (tp_srf_trp - 0.5*tpdiffhor_tropo &
                       - 0.5*ptdiffvert_tropo/kappa * log(pistrat(k))))
             
             thetaStrat(k) = T_bar/pistrat(k)

             pStrat(k) = pistrat(k)**((1.0 - kappa)/kappa)

             rhoStrat(k) = pStrat(k)/thetaStrat(k)
          end do
          
          ! for k > 1:
          ! Exner pressure and potential temperature by upward 
          ! integration of hydrostatic equilibrium, using a trapezoidal 
          ! leapfrog
          ! P and density determined as well

          do k = 2, nz-ceiling((spongeHeight+spongeHeight/3.)*real(nz)) !nz+2 !FS
             pistar = pistrat(k-2) - 2.0*dz * kappa/thetaStrat(k-1)

             T_bar &
             = max( tp_strato, &
                    pistar &
                    * (tp_srf_trp - 0.5*tpdiffhor_tropo &
                       - 0.5*ptdiffvert_tropo/kappa * log(pistar))) 

             thetastar = T_bar/pistar

             pistrat(k) &
             =   pistrat(k-1) &
               - 0.5*dz * (kappa/thetastar + kappa/thetaStrat(k-1))

             T_bar &
             = max( tp_strato, &
                    pistrat(k) &
                    * (tp_srf_trp - 0.5*tpdiffhor_tropo &
                       - 0.5*ptdiffvert_tropo/kappa * log(pistrat(k))))

             thetaStrat(k) = T_bar/pistrat(k)

             pStrat(k) = pistrat(k)**((1.0 - kappa)/kappa)

             rhoStrat(k) = pStrat(k)/thetaStrat(k)
          end do

          tp_sponge = T_bar

          ! close jets below sponge layer !FS
          do k = nz+1-ceiling((spongeHeight+spongeHeight/3.)*real(nz)), nz+2 !FS
             pistar = pistrat(k-2) - 2.0*dz * kappa/thetaStrat(k-1)

             T_bar &
             =  tp_sponge + pistar*(0.5*tpdiffhor_tropo &
                       + 0.5*ptdiffvert_tropo/kappa * log(pistar)) 
      

             thetastar = T_bar/pistar

             pistrat(k) &
             =   pistrat(k-1) &
               - 0.5*dz * (kappa/thetastar + kappa/thetaStrat(k-1))

             T_bar &
             = tp_sponge + pistrat(k)*(0.5*tpdiffhor_tropo &
                       + 0.5*ptdiffvert_tropo/kappa * log(pistrat(k)))

             thetaStrat(k) = T_bar/pistrat(k)

             pStrat(k) = pistrat(k)**((1.0 - kappa)/kappa)

             rhoStrat(k) = pStrat(k)/thetaStrat(k)
          end do
          
          


          ! quantities at half levels

          ! with the exception of the uppermost ghost layer the
          ! values there are obtained by linear interpolation, in order
          ! to be consistent with the handling of the half levels in the
          ! semi-implicit time stepping and the corresponding pressure
          ! solver
          
          do k = -1, nz+2
             if (k < nz+2) then
                PStratTilde(k) = 0.5 * (PStrat(k) + PStrat(k+1))

                rhoStratTilde(k) = 0.5 * (rhoStrat(k) + rhoStrat(k+1))

                thetaStratTilde(k) = PStratTilde(k)/rhoStratTilde(k)
               else
                ! linear extrapolation at uppermost layer

                PStratTilde(k) = 2.0*PStratTilde(k-1) - PStratTilde(k-2)

                rhoStratTilde(k) &
                = 2.0*rhoStratTilde(k-1) - rhoStratTilde(k-2)

                thetaStratTilde(k) &
                = 2.0*thetaStratTilde(k-1) - thetaStratTilde(k-2)
             end if
          enddo 

        
       !------------------------------------------------------------------

       case default
          print*,"background = ", trim(background)
          stop "atmosphere.f90/init_background: background not defined"
       end select

       ! non-dimensional squared Brunt-Vaisala frequency
       ! (this could be done a bit nicer)
          
       bvsStrat(-1) &
            = g_ndim/thetaStrat(0) * (thetaStrat(1) - thetaStrat(0))/dz
       
       bvsStrat(0) &
            = g_ndim/thetaStrat(0) * (thetaStrat(1) - thetaStrat(0))/dz
       
       !UAB
       N2 = max(bvsStrat(-1),bvsStrat(0))
       !UAE
       
       do k = 1,nz
          bvsStrat(k) &
               = g_ndim/thetaStrat(k) &
               * (thetaStrat(k+1) - thetaStrat(k-1))/(2.0 * dz)
          
          !UAB
          N2 = max(N2, bvsStrat(k))
          !UAE
       end do
       
       bvsStrat(nz+1) &
            = g_ndim/thetaStrat(nz+1) * (thetaStrat(nz+1) - thetaStrat(nz))/dz
       
       !UAB
       N2 = max(N2, bvsStrat(nz+1))
       
       if(N2 < 0.) then
          stop 'ERROR: N2 < 0'
       else
          NN = sqrt(N2)
       end if
       !UAE
       

    case( "Boussinesq" )

       !-------------------------------------------------------------------
       !             Atmospheres for Boussinesq
       !-------------------------------------------------------------------

       select case (background)

       case( 'uniform_Boussinesq' )  ! uniform atmosphere for testing

          rho00   = 1.0
          theta00 = theta0_dim / thetaRef
          P00     = rho00 * theta00
          NN      = 0.0
          N2      = NN**2

          ! some stratified fields are still in use...
          rhoStrat = rho00
          rhoStratTilde = rho00       ! background density at half levels
          thetaStrat = theta00
          PStrat = P00


       case( 'stratified_Boussinesq' )          

          rho00   = 1.0
          theta00 = theta0_dim / thetaRef
          P00     = rho00 * theta00
          NN      = N_BruntVaisala_dim*tRef
          N2      = NN**2

          ! some stratified fields are still in use...
          rhoStrat = rho00
          rhoStratTilde = rho00       ! background density at half levels
          thetaStrat = theta00
          PStrat = P00

       case default
          print*,"background = ", trim(background)
          stop "atmosphere.f90/init_background: background not defined"
       end select

       bvsStrat = N2

    case default
       print*,"model = ", model
       stop "init_atmosphere: unknown case model."
    end select


  end subroutine init_atmosphere


  !---------------------------------------------------------------------------


  subroutine terminate_atmosphere
    !---------------------------------
    ! deallocate background varaibles
    !---------------------------------

    ! local variables
    integer :: allocstat


    !---------------- deallocate variables -----------------------

    deallocate(Pstrat,stat=allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not deallocate pStrat"

    !UAB
    deallocate(pistrat,stat=allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not deallocate pistrat"
    !UAE

    deallocate(PstratTilde,stat=allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not deallocate pStratTilde"

    deallocate(rhoStrat,stat=allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not deallocate pStrat"

    deallocate(rhoStratTilde,stat=allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not deallocate pStrat"

    deallocate(thetaStrat,stat=allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not dealloc thetaStrat"

    deallocate(thetaStratTilde,stat=allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not dealloc thetaStratTilde"

    deallocate(Ro,stat=allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not dealloc Ro"

    deallocate(RoInv,stat=allocstat)
    if(allocstat /= 0) stop "atmosphere.f90: could not dealloc RoInv"

  end subroutine terminate_atmosphere


end module atmosphere_module



