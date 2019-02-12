module atmosphere_module

  use type_module

  implicit none

  public   ! all objects are known to programmes using this module

  !--------------
  !    public
  !--------------
  public :: init_atmosphere
  public :: terminate_atmosphere

  real, dimension(:), allocatable :: PStrat, rhoStrat, thetaStrat
  real, dimension(:), allocatable :: PStratTilde, rhoStratTilde, thetaStratTilde


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
  real, parameter  :: Rsp = 287.0   ! spec. gas constant for dry air in J/kg/K

  ! flow parameters
  real :: Re                     ! Reynolds number (calculated from input 1/Re)
  real :: Ma, MaInv2,Ma2         ! Mach number and 1/Ma^2, Ma^2
  real :: Fr, FrInv2,Fr2         ! Froude number Fr and 1/Fr^2, Fr^2
  real :: sig                    ! Ma^2/Fr^2
  real :: Ro, RoInv              ! Rossby number and its inverse

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

    ! debugging
    integer,parameter :: errorlevel = 10 ! 0 -> no output

    ! allocate PStrat
    if( .not. allocated(pStrat) ) then
       allocate( Pstrat(0:nz+1),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate pStrat"
    end if

    ! allocate pStratTilde -> P at half levels
    if( .not. allocated(pStratTilde) ) then
       allocate( PstratTilde(0:nz),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate pStratTilde"
    end if

    ! allocate rhoStrat
    if( .not. allocated(rhoStrat) ) then
       allocate( rhoStrat(0:nz+1),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate rhoStrat"
    end if

    ! allocate rhoStratTilde 
    if( .not. allocated(rhoStratTilde) ) then
       allocate( rhoStratTilde(-1:nz+1),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate rhoStratTilde"
    end if


    ! allocate thetaStrat
    if( .not. allocated(thetaStrat) ) then
       allocate( thetaStrat(0:nz+1),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate thetaStrat"
    end if

    ! allocate thetaStratTilde
    if( .not. allocated(thetaStratTilde) ) then
       allocate( thetaStratTilde(0:nz+1),stat=allocstat)
       if(allocstat /= 0) stop "atmosphere.f90: could not allocate thetaStratTilde"
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
       eps = 0.1                 ! asymptotic parameter ( Ma = Fr for eps = kappa! )
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
       stop "init_atmosphere: Problems with Exner pressure. Use Klein's scaling."
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
       !                          ! but we use Exner-pressure/kappa like R. Klein
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

    ! Rossby numer
!   achatzb correction for zero Coriolis
!   Ro = uRef/f_Coriolis_dim/lRef  
!   RoInv = 1.0/Ro       
    if(f_Coriolis_dim /= 0.0) then
       Ro = uRef/f_Coriolis_dim/lRef  
       RoInv = 1.0/Ro       
      else
       RoInv=0.0
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


    !----------------------------------
    !            setup domain
    !----------------------------------

    ! scale the domain by reference length lRef
    lx = lx_dim / lRef
    ly = ly_dim / lRef
    lz = lz_dim / lRef

    ! init cell size
!    dx = ( lx(1) - lx(0) ) / real(nx)   ! modified by Junhong Wei (20161104)
!    dy = ( ly(1) - ly(0) ) / real(ny)   ! modified by Junhong Wei (20161104)
!    dz = ( lz(1) - lz(0) ) / real(nz)   ! modified by Junhong Wei (20161104)

     dx = ( lx(1) - lx(0) ) / real(sizeX)   ! modified by Junhong Wei (20161104)
     dy = ( ly(1) - ly(0) ) / real(sizeY)   ! modified by Junhong Wei (20161104)
     dz = ( lz(1) - lz(0) ) / real(sizeZ)   ! modified by Junhong Wei (20161104)

    ! init cell coordinates
!    do i = -nbx, nx+nbx   ! modified by Junhong Wei (20161104)
     do i = -nbx, sizeX+nbx   ! modified by Junhong Wei (20161104)
       x(i) = lx(0) + real(i-1)*dx + dx/2.0
    end do

!    do j = -nby, ny+nby   ! modified by Junhong Wei (20161104)
    do j = -nby, sizeY+nby   ! modified by Junhong Wei (20161104)
       y(j) = ly(0) + real(j-1)*dy + dy/2.0
    end do

!    do k = -nbz, nz+nbz   ! modified by Junhong Wei (20161104)
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
             topography_surface(i,j) &
             = 0.5*mountainHeight*(1. + cos(k_mountainw*(x(i)-x_center)))
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
               & stop "atmosphere.f90: referenceQuantities = SI not implemented."
          if( .not. fluctuationMode ) &
               & stop "atmosphere.f90: only fluctuationMode = TRUE implemented."

          ! quantities at tropopause
          z_tr = z_tr_dim / lRef
          theta_tr = theta_tr_dim / thetaRef
          press_tr = p0*(1.0-kappa*sig/theta_tr * z_tr)**(1/kappa)   ! pressure
          T_tr = theta_tr * (press_tr/p0)**kappa                     ! temperature
          
          ! quantities at full levels
          do k = 0,nz+1
             zk = z(k)
             delZ = zk - z_tr
             
             if( zk <= z_tr ) then ! isentropic troposphere

                thetaStrat(k) = theta_tr

                power = 1.0/(gamma-1.0)
                term = kappa*sig/theta_tr * zk
                if (term > 1.0) stop "init_atmosphere: negative term. Stop."
                
                Pstrat(k) = p0*( 1.0 - term)**power
                rhoStrat(k) = PStrat(k) / thetaStrat(k)
                
             else                  ! isothermal stratosphere
                
                thetaStrat(k) = theta_tr * exp(kappa*sig/T_tr*delZ)
                Pstrat(k) = p0**kappa * press_tr**(1/gamma) * exp(-sig/gamma/T_tr*delZ)
                rhoStrat(k) = PStrat(k) / thetaStrat(k)

             end if

          end do ! k loop
          
          
          if( errorlevel < 5 ) then
             ! xxx check atmosphere at full levels
             print"(4a10)","k","z(k)","thetaStrat(k)","PStrat(k)","rhoStrat(k)"
             do k = 0,nz+1
                
                print"(i10,4f10.1)",k,z(k)*lRef,thetaStrat(k)*thetaRef,&
                     & PStrat(k)*pRef,rhoStrat(k)*rhoRef
                
             end do
             
          end if
          
          ! quantities at half levels
          do k = 0,nz
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
                PstratTilde(k) = p0**kappa * press_tr**(1/gamma) * exp(-sig/gamma/T_tr*delZ)
                rhoStratTilde(k) = PStratTilde(k) / thetaStratTilde(k)
                
             end if

          end do ! k loop 
          

          if( errorlevel < 5 ) then
             ! xxx check atmosphere at half levels
             print"(4a10)","k","z(k)+dz","thetaStratTilde(k)",&
                  & "PStratTilde(k)","rhoStratTilde(k)"
             do k = 0,nz
                
                print"(i10,4f10.1)",k,(z(k)+dz/2.)*lRef,thetaStratTilde(k)*thetaRef,&
                     & PStratTilde(k)*pRef,rhoStratTilde(k)*rhoRef
                
             end do
             
          end if

          


          !-----------------------------------------------------------
          !                    Isothermal atmosphere
          !-----------------------------------------------------------

       case( 'isothermal' )     

          if ( referenceQuantities == "SI" ) then  
             !------------------------------------
             !   original equations in SI units
             !------------------------------------

             if( fluctuationMode ) &
                  & stop "init_atmosphere: fluctuationMode not implmented for SI!" 
             T0 = Temp0_dim / thetaRef   ! T0 in K
             N2 = kappa*g**2/Rsp/T0      ! isothermal Brunt-Vaisala frequency^2
             NN = sqrt(N2)               ! 

             hp = Rsp*T0/g               ! pressure scale height
             hTheta = hp/kappa           ! pot. temperature scale height

             do k = 0,nz+1
                PStrat(k) = p0*exp( -z(k)/gamma/hp )
                thetaStrat(k) = T0*exp( z(k) / hTheta )
                rhoStrat(k) = 1.0/Rsp * PStrat(k) / thetaStrat(k)
             end do

          else
             !-----------------------------------------
             !    with reference quantities 
             !-----------------------------------------

             T0 = Temp0_dim / thetaRef      
             N2 = Ma**2/Fr**4*kappa/T0   ! isothermal Brunt-Vaisala frequency^2
             NN = sqrt(N2)                 

             do k = 0,nz+1
                PStrat(k) = p0*exp( -sig*z(k) / gamma / T0 )       
                thetaStrat(k) = T0 * exp( kappa*sig/T0 * z(k) )
                rhoStrat(k) = PStrat(k) / thetaStrat(k)
             end do

             ! rhoStrat and PStrat at half levels
             do k = 0,nz
                zk_half = z(k) + 0.5*dz                           ! half level height
                PStratTilde(k) = p0*exp( -sig*zk_half / gamma / T0 )       
                thetaStratTilde(k) = T0 * exp( kappa*sig/T0 * zk_half )
                rhoStratTilde(k) = PStratTilde(k) / thetaStratTilde(k)
             end do

!xxx need rhoStratTilde(-1) in fluxes.f90, line 1927 \pm
             rhoStratTilde(-1) = rhoStratTilde(0)
             rhoStratTilde(nz+1) = rhoStratTilde(nz)

          end if


          !-------------------------------------------------------------------
          !                        Isentropic atmosphere
          !-------------------------------------------------------------------

       case( 'isentropic' )

          if ( referenceQuantities == "SI" ) then  
             !------------------------------------
             !   original equations in SI units
             !------------------------------------

             if( fluctuationMode ) stop "init_atmosphere: fluctuationMode not implmented for SI!" 
             NN = 0.0
             N2 = 0.0
             theta0 = theta0_dim / thetaRef       ! theta0 in K
             power = 1.0/(gamma-1.0)
             do k = 0,nz+1
                term = kappa*g/Rsp/theta0*z(k)
                if (term > 1.0) stop "init_atmosphere: negative term with power.Stop."
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


             do k = 0,nz+1
                power = 1.0/(gamma-1.0)
                term = kappa*sig/theta0 * z(k)
                if (term > 1.0) then
                   print*, "init_atmosphere: negative term with power.Stopping."
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
             do k = 0,nz
                zk_half = z(k) + dz/2.
                power = 1.0/(gamma-1.0)
                term = kappa*sig/theta0 * zk_half
                if (term > 1.0) then
                   print*, "init_atmosphere: negative term with power.Stopping."
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
             end do


          end if

          !--------------------------------------------------------------------
          !                   Stable atmosphere with constant N
          !--------------------------------------------------------------------

       case( 'const-N' )

          if ( referenceQuantities == "SI" ) then  
             !------------------------------------
             !   original equations in SI units
             !------------------------------------
             if( fluctuationMode ) stop "init_atmosphere: fluctuationMode not implmented for SI!" 

             theta0 = theta0_dim / thetaRef       ! theta0 at z=0 in K
             NN = N_BruntVaisala_dim * tRef
             N2 = NN**2  ! Brunt-Vaisala N0^2 in 1/s^2

             coeff = kappa * g**2 / (Rsp * N2 * theta0)
             power = 1./(gamma-1.)

             do k = 1,nz+1
                term =  exp( -N2/g*z(k) )
                if( term > 1.0 ) stop "init_atmosphere: root of a negative number." 
                PStrat(k) = p0 * (1.0 + coeff*( term - 1.0) )**power
                thetaStrat(k) = theta0 * exp(N2/g*z(k))
                rhoStrat(k) = 1.0/Rsp * PStrat(k) / thetaStrat(k)
             end do

          else
             !-----------------------------------------
             !    with reference quantities 
             !-----------------------------------------

             theta0 = theta0_dim / thetaRef         ! theta0 at z=0
             N2 = (N_BruntVaisala_dim * tRef)**2    ! Brunt-Vaisala frequency^2
             NN = sqrt(N2)

             power = 1.0 / (gamma-1.0)

             ! potential temperature and pressure Pstrat
             do k = 1, nz+1
                thetaStrat(k) = theta0 * exp(Fr**2*N2*z(k))

                term = exp( -Fr2 * N2 * z(k) )

                if( term > 1.0 ) stop "init_atmosphere: root of a negative number." 
                Pstrat(k) = p0 * (1.0 + FrInv2*kappa*sig/N2/theta0 * (term-1.0) )**power
                rhoStrat(k) = pStrat(k) / thetaStrat(k)             
             end do


             !xxx treat exception for k = 0
             thetaStrat(0) = thetaStrat(1)
             Pstrat(0) = Pstrat(1)
             rhoStrat(0) = rhoStrat(1)

             ! rhoStrat at half levels
             do k = 0, nz
                zk_half = z(k) + dz/2.
                thetaStratTilde(k) = theta0 * exp(Fr**2*N2*zk_half)

                term = exp( -Fr2 * N2 * zk_half )

                if( term > 1.0 ) stop "init_atmosphere: root of a negative number." 
                PstratTilde(k) = p0 * (1.0 + FrInv2*kappa*sig/N2/theta0 * (term-1.0) )**power
                rhoStratTilde(k) = PstratTilde(k) / thetaStratTilde(k)             
             end do

          end if

       case default
          print*,"background = ", trim(background)
          stop "atmosphere.f90/init_background: background not defined"
       end select


    case( "Boussinesq" )

       !--------------------------------------------------------------------
       !             Atmospheres for Boussinesq
       !--------------------------------------------------------------------

       select case (background)

       case( 'uniform_Boussinesq' )          ! uniform atmosphere for testing

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

  end subroutine terminate_atmosphere





end module atmosphere_module



