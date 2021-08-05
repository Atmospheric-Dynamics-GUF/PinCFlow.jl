module timeScheme_module

  use type_module

  implicit none

  public

  !------------------------
  !   public subroutines
  !------------------------
  public :: init_timeScheme



  !-------------------------------
  !    public module variables
  !------------------------------

  integer :: RKstage
  real, dimension(3), public :: alpha  
  real, dimension(3), public :: beta
  real, dimension(3), public :: stepFrac
  real, dimension(3,3), public :: rk



contains
  

  subroutine init_timeScheme
    !----------------------------------------------------
    ! initialize parameters for Runge-Kutta time schemes
    !----------------------------------------------------

    ! local variables
    real :: z1,z2,z3,z4,z5,z6,c2       ! TVD-RK3

    ! -------------------------------
    !        setup time scheme 
    ! -------------------------------

    select case (timeScheme)

    case ("semiimplicit")  ! ref. Durran (????), Benaccio & Klein (2019), 
                           ! and Achatz et al (????)
       if(verbose) then
          print*,"update.f90/init_update: timeScheme = semiimplicit"
       end if
       alpha = (/ 0.0,   -5./9., -153./128. /)
       beta =  (/ 1./3., 15./16.,   8./15. /)
       stepFrac = (/ 1./3., 5./12., 1./4. /) 
       ! set phi = t, F(phi) = 1 in time scheme to calc stepFrac
       nStages = 3
       timeSchemeType = "lowStorage"

    case ("LS_Will_RK3")  ! ref. Durran 
       if(verbose) print*,"update.f90/init_update: timeScheme = LS_Will_RK3"
       alpha = (/ 0.0,   -5./9., -153./128. /)
       beta =  (/ 1./3., 15./16.,   8./15. /)
       stepFrac = (/ 1./3., 5./12., 1./4. /) 
       ! set phi = t, F(phi) = 1 in time scheme to calc stepFrac
       nStages = 3
       timeSchemeType = "lowStorage"


    case ("Euler")
       if(verbose) print*,"update.f90/init_update: timeScheme = Euler"
       alpha = (/ 0.0, 0.0, 0.0 /)
       beta =  (/ 1.0, 0.0, 0.0 /)
       stepFrac = (/ 1.0, 0.0, 0.0 /)
       nStages = 1
       timeSchemeType = "lowStorage"


    case ("LS_TVD_RK3") ! ref. Gottlieb, Shu 1998   
       c2 = 0.924574
       z1 = sqrt(36.*c2**4 + 36.*c2**3 - 135.*c2**2 + 84.*c2 - 12.)
       z2 = 2.*c2**2 + c2 - 2.
       z3 = 12.*c2**4 - 18.*c2**3 + 18.*c2**2 - 11.*c2 + 2.
       z4 = 36.*c2**4 - 36.*c2**3 + 13.*c2**2 - 8.*c2 + 4.
       z5 = 69.*c2**3 - 62.*c2**2 + 28.*c2 - 8.
       z6 = 34.*c2**4 - 46.*c2**3 + 34.*c2**2 - 13.*c2 + 2.

       alpha(1) = 0.0
       alpha(2) = (-z1*(6.*c2**2 - 4.*c2 + 1.) + 3.*z3) / &
            & ( (2.*c2 + 1.)*z1 - 3.*(c2 + 2.)*(2.*c2-1.)**2   )
       alpha(3) = (-z4*z1 + 108.*(2.*c2 - 1.)*c2**5 - 3.*(2.*c2 - 1.)*z5) / &
            & ( 24.*z1*c2*(c2-1.)**4 + 72.*c2*z6 + 72.*c2**6*(2.*c2-13.) )

       beta(1) = c2
       beta(2) = ( 12.*c2*(c2-1.)*(3.*z2-z1) - (3.*z2-z1)**2 ) / &
            & ( 144.*c2*(3.*c2-2.)*(c2-1.)**2  )
       beta(3) = ( -24.*(3.*c2-2.)*(c2-1.)**2  ) / &
            & ( (3.*z2 - z1)**2 - 12.*c2*(c2-1.)*(3.*z2-z1)  )

       stepFrac(1) = beta(1)
       stepFrac(2) = beta(2)*(alpha(2) + 1.)
       stepFrac(3) = beta(3)*(alpha(3)*(alpha(2) + 1.) + 1.)

       nStages = 3
       timeSchemeType = "lowStorage"


    case( "CL_TVD_RK3" )
       rk(:,1) = (/ 1.    , 0.    , 1.    /)
       rk(:,2) = (/ 3./4. , 1./4. , 1./4. /)
       rk(:,3) = (/ 1./3. , 2./3. , 2./3. /)
       nStages = 3
       timeSchemeType = "classical"


    case default
       stop "update.f90/init_update: unknown timeScheme."

    end select


  end subroutine init_timeScheme


end module timeScheme_module
