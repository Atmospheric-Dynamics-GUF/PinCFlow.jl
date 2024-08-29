module flux_module

  use type_module
  use xweno_module
  use muscl_module
  use atmosphere_module
  use algebra_module
  use ice_module
  use sizeof_module

  implicit none

  private
  ! all module objects (variables, functions, subroutines)
  ! are internal to the module by default

  ! public routines
  public :: reconstruction
  public :: massFlux
  public :: thetaFlux
  public :: momentumFlux
  public :: iceFlux
  public :: volumeForce
  public :: init_fluxes
  public :: terminate_fluxes
  public :: thetaSource
  public :: massSource
  public :: momentumSource
  public :: iceSource
  public :: tracerFlux
  public :: ice2Flux

  ! TFC FJ
  public :: setHalosOfField
  public :: momentumFluxTestTFC, massFluxTestTFC, reconstructionTestTFC

  ! private routines
  private :: absDiff
  private :: slopeFunction

  ! internal module variables
  real, dimension(:, :, :), allocatable :: rhoBar, rhopBar, rhoOld, rhopOld
  real, dimension(:, :, :), allocatable :: uBar
  real, dimension(:, :, :), allocatable :: vBar
  real, dimension(:, :, :), allocatable :: wBar
  real, dimension(:, :, :), allocatable :: PBar
  real, dimension(:, :, :), allocatable :: thetaBar
  real, dimension(:, :, :), allocatable :: nAerBar
  real, dimension(:, :, :), allocatable :: nIceBar
  real, dimension(:, :, :), allocatable :: qIceBar
  real, dimension(:, :, :), allocatable :: qvBar
  real, dimension(:, :, :), allocatable :: tracerBar
  real, dimension(:, :, :), allocatable :: IceBar

  ! TFC FJ
  ! Needed for semi-implicit time scheme in TFC.
  real, dimension(:, :, :), allocatable :: uOldTFC, vOldTFC, wOldTFC

  ! SK: Needed for compressible explicit Euler
  real, dimension(:, :, :), allocatable :: pinew

  ! if reconstType = MUSCL then uTilde, vTilde, and wTilde are the
  ! reconstructed specific momenta

  real, dimension(:, :, :, :, :), allocatable :: rhoTilde
  real, dimension(:, :, :, :, :), allocatable :: rhoTilde_mom !UA
  real, dimension(:, :, :, :, :), allocatable :: rhopTilde
  real, dimension(:, :, :, :, :), allocatable :: uTilde
  real, dimension(:, :, :, :, :), allocatable :: vTilde
  real, dimension(:, :, :, :, :), allocatable :: wTilde
  real, dimension(:, :, :, :, :), allocatable :: PTilde
  real, dimension(:, :, :, :, :), allocatable :: thetaTilde
  real, dimension(:, :, :, :, :), allocatable :: nAerTilde
  real, dimension(:, :, :, :, :), allocatable :: nIceTilde
  real, dimension(:, :, :, :, :), allocatable :: qIceTilde
  real, dimension(:, :, :, :, :), allocatable :: qvTilde
  real, dimension(:, :, :, :, :), allocatable :: tracerTilde

  ! public variables
  ! needed for
  ! 1) BC correction
  ! 2) explicit boundary setting
  ! 3) update module
  public :: rhoTilde, rhopTilde, thetaTilde, rhoTilde_mom, PTilde
  public :: uTilde, vTilde, wTilde
  public :: rhoOld, rhopOld
  public :: nIceTilde, qIceTilde, qvTilde, nAerTilde
  public :: tracerTilde

  ! TFC FJ
  ! Needed for semi-implicit time scheme in TFC.
  public :: uOldTFC, vOldTFC, wOldTFC

  ! SK Needed for compressible explicit Euler
  public :: pinew

  ! phiTilde(i,j,k,dir,Left/Right) with
  ! dir = 1,2,3 for x,y and z reconstruction directions
  ! 0 = Left and 1 = Right for left and right cell edge

  contains

  real function slopeFunction(xx)
    !
    ! calculates slope of "witch of agnesi" mountain
    !
    real, intent(in) :: xx ! x is reserved for grid

    ! local vars
    real :: l, h, s

    ! init vars

    ! scaled quantities
    h = mountainHeight_dim / lRef
    l = mountainWidth_dim / lRef
    s = xx / l

    if(abs(s) < 5.0) then
      slopeFunction = - 2.0 * h / l * s / (1 + s ** 2) ** 2
    else
      slopeFunction = 0.0
    end if

  end function slopeFunction

  !-----------------------------------------------------------------------

  subroutine reconstruction(var, variable)
    !--------------------------------------------------
    ! reconstructs "variable" with
    ! SALD, ALDM, constant, MUSCL
    ! reconstructed variables: \rho/P, \rho'/P, \rho\vec{v}/P
    !-------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var

    character(len = *), intent(in) :: variable
    integer :: i, j, k
    integer :: ix, jy, kz

    ! locals
    integer :: dir, lr

    ! TFC FJ
    real :: rhoEdgeR, rhoEdgeF, rhoEdgeU
    real :: pEdgeR, pEdgeF, pEdgeU

    ! SD
    integer :: ii, iVar

    select case(reconstType)

      !----------------------------
      !   Constant reconstruction
      !----------------------------

    case('constant2')

      stop 'constant2 reconstruction disabled'

      !       select case( variable )
      !
      !       case( "rho" )
      !
      !          do dir = 1,3
      !             do lr = 0,1
      !                rhoTilde(:,:,:,dir,lr) = var(:,:,:,1)
      !             end do
      !          end do
      !
      !       case( "rhop" )
      !
      !          do dir = 1,3
      !             do lr = 0,1
      !                rhopTilde(:,:,:,dir,lr) = var(:,:,:,6)
      !             end do
      !          end do
      !
      !       case( "uvw" )
      !
      !          do dir = 1,3
      !             do lr = 0,1
      !                uTilde(:,:,:,dir,lr) = var(:,:,:,2)
      !                vTilde(:,:,:,dir,lr) = var(:,:,:,3)
      !                wTilde(:,:,:,dir,lr) = var(:,:,:,4)
      !             end do
      !          end do
      !
      !
      !       case( "theta" )
      !
      !          do dir = 1,3
      !             do lr = 0,1
      !                thetaTilde(:,:,:,dir,lr) = var(:,:,:,6)
      !             end do
      !          end do
      !
      !       case default
      !          stop "reconstruction: unknown case model."
      !       end select

      !-----------------------
      !   MUSCL reconstrcution
      !-----------------------

      ! GBcorr: unified 2 reconstruction options for MUSCL
      ! muscl1: cheap but less accurate (???)
      ! muscl2: accurate (???) but expensive

    case('MUSCL')

      select case(musclType)

        ! muscl1: cheap but less accurate (???)
      case("muscl1")

        select case(variable)

        case("rho")

          !rhoBar(:,:,:) = var(:,:,:,1)
          !call reconstruct_MUSCL(rhoBar,rhoTilde,nxx,nyy,nzz,limiterType1)

          ! GBcorr: to be consistent with current rouine: momentumFlux
          rhoBar = 0.0
          if(topography) then
            ! TFC FJ
            ! Adjust reconstruction for 3D fields.
            do ix = - nbx, nx + nbx
              do jy = - nby, ny + nby
                do kz = 0, nz + 1
                  if(pStratTFC(ix, jy, kz) == 0.0) then
                    print *, "ERROR in rec. rho: pStratTFC = 0 at k = ", kz
                    stop
                  end if
                  rhoBar(ix, jy, kz) = var%rho(ix, jy, kz) / pStratTFC(ix, jy, &
                      &kz)
                end do
              end do
            end do
          else
            do kz = 0, nz + 1
              if(Pstrat(kz) == 0.0) then
                print *, 'ERROR in rec. rho: Pstrat(kz) = 0 at kz =', kz
                stop
              end if
              rhoBar(:, :, kz) = (var%rho(:, :, kz)) / Pstrat(kz)
            end do
          end if
          call reconstruct_MUSCL(rhoBar, rhoTilde, nxx, nyy, nzz, limiterType1)

        case("rhop")

          !rhopBar(:,:,:) = var(:,:,:,6)
          !call reconstruct_MUSCL(rhopBar,rhopTilde,nxx,nyy,nzz,&
          !                     & limiterType1)

          ! GBcorr: to be consistent with current rouine: momentumFlux
          rhopBar = 0.0
          if(topography) then
            ! TFC FJ
            ! Adjust reconstruction for 3D fields.
            do ix = - nbx, nx + nbx
              do jy = - nby, ny + nby
                do kz = 0, nz + 1
                  if(pStratTFC(ix, jy, kz) == 0.0) then
                    print *, "ERROR in rec. rho: pStratTFC = 0 at k = ", kz
                    stop
                  end if
                  rhopBar(ix, jy, kz) = var%rhop(ix, jy, kz) / pStratTFC(ix, &
                      &jy, kz)
                end do
              end do
            end do
          else
            do kz = 0, nz + 1
              if(Pstrat(kz) == 0.0) then
                print *, 'ERROR in rec. rhop: Pstrat(kz) = 0 at kz =', kz
                stop
              end if
              rhopBar(:, :, kz) = (var%rhop(:, :, kz)) / Pstrat(kz)
            end do
          end if
          call reconstruct_MUSCL(rhopBar, rhopTilde, nxx, nyy, nzz, &
              &limiterType1)

        case("uvw")

          ! calculate specific momenta to be reconstructed
          if(fluctuationMode) then
            if(topography) then
              ! TFC FJ
              ! Adjust reconstruction for 3D fields.
              do ix = - nbx, nx + nbx - 1
                do jy = - nby, ny + nby
                  do kz = 0, nz + 1
                    rhoEdgeR = 0.5 * (var%rho(ix, jy, kz) + var%rho(ix + 1, &
                        &jy, kz) + rhoStratTFC(ix, jy, kz) + rhoStratTFC(ix &
                        &+ 1, jy, kz))
                    pEdgeR = 0.5 * (pStratTFC(ix, jy, kz) + pStratTFC(ix + 1, &
                        &jy, kz))
                    uBar(ix, jy, kz) = var%u(ix, jy, kz) * rhoEdgeR / pEdgeR
                  end do
                end do
              end do
            else
              do kz = 0, nz + 1
                do ix = - nbx, nx + nbx - 1
                  uBar(ix, :, kz) = var%u(ix, :, kz) * (0.5 * (var%rho(ix, :, &
                      &kz) + var%rho(ix + 1, :, kz)) + rhoStrat(kz)) &
                      &/ Pstrat(kz)
                  ! GBcorr: to be consistent with current rouine: momentumFlux
                  ! * (  0.5*(var(ix,:,kz,1) + var(ix+1,:,kz,1)) &
                  ! + rhoStrat(kz))
                end do
              end do
            end if
          else
            do kz = 0, nz + 1
              do ix = - nbx, nx + nbx - 1
                uBar(ix, :, kz) = var%u(ix, :, kz) * 0.5 * (var%rho(ix, :, kz) &
                    &+ var%rho(ix + 1, :, kz)) / Pstrat(kz)
                ! GBcorr: to be consistent with current rouine: momentumFlux
                ! * 0.5*(var(ix,:,kz,1) + var(ix+1,:,kz,1))
              end do
            end do
          end if

          if(fluctuationMode) then
            if(topography) then
              ! TFC FJ
              ! Adjust reconstruction for 3D fields.
              do ix = - nbx, nx + nbx
                do jy = - nby, ny + nby - 1
                  do kz = 0, nz + 1
                    rhoEdgeF = 0.5 * (var%rho(ix, jy, kz) + var%rho(ix, jy &
                        &+ 1, kz) + rhoStratTFC(ix, jy, kz) + rhoStratTFC(ix, &
                        &jy + 1, kz))
                    pEdgeF = 0.5 * (pStratTFC(ix, jy, kz) + pStratTFC(ix, jy &
                        &+ 1, kz))
                    vBar(ix, jy, kz) = var%v(ix, jy, kz) * rhoEdgeF / pEdgeF
                  end do
                end do
              end do
            else
              do kz = 0, nz + 1
                do jy = - nby, ny + nby - 1
                  vBar(:, jy, kz) = var%v(:, jy, kz) * (0.5 * (var%rho(:, jy, &
                      &kz) + var%rho(:, jy + 1, kz)) + rhoStrat(kz)) &
                      &/ Pstrat(kz)
                  ! GBcorr: to be consistent with current rouine: momentumFlux
                  ! * (  0.5*(var(:,jy,kz,1) + var(:,jy+1,kz,1)) &
                  ! + rhoStrat(kz))
                end do
              end do
            end if
          else
            do kz = 0, nz + 1
              do jy = - nby, ny + nby - 1
                vBar(:, jy, kz) = var%v(:, jy, kz) * 0.5 * (var%rho(:, jy, kz) &
                    &+ var%rho(:, jy + 1, kz)) / Pstrat(kz)
                ! GBcorr: to be consistent with current rouine: momentumFlux
                ! * 0.5*(var(:,jy,kz,1) + var(:,jy+1,kz,1))
              end do
            end do
          end if

          if(fluctuationMode) then
            if(topography) then
              ! TFC FJ
              ! Reconstruct Cartesian vertical momentum.
              wBar(:, :, 0:(nz + 1)) = var%w(:, :, 0:(nz + 1))
              do ix = 1, nx
                do jy = 1, ny
                  do kz = 0, nz + 1
                    wBar(ix, jy, kz) = vertWindTFC(ix, jy, kz, var)
                  end do
                end do
              end do
              call setHalosOfField(wBar)
              do ix = - nbx, nx + nbx
                do jy = - nby, ny + nby
                  do kz = 0, nz + 1
                    rhoEdgeU = 0.5 * (var%rho(ix, jy, kz) + var%rho(ix, jy, kz &
                        &+ 1) + rhoStratTFC(ix, jy, kz) + rhoStratTFC(ix, jy, &
                        &kz + 1))
                    pEdgeU = 0.5 * (pStratTFC(ix, jy, kz) + pStratTFC(ix, jy, &
                        &kz + 1))
                    wBar(ix, jy, kz) = wBar(ix, jy, kz) * rhoEdgeU / pEdgeU
                  end do
                end do
              end do
            else
              do kz = 0, nz + 1
                wBar(:, :, kz) = var%w(:, :, kz) * (0.5 * (var%rho(:, :, kz) &
                    &+ var%rho(:, :, kz + 1)) + rhoStratTilde(kz)) &
                    &/ PstratTilde(kz)
                ! GBcorr: to be consistent with current rouine: momentumFlux
                ! * (  0.5*(var(:,:,kz,1) + var(:,:,kz+1,1)) &
                ! + rhoStratTilde(kz))
                ! FS:
                ! * (  0.5*((var(:,:,kz,1)+rhoStrat(kz))/PStrat(kz) &
                ! + (var(:,:,kz+1,1)+rhoStrat(kz+1))/PStrat(kz+1)))
              end do
            end if
          else
            do kz = 0, nz + 1
              wBar(:, :, kz) = var%w(:, :, kz) * 0.5 * (var%rho(:, :, kz) &
                  &+ var%rho(:, :, kz + 1)) / PstratTilde(kz)
              ! GBcorr: to be consistent with current rouine: momentumFlux
              ! * 0.5*(var(:,:,kz,1) + var(:,:,kz+1,1))
              ! FS:
              ! * (  0.5*((var(:,:,kz,1))/PStrat(kz) &
              ! + (var(:,:,kz+1,1))/PStrat(kz+1)))
            end do
          end if

          ! reconstruct spcific momenta
          call reconstruct_MUSCL(uBar, uTilde, nxx, nyy, nzz, limiterType1)
          call reconstruct_MUSCL(vBar, vTilde, nxx, nyy, nzz, limiterType1)
          call reconstruct_MUSCL(wBar, wTilde, nxx, nyy, nzz, limiterType1)

        case("theta")

          thetaBar(:, :, :) = var%rhop(:, :, :)
          call reconstruct_MUSCL(thetaBar, thetaTilde, nxx, nyy, nzz, &
              &limiterType1)

        case("P") ! reconstruct 1 = P/P
          PBar(:, :, :) = 1
          call reconstruct_MUSCL(PBar, PTilde, nxx, nyy, nzz, limiterType1)

        case("ice")

          if(include_ice) then

            nAerBar(:, :, :) = var%ICE(:, :, :, 1)
            nIceBar(:, :, :) = var%ICE(:, :, :, 2)
            qIceBar(:, :, :) = var%ICE(:, :, :, 3)
            qvBar(:, :, :) = var%ICE(:, :, :, 4)

            call reconstruct_MUSCL(nIceBar, nIceTilde, nxx, nyy, nzz, &
                &limiterType1)
            call reconstruct_MUSCL(qIceBar, qIceTilde, nxx, nyy, nzz, &
                &limiterType1)
            call reconstruct_MUSCL(qvBar, qvTilde, nxx, nyy, nzz, limiterType1)

          end if

        case("ice2")

          do iVar = 1, nVarIce

            iceBar = 0.0
            if(topography) then
              ! TFC FJ
              ! Adjust reconstruction for 3D fields.
              do ix = - nbx, nx + nbx
                do jy = - nby, ny + nby
                  do kz = 0, nz + 1
                    if(pStratTFC(ix, jy, kz) == 0.0) then
                      print *, "ERROR in rec. rho: pStratTFC = 0 at k = ", kz
                      stop
                    end if
                    iceBar(ix, jy, kz) = var%ICE2(ix, jy, kz, iVar) &
                        &/ pStratTFC(ix, jy, kz)
                  end do
                end do
              end do
            else
              do kz = 0, nz + 1
                if(Pstrat(kz) == 0.0) then
                  print *, 'ERROR in rec. rho: Pstrat(kz) = 0 at kz =', kz
                  stop
                end if
                iceBar(:, :, kz) = (var%ICE2(:, :, kz, iVar)) / Pstrat(kz)
              end do
            end if

            if(iVar .eq. inN) then

              call reconstruct_MUSCL(iceBar, nIceTilde, nxx, nyy, nzz, &
                  &limiterType1)

            elseif(iVar .eq. inQ) then

              call reconstruct_MUSCL(iceBar, qIceTilde, nxx, nyy, nzz, &
                  &limiterType1)

            elseif(iVar .eq. inQv) then

              call reconstruct_MUSCL(iceBar, qvTilde, nxx, nyy, nzz, &
                  &limiterType1)

            end if
          end do !ii

        case("tracer")
          ! calucate tracerBar

          tracerBar = 0.0
          if(topography) then
            do ix = - nbx, nx + nbx
              do jy = - nby, ny + nby
                do kz = 0, nz + 1
                  if(pStratTFC(ix, jy, kz) == 0.0) then
                    print *, "ERROR in rec. rho: pStratTFC = 0 at k =", kz
                    stop
                  end if
                  tracerBar(ix, jy, kz) = var%chi(ix, jy, kz) / pStratTFC(ix, &
                      &jy, kz)
                end do
              end do
            end do
          else
            do kz = 0, nz + 1
              if(Pstrat(kz) == 0.0) then
                print *, "Error in rec. rho: Pstrat(kz) = 0 at k = ", kz
                stop
              end if

              tracerBar(:, :, kz) = var%chi(:, :, kz) / Pstrat(kz)
            end do
          end if

          call reconstruct_MUSCL(tracerBar, tracerTilde, nxx, nyy, nzz, &
              &limiterType1)

        case default
          stop "reconstruction: unknown case variable."
        end select

        ! muscl2: accurate (???) but expensive
      case("muscl2")

        select case(variable)

        case("rho")

          !UAB
          rhoBar = 0.0
          do kz = 0, nz + 1
            if(Pstrat(kz) == 0.0) then
              print *, 'ERROR in rec. rho: Pstrat(kz) = 0 at kz =', kz
              stop
            end if
            rhoBar(:, :, kz) = (var%rho(:, :, kz)) / Pstrat(kz)
          end do
          call reconstruct_MUSCL(rhoBar, rhoTilde, nxx, nyy, nzz, limiterType1)
          !UAE

        case("rhop")

          !UAB
          rhopBar = 0.0
          do kz = 0, nz + 1
            if(Pstrat(kz) == 0.0) then
              print *, 'ERROR in rec. rhop: Pstrat(kz) = 0 at kz =', kz
              stop
            end if
            rhopBar(:, :, kz) = var%rhop(:, :, kz) / Pstrat(kz)
          end do
          call reconstruct_MUSCL(rhopBar, rhopTilde, nxx, nyy, nzz, &
              &limiterType1)
          !UAE

        case("uvw")

          !UAB
          ! reconstruct specific momenta \rho \vec{v}/P
          ! for this \rho/P and \vec v are reconstructed each and then the
          ! product of the two is taken

          ! reconstruct momentum in x direction

          rhoBar = 0.0

          if(fluctuationMode) then
            do kz = 0, nz + 1
              if(Pstrat(kz) == 0.0) then
                print *, 'ERROR in rec. u: Pstrat(kz) = 0 at kz =', kz
                stop
              end if
              do ix = - nbx, nx + nbx - 1
                rhoBar(ix, :, kz) = (0.5 * (var%rho(ix, :, kz) + var%rho(ix &
                    &+ 1, :, kz)) + rhoStrat(kz)) / Pstrat(kz)
                !FS:
                ! = 0.5*&
                !   ((var(ix,:,kz,1)+rhoStrat(kz))/Pstrat(kz) + &
                !   (var(ix+1,:,kz,1)+rhoStrat(kz))/Pstrat(kz))
              end do
            end do
          else
            do kz = 0, nz + 1
              if(Pstrat(kz) == 0.0) then
                print *, 'ERROR in rec. u: Pstrat(kz) = 0 at kz =', kz
                stop
              end if
              do ix = - nbx, nx + nbx - 1
                rhoBar(ix, :, kz) = 0.5 * (var%rho(ix, :, kz) + var%rho(ix &
                    &+ 1, :, kz)) / Pstrat(kz)
                !FS:
                != 0.5*&
                !  ((var(ix,:,kz,1))/Pstrat(kz) + &
                !  (var(ix+1,:,kz,1))/Pstrat(kz))
              end do
            end do
          end if

          call reconstruct_MUSCL(rhoBar, rhoTilde_mom, nxx, nyy, nzz, &
              &limiterType1)

          uBar = 0.0

          do kz = 0, nz + 1
            do ix = - nbx, nx + nbx - 1
              uBar(ix, :, kz) = var%u(ix, :, kz)
            end do
          end do

          call reconstruct_MUSCL(uBar, uTilde, nxx, nyy, nzz, limiterType1)

          uTilde = uTilde * rhoTilde_mom

          ! reconstruct momentum in y direction

          rhoBar = 0.0

          if(fluctuationMode) then
            do kz = 0, nz + 1
              if(Pstrat(kz) == 0.0) then
                print *, 'ERROR in rec. v: Pstrat(kz) = 0 at kz =', kz
                stop
              end if
              do jy = - nby, ny + nby - 1
                rhoBar(:, jy, kz) = (0.5 * (var%rho(:, jy, kz) + var%rho(:, jy &
                    &+ 1, kz)) + rhoStrat(kz)) / Pstrat(kz)
                !FS:
                != 0.5*&
                !  ((var(:,jy,kz,1)+rhoStrat(kz))/Pstrat(kz) + &
                !  (var(:,jy+1,kz,1)+rhoStrat(kz))/Pstrat(kz))
              end do
            end do
          else
            do kz = 0, nz + 1
              if(Pstrat(kz) == 0.0) then
                print *, 'ERROR in rec. v: Pstrat(kz) = 0 at kz =', kz
                stop
              end if
              do jy = - nby, ny + nby - 1
                rhoBar(:, jy, kz) = 0.5 * (var%rho(:, jy, kz) + var%rho(:, jy &
                    &+ 1, kz)) / Pstrat(kz)
                !FS:
                ! = 0.5*&
                !   ((var(:,jy,kz,1))/Pstrat(kz) + &
                !   (var(:,jy+1,kz,1))/Pstrat(kz))
              end do
            end do
          end if

          call reconstruct_MUSCL(rhoBar, rhoTilde_mom, nxx, nyy, nzz, &
              &limiterType1)

          vBar = 0.0

          do kz = 0, nz + 1
            do jy = - nby, ny + nby - 1
              vBar(:, jy, kz) = var%v(:, jy, kz)
            end do
          end do

          call reconstruct_MUSCL(vBar, vTilde, nxx, nyy, nzz, limiterType1)

          vTilde = vTilde * rhoTilde_mom

          ! reconstruct momentum in z direction

          rhoBar = 0.0

          if(fluctuationMode) then
            do kz = 0, nz + 1
              if(Pstrat(kz) == 0.0) then
                print *, 'ERROR in rec. w: Pstrat(kz) = 0 at kz =', kz
                stop
              end if
              rhoBar(:, :, kz) = (0.5 * (var%rho(:, :, kz) + var%rho(:, :, kz &
                  &+ 1)) + rhoStratTilde(kz)) / PstratTilde(kz)
              !FS:
              !   = 0.5*&
              !     ((var(:,:,kz,1)+rhoStrat(kz))/Pstrat(kz) + &
              !     (var(:,:,kz+1,1)+rhoStrat(kz+1))/Pstrat(kz+1))
            end do
          else
            do kz = 0, nz + 1
              if(Pstrat(kz) == 0.0) then
                print *, 'ERROR in rec. w: Pstrat(kz) = 0 at kz =', kz
                stop
              end if
              rhoBar(:, :, kz) = 0.5 * (var%rho(:, :, kz) + var%rho(:, :, kz &
                  &+ 1)) / PstratTilde(kz)
              !FS:
              !   = 0.5*&
              !     ((var(:,:,kz,1))/Pstrat(kz) + &
              !     (var(:,:,kz+1,1))/Pstrat(kz+1))
            end do
          end if

          call reconstruct_MUSCL(rhoBar, rhoTilde_mom, nxx, nyy, nzz, &
              &limiterType1)

          wBar = 0.0

          do kz = 0, nz + 1
            wBar(:, :, kz) = var%w(:, :, kz)
          end do

          call reconstruct_MUSCL(wBar, wTilde, nxx, nyy, nzz, limiterType1)

          wTilde = wTilde * rhoTilde_mom
          !UAE

        case("theta")

          thetaBar(:, :, :) = var%rhop(:, :, :)
          call reconstruct_MUSCL(thetaBar, thetaTilde, nxx, nyy, nzz, &
              &limiterType1)

        case("ice")

          if(include_ice) then

            nAerBar(:, :, :) = var%ICE(:, :, :, 1)
            nIceBar(:, :, :) = var%ICE(:, :, :, 2)
            qIceBar(:, :, :) = var%ICE(:, :, :, 3)
            qvBar(:, :, :) = var%ICE(:, :, :, 4)

            call reconstruct_MUSCL(nIceBar, nIceTilde, nxx, nyy, nzz, &
                &limiterType1)
            call reconstruct_MUSCL(qIceBar, qIceTilde, nxx, nyy, nzz, &
                &limiterType1)
            call reconstruct_MUSCL(qvBar, qvTilde, nxx, nyy, nzz, limiterType1)

          end if

        case("ice2")

          print *, 'reconstruction: muscl2 does not work with ice2'

        case default
          stop "reconstruction: unknown case variable."
        end select

      case default
        stop "reconstruction: unknown case musclType."
      end select

      !---------------------------
      !     no reconstruction
      !---------------------------

    case('constant')

      return

      !---------------------------
      !    SALD reconstruction
      !---------------------------

    case('SALD') ! simplified ALDM using reconstruct_SALD

      select case(variable)

      case("rho")

        rhoBar(:, :, :) = var%rho(:, :, :)
        call reconstruct_SALD(rhoBar, rhoTilde)

      case("rhop")

        rhopBar(:, :, :) = var%rhop(:, :, :)
        call reconstruct_SALD(rhopBar, rhopTilde)

      case("uvw")

        uBar(:, :, :) = var%u(:, :, :)
        vBar(:, :, :) = var%v(:, :, :)
        wBar(:, :, :) = var%w(:, :, :)

        call reconstruct_SALD(uBar, uTilde)
        call reconstruct_SALD(vBar, vTilde)
        call reconstruct_SALD(wBar, wTilde)

      case("theta")

        thetaBar(:, :, :) = var%rhop(:, :, :)
        call reconstruct_SALD(thetaBar, thetaTilde)

      case("ice")

        if(include_ice) then

          nAerBar(:, :, :) = var%ICE(:, :, :, 1)
          nIceBar(:, :, :) = var%ICE(:, :, :, 2)
          qIceBar(:, :, :) = var%ICE(:, :, :, 3)
          qvBar(:, :, :) = var%ICE(:, :, :, 4)

          call reconstruct_SALD(nAerBar, nAerTilde)
          call reconstruct_SALD(nIceBar, nIceTilde)
          call reconstruct_SALD(qIceBar, qIceTilde)
          call reconstruct_SALD(qvBar, qvTilde)

        end if

      case("tracer")

        if(include_tracer) then

          tracerBar(:, :, :) = var%chi(:, :, :)

          call reconstruct_SALD(tracerBar, tracerTilde)

        end if

      case default
        stop "reconstruction: unknown case variable."
      end select

      !---------------------------
      !    ALDM reconstruction
      !---------------------------

    case('ALDM') ! full 3D reconstruction according to ALDM

      select case(variable)

      case("rho")

        rhoBar(:, :, :) = var%rho(:, :, :)
        call reconstruct_ALDM(rhoBar, rhoTilde)

      case("rhop")

        rhopBar(:, :, :) = var%rhop(:, :, :)
        call reconstruct_ALDM(rhopBar, rhopTilde)

      case("uvw")

        uBar(:, :, :) = var%u(:, :, :)
        vBar(:, :, :) = var%v(:, :, :)
        wBar(:, :, :) = var%w(:, :, :)

        call reconstruct_ALDM(uBar, uTilde)
        call reconstruct_ALDM(vBar, vTilde)
        call reconstruct_ALDM(wBar, wTilde)

      case("theta")

        thetaBar(:, :, :) = var%rhop(:, :, :)
        call reconstruct_ALDM(thetaBar, thetaTilde)

      case("ice")

        if(include_ice) then

          nAerBar(:, :, :) = var%ICE(:, :, :, 1)
          nIceBar(:, :, :) = var%ICE(:, :, :, 2)
          qIceBar(:, :, :) = var%ICE(:, :, :, 3)
          qvBar(:, :, :) = var%ICE(:, :, :, 4)

          call reconstruct_ALDM(nAerBar, nAerTilde)
          call reconstruct_ALDM(nIceBar, nIceTilde)
          call reconstruct_ALDM(qIceBar, qIceTilde)
          call reconstruct_ALDM(qvBar, qvTilde)

        end if

      case("tracer")

        if(include_tracer) then

          tracerBar(:, :, :) = var%chi(:, :, :)

          call reconstruct_ALDM(tracerBar, tracerTilde)

        end if

      case default
        stop "reconstruction: unknown case variable."
      end select

    case default
      print *, "reconstruction: unknown reconstruction type."
    end select

  end subroutine reconstruction

  !--------------------------------------------------------------------

  subroutine thetaSource(var, source)
    !---------------------------------------------------------------------
    ! computes theta*div(u) for reconstructed u, which do not satisfy
    ! the divergence constraint -> source term corrects flux difference
    !---------------------------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var

    type(var_type), intent(inout) :: source

    integer :: i, j, k
    real :: uL, uR ! L=Left i-1/2, R=Right i+1/2
    real :: vB, vF ! B=Backward j-1/2, F=Forward j+1/2
    real :: wD, wU ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurfL, vSurfF, wSurfU ! velocities at cell surface
    real :: uSurfR, vSurfB, wSurfD ! velocities at cell surface
    real :: div, theta

    !--------------------------------------------------------
    !                  Divergence correction
    !--------------------------------------------------------
    !   u*grad(theta) = div(theta*u) - theta*div(u)

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

          select case(fluxType)

          case("central")
            return

          case("upwind", "ILES")

            uL = uTilde(i, j, k, 1, 0)
            uR = uTilde(i, j, k, 1, 1)
            uSurfR = 0.5 * (uL + uR)

            uL = uTilde(i - 1, j, k, 1, 0)
            uR = uTilde(i - 1, j, k, 1, 1)
            uSurfL = 0.5 * (uL + uR)

            vB = vTilde(i, j, k, 2, 0)
            vF = vTilde(i, j, k, 2, 1)
            vSurfF = 0.5 * (vB + vF)

            vB = vTilde(i, j - 1, k, 2, 0)
            vF = vTilde(i, j - 1, k, 2, 1)
            vSurfB = 0.5 * (vB + vF)

            wD = wTilde(i, j, k, 3, 0)
            wU = wTilde(i, j, k, 3, 1)
            wSurfU = 0.5 * (wD + wU)

            wD = wTilde(i, j, k - 1, 3, 0)
            wU = wTilde(i, j, k - 1, 3, 1)
            wSurfD = 0.5 * (wD + wU)

            div = (uSurfR - uSurfL) / dx + (vSurfF - vSurfB) / dy + (wSurfU &
                &- wSurfD) / dz

            theta = var%rhop(i, j, k)

            source%rhop(i, j, k) = theta * div

          case default
            stop "thetaFlux: unknown case fluxType"
          end select

        end do
      end do
    end do

  end subroutine thetaSource

  !---------------------------------------------------------------------------

  subroutine thetaFlux(var, flux)
    !---------------------------------------------------------------------
    ! computes the theta flux at all cell edges using reconstructed values
    !---------------------------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var
    type(flux_type), intent(inout) :: flux

    integer :: i, j, k, l
    real :: thetaL, thetaR, uL, uR ! L=Left i-1/2, R=Right i+1/2
    real :: thetaB, thetaF, vB, vF ! B=Backward j-1/2, F=Forward j+1/2
    real :: thetaD, thetaU, wD, wU ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurf, vSurf, wSurf ! velocities at cell surface

    real :: fTheta, gTheta, hTheta

    ! avoid abs() for linerisation
    real :: delta
    real, parameter :: delta0 = 1.0e-6

    !--------------------------------------------------------
    !                          Advection
    !--------------------------------------------------------

    !-----------------------------------------
    !       Zonal theta fluxes in x: f
    !-----------------------------------------

    do k = 1, nz
      do j = 1, ny
        do i = 0, nx

          select case(fluxType)

          case("central")
            thetaL = var%rhop(i, j, k)
            thetaR = var%rhop(i + 1, j, k)
            uSurf = var%u(i, j, k)
            fTheta = uSurf * 0.5 * (thetaL + thetaR)

          case("upwind")
            thetaR = thetaTilde(i + 1, j, k, 1, 0)
            thetaL = thetaTilde(i, j, k, 1, 1)
            uL = uTilde(i, j, k, 1, 0)
            uR = uTilde(i, j, k, 1, 1)
            uSurf = 0.5 * (uL + uR)
            fTheta = flux_muscl(uSurf, thetaL, thetaR)

          case("ILES")
            thetaR = thetaTilde(i + 1, j, k, 1, 0)
            thetaL = thetaTilde(i, j, k, 1, 1)
            uL = uTilde(i, j, k, 1, 0)
            uR = uTilde(i, j, k, 1, 1)
            uSurf = 0.5 * (uL + uR)
            fTheta = flux_aldm(thetaL, thetaR, uSurf, thetaL, thetaR, uL, uR, &
                &sigmaC)
          case default
            stop "thetaFlux: unknown case fluxType"
          end select

          flux%rhop(i, j, k, 1) = fTheta
        end do
      end do
    end do

    !-----------------------------------------
    !    Meridional theta fluxes in y: g
    !-----------------------------------------

    do k = 1, nz
      do j = 0, ny
        do i = 1, nx

          select case(fluxType)

          case("central")
            thetaF = var%rhop(i, j + 1, k)
            thetaB = var%rhop(i, j, k)
            vSurf = var%v(i, j, k)
            gTheta = vSurf * 0.5 * (thetaB + thetaF)

          case("upwind")
            thetaF = thetaTilde(i, j + 1, k, 2, 0)
            thetaB = thetaTilde(i, j, k, 2, 1)
            vB = vTilde(i, j, k, 2, 0)
            vF = vTilde(i, j, k, 2, 1)
            vSurf = 0.5 * (vB + vF)
            gTheta = flux_muscl(vSurf, thetaB, thetaF)

          case("ILES")
            thetaF = thetaTilde(i, j + 1, k, 2, 0)
            thetaB = thetaTilde(i, j, k, 2, 1)
            vB = vTilde(i, j, k, 2, 0)
            vF = vTilde(i, j, k, 2, 1)
            vSurf = 0.5 * (vB + vF)
            gTheta = flux_aldm(thetaB, thetaF, vSurf, thetaB, thetaF, vB, vF, &
                &sigmaC)

          case default
            stop "thetaFlux: unknown case fluxType"
          end select

          flux%rhop(i, j, k, 2) = gTheta

        end do
      end do
    end do

    !-----------------------------------------
    !      Vertical theta fluxes in z: h
    !-----------------------------------------

    do k = 0, nz
      do j = 1, ny
        do i = 1, nx

          select case(fluxType)

          case("central")
            thetaU = var%rhop(i, j, k + 1)
            thetaD = var%rhop(i, j, k)
            wSurf = var%w(i, j, k)
            hTheta = wSurf * 0.5 * (thetaD + thetaU)

          case("upwind")
            thetaU = thetaTilde(i, j, k + 1, 3, 0)
            thetaD = thetaTilde(i, j, k, 3, 1)
            wD = wTilde(i, j, k, 3, 0)
            wU = wTilde(i, j, k, 3, 1)
            wSurf = 0.5 * (wD + wU)
            hTheta = flux_muscl(wSurf, thetaD, thetaU)

          case("ILES")
            thetaU = thetaTilde(i, j, k + 1, 3, 0)
            thetaD = thetaTilde(i, j, k, 3, 1)
            wD = wTilde(i, j, k, 3, 0)
            wU = wTilde(i, j, k, 3, 1)
            wSurf = 0.5 * (wD + wU)
            hTheta = flux_aldm(thetaD, thetaU, wSurf, thetaD, thetaU, wU, wD, &
                &sigmaC)

          case default
            stop "thetaFlux: unknown case fluxType"
          end select

          flux%rhop(i, j, k, 3) = hTheta

        end do
      end do
    end do

    !--------------------------------------------------------------
    !                      Heat conduction
    !--------------------------------------------------------------

    if(mu_conduct > 0.0) then

      !-----------------------------------------
      !       Zonal theta fluxes in x: f
      !-----------------------------------------

      do k = 1, nz
        do j = 1, ny
          do i = 0, nx

            thetaL = var%rhop(i, j, k)
            thetaR = var%rhop(i + 1, j, k)
            fTheta = mu_conduct * (thetaR - thetaL) / dx

            flux%rhop(i, j, k, 1) = flux%rhop(i, j, k, 1) - fTheta
          end do
        end do
      end do

      !-----------------------------------------
      !    Meridional theta fluxes in y: g
      !-----------------------------------------

      do k = 1, nz
        do j = 0, ny
          do i = 1, nx

            thetaF = var%rhop(i, j + 1, k)
            thetaB = var%rhop(i, j, k)
            gTheta = mu_conduct * (thetaF - thetaB) / dy

            flux%rhop(i, j, k, 2) = flux%rhop(i, j, k, 2) - gTheta
          end do
        end do
      end do

      !-----------------------------------------
      !      Vertical theta fluxes in z: h
      !-----------------------------------------

      do k = 0, nz
        do j = 1, ny
          do i = 1, nx

            thetaU = var%rhop(i, j, k + 1)
            thetaD = var%rhop(i, j, k)
            hTheta = mu_conduct * (thetaU - thetaD) / dz

            flux%rhop(i, j, k, 3) = flux%rhop(i, j, k, 3) - hTheta
          end do
        end do
      end do

    end if ! mu_conduct > 0.0

    if(verbose) print *, "thetaFlux: theta fluxes fTheta, gTheta and fTheta &
        &calculated"

  end subroutine thetaFlux

  !---------------------------------------------------------------------------

  subroutine massSource(var, source)
    !---------------------------------------------------------------------
    ! computes source term for the conti equation
    ! 1) divergence term rho*div(u)
    !---------------------------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var

    type(var_type), intent(inout) :: source

    integer :: i, j, k
    real :: uL, uR ! L=Left i-1/2, R=Right i+1/2
    real :: vB, vF ! B=Backward j-1/2, F=Forward j+1/2
    real :: wD, wU ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurfL, vSurfF, wSurfU ! velocities at cell surface
    real :: uSurfR, vSurfB, wSurfD ! velocities at cell surface
    real :: PstratU, PstratD
    real :: divPu

    !--------------------------------------------------------
    !                  Divergence correction
    !--------------------------------------------------------
    !   u*grad(rho) = div(rho*u) - rho*div(u)

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx

          select case(fluxType)

          case("central")
            return

          case("upwind", "ILES")

            uR = var%u(i, j, k); uL = var%u(i - 1, j, k)
            vF = var%v(i, j, k); vB = var%v(i, j - 1, k)
            wU = var%w(i, j, k); wD = var%w(i, j, k - 1)

            PstratU = PstratTilde(k)
            PstratD = PstratTilde(k - 1)

            divPu = Pstrat(k) * ((uR - uL) / dx + (vF - vB) / dy) + (PstratU &
                &* wU - PstratD * wD) / dz

            !!$                uSurfR = var(i,j,k,2)
            !!$                uSurfL = var(i-1,j,k,2)
            !!$
            !!$                vSurfF = var(i,j,k,3)
            !!$                vSurfB = var(i,j-1,k,3)
            !!$
            !!$                wSurfU = var(i,j,k,4)
            !!$                wSurfD = var(i,j,k-1,4)
            !!$
            !!$                div = (uSurfR - uSurfL)/dx &
            !!$                     & + (vSurfF - vSurfB) / dy &
            !!$                     & + (wSurfU - wSurfD) / dz

            source%rho(i, j, k) = divPu / thetaStrat(k)

          case default
            stop "rhoFlux: unknown case fluxType"
          end select

        end do
      end do
    end do

  end subroutine massSource

  !-----------------------------------------------------------------------

  subroutine massFlux(vara, var, flux, fluxmode, Pstrata, PStratTildea)
    !---------------------------------------------------------------------
    ! computes the mass flux at all cell edges using reconstructed values
    ! fluxmode = lin => linear flux, advecting velocities prescribed in
    !                   vara
    !            nln => nonlinear flux, advecting velocities from var
    !
    ! MUSCL assumes that the reconstructed densities are \rho/P
    !---------------------------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: vara, var
    character(len = *), intent(in) :: fluxmode

    type(flux_type), intent(inout) :: flux
    real, dimension(- 1:nz + 2), intent(in) :: PStrata, PStratTildea

    integer :: i, j, k, l
    real :: rhoL, rhoR, uL, uR ! L=Left i-1/2, R=Right i+1/2
    real :: rhoB, rhoF, vB, vF ! B=Backward j-1/2, F=Forward j+1/2
    real :: rhoD, rhoU, wD, wU ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurf, vSurf, wSurf ! velocities at cell surface

    real :: pL, pR, pB, pF, pD, pU
    real :: metD, metU

    ! TFC FJ
    real :: rhoStratEdgeR, rhoStratEdgeF, rhoStratEdgeU
    real :: pEdgeR, pEdgeF, pEdgeU

    real :: fRho, gRho, hRho

    ! avoid abs() for linerisation
    real :: delta
    real, parameter :: delta0 = 1.0e-6

    !   variables for the turbulence scheme
    real :: Pr_turb
    real :: coef_t, drho_dxi, dtht_dxi

    real :: delta_hs, delta_vs

    Pr_turb = 1.0 !0.5   !FS

    ! squared grid scales for the anisotropic turbulence scheme

    if(TurbScheme) then
      if(ny == 1 .and. nx == 1) then
        stop 'ERROR: turbulence scheme assumes either nx > 1 or ny > 1'
      else
        if(nx == 1) then
          delta_hs = dy ** 2 ! 2D problems in y and z
        else if(ny == 1) then
          delta_hs = dx ** 2 ! 2D problems in x and z
        else
          delta_hs = dx * dy ! 3D problems

          if(dx / dy > 10.) then
            print *, 'WARNING: dx/dy > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          elseif(dy / dx > 10.) then
            print *, 'WARNING: dy/dx > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          end if
        end if

        delta_vs = dz ** 2
      end if
    end if

    !-----------------------------------------
    !       Zonal rho fluxes in x: f
    !-----------------------------------------

    do k = 1, nz
      do j = 1, ny
        do i = 0, nx
          select case(fluxType)

          case("central")

            if(fluctuationMode) then
              rhoL = var%rho(i, j, k) + rhoStrat(k)
              rhoR = var%rho(i + 1, j, k) + rhoStrat(k)
            else
              rhoL = var%rho(i, j, k)
              rhoR = var%rho(i + 1, j, k)
            end if

            if(fluxmode == "nln") then
              uSurf = var%u(i, j, k)
            else if(fluxmode == "lin") then
              uSurf = vara%u(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            fRho = uSurf * 0.5 * (rhoL + rhoR)

          case("upwind")

            if(fluctuationMode) then
              if(topography) then
                ! TFC FJ
                ! Adjust for 3D fields.
                rhoStratEdgeR = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i &
                    &+ 1, j, k))
                pEdgeR = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i + 1, j, k))
                rhoR = rhoTilde(i + 1, j, k, 1, 0) + rhoStratEdgeR / pEdgeR
                rhoL = rhoTilde(i, j, k, 1, 1) + rhoStratEdgeR / pEdgeR
              else
                !UAB reference density to be divided by P as well!
                rhoR = rhoTilde(i + 1, j, k, 1, 0) + rhoStrat(k) / Pstrat(k)
                rhoL = rhoTilde(i, j, k, 1, 1) + rhoStrat(k) / Pstrat(k)
                !UAE
              end if
            else
              rhoR = rhoTilde(i + 1, j, k, 1, 0)
              rhoL = rhoTilde(i, j, k, 1, 1)
            end if

            if(topography) then
              ! TFC FJ
              pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, &
                  &j, k) * pStratTFC(i + 1, j, k))
              if(fluxmode == "nln") then
                uSurf = pEdgeR * var%u(i, j, k)
              else if(fluxmode == "lin") then
                uSurf = pEdgeR * vara%u(i, j, k)
              else
                stop "ERROR: wrong fluxmode"
              end if
            else
              if(fluxmode == "nln") then
                uSurf = var%u(i, j, k) * Pstrat(k) !UA
              else if(fluxmode == "lin") then
                uSurf = vara%u(i, j, k) * Pstrata(k) !UA
              else
                stop 'ERROR: worng fluxmode'
              end if
            end if

            fRho = flux_muscl(uSurf, rhoL, rhoR)

          case("ILES")
            rhoR = rhoTilde(i + 1, j, k, 1, 0)
            rhoL = rhoTilde(i, j, k, 1, 1)
            uL = uTilde(i, j, k, 1, 0)
            uR = uTilde(i, j, k, 1, 1)

            if(fluxmode == "nln") then
              uSurf = var%u(i, j, k)
            else if(fluxmode == "lin") then
              uSurf = vara%u(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            if(fluctuationMode) then
              fRho = flux_aldm(rhoL + rhoStrat(k), rhoR + rhoStrat(k), uSurf, &
                  &rhoL, rhoR, uL, uR, sigmaC)
            else
              fRho = flux_aldm(rhoL, rhoR, uSurf, rhoL, rhoR, uL, uR, sigmaC)
            end if

          case default
            stop "rhoFlux: unknown case fluxType"
          end select

          flux%rho(i, j, k, 1) = fRho
        end do
      end do
    end do

    !-----------------------------------------
    !    Meridional rho fluxes in y: g
    !-----------------------------------------

    do k = 1, nz
      do j = 0, ny
        do i = 1, nx
          select case(fluxType)

          case("central")
            if(fluctuationMode) then
              rhoF = var%rho(i, j + 1, k) + rhoStrat(k)
              rhoB = var%rho(i, j, k) + rhoStrat(k)
            else
              rhoF = var%rho(i, j + 1, k)
              rhoB = var%rho(i, j, k)
            end if

            if(fluxmode == "nln") then
              vSurf = var%v(i, j, k)
            else if(fluxmode == "lin") then
              vSurf = vara%v(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            gRho = vSurf * 0.5 * (rhoB + rhoF)

          case("upwind")
            if(fluctuationMode) then
              if(topography) then
                ! TFC FJ
                ! Adjust for 3D fields.
                rhoStratEdgeF = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, j &
                    &+ 1, k))
                pEdgeF = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j + 1, k))
                rhoF = rhoTilde(i, j + 1, k, 2, 0) + rhoStratEdgeF / pEdgeF
                rhoB = rhoTilde(i, j, k, 2, 1) + rhoStratEdgeF / pEdgeF
              else
                !UAB reference density to be divided by P as well!
                rhoF = rhoTilde(i, j + 1, k, 2, 0) + rhoStrat(k) / Pstrat(k)
                rhoB = rhoTilde(i, j, k, 2, 1) + rhoStrat(k) / Pstrat(k)
                !UAE
              end if
            else
              rhoF = rhoTilde(i, j + 1, k, 2, 0)
              rhoB = rhoTilde(i, j, k, 2, 1)
            end if

            if(topography) then
              ! TFC FJ
              pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                  &+ 1, k) * pStratTFC(i, j + 1, k))
              if(fluxmode == "nln") then
                vSurf = pEdgeF * var%v(i, j, k)
              else if(fluxmode == "lin") then
                vSurf = pEdgeF * vara%v(i, j, k)
              else
                stop "ERROR: wrong fluxmode"
              end if
            else
              if(fluxmode == "nln") then
                vSurf = var%v(i, j, k) * Pstrat(k) !UA
              else if(fluxmode == "lin") then
                vSurf = vara%v(i, j, k) * Pstrata(k) !UA
              else
                stop 'ERROR: worng fluxmode'
              end if
            end if

            gRho = flux_muscl(vSurf, rhoB, rhoF)

          case("ILES")
            rhoF = rhoTilde(i, j + 1, k, 2, 0)
            rhoB = rhoTilde(i, j, k, 2, 1)
            vB = vTilde(i, j, k, 2, 0)
            vF = vTilde(i, j, k, 2, 1)

            if(fluxmode == "nln") then
              vSurf = var%v(i, j, k)
            else if(fluxmode == "lin") then
              vSurf = vara%v(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            if(fluctuationMode) then
              gRho = flux_aldm(rhoB + rhoStrat(k), rhoF + rhoStrat(k), vSurf, &
                  &rhoB, rhoF, vB, vF, sigmaC)
            else
              gRho = flux_aldm(rhoB, rhoF, vSurf, rhoB, rhoF, vB, vF, sigmaC)
            end if

          case default
            stop "rhoFlux: unknown case fluxType"
          end select

          flux%rho(i, j, k, 2) = gRho
        end do
      end do
    end do

    !-----------------------------------------
    !      Vertical rho fluxes in z: h
    !-----------------------------------------

    do k = 0, nz
      do j = 1, ny
        do i = 1, nx

          select case(fluxType)

          case("central")
            if(fluctuationMode) then
              rhoU = var%rho(i, j, k + 1) + rhoStratTilde(k)
              ! background rho at half level
              rhoD = var%rho(i, j, k) + rhoStratTilde(k)
            else
              rhoU = var%rho(i, j, k + 1)
              rhoD = var%rho(i, j, k)
            end if

            if(fluxmode == "nln") then
              wSurf = var%w(i, j, k)
            else if(fluxmode == "lin") then
              wSurf = vara%w(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            hRho = wSurf * 0.5 * (rhoD + rhoU)

          case("upwind")

            if(fluctuationMode) then
              if(topography) then
                ! TFC FJ
                ! Adjust for 3D fields.
                rhoStratEdgeU = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, &
                    &j, k + 1))
                pEdgeU = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j, k + 1))
                rhoU = rhoTilde(i, j, k + 1, 3, 0) + rhoStratEdgeU / pEdgeU
                rhoD = rhoTilde(i, j, k, 3, 1) + rhoStratEdgeU / pEdgeU
              else
                ! background at half level
                !UAB reference density to be divided by P as well!
                rhoU = rhoTilde(i, j, k + 1, 3, 0) + rhoStratTilde(k) &
                    &/ PstratTilde(k)
                rhoD = rhoTilde(i, j, k, 3, 1) + rhoStratTilde(k) &
                    &/ PstratTilde(k)
                !UAE
              end if
            else
              rhoU = rhoTilde(i, j, k + 1, 3, 0)
              rhoD = rhoTilde(i, j, k, 3, 1)
            end if

            if(topography) then
              ! TFC FJ
              pEdgeU = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, k &
                  &+ 1) * pStratTFC(i, j, k + 1))
              if(fluxmode == "nln") then
                wSurf = pEdgeU * var%w(i, j, k)
              else if(fluxmode == "lin") then
                wSurf = pEdgeU * vara%w(i, j, k)
              else
                stop "ERROR: wrong fluxmode"
              end if
            else
              if(fluxmode == "nln") then
                wSurf = var%w(i, j, k) * PstratTilde(k) !UA
              else if(fluxmode == "lin") then
                wSurf = vara%w(i, j, k) * PstratTildea(k) !UA
              else
                stop 'ERROR: worng fluxmode'
              end if
            end if

            hRho = flux_muscl(wSurf, rhoD, rhoU)

          case("ILES")
            rhoU = rhoTilde(i, j, k + 1, 3, 0)
            rhoD = rhoTilde(i, j, k, 3, 1)
            wD = wTilde(i, j, k, 3, 0)
            wU = wTilde(i, j, k, 3, 1)

            if(fluxmode == "nln") then
              wSurf = var%w(i, j, k)
            else if(fluxmode == "lin") then
              wSurf = vara%w(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            if(fluctuationMode) then
              hRho = flux_aldm(rhoD + rhoStratTilde(k), rhoU &
                  &+ rhoStratTilde(k), wSurf, rhoD, rhoU, wU, wD, sigmaC)
            else
              hRho = flux_aldm(rhoD, rhoU, wSurf, rhoD, rhoU, wU, wD, sigmaC)
            end if

          case default
            stop "rhoFlux: unknown case fluxType"
          end select

          flux%rho(i, j, k, 3) = hRho
        end do
      end do
    end do

    ! --------------------------------------------
    !      fluxes density fluctuations
    ! --------------------------------------------

    if(timeScheme == "semiimplicit" .or. auxil_equ) then
      ! Zonal rhop fluxes in x: f

      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            select case(fluxType)

            case("central")
              rhoL = var%rhop(i, j, k)
              rhoR = var%rhop(i + 1, j, k)

              if(fluxmode == "nln") then
                uSurf = var%u(i, j, k)
              else if(fluxmode == "lin") then
                uSurf = vara%u(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              fRho = uSurf * 0.5 * (rhoL + rhoR)

            case("upwind")
              rhoR = rhopTilde(i + 1, j, k, 1, 0)
              rhoL = rhopTilde(i, j, k, 1, 1)

              if(topography) then
                ! TFC FJ
                pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, &
                    &j, k) * pStratTFC(i + 1, j, k))
                if(fluxmode == "nln") then
                  uSurf = pEdgeR * var%u(i, j, k)
                else if(fluxmode == "lin") then
                  uSurf = pEdgeR * vara%u(i, j, k)
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                if(fluxmode == "nln") then
                  uSurf = var%u(i, j, k) * Pstrat(k) !UA
                else if(fluxmode == "lin") then
                  uSurf = vara%u(i, j, k) * Pstrata(k) !UA
                else
                  stop 'ERROR: worng fluxmode'
                end if
              end if

              fRho = flux_muscl(uSurf, rhoL, rhoR)

            case("ILES")
              rhoR = rhopTilde(i + 1, j, k, 1, 0)
              rhoL = rhopTilde(i, j, k, 1, 1)
              uL = uTilde(i, j, k, 1, 0)
              uR = uTilde(i, j, k, 1, 1)

              if(fluxmode == "nln") then
                uSurf = var%u(i, j, k)
              else if(fluxmode == "lin") then
                uSurf = vara%u(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              fRho = flux_aldm(rhoL, rhoR, uSurf, rhoL, rhoR, uL, uR, sigmaC)

            case default
              stop "rhopFlux: unknown case fluxType"
            end select

            flux%rhop(i, j, k, 1) = fRho
          end do
        end do
      end do

      ! Meridional rhop fluxes in y: g

      do k = 1, nz
        do j = 0, ny
          do i = 1, nx
            select case(fluxType)

            case("central")
              rhoF = var%rhop(i, j + 1, k)
              rhoB = var%rhop(i, j, k)

              if(fluxmode == "nln") then
                vSurf = var%v(i, j, k)
              else if(fluxmode == "lin") then
                vSurf = vara%v(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              gRho = vSurf * 0.5 * (rhoB + rhoF)

            case("upwind")
              rhoF = rhopTilde(i, j + 1, k, 2, 0)
              rhoB = rhopTilde(i, j, k, 2, 1)

              if(topography) then
                ! TFC FJ
                pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                    &+ 1, k) * pStratTFC(i, j + 1, k))
                if(fluxmode == "nln") then
                  vSurf = pEdgeF * var%v(i, j, k)
                else if(fluxmode == "lin") then
                  vSurf = pEdgeF * vara%v(i, j, k)
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                if(fluxmode == "nln") then
                  vSurf = var%v(i, j, k) * Pstrat(k) !UA
                else if(fluxmode == "lin") then
                  vSurf = vara%v(i, j, k) * Pstrata(k) !UA
                else
                  stop 'ERROR: worng fluxmode'
                end if
              end if

              gRho = flux_muscl(vSurf, rhoB, rhoF)

            case("ILES")
              rhoF = rhopTilde(i, j + 1, k, 2, 0)
              rhoB = rhopTilde(i, j, k, 2, 1)
              vB = vTilde(i, j, k, 2, 0)
              vF = vTilde(i, j, k, 2, 1)

              if(fluxmode == "nln") then
                vSurf = var%v(i, j, k)
              else if(fluxmode == "lin") then
                vSurf = vara%v(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              gRho = flux_aldm(rhoB, rhoF, vSurf, rhoB, rhoF, vB, vF, sigmaC)

            case default
              stop "rhoFlux: unknown case fluxType"
            end select

            flux%rhop(i, j, k, 2) = gRho
          end do
        end do
      end do

      ! Vertical rhop fluxes in z: h

      do k = 0, nz
        do j = 1, ny
          do i = 1, nx

            select case(fluxType)

            case("central")
              rhoU = var%rhop(i, j, k + 1)
              rhoD = var%rhop(i, j, k)

              if(fluxmode == "nln") then
                wSurf = var%w(i, j, k)
              else if(fluxmode == "lin") then
                wSurf = vara%w(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              hRho = wSurf * 0.5 * (rhoD + rhoU)

            case("upwind")
              rhoU = rhopTilde(i, j, k + 1, 3, 0)
              rhoD = rhopTilde(i, j, k, 3, 1)

              if(topography) then
                ! TFC FJ
                pEdgeU = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, &
                    &k + 1) * pStratTFC(i, j, k + 1))
                if(fluxmode == "nln") then
                  wSurf = pEdgeU * var%w(i, j, k)
                else if(fluxmode == "lin") then
                  wSurf = pEdgeU * vara%w(i, j, k)
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                if(fluxmode == "nln") then
                  wSurf = var%w(i, j, k) * PstratTilde(k) !UA
                else if(fluxmode == "lin") then
                  wSurf = vara%w(i, j, k) * PstratTildea(k) !UA
                else
                  stop 'ERROR: worng fluxmode'
                end if
              end if

              hRho = flux_muscl(wSurf, rhoD, rhoU)

            case("ILES")
              rhoU = rhopTilde(i, j, k + 1, 3, 0)
              rhoD = rhopTilde(i, j, k, 3, 1)
              wD = wTilde(i, j, k, 3, 0)
              wU = wTilde(i, j, k, 3, 1)

              if(fluxmode == "nln") then
                wSurf = var%w(i, j, k)
              else if(fluxmode == "lin") then
                wSurf = vara%w(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              hRho = flux_aldm(rhoD, rhoU, wSurf, rhoD, rhoU, wU, wD, sigmaC)

            case default
              stop "rhoFlux: unknown case fluxType"
            end select

            flux%rhop(i, j, k, 3) = hRho
          end do
        end do
      end do
    end if

    ! -------------------------------------------
    ! fluxes mass-weighted potential temperature
    !--------------------------------------------
    if(model == "compressible") then
      ! Zonal fluxes in x: f

      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            select case(fluxType)

            case("upwind")
              rhoR = 1.0
              rhoL = 1.0

              if(topography) then
                ! TFC FJ
                pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, &
                    &j, k) * pStratTFC(i + 1, j, k))
                if(fluxmode == "nln") then
                  uSurf = pEdgeR * var%u(i, j, k)
                else if(fluxmode == "lin") then
                  uSurf = pEdgeR * vara%u(i, j, k)
                else
                  stop "ERROR: wrong fluxmode"
                end if
              end if

              fRho = flux_muscl(uSurf, rhoL, rhoR)

            case default
              stop "PFlux: unknown case fluxType"
            end select

            flux%P(i, j, k, 1) = fRho
          end do
        end do
      end do

      ! Meridional rhop fluxes in y: g

      do k = 1, nz
        do j = 0, ny
          do i = 1, nx
            select case(fluxType)

            case("upwind")
              rhoF = 1.0
              rhoB = 1.0

              if(topography) then
                ! TFC FJ
                pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                    &+ 1, k) * pStratTFC(i, j + 1, k))
                if(fluxmode == "nln") then
                  vSurf = pEdgeF * var%v(i, j, k)
                else if(fluxmode == "lin") then
                  vSurf = pEdgeF * vara%v(i, j, k)
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                if(fluxmode == "nln") then
                  vSurf = var%v(i, j, k) * Pstrat(k) !UA
                else if(fluxmode == "lin") then
                  vSurf = vara%v(i, j, k) * Pstrata(k) !UA
                else
                  stop 'ERROR: worng fluxmode'
                end if
              end if

              gRho = flux_muscl(vSurf, rhoB, rhoF)

            case default
              stop "PFlux: unknown case fluxType"
            end select

            flux%P(i, j, k, 2) = gRho
          end do
        end do
      end do

      ! Vertical rhop fluxes in z: h

      do k = 0, nz
        do j = 1, ny
          do i = 1, nx

            select case(fluxType)

            case("upwind")
              rhoU = 1.0
              rhoD = 1.0

              if(topography) then
                ! TFC FJ
                pEdgeU = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, &
                    &k + 1) * pStratTFC(i, j, k + 1))
                if(fluxmode == "nln") then
                  wSurf = pEdgeU * var%w(i, j, k)
                else if(fluxmode == "lin") then
                  wSurf = pEdgeU * vara%w(i, j, k)
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                if(fluxmode == "nln") then
                  wSurf = var%w(i, j, k) * PstratTilde(k) !UA
                else if(fluxmode == "lin") then
                  wSurf = vara%w(i, j, k) * PstratTildea(k) !UA
                else
                  stop 'ERROR: worng fluxmode'
                end if
              end if

              hRho = flux_muscl(wSurf, rhoD, rhoU)

            case default
              stop "rhoFlux: unknown case fluxType"
            end select

            flux%P(i, j, k, 3) = hRho
          end do
        end do
      end do
    end if

    !--------------------------------------------------------
    !  contributions from molecular and turbulent diffusion
    !  to the potential-temperature fluxes
    !  --> stored in flux(:,:,:,:,5)
    !--------------------------------------------------------

    if(mu_conduct == 0.0 .and. .not. TurbScheme) return

    ! flux in x direction

    if(topography) then
      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            select case(model)
            case("pseudo_incompressible", "compressible")
              coef_t = mu_conduct * 0.5 * (rhoStratTFC(i, j, 1) &
                  &/ rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, j, 1) &
                  &/ rhoStratTFC(i + 1, j, k))
            case("Boussinesq")
              coef_t = mu_conduct
            case default
              stop "diffusivity: unkown case model."
            end select

            if(TurbScheme) then
              coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i + 1, j, &
                  &k)) * delta_hs / Pr_turb
            end if

            rhoL = var%rho(i, j, k) + rhoStratTFC(i, j, k)
            rhoR = var%rho(i + 1, j, k) + rhoStratTFC(i + 1, j, k)
            rhoD = 0.5 * (var%rho(i, j, k - 1) + rhoStratTFC(i, j, k - 1) &
                &+ var%rho(i + 1, j, k - 1) + rhoStratTFC(i + 1, j, k - 1))
            rhoU = 0.5 * (var%rho(i, j, k + 1) + rhoStratTFC(i, j, k + 1) &
                &+ var%rho(i + 1, j, k + 1) + rhoStratTFC(i + 1, j, k + 1))

            pL = pStratTFC(i, j, k)
            pR = pStratTFC(i + 1, j, k)
            pD = 0.5 * (pStratTFC(i, j, k - 1) + pStratTFC(i + 1, j, k - 1))
            pU = 0.5 * (pStratTFC(i, j, k + 1) + pStratTFC(i + 1, j, k + 1))

            dtht_dxi = 0.5 * (jac(i, j, k) + jac(i + 1, j, k)) * (pR / rhoR &
                &- pL / rhoL) / dx + 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) &
                &+ jac(i + 1, j, k) * met(i + 1, j, k, 1, 3)) * (pU / rhoU &
                &- pD / rhoD) / (2.0 * dz)

            flux%theta(i, j, k, 1) = - coef_t * dtht_dxi
          end do
        end do
      end do
    else
      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            ! density dependent diffusivity
            ! turbulence scheme allowing for anisotropic grids

            select case(model)
            case("pseudo_incompressible")
              coef_t = mu_conduct * rhoStrat(1) / rhoStrat(k)
            case("Boussinesq")
              coef_t = mu_conduct
            case default
              stop "diffusivity: unkown case model."
            end select

            if(TurbScheme) then
              coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i + 1, j, &
                  &k)) * delta_hs / Pr_turb
            end if

            if(fluctuationMode) then
              rhoL = var%rho(i, j, k) + rhoStrat(k)
              rhoR = var%rho(i + 1, j, k) + rhoStrat(k)
            else
              rhoL = var%rho(i, j, k)
              rhoR = var%rho(i + 1, j, k)
            end if

            dtht_dxi = (Pstrat(k) / rhoR - Pstrat(k) / rhoL) / dx

            flux%theta(i, j, k, 1) = - coef_t * dtht_dxi

          end do
        end do
      end do
    end if

    ! flux in y direction

    if(topography) then
      do k = 1, nz
        do j = 0, ny
          do i = 1, nx
            select case(model)
            case("pseudo_incompressible", "compressible")
              coef_t = mu_conduct * 0.5 * (rhoStratTFC(i, j, 1) &
                  &/ rhoStratTFC(i, j, k) + rhoStratTFC(i, j + 1, 1) &
                  &/ rhoStratTFC(i, j + 1, k))
            case("Boussinesq")
              coef_t = mu_conduct
            case default
              stop "diffusivity: unkown case model."
            end select

            if(TurbScheme) then
              coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i, j + 1, &
                  &k)) * delta_hs / Pr_turb
            end if

            rhoB = var%rho(i, j, k) + rhoStratTFC(i, j, k)
            rhoF = var%rho(i, j + 1, k) + rhoStratTFC(i, j + 1, k)
            rhoD = 0.5 * (var%rho(i, j, k - 1) + rhoStratTFC(i, j, k - 1) &
                &+ var%rho(i, j + 1, k - 1) + rhoStratTFC(i, j + 1, k - 1))
            rhoU = 0.5 * (var%rho(i, j, k + 1) + rhoStratTFC(i, j, k + 1) &
                &+ var%rho(i, j + 1, k + 1) + rhoStratTFC(i, j + 1, k + 1))

            pB = pStratTFC(i, j, k)
            pF = pStratTFC(i, j + 1, k)
            pD = 0.5 * (pStratTFC(i, j, k - 1) + pStratTFC(i, j + 1, k - 1))
            pU = 0.5 * (pStratTFC(i, j, k + 1) + pStratTFC(i, j + 1, k + 1))

            dtht_dxi = 0.5 * (jac(i, j, k) + jac(i, j + 1, k)) * (pF / rhoF &
                &- pB / rhoB) / dy + 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) &
                &+ jac(i, j + 1, k) * met(i, j + 1, k, 2, 3)) * (pU / rhoU &
                &- pD / rhoD) / (2.0 * dz)

            flux%theta(i, j, k, 2) = - coef_t * dtht_dxi
          end do
        end do
      end do
    else
      do k = 1, nz
        do j = 0, ny
          do i = 1, nx
            ! density dependent diffusivity
            ! turbulence scheme allowing for anisotropic grids

            select case(model)
            case("pseudo_incompressible")
              coef_t = mu_conduct * rhoStrat(1) / rhoStrat(k)
            case("Boussinesq")
              coef_t = mu_conduct
            case default
              stop "diffusivity: unkown case model."
            end select

            if(TurbScheme) then
              coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i, j + 1, &
                  &k)) * delta_hs / Pr_turb !UA
            end if

            if(fluctuationMode) then
              rhoF = var%rho(i, j + 1, k) + rhoStrat(k)
              rhoB = var%rho(i, j, k) + rhoStrat(k)
            else
              rhoF = var%rho(i, j + 1, k)
              rhoB = var%rho(i, j, k)
            end if

            dtht_dxi = (Pstrat(k) / rhoF - Pstrat(k) / rhoB) / dy

            flux%theta(i, j, k, 2) = - coef_t * dtht_dxi

          end do
        end do
      end do
    end if

    ! flux in z direction

    if(topography) then
      do k = 0, nz
        do j = 1, ny
          do i = 1, nx
            select case(model)
            case("pseudo_incompressible", "compressible")
              coef_t = mu_conduct * 0.5 * (rhoStratTFC(i, j, 1) &
                  &/ rhoStratTFC(i, j, k) + rhoStratTFC(i, j, 1) &
                  &/ rhoStratTFC(i, j, k + 1))
            case("Boussinesq")
              coef_t = mu_conduct
            case default
              stop "diffusivity: unkown case model."
            end select

            if(TurbScheme) then
              coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i, j, k &
                  &+ 1)) * jac(i, j, k) ** 2.0 * delta_vs / Pr_turb
            end if

            rhoL = 0.5 * (var%rho(i - 1, j, k) + rhoStratTFC(i - 1, j, k) &
                &+ var%rho(i - 1, j, k + 1) + rhoStratTFC(i - 1, j, k + 1))
            rhoR = 0.5 * (var%rho(i + 1, j, k) + rhoStratTFC(i + 1, j, k) &
                &+ var%rho(i + 1, j, k + 1) + rhoStratTFC(i + 1, j, k + 1))
            rhoB = 0.5 * (var%rho(i, j - 1, k) + rhoStratTFC(i, j - 1, k) &
                &+ var%rho(i, j - 1, k + 1) + rhoStratTFC(i, j - 1, k + 1))
            rhoF = 0.5 * (var%rho(i, j + 1, k) + rhoStratTFC(i, j + 1, k) &
                &+ var%rho(i, j + 1, k + 1) + rhoStratTFC(i, j + 1, k + 1))
            rhoD = var%rho(i, j, k) + rhoStratTFC(i, j, k)
            rhoU = var%rho(i, j, k + 1) + rhoStratTFC(i, j, k + 1)

            pL = 0.5 * (pStratTFC(i - 1, j, k) + pStratTFC(i - 1, j, k + 1))
            pR = 0.5 * (pStratTFC(i + 1, j, k) + pStratTFC(i + 1, j, k + 1))
            pB = 0.5 * (pStratTFC(i, j - 1, k) + pStratTFC(i, j - 1, k + 1))
            pF = 0.5 * (pStratTFC(i, j + 1, k) + pStratTFC(i, j + 1, k + 1))
            pD = pStratTFC(i, j, k)
            pU = pStratTFC(i, j, k + 1)

            dtht_dxi = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) + jac(i, j, k &
                &+ 1) * met(i, j, k + 1, 1, 3)) * (pR / rhoR - pL / rhoL) &
                &/ (2.0 * dx) + 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) &
                &+ jac(i, j, k + 1) * met(i, j, k + 1, 2, 3)) * (pF / rhoF &
                &- pB / rhoB) / (2.0 * dy) + 0.5 * (jac(i, j, k) * met(i, j, &
                &k, 3, 3) + jac(i, j, k + 1) * met(i, j, k + 1, 3, 3)) * (pU &
                &/ rhoU - pD / rhoD) / dz

            flux%theta(i, j, k, 3) = - coef_t * dtht_dxi
          end do
        end do
      end do
    else
      do k = 0, nz
        do j = 1, ny
          do i = 1, nx
            ! density dependent diffusivity
            ! turbulence scheme allowing for anisotropic grids

            select case(model)
            case("pseudo_incompressible")
              coef_t = mu_conduct * rhoStrat(1) / rhoStrat(k)
            case("Boussinesq")
              coef_t = mu_conduct
            case default
              stop "diffusivity: unkown case model."
            end select

            if(TurbScheme) then
              coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i, j, k &
                  &+ 1)) * delta_vs / Pr_turb !UA
            end if

            if(fluctuationMode) then
              rhoU = var%rho(i, j, k + 1) + rhoStrat(k + 1)
              rhoD = var%rho(i, j, k) + rhoStrat(k)
            else
              rhoU = var%rho(i, j, k + 1)
              rhoD = var%rho(i, j, k)
            end if

            dtht_dxi = (Pstrat(k + 1) / rhoU - Pstrat(k) / rhoD) / dz

            flux%theta(i, j, k, 3) = - coef_t * dtht_dxi

          end do
        end do
      end do
    end if

    if(verbose) print *, "rhoFlux: rho fluxes fRho, gRho and fRho calculated"

  end subroutine massFlux

  !-----------------------------------------------------------------------

  subroutine massFlux_0(vara, var, flux, fluxmode)
    !---------------------------------------------------------------------
    ! computes the mass flux at all cell edges using reconstructed values
    ! fluxmode = lin => linear flux, advecting velocities prescribed in
    !                   vara
    !            nln => nonlinear flux, advecting velocities from var
    !
    ! this version assumes that the reconstructed densities are \rho
    !---------------------------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: vara, var
    character(len = *), intent(in) :: fluxmode

    type(flux_type), intent(inout) :: flux

    integer :: i, j, k, l
    real :: rhoL, rhoR, uL, uR ! L=Left i-1/2, R=Right i+1/2
    real :: rhoB, rhoF, vB, vF ! B=Backward j-1/2, F=Forward j+1/2
    real :: rhoD, rhoU, wD, wU ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurf, vSurf, wSurf ! velocities at cell surface

    real :: fRho, gRho, hRho

    ! avoid abs() for linerisation
    real :: delta
    real, parameter :: delta0 = 1.0e-6

    !   variables for the turbulence scheme
    real :: Pr_turb
    real :: coef_t, drho_dxi, dtht_dxi

    real :: delta_hs, delta_vs

    Pr_turb = 0.5

    ! squared grid scales for the anisotropic turbulence scheme

    if(TurbScheme) then
      if(ny == 1 .and. nx == 1) then
        stop 'ERROR: turbulence scheme assumes either nx > 1 or ny > 1'
      else
        if(nx == 1) then
          delta_hs = dy ** 2 ! 2D problems in y and z
        else if(ny == 1) then
          delta_hs = dx ** 2 ! 2D problems in x and z
        else
          delta_hs = dx * dy ! 3D problems

          if(dx / dy > 10.) then
            print *, 'WARNING: dx/dy > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          elseif(dy / dx > 10.) then
            print *, 'WARNING: dy/dx > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          end if
        end if

        delta_vs = dz ** 2
      end if
    end if

    !-----------------------------------------
    !       Zonal rho fluxes in x: f
    !-----------------------------------------

    do k = 1, nz
      do j = 1, ny
        do i = 0, nx
          select case(fluxType)

          case("central")

            if(fluctuationMode) then
              rhoL = var%rho(i, j, k) + rhoStrat(k)
              rhoR = var%rho(i + 1, j, k) + rhoStrat(k)
            else
              rhoL = var%rho(i, j, k)
              rhoR = var%rho(i + 1, j, k)
            end if

            if(fluxmode == "nln") then
              uSurf = var%u(i, j, k)
            else if(fluxmode == "lin") then
              uSurf = vara%u(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            fRho = uSurf * 0.5 * (rhoL + rhoR)

          case("upwind")

            if(fluctuationMode) then
              rhoR = rhoTilde(i + 1, j, k, 1, 0) + rhoStrat(k)
              rhoL = rhoTilde(i, j, k, 1, 1) + rhoStrat(k)
            else
              rhoR = rhoTilde(i + 1, j, k, 1, 0)
              rhoL = rhoTilde(i, j, k, 1, 1)
            end if

            if(fluxmode == "nln") then
              uSurf = var%u(i, j, k)
            else if(fluxmode == "lin") then
              uSurf = vara%u(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            fRho = flux_muscl(uSurf, rhoL, rhoR)

          case("ILES")
            rhoR = rhoTilde(i + 1, j, k, 1, 0)
            rhoL = rhoTilde(i, j, k, 1, 1)
            uL = uTilde(i, j, k, 1, 0)
            uR = uTilde(i, j, k, 1, 1)

            if(fluxmode == "nln") then
              uSurf = var%u(i, j, k)
            else if(fluxmode == "lin") then
              uSurf = vara%u(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            if(fluctuationMode) then
              fRho = flux_aldm(rhoL + rhoStrat(k), rhoR + rhoStrat(k), uSurf, &
                  &rhoL, rhoR, uL, uR, sigmaC)
            else
              fRho = flux_aldm(rhoL, rhoR, uSurf, rhoL, rhoR, uL, uR, sigmaC)
            end if

          case default
            stop "rhoFlux: unknown case fluxType"
          end select

          flux%rho(i, j, k, 1) = fRho
        end do
      end do
    end do

    !-----------------------------------------
    !    Meridional rho fluxes in y: g
    !-----------------------------------------

    do k = 1, nz
      do j = 0, ny
        do i = 1, nx
          select case(fluxType)

          case("central")
            if(fluctuationMode) then
              rhoF = var%rho(i, j + 1, k) + rhoStrat(k)
              rhoB = var%rho(i, j, k) + rhoStrat(k)
            else
              rhoF = var%rho(i, j + 1, k)
              rhoB = var%rho(i, j, k)
            end if

            if(fluxmode == "nln") then
              vSurf = var%v(i, j, k)
            else if(fluxmode == "lin") then
              vSurf = vara%v(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            gRho = vSurf * 0.5 * (rhoB + rhoF)

          case("upwind")
            if(fluctuationMode) then
              rhoF = rhoTilde(i, j + 1, k, 2, 0) + rhoStrat(k)
              rhoB = rhoTilde(i, j, k, 2, 1) + rhoStrat(k)
            else
              rhoF = rhoTilde(i, j + 1, k, 2, 0)
              rhoB = rhoTilde(i, j, k, 2, 1)
            end if

            if(fluxmode == "nln") then
              vSurf = var%v(i, j, k)
            else if(fluxmode == "lin") then
              vSurf = vara%v(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            gRho = flux_muscl(vSurf, rhoB, rhoF)

          case("ILES")
            rhoF = rhoTilde(i, j + 1, k, 2, 0)
            rhoB = rhoTilde(i, j, k, 2, 1)
            vB = vTilde(i, j, k, 2, 0)
            vF = vTilde(i, j, k, 2, 1)

            if(fluxmode == "nln") then
              vSurf = var%v(i, j, k)
            else if(fluxmode == "lin") then
              vSurf = vara%v(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            if(fluctuationMode) then
              gRho = flux_aldm(rhoB + rhoStrat(k), rhoF + rhoStrat(k), vSurf, &
                  &rhoB, rhoF, vB, vF, sigmaC)
            else
              gRho = flux_aldm(rhoB, rhoF, vSurf, rhoB, rhoF, vB, vF, sigmaC)
            end if

          case default
            stop "rhoFlux: unknown case fluxType"
          end select

          flux%rho(i, j, k, 2) = gRho
        end do
      end do
    end do

    !-----------------------------------------
    !      Vertical rho fluxes in z: h
    !-----------------------------------------

    do k = 0, nz
      do j = 1, ny
        do i = 1, nx

          select case(fluxType)

          case("central")
            if(fluctuationMode) then
              rhoU = var%rho(i, j, k + 1) + rhoStratTilde(k)
              ! background rho at half level
              rhoD = var%rho(i, j, k) + rhoStratTilde(k)
            else
              rhoU = var%rho(i, j, k + 1)
              rhoD = var%rho(i, j, k)
            end if

            if(fluxmode == "nln") then
              wSurf = var%w(i, j, k)
            else if(fluxmode == "lin") then
              wSurf = vara%w(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            hRho = wSurf * 0.5 * (rhoD + rhoU)

          case("upwind")

            if(fluctuationMode) then
              ! background at half level
              rhoU = rhoTilde(i, j, k + 1, 3, 0) + rhoStratTilde(k)
              rhoD = rhoTilde(i, j, k, 3, 1) + rhoStratTilde(k)
            else
              rhoU = rhoTilde(i, j, k + 1, 3, 0)
              rhoD = rhoTilde(i, j, k, 3, 1)
            end if

            if(fluxmode == "nln") then
              wSurf = var%w(i, j, k)
            else if(fluxmode == "lin") then
              wSurf = vara%w(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            hRho = flux_muscl(wSurf, rhoD, rhoU)

          case("ILES")
            rhoU = rhoTilde(i, j, k + 1, 3, 0)
            rhoD = rhoTilde(i, j, k, 3, 1)
            wD = wTilde(i, j, k, 3, 0)
            wU = wTilde(i, j, k, 3, 1)

            if(fluxmode == "nln") then
              wSurf = var%w(i, j, k)
            else if(fluxmode == "lin") then
              wSurf = vara%w(i, j, k)
            else
              stop 'ERROR: worng fluxmode'
            end if

            if(fluctuationMode) then
              hRho = flux_aldm(rhoD + rhoStratTilde(k), rhoU &
                  &+ rhoStratTilde(k), wSurf, rhoD, rhoU, wU, wD, sigmaC)
            else
              hRho = flux_aldm(rhoD, rhoU, wSurf, rhoD, rhoU, wU, wD, sigmaC)
            end if

          case default
            stop "rhoFlux: unknown case fluxType"
          end select

          flux%rho(i, j, k, 3) = hRho
        end do
      end do
    end do

    ! --------------------------------------------
    !      fluxes density fluctuations
    ! --------------------------------------------

    if(timeScheme == "semiimplicit" .or. auxil_equ) then
      ! Zonal rhop fluxes in x: f

      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            select case(fluxType)

            case("central")

              rhoL = var%rhop(i, j, k)
              rhoR = var%rhop(i + 1, j, k)

              if(fluxmode == "nln") then
                uSurf = var%u(i, j, k)
              else if(fluxmode == "lin") then
                uSurf = vara%u(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              fRho = uSurf * 0.5 * (rhoL + rhoR)

            case("upwind")

              rhoR = rhopTilde(i + 1, j, k, 1, 0)
              rhoL = rhopTilde(i, j, k, 1, 1)

              if(fluxmode == "nln") then
                uSurf = var%u(i, j, k)
              else if(fluxmode == "lin") then
                uSurf = vara%u(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              fRho = flux_muscl(uSurf, rhoL, rhoR)

            case("ILES")
              rhoR = rhopTilde(i + 1, j, k, 1, 0)
              rhoL = rhopTilde(i, j, k, 1, 1)
              uL = uTilde(i, j, k, 1, 0)
              uR = uTilde(i, j, k, 1, 1)

              if(fluxmode == "nln") then
                uSurf = var%u(i, j, k)
              else if(fluxmode == "lin") then
                uSurf = vara%u(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              fRho = flux_aldm(rhoL, rhoR, uSurf, rhoL, rhoR, uL, uR, sigmaC)

            case default
              stop "rhopFlux: unknown case fluxType"
            end select

            flux%rhop(i, j, k, 1) = fRho
          end do
        end do
      end do

      ! Meridional rhop fluxes in y: g

      do k = 1, nz
        do j = 0, ny
          do i = 1, nx
            select case(fluxType)

            case("central")
              rhoF = var%rhop(i, j + 1, k)
              rhoB = var%rhop(i, j, k)

              if(fluxmode == "nln") then
                vSurf = var%v(i, j, k)
              else if(fluxmode == "lin") then
                vSurf = vara%v(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              gRho = vSurf * 0.5 * (rhoB + rhoF)

            case("upwind")
              rhoF = rhopTilde(i, j + 1, k, 2, 0)
              rhoB = rhopTilde(i, j, k, 2, 1)

              if(fluxmode == "nln") then
                vSurf = var%v(i, j, k)
              else if(fluxmode == "lin") then
                vSurf = vara%v(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              gRho = flux_muscl(vSurf, rhoB, rhoF)

            case("ILES")
              rhoF = rhopTilde(i, j + 1, k, 2, 0)
              rhoB = rhopTilde(i, j, k, 2, 1)
              vB = vTilde(i, j, k, 2, 0)
              vF = vTilde(i, j, k, 2, 1)

              if(fluxmode == "nln") then
                vSurf = var%v(i, j, k)
              else if(fluxmode == "lin") then
                vSurf = vara%v(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              gRho = flux_aldm(rhoB, rhoF, vSurf, rhoB, rhoF, vB, vF, sigmaC)

            case default
              stop "rhoFlux: unknown case fluxType"
            end select

            flux%rhop(i, j, k, 2) = gRho
          end do
        end do
      end do

      ! Vertical rhop fluxes in z: h

      do k = 0, nz
        do j = 1, ny
          do i = 1, nx

            select case(fluxType)

            case("central")
              rhoU = var%rhop(i, j, k + 1)
              rhoD = var%rhop(i, j, k)

              if(fluxmode == "nln") then
                wSurf = var%w(i, j, k)
              else if(fluxmode == "lin") then
                wSurf = vara%w(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              hRho = wSurf * 0.5 * (rhoD + rhoU)

            case("upwind")

              rhoU = rhopTilde(i, j, k + 1, 3, 0)
              rhoD = rhopTilde(i, j, k, 3, 1)

              if(fluxmode == "nln") then
                wSurf = var%w(i, j, k)
              else if(fluxmode == "lin") then
                wSurf = vara%w(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              hRho = flux_muscl(wSurf, rhoD, rhoU)

            case("ILES")
              rhoU = rhopTilde(i, j, k + 1, 3, 0)
              rhoD = rhopTilde(i, j, k, 3, 1)
              wD = wTilde(i, j, k, 3, 0)
              wU = wTilde(i, j, k, 3, 1)

              if(fluxmode == "nln") then
                wSurf = var%w(i, j, k)
              else if(fluxmode == "lin") then
                wSurf = vara%w(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              hRho = flux_aldm(rhoD, rhoU, wSurf, rhoD, rhoU, wU, wD, sigmaC)

            case default
              stop "rhoFlux: unknown case fluxType"
            end select

            flux%rhop(i, j, k, 3) = hRho
          end do
        end do
      end do
    end if

    !--------------------------------------------------------
    !  contributions from molecular and turbulent diffusion
    !  to the potential-temperature fluxes
    !  --> stored in flux(:,:,:,:,5)
    !--------------------------------------------------------

    if(mu_conduct == 0.0 .and. .not. TurbScheme) return

    ! flux in x direction

    do k = 1, nz
      do j = 1, ny
        do i = 0, nx
          ! density dependent diffusivity
          ! turbulence scheme allowing for anisotropic grids

          select case(model)
          case("pseudo_incompressible", "compressible")
            coef_t = mu_conduct * rhoStrat(1) / rhoStrat(k)
          case("Boussinesq")
            coef_t = mu_conduct
          case default
            stop "diffusivity: unkown case model."
          end select

          if(TurbScheme) then
            coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i + 1, j, k)) &
                &* delta_hs / Pr_turb
          end if

          if(fluctuationMode) then
            rhoL = var%rho(i, j, k) + rhoStrat(k)
            rhoR = var%rho(i + 1, j, k) + rhoStrat(k)
          else
            rhoL = var%rho(i, j, k)
            rhoR = var%rho(i + 1, j, k)
          end if

          dtht_dxi = (Pstrat(k) / rhoR - Pstrat(k) / rhoL) / dx

          flux%theta(i, j, k, 1) = - coef_t * dtht_dxi
        end do
      end do
    end do

    ! flux in y direction

    do k = 1, nz
      do j = 0, ny
        do i = 1, nx
          ! density dependent diffusivity
          ! turbulence scheme allowing for anisotropic grids

          select case(model)
          case("pseudo_incompressible", "compressible")
            coef_t = mu_conduct * rhoStrat(1) / rhoStrat(k)
          case("Boussinesq")
            coef_t = mu_conduct
          case default
            stop "diffusivity: unkown case model."
          end select

          if(TurbScheme) then
            coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i, j + 1, k)) &
                &* delta_hs / Pr_turb
          end if

          if(fluctuationMode) then
            rhoF = var%rho(i, j + 1, k) + rhoStrat(k)
            rhoB = var%rho(i, j, k) + rhoStrat(k)
          else
            rhoF = var%rho(i, j + 1, k)
            rhoB = var%rho(i, j, k)
          end if

          dtht_dxi = (Pstrat(k) / rhoF - Pstrat(k) / rhoB) / dy

          flux%theta(i, j, k, 2) = - coef_t * dtht_dxi
        end do
      end do
    end do

    ! flux in z direction

    do k = 0, nz
      do j = 1, ny
        do i = 1, nx
          ! density dependent diffusivity
          ! turbulence scheme allowing for anisotropic grids

          select case(model)
          case("pseudo_incompressible", "compressible")
            coef_t = mu_conduct * rhoStrat(1) / rhoStrat(k)
          case("Boussinesq")
            coef_t = mu_conduct
          case default
            stop "diffusivity: unkown case model."
          end select

          if(TurbScheme) then
            coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i, j, k + 1)) &
                &* delta_vs / Pr_turb
          end if

          if(fluctuationMode) then
            rhoU = var%rho(i, j, k + 1) + rhoStrat(k + 1)
            rhoD = var%rho(i, j, k) + rhoStrat(k)
          else
            rhoU = var%rho(i, j, k + 1)
            rhoD = var%rho(i, j, k)
          end if

          dtht_dxi = (Pstrat(k + 1) / rhoU - Pstrat(k) / rhoD) / dz

          flux%theta(i, j, k, 3) = - coef_t * dtht_dxi
        end do
      end do
    end do

    if(verbose) print *, "rhoFlux: rho fluxes fRho, gRho and fRho calculated"

  end subroutine massFlux_0

  !---------------------------------------------------------------------------
  subroutine tracerFlux(vara, var, flux, fluxmode, Pstrata, PStratTildea)
    ! calculate tracer fluxes
    ! zonal:      u*rho*chi = flux(:, :, :, 1, iVarT)
    ! meridional: v*rho*chi = flux(:, :, :, 2, iVarT)
    ! vertical:   w*rho*chi = flux(:, :, :, 3, iVarT)
    implicit none

    type(var_type), intent(in) :: vara, var
    character(len = *), intent(in) :: fluxmode

    type(flux_type), intent(inout) :: flux

    real, dimension(- 1:nz + 2), intent(in) :: PStrata, PStratTildea

    integer :: i, j, k, l
    real :: tracerL, tracerR, uL, uR ! L=Left i-1/2, R=Right i+1/2
    real :: tracerB, tracerF, vB, vF ! B=Backward j-1/2, F=Forward j+1/2
    real :: tracerD, tracerU, wD, wU ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurf, vSurf, wSurf ! velocities at cell surface
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: rho

    real :: rhoStratEdgeR, rhoStratEdgeF, rhoStratEdgeU
    real :: pEdgeR, pEdgeF, pEdgeU

    real :: fTracer, gTracer, hTracer

    if(TurbScheme) then
      stop "fluxes.f90: TurbScheme not implemented for tracer."
    end if

    !-----------------------------------------
    !       Zonal tracer fluxes in x: f
    !-----------------------------------------
    do k = 1, nz
      do j = 1, ny
        do i = 0, nx
          select case(fluxType)
          case("central")

            tracerL = var%chi(i, j, k)
            tracerR = var%chi(i + 1, j, k)

            if(fluxmode == "nln") then
              uSurf = var%u(i, j, k)
            else if(fluxmode == "lin") then
              uSurf = vara%u(i, j, k)
            else
              stop "fluxes.f90/tracerFlux: wrong fluxmode."
            end if

            fTracer = uSurf * 0.5 * (tracerL + tracerR)

          case("upwind")

            tracerR = tracerTilde(i + 1, j, k, 1, 0)
            tracerL = tracerTilde(i, j, k, 1, 1)

            if(topography) then
              pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, &
                  &j, k) * pStratTFC(i + 1, j, k))
              if(fluxmode == "nln") then
                uSurf = pEdgeR * var%u(i, j, k)
              else if(fluxmode == "lin") then
                uSurf = pEdgeR * vara%u(i, j, k)
              else
                stop "fluxes.f90: wrong fluxmode"
              end if
            else
              if(fluxmode == "nln") then
                uSurf = var%u(i, j, k) * Pstrat(k)
              else if(fluxmode == "lin") then
                uSurf = vara%u(i, j, k) * Pstrata(k)
              else
                stop "fluxes.f90: wrong fluxmode"
              end if
            end if

            fTracer = flux_muscl(uSurf, tracerL, tracerR)

          case default
            stop "fluxes.f90: only fluxType central and upwind implemented for &
                &tracer."
          end select

          flux%chi(i, j, k, 1) = fTracer

        end do
      end do
    end do

    !-----------------------------------
    ! Meridional tracer fluxes in y: g
    !-----------------------------------
    do k = 1, nz
      do j = 0, ny
        do i = 1, nx
          select case(fluxType)

          case("central")

            tracerF = var%chi(i, j + 1, k)
            tracerB = var%chi(i, j, k)

            if(fluxmode == "nln") then
              vSurf = var%v(i, j, k)
            else if(fluxmode == "lin") then
              vSurf = vara%v(i, j, k)
            else
              stop "fluxes.f90/tracerFlux: wrong fluxmode."
            end if

            gTracer = vSurf * 0.5 * (tracerF + tracerB)

          case("upwind")

            tracerF = tracerTilde(i, j + 1, k, 2, 0)
            tracerB = tracerTilde(i, j, k, 2, 1)

            if(topography) then
              pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                  &+ 1, k) * pStratTFC(i, j + 1, k))
              if(fluxmode == "nln") then
                vSurf = pEdgeF * var%v(i, j, k)
              else if(fluxmode == "lin") then
                vSurf = pEdgeF * vara%v(i, j, k)
              else
                stop "fluxes.f90: wrong fluxmode"
              end if
            else
              if(fluxmode == "nln") then
                vSurf = var%v(i, j, k) * Pstrat(k)
              else if(fluxmode == "lin") then
                vSurf = vara%v(i, j, k) * Pstrata(k)
              else
                stop "fluxes.f90: wrong fluxmode"
              end if
            end if

            gTracer = flux_muscl(vSurf, tracerB, tracerF)

          case default
            stop "fluxes.f90: only fluxType central and upwind implemented for &
                &tracer."
          end select
          flux%chi(i, j, k, 2) = gTracer
        end do
      end do
    end do

    !--------------------------------
    ! Vertical tracer fluxes in z: h
    !-------------------------------
    do k = 0, nz
      do j = 1, ny
        do i = 1, nx
          select case(fluxType)

          case("central")

            tracerU = var%chi(i, j, k + 1)
            tracerD = var%chi(i, j, k)

            if(fluxmode == "nln") then
              wSurf = var%w(i, j, k)
            else if(fluxmode == "lin") then
              wSurf = vara%w(i, j, k)
            else
              stop "fluxes.f90/tracerFlux: wrong fluxmode."
            end if

            hTracer = wSurf * 0.5 * (tracerD + tracerU)

          case("upwind")

            tracerU = tracerTilde(i, j, k + 1, 3, 0)
            tracerD = tracerTilde(i, j, k, 3, 1)

            if(topography) then
              pEdgeU = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, k &
                  &+ 1) * pStratTFC(i, j, k + 1))
              if(fluxmode == "nln") then
                wSurf = pEdgeU * var%w(i, j, k)
              else if(fluxmode == "lin") then
                wSurf = pEdgeU * vara%w(i, j, k)
              else
                stop "fluxes.f90: wrong fluxmode"
              end if
            else
              if(fluxmode == "nln") then
                wSurf = var%w(i, j, k) * PstratTilde(k)
              else if(fluxmode == "lin") then
                wSurf = vara%w(i, j, k) * PstratTildea(k)
              else
                stop "fluxes.f90: wrong fluxmode"
              end if
            end if

            hTracer = flux_muscl(wSurf, tracerD, tracerU)

          case default
            stop "fluxes.f90: only fluxType central and upwind implemented for &
                &tracer."
          end select

          flux%chi(i, j, k, 3) = hTracer

        end do
      end do
    end do

  end subroutine tracerFlux

  !---------------------------------------------------------------------------

  subroutine iceFlux(var, flux)
    ! in/out variables
    type(var_type), intent(in) :: var

    type(flux_type), intent(inout) :: flux
    ! flux(i,j,k,dir,iFlux)
    ! dir = 1..3 > f,g,h-flux in x,y,z-direction

    integer :: nqS, k, j, i, iVar

    real :: UpFlux, DownFlux, TotFlux, wU, wD
    real :: delta_w, wSurf, coef_t, d_dxi, m_ice, T, p
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: rho
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, 0:1) &
        &:: rhoTil

    real, parameter :: Pr_t_DySma = 0.5

    ! All ice particles obey the general mass flux in x and y direction

    coef_t = mu_conduct

    if(fluctuationMode) then
      do k = - 1, nz + 1
        rho(:, :, k) = var%rho(:, :, k) + rhoStrat(k)
        rhoTil(:, :, k, :, :) = rhoTilde(:, :, k, :, :) + rhoStrat(k)
      end do
    else
      rho = var%rho(:, :, :)
      rhoTil(:, :, :, :, :) = rhoTilde(:, :, :, :, :)
    end if

    do k = 0, nz
      do j = 0, ny
        do i = 0, nx

          ! find the current temperature in Kelvin inside the grid cell
          call find_temperature(T, i, j, k, var)

          ! find the current pressure in Pascal inside the grid cell
          p = press0_dim * ((PStrat(k) / p0) ** gamma_1 + var%pi(i, j, k)) &
              &** kappaInv

          ! find the current average ice particle mass, if there is no abundant ice
          ! use the initial ice particle mass
          if(var%ICE(i, j, k, 2) == 0.0) then
            m_ice = init_m_ice
          else
            m_ice = var%ICE(i, j, k, 3) / var%ICE(i, j, k, 2) * rhoRef * lRef &
                &** 3
            !print*,"ATTENTION:", m_ice
          end if

          if(model == "pseudo_incompressible") then
            coef_t = mu_conduct * rhoStrat(1) / rhoStrat(k)
          end if

          if(DySmaScheme) then
            coef_t = coef_t + 0.5 * (var%DSC(i, j, k) + var%DSC(i + 1, j, k)) &
                &/ Pr_t_DySma
          end if

          !-----------------------------------------
          !--- horizontal direction ---!
          !-----------------------------------------
          select case(fluxType)

          case("central")
            do iVar = 1, 4
              UpFlux = rho(i + 1, j, k) * var%ICE(i + 1, j, k, iVar)
              DownFLux = rho(i, j, k) * var%ICE(i, j, k, iVar)
              wSurf = var%u(i, j, k)
              TotFlux = wSurf * 0.5 * (DownFlux + UpFlux)
              d_dxi = (UpFlux - DownFlux) / dx
              flux%ICE(i, j, k, 1, iVar) = TotFlux - coef_t * d_dxi
            end do

          case("upwind")
            do iVar = 1, 4
              UpFlux = rho(i + 1, j, k) * var%ICE(i + 1, j, k, iVar)
              DownFLux = rho(i, j, k) * var%ICE(i, j, k, iVar)
              d_dxi = (UpFlux - DownFlux) / dx
              select case(iVar)
              case(4)
                wSurf = var%u(i, j, k)
                UpFlux = rhoTil(i + 1, j, k, 1, 0) * qvTilde(i + 1, j, k, 1, 0)
                DownFlux = rhoTil(i, j, k, 1, 1) * qvTilde(i, j, k, 1, 1)
              case(3)
                wSurf = var%u(i, j, k)
                UpFlux = rhoTil(i + 1, j, k, 1, 0) * qIceTilde(i + 1, j, k, 1, &
                    &0)
                DownFlux = rhoTil(i, j, k, 1, 1) * qIceTilde(i, j, k, 1, 1)
              case(2)
                wSurf = var%u(i, j, k)
                UpFlux = rhoTil(i + 1, j, k, 1, 0) * nIceTilde(i + 1, j, k, 1, &
                    &0)
                DownFlux = rhoTil(i, j, k, 1, 1) * nIceTilde(i, j, k, 1, 1)
              case(1)
                wSurf = var%u(i, j, k)
                UpFlux = rhoTil(i + 1, j, k, 1, 0) * nAerTilde(i + 1, j, k, 1, &
                    &0)
                DownFlux = rhoTil(i, j, k, 1, 1) * nAerTilde(i, j, k, 1, 1)
              case default
                stop "Bug in iceHorizontalFlux"
              end select
              TotFlux = flux_muscl(wSurf, DownFlux, UpFlux)
              flux%ICE(i, j, k, 1, iVar) = TotFlux - coef_t * d_dxi
            end do

          case("ILES")
            do iVar = 1, 4
              UpFlux = rho(i + 1, j, k) * var%ICE(i + 1, j, k, iVar)
              DownFLux = rho(i, j, k) * var%ICE(i, j, k, iVar)
              d_dxi = (UpFlux - DownFlux) / dx
              select case(iVar)
              case(4)
                wD = wTilde(i, j, k, 1, 0)
                wU = wTilde(i, j, k, 1, 1)
                wSurf = var%u(i, j, k)
                UpFlux = rhoTil(i + 1, j, k, 1, 0) * qvTilde(i + 1, j, k, 1, 0)
                DownFlux = rhoTil(i, j, k, 1, 1) * qvTilde(i, j, k, 1, 1)
              case(3)
                wD = wTilde(i, j, k, 1, 0)
                wU = wTilde(i, j, k, 1, 1)
                wSurf = var%u(i, j, k)
                UpFlux = rhoTil(i + 1, j, k, 1, 0) * qIceTilde(i + 1, j, k, 1, &
                    &0)
                DownFlux = rhoTil(i, j, k, 1, 1) * qIceTilde(i, j, k, 1, 1)
              case(2)
                wD = wTilde(i, j, k, 1, 0)
                wU = wTilde(i, j, k, 1, 1)
                wSurf = var%u(i, j, k)
                UpFlux = rhoTil(i + 1, j, k, 1, 0) * nIceTilde(i + 1, j, k, 1, &
                    &0)
                DownFlux = rhoTil(i, j, k, 1, 1) * nIceTilde(i, j, k, 1, 1)
              case(1)
                wD = wTilde(i, j, k, 1, 0)
                wU = wTilde(i, j, k, 1, 1)
                wSurf = var%u(i, j, k)
                UpFlux = rhoTil(i + 1, j, k, 1, 0) * nAerTilde(i + 1, j, k, 1, &
                    &0)
                DownFlux = rhoTil(i, j, k, 1, 1) * nAerTilde(i, j, k, 1, 1)
              case default
                stop "Bug in iceHorizontalFlux"
              end select
              TotFlux = flux_aldm(DownFlux, UpFlux, wSurf, DownFlux, UpFlux, &
                  &wU, wD, sigmaC)
              flux%ICE(i, j, k, 1, iVar) = TotFlux - coef_t * d_dxi
            end do

          case default
            stop "iceSedimentationFlux: unknown case fluxType"

          end select

          !-----------------------------------------
          !--- meridional direction ---!
          !-----------------------------------------
          select case(fluxType)

          case("central")
            do iVar = 1, 4
              UpFlux = rho(i, j + 1, k) * var%ICE(i, j + 1, k, iVar)
              DownFLux = rho(i, j, k) * var%ICE(i, j, k, iVar)
              wSurf = var%v(i, j, k)
              TotFlux = wSurf * 0.5 * (DownFlux + UpFlux)
              d_dxi = (UpFlux - DownFlux) / dy
              flux%ICE(i, j, k, 2, iVar) = TotFlux - coef_t * d_dxi
            end do

          case("upwind")
            do iVar = 1, 4
              UpFlux = rho(i, j + 1, k) * var%ICE(i, j + 1, k, iVar)
              DownFLux = rho(i, j, k) * var%ICE(i, j, k, iVar)
              d_dxi = (UpFlux - DownFlux) / dy
              select case(iVar)
              case(4)
                wSurf = var%v(i, j, k)
                UpFlux = rhoTil(i, j + 1, k, 2, 0) * qvTilde(i, j + 1, k, 2, 0)
                DownFlux = rhoTil(i, j, k, 2, 1) * qvTilde(i, j, k, 2, 1)
              case(3)
                wSurf = var%v(i, j, k)
                UpFlux = rhoTil(i, j + 1, k, 2, 0) * qIceTilde(i, j + 1, k, 2, &
                    &0)
                DownFlux = rhoTil(i, j, k, 2, 1) * qIceTilde(i, j, k, 2, 1)
              case(2)
                wSurf = var%v(i, j, k)
                UpFlux = rhoTil(i, j + 1, k, 2, 0) * nIceTilde(i, j + 1, k, 2, &
                    &0)
                DownFlux = rhoTil(i, j, k, 2, 1) * nIceTilde(i, j, k, 2, 1)
              case(1)
                wSurf = var%v(i, j, k)
                UpFlux = rhoTil(i, j + 1, k, 2, 0) * nAerTilde(i, j + 1, k, 2, &
                    &0)
                DownFlux = rhoTil(i, j, k, 2, 1) * nAerTilde(i, j, k, 2, 1)
              case default
                stop "Bug in iceMeridionalFlux"
              end select
              TotFlux = flux_muscl(wSurf, DownFlux, UpFlux)
              flux%ICE(i, j, k, 2, iVar) = TotFlux - coef_t * d_dxi
            end do

          case("ILES")
            do iVar = 1, 4
              UpFlux = rho(i, j + 1, k) * var%ICE(i, j + 1, k, iVar)
              DownFLux = rho(i, j, k) * var%ICE(i, j, k, iVar)
              d_dxi = (UpFlux - DownFlux) / dy
              select case(iVar)
              case(4)
                wD = wTilde(i, j, k, 2, 0)
                wU = wTilde(i, j, k, 2, 1)
                wSurf = var%v(i, j, k)
                UpFlux = rhoTil(i, j + 1, k, 2, 0) * qvTilde(i, j + 1, k, 2, 0)
                DownFlux = rhoTil(i, j, k, 2, 1) * qvTilde(i, j, k, 2, 1)
              case(3)
                wD = wTilde(i, j, k, 2, 0)
                wU = wTilde(i, j, k, 2, 1)
                wSurf = var%v(i, j, k)
                UpFlux = rhoTil(i, j + 1, k, 2, 0) * qIceTilde(i, j + 1, k, 2, &
                    &0)
                DownFlux = rhoTil(i, j, k, 2, 1) * qIceTilde(i, j, k, 2, 1)
              case(2)
                wD = wTilde(i, j, k, 2, 0)
                wU = wTilde(i, j, k, 2, 1)
                wSurf = var%v(i, j, k)
                UpFlux = rhoTil(i, j + 1, k, 2, 0) * nIceTilde(i, j + 1, k, 2, &
                    &0)
                DownFlux = rhoTil(i, j, k, 2, 1) * nIceTilde(i, j, k, 2, 1)
              case(1)
                wD = wTilde(i, j, k, 2, 0)
                wU = wTilde(i, j, k, 2, 1)
                wSurf = var%v(i, j, k)
                UpFlux = rhoTil(i, j + 1, k, 2, 0) * nAerTilde(i, j + 1, k, 2, &
                    &0)
                DownFlux = rhoTil(i, j, k, 2, 1) * nAerTilde(i, j, k, 2, 1)
              case default
                stop "Bug in iceMeridionalFlux"
              end select
              TotFlux = flux_aldm(DownFlux, UpFlux, wSurf, DownFlux, UpFlux, &
                  &wU, wD, sigmaC)
              flux%ICE(i, j, k, 2, iVar) = TotFlux - coef_t * d_dxi
            end do

          case default
            stop "iceSedimentationFlux: unknown case fluxType"

          end select

          !-----------------------------------------
          !--- vertical direction ---!
          !-----------------------------------------
          ! inspired by massFlux in vertical direction, but here the
          ! vertical wind is superposed with the terminal velocity for ice crystals
          select case(fluxType)

          case("central")
            do iVar = 1, 4
              UpFlux = rho(i, j, k + 1) * var%ICE(i, j, k + 1, iVar)
              DownFLux = rho(i, j, k) * var%ICE(i, j, k, iVar)
              ! to stop ice particles from falling to the ground:
              ! (causes ice piling up near boundary)
              !if (topography_mask(i+is+nbx-1,j+js+nby-1,k)) then
              !  wSurf = 0.0
              !else
              select case(iVar)
              case(3)
                wSurf = var%w(i, j, k) - terminal_v_qIce(m_ice, T, p) / uRef
              case(2)
                wSurf = var%w(i, j, k) - terminal_v_nIce(m_ice, T, p) / uRef
              case default
                wSurf = var%w(i, j, k)
              end select
              !end if
              TotFlux = wSurf * 0.5 * (DownFlux + UpFlux)
              d_dxi = (UpFlux - DownFlux) / dz
              flux%ICE(i, j, k, 3, iVar) = TotFlux - coef_t * d_dxi
            end do

          case("upwind")
            do iVar = 1, 4
              UpFlux = rho(i, j, k + 1) * var%ICE(i, j, k + 1, iVar)
              DownFLux = rho(i, j, k) * var%ICE(i, j, k, iVar)
              d_dxi = (UpFlux - DownFlux) / dz

              select case(iVar)
              case(4)
                wSurf = var%w(i, j, k)
                UpFlux = rhoTil(i, j, k + 1, 3, 0) * qvTilde(i, j, k + 1, 3, 0)
                DownFlux = rhoTil(i, j, k, 3, 1) * qvTilde(i, j, k, 3, 1)
              case(3)
                !if (topography_mask(i+is+nbx-1,j+js+nby-1,k)) then
                ! wSurf = 0.0
                !else
                wSurf = var%w(i, j, k) - terminal_v_qIce(m_ice, T, p) / uRef
                !end if
                UpFlux = rhoTil(i, j, k + 1, 3, 0) * qIceTilde(i, j, k + 1, 3, &
                    &0)
                DownFlux = rhoTil(i, j, k, 3, 1) * qIceTilde(i, j, k, 3, 1)
              case(2)
                !if (topography_mask(i+is+nbx-1,j+js+nby-1,k)) then
                ! wSurf = 0.0
                !else
                wSurf = var%w(i, j, k) - terminal_v_nIce(m_ice, T, p) / uRef
                !end if
                UpFlux = rhoTil(i, j, k + 1, 3, 0) * nIceTilde(i, j, k + 1, 3, &
                    &0)
                DownFlux = rhoTil(i, j, k, 3, 1) * nIceTilde(i, j, k, 3, 1)
              case(1)
                wSurf = var%w(i, j, k)
                UpFlux = rhoTil(i, j, k + 1, 3, 0) * nAerTilde(i, j, k + 1, 3, &
                    &0)
                DownFlux = rhoTil(i, j, k, 3, 1) * nAerTilde(i, j, k, 3, 1)
              case default
                stop "Bug in iceSedimentationFlux"
              end select
              TotFlux = flux_muscl(wSurf, DownFlux, UpFlux)
              flux%ICE(i, j, k, 3, iVar) = TotFlux - coef_t * d_dxi
            end do

          case("ILES")
            do iVar = 1, 4
              UpFlux = rho(i, j, k + 1) * var%ICE(i, j, k + 1, iVar)
              DownFLux = rho(i, j, k) * var%ICE(i, j, k, iVar)
              d_dxi = (UpFlux - DownFlux) / dz
              select case(iVar)
              case(4)
                wD = wTilde(i, j, k, 3, 0)
                wU = wTilde(i, j, k, 3, 1)
                wSurf = var%w(i, j, k)
                UpFlux = rhoTil(i, j, k + 1, 3, 0) * qvTilde(i, j, k + 1, 3, 0)
                DownFlux = rhoTil(i, j, k, 3, 1) * qvTilde(i, j, k, 3, 1)
              case(3)
                !if (topography_mask(i+is+nbx-1,j+js+nby-1,k)) then
                ! delta_w = 0.0
                !else
                delta_w = - terminal_v_qIce(m_ice, T, p) / uRef
                !end if
                wD = wTilde(i, j, k, 3, 0) + delta_w
                wU = wTilde(i, j, k, 3, 1) + delta_w
                wSurf = var%w(i, j, k) + delta_w
                UpFlux = rhoTil(i, j, k + 1, 3, 0) * qIceTilde(i, j, k + 1, 3, &
                    &0)
                DownFlux = rhoTil(i, j, k, 3, 1) * qIceTilde(i, j, k, 3, 1)
              case(2)
                !if (topography_mask(i+is+nbx-1,j+js+nby-1,k)) then
                !  delta_w = 0.0
                !else
                delta_w = - terminal_v_nIce(m_ice, T, p) / uRef
                !end if
                wD = wTilde(i, j, k, 3, 0) + delta_w
                wU = wTilde(i, j, k, 3, 1) + delta_w
                wSurf = var%w(i, j, k) + delta_w
                UpFlux = rhoTil(i, j, k + 1, 3, 0) * nIceTilde(i, j, k + 1, 3, &
                    &0)
                DownFlux = rhoTil(i, j, k, 3, 1) * nIceTilde(i, j, k, 3, 1)
              case(1)
                wD = wTilde(i, j, k, 3, 0)
                wU = wTilde(i, j, k, 3, 1)
                wSurf = var%w(i, j, k)
                UpFlux = rhoTil(i, j, k + 1, 3, 0) * nAerTilde(i, j, k + 1, 3, &
                    &0)
                DownFlux = rhoTil(i, j, k, 3, 1) * nAerTilde(i, j, k, 3, 1)
              case default
                stop "Bug in iceSedimentationFlux"
              end select
              TotFlux = flux_aldm(DownFlux, UpFlux, wSurf, DownFlux, UpFlux, &
                  &wU, wD, sigmaC)
              flux%ICE(i, j, k, 3, iVar) = TotFlux - coef_t * d_dxi
            end do

          case default
            stop "iceSedimentationFlux: unknown case fluxType"

          end select

        end do
      end do
    end do

  end subroutine iceFlux

  ! --------------------------------------------------------------------

  subroutine iceSource(var, source)
    ! in/out variables
    type(var_type), intent(in) :: var

    type(var_type), intent(inout) :: source

    real :: SIce

    integer :: k, j, i

    real :: T ! current temperature in Kelvin
    real :: p ! current pressure in Pascal
    real :: m_ice ! mean ice crystal mass in kg
    real :: nucleation, deposition, evaporation

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ! find the current temperature in Kelvin inside the grid cell
          call find_temperature(T, i, j, k, var)

          ! find the current pressure in Pascal inside the grid cell
          p = press0_dim * ((PStrat(k) / p0) ** gamma_1 + var%pi(i, j, k)) &
              &** kappaInv

          ! find the current saturation with respect to ice inside the grid cell
          SIce = var%ICE(i, j, k, 4) * p / (epsilon0 * pIce(T))
          ! print a warning in case SIce takes unreasonable values
          if((SIce < 0) .or. (SIce > 2)) then
            print *, "#+#+#+#+#+#+#+#+#"
            print *, "SIce=", SIce, ", k = ", k, "i = ", i + is + nbx - 1
            print *, "qv = ", var%ICE(i, j, k, 4)
            print *, "T = ", T
            print *, "#+#+#+#+#+#+#+#+#"
          end if

          ! find the current average ice particle mass, if there is no abundant ice
          ! use the initial ice particle mass
          if(var%ICE(i, j, k, 2) .le. 0.0) then
            m_ice = init_m_ice
          else
            m_ice = abs(var%ICE(i, j, k, 3) / var%ICE(i, j, k, 2) * rhoRef &
                &* lRef ** 3)
          end if

          nucleation = NUCn(i, j, k, var, SIce, T, p, m_ice) ! nucleation of ice crystals by aerosols
          deposition = DEPq(i, j, k, var, SIce, T, p, m_ice) ! depositional growth of ice crystals
          !print*, "nucleation = ", nucleation
          !print*, "deposition = ", deposition

          if(evaporation_on) then
            if(deposition .lt. 0.0) then
              evaporation = min(- 1 / m_ice * rhoRef * lRef ** 3 * deposition, &
                  &var%ICE(i, j, k, 2) / dt_ice * tRef)
              deposition = - m_ice / rhoRef / lRef ** 3 * evaporation
            else
              evaporation = 0.0
            end if
          else
            evaporation = 0.0
            if(deposition .lt. 0.0) then
              deposition = 0.0
            end if
          end if

          ! nAerosol equation
          source%ICE(i, j, k, 1) = - nucleation + evaporation

          ! nIce equation
          source%ICE(i, j, k, 2) = nucleation - evaporation

          ! qIce equation
          source%ICE(i, j, k, 3) = init_m_ice / (rhoRef * lRef ** 3) &
              &* nucleation + deposition

          ! qv equation
          source%ICE(i, j, k, 4) = - init_m_ice / (rhoRef * lRef ** 3) &
              &* nucleation - deposition

        end do
      end do
    end do

  end subroutine iceSource

  !-----------------------------------------------------------------------!

  subroutine volumeForce(var, time, force)

    !------------------------------------------------------------
    ! supplememts volume forces on a grid cell by gravitation and
    ! Coriolis force
    !------------------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var

    ! volume forces
    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1, 3), intent(inout) :: force

    ! local variables
    integer :: i, j, k, l

    ! gravitational forces
    real :: dRho, theta
    real, dimension(3) :: gForce

    ! Coriolis force
    real :: n1, n2, n3
    real :: u1, u2, u3
    real :: f1, f2, f3
    real :: rho

    ! zonal wind relaxation
    real :: ft_relax
    real :: time
    real, dimension(0:nx + 1) :: fx_relax
    real, dimension(0:ny + 1) :: fy_relax
    integer :: i0, j0

    if(timeScheme /= 'semiimplicit') then
      !--------------------------------------------
      !             Gravitational force
      !--------------------------------------------

      do k = 0, nz + 1
        do j = 0, ny + 1
          do i = 0, nx + 1

            select case(model)

              ! case( "Boussinesq" )

              ! theta = var(i,j,k,6)
              ! gForce = FrInv2*rho00/theta00 * theta * vertical

              ! TFC FJ
              ! No changes for Boussinesq model required.
            case("pseudo_incompressible", "Boussinesq", "compressible")

              if(auxil_equ) then
                dRho = var%rhop(i, j, k)
              else
                if(fluctuationMode) then
                  dRho = var%rho(i, j, k)
                else
                  dRho = var%rho(i, j, k) - rhoStrat(k)
                end if
              end if

              gForce = - FrInv2 * dRho * vertical

            case default
              stop "volumeForce: unknown case model."
            end select

            ! TFC FJ
            ! Adjust gravitational force.
            if(topography) then
              gForce = gForce / jac(i, j, k)
            end if

            force(i, j, k, 1:3) = force(i, j, k, 1:3) + gForce

          end do
        end do
      end do

      !--------------------------------------------
      !             Coriolis force
      !--------------------------------------------
      ! Coriolis force defined in scalar volume cells
      ! -> needs interpolation

      !FS  if( RoInv > 0.0 ) then

      do k = 0, nz + 1
        do j = 0, ny + 1
          do i = 0, nx + 1
            u1 = 0.5 * (var%u(i, j, k) + var%u(i - 1, j, k))
            if(testCase == "SkamarockKlemp94") then
              u1 = u1 - backgroundFlow_dim(1) / uRef
            end if

            u2 = 0.5 * (var%v(i, j, k) + var%v(i, j - 1, k))
            u3 = 0.5 * (var%w(i, j, k) + var%w(i, j, k - 1))

            select case(model)

            case("Boussinesq")
              rho = rho00

            case("pseudo_incompressible", "compressible")
              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  ! Adjust for 3D fields.
                  rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)
                else
                  rho = var%rho(i, j, k) + rhoStrat(k)
                end if
              else
                rho = var%rho(i, j, k)
              end if

            case default
              stop "volumeForce: unknown case model."
            end select

            n1 = vertical(1)
            n2 = vertical(2)
            n3 = vertical(3)

            f1 = n2 * u3 - n3 * u2
            f2 = n3 * u1 - n1 * u3
            f3 = n1 * u2 - n2 * u1

            ! Coriolis force normally written in LHS with "+"
            ! gets now a "-" since force is assumed on the RHS

            force(i, j, k, 1:3) = force(i, j, k, 1:3) - (rho * RoInv(j) * &
                &(/f1, f2, f3/))

            ! TFC FJ
            ! Add vertical Coriolis force component in TFC.
            if(topography) then
              force(i, j, k, 1:3) = force(i, j, k, 1:3) + rho * RoInv(j) &
                  &* vertical * (met(i, j, k, 3, 1) * u2 - met(i, j, k, 3, 2) &
                  &* u1)
            end if

          end do
        end do
      end do

      !FS end if ! RoInv > 0.0
    end if ! not semiimplicit

    !--------------------------------------------
    !             topography growth
    !--------------------------------------------

    if(topography .and. topographyTime > 0.0) then
      do k = 0, nz + 1
        do j = 0, ny + 1
          do i = 0, nx + 1
            ! Zonal part
            force(i, j, k, 3) = force(i, j, k, 3) + (var%rho(i, j, k) &
                &+ rhoStratTFC(i, j, k)) * 0.5 * (var%u(i, j, k) + var%u(i &
                &- 1, j, k)) * tRef / topographyTime &
                &* ((final_topography_surface(i + 1, j) &
                &- final_topography_surface(i - 1, j)) / (2.0 * dx) * (z(k) &
                &- (lz(1) - lz(0))) / (lz(1) - lz(0) - topography_surface(i, &
                &j)) + met(i, j, k, 1, 3) * final_topography_surface(i, j) &
                &/ (lz(1) - lz(0) - topography_surface(i, j)))
            ! Meridional part
            force(i, j, k, 3) = force(i, j, k, 3) + (var%rho(i, j, k) &
                &+ rhoStratTFC(i, j, k)) * 0.5 * (var%v(i, j, k) + var%v(i, j &
                &- 1, k)) * tRef / topographyTime &
                &* ((final_topography_surface(i, j + 1) &
                &- final_topography_surface(i, j - 1)) / (2.0 * dy) * (z(k) &
                &- (lz(1) - lz(0))) / (lz(1) - lz(0) - topography_surface(i, &
                &j)) + met(i, j, k, 2, 3) * final_topography_surface(i, j) &
                &/ (lz(1) - lz(0) - topography_surface(i, j)))
            ! Vertical part
            force(i, j, k, 3) = force(i, j, k, 3) + (var%rho(i, j, k) &
                &+ rhoStratTFC(i, j, k)) * 0.5 * (vertWindTFC(i, j, k, var) &
                &+ vertWindTFC(i, j, k - 1, var)) * tRef / topographyTime &
                &/ jac(i, j, k) * final_topography_surface(i, j) / (lz(1) &
                &- lz(0) - topography_surface(i, j))
          end do
        end do
      end do
    end if

    !--------------------------------------------
    !               wind relaxation
    !--------------------------------------------

    ! TFC FJ
    if(.not. wind_relaxation) return

    i0 = is + nbx - 1
    j0 = js + nby - 1

    if((testCase == "mountainwave") .or. (raytracer .and. case_wkb == 3)) then

      ! if(time < t_ramp) then
      !   ft_relax = (1.0 - cos(time * pi / (t_ramp * 2.0))) / t_relax
      ! else if(time < t_relax - t_ramp) then
      !   ft_relax = 1.0 / t_relax
      ! else if(time < t_relax) then
      !   ft_relax = (1.0 - cos((t_relax - time) * pi / (t_ramp * 2.0))) / t_relax
      ! else
      !   ft_relax = 0.0
      ! end if

      if(xextent_relax > 0.0) then
        do i = 0, nx + 1
          if(x(i0 + i) < lx(0) + 0.5 * xextent_relax) then
            fx_relax(i) = cos((x(i0 + i) - lx(0)) * pi / xextent_relax)
          else if(x(i0 + i) < lx(1) - 0.5 * xextent_relax) then
            fx_relax(i) = 0.0
          else
            fx_relax(i) = cos((lx(1) - x(i0 + i)) * pi / xextent_relax)
          end if
        end do
      else
        do i = 0, nx + 1
          fx_relax(i) = 0.0
        end do
      end if

      if(yextent_relax > 0.0) then
        do j = 0, ny + 1
          if(y(j0 + j) < ly(0) + 0.5 * yextent_relax) then
            fy_relax(j) = cos((y(j0 + j) - ly(0)) * pi / yextent_relax)
          else if(y(j0 + j) < ly(1) - 0.5 * yextent_relax) then
            fy_relax(j) = 0.0
          else
            fy_relax(j) = cos((ly(1) - y(j0 + j)) * pi / yextent_relax)
          end if
        end do
      else
        do j = 0, ny + 1
          fy_relax(j) = 0.0
        end do
      end if

      do k = 0, nz + 1
        do j = 0, ny + 1
          do i = 0, nx + 1

            select case(model)

            case("Boussinesq")
              rho = rho00

            case("pseudo_incompressible", "compressible")

              if(fluctuationMode) then
                if(topography) then
                  ! TFC FJ
                  rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)
                else
                  rho = var%rho(i, j, k) + rhoStrat(k)
                end if
              else
                rho = var%rho(i, j, k)
              end if

            case default
              stop "volumeForce: unknown case model."
            end select

            force(i, j, k, 1) = force(i, j, k, 1) - rho * (0.5 * (var%u(i, j, &
                &k) + var%u(i - 1, j, k)) - u_relax) * fx_relax(i) &
                &* fy_relax(j) / t_relax

            force(i, j, k, 2) = force(i, j, k, 2) - rho * (0.5 * (var%v(i, j, &
                &k) + var%v(i, j - 1, k)) - v_relax) * fx_relax(i) &
                &* fy_relax(j) / t_relax

            force(i, j, k, 3) = force(i, j, k, 3) - rho * (0.5 * (var%w(i, j, &
                &k) + var%w(i, j, k - 1)) - w_relax) * fx_relax(i) &
                &* fy_relax(j) / t_relax

          end do
        end do
      end do

    end if ! mountainwave

  end subroutine volumeForce

  !----------------------------------------------------------------------

  function flux_aldm(uB, uF, vSurf, vL, vR, vBarL, vBarR, sigma)
    !--------------------------------------------
    !   ALDM flux function (cf. Adams, Hickel)
    !--------------------------------------------

    ! in/out arguments
    real, intent(in) :: uB, uF ! velocity to be transported
    real, intent(in) :: vSurf ! averaged cell face transport velocity
    real, intent(in) :: vL, vR ! reconstructed transport velocity
    real, intent(in) :: vBarL, vBarR ! filtered transport velocities
    real, intent(in) :: sigma ! ILES parameter
    real :: flux_aldm

    flux_aldm = 0.5 * (uB + uF) * vSurf - sigma * (vR - vL) * abs(vBarR - vBarL)

  end function flux_aldm

  !---------------------------------------------------------------------------

  function flux_muscl(uSurf, phiUp, phiDown)
    !----------------------------
    !   upwind flux function
    !----------------------------

    ! in/out arguments
    real, intent(in) :: uSurf ! cell face value
    real, intent(in) :: phiUp, phiDown ! upwind, downwind values
    real :: flux_muscl

    !    flux_muscl = uSurf*0.5*(phiUp+phiDown) &
    !         &     + abs(uSurf)*0.5*(phiUp-phiDown)

    if(uSurf > 0.0) then
      flux_muscl = uSurf * phiUp
    else
      flux_muscl = uSurf * phiDown
    end if

  end function flux_muscl

  !-----------------------------------------------------------------------

  subroutine momentumFlux(vara, var, flux, fluxmode, PStrata, PStratTildea)

    !----------------------------------------------------------------------
    ! computes the momentum fluxes at the cell edges using reconstr. values
    ! fluxmode = lin => linear flux, advecting velocities prescribed in
    !                   vara
    !            nln => nonlinear flux, advecting velocities from var
    !
    ! MUSCL assumes that the reconstructed momenta are \rho/P * \vec v
    !----------------------------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: vara, var
    character(len = *), intent(in) :: fluxmode

    type(flux_type), intent(inout) :: flux
    ! flux(i,j,k,dir,iFlux)
    ! dir = 1..3 > f,g,h-flux in x,y,z-direction
    ! iFlux = 1..4 > fRho, fRhoU, fRhoV, fRhoW

    real, dimension(- 1:nz + 2), intent(in) :: PStrata, PStratTildea

    ! local variables
    integer :: i, j, k, l

    ! density at cell edge
    real :: rhoEdge

    ! uTilde at cell edges
    real :: uL, uR, vL, vR, wL, wR ! L=Left at i+1, R=Right at i
    real :: uB, uF, vB, vF, wB, wF ! B=Backward at j+1, F=Forward at j
    real :: uD, uU, vD, vU, wD, wU ! D=Downward at k+1, U=Upward at k

    ! advecting velocities
    real :: uL0, uR0, vL0, vR0, wL0, wR0
    real :: uB0, uF0, vB0, vF0, wB0, wF0
    real :: uD0, uU0, vD0, vU0, wD0, wU0

    ! cell averaged values at cell centres
    real :: uBarL, uBarR, vBarL, vBarR, wBarL, wBarR ! at i and i+1
    real :: uBarB, uBarF, vBarB, vBarF, wBarB, wBarF ! at j and j+1
    real :: uBarD, uBarU, vBarD, vBarU, wBarD, wBarU ! at k and k+1

    ! TFC FJ
    real :: pEdgeR, pREdgeR, pEdgeF, pREdgeF, pEdgeU, pREdgeU, pFEdgeR, &
        &pFEdgeF, pFedgeU, pUEdgeR, pUEdgeF, pUEdgeU
    real :: stressTens13, stressTens13R, stressTens13U, stressTens13RU, &
        &stressTens23, stressTens23F, stressTens23U, stressTens23FU

    ! upwinding
    real :: uSurf, vSurf, wSurf

    ! local flux variables
    real :: fRhoU, gRhoU, hRhoU ! rho*U momentum fluxes
    real :: fRhoV, gRhoV, hRhoV ! rho*V momentum fluxes
    real :: fRhoW, gRhoW, hRhoW ! rho*W momentum fluxes

    ! viscous fluxes
    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1) :: divU ! div(u) field of visc.
    real :: div
    real :: du_dx, du_dy, du_dz ! partial derivatives
    real :: dv_dx, dv_dy, dv_dz
    real :: dw_dx, dw_dy, dw_dz
    real :: fRhoU_visc, gRhoU_visc, hRhoU_visc ! viscous mom. fluxes
    real :: fRhoV_visc, gRhoV_visc, hRhoV_visc
    real :: fRhoW_visc, gRhoW_visc, hRhoW_visc

    ! debugging
    real :: rhoEdge2

    ! avoid abs() for linerisation
    real :: delta

    !for basic-state density
    real :: rhos

    !achatzb for putting together molecular and turbulent vioscosity
    real :: coef_v
    !achatze

    real :: delta_hs, delta_vs

    ! squared grid scales for the anisotropic turbulence scheme

    if(TurbScheme) then
      if(ny == 1 .and. nx == 1) then
        stop 'ERROR: turbulence scheme assumes either nx > 1 or ny > 1'
      else
        if(nx == 1) then
          delta_hs = dy ** 2 ! 2D problems in y and z
        else if(ny == 1) then
          delta_hs = dx ** 2 ! 2D problems in x and z
        else
          delta_hs = dx * dy ! 3D problems

          if(dx / dy > 10.) then
            print *, 'WARNING: dx/dy > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          elseif(dy / dx > 10.) then
            print *, 'WARNING: dy/dx > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          end if
        end if

        delta_vs = dz ** 2
      end if
    end if

    if(verbose) print *, "fluxes.f90/momentumFlux: Entering subroutine..."

    !------------------------------
    !     flux for rho*u
    !------------------------------

    ! flux fRhoU
    do k = 1, nz
      do j = 1, ny
        do i = - 1, nx

          rhoedge = 0.0

          ! density at flux point
          select case(model)

          case("Boussinesq")
            rhoEdge = rho00

          case default

            select case(fluxType)

            case("central")
              ! density interpolation consistent with conti eq

              rhoEdge = 0.25 * (var%rho(i, j, k) + var%rho(i + 1, j, k) &
                  &+ var%rho(i + 1, j, k) + var%rho(i + 2, j, k))

            case("upwind", "ILES")
              ! density interpolation consistent with conti eq
              ! not used in MUSCL case

              if(reconstType /= "MUSCL") then
                rhoEdge = 0.25 * (rhoTilde(i, j, k, 1, 1) + rhoTilde(i + 1, j, &
                    &k, 1, 0) + rhoTilde(i + 1, j, k, 1, 1) + rhoTilde(i + 2, &
                    &j, k, 1, 0))
              end if

            case default
              stop "momentumFlux: unknown fluxType."
            end select

            if(fluctuationMode) rhoEdge = rhoEdge + rhoStrat(k)

          end select ! model

          ! velocity
          select case(fluxType)

          case("central")
            uL = var%u(i, j, k)
            uR = var%u(i + 1, j, k)

            if(fluxmode == "nln") then
              uL0 = var%u(i, j, k)
              uR0 = var%u(i + 1, j, k)
            else if(fluxmode == "lin") then
              uL0 = vara%u(i, j, k)
              uR0 = vara%u(i + 1, j, k)
            else
              stop 'ERROR; wrong fluxmode'
            end if

            fRhoU = 0.25 * (uL + uR) * (uL0 + uR0)

          case("upwind", "ILES")
            ! in MUSCL case the uTilde are the reconstructed
            ! specific momenta, divided by P
            ! in an upwinding scheme these are to be multiplied by
            ! the linearly interpolated velocities (times P) in order
            ! to obtain the desired momentum fluxes

            uR = uTilde(i + 1, j, k, 1, 0)
            uL = uTilde(i, j, k, 1, 1)

            select case(fluxType)

            case("upwind")
              if(topography) then
                ! TFC FJ
                pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, &
                    &j, k) * pStratTFC(i + 1, j, k))
                pREdgeR = 0.5 * (jac(i + 1, j, k) * pStratTFC(i + 1, j, k) &
                    &+ jac(i + 2, j, k) * pStratTFC(i + 2, j, k))
                if(fluxmode == "nln") then
                  uSurf = 0.5 * (pEdgeR * var%u(i, j, k) + pREdgeR * var%u(i &
                      &+ 1, j, k))
                else if(fluxmode == "lin") then
                  uSurf = 0.5 * (pEdgeR * vara%u(i, j, k) + pREdgeR * vara%u(i &
                      &+ 1, j, k))
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                !UAB
                if(fluxmode == "nln") then
                  usurf = 0.5 * (var%u(i, j, k) + var%u(i + 1, j, k)) &
                      &* Pstrat(k)
                else if(fluxmode == "lin") then
                  usurf = 0.5 * (vara%u(i, j, k) + vara%u(i + 1, j, k)) &
                      &* Pstrata(k)
                else
                  stop 'ERROR; wrong fluxmode'
                end if
                !UAE
              end if

              fRhoU = flux_muscl(uSurf, uL, uR)
            case("ILES")
              uBarL = uBar(i, j, k)
              uBarR = uBar(i + 1, j, k)

              if(fluxmode == "nln") then
                uSurf = 0.5 * (uL + uR)
              else if(fluxmode == "lin") then
                usurf = 0.5 * (vara%u(i, j, k) + vara%u(i + 1, j, k))
              else
                stop 'ERROR; wrong fluxmode'
              end if

              fRhoU = flux_aldm(uL, uR, uSurf, uL, uR, uBarL, uBarR, sigmaX)
            case default
              stop "momentumFlux: unknown fluxType."
            end select

          case default
            stop "momentumFlux: unknown fluxType."
          end select

          if(reconstType /= "MUSCL") then
            flux%u(i, j, k, 1) = rhoEdge * fRhoU
          else
            ! for MUSCL case see comment above

            flux%u(i, j, k, 1) = fRhoU
          end if

        end do
      end do
    end do

    !----------------------------------------------------------

    !  flux gRhoU
    do k = 1, nz
      do j = 0, ny
        do i = 0, nx

          rhoedge = 0.0

          ! density at flux point
          select case(model)

          case("Boussinesq")
            rhoEdge = rho00

          case default

            select case(fluxType)

            case("central")
              ! density interpolation consistent with conti eq
              rhoEdge = 0.25 * (var%rho(i + 1, j, k) + var%rho(i + 1, j + 1, &
                  &k) + var%rho(i, j, k) + var%rho(i, j + 1, k))

            case("upwind", "ILES")
              ! density interpolation consistent with conti eq
              ! not used in MUSCL case

              if(reconstType /= "MUSCL") then
                rhoEdge = 0.25 * (rhoTilde(i + 1, j, k, 2, 1) + rhoTilde(i &
                    &+ 1, j + 1, k, 2, 0) + rhoTilde(i, j, k, 2, 1) &
                    &+ rhoTilde(i, j + 1, k, 2, 0))
              end if
            case default
              stop "momentumFlux: unknown fluxType."
            end select

            if(fluctuationMode) rhoEdge = rhoEdge + rhoStrat(k)

          end select ! model

          ! velocity
          select case(fluxType)

          case("central")
            uF = var%u(i, j + 1, k)
            uB = var%u(i, j, k)

            if(fluxmode == "nln") then
              vR = var%v(i + 1, j, k)
              vL = var%v(i, j, k)
            else if(fluxmode == "lin") then
              vR = vara%v(i + 1, j, k)
              vL = vara%v(i, j, k)
            else
              stop 'ERROR; wrong fluxmode'
            end if

            gRhoU = 0.25 * (uB + uF) * (vL + vR)

          case("upwind", "ILES")
            ! in MUSCL case the uTilde are the reconstructed
            ! specific momenta, divided by P
            ! in an upwinding scheme these are to be multiplied by
            ! the linearly interpolated velocities (times P) in order
            ! to obtain the desired momentum fluxes

            uF = uTilde(i, j + 1, k, 2, 0)
            uB = uTilde(i, j, k, 2, 1)

            select case(fluxType)

            case("upwind")
              if(topography) then
                ! TFC FJ
                pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                    &+ 1, k) * pStratTFC(i, j + 1, k))
                pREdgeF = 0.5 * (jac(i + 1, j, k) * pStratTFC(i + 1, j, k) &
                    &+ jac(i + 1, j + 1, k) * pStratTFC(i + 1, j + 1, k))
                if(fluxmode == "nln") then
                  vSurf = 0.5 * (pEdgeF * var%v(i, j, k) + pREdgeF * var%v(i &
                      &+ 1, j, k))
                else if(fluxmode == "lin") then
                  vSurf = 0.5 * (pEdgeF * vara%v(i, j, k) + pREdgeF * vara%v(i &
                      &+ 1, j, k))
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                !UAB
                if(fluxmode == "nln") then
                  vsurf = 0.5 * (var%v(i, j, k) + var%v(i + 1, j, k)) &
                      &* Pstrat(k)
                else if(fluxmode == "lin") then
                  vsurf = 0.5 * (vara%v(i, j, k) + vara%v(i + 1, j, k)) &
                      &* Pstrata(k)
                else
                  stop 'ERROR; wrong fluxmode'
                end if
                !UAE
              end if

              gRhoU = flux_muscl(vSurf, uB, uF)
            case("ILES")
              if(fluxmode == "nln") then
                vR = vTilde(i + 1, j, k, 1, 0)
                vL = vTilde(i, j, k, 1, 1)
              else if(fluxmode == "lin") then
                vR = vara%v(i + 1, j, k)
                vL = vara%v(i, j, k)
              else
                stop 'ERROR; wrong fluxmode'
              end if

              vSurf = 0.5 * (vR + vL)

              uBarF = uBar(i, j + 1, k)
              uBarB = uBar(i, j, k)

              gRhoU = flux_aldm(uB, uF, vSurf, uB, uF, uBarB, uBarF, sigmaX)
            case default
              stop "momentumFlux: unknown fluxType."
            end select

          case default
            stop "momentumFlux: unknown fluxType."
          end select

          if(reconstType /= "MUSCL") then
            flux%u(i, j, k, 2) = rhoEdge * gRhoU
          else
            ! for MUSCL case see comment above

            flux%u(i, j, k, 2) = gRhoU
          end if
        end do
      end do
    end do

    !---------------------------------------------------------------------

    ! flux hRhoU
    do k = 0, nz
      do j = 1, ny
        do i = 0, nx

          rhoedge = 0.0

          ! density at flux point
          select case(model)

          case("Boussinesq")
            rhoEdge = rho00

          case default

            select case(fluxType)

            case("central")
              ! density interpolation consistent with conti eq
              rhoEdge = 0.25 * (var%rho(i, j, k) + var%rho(i, j, k + 1) &
                  &+ var%rho(i + 1, j, k) + var%rho(i + 1, j, k + 1))

            case("upwind", "ILES")
              ! density interpolation consistent with conti eq
              ! not used in MUSCL case

              if(reconstType /= "MUSCL") then
                rhoEdge = 0.25 * (rhoTilde(i, j, k, 3, 1) + rhoTilde(i, j, k &
                    &+ 1, 3, 0) + rhoTilde(i + 1, j, k, 3, 1) + rhoTilde(i &
                    &+ 1, j, k + 1, 3, 0))
              end if

            case default
              stop "momentumFlux: unknown fluxType."
            end select

            ! comment: for CDS rhoEdge should add rhoStrat
            ! for each var(...1) individually for 100% correctness

            if(fluctuationMode) rhoEdge = rhoEdge + rhoStratTilde(k)

          end select ! model

          ! velocity
          select case(fluxType)

          case("central")
            uU = var%u(i, j, k + 1)
            uD = var%u(i, j, k)

            if(fluxmode == "nln") then
              wR = var%w(i + 1, j, k)
              wL = var%w(i, j, k)
            else if(fluxmode == "lin") then
              wR = vara%w(i + 1, j, k)
              wL = vara%w(i, j, k)
            else
              stop 'ERROR; wrong fluxmode'
            end if

            hRhoU = 0.25 * (uD + uU) * (wL + wR)

          case("upwind", "ILES")
            ! in MUSCL case the uTilde are the reconstructed
            ! specific momenta, divided by P
            ! in an upwinding scheme these are to be multiplied by
            ! the linearly interpolated velocities (times P) in order
            ! to obtain the desired momentum fluxes

            uU = uTilde(i, j, k + 1, 3, 0)
            uD = uTilde(i, j, k, 3, 1)

            select case(fluxType)

            case("upwind")
              if(topography) then
                ! TFC FJ
                pEdgeU = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, &
                    &k + 1) * pStratTFC(i, j, k + 1))
                pREdgeU = 0.5 * (jac(i + 1, j, k) * pStratTFC(i + 1, j, k) &
                    &+ jac(i + 1, j, k + 1) * pStratTFC(i + 1, j, k + 1))
                if(fluxmode == "nln") then
                  wSurf = 0.5 * (pEdgeU * var%w(i, j, k) + pREdgeU * var%w(i &
                      &+ 1, j, k))
                else if(fluxmode == "lin") then
                  wSurf = 0.5 * (pEdgeU * vara%w(i, j, k) + pREdgeU * vara%w(i &
                      &+ 1, j, k))
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                !UAB
                if(fluxmode == "nln") then
                  wsurf = 0.5 * (var%w(i, j, k) + var%w(i + 1, j, k)) &
                      &* PstratTilde(k)
                else if(fluxmode == "lin") then
                  wsurf = 0.5 * (vara%w(i, j, k) + vara%w(i + 1, j, k)) &
                      &* PstratTildea(k)
                else
                  stop 'ERROR; wrong fluxmode'
                end if
                !UAE
              end if

              hRhoU = flux_muscl(wSurf, uD, uU)
            case("ILES")
              if(fluxmode == "nln") then
                wR = wTilde(i + 1, j, k, 1, 0)
                wL = wTilde(i, j, k, 1, 1)
              else if(fluxmode == "lin") then
                wR = vara%w(i + 1, j, k)
                wL = vara%w(i, j, k)
              else
                stop 'ERROR; wrong fluxmode'
              end if

              wSurf = 0.5 * (wL + wR)

              uBarU = uBar(i, j, k + 1)
              uBarD = uBar(i, j, k)

              hRhoU = flux_aldm(uU, uD, wSurf, uD, uU, uBarD, uBarU, sigmaX)
            case default
              stop "momentumFlux: unknown fluxType."
            end select

          case default
            stop "momentumFlux: unknown fluxType."
          end select

          if(reconstType /= "MUSCL") then
            flux%u(i, j, k, 3) = rhoEdge * hRhoU
          else
            ! for MUSCL case see comment above

            flux%u(i, j, k, 3) = hRhoU
          end if
        end do
      end do
    end do

    !------------------------------
    !     flux for rho*v
    !------------------------------

    !  flux fRhoV
    do k = 1, nz
      do j = 0, ny
        do i = 0, nx

          rhoedge = 0.0

          ! density at flux point
          select case(model)

          case("Boussinesq")
            rhoEdge = rho00

          case default

            select case(fluxType)

            case("central")
              ! density interpolation consistent with conti eq
              rhoEdge = 0.25 * (var%rho(i, j, k) + var%rho(i + 1, j, k) &
                  &+ var%rho(i, j + 1, k) + var%rho(i + 1, j + 1, k))

            case("upwind", "ILES")
              ! density interpolation consistent with conti eq
              ! not used in MUSCL case

              if(reconstType /= "MUSCL") then
                rhoEdge = 0.25 * (rhoTilde(i, j, k, 1, 1) + rhoTilde(i + 1, j, &
                    &k, 1, 0) + rhoTilde(i, j + 1, k, 1, 1) + rhoTilde(i + 1, &
                    &j + 1, k, 1, 0))
              end if

            case default
              stop "momentumFlux: unknown fluxType."
            end select

            if(fluctuationMode) rhoEdge = rhoEdge + rhoStrat(k)

          end select ! model

          ! velocity
          select case(fluxType)

          case("central")
            vR = var%v(i + 1, j, k)
            vL = var%v(i, j, k)

            if(fluxmode == "nln") then
              uF = var%u(i, j + 1, k)
              uB = var%u(i, j, k)
            else if(fluxmode == "lin") then
              uF = vara%u(i, j + 1, k)
              uB = vara%u(i, j, k)
            else
              stop 'ERROR; wrong fluxmode'
            end if

            fRhoV = 0.25 * (vL + vR) * (uB + uF)

          case("upwind", "ILES")
            ! in MUSCL case the vTilde are the reconstructed
            ! specific momenta, divided by P
            ! in an upwinding scheme these are to be multiplied by
            ! the linearly interpolated velocities (times P) in order
            ! to obtain the desired momentum fluxes

            vR = vTilde(i + 1, j, k, 1, 0)
            vL = vTilde(i, j, k, 1, 1)

            select case(fluxType)

            case("upwind")
              if(topography) then
                ! TFC FJ
                pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, &
                    &j, k) * pStratTFC(i + 1, j, k))
                pFEdgeR = 0.5 * (jac(i, j + 1, k) * pStratTFC(i, j + 1, k) &
                    &+ jac(i + 1, j + 1, k) * pStratTFC(i + 1, j + 1, k))
                if(fluxmode == "nln") then
                  uSurf = 0.5 * (pEdgeR * var%u(i, j, k) + pFEdgeR * var%u(i, &
                      &j + 1, k))
                else if(fluxmode == "lin") then
                  uSurf = 0.5 * (pEdgeR * vara%u(i, j, k) + pFEdgeR &
                      &* vara%u(i, j + 1, k))
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                !UAB
                if(fluxmode == "nln") then
                  usurf = 0.5 * (var%u(i, j, k) + var%u(i, j + 1, k)) &
                      &* Pstrat(k)
                else if(fluxmode == "lin") then
                  usurf = 0.5 * (vara%u(i, j, k) + vara%u(i, j + 1, k)) &
                      &* Pstrata(k)
                else
                  stop 'ERROR; wrong fluxmode'
                end if
                !UAE
              end if

              fRhoV = flux_muscl(uSurf, vL, vR)
            case("ILES")
              if(fluxmode == "nln") then
                uF = uTilde(i, j + 1, k, 2, 0)
                uB = uTilde(i, j, k, 2, 1)
              else if(fluxmode == "lin") then
                uF = vara%u(i, j + 1, k)
                uB = vara%u(i, j, k)
              else
                stop 'ERROR; wrong fluxmode'
              end if

              uSurf = 0.5 * (uB + uF)

              vBarR = vBar(i + 1, j, k)
              vBarL = vBar(i, j, k)

              fRhoV = flux_aldm(vL, vR, uSurf, vL, vR, vBarL, vBarR, sigmaY)
            case default
              stop "momentumFlux: unknown fluxType."
            end select

          case default
            stop "momentumFlux: unknown fluxType."
          end select

          if(reconstType /= "MUSCL") then
            flux%v(i, j, k, 1) = rhoEdge * fRhoV
          else
            ! for MUSCL case see comment above

            flux%v(i, j, k, 1) = fRhoV
          end if
        end do
      end do
    end do

    !---------------------------------------------------------------------

    ! flux gRhoV
    do k = 1, nz
      do j = - 1, ny
        do i = 1, nx

          rhoedge = 0.0

          ! density at flux point
          select case(model)

          case("Boussinesq")
            rhoEdge = rho00

          case default

            select case(fluxType)

            case("central")
              ! density interpolation consistent with conti eq
              rhoEdge = 0.25 * (var%rho(i, j, k) + var%rho(i, j + 1, k) &
                  &+ var%rho(i, j + 1, k) + var%rho(i, j + 2, k))

            case("upwind", "ILES")
              ! density interpolation consistent with conti eq
              ! not used in MUSCL case

              if(reconstType /= "MUSCL") then
                rhoEdge = 0.25 * (rhoTilde(i, j, k, 2, 1) + rhoTilde(i, j + 1, &
                    &k, 2, 0) + rhoTilde(i, j + 1, k, 2, 1) + rhoTilde(i, j &
                    &+ 2, k, 2, 0))
              end if

            case default
              stop "momentumFlux: unknown fluxType."
            end select

            if(fluctuationMode) rhoEdge = rhoEdge + rhoStrat(k)

          end select ! model

          ! velocity
          select case(fluxType)

          case("central")
            vF = var%v(i, j + 1, k)
            vB = var%v(i, j, k)

            if(fluxmode == "nln") then
              vF0 = var%v(i, j + 1, k)
              vB0 = var%v(i, j, k)
            else if(fluxmode == "lin") then
              vF0 = vara%v(i, j + 1, k)
              vB0 = vara%v(i, j, k)
            else
              stop 'ERROR; wrong fluxmode'
            end if

            gRhoV = 0.25 * (vB + vF) * (vB0 + vF0)

          case("upwind", "ILES")
            ! in MUSCL case the vTilde are the reconstructed
            ! specific momenta, divided by P
            ! in an upwinding scheme these are to be multiplied by
            ! the linearly interpolated velocities (times P) in order
            ! to obtain the desired momentum fluxes

            vF = vTilde(i, j + 1, k, 2, 0)
            vB = vTilde(i, j, k, 2, 1)

            select case(fluxType)

            case("upwind")
              if(topography) then
                ! TFC FJ
                pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                    &+ 1, k) * pStratTFC(i, j + 1, k))
                pFEdgeF = 0.5 * (jac(i, j + 1, k) * pStratTFC(i, j + 1, k) &
                    &+ jac(i, j + 2, k) * pStratTFC(i, j + 2, k))
                if(fluxmode == "nln") then
                  vSurf = 0.5 * (pEdgeF * var%v(i, j, k) + pFEdgeF * var%v(i, &
                      &j + 1, k))
                else if(fluxmode == "lin") then
                  vSurf = 0.5 * (pEdgeF * vara%v(i, j, k) + pFEdgeF &
                      &* vara%v(i, j + 1, k))
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                !UAB
                if(fluxmode == "nln") then
                  vsurf = 0.5 * (var%v(i, j, k) + var%v(i, j + 1, k)) &
                      &* Pstrat(k)
                else if(fluxmode == "lin") then
                  vsurf = 0.5 * (vara%v(i, j, k) + vara%v(i, j + 1, k)) &
                      &* Pstrata(k)
                else
                  stop 'ERROR; wrong fluxmode'
                end if
                !UAE
              end if

              gRhoV = flux_muscl(vSurf, vB, vF)
            case("ILES")
              vBarF = vBar(i, j + 1, k)
              vBarB = vBar(i, j, k)

              if(fluxmode == "nln") then
                vSurf = 0.5 * (vB + vF)
              else if(fluxmode == "lin") then
                vsurf = 0.5 * (vara%v(i, j, k) + vara%v(i, j + 1, k))
              else
                stop 'ERROR; wrong fluxmode'
              end if

              gRhoV = flux_aldm(vB, vF, vSurf, vB, vF, vBarB, vBarF, sigmaY)
            case default
              stop "momentumFlux: unknown fluxType."
            end select

          case default
            stop "momentumFlux: unknown fluxType."
          end select

          if(reconstType /= "MUSCL") then
            flux%v(i, j, k, 2) = rhoEdge * gRhoV
          else
            ! for MUSCL case see comment above

            flux%v(i, j, k, 2) = gRhoV
          end if
        end do
      end do
    end do

    !----------------------------------------------------------------------

    ! vertical flux hRhoV
    do k = 0, nz
      do j = 0, ny
        do i = 1, nx

          rhoedge = 0.0

          ! density at flux point
          select case(model)

          case("Boussinesq")
            rhoEdge = rho00

          case default

            select case(fluxType)

            case("central")
              ! density interpolation consistent with conti eq
              rhoEdge = 0.25 * (var%rho(i, j, k) + var%rho(i, j, k + 1) &
                  &+ var%rho(i, j + 1, k) + var%rho(i, j + 1, k + 1))

            case("upwind", "ILES")
              ! density interpolation consistent with conti eq
              ! not used in MUSCL case

              if(reconstType /= "MUSCL") then
                rhoEdge = 0.25 * (rhoTilde(i, j, k, 3, 1) + rhoTilde(i, j, k &
                    &+ 1, 3, 0) + rhoTilde(i, j + 1, k, 3, 1) + rhoTilde(i, j &
                    &+ 1, k + 1, 3, 0))
              end if
            case default
              stop "momentumFlux: unknown fluxType."
            end select

            ! comment: for CDS rhoEdge should add rhoStrat for each
            ! var(...1) individually for 100% consist. with conti eq.

            if(fluctuationMode) rhoEdge = rhoEdge + rhoStratTilde(k)

          end select ! model

          ! velocity
          select case(fluxType)

          case("central")
            vU = var%v(i, j, k + 1)
            vD = var%v(i, j, k)

            if(fluxmode == "nln") then
              wF = var%w(i, j + 1, k)
              wB = var%w(i, j, k)
            else if(fluxmode == "lin") then
              wF = vara%w(i, j + 1, k)
              wB = vara%w(i, j, k)
            else
              stop 'ERROR; wrong fluxmode'
            end if

            hRhoV = 0.25 * (vD + vU) * (wB + wF)

          case("upwind", "ILES")
            ! in MUSCL case the vTilde are the reconstructed
            ! specific momenta, divided by P
            ! in an upwinding scheme these are to be multiplied by
            ! the linearly interpolated velocities (times P) in order
            ! to obtain the desired momentum fluxes

            vU = vTilde(i, j, k + 1, 3, 0)
            vD = vTilde(i, j, k, 3, 1)

            select case(fluxType)

            case("upwind")
              if(topography) then
                ! TFC FJ
                pEdgeU = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, &
                    &k + 1) * pStratTFC(i, j, k + 1))
                pFEdgeU = 0.5 * (jac(i, j + 1, k) * pStratTFC(i, j + 1, k) &
                    &+ jac(i, j + 1, k + 1) * pStratTFC(i, j + 1, k + 1))
                if(fluxmode == "nln") then
                  wSurf = 0.5 * (pEdgeU * var%w(i, j, k) + pFEdgeU * var%w(i, &
                      &j + 1, k))
                else if(fluxmode == "lin") then
                  wSurf = 0.5 * (pEdgeU * vara%w(i, j, k) + pFEdgeU &
                      &* vara%w(i, j + 1, k))
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                !UAB
                if(fluxmode == "nln") then
                  wsurf = 0.5 * (var%w(i, j, k) + var%w(i, j + 1, k)) &
                      &* PstratTilde(k)
                else if(fluxmode == "lin") then
                  wsurf = 0.5 * (vara%w(i, j, k) + vara%w(i, j + 1, k)) &
                      &* PstratTildea(k)
                else
                  stop 'ERROR; wrong fluxmode'
                end if
                !UAE
              end if

              hRhoV = flux_muscl(wSurf, vD, vU)
            case("ILES")
              if(fluxmode == "nln") then
                wF = wTilde(i, j + 1, k, 2, 0)
                wB = wTilde(i, j, k, 2, 1)
              else if(fluxmode == "lin") then
                wF = vara%w(i, j + 1, k)
                wB = vara%w(i, j, k)
              else
                stop 'ERROR; wrong fluxmode'
              end if

              wSurf = 0.5 * (wB + wF)

              vBarU = vBar(i, j, k + 1)
              vBarD = vBar(i, j, k)

              hRhoV = flux_aldm(vD, vU, wSurf, vD, vU, vBarD, vBarU, sigmaY)
            case default
              stop "momentumFlux: unknown fluxType."
            end select

          case default
            stop "momentumFlux: unknown fluxType."
          end select

          if(reconstType /= "MUSCL") then
            flux%v(i, j, k, 3) = rhoEdge * hRhoV
          else
            ! for MUSCL case see comment above

            flux%v(i, j, k, 3) = hRhoV
          end if
        end do
      end do
    end do

    !------------------------------
    !     flux for rho*w
    !------------------------------

    ! flux fRhoW
    do k = 0, nz
      do j = 1, ny
        do i = 0, nx

          rhoedge = 0.0

          ! density at flux point
          select case(model)

          case("Boussinesq")
            rhoEdge = rho00

          case default

            select case(fluxType)

            case("central")

              if(fluctuationMode) then
                rhoEdge = 0.25 * (var%rho(i, j, k) + var%rho(i + 1, j, k) &
                    &+ var%rho(i, j, k + 1) + var%rho(i + 1, j, k + 1)) + 0.5 &
                    &* (rhoStrat(k) + rhoStrat(k + 1))
              else
                rhoEdge = 0.25 * (var%rho(i, j, k) + var%rho(i + 1, j, k) &
                    &+ var%rho(i, j, k + 1) + var%rho(i + 1, j, k + 1))
              end if

            case("upwind", "ILES")
              ! density interpolation consistent with conti eq
              ! not used in MUSCL case

              if(reconstType /= "MUSCL") then
                if(fluctuationMode) then
                  rhoEdge = 0.25 * (rhoTilde(i, j, k, 1, 1) + rhoTilde(i + 1, &
                      &j, k, 1, 0) + rhoTilde(i, j, k + 1, 1, 1) + rhoTilde(i &
                      &+ 1, j, k + 1, 1, 0)) + 0.5 * (rhoStrat(k) + rhoStrat(k &
                      &+ 1))
                else
                  rhoEdge = 0.25 * (rhoTilde(i, j, k, 1, 1) + rhoTilde(i + 1, &
                      &j, k, 1, 0) + rhoTilde(i, j, k + 1, 1, 1) + rhoTilde(i &
                      &+ 1, j, k + 1, 1, 0))
                end if
              end if

            case default
              stop "momentumFlux: unknown fluxType."
            end select

          end select ! model

          ! velocity
          select case(fluxType)

          case("central")
            wR = var%w(i + 1, j, k)
            wL = var%w(i, j, k)

            if(fluxmode == "nln") then
              uU = var%u(i, j, k + 1)
              uD = var%u(i, j, k)
            else if(fluxmode == "lin") then
              uU = vara%u(i, j, k + 1)
              uD = vara%u(i, j, k)
            else
              stop 'ERROR; wrong fluxmode'
            end if

            fRhoW = 0.25 * (wL + wR) * (uD + uU)

          case("upwind", "ILES")
            ! in MUSCL case the wTilde are the reconstructed
            ! specific momenta, divided by P
            ! in an upwinding scheme these are to be multiplied by
            ! the linearly interpolated velocities (times P) in order �
            ! to obtain the desired momentum fluxes

            wR = wTilde(i + 1, j, k, 1, 0)
            wL = wTilde(i, j, k, 1, 1)

            select case(fluxType)

            case("upwind")
              if(topography) then
                ! TFC FJ
                pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, &
                    &j, k) * pStratTFC(i + 1, j, k))
                pUEdgeR = 0.5 * (jac(i, j, k + 1) * pStratTFC(i, j, k + 1) &
                    &+ jac(i + 1, j, k + 1) * pStratTFC(i + 1, j, k + 1))
                if(fluxmode == "nln") then
                  uSurf = 0.5 * (pEdgeR * var%u(i, j, k) + pUEdgeR * var%u(i, &
                      &j, k + 1))
                else if(fluxmode == "lin") then
                  uSurf = 0.5 * (pEdgeR * vara%u(i, j, k) + pUEdgeR &
                      &* vara%u(i, j, k + 1))
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                !UAB
                if(fluxmode == "nln") then
                  usurf = 0.5 * (var%u(i, j, k) * Pstrat(k) + var%u(i, j, k &
                      &+ 1) * Pstrat(k + 1))
                else if(fluxmode == "lin") then
                  usurf = 0.5 * (vara%u(i, j, k) * Pstrata(k) + vara%u(i, j, k &
                      &+ 1) * Pstrata(k + 1))
                else
                  stop 'ERROR; wrong fluxmode'
                end if
                !UAE
              end if

              fRhoW = flux_muscl(uSurf, wL, wR)
            case("ILES")
              if(fluxmode == "nln") then
                uU = uTilde(i, j, k + 1, 3, 0)
                uD = uTilde(i, j, k, 3, 1)
              else if(fluxmode == "lin") then
                uU = vara%u(i, j, k + 1)
                uD = vara%u(i, j, k)
              else
                stop 'ERROR; wrong fluxmode'
              end if

              uSurf = 0.5 * (uD + uU)

              wBarR = wBar(i + 1, j, k)
              wBarL = wBar(i, j, k)

              fRhoW = flux_aldm(wL, wR, uSurf, wL, wR, wBarR, wBarL, sigmaZ)
            case default
              stop "momentumFlux: unknown fluxType."
            end select

          case default
            stop "momentumFlux: unknown fluxType."
          end select

          if(reconstType /= "MUSCL") then
            flux%w(i, j, k, 1) = rhoEdge * fRhoW
          else
            ! for MUSCL case see comment above

            flux%w(i, j, k, 1) = fRhoW
          end if
        end do
      end do
    end do

    !-------------------------------------------------------------------

    ! flux gRhoW
    do k = 0, nz
      do j = 0, ny
        do i = 1, nx

          rhoedge = 0.0

          ! density at flux point
          select case(model)

          case("Boussinesq")
            rhoEdge = rho00

          case default

            select case(fluxType)

            case("central")

              if(fluctuationMode) then
                rhoEdge = 0.25 * (var%rho(i, j, k) + var%rho(i, j + 1, k) &
                    &+ var%rho(i, j, k + 1) + var%rho(i, j + 1, k + 1)) + 0.5 &
                    &* (rhoStrat(k) + rhoStrat(k + 1))
              else
                rhoEdge = 0.25 * (var%rho(i, j, k) + var%rho(i, j + 1, k) &
                    &+ var%rho(i, j, k + 1) + var%rho(i, j + 1, k + 1))
              end if

            case("upwind", "ILES")
              ! density interpolation consistent with conti eq
              ! not used in MUSCL case

              if(reconstType /= "MUSCL") then
                if(fluctuationMode) then
                  rhoEdge = 0.25 * (rhoTilde(i, j, k, 2, 1) + rhoTilde(i, j &
                      &+ 1, k, 2, 0) + rhoTilde(i, j, k + 1, 2, 1) &
                      &+ rhoTilde(i, j + 1, k + 1, 2, 0)) + 0.5 * (rhoStrat(k) &
                      &+ rhoStrat(k + 1))
                else
                  rhoEdge = 0.25 * (rhoTilde(i, j, k, 2, 1) + rhoTilde(i, j &
                      &+ 1, k, 2, 0) + rhoTilde(i, j, k + 1, 2, 1) &
                      &+ rhoTilde(i, j + 1, k + 1, 2, 0))
                end if
              end if
            case default
              stop "momentumFlux: unknown fluxType."
            end select

          end select ! model

          ! velocity
          select case(fluxType)

          case("central")
            wF = var%w(i, j + 1, k)
            wB = var%w(i, j, k)

            if(fluxmode == "nln") then
              vU = var%v(i, j, k + 1)
              vD = var%v(i, j, k)
            else if(fluxmode == "lin") then
              vU = vara%v(i, j, k + 1)
              vD = vara%v(i, j, k)
            else
              stop 'ERROR; wrong fluxmode'
            end if

            gRhoW = 0.25 * (wB + wF) * (vD + vU)

          case("upwind", "ILES")
            ! in MUSCL case the wTilde are the reconstructed
            ! specific momenta, divided by P
            ! in an upwinding scheme these are to be multiplied by
            ! the linearly interpolated velocities (times P) in order
            ! to obtain the desired momentum fluxes

            wF = wTilde(i, j + 1, k, 2, 0)
            wB = wTilde(i, j, k, 2, 1)

            select case(fluxType)

            case("upwind")
              if(topography) then
                ! TFC FJ
                pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                    &+ 1, k) * pStratTFC(i, j + 1, k))
                pUEdgeF = 0.5 * (jac(i, j, k + 1) * pStratTFC(i, j, k + 1) &
                    &+ jac(i, j + 1, k + 1) * pStratTFC(i, j + 1, k + 1))
                if(fluxmode == "nln") then
                  vSurf = 0.5 * (pEdgeF * var%v(i, j, k) + pUEdgeF * var%v(i, &
                      &j, k + 1))
                else if(fluxmode == "lin") then
                  vSurf = 0.5 * (pEdgeF * vara%v(i, j, k) + pUEdgeF &
                      &* vara%v(i, j, k + 1))
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                !UAB
                if(fluxmode == "nln") then
                  vsurf = 0.5 * (var%v(i, j, k) * Pstrat(k) + var%v(i, j, k &
                      &+ 1) * Pstrat(k + 1))
                else if(fluxmode == "lin") then
                  vsurf = 0.5 * (vara%v(i, j, k) * Pstrata(k) + vara%v(i, j, k &
                      &+ 1) * Pstrata(k + 1))
                else
                  stop 'ERROR; wrong fluxmode'
                end if
                !UAE
              end if

              gRhoW = flux_muscl(vSurf, wB, wF)
            case("ILES")
              if(fluxmode == "nln") then
                vU = vTilde(i, j, k + 1, 3, 0)
                vD = vTilde(i, j, k, 3, 1)
              else if(fluxmode == "lin") then
                vU = vara%v(i, j, k + 1)
                vD = vara%v(i, j, k)
              else
                stop 'ERROR; wrong fluxmode'
              end if

              vSurf = 0.5 * (vU + vD)

              wBarF = wBar(i, j + 1, k)
              wBarB = wBar(i, j, k)

              gRhoW = flux_aldm(wB, wF, vSurf, wB, wF, wBarB, wBarF, sigmaZ)
            case default
              stop "momentumFlux: unknown fluxType."
            end select

          case default
            stop "momentumFlux: unknown fluxType."
          end select

          if(reconstType /= "MUSCL") then
            flux%w(i, j, k, 2) = rhoEdge * gRhoW
          else
            ! for MUSCL case see comment above

            flux%w(i, j, k, 2) = gRhoW
          end if
        end do
      end do
    end do

    !---------------------------------------------------------------------

    ! flux hRhoW
    do k = - 1, nz
      do j = 1, ny
        do i = 1, nx

          rhoedge = 0.0

          ! density at flux point
          select case(model)

          case("Boussinesq")
            rhoEdge = rho00

          case default

            select case(fluxType)

            case("central")

              if(fluctuationMode) then
                rhoEdge = 0.25 * (var%rho(i, j, k) + rhoStrat(k) + var%rho(i, &
                    &j, k + 1) + rhoStrat(k + 1) + var%rho(i, j, k + 1) &
                    &+ rhoStrat(k + 1) + var%rho(i, j, k + 2) + rhoStrat(k + 2))

              else
                rhoEdge = 0.25 * (var%rho(i, j, k) + var%rho(i, j, k + 1) &
                    &+ var%rho(i, j, k + 1) + var%rho(i, j, k + 2))
              end if

            case("upwind", "ILES")
              ! density interpolation consistent with conti eq
              ! not used in MUSCL case

              if(reconstType /= "MUSCL") then
                if(fluctuationMode) then
                  rhoEdge = 0.25 * (rhoTilde(i, j, k, 3, 1) + rhoTilde(i, j, k &
                      &+ 1, 3, 0) + rhoTilde(i, j, k + 1, 3, 1) + rhoTilde(i, &
                      &j, k + 2, 3, 0)) + 0.5 * (rhoStratTilde(k) &
                      &+ rhoStratTilde(k + 1))
                else
                  rhoEdge = 0.25 * (rhoTilde(i, j, k, 3, 1) + rhoTilde(i, j, k &
                      &+ 1, 3, 0) + rhoTilde(i, j, k + 1, 3, 1) + rhoTilde(i, &
                      &j, k + 2, 3, 0))
                end if
              end if

            case default
              stop "momentumFlux: unknown fluxType."
            end select

          end select ! model

          ! velocity
          select case(fluxType)

          case("central")
            wU = var%w(i, j, k + 1)
            wD = var%w(i, j, k)

            if(fluxmode == "nln") then
              wU0 = var%w(i, j, k + 1)
              wD0 = var%w(i, j, k)
            else if(fluxmode == "lin") then
              wU0 = vara%w(i, j, k + 1)
              wD0 = vara%w(i, j, k)
            else
              stop 'ERROR; wrong fluxmode'
            end if

            hRhoW = 0.25 * (wD + wU) * (wD0 + wU0)

          case("upwind", "ILES")
            ! in MUSCL case the wTilde are the reconstructed
            ! specific momenta, divided by P
            ! in an upwinding scheme these are to be multiplied by
            ! the linearly interpolated velocities (times P) in order
            ! to obtain the desired momentum fluxes

            wU = wTilde(i, j, k + 1, 3, 0)
            wD = wTilde(i, j, k, 3, 1)

            select case(fluxType)

            case("upwind")
              if(topography) then
                ! TFC FJ
                pEdgeU = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, &
                    &k + 1) * pStratTFC(i, j, k + 1))
                pUEdgeU = 0.5 * (jac(i, j, k + 1) * pStratTFC(i, j, k + 1) &
                    &+ jac(i, j, k + 2) * pStratTFC(i, j, k + 2))
                if(fluxmode == "nln") then
                  wSurf = 0.5 * (pEdgeU * var%w(i, j, k) + pUEdgeU * var%w(i, &
                      &j, k + 1))
                else if(fluxmode == "lin") then
                  wSurf = 0.5 * (pEdgeU * vara%w(i, j, k) + pUEdgeU &
                      &* vara%w(i, j, k + 1))
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                !UAB
                if(fluxmode == "nln") then
                  wsurf = 0.5 * (var%w(i, j, k) * PstratTilde(k) + var%w(i, j, &
                      &k + 1) * PstratTilde(k + 1))
                else if(fluxmode == "lin") then
                  wsurf = 0.5 * (vara%w(i, j, k) * PstratTildea(k) + vara%w(i, &
                      &j, k + 1) * PstratTildea(k + 1))
                else
                  stop 'ERROR; wrong fluxmode'
                end if
                !UAE
              end if

              hRhoW = flux_muscl(wSurf, wD, wU)
            case("ILES")
              wBarU = wBar(i, j, k + 1)
              wBarD = wBar(i, j, k)

              if(fluxmode == "nln") then
                wSurf = 0.5 * (wD + wU)
              else if(fluxmode == "lin") then
                wsurf = 0.5 * (vara%w(i, j, k) + vara%w(i, j, k + 1))
              else
                stop 'ERROR; wrong fluxmode'
              end if

              hRhoW = flux_aldm(wD, wU, wSurf, wD, wU, wBarD, wBarU, sigmaZ)
            case default
              stop "momentumFlux: unknown fluxType."
            end select

          case default
            stop "momentumFlux: unknown fluxType."
          end select

          if(reconstType /= "MUSCL") then
            flux%w(i, j, k, 3) = rhoEdge * hRhoW
          else
            ! for MUSCL case see comment above

            flux%w(i, j, k, 3) = hRhoW
          end if
        end do
      end do
    end do

    !-------------------------------------------------------------------
    !                          Viscous Fluxes
    !-------------------------------------------------------------------

    ! TFC FJ
    ! if(topography .and. mu_viscous_dim == 0.0) then
    !   return
    ! end if

    select case(model)

    case("Boussinesq")

      !------------------------------
      !     flux for rho*u
      !------------------------------

      ! horizontal flux fRhoU
      do k = 1, nz
        do j = 1, ny
          do i = - 1, nx

            ! du/dx at i+1/2
            du_dx = (var%u(i + 1, j, k) - var%u(i, j, k)) / dx

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids
              if(topography) then
                fRhoU_visc = (ReInv + (rho00 * delta_hs * var%DSC(i + 1, j, &
                    &k))) * jac(i + 1, j, k) * stressTensTFC(i + 1, j, k, 1, &
                    &1, var)
              else
                ! fRhoU_visc = (ReInv + (rho00 * delta_hs * var(i, j, k, 7))) &
                !     * (du_dx + du_dx)
                fRhoU_visc = (ReInv + (rho00 * delta_hs * var%DSC(i + 1, j, &
                    &k))) * (du_dx + du_dx)
              end if
            else
              if(topography) then
                fRhoU_visc = ReInv * jac(i + 1, j, k) * stressTensTFC(i + 1, &
                    &j, k, 1, 1, var)
              else
                fRhoU_visc = ReInv * (du_dx + du_dx)
              end if
            end if

            flux%u(i, j, k, 1) = flux%u(i, j, k, 1) - fRhoU_visc
          end do
        end do
      end do

      ! horizontal flux gRhoU
      do k = 1, nz
        do j = 0, ny
          do i = 0, nx
            ! du/dy at j+1/2
            du_dy = (var%u(i, j + 1, k) - var%u(i, j, k)) / dy
            dv_dx = (var%v(i + 1, j, k) - var%v(i, j, k)) / dx ! dv/dx

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids
              if(topography) then
                gRhoU_visc = (ReInv + (rho00 * 0.25 * delta_hs * (var%DSC(i, &
                    &j, k) + var%DSC(i + 1, j, k) + var%DSC(i, j + 1, k) &
                    &+ var%DSC(i + 1, j + 1, k)))) * 0.25 * (jac(i, j, k) &
                    &* stressTensTFC(i, j, k, 1, 2, var) + jac(i + 1, j, k) &
                    &* stressTensTFC(i + 1, j, k, 1, 2, var) + jac(i, j + 1, &
                    &k) * stressTensTFC(i, j + 1, k, 1, 2, var) + jac(i + 1, j &
                    &+ 1, k) * stressTensTFC(i + 1, j + 1, k, 1, 2, var))
              else
                ! gRhoU_visc = (ReInv + (rho00 * delta_hs * var(i, j, k, 7))) &
                !     * (du_dy + dv_dx)
                gRhoU_visc = (ReInv + (rho00 * delta_hs * 0.25 * (var%DSC(i, &
                    &j, k) + var%DSC(i + 1, j, k) + var%DSC(i, j + 1, k) &
                    &+ var%DSC(i + 1, j + 1, k)))) * (du_dy + dv_dx)
              end if
            else
              if(topography) then
                gRhoU_visc = ReInv * 0.25 * (jac(i, j, k) * stressTensTFC(i, &
                    &j, k, 1, 2, var) + jac(i + 1, j, k) * stressTensTFC(i &
                    &+ 1, j, k, 1, 2, var) + jac(i, j + 1, k) &
                    &* stressTensTFC(i, j + 1, k, 1, 2, var) + jac(i + 1, j &
                    &+ 1, k) * stressTensTFC(i + 1, j + 1, k, 1, 2, var))
              else
                gRhoU_visc = ReInv * (du_dy + dv_dx)
              end if
            end if

            flux%u(i, j, k, 2) = flux%u(i, j, k, 2) - gRhoU_visc
          end do
        end do
      end do

      ! vertical flux hRhoU
      do k = 0, nz
        do j = 1, ny
          do i = 0, nx
            ! du/dz  at k+1/2
            du_dz = (var%u(i, j, k + 1) - var%u(i, j, k)) / dz
            dw_dx = (var%w(i + 1, j, k) - var%w(i, j, k)) / dx ! dw/dx

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids
              if(topography) then
                stressTens13 = met(i, j, k, 1, 3) * stressTensTFC(i, j, k, 1, &
                    &1, var) + met(i, j, k, 2, 3) * stressTensTFC(i, j, k, 1, &
                    &2, var) + stressTensTFC(i, j, k, 1, 3, var) / jac(i, j, k)
                stressTens13R = met(i + 1, j, k, 1, 3) * stressTensTFC(i + 1, &
                    &j, k, 1, 1, var) + met(i + 1, j, k, 2, 3) &
                    &* stressTensTFC(i + 1, j, k, 1, 2, var) + stressTensTFC(i &
                    &+ 1, j, k, 1, 3, var) / jac(i + 1, j, k)
                stressTens13U = met(i, j, k + 1, 1, 3) * stressTensTFC(i, j, k &
                    &+ 1, 1, 1, var) + met(i, j, k + 1, 2, 3) &
                    &* stressTensTFC(i, j, k + 1, 1, 2, var) &
                    &+ stressTensTFC(i, j, k + 1, 1, 3, var) / jac(i, j, k + 1)
                stressTens13RU = met(i + 1, j, k + 1, 1, 3) * stressTensTFC(i &
                    &+ 1, j, k + 1, 1, 1, var) + met(i + 1, j, k + 1, 2, 3) &
                    &* stressTensTFC(i + 1, j, k + 1, 1, 2, var) &
                    &+ stressTensTFC(i + 1, j, k + 1, 1, 3, var) / jac(i + 1, &
                    &j, k + 1)
                hRhoU_visc = (ReInv + (rho00 * 0.25 * delta_vs * (jac(i, j, k) &
                    &** 2.0 * var%DSC(i, j, k) + jac(i + 1, j, k) ** 2.0 &
                    &* var%DSC(i + 1, j, k) + jac(i, j, k + 1) ** 2.0 &
                    &* var%DSC(i, j, k + 1) + jac(i + 1, j, k + 1) ** 2.0 &
                    &* var%DSC(i + 1, j, k + 1)))) * 0.25 * (jac(i, j, k) &
                    &* stressTens13 + jac(i + 1, j, k) * stressTens13R &
                    &+ jac(i, j, k + 1) * stressTens13U + jac(i + 1, j, k + 1) &
                    &* stressTens13RU)
              else
                ! hRhoU_visc = (ReInv + (rho00 * delta_vs * var(i, j, k, 7))) &
                !     * (du_dz + dw_dx)
                hRhoU_visc = (ReInv + (rho00 * delta_vs * 0.25 * (var%DSC(i, &
                    &j, k) + var%DSC(i + 1, j, k) + var%DSC(i, j, k + 1) &
                    &+ var%DSC(i + 1, j, k + 1)))) * (du_dz + dw_dx)
              end if
            else
              if(topography) then
                stressTens13 = met(i, j, k, 1, 3) * stressTensTFC(i, j, k, 1, &
                    &1, var) + met(i, j, k, 2, 3) * stressTensTFC(i, j, k, 1, &
                    &2, var) + stressTensTFC(i, j, k, 1, 3, var) / jac(i, j, k)
                stressTens13R = met(i + 1, j, k, 1, 3) * stressTensTFC(i + 1, &
                    &j, k, 1, 1, var) + met(i + 1, j, k, 2, 3) &
                    &* stressTensTFC(i + 1, j, k, 1, 2, var) + stressTensTFC(i &
                    &+ 1, j, k, 1, 3, var) / jac(i + 1, j, k)
                stressTens13U = met(i, j, k + 1, 1, 3) * stressTensTFC(i, j, k &
                    &+ 1, 1, 1, var) + met(i, j, k + 1, 2, 3) &
                    &* stressTensTFC(i, j, k + 1, 1, 2, var) &
                    &+ stressTensTFC(i, j, k + 1, 1, 3, var) / jac(i, j, k + 1)
                stressTens13RU = met(i + 1, j, k + 1, 1, 3) * stressTensTFC(i &
                    &+ 1, j, k + 1, 1, 1, var) + met(i + 1, j, k + 1, 2, 3) &
                    &* stressTensTFC(i + 1, j, k + 1, 1, 2, var) &
                    &+ stressTensTFC(i + 1, j, k + 1, 1, 3, var) / jac(i + 1, &
                    &j, k + 1)
                hRhoU_visc = ReInv * 0.25 * (jac(i, j, k) * stressTens13 &
                    &+ jac(i + 1, j, k) * stressTens13R + jac(i, j, k + 1) &
                    &* stressTens13U + jac(i + 1, j, k + 1) * stressTens13RU)
              else
                hRhoU_visc = ReInv * (du_dz + dw_dx)
              end if
            end if

            flux%u(i, j, k, 3) = flux%u(i, j, k, 3) - hRhoU_visc
          end do
        end do
      end do

      !------------------------------
      !     flux for rho*v
      !------------------------------

      ! horizontal flux fRhoV
      do k = 1, nz
        do j = 0, ny
          do i = 0, nx
            ! dv/dx at i+1/2
            dv_dx = (var%v(i + 1, j, k) - var%v(i, j, k)) / dx
            du_dy = (var%u(i, j + 1, k) - var%u(i, j, k)) / dy ! du/dy

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids
              if(topography) then
                fRhoV_visc = (ReInv + (rho00 * 0.25 * delta_hs * (var%DSC(i, &
                    &j, k) + var%DSC(i + 1, j, k) + var%DSC(i, j + 1, k) &
                    &+ var%DSC(i + 1, j + 1, k)))) * 0.25 * (jac(i, j, k) &
                    &* stressTensTFC(i, j, k, 2, 1, var) + jac(i + 1, j, k) &
                    &* stressTensTFC(i + 1, j, k, 2, 1, var) + jac(i, j + 1, &
                    &k) * stressTensTFC(i, j + 1, k, 2, 1, var) + jac(i + 1, j &
                    &+ 1, k) * stressTensTFC(i + 1, j + 1, k, 2, 1, var))
              else
                ! fRhoV_visc = (ReInv + (rho00 * delta_hs * var(i, j, k, 7))) &
                !     * (dv_dx + du_dy)
                fRhoV_visc = (ReInv + (rho00 * delta_hs * 0.25 * (var%DSC(i, &
                    &j, k) + var%DSC(i + 1, j, k) + var%DSC(i, j + 1, k) &
                    &+ var%DSC(i + 1, j + 1, k)))) * (dv_dx + du_dy)
              end if
            else
              if(topography) then
                fRhoV_visc = ReInv * 0.25 * (jac(i, j, k) * stressTensTFC(i, &
                    &j, k, 2, 1, var) + jac(i + 1, j, k) * stressTensTFC(i &
                    &+ 1, j, k, 2, 1, var) + jac(i, j + 1, k) &
                    &* stressTensTFC(i, j + 1, k, 2, 1, var) + jac(i + 1, j &
                    &+ 1, k) * stressTensTFC(i + 1, j + 1, k, 2, 1, var))
              else
                fRhoV_visc = ReInv * (dv_dx + du_dy)
              end if
            end if

            flux%v(i, j, k, 1) = flux%v(i, j, k, 1) - fRhoV_visc
          end do
        end do
      end do

      ! horizontal flux gRhoV
      do k = 1, nz
        do j = - 1, ny
          do i = 1, nx
            ! dv/dy at j+1/2
            dv_dy = (var%v(i, j + 1, k) - var%v(i, j, k)) / dy

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids
              if(topography) then
                gRhoV_visc = (ReInv + (rho00 * delta_hs * var%DSC(i, j + 1, &
                    &k))) * jac(i, j + 1, k) * stressTensTFC(i, j + 1, k, 2, &
                    &2, var)
              else
                ! gRhoV_visc = (ReInv + (rho00 * delta_hs * var(i, j, k, 7))) &
                !     * (dv_dy + dv_dy)
                gRhoV_visc = (ReInv + (rho00 * delta_hs * var%DSC(i, j + 1, &
                    &k))) * (dv_dy + dv_dy)
              end if
            else
              if(topography) then
                gRhoV_visc = ReInv * jac(i, j + 1, k) * stressTensTFC(i, j &
                    &+ 1, k, 2, 2, var)
              else
                gRhoV_visc = ReInv * (dv_dy + dv_dy)
              end if
            end if

            flux%v(i, j, k, 2) = flux%v(i, j, k, 2) - gRhoV_visc
          end do
        end do
      end do

      ! vertical flux hRhoV
      do k = 0, nz
        do j = 0, ny
          do i = 1, nx
            ! dv/dz at k+1/2
            dv_dz = (var%v(i, j, k + 1) - var%v(i, j, k)) / dz
            dw_dy = (var%w(i, j + 1, k) - var%w(i, j, k)) / dy ! dw/dy

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids
              if(topography) then
                stressTens23 = met(i, j, k, 1, 3) * stressTensTFC(i, j, k, 2, &
                    &1, var) + met(i, j, k, 2, 3) * stressTensTFC(i, j, k, 2, &
                    &2, var) + stressTensTFC(i, j, k, 2, 3, var) / jac(i, j, k)
                stressTens23F = met(i, j + 1, k, 1, 3) * stressTensTFC(i, j &
                    &+ 1, k, 2, 1, var) + met(i, j + 1, k, 2, 3) &
                    &* stressTensTFC(i, j + 1, k, 2, 2, var) &
                    &+ stressTensTFC(i, j + 1, k, 2, 3, var) / jac(i, j + 1, k)
                stressTens23U = met(i, j, k + 1, 1, 3) * stressTensTFC(i, j, k &
                    &+ 1, 2, 1, var) + met(i, j, k + 1, 2, 3) &
                    &* stressTensTFC(i, j, k + 1, 2, 2, var) &
                    &+ stressTensTFC(i, j, k + 1, 2, 3, var) / jac(i, j, k + 1)
                stressTens23FU = met(i, j + 1, k + 1, 1, 3) * stressTensTFC(i, &
                    &j + 1, k + 1, 2, 1, var) + met(i, j + 1, k + 1, 2, 3) &
                    &* stressTensTFC(i, j + 1, k + 1, 2, 2, var) &
                    &+ stressTensTFC(i, j + 1, k + 1, 2, 3, var) / jac(i, j &
                    &+ 1, k + 1)
                hRhoV_visc = (ReInv + (rho00 * 0.25 * delta_vs * (jac(i, j, k) &
                    &** 2.0 * var%DSC(i, j, k) + jac(i, j + 1, k) ** 2.0 &
                    &* var%DSC(i, j + 1, k) + jac(i, j, k + 1) ** 2.0 &
                    &* var%DSC(i, j, k + 1) + jac(i, j + 1, k + 1) ** 2.0 &
                    &* var%DSC(i, j + 1, k + 1)))) * 0.25 * (jac(i, j, k) &
                    &* stressTens23 + jac(i, j + 1, k) * stressTens23F &
                    &+ jac(i, j, k + 1) * stressTens23U + jac(i, j + 1, k + 1) &
                    &* stressTens23FU)
              else
                ! hRhoV_visc = (ReInv + (rho00 * delta_vs * var(i, j, k, 7))) &
                !     * (dv_dz + dw_dy)
                hRhoV_visc = (ReInv + (rho00 * delta_vs * 0.25 * (var%DSC(i, &
                    &j, k) + var%DSC(i, j + 1, k) + var%DSC(i, j, k + 1) &
                    &+ var%DSC(i, j + 1, k + 1)))) * (dv_dz + dw_dy)
              end if
            else
              if(topography) then
                hRhoV_visc = ReInv * 0.25 * (jac(i, j, k) * stressTens23 &
                    &+ jac(i, j + 1, k) * stressTens23F + jac(i, j, k + 1) &
                    &* stressTens23U + jac(i, j + 1, k + 1) * stressTens23FU)
              else
                hRhoV_visc = ReInv * (dv_dz + dw_dy)
              end if
            end if

            flux%v(i, j, k, 3) = flux%v(i, j, k, 3) - hRhoV_visc
          end do
        end do
      end do

      !------------------------------
      !     flux for rho*w
      !------------------------------

      ! horizontal flux fRhoW
      do k = 0, nz
        do j = 1, ny
          do i = 0, nx
            ! dw/dx at i+1/2
            dw_dx = (var%w(i + 1, j, k) - var%w(i, j, k)) / dx
            du_dz = (var%u(i, j, k + 1) - var%u(i, j, k)) / dz ! du/dz

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids
              if(topography) then
                fRhoW_visc = (ReInv + (rho00 * 0.25 * delta_vs * (jac(i, j, k) &
                    &** 2.0 * var%DSC(i, j, k) + jac(i + 1, j, k) ** 2.0 &
                    &* var%DSC(i + 1, j, k) + jac(i, j, k + 1) ** 2.0 &
                    &* var%DSC(i, j, k + 1) + jac(i + 1, j, k + 1) ** 2.0 &
                    &* var%DSC(i + 1, j, k + 1)))) * 0.25 * (jac(i, j, k) &
                    &* stressTensTFC(i, j, k, 3, 1, var) + jac(i + 1, j, k) &
                    &* stressTensTFC(i + 1, j, k, 3, 1, var) + jac(i, j, k &
                    &+ 1) * stressTensTFC(i, j, k + 1, 3, 1, var) + jac(i + 1, &
                    &j, k + 1) * stressTensTFC(i + 1, j, k + 1, 3, 1, var))
              else
                ! fRhoW_visc = (ReInv + (rho00 * delta_vs * var(i, j, k, 7))) &
                !     * (dw_dx + du_dz)
                fRhoW_visc = (ReInv + (rho00 * delta_vs * 0.25 * (var%DSC(i, &
                    &j, k) + var%DSC(i + 1, j, k) + var%DSC(i, j, k + 1) &
                    &+ var%DSC(i + 1, j, k + 1)))) * (dw_dx + du_dz)
              end if
            else
              if(topography) then
                fRhoW_visc = ReInv * 0.25 * (jac(i, j, k) * stressTensTFC(i, &
                    &j, k, 3, 1, var) + jac(i + 1, j, k) * stressTensTFC(i &
                    &+ 1, j, k, 3, 1, var) + jac(i, j, k + 1) &
                    &* stressTensTFC(i, j, k + 1, 3, 1, var) + jac(i + 1, j, k &
                    &+ 1) * stressTensTFC(i + 1, j, k + 1, 3, 1, var))
              else
                fRhoW_visc = ReInv * (dw_dx + du_dz)
              end if
            end if

            flux%w(i, j, k, 1) = flux%w(i, j, k, 1) - fRhoW_visc
          end do
        end do
      end do

      ! horizontal flux gRhoW
      do k = 0, nz
        do j = 0, ny
          do i = 1, nx
            ! dw/dy at j+1/2
            dw_dy = (var%w(i, j + 1, k) - var%w(i, j, k)) / dy
            dv_dz = (var%u(i, j, k + 1) - var%u(i, j, k)) / dz ! dv/dz

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids
              if(topography) then
                gRhoW_visc = (ReInv + (rho00 * 0.25 * delta_vs * (jac(i, j, k) &
                    &** 2.0 * var%DSC(i, j, k) + jac(i, j + 1, k) ** 2.0 &
                    &* var%DSC(i, j + 1, k) + jac(i, j, k + 1) ** 2.0 &
                    &* var%DSC(i, j, k + 1) + jac(i, j + 1, k + 1) ** 2.0 &
                    &* var%DSC(i, j + 1, k + 1)))) * 0.25 * (jac(i, j, k) &
                    &* stressTensTFC(i, j, k, 3, 2, var) + jac(i, j + 1, k) &
                    &* stressTensTFC(i, j + 1, k, 3, 2, var) + jac(i, j, k &
                    &+ 1) * stressTensTFC(i, j, k + 1, 3, 2, var) + jac(i, j &
                    &+ 1, k + 1) * stressTensTFC(i, j + 1, k + 1, 3, 2, var))
              else
                ! gRhoW_visc = (ReInv + (rho00 * delta_vs * var(i, j, k, 7))) &
                !     * (dw_dy + dv_dz)
                gRhoW_visc = (ReInv + (rho00 * delta_vs * 0.25 * (var%DSC(i, &
                    &j, k) + var%DSC(i, j + 1, k) + var%DSC(i, j, k + 1) &
                    &+ var%DSC(i, j + 1, k + 1)))) * (dw_dy + dv_dz)
              end if
            else
              if(topography) then
                gRhoW_visc = ReInv * 0.25 * (jac(i, j, k) * stressTensTFC(i, &
                    &j, k, 3, 2, var) + jac(i, j + 1, k) * stressTensTFC(i, j &
                    &+ 1, k, 3, 2, var) + jac(i, j, k + 1) * stressTensTFC(i, &
                    &j, k + 1, 3, 2, var) + jac(i, j + 1, k + 1) &
                    &* stressTensTFC(i, j + 1, k + 1, 3, 2, var))
              else
                gRhoW_visc = ReInv * (dw_dy + dv_dz)
              end if
            end if

            flux%w(i, j, k, 2) = flux%w(i, j, k, 2) - gRhoW_visc
          end do
        end do
      end do

      ! vertical flux hRhoW
      do k = - 1, nz
        do j = 1, ny
          do i = 1, nx
            ! dw/dz at k+1/2
            dw_dz = (var%w(i, j, k + 1) - var%w(i, j, k)) / dz

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids
              if(topography) then
                hRhoW_visc = (ReInv + (rho00 * delta_vs * jac(i, j, k + 1) &
                    &** 2.0 * var%DSC(i, j, k + 1))) * (jac(i, j, k + 1) &
                    &* met(i, j, k + 1, 1, 3) * stressTensTFC(i, j, k + 1, 3, &
                    &1, var) + jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                    &* stressTensTFC(i, j, k + 1, 3, 2, var) &
                    &+ stressTensTFC(i, j, k + 1, 3, 3, var))
              else
                ! hRhoW_visc = (ReInv + (rho00 * delta_vs * var(i, j, k, 7))) &
                !     * (dw_dz + dw_dz)
                hRhoW_visc = (ReInv + (rho00 * delta_vs * var%DSC(i, j, k &
                    &+ 1))) * (dw_dz + dw_dz)
              end if
            else
              if(topography) then
                hRhoW_visc = ReInv * (jac(i, j, k + 1) * met(i, j, k + 1, 1, &
                    &3) * stressTensTFC(i, j, k + 1, 3, 1, var) + jac(i, j, k &
                    &+ 1) * met(i, j, k + 1, 2, 3) * stressTensTFC(i, j, k &
                    &+ 1, 3, 2, var) + stressTensTFC(i, j, k + 1, 3, 3, var))
              else
                hRhoW_visc = ReInv * (dw_dz + dw_dz)
              end if
            end if

            flux%w(i, j, k, 3) = flux%w(i, j, k, 3) - hRhoW_visc
          end do
        end do
      end do

      !------------------------------------------------------------------

    case("pseudo_incompressible", "compressible")

      !------------------------------------
      !      calc div(u) for viscosity
      !------------------------------------

      do k = 1, nz
        do j = 0, ny + 1 ! j = 0 and ny+1 are div's at ghost cells
          do i = 0, nx + 1 ! same for i
            uR = var%u(i, j, k)
            uL = var%u(i - 1, j, k)
            vF = var%v(i, j, k)
            vB = var%v(i, j - 1, k)
            wU = var%w(i, j, k)
            wD = var%w(i, j, k - 1)

            ! divergence at (i,j,k):
            divU(i, j, k) = (uR - uL) / dx + (vF - vB) / dy + (wU - wD) / dz
          end do
        end do
      end do
      divU(:, :, 0) = divU(:, :, 1) ! set div's in bottom ghost cells
      divU(:, :, nz + 1) = divU(:, :, nz) ! set div's in top ghost cells

      !------------------------------
      !     flux for rho*u
      !------------------------------

      ! horizontal flux fRhoU
      do k = 1, nz
        do j = 1, ny
          do i = - 1, nx
            ! du/dx at i+1/2
            du_dx = (var%u(i + 1, j, k) - var%u(i, j, k)) / dx

            div = divU(i + 1, j, k)

            if(topography) then
              coef_v = ReInv * rhoStratTFC(i + 1, j, 1)
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids

              if(fluctuationMode) then
                if(topography) then
                  coef_v = coef_v + (var%rho(i + 1, j, k) + rhoStratTFC(i + 1, &
                      &j, k)) * delta_hs * var%DSC(i + 1, j, k)
                else
                  ! coef_v = coef_v + (var(i, j, k, 1) + rhoStrat(k)) * delta_hs &
                  !     * var(i, j, k, 7)
                  coef_v = coef_v + (var%rho(i + 1, j, k) + rhoStrat(k)) &
                      &* delta_hs * var%DSC(i + 1, j, k)
                end if
              else
                ! coef_v = coef_v + var(i, j, k, 1) * delta_hs * var(i, j, k, 7)
                coef_v = coef_v + var%rho(i + 1, j, k) * delta_hs * var%DSC(i &
                    &+ 1, j, k)
              end if
            end if

            if(topography) then
              ! TFC FJ
              fRhoU_visc = coef_v * jac(i + 1, j, k) * stressTensTFC(i + 1, j, &
                  &k, 1, 1, var)
            else
              fRhoU_visc = coef_v * (du_dx + du_dx - 2. / 3. * div)
            end if

            flux%u(i, j, k, 1) = flux%u(i, j, k, 1) - fRhoU_visc
          end do
        end do
      end do

      ! horizontal flux gRhoU
      do k = 1, nz
        do j = 0, ny
          do i = 0, nx
            ! du/dy at j+1/2
            du_dy = (var%u(i, j + 1, k) - var%u(i, j, k)) / dy
            dv_dx = (var%v(i + 1, j, k) - var%v(i, j, k)) / dx ! dv/dx

            if(topography) then
              coef_v = ReInv * 0.25 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i &
                  &+ 1, j, 1) + rhoStratTFC(i, j + 1, 1) + rhoStratTFC(i + 1, &
                  &j + 1, 1))
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids

              if(fluctuationMode) then
                if(topography) then
                  coef_v = coef_v + 0.25 * delta_hs * ((var%rho(i, j, k) &
                      &+ rhoStratTFC(i, j, k)) * var%DSC(i, j, k) + (var%rho(i &
                      &+ 1, j, k) + rhoStratTFC(i + 1, j, k)) * var%DSC(i + 1, &
                      &j, k) + (var%rho(i, j + 1, k) + rhoStratTFC(i, j + 1, &
                      &k)) * var%DSC(i, j + 1, k) + (var%rho(i + 1, j + 1, k) &
                      &+ rhoStratTFC(i + 1, j + 1, k)) * var%DSC(i + 1, j + 1, &
                      &k))
                else
                  ! coef_v = coef_v + (var(i, j, k, 1) + rhoStrat(k)) * delta_hs &
                  !     * var(i, j, k, 7)
                  coef_v = coef_v + (0.25 * (var%rho(i, j, k) + var%rho(i + 1, &
                      &j, k) + var%rho(i, j + 1, k) + var%rho(i + 1, j + 1, &
                      &k)) + rhoStrat(k)) * delta_hs * 0.25 * (var%DSC(i, j, &
                      &k) + var%DSC(i + 1, j, k) + var%DSC(i, j + 1, k) &
                      &+ var%DSC(i + 1, j + 1, k))
                end if
              else
                ! coef_v = coef_v + var(i, j, k, 1) * delta_hs * var(i, j, k, 7)
                coef_v = coef_v + 0.25 * (var%rho(i, j, k) + var%rho(i + 1, j, &
                    &k) + var%rho(i, j + 1, k) + var%rho(i + 1, j + 1, k)) &
                    &* delta_hs * 0.25 * (var%DSC(i, j, k) + var%DSC(i + 1, j, &
                    &k) + var%DSC(i, j + 1, k) + var%DSC(i + 1, j + 1, k))
              end if
            end if

            if(topography) then
              ! TFC FJ
              gRhoU_visc = coef_v * 0.25 * (jac(i, j, k) * stressTensTFC(i, j, &
                  &k, 1, 2, var) + jac(i + 1, j, k) * stressTensTFC(i + 1, j, &
                  &k, 1, 2, var) + jac(i, j + 1, k) * stressTensTFC(i, j + 1, &
                  &k, 1, 2, var) + jac(i + 1, j + 1, k) * stressTensTFC(i + 1, &
                  &j + 1, k, 1, 2, var))
            else
              gRhoU_visc = coef_v * (du_dy + dv_dx)
            end if

            flux%u(i, j, k, 2) = flux%u(i, j, k, 2) - gRhoU_visc
          end do
        end do
      end do

      ! vertical flux hRhoU
      do k = 0, nz
        do j = 1, ny
          do i = 0, nx
            ! du/dz at k+1/2
            du_dz = (var%u(i, j, k + 1) - var%u(i, j, k)) / dz
            dw_dx = (var%w(i + 1, j, k) - var%w(i, j, k)) / dx ! dw/dx

            if(topography) then
              coef_v = ReInv * 0.5 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i &
                  &+ 1, j, 1))
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids

              if(fluctuationMode) then
                if(topography) then
                  coef_v = coef_v + 0.25 * delta_vs * ((var%rho(i, j, k) &
                      &+ rhoStratTFC(i, j, k)) * jac(i, j, k) ** 2.0 &
                      &* var%DSC(i, j, k) + (var%rho(i + 1, j, k) &
                      &+ rhoStratTFC(i + 1, j, k)) * jac(i + 1, j, k) ** 2.0 &
                      &* var%DSC(i + 1, j, k) + (var%rho(i, j, k + 1) &
                      &+ rhoStratTFC(i, j, k + 1)) * jac(i, j, k + 1) ** 2.0 &
                      &* var%DSC(i, j, k + 1) + (var%rho(i + 1, j, k + 1) &
                      &+ rhoStratTFC(i + 1, j, k + 1)) * jac(i + 1, j, k + 1) &
                      &** 2.0 * var%DSC(i + 1, j, k + 1))
                else
                  ! coef_v = coef_v + (var(i, j, k, 1) + rhoStrat(k)) * delta_vs &
                  !     * var(i, j, k, 7)
                  coef_v = coef_v + (0.25 * (var%rho(i, j, k) + var%rho(i + 1, &
                      &j, k) + var%rho(i, j, k + 1) + var%rho(i + 1, j, k &
                      &+ 1)) + rhoStratTilde(k)) * delta_vs * 0.25 &
                      &* (var%DSC(i, j, k) + var%DSC(i + 1, j, k) + var%DSC(i, &
                      &j, k + 1) + var%DSC(i + 1, j, k + 1))
                end if
              else
                ! coef_v = coef_v + var(i, j, k, 1) * delta_vs * var(i, j, k, 7)
                coef_v = coef_v + 0.25 * (var%rho(i, j, k) + var%rho(i + 1, j, &
                    &k) + var%rho(i, j, k + 1) + var%rho(i + 1, j, k + 1)) &
                    &* delta_vs * 0.25 * (var%DSC(i, j, k) + var%DSC(i + 1, j, &
                    &k) + var%DSC(i, j, k + 1) + var%DSC(i + 1, j, k + 1))
              end if
            end if

            ! vertical flux hRhoU
            if(topography) then
              ! TFC FJ
              stressTens13 = met(i, j, k, 1, 3) * stressTensTFC(i, j, k, 1, 1, &
                  &var) + met(i, j, k, 2, 3) * stressTensTFC(i, j, k, 1, 2, &
                  &var) + stressTensTFC(i, j, k, 1, 3, var) / jac(i, j, k)
              stressTens13R = met(i + 1, j, k, 1, 3) * stressTensTFC(i + 1, j, &
                  &k, 1, 1, var) + met(i + 1, j, k, 2, 3) * stressTensTFC(i &
                  &+ 1, j, k, 1, 2, var) + stressTensTFC(i + 1, j, k, 1, 3, &
                  &var) / jac(i + 1, j, k)
              stressTens13U = met(i, j, k + 1, 1, 3) * stressTensTFC(i, j, k &
                  &+ 1, 1, 1, var) + met(i, j, k + 1, 2, 3) * stressTensTFC(i, &
                  &j, k + 1, 1, 2, var) + stressTensTFC(i, j, k + 1, 1, 3, &
                  &var) / jac(i, j, k + 1)
              stressTens13RU = met(i + 1, j, k + 1, 1, 3) * stressTensTFC(i &
                  &+ 1, j, k + 1, 1, 1, var) + met(i + 1, j, k + 1, 2, 3) &
                  &* stressTensTFC(i + 1, j, k + 1, 1, 2, var) &
                  &+ stressTensTFC(i + 1, j, k + 1, 1, 3, var) / jac(i + 1, j, &
                  &k + 1)
              hRhoU_visc = coef_v * 0.25 * (jac(i, j, k) * stressTens13 &
                  &+ jac(i + 1, j, k) * stressTens13R + jac(i, j, k + 1) &
                  &* stressTens13U + jac(i + 1, j, k + 1) * stressTens13RU)
            else
              hRhoU_visc = coef_v * (du_dz + dw_dx)
            end if

            flux%u(i, j, k, 3) = flux%u(i, j, k, 3) - hRhoU_visc
          end do
        end do
      end do

      !------------------------------
      !     flux for rho*v
      !------------------------------

      ! horizontal flux fRhoV
      do k = 1, nz
        do j = 0, ny
          do i = 0, nx
            ! dv/dx at i+1/2
            dv_dx = (var%v(i + 1, j, k) - var%v(i, j, k)) / dx
            du_dy = (var%u(i, j + 1, k) - var%u(i, j, k)) / dy ! dv/dy

            if(topography) then
              coef_v = ReInv * 0.25 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i &
                  &+ 1, j, 1) + rhoStratTFC(i, j + 1, 1) + rhoStratTFC(i + 1, &
                  &j + 1, 1))
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids

              if(fluctuationMode) then
                if(topography) then
                  coef_v = coef_v + 0.25 * delta_hs * ((var%rho(i, j, k) &
                      &+ rhoStratTFC(i, j, k)) * var%DSC(i, j, k) + (var%rho(i &
                      &+ 1, j, k) + rhoStratTFC(i + 1, j, k)) * var%DSC(i + 1, &
                      &j, k) + (var%rho(i, j + 1, k) + rhoStratTFC(i, j + 1, &
                      &k)) * var%DSC(i, j + 1, k) + (var%rho(i + 1, j + 1, k) &
                      &+ rhoStratTFC(i + 1, j + 1, k)) * var%DSC(i + 1, j + 1, &
                      &k))
                else
                  ! coef_v = coef_v + (var(i, j, k, 1) + rhoStrat(k)) * delta_hs &
                  !     * var(i, j, k, 7)
                  coef_v = coef_v + (0.25 * (var%rho(i, j, k) + var%rho(i + 1, &
                      &j, k) + var%rho(i, j + 1, k) + var%rho(i + 1, j + 1, &
                      &k)) + rhoStrat(k)) * delta_hs * 0.25 * (var%DSC(i, j, &
                      &k) + var%DSC(i + 1, j, k) + var%DSC(i, j + 1, k) &
                      &+ var%DSC(i + 1, j + 1, k))
                end if
              else
                ! coef_v = coef_v + var(i, j, k, 1) * delta_hs * var(i, j, k, 7)
                coef_v = coef_v + 0.25 * (var%rho(i, j, k) + var%rho(i + 1, j, &
                    &k) + var%rho(i, j + 1, k) + var%rho(i + 1, j + 1, k)) &
                    &* delta_hs * 0.25 * (var%DSC(i, j, k) + var%DSC(i + 1, j, &
                    &k) + var%DSC(i, j + 1, k) + var%DSC(i + 1, j + 1, k))
              end if
            end if

            if(topography) then
              ! TFC FJ
              fRhoV_visc = coef_v * 0.25 * (jac(i, j, k) * stressTensTFC(i, j, &
                  &k, 2, 1, var) + jac(i + 1, j, k) * stressTensTFC(i + 1, j, &
                  &k, 2, 1, var) + jac(i, j + 1, k) * stressTensTFC(i, j + 1, &
                  &k, 2, 1, var) + jac(i + 1, j + 1, k) * stressTensTFC(i + 1, &
                  &j + 1, k, 2, 1, var))
            else
              fRhoV_visc = coef_v * (dv_dx + du_dy)
            end if

            flux%v(i, j, k, 1) = flux%v(i, j, k, 1) - fRhoV_visc
          end do
        end do
      end do

      ! horizontal flux gRhoV
      do k = 1, nz
        do j = - 1, ny
          do i = 1, nx
            ! dv/dy at j+1/2
            dv_dy = (var%v(i, j + 1, k) - var%v(i, j, k)) / dy

            div = divU(i, j + 1, k)

            if(topography) then
              coef_v = ReInv * rhoStratTFC(i, j + 1, 1)
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids

              if(fluctuationMode) then
                if(topography) then
                  coef_v = coef_v + (var%rho(i, j + 1, k) + rhoStratTFC(i, j &
                      &+ 1, k)) * delta_hs * var%DSC(i, j + 1, k)
                else
                  ! coef_v = coef_v + (var(i, j, k, 1) + rhoStrat(k)) * delta_hs &
                  !     * var(i, j, k, 7)
                  coef_v = coef_v + (var%rho(i, j + 1, k) + rhoStrat(k)) &
                      &* delta_hs * var%DSC(i, j + 1, k)
                end if
              else
                ! coef_v = coef_v + var(i, j, k, 1) * delta_hs * var(i, j, k, 7)
                coef_v = coef_v + var%rho(i, j + 1, k) * delta_hs * var%DSC(i, &
                    &j + 1, k)
              end if
            end if

            if(topography) then
              ! TFC FJ
              gRhoV_visc = coef_v * jac(i, j + 1, k) * stressTensTFC(i, j + 1, &
                  &k, 2, 2, var)
            else
              gRhoV_visc = coef_v * (dv_dy + dv_dy - 2. / 3. * div)
            end if

            flux%v(i, j, k, 2) = flux%v(i, j, k, 2) - gRhoV_visc
          end do
        end do
      end do

      ! vertical flux hRhoV
      do k = 0, nz
        do j = 0, ny
          do i = 1, nx
            ! dv/dz at k+1/2
            dv_dz = (var%v(i, j, k + 1) - var%v(i, j, k)) / dz
            dw_dy = (var%w(i, j + 1, k) - var%w(i, j, k)) / dy ! dw/dy

            if(topography) then
              coef_v = ReInv * 0.5 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i, j &
                  &+ 1, 1))
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids

              if(fluctuationMode) then
                if(topography) then
                  coef_v = coef_v + 0.25 * delta_vs * ((var%rho(i, j, k) &
                      &+ rhoStratTFC(i, j, k)) * jac(i, j, k) ** 2.0 &
                      &* var%DSC(i, j, k) + (var%rho(i, j + 1, k) &
                      &+ rhoStratTFC(i, j + 1, k)) * jac(i, j + 1, k) ** 2.0 &
                      &* var%DSC(i, j + 1, k) + (var%rho(i, j, k + 1) &
                      &+ rhoStratTFC(i, j, k + 1)) * jac(i, j, k + 1) ** 2.0 &
                      &* var%DSC(i, j, k + 1) + (var%rho(i, j + 1, k + 1) &
                      &+ rhoStratTFC(i, j + 1, k + 1)) * jac(i, j + 1, k + 1) &
                      &** 2.0 * var%DSC(i, j + 1, k + 1))
                else
                  ! coef_v = coef_v + (var(i, j, k, 1) + rhoStrat(k)) * delta_vs &
                  !     * var(i, j, k, 7)
                  coef_v = coef_v + (0.25 * (var%rho(i, j, k) + var%rho(i, j &
                      &+ 1, k) + var%rho(i, j, k + 1) + var%rho(i, j + 1, k &
                      &+ 1)) + rhoStratTilde(k)) * delta_vs * 0.25 &
                      &* (var%DSC(i, j, k) + var%DSC(i + 1, j, k) + var%DSC(i, &
                      &j, k + 1) + var%DSC(i + 1, j, k + 1))
                end if
              else
                ! coef_v = coef_v + var(i, j, k, 1) * delta_vs * var(i, j, k, 7)
                coef_v = coef_v + 0.25 * (var%rho(i, j, k) + var%rho(i, j + 1, &
                    &k) + var%rho(i, j, k + 1) + var%rho(i, j + 1, k + 1)) &
                    &* delta_vs * 0.25 * (var%DSC(i, j, k) + var%DSC(i, j + 1, &
                    &k) + var%DSC(i, j, k + 1) + var%DSC(i, j + 1, k + 1))
              end if
            end if

            if(topography) then
              ! TFC FJ
              stressTens23 = met(i, j, k, 1, 3) * stressTensTFC(i, j, k, 2, 1, &
                  &var) + met(i, j, k, 2, 3) * stressTensTFC(i, j, k, 2, 2, &
                  &var) + stressTensTFC(i, j, k, 2, 3, var) / jac(i, j, k)
              stressTens23F = met(i, j + 1, k, 1, 3) * stressTensTFC(i, j + 1, &
                  &k, 2, 1, var) + met(i, j + 1, k, 2, 3) * stressTensTFC(i, j &
                  &+ 1, k, 2, 2, var) + stressTensTFC(i, j + 1, k, 2, 3, var) &
                  &/ jac(i, j + 1, k)
              stressTens23U = met(i, j, k + 1, 1, 3) * stressTensTFC(i, j, k &
                  &+ 1, 2, 1, var) + met(i, j, k + 1, 2, 3) * stressTensTFC(i, &
                  &j, k + 1, 2, 2, var) + stressTensTFC(i, j, k + 1, 2, 3, &
                  &var) / jac(i, j, k + 1)
              stressTens23FU = met(i, j + 1, k + 1, 1, 3) * stressTensTFC(i, j &
                  &+ 1, k + 1, 2, 1, var) + met(i, j + 1, k + 1, 2, 3) &
                  &* stressTensTFC(i, j + 1, k + 1, 2, 2, var) &
                  &+ stressTensTFC(i, j + 1, k + 1, 2, 3, var) / jac(i, j + 1, &
                  &k + 1)
              hRhoV_visc = coef_v * 0.25 * (jac(i, j, k) * stressTens23 &
                  &+ jac(i, j + 1, k) * stressTens23F + jac(i, j, k + 1) &
                  &* stressTens23U + jac(i, j + 1, k + 1) * stressTens23FU)
            else
              hRhoV_visc = coef_v * (dv_dz + dw_dy)
            end if

            flux%v(i, j, k, 3) = flux%v(i, j, k, 3) - hRhoV_visc
          end do
        end do
      end do

      !------------------------------
      !     flux for rho*w
      !------------------------------

      ! horizontal flux fRhoW
      do k = 0, nz
        do j = 1, ny
          do i = 0, nx
            ! dw/dx at i     +1/2
            dw_dx = (var%w(i + 1, j, k) - var%w(i, j, k)) / dx
            du_dz = (var%u(i, j, k + 1) - var%u(i, j, k)) / dz ! du/dz

            if(topography) then
              coef_v = ReInv * 0.5 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i &
                  &+ 1, j, 1))
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids

              if(fluctuationMode) then
                if(topography) then
                  coef_v = coef_v + 0.25 * delta_vs * ((var%rho(i, j, k) &
                      &+ rhoStratTFC(i, j, k)) * jac(i, j, k) ** 2.0 &
                      &* var%DSC(i, j, k) + (var%rho(i + 1, j, k) &
                      &+ rhoStratTFC(i + 1, j, k)) * jac(i + 1, j, k) ** 2.0 &
                      &* var%DSC(i + 1, j, k) + (var%rho(i, j, k + 1) &
                      &+ rhoStratTFC(i, j, k + 1)) * jac(i, j, k + 1) ** 2.0 &
                      &* var%DSC(i, j, k + 1) + (var%rho(i + 1, j, k + 1) &
                      &+ rhoStratTFC(i + 1, j, k + 1)) * jac(i + 1, j, k + 1) &
                      &** 2.0 * var%DSC(i + 1, j, k + 1))
                else
                  ! coef_v = coef_v + (var(i, j, k, 1) + rhoStrat(k)) * delta_vs &
                  !     * var(i, j, k, 7)
                  coef_v = coef_v + (0.25 * (var%rho(i, j, k) + var%rho(i + 1, &
                      &j, k) + var%rho(i, j, k + 1) + var%rho(i + 1, j, k &
                      &+ 1)) + rhoStratTilde(k)) * delta_vs * 0.25 &
                      &* (var%DSC(i, j, k) + var%DSC(i + 1, j, k) + var%DSC(i, &
                      &j, k + 1) + var%DSC(i + 1, j, k + 1))
                end if
              else
                ! coef_v = coef_v + var(i, j, k, 1) * delta_vs * var(i, j, k, 7)
                coef_v = coef_v + 0.25 * (var%rho(i, j, k) + var%rho(i + 1, j, &
                    &k) + var%rho(i, j, k + 1) + var%rho(i + 1, j, k + 1)) &
                    &* delta_vs * 0.25 * (var%DSC(i, j, k) + var%DSC(i + 1, j, &
                    &k) + var%DSC(i, j, k + 1) + var%DSC(i + 1, j, k + 1))
              end if
            end if

            if(topography) then
              ! TFC FJ
              fRhoW_visc = coef_v * 0.25 * (jac(i, j, k) * stressTensTFC(i, j, &
                  &k, 3, 1, var) + jac(i + 1, j, k) * stressTensTFC(i + 1, j, &
                  &k, 3, 1, var) + jac(i, j, k + 1) * stressTensTFC(i, j, k &
                  &+ 1, 3, 1, var) + jac(i + 1, j, k + 1) * stressTensTFC(i &
                  &+ 1, j, k + 1, 3, 1, var))
            else
              fRhoW_visc = coef_v * (dw_dx + du_dz)
            end if

            flux%w(i, j, k, 1) = flux%w(i, j, k, 1) - fRhoW_visc
          end do
        end do
      end do

      ! horizontal flux gRhoW
      do k = 0, nz
        do j = 0, ny
          do i = 1, nx
            ! dw/dy at j+1/2
            dw_dy = (var%w(i, j + 1, k) - var%w(i, j, k)) / dy
            dv_dz = (var%u(i, j, k + 1) - var%u(i, j, k)) / dz ! dv/dz

            if(topography) then
              coef_v = ReInv * 0.5 * (rhoStratTFC(i, j, 1) + rhoStratTFC(i, j &
                  &+ 1, 1))
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids

              if(fluctuationMode) then
                if(topography) then
                  coef_v = coef_v + 0.25 * delta_vs * ((var%rho(i, j, k) &
                      &+ rhoStratTFC(i, j, k)) * jac(i, j, k) ** 2.0 &
                      &* var%DSC(i, j, k) + (var%rho(i, j + 1, k) &
                      &+ rhoStratTFC(i, j + 1, k)) * jac(i, j + 1, k) ** 2.0 &
                      &* var%DSC(i, j + 1, k) + (var%rho(i, j, k + 1) &
                      &+ rhoStratTFC(i, j, k + 1)) * jac(i, j, k + 1) ** 2.0 &
                      &* var%DSC(i, j, k + 1) + (var%rho(i, j + 1, k + 1) &
                      &+ rhoStratTFC(i, j + 1, k + 1)) * jac(i, j + 1, k + 1) &
                      &** 2.0 * var%DSC(i, j + 1, k + 1))
                else
                  ! coef_v = coef_v + (var(i, j, k, 1) + rhoStrat(k)) * delta_vs &
                  !     * var(i, j, k, 7)
                  coef_v = coef_v + (0.25 * (var%rho(i, j, k) + var%rho(i, j &
                      &+ 1, k) + var%rho(i, j, k + 1) + var%rho(i, j + 1, k &
                      &+ 1)) + rhoStratTilde(k)) * delta_vs * 0.25 &
                      &* (var%DSC(i, j, k) + var%DSC(i, j + 1, k) + var%DSC(i, &
                      &j, k + 1) + var%DSC(i, j + 1, k + 1))
                end if
              else
                ! coef_v = coef_v + var(i, j, k, 1) * delta_vs * var(i, j, k, 7)
                coef_v = coef_v + 0.25 * (var%rho(i, j, k) + var%rho(i, j + 1, &
                    &k) + var%rho(i, j, k + 1) + var%rho(i, j + 1, k + 1)) &
                    &* delta_vs * 0.25 * (var%DSC(i, j, k) + var%DSC(i, j + 1, &
                    &k) + var%DSC(i, j, k + 1) + var%DSC(i, j + 1, k + 1))
              end if
            end if

            if(topography) then
              ! TFC FJ
              gRhoW_visc = coef_v * 0.25 * (jac(i, j, k) * stressTensTFC(i, j, &
                  &k, 3, 2, var) + jac(i, j + 1, k) * stressTensTFC(i, j + 1, &
                  &k, 3, 2, var) + jac(i, j, k + 1) * stressTensTFC(i, j, k &
                  &+ 1, 3, 2, var) + jac(i, j + 1, k + 1) * stressTensTFC(i, j &
                  &+ 1, k + 1, 3, 2, var))
            else
              gRhoW_visc = coef_v * (dw_dy + dv_dz)
            end if

            flux%w(i, j, k, 2) = flux%w(i, j, k, 2) - gRhoW_visc
          end do
        end do
      end do

      ! vertical flux hRhoW
      do k = - 1, nz
        do j = 1, ny
          do i = 1, nx
            ! dw/dz at k+1/2
            dw_dz = (var%w(i, j, k + 1) - var%w(i, j, k)) / dz

            div = divU(i, j, k + 1)

            if(topography) then
              coef_v = ReInv * rhoStratTFC(i, j, 1)
            else
              coef_v = ReInv * rhoStrat(1)
            end if

            if(TurbScheme) then
              ! turbulence scheme allowing for anisotropic grids

              if(fluctuationMode) then
                if(topography) then
                  coef_v = coef_v + (var%rho(i, j, k + 1) + rhoStratTFC(i, j, &
                      &k + 1)) * delta_vs * jac(i, j, k + 1) ** 2.0 &
                      &* var%DSC(i, j, k + 1)
                else
                  ! if(k == - 1) then
                  !   rhos = rhoStrat(0)
                  ! else
                  !   rhos = rhoStrat(k)
                  ! end if
                  ! coef_v = coef_v + (var(i, j, k, 1) + rhos) * delta_vs &
                  !     * var(i, j, k, 7)
                  coef_v = coef_v + (var%rho(i, j, k + 1) + rhoStrat(k + 1)) &
                      &* delta_vs * var%DSC(i, j, k + 1)
                end if
              else
                ! coef_v = coef_v + var(i, j, k, 1) * delta_vs * var(i, j, k, 7)
                coef_v = coef_v + var%rho(i, j, k + 1) * delta_vs * var%DSC(i, &
                    &j, k + 1)
              end if
            end if

            if(topography) then
              ! TFC FJ
              hRhoW_visc = coef_v * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
                  &* stressTensTFC(i, j, k + 1, 3, 1, var) + jac(i, j, k + 1) &
                  &* met(i, j, k + 1, 2, 3) * stressTensTFC(i, j, k + 1, 3, 2, &
                  &var) + stressTensTFC(i, j, k + 1, 3, 3, var))
            else
              hRhoW_visc = coef_v * (dw_dz + dw_dz - 2. / 3. * div)
            end if

            flux%w(i, j, k, 3) = flux%w(i, j, k, 3) - hRhoW_visc
          end do
        end do
      end do

    case default
      stop "momentumFlux: unknown case model"
    end select

    if(verbose) print *, "fluxes.f90/momentumFlux:  momentum fluxes fRhoU, &
        &fRhoV, fRhoW,  gRhoU, gRhoV, gRhoW hRhoU, hRhoV, hRhoW calculated."

  end subroutine momentumFlux

  !-----------------------------------------------------------------------

  subroutine momentumSource(var, source)
    !---------------------------------------------------------------------
    ! computes the source terms in the momentum equation
    ! 1) div error correction: rhoU * div(u)
    !---------------------------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var

    type(var_type), intent(inout) :: source

    integer :: i, j, k
    real :: uL, uR, uB, uF, uD, uU
    real :: vL, vR, vB, vF, vD, vU
    real :: wL, wR, wB, wF, wD, wU
    real :: uSurfL, vSurfF, wSurfU ! velocities at cell surface
    real :: uSurfR, vSurfB, wSurfD ! velocities at cell surface
    real :: divPu, u, v, w, theta
    real :: PstratU, PstratD, PstratC

    !--------------------------------------------------------
    !              Divergence correction for rhoU
    !--------------------------------------------------------

    do k = 1, nz
      do j = 1, ny
        do i = 0, nx

          select case(fluxType)

          case("central")

            ! no correction implemented
            return

          case("upwind", "ILES")

            uR = uTilde(i + 1, j, k, 1, 0)
            uL = uTilde(i, j, k, 1, 1)
            uSurfR = 0.5 * (uL + uR)

            uR = uTilde(i, j, k, 1, 0)
            uL = uTilde(i - 1, j, k, 1, 1)
            uSurfL = 0.5 * (uL + uR)

            vR = vTilde(i + 1, j, k, 1, 0)
            vL = vTilde(i, j, k, 1, 1)
            vSurfF = 0.5 * (vR + vL)

            vR = vTilde(i + 1, j - 1, k, 1, 0)
            vL = vTilde(i, j - 1, k, 1, 1)
            vSurfB = 0.5 * (vR + vL)

            wR = wTilde(i + 1, j, k, 1, 0)
            wL = wTilde(i, j, k, 1, 1)
            wSurfU = 0.5 * (wL + wR)

            wR = wTilde(i + 1, j, k - 1, 1, 0)
            wL = wTilde(i, j, k - 1, 1, 1)
            wSurfD = 0.5 * (wL + wR)

            PstratU = PstratTilde(k)
            PstratD = PstratTilde(k - 1)
            PstratC = Pstrat(k)

            divPu = PstratC * ((uSurfR - uSurfL) / dx + (vSurfF - vSurfB) &
                &/ dy) + (PstratU * wSurfU - PstratD * wSurfD) / dz

            theta = thetaStrat(k)

            u = var%u(i, j, k)

            source%u(i, j, k) = u * divPu / theta

          case default
            stop "thetaFlux: unknown case fluxType"
          end select

        end do
      end do
    end do

    !--------------------------------------------------------
    !              Divergence correction for rhoV
    !--------------------------------------------------------

    do k = 1, nz
      do j = 0, ny
        do i = 1, nx

          select case(fluxType)

          case("central")

            ! no correction implemented
            return

          case("upwind", "ILES")

            uF = uTilde(i, j + 1, k, 2, 0)
            uB = uTilde(i, j, k, 2, 1)
            uSurfR = 0.5 * (uB + uF)

            uF = uTilde(i - 1, j + 1, k, 2, 0)
            uB = uTilde(i - 1, j, k, 2, 1)
            uSurfL = 0.5 * (uB + uF)

            vF = vTilde(i, j + 1, k, 2, 0)
            vB = vTilde(i, j, k, 2, 1)
            vSurfF = 0.5 * (vB + vF)

            vF = vTilde(i, j, k, 2, 0)
            vB = vTilde(i, j - 1, k, 2, 1)
            vSurfB = 0.5 * (vB + vF)

            wF = wTilde(i, j + 1, k, 2, 0)
            wB = wTilde(i, j, k, 2, 1)
            wSurfU = 0.5 * (wB + wF)

            wF = wTilde(i, j + 1, k - 1, 2, 0)
            wB = wTilde(i, j, k - 1, 2, 1)
            wSurfD = 0.5 * (wB + wF)

            PstratU = PstratTilde(k)
            PstratD = PstratTilde(k - 1)
            PstratC = Pstrat(k)

            divPu = PstratC * ((uSurfR - uSurfL) / dx + (vSurfF - vSurfB) &
                &/ dy) + (PstratU * wSurfU - PstratD * wSurfD) / dz

            theta = thetaStrat(k)

            v = var%v(i, j, k)

            source%v(i, j, k) = v * divPu / theta

          case default
            stop "thetaFlux: unknown case fluxType"
          end select

        end do
      end do
    end do

    !--------------------------------------------------------
    !               Divergence correction for rhoW
    !--------------------------------------------------------

    do k = 1, nz - 1 !xxxx assume solid wall here
      do j = 1, ny
        do i = 0, nx

          select case(fluxType)

          case("central")

            ! no correction implemented
            return

          case("upwind", "ILES")

            uU = uTilde(i, j, k + 1, 3, 0)
            uD = uTilde(i, j, k, 3, 1)
            uSurfR = 0.5 * (uD + uU)

            uU = uTilde(i - 1, j, k + 1, 3, 0)
            uD = uTilde(i - 1, j, k, 3, 1)
            uSurfL = 0.5 * (uD + uU)

            vU = vTilde(i, j, k + 1, 3, 0)
            vD = vTilde(i, j, k, 3, 1)
            vSurfF = 0.5 * (vU + vD)

            vU = vTilde(i, j - 1, k + 1, 3, 0)
            vD = vTilde(i, j - 1, k, 3, 1)
            vSurfB = 0.5 * (vU + vD)

            wU = wTilde(i, j, k + 1, 3, 0)
            wD = wTilde(i, j, k, 3, 1)
            wSurfU = 0.5 * (wD + wU)

            wU = wTilde(i, j, k, 3, 0)
            wD = wTilde(i, j, k - 1, 3, 1)
            wSurfD = 0.5 * (wD + wU)

            PstratU = Pstrat(k + 1)
            PstratD = Pstrat(k)
            PstratC = PstratTilde(k)

            divPu = PstratC * ((uSurfR - uSurfL) / dx + (vSurfF - vSurfB) &
                &/ dy) + (PstratU * wSurfU - PstratD * wSurfD) / dz

            theta = thetaStratTilde(k)

            w = var%w(i, j, k)

            source%w(i, j, k) = w * divPu / theta

          case default
            stop "thetaFlux: unknown case fluxType"
          end select

        end do
      end do
    end do

  end subroutine momentumSource

  !---------------------------------------------------------------------------

  subroutine init_fluxes
    !---------------------------------------------
    ! 1) set parameter for central or ILES flux
    ! 2) allocate flux module variables
    !---------------------------------------------

    ! local variables
    integer :: allocstat

    ! module variables

    ! rhoBar
    allocate(rhoBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhoBar"

    ! rhopBar
    allocate(rhopBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhopBar"

    ! rhoOld
    allocate(rhoOld(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "init_fluxes: alloc of rhoOld failed"

    ! rhopOld
    allocate(rhopOld(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "init_fluxes: alloc of rhopOld failed"

    ! uBar
    allocate(uBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "init_fluxes: could not allocate uBar"

    ! vBar
    allocate(vBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate vBar"

    ! wBar
    allocate(wBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate wBar"

    ! thetaBar
    allocate(thetaBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate thetaBar"

    ! rhoTilde
    allocate(rhoTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
        &0:1), stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhoTilde"

    !UAB
    ! rhoTilde_mom
    allocate(rhoTilde_mom(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
        &0:1), stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhoTilde"
    !UAE

    ! rhopTilde
    allocate(rhopTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
        &0:1), stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate rhopTilde"

    ! uTilde
    allocate(uTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, 0:1), &
        &stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate uTilde"

    ! vTilde
    allocate(vTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, 0:1), &
        &stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate vTilde"

    ! wTilde
    allocate(wTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, 0:1), &
        &stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate wTilde"

    ! thetaTilde
    allocate(thetaTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
        &0:1), stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate thetaTilde"

    if(model == "compressible") then
      allocate(PBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate PBar"

      allocate(PTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
          &0:1), stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate PTilde"

      allocate(pinew(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate pinew"

    end if

    !SD added if
    if(include_ice .or. include_ice2) then
      !  nIceTilde
      allocate(nIceTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
          &0:1), stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate nIceTilde"

      ! qIceTilde
      allocate(qIceTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
          &0:1), stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate qIceTilde"

      ! qvTilde
      allocate(qvTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
          &0:1), stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate qvTilde"

      ! nAerTilde
      allocate(nAerTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
          &0:1), stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate nAerTilde"
    end if

    !SD
    if(include_ice) then
      ! nIceBar
      allocate(nIceBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate nIceBar"

      ! qIceBar
      allocate(qIceBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate qIceBar"

      ! qvBar
      allocate(qvBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate qvBar"

      ! nAerBar
      allocate(nAerBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate nAerBar"
    end if

    if(include_ice2) then

      allocate(IceBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
          &= allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not allocate IceBar"

    end if

    allocate(tracerBar(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), stat &
        &= allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate tracerBar"

    allocate(tracerTilde(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 1:3, &
        &0:1), stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not allocate tracerTilde"

    ! TFC FJ
    if(topography) then

      allocate(uOldTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
          &+ nbz)), stat = allocstat)
      if(allocstat /= 0) stop "init_fluxes: alloc of uOldTFC failed"

      allocate(vOldTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
          &+ nbz)), stat = allocstat)
      if(allocstat /= 0) stop "init_fluxes: alloc of vOldTFC failed"

      allocate(wOldTFC((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
          &+ nbz)), stat = allocstat)
      if(allocstat /= 0) stop "init_fluxes: alloc of wOldTFC failed"

    end if

    if(verbose) print *, "init_fluxes: rhoBar, uBar, vBar, wBar, thetaBar,  &
        &rhoTilde, uTilde, vTilde, wTilde allocated."
  end subroutine init_fluxes

  !---------------------------------------------------------------------------

  subroutine terminate_fluxes
    !-----------------------------------
    ! deallocates flux module variables
    !-----------------------------------

    ! local variables
    integer :: allocstat

    !---------------- deallocate variables -----------------------

    deallocate(rhoBar, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate rhoBar"

    deallocate(rhopBar, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate rhopBar"

    deallocate(rhoOld, stat = allocstat)
    if(allocstat /= 0) stop "terminate_fluxes: dealloc of rhoOld failed"

    deallocate(rhopOld, stat = allocstat)
    if(allocstat /= 0) stop "terminate_fluxes: dealloc of rhopOld failed"

    deallocate(uBar, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate uBar"

    deallocate(vBar, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate vBar"

    deallocate(wBar, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate wBar"

    deallocate(thetaBar, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate thetaBar"

    deallocate(rhoTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate rhoTilde"

    !UAB
    deallocate(rhoTilde_mom, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate rhoTilde"
    !UAE

    deallocate(rhopTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate rhopTilde"

    deallocate(uTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate uTilde"

    deallocate(vTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate vTilde"

    deallocate(wTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate wTilde"

    deallocate(thetaTilde, stat = allocstat)
    if(allocstat /= 0) stop "fluxes.f90: could not deallocate thetaTilde"

    ! SK compressible model
    if(model == "compressible") then
      deallocate(PBar, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate PBar"

      deallocate(PTilde, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate PTilde"

      deallocate(pinew, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate pinew"

    end if

    !SD
    if(include_ice .or. include_ice2) then

      deallocate(nIceTilde, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate nIceTilde"

      deallocate(qIceTilde, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate qIceTilde"

      deallocate(qvTilde, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate qvTilde"

      deallocate(nAerTilde, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate nAerTilde"

    end if

    !SD
    if(include_ice) then

      deallocate(nIceBar, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate nIceBar"

      deallocate(qIceBar, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate qIceBar"

      deallocate(qvBar, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate qvBar"

      deallocate(nAerBar, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate nAerBar"

    endif

    if(include_ice2) then

      deallocate(IceBar, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate IceBar"

    end if

    if(include_tracer) then
      deallocate(tracerBar, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate tracerBar"

      deallocate(tracerTilde, stat = allocstat)
      if(allocstat /= 0) stop "fluxes.f90: could not deallocate tracerTilde"
    end if

    ! TFC FJ
    if(topography) then

      deallocate(uOldTFC, stat = allocstat)
      if(allocstat /= 0) stop "terminate_fluxes: dealloc of uOldTFC failed"

      deallocate(vOldTFC, stat = allocstat)
      if(allocstat /= 0) stop "terminate_fluxes: dealloc of vOldTFC failed"

      deallocate(wOldTFC, stat = allocstat)
      if(allocstat /= 0) stop "terminate_fluxes: dealloc of wOldTFC failed"

    end if

  end subroutine terminate_fluxes

  ! ----------------------------------------------------

  subroutine absDiff(x, absX)
    !----------------------------------------------------------------------
    ! differentiable approx. of absolute value abs()-function for linFloit
    !----------------------------------------------------------------------

    ! in/out variables
    real :: x
    intent(in) :: x
    real :: absX
    intent(out) :: absX

    ! local vars
    real, parameter :: delta0 = 1.0e-18

    if(abs(x) > delta0) then
      absX = abs(x)
    else
      absX = (x ** 2 + delta0 ** 2) / 2.0 / delta0
    end if

  end subroutine absDiff

  ! ----------------------------------------------------

  subroutine setHalosOfField(field)

    !-------------------------------
    !  set values in halo cells
    !-------------------------------

    ! in/out variables
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(inout) :: field

    ! auxiliary fields
    real, dimension(nbx, - nby:ny + nby, - nbz:nz + nbz) :: xSliceLeft_send, &
        &xSliceRight_send
    real, dimension(nbx, - nby:ny + nby, - nbz:nz + nbz) :: xSliceLeft_recv, &
        &xSliceRight_recv

    real, dimension(- nbx:nx + nbx, nby, - nbz:nz + nbz) :: ySliceBack_send, &
        &ySliceForw_send
    real, dimension(- nbx:nx + nbx, nby, - nbz:nz + nbz) :: ySliceBack_recv, &
        &ySliceForw_recv

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! locals
    integer :: i, j, k
    integer :: i0, j0, k0

    !-----------------------------
    !     find neighbour procs
    !-----------------------------

    if(idim > 1) call mpi_cart_shift(comm, 0, 1, left, right, ierror)
    if(jdim > 1) call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

    !------------------------------
    !          x-direction
    !------------------------------

    if(idim > 1) then

      ! slice size
      sendcount = nbx * (ny + 2 * nby + 1) * (nz + 2 * nbz + 1)
      recvcount = sendcount

      ! read slice into contiguous array
      do i = 1, nbx
        xSliceLeft_send(i, :, :) = field(i, :, :)
        xSliceRight_send(i, :, :) = field(nx - nbx + i, :, :)
      end do

      ! left -> right
      source = left
      dest = right
      tag = 100

      i0 = 1; j0 = - nby; k0 = - nbz

      call mpi_sendrecv(xSliceRight_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, xSliceLeft_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xSliceLeft_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, xSliceRight_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_right, ierror)

      ! write auxiliary slice to var field
      do i = 1, nbx
        ! right halos
        field(nx + i, :, :) = xSliceRight_recv(i, :, :)
        ! left halos
        field(- nbx + i, :, :) = xSliceLeft_recv(i, :, :)
      end do

    else

      do i = 1, nbx
        field(nx + i, :, :) = field(i, :, :)
        field(- i + 1, :, :) = field(nx - i + 1, :, :)
      end do

    end if

    !------------------------------
    !          y-direction
    !------------------------------

    if(jdim > 1) then

      ! slice size
      sendcount = nby * (nx + 2 * nbx + 1) * (nz + 2 * nbz + 1)
      recvcount = sendcount

      ! read slice into contiguous array
      do j = 1, nby
        ySliceBack_send(:, j, :) = field(:, j, :)
        ySliceForw_send(:, j, :) = field(:, ny - nby + j, :)
      end do

      ! back -> forw
      source = back
      dest = forw
      tag = 100

      i0 = - nbx; j0 = 1; k0 = - nbz

      call mpi_sendrecv(ySliceForw_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, ySliceBack_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(ySliceBack_send(i0, j0, k0), sendcount, &
          &mpi_double_precision, dest, tag, ySliceForw_recv(i0, j0, k0), &
          &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
          &sts_forw, ierror)

      ! write auxiliary slice to var field
      do j = 1, nby
        ! right halos
        field(:, ny + j, :) = ySliceForw_recv(:, j, :)
        ! left halos
        field(:, - nby + j, :) = ySliceBack_recv(:, j, :)
      end do

    else

      do j = 1, nby
        field(:, ny + j, :) = field(:, j, :)
        field(:, - j + 1, :) = field(:, ny - j + 1, :)
      end do

    end if

  end subroutine setHalosOfField

  ! ----------------------------------------------------

  ! TFC FJ
  subroutine reconstructionTestTFC(var)

    type(var_type), intent(inout) :: var

    real, dimension((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        &+ nbz), 1:3, 0:1) :: rhoTilde_tfc, rhopTilde_tfc, uTilde_tfc, &
        &vTilde_tfc, wTilde_tfc

    call reconstruction(var, "rho")
    rhoTilde_tfc = rhoTilde
    topography = .false.
    call reconstruction(var, "rho")
    topography = .true.
    print *, "Reconstruction difference (rho):", maxval(abs(rhoTilde &
        &- rhoTilde_tfc))

    call reconstruction(var, "rhop")
    rhopTilde_tfc = rhopTilde
    topography = .false.
    call reconstruction(var, "rhop")
    topography = .true.
    print *, "Reconstruction difference (rhop):", maxval(abs(rhopTilde &
        &- rhopTilde_tfc))

    call reconstruction(var, "uvw")
    uTilde_tfc = uTilde
    vTilde_tfc = vTilde
    wTilde_tfc = wTilde
    topography = .false.
    call reconstruction(var, "uvw")
    topography = .true.
    print *, "Reconstruction difference (u):", maxval(abs(uTilde - uTilde_tfc))
    print *, "Reconstruction difference (v):", maxval(abs(vTilde - vTilde_tfc))
    print *, "Reconstruction difference (w):", maxval(abs(wTilde - wTilde_tfc))

  end subroutine reconstructionTestTFC

  ! ----------------------------------------------------

  ! TFC FJ
  subroutine momentumFluxTestTFC(var, flux, RKStage)

    type(var_type), intent(inout) :: var

    type(flux_type), intent(inout) :: flux
    integer :: RKStage

    type(flux_type) :: flux_tfc

    integer :: i, j, k

    real :: fL, fR, gB, gF, hD, hU, fluxDiff
    real :: jacEdgeU
    real, dimension(1:nx, 1:ny, 0:nz) :: fluxDiffW
    real :: sumLoc, sumLocSquared, sumGlob, sumGlobSquared
    real * 4 :: sumOut1, sumOut2
    integer :: recTFC

    logical :: divergenceTest

    divergenceTest = .false.

    if(.not. divergenceTest) then

      call momentumFlux(var, var, flux, "nln", pStrat, pStratTilde)
      flux_tfc = flux
      topography = .false.
      call momentumFlux(var, var, flux, "nln", pStrat, pStratTilde)
      topography = .true.
      print *, "Zonal momentum flux difference (nln):", maxval(abs(flux_tfc%u &
          &- flux%u))
      print *, "Meridional momentum flux difference (nln):", &
          &maxval(abs(flux_tfc%v - flux%v))
      print *, "Vertical momentum flux difference (nln):", &
          &maxval(abs(flux_tfc%w - flux%w))

      call momentumFlux(var, var, flux, "lin", pStrat, pStratTilde)
      flux_tfc = flux
      topography = .false.
      call momentumFlux(var, var, flux, "lin", pStrat, pStratTilde)
      topography = .true.
      print *, "Zonal momentum flux difference (lin):", maxval(abs(flux_tfc%u &
          &- flux%u))
      print *, "Meridional momentum flux difference (lin):", &
          &maxval(abs(flux_tfc%v - flux%v))
      print *, "Vertical momentum flux difference (lin):", &
          &maxval(abs(flux_tfc%w - flux%w))

    else

      do k = 0, nz
        do j = 1, ny
          do i = 1, nx
            ! Compute Cartesian vertical momentum flux divergence.
            fR = flux%w(i, j, k, 1)
            fL = flux%w(i - 1, j, k, 1)
            gF = flux%w(i, j, k, 2)
            gB = flux%w(i, j - 1, k, 2)
            hU = flux%w(i, j, k, 3)
            hD = flux%w(i, j, k - 1, 3)
            fluxDiff = (fR - fL) / dx + (gF - gB) / dy + (hU - hD) / dz
            ! fluxDiff = (hU - hD) / dz
            ! Adjust Cartesian vertical momentum flux divergence.
            jacEdgeU = 0.5 * (jac(i, j, k) + jac(i, j, k + 1))
            ! jacEdgeU = 0.5 * (jac(i, j, k) + jac(i + 1, j, k))
            fluxDiff = fluxDiff / jacEdgeU
            fluxdiffW(i, j, k) = fluxDiff
          end do
        end do
      end do

      sumLoc = 0.0
      sumLocSquared = 0.0
      sumGlob = 0.0
      sumGlobSquared = 0.0

      do k = 0, nz
        do j = 1, ny
          do i = 1, nx
            jacEdgeU = 0.5 * (jac(i, j, k) + jac(i, j, k + 1))
            ! jacEdgeU = 0.5 * (jac(i, j, k) + jac(i + 1, j, k))
            sumLoc = sumLoc + fluxDiffW(i, j, k) * jacEdgeU * dx * dy * dz
            sumLocSquared = sumLocSquared + (fluxDiffW(i, j, k) * jacEdgeU &
                &* dx * dy * dz) ** 2.0
          end do
        end do
      end do

      call mpi_allreduce(sumLoc, sumGlob, 1, mpi_double_precision, mpi_sum, &
          &comm, ierror)
      call mpi_allreduce(sumLocSquared, sumGlobSquared, 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)

      sumOut1 = sumGlob
      if(sumGlobSquared /= 0.0) then
        sumOut2 = sumGlob / sqrt(sumGlobSquared)
      else
        sumOut2 = 0.0
      end if
      if(master) then
        print *, "Absolute error in Cartesian vertical momentum"
        print *, "flux divergence = ", sumOut1
        print *, "Relative error in Cartesian vertical momentum"
        print *, "flux divergence = ", sumOut2
        if(.not. testTFC) then
          recTFC = (iOut - 2) * nStages + RKStage
          open(42, file = "momentum_flux_divergence_error.dat", form &
              &= "unformatted", access = "direct", recl = 2 * sizeofreal4)
          write(42, rec = recTFC) sumOut1, sumOut2
          close(42)
        end if
      end if

    end if

  end subroutine momentumFluxTestTFC

  ! ----------------------------------------------------

  ! TFC FJ
  subroutine massFluxTestTFC(var, flux)

    type(var_type), intent(inout) :: var
    type(flux_type), intent(inout) :: flux

    type(flux_type) :: flux_tfc

    mu_conduct = 0.0

    call massFlux(var, var, flux, "nln", pStrat, pStratTilde)
    flux_tfc = flux
    topography = .false.
    call massFlux(var, var, flux, "nln", pStrat, pStratTilde)
    topography = .true.
    print *, "Mass flux difference (nln):", maxval(abs(flux_tfc%rho - flux%rho))

    call massFlux(var, var, flux, "lin", pStrat, pStratTilde)
    flux_tfc = flux
    topography = .false.
    call massFlux(var, var, flux, "lin", pStrat, pStratTilde)
    topography = .true.
    print *, "Mass flux difference (lin):", maxval(abs(flux_tfc%rho - flux%rho))

  end subroutine massFluxTestTFC

  subroutine ice2Flux(vara, var, flux, fluxmode, Pstrata, PStratTildea)
    !---------------------------------------------------------------------
    ! computes the tracer flux at all cell edges using reconstructed values
    ! fluxmode = lin => linear flux, advecting velocities prescribed in
    !                   vara
    !            nln => nonlinear flux, advecting velocities from var
    !
    ! MUSCL assumes that the reconstructed tracer is tracer/P
    !---------------------------------------------------------------------
    implicit none

    ! in/out variables
    type(var_type), intent(in) :: vara, var
    character(len = *), intent(in) :: fluxmode

    type(flux_type), intent(inout) :: flux
    real, dimension(- 1:nz + 2), intent(in) :: PStrata, PStratTildea

    integer :: i, j, k !,l
    !    real :: rhoL,rhoR, uL,uR       ! L=Left i-1/2, R=Right i+1/2
    !    real :: rhoB,rhoF, vB,vF       ! B=Backward j-1/2, F=Forward j+1/2
    !    real :: rhoD,rhoU, wD,wU       ! D=Downward k-1/2, U=Upward k+1/2
    real :: uSurf, vSurf, wSurf ! velocities at cell surface

    real :: iceL, iceR ! L=Left i-1/2, R=Right i+1/2
    real :: fIce, gIce, hIce

    ! TFC FJ
    real :: rhoStratEdgeR, rhoStratEdgeF, rhoStratEdgeU
    real :: pEdgeR, pEdgeF, pEdgeU

    !    real  :: fRho, gRho, hRho

    ! avoid abs() for linerisation
    real :: delta
    real, parameter :: delta0 = 1.0e-6

    !   variables for the turbulence scheme
    real :: Pr_turb
    real :: coef_t, drho_dxi, dtht_dxi

    real :: delta_hs, delta_vs

    integer :: iVar, ii

    Pr_turb = 1.0 !0.5   !FS

    ! squared grid scales for the anisotropic turbulence scheme

    if(TurbScheme) then
      if(ny == 1 .and. nx == 1) then
        stop 'ERROR: turbulence scheme assumes either nx > 1 or ny > 1'
      else
        if(nx == 1) then
          delta_hs = dy ** 2 ! 2D problems in y and z
        else if(ny == 1) then
          delta_hs = dx ** 2 ! 2D problems in x and z
        else
          delta_hs = dx * dy ! 3D problems

          if(dx / dy > 10.) then
            print *, 'WARNING: dx/dy > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          elseif(dy / dx > 10.) then
            print *, 'WARNING: dy/dx > 10!'
            print *, 'The turbulence scheme is not ready for such  horizontal &
                &grid anisotropies!'
          end if
        end if

        delta_vs = dz ** 2
      end if
    end if

    do iVar = 1, nVarIce

      !-----------------------------------------
      !       Zonal rho fluxes in x: f
      !-----------------------------------------

      do k = 1, nz
        do j = 1, ny
          do i = 0, nx

            select case(fluxType)

            case("central")

              !!$                if( fluctuationMode ) then
              !!$                   rhoL = var(i,j,k,1) + rhoStrat(k)
              !!$                   rhoR = var(i+1,j,k,1) + rhoStrat(k)
              !!$                else
              !!$                   rhoL = var(i,j,k,1)
              !!$                   rhoR = var(i+1,j,k,1)
              !!$                end if

              iceR = var%ICE2(i + 1, j, k, iVar)
              iceL = var%ICE2(i, j, k, iVar)

              if(fluxmode == "nln") then
                uSurf = var%u(i, j, k)
              else if(fluxmode == "lin") then
                uSurf = vara%u(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              fIce = uSurf * 0.5 * (iceL + iceR)

            case("upwind")

              !!$                if( fluctuationMode ) then
              !!$                    if(topography) then
              !!$                        ! TFC FJ
              !!$                        ! Adjust for 3D fields.
              !!$                        rhoStratEdgeR = 0.5 * (rhoStratTFC(i, j, k) &
              !!$                                        + rhoStratTFC(i + 1, j, k))
              !!$                        pEdgeR = 0.5 * (pStratTFC(i, j, k) &
              !!$                                 + pStratTFC(i + 1, j, k))
              !!$                        rhoR = rhoTilde(i + 1, j, k, 1, 0) &
              !!$                               + rhoStratEdgeR / pEdgeR
              !!$                        rhoL = rhoTilde(i, j, k, 1, 1) &
              !!$                               + rhoStratEdgeR / pEdgeR
              !!$                    else
              !!$                        !UAB reference density to be divided by P as well!
              !!$                        rhoR = rhoTilde(i+1,j,k,1,0)+ rhoStrat(k)/Pstrat(k)
              !!$                        rhoL = rhoTilde(i,j,k,1,1)+ rhoStrat(k)/Pstrat(k)
              !!$                        !UAE
              !!$                    end if
              !!$                else
              !!$                   rhoR = rhoTilde(i+1,j,k,1,0)
              !!$                    rhoL = rhoTilde(i,j,k,1,1)
              !!$                end if

              if(iVar .eq. inN) then

                iceR = nIceTilde(i + 1, j, k, 1, 0)
                iceL = nIceTilde(i, j, k, 1, 1)

              elseif(iVar .eq. inQ) then

                iceR = qIceTilde(i + 1, j, k, 1, 0)
                iceL = qIceTilde(i, j, k, 1, 1)

              elseif(iVar .eq. inQv) then

                iceR = qvTilde(i + 1, j, k, 1, 0)
                iceL = qvTilde(i, j, k, 1, 1)

              end if

              if(topography) then
                ! TFC FJ
                pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, &
                    &j, k) * pStratTFC(i + 1, j, k))
                if(fluxmode == "nln") then
                  uSurf = pEdgeR * var%u(i, j, k)
                else if(fluxmode == "lin") then
                  uSurf = pEdgeR * vara%u(i, j, k)
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                if(fluxmode == "nln") then
                  uSurf = var%u(i, j, k) * Pstrat(k) !UA
                else if(fluxmode == "lin") then
                  uSurf = vara%u(i, j, k) * Pstrata(k) !UA
                else
                  stop 'ERROR: worng fluxmode'
                end if
              end if

              fIce = flux_muscl(uSurf, iceL, iceR)

            case("ILES")
              print *, 'ice2Flux does not work with ILES'
              stop

            case default
              stop "ice2Flux: unknown case fluxType"
            end select

            flux%ICE2(i, j, k, 1, iVar) = fIce

          end do
        end do
      end do

      !-----------------------------------------
      !    Meridional rho fluxes in y: g
      !-----------------------------------------

      do k = 1, nz
        do j = 0, ny
          do i = 1, nx

            select case(fluxType)

            case("central")
              !!$                   if( fluctuationMode ) then
              !!$                      rhoF = var(i,j+1,k,1) + rhoStrat(k)
              !!$                      rhoB = var(i,j,k,1)   + rhoStrat(k)
              !!$                   else
              !!$                      rhoF = var(i,j+1,k,1)
              !!$                      rhoB = var(i,j,k,1)
              !!$                   end if

              iceR = var%ICE2(i, j + 1, k, iVar)
              iceL = var%ICE2(i, j, k, iVar)

              if(fluxmode == "nln") then
                vSurf = var%v(i, j, k)
              else if(fluxmode == "lin") then
                vSurf = vara%v(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              gIce = vSurf * 0.5 * (iceL + iceR)

            case("upwind")
              !!$                   if( fluctuationMode ) then
              !!$                      if(topography) then
              !!$                         ! TFC FJ
              !!$                         ! Adjust for 3D fields.
              !!$                         rhoStratEdgeF = 0.5 * (rhoStratTFC(i, j, k) &
              !!$                              + rhoStratTFC(i, j + 1, k))
              !!$                         pEdgeF = 0.5 * (pStratTFC(i, j, k) &
              !!$                              + pStratTFC(i, j + 1, k))
              !!$                         rhoF = rhoTilde(i, j + 1, k, 2, 0) &
              !!$                              + rhoStratEdgeF / pEdgeF
              !!$                         rhoB = rhoTilde(i, j, k, 2, 1) &
              !!$                              + rhoStratEdgeF / pEdgeF
              !!$                      else
              !!$                         !UAB reference density to be divided by P as well!
              !!$                         rhoF = rhoTilde(i,j+1,k,2,0) + rhoStrat(k)/Pstrat(k)
              !!$                         rhoB = rhoTilde(i,j,k,2,1)   + rhoStrat(k)/Pstrat(k)
              !!$                         !UAE
              !!$                      end if
              !!$                   else
              !!$                      rhoF = rhoTilde(i,j+1,k,2,0)
              !!$                      rhoB = rhoTilde(i,j,k,2,1)
              !!$                   end if

              if(iVar .eq. inN) then

                iceR = nIceTilde(i, j + 1, k, 2, 0)
                iceL = nIceTilde(i, j, k, 2, 1)

              elseif(iVar .eq. inQ) then

                iceR = qIceTilde(i, j + 1, k, 2, 0)
                iceL = qIceTilde(i, j, k, 2, 1)

              elseif(iVar .eq. inQv) then

                iceR = qvTilde(i, j + 1, k, 2, 0)
                iceL = qvTilde(i, j, k, 2, 1)

              end if

              if(topography) then
                ! TFC FJ
                pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                    &+ 1, k) * pStratTFC(i, j + 1, k))
                if(fluxmode == "nln") then
                  vSurf = pEdgeF * var%v(i, j, k)
                else if(fluxmode == "lin") then
                  vSurf = pEdgeF * vara%v(i, j, k)
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                if(fluxmode == "nln") then
                  vSurf = var%v(i, j, k) * Pstrat(k) !UA
                else if(fluxmode == "lin") then
                  vSurf = vara%v(i, j, k) * Pstrata(k) !UA
                else
                  stop 'ERROR: worng fluxmode'
                end if
              end if

              gIce = flux_muscl(vSurf, iceL, iceR)

            case("ILES")

              print *, 'ice2Flux does not work with ILES'
              stop

            case default
              stop "ice2Flux: unknown case fluxType"
            end select

            flux%ICE2(i, j, k, 2, iVar) = gIce
          end do
        end do
      end do

      !-----------------------------------------
      !      Vertical rho fluxes in z: h
      !-----------------------------------------

      do k = 0, nz
        do j = 1, ny
          do i = 1, nx

            select case(fluxType)

            case("central")

              !!$                if( fluctuationMode ) then
              !!$                   rhoU = var(i,j,k+1,1) + rhoStratTilde(k)
              !!$                   ! background rho at half level
              !!$                   rhoD = var(i,j,k,1)   + rhoStratTilde(k)
              !!$                else
              !!$                   rhoU = var(i,j,k+1,1)
              !!$                   rhoD = var(i,j,k,1)
              !!$                end if

              iceR = var%ICE2(i, j, k + 1, iVar)
              iceL = var%ICE2(i, j, k, iVar)

              if(fluxmode == "nln") then
                wSurf = var%w(i, j, k)
              else if(fluxmode == "lin") then
                wSurf = vara%w(i, j, k)
              else
                stop 'ERROR: worng fluxmode'
              end if

              hIce = wSurf * 0.5 * (iceL + iceR)

            case("upwind")

              !!$                if( fluctuationMode ) then
              !!$                    if(topography) then
              !!$                        ! TFC FJ
              !!$                        ! Adjust for 3D fields.
              !!$                        rhoStratEdgeU = 0.5 * (rhoStratTFC(i, j, k) &
              !!$                                        + rhoStratTFC(i, j, k + 1))
              !!$                        pEdgeU = 0.5 * (pStratTFC(i, j, k) &
              !!$                                 + pStratTFC(i, j, k + 1))
              !!$                        rhoU = rhoTilde(i, j, k + 1, 3, 0) &
              !!$                               + rhoStratEdgeU / pEdgeU
              !!$                        rhoD = rhoTilde(i, j, k, 3, 1) &
              !!$                               + rhoStratEdgeU / pEdgeU
              !!$                    else
              !!$                        ! background at half level
              !!$                        !UAB reference density to be divided by P as well!
              !!$                        rhoU = rhoTilde(i,j,k+1,3,0) &
              !!$                               + rhoStratTilde(k)/PstratTilde(k)
              !!$                        rhoD = rhoTilde(i,j,k  ,3,1)&
              !!$                               + rhoStratTilde(k)/PstratTilde(k)
              !!$                        !UAE
              !!$                    end if
              !!$                else
              !!$                    rhoU = rhoTilde(i,j,k+1,3,0)
              !!$                    rhoD = rhoTilde(i,j,k,3,1)
              !!$                end if
              !!$
              if(iVar .eq. inN) then

                iceR = nIceTilde(i, j, k + 1, 3, 0)
                iceL = nIceTilde(i, j, k, 3, 1)

              elseif(iVar .eq. inQ) then

                iceR = qIceTilde(i, j, k + 1, 3, 0)
                iceL = qIceTilde(i, j, k, 3, 1)

              elseif(iVar .eq. inQv) then

                iceR = qvTilde(i, j, k + 1, 3, 0)
                iceL = qvTilde(i, j, k, 3, 1)

              end if

              if(topography) then
                ! TFC FJ
                pEdgeU = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, &
                    &k + 1) * pStratTFC(i, j, k + 1))
                if(fluxmode == "nln") then
                  wSurf = pEdgeU * var%w(i, j, k)
                else if(fluxmode == "lin") then
                  wSurf = pEdgeU * vara%w(i, j, k)
                else
                  stop "ERROR: wrong fluxmode"
                end if
              else
                if(fluxmode == "nln") then
                  wSurf = var%w(i, j, k) * PstratTilde(k) !UA
                else if(fluxmode == "lin") then
                  wSurf = vara%w(i, j, k) * PstratTildea(k) !UA
                else
                  stop 'ERROR: worng fluxmode'
                end if
              end if

              hIce = flux_muscl(wSurf, iceL, iceR)

            case("ILES")
              print *, 'ice2Flux does not work with ILES'
              stop

            case default
              stop "rhoFlux: unknown case fluxType"
            end select

            flux%ICE2(i, j, k, 3, iVar) = hIce

          end do
        end do
      end do

    end do ! ii

  end subroutine ice2Flux

end module flux_module
