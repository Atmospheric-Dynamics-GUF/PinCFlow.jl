module poisson_module

  use type_module
  use mpi_module ! modified by Junhong Wei (20161106)
  use timeScheme_module
  use atmosphere_module
  use bicgstab_tools_module
  use sizeof_module
  use mpi

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

  ! TFC
  public :: linearOperatorTestTFC, correctorStepTestTFC

  !------------------------
  !   private subroutines
  !------------------------
  private :: getIndex
  private :: pressureBoundaryCondition
  private :: correctorStep
  private :: linOpr
  private :: bicgstab
  private :: poissonSolver
  private :: thomas ! tridiagonal matrix algoritm (TDMA) / Thomas algorithm

  !-------------------------------
  !    private module variables
  !------------------------------

  ! pressure correction
  real, dimension(:, :, :), allocatable :: dp

  ! solution to Poisson problem
  real, dimension(:, :, :), allocatable :: sol_old1, sol_old2

  ! predicted pressure
  real, dimension(:, :, :), allocatable :: p_pred

  ! tolerance for initial divergence cleaning
  real, parameter :: tolInitial = 1.0e-9

  real :: tol

  ! TFC
  ! Status of Boussinesq tensor elements.
  logical :: expEle, impEle

  contains

  subroutine Corrector(var, flux, dMom, dt, errFlagBicg, nIter, m, opt, &
      &facray, facprs)
    ! -------------------------------------------------
    !              correct uStar, bStar, and p
    ! -------------------------------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var
    type(flux_type), intent(in) :: flux
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 3), &
        &intent(inout) :: dMom

    real, intent(in) :: dt, facray, facprs
    logical, intent(out) :: errFlagBicg
    integer, intent(out) :: nIter
    integer, intent(in) :: m

    ! facray multiplies the Rayleigh-damping terms so that they are only
    ! handled in the implicit time stepping (sponge and immersed boundary)

    ! facprs multiplies the time step so that the routine can be used
    ! properly also in the implicit mode (where in sub-step 5 of the
    ! semi-implicit scheme the pressure correction is over a full
    ! time step, instead of half a time step)

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction
    ! of the winds and density fluctuations
    character(len = *), intent(in) :: opt

    ! local variables
    real, dimension(1:nx, 1:ny, 1:nz) :: rhs ! RHS
    logical :: onlyinfo

    ! Note: dp is a module variable
    ! calc dp

    ! Calc RHS of Poisson problem
    onlyinfo = .false.

    call calc_RHS(rhs, var, flux, dt, onlyinfo)

    call poissonSolver(rhs, var, dt, errFlagBicg, nIter, m, opt, facray, facprs)

    if(errFlagBicg) return

    ! set horizontal and vertical BC for dp
    call pressureBoundaryCondition

    ! correct p, rhopStar, and uStar with dp
    call correctorStep(var, dMom, dt, m, opt, facray, facprs)

    if(detailedinfo) call calc_RHS(rhs, var, flux, dt, detailedinfo)

  end subroutine Corrector

  !----------------------------------------------------------------------

  subroutine preCond(sIn, sOut, opt)
    ! --------------------------------------
    !   preconditioner for BiCGStab
    !   solves vertical problem exploiting its tri-diagonal character
    !   (Isaacson & Keller 1966, see also Durran's book appendix A.2)
    ! --------------------------------------

    ! in/out variables
    real, dimension(1:nx, 1:ny, 1:nz), intent(out) :: sOut
    real, dimension(1:nx, 1:ny, 1:nz), intent(in) :: sIn

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction
    ! of the winds and density fluctuations
    character(len = *), intent(in) :: opt

    ! local field
    real, dimension(1:nx, 1:ny, 1:nz) :: s_pc, q_pc
    real, dimension(1:nx, 1:ny) :: p_pc

    ! local variables
    integer :: k
    integer :: i, j
    integer :: niter

    real :: deta

    ! pseudo timestep

    deta = dtau / (2. * (1. / dx ** 2 + 1. / dy ** 2))

    ! work with auxiliary field s_pc

    s_pc = 0.

    do niter = 1, maxIterADI
      if(niter == 0) then
        s_pc = sIn
      else
        ! Treat all diagonal elements implicitly.
        call linOpr(s_pc, q_pc, opt, 'hnd')

        s_pc = s_pc + deta * (q_pc - sIn)
      end if

      ! upward sweep

      do j = 1, ny
        do i = 1, nx
          au_b(i, j, nz) = 0.0
        end do
      end do

      do j = 1, ny
        do i = 1, nx
          if(niter == 0) then
            q_pc(i, j, 1) = - au_b(i, j, 1) / ac_b(i, j, 1)
            s_pc(i, j, 1) = s_pc(i, j, 1) / ac_b(i, j, 1)
          else
            ! Treat all diagonal elements implicity.
            q_pc(i, j, 1) = deta * au_b(i, j, 1) / (1. - deta * ac_b(i, j, 1))
            s_pc(i, j, 1) = s_pc(i, j, 1) / (1. - deta * ac_b(i, j, 1))
          end if
        end do
      end do

      do k = 2, nz
        do j = 1, ny
          do i = 1, nx
            if(niter == 0) then
              p_pc(i, j) = 1.0 / (ac_b(i, j, k) + ad_b(i, j, k) * q_pc(i, j, k &
                  &- 1))

              q_pc(i, j, k) = - au_b(i, j, k) * p_pc(i, j)

              s_pc(i, j, k) = (s_pc(i, j, k) - ad_b(i, j, k) * s_pc(i, j, k &
                  &- 1)) * p_pc(i, j)
            else
              ! Treat all diagonal elements implicitly.
              p_pc(i, j) = 1.0 / (1. - deta * ac_b(i, j, k) - deta * ad_b(i, &
                  &j, k) * q_pc(i, j, k - 1))

              q_pc(i, j, k) = deta * au_b(i, j, k) * p_pc(i, j)

              s_pc(i, j, k) = (s_pc(i, j, k) + deta * ad_b(i, j, k) * s_pc(i, &
                  &j, k - 1)) * p_pc(i, j)
            end if
          end do
        end do
      end do

      ! backward pass

      do k = nz - 1, 1, - 1
        do j = 1, ny
          do i = 1, nx
            s_pc(i, j, k) = s_pc(i, j, k) + q_pc(i, j, k) * s_pc(i, j, k + 1)
          end do
        end do
      end do
    end do

    ! final result

    sOut = s_pc

    return

  end subroutine preCond

  !----------------------------------------------------------------------

  subroutine preCondExperiment(sIn, sOut, opt)
    ! --------------------------------------
    !   preconditioner for BiCGStab
    !   solves vertical problem exploiting its tri-diagonal character
    !   (Isaacson & Keller 1966, see also Durran's book appendix A.2)
    ! --------------------------------------

    ! in/out variables
    real, dimension(1:nx, 1:ny, 1:nz), intent(out) :: sOut
    real, dimension(1:nx, 1:ny, 1:nz), intent(in) :: sIn

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction
    ! of the winds and density fluctuations
    character(len = *), intent(in) :: opt

    ! local field
    real, dimension(1:nx, 1:ny, 1:nz) :: s_pc, q_pc
    real, dimension(1:nx, 1:ny) :: p_pc

    ! local variables
    integer :: k
    integer :: i, j
    integer :: niter

    real :: deta

    ! pseudo timestep

    deta = dtau / (2. * (1. / dx ** 2 + 1. / dy ** 2))

    s_pc = 0.

    do niter = 1, maxIterADI
      if(niter == 0) then
        s_pc = sIn
      else
        call linOpr(s_pc, q_pc, opt, 'hnd')

        s_pc = ach_b * s_pc + deta * (q_pc - sIn)
      end if

      ! upward sweep

      do j = 1, ny
        do i = 1, nx
          au_b(i, j, nz) = 0.0
        end do
      end do

      do j = 1, ny
        do i = 1, nx
          if(niter == 0) then
            q_pc(i, j, 1) = - au_b(i, j, 1) / ac_b(i, j, 1)
            s_pc(i, j, 1) = s_pc(i, j, 1) / ac_b(i, j, 1)
          else
            q_pc(i, j, 1) = deta * au_b(i, j, 1) / (ach_b(i, j, 1) - deta &
                &* ac_b(i, j, 1))
            s_pc(i, j, 1) = s_pc(i, j, 1) / (ach_b(i, j, 1) - deta * ac_b(i, &
                &j, 1))
          end if
        end do
      end do

      do k = 2, nz
        do j = 1, ny
          do i = 1, nx
            if(niter == 0) then
              p_pc(i, j) = 1.0 / (ac_b(i, j, k) + ad_b(i, j, k) * q_pc(i, j, k &
                  &- 1))

              q_pc(i, j, k) = - au_b(i, j, k) * p_pc(i, j)

              s_pc(i, j, k) = (s_pc(i, j, k) - ad_b(i, j, k) * s_pc(i, j, k &
                  &- 1)) * p_pc(i, j)
            else
              p_pc(i, j) = 1.0 / (ach_b(i, j, k) - deta * ac_b(i, j, k) - deta &
                  &* ad_b(i, j, k) * q_pc(i, j, k - 1))

              q_pc(i, j, k) = deta * au_b(i, j, k) * p_pc(i, j)

              s_pc(i, j, k) = (s_pc(i, j, k) + deta * ad_b(i, j, k) * s_pc(i, &
                  &j, k - 1)) * p_pc(i, j)
            end if
          end do
        end do
      end do

      ! backward pass

      do k = nz - 1, 1, - 1
        do j = 1, ny
          do i = 1, nx
            s_pc(i, j, k) = s_pc(i, j, k) + q_pc(i, j, k) * s_pc(i, j, k + 1)
          end do
        end do
      end do
    end do

    ! final result

    sOut = s_pc

    return

  end subroutine preCondExperiment

  !----------------------------------------------------------------------

  subroutine linOpr(sIn, Ls, opt, hortot)
    ! --------------------------------------
    !   Linear Operator in Poisson problem
    !   Functions as A*x
    ! --------------------------------------

    ! in/out variables
    real, dimension(1:nx, 1:ny, 1:nz), intent(out) :: Ls
    real, dimension(1:nx, 1:ny, 1:nz), intent(in) :: sIn

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction
    ! of the winds and density fluctuations
    character(len = *), intent(in) :: opt

    ! hortot = tot =>
    ! linear operator for total problem
    ! hortot = hor =>
    ! linear operator for horizontal problem
    ! hortot = hnd =>
    ! linear operator for horizontal problem without diagonal term
    character(len = *), intent(in) :: hortot

    ! local field (extended by ghost cells)
    real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1) :: s

    ! auxiliary fields for "dp"
    real, dimension(0:ny + 1, 0:nz + 1) :: xSliceLeft_send, xSliceRight_send
    real, dimension(0:ny + 1, 0:nz + 1) :: xSliceLeft_recv, xSliceRight_recv

    real, dimension(0:nx + 1, 0:nz + 1) :: ySliceBack_send, ySliceForw_send
    real, dimension(0:nx + 1, 0:nz + 1) :: ySliceBack_recv, ySliceForw_recv

    ! local variables
    integer :: i, j, k
    real :: AL, AR, AB, AF, AD, AU, AC, ACH, ACV, ALB, ALF, ARB, ARF
    real :: sL, sR, sB, sF, sD, sU, sC, sLB, sLF, sRB, sRF

    ! TFC
    real :: ARU, ARD, ALU, ALD, AFU, AFD, ABU, ABD
    real :: AUU, ADD, ARUU, ARDD, ALUU, ALDD, AFUU, AFDD, ABUU, ABDD
    real :: sRU, sRD, sLU, sLD, sFU, sFD, sBU, sBD
    real :: sUU, sDD, sRUU, sRDD, sLUU, sLDD, sFUU, sFDD, sBUU, sBDD

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    ! work with auxiliary field s
    s(1:nx, 1:ny, 1:nz) = sIn

    ! Find neighbour procs
    if(idim > 1) call mpi_cart_shift(comm, 0, 1, left, right, ierror)
    if(jdim > 1) call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

    select case(model)

      !----------------------------------------
      !       Pseudo-incompressible model
      !----------------------------------------

      ! TFC
      ! No changes for Boussinesq model required.
    case("pseudo_incompressible", "Boussinesq", "compressible")

      !----------------------------
      !   set Halo cells: xSlice
      !----------------------------

      if(xBoundary == "periodic") then
        if(idim > 1) then
          ! slice size
          sendcount = (ny + 2) * (nz + 2)
          recvcount = sendcount

          ! read slice into contiguous array
          xSliceLeft_send(:, :) = s(1, :, :)
          xSliceRight_send(:, :) = s(nx, :, :)

          ! left -> right
          source = left
          dest = right
          tag = 100

          call mpi_sendrecv(xSliceRight_send(0, 0), sendcount, &
              &mpi_double_precision, dest, tag, xSliceLeft_recv(0, 0), &
              &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              &sts_left, ierror)

          ! right -> left
          source = right
          dest = left
          tag = 100

          call mpi_sendrecv(xSliceLeft_send(0, 0), sendcount, &
              &mpi_double_precision, dest, tag, xSliceRight_recv(0, 0), &
              &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              &sts_right, ierror)

          ! right halos
          s(nx + 1, :, :) = xSliceRight_recv(:, :)

          ! left halos
          s(0, :, :) = xSliceLeft_recv(:, :)
        else
          s(0, :, :) = s(nx, :, :)
          s(nx + 1, :, :) = s(1, :, :)
        end if
      else
        stop "Poisson: unknown case xBoundary"
      endif

      if(verbose .and. master) print *, "horizontalHalos:  x-horizontal halos &
          &copied."

      !------------------------------
      !   set Halo cells: ySlice
      !------------------------------

      if(yBoundary == "periodic") then
        if(jdim > 1) then
          ! slice size
          sendcount = (nx + 2) * (nz + 2)
          recvcount = sendcount

          ! read slice into contiguous array
          ySliceBack_send(:, :) = s(:, 1, :)
          ySliceForw_send(:, :) = s(:, ny, :)

          ! back -> forw
          source = back
          dest = forw
          tag = 100

          call mpi_sendrecv(ySliceForw_send(0, 0), sendcount, &
              &mpi_double_precision, dest, tag, ySliceBack_recv(0, 0), &
              &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              &sts_back, ierror)

          ! forw -> back
          source = forw
          dest = back
          tag = 100

          call mpi_sendrecv(ySliceBack_send(0, 0), sendcount, &
              &mpi_double_precision, dest, tag, ySliceForw_recv(0, 0), &
              &recvcount, mpi_double_precision, source, mpi_any_tag, comm, &
              &sts_right, ierror)

          ! forward halos
          s(:, ny + 1, :) = ySliceForw_recv(:, :)

          ! backward halos
          s(:, 0, :) = ySliceBack_recv(:, :)
        else
          s(:, 0, :) = s(:, ny, :)
          s(:, ny + 1, :) = s(:, 1, :)
        end if
      else
        stop "Poisson: unknown case xBoundary"
      end if

      ! set vertical boundary conditions
      if(zBoundary == "periodic") then
        s(:, :, 0) = s(:, :, nz)
        s(:, :, nz + 1) = s(:, :, 1)
      end if

      if(verbose .and. master) print *, "horizontalHalos:  x-horizontal halos &
          &copied."

      !---------------------------------
      !         Loop over field
      !---------------------------------

      k_loop: do k = 1, nz
        j_loop: do j = 1, ny
          i_loop: do i = 1, nx

            ! ------------------ A(i+1,j,k) ------------------

            AR = ar_b(i, j, k)
            sR = s(i + 1, j, k)

            ! ------------------- A(i-1,j,k) --------------------

            AL = al_b(i, j, k)
            sL = s(i - 1, j, k)

            ! -------------------- A(i,j+1,k) ----------------------

            AF = af_b(i, j, k)
            sF = s(i, j + 1, k)

            ! --------------------- A(i,j-1,k) -----------------------

            AB = ab_b(i, j, k)
            sB = s(i, j - 1, k)

            ! --------------------- A(i,j,k+1) ------------------------

            ! TFC
            if(k < nz .or. zBoundary == "periodic") then
              AU = au_b(i, j, k)
              sU = s(i, j, k + 1)
            else ! k = nz -> upwad boundary (solid wall)
              ! A(i,j,nz+1) = 0
              AU = 0.0
              sU = 0.0
            end if

            ! --------------------- A(i,j,k-1) ------------------------

            ! TFC
            if(k > 1 .or. zBoundary == "periodic") then
              AD = ad_b(i, j, k)
              sD = s(i, j, k - 1)
            else ! k = 1 -> downward boundary (solid wall)
              ! A(i,j,0) = 0
              AD = 0.0
              sD = 0.0
            end if

            ! -------------------- A(i,j,k) --------------------------

            ACH = ach_b(i, j, k)
            ACV = acv_b(i, j, k)

            AC = ac_b(i, j, k)
            sC = s(i, j, k)

            ! -------------------- apply Operator ---------------------

            if(hortot == 'tot') then
              Ls(i, j, k) = AL * sL + AR * sR + AF * sF + AB * sB + AU * sU &
                  &+ AD * sD + AC * sC
            else if(hortot == 'hor') then
              Ls(i, j, k) = AL * sL + AR * sR + AF * sF + AB * sB + ACH * sC
            else if(hortot == 'hnd') then
              Ls(i, j, k) = AL * sL + AR * sR + AF * sF + AB * sB
            else
              stop "wrong hortot in linOpr"
            end if

            ! TFC
            ! Set additional matrix elements for TFC.
            if(topography) then
              ! ----------------- A(i+1,j,k+1) -----------------

              if(k < nz .or. zBoundary == "periodic") then
                ARU = aru_b(i, j, k)
                sRU = s(i + 1, j, k + 1)
              else
                ARU = 0.0
                sRU = 0.0
              end if

              ! ----------------- A(i+1,j,k-1) -----------------

              if(k > 1 .or. zBoundary == "periodic") then
                ARD = ard_b(i, j, k)
                sRD = s(i + 1, j, k - 1)
              else
                ARD = 0.0
                sRD = 0.0
              end if

              ! ----------------- A(i-1,j,k+1) -----------------

              if(k < nz .or. zBoundary == "periodic") then
                ALU = alu_b(i, j, k)
                sLU = s(i - 1, j, k + 1)
              else
                ALU = 0.0
                sLU = 0.0
              end if

              ! ----------------- A(i-1,j,k-1) -----------------

              if(k > 1 .or. zBoundary == "periodic") then
                ALD = ald_b(i, j, k)
                sLD = s(i - 1, j, k - 1)
              else
                ALD = 0.0
                sLD = 0.0
              end if

              ! ----------------- A(i,j+1,k+1) -----------------

              if(k < nz .or. zBoundary == "periodic") then
                AFU = afu_b(i, j, k)
                sFU = s(i, j + 1, k + 1)
              else
                AFU = 0.0
                sFU = 0.0
              end if

              ! ----------------- A(i,j+1,k-1) -----------------

              if(k > 1 .or. zBoundary == "periodic") then
                AFD = afd_b(i, j, k)
                sFD = s(i, j + 1, k - 1)
              else
                AFD = 0.0
                sFD = 0.0
              end if

              ! ----------------- A(i,j-1,k+1) -----------------

              if(k < nz .or. zBoundary == "periodic") then
                ABU = abu_b(i, j, k)
                sBU = s(i, j - 1, k + 1)
              else
                ABU = 0.0
                sBU = 0.0
              end if

              ! ----------------- A(i,j-1,k-1) -----------------

              if(k > 1 .or. zBoundary == "periodic") then
                ABD = abd_b(i, j, k)
                sBD = s(i, j - 1, k - 1)
              else
                ABD = 0.0
                sBD = 0.0
              end if

              ! ------------------ A(i,j,k+2) -----------------

              if(k < nz - 1 .or. zBoundary == "periodic") then
                AUU = auu_b(i, j, k)
                sUU = s(i, j, k + 2)
              else
                AUU = 0.0
                sUU = 0.0
              end if

              ! ------------------ A(i,j,k-2) -----------------

              if(k > 2 .or. zBoundary == "periodic") then
                ADD = add_b(i, j, k)
                sDD = s(i, j, k - 2)
              else
                ADD = 0.0
                sDD = 0.0
              end if

              ! ----------------- A(i+1,j,k+2) -----------------

              if(k < nz - 1 .or. zBoundary == "periodic") then
                ARUU = aruu_b(i, j, k)
                sRUU = s(i + 1, j, k + 2)
              else
                ARUU = 0.0
                sRUU = 0.0
              end if

              ! ----------------- A(i+1,j,k-2) -----------------

              if(k > 2 .or. zBoundary == "periodic") then
                ARDD = ardd_b(i, j, k)
                sRDD = s(i + 1, j, k - 2)
              else
                ARDD = 0.0
                sRDD = 0.0
              end if

              ! ----------------- A(i-1,j,k+2) -----------------

              if(k < nz - 1 .or. zBoundary == "periodic") then
                ALUU = aluu_b(i, j, k)
                sLUU = s(i - 1, j, k + 2)
              else
                ALUU = 0.0
                sLUU = 0.0
              end if

              ! ----------------- A(i-1,j,k-2) -----------------

              if(k > 2 .or. zBoundary == "periodic") then
                ALDD = aldd_b(i, j, k)
                sLDD = s(i - 1, j, k - 2)
              else
                ALDD = 0.0
                sLDD = 0.0
              end if

              ! ----------------- A(i,j+1,k+2) -----------------

              if(k < nz - 1 .or. zBoundary == "periodic") then
                AFUU = afuu_b(i, j, k)
                sFUU = s(i, j + 1, k + 2)
              else
                AFUU = 0.0
                sFUU = 0.0
              end if

              ! ----------------- A(i,j+1,k-2) -----------------

              if(k > 2 .or. zBoundary == "periodic") then
                AFDD = afdd_b(i, j, k)
                sFDD = s(i, j + 1, k - 2)
              else
                AFDD = 0.0
                sFDD = 0.0
              end if

              ! ----------------- A(i,j-1,k+2) -----------------

              if(k < nz - 1 .or. zBoundary == "periodic") then
                ABUU = abuu_b(i, j, k)
                sBUU = s(i, j - 1, k + 2)
              else
                ABUU = 0.0
                sBUU = 0.0
              end if

              ! ----------------- A(i,j-1,k-2) -----------------

              if(k > 2 .or. zBoundary == "periodic") then
                ABDD = abdd_b(i, j, k)
                sBDD = s(i, j - 1, k - 2)
              else
                ABDD = 0.0
                sBDD = 0.0
              end if

              ! Correct operator.
              Ls(i, j, k) = Ls(i, j, k) + ARU * sRU + ARD * sRD + ALU * sLU &
                  &+ ALD * sLD + AFU * sFU + AFD * sFD + ABU * sBU + ABD * sBD &
                  &+ AUU * sUU + ADD * sDD + ARUU * sRUU + ARDD * sRDD + ALUU &
                  &* sLUU + ALDD * sLDD + AFUU * sFUU + AFDD * sFDD + ABUU &
                  &* sBUU + ABDD * sBDD
            end if

            ! TFC FJ
            if(timeScheme == "semiimplicit" .and. .not. topography) then
              if(opt == 'impl') then
                ! -------------------- A(i,j,k) ---------------------

                ALB = alb_b(i, j, k)
                sLB = s(i - 1, j - 1, k)

                ! -------------------- A(i,j,k) ---------------------

                ALF = alf_b(i, j, k)
                sLF = s(i - 1, j + 1, k)

                ! -------------------- A(i,j,k) ---------------------

                ARB = arb_b(i, j, k)
                sRB = s(i + 1, j - 1, k)

                ! -------------------- A(i,j,k) ---------------------

                ARF = arf_b(i, j, k)
                sRF = s(i + 1, j + 1, k)

                Ls(i, j, k) = Ls(i, j, k) + ALB * sLB + ALF * sLF + ARB * sRB &
                    &+ ARF * sRF
              else if(opt /= 'expl') then
                stop 'ERROR: linOpr expects opt = expl or opt = impl'
              end if
            end if

            ! ---------------- scale with thetaStrat ------------------
            if(pressureScaling) then
              Ls(i, j, k) = Ls(i, j, k) / Pstrat(k)
            end if
          end do i_loop
        end do j_loop
      end do k_loop

    case default
      stop "linOpr: unknown case model"
    end select

  end subroutine linOpr

  !----------------------------------------------------------------------------

  subroutine linearOperatorTestTFC(var, dt, opt, facray)

    ! TFC
    ! Linear operator test for tensor elements in TFC.

    type(var_type), intent(in) :: var
    real, intent(in) :: dt, facray
    character(len = *), intent(in) :: opt

    integer :: i, j, k

    real :: jacInv
    real :: pEdgeRDiv, pEdgeLDiv, pEdgeFDiv, pEdgeBDiv, pEdgeUDiv, pEdgeDDiv
    real :: pEdgeRGra, pEdgeLGra, pEdgeFGra, pEdgeBGra, pEdgeUGra, pEdgeDGra
    real :: rhoEdgeR, rhoEdgeL, rhoEdgeF, rhoEdgeB, rhoEdgeU, rhoEdgeD

    real, dimension(0:(nx + 1), 0:(ny + 1), 0:(nz + 1)) :: s
    real, dimension(1:nx, 1:ny, 1:nz) :: sIn, Ls, LsTest

    real :: AR, AL, AF, AB, AU, AD, AC
    real :: ARU, ARD, ALU, ALD, AFU, AFD, ABU, ABD
    real :: AUU, ADD, ARUU, ARDD, ALUU, ALDD, AFUU, AFDD, ABUU, ABDD

    real :: sR, sL, sF, sB, sU, sD, sC
    real :: sRU, sRD, sLU, sLD, sFU, sFD, sBU, sBD
    real :: sUU, sDD, sRUU, sRDD, sLUU, sLDD, sFUU, sFDD, sBUU, sBDD

    real, dimension(1:nx, 1:ny, 1:nz) :: ac_tfc, ar_tfc, al_tfc, af_tfc, &
        &ab_tfc, au_tfc, ad_tfc, aru_tfc, ard_tfc, alu_tfc, ald_tfc, afu_tfc, &
        &afd_tfc, abu_tfc, abd_tfc, auu_tfc, add_tfc, aruu_tfc, ardd_tfc, &
        &aluu_tfc, aldd_tfc, afuu_tfc, afdd_tfc, abuu_tfc, abdd_tfc

    real, dimension(1:nx, 1:ny, 1:nz) :: ach_tfc

    real :: sumLoc, sumLocSquared, sumGlob, sumGlobSquared
    real * 4 :: sumOut1, sumOut2
    real * 4, dimension(1:nx, 1:ny, 1:nz) :: diffOut
    integer :: recTFC

    logical :: directTensorElements, directElementAssignment
    character(len = 50) :: linearOperatorTestCase

    preconditioner = "yes"

    directTensorElements = .false.
    directElementAssignment = .false.
    linearOperatorTestCase = "noTopography"

    s = var%pi(0:(nx + 1), 0:(ny + 1), 0:(nz + 1))
    sIn = var%pi(1:nx, 1:ny, 1:nz)

    if(directTensorElements) then

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx

            ! Compute inverse Jacobian.
            jacInv = 1.0 / jac(i, j, k)

            ! Compute P coefficients (divergence).
            pEdgeRDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, &
                &j, k) * pStratTFC(i + 1, j, k))
            pEdgeLDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i - 1, &
                &j, k) * pStratTFC(i - 1, j, k))
            pEdgeFDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                &+ 1, k) * pStratTFC(i, j + 1, k))
            pEdgeBDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                &- 1, k) * pStratTFC(i, j - 1, k))
            pEdgeUDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, k &
                &+ 1) * pStratTFC(i, j, k + 1))
            pEdgeDDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, k &
                &- 1) * pStratTFC(i, j, k - 1))

            ! Compute P coefficients (pressure gradient).
            pEdgeRGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) + pStratTFC(i &
                &+ 1, j, k) / jac(i + 1, j, k))
            pEdgeLGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) + pStratTFC(i &
                &- 1, j, k) / jac(i - 1, j, k))
            pEdgeFGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                &+ pStratTFC(i, j + 1, k) / jac(i, j + 1, k))
            pEdgeBGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                &+ pStratTFC(i, j - 1, k) / jac(i, j - 1, k))
            pEdgeUGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                &+ pStratTFC(i, j, k + 1) / jac(i, j, k + 1))
            pEdgeDGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                &+ pStratTFC(i, j, k - 1) / jac(i, j, k - 1))

            ! Compute rho coefficients.
            rhoEdgeR = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, j, k))
            rhoEdgeL = 0.5 * (var%rho(i, j, k) + var%rho(i - 1, j, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i - 1, j, k))
            rhoEdgeF = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j + 1, k))
            rhoEdgeB = 0.5 * (var%rho(i, j, k) + var%rho(i, j - 1, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j - 1, k))
            rhoEdgeU = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j, k + 1))
            rhoEdgeD = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k - 1) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j, k - 1))

            ! --------------------- A(i,j,k) ---------------------

            if(k == 1) then
              AC = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * (1.0 &
                  &/ dx + 0.75 * met(i, j, k, 1, 3) / dz) + pEdgeLDiv &
                  &/ rhoEdgeL * pEdgeLGra * (1.0 / dx - 0.75 * met(i, j, k, 1, &
                  &3) / dz)) - jacInv / dy * (pEdgeFDiv / rhoEdgeF * pEdgeFGra &
                  &* (1.0 / dy + 0.75 * met(i, j, k, 2, 3) / dz) + pEdgeBDiv &
                  &/ rhoEdgeB * pEdgeBGra * (1.0 / dy - 0.75 * met(i, j, k, 2, &
                  &3) / dz)) - jacInv / dz * pEdgeUDiv / rhoEdgeU * pEdgeUGra &
                  &* met(i, j, k, 3, 3) / dz + jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * 0.5 * pStratTFC(i, j, k) / jac(i, j, k) &
                  &* (chris(i, j, k, 1, 1) + chris(i, j, k, 2, 2) + 2.0 &
                  &* chris(i, j, k, 1, 3) * met(i, j, k, 1, 3) + 2.0 &
                  &* chris(i, j, k, 2, 3) * met(i, j, k, 2, 3))
            else if(k == nz) then
              AC = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * (1.0 &
                  &/ dx - 0.75 * met(i, j, k, 1, 3) / dz) + pEdgeLDiv &
                  &/ rhoEdgeL * pEdgeLGra * (1.0 / dx + 0.75 * met(i, j, k, 1, &
                  &3) / dz)) - jacInv / dy * (pEdgeFDiv / rhoEdgeF * pEdgeFGra &
                  &* (1.0 / dy - 0.75 * met(i, j, k, 2, 3) / dz) + pEdgeBDiv &
                  &/ rhoEdgeB * pEdgeBGra * (1.0 / dy + 0.75 * met(i, j, k, 2, &
                  &3) / dz)) - jacInv / dz * pEdgeDDiv / rhoEdgeD * pEdgeDGra &
                  &* met(i, j, k, 3, 3) / dz - jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * 0.5 * pStratTFC(i, j, k) / jac(i, j, k) &
                  &* (chris(i, j, k, 1, 1) + chris(i, j, k, 2, 2) + 2.0 &
                  &* chris(i, j, k, 1, 3) * met(i, j, k, 1, 3) + 2.0 &
                  &* chris(i, j, k, 2, 3) * met(i, j, k, 2, 3))
            else
              AC = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra / dx &
                  &+ pEdgeLDiv / rhoEdgeL * pEdgeLGra / dx) - jacInv / dy &
                  &* (pEdgeFDiv / rhoEdgeF * pEdgeFGra / dy + pEdgeBDiv &
                  &/ rhoEdgeB * pEdgeBGra / dy) - jacInv / dz * (pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met(i, j, k, 3, 3) / dz &
                  &+ pEdgeDDiv / rhoEdgeD * pEdgeDGra * met(i, j, k, 3, 3) &
                  &/ dz) + jacInv / dz * pEdgeUDiv / rhoEdgeU * 0.5 &
                  &* pStratTFC(i, j, k) / jac(i, j, k) * (chris(i, j, k, 1, 1) &
                  &+ chris(i, j, k, 2, 2) + 2.0 * chris(i, j, k, 1, 3) &
                  &* met(i, j, k, 1, 3) + 2.0 * chris(i, j, k, 2, 3) * met(i, &
                  &j, k, 2, 3)) - jacInv / dz * pEdgeDDiv / rhoEdgeD * 0.5 &
                  &* pStratTFC(i, j, k) / jac(i, j, k) * (chris(i, j, k, 1, 1) &
                  &+ chris(i, j, k, 2, 2) + 2.0 * chris(i, j, k, 1, 3) &
                  &* met(i, j, k, 1, 3) + 2.0 * chris(i, j, k, 2, 3) * met(i, &
                  &j, k, 2, 3))
            end if

            ! -------------------- A(i+1,j,k) --------------------

            if(k == 1) then
              AR = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * (1.0 / dx &
                  &- 0.75 * met(i + 1, j, k, 1, 3) / dz) + jacInv / dz &
                  &* pEdgeUDiv / rhoEdgeU * pEdgeUGra * 0.25 * met(i + 1, j, &
                  &k, 1, 3) / dx
            else if(k == nz) then
              AR = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * (1.0 / dx &
                  &+ 0.75 * met(i + 1, j, k, 1, 3) / dz) - jacInv / dz &
                  &* pEdgeDDiv / rhoEdgeD * pEdgeDGra * 0.25 * met(i + 1, j, &
                  &k, 1, 3) / dx
            else
              AR = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra / dx &
                  &+ jacInv / dz * (pEdgeUDiv / rhoEdgeU * pEdgeUGra * met(i &
                  &+ 1, j, k, 1, 3) * 0.25 / dx - pEdgeDDiv / rhoEdgeD &
                  &* pEdgeDGra * met(i + 1, j, k, 1, 3) * 0.25 / dx)
            end if

            ! -------------------- A(i-1,j,k) --------------------

            if(k == 1) then
              AL = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * (1.0 / dx &
                  &+ 0.75 * met(i - 1, j, k, 1, 3) / dz) - jacInv / dz &
                  &* pEdgeUDiv / rhoEdgeU * pEdgeUGra * 0.25 * met(i - 1, j, &
                  &k, 1, 3) / dx
            else if(k == nz) then
              AL = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * (1.0 / dx &
                  &- 0.75 * met(i - 1, j, k, 1, 3) / dz) + jacInv / dz &
                  &* pEdgeDDiv / rhoEdgeD * pEdgeDGra * 0.25 * met(i - 1, j, &
                  &k, 1, 3) / dx
            else
              AL = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra / dx &
                  &- jacInv / dz * (pEdgeUDiv / rhoEdgeU * pEdgeUGra * met(i &
                  &- 1, j, k, 1, 3) * 0.25 / dx - pEdgeDDiv / rhoEdgeD &
                  &* pEdgeDGra * met(i - 1, j, k, 1, 3) * 0.25 / dx)
            end if

            ! -------------------- A(i,j+1,k) --------------------

            if(k == 1) then
              AF = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * (1.0 / dy &
                  &- 0.75 * met(i, j + 1, k, 2, 3) / dz) + jacInv / dz &
                  &* pEdgeUDiv / rhoEdgeU * pEdgeUGra * 0.25 * met(i, j + 1, &
                  &k, 2, 3) / dy
            else if(k == nz) then
              AF = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * (1.0 / dy &
                  &+ 0.75 * met(i, j + 1, k, 2, 3) / dz) - jacInv / dz &
                  &* pEdgeDDiv / rhoEdgeD * pEdgeDGra * 0.25 * met(i, j + 1, &
                  &k, 2, 3) / dy
            else
              AF = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra / dy &
                  &+ jacInv / dz * (pEdgeUDiv / rhoEdgeU * pEdgeUGra * met(i, &
                  &j + 1, k, 2, 3) * 0.25 / dy - pEdgeDDiv / rhoEdgeD &
                  &* pEdgeDGra * met(i, j + 1, k, 2, 3) * 0.25 / dy)
            end if

            ! -------------------- A(i,j-1,k) --------------------

            if(k == 1) then
              AB = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * (1.0 / dy &
                  &+ 0.75 * met(i, j - 1, k, 2, 3) / dz) - jacInv / dz &
                  &* pEdgeUDiv / rhoEdgeU * pEdgeUGra * 0.25 * met(i, j - 1, &
                  &k, 2, 3) / dy
            else if(k == nz) then
              AB = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * (1.0 / dy &
                  &- 0.75 * met(i, j - 1, k, 2, 3) / dz) + jacInv / dz &
                  &* pEdgeDDiv / rhoEdgeD * pEdgeDGra * 0.25 * met(i, j - 1, &
                  &k, 2, 3) / dy
            else
              AB = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra / dy &
                  &- jacInv / dz * (pEdgeUDiv / rhoEdgeU * pEdgeUGra * met(i, &
                  &j - 1, k, 2, 3) * 0.25 / dy - pEdgeDDiv / rhoEdgeD &
                  &* pEdgeDGra * met(i, j - 1, k, 2, 3) * 0.25 / dy)
            end if

            ! -------------------- A(i,j,k+1) --------------------

            if(k == 1) then
              AU = jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * met(i, j, &
                  &k + 1, 1, 3) / dz - pEdgeLDiv / rhoEdgeL * pEdgeLGra &
                  &* met(i, j, k + 1, 1, 3) / dz) + jacInv / dy * (pEdgeFDiv &
                  &/ rhoEdgeF * pEdgeFGra * met(i, j, k + 1, 2, 3) / dz &
                  &- pEdgeBDiv / rhoEdgeB * pEdgeBGra * met(i, j, k + 1, 2, 3) &
                  &/ dz) + jacInv / dz * pEdgeUDiv / rhoEdgeU * pEdgeUGra &
                  &* met(i, j, k + 1, 3, 3) / dz + jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * 0.5 * pStratTFC(i, j, k + 1) / jac(i, j, k &
                  &+ 1) * (chris(i, j, k + 1, 1, 1) + chris(i, j, k + 1, 2, 2) &
                  &+ 2.0 * chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, 1, 3) &
                  &+ 2.0 * chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, 2, 3))
            else if(k == nz) then
              AU = 0.0
            else
              AU = jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * met(i, j, &
                  &k + 1, 1, 3) * 0.25 / dz - pEdgeLDiv / rhoEdgeL * pEdgeLGra &
                  &* met(i, j, k + 1, 1, 3) * 0.25 / dz) + jacInv / dy &
                  &* (pEdgeFDiv / rhoEdgeF * pEdgeFGra * met(i, j, k + 1, 2, &
                  &3) * 0.25 / dz - pEdgeBDiv / rhoEdgeB * pEdgeBGra * met(i, &
                  &j, k + 1, 2, 3) * 0.25 / dz) + jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met(i, j, k + 1, 3, 3) / dz &
                  &+ jacInv / dz * pEdgeUDiv / rhoEdgeU * 0.5 * pStratTFC(i, &
                  &j, k + 1) / jac(i, j, k + 1) * (chris(i, j, k + 1, 1, 1) &
                  &+ chris(i, j, k + 1, 2, 2) + 2.0 * chris(i, j, k + 1, 1, 3) &
                  &* met(i, j, k + 1, 1, 3) + 2.0 * chris(i, j, k + 1, 2, 3) &
                  &* met(i, j, k + 1, 2, 3))
            end if

            ! -------------------- A(i,j,k-1) --------------------

            if(k == 1) then
              AD = 0.0
            else if(k == nz) then
              AD = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * met(i, &
                  &j, k - 1, 1, 3) / dz - pEdgeLDiv / rhoEdgeL * pEdgeLGra &
                  &* met(i, j, k - 1, 1, 3) / dz) - jacInv / dy * (pEdgeFDiv &
                  &/ rhoEdgeF * pEdgeFGra * met(i, j, k - 1, 2, 3) / dz &
                  &- pEdgeBDiv / rhoEdgeB * pEdgeBGra * met(i, j, k - 1, 2, 3) &
                  &/ dz) + jacInv / dz * pEdgeDDiv / rhoEdgeD * pEdgeDGra &
                  &* met(i, j, k - 1, 3, 3) / dz - jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * 0.5 * pStratTFC(i, j, k - 1) / jac(i, j, k &
                  &- 1) * (chris(i, j, k - 1, 1, 1) + chris(i, j, k - 1, 2, 2) &
                  &+ 2.0 * chris(i, j, k - 1, 1, 3) * met(i, j, k - 1, 1, 3) &
                  &+ 2.0 * chris(i, j, k - 1, 2, 3) * met(i, j, k - 1, 2, 3))
            else
              AD = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * met(i, &
                  &j, k - 1, 1, 3) * 0.25 / dz - pEdgeLDiv / rhoEdgeL &
                  &* pEdgeLGra * met(i, j, k - 1, 1, 3) * 0.25 / dz) - jacInv &
                  &/ dy * (pEdgeFDiv / rhoEdgeF * pEdgeFGra * met(i, j, k - 1, &
                  &2, 3) * 0.25 / dz - pEdgeBDiv / rhoEdgeB * pEdgeBGra &
                  &* met(i, j, k - 1, 2, 3) * 0.25 / dz) + jacInv / dz &
                  &* pEdgeDDiv / rhoEdgeD * pEdgeDGra * met(i, j, k - 1, 3, 3) &
                  &/ dz - jacInv / dz * pEdgeDDiv / rhoEdgeD * 0.5 &
                  &* pStratTFC(i, j, k - 1) / jac(i, j, k - 1) * (chris(i, j, &
                  &k - 1, 1, 1) + chris(i, j, k - 1, 2, 2) + 2.0 * chris(i, j, &
                  &k - 1, 1, 3) * met(i, j, k - 1, 1, 3) + 2.0 * chris(i, j, k &
                  &- 1, 2, 3) * met(i, j, k - 1, 2, 3))
            end if

            ! ------------------- A(i+1,j,k+1) -------------------

            if(k == 1) then
              ARU = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * met(i &
                  &+ 1, j, k + 1, 1, 3) / dz + jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met(i + 1, j, k + 1, 1, 3) * 0.25 &
                  &/ dx
            else if(k == nz) then
              ARU = 0.0
            else
              ARU = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * met(i &
                  &+ 1, j, k + 1, 1, 3) * 0.25 / dz + jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met(i + 1, j, k + 1, 1, 3) * 0.25 &
                  &/ dx
            end if

            ! ------------------- A(i+1,j,k-1) -------------------

            if(k == 1) then
              ARD = 0.0
            else if(k == nz) then
              ARD = - jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * met(i &
                  &+ 1, j, k - 1, 1, 3) / dz - jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * met(i + 1, j, k - 1, 1, 3) * 0.25 &
                  &/ dx
            else
              ARD = - jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * met(i &
                  &+ 1, j, k - 1, 1, 3) * 0.25 / dz - jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * met(i + 1, j, k - 1, 1, 3) * 0.25 &
                  &/ dx
            end if

            ! ------------------- A(i-1,j,k+1) -------------------

            if(k == 1) then
              ALU = - jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * met(i &
                  &- 1, j, k + 1, 1, 3) / dz - jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met(i - 1, j, k + 1, 1, 3) * 0.25 &
                  &/ dx
            else if(k == nz) then
              ALU = 0.0
            else
              ALU = - jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * met(i &
                  &- 1, j, k + 1, 1, 3) * 0.25 / dz - jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met(i - 1, j, k + 1, 1, 3) * 0.25 &
                  &/ dx
            end if

            ! ------------------- A(i-1,j,k-1) -------------------

            if(k == 1) then
              ALD = 0.0
            else if(k == nz) then
              ALD = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * met(i &
                  &- 1, j, k - 1, 1, 3) / dz + jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * met(i - 1, j, k - 1, 1, 3) * 0.25 &
                  &/ dx
            else
              ALD = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * met(i &
                  &- 1, j, k - 1, 1, 3) * 0.25 / dz + jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * met(i - 1, j, k - 1, 1, 3) * 0.25 &
                  &/ dx
            end if

            ! ------------------- A(i,j+1,k+1) -------------------

            if(k == 1) then
              AFU = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * met(i, j &
                  &+ 1, k + 1, 2, 3) / dz + jacInv / dz * pEdgeUDiv / rhoEdgeU &
                  &* pEdgeUGra * met(i, j + 1, k + 1, 2, 3) * 0.25 / dy
            else if(k == nz) then
              AFU = 0.0
            else
              AFU = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * met(i, j &
                  &+ 1, k + 1, 2, 3) * 0.25 / dz + jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met(i, j + 1, k + 1, 2, 3) * 0.25 &
                  &/ dy
            end if

            ! ------------------- A(i,j+1,k-1) -------------------

            if(k == 1) then
              AFD = 0.0
            else if(k == nz) then
              AFD = - jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * met(i, &
                  &j + 1, k - 1, 2, 3) / dz - jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * met(i, j + 1, k - 1, 2, 3) * 0.25 &
                  &/ dy
            else
              AFD = - jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * met(i, &
                  &j + 1, k - 1, 2, 3) * 0.25 / dz - jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * met(i, j + 1, k - 1, 2, 3) * 0.25 &
                  &/ dy
            end if

            ! ------------------- A(i,j-1,k+1) -------------------

            if(k == 1) then
              ABU = - jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * met(i, &
                  &j - 1, k + 1, 2, 3) / dz - jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met(i, j - 1, k + 1, 2, 3) * 0.25 &
                  &/ dy
            else if(k == nz) then
              ABU = 0.0
            else
              ABU = - jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * met(i, &
                  &j - 1, k + 1, 2, 3) * 0.25 / dz - jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met(i, j - 1, k + 1, 2, 3) * 0.25 &
                  &/ dy
            end if

            ! ------------------- A(i,j-1,k-1) -------------------

            if(k == 1) then
              ABD = 0.0
            else if(k == nz) then
              ABD = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * met(i, j &
                  &- 1, k - 1, 2, 3) / dz + jacInv / dz * pEdgeDDiv / rhoEdgeD &
                  &* pEdgeDGra * met(i, j - 1, k - 1, 2, 3) * 0.25 / dy
            else
              ABD = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * met(i, j &
                  &- 1, k - 1, 2, 3) * 0.25 / dz + jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * met(i, j - 1, k - 1, 2, 3) * 0.25 &
                  &/ dy
            end if

            ! ------------------- A(i,j,k+2) ---------------------

            if(k == 1) then
              AUU = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * 0.25 &
                  &* met(i, j, k + 2, 1, 3) / dz - pEdgeLDiv / rhoEdgeL &
                  &* pEdgeLGra * 0.25 * met(i, j, k + 2, 1, 3) / dz) - jacInv &
                  &/ dy * (pEdgeFDiv / rhoEdgeF * pEdgeFGra * 0.25 * met(i, j, &
                  &k + 2, 2, 3) / dz - pEdgeBDiv / rhoEdgeB * pEdgeBGra * 0.25 &
                  &* met(i, j, k + 2, 2, 3) / dz)
            else
              AUU = 0.0
            end if

            ! ------------------- A(i,j,k-2) ---------------------

            if(k == nz) then
              ADD = jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * 0.25 &
                  &* met(i, j, k - 2, 1, 3) / dz - pEdgeLDiv / rhoEdgeL &
                  &* pEdgeLGra * 0.25 * met(i, j, k - 2, 1, 3) / dz) + jacInv &
                  &/ dy * (pEdgeFDiv / rhoEdgeF * pEdgeFGra * 0.25 * met(i, j, &
                  &k - 2, 2, 3) / dz - pEdgeBDiv / rhoEdgeB * pEdgeBGra * 0.25 &
                  &* met(i, j, k - 2, 2, 3) / dz)
            else
              ADD = 0.0
            end if

            ! ------------------ A(i+1,j,k+2) --------------------

            if(k == 1) then
              ARUU = - jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * 0.25 &
                  &* met(i + 1, j, k + 2, 1, 3) / dz
            else
              ARUU = 0.0
            end if

            ! ------------------ A(i+1,j,k-2) --------------------

            if(k == nz) then
              ARDD = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * 0.25 &
                  &* met(i + 1, j, k - 2, 1, 3) / dz
            else
              ARDD = 0.0
            end if

            ! ------------------ A(i-1,j,k+2) --------------------

            if(k == 1) then
              ALUU = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * 0.25 &
                  &* met(i - 1, j, k + 2, 1, 3) / dz
            else
              ALUU = 0.0
            end if

            ! ------------------ A(i-1,j,k-2) --------------------

            if(k == nz) then
              ALDD = - jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * 0.25 &
                  &* met(i - 1, j, k - 2, 1, 3) / dz
            else
              ALDD = 0.0
            end if

            ! ------------------ A(i,j+1,k+2) --------------------

            if(k == 1) then
              AFUU = - jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * 0.25 &
                  &* met(i, j + 1, k + 2, 2, 3) / dz
            else
              AFUU = 0.0
            end if

            ! ------------------ A(i,j+1,k-2) --------------------

            if(k == nz) then
              AFDD = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * 0.25 &
                  &* met(i, j + 1, k - 2, 2, 3) / dz
            else
              AFDD = 0.0
            end if

            ! ------------------ A(i,j-1,k+2) --------------------

            if(k == 1) then
              ABUU = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * 0.25 &
                  &* met(i, j - 1, k + 2, 2, 3) / dz
            else
              ABUU = 0.0
            end if

            ! ------------------ A(i,j-1,k-2) --------------------

            if(k == nz) then
              ABDD = - jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * 0.25 &
                  &* met(i, j - 1, k - 2, 2, 3) / dz
            else
              ABDD = 0.0
            end if

            ! Set matrix elements for bicgstab.
            ac_b(i, j, k) = AC
            ar_b(i, j, k) = AR
            al_b(i, j, k) = AL
            af_b(i, j, k) = AF
            ab_b(i, j, k) = AB
            au_b(i, j, k) = AU
            ad_b(i, j, k) = AD
            aru_b(i, j, k) = ARU
            ard_b(i, j, k) = ARD
            alu_b(i, j, k) = ALU
            ald_b(i, j, k) = ALD
            afu_b(i, j, k) = AFU
            afd_b(i, j, k) = AFD
            abu_b(i, j, k) = ABU
            abd_b(i, j, k) = ABD
            auu_b(i, j, k) = AUU
            add_b(i, j, k) = ADD
            aruu_b(i, j, k) = ARUU
            ardd_b(i, j, k) = ARDD
            aluu_b(i, j, k) = ALUU
            aldd_b(i, j, k) = ALDD
            afuu_b(i, j, k) = AFUU
            afdd_b(i, j, k) = AFDD
            abuu_b(i, j, k) = ABUU
            abdd_b(i, j, k) = ABDD

          end do
        end do
      end do

    else

      call val_PsIn(var, dt, opt, facray)

    end if

    if(directElementAssignment) then

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx

            ! ------------------ A(i+1,j,k) ------------------

            AR = ar_b(i, j, k)
            sR = s(i + 1, j, k)

            ! ------------------- A(i-1,j,k) --------------------

            AL = al_b(i, j, k)
            sL = s(i - 1, j, k)

            ! -------------------- A(i,j+1,k) ----------------------

            AF = af_b(i, j, k)
            sF = s(i, j + 1, k)

            ! --------------------- A(i,j-1,k) -----------------------

            AB = ab_b(i, j, k)
            sB = s(i, j - 1, k)

            ! --------------------- A(i,j,k+1) ------------------------

            if(k < nz) then
              AU = au_b(i, j, k)
              sU = s(i, j, k + 1)
            else ! k = nz -> upwad boundary (solid wall)
              ! A(i,j,nz+1) = 0
              AU = 0.0
              sU = 0.0
            end if

            ! --------------------- A(i,j,k-1) ------------------------

            if(k > 1) then
              AD = ad_b(i, j, k)
              sD = s(i, j, k - 1)
            else ! k = 1 -> downward boundary (solid wall)
              ! A(i,j,0) = 0
              AD = 0.0
              sD = 0.0
            end if

            ! -------------------- A(i,j,k) --------------------------

            AC = ac_b(i, j, k)
            sC = s(i, j, k)

            ! -------------------- Apply operator ---------------------

            LsTest(i, j, k) = AL * sL + AR * sR + AF * sF + AB * sB + AU * sU &
                &+ AD * sD + AC * sC

            ! ----------------- A(i+1,j,k+1) -----------------

            if(k < nz) then
              ARU = aru_b(i, j, k)
              sRU = s(i + 1, j, k + 1)
            else
              ARU = 0.0
              sRU = 0.0
            end if

            ! ----------------- A(i+1,j,k-1) -----------------

            if(k > 1) then
              ARD = ard_b(i, j, k)
              sRD = s(i + 1, j, k - 1)
            else
              ARD = 0.0
              sRD = 0.0
            end if

            ! ----------------- A(i-1,j,k+1) -----------------

            if(k < nz) then
              ALU = alu_b(i, j, k)
              sLU = s(i - 1, j, k + 1)
            else
              ALU = 0.0
              sLU = 0.0
            end if

            ! ----------------- A(i-1,j,k-1) -----------------

            if(k > 1) then
              ALD = ald_b(i, j, k)
              sLD = s(i - 1, j, k - 1)
            else
              ALD = 0.0
              sLD = 0.0
            end if

            ! ----------------- A(i,j+1,k+1) -----------------

            if(k < nz) then
              AFU = afu_b(i, j, k)
              sFU = s(i, j + 1, k + 1)
            else
              AFU = 0.0
              sFU = 0.0
            end if

            ! ----------------- A(i,j+1,k-1) -----------------

            if(k > 1) then
              AFD = afd_b(i, j, k)
              sFD = s(i, j + 1, k - 1)
            else
              AFD = 0.0
              sFD = 0.0
            end if

            ! ----------------- A(i,j-1,k+1) -----------------

            if(k < nz) then
              ABU = abu_b(i, j, k)
              sBU = s(i, j - 1, k + 1)
            else
              ABU = 0.0
              sBU = 0.0
            end if

            ! ----------------- A(i,j-1,k-1) -----------------

            if(k > 1) then
              ABD = abd_b(i, j, k)
              sBD = s(i, j - 1, k - 1)
            else
              ABD = 0.0
              sBD = 0.0
            end if

            ! ------------------ A(i,j,k+2) -----------------

            if(k == 1) then
              AUU = auu_b(i, j, k)
              sUU = s(i, j, k + 2)
            else
              AUU = 0.0
              sUU = 0.0
            end if

            ! ------------------ A(i,j,k-2) -----------------

            if(k == nz) then
              ADD = add_b(i, j, k)
              sDD = s(i, j, k - 2)
            else
              ADD = 0.0
              sDD = 0.0
            end if

            ! ----------------- A(i+1,j,k+2) -----------------

            if(k == 1) then
              ARUU = aruu_b(i, j, k)
              sRUU = s(i + 1, j, k + 2)
            else
              ARUU = 0.0
              sRUU = 0.0
            end if

            ! ----------------- A(i+1,j,k-2) -----------------

            if(k == nz) then
              ARDD = ardd_b(i, j, k)
              sRDD = s(i + 1, j, k - 2)
            else
              ARDD = 0.0
              sRDD = 0.0
            end if

            ! ----------------- A(i-1,j,k+2) -----------------

            if(k == 1) then
              ALUU = aluu_b(i, j, k)
              sLUU = s(i - 1, j, k + 2)
            else
              ALUU = 0.0
              sLUU = 0.0
            end if

            ! ----------------- A(i-1,j,k-2) -----------------

            if(k == nz) then
              ALDD = aldd_b(i, j, k)
              sLDD = s(i - 1, j, k - 2)
            else
              ALDD = 0.0
              sLDD = 0.0
            end if

            ! ----------------- A(i,j+1,k+2) -----------------

            if(k == 1) then
              AFUU = afuu_b(i, j, k)
              sFUU = s(i, j + 1, k + 2)
            else
              AFUU = 0.0
              sFUU = 0.0
            end if

            ! ----------------- A(i,j+1,k-2) -----------------

            if(k == nz) then
              AFDD = afdd_b(i, j, k)
              sFDD = s(i, j + 1, k - 2)
            else
              AFDD = 0.0
              sFDD = 0.0
            end if

            ! ----------------- A(i,j-1,k+2) -----------------

            if(k == 1) then
              ABUU = abuu_b(i, j, k)
              sBUU = s(i, j - 1, k + 2)
            else
              ABUU = 0.0
              sBUU = 0.0
            end if

            ! ----------------- A(i,j-1,k-2) -----------------

            if(k == nz) then
              ABDD = abdd_b(i, j, k)
              sBDD = s(i, j - 1, k - 2)
            else
              ABDD = 0.0
              sBDD = 0.0
            end if

            ! Correct operator.
            LsTest(i, j, k) = LsTest(i, j, k) + ARU * sRU + ARD * sRD + ALU &
                &* sLU + ALD * sLD + AFU * sFU + AFD * sFD + ABU * sBU + ABD &
                &* sBD + AUU * sUU + ADD * sDD + ARUU * sRUU + ARDD * sRDD &
                &+ ALUU * sLUU + ALDD * sLDD + AFUU * sFUU + AFDD * sFDD &
                &+ ABUU * sBUU + ABDD * sBDD

          end do
        end do
      end do

    else

      call linOpr(sIn, LsTest, opt, "tot")

    end if

    select case(linearOperatorTestCase)

    case("directLinearOperator")

      k_loop: do k = 1, nz
        j_loop: do j = 1, ny
          i_loop: do i = 1, nx

            ! Compute inverse Jacobian.
            jacInv = 1.0 / jac(i, j, k)

            ! Compute P coefficients (divergence).
            pEdgeRDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, &
                &j, k) * pStratTFC(i + 1, j, k))
            pEdgeLDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i - 1, &
                &j, k) * pStratTFC(i - 1, j, k))
            pEdgeFDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                &+ 1, k) * pStratTFC(i, j + 1, k))
            pEdgeBDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                &- 1, k) * pStratTFC(i, j - 1, k))
            pEdgeUDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, k &
                &+ 1) * pStratTFC(i, j, k + 1))
            pEdgeDDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, k &
                &- 1) * pStratTFC(i, j, k - 1))

            ! Compute P coefficients (pressure gradient).
            pEdgeRGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) + pStratTFC(i &
                &+ 1, j, k) / jac(i + 1, j, k))
            pEdgeLGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) + pStratTFC(i &
                &- 1, j, k) / jac(i - 1, j, k))
            pEdgeFGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                &+ pStratTFC(i, j + 1, k) / jac(i, j + 1, k))
            pEdgeBGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                &+ pStratTFC(i, j - 1, k) / jac(i, j - 1, k))
            pEdgeUGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                &+ pStratTFC(i, j, k + 1) / jac(i, j, k + 1))
            pEdgeDGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                &+ pStratTFC(i, j, k - 1) / jac(i, j, k - 1))

            ! Compute rho coefficients.
            rhoEdgeR = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, j, k))
            rhoEdgeL = 0.5 * (var%rho(i, j, k) + var%rho(i - 1, j, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i - 1, j, k))
            rhoEdgeF = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j + 1, k))
            rhoEdgeB = 0.5 * (var%rho(i, j, k) + var%rho(i, j - 1, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j - 1, k))
            rhoEdgeU = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j, k + 1))
            rhoEdgeD = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k - 1) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j, k - 1))

            ! Direct computation.

            if(k == 1) then

              Ls(i, j, k) = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra &
                  &* ((s(i + 1, j, k) - s(i, j, k)) / dx + 0.25 * (- met(i, j, &
                  &k + 2, 1, 3) * s(i, j, k + 2) - met(i + 1, j, k + 2, 1, 3) &
                  &* s(i + 1, j, k + 2) + 4.0 * met(i, j, k + 1, 1, 3) * s(i, &
                  &j, k + 1) + 4.0 * met(i + 1, j, k + 1, 1, 3) * s(i + 1, j, &
                  &k + 1) - 3.0 * met(i, j, k, 1, 3) * s(i, j, k) - 3.0 &
                  &* met(i + 1, j, k, 1, 3) * s(i + 1, j, k)) / dz) - jacInv &
                  &/ dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * ((s(i, j, k) &
                  &- s(i - 1, j, k)) / dx + 0.25 * (- met(i, j, k + 2, 1, 3) &
                  &* s(i, j, k + 2) - met(i - 1, j, k + 2, 1, 3) * s(i - 1, j, &
                  &k + 2) + 4.0 * met(i, j, k + 1, 1, 3) * s(i, j, k + 1) &
                  &+ 4.0 * met(i - 1, j, k + 1, 1, 3) * s(i - 1, j, k + 1) &
                  &- 3.0 * met(i, j, k, 1, 3) * s(i, j, k) - 3.0 * met(i - 1, &
                  &j, k, 1, 3) * s(i - 1, j, k)) / dz) + jacInv / dy &
                  &* pEdgeFDiv / rhoEdgeF * pEdgeFGra * ((s(i, j + 1, k) &
                  &- s(i, j, k)) / dy + 0.25 * (- met(i, j, k + 2, 2, 3) &
                  &* s(i, j, k + 2) - met(i, j + 1, k + 2, 2, 3) * s(i, j + 1, &
                  &k + 2) + 4.0 * met(i, j, k + 1, 2, 3) * s(i, j, k + 1) &
                  &+ 4.0 * met(i, j + 1, k + 1, 2, 3) * s(i, j + 1, k + 1) &
                  &- 3.0 * met(i, j, k, 2, 3) * s(i, j, k) - 3.0 * met(i, j &
                  &+ 1, k, 2, 3) * s(i, j + 1, k)) / dz) - jacInv / dy &
                  &* pEdgeBDiv / rhoEdgeB * pEdgeBGra * ((s(i, j, k) - s(i, j &
                  &- 1, k)) / dy + 0.25 * (- met(i, j, k + 2, 2, 3) * s(i, j, &
                  &k + 2) - met(i, j - 1, k + 2, 2, 3) * s(i, j - 1, k + 2) &
                  &+ 4.0 * met(i, j, k + 1, 2, 3) * s(i, j, k + 1) + 4.0 &
                  &* met(i, j - 1, k + 1, 2, 3) * s(i, j - 1, k + 1) - 3.0 &
                  &* met(i, j, k, 2, 3) * s(i, j, k) - 3.0 * met(i, j - 1, k, &
                  &2, 3) * s(i, j - 1, k)) / dz) + jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * (0.25 * (met(i + 1, j, k, 3, 1) &
                  &* s(i + 1, j, k) + met(i + 1, j, k + 1, 3, 1) * s(i + 1, j, &
                  &k + 1) - met(i - 1, j, k, 3, 1) * s(i - 1, j, k) - met(i &
                  &- 1, j, k + 1, 3, 1) * s(i - 1, j, k + 1)) / dx + 0.25 &
                  &* (met(i, j + 1, k, 3, 2) * s(i, j + 1, k) + met(i, j + 1, &
                  &k + 1, 3, 2) * s(i, j + 1, k + 1) - met(i, j - 1, k, 3, 2) &
                  &* s(i, j - 1, k) - met(i, j - 1, k + 1, 3, 2) * s(i, j - 1, &
                  &k + 1)) / dy + (met(i, j, k + 1, 3, 3) * s(i, j, k + 1) &
                  &- met(i, j, k, 3, 3) * s(i, j, k)) / dz) + jacInv / dz &
                  &* pEdgeUDiv / rhoEdgeU * 0.5 * (pStratTFC(i, j, k) / jac(i, &
                  &j, k) * (chris(i, j, k, 1, 1) + chris(i, j, k, 2, 2) + 2.0 &
                  &* chris(i, j, k, 1, 3) * met(i, j, k, 1, 3) + 2.0 &
                  &* chris(i, j, k, 2, 3) * met(i, j, k, 2, 3)) * s(i, j, k) &
                  &+ pStratTFC(i, j, k + 1) / jac(i, j, k + 1) * (chris(i, j, &
                  &k + 1, 1, 1) + chris(i, j, k + 1, 2, 2) + 2.0 * chris(i, j, &
                  &k + 1, 1, 3) * met(i, j, k + 1, 1, 3) + 2.0 * chris(i, j, k &
                  &+ 1, 2, 3) * met(i, j, k + 1, 2, 3)) * s(i, j, k + 1))

            else if(k == nz) then

              Ls(i, j, k) = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra &
                  &* ((s(i + 1, j, k) - s(i, j, k)) / dx + 0.25 * (met(i, j, k &
                  &- 2, 1, 3) * s(i, j, k - 2) + met(i + 1, j, k - 2, 1, 3) &
                  &* s(i + 1, j, k - 2) - 4.0 * met(i, j, k - 1, 1, 3) * s(i, &
                  &j, k - 1) - 4.0 * met(i + 1, j, k - 1, 1, 3) * s(i + 1, j, &
                  &k - 1) + 3.0 * met(i, j, k, 1, 3) * s(i, j, k) + 3.0 &
                  &* met(i + 1, j, k, 1, 3) * s(i + 1, j, k)) / dz) - jacInv &
                  &/ dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * ((s(i, j, k) &
                  &- s(i - 1, j, k)) / dx + 0.25 * (met(i, j, k - 2, 1, 3) &
                  &* s(i, j, k - 2) + met(i - 1, j, k - 2, 1, 3) * s(i - 1, j, &
                  &k - 2) - 4.0 * met(i, j, k - 1, 1, 3) * s(i, j, k - 1) &
                  &- 4.0 * met(i - 1, j, k - 1, 1, 3) * s(i - 1, j, k - 1) &
                  &+ 3.0 * met(i, j, k, 1, 3) * s(i, j, k) + 3.0 * met(i - 1, &
                  &j, k, 1, 3) * s(i - 1, j, k)) / dz) + jacInv / dy &
                  &* pEdgeFDiv / rhoEdgeF * pEdgeFGra * ((s(i, j + 1, k) &
                  &- s(i, j, k)) / dy + 0.25 * (met(i, j, k - 2, 2, 3) * s(i, &
                  &j, k - 2) + met(i, j + 1, k - 2, 2, 3) * s(i, j + 1, k - 2) &
                  &- 4.0 * met(i, j, k - 1, 2, 3) * s(i, j, k - 1) - 4.0 &
                  &* met(i, j + 1, k - 1, 2, 3) * s(i, j + 1, k - 1) + 3.0 &
                  &* met(i, j, k, 2, 3) * s(i, j, k) + 3.0 * met(i, j + 1, k, &
                  &2, 3) * s(i, j + 1, k)) / dz) - jacInv / dy * pEdgeBDiv &
                  &/ rhoEdgeB * pEdgeBGra * ((s(i, j, k) - s(i, j - 1, k)) &
                  &/ dy + 0.25 * (met(i, j, k - 2, 2, 3) * s(i, j, k - 2) &
                  &+ met(i, j - 1, k - 2, 2, 3) * s(i, j - 1, k - 2) - 4.0 &
                  &* met(i, j, k - 1, 2, 3) * s(i, j, k - 1) - 4.0 * met(i, j &
                  &- 1, k - 1, 2, 3) * s(i, j - 1, k - 1) + 3.0 * met(i, j, k, &
                  &2, 3) * s(i, j, k) + 3.0 * met(i, j - 1, k, 2, 3) * s(i, j &
                  &- 1, k)) / dz) - jacInv / dz * pEdgeDDiv / rhoEdgeD &
                  &* pEdgeDGra * (0.25 * (met(i + 1, j, k, 3, 1) * s(i + 1, j, &
                  &k) + met(i + 1, j, k - 1, 3, 1) * s(i + 1, j, k - 1) &
                  &- met(i - 1, j, k, 3, 1) * s(i - 1, j, k) - met(i - 1, j, k &
                  &- 1, 3, 1) * s(i - 1, j, k - 1)) / dx + 0.25 * (met(i, j &
                  &+ 1, k, 3, 2) * s(i, j + 1, k) + met(i, j + 1, k - 1, 3, 2) &
                  &* s(i, j + 1, k - 1) - met(i, j - 1, k, 3, 2) * s(i, j - 1, &
                  &k) - met(i, j - 1, k - 1, 3, 2) * s(i, j - 1, k - 1)) / dy &
                  &+ (met(i, j, k, 3, 3) * s(i, j, k) - met(i, j, k - 1, 3, 3) &
                  &* s(i, j, k - 1)) / dz) - jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                  &* (chris(i, j, k, 1, 1) + chris(i, j, k, 2, 2) + 2.0 &
                  &* chris(i, j, k, 1, 3) * met(i, j, k, 1, 3) + 2.0 &
                  &* chris(i, j, k, 2, 3) * met(i, j, k, 2, 3)) * s(i, j, k) &
                  &+ pStratTFC(i, j, k - 1) / jac(i, j, k - 1) * (chris(i, j, &
                  &k - 1, 1, 1) + chris(i, j, k - 1, 2, 2) + 2.0 * chris(i, j, &
                  &k - 1, 1, 3) * met(i, j, k - 1, 1, 3) + 2.0 * chris(i, j, k &
                  &- 1, 2, 3) * met(i, j, k - 1, 2, 3)) * s(i, j, k - 1))

            else

              Ls(i, j, k) = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra &
                  &* ((s(i + 1, j, k) - s(i, j, k)) / dx + 0.25 * (met(i, j, k &
                  &+ 1, 1, 3) * s(i, j, k + 1) + met(i + 1, j, k + 1, 1, 3) &
                  &* s(i + 1, j, k + 1) - met(i, j, k - 1, 1, 3) * s(i, j, k &
                  &- 1) - met(i + 1, j, k - 1, 1, 3) * s(i + 1, j, k - 1)) &
                  &/ dz) - jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra &
                  &* ((s(i, j, k) - s(i - 1, j, k)) / dx + 0.25 * (met(i, j, k &
                  &+ 1, 1, 3) * s(i, j, k + 1) + met(i - 1, j, k + 1, 1, 3) &
                  &* s(i - 1, j, k + 1) - met(i, j, k - 1, 1, 3) * s(i, j, k &
                  &- 1) - met(i - 1, j, k - 1, 1, 3) * s(i - 1, j, k - 1)) &
                  &/ dz) + jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra &
                  &* ((s(i, j + 1, k) - s(i, j, k)) / dy + 0.25 * (met(i, j, k &
                  &+ 1, 2, 3) * s(i, j, k + 1) + met(i, j + 1, k + 1, 2, 3) &
                  &* s(i, j + 1, k + 1) - met(i, j, k - 1, 2, 3) * s(i, j, k &
                  &- 1) - met(i, j + 1, k - 1, 2, 3) * s(i, j + 1, k - 1)) &
                  &/ dz) - jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra &
                  &* ((s(i, j, k) - s(i, j - 1, k)) / dy + 0.25 * (met(i, j, k &
                  &+ 1, 2, 3) * s(i, j, k + 1) + met(i, j - 1, k + 1, 2, 3) &
                  &* s(i, j - 1, k + 1) - met(i, j, k - 1, 2, 3) * s(i, j, k &
                  &- 1) - met(i, j - 1, k - 1, 2, 3) * s(i, j - 1, k - 1)) &
                  &/ dz) + jacInv / dz * pEdgeUDiv / rhoEdgeU * pEdgeUGra &
                  &* (0.25 * (met(i + 1, j, k, 3, 1) * s(i + 1, j, k) + met(i &
                  &+ 1, j, k + 1, 3, 1) * s(i + 1, j, k + 1) - met(i - 1, j, &
                  &k, 3, 1) * s(i - 1, j, k) - met(i - 1, j, k + 1, 3, 1) &
                  &* s(i - 1, j, k + 1)) / dx + 0.25 * (met(i, j + 1, k, 3, 2) &
                  &* s(i, j + 1, k) + met(i, j + 1, k + 1, 3, 2) * s(i, j + 1, &
                  &k + 1) - met(i, j - 1, k, 3, 2) * s(i, j - 1, k) - met(i, j &
                  &- 1, k + 1, 3, 2) * s(i, j - 1, k + 1)) / dy + (met(i, j, k &
                  &+ 1, 3, 3) * s(i, j, k + 1) - met(i, j, k, 3, 3) * s(i, j, &
                  &k)) / dz) + jacInv / dz * pEdgeUDiv / rhoEdgeU * 0.5 &
                  &* (pStratTFC(i, j, k) / jac(i, j, k) * (chris(i, j, k, 1, &
                  &1) + chris(i, j, k, 2, 2) + 2.0 * chris(i, j, k, 1, 3) &
                  &* met(i, j, k, 1, 3) + 2.0 * chris(i, j, k, 2, 3) * met(i, &
                  &j, k, 2, 3)) * s(i, j, k) + pStratTFC(i, j, k + 1) / jac(i, &
                  &j, k + 1) * (chris(i, j, k + 1, 1, 1) + chris(i, j, k + 1, &
                  &2, 2) + 2.0 * chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, &
                  &1, 3) + 2.0 * chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, &
                  &2, 3)) * s(i, j, k + 1)) - jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * (0.25 * (met(i + 1, j, k, 3, 1) &
                  &* s(i + 1, j, k) + met(i + 1, j, k - 1, 3, 1) * s(i + 1, j, &
                  &k - 1) - met(i - 1, j, k, 3, 1) * s(i - 1, j, k) - met(i &
                  &- 1, j, k - 1, 3, 1) * s(i - 1, j, k - 1)) / dx + 0.25 &
                  &* (met(i, j + 1, k, 3, 2) * s(i, j + 1, k) + met(i, j + 1, &
                  &k - 1, 3, 2) * s(i, j + 1, k - 1) - met(i, j - 1, k, 3, 2) &
                  &* s(i, j - 1, k) - met(i, j - 1, k - 1, 3, 2) * s(i, j - 1, &
                  &k - 1)) / dy + (met(i, j, k, 3, 3) * s(i, j, k) - met(i, j, &
                  &k - 1, 3, 3) * s(i, j, k - 1)) / dz) - jacInv / dz &
                  &* pEdgeDDiv / rhoEdgeD * 0.5 * (pStratTFC(i, j, k) / jac(i, &
                  &j, k) * (chris(i, j, k, 1, 1) + chris(i, j, k, 2, 2) + 2.0 &
                  &* chris(i, j, k, 1, 3) * met(i, j, k, 1, 3) + 2.0 &
                  &* chris(i, j, k, 2, 3) * met(i, j, k, 2, 3)) * s(i, j, k) &
                  &+ pStratTFC(i, j, k - 1) / jac(i, j, k - 1) * (chris(i, j, &
                  &k - 1, 1, 1) + chris(i, j, k - 1, 2, 2) + 2.0 * chris(i, j, &
                  &k - 1, 1, 3) * met(i, j, k - 1, 1, 3) + 2.0 * chris(i, j, k &
                  &- 1, 2, 3) * met(i, j, k - 1, 2, 3)) * s(i, j, k - 1))

            end if

          end do i_loop
        end do j_loop
      end do k_loop

    case("noTopography")

      ac_tfc = ac_b
      ar_tfc = ar_b
      al_tfc = al_b
      af_tfc = af_b
      ab_tfc = ab_b
      au_tfc = au_b
      ad_tfc = ad_b
      aru_tfc = aru_b
      ard_tfc = ard_b
      alu_tfc = alu_b
      ald_tfc = ald_b
      afu_tfc = afu_b
      afd_tfc = afd_b
      abu_tfc = abu_b
      abd_tfc = abd_b
      auu_tfc = auu_b
      add_tfc = add_b
      aruu_tfc = aruu_b
      ardd_tfc = ardd_b
      aluu_tfc = aluu_b
      aldd_tfc = aldd_b
      afuu_tfc = afuu_b
      afdd_tfc = afdd_b
      abuu_tfc = abuu_b
      abdd_tfc = abdd_b

      ach_tfc = ach_b

      topography = .false.
      call val_PsIn(var, dt, opt, facray)
      call linOpr(sIn, Ls, opt, "tot")
      topography = .true.

      if(master) then
        print *, "AC difference: ", maxval(abs(ac_tfc - ac_b))
        print *, "AR difference: ", maxval(abs(ar_tfc - ar_b))
        print *, "AL difference: ", maxval(abs(al_tfc - al_b))
        print *, "AF difference: ", maxval(abs(af_tfc - af_b))
        print *, "AB difference: ", maxval(abs(ab_tfc - ab_b))
        print *, "AU difference: ", maxval(abs(au_tfc - au_b))
        print *, "AD difference: ", maxval(abs(ad_tfc - ad_b))
        print *, "AC: ", maxval(abs(ac_tfc))
        print *, "AR: ", maxval(abs(ar_tfc))
        print *, "AL: ", maxval(abs(al_tfc))
        print *, "AF: ", maxval(abs(af_tfc))
        print *, "AB: ", maxval(abs(ab_tfc))
        print *, "AU: ", maxval(abs(au_tfc))
        print *, "AD: ", maxval(abs(ad_tfc))
        print *, "ARU: ", maxval(abs(aru_tfc))
        print *, "ARD: ", maxval(abs(ard_tfc))
        print *, "ALU: ", maxval(abs(alu_tfc))
        print *, "ALD: ", maxval(abs(ald_tfc))
        print *, "AFU: ", maxval(abs(afu_tfc))
        print *, "AFD: ", maxval(abs(afd_tfc))
        print *, "ABU: ", maxval(abs(abu_tfc))
        print *, "ABD: ", maxval(abs(abd_tfc))
        print *, "AUU: ", maxval(abs(auu_tfc))
        print *, "ADD: ", maxval(abs(add_tfc))
        print *, "ARUU: ", maxval(abs(aruu_tfc))
        print *, "ARDD: ", maxval(abs(ardd_tfc))
        print *, "ALUU: ", maxval(abs(aluu_tfc))
        print *, "ALDD: ", maxval(abs(aldd_tfc))
        print *, "AFUU: ", maxval(abs(afuu_tfc))
        print *, "AFDD: ", maxval(abs(afdd_tfc))
        print *, "ABUU: ", maxval(abs(abuu_tfc))
        print *, "ABDD: ", maxval(abs(abdd_tfc))

        print *, "ACH difference: ", maxval(abs(ach_tfc - ach_b))
        print *, "ACH: ", maxval(abs(ach_tfc))
      end if

      if(master) then
        print *, "Left-hand side difference", maxval(abs(Ls - LsTest))
      end if

    end select

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          diffOut(i, j, k) = abs(Ls(i, j, k) - LsTest(i, j, k)) / &
              &maxval((/abs(Ls(i, j, k)), abs(LsTest(i, j, k))/))
        end do
      end do
    end do

    if(master) then
      open(42, file = "linear_operator_test.dat", form = "unformatted", access &
          &= "direct", recl = sizeX * sizeY * sizeZ * sizeofreal4)
      write(42, rec = 1) diffOut
      close(42)
    end if

  end subroutine linearOperatorTestTFC

  !----------------------------------------------------------------------------

  subroutine calc_RHS(b, var, flux, dt, onlyinfo)

    !----------------------------------------
    !   calculates the RHS of the
    !   Poisson problem
    !----------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var
    type(flux_type), intent(in) :: flux
    real, intent(in) :: dt
    real, dimension(1:nx, 1:ny, 1:nz), intent(out) :: b ! RHS
    logical, intent(in) :: onlyinfo ! give info in div

    ! local vars
    real :: uR, uL, vF, vB, wU, wD
    real :: pStratU, pStratD
    real, dimension(1:nz) :: sum_local, sum_global

    ! TFC
    real :: pEdgeR, pEdgeL, pEdgeF, pEdgeB, pEdgeU, pEdgeD

    integer :: i, j, k
    real :: div, divSum, divSumScaled

    real :: bu, bv, bw, bl2loc, divL2_norm, divL2_norm_local

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

    integer :: i0, j0

    real :: rho, the
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz) :: heat
    real, dimension(- nbz:nz + nbz) :: w_0
    real, dimension(- nbz:nz + nbz) :: S_bar

    real :: fcscal

    if(giveInfo .and. master) then
      print *, ""
      print "(a)", repeat("-", 80)
      print *, "   calc_RHS: computing RHS... "
      print "(a)", repeat("-", 80)
      print *, ""
    end if

    ! ----------------------------------
    !          Poisson Problem
    ! ----------------------------------

    if(master .and. verbose) print *, "update.f90/poissonSolver:  Setting up &
        &Poisson problem."

    i0 = is + nbx - 1
    j0 = js + nby - 1

    ! heating due to thermal relaxation, molecular or turbulent diffusion,
    ! or gravity waves, and the horizontal-mean vertical wind due to it

    if(model == "compressible") then
      heat = 0.
      S_bar = 0.
      w_0 = 0.
    elseif(heatingONK14 .or. TurbScheme .or. rayTracer) then
      call heat_w0(var, flux, dt, heat, S_bar, w_0)
    else
      heat = 0.
      S_bar = 0.
      w_0 = 0.
    end if

    ! subtreact horizontal mean of heat(:,:,:)

    if(model /= "compressible") then
      do i = 1, nx
        do j = 1, ny
          heat(i, j, 1:nz) = heat(i, j, 1:nz) - S_bar(1:nz)
        end do
      end do
    end if

    ! scale RHS with Ma^2 * kappa, hence ...
    heat(:, :, :) = heat(:, :, :) * Ma ** 2 * kappa

    !--------------------------------------------------
    !    setup b = Ma^2 * P * u^*  (right hand side)
    !--------------------------------------------------

    divSum = 0.0
    divSumScaled = 0.0
    divL2 = 0.0
    divMax = 0.0

    divSum_local = 0.0
    divL2_local = 0.0

    divL2_norm = 0.0
    divL2_norm_local = 0.0

    if(RHS_diagnostics) then
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

    select case(model)

    case("pseudo_incompressible", "compressible")

      if(topography) then
        ! TFC
        ! Calculate RHS for TFC.
        do k = 1, nz
          ! Calculate scaling factor.
          fcscal = sqrt(pStrat(k) ** 2.0 / rhoStrat(k))
          do j = 1, ny
            do i = 1, nx
              ! Store velocities at cell edges.
              uR = var%u(i, j, k)
              uL = var%u(i - 1, j, k)
              vF = var%v(i, j, k)
              vB = var%v(i, j - 1, k)
              wU = var%w(i, j, k) - w_0(k)
              wD = var%w(i, j, k - 1) - w_0(k - 1)

              if(model == "compressible") then ! JP already in var
                pEdgeR = 1.0
                pEdgeL = 1.0
                pEdgeF = 1.0
                pEdgeB = 1.0
                pEdgeU = 1.0
                pEdgeD = 1.0
              else ! Compute pStrat at cell edges.
                pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, &
                    &j, k) * pStratTFC(i + 1, j, k))
                pEdgeL = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i - 1, &
                    &j, k) * pStratTFC(i - 1, j, k))
                pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                    &+ 1, k) * pStratTFC(i, j + 1, k))
                pEdgeB = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                    &- 1, k) * pStratTFC(i, j - 1, k))
                pEdgeU = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, &
                    &k + 1) * pStratTFC(i, j, k + 1))
                pEdgeD = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, &
                    &k - 1) * pStratTFC(i, j, k - 1))
              end if
              ! Compute RHS.
              bu = (pEdgeR * uR - pEdgeL * uL) / dx / jac(i, j, k) * Ma ** 2.0 &
                  &* kappa
              bv = (pEdgeF * vF - pEdgeB * vB) / dy / jac(i, j, k) * Ma ** 2.0 &
                  &* kappa
              bw = (pEdgeU * wU - pEdgeD * wD) / dz / jac(i, j, k) * Ma ** 2.0 &
                  &* kappa
              divSum_local = divSum_local + bu + bv + bw + heat(i, j, k)
              bu = bu / fcscal
              bv = bv / fcscal
              bw = bw / fcscal
              heat(i, j, k) = heat(i, j, k) / fcscal
              b(i, j, k) = bu + bv + bw + heat(i, j, k)
              ! Compute check sum for solvability criterion.
              divL2_local = divL2_local + b(i, j, k) ** 2.0
              bl2loc = bu ** 2.0 + bv ** 2.0 + bw ** 2.0 + (heat(i, j, k)) ** 2
              divL2_norm_local = divL2_norm_local + bl2loc
              if(RHS_diagnostics) then
                dPudx_norm_local = dPudx_norm_local + bu ** 2.0
                dPvdy_norm_local = dPvdy_norm_local + bv ** 2.0
                dPwdz_norm_local = dPwdz_norm_local + bw ** 2.0
                Q_norm_local = Q_norm_local + (heat(i, j, k)) ** 2
              end if
              if(abs(b(i, j, k)) > divMax) then
                divMax = abs(b(i, j, k))
              end if
            end do
          end do
        end do
      else
        do k = 1, nz
          fcscal = sqrt(Pstrat(k) ** 2 / rhoStrat(k))

          do j = 1, ny
            do i = 1, nx

              uR = var%u(i, j, k); uL = var%u(i - 1, j, k)
              vF = var%v(i, j, k); vB = var%v(i, j - 1, k)
              wU = var%w(i, j, k) - w_0(k); wD = var%w(i, j, k - 1) - w_0(k - 1)

              PstratU = PstratTilde(k)
              PstratD = PstratTilde(k - 1)

              ! scale RHS with Ma^2 * kappa, hence ...

              bu = Pstrat(k) * (uR - uL) / dx * Ma ** 2 * kappa
              bv = Pstrat(k) * (vF - vB) / dy * Ma ** 2 * kappa
              bw = (PstratU * wU - PstratD * wD) / dz * Ma ** 2 * kappa

              ! check sum for solvability criterion
              divSum_local = divSum_local + bu + bv + bw + heat(i, j, k)

              bu = bu / fcscal
              bv = bv / fcscal
              bw = bw / fcscal
              heat(i, j, k) = heat(i, j, k) / fcscal

              b(i, j, k) = bu + bv + bw + heat(i, j, k)

              ! L2-norm of the divergence div(Pu)
              divL2_local = divL2_local + b(i, j, k) ** 2

              ! introduce a reference norm for the RHS

              bl2loc = bu ** 2 + bv ** 2 + bw ** 2 + (heat(i, j, k)) ** 2
              divL2_norm_local = divL2_norm_local + bl2loc

              if(RHS_diagnostics) then
                dPudx_norm_local = dPudx_norm_local + bu ** 2
                dPvdy_norm_local = dPvdy_norm_local + bv ** 2
                dPwdz_norm_local = dPwdz_norm_local + bw ** 2
                Q_norm_local = Q_norm_local + (heat(i, j, k)) ** 2
              end if

              ! max norm of div(Pu)
              if(abs(b(i, j, k)) > divMax) then
                divMax = abs(b(i, j, k))
              end if

              ! Skalierung mit thetaStrat
              if(pressureScaling) then
                b(i, j, k) = b(i, j, k) / PStrat(k)
              end if

            end do
          end do
        end do
      end if

      !MPI: sum divSum_local over all procs
      root = 0
      call mpi_reduce(divSum_local, divSum, 1, mpi_double_precision, mpi_sum, &
          &root, comm, ierror)

      call mpi_bcast(divSum, 1, mpi_double_precision, root, comm, ierror)

      !MPI: sum divL2_local over all procs
      root = 0
      call mpi_reduce(divL2_local, divL2, 1, mpi_double_precision, mpi_sum, &
          &root, comm, ierror)

      call mpi_bcast(divL2, 1, mpi_double_precision, root, comm, ierror)

      !MPI: sum divL2_norm_local over all procs
      root = 0
      call mpi_reduce(divL2_norm_local, divL2_norm, 1, mpi_double_precision, &
          &mpi_sum, root, comm, ierror)

      call mpi_bcast(divL2_norm, 1, mpi_double_precision, root, comm, ierror)

      if(RHS_diagnostics) then
        call mpi_allreduce(dPudx_norm_local, dPudx_norm, 1, &
            &mpi_double_precision, mpi_sum, comm, ierror)
        call mpi_allreduce(dPvdy_norm_local, dPvdy_norm, 1, &
            &mpi_double_precision, mpi_sum, comm, ierror)
        call mpi_allreduce(dPwdz_norm_local, dPwdz_norm, 1, &
            &mpi_double_precision, mpi_sum, comm, ierror)
        call mpi_allreduce(Q_norm_local, Q_norm, 1, mpi_double_precision, &
            &mpi_sum, comm, ierror)

        dPudx_norm = sqrt(dPudx_norm / sizeX / sizeY / sizeZ)
        dPvdy_norm = sqrt(dPvdy_norm / sizeX / sizeY / sizeZ)
        dPwdz_norm = sqrt(dPwdz_norm / sizeX / sizeY / sizeZ)
        Q_norm = sqrt(Q_norm / sizeX / sizeY / sizeZ)
      end if

      ! scale div
      divL2_local = sqrt(divL2_local / nx / ny / nz)
      divL2 = sqrt(divL2 / sizeX / sizeY / sizeZ)

      divL2_norm_local = sqrt(divL2_norm_local / nx / ny / nz)
      divL2_norm = sqrt(divL2_norm / sizeX / sizeY / sizeZ)

      alpha_tol = divL2_norm
      b_norm = divL2

      if(divL2_norm /= 0.0) then
        tolref = divL2 / divL2_norm
      else
        if(divL2 == 0.0) then
          tolref = 1.0
        else
          stop 'ERROR: divL2_norm = 0 while divL2 /= 0'
        end if
      end if

      if(master) then
        print *, "tolref = ", tolref

        if(RHS_diagnostics) then
          print *, "RHS diagnostics:"
          print *, "dPudx_norm = ", dPudx_norm
          print *, "dPvdy_norm = ", dPvdy_norm
          print *, "dPwdz_norm = ", dPwdz_norm
          print *, "Q_norm = ", Q_norm
          print *, "alpha = ", divL2_norm
          print *, "unscaled |b|=", divL2
        end if
      end if

    case("Boussinesq")

      if(topography) then
        ! TFC
        ! Calculate RHS for TFC.
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              ! Compute values at cell edges.
              uR = 0.5 * (jac(i, j, k) + jac(i + 1, j, k)) * var%u(i, j, k)
              uL = 0.5 * (jac(i, j, k) + jac(i - 1, j, k)) * var%u(i - 1, j, k)
              vF = 0.5 * (jac(i, j, k) + jac(i, j + 1, k)) * var%v(i, j, k)
              vB = 0.5 * (jac(i, j, k) + jac(i, j - 1, k)) * var%v(i, j - 1, k)
              wU = 0.5 * (jac(i, j, k) + jac(i, j, k + 1)) * var%w(i, j, k)
              wD = 0.5 * (jac(i, j, k) + jac(i, j, k - 1)) * var%w(i, j, k - 1)
              ! Compute RHS.
              bu = Ma ** 2.0 * kappa / theta00 * (uR - uL) / dx / jac(i, j, k)
              bv = Ma ** 2.0 * kappa / theta00 * (vF - vB) / dy / jac(i, j, k)
              bw = Ma ** 2.0 * kappa / theta00 * (wU - wD) / dz / jac(i, j, k)
              b(i, j, k) = bu + bv + bw
              bl2loc = bu ** 2.0 + bv ** 2.0 + bw ** 2.0
              ! Compute check sum for solvability criterion.
              divSum_local = divSum_local + b(i, j, k)
              divL2_local = divL2_local + b(i, j, k) ** 2.0
              divL2_norm_local = divL2_norm_local + bl2loc
              if(abs(b(i, j, k)) > divMax) then
                divMax = abs(b(i, j, k))
              end if
            end do
          end do
        end do
      else
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx

              uR = var%u(i, j, k); uL = var%u(i - 1, j, k)
              vF = var%v(i, j, k); vB = var%v(i, j - 1, k)
              wU = var%w(i, j, k); wD = var%w(i, j, k - 1)

              bu = Ma ** 2 * kappa / theta00 * (uR - uL) / dx
              bv = Ma ** 2 * kappa / theta00 * (vF - vB) / dy
              bw = Ma ** 2 * kappa / theta00 * (wU - wD) / dz

              ! TFC
              b(i, j, k) = bu + bv + bw

              !               introduce a reference norm for the RHS
              !               div = (uR-uL)/dx + (vF-vB)/dy + (wU-wD)/dz
              bl2loc = bu ** 2 + bv ** 2 + bw ** 2

              ! check sum for solvability criterion (shoud be zero)
              ! divSum = divSum + b(i,j,k)
              divSum_local = divSum_local + b(i, j, k)

              ! TFC
              divL2_local = divL2_local + b(i, j, k) ** 2

              divL2_norm_local = divL2_norm_local + bl2loc

              ! TFC
              if(abs(b(i, j, k)) > divMax) then
                divMax = abs(b(i, j, k))
              end if

            end do
          end do
        end do
      end if

      !MPI: sum divSum_local over all procs
      root = 0
      call mpi_reduce(divSum_local, divSum, 1, mpi_double_precision, mpi_sum, &
          &root, comm, ierror)

      call mpi_bcast(divSum, 1, mpi_double_precision, root, comm, ierror)

      !MPI: sum divL2_local over all procs
      root = 0
      call mpi_reduce(divL2_local, divL2, 1, mpi_double_precision, mpi_sum, &
          &root, comm, ierror)

      call mpi_bcast(divL2, 1, mpi_double_precision, root, comm, ierror)

      !MPI: sum divL2_norm_local over all procs
      root = 0
      call mpi_reduce(divL2_norm_local, divL2_norm, 1, mpi_double_precision, &
          &mpi_sum, root, comm, ierror)

      call mpi_bcast(divL2_norm, 1, mpi_double_precision, root, comm, ierror)

      ! scale div
      divL2_local = sqrt(divL2_local / nx / ny / nz)
      divL2 = sqrt(divL2 / sizeX / sizeY / sizeZ)

      ! scale by number of cells
      divL2_local = sqrt(divL2_local / nx / ny / nz)
      divL2 = sqrt(divL2 / sizeX / sizeY / sizeZ)

      divL2_norm_local = sqrt(divL2_norm_local / nx / ny / nz)
      divL2_norm = sqrt(divL2_norm / sizeX / sizeY / sizeZ)

      alpha_tol = divL2_norm
      b_norm = divL2

      if(divL2_norm /= 0.0) then
        tolref = divL2 / divL2_norm
      else
        if(divL2 == 0.0) then
          tolref = 1.0
        else
          stop 'ERROR: divL2_norm = 0 while divL2 /= 0'
        end if
      end if

      if(master) print *, "tolref = ", tolref

    case default
      stop "poissonSolver: unknown case model."
    end select

    !-------------------------------
    !     Display info on screen
    !-------------------------------

    if(master .and. onlyinfo) then
      ! Information on divergence

      print *, ""
      print *, " Poisson Solver ", trim(poissonSolverType), ": Final state"
      print *, ""

      select case(model)
      case("Boussinesq")

        write(*, fmt = "(a25,es17.6)") "L2(div(u)) [1/s] = ", divL2 * uRef &
            &/ lRef

        write(*, fmt = "(a25,es17.6)") "max(div(u)) [1/s] = ", divMax * uRef &
            &/ lRef

        write(*, fmt = "(a25,es17.6)") "rms terms (div(u)) [1/s] = ", &
            &divL2_norm * uRef / lRef

        write(*, fmt = "(a25,es17.6)") "normalized L2(div(u)) = ", divL2 &
            &/ divL2_norm

      case("pseudo_incompressible", "compressible")

        write(*, fmt = "(a25,es17.6)") "L2(div(Pu)) [Pa/s] = ", divL2 * rhoRef &
            &* thetaRef * uRef / lRef

        write(*, fmt = "(a25,es17.6)") "max(div(Pu)) [Pa/s] = ", divMax &
            &* rhoRef * thetaRef * uRef / lRef

        write(*, fmt = "(a25,es17.6)") "rms terms (div(Pu)) [Pa/s] = ", &
            &divL2_norm * rhoRef * thetaRef * uRef / lRef

        write(*, fmt = "(a25,es17.6)") "normalized L2(div(Pu)) = ", tolref
      case default
        stop "calc_RHS: unkown case model"
      end select

      print *, ""
      print "(a)", repeat("-", 80)
      print *, ""
    else
      if(master .and. giveInfo) then
        print *, ""
        print *, " Poisson Solver ", trim(poissonSolverType), ": Initial state"
        print *, ""
        write(*, fmt = "(a25,es17.6)") "Sum over RHS = ", divSum
        write(*, fmt = "(a25,es17.6)") "Sum over scaled RHS  = ", divSumScaled

        ! Information on divergence
        select case(model)
        case("Boussinesq")

          write(*, fmt = "(a25,es17.6)") "L2(div(u)) [1/s] = ", divL2 * uRef &
              &/ lRef

          write(*, fmt = "(a25,es17.6)") "max(div(u)) [1/s] = ", divMax * uRef &
              &/ lRef

          write(*, fmt = "(a25,es17.6)") "rms terms (div(u)) [1/s] = ", &
              &divL2_norm * uRef / lRef

          write(*, fmt = "(a25,es17.6)") "normalized L2(div(u)) = ", divL2 &
              &/ divL2_norm

        case("pseudo_incompressible", "compressible")

          write(*, fmt = "(a25,es17.6)") "L2(div(Pu)) [Pa/s] = ", divL2 &
              &* rhoRef * thetaRef * uRef / lRef

          write(*, fmt = "(a25,es17.6)") "max(div(Pu)) [Pa/s] = ", divMax &
              &* rhoRef * thetaRef * uRef / lRef

          write(*, fmt = "(a25,es17.6)") "rms terms (div(Pu))  [Pa/s] = ", &
              &divL2_norm * rhoRef * thetaRef * uRef / lRef

          write(*, fmt = "(a25,es17.6)") "normalized L2(div(Pu)) = ", tolref

        case default
          stop "calc_RHS: unkown case model"
        end select
      end if
    end if

  end subroutine calc_RHS

  !----------------------------------------------------------------------

  subroutine poissonSolver(b, var, dt, errFlagBicg, nIter, m, opt, facray, &
      &facprs)
    ! -------------------------------------------------
    ! solves the Poisson problem with
    ! application of linear operator L
    ! -------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var
    real, intent(in) :: dt, facray, facprs

    logical, intent(out) :: errFlagBicg
    integer, intent(out) :: nIter

    integer, intent(in) :: m
    real, dimension(1:nx, 1:ny, 1:nz), intent(in) :: b ! RHS

    ! facray multiplies the Rayleigh-damping terms so that they are only
    ! handled in the implicit time stepping (sponge and immersed boundary)

    ! facprs multiplies the time step so that the routine can be used
    ! properly also in the implicit mode (where in sub-step 5 of the
    ! semi-implicit scheme the pressure correction is over a full
    ! time step, instead of half a time step)

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction
    ! of the winds and density fluctuations
    character(len = *), intent(in) :: opt

    ! local vars
    real, dimension(1:nx, 1:ny, 1:nz) :: sol ! solution of Poisson problem
    real :: res

    real :: dtInv

    ! verbose
    logical, parameter :: giveInfo = .true.

    integer :: k
    real :: fcscal

    ! TFC
    integer :: i, j

    ! Init
    if(dt == 0.0) stop "poissonSolver: dt = 0.0. Stopping."
    dtInv = 1.0 / dt

    !--------------------------------
    !     Linear equation solver
    !     solve for dt * dp ...
    !--------------------------------

    sol = 0.0

    select case(poissonSolverType)

      ! bicgstab solver
    case("bicgstab")

      select case(model)

      case("pseudo_incompressible", "compressible")

        call val_PsIn(var, dt, opt, facray)

      case("Boussinesq")

        ! Update tensor elements when topography is growing.
        if(topography .and. topographyTime > 0.0) then
          call val_PsIn(var, dt, opt, facray)
          if(opt == "expl") then
            expEle = .true.
            impEle = .false.
          else if(opt == "impl") then
            impEle = .true.
            expEle = .false.
          end if
        end if

        ! TFC
        ! Tensor elements are constant. Due to the scaling, the subroutine also
        ! works for the Boussinesq model.
        if(opt == "expl" .and. .not. expEle) then
          call val_PsIn(var, dt, opt, facray)
          expEle = .true.
          impEle = .false.
        else if(opt == "impl" .and. .not. impEle) then
          call val_PsIn(var, dt, opt, facray)
          impEle = .true.
          expEle = .false.
        end if

        ! stop'ERROR: BiCGStab still to be made ready for Boussinesq'

      case default
        stop "linOpr: unknown case model"
      end select

      call bicgstab(b, dt, sol, res, nIter, errFlagBicg, opt)

    case("gcr")

      stop 'ERROR: no gcr provided anymore'

    case("adi")

      stop 'ERROR: no adi provided anymore'

    case("hypre")

      stop "ERROR: no hypre provided anymore"

    case default
      stop "Unknown PoissonSolver. Stop"
    end select

    if(errFlagBicg) return

    select case(model)
    case("pseudo_incompressible", "compressible")
      do k = 1, nz
        fcscal = sqrt(Pstrat(k) ** 2 / rhoStrat(k))
        sol(:, :, k) = sol(:, :, k) / fcscal
      end do
    case default
    end select

    ! now get dp from dt * dp ...
    if(topography) then
      ! TFC
      ! Solution must be divided by Jacobian.
      do i = 1, nx
        do j = 1, ny
          do k = 1, nz
            dp(i, j, k) = dtInv / facprs / jac(i, j, k) * sol(i, j, k)
          end do
        end do
      end do
    else
      ! pass solution to pressure corrector
      dp(1:nx, 1:ny, 1:nz) = dtInv / facprs * sol
    end if

  end subroutine poissonSolver

  !----------------------------------------------------------------------

  subroutine thomas(l, c, r, sol, b)
    ! -------------------------------------
    !      Thomas algorithmus
    !      using diagonals l,c,r
    ! -------------------------------------

    ! In/Out variables
    real, dimension(:), intent(inout) :: l, c, r ! matrix diagonals:
    ! left, center, right
    real, dimension(:), intent(inout) :: b ! rhs
    real, dimension(:), intent(out) :: sol ! solution

    ! Local vars
    integer :: i, n

    n = size(b)

    ! Forward sweep
    do i = 2, n
      c(i) = c(i) - l(i) / c(i - 1) * r(i - 1)
      b(i) = b(i) - l(i) / c(i - 1) * b(i - 1)
    end do

    ! Backward sweep
    do i = n - 1, 1, - 1
      b(i) = b(i) - r(i) / c(i + 1) * b(i + 1)
    end do

    ! Solve
    do i = 1, n
      sol(i) = b(i) / c(i)
    end do

  end subroutine thomas

  !----------------------------------------------------------------------

  subroutine bicgstab(b_in, dt, sol, res, nIter, errFlag, opt)
    ! --------------------------------------
    !    BiCGStab using linear operator
    !    preconditioner applied via A M^-1 M x = b
    !---------------------------------------

    ! in/out variables
    real, dimension(1:nx, 1:ny, 1:nz), intent(in) :: b_in ! RHS
    real, intent(in) :: dt
    real, dimension(1:nx, 1:ny, 1:nz), intent(inout) :: sol
    real, intent(out) :: res ! residual
    integer, intent(out) :: nIter
    logical, intent(out) :: errFlag

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction
    ! of the winds and density fluctuations
    character(len = *), intent(in) :: opt

    ! Local parameters
    integer :: maxIt

    ! local variables
    integer :: i, j, k, allocstat
    integer :: j_b
    real, dimension(:, :, :), allocatable :: p, r0, rOld, r, s, t, v, matVec, &
        &v_pc
    real :: alpha, beta, omega

    ! verbose
    logical, parameter :: giveInfo = .true.

    ! MPI stuff
    integer :: root
    real :: res_local

    real, dimension(1:nx, 1:ny, 1:nz) :: b ! RHS
    real, dimension(1:nx, 1:ny) :: r_vm, b_vm
    real :: b_vm_norm, res_vm

    if(giveInfo .and. master) then
      print *, ""
      print "(a)", repeat("-", 80)
      print *, "BICGSTAB: solving linear system... "
      print "(a)", repeat("-", 80)
      print *, ""
    end if

    sol = 0.0

    ! Set parameters
    maxIt = maxIterPoisson

    ! modified convergence criterion so that iterations stop when either
    ! (a) tolcrit = abs  =>  |Ax - b| < eps b_*
    !     with b_* a suitable norm deciding whether b (= the divergence
    !     criterion for the winds from the predictor) is small or not
    ! (b) tolcrit = rel  =>  |Ax - b| < eps |b|
    ! here eps = tolPoisson is the user-set convergence criterion
    ! hypre has the criterion |Ax - b| < tol * |b|, hence, with
    ! tolref = divL2/divL2_norm = |b|/b_*

    if(tolcrit == "abs") then
      tol = tolPoisson / tolref
    else if(tolcrit == "rel") then
      tol = tolPoisson
    end if

    ! Allocate local fields
    allocate(p(1:nx, 1:ny, 1:nz), stat = allocstat); if(allocstat /= 0) stop &
        &"bicgstab:alloc failed"
    allocate(r0(1:nx, 1:ny, 1:nz), stat = allocstat); if(allocstat /= 0) stop &
        &"bicgstab:alloc failed"
    allocate(rOld(1:nx, 1:ny, 1:nz), stat = allocstat); if(allocstat /= 0) &
        &stop "bicgstab:alloc failed"
    allocate(r(1:nx, 1:ny, 1:nz), stat = allocstat); if(allocstat /= 0) stop &
        &"bicgstab:alloc failed"
    allocate(s(1:nx, 1:ny, 1:nz), stat = allocstat); if(allocstat /= 0) stop &
        &"bicgstab:alloc failed"
    allocate(t(1:nx, 1:ny, 1:nz), stat = allocstat); if(allocstat /= 0) stop &
        &"bicgstab:alloc failed"
    allocate(v(1:nx, 1:ny, 1:nz), stat = allocstat); if(allocstat /= 0) stop &
        &"bicgstab:alloc failed"
    allocate(matVec(1:nx, 1:ny, 1:nz), stat = allocstat); if(allocstat /= 0) &
        &stop "bicgstab:alloc failed"
    allocate(v_pc(1:nx, 1:ny, 1:nz), stat = allocstat); if(allocstat /= 0) &
        &stop "bicgstab:alloc failed"

    ! error flag
    errFlag = .false.

    b = b_in

    call linOpr(sol, matVec, opt, 'tot')
    r0 = b - matVec
    p = r0
    r = r0

    res_local = 0.0
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          res_local = res_local + r(i, j, k) ** 2
        end do
      end do
    end do

    !MPI find global residual
    root = 0
    call mpi_reduce(res_local, res, 1, mpi_double_precision, mpi_sum, root, &
        &comm, ierror)

    call mpi_bcast(res, 1, mpi_double_precision, root, comm, ierror)

    res = sqrt(res / sizeX / sizeY / sizeZ)

    b_norm = res

    r_vm = 0.0
    do k = 1, nz
      r_vm(:, :) = r_vm(:, :) + r(:, :, k)
    end do
    r_vm = r_vm / sizeZ

    res_local = 0.0
    do j = 1, ny
      do i = 1, nx
        res_local = res_local + r_vm(i, j) ** 2
      end do
    end do

    root = 0
    call mpi_reduce(res_local, res_vm, 1, mpi_double_precision, mpi_sum, root, &
        &comm, ierror)

    call mpi_bcast(res_vm, 1, mpi_double_precision, root, comm, ierror)

    res_vm = sqrt(res_vm / sizeX / sizeY)

    b_vm_norm = res_vm

    if(master) then
      print *, ""
      print *, " BiCGStab solver: "
      if(res == 0.0) then
        if(giveInfo) write(*, fmt = "(a25,es17.6)") " Initial residual: res &
            &= ", res
      else
        if(giveInfo) write(*, fmt = "(a25,es17.6)") " Initial residual: res &
            &= ", res / b_norm
      end if
      if(res_vm == 0.0) then
        if(giveInfo) write(*, fmt = "(a25,es17.6)") " Initial residual: res_vm &
            &= ", res_vm
      else
        if(giveInfo) write(*, fmt = "(a25,es17.6)") " Initial residual: res_vm &
            &= ", res_vm / b_vm_norm
      end if
      if(giveInfo) write(*, fmt = "(a25,es17.6)") " tol = ", tol
    end if

    if(res == 0.0 .or. res / b_norm <= tol) then
      if(master .and. giveInfo) print *, " ==> no iteration needed."
      nIter = 0
      return
    end if

    ! Loop

    iteration: do j_b = 1, maxIt

      ! v = A*p
      if(preconditioner == 'yes') then
        call preCond(p, v_pc, opt)
      else
        v_pc = p
      end if
      call linOpr(v_pc, matVec, opt, 'tot')
      v = matVec

      alpha = dot_product3D_glob(r, r0) / dot_product3D_glob(v, r0)
      s = r - alpha * v

      ! t = A*s
      if(preconditioner == 'yes') then
        call preCond(s, v_pc, opt)
      else
        v_pc = s
      end if
      call linOpr(v_pc, matVec, opt, 'tot')
      t = matVec

      omega = dot_product3D_glob(t, s) / dot_product3D_glob(t, t)
      sol = sol + alpha * p + omega * s

      rOld = r
      r = s - omega * t

      !-----------------------
      !   Abort criterion
      !-----------------------

      res_local = 0.0
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            res_local = res_local + r(i, j, k) ** 2
          end do
        end do
      end do

      !MPI find global residual
      root = 0
      call mpi_reduce(res_local, res, 1, mpi_double_precision, mpi_sum, root, &
          &comm, ierror)

      call mpi_bcast(res, 1, mpi_double_precision, root, comm, ierror)

      res = sqrt(res / sizeX / sizeY / sizeZ)

      r_vm = 0.0
      do k = 1, nz
        r_vm(:, :) = r_vm(:, :) + r(:, :, k)
      end do
      r_vm = r_vm / sizeZ

      res_local = 0.0
      do j = 1, ny
        do i = 1, nx
          res_local = res_local + r_vm(i, j) ** 2
        end do
      end do

      root = 0
      call mpi_reduce(res_local, res_vm, 1, mpi_double_precision, mpi_sum, &
          &root, comm, ierror)

      call mpi_bcast(res_vm, 1, mpi_double_precision, root, comm, ierror)

      res_vm = sqrt(res_vm / sizeX / sizeY)

      if(max(res / b_norm, res_vm / b_vm_norm) <= tol) then
        if(master .and. giveInfo) then
          write(*, fmt = "(a25,i25)") " Nb.of iterations: j = ", j_b
          write(*, fmt = "(a25,es17.6)") " Final residual: res = ", res / b_norm
          write(*, fmt = "(a25,es17.6)") " Final residual v.m. = ", res_vm &
              &/ b_vm_norm
          print *, ""
        end if

        nIter = j_b

        if(preconditioner == 'yes') then
          s = sol
          call preCond(s, sol, opt)
        end if

        return
      end if

      beta = alpha / omega * dot_product3D_glob(r, r0) &
          &/ dot_product3D_glob(rOld, r0)
      p = r + beta * (p - omega * v)

    end do iteration

    ! max iteration

    if(master) then ! modified by Junhong Wei (20161107)
      write(*, fmt = "(a25,i25)") " Bicgstab: max iterations!!!", maxIt
      write(*, fmt = "(a25,es17.6)") " Final BICGSTAB residual = ", res / b_norm
      print "(a)", repeat("-", 80)
      print *, ""

    end if ! modified by Junhong Wei (20161107)

    errFlag = .true.
    nIter = j_b

    ! deallocate local fields
    deallocate(p, stat = allocstat); if(allocstat /= 0) stop "bicgstab:dealloc &
        &p failed"
    deallocate(r0, stat = allocstat); if(allocstat /= 0) stop &
        &"bicgstab:dealloc r0 failed"
    deallocate(rOld, stat = allocstat); if(allocstat /= 0) stop &
        &"bicgstab:dealloc rOld failed"
    deallocate(r, stat = allocstat); if(allocstat /= 0) stop "bicgstab:dealloc &
        &r failed"
    deallocate(s, stat = allocstat); if(allocstat /= 0) stop "bicgstab:dealloc &
        &s failed"
    deallocate(t, stat = allocstat); if(allocstat /= 0) stop "bicgstab:dealloc &
        &t failed"
    deallocate(v, stat = allocstat); if(allocstat /= 0) stop "bicgstab:dealloc &
        &v failed"
    deallocate(v_pc, stat = allocstat); if(allocstat /= 0) stop &
        &"bicgstab:dealloc v_pcfailed"
    deallocate(matVec, stat = allocstat); if(allocstat /= 0) stop &
        &"bicgstab:dealloc matvec failed"

  end subroutine bicgstab

  !-----------------------------------------------------------------------

  subroutine pressureBoundaryCondition
    !--------------------------------------------------
    ! set pressure correction dp in ghost cells for BC
    !--------------------------------------------------

    ! auxiliary fields for "dp"
    real, dimension(0:ny + 1, 0:nz + 1) :: xSliceLeft_send, xSliceRight_send
    real, dimension(0:ny + 1, 0:nz + 1) :: xSliceLeft_recv, xSliceRight_recv

    real, dimension(0:nx + 1, 0:nz + 1) :: ySliceBack_send, ySliceForw_send
    real, dimension(0:nx + 1, 0:nz + 1) :: ySliceBack_recv, ySliceForw_recv

    ! MPI variables
    integer :: dest, source, tag
    integer :: sendcount, recvcount

    logical, parameter :: giveInfo = .true.

    ! Find neighbour procs
    if(idim > 1) call mpi_cart_shift(comm, 0, 1, left, right, ierror)
    if(jdim > 1) call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

    if(giveInfo .and. master) then
      print *, ""
      print "(a)", repeat("-", 80)
      print *, "pressureBoundaryCondition: setting dp in Halos..."
      print "(a)", repeat("-", 80)
      print *, ""
    end if

    !----------------------------
    !   set Halo cells: xSlice
    !----------------------------

    if(idim > 1) then
      ! slice size
      sendcount = (ny + 2) * (nz + 2)
      recvcount = sendcount

      ! read slice into contiguous array
      xSliceLeft_send(:, :) = dp(1, 0:ny + 1, :)
      xSliceRight_send(:, :) = dp(nx, 0:ny + 1, :)

      ! left -> right
      source = left
      dest = right
      tag = 100

      call mpi_sendrecv(xSliceRight_send(0, 0), sendcount, &
          &mpi_double_precision, dest, tag, xSliceLeft_recv(0, 0), recvcount, &
          &mpi_double_precision, source, mpi_any_tag, comm, sts_left, ierror)

      ! right -> left
      source = right
      dest = left
      tag = 100

      call mpi_sendrecv(xSliceLeft_send(0, 0), sendcount, &
          &mpi_double_precision, dest, tag, xSliceRight_recv(0, 0), recvcount, &
          &mpi_double_precision, source, mpi_any_tag, comm, sts_right, ierror)

      ! right halos
      dp(nx + 1, 0:ny + 1, :) = xSliceRight_recv(0:ny + 1, :)

      ! left halos
      dp(0, 0:ny + 1, :) = xSliceLeft_recv(0:ny + 1, :)

    else

      dp(0, :, :) = dp(nx, :, :)
      dp(nx + 1, :, :) = dp(1, :, :)

    end if

    if(verbose .and. master) print *, "horizontalHalos:  x-horizontal halos &
        &copied."

    !------------------------------
    !   set Halo cells: ySlice
    !------------------------------

    if(jdim > 1) then
      ! slice size
      sendcount = (nx + 2) * (nz + 2)

      recvcount = sendcount

      ! read slice into contiguous array
      ySliceBack_send(:, :) = dp(0:nx + 1, 1, :)

      ySliceForw_send(:, :) = dp(0:nx + 1, ny, :)

      ! back -> forw
      source = back
      dest = forw
      tag = 100

      call mpi_sendrecv(ySliceForw_send(0, 0), sendcount, &
          &mpi_double_precision, dest, tag, ySliceBack_recv(0, 0), recvcount, &
          &mpi_double_precision, source, mpi_any_tag, comm, sts_back, ierror)

      ! forw -> back
      source = forw
      dest = back
      tag = 100

      call mpi_sendrecv(ySliceBack_send(0, 0), sendcount, &
          &mpi_double_precision, dest, tag, ySliceForw_recv(0, 0), recvcount, &
          &mpi_double_precision, source, mpi_any_tag, comm, sts_right, ierror)

      ! forward halos
      dp(0:nx + 1, ny + 1, :) = ySliceForw_recv(0:nx + 1, :)

      ! backward halos
      dp(0:nx + 1, 0, :) = ySliceBack_recv(0:nx + 1, :)

    else

      dp(:, 0, :) = dp(:, ny, :)
      dp(:, ny + 1, :) = dp(:, 1, :)

    end if

    if(verbose .and. master) print *, "horizontalHalos:  x-horizontal halos &
        &copied."

    !----------------
    !   z-Boundary
    !----------------

    select case(zBoundary)

    case("periodic")

      dp(:, :, 0) = dp(:, :, nz)
      dp(:, :, nz + 1) = dp(:, :, 1)

    case("solid_wall")

      ! zero gradient
      dp(:, :, 0) = dp(:, :, 1)
      dp(:, :, nz + 1) = dp(:, :, nz)

    case default
      stop "pressureBoundaryCondition: unknown case zBoundary."
    end select

  end subroutine pressureBoundaryCondition

  !-----------------------------------------------------------------------

  subroutine correctorStep(var, dMom, dt, RKstage, opt, facray, facprs)

    !-------------------------------------------------
    !         correct pressure & velocity
    !-------------------------------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz, 3), &
        &intent(inout) :: dMom
    real, intent(in) :: dt, facray, facprs
    integer, intent(in) :: RKstage

    ! facray multiplies the Rayleigh-damping terms so that they are only
    ! handled in the implicit time stepping (sponge and immersed boundary)

    ! facprs multiplies the time step so that the routine can be used
    ! properly also in the implicit mode (where in sub-step 5 of the
    ! semi-implicit scheme the pressure correction is over a full
    ! time step, instead of half a time step)

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction
    ! of the winds and density fluctuations
    character(len = *), intent(in) :: opt

    ! local variables
    integer :: i, j, k
    real :: pEdge, pEdge_0, rhoEdge, rhou, rhov, rho, rho10, rho01
    real :: pGradX, pGradY, pGradZ
    real :: du, dv, dw, db
    real :: facu, facv, facw, facr
    real, dimension(0:ny + 1) :: f_cor_nd
    real :: bvsstw

    real :: rhov0m, rhov00, rhov1m, rhov10
    real :: rhoum0, rhou00, rhoum1, rhou01
    real :: rhow0, rhowm

    ! divergence check
    real :: maxDivPu, divPu
    real :: uR, uL, vF, vB, wU, wD
    real :: pStratU, pStratD

    ! TFC
    real :: rhoStratEdgeU
    real :: pEdgeR, pEdgeF, pEdgeU, pEdgeD
    real :: dpEdgeR, dpUEdgeR, dpUUEdgeR, dpDEdgeR, dpDDEdgeR, dpEdgeF, &
        &dpUEdgeF, dpUUEdgeF, dpDEdgeF, dpDDEdgeF, dpREdgeU, dpLEdgeU, &
        &dpFEdgeU, dpBEdgeU, dpREdgeD, dpLEdgeD, dpFEdgeD, dpBEdgeD
    real :: chris11EdgeU, chris22EdgeU, chris13EdgeU, chris23EdgeU, &
        &chris11EdgeD, chris22EdgeD, chris13EdgeD, chris23EdgeD
    real :: pGradZEdgeU, pGradZEdgeD
    real, dimension((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        &+ nbz)) :: corX, corY

    ! TFC
    integer :: k0, k1

    ! SK compressible
    real :: JPR, JPF, JPU

    ! verbose
    logical, parameter :: giveInfo = .true.

    integer :: i0, j0

    real :: ymin, ymax, yloc

    real :: f_cor_v

    i0 = is + nbx - 1
    j0 = js + nby - 1

    if(corset == 'periodic') then
      ymax = ly_dim(1) / lRef
      ymin = ly_dim(0) / lRef

      j0 = js + nby - 1

      do j = 0, ny + 1
        yloc = y(j + j0)

        f_cor_nd(j) = - 4. * pi / 8.64e4 * tRef * cos(2. * pi * (yloc - ymin) &
            &/ (ymax - ymin))
      end do
    else if(corset == 'constant') then
      f_cor_nd(0:ny + 1) = f_Coriolis_dim * tRef
    else
      stop 'ERROR: wrong corset'
    end if

    ! --------------------------------------
    !             calc p + dp
    ! --------------------------------------

    var%pi(0:nx + 1, 0:ny + 1, 0:nz + 1) = var%pi(0:nx + 1, 0:ny + 1, 0:nz &
        &+ 1) + dp(0:nx + 1, 0:ny + 1, 0:nz + 1)

    if(timeScheme == "semiimplicit") then
      if(opt == "impl") then
        kr_sp = kr_sp * facray
        alprlx = alprlx * facray
        if(topography) then
          kr_sp_tfc = kr_sp_tfc * facray
          kr_sp_w_tfc = kr_sp_w_tfc * facray
        end if
      end if
    end if

    ! --------------------------------------
    !           calc du and u + du
    ! --------------------------------------

    if(timeScheme == "semiimplicit") then
      if(opt == "impl") then
        do k = 1, nz
          do j = 1, ny
            do i = 0, nx
              facu = 1.0

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  facu = facu + dt * kv_hs(j, k)
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                if(topography) then
                  facu = facu + dt * 0.5 * (kr_sp_tfc(i, j, k) + kr_sp_tfc(i &
                      &+ 1, j, k))
                else
                  facu = facu + dt * kr_sp(j, k)
                end if
              end if

              facv = facu

              if(topography) then
                ! TFC
                ! Compute values at cell edges.
                rhou = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k) &
                    &+ rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, j, k))
                pEdgeR = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    &+ pStratTFC(i + 1, j, k) / jac(i + 1, j, k))
                ! Compute pressure difference gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  dpUUEdgeR = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 1, 3) &
                      &* dp(i, j, k + 2) + jac(i + 1, j, k + 2) * met(i + 1, &
                      &j, k + 2, 1, 3) * dp(i + 1, j, k + 2))
                  dpUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
                      &* dp(i, j, k + 1) + jac(i + 1, j, k + 1) * met(i + 1, &
                      &j, k + 1, 1, 3) * dp(i + 1, j, k + 1))
                  dpEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * dp(i, &
                      &j, k) + jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) &
                      &* dp(i + 1, j, k))
                  pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((jac(i + 1, j, &
                      &k) * dp(i + 1, j, k) - jac(i, j, k) * dp(i, j, k)) / dx &
                      &+ (- dpUUEdgeR + 4.0 * dpUEdgeR - 3.0 * dpEdgeR) * 0.5 &
                      &/ dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  dpDDEdgeR = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 1, 3) &
                      &* dp(i, j, k - 2) + jac(i + 1, j, k - 2) * met(i + 1, &
                      &j, k - 2, 1, 3) * dp(i + 1, j, k - 2))
                  dpDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, 3) &
                      &* dp(i, j, k - 1) + jac(i + 1, j, k - 1) * met(i + 1, &
                      &j, k - 1, 1, 3) * dp(i + 1, j, k - 1))
                  dpEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * dp(i, &
                      &j, k) + jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) &
                      &* dp(i + 1, j, k))
                  pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((jac(i + 1, j, &
                      &k) * dp(i + 1, j, k) - jac(i, j, k) * dp(i, j, k)) / dx &
                      &+ (dpDDEdgeR - 4.0 * dpDEdgeR + 3.0 * dpEdgeR) * 0.5 &
                      &/ dz)
                else
                  dpUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
                      &* dp(i, j, k + 1) + jac(i + 1, j, k + 1) * met(i + 1, &
                      &j, k + 1, 1, 3) * dp(i + 1, j, k + 1))
                  dpDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, 3) &
                      &* dp(i, j, k - 1) + jac(i + 1, j, k - 1) * met(i + 1, &
                      &j, k - 1, 1, 3) * dp(i + 1, j, k - 1))
                  pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((jac(i + 1, j, &
                      &k) * dp(i + 1, j, k) - jac(i, j, k) * dp(i, j, k)) / dx &
                      &+ (dpUEdgeR - dpDEdgeR) * 0.5 / dz)
                end if
                ! Compute velocity correction.
                corX(i, j, k) = facprs * dt / facu * pGradX
                if(model == "compressible") then ! multiply with JP
                  JPR = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i + 1, j, &
                      &k) * var%P(i + 1, j, k))

                  du = - JPR * corX(i, j, k)
                else
                  du = - corX(i, j, k)
                end if
              else
                rhov0m = 0.5 * (var%rho(i, j - 1, k) + var%rho(i, j, k)) &
                    &+ rhoStrat(k)

                rhov00 = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k)) &
                    &+ rhoStrat(k)

                rhov1m = 0.5 * (var%rho(i + 1, j - 1, k) + var%rho(i + 1, j, &
                    &k)) + rhoStrat(k)

                rhov10 = 0.5 * (var%rho(i + 1, j, k) + var%rho(i + 1, j + 1, &
                    &k)) + rhoStrat(k)

                rhou = 0.5 * (var%rho(i + 1, j, k) + var%rho(i, j, k)) &
                    &+ rhoStrat(k)

                pGradX = kappaInv * MaInv2 * pStrat(k) / rhou * (dp(i + 1, j, &
                    &k) - dp(i, j, k)) / dx

                pGradY = kappaInv * MaInv2 * 0.25 * (pStrat(k) / rhov0m &
                    &* (dp(i, j, k) - dp(i, j - 1, k)) / dy + pStrat(k) &
                    &/ rhov00 * (dp(i, j + 1, k) - dp(i, j, k)) / dy &
                    &+ pStrat(k) / rhov1m * (dp(i + 1, j, k) - dp(i + 1, j &
                    &- 1, k)) / dy + pStrat(k) / rhov10 * (dp(i + 1, j + 1, k) &
                    &- dp(i + 1, j, k)) / dy)

                du = - facprs * dt / (facu * facv + (f_cor_nd(j) * dt) ** 2) &
                    &* (facv * pGradX + f_cor_nd(j) * dt * pGradY)
              end if

              var%u(i, j, k) = var%u(i, j, k) + du
            end do
          end do
        end do
      else if(opt == "expl") then
        if(facprs /= 1.) stop 'ERROR: wrong facprs in explicit sub-step'
        do k = 1, nz
          do j = 1, ny
            do i = 0, nx
              if(topography) then
                ! TFC
                ! Compute values at cell edges.
                rhou = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k) &
                    &+ rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, j, k))
                pEdgeR = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    &+ pStratTFC(i + 1, j, k) / jac(i + 1, j, k))
                ! Compute pressure difference gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  dpUUEdgeR = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 1, 3) &
                      &* dp(i, j, k + 2) + jac(i + 1, j, k + 2) * met(i + 1, &
                      &j, k + 2, 1, 3) * dp(i + 1, j, k + 2))
                  dpUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
                      &* dp(i, j, k + 1) + jac(i + 1, j, k + 1) * met(i + 1, &
                      &j, k + 1, 1, 3) * dp(i + 1, j, k + 1))
                  dpEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * dp(i, &
                      &j, k) + jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) &
                      &* dp(i + 1, j, k))
                  pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((jac(i + 1, j, &
                      &k) * dp(i + 1, j, k) - jac(i, j, k) * dp(i, j, k)) / dx &
                      &+ (- dpUUEdgeR + 4.0 * dpUEdgeR - 3.0 * dpEdgeR) * 0.5 &
                      &/ dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  dpDDEdgeR = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 1, 3) &
                      &* dp(i, j, k - 2) + jac(i + 1, j, k - 2) * met(i + 1, &
                      &j, k - 2, 1, 3) * dp(i + 1, j, k - 2))
                  dpDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, 3) &
                      &* dp(i, j, k - 1) + jac(i + 1, j, k - 1) * met(i + 1, &
                      &j, k - 1, 1, 3) * dp(i + 1, j, k - 1))
                  dpEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * dp(i, &
                      &j, k) + jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) &
                      &* dp(i + 1, j, k))
                  pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((jac(i + 1, j, &
                      &k) * dp(i + 1, j, k) - jac(i, j, k) * dp(i, j, k)) / dx &
                      &+ (dpDDEdgeR - 4.0 * dpDEdgeR + 3.0 * dpEdgeR) * 0.5 &
                      &/ dz)
                else
                  dpUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
                      &* dp(i, j, k + 1) + jac(i + 1, j, k + 1) * met(i + 1, &
                      &j, k + 1, 1, 3) * dp(i + 1, j, k + 1))
                  dpDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, 3) &
                      &* dp(i, j, k - 1) + jac(i + 1, j, k - 1) * met(i + 1, &
                      &j, k - 1, 1, 3) * dp(i + 1, j, k - 1))
                  pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((jac(i + 1, j, &
                      &k) * dp(i + 1, j, k) - jac(i, j, k) * dp(i, j, k)) / dx &
                      &+ (dpUEdgeR - dpDEdgeR) * 0.5 / dz)
                end if
              else
                rhou = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k)) &
                    &+ rhoStrat(k)

                pGradX = kappaInv * MaInv2 * pStrat(k) / rhou * (dp(i + 1, j, &
                    &k) - dp(i, j, k)) / dx
              end if

              if(model == "compressible") then ! multiply with JP
                JPR = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i + 1, j, k) &
                    &* var%P(i + 1, j, k))

                du = - dt * pGradX * JPR
              else
                du = - dt * pGradX
              end if

              var%u(i, j, k) = var%u(i, j, k) + du
            end do
          end do
        end do
      else
        stop 'ERROR: wrong opt in correctorStep'
      end if
    else
      if(facprs /= 1.) stop 'ERROR: wrong facprs in explicit sub-step'
      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            if(topography) then
              ! TFC
              ! Compute values at cell edges.
              rhou = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k) &
                  &+ rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, j, k))
              pEdgeR = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) + pStratTFC(i &
                  &+ 1, j, k) / jac(i + 1, j, k))
              ! Compute pressure difference gradient component.
              if(k == 1 .and. zBoundary == "solid_wall") then
                dpUUEdgeR = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 1, 3) &
                    &* dp(i, j, k + 2) + jac(i + 1, j, k + 2) * met(i + 1, j, &
                    &k + 2, 1, 3) * dp(i + 1, j, k + 2))
                dpUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
                    &* dp(i, j, k + 1) + jac(i + 1, j, k + 1) * met(i + 1, j, &
                    &k + 1, 1, 3) * dp(i + 1, j, k + 1))
                dpEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * dp(i, j, &
                    &k) + jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) * dp(i &
                    &+ 1, j, k))
                pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((jac(i + 1, j, &
                    &k) * dp(i + 1, j, k) - jac(i, j, k) * dp(i, j, k)) / dx &
                    &+ (- dpUUEdgeR + 4.0 * dpUEdgeR - 3.0 * dpEdgeR) * 0.5 &
                    &/ dz)
              else if(k == nz .and. zBoundary == "solid_wall") then
                dpDDEdgeR = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 1, 3) &
                    &* dp(i, j, k - 2) + jac(i + 1, j, k - 2) * met(i + 1, j, &
                    &k - 2, 1, 3) * dp(i + 1, j, k - 2))
                dpDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, 3) &
                    &* dp(i, j, k - 1) + jac(i + 1, j, k - 1) * met(i + 1, j, &
                    &k - 1, 1, 3) * dp(i + 1, j, k - 1))
                dpEdgeR = 0.5 * (jac(i, j, k) * met(i, j, k, 1, 3) * dp(i, j, &
                    &k) + jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) * dp(i &
                    &+ 1, j, k))
                pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((jac(i + 1, j, &
                    &k) * dp(i + 1, j, k) - jac(i, j, k) * dp(i, j, k)) / dx &
                    &+ (dpDDEdgeR - 4.0 * dpDEdgeR + 3.0 * dpEdgeR) * 0.5 / dz)
              else
                dpUEdgeR = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 1, 3) &
                    &* dp(i, j, k + 1) + jac(i + 1, j, k + 1) * met(i + 1, j, &
                    &k + 1, 1, 3) * dp(i + 1, j, k + 1))
                dpDEdgeR = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 1, 3) &
                    &* dp(i, j, k - 1) + jac(i + 1, j, k - 1) * met(i + 1, j, &
                    &k - 1, 1, 3) * dp(i + 1, j, k - 1))
                pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((jac(i + 1, j, &
                    &k) * dp(i + 1, j, k) - jac(i, j, k) * dp(i, j, k)) / dx &
                    &+ (dpUEdgeR - dpDEdgeR) * 0.5 / dz)
              end if
            else
              rhou = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k)) &
                  &+ rhoStrat(k)

              pGradX = kappaInv * MaInv2 * pStrat(k) / rhou * (dp(i + 1, j, k) &
                  &- dp(i, j, k)) / dx
            end if

            du = - dt * pGradX

            var%u(i, j, k) = var%u(i, j, k) + du

            ! correct x-momentum tendency
            dMom(i, j, k, 1) = dMom(i, j, k, 1) + rhou * du / betaRK(RKstage)
          end do
        end do
      end do
    end if

    !--------------------------------------
    !         calc dv and v + dv
    !--------------------------------------

    if(timeScheme == "semiimplicit") then
      if(opt == "impl") then
        do k = 1, nz
          do j = 0, ny
            do i = 1, nx
              facv = 1.0

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  facv = facv + dt * 0.5 * (kv_hs(j, k) + kv_hs(j + 1, k))
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                if(topography) then
                  facv = facv + dt * 0.5 * (kr_sp_tfc(i, j, k) + kr_sp_tfc(i, &
                      &j + 1, k))
                else
                  facv = facv + dt * 0.5 * (kr_sp(j, k) + kr_sp(j + 1, k))
                end if
              end if

              facu = facv

              if(topography) then
                ! TFC
                ! Compute values at cell edges.
                rhov = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k) &
                    &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j + 1, k))
                pEdgeF = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    &+ pStratTFC(i, j + 1, k) / jac(i, j + 1, k))
                ! Compute pressure difference gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  dpUUEdgeF = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 2, 3) &
                      &* dp(i, j, k + 2) + jac(i, j + 1, k + 2) * met(i, j &
                      &+ 1, k + 2, 2, 3) * dp(i, j + 1, k + 2))
                  dpUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                      &* dp(i, j, k + 1) + jac(i, j + 1, k + 1) * met(i, j &
                      &+ 1, k + 1, 2, 3) * dp(i, j + 1, k + 1))
                  dpEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) * dp(i, &
                      &j, k) + jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) &
                      &* dp(i, j + 1, k))
                  pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((jac(i, j + 1, &
                      &k) * dp(i, j + 1, k) - jac(i, j, k) * dp(i, j, k)) / dy &
                      &+ (- dpUUEdgeF + 4.0 * dpUEdgeF - 3.0 * dpEdgeF) * 0.5 &
                      &/ dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  dpDDEdgeF = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 2, 3) &
                      &* dp(i, j, k - 2) + jac(i, j + 1, k - 2) * met(i, j &
                      &+ 1, k - 2, 2, 3) * dp(i, j + 1, k - 2))
                  dpDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) &
                      &* dp(i, j, k - 1) + jac(i, j + 1, k - 1) * met(i, j &
                      &+ 1, k - 1, 2, 3) * dp(i, j + 1, k - 1))
                  dpEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) * dp(i, &
                      &j, k) + jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) &
                      &* dp(i, j + 1, k))
                  pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((jac(i, j + 1, &
                      &k) * dp(i, j + 1, k) - jac(i, j, k) * dp(i, j, k)) / dy &
                      &+ (dpDDEdgeF - 4.0 * dpDEdgeF + 3.0 * dpEdgeF) * 0.5 &
                      &/ dz)
                else
                  dpUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                      &* dp(i, j, k + 1) + jac(i, j + 1, k + 1) * met(i, j &
                      &+ 1, k + 1, 2, 3) * dp(i, j + 1, k + 1))
                  dpDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) &
                      &* dp(i, j, k - 1) + jac(i, j + 1, k - 1) * met(i, j &
                      &+ 1, k - 1, 2, 3) * dp(i, j + 1, k - 1))
                  pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((jac(i, j + 1, &
                      &k) * dp(i, j + 1, k) - jac(i, j, k) * dp(i, j, k)) / dy &
                      &+ (dpUEdgeF - dpDEdgeF) * 0.5 / dz)
                end if
                ! Compute velocity correction.
                corY(i, j, k) = facprs * dt / facv * pGradY
                if(model == "compressible") then ! multiply with JP
                  JPF = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i, j + 1, &
                      &k) * var%P(i, j + 1, k))

                  dv = - JPF * corY(i, j, k)
                else
                  dv = - corY(i, j, k)
                end if
              else
                rhou00 = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k)) &
                    &+ rhoStrat(k)

                rhoum0 = 0.5 * (var%rho(i - 1, j, k) + var%rho(i, j, k)) &
                    &+ rhoStrat(k)

                rhou01 = 0.5 * (var%rho(i, j + 1, k) + var%rho(i + 1, j + 1, &
                    &k)) + rhoStrat(k)

                rhoum1 = 0.5 * (var%rho(i - 1, j + 1, k) + var%rho(i, j + 1, &
                    &k)) + rhoStrat(k)

                rhov = 0.5 * (var%rho(i, j + 1, k) + var%rho(i, j, k)) &
                    &+ rhoStrat(k)

                pGradX = kappaInv * MaInv2 * 0.25 * (pStrat(k) / rhou00 &
                    &* (dp(i + 1, j, k) - dp(i, j, k)) / dx + pStrat(k) &
                    &/ rhoum0 * (dp(i, j, k) - dp(i - 1, j, k)) / dx &
                    &+ pStrat(k) / rhou01 * (dp(i + 1, j + 1, k) - dp(i, j &
                    &+ 1, k)) / dx + pStrat(k) / rhoum1 * (dp(i, j + 1, k) &
                    &- dp(i - 1, j + 1, k)) / dx)

                pGradY = kappaInv * MaInv2 * pStrat(k) / rhov * (dp(i, j + 1, &
                    &k) - dp(i, j, k)) / dy

                f_cor_v = 0.5 * (f_cor_nd(j) + f_cor_nd(j + 1))

                dv = - facprs * dt / (facu * facv + (f_cor_v * dt) ** 2) * (- &
                    &f_cor_v * dt * pGradX + facu * pGradY)
              end if

              var%v(i, j, k) = var%v(i, j, k) + dv
            end do
          end do
        end do
      else if(opt == "expl") then
        if(facprs /= 1.) stop 'ERROR: wrong facprs in explicit sub-step'
        do k = 1, nz
          do j = 0, ny
            do i = 1, nx
              if(topography) then
                ! TFC
                ! Compute values at cell edges.
                rhov = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k) &
                    &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j + 1, k))
                pEdgeF = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    &+ pStratTFC(i, j + 1, k) / jac(i, j + 1, k))
                ! Compute pressure difference gradient component.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  dpUUEdgeF = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 2, 3) &
                      &* dp(i, j, k + 2) + jac(i, j + 1, k + 2) * met(i, j &
                      &+ 1, k + 2, 2, 3) * dp(i, j + 1, k + 2))
                  dpUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                      &* dp(i, j, k + 1) + jac(i, j + 1, k + 1) * met(i, j &
                      &+ 1, k + 1, 2, 3) * dp(i, j + 1, k + 1))
                  dpEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) * dp(i, &
                      &j, k) + jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) &
                      &* dp(i, j + 1, k))
                  pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((jac(i, j + 1, &
                      &k) * dp(i, j + 1, k) - jac(i, j, k) * dp(i, j, k)) / dy &
                      &+ (- dpUUEdgeF + 4.0 * dpUEdgeF - 3.0 * dpEdgeF) * 0.5 &
                      &/ dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  dpDDEdgeF = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 2, 3) &
                      &* dp(i, j, k - 2) + jac(i, j + 1, k - 2) * met(i, j &
                      &+ 1, k - 2, 2, 3) * dp(i, j + 1, k - 2))
                  dpDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) &
                      &* dp(i, j, k - 1) + jac(i, j + 1, k - 1) * met(i, j &
                      &+ 1, k - 1, 2, 3) * dp(i, j + 1, k - 1))
                  dpEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) * dp(i, &
                      &j, k) + jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) &
                      &* dp(i, j + 1, k))
                  pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((jac(i, j + 1, &
                      &k) * dp(i, j + 1, k) - jac(i, j, k) * dp(i, j, k)) / dy &
                      &+ (dpDDEdgeF - 4.0 * dpDEdgeF + 3.0 * dpEdgeF) * 0.5 &
                      &/ dz)
                else
                  dpUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                      &* dp(i, j, k + 1) + jac(i, j + 1, k + 1) * met(i, j &
                      &+ 1, k + 1, 2, 3) * dp(i, j + 1, k + 1))
                  dpDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) &
                      &* dp(i, j, k - 1) + jac(i, j + 1, k - 1) * met(i, j &
                      &+ 1, k - 1, 2, 3) * dp(i, j + 1, k - 1))
                  pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((jac(i, j + 1, &
                      &k) * dp(i, j + 1, k) - jac(i, j, k) * dp(i, j, k)) / dy &
                      &+ (dpUEdgeF - dpDEdgeF) * 0.5 / dz)
                end if
              else
                rhov = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k)) &
                    &+ rhoStrat(k)

                pGradY = kappaInv * MaInv2 * pStrat(k) / rhov * (dp(i, j + 1, &
                    &k) - dp(i, j, k)) / dy
              end if

              if(model == "compressible") then ! multiply with JP
                JPF = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i, j + 1, k) &
                    &* var%P(i, j + 1, k))

                dv = - dt * pGradY * JPF
              else
                dv = - dt * pGradY
              end if

              var%v(i, j, k) = var%v(i, j, k) + dv
            end do
          end do
        end do
      else
        stop 'ERROR: wrong opt in correctorStep'
      end if
    else
      if(facprs /= 1.) stop 'ERROR: wrong facprs in explicit sub-step'
      do k = 1, nz
        do j = 0, ny
          do i = 1, nx
            if(topography) then
              ! TFC
              ! Compute values at cell edges.
              rhov = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k) &
                  &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j + 1, k))
              pEdgeF = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) + pStratTFC(i, &
                  &j + 1, k) / jac(i, j + 1, k))
              ! Compute pressure difference gradient component.
              if(k == 1 .and. zBoundary == "solid_wall") then
                dpUUEdgeF = 0.5 * (jac(i, j, k + 2) * met(i, j, k + 2, 2, 3) &
                    &* dp(i, j, k + 2) + jac(i, j + 1, k + 2) * met(i, j + 1, &
                    &k + 2, 2, 3) * dp(i, j + 1, k + 2))
                dpUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                    &* dp(i, j, k + 1) + jac(i, j + 1, k + 1) * met(i, j + 1, &
                    &k + 1, 2, 3) * dp(i, j + 1, k + 1))
                dpEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) * dp(i, j, &
                    &k) + jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) * dp(i, j &
                    &+ 1, k))
                pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((jac(i, j + 1, &
                    &k) * dp(i, j + 1, k) - jac(i, j, k) * dp(i, j, k)) / dy &
                    &+ (- dpUUEdgeF + 4.0 * dpUEdgeF - 3.0 * dpEdgeF) * 0.5 &
                    &/ dz)
              else if(k == nz .and. zBoundary == "solid_wall") then
                dpDDEdgeF = 0.5 * (jac(i, j, k - 2) * met(i, j, k - 2, 2, 3) &
                    &* dp(i, j, k - 2) + jac(i, j + 1, k - 2) * met(i, j + 1, &
                    &k - 2, 2, 3) * dp(i, j + 1, k - 2))
                dpDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) &
                    &* dp(i, j, k - 1) + jac(i, j + 1, k - 1) * met(i, j + 1, &
                    &k - 1, 2, 3) * dp(i, j + 1, k - 1))
                dpEdgeF = 0.5 * (jac(i, j, k) * met(i, j, k, 2, 3) * dp(i, j, &
                    &k) + jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) * dp(i, j &
                    &+ 1, k))
                pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((jac(i, j + 1, &
                    &k) * dp(i, j + 1, k) - jac(i, j, k) * dp(i, j, k)) / dy &
                    &+ (dpDDEdgeF - 4.0 * dpDEdgeF + 3.0 * dpEdgeF) * 0.5 / dz)
              else
                dpUEdgeF = 0.5 * (jac(i, j, k + 1) * met(i, j, k + 1, 2, 3) &
                    &* dp(i, j, k + 1) + jac(i, j + 1, k + 1) * met(i, j + 1, &
                    &k + 1, 2, 3) * dp(i, j + 1, k + 1))
                dpDEdgeF = 0.5 * (jac(i, j, k - 1) * met(i, j, k - 1, 2, 3) &
                    &* dp(i, j, k - 1) + jac(i, j + 1, k - 1) * met(i, j + 1, &
                    &k - 1, 2, 3) * dp(i, j + 1, k - 1))
                pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((jac(i, j + 1, &
                    &k) * dp(i, j + 1, k) - jac(i, j, k) * dp(i, j, k)) / dy &
                    &+ (dpUEdgeF - dpDEdgeF) * 0.5 / dz)
              end if
            else
              rhov = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k)) &
                  &+ rhoStrat(k)

              pGradY = kappaInv * MaInv2 * pStrat(k) / rhov * (dp(i, j + 1, k) &
                  &- dp(i, j, k)) / dy
            end if

            dv = - dt * pGradY

            var%v(i, j, k) = var%v(i, j, k) + dv

            ! correct y-momentum tendency
            dMom(i, j, k, 2) = dMom(i, j, k, 2) + rhov * dv / betaRK(RKstage)
          end do
        end do
      end do
    end if

    !--------------------------------------
    !         calc w and  w + dw
    !--------------------------------------

    ! TFC
    select case(zBoundary)
    case("solid_wall")
      k0 = 1
      k1 = nz - 1
    case("periodic")
      k0 = 0
      k1 = nz
    case default
      stop "correctorStep: unknown case zBoundary."
    end select

    if(timeScheme == "semiimplicit") then
      if(opt == "impl") then
        ! solid wall implies zero change of w at the bottom and top
        ! TFC
        do k = k0, k1
          do j = 1, ny
            do i = 1, nx
              facw = 1.0

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  facw = facw + dt * 0.5 * (kw_hs(k) + kw_hs(k + 1))
                end if
              end if

              if(spongeLayer) then
                if(topography) then
                  facw = facw + dt * 0.5 * (kr_sp_w_tfc(i, j, k) &
                      &+ kr_sp_w_tfc(i, j, k + 1))
                else
                  facw = facw + dt * 0.5 * (kr_sp(j, k) + kr_sp(j, k + 1))
                end if
              end if

              if(topography) then
                ! TFC
                ! Compute values at cell edges.
                rhoStratEdgeU = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, &
                    &j, k + 1))
                rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1)) &
                    &+ rhoStratEdgeU
                pEdgeU = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    &+ pStratTFC(i, j, k + 1) / jac(i, j, k + 1))
                bvsstw = 0.5 * (bvsStratTFC(i, j, k) + bvsStratTFC(i, j, k + 1))
                dpREdgeU = 0.5 * (jac(i + 1, j, k) * met(i + 1, j, k, 3, 1) &
                    &* dp(i + 1, j, k) + jac(i + 1, j, k + 1) * met(i + 1, j, &
                    &k + 1, 3, 1) * dp(i + 1, j, k + 1))
                dpLEdgeU = 0.5 * (jac(i - 1, j, k) * met(i - 1, j, k, 3, 1) &
                    &* dp(i - 1, j, k) + jac(i - 1, j, k + 1) * met(i - 1, j, &
                    &k + 1, 3, 1) * dp(i - 1, j, k + 1))
                dpFEdgeU = 0.5 * (jac(i, j + 1, k) * met(i, j + 1, k, 3, 2) &
                    &* dp(i, j + 1, k) + jac(i, j + 1, k + 1) * met(i, j + 1, &
                    &k + 1, 3, 2) * dp(i, j + 1, k + 1))
                dpBEdgeU = 0.5 * (jac(i, j - 1, k) * met(i, j - 1, k, 3, 2) &
                    &* dp(i, j - 1, k) + jac(i, j - 1, k + 1) * met(i, j - 1, &
                    &k + 1, 3, 2) * dp(i, j - 1, k + 1))
                chris11EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                    &1) * dp(i, j, k) + pStratTFC(i, j, k + 1) * chris(i, j, k &
                    &+ 1, 1, 1) * dp(i, j, k + 1))
                chris22EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                    &2) * dp(i, j, k) + pStratTFC(i, j, k + 1) * chris(i, j, k &
                    &+ 1, 2, 2) * dp(i, j, k + 1))
                chris13EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                    &3) * met(i, j, k, 1, 3) * dp(i, j, k) + pStratTFC(i, j, k &
                    &+ 1) * chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, 1, 3) &
                    &* dp(i, j, k + 1))
                chris23EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                    &3) * met(i, j, k, 2, 3) * dp(i, j, k) + pStratTFC(i, j, k &
                    &+ 1) * chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, 2, 3) &
                    &* dp(i, j, k + 1))
                ! Compute pressure difference gradient component.
                pGradZ = kappaInv * MaInv2 / rhoEdge * pEdgeU * ((dpREdgeU &
                    &- dpLEdgeU) * 0.5 / dx + (dpFEdgeU - dpBEdgeU) * 0.5 / dy &
                    &+ (jac(i, j, k + 1) * met(i, j, k + 1, 3, 3) * dp(i, j, k &
                    &+ 1) - jac(i, j, k) * met(i, j, k, 3, 3) * dp(i, j, k)) &
                    &/ dz) + kappaInv * MaInv2 / rhoEdge * (chris11EdgeU &
                    &+ chris22EdgeU + 2.0 * (chris13EdgeU + chris23EdgeU))
                ! Compute velocity correction.
                if(model == "compressible") then
                  JPU = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i, j, k &
                      &+ 1) * var%P(i, j, k + 1))

                  dw = - facprs * dt / (facw + bvsstw * dt ** 2.0) * JPU &
                      &* pGradZ - 1.0 / (facw + bvsstw * dt ** 2.0) * bvsstw &
                      &* dt ** 2.0 * JPU * (0.25 * (met(i, j, k, 1, 3) &
                      &* (corX(i, j, k) + corX(i - 1, j, k)) + met(i, j, k &
                      &+ 1, 1, 3) * (corX(i, j, k + 1) + corX(i - 1, j, k &
                      &+ 1))) + 0.25 * (met(i, j, k, 2, 3) * (corY(i, j, k) &
                      &+ corY(i, j - 1, k)) + met(i, j, k + 1, 2, 3) &
                      &* (corY(i, j, k + 1) + corY(i, j - 1, k + 1))))
                else
                  dw = - facprs * dt / (facw + rhoStratEdgeU / rhoEdge &
                      &* bvsstw * dt ** 2.0) * pGradZ - 1.0 / (facw &
                      &+ rhoStratEdgeU / rhoEdge * bvsstw * dt ** 2.0) &
                      &* rhoStratEdgeU / rhoEdge * bvsstw * dt ** 2.0 * (0.25 &
                      &* (met(i, j, k, 1, 3) * (corX(i, j, k) + corX(i - 1, j, &
                      &k)) + met(i, j, k + 1, 1, 3) * (corX(i, j, k + 1) &
                      &+ corX(i - 1, j, k + 1))) + 0.25 * (met(i, j, k, 2, 3) &
                      &* (corY(i, j, k) + corY(i, j - 1, k)) + met(i, j, k &
                      &+ 1, 2, 3) * (corY(i, j, k + 1) + corY(i, j - 1, k &
                      &+ 1))))
                end if
              else
                pEdge = 0.5 * (pStrat(k + 1) + pStrat(k))
                pEdge_0 = 0.5 * (pStrat_0(k + 1) + pStrat_0(k))

                rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1)) &
                    &+ rhoStratTilde(k)

                bvsstw = 0.5 * (bvsStrat(k) + bvsStrat(k + 1))

                pGradZ = (dp(i, j, k + 1) - dp(i, j, k)) / dz

                dw = - facprs * dt * kappaInv * MaInv2 / (facw &
                    &+ rhoStratTilde(k) / rhoEdge * pEdge / pEdge_0 * bvsstw &
                    &* dt ** 2) * pEdge / rhoEdge * pGradz
              end if

              var%w(i, j, k) = var%w(i, j, k) + dw
            end do
          end do
        end do

      else if(opt == "expl") then
        ! solid wall implies zero change of w at the bottom and top
        if(facprs /= 1.) stop 'ERROR: wrong facprs in explicit sub-step'
        ! TFC
        do k = k0, k1
          do j = 1, ny
            do i = 1, nx
              if(topography) then
                ! TFC
                ! Compute values at cell edges.
                rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1) &
                    &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j, k + 1))
                pEdgeU = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    &+ pStratTFC(i, j, k + 1) / jac(i, j, k + 1))
                dpREdgeU = 0.5 * (jac(i + 1, j, k) * met(i + 1, j, k, 3, 1) &
                    &* dp(i + 1, j, k) + jac(i + 1, j, k + 1) * met(i + 1, j, &
                    &k + 1, 3, 1) * dp(i + 1, j, k + 1))
                dpLEdgeU = 0.5 * (jac(i - 1, j, k) * met(i - 1, j, k, 3, 1) &
                    &* dp(i - 1, j, k) + jac(i - 1, j, k + 1) * met(i - 1, j, &
                    &k + 1, 3, 1) * dp(i - 1, j, k + 1))
                dpFEdgeU = 0.5 * (jac(i, j + 1, k) * met(i, j + 1, k, 3, 2) &
                    &* dp(i, j + 1, k) + jac(i, j + 1, k + 1) * met(i, j + 1, &
                    &k + 1, 3, 2) * dp(i, j + 1, k + 1))
                dpBEdgeU = 0.5 * (jac(i, j - 1, k) * met(i, j - 1, k, 3, 2) &
                    &* dp(i, j - 1, k) + jac(i, j - 1, k + 1) * met(i, j - 1, &
                    &k + 1, 3, 2) * dp(i, j - 1, k + 1))
                chris11EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                    &1) * dp(i, j, k) + pStratTFC(i, j, k + 1) * chris(i, j, k &
                    &+ 1, 1, 1) * dp(i, j, k + 1))
                chris22EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                    &2) * dp(i, j, k) + pStratTFC(i, j, k + 1) * chris(i, j, k &
                    &+ 1, 2, 2) * dp(i, j, k + 1))
                chris13EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                    &3) * met(i, j, k, 1, 3) * dp(i, j, k) + pStratTFC(i, j, k &
                    &+ 1) * chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, 1, 3) &
                    &* dp(i, j, k + 1))
                chris23EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                    &3) * met(i, j, k, 2, 3) * dp(i, j, k) + pStratTFC(i, j, k &
                    &+ 1) * chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, 2, 3) &
                    &* dp(i, j, k + 1))
                ! Compute pressure difference gradient component.
                pGradZ = kappaInv * MaInv2 / rhoEdge * pEdgeU * ((dpREdgeU &
                    &- dpLEdgeU) * 0.5 / dx + (dpFEdgeU - dpBEdgeU) * 0.5 / dy &
                    &+ (jac(i, j, k + 1) * met(i, j, k + 1, 3, 3) * dp(i, j, k &
                    &+ 1) - jac(i, j, k) * met(i, j, k, 3, 3) * dp(i, j, k)) &
                    &/ dz) + kappaInv * MaInv2 / rhoEdge * (chris11EdgeU &
                    &+ chris22EdgeU + 2.0 * (chris13EdgeU + chris23EdgeU))
                ! Correct vertical velocity.
                if(model == "compressible") then
                  JPU = 0.5 * (jac(i, j, k) * var%P(i, j, k) + jac(i, j, k &
                      &+ 1) * var%P(i, j, k + 1))

                  dw = - dt * pGradZ * JPU
                else
                  dw = - dt * pGradZ
                end if
              else
                pEdge = pStratTilde(k)

                rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1)) &
                    &+ rhoStratTilde(k)

                pGradZ = (dp(i, j, k + 1) - dp(i, j, k)) / dz

                dw = - dt * kappaInv * MaInv2 * pEdge / rhoEdge * pGradz
              end if

              var%w(i, j, k) = var%w(i, j, k) + dw
            end do
          end do
        end do

      else
        stop 'ERROR: wrong opt in correctorStep'
      end if
    else
      if(facprs /= 1.) stop 'ERROR: wrong facprs in explicit sub-step'
      ! solid wall implies zero change of w at the bottom and top
      do k = k0, k1
        do j = 1, ny
          do i = 1, nx
            if(topography) then
              ! Compute values at cell edges.
              rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1) &
                  &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j, k + 1))
              pEdgeU = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) + pStratTFC(i, &
                  &j, k + 1) / jac(i, j, k + 1))
              dpREdgeU = 0.5 * (jac(i + 1, j, k) * met(i + 1, j, k, 3, 1) &
                  &* dp(i + 1, j, k) + jac(i + 1, j, k + 1) * met(i + 1, j, k &
                  &+ 1, 3, 1) * dp(i + 1, j, k + 1))
              dpLEdgeU = 0.5 * (jac(i - 1, j, k) * met(i - 1, j, k, 3, 1) &
                  &* dp(i - 1, j, k) + jac(i - 1, j, k + 1) * met(i - 1, j, k &
                  &+ 1, 3, 1) * dp(i - 1, j, k + 1))
              dpFEdgeU = 0.5 * (jac(i, j + 1, k) * met(i, j + 1, k, 3, 2) &
                  &* dp(i, j + 1, k) + jac(i, j + 1, k + 1) * met(i, j + 1, k &
                  &+ 1, 3, 2) * dp(i, j + 1, k + 1))
              dpBEdgeU = 0.5 * (jac(i, j - 1, k) * met(i, j - 1, k, 3, 2) &
                  &* dp(i, j - 1, k) + jac(i, j - 1, k + 1) * met(i, j - 1, k &
                  &+ 1, 3, 2) * dp(i, j - 1, k + 1))
              chris11EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, 1) &
                  &* dp(i, j, k) + pStratTFC(i, j, k + 1) * chris(i, j, k + 1, &
                  &1, 1) * dp(i, j, k + 1))
              chris22EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, 2) &
                  &* dp(i, j, k) + pStratTFC(i, j, k + 1) * chris(i, j, k + 1, &
                  &2, 2) * dp(i, j, k + 1))
              chris13EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, 3) &
                  &* met(i, j, k, 1, 3) * dp(i, j, k) + pStratTFC(i, j, k + 1) &
                  &* chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, 1, 3) * dp(i, &
                  &j, k + 1))
              chris23EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, 3) &
                  &* met(i, j, k, 2, 3) * dp(i, j, k) + pStratTFC(i, j, k + 1) &
                  &* chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, 2, 3) * dp(i, &
                  &j, k + 1))
              ! Compute pressure difference gradient component.
              pGradZ = kappaInv * MaInv2 / rhoEdge * pEdgeU * ((dpREdgeU &
                  &- dpLEdgeU) * 0.5 / dx + (dpFEdgeU - dpBEdgeU) * 0.5 / dy &
                  &+ (jac(i, j, k + 1) * met(i, j, k + 1, 3, 3) * dp(i, j, k &
                  &+ 1) - jac(i, j, k) * met(i, j, k, 3, 3) * dp(i, j, k)) &
                  &/ dz) + kappaInv * MaInv2 / rhoEdge * (chris11EdgeU &
                  &+ chris22EdgeU + 2.0 * (chris13EdgeU + chris23EdgeU))
            else
              pEdge = pStratTilde(k)

              rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1)) &
                  &+ rhoStratTilde(k)

              pGradZ = kappaInv * MaInv2 * pEdge / rhoEdge * (dp(i, j, k + 1) &
                  &- dp(i, j, k)) / dz
            end if

            dw = - dt * pGradz

            var%w(i, j, k) = var%w(i, j, k) + dw

            ! correct z-momentum tendency
            dMom(i, j, k, 3) = dMom(i, j, k, 3) + rhoEdge * dw / betaRK(RKstage)
          end do
        end do
      end do

    end if

    !------------------------------------------------------------------
    !         calc rhop and rhop + drhop (only for implicit time step)
    !------------------------------------------------------------------

    if(timeScheme == "semiimplicit") then
      if(opt == "impl") then
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              facw = 1.0

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  facw = facw + dt * kw_hs(k)
                end if
              end if

              if(spongeLayer) then
                if(topography) then
                  facw = facw + dt * kr_sp_w_tfc(i, j, k)
                else
                  facw = facw + dt * kr_sp(j, k)
                end if
              end if

              if(topography) then
                ! Compute P coefficients.
                pEdgeU = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    &+ pStratTFC(i, j, k + 1) / jac(i, j, k + 1))
                pEdgeD = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                    &+ pStratTFC(i, j, k - 1) / jac(i, j, k - 1))
                ! Compute density coefficients.
                rhow0 = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1) &
                    &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j, k + 1))
                rhowm = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k - 1) &
                    &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j, k - 1))
                rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)
                ! Interpolate pressure differences.
                dpREdgeU = 0.5 * (jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) &
                    &* dp(i + 1, j, k) + jac(i + 1, j, k + 1) * met(i + 1, j, &
                    &k + 1, 1, 3) * dp(i + 1, j, k + 1))
                dpLEdgeU = 0.5 * (jac(i - 1, j, k) * met(i - 1, j, k, 1, 3) &
                    &* dp(i - 1, j, k) + jac(i - 1, j, k + 1) * met(i - 1, j, &
                    &k + 1, 1, 3) * dp(i - 1, j, k + 1))
                dpREdgeD = 0.5 * (jac(i + 1, j, k) * met(i + 1, j, k, 1, 3) &
                    &* dp(i + 1, j, k) + jac(i + 1, j, k - 1) * met(i + 1, j, &
                    &k - 1, 1, 3) * dp(i + 1, j, k - 1))
                dpLEdgeD = 0.5 * (jac(i - 1, j, k) * met(i - 1, j, k, 1, 3) &
                    &* dp(i - 1, j, k) + jac(i - 1, j, k - 1) * met(i - 1, j, &
                    &k - 1, 1, 3) * dp(i - 1, j, k - 1))
                dpFEdgeU = 0.5 * (jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) &
                    &* dp(i, j + 1, k) + jac(i, j + 1, k + 1) * met(i, j + 1, &
                    &k + 1, 2, 3) * dp(i, j + 1, k + 1))
                dpBEdgeU = 0.5 * (jac(i, j - 1, k) * met(i, j - 1, k, 2, 3) &
                    &* dp(i, j - 1, k) + jac(i, j - 1, k + 1) * met(i, j - 1, &
                    &k + 1, 2, 3) * dp(i, j - 1, k + 1))
                dpFEdgeD = 0.5 * (jac(i, j + 1, k) * met(i, j + 1, k, 2, 3) &
                    &* dp(i, j + 1, k) + jac(i, j + 1, k - 1) * met(i, j + 1, &
                    &k - 1, 2, 3) * dp(i, j + 1, k - 1))
                dpBEdgeD = 0.5 * (jac(i, j - 1, k) * met(i, j - 1, k, 2, 3) &
                    &* dp(i, j - 1, k) + jac(i, j - 1, k - 1) * met(i, j - 1, &
                    &k - 1, 2, 3) * dp(i, j - 1, k - 1))
                ! Interpolate Christoffel symbols.
                chris11EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                    &1) * dp(i, j, k) + pStratTFC(i, j, k + 1) * chris(i, j, k &
                    &+ 1, 1, 1) * dp(i, j, k + 1))
                chris11EdgeD = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                    &1) * dp(i, j, k) + pStratTFC(i, j, k - 1) * chris(i, j, k &
                    &- 1, 1, 1) * dp(i, j, k - 1))
                chris22EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                    &2) * dp(i, j, k) + pStratTFC(i, j, k + 1) * chris(i, j, k &
                    &+ 1, 2, 2) * dp(i, j, k + 1))
                chris22EdgeD = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                    &2) * dp(i, j, k) + pStratTFC(i, j, k - 1) * chris(i, j, k &
                    &- 1, 2, 2) * dp(i, j, k - 1))
                chris13EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                    &3) * met(i, j, k, 1, 3) * dp(i, j, k) + pStratTFC(i, j, k &
                    &+ 1) * chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, 1, 3) &
                    &* dp(i, j, k + 1))
                chris13EdgeD = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 1, &
                    &3) * met(i, j, k, 1, 3) * dp(i, j, k) + pStratTFC(i, j, k &
                    &- 1) * chris(i, j, k - 1, 1, 3) * met(i, j, k - 1, 1, 3) &
                    &* dp(i, j, k - 1))
                chris23EdgeU = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                    &3) * met(i, j, k, 2, 3) * dp(i, j, k) + pStratTFC(i, j, k &
                    &+ 1) * chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, 2, 3) &
                    &* dp(i, j, k + 1))
                chris23EdgeD = 0.5 * (pStratTFC(i, j, k) * chris(i, j, k, 2, &
                    &3) * met(i, j, k, 2, 3) * dp(i, j, k) + pStratTFC(i, j, k &
                    &- 1) * chris(i, j, k - 1, 2, 3) * met(i, j, k - 1, 2, 3) &
                    &* dp(i, j, k - 1))
                ! Compute pressure difference gradients.
                pGradZEdgeU = kappaInv * MaInv2 * pEdgeU / rhow0 * (0.5 &
                    &* (dpREdgeU - dpLEdgeU) / dx + 0.5 * (dpFEdgeU &
                    &- dpBEdgeU) / dy + (jac(i, j, k + 1) * met(i, j, k + 1, &
                    &3, 3) * dp(i, j, k + 1) - jac(i, j, k) * met(i, j, k, 3, &
                    &3) * dp(i, j, k)) / dz) + kappaInv * MaInv2 / rhow0 &
                    &* (chris11EdgeU + chris22EdgeU + 2.0 * (chris13EdgeU &
                    &+ chris23EdgeU))
                pGradZEdgeD = kappaInv * MaInv2 * pEdgeD / rhowm * (0.5 &
                    &* (dpREdgeD - dpLEdgeD) / dx + 0.5 * (dpFEdgeD &
                    &- dpBEdgeD) / dy + (jac(i, j, k) * met(i, j, k, 3, 3) &
                    &* dp(i, j, k) - jac(i, j, k - 1) * met(i, j, k - 1, 3, 3) &
                    &* dp(i, j, k - 1)) / dz) + kappaInv * MaInv2 / rhowm &
                    &* (chris11EdgeD + chris22EdgeD + 2.0 * (chris13EdgeD &
                    &+ chris23EdgeD))
                ! Adjust at boundaries.
                if(k == 1 .and. zBoundary == "solid_wall") then
                  pGradZEdgeD = 0.0
                else if(k == nz .and. zBoundary == "solid_wall") then
                  pGradZEdgeU = 0.0
                end if
                ! Interpolate.
                pGradZ = 0.5 * (pGradZEdgeU + pGradZEdgeD)
                ! Compute buoyancy correction.
                if(model == "compressible") then
                  db = - 1.0 / (facw + bvsStratTFC(i, j, k) * dt ** 2.0) * (- &
                      &bvsStratTFC(i, j, k) * facprs * dt ** 2.0 * jac(i, j, &
                      &k) * pGradZ + bvsStratTFC(i, j, k) * dt * jac(i, j, k) &
                      &* facw * 0.5 * (met(i, j, k, 1, 3) * (corX(i, j, k) &
                      &+ corX(i - 1, j, k)) + met(i, j, k, 2, 3) * (corY(i, j, &
                      &k) + corY(i, j - 1, k))))
                else
                  db = - 1.0 / (facw + rhoStratTFC(i, j, k) / rho &
                      &* bvsStratTFC(i, j, k) * dt ** 2.0) * (- rhoStratTFC(i, &
                      &j, k) / rho * bvsStratTFC(i, j, k) * facprs * dt ** 2.0 &
                      &* jac(i, j, k) * pGradZ + rhoStratTFC(i, j, k) / rho &
                      &* bvsStratTFC(i, j, k) * dt * jac(i, j, k) * facw * 0.5 &
                      &* (met(i, j, k, 1, 3) * (corX(i, j, k) + corX(i - 1, j, &
                      &k)) + met(i, j, k, 2, 3) * (corY(i, j, k) + corY(i, j &
                      &- 1, k))))
                end if
              else
                rho = var%rho(i, j, k) + rhoStrat(k)

                rhowm = 0.5 * (var%rho(i, j, k - 1) + var%rho(i, j, k)) &
                    &+ rhoStratTilde(k - 1)

                rhow0 = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1)) &
                    &+ rhoStratTilde(k)

                if(k == 1 .and. zBoundary == "solid_wall") then
                  pGradZ = kappaInv * MaInv2 * 0.5 * (pStratTilde(k) / rhow0 &
                      &* (dp(i, j, k + 1) - dp(i, j, k)) / dz)
                else if(k == nz .and. zBoundary == "solid_wall") then
                  pGradZ = kappaInv * MaInv2 * 0.5 * (pStratTilde(k - 1) &
                      &/ rhowm * (dp(i, j, k) - dp(i, j, k - 1)) / dz)
                else
                  pGradZ = kappaInv * MaInv2 * 0.5 * (pStratTilde(k) / rhow0 &
                      &* (dp(i, j, k + 1) - dp(i, j, k)) / dz + pStratTilde(k &
                      &- 1) / rhowm * (dp(i, j, k) - dp(i, j, k - 1)) / dz)
                end if

                db = rhoStrat(k) / rho * pStrat(k) / pStrat_0(k) * bvsStrat(k) &
                    &* facprs * dt ** 2 / (facw + rhoStrat(k) / rho &
                    &* pStrat(k) / pStrat_0(k) * bvsStrat(k) * dt ** 2) * pGradz
              end if

              var%rhop(i, j, k) = var%rhop(i, j, k) - rho / g_ndim * db
            end do
          end do
        end do

      end if
    end if

    if(timeScheme == "semiimplicit") then
      if(opt == "impl") then
        kr_sp = kr_sp / facray
        alprlx = alprlx / facray
        if(topography) then
          kr_sp_tfc = kr_sp_tfc / facray
          kr_sp_w_tfc = kr_sp_w_tfc / facray
        end if
      end if
    end if

  end subroutine correctorStep

  !----------------------------------------------------------------------

  function getIndex(i, j, k)
    !----------------------------------------------------
    ! compute matrix index l from spatial indicies i,j,k
    !----------------------------------------------------

    ! in/out variables
    integer :: getIndex
    integer, intent(in) :: i, j, k

    getIndex = (k - 1) * nx * ny + (j - 1) * nx + i

  end function getIndex

  !--------------------------------------------------------------------

  subroutine init_poisson

    !-----------------------------------
    ! allocate poisson module variables
    !-----------------------------------

    ! local variables
    integer :: allocstat

    ! allocate fields
    allocate(dp(0:nx + 1, 0:ny + 1, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "init_poisson: alloc failed"

    allocate(sol_old1(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "init_poisson: alloc failed"

    allocate(sol_old2(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "init_poisson: alloc failed"

    allocate(p_pred(1:nx, 1:ny, 1:nz), stat = allocstat)
    if(allocstat /= 0) stop "init_poisson: alloc failed"

    ! Initial status of Boussinesq tensor elements.
    if(model == "Boussinesq") then
      expEle = .false.
      impEle = .false.
    end if

  end subroutine init_poisson

  !----------------------------------------------------------------------

  subroutine terminate_poisson
    !-----------------------------------
    ! deallocate poisson module variables
    !-----------------------------------

    ! local variables
    integer :: allocstat

    ! deallocate fields
    deallocate(dp, stat = allocstat)
    if(allocstat /= 0) stop "init_poisson: dealloc failed"

    deallocate(sol_old1, stat = allocstat)
    if(allocstat /= 0) stop "init_poisson: dealloc failed"

    deallocate(sol_old2, stat = allocstat)
    if(allocstat /= 0) stop "init_poisson: dealloc failed"

    deallocate(p_pred, stat = allocstat)
    if(allocstat /= 0) stop "init_poisson: dealloc failed"

  end subroutine terminate_poisson

  !----------------------------------------------------------------------

  subroutine calculate_heating(var, flux, heat)

    !-----------------------------------------------------------------
    ! calculate heating in the divergence constraint
    ! supplemented by 'heating' due to turbulent and GW entropy fluxes
    !-----------------------------------------------------------------

    type(var_type), intent(in) :: var
    type(flux_type), intent(in) :: flux

    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(out) :: heat !, term1, term2

    ! local variables
    integer :: i, j, k, khmax
    real :: rho, the, theta_bar_0, rho_e
    real, dimension(1:ny, 0:nz) :: tau_relax_i ! inverse scaled relaxation
    !FS                 ! time for barocl. l. c.
    real :: tau_jet_sc, tau_relax_sc ! Klein scaling for relax param.

    real, dimension(1:nx, 1:ny, 1:nz) :: tau_relax_i_tfc

    real :: ymax, ymin, yloc, r, theta, dTheta_dim, delX, delZ, x_dim, z_dim
    integer :: j00, i00

    real, dimension(1:nz) :: sum_local, sum_global

    heat = 0.0

    ymin = ly_dim(0) / lRef
    ymax = ly_dim(1) / lRef

    i00 = is + nbx - 1
    j00 = js + nby - 1

    ! No heating in TFC (FJApr2023)
    if(topography .and. model /= "compressible") return

    !-------------------------------------------------------
    ! calculate the environment-induced negative (!) heating
    !-------------------------------------------------------

    if(testCase == "hotBubble_heat") then
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            x_dim = x(i + i00) * lRef ! dimensional lenghts
            z_dim = z(k) * lRef

            delX = (x_dim) / 1000.
            delZ = (z_dim - 3000.) / 1000.

            r = sqrt(delX ** 2 + delZ ** 2) ! scaled radius

            if(r <= 1.0) then ! inside bubble

              dTheta_dim = 0.5 * (cos(0.5 * pi * r)) ** 2
              theta = dTheta_dim / thetaRef

              heat(i, j, k) = - theta

            else
              heat(i, j, k) = 0.
            end if

          end do
        end do
      end do

    end if

    if(TestCase == "heatedLayer") then

      do k = 1, nz
        heat(:, :, k) = - 0.5 / thetaRef * exp(- (z(k) - 3000. / lRef) ** 2 &
            &/ (1000. / lRef) ** 2)
      end do
    end if

    if(testCase == "hotBubble_heatedLayer") then
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            x_dim = x(i + i00) * lRef ! dimensional lenghts
            z_dim = z(k) * lRef

            delX = (x_dim) / 1000.
            delZ = (z_dim - 3000.) / 1000.

            r = sqrt(delX ** 2 + delZ ** 2) ! scaled radius

            if(r <= 1.0) then ! inside bubble

              dTheta_dim = 0.5 * ((cos(0.5 * pi * r)) ** 2 + exp(- (z(k) &
                  &- 3000. / lRef) ** 2 / (1000. / lRef) ** 2))
              theta = dTheta_dim / thetaRef

              heat(i, j, k) = - theta

            else
              heat(i, j, k) = - 0.5 / thetaRef * exp(- (z(k) - 3000. / lRef) &
                  &** 2 / (1000. / lRef) ** 2)
            end if

          end do
        end do
      end do

    end if

    if((TestCase == "baroclinic_LC") .or. (TestCase == "baroclinic_ID")) then

      if(background == "HeldSuarez") then
        if(topography) then
          tau_relax_i_tfc = kt_hs_tfc
        else
          do k = 1, nz
            do j = 1, ny
              tau_relax_i(j, k) = kt_hs(j, k)
            end do
          end do
        end if
      else
        do k = 1, nz
          if(referenceQuantities == "SI") then
            tau_relax_sc = tau_relax !tau_z(k)  !tau_relax
            tau_jet_sc = tau_jet
          else
            tau_relax_sc = tau_relax / tref !tau_z(k)  !tau_relax/tref
            tau_jet_sc = tau_jet / tref
          end if

          do j = 1, ny
            yloc = y(j + j00)

            if(yloc > 0.5 * (ymax + ymin)) then
              ! meridionally dependent tau_sc

              tau_relax_sc = tau_relax / tref + (tau_relax_low / tref &
                  &- tau_relax / tref) * (1.0 - 0.5 * (tanh((yloc / ymax &
                  &- 0.1) / 0.02) - tanh((yloc / ymax - 0.9) / 0.02)))
            else
              tau_relax_sc = tau_relax / tref + (tau_relax_low / tref &
                  &- tau_relax / tref) * (1.0 - 0.5 * (tanh(((- 1) * yloc &
                  &/ ymax - 0.1) / 0.02) - tanh(((- 1) * yloc / ymax - 0.9) &
                  &/ 0.02)))
            end if

            if(topography) then
              tau_relax_i_tfc(:, j, k) = 1.0 / tau_relax_sc
            else
              tau_relax_i(j, k) = 1. / (tau_relax_sc)
            end if
          end do
        end do
      end if

      if(dens_relax) tau_relax_i = 0.0

      if(master) then
        if(topography) then
          print *, ""
          print *, " Poisson Solver, Thermal Relaxation is on: "
          print *, " Relaxation factor: div = - rho * (theta - theta_e) / tau: "
          print *, " tau = ", tref / tau_relax_i_tfc(1, 1, 1), " s"
        else
          print *, ""
          print *, " Poisson Solver, Thermal Relaxation is on: "
          print *, " Relaxation factor: Div = - rho(Th - Th_e)/tau: "
          print *, "tau = ", tref / tau_relax_i(1, 1), " s"
        end if
      end if

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            if(topography) then
              rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)
            else
              rho = var%rho(i, j, k) + rhoStrat(k)
            end if

            if(topography) then
              the = pStratTFC(i, j, k) / rho
            else
              the = Pstrat(k) / rho
            end if

            heat(i, j, k) = the - the_env_pp(i, j, k)

          end do
        end do
      end do

      ! heating

      if(topography) then
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)
              heat(i, j, k) = rho * heat(i, j, k) * tau_relax_i_tfc(i, j, k)
            end do
          end do
        end do
      else
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              rho = var%rho(i, j, k) + rhoStrat(k)

              heat(i, j, k) = rho * heat(i, j, k) * tau_relax_i(j, k)

            end do
          end do
        end do
      end if

    end if

    !------------------------------------------------------------------
    ! supplement by negative (!) heating due to molecular and turbulent
    ! diffusivity
    !------------------------------------------------------------------

    if(TurbScheme) then
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            if(topography) then
              ! Divide by Jacobian for divergence below.
              rho = (var%rho(i, j, k) + rhoStratTFC(i, j, k)) / jac(i, j, k)
            else
              rho = var%rho(i, j, k) + rhoStrat(k)
            end if

            heat(i, j, k) = heat(i, j, k) + rho * ((flux%theta(i, j, k, 1) &
                &- flux%theta(i - 1, j, k, 1)) / dx + (flux%theta(i, j, k, 2) &
                &- flux%theta(i, j - 1, k, 2)) / dy + (flux%theta(i, j, k, 3) &
                &- flux%theta(i, j, k - 1, 3)) / dz)
          end do
        end do
      end do
    end if

    !-------------------------------------------------------------------
    ! supplement by negative (!) heating due GW entropy-flux convergence
    !-------------------------------------------------------------------

    if(raytracer) heat(:, :, :) = heat(:, :, :) + var%GWH(:, :, :)

    if(spongeLayer .and. .not. unifiedSponge) then
      if(topography) then
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              if(heightTFC(i, j, k) > zSponge + 0.5 * (lz(1) - zSponge)) then
                heat(i, j, k) = 0.0
              else if(heightTFC(i, j, k) >= zSponge) then
                heat(i, j, k) = heat(i, j, k) * cos(0.5 * pi * (heightTFC(i, &
                    &j, k) - zSponge) / (0.5 * (lz(1) - zSponge))) ** 2.0
              end if
            end do
          end do
        end do
      else
        khmax = ksponge + int((nz - kSponge) / 2)

        do k = kSponge, khmax
          heat(:, :, k) = heat(:, :, k) * cos(0.5 * pi * (z(k) - zSponge) &
              &/ (z(khmax) - zSponge)) ** 2
        end do

        do k = khmax + 1, nz
          heat(:, :, k) = 0.

        end do
      end if
    end if

  end subroutine calculate_heating

  !==============================================================

  subroutine val_PsIn(var, dt, opt, facray)

    ! calculates the matrix values for the pressure solver
    ! the solver solves for dt * dp, hence no dt in the matrix elements

    type(var_type), intent(in) :: var
    real, intent(in) :: dt, facray

    ! facray multiplies the Rayleigh-damping terms so that they are only
    ! handled in the implicit time stepping (sponge and immersed boundary)

    ! opt = expl =>
    ! pressure solver for explicit problem and corresponding correction
    ! of the winds
    ! opt = impl =>
    ! pressure solver for implicit problem and corresponding correction
    ! of the winds and density fluctuations
    character(len = *), intent(in) :: opt

    ! local variables
    real :: dx2, dy2, dz2, dxy
    real :: pStratU, pStratD, pStratU_0, pStratD_0, rhoEdge
    real :: AL, AR, AB, AF, AD, AU, ALB, ALF, ARB, ARF, AC, ACH, ACV
    real :: facu, facv, facw, facr
    real :: acontr
    real :: bvsstw

    ! non-dimensional Corilois parameter (= inverse Rossby number)
    real, dimension(0:ny + 1) :: f_cor_nd

    integer :: i0, j0, i, j, k
    integer :: index_count_hypre

    real :: ymin, ymax, yloc

    real :: f_cor_v
    real :: fcscal, fcscal_u, fcscal_d

    ! SK compressible
    real :: dPdPi

    ! TFC FJ
    real :: ARU, ARD, ALU, ALD, AFU, AFD, ABU, ABD
    real :: AUU, ADD, ARUU, ARDD, ALUU, ALDD, AFUU, AFDD, ABUU, ABDD
    real :: fcscal_uu, fcscal_dd
    real :: jacInv
    real :: pEdgeRDiv, pEdgeLDiv, pEdgeFDiv, pEdgeBDiv, pEdgeUDiv, pEdgeDDiv
    real :: pEdgeRGra, pEdgeLGra, pEdgeFGra, pEdgeBGra, pEdgeUGra, pEdgeDGra
    real :: rhoEdgeR, rhoEdgeL, rhoEdgeF, rhoEdgeB, rhoEdgeU, rhoEdgeD
    real :: pUEdgeRGra, pUEdgeLGra, pUEdgeFGra, pUEdgeBGra, pDEdgeRGra, &
        &pDEdgeLGra, pDEdgeFGra, pDEdgeBGra
    real :: rhoStratEdgeR, rhoStratEdgeL, rhoStratEdgeF, rhoStratEdgeB, &
        &rhoStratEdgeU, rhoStratEdgeD, rhoUEdgeR, rhoUEdgeL, rhoUEdgeF, &
        &rhoUEdgeB, rhoDEdgeR, rhoDEdgeL, rhoDEdgeF, rhoDEdgeB
    real :: bvsStratEdgeU, bvsStratEdgeD
    real :: facEdgeR, facEdgeL, facEdgeF, facEdgeB, facEdgeU, facEdgeD, &
        &facUEdgeR, facUEdgeL, facUEdgeF, facUEdgeB, facDEdgeR, facDEdgeL, &
        &facDEdgeF, facDEdgeB
    real :: impHorEdgeR, impHorEdgeL, impHorEdgeF, impHorEdgeB, impHorUEdgeR, &
        &impHorUEdgeL, impHorUEdgeF, impHorUEdgeB, impHorDEdgeR, impHorDEdgeL, &
        &impHorDEdgeF, impHorDEdgeB, impVerEdgeU, impVerEdgeD
    real :: gEdgeR, gEdgeL, gEdgeF, gEdgeB, gEdgeU, gEdgeD, gUEdgeR, gUEdgeL, &
        &gUEdgeF, gUEdgeB, gDEdgeR, gDEdgeL, gDEdgeF, gDEdgeB

    if(corset == 'periodic') then
      ymax = ly_dim(1) / lRef
      ymin = ly_dim(0) / lRef

      j0 = js + nby - 1

      do j = 0, ny + 1
        yloc = y(j + j0)

        f_cor_nd(j) = - 4. * pi / 8.64e4 * tRef * cos(2. * pi * (yloc - ymin) &
            &/ (ymax - ymin))
      end do
    else if(corset == 'constant') then
      f_cor_nd(0:ny + 1) = f_Coriolis_dim * tRef
    else
      stop 'ERROR: wrong corset'
    end if

    ! auxiliary variables
    dx2 = 1.0 / dx ** 2
    dy2 = 1.0 / dy ** 2
    dz2 = 1.0 / dz ** 2
    dxy = 1.0 / (dx * dy)

    i0 = is + nbx - 1
    j0 = js + nby - 1

    !---------------------------------
    !         Loop over field
    !---------------------------------

    if(opt == "expl") then
      ! TFC FJ
      if(timeScheme == "semiimplicit" .and. .not. topography) then
        do k = 1, nz
          fcscal = sqrt(Pstrat(k) ** 2 / rhoStrat(k))
          fcscal_u = sqrt(Pstrat(k + 1) ** 2 / rhoStrat(k + 1))
          fcscal_d = sqrt(Pstrat(k - 1) ** 2 / rhoStrat(k - 1))

          do j = 1, ny
            do i = 1, nx
              ! ------------------ A(i+1,j,k) ------------------

              rhoEdge = 0.5 * (var%rho(i + 1, j, k) + var%rho(i, j, k)) &
                  &+ rhoStrat(k)

              AR = dx2 * pStrat(k) ** 2 / rhoEdge
              if(pressureScaling) then
                AR = AR / Pstrat(k)
              end if

              ! ------------------- A(i-1,j,k) --------------------

              rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i - 1, j, k)) &
                  &+ rhoStrat(k)

              AL = dx2 * pStrat(k) ** 2 / rhoEdge
              if(pressureScaling) then
                AL = AL / Pstrat(k)
              end if

              ! -------------------- A(i,j+1,k) ----------------------

              rhoEdge = 0.5 * (var%rho(i, j + 1, k) + var%rho(i, j, k)) &
                  &+ rhoStrat(k)

              AF = dy2 * pStrat(k) ** 2 / rhoEdge

              if(pressureScaling) then
                AF = AF / Pstrat(k)
              end if

              ! --------------------- A(i,j-1,k) -------------------

              rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i, j - 1, k)) &
                  &+ rhoStrat(k)

              AB = dy2 * pStrat(k) ** 2 / rhoEdge
              if(pressureScaling) then
                AB = AB / Pstrat(k)
              end if

              ! ---------------------- A(i,j,k+1) ------------------

              ! TFC FJ
              if(k == nz .and. zBoundary == "solid_wall") then
                AU = 0.0
              else
                rhoEdge = 0.5 * (var%rho(i, j, k + 1) + var%rho(i, j, k)) &
                    &+ rhoStratTilde(k)

                pStratU = 0.5 * (pStrat(k + 1) + pStrat(k))
                AU = dz2 * pStratU ** 2 / rhoEdge
              end if
              if(pressureScaling) then
                AU = AU / PstratTilde(k)
              end if

              ! ----------------------- A(i,j,k-1) -----------------

              ! TFC FJ
              if(k == 1 .and. zBoundary == "solid_wall") then
                AD = 0.0
              else
                rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k - 1)) &
                    &+ rhoStratTilde(k - 1)

                pStratD = 0.5 * (pStrat(k) + pStrat(k - 1))
                AD = dz2 * pStratD ** 2 / rhoEdge
              end if
              if(pressureScaling) then
                AD = AD / PstratTilde(k - 1)
              end if

              ! ----------------------- A(i,j,k) -------------------

              AC = - AR - AL - AF - AB - AU - AD

              ACH = - AR - AL - AF - AB
              ACV = - AU - AD

              AC = AC / fcscal ** 2

              ACH = ACH / fcscal ** 2
              ACV = ACV / fcscal ** 2

              AL = AL / fcscal ** 2
              AR = AR / fcscal ** 2
              AB = AB / fcscal ** 2
              AF = AF / fcscal ** 2
              AD = AD / (fcscal * fcscal_d)
              AU = AU / (fcscal * fcscal_u)

              ! ------------------ define matrix A -------------------

              if(poissonSolverType == 'bicgstab') then
                ac_b(i, j, k) = AC

                ach_b(i, j, k) = ACH
                acv_b(i, j, k) = ACV

                al_b(i, j, k) = AL
                ar_b(i, j, k) = AR

                ab_b(i, j, k) = AB
                af_b(i, j, k) = AF

                ad_b(i, j, k) = AD
                au_b(i, j, k) = AU

                alb_b(i, j, k) = 0.0
                alf_b(i, j, k) = 0.0
                arb_b(i, j, k) = 0.0
                arf_b(i, j, k) = 0.0
              else
                stop 'ERROR: val_PsIn expects bicgstab'
              end if
            end do ! i_loop
          end do ! j_loop
        end do ! k_loop
      else if(topography) then
        ! Compute tensor elements for TFC.
        do k = 1, nz
          ! Compute scaling factors.
          fcscal = sqrt(pStrat(k) ** 2.0 / rhoStrat(k))
          fcscal_u = sqrt(pStrat(k + 1) ** 2.0 / rhoStrat(k + 1))
          fcscal_d = sqrt(pStrat(k - 1) ** 2.0 / rhoStrat(k - 1))
          fcscal_uu = sqrt(pStrat(k + 2) ** 2.0 / rhoStrat(k + 2))
          fcscal_dd = sqrt(pStrat(k - 2) ** 2.0 / rhoStrat(k - 2))
          do j = 1, ny
            do i = 1, nx
              ! Compute inverse Jacobian.
              jacInv = 1.0 / jac(i, j, k)

              ! Compute P coefficients (divergence).
              pEdgeRDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i &
                  &+ 1, j, k) * pStratTFC(i + 1, j, k))
              pEdgeLDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i &
                  &- 1, j, k) * pStratTFC(i - 1, j, k))
              pEdgeFDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                  &+ 1, k) * pStratTFC(i, j + 1, k))
              pEdgeBDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                  &- 1, k) * pStratTFC(i, j - 1, k))
              pEdgeUDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, &
                  &k + 1) * pStratTFC(i, j, k + 1))
              pEdgeDDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, &
                  &k - 1) * pStratTFC(i, j, k - 1))

              ! Compute P coefficients (pressure gradient).
              pEdgeRGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                  &+ pStratTFC(i + 1, j, k) / jac(i + 1, j, k))
              pEdgeLGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                  &+ pStratTFC(i - 1, j, k) / jac(i - 1, j, k))
              pEdgeFGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                  &+ pStratTFC(i, j + 1, k) / jac(i, j + 1, k))
              pEdgeBGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                  &+ pStratTFC(i, j - 1, k) / jac(i, j - 1, k))
              pEdgeUGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                  &+ pStratTFC(i, j, k + 1) / jac(i, j, k + 1))
              pEdgeDGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                  &+ pStratTFC(i, j, k - 1) / jac(i, j, k - 1))

              ! Compute density coefficients.
              rhoEdgeR = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k) &
                  &+ rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, j, k))
              rhoEdgeL = 0.5 * (var%rho(i, j, k) + var%rho(i - 1, j, k) &
                  &+ rhoStratTFC(i, j, k) + rhoStratTFC(i - 1, j, k))
              rhoEdgeF = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k) &
                  &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j + 1, k))
              rhoEdgeB = 0.5 * (var%rho(i, j, k) + var%rho(i, j - 1, k) &
                  &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j - 1, k))
              rhoEdgeU = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1) &
                  &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j, k + 1))
              rhoEdgeD = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k - 1) &
                  &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j, k - 1))

              ! Compute (dP / dpi) for compressible model for AC
              if(model == "compressible") then
                dPdPi = 1 / (gamma - 1) * (Rsp / pref) ** (1 - gamma) &
                    &* var%P(i, j, k) ** (2 - gamma)
              end if

              ! --------------------- A(i,j,k) ---------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                AC = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * (1.0 &
                    &/ dx + 0.75 * met(i, j, k, 1, 3) / dz) + pEdgeLDiv &
                    &/ rhoEdgeL * pEdgeLGra * (1.0 / dx - 0.75 * met(i, j, k, &
                    &1, 3) / dz)) - jacInv / dy * (pEdgeFDiv / rhoEdgeF &
                    &* pEdgeFGra * (1.0 / dy + 0.75 * met(i, j, k, 2, 3) / dz) &
                    &+ pEdgeBDiv / rhoEdgeB * pEdgeBGra * (1.0 / dy - 0.75 &
                    &* met(i, j, k, 2, 3) / dz)) - jacInv / dz * pEdgeUDiv &
                    &/ rhoEdgeU * pEdgeUGra * met(i, j, k, 3, 3) / dz + jacInv &
                    &/ dz * pEdgeUDiv / rhoEdgeU * 0.5 * pStratTFC(i, j, k) &
                    &/ jac(i, j, k) * (chris(i, j, k, 1, 1) + chris(i, j, k, &
                    &2, 2) + 2.0 * chris(i, j, k, 1, 3) * met(i, j, k, 1, 3) &
                    &+ 2.0 * chris(i, j, k, 2, 3) * met(i, j, k, 2, 3))
              else if(k == nz .and. zBoundary == "solid_wall") then
                AC = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * (1.0 &
                    &/ dx - 0.75 * met(i, j, k, 1, 3) / dz) + pEdgeLDiv &
                    &/ rhoEdgeL * pEdgeLGra * (1.0 / dx + 0.75 * met(i, j, k, &
                    &1, 3) / dz)) - jacInv / dy * (pEdgeFDiv / rhoEdgeF &
                    &* pEdgeFGra * (1.0 / dy - 0.75 * met(i, j, k, 2, 3) / dz) &
                    &+ pEdgeBDiv / rhoEdgeB * pEdgeBGra * (1.0 / dy + 0.75 &
                    &* met(i, j, k, 2, 3) / dz)) - jacInv / dz * pEdgeDDiv &
                    &/ rhoEdgeD * pEdgeDGra * met(i, j, k, 3, 3) / dz - jacInv &
                    &/ dz * pEdgeDDiv / rhoEdgeD * 0.5 * pStratTFC(i, j, k) &
                    &/ jac(i, j, k) * (chris(i, j, k, 1, 1) + chris(i, j, k, &
                    &2, 2) + 2.0 * chris(i, j, k, 1, 3) * met(i, j, k, 1, 3) &
                    &+ 2.0 * chris(i, j, k, 2, 3) * met(i, j, k, 2, 3))
              else
                AC = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra / dx &
                    &+ pEdgeLDiv / rhoEdgeL * pEdgeLGra / dx) - jacInv / dy &
                    &* (pEdgeFDiv / rhoEdgeF * pEdgeFGra / dy + pEdgeBDiv &
                    &/ rhoEdgeB * pEdgeBGra / dy) - jacInv / dz * (pEdgeUDiv &
                    &/ rhoEdgeU * pEdgeUGra * met(i, j, k, 3, 3) / dz &
                    &+ pEdgeDDiv / rhoEdgeD * pEdgeDGra * met(i, j, k, 3, 3) &
                    &/ dz) + jacInv / dz * pEdgeUDiv / rhoEdgeU * 0.5 &
                    &* pStratTFC(i, j, k) / jac(i, j, k) * (chris(i, j, k, 1, &
                    &1) + chris(i, j, k, 2, 2) + 2.0 * chris(i, j, k, 1, 3) &
                    &* met(i, j, k, 1, 3) + 2.0 * chris(i, j, k, 2, 3) &
                    &* met(i, j, k, 2, 3)) - jacInv / dz * pEdgeDDiv &
                    &/ rhoEdgeD * 0.5 * pStratTFC(i, j, k) / jac(i, j, k) &
                    &* (chris(i, j, k, 1, 1) + chris(i, j, k, 2, 2) + 2.0 &
                    &* chris(i, j, k, 1, 3) * met(i, j, k, 1, 3) + 2.0 &
                    &* chris(i, j, k, 2, 3) * met(i, j, k, 2, 3))
              end if

              if(model == "compressible") then
                AC = AC - (dPdPi / dt) / (dt * Rsp / kappa)
              end if

              ! -------------------- A(i+1,j,k) --------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                AR = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * (1.0 &
                    &/ dx - 0.75 * met(i + 1, j, k, 1, 3) / dz) + jacInv / dz &
                    &* pEdgeUDiv / rhoEdgeU * pEdgeUGra * 0.25 * met(i + 1, j, &
                    &k, 1, 3) / dx
              else if(k == nz .and. zBoundary == "solid_wall") then
                AR = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * (1.0 &
                    &/ dx + 0.75 * met(i + 1, j, k, 1, 3) / dz) - jacInv / dz &
                    &* pEdgeDDiv / rhoEdgeD * pEdgeDGra * 0.25 * met(i + 1, j, &
                    &k, 1, 3) / dx
              else
                AR = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra / dx &
                    &+ jacInv / dz * (pEdgeUDiv / rhoEdgeU * pEdgeUGra * met(i &
                    &+ 1, j, k, 1, 3) * 0.25 / dx - pEdgeDDiv / rhoEdgeD &
                    &* pEdgeDGra * met(i + 1, j, k, 1, 3) * 0.25 / dx)
              end if

              ! -------------------- A(i-1,j,k) --------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                AL = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * (1.0 &
                    &/ dx + 0.75 * met(i - 1, j, k, 1, 3) / dz) - jacInv / dz &
                    &* pEdgeUDiv / rhoEdgeU * pEdgeUGra * 0.25 * met(i - 1, j, &
                    &k, 1, 3) / dx
              else if(k == nz .and. zBoundary == "solid_wall") then
                AL = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * (1.0 &
                    &/ dx - 0.75 * met(i - 1, j, k, 1, 3) / dz) + jacInv / dz &
                    &* pEdgeDDiv / rhoEdgeD * pEdgeDGra * 0.25 * met(i - 1, j, &
                    &k, 1, 3) / dx
              else
                AL = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra / dx &
                    &- jacInv / dz * (pEdgeUDiv / rhoEdgeU * pEdgeUGra * met(i &
                    &- 1, j, k, 1, 3) * 0.25 / dx - pEdgeDDiv / rhoEdgeD &
                    &* pEdgeDGra * met(i - 1, j, k, 1, 3) * 0.25 / dx)
              end if

              ! -------------------- A(i,j+1,k) --------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                AF = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * (1.0 &
                    &/ dy - 0.75 * met(i, j + 1, k, 2, 3) / dz) + jacInv / dz &
                    &* pEdgeUDiv / rhoEdgeU * pEdgeUGra * 0.25 * met(i, j + 1, &
                    &k, 2, 3) / dy
              else if(k == nz .and. zBoundary == "solid_wall") then
                AF = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * (1.0 &
                    &/ dy + 0.75 * met(i, j + 1, k, 2, 3) / dz) - jacInv / dz &
                    &* pEdgeDDiv / rhoEdgeD * pEdgeDGra * 0.25 * met(i, j + 1, &
                    &k, 2, 3) / dy
              else
                AF = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra / dy &
                    &+ jacInv / dz * (pEdgeUDiv / rhoEdgeU * pEdgeUGra &
                    &* met(i, j + 1, k, 2, 3) * 0.25 / dy - pEdgeDDiv &
                    &/ rhoEdgeD * pEdgeDGra * met(i, j + 1, k, 2, 3) * 0.25 &
                    &/ dy)
              end if

              ! -------------------- A(i,j-1,k) --------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                AB = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * (1.0 &
                    &/ dy + 0.75 * met(i, j - 1, k, 2, 3) / dz) - jacInv / dz &
                    &* pEdgeUDiv / rhoEdgeU * pEdgeUGra * 0.25 * met(i, j - 1, &
                    &k, 2, 3) / dy
              else if(k == nz .and. zBoundary == "solid_wall") then
                AB = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * (1.0 &
                    &/ dy - 0.75 * met(i, j - 1, k, 2, 3) / dz) + jacInv / dz &
                    &* pEdgeDDiv / rhoEdgeD * pEdgeDGra * 0.25 * met(i, j - 1, &
                    &k, 2, 3) / dy
              else
                AB = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra / dy &
                    &- jacInv / dz * (pEdgeUDiv / rhoEdgeU * pEdgeUGra &
                    &* met(i, j - 1, k, 2, 3) * 0.25 / dy - pEdgeDDiv &
                    &/ rhoEdgeD * pEdgeDGra * met(i, j - 1, k, 2, 3) * 0.25 &
                    &/ dy)
              end if

              ! -------------------- A(i,j,k+1) --------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                AU = jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * met(i, &
                    &j, k + 1, 1, 3) / dz - pEdgeLDiv / rhoEdgeL * pEdgeLGra &
                    &* met(i, j, k + 1, 1, 3) / dz) + jacInv / dy * (pEdgeFDiv &
                    &/ rhoEdgeF * pEdgeFGra * met(i, j, k + 1, 2, 3) / dz &
                    &- pEdgeBDiv / rhoEdgeB * pEdgeBGra * met(i, j, k + 1, 2, &
                    &3) / dz) + jacInv / dz * pEdgeUDiv / rhoEdgeU * pEdgeUGra &
                    &* met(i, j, k + 1, 3, 3) / dz + jacInv / dz * pEdgeUDiv &
                    &/ rhoEdgeU * 0.5 * pStratTFC(i, j, k + 1) / jac(i, j, k &
                    &+ 1) * (chris(i, j, k + 1, 1, 1) + chris(i, j, k + 1, 2, &
                    &2) + 2.0 * chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, 1, &
                    &3) + 2.0 * chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, 2, &
                    &3))
              else if(k == nz .and. zBoundary == "solid_wall") then
                AU = 0.0
              else
                AU = jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * met(i, &
                    &j, k + 1, 1, 3) * 0.25 / dz - pEdgeLDiv / rhoEdgeL &
                    &* pEdgeLGra * met(i, j, k + 1, 1, 3) * 0.25 / dz) &
                    &+ jacInv / dy * (pEdgeFDiv / rhoEdgeF * pEdgeFGra &
                    &* met(i, j, k + 1, 2, 3) * 0.25 / dz - pEdgeBDiv &
                    &/ rhoEdgeB * pEdgeBGra * met(i, j, k + 1, 2, 3) * 0.25 &
                    &/ dz) + jacInv / dz * pEdgeUDiv / rhoEdgeU * pEdgeUGra &
                    &* met(i, j, k + 1, 3, 3) / dz + jacInv / dz * pEdgeUDiv &
                    &/ rhoEdgeU * 0.5 * pStratTFC(i, j, k + 1) / jac(i, j, k &
                    &+ 1) * (chris(i, j, k + 1, 1, 1) + chris(i, j, k + 1, 2, &
                    &2) + 2.0 * chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, 1, &
                    &3) + 2.0 * chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, 2, &
                    &3))
              end if

              ! -------------------- A(i,j,k-1) --------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                AD = 0.0
              else if(k == nz .and. zBoundary == "solid_wall") then
                AD = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra &
                    &* met(i, j, k - 1, 1, 3) / dz - pEdgeLDiv / rhoEdgeL &
                    &* pEdgeLGra * met(i, j, k - 1, 1, 3) / dz) - jacInv / dy &
                    &* (pEdgeFDiv / rhoEdgeF * pEdgeFGra * met(i, j, k - 1, 2, &
                    &3) / dz - pEdgeBDiv / rhoEdgeB * pEdgeBGra * met(i, j, k &
                    &- 1, 2, 3) / dz) + jacInv / dz * pEdgeDDiv / rhoEdgeD &
                    &* pEdgeDGra * met(i, j, k - 1, 3, 3) / dz - jacInv / dz &
                    &* pEdgeDDiv / rhoEdgeD * 0.5 * pStratTFC(i, j, k - 1) &
                    &/ jac(i, j, k - 1) * (chris(i, j, k - 1, 1, 1) + chris(i, &
                    &j, k - 1, 2, 2) + 2.0 * chris(i, j, k - 1, 1, 3) * met(i, &
                    &j, k - 1, 1, 3) + 2.0 * chris(i, j, k - 1, 2, 3) * met(i, &
                    &j, k - 1, 2, 3))
              else
                AD = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra &
                    &* met(i, j, k - 1, 1, 3) * 0.25 / dz - pEdgeLDiv &
                    &/ rhoEdgeL * pEdgeLGra * met(i, j, k - 1, 1, 3) * 0.25 &
                    &/ dz) - jacInv / dy * (pEdgeFDiv / rhoEdgeF * pEdgeFGra &
                    &* met(i, j, k - 1, 2, 3) * 0.25 / dz - pEdgeBDiv &
                    &/ rhoEdgeB * pEdgeBGra * met(i, j, k - 1, 2, 3) * 0.25 &
                    &/ dz) + jacInv / dz * pEdgeDDiv / rhoEdgeD * pEdgeDGra &
                    &* met(i, j, k - 1, 3, 3) / dz - jacInv / dz * pEdgeDDiv &
                    &/ rhoEdgeD * 0.5 * pStratTFC(i, j, k - 1) / jac(i, j, k &
                    &- 1) * (chris(i, j, k - 1, 1, 1) + chris(i, j, k - 1, 2, &
                    &2) + 2.0 * chris(i, j, k - 1, 1, 3) * met(i, j, k - 1, 1, &
                    &3) + 2.0 * chris(i, j, k - 1, 2, 3) * met(i, j, k - 1, 2, &
                    &3))
              end if

              ! ------------------- A(i+1,j,k+1) -------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                ARU = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * met(i &
                    &+ 1, j, k + 1, 1, 3) / dz + jacInv / dz * pEdgeUDiv &
                    &/ rhoEdgeU * pEdgeUGra * met(i + 1, j, k + 1, 1, 3) &
                    &* 0.25 / dx
              else if(k == nz .and. zBoundary == "solid_wall") then
                ARU = 0.0
              else
                ARU = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * met(i &
                    &+ 1, j, k + 1, 1, 3) * 0.25 / dz + jacInv / dz &
                    &* pEdgeUDiv / rhoEdgeU * pEdgeUGra * met(i + 1, j, k + 1, &
                    &1, 3) * 0.25 / dx
              end if

              ! ------------------- A(i+1,j,k-1) -------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                ARD = 0.0
              else if(k == nz .and. zBoundary == "solid_wall") then
                ARD = - jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * met(i &
                    &+ 1, j, k - 1, 1, 3) / dz - jacInv / dz * pEdgeDDiv &
                    &/ rhoEdgeD * pEdgeDGra * met(i + 1, j, k - 1, 1, 3) &
                    &* 0.25 / dx
              else
                ARD = - jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * met(i &
                    &+ 1, j, k - 1, 1, 3) * 0.25 / dz - jacInv / dz &
                    &* pEdgeDDiv / rhoEdgeD * pEdgeDGra * met(i + 1, j, k - 1, &
                    &1, 3) * 0.25 / dx
              end if

              ! ------------------- A(i-1,j,k+1) -------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                ALU = - jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * met(i &
                    &- 1, j, k + 1, 1, 3) / dz - jacInv / dz * pEdgeUDiv &
                    &/ rhoEdgeU * pEdgeUGra * met(i - 1, j, k + 1, 1, 3) &
                    &* 0.25 / dx
              else if(k == nz .and. zBoundary == "solid_wall") then
                ALU = 0.0
              else
                ALU = - jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * met(i &
                    &- 1, j, k + 1, 1, 3) * 0.25 / dz - jacInv / dz &
                    &* pEdgeUDiv / rhoEdgeU * pEdgeUGra * met(i - 1, j, k + 1, &
                    &1, 3) * 0.25 / dx
              end if

              ! ------------------- A(i-1,j,k-1) -------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                ALD = 0.0
              else if(k == nz .and. zBoundary == "solid_wall") then
                ALD = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * met(i &
                    &- 1, j, k - 1, 1, 3) / dz + jacInv / dz * pEdgeDDiv &
                    &/ rhoEdgeD * pEdgeDGra * met(i - 1, j, k - 1, 1, 3) &
                    &* 0.25 / dx
              else
                ALD = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * met(i &
                    &- 1, j, k - 1, 1, 3) * 0.25 / dz + jacInv / dz &
                    &* pEdgeDDiv / rhoEdgeD * pEdgeDGra * met(i - 1, j, k - 1, &
                    &1, 3) * 0.25 / dx
              end if

              ! ------------------- A(i,j+1,k+1) -------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                AFU = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * met(i, &
                    &j + 1, k + 1, 2, 3) / dz + jacInv / dz * pEdgeUDiv &
                    &/ rhoEdgeU * pEdgeUGra * met(i, j + 1, k + 1, 2, 3) &
                    &* 0.25 / dy
              else if(k == nz .and. zBoundary == "solid_wall") then
                AFU = 0.0
              else
                AFU = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * met(i, &
                    &j + 1, k + 1, 2, 3) * 0.25 / dz + jacInv / dz * pEdgeUDiv &
                    &/ rhoEdgeU * pEdgeUGra * met(i, j + 1, k + 1, 2, 3) &
                    &* 0.25 / dy
              end if

              ! ------------------- A(i,j+1,k-1) -------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                AFD = 0.0
              else if(k == nz .and. zBoundary == "solid_wall") then
                AFD = - jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra &
                    &* met(i, j + 1, k - 1, 2, 3) / dz - jacInv / dz &
                    &* pEdgeDDiv / rhoEdgeD * pEdgeDGra * met(i, j + 1, k - 1, &
                    &2, 3) * 0.25 / dy
              else
                AFD = - jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra &
                    &* met(i, j + 1, k - 1, 2, 3) * 0.25 / dz - jacInv / dz &
                    &* pEdgeDDiv / rhoEdgeD * pEdgeDGra * met(i, j + 1, k - 1, &
                    &2, 3) * 0.25 / dy
              end if

              ! ------------------- A(i,j-1,k+1) -------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                ABU = - jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra &
                    &* met(i, j - 1, k + 1, 2, 3) / dz - jacInv / dz &
                    &* pEdgeUDiv / rhoEdgeU * pEdgeUGra * met(i, j - 1, k + 1, &
                    &2, 3) * 0.25 / dy
              else if(k == nz .and. zBoundary == "solid_wall") then
                ABU = 0.0
              else
                ABU = - jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra &
                    &* met(i, j - 1, k + 1, 2, 3) * 0.25 / dz - jacInv / dz &
                    &* pEdgeUDiv / rhoEdgeU * pEdgeUGra * met(i, j - 1, k + 1, &
                    &2, 3) * 0.25 / dy
              end if

              ! ------------------- A(i,j-1,k-1) -------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                ABD = 0.0
              else if(k == nz .and. zBoundary == "solid_wall") then
                ABD = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * met(i, &
                    &j - 1, k - 1, 2, 3) / dz + jacInv / dz * pEdgeDDiv &
                    &/ rhoEdgeD * pEdgeDGra * met(i, j - 1, k - 1, 2, 3) &
                    &* 0.25 / dy
              else
                ABD = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * met(i, &
                    &j - 1, k - 1, 2, 3) * 0.25 / dz + jacInv / dz * pEdgeDDiv &
                    &/ rhoEdgeD * pEdgeDGra * met(i, j - 1, k - 1, 2, 3) &
                    &* 0.25 / dy
              end if

              ! ------------------- A(i,j,k+2) ---------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                AUU = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * 0.25 &
                    &* met(i, j, k + 2, 1, 3) / dz - pEdgeLDiv / rhoEdgeL &
                    &* pEdgeLGra * 0.25 * met(i, j, k + 2, 1, 3) / dz) &
                    &- jacInv / dy * (pEdgeFDiv / rhoEdgeF * pEdgeFGra * 0.25 &
                    &* met(i, j, k + 2, 2, 3) / dz - pEdgeBDiv / rhoEdgeB &
                    &* pEdgeBGra * 0.25 * met(i, j, k + 2, 2, 3) / dz)
              else
                AUU = 0.0
              end if

              ! ------------------- A(i,j,k-2) ---------------------

              if(k == nz .and. zBoundary == "solid_wall") then
                ADD = jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * 0.25 &
                    &* met(i, j, k - 2, 1, 3) / dz - pEdgeLDiv / rhoEdgeL &
                    &* pEdgeLGra * 0.25 * met(i, j, k - 2, 1, 3) / dz) &
                    &+ jacInv / dy * (pEdgeFDiv / rhoEdgeF * pEdgeFGra * 0.25 &
                    &* met(i, j, k - 2, 2, 3) / dz - pEdgeBDiv / rhoEdgeB &
                    &* pEdgeBGra * 0.25 * met(i, j, k - 2, 2, 3) / dz)
              else
                ADD = 0.0
              end if

              ! ------------------ A(i+1,j,k+2) --------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                ARUU = - jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * 0.25 &
                    &* met(i + 1, j, k + 2, 1, 3) / dz
              else
                ARUU = 0.0
              end if

              ! ------------------ A(i+1,j,k-2) --------------------

              if(k == nz .and. zBoundary == "solid_wall") then
                ARDD = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * 0.25 &
                    &* met(i + 1, j, k - 2, 1, 3) / dz
              else
                ARDD = 0.0
              end if

              ! ------------------ A(i-1,j,k+2) --------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                ALUU = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * 0.25 &
                    &* met(i - 1, j, k + 2, 1, 3) / dz
              else
                ALUU = 0.0
              end if

              ! ------------------ A(i-1,j,k-2) --------------------

              if(k == nz .and. zBoundary == "solid_wall") then
                ALDD = - jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * 0.25 &
                    &* met(i - 1, j, k - 2, 1, 3) / dz
              else
                ALDD = 0.0
              end if

              ! ------------------ A(i,j+1,k+2) --------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                AFUU = - jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * 0.25 &
                    &* met(i, j + 1, k + 2, 2, 3) / dz
              else
                AFUU = 0.0
              end if

              ! ------------------ A(i,j+1,k-2) --------------------

              if(k == nz .and. zBoundary == "solid_wall") then
                AFDD = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * 0.25 &
                    &* met(i, j + 1, k - 2, 2, 3) / dz
              else
                AFDD = 0.0
              end if

              ! ------------------ A(i,j-1,k+2) --------------------

              if(k == 1 .and. zBoundary == "solid_wall") then
                ABUU = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * 0.25 &
                    &* met(i, j - 1, k + 2, 2, 3) / dz
              else
                ABUU = 0.0
              end if

              ! ------------------ A(i,j-1,k-2) --------------------

              if(k == nz .and. zBoundary == "solid_wall") then
                ABDD = - jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * 0.25 &
                    &* met(i, j - 1, k - 2, 2, 3) / dz
              else
                ABDD = 0.0
              end if

              ! Scale the tensor elements.
              AC = AC / (fcscal ** 2.0)
              AR = AR / (fcscal ** 2.0)
              AL = AL / (fcscal ** 2.0)
              AF = AF / (fcscal ** 2.0)
              AB = AB / (fcscal ** 2.0)
              AU = AU / fcscal / fcscal_u
              AD = AD / fcscal / fcscal_d
              ARU = ARU / fcscal / fcscal_u
              ARD = ARD / fcscal / fcscal_d
              ALU = ALU / fcscal / fcscal_u
              ALD = ALD / fcscal / fcscal_d
              AFU = AFU / fcscal / fcscal_u
              AFD = AFD / fcscal / fcscal_d
              ABU = ABU / fcscal / fcscal_u
              ABD = ABD / fcscal / fcscal_d
              AUU = AUU / fcscal / fcscal_uu
              ADD = ADD / fcscal / fcscal_dd
              ARUU = ARUU / fcscal / fcscal_uu
              ARDD = ARDD / fcscal / fcscal_dd
              ALUU = ALUU / fcscal / fcscal_uu
              ALDD = ALDD / fcscal / fcscal_dd
              AFUU = AFUU / fcscal / fcscal_uu
              AFDD = AFDD / fcscal / fcscal_dd
              ABUU = ABUU / fcscal / fcscal_uu
              ABDD = ABDD / fcscal / fcscal_dd

              ! Set tensor elements for bicgstab.
              ac_b(i, j, k) = AC
              ar_b(i, j, k) = AR
              al_b(i, j, k) = AL
              af_b(i, j, k) = AF
              ab_b(i, j, k) = AB
              au_b(i, j, k) = AU
              ad_b(i, j, k) = AD
              aru_b(i, j, k) = ARU
              ard_b(i, j, k) = ARD
              alu_b(i, j, k) = ALU
              ald_b(i, j, k) = ALD
              afu_b(i, j, k) = AFU
              afd_b(i, j, k) = AFD
              abu_b(i, j, k) = ABU
              abd_b(i, j, k) = ABD
              auu_b(i, j, k) = AUU
              add_b(i, j, k) = ADD
              aruu_b(i, j, k) = ARUU
              ardd_b(i, j, k) = ARDD
              aluu_b(i, j, k) = ALUU
              aldd_b(i, j, k) = ALDD
              afuu_b(i, j, k) = AFUU
              afdd_b(i, j, k) = AFDD
              abuu_b(i, j, k) = ABUU
              abdd_b(i, j, k) = ABDD

              ! Store horizontal and vertical components of AC (for
              ! preconditioner).
              if(preconditioner == "yes") then
                ach_b(i, j, k) = - AR - AL - AF - AB
                acv_b(i, j, k) = - AU - AD
              end if
            end do
          end do
        end do
      else
        do k = 1, nz
          fcscal = sqrt(Pstrat(k) ** 2 / rhoStrat(k))
          fcscal_u = sqrt(Pstrat(k + 1) ** 2 / rhoStrat(k + 1))
          fcscal_d = sqrt(Pstrat(k - 1) ** 2 / rhoStrat(k - 1))
          do j = 1, ny
            do i = 1, nx

              ! ------------------ A(i+1,j,k) ------------------

              rhoEdge = 0.5 * (var%rho(i + 1, j, k) + var%rho(i, j, k)) &
                  &+ rhoStrat(k)

              AR = dx2 * pStrat(k) ** 2 / rhoEdge
              if(pressureScaling) then
                AR = AR / Pstrat(k)
              end if

              ! ------------------- A(i-1,j,k) --------------------

              rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i - 1, j, k)) &
                  &+ rhoStrat(k)

              AL = dx2 * pStrat(k) ** 2 / rhoEdge
              if(pressureScaling) then
                AL = AL / Pstrat(k)
              end if

              ! -------------------- A(i,j+1,k) ----------------------

              rhoEdge = 0.5 * (var%rho(i, j + 1, k) + var%rho(i, j, k)) &
                  &+ rhoStrat(k)

              AF = dy2 * pStrat(k) ** 2 / rhoEdge
              if(pressureScaling) then
                AF = AF / Pstrat(k)
              end if

              ! --------------------- A(i,j-1,k) -------------------

              rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i, j - 1, k)) &
                  &+ rhoStrat(k)

              AB = dy2 * pStrat(k) ** 2 / rhoEdge
              if(pressureScaling) then
                AB = AB / Pstrat(k)
              end if

              ! ---------------------- A(i,j,k+1) ------------------

              ! TFC FJ
              if(k == nz .and. zBoundary == "solid_wall") then
                AU = 0.0
              else
                rhoEdge = 0.5 * (var%rho(i, j, k + 1) + var%rho(i, j, k)) &
                    &+ rhoStratTilde(k)

                pStratU = 0.5 * (pStrat(k + 1) + pStrat(k))
                AU = dz2 * pStratU ** 2 / rhoEdge
              end if
              if(pressureScaling) then
                AU = AU / PstratTilde(k)
              end if

              ! ----------------------- A(i,j,k-1) -----------------

              ! TFC
              if(k == 1 .and. zBoundary == "solid_wall") then
                AD = 0.0
              else
                rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k - 1)) &
                    &+ rhoStratTilde(k - 1)

                pStratD = 0.5 * (pStrat(k) + pStrat(k - 1))
                AD = dz2 * pStratD ** 2 / rhoEdge
              end if
              if(pressureScaling) then
                AD = AD / PstratTilde(k - 1)
              end if

              if(pressureScaling) then
                ! stop'ERROR: no pressure scaling allowed'
              end if

              ! ----------------------- A(i,j,k) -------------------

              AC = - AR - AL - AF - AB - AU - AD

              ACH = - AR - AL - AF - AB
              ACV = - AU - AD

              ! ------------------ define matrix A -------------------

              AC = AC / fcscal ** 2

              ACV = ACV / fcscal ** 2
              ACH = ACH / fcscal ** 2

              AL = AL / fcscal ** 2
              AR = AR / fcscal ** 2
              AB = AB / fcscal ** 2
              AF = AF / fcscal ** 2
              AD = AD / (fcscal * fcscal_d)
              AU = AU / (fcscal * fcscal_u)

              if(poissonSolverType == 'bicgstab') then
                ac_b(i, j, k) = AC

                ach_b(i, j, k) = ACH
                acv_b(i, j, k) = ACV

                al_b(i, j, k) = AL
                ar_b(i, j, k) = AR

                ab_b(i, j, k) = AB
                af_b(i, j, k) = AF

                ad_b(i, j, k) = AD
                au_b(i, j, k) = AU
              else
                stop 'ERROR: val_PsIn expects bicgstab'
              end if

            end do ! i_loop
          end do ! j_loop
        end do ! k_loop
      end if
    else if(opt == "impl") then
      if(timeScheme /= "semiimplicit") then
        stop 'ERROR: for opt = impl must have timeScheme = semiimplicit'
      end if

      kr_sp = kr_sp * facray
      alprlx = alprlx * facray

      if(topography) then
        kr_sp_tfc = kr_sp_tfc * facray
        kr_sp_w_tfc = kr_sp_w_tfc * facray
        do k = 1, nz
          ! Compute scaling factors
          fcscal = sqrt(pStrat(k) ** 2.0 / rhoStrat(k))
          fcscal_u = sqrt(pStrat(k + 1) ** 2.0 / rhoStrat(k + 1))
          fcscal_d = sqrt(pStrat(k - 1) ** 2.0 / rhoStrat(k - 1))
          fcscal_uu = sqrt(pStrat(k + 2) ** 2.0 / rhoStrat(k + 2))
          fcscal_dd = sqrt(pStrat(k - 2) ** 2.0 / rhoStrat(k - 2))
          do j = 1, ny
            do i = 1, nx
              ! Compute tensor elements for TFC.

              ! Compute inverse Jacobian.
              jacInv = 1.0 / jac(i, j, k)

              ! Compute P coefficients (divergence).
              pEdgeRDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i &
                  &+ 1, j, k) * pStratTFC(i + 1, j, k))
              pEdgeLDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i &
                  &- 1, j, k) * pStratTFC(i - 1, j, k))
              pEdgeFDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                  &+ 1, k) * pStratTFC(i, j + 1, k))
              pEdgeBDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j &
                  &- 1, k) * pStratTFC(i, j - 1, k))
              pEdgeUDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, &
                  &k + 1) * pStratTFC(i, j, k + 1))
              pEdgeDDiv = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j, &
                  &k - 1) * pStratTFC(i, j, k - 1))

              ! Compute P coefficients (pressure gradient).
              pEdgeRGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                  &+ pStratTFC(i + 1, j, k) / jac(i + 1, j, k))
              pEdgeLGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                  &+ pStratTFC(i - 1, j, k) / jac(i - 1, j, k))
              pEdgeFGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                  &+ pStratTFC(i, j + 1, k) / jac(i, j + 1, k))
              pEdgeBGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                  &+ pStratTFC(i, j - 1, k) / jac(i, j - 1, k))
              pEdgeUGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                  &+ pStratTFC(i, j, k + 1) / jac(i, j, k + 1))
              pEdgeDGra = 0.5 * (pStratTFC(i, j, k) / jac(i, j, k) &
                  &+ pStratTFC(i, j, k - 1) / jac(i, j, k - 1))
              pUEdgeRGra = 0.5 * (pStratTFC(i, j, k + 1) / jac(i, j, k + 1) &
                  &+ pStratTFC(i + 1, j, k + 1) / jac(i + 1, j, k + 1))
              pUEdgeLGra = 0.5 * (pStratTFC(i, j, k + 1) / jac(i, j, k + 1) &
                  &+ pStratTFC(i - 1, j, k + 1) / jac(i - 1, j, k + 1))
              pUEdgeFGra = 0.5 * (pStratTFC(i, j, k + 1) / jac(i, j, k + 1) &
                  &+ pStratTFC(i, j + 1, k + 1) / jac(i, j + 1, k + 1))
              pUEdgeBGra = 0.5 * (pStratTFC(i, j, k + 1) / jac(i, j, k + 1) &
                  &+ pStratTFC(i, j - 1, k + 1) / jac(i, j - 1, k + 1))
              pDEdgeRGra = 0.5 * (pStratTFC(i, j, k - 1) / jac(i, j, k - 1) &
                  &+ pStratTFC(i + 1, j, k - 1) / jac(i + 1, j, k - 1))
              pDEdgeLGra = 0.5 * (pStratTFC(i, j, k - 1) / jac(i, j, k - 1) &
                  &+ pStratTFC(i - 1, j, k - 1) / jac(i - 1, j, k - 1))
              pDEdgeFGra = 0.5 * (pStratTFC(i, j, k - 1) / jac(i, j, k - 1) &
                  &+ pStratTFC(i, j + 1, k - 1) / jac(i, j + 1, k - 1))
              pDEdgeBGra = 0.5 * (pStratTFC(i, j, k - 1) / jac(i, j, k - 1) &
                  &+ pStratTFC(i, j - 1, k - 1) / jac(i, j - 1, k - 1))

              ! Compute density coefficients.
              rhoStratEdgeR = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, &
                  &j, k))
              rhoStratEdgeL = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i - 1, &
                  &j, k))
              rhoStratEdgeF = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, j &
                  &+ 1, k))
              rhoStratEdgeB = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, j &
                  &- 1, k))
              rhoStratEdgeU = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, j, &
                  &k + 1))
              rhoStratEdgeD = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, j, &
                  &k - 1))
              rhoEdgeR = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k)) &
                  &+ rhoStratEdgeR
              rhoEdgeL = 0.5 * (var%rho(i, j, k) + var%rho(i - 1, j, k)) &
                  &+ rhoStratEdgeL
              rhoEdgeF = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k)) &
                  &+ rhoStratEdgeF
              rhoEdgeB = 0.5 * (var%rho(i, j, k) + var%rho(i, j - 1, k)) &
                  &+ rhoStratEdgeB
              rhoEdgeU = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k + 1)) &
                  &+ rhoStratEdgeU
              rhoEdgeD = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k - 1)) &
                  &+ rhoStratEdgeD

              ! In the compressible model, N2 depends on the full density,
              ! so that the factor rhoStrat / rho is replaced by 1.
              if(model == "compressible") then
                rhoStratEdgeR = rhoEdgeR
                rhoStratEdgeL = rhoEdgeL
                rhoStratEdgeF = rhoEdgeF
                rhoStratEdgeB = rhoEdgeB
                rhoStratEdgeU = rhoEdgeU
                rhoStratEdgeD = rhoEdgeD
              end if

              rhoUEdgeR = 0.5 * (var%rho(i, j, k + 1) + var%rho(i + 1, j, k &
                  &+ 1) + rhoStratTFC(i, j, k + 1) + rhoStratTFC(i + 1, j, k &
                  &+ 1))
              rhoUEdgeL = 0.5 * (var%rho(i, j, k + 1) + var%rho(i - 1, j, k &
                  &+ 1) + rhoStratTFC(i, j, k + 1) + rhoStratTFC(i - 1, j, k &
                  &+ 1))
              rhoUEdgeF = 0.5 * (var%rho(i, j, k + 1) + var%rho(i, j + 1, k &
                  &+ 1) + rhoStratTFC(i, j, k + 1) + rhoStratTFC(i, j + 1, k &
                  &+ 1))
              rhoUEdgeB = 0.5 * (var%rho(i, j, k + 1) + var%rho(i, j - 1, k &
                  &+ 1) + rhoStratTFC(i, j, k + 1) + rhoStratTFC(i, j - 1, k &
                  &+ 1))
              rhoDEdgeR = 0.5 * (var%rho(i, j, k - 1) + var%rho(i + 1, j, k &
                  &- 1) + rhoStratTFC(i, j, k - 1) + rhoStratTFC(i + 1, j, k &
                  &- 1))
              rhoDEdgeL = 0.5 * (var%rho(i, j, k - 1) + var%rho(i - 1, j, k &
                  &- 1) + rhoStratTFC(i, j, k - 1) + rhoStratTFC(i - 1, j, k &
                  &- 1))
              rhoDEdgeF = 0.5 * (var%rho(i, j, k - 1) + var%rho(i, j + 1, k &
                  &- 1) + rhoStratTFC(i, j, k - 1) + rhoStratTFC(i, j + 1, k &
                  &- 1))
              rhoDEdgeB = 0.5 * (var%rho(i, j, k - 1) + var%rho(i, j - 1, k &
                  &- 1) + rhoStratTFC(i, j, k - 1) + rhoStratTFC(i, j - 1, k &
                  &- 1))

              ! Compute squared buoyancy frequency at edges.
              bvsStratEdgeU = 0.5 * (bvsStratTFC(i, j, k) + bvsStratTFC(i, j, &
                  &k + 1))
              bvsStratEdgeD = 0.5 * (bvsStratTFC(i, j, k) + bvsStratTFC(i, j, &
                  &k - 1))

              ! Compute (dP / dpi) for compressible model for AC
              if(model == "compressible") then
                dPdPi = 1 / (gamma - 1) * (Rsp / pref) ** (1 - gamma) &
                    &* var%P(i, j, k) ** (2 - gamma)
              end if

              ! Compute Rayleigh damping terms.
              facEdgeR = 1.0
              facEdgeL = 1.0
              facEdgeF = 1.0
              facEdgeB = 1.0
              facUEdgeR = 1.0
              facUEdgeL = 1.0
              facUEdgeF = 1.0
              facUEdgeB = 1.0
              facDEdgeR = 1.0
              facDEdgeL = 1.0
              facDEdgeF = 1.0
              facDEdgeB = 1.0
              facEdgeU = 1.0
              facEdgeD = 1.0
              if(spongeLayer) then
                if(sponge_uv) then
                  facEdgeR = facEdgeR + dt * 0.5 * (kr_sp_tfc(i, j, k) &
                      &+ kr_sp_tfc(i + 1, j, k))
                  facEdgeL = facEdgeL + dt * 0.5 * (kr_sp_tfc(i, j, k) &
                      &+ kr_sp_tfc(i - 1, j, k))
                  facEdgeF = facEdgeF + dt * 0.5 * (kr_sp_tfc(i, j, k) &
                      &+ kr_sp_tfc(i, j + 1, k))
                  facEdgeB = facEdgeB + dt * 0.5 * (kr_sp_tfc(i, j, k) &
                      &+ kr_sp_tfc(i, j - 1, k))
                  facUEdgeR = facUEdgeR + dt * 0.5 * (kr_sp_tfc(i, j, k + 1) &
                      &+ kr_sp_tfc(i + 1, j, k + 1))
                  facUEdgeL = facUEdgeL + dt * 0.5 * (kr_sp_tfc(i, j, k + 1) &
                      &+ kr_sp_tfc(i - 1, j, k + 1))
                  facUEdgeF = facUEdgeF + dt * 0.5 * (kr_sp_tfc(i, j, k + 1) &
                      &+ kr_sp_tfc(i, j + 1, k + 1))
                  facUEdgeB = facUEdgeB + dt * 0.5 * (kr_sp_tfc(i, j, k + 1) &
                      &+ kr_sp_tfc(i, j - 1, k + 1))
                  facDEdgeR = facDEdgeR + dt * 0.5 * (kr_sp_tfc(i, j, k - 1) &
                      &+ kr_sp_tfc(i + 1, j, k - 1))
                  facDEdgeL = facDEdgeL + dt * 0.5 * (kr_sp_tfc(i, j, k - 1) &
                      &+ kr_sp_tfc(i - 1, j, k - 1))
                  facDEdgeF = facDEdgeF + dt * 0.5 * (kr_sp_tfc(i, j, k - 1) &
                      &+ kr_sp_tfc(i, j + 1, k - 1))
                  facDEdgeB = facDEdgeB + dt * 0.5 * (kr_sp_tfc(i, j, k - 1) &
                      &+ kr_sp_tfc(i, j - 1, k - 1))
                end if
                facEdgeU = facEdgeU + dt * 0.5 * (kr_sp_w_tfc(i, j, k) &
                    &+ kr_sp_w_tfc(i, j, k + 1))
                facEdgeD = facEdgeD + dt * 0.5 * (kr_sp_w_tfc(i, j, k) &
                    &+ kr_sp_w_tfc(i, j, k - 1))
              end if

              ! Compute implicit coefficients.
              impHorEdgeR = 1.0 / (facEdgeR ** 2.0)
              impHorEdgeL = 1.0 / (facEdgeL ** 2.0)
              impHorEdgeF = 1.0 / (facEdgeF ** 2.0)
              impHorEdgeB = 1.0 / (facEdgeB ** 2.0)
              impHorUEdgeR = 1.0 / (facUEdgeR ** 2.0)
              impHorUEdgeL = 1.0 / (facUEdgeL ** 2.0)
              impHorUEdgeF = 1.0 / (facUEdgeF ** 2.0)
              impHorUEdgeB = 1.0 / (facUEdgeB ** 2.0)
              impHorDEdgeR = 1.0 / (facDEdgeR ** 2.0)
              impHorDEdgeL = 1.0 / (facDEdgeL ** 2.0)
              impHorDEdgeF = 1.0 / (facDEdgeF ** 2.0)
              impHorDEdgeB = 1.0 / (facDEdgeB ** 2.0)
              impVerEdgeU = 1.0 / (facEdgeU + rhoStratEdgeU / rhoEdgeU &
                  &* bvsStratEdgeU * dt ** 2.0)
              impVerEdgeD = 1.0 / (facEdgeD + rhoStratEdgeD / rhoEdgeD &
                  &* bvsStratEdgeD * dt ** 2.0)

              ! Compute gradient coefficients

              ! G(i + 1 / 2)
              if(k == 1 .and. zBoundary == "solid_wall") then
                gEdgeR = jacInv / dx * pEdgeRDiv * impHorEdgeR * facEdgeR &
                    &/ rhoEdgeR + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                    &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 1, 3) * impHorEdgeR * facEdgeR &
                    &/ rhoEdgeR
              else if(k == nz .and. zBoundary == "solid_wall") then
                gEdgeR = jacInv / dx * pEdgeRDiv * impHorEdgeR * facEdgeR &
                    &/ rhoEdgeR - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                    &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 1, 3) * impHorEdgeR * facEdgeR &
                    &/ rhoEdgeR
              else
                gEdgeR = jacInv / dx * pEdgeRDiv * impHorEdgeR * facEdgeR &
                    &/ rhoEdgeR + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                    &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 1, 3) * impHorEdgeR * facEdgeR &
                    &/ rhoEdgeR - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                    &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 1, 3) * impHorEdgeR * facEdgeR &
                    &/ rhoEdgeR
              end if

              ! G(i - 1 / 2)
              if(k == 1 .and. zBoundary == "solid_wall") then
                gEdgeL = - jacInv / dx * pEdgeLDiv * impHorEdgeL * facEdgeL &
                    &/ rhoEdgeL + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                    &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 1, 3) * impHorEdgeL * facEdgeL &
                    &/ rhoEdgeL
              else if(k == nz .and. zBoundary == "solid_wall") then
                gEdgeL = - jacInv / dx * pEdgeLDiv * impHorEdgeL * facEdgeL &
                    &/ rhoEdgeL - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                    &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 1, 3) * impHorEdgeL * facEdgeL &
                    &/ rhoEdgeL
              else
                gEdgeL = - jacInv / dx * pEdgeLDiv * impHorEdgeL * facEdgeL &
                    &/ rhoEdgeL + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                    &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 1, 3) * impHorEdgeL * facEdgeL &
                    &/ rhoEdgeL - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                    &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 1, 3) * impHorEdgeL * facEdgeL &
                    &/ rhoEdgeL
              end if

              ! G(j + 1 / 2)
              if(k == 1 .and. zBoundary == "solid_wall") then
                gEdgeF = jacInv / dy * pEdgeFDiv * impHorEdgeF * facEdgeF &
                    &/ rhoEdgeF + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                    &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 2, 3) * impHorEdgeF * facEdgeF &
                    &/ rhoEdgeF
              else if(k == nz .and. zBoundary == "solid_wall") then
                gEdgeF = jacInv / dy * pEdgeFDiv * impHorEdgeF * facEdgeF &
                    &/ rhoEdgeF - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                    &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 2, 3) * impHorEdgeF * facEdgeF &
                    &/ rhoEdgeF
              else
                gEdgeF = jacInv / dy * pEdgeFDiv * impHorEdgeF * facEdgeF &
                    &/ rhoEdgeF + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                    &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 2, 3) * impHorEdgeF * facEdgeF &
                    &/ rhoEdgeF - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                    &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 2, 3) * impHorEdgeF * facEdgeF &
                    &/ rhoEdgeF
              end if

              ! G(j - 1 / 2)
              if(k == 1 .and. zBoundary == "solid_wall") then
                gEdgeB = - jacInv / dy * pEdgeBDiv * impHorEdgeB * facEdgeB &
                    &/ rhoEdgeB + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                    &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 2, 3) * impHorEdgeB * facEdgeB &
                    &/ rhoEdgeB
              else if(k == nz .and. zBoundary == "solid_wall") then
                gEdgeB = - jacInv / dy * pEdgeBDiv * impHorEdgeB * facEdgeB &
                    &/ rhoEdgeB - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                    &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 2, 3) * impHorEdgeB * facEdgeB &
                    &/ rhoEdgeB
              else
                gEdgeB = - jacInv / dy * pEdgeBDiv * impHorEdgeB * facEdgeB &
                    &/ rhoEdgeB + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                    &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 2, 3) * impHorEdgeB * facEdgeB &
                    &/ rhoEdgeB - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                    &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                    &* 0.25 * met(i, j, k, 2, 3) * impHorEdgeB * facEdgeB &
                    &/ rhoEdgeB
              end if

              ! G(k + 1 / 2)
              if(k == nz .and. zBoundary == "solid_wall") then
                gEdgeU = 0.0
              else
                gEdgeU = jacInv / dz * pEdgeUDiv * impVerEdgeU / rhoEdgeU
              end if

              ! G(k - 1 / 2)
              if(k == 1 .and. zBoundary == "solid_wall") then
                gEdgeD = 0.0
              else
                gEdgeD = - jacInv / dz * pEdgeDDiv * impVerEdgeD / rhoEdgeD
              end if

              ! G(i + 1 / 2, k + 1)
              if(k == nz .and. zBoundary == "solid_wall") then
                gUEdgeR = 0.0
              else
                gUEdgeR = jacInv / dz * pEdgeUDiv * impVerEdgeU &
                    &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                    &* 0.25 * met(i, j, k + 1, 1, 3) * impHorUEdgeR &
                    &* facUEdgeR / rhoUEdgeR
              end if

              ! G(i - 1 / 2, k + 1)
              if(k == nz .and. zBoundary == "solid_wall") then
                gUEdgeL = 0.0
              else
                gUEdgeL = jacInv / dz * pEdgeUDiv * impVerEdgeU &
                    &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                    &* 0.25 * met(i, j, k + 1, 1, 3) * impHorUEdgeL &
                    &* facUEdgeL / rhoUEdgeL
              end if

              ! G(j + 1 / 2, k + 1)
              if(k == nz .and. zBoundary == "solid_wall") then
                gUEdgeF = 0.0
              else
                gUEdgeF = jacInv / dz * pEdgeUDiv * impVerEdgeU &
                    &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                    &* 0.25 * met(i, j, k + 1, 2, 3) * impHorUEdgeF &
                    &* facUEdgeF / rhoUEdgeF
              end if

              ! G(j - 1 / 2, k + 1)
              if(k == nz .and. zBoundary == "solid_wall") then
                gUEdgeB = 0.0
              else
                gUEdgeB = jacInv / dz * pEdgeUDiv * impVerEdgeU &
                    &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                    &* 0.25 * met(i, j, k + 1, 2, 3) * impHorUEdgeB &
                    &* facUEdgeB / rhoUEdgeB
              end if

              ! G(i + 1 / 2, k - 1)
              if(k == 1 .and. zBoundary == "solid_wall") then
                gDEdgeR = 0.0
              else
                gDEdgeR = - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                    &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                    &* 0.25 * met(i, j, k - 1, 1, 3) * impHorDEdgeR &
                    &* facDEdgeR / rhoDEdgeR
              end if

              ! G(i - 1 / 2, k - 1)
              if(k == 1 .and. zBoundary == "solid_wall") then
                gDEdgeL = 0.0
              else
                gDEdgeL = - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                    &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                    &* 0.25 * met(i, j, k - 1, 1, 3) * impHorDEdgeL &
                    &* facDEdgeL / rhoDEdgeL
              end if

              ! G(j + 1 / 2, k - 1)
              if(k == 1 .and. zBoundary == "solid_wall") then
                gDEdgeF = 0.0
              else
                gDEdgeF = - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                    &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                    &* 0.25 * met(i, j, k - 1, 2, 3) * impHorDEdgeF &
                    &* facDEdgeF / rhoDEdgeF
              end if

              ! G(j - 1 / 2, k - 1)
              if(k == 1 .and. zBoundary == "solid_wall") then
                gDEdgeB = 0.0
              else
                gDEdgeB = - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                    &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                    &* 0.25 * met(i, j, k - 1, 2, 3) * impHorDEdgeB &
                    &* facDEdgeB / rhoDEdgeB
              end if

              ! Compute tensor elements

              ! ------------------- A(i,j,k) --------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                AC = - gEdgeR * pEdgeRGra * (1.0 / dx + 0.75 * met(i, j, k, 1, &
                    &3) / dz) + gEdgeL * pEdgeLGra * (1.0 / dx - 0.75 * met(i, &
                    &j, k, 1, 3) / dz) - gEdgeF * pEdgeFGra * (1.0 / dy + 0.75 &
                    &* met(i, j, k, 2, 3) / dz) + gEdgeB * pEdgeBGra * (1.0 &
                    &/ dy - 0.75 * met(i, j, k, 2, 3) / dz) + gEdgeU * (- &
                    &pEdgeUGra * met(i, j, k, 3, 3) / dz + 0.5 * pStratTFC(i, &
                    &j, k) * jacInv * (chris(i, j, k, 1, 1) + chris(i, j, k, &
                    &2, 2) + 2.0 * chris(i, j, k, 1, 3) * met(i, j, k, 1, 3) &
                    &+ 2.0 * chris(i, j, k, 2, 3) * met(i, j, k, 2, 3))) &
                    &- gUEdgeR * pUEdgeRGra * 0.25 * met(i, j, k, 1, 3) / dz &
                    &- gUEdgeL * pUEdgeLGra * 0.25 * met(i, j, k, 1, 3) / dz &
                    &- gUEdgeF * pUEdgeFGra * 0.25 * met(i, j, k, 2, 3) / dz &
                    &- gUEdgeB * pUEdgeBGra * 0.25 * met(i, j, k, 2, 3) / dz
              else if(k == 2 .and. zBoundary == "solid_wall") then
                AC = - gEdgeR * pEdgeRGra / dx + gEdgeL * pEdgeLGra / dx &
                    &- gEdgeF * pEdgeFGra / dy + gEdgeB * pEdgeBGra / dy &
                    &+ gEdgeU * (- pEdgeUGra * met(i, j, k, 3, 3) / dz + 0.5 &
                    &* pStratTFC(i, j, k) * jacInv * (chris(i, j, k, 1, 1) &
                    &+ chris(i, j, k, 2, 2) + 2.0 * chris(i, j, k, 1, 3) &
                    &* met(i, j, k, 1, 3) + 2.0 * chris(i, j, k, 2, 3) &
                    &* met(i, j, k, 2, 3))) + gEdgeD * (pEdgeDGra * met(i, j, &
                    &k, 3, 3) / dz + 0.5 * pStratTFC(i, j, k) * jacInv &
                    &* (chris(i, j, k, 1, 1) + chris(i, j, k, 2, 2) + 2.0 &
                    &* chris(i, j, k, 1, 3) * met(i, j, k, 1, 3) + 2.0 &
                    &* chris(i, j, k, 2, 3) * met(i, j, k, 2, 3))) - gUEdgeR &
                    &* pUEdgeRGra * 0.25 * met(i, j, k, 1, 3) / dz - gUEdgeL &
                    &* pUEdgeLGra * 0.25 * met(i, j, k, 1, 3) / dz - gUEdgeF &
                    &* pUEdgeFGra * 0.25 * met(i, j, k, 2, 3) / dz - gUEdgeB &
                    &* pUEdgeBGra * 0.25 * met(i, j, k, 2, 3) / dz + gDEdgeR &
                    &* pDEdgeRGra * met(i, j, k, 1, 3) / dz + gDEdgeL &
                    &* pDEdgeLGra * met(i, j, k, 1, 3) / dz + gDEdgeF &
                    &* pDEdgeFGra * met(i, j, k, 2, 3) / dz + gDEdgeB &
                    &* pDEdgeBGra * met(i, j, k, 2, 3) / dz
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                AC = - gEdgeR * pEdgeRGra / dx + gEdgeL * pEdgeLGra / dx &
                    &- gEdgeF * pEdgeFGra / dy + gEdgeB * pEdgeBGra / dy &
                    &+ gEdgeU * (- pEdgeUGra * met(i, j, k, 3, 3) / dz + 0.5 &
                    &* pStratTFC(i, j, k) * jacInv * (chris(i, j, k, 1, 1) &
                    &+ chris(i, j, k, 2, 2) + 2.0 * chris(i, j, k, 1, 3) &
                    &* met(i, j, k, 1, 3) + 2.0 * chris(i, j, k, 2, 3) &
                    &* met(i, j, k, 2, 3))) + gEdgeD * (pEdgeDGra * met(i, j, &
                    &k, 3, 3) / dz + 0.5 * pStratTFC(i, j, k) * jacInv &
                    &* (chris(i, j, k, 1, 1) + chris(i, j, k, 2, 2) + 2.0 &
                    &* chris(i, j, k, 1, 3) * met(i, j, k, 1, 3) + 2.0 &
                    &* chris(i, j, k, 2, 3) * met(i, j, k, 2, 3))) - gUEdgeR &
                    &* pUEdgeRGra * met(i, j, k, 1, 3) / dz - gUEdgeL &
                    &* pUEdgeLGra * met(i, j, k, 1, 3) / dz - gUEdgeF &
                    &* pUEdgeFGra * met(i, j, k, 2, 3) / dz - gUEdgeB &
                    &* pUEdgeBGra * met(i, j, k, 2, 3) / dz + gDEdgeR &
                    &* pDEdgeRGra * 0.25 * met(i, j, k, 1, 3) / dz + gDEdgeL &
                    &* pDEdgeLGra * 0.25 * met(i, j, k, 1, 3) / dz + gDEdgeF &
                    &* pDEdgeFGra * 0.25 * met(i, j, k, 2, 3) / dz + gDEdgeB &
                    &* pDEdgeBGra * 0.25 * met(i, j, k, 2, 3) / dz
              else if(k == nz .and. zBoundary == "solid_wall") then
                AC = - gEdgeR * pEdgeRGra * (1.0 / dx - 0.75 * met(i, j, k, 1, &
                    &3) / dz) + gEdgeL * pEdgeLGra * (1.0 / dx + 0.75 * met(i, &
                    &j, k, 1, 3) / dz) - gEdgeF * pEdgeFGra * (1.0 / dy - 0.75 &
                    &* met(i, j, k, 2, 3) / dz) + gEdgeB * pEdgeBGra * (1.0 &
                    &/ dy + 0.75 * met(i, j, k, 2, 3) / dz) + gEdgeD &
                    &* (pEdgeDGra * met(i, j, k, 3, 3) / dz + 0.5 &
                    &* pStratTFC(i, j, k) * jacInv * (chris(i, j, k, 1, 1) &
                    &+ chris(i, j, k, 2, 2) + 2.0 * chris(i, j, k, 1, 3) &
                    &* met(i, j, k, 1, 3) + 2.0 * chris(i, j, k, 2, 3) &
                    &* met(i, j, k, 2, 3))) + gDEdgeR * pDEdgeRGra * 0.25 &
                    &* met(i, j, k, 1, 3) / dz + gDEdgeL * pDEdgeLGra * 0.25 &
                    &* met(i, j, k, 1, 3) / dz + gDEdgeF * pDEdgeFGra * 0.25 &
                    &* met(i, j, k, 2, 3) / dz + gDEdgeB * pDEdgeBGra * 0.25 &
                    &* met(i, j, k, 2, 3) / dz
              else
                AC = - gEdgeR * pEdgeRGra / dx + gEdgeL * pEdgeLGra / dx &
                    &- gEdgeF * pEdgeFGra / dy + gEdgeB * pEdgeBGra / dy &
                    &+ gEdgeU * (- pEdgeUGra * met(i, j, k, 3, 3) / dz + 0.5 &
                    &* pStratTFC(i, j, k) * jacInv * (chris(i, j, k, 1, 1) &
                    &+ chris(i, j, k, 2, 2) + 2.0 * chris(i, j, k, 1, 3) &
                    &* met(i, j, k, 1, 3) + 2.0 * chris(i, j, k, 2, 3) &
                    &* met(i, j, k, 2, 3))) + gEdgeD * (pEdgeDGra * met(i, j, &
                    &k, 3, 3) / dz + 0.5 * pStratTFC(i, j, k) * jacInv &
                    &* (chris(i, j, k, 1, 1) + chris(i, j, k, 2, 2) + 2.0 &
                    &* chris(i, j, k, 1, 3) * met(i, j, k, 1, 3) + 2.0 &
                    &* chris(i, j, k, 2, 3) * met(i, j, k, 2, 3))) - gUEdgeR &
                    &* pUEdgeRGra * 0.25 * met(i, j, k, 1, 3) / dz - gUEdgeL &
                    &* pUEdgeLGra * 0.25 * met(i, j, k, 1, 3) / dz - gUEdgeF &
                    &* pUEdgeFGra * 0.25 * met(i, j, k, 2, 3) / dz - gUEdgeB &
                    &* pUEdgeBGra * 0.25 * met(i, j, k, 2, 3) / dz + gDEdgeR &
                    &* pDEdgeRGra * 0.25 * met(i, j, k, 1, 3) / dz + gDEdgeL &
                    &* pDEdgeLGra * 0.25 * met(i, j, k, 1, 3) / dz + gDEdgeF &
                    &* pDEdgeFGra * 0.25 * met(i, j, k, 2, 3) / dz + gDEdgeB &
                    &* pDEdgeBGra * 0.25 * met(i, j, k, 2, 3) / dz
              end if

              if(model == "compressible") then
                AC = AC - (dPdPi / dt) / (dt * Rsp / kappa)
              end if

              ! ------------------ A(i+1,j,k) -------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                AR = gEdgeR * pEdgeRGra * (1.0 / dx - 0.75 * met(i + 1, j, k, &
                    &1, 3) / dz) + gEdgeU * pEdgeUGra * 0.25 * met(i + 1, j, &
                    &k, 1, 3) / dx - gUEdgeR * pUEdgeRGra * 0.25 * met(i + 1, &
                    &j, k, 1, 3) / dz
              else if(k == 2 .and. zBoundary == "solid_wall") then
                AR = gEdgeR * pEdgeRGra / dx + gEdgeU * pEdgeUGra * 0.25 &
                    &* met(i + 1, j, k, 1, 3) / dx + gEdgeD * pEdgeDGra * 0.25 &
                    &* met(i + 1, j, k, 1, 3) / dx - gUEdgeR * pUEdgeRGra &
                    &* 0.25 * met(i + 1, j, k, 1, 3) / dz + gDEdgeR &
                    &* pDEdgeRGra * met(i + 1, j, k, 1, 3) / dz
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                AR = gEdgeR * pEdgeRGra / dx + gEdgeU * pEdgeUGra * 0.25 &
                    &* met(i + 1, j, k, 1, 3) / dx + gEdgeD * pEdgeDGra * 0.25 &
                    &* met(i + 1, j, k, 1, 3) / dx - gUEdgeR * pUEdgeRGra &
                    &* met(i + 1, j, k, 1, 3) / dz + gDEdgeR * pDEdgeRGra &
                    &* 0.25 * met(i + 1, j, k, 1, 3) / dz
              else if(k == nz .and. zBoundary == "solid_wall") then
                AR = gEdgeR * pEdgeRGra * (1.0 / dx + 0.75 * met(i + 1, j, k, &
                    &1, 3) / dz) + gEdgeD * pEdgeDGra * 0.25 * met(i + 1, j, &
                    &k, 1, 3) / dx + gDEdgeR * pDEdgeRGra * 0.25 * met(i + 1, &
                    &j, k, 1, 3) / dz
              else
                AR = gEdgeR * pEdgeRGra / dx + gEdgeU * pEdgeUGra * 0.25 &
                    &* met(i + 1, j, k, 1, 3) / dx + gEdgeD * pEdgeDGra * 0.25 &
                    &* met(i + 1, j, k, 1, 3) / dx - gUEdgeR * pUEdgeRGra &
                    &* 0.25 * met(i + 1, j, k, 1, 3) / dz + gDEdgeR &
                    &* pDEdgeRGra * 0.25 * met(i + 1, j, k, 1, 3) / dz
              end if

              ! ------------------ A(i-1,j,k) -------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                AL = - gEdgeL * pEdgeLGra * (1.0 / dx + 0.75 * met(i - 1, j, &
                    &k, 1, 3) / dz) - gEdgeU * pEdgeUGra * 0.25 * met(i - 1, &
                    &j, k, 1, 3) / dx - gUEdgeL * pUEdgeLGra * 0.25 * met(i &
                    &- 1, j, k, 1, 3) / dz
              else if(k == 2 .and. zBoundary == "solid_wall") then
                AL = - gEdgeL * pEdgeLGra / dx - gEdgeU * pEdgeUGra * 0.25 &
                    &* met(i - 1, j, k, 1, 3) / dx - gEdgeD * pEdgeDGra * 0.25 &
                    &* met(i - 1, j, k, 1, 3) / dx - gUEdgeL * pUEdgeLGra &
                    &* 0.25 * met(i - 1, j, k, 1, 3) / dz + gDEdgeL &
                    &* pDEdgeLGra * met(i - 1, j, k, 1, 3) / dz
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                AL = - gEdgeL * pEdgeLGra / dx - gEdgeU * pEdgeUGra * 0.25 &
                    &* met(i - 1, j, k, 1, 3) / dx - gEdgeD * pEdgeDGra * 0.25 &
                    &* met(i - 1, j, k, 1, 3) / dx - gUEdgeL * pUEdgeLGra &
                    &* met(i - 1, j, k, 1, 3) / dz + gDEdgeL * pDEdgeLGra &
                    &* 0.25 * met(i - 1, j, k, 1, 3) / dz
              else if(k == nz .and. zBoundary == "solid_wall") then
                AL = - gEdgeL * pEdgeLGra * (1.0 / dx - 0.75 * met(i - 1, j, &
                    &k, 1, 3) / dz) - gEdgeD * pEdgeDGra * 0.25 * met(i - 1, &
                    &j, k, 1, 3) / dx + gDEdgeL * pDEdgeLGra * 0.25 * met(i &
                    &- 1, j, k, 1, 3) / dz
              else
                AL = - gEdgeL * pEdgeLGra / dx - gEdgeU * pEdgeUGra * 0.25 &
                    &* met(i - 1, j, k, 1, 3) / dx - gEdgeD * pEdgeDGra * 0.25 &
                    &* met(i - 1, j, k, 1, 3) / dx - gUEdgeL * pUEdgeLGra &
                    &* 0.25 * met(i - 1, j, k, 1, 3) / dz + gDEdgeL &
                    &* pDEdgeLGra * 0.25 * met(i - 1, j, k, 1, 3) / dz
              end if

              ! ------------------ A(i,j+1,k) -------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                AF = gEdgeF * pEdgeFGra * (1.0 / dy - 0.75 * met(i, j + 1, k, &
                    &2, 3) / dz) + gEdgeU * pEdgeUGra * 0.25 * met(i, j + 1, &
                    &k, 2, 3) / dy - gUEdgeF * pUEdgeFGra * 0.25 * met(i, j &
                    &+ 1, k, 2, 3) / dz
              else if(k == 2 .and. zBoundary == "solid_wall") then
                AF = gEdgeF * pEdgeFGra / dy + gEdgeU * pEdgeUGra * 0.25 &
                    &* met(i, j + 1, k, 2, 3) / dy + gEdgeD * pEdgeDGra * 0.25 &
                    &* met(i, j + 1, k, 2, 3) / dy - gUEdgeF * pUEdgeFGra &
                    &* 0.25 * met(i, j + 1, k, 2, 3) / dz + gDEdgeF &
                    &* pDEdgeFGra * met(i, j + 1, k, 2, 3) / dz
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                AF = gEdgeF * pEdgeFGra / dy + gEdgeU * pEdgeUGra * 0.25 &
                    &* met(i, j + 1, k, 2, 3) / dy + gEdgeD * pEdgeDGra * 0.25 &
                    &* met(i, j + 1, k, 2, 3) / dy - gUEdgeF * pUEdgeFGra &
                    &* met(i, j + 1, k, 2, 3) / dz + gDEdgeF * pDEdgeFGra &
                    &* 0.25 * met(i, j + 1, k, 2, 3) / dz
              else if(k == nz .and. zBoundary == "solid_wall") then
                AF = gEdgeF * pEdgeFGra * (1.0 / dy + 0.75 * met(i, j + 1, k, &
                    &2, 3) / dz) + gEdgeD * pEdgeDGra * 0.25 * met(i, j + 1, &
                    &k, 2, 3) / dy + gDEdgeF * pDEdgeFGra * 0.25 * met(i, j &
                    &+ 1, k, 2, 3) / dz
              else
                AF = gEdgeF * pEdgeFGra / dy + gEdgeU * pEdgeUGra * 0.25 &
                    &* met(i, j + 1, k, 2, 3) / dy + gEdgeD * pEdgeDGra * 0.25 &
                    &* met(i, j + 1, k, 2, 3) / dy - gUEdgeF * pUEdgeFGra &
                    &* 0.25 * met(i, j + 1, k, 2, 3) / dz + gDEdgeF &
                    &* pDEdgeFGra * 0.25 * met(i, j + 1, k, 2, 3) / dz
              end if

              ! ------------------ A(i,j-1,k) -------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                AB = - gEdgeB * pEdgeBGra * (1.0 / dy + 0.75 * met(i, j - 1, &
                    &k, 2, 3) / dz) - gEdgeU * pEdgeUGra * 0.25 * met(i, j &
                    &- 1, k, 2, 3) / dy - gUEdgeB * pUEdgeBGra * 0.25 * met(i, &
                    &j - 1, k, 2, 3) / dz
              else if(k == 2 .and. zBoundary == "solid_wall") then
                AB = - gEdgeB * pEdgeBGra / dy - gEdgeU * pEdgeUGra * 0.25 &
                    &* met(i, j - 1, k, 2, 3) / dy - gEdgeD * pEdgeDGra * 0.25 &
                    &* met(i, j - 1, k, 2, 3) / dy - gUEdgeB * pUEdgeBGra &
                    &* 0.25 * met(i, j - 1, k, 2, 3) / dz + gDEdgeB &
                    &* pDEdgeBGra * met(i, j - 1, k, 2, 3) / dz
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                AB = - gEdgeB * pEdgeBGra / dy - gEdgeU * pEdgeUGra * 0.25 &
                    &* met(i, j - 1, k, 2, 3) / dy - gEdgeD * pEdgeDGra * 0.25 &
                    &* met(i, j - 1, k, 2, 3) / dy - gUEdgeB * pUEdgeBGra &
                    &* met(i, j - 1, k, 2, 3) / dz + gDEdgeB * pDEdgeBGra &
                    &* 0.25 * met(i, j - 1, k, 2, 3) / dz
              else if(k == nz .and. zBoundary == "solid_wall") then
                AB = - gEdgeB * pEdgeBGra * (1.0 / dy - 0.75 * met(i, j - 1, &
                    &k, 2, 3) / dz) - gEdgeD * pEdgeDGra * 0.25 * met(i, j &
                    &- 1, k, 2, 3) / dy + gDEdgeB * pDEdgeBGra * 0.25 * met(i, &
                    &j - 1, k, 2, 3) / dz
              else
                AB = - gEdgeB * pEdgeBGra / dy - gEdgeU * pEdgeUGra * 0.25 &
                    &* met(i, j - 1, k, 2, 3) / dy - gEdgeD * pEdgeDGra * 0.25 &
                    &* met(i, j - 1, k, 2, 3) / dy - gUEdgeB * pUEdgeBGra &
                    &* 0.25 * met(i, j - 1, k, 2, 3) / dz + gDEdgeB &
                    &* pDEdgeBGra * 0.25 * met(i, j - 1, k, 2, 3) / dz
              end if

              ! ------------------ A(i,j,k+1) -------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                AU = gEdgeR * pEdgeRGra * met(i, j, k + 1, 1, 3) / dz + gEdgeL &
                    &* pEdgeLGra * met(i, j, k + 1, 1, 3) / dz + gEdgeF &
                    &* pEdgeFGra * met(i, j, k + 1, 2, 3) / dz + gEdgeB &
                    &* pEdgeBGra * met(i, j, k + 1, 2, 3) / dz + gEdgeU &
                    &* (pEdgeUGra * met(i, j, k + 1, 3, 3) / dz + 0.5 &
                    &* pStratTFC(i, j, k + 1) / jac(i, j, k + 1) * (chris(i, &
                    &j, k + 1, 1, 1) + chris(i, j, k + 1, 2, 2) + 2.0 &
                    &* chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, 1, 3) + 2.0 &
                    &* chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, 2, 3))) &
                    &- gUEdgeR * pUEdgeRGra / dx + gUEdgeL * pUEdgeLGra / dx &
                    &- gUEdgeF * pUEdgeFGra / dy + gUEdgeB * pUEdgeBGra / dy
              else if(k == 2 .and. zBoundary == "solid_wall") then
                AU = gEdgeR * pEdgeRGra * 0.25 * met(i, j, k + 1, 1, 3) / dz &
                    &+ gEdgeL * pEdgeLGra * 0.25 * met(i, j, k + 1, 1, 3) / dz &
                    &+ gEdgeF * pEdgeFGra * 0.25 * met(i, j, k + 1, 2, 3) / dz &
                    &+ gEdgeB * pEdgeBGra * 0.25 * met(i, j, k + 1, 2, 3) / dz &
                    &+ gEdgeU * (pEdgeUGra * met(i, j, k + 1, 3, 3) / dz + 0.5 &
                    &* pStratTFC(i, j, k + 1) / jac(i, j, k + 1) * (chris(i, &
                    &j, k + 1, 1, 1) + chris(i, j, k + 1, 2, 2) + 2.0 &
                    &* chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, 1, 3) + 2.0 &
                    &* chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, 2, 3))) &
                    &- gUEdgeR * pUEdgeRGra / dx + gUEdgeL * pUEdgeLGra / dx &
                    &- gUEdgeF * pUEdgeFGra / dy + gUEdgeB * pUEdgeBGra / dy &
                    &- gDEdgeR * pDEdgeRGra * 0.25 * met(i, j, k + 1, 1, 3) &
                    &/ dz - gDEdgeL * pDEdgeLGra * 0.25 * met(i, j, k + 1, 1, &
                    &3) / dz - gDEdgeF * pDEdgeFGra * 0.25 * met(i, j, k + 1, &
                    &2, 3) / dz - gDEdgeB * pDEdgeBGra * 0.25 * met(i, j, k &
                    &+ 1, 2, 3)
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                AU = gEdgeR * pEdgeRGra * 0.25 * met(i, j, k + 1, 1, 3) / dz &
                    &+ gEdgeL * pEdgeLGra * 0.25 * met(i, j, k + 1, 1, 3) / dz &
                    &+ gEdgeF * pEdgeFGra * 0.25 * met(i, j, k + 1, 2, 3) / dz &
                    &+ gEdgeB * pEdgeBGra * 0.25 * met(i, j, k + 1, 2, 3) / dz &
                    &+ gEdgeU * (pEdgeUGra * met(i, j, k + 1, 3, 3) / dz + 0.5 &
                    &* pStratTFC(i, j, k + 1) / jac(i, j, k + 1) * (chris(i, &
                    &j, k + 1, 1, 1) + chris(i, j, k + 1, 2, 2) + 2.0 &
                    &* chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, 1, 3) + 2.0 &
                    &* chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, 2, 3))) &
                    &- gUEdgeR * pUEdgeRGra * (1.0 / dx - 0.75 * met(i, j, k &
                    &+ 1, 1, 3) / dz) + gUEdgeL * pUEdgeLGra * (1.0 / dx &
                    &+ 0.75 * met(i, j, k + 1, 1, 3) / dz) - gUEdgeF &
                    &* pUEdgeFGra * (1.0 / dy - 0.75 * met(i, j, k + 1, 2, 3) &
                    &/ dz) + gUEdgeB * pUEdgeBGra * (1.0 / dy + 0.75 * met(i, &
                    &j, k + 1, 2, 3) / dz)
              else if(k == nz .and. zBoundary == "solid_wall") then
                AU = 0.0
              else
                AU = gEdgeR * pEdgeRGra * 0.25 * met(i, j, k + 1, 1, 3) / dz &
                    &+ gEdgeL * pEdgeLGra * 0.25 * met(i, j, k + 1, 1, 3) / dz &
                    &+ gEdgeF * pEdgeFGra * 0.25 * met(i, j, k + 1, 2, 3) / dz &
                    &+ gEdgeB * pEdgeBGra * 0.25 * met(i, j, k + 1, 2, 3) / dz &
                    &+ gEdgeU * (pEdgeUGra * met(i, j, k + 1, 3, 3) / dz + 0.5 &
                    &* pStratTFC(i, j, k + 1) / jac(i, j, k + 1) * (chris(i, &
                    &j, k + 1, 1, 1) + chris(i, j, k + 1, 2, 2) + 2.0 &
                    &* chris(i, j, k + 1, 1, 3) * met(i, j, k + 1, 1, 3) + 2.0 &
                    &* chris(i, j, k + 1, 2, 3) * met(i, j, k + 1, 2, 3))) &
                    &- gUEdgeR * pUEdgeRGra / dx + gUEdgeL * pUEdgeLGra / dx &
                    &- gUEdgeF * pUEdgeFGra / dy + gUEdgeB * pUEdgeBGra / dy
              end if

              ! ------------------ A(i,j,k-1) -------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                AD = 0.0
              else if(k == 2 .and. zBoundary == "solid_wall") then
                AD = - gEdgeR * pEdgeRGra * 0.25 * met(i, j, k - 1, 1, 3) / dz &
                    &- gEdgeL * pEdgeLGra * 0.25 * met(i, j, k - 1, 1, 3) / dz &
                    &- gEdgeF * pEdgeFGra * 0.25 * met(i, j, k - 1, 2, 3) / dz &
                    &- gEdgeB * pEdgeBGra * 0.25 * met(i, j, k - 1, 2, 3) / dz &
                    &+ gEdgeD * (- pEdgeDGra * met(i, j, k - 1, 3, 3) / dz &
                    &+ 0.5 * pStratTFC(i, j, k - 1) / jac(i, j, k - 1) &
                    &* (chris(i, j, k - 1, 1, 1) + chris(i, j, k - 1, 2, 2) &
                    &+ 2.0 * chris(i, j, k - 1, 1, 3) * met(i, j, k - 1, 1, 3) &
                    &+ 2.0 * chris(i, j, k - 1, 2, 3) * met(i, j, k - 1, 2, &
                    &3))) - gDEdgeR * pDEdgeRGra * (1.0 / dx + 0.75 * met(i, &
                    &j, k - 1, 1, 3) / dz) + gDEdgeL * pDEdgeLGra * (1.0 / dx &
                    &- 0.75 * met(i, j, k - 1, 1, 3) / dz) - gDEdgeF &
                    &* pDEdgeFGra * (1.0 / dy + 0.75 * met(i, j, k - 1, 2, 3) &
                    &/ dz) + gDEdgeB * pDEdgeBGra * (1.0 / dy - 0.75 * met(i, &
                    &j, k - 1, 2, 3) / dz)
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                AD = - gEdgeR * pEdgeRGra * 0.25 * met(i, j, k - 1, 1, 3) / dz &
                    &- gEdgeL * pEdgeLGra * 0.25 * met(i, j, k - 1, 1, 3) / dz &
                    &- gEdgeF * pEdgeFGra * 0.25 * met(i, j, k - 1, 2, 3) / dz &
                    &- gEdgeB * pEdgeBGra * 0.25 * met(i, j, k - 1, 2, 3) / dz &
                    &+ gEdgeD * (- pEdgeDGra * met(i, j, k - 1, 3, 3) / dz &
                    &+ 0.5 * pStratTFC(i, j, k - 1) / jac(i, j, k - 1) &
                    &* (chris(i, j, k - 1, 1, 1) + chris(i, j, k - 1, 2, 2) &
                    &+ 2.0 * chris(i, j, k - 1, 1, 3) * met(i, j, k - 1, 1, 3) &
                    &+ 2.0 * chris(i, j, k - 1, 2, 3) * met(i, j, k - 1, 2, &
                    &3))) - gDEdgeR * pDEdgeRGra / dx + gDEdgeL * pDEdgeLGra &
                    &/ dx - gDEdgeF * pDEdgeFGra / dy + gDEdgeB * pDEdgeBGra &
                    &/ dy + gUEdgeR * pUEdgeRGra * 0.25 * met(i, j, k - 1, 1, &
                    &3) / dz + gUEdgeL * pUEdgeLGra * 0.25 * met(i, j, k - 1, &
                    &1, 3) / dz + gUEdgeF * pUEdgeFGra * 0.25 * met(i, j, k &
                    &- 1, 2, 3) / dz + gUEdgeB * pUEdgeBGra * 0.25 * met(i, j, &
                    &k - 1, 2, 3) / dz
              else if(k == nz .and. zBoundary == "solid_wall") then
                AD = - gEdgeR * pEdgeRGra * met(i, j, k - 1, 1, 3) / dz &
                    &- gEdgeL * pEdgeLGra * met(i, j, k - 1, 1, 3) / dz &
                    &- gEdgeF * pEdgeFGra * met(i, j, k - 1, 2, 3) / dz &
                    &- gEdgeB * pEdgeBGra * met(i, j, k - 1, 2, 3) / dz &
                    &+ gEdgeD * (- pEdgeDGra * met(i, j, k - 1, 3, 3) / dz &
                    &+ 0.5 * pStratTFC(i, j, k - 1) / jac(i, j, k - 1) &
                    &* (chris(i, j, k - 1, 1, 1) + chris(i, j, k - 1, 2, 2) &
                    &+ 2.0 * chris(i, j, k - 1, 1, 3) * met(i, j, k - 1, 1, 3) &
                    &+ 2.0 * chris(i, j, k - 1, 2, 3) * met(i, j, k - 1, 2, &
                    &3))) - gDEdgeR * pDEdgeRGra / dx + gDEdgeL * pDEdgeLGra &
                    &/ dx - gDEdgeF * pDEdgeFGra / dy + gDEdgeB * pDEdgeBGra &
                    &/ dy
              else
                AD = - gEdgeR * pEdgeRGra * 0.25 * met(i, j, k - 1, 1, 3) / dz &
                    &- gEdgeL * pEdgeLGra * 0.25 * met(i, j, k - 1, 1, 3) / dz &
                    &- gEdgeF * pEdgeFGra * 0.25 * met(i, j, k - 1, 2, 3) / dz &
                    &- gEdgeB * pEdgeBGra * 0.25 * met(i, j, k - 1, 2, 3) / dz &
                    &+ gEdgeD * (- pEdgeDGra * met(i, j, k - 1, 3, 3) / dz &
                    &+ 0.5 * pStratTFC(i, j, k - 1) / jac(i, j, k - 1) &
                    &* (chris(i, j, k - 1, 1, 1) + chris(i, j, k - 1, 2, 2) &
                    &+ 2.0 * chris(i, j, k - 1, 1, 3) * met(i, j, k - 1, 1, 3) &
                    &+ 2.0 * chris(i, j, k - 1, 2, 3) * met(i, j, k - 1, 2, &
                    &3))) - gDEdgeR * pDEdgeRGra / dx + gDEdgeL * pDEdgeLGra &
                    &/ dx - gDEdgeF * pDEdgeFGra / dy + gDEdgeB * pDEdgeBGra &
                    &/ dy
              end if

              ! ----------------- A(i+1,j,k+1) ------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                ARU = gEdgeR * pEdgeRGra * met(i + 1, j, k + 1, 1, 3) / dz &
                    &+ gEdgeU * pEdgeUGra * 0.25 * met(i + 1, j, k + 1, 1, 3) &
                    &/ dx + gUEdgeR * pUEdgeRGra / dx
              else if(k == 2 .and. zBoundary == "solid_wall") then
                ARU = gEdgeR * pEdgeRGra * 0.25 * met(i + 1, j, k + 1, 1, 3) &
                    &/ dz + gEdgeU * pEdgeUGra * 0.25 * met(i + 1, j, k + 1, &
                    &1, 3) / dx + gUEdgeR * pUEdgeRGra / dx - gDEdgeR &
                    &* pDEdgeRGra * 0.25 * met(i + 1, j, k + 1, 1, 3) / dz
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                ARU = gEdgeR * pEdgeRGra * 0.25 * met(i + 1, j, k + 1, 1, 3) &
                    &/ dz + gEdgeU * pEdgeUGra * 0.25 * met(i + 1, j, k + 1, &
                    &1, 3) / dx + gUEdgeR * pUEdgeRGra * (1.0 / dx + 0.75 &
                    &* met(i + 1, j, k + 1, 1, 3) / dz)
              else if(k == nz .and. zBoundary == "solid_wall") then
                ARU = 0.0
              else
                ARU = gEdgeR * pEdgeRGra * 0.25 * met(i + 1, j, k + 1, 1, 3) &
                    &/ dz + gEdgeU * pEdgeUGra * 0.25 * met(i + 1, j, k + 1, &
                    &1, 3) / dx + gUEdgeR * pUEdgeRGra / dx
              end if

              ! ----------------- A(i+1,j,k-1) ------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                ARD = 0.0
              else if(k == 2 .and. zBoundary == "solid_wall") then
                ARD = - gEdgeR * pEdgeRGra * 0.25 * met(i + 1, j, k - 1, 1, 3) &
                    &/ dz + gEdgeD * pEdgeDGra * 0.25 * met(i + 1, j, k - 1, &
                    &1, 3) / dx + gDEdgeR * pDEdgeRGra * (1.0 / dx - 0.75 &
                    &* met(i + 1, j, k - 1, 1, 3) / dz)
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                ARD = - gEdgeR * pEdgeRGra * 0.25 * met(i + 1, j, k - 1, 1, 3) &
                    &/ dz + gEdgeD * pEdgeDGra * 0.25 * met(i + 1, j, k - 1, &
                    &1, 3) / dx + gDEdgeR * pDEdgeRGra / dx + gUEdgeR &
                    &* pUEdgeRGra * 0.25 * met(i + 1, j, k - 1, 1, 3) / dz
              else if(k == nz .and. zBoundary == "solid_wall") then
                ARD = - gEdgeR * pEdgeRGra * met(i + 1, j, k - 1, 1, 3) / dz &
                    &+ gEdgeD * pEdgeDGra * 0.25 * met(i + 1, j, k - 1, 1, 3) &
                    &/ dx + gDEdgeR * pDEdgeRGra / dx
              else
                ARD = - gEdgeR * pEdgeRGra * 0.25 * met(i + 1, j, k - 1, 1, 3) &
                    &/ dz + gEdgeD * pEdgeDGra * 0.25 * met(i + 1, j, k - 1, &
                    &1, 3) / dx + gDEdgeR * pDEdgeRGra / dx
              end if

              ! ----------------- A(i-1,j,k+1) ------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                ALU = gEdgeL * pEdgeLGra * met(i - 1, j, k + 1, 1, 3) / dz &
                    &- gEdgeU * pEdgeUGra * 0.25 * met(i - 1, j, k + 1, 1, 3) &
                    &/ dx - gUEdgeL * pUEdgeLGra / dx
              else if(k == 2 .and. zBoundary == "solid_wall") then
                ALU = gEdgeL * pEdgeLGra * 0.25 * met(i - 1, j, k + 1, 1, 3) &
                    &/ dz - gEdgeU * pEdgeUGra * 0.25 * met(i - 1, j, k + 1, &
                    &1, 3) / dx - gUEdgeL * pUEdgeLGra / dx - gDEdgeL &
                    &* pDEdgeLGra * 0.25 * met(i - 1, j, k + 1, 1, 3) / dz
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                ALU = gEdgeL * pEdgeLGra * 0.25 * met(i - 1, j, k + 1, 1, 3) &
                    &/ dz - gEdgeU * pEdgeUGra * 0.25 * met(i - 1, j, k + 1, &
                    &1, 3) / dx - gUEdgeL * pUEdgeLGra * (1.0 / dx - 0.75 &
                    &* met(i - 1, j, k + 1, 1, 3) / dz)
              else if(k == nz .and. zBoundary == "solid_wall") then
                ALU = 0.0
              else
                ALU = gEdgeL * pEdgeLGra * 0.25 * met(i - 1, j, k + 1, 1, 3) &
                    &/ dz - gEdgeU * pEdgeUGra * 0.25 * met(i - 1, j, k + 1, &
                    &1, 3) / dx - gUEdgeL * pUEdgeLGra / dx
              end if

              ! ----------------- A(i-1,j,k-1) ------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                ALD = 0.0
              else if(k == 2 .and. zBoundary == "solid_wall") then
                ALD = - gEdgeL * pEdgeLGra * 0.25 * met(i - 1, j, k - 1, 1, 3) &
                    &/ dz - gEdgeD * pEdgeDGra * 0.25 * met(i - 1, j, k - 1, &
                    &1, 3) / dx - gDEdgeL * pDEdgeLGra * (1.0 / dx + 0.75 &
                    &* met(i - 1, j, k - 1, 1, 3) / dz)
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                ALD = - gEdgeL * pEdgeLGra * 0.25 * met(i - 1, j, k - 1, 1, 3) &
                    &/ dz - gEdgeD * pEdgeDGra * 0.25 * met(i - 1, j, k - 1, &
                    &1, 3) / dx - gDEdgeL * pDEdgeLGra / dx + gUEdgeL &
                    &* pUEdgeLGra * 0.25 * met(i - 1, j, k - 1, 1, 3) / dz
              else if(k == nz .and. zBoundary == "solid_wall") then
                ALD = - gEdgeL * pEdgeLGra * met(i - 1, j, k - 1, 1, 3) / dz &
                    &- gEdgeD * pEdgeDGra * 0.25 * met(i - 1, j, k - 1, 1, 3) &
                    &/ dx - gDEdgeL * pDEdgeLGra / dx
              else
                ALD = - gEdgeL * pEdgeLGra * 0.25 * met(i - 1, j, k - 1, 1, 3) &
                    &/ dz - gEdgeD * pEdgeDGra * 0.25 * met(i - 1, j, k - 1, &
                    &1, 3) / dx - gDEdgeL * pDEdgeLGra / dx
              end if

              ! ----------------- A(i,j+1,k+1) ------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                AFU = gEdgeF * pEdgeFGra * met(i, j + 1, k + 1, 2, 3) / dz &
                    &+ gEdgeU * pEdgeUGra * 0.25 * met(i, j + 1, k + 1, 2, 3) &
                    &/ dy + gUEdgeF * pUEdgeFGra / dy
              else if(k == 2 .and. zBoundary == "solid_wall") then
                AFU = gEdgeF * pEdgeFGra * 0.25 * met(i, j + 1, k + 1, 2, 3) &
                    &/ dz + gEdgeU * pEdgeUGra * 0.25 * met(i, j + 1, k + 1, &
                    &2, 3) / dy + gUEdgeF * pUEdgeFGra / dy - gDEdgeF &
                    &* pDEdgeFGra * 0.25 * met(i, j + 1, k + 1, 2, 3) / dz
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                AFU = gEdgeF * pEdgeFGra * 0.25 * met(i, j + 1, k + 1, 2, 3) &
                    &/ dz + gEdgeU * pEdgeUGra * 0.25 * met(i, j + 1, k + 1, &
                    &2, 3) / dy + gUEdgeF * pUEdgeFGra * (1.0 / dy + 0.75 &
                    &* met(i, j + 1, k + 1, 2, 3) / dz)
              else if(k == nz .and. zBoundary == "solid_wall") then
                AFU = 0.0
              else
                AFU = gEdgeF * pEdgeFGra * 0.25 * met(i, j + 1, k + 1, 2, 3) &
                    &/ dz + gEdgeU * pEdgeUGra * 0.25 * met(i, j + 1, k + 1, &
                    &2, 3) / dy + gUEdgeF * pUEdgeFGra / dy
              end if

              ! ----------------- A(i,j+1,k-1) ------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                AFD = 0.0
              else if(k == 2 .and. zBoundary == "solid_wall") then
                AFD = - gEdgeF * pEdgeFGra * 0.25 * met(i, j + 1, k - 1, 2, 3) &
                    &/ dz + gEdgeD * pEdgeDGra * 0.25 * met(i, j + 1, k - 1, &
                    &2, 3) / dy + gDEdgeF * pDEdgeFGra * (1.0 / dy - 0.75 &
                    &* met(i, j + 1, k - 1, 2, 3) / dz)
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                AFD = - gEdgeF * pEdgeFGra * 0.25 * met(i, j + 1, k - 1, 2, 3) &
                    &/ dz + gEdgeD * pEdgeDGra * 0.25 * met(i, j + 1, k - 1, &
                    &2, 3) / dy + gDEdgeF * pDEdgeFGra / dy + gUEdgeF &
                    &* pUEdgeFGra * 0.25 * met(i, j + 1, k - 1, 2, 3) / dz
              else if(k == nz .and. zBoundary == "solid_wall") then
                AFD = - gEdgeF * pEdgeFGra * met(i, j + 1, k - 1, 2, 3) / dz &
                    &+ gEdgeD * pEdgeDGra * 0.25 * met(i, j + 1, k - 1, 2, 3) &
                    &/ dy + gDEdgeF * pDEdgeFGra / dy
              else
                AFD = - gEdgeF * pEdgeFGra * 0.25 * met(i, j + 1, k - 1, 2, 3) &
                    &/ dz + gEdgeD * pEdgeDGra * 0.25 * met(i, j + 1, k - 1, &
                    &2, 3) / dy + gDEdgeF * pDEdgeFGra / dy
              end if

              ! ----------------- A(i,j-1,k+1) ------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                ABU = gEdgeB * pEdgeBGra * met(i, j - 1, k + 1, 2, 3) / dz &
                    &- gEdgeU * pEdgeUGra * 0.25 * met(i, j - 1, k + 1, 2, 3) &
                    &/ dy - gUEdgeB * pUEdgeBGra / dy
              else if(k == 2 .and. zBoundary == "solid_wall") then
                ABU = gEdgeB * pEdgeBGra * 0.25 * met(i, j - 1, k + 1, 2, 3) &
                    &/ dz - gEdgeU * pEdgeUGra * 0.25 * met(i, j - 1, k + 1, &
                    &2, 3) / dy - gUEdgeB * pUEdgeBGra / dy - gDEdgeB &
                    &* pDEdgeBGra * 0.25 * met(i, j - 1, k + 1, 2, 3) / dz
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                ABU = gEdgeB * pEdgeBGra * 0.25 * met(i, j - 1, k + 1, 2, 3) &
                    &/ dz - gEdgeU * pEdgeUGra * 0.25 * met(i, j - 1, k + 1, &
                    &2, 3) / dy - gUEdgeB * pUEdgeBGra * (1.0 / dy - 0.75 &
                    &* met(i, j - 1, k + 1, 2, 3) / dz)
              else if(k == nz .and. zBoundary == "solid_wall") then
                ABU = 0.0
              else
                ABU = gEdgeB * pEdgeBGra * 0.25 * met(i, j - 1, k + 1, 2, 3) &
                    &/ dz - gEdgeU * pEdgeUGra * 0.25 * met(i, j - 1, k + 1, &
                    &2, 3) / dy - gUEdgeB * pUEdgeBGra / dy
              end if

              ! ----------------- A(i,j-1,k-1) ------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                ABD = 0.0
              else if(k == 2 .and. zBoundary == "solid_wall") then
                ABD = - gEdgeB * pEdgeBGra * 0.25 * met(i, j - 1, k - 1, 2, 3) &
                    &/ dz - gEdgeD * pEdgeDGra * 0.25 * met(i, j - 1, k - 1, &
                    &2, 3) / dy - gDEdgeB * pDEdgeBGra * (1.0 / dy + 0.75 &
                    &* met(i, j - 1, k - 1, 2, 3) / dz)
              else if(k == nz - 1 .and. zBoundary == "solid_wall") then
                ABD = - gEdgeB * pEdgeBGra * 0.25 * met(i, j - 1, k - 1, 2, 3) &
                    &/ dz - gEdgeD * pEdgeDGra * 0.25 * met(i, j - 1, k - 1, &
                    &2, 3) / dy - gDEdgeB * pDEdgeBGra / dy + gUEdgeB &
                    &* pUEdgeBGra * 0.25 * met(i, j - 1, k - 1, 2, 3) / dz
              else if(k == nz .and. zBoundary == "solid_wall") then
                ABD = - gEdgeB * pEdgeBGra * met(i, j - 1, k - 1, 2, 3) / dz &
                    &- gEdgeD * pEdgeDGra * 0.25 * met(i, j - 1, k - 1, 2, 3) &
                    &/ dy - gDEdgeB * pDEdgeBGra / dy
              else
                ABD = - gEdgeB * pEdgeBGra * 0.25 * met(i, j - 1, k - 1, 2, 3) &
                    &/ dz - gEdgeD * pEdgeDGra * 0.25 * met(i, j - 1, k - 1, &
                    &2, 3) / dy - gDEdgeB * pDEdgeBGra / dy
              end if

              ! ------------------ A(i,j,k+2) -------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                AUU = - gEdgeR * pEdgeRGra * 0.25 * met(i, j, k + 2, 1, 3) &
                    &/ dz - gEdgeL * pEdgeLGra * 0.25 * met(i, j, k + 2, 1, 3) &
                    &/ dz - gEdgeF * pEdgeFGra * 0.25 * met(i, j, k + 2, 2, 3) &
                    &/ dz - gEdgeB * pEdgeBGra * 0.25 * met(i, j, k + 2, 2, 3) &
                    &/ dz + gUEdgeR * pUEdgeRGra * 0.25 * met(i, j, k + 2, 1, &
                    &3) / dz + gUEdgeL * pUEdgeLGra * 0.25 * met(i, j, k + 2, &
                    &1, 3) / dz + gUEdgeF * pUEdgeFGra * 0.25 * met(i, j, k &
                    &+ 2, 2, 3) / dz + gUEdgeB * pUEdgeBGra * 0.25 * met(i, j, &
                    &k + 2, 2, 3) / dz
              else if((k == nz - 1 .or. k == nz) .and. zBoundary &
                  &== "solid_wall") then
                AUU = 0.0
              else
                AUU = gUEdgeR * pUEdgeRGra * 0.25 * met(i, j, k + 2, 1, 3) &
                    &/ dz + gUEdgeL * pUEdgeLGra * 0.25 * met(i, j, k + 2, 1, &
                    &3) / dz + gUEdgeF * pUEdgeFGra * 0.25 * met(i, j, k + 2, &
                    &2, 3) / dz + gUEdgeB * pUEdgeBGra * 0.25 * met(i, j, k &
                    &+ 2, 2, 3) / dz
              end if

              ! ------------------ A(i,j,k-2) -------------------!

              if((k == 1 .or. k == 2) .and. zBoundary == "solid_wall") then
                ADD = 0.0
              else if(k == nz .and. zBoundary == "solid_wall") then
                ADD = gEdgeR * pEdgeRGra * 0.25 * met(i, j, k - 2, 1, 3) / dz &
                    &+ gEdgeL * pEdgeLGra * 0.25 * met(i, j, k - 2, 1, 3) / dz &
                    &+ gEdgeF * pEdgeFGra * 0.25 * met(i, j, k - 2, 2, 3) / dz &
                    &+ gEdgeB * pEdgeBGra * 0.25 * met(i, j, k - 2, 2, 3) / dz &
                    &- gDEdgeR * pDEdgeRGra * 0.25 * met(i, j, k - 2, 1, 3) &
                    &/ dz - gDEdgeL * pDEdgeLGra * 0.25 * met(i, j, k - 2, 1, &
                    &3) / dz - gDEdgeF * pDEdgeFGra * 0.25 * met(i, j, k - 2, &
                    &2, 3) / dz - gDEdgeB * pDEdgeBGra * 0.25 * met(i, j, k &
                    &- 2, 2, 3) / dz
              else
                ADD = - gDEdgeR * pDEdgeRGra * 0.25 * met(i, j, k - 2, 1, 3) &
                    &/ dz - gDEdgeL * pDEdgeLGra * 0.25 * met(i, j, k - 2, 1, &
                    &3) / dz - gDEdgeF * pDEdgeFGra * 0.25 * met(i, j, k - 2, &
                    &2, 3) / dz - gDEdgeB * pDEdgeBGra * 0.25 * met(i, j, k &
                    &- 2, 2, 3) / dz
              end if

              ! ----------------- A(i+1,j,k+2) ------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                ARUU = - gEdgeR * pEdgeRGra * 0.25 * met(i + 1, j, k + 2, 1, &
                    &3) / dz + gUEdgeR * pUEdgeRGra * 0.25 * met(i + 1, j, k &
                    &+ 2, 1, 3) / dz
              else if((k == nz - 1 .or. k == nz) .and. zBoundary &
                  &== "solid_wall") then
                ARUU = 0.0
              else
                ARUU = gUEdgeR * pUEdgeRGra * 0.25 * met(i + 1, j, k + 2, 1, &
                    &3) / dz
              end if

              ! ----------------- A(i+1,j,k-2) ------------------!

              if((k == 1 .or. k == 2) .and. zBoundary == "solid_wall") then
                ARDD = 0.0
              else if(k == nz .and. zBoundary == "solid_wall") then
                ARDD = gEdgeR * pEdgeRGra * 0.25 * met(i + 1, j, k - 2, 1, 3) &
                    &/ dz - gDEdgeR * pDEdgeRGra * 0.25 * met(i + 1, j, k - 2, &
                    &1, 3) / dz
              else
                ARDD = - gDEdgeR * pDEdgeRGra * 0.25 * met(i + 1, j, k - 2, 1, &
                    &3) / dz
              end if

              ! ----------------- A(i-1,j,k+2) ------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                ALUU = - gEdgeL * pEdgeLGra * 0.25 * met(i - 1, j, k + 2, 1, &
                    &3) / dz + gUEdgeL * pUEdgeLGra * 0.25 * met(i - 1, j, k &
                    &+ 2, 1, 3) / dz
              else if((k == nz - 1 .or. k == nz) .and. zBoundary &
                  &== "solid_wall") then
                ALUU = 0.0
              else
                ALUU = gUEdgeL * pUEdgeLGra * 0.25 * met(i - 1, j, k + 2, 1, &
                    &3) / dz
              end if

              ! ----------------- A(i-1,j,k-2) ------------------!

              if((k == 1 .or. k == 2) .and. zBoundary == "solid_wall") then
                ALDD = 0.0
              else if(k == nz .and. zBoundary == "solid_wall") then
                ALDD = gEdgeL * pEdgeLGra * 0.25 * met(i - 1, j, k - 2, 1, 3) &
                    &/ dz - gDEdgeL * pDEdgeLGra * 0.25 * met(i - 1, j, k - 2, &
                    &1, 3) / dz
              else
                ALDD = - gDEdgeL * pDEdgeLGra * 0.25 * met(i - 1, j, k - 2, 1, &
                    &3) / dz
              end if

              ! ----------------- A(i,j+1,k+2) ------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                AFUU = - gEdgeF * pEdgeFGra * 0.25 * met(i, j + 1, k + 2, 2, &
                    &3) / dz + gUEdgeF * pUEdgeFGra * 0.25 * met(i, j + 1, k &
                    &+ 2, 2, 3) / dz
              else if((k == nz - 1 .or. k == nz) .and. zBoundary &
                  &== "solid_wall") then
                AFUU = 0.0
              else
                AFUU = gUEdgeF * pUEdgeFGra * 0.25 * met(i, j + 1, k + 2, 2, &
                    &3) / dz
              end if

              ! ----------------- A(i,j+1,k-2) ------------------!

              if((k == 1 .or. k == 2) .and. zBoundary == "solid_wall") then
                AFDD = 0.0
              else if(k == nz .and. zBoundary == "solid_wall") then
                AFDD = gEdgeF * pEdgeFGra * 0.25 * met(i, j + 1, k - 2, 2, 3) &
                    &/ dz - gDEdgeF * pDEdgeFGra * 0.25 * met(i, j + 1, k - 2, &
                    &2, 3) / dz
              else
                AFDD = - gDEdgeF * pDEdgeFGra * 0.25 * met(i, j + 1, k - 2, 2, &
                    &3) / dz
              end if

              ! ----------------- A(i,j-1,k+2) ------------------!

              if(k == 1 .and. zBoundary == "solid_wall") then
                ABUU = - gEdgeB * pEdgeBGra * 0.25 * met(i, j - 1, k + 2, 2, &
                    &3) / dz + gUEdgeB * pUEdgeBGra * 0.25 * met(i, j - 1, k &
                    &+ 2, 2, 3) / dz
              else if((k == nz - 1 .or. k == nz) .and. zBoundary &
                  &== "solid_wall") then
                ABUU = 0.0
              else
                ABUU = gUEdgeB * pUEdgeBGra * 0.25 * met(i, j - 1, k + 2, 2, &
                    &3) / dz
              end if

              ! ----------------- A(i,j-1,k-2) ------------------!

              if((k == 1 .or. k == 2) .and. zBoundary == "solid_wall") then
                ABDD = 0.0
              else if(k == nz .and. zBoundary == "solid_wall") then
                ABDD = gEdgeB * pEdgeBGra * 0.25 * met(i, j - 1, k - 2, 2, 3) &
                    &/ dz - gDEdgeB * pDEdgeBGra * 0.25 * met(i, j - 1, k - 2, &
                    &2, 3) / dz
              else
                ABDD = - gDEdgeB * pDEdgeBGra * 0.25 * met(i, j - 1, k - 2, 2, &
                    &3) / dz
              end if

              ! Scale the matrix elements.
              AC = AC / (fcscal ** 2.0)
              AR = AR / (fcscal ** 2.0)
              AL = AL / (fcscal ** 2.0)
              AF = AF / (fcscal ** 2.0)
              AB = AB / (fcscal ** 2.0)
              AU = AU / fcscal / fcscal_u
              AD = AD / fcscal / fcscal_d
              ARU = ARU / fcscal / fcscal_u
              ARD = ARD / fcscal / fcscal_d
              ALU = ALU / fcscal / fcscal_u
              ALD = ALD / fcscal / fcscal_d
              AFU = AFU / fcscal / fcscal_u
              AFD = AFD / fcscal / fcscal_d
              ABU = ABU / fcscal / fcscal_u
              ABD = ABD / fcscal / fcscal_d
              AUU = AUU / fcscal / fcscal_uu
              ADD = ADD / fcscal / fcscal_dd
              ARUU = ARUU / fcscal / fcscal_uu
              ARDD = ARDD / fcscal / fcscal_dd
              ALUU = ALUU / fcscal / fcscal_uu
              ALDD = ALDD / fcscal / fcscal_dd
              AFUU = AFUU / fcscal / fcscal_uu
              AFDD = AFDD / fcscal / fcscal_dd
              ABUU = ABUU / fcscal / fcscal_uu
              ABDD = ABDD / fcscal / fcscal_dd

              ! Set matrix elements for bicgstab.
              ac_b(i, j, k) = AC
              ar_b(i, j, k) = AR
              al_b(i, j, k) = AL
              af_b(i, j, k) = AF
              ab_b(i, j, k) = AB
              au_b(i, j, k) = AU
              ad_b(i, j, k) = AD
              aru_b(i, j, k) = ARU
              ard_b(i, j, k) = ARD
              alu_b(i, j, k) = ALU
              ald_b(i, j, k) = ALD
              afu_b(i, j, k) = AFU
              afd_b(i, j, k) = AFD
              abu_b(i, j, k) = ABU
              abd_b(i, j, k) = ABD
              auu_b(i, j, k) = AUU
              add_b(i, j, k) = ADD
              aruu_b(i, j, k) = ARUU
              ardd_b(i, j, k) = ARDD
              aluu_b(i, j, k) = ALUU
              aldd_b(i, j, k) = ALDD
              afuu_b(i, j, k) = AFUU
              afdd_b(i, j, k) = AFDD
              abuu_b(i, j, k) = ABUU
              abdd_b(i, j, k) = ABDD

              ! Store horizontal and vertical components of AC (for
              ! preconditioner).
              if(preconditioner == "yes") then
                ach_b(i, j, k) = - AR - AL - AF - AB
                acv_b(i, j, k) = - AU - AD
              end if
            end do
          end do
        end do
        kr_sp_tfc = kr_sp_tfc / facray
        kr_sp_w_tfc = kr_sp_w_tfc / facray
      else
        do k = 1, nz
          fcscal = sqrt(Pstrat(k) ** 2 / rhoStrat(k))
          fcscal_u = sqrt(Pstrat(k + 1) ** 2 / rhoStrat(k + 1))
          fcscal_d = sqrt(Pstrat(k - 1) ** 2 / rhoStrat(k - 1))
          do j = 1, ny
            do i = 1, nx
              AL = 0.0
              AR = 0.0
              AB = 0.0
              AF = 0.0
              AD = 0.0
              AU = 0.0
              AC = 0.0

              ACH = 0.0
              ACV = 0.0

              ALB = 0.0
              ALF = 0.0
              ARB = 0.0
              ARF = 0.0

              ! ------------------- from P UR/dx ------------------------

              facu = 1.0

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  facu = facu + dt * kv_hs(j, k)
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                facu = facu + dt * kr_sp(j, k)
              end if

              facv = facu

              ! A(i+1,j,k) and A(i,j,k)

              rhoEdge = 0.5 * (var%rho(i + 1, j, k) + var%rho(i, j, k)) &
                  &+ rhoStrat(k)

              acontr = dx2 * pStrat(k) ** 2 / rhoEdge * facv / (facu * facv &
                  &+ (f_cor_nd(j) * dt) ** 2)

              AR = AR + acontr
              AC = AC - acontr

              ACH = ACH - acontr

              ! A(i,j,k) and A(i,j-1,k)

              rhoEdge = 0.5 * (var%rho(i, j - 1, k) + var%rho(i, j, k)) &
                  &+ rhoStrat(k)

              acontr = 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_nd(j) &
                  &* dt / (facu * facv + (f_cor_nd(j) * dt) ** 2)

              ACH = ACH + acontr

              AC = AC + acontr
              AB = AB - acontr

              ! A(i,j,k) and A(i,j+1,k)

              rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k))
              rhoEdge = rhoEdge + rhoStrat(k)

              acontr = 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_nd(j) &
                  &* dt / (facu * facv + (f_cor_nd(j) * dt) ** 2)

              AF = AF + acontr
              AC = AC - acontr

              ACH = ACH - acontr

              ! A(i+1,j,k) and A(i+1,j-1,k)

              rhoEdge = 0.5 * (var%rho(i + 1, j - 1, k) + var%rho(i + 1, j, &
                  &k)) + rhoStrat(k)

              acontr = 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_nd(j) &
                  &* dt / (facu * facv + (f_cor_nd(j) * dt) ** 2)

              AR = AR + acontr
              ARB = ARB - acontr

              ! A(i+1,j,k) and A(i+1,j+1,k)

              rhoEdge = 0.5 * (var%rho(i + 1, j, k) + var%rho(i + 1, j + 1, &
                  &k)) + rhoStrat(k)

              acontr = 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_nd(j) &
                  &* dt / (facu * facv + (f_cor_nd(j) * dt) ** 2)

              ARF = ARF + acontr
              AR = AR - acontr

              ! ------------------- from - P UL/dx --------------------

              facu = 1.0

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  facu = facu + dt * kv_hs(j, k)
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                facu = facu + dt * kr_sp(j, k)
              end if

              facv = facu

              ! A(i,j,k) and A(i-1,j,k)

              rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i - 1, j, k)) &
                  &+ rhoStrat(k)

              acontr = - dx2 * pStrat(k) ** 2 / rhoEdge * facv / (facu * facv &
                  &+ (f_cor_nd(j) * dt) ** 2)

              ACH = ACH + acontr

              AC = AC + acontr
              AL = AL - acontr

              ! A(i-1,j,k) and A(i-1,j-1,k)

              rhoEdge = 0.5 * (var%rho(i - 1, j - 1, k) + var%rho(i - 1, j, &
                  &k)) + rhoStrat(k)

              acontr = - 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_nd(j) &
                  &* dt / (facu * facv + (f_cor_nd(j) * dt) ** 2)

              AL = AL + acontr
              ALB = ALB - acontr

              ! A(i-1,j,k) and A(i-1,j+1,k)

              rhoEdge = 0.5 * (var%rho(i - 1, j, k) + var%rho(i - 1, j + 1, &
                  &k)) + rhoStrat(k)

              acontr = - 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_nd(j) &
                  &* dt / (facu * facv + (f_cor_nd(j) * dt) ** 2)

              ALF = ALF + acontr
              AL = AL - acontr

              ! A(i,j,k) and A(i,j-1,k)

              rhoEdge = 0.5 * (var%rho(i, j - 1, k) + var%rho(i, j, k)) &
                  &+ rhoStrat(k)

              acontr = - 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_nd(j) &
                  &* dt / (facu * facv + (f_cor_nd(j) * dt) ** 2)

              ACH = ACH + acontr

              AC = AC + acontr
              AB = AB - acontr

              ! A(i,j,k) and A(i,j+1,k)

              rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k)) &
                  &+ rhoStrat(k)

              acontr = - 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_nd(j) &
                  &* dt / (facu * facv + (f_cor_nd(j) * dt) ** 2)

              AF = AF + acontr
              AC = AC - acontr

              ACH = ACH - acontr

              ! ------------------- from P VF/dy ------------------------

              facv = 1.0

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  facv = facv + dt * 0.5 * (kv_hs(j, k) + kv_hs(j + 1, k))
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                facv = facv + dt * 0.5 * (kr_sp(j, k) + kr_sp(j + 1, k))
              end if

              ! Coriolis parameter interpolated to v point
              f_cor_v = 0.5 * (f_cor_nd(j) + f_cor_nd(j + 1))

              facu = facv

              ! A(i+1,j,k) and A(i,j,k)

              rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k)) &
                  &+ rhoStrat(k)

              acontr = - 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_v * dt &
                  &/ (facu * facv + (f_cor_v * dt) ** 2)

              AR = AR + acontr
              AC = AC - acontr

              ACH = ACH - acontr

              ! A(i,j,k) and A(i-1,j,k)

              rhoEdge = 0.5 * (var%rho(i - 1, j, k) + var%rho(i, j, k)) &
                  &+ rhoStrat(k)

              acontr = - 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_v * dt &
                  &/ (facu * facv + (f_cor_v * dt) ** 2)

              ACH = ACH + acontr

              AC = AC + acontr
              AL = AL - acontr

              ! A(i+1,j+1,k) and A(i,j+1,k)

              rhoEdge = 0.5 * (var%rho(i, j + 1, k) + var%rho(i + 1, j + 1, &
                  &k)) + rhoStrat(k)

              acontr = - 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_v * dt &
                  &/ (facu * facv + (f_cor_v * dt) ** 2)

              ARF = ARF + acontr
              AF = AF - acontr

              ! A(i,j+1,k) and A(i-1,j+1,k)

              rhoEdge = 0.5 * (var%rho(i - 1, j + 1, k) + var%rho(i, j + 1, &
                  &k)) + rhoStrat(k)

              acontr = - 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_v * dt &
                  &/ (facu * facv + (f_cor_v * dt) ** 2)

              AF = AF + acontr
              ALF = ALF - acontr

              ! A(i,j+1,k) and A(i,j,k)

              rhoEdge = 0.5 * (var%rho(i, j + 1, k) + var%rho(i, j, k)) &
                  &+ rhoStrat(k)

              acontr = dy2 * pStrat(k) ** 2 / rhoEdge * facu / (facu * facv &
                  &+ (f_cor_v * dt) ** 2)

              AF = AF + acontr
              AC = AC - acontr

              ACH = ACH - acontr

              ! ------------------- from - P VB/dy ---------------------

              facv = 1.0

              if(TestCase == "baroclinic_LC") then
                if(background == "HeldSuarez") then
                  ! Rayleigh damping
                  facv = facv + dt * 0.5 * (kv_hs(j, k) + kv_hs(j - 1, k))
                end if
              end if

              if(spongeLayer .and. sponge_uv) then
                facv = facv + dt * 0.5 * (kr_sp(j, k) + kr_sp(j - 1, k))
              end if

              ! Coriolis parameter interpolated to v point
              f_cor_v = 0.5 * (f_cor_nd(j) + f_cor_nd(j - 1))

              facu = facv

              ! A(i+1,j-1,k) and A(i,j-1,k)

              rhoEdge = 0.5 * (var%rho(i, j - 1, k) + var%rho(i + 1, j - 1, &
                  &k)) + rhoStrat(k)

              acontr = 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_v * dt &
                  &/ (facu * facv + (f_cor_v * dt) ** 2)

              ARB = ARB + acontr
              AB = AB - acontr

              ! A(i,j-1,k) and A(i-1,j-1,k)

              rhoEdge = 0.5 * (var%rho(i - 1, j - 1, k) + var%rho(i, j - 1, &
                  &k)) + rhoStrat(k)

              acontr = 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_v * dt &
                  &/ (facu * facv + (f_cor_v * dt) ** 2)

              AB = AB + acontr
              ALB = ALB - acontr

              ! A(i+1,j,k) and A(i,j,k)

              rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k)) &
                  &+ rhoStrat(k)

              acontr = 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_v * dt &
                  &/ (facu * facv + (f_cor_v * dt) ** 2)

              AR = AR + acontr
              AC = AC - acontr

              ACH = ACH - acontr

              ! A(i,j,k) and A(i-1,j,k)

              rhoEdge = 0.5 * (var%rho(i - 1, j, k) + var%rho(i, j, k)) &
                  &+ rhoStrat(k)

              acontr = 0.25 * dxy * pStrat(k) ** 2 / rhoEdge * f_cor_v * dt &
                  &/ (facu * facv + (f_cor_v * dt) ** 2)

              ACH = ACH + acontr

              AC = AC + acontr
              AL = AL - acontr

              ! A(i,j,k) and A(i,j-1,k)

              rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i, j - 1, k)) &
                  &+ rhoStrat(k)

              acontr = - dy2 * pStrat(k) ** 2 / rhoEdge * facu / (facu * facv &
                  &+ (f_cor_v * dt) ** 2)

              ACH = ACH + acontr

              AC = AC + acontr
              AB = AB - acontr

              ! ------------------- from PU WU/dz ---------------------

              ! TFC
              if(k == nz .and. zBoundary == "solid_wall") then
                AU = 0.0
              else
                facw = 1.0

                if(TestCase == "baroclinic_LC") then
                  if(background == "HeldSuarez") then
                    ! Rayleigh damping

                    facw = facw + dt * 0.5 * (kw_hs(k) + kw_hs(k + 1))
                  end if
                end if

                if(spongeLayer) then
                  facw = facw + dt * 0.5 * (kr_sp(j, k) + kr_sp(j, k + 1))
                end if

                bvsstw = 0.5 * (bvsStrat(k) + bvsStrat(k + 1))

                ! A(i,j,k+1) and A(i,j,k)

                rhoEdge = 0.5 * (var%rho(i, j, k + 1) + var%rho(i, j, k)) &
                    &+ rhoStratTilde(k)

                pStratU = 0.5 * (pStrat(k + 1) + pStrat(k))
                pStratU_0 = 0.5 * (pStrat_0(k + 1) + pStrat_0(k))

                AU = dz2 * pStratU ** 2 / rhoEdge / (facw + rhoStratTilde(k) &
                    &/ rhoEdge * pStratU / pStratU_0 * bvsstw * dt ** 2)

                AC = AC - AU

                ACV = ACV - AU
              end if

              ! ------------------- from - PD WD/dz ---------------------

              ! TFC
              if(k == 1 .and. zBoundary == "solid_wall") then
                AD = 0.0
              else
                facw = 1.0

                if(TestCase == "baroclinic_LC") then
                  if(background == "HeldSuarez") then
                    ! Rayleigh damping

                    facw = facw + dt * 0.5 * (kw_hs(k) + kw_hs(k - 1))
                  end if
                end if

                if(spongeLayer) then
                  facw = facw + dt * 0.5 * (kr_sp(j, k) + kr_sp(j, k - 1))
                end if

                bvsstw = 0.5 * (bvsStrat(k - 1) + bvsStrat(k))

                ! A(i,j,k) and A(i,j,k-1)

                rhoEdge = 0.5 * (var%rho(i, j, k) + var%rho(i, j, k - 1)) &
                    &+ rhoStratTilde(k - 1)

                pStratD = 0.5 * (pStrat(k) + pStrat(k - 1))
                pStratD_0 = 0.5 * (pStrat_0(k) + pStrat_0(k - 1))

                AD = dz2 * pStratD ** 2 / rhoEdge / (facw + rhoStratTilde(k &
                    &- 1) / rhoEdge * pStratD / pStratD_0 * bvsstw * dt ** 2)

                AC = AC - AD

                ACV = ACV - AD
              end if

              AC = AC / fcscal ** 2

              ACH = ACH / fcscal ** 2
              ACV = ACV / fcscal ** 2

              AL = AL / fcscal ** 2
              AR = AR / fcscal ** 2
              AB = AB / fcscal ** 2
              AF = AF / fcscal ** 2
              AD = AD / (fcscal * fcscal_d)
              AU = AU / (fcscal * fcscal_u)
              ALB = ALB / fcscal ** 2
              ALF = ALF / fcscal ** 2
              ARF = ARF / fcscal ** 2
              ARB = ARB / fcscal ** 2

              if(pressureScaling) then
                AC = AC / Pstrat(k)

                ACH = ACH / Pstrat(k)
                ACV = ACV / Pstrat(k)

                AL = AL / Pstrat(k)
                AR = AR / Pstrat(k)
                AB = AB / Pstrat(k)
                AF = AF / Pstrat(k)
                AD = AD / PstratTilde(k - 1)
                AU = AU / PstratTilde(k)
                ALB = ALB / Pstrat(k)
                ALF = ALF / Pstrat(k)
                ARF = ARF / Pstrat(k)
                ARB = ARB / Pstrat(k)
              end if

              ! ------------------- define matrix A -------------------

              if(poissonSolverType == 'bicgstab') then
                ac_b(i, j, k) = AC

                ach_b(i, j, k) = ACH
                acv_b(i, j, k) = ACV

                al_b(i, j, k) = AL
                ar_b(i, j, k) = AR

                ab_b(i, j, k) = AB
                af_b(i, j, k) = AF

                ad_b(i, j, k) = AD
                au_b(i, j, k) = AU

                alb_b(i, j, k) = ALB
                alf_b(i, j, k) = ALF
                arb_b(i, j, k) = ARB
                arf_b(i, j, k) = ARF
              else
                stop 'ERROR: val_PsIn expects bicgstab'
              end if
            end do ! i_loop
          end do ! j_loop
        end do ! k_loop
      end if

      kr_sp = kr_sp / facray
      alprlx = alprlx / facray

    else
      stop 'ERROR: wrong opt'
    end if

    return

  end subroutine val_PsIn

  !------------------------------------------------------------------

  subroutine heat_w0(var, flux, dt, heat, S_bar, w_0)

    ! negative (!) heating, its horizontal mean,
    ! and the thereby induced vertical wind

    ! in/out variables
    type(var_type), intent(in) :: var
    type(flux_type), intent(in) :: flux
    real, intent(in) :: dt
    real, dimension(- nbx:nx + nbx, - nby:ny + nby, - nbz:nz + nbz), &
        &intent(out) :: heat
    real, dimension(- nbz:nz + nbz), intent(out) :: w_0
    real, dimension(- nbz:nz + nbz), intent(out) :: S_bar

    character(len = 40) :: w0_mod

    integer :: i, j, k
    real, dimension(1:nz) :: sum_local, sum_global
    real, dimension(- nbz:nz + nbz) :: press0
    real, dimension(- nbz:nz + nbz) :: rhow_bar

    real :: dptopdt
    real :: expo
    real :: rho, wvert

    real :: sum_d, sum_n

    w0_mod = 'ONK14'

    w_0 = 0.
    S_bar = 0.
    heat = 0.

    ! No heating in TFC
    if(topography) return

    ! negative (!) heating, i.e. -S eq(9)  ONeill+Klein2014
    call calculate_heating(var, flux, heat)

    ! calculate horizontal mean of heat(:,:,:)
    do k = 1, nz
      sum_local(k) = sum(heat(1:nx, 1:ny, k))
    end do
    !global sum and average
    call mpi_allreduce(sum_local(1), sum_global(1), nz - 1 + 1, &
        &mpi_double_precision, mpi_sum, comm, ierror)
    sum_global = sum_global / (sizeX * sizeY)

    S_bar(1:nz) = sum_global(1:nz)

    if(w0_mod == 'Almgrenetal08') then
      ! horizontal mean of the vertical density flux

      rhow_bar = 0.

      sum_local = 0.
      sum_global = 0.

      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            rho = rhoStrat(k) + var%rho(i, j, k)

            if(k == 1) then
              wvert = 0.5 * var%w(i, j, k)
            else if(k == nz) then
              wvert = 0.5 * var%w(i, j, k - 1)
            else
              wvert = 0.5 * (var%w(i, j, k - 1) + var%w(i, j, k))
            end if

            sum_local(k) = sum_local(k) + rho * wvert
          end do
        end do
      end do
      call mpi_allreduce(sum_local(1), sum_global(1), nz - 1 + 1, &
          &mpi_double_precision, mpi_sum, comm, ierror)
      sum_global = sum_global / (sizeX * sizeY)

      rhow_bar(1:nz) = sum_global(1:nz)

    end if

    !non_dim. pressure; eq(12) ONeill+Klein2014

    do k = 1, nz
      press0(k) = PStrat(k) ** gamma
    end do

    !  time derivative of reference pressure at the model top

    sum_d = 0.0
    sum_n = 0.0

    if(w0_mod == 'Almgrenetal08') then
      do k = 1, nz
        expo = exp(- g_ndim * rhoStrat(k) / (gamma * press0(k)) * z(k))

        sum_n = sum_n + expo * (- S_bar(k) / PStrat(k) - g_ndim * rhow_bar(k) &
            &/ (gamma * press0(k)))

        sum_d = sum_d + expo / (gamma * press0(k))
      end do
    else if(w0_mod == 'ONK14') then
      do k = 1, nz
        sum_n = sum_n - S_bar(k) / PStrat(k)
        sum_d = sum_d + 1. / (gamma * press0(k))
      end do
    else
      stop 'ERROR: wrong w0_mod'
    end if

    dptopdt = sum_n / sum_d

    ! horizontal-mean vertical wind

    w_0 = 0.

    if(w0_mod == 'Almgrenetal08') then
      do k = 1, nz - 1
        expo = exp(- g_ndim * rhoStrat(k) / (gamma * press0(k)) * z(k))

        w_0(k) = w_0(k - 1) + dz * expo * (- S_bar(k) / Pstrat(k) - g_ndim &
            &* rhow_bar(k) / (gamma * press0(k)) - dptopdt / (gamma &
            &* press0(k)))
      end do

      do k = 1, nz - 1
        expo = exp(g_ndim * rhoStrat(k) / (gamma * press0(k)) * 0.5 * (z(k) &
            &+ z(k + 1)))

        w_0(k) = expo * w_0(k)
      end do
    else if(w0_mod == 'ONK14') then
      w_0(1) = dz * (- S_bar(1) / Pstrat(1) - (1. / (gamma * press0(1))) &
          &* dptopdt)

      do k = 2, nz - 1
        w_0(k) = w_0(k - 1) + dz * (- S_bar(k) / Pstrat(k) - (1. / (gamma &
            &* press0(k))) * dptopdt)
      end do
    else
      stop 'ERROR: wrong w0_mod'
    end if

  end subroutine heat_w0

  !------------------------------------------------------------------
  subroutine correctorStepTestTFC(var, dMom, int_mod)

    type(var_type), intent(inout) :: var
    real, dimension((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        &+ nbz), 3), intent(inout) :: dMom
    character(len = *), intent(in) :: int_mod

    type(var_type) :: var_tfc

    call random_number(dp)

    var_tfc = var

    if(int_mod == "expl") then
      call correctorStep(var_tfc, dMom, 1.0, 1, "expl", 1.0, 1.0)
      topography = .false.
      call correctorStep(var, dMom, 1.0, 1, "expl", 1.0, 1.0)
      topography = .true.
      print *, "Zonal-wind corrector-step difference: ", maxval(abs(var_tfc%u &
          &- var%u))
      print *, "Meridional-wind corrector-step difference: ", &
          &maxval(abs(var_tfc%v - var%v))
      print *, "Vertical-wind corrector-step difference: ", &
          &maxval(abs(var_tfc%w - var%w))
      print *, "Density-fluctuation corrector-step difference: ", &
          &maxval(abs(var_tfc%rhop - var%rhop))
    else if(int_mod == "impl") then
      call correctorStep(var_tfc, dMom, 1.0, 1, "impl", 1.0, 1.0)
      topography = .false.
      call correctorStep(var, dMom, 1.0, 1, "impl", 1.0, 1.0)
      topography = .true.
      print *, "Zonal-wind corrector-step difference (impl): ", &
          &maxval(abs(var_tfc%u - var%u))
      print *, "Meridional-wind corrector-step difference (impl): ", &
          &maxval(abs(var_tfc%v - var%v))
      print *, "Vertical-wind corrector-step difference (impl): ", &
          &maxval(abs(var_tfc%w - var%w))
      print *, "Density-fluctuation corrector-step difference (impl): ", &
          &maxval(abs(var_tfc%rhop - var%rhop))

      var_tfc = var
      call correctorStep(var_tfc, dMom, 1.0, 1, "expl", 1.0, 1.0)
      topography = .false.
      call correctorStep(var, dMom, 1.0, 1, "expl", 1.0, 1.0)
      topography = .true.
      print *, "Zonal-wind corrector-step difference (expl): ", &
          &maxval(abs(var_tfc%u - var%u))
      print *, "Meridional-wind corrector-step difference (expl): ", &
          &maxval(abs(var_tfc%v - var%v))
      print *, "Vertical-wind corrector-step difference (expl): ", &
          &maxval(abs(var_tfc%w - var%w))
      print *, "Density-fluctuation corrector-step difference (expl): ", &
          &maxval(abs(var_tfc%rhop - var%rhop))
    end if

  end subroutine correctorStepTestTFC

end module poisson_module
