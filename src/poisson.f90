module poisson_module

  use type_module
  use mpi_module
  use timeScheme_module
  use atmosphere_module
  use bicgstab_tools_module
  use mpi

  implicit none

  private

  !------------------------
  !   public subroutines
  !------------------------

  public :: Corrector
  public :: init_poisson
  public :: terminate_poisson

  !------------------------
  !   private subroutines
  !------------------------

  private :: pressureBoundaryCondition
  private :: correctorStep
  private :: linOpr
  private :: bicgstab
  private :: poissonSolver

  !-------------------------------
  !    private module variables
  !------------------------------

  ! pressure correction
  real, dimension(:, :, :), allocatable :: dp

  real :: tol

  contains

  subroutine Corrector(var, flux, dt, errFlagBicg, nIter, opt, facray, facprs)

    ! -------------------------------------------------
    !              Correct uStar, bStar, and p
    ! -------------------------------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var
    type(flux_type), intent(in) :: flux

    real, intent(in) :: dt, facray, facprs
    logical, intent(out) :: errFlagBicg
    integer, intent(out) :: nIter

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

    ! calculate RHS
    call calc_RHS(rhs, var, flux, dt)

    ! calculate dp
    call poissonSolver(rhs, var, dt, errFlagBicg, nIter, opt, facray, facprs)

    if(errFlagBicg) return

    ! set horizontal and vertical BC for dp
    call pressureBoundaryCondition

    ! correct p, rhopStar, and uStar with dp
    call correctorStep(var, dt, opt, facray, facprs)

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
    real :: AL, AR, AB, AF, AD, AU, AC, ACH, ACV
    real :: sL, sR, sB, sF, sD, sU, sC

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
    case("pseudo_incompressible")

      !----------------------------
      !   set Halo cells: xSlice
      !----------------------------

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

      !------------------------------
      !   set Halo cells: ySlice
      !------------------------------

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

            if(k < nz - 1) then
              AUU = auu_b(i, j, k)
              sUU = s(i, j, k + 2)
            else
              AUU = 0.0
              sUU = 0.0
            end if

            ! ------------------ A(i,j,k-2) -----------------

            if(k > 2) then
              ADD = add_b(i, j, k)
              sDD = s(i, j, k - 2)
            else
              ADD = 0.0
              sDD = 0.0
            end if

            ! ----------------- A(i+1,j,k+2) -----------------

            if(k < nz - 1) then
              ARUU = aruu_b(i, j, k)
              sRUU = s(i + 1, j, k + 2)
            else
              ARUU = 0.0
              sRUU = 0.0
            end if

            ! ----------------- A(i+1,j,k-2) -----------------

            if(k > 2) then
              ARDD = ardd_b(i, j, k)
              sRDD = s(i + 1, j, k - 2)
            else
              ARDD = 0.0
              sRDD = 0.0
            end if

            ! ----------------- A(i-1,j,k+2) -----------------

            if(k < nz - 1) then
              ALUU = aluu_b(i, j, k)
              sLUU = s(i - 1, j, k + 2)
            else
              ALUU = 0.0
              sLUU = 0.0
            end if

            ! ----------------- A(i-1,j,k-2) -----------------

            if(k > 2) then
              ALDD = aldd_b(i, j, k)
              sLDD = s(i - 1, j, k - 2)
            else
              ALDD = 0.0
              sLDD = 0.0
            end if

            ! ----------------- A(i,j+1,k+2) -----------------

            if(k < nz - 1) then
              AFUU = afuu_b(i, j, k)
              sFUU = s(i, j + 1, k + 2)
            else
              AFUU = 0.0
              sFUU = 0.0
            end if

            ! ----------------- A(i,j+1,k-2) -----------------

            if(k > 2) then
              AFDD = afdd_b(i, j, k)
              sFDD = s(i, j + 1, k - 2)
            else
              AFDD = 0.0
              sFDD = 0.0
            end if

            ! ----------------- A(i,j-1,k+2) -----------------

            if(k < nz - 1) then
              ABUU = abuu_b(i, j, k)
              sBUU = s(i, j - 1, k + 2)
            else
              ABUU = 0.0
              sBUU = 0.0
            end if

            ! ----------------- A(i,j-1,k-2) -----------------

            if(k > 2) then
              ABDD = abdd_b(i, j, k)
              sBDD = s(i, j - 1, k - 2)
            else
              ABDD = 0.0
              sBDD = 0.0
            end if

            ! Update operator.
            Ls(i, j, k) = Ls(i, j, k) + ARU * sRU + ARD * sRD + ALU * sLU &
                &+ ALD * sLD + AFU * sFU + AFD * sFD + ABU * sBU + ABD * sBD &
                &+ AUU * sUU + ADD * sDD + ARUU * sRUU + ARDD * sRDD + ALUU &
                &* sLUU + ALDD * sLDD + AFUU * sFUU + AFDD * sFDD + ABUU &
                &* sBUU + ABDD * sBDD
          end do i_loop
        end do j_loop
      end do k_loop

    case default
      stop "linOpr: unknown case model"
    end select

  end subroutine linOpr

  !----------------------------------------------------------------------------

  subroutine calc_RHS(b, var, flux, dt)

    !----------------------------------------
    !   calculates the RHS of the
    !   Poisson problem
    !----------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var
    type(flux_type), intent(in) :: flux
    real, intent(in) :: dt
    real, dimension(1:nx, 1:ny, 1:nz), intent(out) :: b ! RHS

    ! local vars
    real :: uR, uL, vF, vB, wU, wD
    real, dimension(1:nz) :: sum_local, sum_global

    real :: pEdgeR, pEdgeL, pEdgeF, pEdgeB, pEdgeU, pEdgeD

    integer :: i, j, k
    real :: divSum

    real :: bu, bv, bw, bl2loc, divL2_norm, divL2_norm_local

    ! check L2-norm of divergence
    real :: divL2, divMax

    ! MPI stuff
    real :: divL2_local, divSum_local
    integer :: root

    real :: fcscal

    !--------------------------------------------------
    !    Setup b = Ma^2 * P * u^*  (right hand side)
    !--------------------------------------------------

    divSum = 0.0
    divL2 = 0.0
    divMax = 0.0

    divSum_local = 0.0
    divL2_local = 0.0

    divL2_norm = 0.0
    divL2_norm_local = 0.0

    select case(model)

    case("pseudo_incompressible")

      ! Calculate RHS for TFC.
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            ! Calculate scaling factor.
            fcscal = sqrt(pStratTFC(i, j, k) ** 2.0 / rhoStratTFC(i, j, k))
            ! Store velocities at cell edges.
            uR = var%u(i, j, k)
            uL = var%u(i - 1, j, k)
            vF = var%v(i, j, k)
            vB = var%v(i, j - 1, k)
            wU = var%w(i, j, k)
            wD = var%w(i, j, k - 1)
            ! Calculate P at cell edges.
            pEdgeR = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i + 1, j, &
                &k) * pStratTFC(i + 1, j, k))
            pEdgeL = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i - 1, j, &
                &k) * pStratTFC(i - 1, j, k))
            pEdgeF = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j + 1, &
                &k) * pStratTFC(i, j + 1, k))
            pEdgeB = 0.5 * (jac(i, j, k) * pStratTFC(i, j, k) + jac(i, j - 1, &
                &k) * pStratTFC(i, j - 1, k))
            pEdgeU = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
                &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            pEdgeD = jac(i, j, k) * jac(i, j, k - 1) * (pStratTFC(i, j, k) &
                &+ pStratTFC(i, j, k - 1)) / (jac(i, j, k) + jac(i, j, k - 1))
            ! Compute RHS.
            bu = (pEdgeR * uR - pEdgeL * uL) / dx / jac(i, j, k) * Ma ** 2.0 &
                &* kappa
            bv = (pEdgeF * vF - pEdgeB * vB) / dy / jac(i, j, k) * Ma ** 2.0 &
                &* kappa
            bw = (pEdgeU * wU - pEdgeD * wD) / dz / jac(i, j, k) * Ma ** 2.0 &
                &* kappa
            divSum_local = divSum_local + bu + bv + bw
            bu = bu / fcscal
            bv = bv / fcscal
            bw = bw / fcscal
            b(i, j, k) = bu + bv + bw
            ! Compute check sum for solvability criterion.
            divL2_local = divL2_local + b(i, j, k) ** 2.0
            bl2loc = bu ** 2.0 + bv ** 2.0 + bw ** 2.0
            divL2_norm_local = divL2_norm_local + bl2loc
            if(abs(b(i, j, k)) > divMax) then
              divMax = abs(b(i, j, k))
            end if
          end do
        end do
      end do

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

      divL2_norm_local = sqrt(divL2_norm_local / nx / ny / nz)
      divL2_norm = sqrt(divL2_norm / sizeX / sizeY / sizeZ)

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

    case default
      stop "poissonSolver: unknown case model."
    end select

  end subroutine calc_RHS

  !----------------------------------------------------------------------

  subroutine poissonSolver(b, var, dt, errFlagBicg, nIter, opt, facray, facprs)

    ! -------------------------------------------------
    ! solves the Poisson problem with
    ! application of linear operator L
    ! -------------------------------------------------

    ! in/out variables
    type(var_type), intent(in) :: var
    real, intent(in) :: dt, facray, facprs

    logical, intent(out) :: errFlagBicg
    integer, intent(out) :: nIter

    real, dimension(1:nx, 1:ny, 1:nz), intent(in) :: b

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

    integer :: i, j, k
    real :: fcscal

    ! Init
    if(dt == 0.0) stop "poissonSolver: dt = 0.0. Stopping."
    dtInv = 1.0 / dt

    !--------------------------------
    !     Linear equation solver
    !     solve for dt * dp ...
    !--------------------------------

    sol = 0.0

    select case(model)

    case("pseudo_incompressible")
      call val_PsIn(var, dt, opt, facray)
    case default
      stop "linOpr: unknown case model"
    end select

    call bicgstab(b, dt, sol, res, nIter, errFlagBicg, opt)

    if(errFlagBicg) return

    select case(model)
    case("pseudo_incompressible")
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            fcscal = sqrt(pStratTFC(i, j, k) ** 2 / rhoStratTFC(i, j, k))
            sol(i, j, k) = sol(i, j, k) / fcscal
          end do
        end do
      end do
    case default
    end select

    ! now get dp from dt * dp ...
    ! pass solution to pressure corrector
    dp(1:nx, 1:ny, 1:nz) = dtInv / facprs * sol

  end subroutine poissonSolver

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

    ! Find neighbour procs
    if(idim > 1) call mpi_cart_shift(comm, 0, 1, left, right, ierror)
    if(jdim > 1) call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

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

    !----------------
    !   z-Boundary
    !----------------

    select case(zBoundary)
    case("solid_wall")
      dp(:, :, 0) = dp(:, :, 1)
      dp(:, :, nz + 1) = dp(:, :, nz)
    case default
      stop "pressureBoundaryCondition: unknown case zBoundary."
    end select

  end subroutine pressureBoundaryCondition

  !-----------------------------------------------------------------------

  subroutine correctorStep(var, dt, opt, facray, facprs)

    !-------------------------------------------------
    !         correct pressure & velocity
    !-------------------------------------------------

    ! in/out variables
    type(var_type), intent(inout) :: var
    real, intent(in) :: dt, facray, facprs

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
    real :: rhoEdge, rhou, rhov, rho
    real :: pGradX, pGradY, pGradZ
    real :: du, dv, dw, db
    real :: facu, facv, facw
    real :: bvsstw

    real :: rhow0, rhowm

    real :: rhoStratEdgeU
    real :: pEdgeR, pEdgeF, pEdgeU, pEdgeD
    real :: dpEdgeR, dpUEdgeR, dpUUEdgeR, dpDEdgeR, dpDDEdgeR, dpEdgeF, &
        &dpUEdgeF, dpUUEdgeF, dpDEdgeF, dpDDEdgeF, dpREdgeU, dpLEdgeU, &
        &dpFEdgeU, dpBEdgeU, dpREdgeD, dpLEdgeD, dpFEdgeD, dpBEdgeD
    real :: met13EdgeR, met23EdgeF, met13EdgeU, met23EdgeU, met33EdgeU, &
        &met13EdgeD, met23EdgeD, met33EdgeD
    real :: pGradZEdgeU, pGradZEdgeD
    real, dimension((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz &
        &+ nbz)) :: corX, corY

    integer :: k0, k1

    ! --------------------------------------
    !             calc p + dp
    ! --------------------------------------

    var%pi(0:nx + 1, 0:ny + 1, 0:nz + 1) = var%pi(0:nx + 1, 0:ny + 1, 0:nz &
        &+ 1) + dp(0:nx + 1, 0:ny + 1, 0:nz + 1)

    if(opt == "impl") then
      kr_sp_tfc = kr_sp_tfc * facray
      kr_sp_w_tfc = kr_sp_w_tfc * facray
    end if

    ! --------------------------------------
    !           calc du and u + du
    ! --------------------------------------

    if(opt == "impl") then
      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            facu = 1.0

            if(spongeLayer .and. sponge_uv) then
              facu = facu + dt * 0.5 * (kr_sp_tfc(i, j, k) + kr_sp_tfc(i + 1, &
                  &j, k))
            end if

            facv = facu

            ! Compute values at cell edges.
            rhou = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, j, k))
            pEdgeR = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i + 1, j, k))
            met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))
            ! Compute pressure difference gradient component.
            if(k == 1 .and. zBoundary == "solid_wall") then
              dpUUEdgeR = 0.5 * (dp(i, j, k + 2) + dp(i + 1, j, k + 2))
              dpUEdgeR = 0.5 * (dp(i, j, k + 1) + dp(i + 1, j, k + 1))
              dpEdgeR = 0.5 * (dp(i, j, k) + dp(i + 1, j, k))
              pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((dp(i + 1, j, k) &
                  &- dp(i, j, k)) / dx + met13EdgeR * (- dpUUEdgeR + 4.0 &
                  &* dpUEdgeR - 3.0 * dpEdgeR) * 0.5 / dz)
            else if(k == nz .and. zBoundary == "solid_wall") then
              dpDDEdgeR = 0.5 * (dp(i, j, k - 2) + dp(i + 1, j, k - 2))
              dpDEdgeR = 0.5 * (dp(i, j, k - 1) + dp(i + 1, j, k - 1))
              dpEdgeR = 0.5 * (dp(i, j, k) + dp(i + 1, j, k))
              pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((dp(i + 1, j, k) &
                  &- dp(i, j, k)) / dx + met13EdgeR * (dpDDEdgeR - 4.0 &
                  &* dpDEdgeR + 3.0 * dpEdgeR) * 0.5 / dz)
            else
              dpUEdgeR = 0.5 * (dp(i, j, k + 1) + dp(i + 1, j, k + 1))
              dpDEdgeR = 0.5 * (dp(i, j, k - 1) + dp(i + 1, j, k - 1))
              pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((dp(i + 1, j, k) &
                  &- dp(i, j, k)) / dx + met13EdgeR * (dpUEdgeR - dpDEdgeR) &
                  &* 0.5 / dz)
            end if
            ! Compute velocity correction.
            corX(i, j, k) = facprs * dt / facu * pGradX
            du = - corX(i, j, k)

            var%u(i, j, k) = var%u(i, j, k) + du
          end do
        end do
      end do
    else if(opt == "expl") then
      if(facprs /= 1.) stop 'ERROR: wrong facprs in explicit sub-step'
      do k = 1, nz
        do j = 1, ny
          do i = 0, nx
            ! Compute values at cell edges.
            rhou = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, j, k))
            pEdgeR = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i + 1, j, k))
            met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))
            ! Compute pressure difference gradient component.
            if(k == 1 .and. zBoundary == "solid_wall") then
              dpUUEdgeR = 0.5 * (dp(i, j, k + 2) + dp(i + 1, j, k + 2))
              dpUEdgeR = 0.5 * (dp(i, j, k + 1) + dp(i + 1, j, k + 1))
              dpEdgeR = 0.5 * (dp(i, j, k) + dp(i + 1, j, k))
              pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((dp(i + 1, j, k) &
                  &- dp(i, j, k)) / dx + met13EdgeR * (- dpUUEdgeR + 4.0 &
                  &* dpUEdgeR - 3.0 * dpEdgeR) * 0.5 / dz)
            else if(k == nz .and. zBoundary == "solid_wall") then
              dpDDEdgeR = 0.5 * (dp(i, j, k - 2) + dp(i + 1, j, k - 2))
              dpDEdgeR = 0.5 * (dp(i, j, k - 1) + dp(i + 1, j, k - 1))
              dpEdgeR = 0.5 * (dp(i, j, k) + dp(i + 1, j, k))
              pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((dp(i + 1, j, k) &
                  &- dp(i, j, k)) / dx + met13EdgeR * (dpDDEdgeR - 4.0 &
                  &* dpDEdgeR + 3.0 * dpEdgeR) * 0.5 / dz)
            else
              dpUEdgeR = 0.5 * (dp(i, j, k + 1) + dp(i + 1, j, k + 1))
              dpDEdgeR = 0.5 * (dp(i, j, k - 1) + dp(i + 1, j, k - 1))
              pGradX = kappaInv * MaInv2 / rhou * pEdgeR * ((dp(i + 1, j, k) &
                  &- dp(i, j, k)) / dx + met13EdgeR * (dpUEdgeR - dpDEdgeR) &
                  &* 0.5 / dz)
            end if

            du = - dt * pGradX

            var%u(i, j, k) = var%u(i, j, k) + du
          end do
        end do
      end do
    else
      stop 'ERROR: wrong opt in correctorStep'
    end if

    !--------------------------------------
    !         calc dv and v + dv
    !--------------------------------------

    if(opt == "impl") then
      do k = 1, nz
        do j = 0, ny
          do i = 1, nx
            facv = 1.0

            if(spongeLayer .and. sponge_uv) then
              facv = facv + dt * 0.5 * (kr_sp_tfc(i, j, k) + kr_sp_tfc(i, j &
                  &+ 1, k))
            end if

            facu = facv

            ! Compute values at cell edges.
            rhov = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j + 1, k))
            pEdgeF = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j + 1, k))
            met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3))
            ! Compute pressure difference gradient component.
            if(k == 1 .and. zBoundary == "solid_wall") then
              dpUUEdgeF = 0.5 * (dp(i, j, k + 2) + dp(i, j + 1, k + 2))
              dpUEdgeF = 0.5 * (dp(i, j, k + 1) + dp(i, j + 1, k + 1))
              dpEdgeF = 0.5 * (dp(i, j, k) + dp(i, j + 1, k))
              pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((dp(i, j + 1, k) &
                  &- dp(i, j, k)) / dy + met23EdgeF * (- dpUUEdgeF + 4.0 &
                  &* dpUEdgeF - 3.0 * dpEdgeF) * 0.5 / dz)
            else if(k == nz .and. zBoundary == "solid_wall") then
              dpDDEdgeF = 0.5 * (dp(i, j, k - 2) + dp(i, j + 1, k - 2))
              dpDEdgeF = 0.5 * (dp(i, j, k - 1) + dp(i, j + 1, k - 1))
              dpEdgeF = 0.5 * (dp(i, j, k) + dp(i, j + 1, k))
              pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((dp(i, j + 1, k) &
                  &- dp(i, j, k)) / dy + met23EdgeF * (dpDDEdgeF - 4.0 &
                  &* dpDEdgeF + 3.0 * dpEdgeF) * 0.5 / dz)
            else
              dpUEdgeF = 0.5 * (dp(i, j, k + 1) + dp(i, j + 1, k + 1))
              dpDEdgeF = 0.5 * (dp(i, j, k - 1) + dp(i, j + 1, k - 1))
              pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((dp(i, j + 1, k) &
                  &- dp(i, j, k)) / dy + met23EdgeF * (dpUEdgeF - dpDEdgeF) &
                  &* 0.5 / dz)
            end if
            ! Compute velocity correction.
            corY(i, j, k) = facprs * dt / facv * pGradY
            dv = - corY(i, j, k)

            var%v(i, j, k) = var%v(i, j, k) + dv
          end do
        end do
      end do
    else if(opt == "expl") then
      if(facprs /= 1.) stop 'ERROR: wrong facprs in explicit sub-step'
      do k = 1, nz
        do j = 0, ny
          do i = 1, nx
            ! Compute values at cell edges.
            rhov = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j + 1, k))
            pEdgeF = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j + 1, k))
            met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3))
            ! Compute pressure difference gradient component.
            if(k == 1 .and. zBoundary == "solid_wall") then
              dpUUEdgeF = 0.5 * (dp(i, j, k + 2) + dp(i, j + 1, k + 2))
              dpUEdgeF = 0.5 * (dp(i, j, k + 1) + dp(i, j + 1, k + 1))
              dpEdgeF = 0.5 * (dp(i, j, k) + dp(i, j + 1, k))
              pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((dp(i, j + 1, k) &
                  &- dp(i, j, k)) / dy + met23EdgeF * (- dpUUEdgeF + 4.0 &
                  &* dpUEdgeF - 3.0 * dpEdgeF) * 0.5 / dz)
            else if(k == nz .and. zBoundary == "solid_wall") then
              dpDDEdgeF = 0.5 * (dp(i, j, k - 2) + dp(i, j + 1, k - 2))
              dpDEdgeF = 0.5 * (dp(i, j, k - 1) + dp(i, j + 1, k - 1))
              dpEdgeF = 0.5 * (dp(i, j, k) + dp(i, j + 1, k))
              pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((dp(i, j + 1, k) &
                  &- dp(i, j, k)) / dy + met23EdgeF * (dpDDEdgeF - 4.0 &
                  &* dpDEdgeF + 3.0 * dpEdgeF) * 0.5 / dz)
            else
              dpUEdgeF = 0.5 * (dp(i, j, k + 1) + dp(i, j + 1, k + 1))
              dpDEdgeF = 0.5 * (dp(i, j, k - 1) + dp(i, j + 1, k - 1))
              pGradY = kappaInv * MaInv2 / rhov * pEdgeF * ((dp(i, j + 1, k) &
                  &- dp(i, j, k)) / dy + met23EdgeF * (dpUEdgeF - dpDEdgeF) &
                  &* 0.5 / dz)
            end if

            dv = - dt * pGradY

            var%v(i, j, k) = var%v(i, j, k) + dv
          end do
        end do
      end do
    else
      stop 'ERROR: wrong opt in correctorStep'
    end if

    !--------------------------------------
    !         calc w and  w + dw
    !--------------------------------------

    select case(zBoundary)
    case("solid_wall")
      k0 = 1
      k1 = nz - 1
    case default
      stop "correctorStep: unknown case zBoundary."
    end select

    if(opt == "impl") then
      do k = k0, k1
        do j = 1, ny
          do i = 1, nx
            facw = 1.0

            if(spongeLayer) then
              facw = facw + dt * (jac(i, j, k + 1) * kr_sp_w_tfc(i, j, k) &
                  &+ jac(i, j, k) * kr_sp_w_tfc(i, j, k + 1)) / (jac(i, j, k) &
                  &+ jac(i, j, k + 1))
            end if

            ! Compute values at cell edges.
            rhoStratEdgeU = (jac(i, j, k + 1) * rhoStratTFC(i, j, k) + jac(i, &
                &j, k) * rhoStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, &
                &k + 1))
            rhoEdge = (jac(i, j, k + 1) * var%rho(i, j, k) + jac(i, j, k) &
                &* var%rho(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1)) &
                &+ rhoStratEdgeU
            pEdgeU = (jac(i, j, k + 1) * pStratTFC(i, j, k) + jac(i, j, k) &
                &* pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            bvsstw = (jac(i, j, k + 1) * bvsStratTFC(i, j, k) + jac(i, j, k) &
                &* bvsStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 1, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 2, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 3, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            dpREdgeU = (jac(i + 1, j, k + 1) * dp(i + 1, j, k) + jac(i + 1, j, &
                &k) * dp(i + 1, j, k + 1)) / (jac(i + 1, j, k) + jac(i + 1, j, &
                &k + 1))
            dpLEdgeU = (jac(i - 1, j, k + 1) * dp(i - 1, j, k) + jac(i - 1, j, &
                &k) * dp(i - 1, j, k + 1)) / (jac(i - 1, j, k) + jac(i - 1, j, &
                &k + 1))
            dpFEdgeU = (jac(i, j + 1, k + 1) * dp(i, j + 1, k) + jac(i, j + 1, &
                &k) * dp(i, j + 1, k + 1)) / (jac(i, j + 1, k) + jac(i, j + 1, &
                &k + 1))
            dpBEdgeU = (jac(i, j - 1, k + 1) * dp(i, j - 1, k) + jac(i, j - 1, &
                &k) * dp(i, j - 1, k + 1)) / (jac(i, j - 1, k) + jac(i, j - 1, &
                &k + 1))
            ! Compute pressure difference gradient component.
            pGradZ = kappaInv * MaInv2 / rhoEdge * pEdgeU * (met13EdgeU &
                &* (dpREdgeU - dpLEdgeU) * 0.5 / dx + met23EdgeU * (dpFEdgeU &
                &- dpBEdgeU) * 0.5 / dy + met33EdgeU * (dp(i, j, k + 1) &
                &- dp(i, j, k)) / dz)
            ! Compute velocity correction.
            dw = - facprs * dt / (facw + rhoStratEdgeU / rhoEdge * bvsstw * dt &
                &** 2.0) * pGradZ - 1.0 / (facw + rhoStratEdgeU / rhoEdge &
                &* bvsstw * dt ** 2.0) * rhoStratEdgeU / rhoEdge * bvsstw * dt &
                &** 2.0 * 0.5 * (jac(i, j, k + 1) * (met(i, j, k, 1, 3) &
                &* (corX(i, j, k) + corX(i - 1, j, k)) + met(i, j, k, 2, 3) &
                &* (corY(i, j, k) + corY(i, j - 1, k))) + jac(i, j, k) &
                &* (met(i, j, k + 1, 1, 3) * (corX(i, j, k + 1) + corX(i - 1, &
                &j, k + 1)) + met(i, j, k + 1, 2, 3) * (corY(i, j, k + 1) &
                &+ corY(i, j - 1, k + 1)))) / (jac(i, j, k) + jac(i, j, k + 1))

            var%w(i, j, k) = var%w(i, j, k) + dw
          end do
        end do
      end do

    else if(opt == "expl") then
      if(facprs /= 1.) stop 'ERROR: wrong facprs in explicit sub-step'
      do k = k0, k1
        do j = 1, ny
          do i = 1, nx
            ! Compute values at cell edges.
            rhoStratEdgeU = (jac(i, j, k + 1) * rhoStratTFC(i, j, k) + jac(i, &
                &j, k) * rhoStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, &
                &k + 1))
            rhoEdge = (jac(i, j, k + 1) * var%rho(i, j, k) + jac(i, j, k) &
                &* var%rho(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1)) &
                &+ rhoStratEdgeU
            pEdgeU = (jac(i, j, k + 1) * pStratTFC(i, j, k) + jac(i, j, k) &
                &* pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            bvsstw = (jac(i, j, k + 1) * bvsStratTFC(i, j, k) + jac(i, j, k) &
                &* bvsStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 1, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 2, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 3, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            dpREdgeU = (jac(i + 1, j, k + 1) * dp(i + 1, j, k) + jac(i + 1, j, &
                &k) * dp(i + 1, j, k + 1)) / (jac(i + 1, j, k) + jac(i + 1, j, &
                &k + 1))
            dpLEdgeU = (jac(i - 1, j, k + 1) * dp(i - 1, j, k) + jac(i - 1, j, &
                &k) * dp(i - 1, j, k + 1)) / (jac(i - 1, j, k) + jac(i - 1, j, &
                &k + 1))
            dpFEdgeU = (jac(i, j + 1, k + 1) * dp(i, j + 1, k) + jac(i, j + 1, &
                &k) * dp(i, j + 1, k + 1)) / (jac(i, j + 1, k) + jac(i, j + 1, &
                &k + 1))
            dpBEdgeU = (jac(i, j - 1, k + 1) * dp(i, j - 1, k) + jac(i, j - 1, &
                &k) * dp(i, j - 1, k + 1)) / (jac(i, j - 1, k) + jac(i, j - 1, &
                &k + 1))
            ! Compute pressure difference gradient component.
            pGradZ = kappaInv * MaInv2 / rhoEdge * pEdgeU * (met13EdgeU &
                &* (dpREdgeU - dpLEdgeU) * 0.5 / dx + met23EdgeU * (dpFEdgeU &
                &- dpBEdgeU) * 0.5 / dy + met33EdgeU * (dp(i, j, k + 1) &
                &- dp(i, j, k)) / dz)
            ! Correct vertical velocity.
            dw = - dt * pGradZ

            var%w(i, j, k) = var%w(i, j, k) + dw
          end do
        end do
      end do

    else
      stop 'ERROR: wrong opt in correctorStep'
    end if

    !------------------------------------------------------------------
    !         calc rhop and rhop + drhop (only for implicit time step)
    !------------------------------------------------------------------

    if(opt == "impl") then
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            facw = 1.0

            if(spongeLayer) then
              facw = facw + dt * kr_sp_w_tfc(i, j, k)
            end if

            ! Compute P coefficients.
            pEdgeU = (jac(i, j, k + 1) * pStratTFC(i, j, k) + jac(i, j, k) &
                &* pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            pEdgeD = (jac(i, j, k - 1) * pStratTFC(i, j, k) + jac(i, j, k) &
                &* pStratTFC(i, j, k - 1)) / (jac(i, j, k) + jac(i, j, k - 1))
            ! Compute density coefficients.
            rhow0 = (jac(i, j, k + 1) * (var%rho(i, j, k) + rhoStratTFC(i, j, &
                &k)) + jac(i, j, k) * (var%rho(i, j, k + 1) + rhoStratTFC(i, &
                &j, k + 1))) / (jac(i, j, k) + jac(i, j, k + 1))
            rhowm = (jac(i, j, k - 1) * (var%rho(i, j, k) + rhoStratTFC(i, j, &
                &k)) + jac(i, j, k) * (var%rho(i, j, k - 1) + rhoStratTFC(i, &
                &j, k - 1))) / (jac(i, j, k) + jac(i, j, k - 1))
            rho = var%rho(i, j, k) + rhoStratTFC(i, j, k)
            ! Interpolate metric tensor elements.
            met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 1, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            met13EdgeD = (jac(i, j, k - 1) * met(i, j, k, 1, 3) + jac(i, j, k) &
                &* met(i, j, k - 1, 1, 3)) / (jac(i, j, k) + jac(i, j, k - 1))
            met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 2, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            met23EdgeD = (jac(i, j, k - 1) * met(i, j, k, 2, 3) + jac(i, j, k) &
                &* met(i, j, k - 1, 2, 3)) / (jac(i, j, k) + jac(i, j, k - 1))
            met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 3, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            met33EdgeD = (jac(i, j, k - 1) * met(i, j, k, 3, 3) + jac(i, j, k) &
                &* met(i, j, k - 1, 3, 3)) / (jac(i, j, k) + jac(i, j, k - 1))
            ! Interpolate pressure differences.
            dpREdgeU = (jac(i + 1, j, k + 1) * dp(i + 1, j, k) + jac(i + 1, j, &
                &k) * dp(i + 1, j, k + 1)) / (jac(i + 1, j, k) + jac(i + 1, j, &
                &k + 1))
            dpLEdgeU = (jac(i - 1, j, k + 1) * dp(i - 1, j, k) + jac(i - 1, j, &
                &k) * dp(i - 1, j, k + 1)) / (jac(i - 1, j, k) + jac(i - 1, j, &
                &k + 1))
            dpREdgeD = (jac(i + 1, j, k - 1) * dp(i + 1, j, k) + jac(i + 1, j, &
                &k) * dp(i + 1, j, k - 1)) / (jac(i + 1, j, k) + jac(i + 1, j, &
                &k - 1))
            dpLEdgeD = (jac(i - 1, j, k - 1) * dp(i - 1, j, k) + jac(i - 1, j, &
                &k) * dp(i - 1, j, k - 1)) / (jac(i - 1, j, k) + jac(i - 1, j, &
                &k - 1))
            dpFEdgeU = (jac(i, j + 1, k + 1) * dp(i, j + 1, k) + jac(i, j + 1, &
                &k) * dp(i, j + 1, k + 1)) / (jac(i, j + 1, k) + jac(i, j + 1, &
                &k + 1))
            dpBEdgeU = (jac(i, j - 1, k + 1) * dp(i, j - 1, k) + jac(i, j - 1, &
                &k) * dp(i, j - 1, k + 1)) / (jac(i, j - 1, k) + jac(i, j - 1, &
                &k + 1))
            dpFEdgeD = (jac(i, j + 1, k - 1) * dp(i, j + 1, k) + jac(i, j + 1, &
                &k) * dp(i, j + 1, k - 1)) / (jac(i, j + 1, k) + jac(i, j + 1, &
                &k - 1))
            dpBEdgeD = (jac(i, j - 1, k - 1) * dp(i, j - 1, k) + jac(i, j - 1, &
                &k) * dp(i, j - 1, k - 1)) / (jac(i, j - 1, k) + jac(i, j - 1, &
                &k - 1))
            ! Compute pressure difference gradients.
            pGradZEdgeU = kappaInv * MaInv2 * pEdgeU / rhow0 * (0.5 &
                &* met13EdgeU * (dpREdgeU - dpLEdgeU) / dx + 0.5 * met23EdgeU &
                &* (dpFEdgeU - dpBEdgeU) / dy + met33EdgeU * (dp(i, j, k + 1) &
                &- dp(i, j, k)) / dz)
            pGradZEdgeD = kappaInv * MaInv2 * pEdgeD / rhowm * (0.5 &
                &* met13EdgeD * (dpREdgeD - dpLEdgeD) / dx + 0.5 * met23EdgeD &
                &* (dpFEdgeD - dpBEdgeD) / dy + met33EdgeD * (dp(i, j, k) &
                &- dp(i, j, k - 1)) / dz)
            ! Adjust at boundaries.
            if(k == 1 .and. zBoundary == "solid_wall") then
              pGradZEdgeD = 0.0
            else if(k == nz .and. zBoundary == "solid_wall") then
              pGradZEdgeU = 0.0
            end if
            ! Interpolate.
            pGradZ = 0.5 * (pGradZEdgeU + pGradZEdgeD)
            ! Compute buoyancy correction.
            db = - 1.0 / (facw + rhoStratTFC(i, j, k) / rho * bvsStratTFC(i, &
                &j, k) * dt ** 2.0) * (- rhoStratTFC(i, j, k) / rho &
                &* bvsStratTFC(i, j, k) * facprs * dt ** 2.0 * jac(i, j, k) &
                &* pGradZ + rhoStratTFC(i, j, k) / rho * bvsStratTFC(i, j, k) &
                &* dt * jac(i, j, k) * facw * 0.5 * (met(i, j, k, 1, 3) &
                &* (corX(i, j, k) + corX(i - 1, j, k)) + met(i, j, k, 2, 3) &
                &* (corY(i, j, k) + corY(i, j - 1, k))))

            var%rhop(i, j, k) = var%rhop(i, j, k) - rho / g_ndim * db
          end do
        end do
      end do

    end if

    if(opt == "impl") then
      kr_sp_tfc = kr_sp_tfc / facray
      kr_sp_w_tfc = kr_sp_w_tfc / facray
    end if

  end subroutine correctorStep

  !----------------------------------------------------------------------

  subroutine init_poisson

    !-----------------------------------
    ! allocate poisson module variables
    !-----------------------------------

    ! local variables
    integer :: allocstat

    ! allocate fields
    allocate(dp(0:nx + 1, 0:ny + 1, 0:nz + 1), stat = allocstat)
    if(allocstat /= 0) stop "init_poisson: alloc failed"

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

  end subroutine terminate_poisson

  !----------------------------------------------------------------------

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
    real :: AL, AR, AB, AF, AD, AU, AC, ACH, ACV

    integer :: i, j, k

    real :: fcscal, fcscal_u, fcscal_d

    real :: ARU, ARD, ALU, ALD, AFU, AFD, ABU, ABD
    real :: AUU, ADD, ARUU, ARDD, ALUU, ALDD, AFUU, AFDD, ABUU, ABDD
    real :: fcscal_r, fcscal_l, fcscal_f, fcscal_b, fcscal_ru, fcscal_rd, &
        &fcscal_lu, fcscal_ld, fcscal_fu, fcscal_fd, fcscal_bu, fcscal_bd
    real :: fcscal_uu, fcscal_dd, fcscal_ruu, fcscal_rdd, fcscal_luu, &
        &fcscal_ldd, fcscal_fuu, fcscal_fdd, fcscal_buu, fcscal_bdd
    real :: jacInv
    real :: pEdgeRDiv, pEdgeLDiv, pEdgeFDiv, pEdgeBDiv, pEdgeUDiv, pEdgeDDiv
    real :: pEdgeRGra, pEdgeLGra, pEdgeFGra, pEdgeBGra, pEdgeUGra, pEdgeDGra
    real :: rhoEdgeR, rhoEdgeL, rhoEdgeF, rhoEdgeB, rhoEdgeU, rhoEdgeD
    real :: met13EdgeR, met13EdgeL, met23EdgeF, met23EdgeB, met13EdgeU, &
        &met23EdgeU, met33EdgeU, met13EdgeD, met23EdgeD, met33EdgeD, &
        &met13UEdgeR, met13UEdgeL, met23UEdgeF, met23UEdgeB, met13DEdgeR, &
        &met13DEdgeL, met23DEdgeF, met23DEdgeB
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

    !---------------------------------
    !         Loop over field
    !---------------------------------

    if(opt == "expl") then
      ! Compute tensor elements for TFC.
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            ! Compute scaling factors.
            fcscal = sqrt(pStratTFC(i, j, k) ** 2.0 / rhoStratTFC(i, j, k))
            fcscal_r = sqrt(pStratTFC(i + 1, j, k) ** 2.0 / rhoStratTFC(i + 1, &
                &j, k))
            fcscal_l = sqrt(pStratTFC(i - 1, j, k) ** 2.0 / rhoStratTFC(i - 1, &
                &j, k))
            fcscal_f = sqrt(pStratTFC(i, j + 1, k) ** 2.0 / rhoStratTFC(i, j &
                &+ 1, k))
            fcscal_b = sqrt(pStratTFC(i, j - 1, k) ** 2.0 / rhoStratTFC(i, j &
                &- 1, k))
            fcscal_u = sqrt(pStratTFC(i, j, k + 1) ** 2.0 / rhoStratTFC(i, j, &
                &k + 1))
            fcscal_d = sqrt(pStratTFC(i, j, k - 1) ** 2.0 / rhoStratTFC(i, j, &
                &k - 1))
            fcscal_ru = sqrt(pStratTFC(i + 1, j, k + 1) ** 2.0 / rhoStratTFC(i &
                &+ 1, j, k + 1))
            fcscal_rd = sqrt(pStratTFC(i + 1, j, k - 1) ** 2.0 / rhoStratTFC(i &
                &+ 1, j, k - 1))
            fcscal_lu = sqrt(pStratTFC(i - 1, j, k + 1) ** 2.0 / rhoStratTFC(i &
                &- 1, j, k + 1))
            fcscal_ld = sqrt(pStratTFC(i - 1, j, k - 1) ** 2.0 / rhoStratTFC(i &
                &- 1, j, k - 1))
            fcscal_fu = sqrt(pStratTFC(i, j + 1, k + 1) ** 2.0 &
                &/ rhoStratTFC(i, j + 1, k + 1))
            fcscal_fd = sqrt(pStratTFC(i, j + 1, k - 1) ** 2.0 &
                &/ rhoStratTFC(i, j + 1, k - 1))
            fcscal_bu = sqrt(pStratTFC(i, j - 1, k + 1) ** 2.0 &
                &/ rhoStratTFC(i, j - 1, k + 1))
            fcscal_bd = sqrt(pStratTFC(i, j - 1, k - 1) ** 2.0 &
                &/ rhoStratTFC(i, j - 1, k - 1))
            fcscal_uu = sqrt(pStratTFC(i, j, k + 2) ** 2.0 / rhoStratTFC(i, j, &
                &k + 2))
            fcscal_dd = sqrt(pStratTFC(i, j, k - 2) ** 2.0 / rhoStratTFC(i, j, &
                &k - 2))
            fcscal_ruu = sqrt(pStratTFC(i + 1, j, k + 2) ** 2.0 &
                &/ rhoStratTFC(i + 1, j, k + 2))
            fcscal_rdd = sqrt(pStratTFC(i + 1, j, k - 2) ** 2.0 &
                &/ rhoStratTFC(i + 1, j, k - 2))
            fcscal_luu = sqrt(pStratTFC(i - 1, j, k + 2) ** 2.0 &
                &/ rhoStratTFC(i - 1, j, k + 2))
            fcscal_ldd = sqrt(pStratTFC(i - 1, j, k - 2) ** 2.0 &
                &/ rhoStratTFC(i - 1, j, k - 2))
            fcscal_fuu = sqrt(pStratTFC(i, j + 1, k + 2) ** 2.0 &
                &/ rhoStratTFC(i, j + 1, k + 2))
            fcscal_fdd = sqrt(pStratTFC(i, j + 1, k - 2) ** 2.0 &
                &/ rhoStratTFC(i, j + 1, k - 2))
            fcscal_buu = sqrt(pStratTFC(i, j - 1, k + 2) ** 2.0 &
                &/ rhoStratTFC(i, j - 1, k + 2))
            fcscal_bdd = sqrt(pStratTFC(i, j - 1, k - 2) ** 2.0 &
                &/ rhoStratTFC(i, j - 1, k - 2))

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
            pEdgeUDiv = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
                &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            pEdgeDDiv = jac(i, j, k) * jac(i, j, k - 1) * (pStratTFC(i, j, k) &
                &+ pStratTFC(i, j, k - 1)) / (jac(i, j, k) + jac(i, j, k - 1))

            ! Compute P coefficients (pressure gradient).
            pEdgeRGra = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i + 1, j, k))
            pEdgeLGra = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i - 1, j, k))
            pEdgeFGra = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j + 1, k))
            pEdgeBGra = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j - 1, k))
            pEdgeUGra = (jac(i, j, k + 1) * pStratTFC(i, j, k) + jac(i, j, k) &
                &* pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            pEdgeDGra = (jac(i, j, k - 1) * pStratTFC(i, j, k) + jac(i, j, k) &
                &* pStratTFC(i, j, k - 1)) / (jac(i, j, k) + jac(i, j, k - 1))

            ! Compute density coefficients.
            rhoEdgeR = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, j, k))
            rhoEdgeL = 0.5 * (var%rho(i, j, k) + var%rho(i - 1, j, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i - 1, j, k))
            rhoEdgeF = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j + 1, k))
            rhoEdgeB = 0.5 * (var%rho(i, j, k) + var%rho(i, j - 1, k) &
                &+ rhoStratTFC(i, j, k) + rhoStratTFC(i, j - 1, k))
            rhoEdgeU = (jac(i, j, k + 1) * (var%rho(i, j, k) + rhoStratTFC(i, &
                &j, k)) + jac(i, j, k) * (var%rho(i, j, k + 1) &
                &+ rhoStratTFC(i, j, k + 1))) / (jac(i, j, k) + jac(i, j, k &
                &+ 1))
            rhoEdgeD = (jac(i, j, k - 1) * (var%rho(i, j, k) + rhoStratTFC(i, &
                &j, k)) + jac(i, j, k) * (var%rho(i, j, k - 1) &
                &+ rhoStratTFC(i, j, k - 1))) / (jac(i, j, k) + jac(i, j, k &
                &- 1))

            ! Interpolate metric-tensor elements.
            met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))
            met13EdgeL = 0.5 * (met(i, j, k, 1, 3) + met(i - 1, j, k, 1, 3))
            met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3))
            met23EdgeB = 0.5 * (met(i, j, k, 2, 3) + met(i, j - 1, k, 2, 3))
            met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 1, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 2, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 3, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            met13EdgeD = (jac(i, j, k - 1) * met(i, j, k, 1, 3) + jac(i, j, k) &
                &* met(i, j, k - 1, 1, 3)) / (jac(i, j, k) + jac(i, j, k - 1))
            met23EdgeD = (jac(i, j, k - 1) * met(i, j, k, 2, 3) + jac(i, j, k) &
                &* met(i, j, k - 1, 2, 3)) / (jac(i, j, k) + jac(i, j, k - 1))
            met33EdgeD = (jac(i, j, k - 1) * met(i, j, k, 3, 3) + jac(i, j, k) &
                &* met(i, j, k - 1, 3, 3)) / (jac(i, j, k) + jac(i, j, k - 1))

            ! --------------------- A(i,j,k) ---------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              AC = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * (1.0 &
                  &/ dx + 0.75 * met13EdgeR / dz) + pEdgeLDiv / rhoEdgeL &
                  &* pEdgeLGra * (1.0 / dx - 0.75 * met13EdgeL / dz)) - jacInv &
                  &/ dy * (pEdgeFDiv / rhoEdgeF * pEdgeFGra * (1.0 / dy + 0.75 &
                  &* met23EdgeF / dz) + pEdgeBDiv / rhoEdgeB * pEdgeBGra &
                  &* (1.0 / dy - 0.75 * met23EdgeB / dz)) - jacInv / dz &
                  &* pEdgeUDiv / rhoEdgeU * pEdgeUGra * met33EdgeU / dz
            else if(k == nz .and. zBoundary == "solid_wall") then
              AC = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * (1.0 &
                  &/ dx - 0.75 * met13EdgeR / dz) + pEdgeLDiv / rhoEdgeL &
                  &* pEdgeLGra * (1.0 / dx + 0.75 * met13EdgeL / dz)) - jacInv &
                  &/ dy * (pEdgeFDiv / rhoEdgeF * pEdgeFGra * (1.0 / dy - 0.75 &
                  &* met23EdgeF / dz) + pEdgeBDiv / rhoEdgeB * pEdgeBGra &
                  &* (1.0 / dy + 0.75 * met23EdgeB / dz)) - jacInv / dz &
                  &* pEdgeDDiv / rhoEdgeD * pEdgeDGra * met33EdgeD / dz
            else
              AC = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra / dx &
                  &+ pEdgeLDiv / rhoEdgeL * pEdgeLGra / dx) - jacInv / dy &
                  &* (pEdgeFDiv / rhoEdgeF * pEdgeFGra / dy + pEdgeBDiv &
                  &/ rhoEdgeB * pEdgeBGra / dy) - jacInv / dz * (pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met33EdgeU / dz + pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * met33EdgeD / dz)
            end if
            ! -------------------- A(i+1,j,k) --------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              AR = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * (1.0 / dx &
                  &- 0.75 * met13EdgeR / dz) + jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx * jac(i + 1, &
                  &j, k + 1) / (jac(i + 1, j, k) + jac(i + 1, j, k + 1))
            else if(k == nz .and. zBoundary == "solid_wall") then
              AR = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * (1.0 / dx &
                  &+ 0.75 * met13EdgeR / dz) - jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx * jac(i + 1, &
                  &j, k - 1) / (jac(i + 1, j, k) + jac(i + 1, j, k - 1))
            else
              AR = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra / dx &
                  &+ jacInv / dz * (pEdgeUDiv / rhoEdgeU * pEdgeUGra &
                  &* met13EdgeU * 0.5 / dx * jac(i + 1, j, k + 1) / (jac(i &
                  &+ 1, j, k) + jac(i + 1, j, k + 1)) - pEdgeDDiv / rhoEdgeD &
                  &* pEdgeDGra * met13EdgeD * 0.5 / dx * jac(i + 1, j, k - 1) &
                  &/ (jac(i + 1, j, k) + jac(i + 1, j, k - 1)))
            end if

            ! -------------------- A(i-1,j,k) --------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              AL = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * (1.0 / dx &
                  &+ 0.75 * met13EdgeL / dz) - jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx * jac(i - 1, &
                  &j, k + 1) / (jac(i - 1, j, k) + jac(i - 1, j, k + 1))
            else if(k == nz .and. zBoundary == "solid_wall") then
              AL = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * (1.0 / dx &
                  &- 0.75 * met13EdgeL / dz) + jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx * jac(i - 1, &
                  &j, k - 1) / (jac(i - 1, j, k) + jac(i - 1, j, k - 1))
            else
              AL = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra / dx &
                  &- jacInv / dz * (pEdgeUDiv / rhoEdgeU * pEdgeUGra &
                  &* met13EdgeU * 0.5 / dx * jac(i - 1, j, k + 1) / (jac(i &
                  &- 1, j, k) + jac(i - 1, j, k + 1)) - pEdgeDDiv / rhoEdgeD &
                  &* pEdgeDGra * met13EdgeD * 0.5 / dx * jac(i - 1, j, k - 1) &
                  &/ (jac(i - 1, j, k) + jac(i - 1, j, k - 1)))
            end if

            ! -------------------- A(i,j+1,k) --------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              AF = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * (1.0 / dy &
                  &- 0.75 * met23EdgeF / dz) + jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy * jac(i, j &
                  &+ 1, k + 1) / (jac(i, j + 1, k) + jac(i, j + 1, k + 1))
            else if(k == nz .and. zBoundary == "solid_wall") then
              AF = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * (1.0 / dy &
                  &+ 0.75 * met23EdgeF / dz) - jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy * jac(i, j &
                  &+ 1, k - 1) / (jac(i, j + 1, k) + jac(i, j + 1, k - 1))
            else
              AF = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra / dy &
                  &+ jacInv / dz * (pEdgeUDiv / rhoEdgeU * pEdgeUGra &
                  &* met23EdgeU * 0.5 / dy * jac(i, j + 1, k + 1) / (jac(i, j &
                  &+ 1, k) + jac(i, j + 1, k + 1)) - pEdgeDDiv / rhoEdgeD &
                  &* pEdgeDGra * met23EdgeD * 0.5 / dy * jac(i, j + 1, k - 1) &
                  &/ (jac(i, j + 1, k) + jac(i, j + 1, k - 1)))
            end if

            ! -------------------- A(i,j-1,k) --------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              AB = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * (1.0 / dy &
                  &+ 0.75 * met23EdgeB / dz) - jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy * jac(i, j &
                  &- 1, k + 1) / (jac(i, j - 1, k) + jac(i, j - 1, k + 1))
            else if(k == nz .and. zBoundary == "solid_wall") then
              AB = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * (1.0 / dy &
                  &- 0.75 * met23EdgeB / dz) + jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy * jac(i, j &
                  &- 1, k - 1) / (jac(i, j - 1, k) + jac(i, j - 1, k - 1))
            else
              AB = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra / dy &
                  &- jacInv / dz * (pEdgeUDiv / rhoEdgeU * pEdgeUGra &
                  &* met23EdgeU * 0.5 / dy * jac(i, j - 1, k + 1) / (jac(i, j &
                  &- 1, k) + jac(i, j - 1, k + 1)) - pEdgeDDiv / rhoEdgeD &
                  &* pEdgeDGra * met23EdgeD * 0.5 / dy * jac(i, j - 1, k - 1) &
                  &/ (jac(i, j - 1, k) + jac(i, j - 1, k - 1)))
            end if

            ! -------------------- A(i,j,k+1) --------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              AU = jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra &
                  &* met13EdgeR / dz - pEdgeLDiv / rhoEdgeL * pEdgeLGra &
                  &* met13EdgeL / dz) + jacInv / dy * (pEdgeFDiv / rhoEdgeF &
                  &* pEdgeFGra * met23EdgeF / dz - pEdgeBDiv / rhoEdgeB &
                  &* pEdgeBGra * met23EdgeB / dz) + jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met33EdgeU / dz
            else if(k == nz .and. zBoundary == "solid_wall") then
              AU = 0.0
            else
              AU = jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra &
                  &* met13EdgeR * 0.25 / dz - pEdgeLDiv / rhoEdgeL * pEdgeLGra &
                  &* met13EdgeL * 0.25 / dz) + jacInv / dy * (pEdgeFDiv &
                  &/ rhoEdgeF * pEdgeFGra * met23EdgeF * 0.25 / dz - pEdgeBDiv &
                  &/ rhoEdgeB * pEdgeBGra * met23EdgeB * 0.25 / dz) + jacInv &
                  &/ dz * pEdgeUDiv / rhoEdgeU * pEdgeUGra * met33EdgeU / dz
            end if

            ! -------------------- A(i,j,k-1) --------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              AD = 0.0
            else if(k == nz .and. zBoundary == "solid_wall") then
              AD = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra &
                  &* met13EdgeR / dz - pEdgeLDiv / rhoEdgeL * pEdgeLGra &
                  &* met13EdgeL / dz) - jacInv / dy * (pEdgeFDiv / rhoEdgeF &
                  &* pEdgeFGra * met23EdgeF / dz - pEdgeBDiv / rhoEdgeB &
                  &* pEdgeBGra * met23EdgeB / dz) + jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * met33EdgeD / dz
            else
              AD = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra &
                  &* met13EdgeR * 0.25 / dz - pEdgeLDiv / rhoEdgeL * pEdgeLGra &
                  &* met13EdgeL * 0.25 / dz) - jacInv / dy * (pEdgeFDiv &
                  &/ rhoEdgeF * pEdgeFGra * met23EdgeF * 0.25 / dz - pEdgeBDiv &
                  &/ rhoEdgeB * pEdgeBGra * met23EdgeB * 0.25 / dz) + jacInv &
                  &/ dz * pEdgeDDiv / rhoEdgeD * pEdgeDGra * met33EdgeD / dz
            end if

            ! ------------------- A(i+1,j,k+1) -------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              ARU = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra &
                  &* met13EdgeR / dz + jacInv / dz * pEdgeUDiv / rhoEdgeU &
                  &* pEdgeUGra * met13EdgeU * 0.5 / dx * jac(i + 1, j, k) &
                  &/ (jac(i + 1, j, k) + jac(i + 1, j, k + 1))
            else if(k == nz .and. zBoundary == "solid_wall") then
              ARU = 0.0
            else
              ARU = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra &
                  &* met13EdgeR * 0.25 / dz + jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met13EdgeU * 0.5 / dx * jac(i + 1, &
                  &j, k) / (jac(i + 1, j, k) + jac(i + 1, j, k + 1))
            end if

            ! ------------------- A(i+1,j,k-1) -------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              ARD = 0.0
            else if(k == nz .and. zBoundary == "solid_wall") then
              ARD = - jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra &
                  &* met13EdgeR / dz - jacInv / dz * pEdgeDDiv / rhoEdgeD &
                  &* pEdgeDGra * met13EdgeD * 0.5 / dx * jac(i + 1, j, k) &
                  &/ (jac(i + 1, j, k) + jac(i + 1, j, k - 1))
            else
              ARD = - jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra &
                  &* met13EdgeR * 0.25 / dz - jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * met13EdgeD * 0.5 / dx * jac(i + 1, &
                  &j, k) / (jac(i + 1, j, k) + jac(i + 1, j, k - 1))
            end if

            ! ------------------- A(i-1,j,k+1) -------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              ALU = - jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra &
                  &* met13EdgeL / dz - jacInv / dz * pEdgeUDiv / rhoEdgeU &
                  &* pEdgeUGra * met13EdgeU * 0.5 / dx * jac(i - 1, j, k) &
                  &/ (jac(i - 1, j, k) + jac(i - 1, j, k + 1))
            else if(k == nz .and. zBoundary == "solid_wall") then
              ALU = 0.0
            else
              ALU = - jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra &
                  &* met13EdgeL * 0.25 / dz - jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met13EdgeU * 0.5 / dx * jac(i - 1, &
                  &j, k) / (jac(i - 1, j, k) + jac(i - 1, j, k + 1))
            end if

            ! ------------------- A(i-1,j,k-1) -------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              ALD = 0.0
            else if(k == nz .and. zBoundary == "solid_wall") then
              ALD = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra &
                  &* met13EdgeL / dz + jacInv / dz * pEdgeDDiv / rhoEdgeD &
                  &* pEdgeDGra * met13EdgeD * 0.5 / dx * jac(i - 1, j, k) &
                  &/ (jac(i - 1, j, k) + jac(i - 1, j, k - 1))
            else
              ALD = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra &
                  &* met13EdgeL * 0.25 / dz + jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * met13EdgeD * 0.5 / dx * jac(i - 1, &
                  &j, k) / (jac(i - 1, j, k) + jac(i - 1, j, k - 1))
            end if

            ! ------------------- A(i,j+1,k+1) -------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              AFU = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra &
                  &* met23EdgeF / dz + jacInv / dz * pEdgeUDiv / rhoEdgeU &
                  &* pEdgeUGra * met23EdgeU * 0.5 / dy * jac(i, j + 1, k) &
                  &/ (jac(i, j + 1, k) + jac(i, j + 1, k + 1))
            else if(k == nz .and. zBoundary == "solid_wall") then
              AFU = 0.0
            else
              AFU = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra &
                  &* met23EdgeF * 0.25 / dz + jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met23EdgeU * 0.5 / dy * jac(i, j &
                  &+ 1, k) / (jac(i, j + 1, k) + jac(i, j + 1, k + 1))
            end if

            ! ------------------- A(i,j+1,k-1) -------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              AFD = 0.0
            else if(k == nz .and. zBoundary == "solid_wall") then
              AFD = - jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra &
                  &* met23EdgeF / dz - jacInv / dz * pEdgeDDiv / rhoEdgeD &
                  &* pEdgeDGra * met23EdgeD * 0.5 / dy * jac(i, j + 1, k) &
                  &/ (jac(i, j + 1, k) + jac(i, j + 1, k - 1))
            else
              AFD = - jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra &
                  &* met23EdgeF * 0.25 / dz - jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * met23EdgeD * 0.5 / dy * jac(i, j &
                  &+ 1, k) / (jac(i, j + 1, k) + jac(i, j + 1, k - 1))
            end if

            ! ------------------- A(i,j-1,k+1) -------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              ABU = - jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra &
                  &* met23EdgeB / dz - jacInv / dz * pEdgeUDiv / rhoEdgeU &
                  &* pEdgeUGra * met23EdgeU * 0.5 / dy * jac(i, j - 1, k) &
                  &/ (jac(i, j - 1, k) + jac(i, j - 1, k + 1))
            else if(k == nz .and. zBoundary == "solid_wall") then
              ABU = 0.0
            else
              ABU = - jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra &
                  &* met23EdgeB * 0.25 / dz - jacInv / dz * pEdgeUDiv &
                  &/ rhoEdgeU * pEdgeUGra * met23EdgeU * 0.5 / dy * jac(i, j &
                  &- 1, k) / (jac(i, j - 1, k) + jac(i, j - 1, k + 1))
            end if

            ! ------------------- A(i,j-1,k-1) -------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              ABD = 0.0
            else if(k == nz .and. zBoundary == "solid_wall") then
              ABD = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra &
                  &* met23EdgeB / dz + jacInv / dz * pEdgeDDiv / rhoEdgeD &
                  &* pEdgeDGra * met23EdgeD * 0.5 / dy * jac(i, j - 1, k) &
                  &/ (jac(i, j - 1, k) + jac(i, j - 1, k - 1))
            else
              ABD = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra &
                  &* met23EdgeB * 0.25 / dz + jacInv / dz * pEdgeDDiv &
                  &/ rhoEdgeD * pEdgeDGra * met23EdgeD * 0.5 / dy * jac(i, j &
                  &- 1, k) / (jac(i, j - 1, k) + jac(i, j - 1, k - 1))
            end if

            ! ------------------- A(i,j,k+2) ---------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              AUU = - jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * 0.25 &
                  &* met13EdgeR / dz - pEdgeLDiv / rhoEdgeL * pEdgeLGra * 0.25 &
                  &* met13EdgeL / dz) - jacInv / dy * (pEdgeFDiv / rhoEdgeF &
                  &* pEdgeFGra * 0.25 * met23EdgeF / dz - pEdgeBDiv / rhoEdgeB &
                  &* pEdgeBGra * 0.25 * met23EdgeB / dz)
            else
              AUU = 0.0
            end if

            ! ------------------- A(i,j,k-2) ---------------------

            if(k == nz .and. zBoundary == "solid_wall") then
              ADD = jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra * 0.25 &
                  &* met13EdgeR / dz - pEdgeLDiv / rhoEdgeL * pEdgeLGra * 0.25 &
                  &* met13EdgeL / dz) + jacInv / dy * (pEdgeFDiv / rhoEdgeF &
                  &* pEdgeFGra * 0.25 * met23EdgeF / dz - pEdgeBDiv / rhoEdgeB &
                  &* pEdgeBGra * 0.25 * met23EdgeB / dz)
            else
              ADD = 0.0
            end if

            ! ------------------ A(i+1,j,k+2) --------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              ARUU = - jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * 0.25 &
                  &* met13EdgeR / dz
            else
              ARUU = 0.0
            end if

            ! ------------------ A(i+1,j,k-2) --------------------

            if(k == nz .and. zBoundary == "solid_wall") then
              ARDD = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * 0.25 &
                  &* met13EdgeR / dz
            else
              ARDD = 0.0
            end if

            ! ------------------ A(i-1,j,k+2) --------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              ALUU = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * 0.25 &
                  &* met13EdgeL / dz
            else
              ALUU = 0.0
            end if

            ! ------------------ A(i-1,j,k-2) --------------------

            if(k == nz .and. zBoundary == "solid_wall") then
              ALDD = - jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * 0.25 &
                  &* met13EdgeL / dz
            else
              ALDD = 0.0
            end if

            ! ------------------ A(i,j+1,k+2) --------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              AFUU = - jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * 0.25 &
                  &* met23EdgeF / dz
            else
              AFUU = 0.0
            end if

            ! ------------------ A(i,j+1,k-2) --------------------

            if(k == nz .and. zBoundary == "solid_wall") then
              AFDD = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * 0.25 &
                  &* met23EdgeF / dz
            else
              AFDD = 0.0
            end if

            ! ------------------ A(i,j-1,k+2) --------------------

            if(k == 1 .and. zBoundary == "solid_wall") then
              ABUU = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * 0.25 &
                  &* met23EdgeB / dz
            else
              ABUU = 0.0
            end if

            ! ------------------ A(i,j-1,k-2) --------------------

            if(k == nz .and. zBoundary == "solid_wall") then
              ABDD = - jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * 0.25 &
                  &* met23EdgeB / dz
            else
              ABDD = 0.0
            end if

            ! Scale the tensor elements.
            AC = AC / (fcscal ** 2.0)
            AR = AR / fcscal / fcscal_r
            AL = AL / fcscal / fcscal_l
            AF = AF / fcscal / fcscal_f
            AB = AB / fcscal / fcscal_b
            AU = AU / fcscal / fcscal_u
            AD = AD / fcscal / fcscal_d
            ARU = ARU / fcscal / fcscal_ru
            ARD = ARD / fcscal / fcscal_rd
            ALU = ALU / fcscal / fcscal_lu
            ALD = ALD / fcscal / fcscal_ld
            AFU = AFU / fcscal / fcscal_fu
            AFD = AFD / fcscal / fcscal_fd
            ABU = ABU / fcscal / fcscal_bu
            ABD = ABD / fcscal / fcscal_bd
            AUU = AUU / fcscal / fcscal_uu
            ADD = ADD / fcscal / fcscal_dd
            ARUU = ARUU / fcscal / fcscal_ruu
            ARDD = ARDD / fcscal / fcscal_rdd
            ALUU = ALUU / fcscal / fcscal_luu
            ALDD = ALDD / fcscal / fcscal_ldd
            AFUU = AFUU / fcscal / fcscal_fuu
            AFDD = AFDD / fcscal / fcscal_fdd
            ABUU = ABUU / fcscal / fcscal_buu
            ABDD = ABDD / fcscal / fcscal_bdd

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
    else if(opt == "impl") then

      ! Compute tensor elements for TFC.
      kr_sp_tfc = kr_sp_tfc * facray
      kr_sp_w_tfc = kr_sp_w_tfc * facray
      do k = 1, nz
        do j = 1, ny
          do i = 1, nx
            ! Compute scaling factors.
            fcscal = sqrt(pStratTFC(i, j, k) ** 2.0 / rhoStratTFC(i, j, k))
            fcscal_r = sqrt(pStratTFC(i + 1, j, k) ** 2.0 / rhoStratTFC(i + 1, &
                &j, k))
            fcscal_l = sqrt(pStratTFC(i - 1, j, k) ** 2.0 / rhoStratTFC(i - 1, &
                &j, k))
            fcscal_f = sqrt(pStratTFC(i, j + 1, k) ** 2.0 / rhoStratTFC(i, j &
                &+ 1, k))
            fcscal_b = sqrt(pStratTFC(i, j - 1, k) ** 2.0 / rhoStratTFC(i, j &
                &- 1, k))
            fcscal_u = sqrt(pStratTFC(i, j, k + 1) ** 2.0 / rhoStratTFC(i, j, &
                &k + 1))
            fcscal_d = sqrt(pStratTFC(i, j, k - 1) ** 2.0 / rhoStratTFC(i, j, &
                &k - 1))
            fcscal_ru = sqrt(pStratTFC(i + 1, j, k + 1) ** 2.0 / rhoStratTFC(i &
                &+ 1, j, k + 1))
            fcscal_rd = sqrt(pStratTFC(i + 1, j, k - 1) ** 2.0 / rhoStratTFC(i &
                &+ 1, j, k - 1))
            fcscal_lu = sqrt(pStratTFC(i - 1, j, k + 1) ** 2.0 / rhoStratTFC(i &
                &- 1, j, k + 1))
            fcscal_ld = sqrt(pStratTFC(i - 1, j, k - 1) ** 2.0 / rhoStratTFC(i &
                &- 1, j, k - 1))
            fcscal_fu = sqrt(pStratTFC(i, j + 1, k + 1) ** 2.0 &
                &/ rhoStratTFC(i, j + 1, k + 1))
            fcscal_fd = sqrt(pStratTFC(i, j + 1, k - 1) ** 2.0 &
                &/ rhoStratTFC(i, j + 1, k - 1))
            fcscal_bu = sqrt(pStratTFC(i, j - 1, k + 1) ** 2.0 &
                &/ rhoStratTFC(i, j - 1, k + 1))
            fcscal_bd = sqrt(pStratTFC(i, j - 1, k - 1) ** 2.0 &
                &/ rhoStratTFC(i, j - 1, k - 1))
            fcscal_uu = sqrt(pStratTFC(i, j, k + 2) ** 2.0 / rhoStratTFC(i, j, &
                &k + 2))
            fcscal_dd = sqrt(pStratTFC(i, j, k - 2) ** 2.0 / rhoStratTFC(i, j, &
                &k - 2))
            fcscal_ruu = sqrt(pStratTFC(i + 1, j, k + 2) ** 2.0 &
                &/ rhoStratTFC(i + 1, j, k + 2))
            fcscal_rdd = sqrt(pStratTFC(i + 1, j, k - 2) ** 2.0 &
                &/ rhoStratTFC(i + 1, j, k - 2))
            fcscal_luu = sqrt(pStratTFC(i - 1, j, k + 2) ** 2.0 &
                &/ rhoStratTFC(i - 1, j, k + 2))
            fcscal_ldd = sqrt(pStratTFC(i - 1, j, k - 2) ** 2.0 &
                &/ rhoStratTFC(i - 1, j, k - 2))
            fcscal_fuu = sqrt(pStratTFC(i, j + 1, k + 2) ** 2.0 &
                &/ rhoStratTFC(i, j + 1, k + 2))
            fcscal_fdd = sqrt(pStratTFC(i, j + 1, k - 2) ** 2.0 &
                &/ rhoStratTFC(i, j + 1, k - 2))
            fcscal_buu = sqrt(pStratTFC(i, j - 1, k + 2) ** 2.0 &
                &/ rhoStratTFC(i, j - 1, k + 2))
            fcscal_bdd = sqrt(pStratTFC(i, j - 1, k - 2) ** 2.0 &
                &/ rhoStratTFC(i, j - 1, k - 2))

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
            pEdgeUDiv = jac(i, j, k) * jac(i, j, k + 1) * (pStratTFC(i, j, k) &
                &+ pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            pEdgeDDiv = jac(i, j, k) * jac(i, j, k - 1) * (pStratTFC(i, j, k) &
                &+ pStratTFC(i, j, k - 1)) / (jac(i, j, k) + jac(i, j, k - 1))

            ! Compute P coefficients (pressure gradient).
            pEdgeRGra = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i + 1, j, k))
            pEdgeLGra = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i - 1, j, k))
            pEdgeFGra = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j + 1, k))
            pEdgeBGra = 0.5 * (pStratTFC(i, j, k) + pStratTFC(i, j - 1, k))
            pEdgeUGra = (jac(i, j, k + 1) * pStratTFC(i, j, k) + jac(i, j, k) &
                &* pStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1))
            pEdgeDGra = (jac(i, j, k - 1) * pStratTFC(i, j, k) + jac(i, j, k) &
                &* pStratTFC(i, j, k - 1)) / (jac(i, j, k) + jac(i, j, k - 1))
            pUEdgeRGra = 0.5 * (pStratTFC(i, j, k + 1) + pStratTFC(i + 1, j, k &
                &+ 1))
            pUEdgeLGra = 0.5 * (pStratTFC(i, j, k + 1) + pStratTFC(i - 1, j, k &
                &+ 1))
            pUEdgeFGra = 0.5 * (pStratTFC(i, j, k + 1) + pStratTFC(i, j + 1, k &
                &+ 1))
            pUEdgeBGra = 0.5 * (pStratTFC(i, j, k + 1) + pStratTFC(i, j - 1, k &
                &+ 1))
            pDEdgeRGra = 0.5 * (pStratTFC(i, j, k - 1) + pStratTFC(i + 1, j, k &
                &- 1))
            pDEdgeLGra = 0.5 * (pStratTFC(i, j, k - 1) + pStratTFC(i - 1, j, k &
                &- 1))
            pDEdgeFGra = 0.5 * (pStratTFC(i, j, k - 1) + pStratTFC(i, j + 1, k &
                &- 1))
            pDEdgeBGra = 0.5 * (pStratTFC(i, j, k - 1) + pStratTFC(i, j - 1, k &
                &- 1))

            ! Compute density coefficients.
            rhoStratEdgeR = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i + 1, &
                &j, k))
            rhoStratEdgeL = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i - 1, &
                &j, k))
            rhoStratEdgeF = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, j &
                &+ 1, k))
            rhoStratEdgeB = 0.5 * (rhoStratTFC(i, j, k) + rhoStratTFC(i, j &
                &- 1, k))
            rhoStratEdgeU = (jac(i, j, k + 1) * rhoStratTFC(i, j, k) + jac(i, &
                &j, k) * rhoStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, &
                &k + 1))
            rhoStratEdgeD = (jac(i, j, k - 1) * rhoStratTFC(i, j, k) + jac(i, &
                &j, k) * rhoStratTFC(i, j, k - 1)) / (jac(i, j, k) + jac(i, j, &
                &k - 1))
            rhoEdgeR = 0.5 * (var%rho(i, j, k) + var%rho(i + 1, j, k)) &
                &+ rhoStratEdgeR
            rhoEdgeL = 0.5 * (var%rho(i, j, k) + var%rho(i - 1, j, k)) &
                &+ rhoStratEdgeL
            rhoEdgeF = 0.5 * (var%rho(i, j, k) + var%rho(i, j + 1, k)) &
                &+ rhoStratEdgeF
            rhoEdgeB = 0.5 * (var%rho(i, j, k) + var%rho(i, j - 1, k)) &
                &+ rhoStratEdgeB
            rhoEdgeU = (jac(i, j, k + 1) * var%rho(i, j, k) + jac(i, j, k) &
                &* var%rho(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, k + 1)) &
                &+ rhoStratEdgeU
            rhoEdgeD = (jac(i, j, k - 1) * var%rho(i, j, k) + jac(i, j, k) &
                &* var%rho(i, j, k - 1)) / (jac(i, j, k) + jac(i, j, k - 1)) &
                &+ rhoStratEdgeD

            rhoUEdgeR = 0.5 * (var%rho(i, j, k + 1) + var%rho(i + 1, j, k + 1) &
                &+ rhoStratTFC(i, j, k + 1) + rhoStratTFC(i + 1, j, k + 1))
            rhoUEdgeL = 0.5 * (var%rho(i, j, k + 1) + var%rho(i - 1, j, k + 1) &
                &+ rhoStratTFC(i, j, k + 1) + rhoStratTFC(i - 1, j, k + 1))
            rhoUEdgeF = 0.5 * (var%rho(i, j, k + 1) + var%rho(i, j + 1, k + 1) &
                &+ rhoStratTFC(i, j, k + 1) + rhoStratTFC(i, j + 1, k + 1))
            rhoUEdgeB = 0.5 * (var%rho(i, j, k + 1) + var%rho(i, j - 1, k + 1) &
                &+ rhoStratTFC(i, j, k + 1) + rhoStratTFC(i, j - 1, k + 1))
            rhoDEdgeR = 0.5 * (var%rho(i, j, k - 1) + var%rho(i + 1, j, k - 1) &
                &+ rhoStratTFC(i, j, k - 1) + rhoStratTFC(i + 1, j, k - 1))
            rhoDEdgeL = 0.5 * (var%rho(i, j, k - 1) + var%rho(i - 1, j, k - 1) &
                &+ rhoStratTFC(i, j, k - 1) + rhoStratTFC(i - 1, j, k - 1))
            rhoDEdgeF = 0.5 * (var%rho(i, j, k - 1) + var%rho(i, j + 1, k - 1) &
                &+ rhoStratTFC(i, j, k - 1) + rhoStratTFC(i, j + 1, k - 1))
            rhoDEdgeB = 0.5 * (var%rho(i, j, k - 1) + var%rho(i, j - 1, k - 1) &
                &+ rhoStratTFC(i, j, k - 1) + rhoStratTFC(i, j - 1, k - 1))

            ! Compute squared buoyancy frequency at edges.
            bvsStratEdgeU = (jac(i, j, k + 1) * bvsStratTFC(i, j, k) + jac(i, &
                &j, k) * bvsStratTFC(i, j, k + 1)) / (jac(i, j, k) + jac(i, j, &
                &k + 1))
            bvsStratEdgeD = (jac(i, j, k - 1) * bvsStratTFC(i, j, k) + jac(i, &
                &j, k) * bvsStratTFC(i, j, k - 1)) / (jac(i, j, k) + jac(i, j, &
                &k - 1))

            ! Interpolate metric-tensor elements.
            met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))
            met13EdgeL = 0.5 * (met(i, j, k, 1, 3) + met(i - 1, j, k, 1, 3))
            met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3))
            met23EdgeB = 0.5 * (met(i, j, k, 2, 3) + met(i, j - 1, k, 2, 3))
            met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 1, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 2, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) + jac(i, j, k) &
                &* met(i, j, k + 1, 3, 3)) / (jac(i, j, k) + jac(i, j, k + 1))
            met13EdgeD = (jac(i, j, k - 1) * met(i, j, k, 1, 3) + jac(i, j, k) &
                &* met(i, j, k - 1, 1, 3)) / (jac(i, j, k) + jac(i, j, k - 1))
            met23EdgeD = (jac(i, j, k - 1) * met(i, j, k, 2, 3) + jac(i, j, k) &
                &* met(i, j, k - 1, 2, 3)) / (jac(i, j, k) + jac(i, j, k - 1))
            met33EdgeD = (jac(i, j, k - 1) * met(i, j, k, 3, 3) + jac(i, j, k) &
                &* met(i, j, k - 1, 3, 3)) / (jac(i, j, k) + jac(i, j, k - 1))
            met13UEdgeR = 0.5 * (met(i, j, k + 1, 1, 3) + met(i + 1, j, k + 1, &
                &1, 3))
            met13UEdgeL = 0.5 * (met(i, j, k + 1, 1, 3) + met(i - 1, j, k + 1, &
                &1, 3))
            met23UEdgeF = 0.5 * (met(i, j, k + 1, 2, 3) + met(i, j + 1, k + 1, &
                &2, 3))
            met23UEdgeB = 0.5 * (met(i, j, k + 1, 2, 3) + met(i, j - 1, k + 1, &
                &2, 3))
            met13DEdgeR = 0.5 * (met(i, j, k - 1, 1, 3) + met(i + 1, j, k - 1, &
                &1, 3))
            met13DEdgeL = 0.5 * (met(i, j, k - 1, 1, 3) + met(i - 1, j, k - 1, &
                &1, 3))
            met23DEdgeF = 0.5 * (met(i, j, k - 1, 2, 3) + met(i, j + 1, k - 1, &
                &2, 3))
            met23DEdgeB = 0.5 * (met(i, j, k - 1, 2, 3) + met(i, j - 1, k - 1, &
                &2, 3))

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
              facEdgeU = facEdgeU + dt * (jac(i, j, k + 1) * kr_sp_w_tfc(i, j, &
                  &k) + jac(i, j, k) * kr_sp_w_tfc(i, j, k + 1)) / (jac(i, j, &
                  &k) + jac(i, j, k + 1))
              facEdgeD = facEdgeD + dt * (jac(i, j, k - 1) * kr_sp_w_tfc(i, j, &
                  &k) + jac(i, j, k) * kr_sp_w_tfc(i, j, k - 1)) / (jac(i, j, &
                  &k) + jac(i, j, k - 1))
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
                  &* 0.5 * met(i, j, k, 1, 3) * jac(i, j, k + 1) / (jac(i, j, &
                  &k) + jac(i, j, k + 1)) * impHorEdgeR * facEdgeR / rhoEdgeR
            else if(k == nz .and. zBoundary == "solid_wall") then
              gEdgeR = jacInv / dx * pEdgeRDiv * impHorEdgeR * facEdgeR &
                  &/ rhoEdgeR - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                  &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                  &* 0.5 * met(i, j, k, 1, 3) * jac(i, j, k - 1) / (jac(i, j, &
                  &k) + jac(i, j, k - 1)) * impHorEdgeR * facEdgeR / rhoEdgeR
            else
              gEdgeR = jacInv / dx * pEdgeRDiv * impHorEdgeR * facEdgeR &
                  &/ rhoEdgeR + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                  &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                  &* 0.5 * met(i, j, k, 1, 3) * jac(i, j, k + 1) / (jac(i, j, &
                  &k) + jac(i, j, k + 1)) * impHorEdgeR * facEdgeR / rhoEdgeR &
                  &- jacInv / dz * pEdgeDDiv * impVerEdgeD * rhoStratEdgeD &
                  &/ rhoEdgeD * bvsStratEdgeD * dt ** 2.0 * 0.5 * met(i, j, k, &
                  &1, 3) * jac(i, j, k - 1) / (jac(i, j, k) + jac(i, j, k &
                  &- 1)) * impHorEdgeR * facEdgeR / rhoEdgeR
            end if

            ! G(i - 1 / 2)
            if(k == 1 .and. zBoundary == "solid_wall") then
              gEdgeL = - jacInv / dx * pEdgeLDiv * impHorEdgeL * facEdgeL &
                  &/ rhoEdgeL + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                  &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                  &* 0.5 * met(i, j, k, 1, 3) * jac(i, j, k + 1) / (jac(i, j, &
                  &k) + jac(i, j, k + 1)) * impHorEdgeL * facEdgeL / rhoEdgeL
            else if(k == nz .and. zBoundary == "solid_wall") then
              gEdgeL = - jacInv / dx * pEdgeLDiv * impHorEdgeL * facEdgeL &
                  &/ rhoEdgeL - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                  &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                  &* 0.5 * met(i, j, k, 1, 3) * jac(i, j, k - 1) / (jac(i, j, &
                  &k) + jac(i, j, k - 1)) * impHorEdgeL * facEdgeL / rhoEdgeL
            else
              gEdgeL = - jacInv / dx * pEdgeLDiv * impHorEdgeL * facEdgeL &
                  &/ rhoEdgeL + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                  &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                  &* 0.5 * met(i, j, k, 1, 3) * jac(i, j, k + 1) / (jac(i, j, &
                  &k) + jac(i, j, k + 1)) * impHorEdgeL * facEdgeL / rhoEdgeL &
                  &- jacInv / dz * pEdgeDDiv * impVerEdgeD * rhoStratEdgeD &
                  &/ rhoEdgeD * bvsStratEdgeD * dt ** 2.0 * 0.5 * met(i, j, k, &
                  &1, 3) * jac(i, j, k - 1) / (jac(i, j, k) + jac(i, j, k &
                  &- 1)) * impHorEdgeL * facEdgeL / rhoEdgeL
            end if

            ! G(j + 1 / 2)
            if(k == 1 .and. zBoundary == "solid_wall") then
              gEdgeF = jacInv / dy * pEdgeFDiv * impHorEdgeF * facEdgeF &
                  &/ rhoEdgeF + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                  &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                  &* 0.5 * met(i, j, k, 2, 3) * jac(i, j, k + 1) / (jac(i, j, &
                  &k) + jac(i, j, k + 1)) * impHorEdgeF * facEdgeF / rhoEdgeF
            else if(k == nz .and. zBoundary == "solid_wall") then
              gEdgeF = jacInv / dy * pEdgeFDiv * impHorEdgeF * facEdgeF &
                  &/ rhoEdgeF - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                  &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                  &* 0.5 * met(i, j, k, 2, 3) * jac(i, j, k - 1) / (jac(i, j, &
                  &k) + jac(i, j, k - 1)) * impHorEdgeF * facEdgeF / rhoEdgeF
            else
              gEdgeF = jacInv / dy * pEdgeFDiv * impHorEdgeF * facEdgeF &
                  &/ rhoEdgeF + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                  &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                  &* 0.5 * met(i, j, k, 2, 3) * jac(i, j, k + 1) / (jac(i, j, &
                  &k) + jac(i, j, k + 1)) * impHorEdgeF * facEdgeF / rhoEdgeF &
                  &- jacInv / dz * pEdgeDDiv * impVerEdgeD * rhoStratEdgeD &
                  &/ rhoEdgeD * bvsStratEdgeD * dt ** 2.0 * 0.5 * met(i, j, k, &
                  &2, 3) * jac(i, j, k - 1) / (jac(i, j, k) + jac(i, j, k &
                  &- 1)) * impHorEdgeF * facEdgeF / rhoEdgeF
            end if

            ! G(j - 1 / 2)
            if(k == 1 .and. zBoundary == "solid_wall") then
              gEdgeB = - jacInv / dy * pEdgeBDiv * impHorEdgeB * facEdgeB &
                  &/ rhoEdgeB + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                  &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                  &* 0.5 * met(i, j, k, 2, 3) * jac(i, j, k + 1) / (jac(i, j, &
                  &k) + jac(i, j, k + 1)) * impHorEdgeB * facEdgeB / rhoEdgeB
            else if(k == nz .and. zBoundary == "solid_wall") then
              gEdgeB = - jacInv / dy * pEdgeBDiv * impHorEdgeB * facEdgeB &
                  &/ rhoEdgeB - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                  &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                  &* 0.5 * met(i, j, k, 2, 3) * jac(i, j, k - 1) / (jac(i, j, &
                  &k) + jac(i, j, k - 1)) * impHorEdgeB * facEdgeB / rhoEdgeB
            else
              gEdgeB = - jacInv / dy * pEdgeBDiv * impHorEdgeB * facEdgeB &
                  &/ rhoEdgeB + jacInv / dz * pEdgeUDiv * impVerEdgeU &
                  &* rhoStratEdgeU / rhoEdgeU * bvsStratEdgeU * dt ** 2.0 &
                  &* 0.5 * met(i, j, k, 2, 3) * jac(i, j, k + 1) / (jac(i, j, &
                  &k) + jac(i, j, k + 1)) * impHorEdgeB * facEdgeB / rhoEdgeB &
                  &- jacInv / dz * pEdgeDDiv * impVerEdgeD * rhoStratEdgeD &
                  &/ rhoEdgeD * bvsStratEdgeD * dt ** 2.0 * 0.5 * met(i, j, k, &
                  &2, 3) * jac(i, j, k - 1) / (jac(i, j, k) + jac(i, j, k &
                  &- 1)) * impHorEdgeB * facEdgeB / rhoEdgeB
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
              gUEdgeR = jacInv / dz * pEdgeUDiv * impVerEdgeU * rhoStratEdgeU &
                  &/ rhoEdgeU * bvsStratEdgeU * dt ** 2.0 * 0.5 * met(i, j, k &
                  &+ 1, 1, 3) * jac(i, j, k) / (jac(i, j, k) + jac(i, j, k &
                  &+ 1)) * impHorUEdgeR * facUEdgeR / rhoUEdgeR
            end if

            ! G(i - 1 / 2, k + 1)
            if(k == nz .and. zBoundary == "solid_wall") then
              gUEdgeL = 0.0
            else
              gUEdgeL = jacInv / dz * pEdgeUDiv * impVerEdgeU * rhoStratEdgeU &
                  &/ rhoEdgeU * bvsStratEdgeU * dt ** 2.0 * 0.5 * met(i, j, k &
                  &+ 1, 1, 3) * jac(i, j, k) / (jac(i, j, k) + jac(i, j, k &
                  &+ 1)) * impHorUEdgeL * facUEdgeL / rhoUEdgeL
            end if

            ! G(j + 1 / 2, k + 1)
            if(k == nz .and. zBoundary == "solid_wall") then
              gUEdgeF = 0.0
            else
              gUEdgeF = jacInv / dz * pEdgeUDiv * impVerEdgeU * rhoStratEdgeU &
                  &/ rhoEdgeU * bvsStratEdgeU * dt ** 2.0 * 0.5 * met(i, j, k &
                  &+ 1, 2, 3) * jac(i, j, k) / (jac(i, j, k) + jac(i, j, k &
                  &+ 1)) * impHorUEdgeF * facUEdgeF / rhoUEdgeF
            end if

            ! G(j - 1 / 2, k + 1)
            if(k == nz .and. zBoundary == "solid_wall") then
              gUEdgeB = 0.0
            else
              gUEdgeB = jacInv / dz * pEdgeUDiv * impVerEdgeU * rhoStratEdgeU &
                  &/ rhoEdgeU * bvsStratEdgeU * dt ** 2.0 * 0.5 * met(i, j, k &
                  &+ 1, 2, 3) * jac(i, j, k) / (jac(i, j, k) + jac(i, j, k &
                  &+ 1)) * impHorUEdgeB * facUEdgeB / rhoUEdgeB
            end if

            ! G(i + 1 / 2, k - 1)
            if(k == 1 .and. zBoundary == "solid_wall") then
              gDEdgeR = 0.0
            else
              gDEdgeR = - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                  &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                  &* 0.5 * met(i, j, k - 1, 1, 3) * jac(i, j, k) / (jac(i, j, &
                  &k) + jac(i, j, k - 1)) * impHorDEdgeR * facDEdgeR / rhoDEdgeR
            end if

            ! G(i - 1 / 2, k - 1)
            if(k == 1 .and. zBoundary == "solid_wall") then
              gDEdgeL = 0.0
            else
              gDEdgeL = - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                  &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                  &* 0.5 * met(i, j, k - 1, 1, 3) * jac(i, j, k) / (jac(i, j, &
                  &k) + jac(i, j, k - 1)) * impHorDEdgeL * facDEdgeL / rhoDEdgeL
            end if

            ! G(j + 1 / 2, k - 1)
            if(k == 1 .and. zBoundary == "solid_wall") then
              gDEdgeF = 0.0
            else
              gDEdgeF = - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                  &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                  &* 0.5 * met(i, j, k - 1, 2, 3) * jac(i, j, k) / (jac(i, j, &
                  &k) + jac(i, j, k - 1)) * impHorDEdgeF * facDEdgeF / rhoDEdgeF
            end if

            ! G(j - 1 / 2, k - 1)
            if(k == 1 .and. zBoundary == "solid_wall") then
              gDEdgeB = 0.0
            else
              gDEdgeB = - jacInv / dz * pEdgeDDiv * impVerEdgeD &
                  &* rhoStratEdgeD / rhoEdgeD * bvsStratEdgeD * dt ** 2.0 &
                  &* 0.5 * met(i, j, k - 1, 2, 3) * jac(i, j, k) / (jac(i, j, &
                  &k) + jac(i, j, k - 1)) * impHorDEdgeB * facDEdgeB / rhoDEdgeB
            end if

            ! Compute tensor elements

            ! ------------------- A(i,j,k) --------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              AC = - gEdgeR * pEdgeRGra * (1.0 / dx + 0.75 * met13EdgeR / dz) &
                  &+ gEdgeL * pEdgeLGra * (1.0 / dx - 0.75 * met13EdgeL / dz) &
                  &- gEdgeF * pEdgeFGra * (1.0 / dy + 0.75 * met23EdgeF / dz) &
                  &+ gEdgeB * pEdgeBGra * (1.0 / dy - 0.75 * met23EdgeB / dz) &
                  &- gEdgeU * pEdgeUGra * met33EdgeU / dz - gUEdgeR &
                  &* pUEdgeRGra * 0.25 * met13UEdgeR / dz - gUEdgeL &
                  &* pUEdgeLGra * 0.25 * met13UEdgeL / dz - gUEdgeF &
                  &* pUEdgeFGra * 0.25 * met23UEdgeF / dz - gUEdgeB &
                  &* pUEdgeBGra * 0.25 * met23UEdgeB / dz
            else if(k == 2 .and. zBoundary == "solid_wall") then
              AC = - gEdgeR * pEdgeRGra / dx + gEdgeL * pEdgeLGra / dx &
                  &- gEdgeF * pEdgeFGra / dy + gEdgeB * pEdgeBGra / dy &
                  &- gEdgeU * pEdgeUGra * met33EdgeU / dz + gEdgeD * pEdgeDGra &
                  &* met33EdgeD / dz - gUEdgeR * pUEdgeRGra * 0.25 &
                  &* met13UEdgeR / dz - gUEdgeL * pUEdgeLGra * 0.25 &
                  &* met13UEdgeL / dz - gUEdgeF * pUEdgeFGra * 0.25 &
                  &* met23UEdgeF / dz - gUEdgeB * pUEdgeBGra * 0.25 &
                  &* met23UEdgeB / dz + gDEdgeR * pDEdgeRGra * met13DEdgeR &
                  &/ dz + gDEdgeL * pDEdgeLGra * met13DEdgeL / dz + gDEdgeF &
                  &* pDEdgeFGra * met23DEdgeF / dz + gDEdgeB * pDEdgeBGra &
                  &* met23DEdgeB / dz
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              AC = - gEdgeR * pEdgeRGra / dx + gEdgeL * pEdgeLGra / dx &
                  &- gEdgeF * pEdgeFGra / dy + gEdgeB * pEdgeBGra / dy &
                  &- gEdgeU * pEdgeUGra * met33EdgeU / dz + gEdgeD * pEdgeDGra &
                  &* met33EdgeD / dz - gUEdgeR * pUEdgeRGra * met13UEdgeR / dz &
                  &- gUEdgeL * pUEdgeLGra * met13UEdgeL / dz - gUEdgeF &
                  &* pUEdgeFGra * met23UEdgeF / dz - gUEdgeB * pUEdgeBGra &
                  &* met23UEdgeB / dz + gDEdgeR * pDEdgeRGra * 0.25 &
                  &* met13DEdgeR / dz + gDEdgeL * pDEdgeLGra * 0.25 &
                  &* met13DEdgeL / dz + gDEdgeF * pDEdgeFGra * 0.25 &
                  &* met23DEdgeF / dz + gDEdgeB * pDEdgeBGra * 0.25 &
                  &* met23DEdgeB / dz
            else if(k == nz .and. zBoundary == "solid_wall") then
              AC = - gEdgeR * pEdgeRGra * (1.0 / dx - 0.75 * met13EdgeR / dz) &
                  &+ gEdgeL * pEdgeLGra * (1.0 / dx + 0.75 * met13EdgeL / dz) &
                  &- gEdgeF * pEdgeFGra * (1.0 / dy - 0.75 * met23EdgeF / dz) &
                  &+ gEdgeB * pEdgeBGra * (1.0 / dy + 0.75 * met23EdgeB / dz) &
                  &+ gEdgeD * pEdgeDGra * met33EdgeD / dz + gDEdgeR &
                  &* pDEdgeRGra * 0.25 * met13DEdgeR / dz + gDEdgeL &
                  &* pDEdgeLGra * 0.25 * met13DEdgeL / dz + gDEdgeF &
                  &* pDEdgeFGra * 0.25 * met23DEdgeF / dz + gDEdgeB &
                  &* pDEdgeBGra * 0.25 * met23DEdgeB / dz
            else
              AC = - gEdgeR * pEdgeRGra / dx + gEdgeL * pEdgeLGra / dx &
                  &- gEdgeF * pEdgeFGra / dy + gEdgeB * pEdgeBGra / dy &
                  &- gEdgeU * pEdgeUGra * met33EdgeU / dz + gEdgeD * pEdgeDGra &
                  &* met33EdgeD / dz - gUEdgeR * pUEdgeRGra * 0.25 &
                  &* met13UEdgeR / dz - gUEdgeL * pUEdgeLGra * 0.25 &
                  &* met13UEdgeL / dz - gUEdgeF * pUEdgeFGra * 0.25 &
                  &* met23UEdgeF / dz - gUEdgeB * pUEdgeBGra * 0.25 &
                  &* met23UEdgeB / dz + gDEdgeR * pDEdgeRGra * 0.25 &
                  &* met13DEdgeR / dz + gDEdgeL * pDEdgeLGra * 0.25 &
                  &* met13DEdgeL / dz + gDEdgeF * pDEdgeFGra * 0.25 &
                  &* met23DEdgeF / dz + gDEdgeB * pDEdgeBGra * 0.25 &
                  &* met23DEdgeB / dz
            end if

            ! ------------------ A(i+1,j,k) -------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              AR = gEdgeR * pEdgeRGra * (1.0 / dx - 0.75 * met13EdgeR / dz) &
                  &+ gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx * jac(i + 1, &
                  &j, k + 1) / (jac(i + 1, j, k) + jac(i + 1, j, k + 1)) &
                  &- gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz
            else if(k == 2 .and. zBoundary == "solid_wall") then
              AR = gEdgeR * pEdgeRGra / dx + gEdgeU * pEdgeUGra * 0.5 &
                  &* met13EdgeU / dx * jac(i + 1, j, k + 1) / (jac(i + 1, j, &
                  &k) + jac(i + 1, j, k + 1)) + gEdgeD * pEdgeDGra * 0.5 &
                  &* met13EdgeD / dx * jac(i + 1, j, k - 1) / (jac(i + 1, j, &
                  &k) + jac(i + 1, j, k - 1)) - gUEdgeR * pUEdgeRGra * 0.25 &
                  &* met13UEdgeR / dz + gDEdgeR * pDEdgeRGra * met13DEdgeR / dz
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              AR = gEdgeR * pEdgeRGra / dx + gEdgeU * pEdgeUGra * 0.5 &
                  &* met13EdgeU / dx * jac(i + 1, j, k + 1) / (jac(i + 1, j, &
                  &k) + jac(i + 1, j, k + 1)) + gEdgeD * pEdgeDGra * 0.5 &
                  &* met13EdgeD / dx * jac(i + 1, j, k - 1) / (jac(i + 1, j, &
                  &k) + jac(i + 1, j, k - 1)) - gUEdgeR * pUEdgeRGra &
                  &* met13UEdgeR / dz + gDEdgeR * pDEdgeRGra * 0.25 &
                  &* met13DEdgeR / dz
            else if(k == nz .and. zBoundary == "solid_wall") then
              AR = gEdgeR * pEdgeRGra * (1.0 / dx + 0.75 * met13EdgeR / dz) &
                  &+ gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx * jac(i + 1, &
                  &j, k - 1) / (jac(i + 1, j, k) + jac(i + 1, j, k - 1)) &
                  &+ gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz
            else
              AR = gEdgeR * pEdgeRGra / dx + gEdgeU * pEdgeUGra * 0.5 &
                  &* met13EdgeU / dx * jac(i + 1, j, k + 1) / (jac(i + 1, j, &
                  &k) + jac(i + 1, j, k + 1)) + gEdgeD * pEdgeDGra * 0.5 &
                  &* met13EdgeD / dx * jac(i + 1, j, k - 1) / (jac(i + 1, j, &
                  &k) + jac(i + 1, j, k - 1)) - gUEdgeR * pUEdgeRGra * 0.25 &
                  &* met13UEdgeR / dz + gDEdgeR * pDEdgeRGra * 0.25 &
                  &* met13DEdgeR / dz
            end if

            ! ------------------ A(i-1,j,k) -------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              AL = - gEdgeL * pEdgeLGra * (1.0 / dx + 0.75 * met13EdgeL / dz) &
                  &- gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx * jac(i - 1, &
                  &j, k + 1) / (jac(i - 1, j, k) + jac(i - 1, j, k + 1)) &
                  &- gUEdgeL * pUEdgeLGra * 0.25 * met13UEdgeL / dz
            else if(k == 2 .and. zBoundary == "solid_wall") then
              AL = - gEdgeL * pEdgeLGra / dx - gEdgeU * pEdgeUGra * 0.5 &
                  &* met13EdgeU / dx * jac(i - 1, j, k + 1) / (jac(i - 1, j, &
                  &k) + jac(i - 1, j, k + 1)) - gEdgeD * pEdgeDGra * 0.5 &
                  &* met13EdgeD / dx * jac(i - 1, j, k - 1) / (jac(i - 1, j, &
                  &k) + jac(i - 1, j, k - 1)) - gUEdgeL * pUEdgeLGra * 0.25 &
                  &* met13UEdgeL / dz + gDEdgeL * pDEdgeLGra * met13DEdgeL / dz
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              AL = - gEdgeL * pEdgeLGra / dx - gEdgeU * pEdgeUGra * 0.5 &
                  &* met13EdgeU / dx * jac(i - 1, j, k + 1) / (jac(i - 1, j, &
                  &k) + jac(i - 1, j, k + 1)) - gEdgeD * pEdgeDGra * 0.5 &
                  &* met13EdgeD / dx * jac(i - 1, j, k - 1) / (jac(i - 1, j, &
                  &k) + jac(i - 1, j, k - 1)) - gUEdgeL * pUEdgeLGra &
                  &* met13UEdgeL / dz + gDEdgeL * pDEdgeLGra * 0.25 &
                  &* met13DEdgeL / dz
            else if(k == nz .and. zBoundary == "solid_wall") then
              AL = - gEdgeL * pEdgeLGra * (1.0 / dx - 0.75 * met13EdgeL / dz) &
                  &- gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx * jac(i - 1, &
                  &j, k - 1) / (jac(i - 1, j, k) + jac(i - 1, j, k - 1)) &
                  &+ gDEdgeL * pDEdgeLGra * 0.25 * met13DEdgeL / dz
            else
              AL = - gEdgeL * pEdgeLGra / dx - gEdgeU * pEdgeUGra * 0.5 &
                  &* met13EdgeU / dx * jac(i - 1, j, k + 1) / (jac(i - 1, j, &
                  &k) + jac(i - 1, j, k + 1)) - gEdgeD * pEdgeDGra * 0.5 &
                  &* met13EdgeD / dx * jac(i - 1, j, k - 1) / (jac(i - 1, j, &
                  &k) + jac(i - 1, j, k - 1)) - gUEdgeL * pUEdgeLGra * 0.25 &
                  &* met13UEdgeL / dz + gDEdgeL * pDEdgeLGra * 0.25 &
                  &* met13DEdgeL / dz
            end if

            ! ------------------ A(i,j+1,k) -------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              AF = gEdgeF * pEdgeFGra * (1.0 / dy - 0.75 * met23EdgeF / dz) &
                  &+ gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy * jac(i, j &
                  &+ 1, k + 1) / (jac(i, j + 1, k) + jac(i, j + 1, k + 1)) &
                  &- gUEdgeF * pUEdgeFGra * 0.25 * met23UEdgeF / dz
            else if(k == 2 .and. zBoundary == "solid_wall") then
              AF = gEdgeF * pEdgeFGra / dy + gEdgeU * pEdgeUGra * 0.5 &
                  &* met23EdgeU / dy * jac(i, j + 1, k + 1) / (jac(i, j + 1, &
                  &k) + jac(i, j + 1, k + 1)) + gEdgeD * pEdgeDGra * 0.5 &
                  &* met23EdgeD / dy * jac(i, j + 1, k - 1) / (jac(i, j + 1, &
                  &k) + jac(i, j + 1, k - 1)) - gUEdgeF * pUEdgeFGra * 0.25 &
                  &* met23UEdgeF / dz + gDEdgeF * pDEdgeFGra * met23DEdgeF / dz
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              AF = gEdgeF * pEdgeFGra / dy + gEdgeU * pEdgeUGra * 0.5 &
                  &* met23EdgeU / dy * jac(i, j + 1, k + 1) / (jac(i, j + 1, &
                  &k) + jac(i, j + 1, k + 1)) + gEdgeD * pEdgeDGra * 0.5 &
                  &* met23EdgeD / dy * jac(i, j + 1, k - 1) / (jac(i, j + 1, &
                  &k) + jac(i, j + 1, k - 1)) - gUEdgeF * pUEdgeFGra &
                  &* met23UEdgeF / dz + gDEdgeF * pDEdgeFGra * 0.25 &
                  &* met23DEdgeF / dz
            else if(k == nz .and. zBoundary == "solid_wall") then
              AF = gEdgeF * pEdgeFGra * (1.0 / dy + 0.75 * met23EdgeF / dz) &
                  &+ gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy * jac(i, j &
                  &+ 1, k - 1) / (jac(i, j + 1, k) + jac(i, j + 1, k - 1)) &
                  &+ gDEdgeF * pDEdgeFGra * 0.25 * met23DEdgeF / dz
            else
              AF = gEdgeF * pEdgeFGra / dy + gEdgeU * pEdgeUGra * 0.5 &
                  &* met23EdgeU / dy * jac(i, j + 1, k + 1) / (jac(i, j + 1, &
                  &k) + jac(i, j + 1, k + 1)) + gEdgeD * pEdgeDGra * 0.5 &
                  &* met23EdgeD / dy * jac(i, j + 1, k - 1) / (jac(i, j + 1, &
                  &k) + jac(i, j + 1, k - 1)) - gUEdgeF * pUEdgeFGra * 0.25 &
                  &* met23UEdgeF / dz + gDEdgeF * pDEdgeFGra * 0.25 &
                  &* met23DEdgeF / dz
            end if

            ! ------------------ A(i,j-1,k) -------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              AB = - gEdgeB * pEdgeBGra * (1.0 / dy + 0.75 * met23EdgeB / dz) &
                  &- gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy * jac(i, j &
                  &- 1, k + 1) / (jac(i, j - 1, k) + jac(i, j - 1, k + 1)) &
                  &- gUEdgeB * pUEdgeBGra * 0.25 * met23UEdgeB / dz
            else if(k == 2 .and. zBoundary == "solid_wall") then
              AB = - gEdgeB * pEdgeBGra / dy - gEdgeU * pEdgeUGra * 0.5 &
                  &* met23EdgeU / dy * jac(i, j - 1, k + 1) / (jac(i, j - 1, &
                  &k) + jac(i, j - 1, k + 1)) - gEdgeD * pEdgeDGra * 0.5 &
                  &* met23EdgeD / dy * jac(i, j - 1, k - 1) / (jac(i, j - 1, &
                  &k) + jac(i, j - 1, k - 1)) - gUEdgeB * pUEdgeBGra * 0.25 &
                  &* met23UEdgeB / dz + gDEdgeB * pDEdgeBGra * met23DEdgeB / dz
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              AB = - gEdgeB * pEdgeBGra / dy - gEdgeU * pEdgeUGra * 0.5 &
                  &* met23EdgeU / dy * jac(i, j - 1, k + 1) / (jac(i, j - 1, &
                  &k) + jac(i, j - 1, k + 1)) - gEdgeD * pEdgeDGra * 0.5 &
                  &* met23EdgeD / dy * jac(i, j - 1, k - 1) / (jac(i, j - 1, &
                  &k) + jac(i, j - 1, k - 1)) - gUEdgeB * pUEdgeBGra &
                  &* met23UEdgeB / dz + gDEdgeB * pDEdgeBGra * 0.25 &
                  &* met23DEdgeB / dz
            else if(k == nz .and. zBoundary == "solid_wall") then
              AB = - gEdgeB * pEdgeBGra * (1.0 / dy - 0.75 * met23EdgeB / dz) &
                  &- gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy * jac(i, j &
                  &- 1, k - 1) / (jac(i, j - 1, k) + jac(i, j - 1, k - 1)) &
                  &+ gDEdgeB * pDEdgeBGra * 0.25 * met23DEdgeB / dz
            else
              AB = - gEdgeB * pEdgeBGra / dy - gEdgeU * pEdgeUGra * 0.5 &
                  &* met23EdgeU / dy * jac(i, j - 1, k + 1) / (jac(i, j - 1, &
                  &k) + jac(i, j - 1, k + 1)) - gEdgeD * pEdgeDGra * 0.5 &
                  &* met23EdgeD / dy * jac(i, j - 1, k - 1) / (jac(i, j - 1, &
                  &k) + jac(i, j - 1, k - 1)) - gUEdgeB * pUEdgeBGra * 0.25 &
                  &* met23UEdgeB / dz + gDEdgeB * pDEdgeBGra * 0.25 &
                  &* met23DEdgeB / dz
            end if

            ! ------------------ A(i,j,k+1) -------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              AU = gEdgeR * pEdgeRGra * met13EdgeR / dz + gEdgeL * pEdgeLGra &
                  &* met13EdgeL / dz + gEdgeF * pEdgeFGra * met23EdgeF / dz &
                  &+ gEdgeB * pEdgeBGra * met23EdgeB / dz + gEdgeU * pEdgeUGra &
                  &* met33EdgeU / dz - gUEdgeR * pUEdgeRGra / dx + gUEdgeL &
                  &* pUEdgeLGra / dx - gUEdgeF * pUEdgeFGra / dy + gUEdgeB &
                  &* pUEdgeBGra / dy
            else if(k == 2 .and. zBoundary == "solid_wall") then
              AU = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz + gEdgeL &
                  &* pEdgeLGra * 0.25 * met13EdgeL / dz + gEdgeF * pEdgeFGra &
                  &* 0.25 * met23EdgeF / dz + gEdgeB * pEdgeBGra * 0.25 &
                  &* met23EdgeB / dz + gEdgeU * pEdgeUGra * met33EdgeU / dz &
                  &- gUEdgeR * pUEdgeRGra / dx + gUEdgeL * pUEdgeLGra / dx &
                  &- gUEdgeF * pUEdgeFGra / dy + gUEdgeB * pUEdgeBGra / dy &
                  &- gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz - gDEdgeL &
                  &* pDEdgeLGra * 0.25 * met13DEdgeL / dz - gDEdgeF &
                  &* pDEdgeFGra * 0.25 * met23DEdgeF / dz - gDEdgeB &
                  &* pDEdgeBGra * 0.25 * met23DEdgeB / dz
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              AU = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz + gEdgeL &
                  &* pEdgeLGra * 0.25 * met13EdgeL / dz + gEdgeF * pEdgeFGra &
                  &* 0.25 * met23EdgeF / dz + gEdgeB * pEdgeBGra * 0.25 &
                  &* met23EdgeB / dz + gEdgeU * pEdgeUGra * met33EdgeU / dz &
                  &- gUEdgeR * pUEdgeRGra * (1.0 / dx - 0.75 * met13UEdgeR &
                  &/ dz) + gUEdgeL * pUEdgeLGra * (1.0 / dx + 0.75 &
                  &* met13UEdgeL / dz) - gUEdgeF * pUEdgeFGra * (1.0 / dy &
                  &- 0.75 * met23UEdgeF / dz) + gUEdgeB * pUEdgeBGra * (1.0 &
                  &/ dy + 0.75 * met23UEdgeB / dz)
            else if(k == nz .and. zBoundary == "solid_wall") then
              AU = 0.0
            else
              AU = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz + gEdgeL &
                  &* pEdgeLGra * 0.25 * met13EdgeL / dz + gEdgeF * pEdgeFGra &
                  &* 0.25 * met23EdgeF / dz + gEdgeB * pEdgeBGra * 0.25 &
                  &* met23EdgeB / dz + gEdgeU * pEdgeUGra * met33EdgeU / dz &
                  &- gUEdgeR * pUEdgeRGra / dx + gUEdgeL * pUEdgeLGra / dx &
                  &- gUEdgeF * pUEdgeFGra / dy + gUEdgeB * pUEdgeBGra / dy
            end if

            ! ------------------ A(i,j,k-1) -------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              AD = 0.0
            else if(k == 2 .and. zBoundary == "solid_wall") then
              AD = - gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz - gEdgeL &
                  &* pEdgeLGra * 0.25 * met13EdgeL / dz - gEdgeF * pEdgeFGra &
                  &* 0.25 * met23EdgeF / dz - gEdgeB * pEdgeBGra * 0.25 &
                  &* met23EdgeB / dz - gEdgeD * pEdgeDGra * met33EdgeD / dz &
                  &- gDEdgeR * pDEdgeRGra * (1.0 / dx + 0.75 * met13DEdgeR &
                  &/ dz) + gDEdgeL * pDEdgeLGra * (1.0 / dx - 0.75 &
                  &* met13DEdgeL / dz) - gDEdgeF * pDEdgeFGra * (1.0 / dy &
                  &+ 0.75 * met23DEdgeF / dz) + gDEdgeB * pDEdgeBGra * (1.0 &
                  &/ dy - 0.75 * met23DEdgeB / dz)
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              AD = - gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz - gEdgeL &
                  &* pEdgeLGra * 0.25 * met13EdgeL / dz - gEdgeF * pEdgeFGra &
                  &* 0.25 * met23EdgeF / dz - gEdgeB * pEdgeBGra * 0.25 &
                  &* met23EdgeB / dz - gEdgeD * pEdgeDGra * met33EdgeD / dz &
                  &- gDEdgeR * pDEdgeRGra / dx + gDEdgeL * pDEdgeLGra / dx &
                  &- gDEdgeF * pDEdgeFGra / dy + gDEdgeB * pDEdgeBGra / dy &
                  &+ gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz + gUEdgeL &
                  &* pUEdgeLGra * 0.25 * met13UEdgeL / dz + gUEdgeF &
                  &* pUEdgeFGra * 0.25 * met23UEdgeF / dz + gUEdgeB &
                  &* pUEdgeBGra * 0.25 * met23UEdgeB / dz
            else if(k == nz .and. zBoundary == "solid_wall") then
              AD = - gEdgeR * pEdgeRGra * met13EdgeR / dz - gEdgeL * pEdgeLGra &
                  &* met13EdgeL / dz - gEdgeF * pEdgeFGra * met23EdgeF / dz &
                  &- gEdgeB * pEdgeBGra * met23EdgeB / dz - gEdgeD * pEdgeDGra &
                  &* met33EdgeD / dz - gDEdgeR * pDEdgeRGra / dx + gDEdgeL &
                  &* pDEdgeLGra / dx - gDEdgeF * pDEdgeFGra / dy + gDEdgeB &
                  &* pDEdgeBGra / dy
            else
              AD = - gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz - gEdgeL &
                  &* pEdgeLGra * 0.25 * met13EdgeL / dz - gEdgeF * pEdgeFGra &
                  &* 0.25 * met23EdgeF / dz - gEdgeB * pEdgeBGra * 0.25 &
                  &* met23EdgeB / dz - gEdgeD * pEdgeDGra * met33EdgeD / dz &
                  &- gDEdgeR * pDEdgeRGra / dx + gDEdgeL * pDEdgeLGra / dx &
                  &- gDEdgeF * pDEdgeFGra / dy + gDEdgeB * pDEdgeBGra / dy
            end if

            ! ----------------- A(i+1,j,k+1) ------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              ARU = gEdgeR * pEdgeRGra * met13EdgeR / dz + gEdgeU * pEdgeUGra &
                  &* 0.5 * met13EdgeU / dx * jac(i + 1, j, k) / (jac(i + 1, j, &
                  &k) + jac(i + 1, j, k + 1)) + gUEdgeR * pUEdgeRGra / dx
            else if(k == 2 .and. zBoundary == "solid_wall") then
              ARU = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz + gEdgeU &
                  &* pEdgeUGra * 0.5 * met13EdgeU / dx * jac(i + 1, j, k) &
                  &/ (jac(i + 1, j, k) + jac(i + 1, j, k + 1)) + gUEdgeR &
                  &* pUEdgeRGra / dx - gDEdgeR * pDEdgeRGra * 0.25 &
                  &* met13DEdgeR / dz
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              ARU = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz + gEdgeU &
                  &* pEdgeUGra * 0.5 * met13EdgeU / dx * jac(i + 1, j, k) &
                  &/ (jac(i + 1, j, k) + jac(i + 1, j, k + 1)) + gUEdgeR &
                  &* pUEdgeRGra * (1.0 / dx + 0.75 * met13UEdgeR / dz)
            else if(k == nz .and. zBoundary == "solid_wall") then
              ARU = 0.0
            else
              ARU = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz + gEdgeU &
                  &* pEdgeUGra * 0.5 * met13EdgeU / dx * jac(i + 1, j, k) &
                  &/ (jac(i + 1, j, k) + jac(i + 1, j, k + 1)) + gUEdgeR &
                  &* pUEdgeRGra / dx
            end if

            ! ----------------- A(i+1,j,k-1) ------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              ARD = 0.0
            else if(k == 2 .and. zBoundary == "solid_wall") then
              ARD = - gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz + gEdgeD &
                  &* pEdgeDGra * 0.5 * met13EdgeD / dx * jac(i + 1, j, k) &
                  &/ (jac(i + 1, j, k) + jac(i + 1, j, k - 1)) + gDEdgeR &
                  &* pDEdgeRGra * (1.0 / dx - 0.75 * met13DEdgeR / dz)
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              ARD = - gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz + gEdgeD &
                  &* pEdgeDGra * 0.5 * met13EdgeD / dx * jac(i + 1, j, k) &
                  &/ (jac(i + 1, j, k) + jac(i + 1, j, k - 1)) + gDEdgeR &
                  &* pDEdgeRGra / dx + gUEdgeR * pUEdgeRGra * 0.25 &
                  &* met13UEdgeR / dz
            else if(k == nz .and. zBoundary == "solid_wall") then
              ARD = - gEdgeR * pEdgeRGra * met13EdgeR / dz + gEdgeD &
                  &* pEdgeDGra * 0.5 * met13EdgeD / dx * jac(i + 1, j, k) &
                  &/ (jac(i + 1, j, k) + jac(i + 1, j, k - 1)) + gDEdgeR &
                  &* pDEdgeRGra / dx
            else
              ARD = - gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz + gEdgeD &
                  &* pEdgeDGra * 0.5 * met13EdgeD / dx * jac(i + 1, j, k) &
                  &/ (jac(i + 1, j, k) + jac(i + 1, j, k - 1)) + gDEdgeR &
                  &* pDEdgeRGra / dx
            end if

            ! ----------------- A(i-1,j,k+1) ------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              ALU = gEdgeL * pEdgeLGra * met13EdgeL / dz - gEdgeU * pEdgeUGra &
                  &* 0.5 * met13EdgeU / dx * jac(i - 1, j, k) / (jac(i - 1, j, &
                  &k) + jac(i - 1, j, k + 1)) - gUEdgeL * pUEdgeLGra / dx
            else if(k == 2 .and. zBoundary == "solid_wall") then
              ALU = gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz - gEdgeU &
                  &* pEdgeUGra * 0.5 * met13EdgeU / dx * jac(i - 1, j, k) &
                  &/ (jac(i - 1, j, k) + jac(i - 1, j, k + 1)) - gUEdgeL &
                  &* pUEdgeLGra / dx - gDEdgeL * pDEdgeLGra * 0.25 &
                  &* met13DEdgeL / dz
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              ALU = gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz - gEdgeU &
                  &* pEdgeUGra * 0.5 * met13EdgeU / dx * jac(i - 1, j, k) &
                  &/ (jac(i - 1, j, k) + jac(i - 1, j, k + 1)) - gUEdgeL &
                  &* pUEdgeLGra * (1.0 / dx - 0.75 * met13UEdgeL / dz)
            else if(k == nz .and. zBoundary == "solid_wall") then
              ALU = 0.0
            else
              ALU = gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz - gEdgeU &
                  &* pEdgeUGra * 0.5 * met13EdgeU / dx * jac(i - 1, j, k) &
                  &/ (jac(i - 1, j, k) + jac(i - 1, j, k + 1)) - gUEdgeL &
                  &* pUEdgeLGra / dx
            end if

            ! ----------------- A(i-1,j,k-1) ------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              ALD = 0.0
            else if(k == 2 .and. zBoundary == "solid_wall") then
              ALD = - gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz - gEdgeD &
                  &* pEdgeDGra * 0.5 * met13EdgeD / dx * jac(i - 1, j, k) &
                  &/ (jac(i - 1, j, k) + jac(i - 1, j, k - 1)) - gDEdgeL &
                  &* pDEdgeLGra * (1.0 / dx + 0.75 * met13DEdgeL / dz)
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              ALD = - gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz - gEdgeD &
                  &* pEdgeDGra * 0.5 * met13EdgeD / dx * jac(i - 1, j, k) &
                  &/ (jac(i - 1, j, k) + jac(i - 1, j, k - 1)) - gDEdgeL &
                  &* pDEdgeLGra / dx + gUEdgeL * pUEdgeLGra * 0.25 &
                  &* met13UEdgeL / dz
            else if(k == nz .and. zBoundary == "solid_wall") then
              ALD = - gEdgeL * pEdgeLGra * met13EdgeL / dz - gEdgeD &
                  &* pEdgeDGra * 0.5 * met13EdgeD / dx * jac(i - 1, j, k) &
                  &/ (jac(i - 1, j, k) + jac(i - 1, j, k - 1)) - gDEdgeL &
                  &* pDEdgeLGra / dx
            else
              ALD = - gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz - gEdgeD &
                  &* pEdgeDGra * 0.5 * met13EdgeD / dx * jac(i - 1, j, k) &
                  &/ (jac(i - 1, j, k) + jac(i - 1, j, k - 1)) - gDEdgeL &
                  &* pDEdgeLGra / dx
            end if

            ! ----------------- A(i,j+1,k+1) ------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              AFU = gEdgeF * pEdgeFGra * met23EdgeF / dz + gEdgeU * pEdgeUGra &
                  &* 0.5 * met23EdgeU / dy * jac(i, j + 1, k) / (jac(i, j + 1, &
                  &k) + jac(i, j + 1, k + 1)) + gUEdgeF * pUEdgeFGra / dy
            else if(k == 2 .and. zBoundary == "solid_wall") then
              AFU = gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz + gEdgeU &
                  &* pEdgeUGra * 0.5 * met23EdgeU / dy * jac(i, j + 1, k) &
                  &/ (jac(i, j + 1, k) + jac(i, j + 1, k + 1)) + gUEdgeF &
                  &* pUEdgeFGra / dy - gDEdgeF * pDEdgeFGra * 0.25 &
                  &* met23DEdgeF / dz
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              AFU = gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz + gEdgeU &
                  &* pEdgeUGra * 0.5 * met23EdgeU / dy * jac(i, j + 1, k) &
                  &/ (jac(i, j + 1, k) + jac(i, j + 1, k + 1)) + gUEdgeF &
                  &* pUEdgeFGra * (1.0 / dy + 0.75 * met23UEdgeF / dz)
            else if(k == nz .and. zBoundary == "solid_wall") then
              AFU = 0.0
            else
              AFU = gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz + gEdgeU &
                  &* pEdgeUGra * 0.5 * met23EdgeU / dy * jac(i, j + 1, k) &
                  &/ (jac(i, j + 1, k) + jac(i, j + 1, k + 1)) + gUEdgeF &
                  &* pUEdgeFGra / dy
            end if

            ! ----------------- A(i,j+1,k-1) ------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              AFD = 0.0
            else if(k == 2 .and. zBoundary == "solid_wall") then
              AFD = - gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz + gEdgeD &
                  &* pEdgeDGra * 0.5 * met23EdgeD / dy * jac(i, j + 1, k) &
                  &/ (jac(i, j + 1, k) + jac(i, j + 1, k - 1)) + gDEdgeF &
                  &* pDEdgeFGra * (1.0 / dy - 0.75 * met23DEdgeF / dz)
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              AFD = - gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz + gEdgeD &
                  &* pEdgeDGra * 0.5 * met23EdgeD / dy * jac(i, j + 1, k) &
                  &/ (jac(i, j + 1, k) + jac(i, j + 1, k - 1)) + gDEdgeF &
                  &* pDEdgeFGra / dy + gUEdgeF * pUEdgeFGra * 0.25 &
                  &* met23UEdgeF / dz
            else if(k == nz .and. zBoundary == "solid_wall") then
              AFD = - gEdgeF * pEdgeFGra * met23EdgeF / dz + gEdgeD &
                  &* pEdgeDGra * 0.5 * met23EdgeD / dy * jac(i, j + 1, k) &
                  &/ (jac(i, j + 1, k) + jac(i, j + 1, k - 1)) + gDEdgeF &
                  &* pDEdgeFGra / dy
            else
              AFD = - gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz + gEdgeD &
                  &* pEdgeDGra * 0.5 * met23EdgeD / dy * jac(i, j + 1, k) &
                  &/ (jac(i, j + 1, k) + jac(i, j + 1, k - 1)) + gDEdgeF &
                  &* pDEdgeFGra / dy
            end if

            ! ----------------- A(i,j-1,k+1) ------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              ABU = gEdgeB * pEdgeBGra * met23EdgeB / dz - gEdgeU * pEdgeUGra &
                  &* 0.5 * met23EdgeU / dy * jac(i, j - 1, k) / (jac(i, j - 1, &
                  &k) + jac(i, j - 1, k + 1)) - gUEdgeB * pUEdgeBGra / dy
            else if(k == 2 .and. zBoundary == "solid_wall") then
              ABU = gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz - gEdgeU &
                  &* pEdgeUGra * 0.5 * met23EdgeU / dy * jac(i, j - 1, k) &
                  &/ (jac(i, j - 1, k) + jac(i, j - 1, k + 1)) - gUEdgeB &
                  &* pUEdgeBGra / dy - gDEdgeB * pDEdgeBGra * 0.25 &
                  &* met23DEdgeB / dz
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              ABU = gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz - gEdgeU &
                  &* pEdgeUGra * 0.5 * met23EdgeU / dy * jac(i, j - 1, k) &
                  &/ (jac(i, j - 1, k) + jac(i, j - 1, k + 1)) - gUEdgeB &
                  &* pUEdgeBGra * (1.0 / dy - 0.75 * met23UEdgeB / dz)
            else if(k == nz .and. zBoundary == "solid_wall") then
              ABU = 0.0
            else
              ABU = gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz - gEdgeU &
                  &* pEdgeUGra * 0.5 * met23EdgeU / dy * jac(i, j - 1, k) &
                  &/ (jac(i, j - 1, k) + jac(i, j - 1, k + 1)) - gUEdgeB &
                  &* pUEdgeBGra / dy
            end if

            ! ----------------- A(i,j-1,k-1) ------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              ABD = 0.0
            else if(k == 2 .and. zBoundary == "solid_wall") then
              ABD = - gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz - gEdgeD &
                  &* pEdgeDGra * 0.5 * met23EdgeD / dy * jac(i, j - 1, k) &
                  &/ (jac(i, j - 1, k) + jac(i, j - 1, k - 1)) - gDEdgeB &
                  &* pDEdgeBGra * (1.0 / dy + 0.75 * met23DEdgeB / dz)
            else if(k == nz - 1 .and. zBoundary == "solid_wall") then
              ABD = - gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz - gEdgeD &
                  &* pEdgeDGra * 0.5 * met23EdgeD / dy * jac(i, j - 1, k) &
                  &/ (jac(i, j - 1, k) + jac(i, j - 1, k - 1)) - gDEdgeB &
                  &* pDEdgeBGra / dy + gUEdgeB * pUEdgeBGra * 0.25 &
                  &* met23UEdgeB / dz
            else if(k == nz .and. zBoundary == "solid_wall") then
              ABD = - gEdgeB * pEdgeBGra * met23EdgeB / dz - gEdgeD &
                  &* pEdgeDGra * 0.5 * met23EdgeD / dy * jac(i, j - 1, k) &
                  &/ (jac(i, j - 1, k) + jac(i, j - 1, k - 1)) - gDEdgeB &
                  &* pDEdgeBGra / dy
            else
              ABD = - gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz - gEdgeD &
                  &* pEdgeDGra * 0.5 * met23EdgeD / dy * jac(i, j - 1, k) &
                  &/ (jac(i, j - 1, k) + jac(i, j - 1, k - 1)) - gDEdgeB &
                  &* pDEdgeBGra / dy
            end if

            ! ------------------ A(i,j,k+2) -------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              AUU = - gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz - gEdgeL &
                  &* pEdgeLGra * 0.25 * met13EdgeL / dz - gEdgeF * pEdgeFGra &
                  &* 0.25 * met23EdgeF / dz - gEdgeB * pEdgeBGra * 0.25 &
                  &* met23EdgeB / dz + gUEdgeR * pUEdgeRGra * 0.25 &
                  &* met13UEdgeR / dz + gUEdgeL * pUEdgeLGra * 0.25 &
                  &* met13UEdgeL / dz + gUEdgeF * pUEdgeFGra * 0.25 &
                  &* met23UEdgeF / dz + gUEdgeB * pUEdgeBGra * 0.25 &
                  &* met23UEdgeB / dz
            else if((k == nz - 1 .or. k == nz) .and. zBoundary &
                &== "solid_wall") then
              AUU = 0.0
            else
              AUU = gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz + gUEdgeL &
                  &* pUEdgeLGra * 0.25 * met13UEdgeL / dz + gUEdgeF &
                  &* pUEdgeFGra * 0.25 * met23UEdgeF / dz + gUEdgeB &
                  &* pUEdgeBGra * 0.25 * met23UEdgeB / dz
            end if

            ! ------------------ A(i,j,k-2) -------------------!

            if((k == 1 .or. k == 2) .and. zBoundary == "solid_wall") then
              ADD = 0.0
            else if(k == nz .and. zBoundary == "solid_wall") then
              ADD = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz + gEdgeL &
                  &* pEdgeLGra * 0.25 * met13EdgeL / dz + gEdgeF * pEdgeFGra &
                  &* 0.25 * met23EdgeF / dz + gEdgeB * pEdgeBGra * 0.25 &
                  &* met23EdgeB / dz - gDEdgeR * pDEdgeRGra * 0.25 &
                  &* met13DEdgeR / dz - gDEdgeL * pDEdgeLGra * 0.25 &
                  &* met13DEdgeL / dz - gDEdgeF * pDEdgeFGra * 0.25 &
                  &* met23DEdgeF / dz - gDEdgeB * pDEdgeBGra * 0.25 &
                  &* met23DEdgeB / dz
            else
              ADD = - gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz - gDEdgeL &
                  &* pDEdgeLGra * 0.25 * met13DEdgeL / dz - gDEdgeF &
                  &* pDEdgeFGra * 0.25 * met23DEdgeF / dz - gDEdgeB &
                  &* pDEdgeBGra * 0.25 * met23DEdgeB / dz
            end if

            ! ----------------- A(i+1,j,k+2) ------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              ARUU = - gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz + gUEdgeR &
                  &* pUEdgeRGra * 0.25 * met13UEdgeR / dz
            else if((k == nz - 1 .or. k == nz) .and. zBoundary &
                &== "solid_wall") then
              ARUU = 0.0
            else
              ARUU = gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz
            end if

            ! ----------------- A(i+1,j,k-2) ------------------!

            if((k == 1 .or. k == 2) .and. zBoundary == "solid_wall") then
              ARDD = 0.0
            else if(k == nz .and. zBoundary == "solid_wall") then
              ARDD = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz - gDEdgeR &
                  &* pDEdgeRGra * 0.25 * met13DEdgeR / dz
            else
              ARDD = - gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz
            end if

            ! ----------------- A(i-1,j,k+2) ------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              ALUU = - gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz + gUEdgeL &
                  &* pUEdgeLGra * 0.25 * met13UEdgeL / dz
            else if((k == nz - 1 .or. k == nz) .and. zBoundary &
                &== "solid_wall") then
              ALUU = 0.0
            else
              ALUU = gUEdgeL * pUEdgeLGra * 0.25 * met13UEdgeL / dz
            end if

            ! ----------------- A(i-1,j,k-2) ------------------!

            if((k == 1 .or. k == 2) .and. zBoundary == "solid_wall") then
              ALDD = 0.0
            else if(k == nz .and. zBoundary == "solid_wall") then
              ALDD = gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz - gDEdgeL &
                  &* pDEdgeLGra * 0.25 * met13DEdgeL / dz
            else
              ALDD = - gDEdgeL * pDEdgeLGra * 0.25 * met13DEdgeL / dz
            end if

            ! ----------------- A(i,j+1,k+2) ------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              AFUU = - gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz + gUEdgeF &
                  &* pUEdgeFGra * 0.25 * met23UEdgeF / dz
            else if((k == nz - 1 .or. k == nz) .and. zBoundary &
                &== "solid_wall") then
              AFUU = 0.0
            else
              AFUU = gUEdgeF * pUEdgeFGra * 0.25 * met23UEdgeF / dz
            end if

            ! ----------------- A(i,j+1,k-2) ------------------!

            if((k == 1 .or. k == 2) .and. zBoundary == "solid_wall") then
              AFDD = 0.0
            else if(k == nz .and. zBoundary == "solid_wall") then
              AFDD = gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz - gDEdgeF &
                  &* pDEdgeFGra * 0.25 * met23DEdgeF / dz
            else
              AFDD = - gDEdgeF * pDEdgeFGra * 0.25 * met23DEdgeF / dz
            end if

            ! ----------------- A(i,j-1,k+2) ------------------!

            if(k == 1 .and. zBoundary == "solid_wall") then
              ABUU = - gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz + gUEdgeB &
                  &* pUEdgeBGra * 0.25 * met23UEdgeB / dz
            else if((k == nz - 1 .or. k == nz) .and. zBoundary &
                &== "solid_wall") then
              ABUU = 0.0
            else
              ABUU = gUEdgeB * pUEdgeBGra * 0.25 * met23UEdgeB / dz
            end if

            ! ----------------- A(i,j-1,k-2) ------------------!

            if((k == 1 .or. k == 2) .and. zBoundary == "solid_wall") then
              ABDD = 0.0
            else if(k == nz .and. zBoundary == "solid_wall") then
              ABDD = gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz - gDEdgeB &
                  &* pDEdgeBGra * 0.25 * met23DEdgeB / dz
            else
              ABDD = - gDEdgeB * pDEdgeBGra * 0.25 * met23DEdgeB / dz
            end if

            ! Scale the tensor elements.
            AC = AC / (fcscal ** 2.0)
            AR = AR / fcscal / fcscal_r
            AL = AL / fcscal / fcscal_l
            AF = AF / fcscal / fcscal_f
            AB = AB / fcscal / fcscal_b
            AU = AU / fcscal / fcscal_u
            AD = AD / fcscal / fcscal_d
            ARU = ARU / fcscal / fcscal_ru
            ARD = ARD / fcscal / fcscal_rd
            ALU = ALU / fcscal / fcscal_lu
            ALD = ALD / fcscal / fcscal_ld
            AFU = AFU / fcscal / fcscal_fu
            AFD = AFD / fcscal / fcscal_fd
            ABU = ABU / fcscal / fcscal_bu
            ABD = ABD / fcscal / fcscal_bd
            AUU = AUU / fcscal / fcscal_uu
            ADD = ADD / fcscal / fcscal_dd
            ARUU = ARUU / fcscal / fcscal_ruu
            ARDD = ARDD / fcscal / fcscal_rdd
            ALUU = ALUU / fcscal / fcscal_luu
            ALDD = ALDD / fcscal / fcscal_ldd
            AFUU = AFUU / fcscal / fcscal_fuu
            AFDD = AFDD / fcscal / fcscal_fdd
            ABUU = ABUU / fcscal / fcscal_buu
            ABDD = ABDD / fcscal / fcscal_bdd

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
      stop 'ERROR: wrong opt'
    end if

    return

  end subroutine val_PsIn

end module poisson_module
