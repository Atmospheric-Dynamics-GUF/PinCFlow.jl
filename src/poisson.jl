global zBoundary = "solid_wall" # TODO - Fix this

# VERY HACKY THING - TO BE FIXED BY DEVELOPING MY OWN Array type
# JR: Why not use arr[]?
using OffsetArrays
using LinearAlgebra
@inline (arr::Base.Array)(indices...) = arr[indices...]
@inline (arr::OffsetArrays.OffsetArray)(indices...) = arr[indices...]


"""
Linear operator for the Poisson equation.
"""
struct PoissonOperator
    # TODO generic over field type
    op_fields::NamedTuple
    cache::NamedTuple
end

function PoissonOperator(pars::Parameters)
    nx, ny, nz = pars.domain.sizex, pars.domain.sizey, pars.domain.sizez
    # TODO: good name for these, clean up
    el_names = (:ac_b, :acv_b, :ach_b, :al_b, :ar_b, :ab_b, :af_b,
        :ad_b, :au_b, :aru_b, :ard_b, :alu_b, :ald_b, :afu_b,
        :afd_b, :abu_b, :abd_b, :auu_b, :add_b, :aruu_b, :ardd_b,
        :aluu_b, :aldd_b, :afuu_b, :afdd_b, :abuu_b, :abdd_b,
    )
    # elements = NamedTuple{el_names}(Array{Float64,3}(undef, nx, ny, nz) for _ in el_names)
    elements = NamedTuple{el_names}(zeros(nx, ny, nz) for _ in el_names)
    cache_vars = (:p, :r0, :rOld, :r, :s, :b, :t, :v, :rhs_bigc, :sol, :v_pc, :matVec, :s_pc, :q_pc)
    # cache = NamedTuple{cache_vars}(Array{Float64,3}(undef, nx, ny, nz) for _ in cache_vars)
    cache = NamedTuple{cache_vars}(zeros(nx, ny, nz) for _ in cache_vars)
    r_vm = zeros(nx, ny)

    s_aux_field_lin_opr_ = zeros(nx + 2, ny + 2, nz + 2)
    s_aux_field_lin_opr = OffsetArray(s_aux_field_lin_opr_, 0:(nx+1), 0:(ny+1),
        0:(nz+1))
    cache = merge(cache, [:p_pc => copy(r_vm), :r_vm => r_vm, :s_aux_field_lin_opr => s_aux_field_lin_opr])
    return PoissonOperator(elements, cache)
end

function Corrector(model, dt, errFlagBicg, nIter, opt, facray, facprs)
    @trixi_timeit timer() "Corrector" begin
    #! format: noindent
    (; variables, grid) = model
    (; nx, ny, nz) = model.domain
    jac = grid.jac
    flux = model.fluxes
    var = model.variables
    # -------------------------------------------------
    #              Correct uStar, bStar, and p
    # -------------------------------------------------

    # in/out variables
    # type(var_type), intent(inout) :: var
    # type(flux_type), intent(in) :: flux

    # real, intent(in) :: dt, facray, facprs
    # logical, intent(out) :: errFlagBicg
    # integer, intent(out) :: nIter

    # facray multiplies the Rayleigh-damping terms so that they are only
    # handled in the implicit time stepping (sponge and immersed boundary)

    # facprs multiplies the time step so that the routine can be used
    # properly also in the implicit mode (where in sub-step 5 of the
    # semi-implicit scheme the pressure correction is over a full
    # time step, instead of half a time step)

    # opt = expl =>
    # pressure solver for explicit problem and corresponding correction
    # of the winds
    # opt = impl =>
    # pressure solver for implicit problem and corresponding correction
    # of the winds and density fluctuations
    # character(len = *), intent(in)::opt

    # local variables
    # real, dimension[1:nx, 1:ny, 1:nz]::rhs # RHS

    # calculate RHS
    @show model.operator.cache.rhs_bigc[1]
    calc_RHS(model.operator.cache.rhs_bigc, model, dt)
    @show model.operator.cache.rhs_bigc[1]

    # @assert false rhs

    # calculate dp
    poissonSolver(model.operator.cache.rhs_bigc, model, dt, errFlagBicg, nIter, opt, facray, facprs)

    if (errFlagBicg)
        return
    end

    # set horizontal and vertical BC for dp
    pressureBoundaryCondition(model)

    # correct p, rhopStar, and uStar with dp
    correctorStep(model, dt, opt, facray, facprs)
    end # timer
end

function preCond(sIn, sOut, opt, model)
    (; nx, ny, nz) = model.domain
    (; dx, dy) = model.grid
    (; au_b, ac_b, ad_b) = model.operator.op_fields
    (; s_pc, q_pc, p_pc) = model.operator.cache

    # --------------------------------------
    #   preconditioner for BiCGStab
    #   solves vertical problem exploiting its tri-diagonal character
    #   (Isaacson & Keller 1966, see also Durran's book appendix A.2)
    # --------------------------------------

    # in/out variables
    # real, dimension[1:nx, 1:ny, 1:nz], intent(out) :: sOut
    # real, dimension[1:nx, 1:ny, 1:nz], intent(in) :: sIn

    # opt = expl =>
    # pressure solver for explicit problem and corresponding correction
    # of the winds
    # opt = impl =>
    # pressure solver for implicit problem and corresponding correction
    # of the winds and density fluctuations
    # character(len = *), intent(in) :: opt

    # local field
    # real, dimension[1:nx, 1:ny, 1:nz] :: s_pc, q_pc
    # real, dimension(1:nx, 1:ny) :: p_pc

    # local variables
    # integer :: k
    # integer :: i, j
    # integer :: niter

    # real :: deta

    # pseudo timestep

    dtau = 4.0 # NAMELIST PARAMETER (TODO - PLEASE)
    deta = dtau / (2.0 * (1.0 / dx^2 + 1.0 / dy^2))

    # work with auxiliary field s_pc

    s_pc .= 0.0

    maxIterADI = 2 # NAMELIST PARAMETER (TODO - PLEASE)
    for niter in 1:maxIterADI
        if (niter == 0)
            s_pc .= sIn
        else
            # Treat all diagonal elements implicitly.
            linOpr(s_pc, q_pc, opt, "hnd", model)

            @. s_pc = s_pc + deta * (q_pc - sIn)
        end

        # upward sweep

        for j in 1:ny
            for i in 1:nx
                au_b[i, j, nz] = 0.0
            end
        end

        for j in 1:ny
            for i in 1:nx
                if (niter == 0)
                    q_pc[i, j, 1] = -au_b[i, j, 1] / ac_b[i, j, 1]
                    s_pc[i, j, 1] = s_pc[i, j, 1] / ac_b[i, j, 1]
                else
                    # Treat all diagonal elements implicity.
                    q_pc[i, j, 1] = deta * au_b[i, j, 1] / (1.0 - deta * ac_b[i, j, 1])
                    s_pc[i, j, 1] = s_pc[i, j, 1] / (1.0 - deta * ac_b[i, j, 1])
                end
            end
        end

        for k in 2:nz
            for j in 1:ny
                for i in 1:nx
                    if (niter == 0)
                        p_pc[i, j] = 1.0 /
                                     (ac_b[i, j, k] + ad_b[i, j, k] * q_pc(i, j, k - 1))

                        q_pc[i, j, k] = -au_b[i, j, k] * p_pc[i, j]

                        s_pc[i, j, k] = (s_pc[i, j, k] - ad_b[i, j, k] * s_pc(i, j, k - 1)) *
                                        p_pc[i, j]
                    else
                        # Treat all diagonal elements implicitly.
                        p_pc[i, j] = 1.0 / (1.0 - deta * ac_b[i, j, k] -
                                            deta * ad_b[i, j, k] * q_pc(i, j, k - 1))

                        q_pc[i, j, k] = deta * au_b[i, j, k] * p_pc[i, j]

                        s_pc[i, j, k] = (s_pc[i, j, k] +
                                         deta * ad_b[i, j, k] * s_pc(i, j, k - 1)) *
                                        p_pc[i, j]
                    end
                end
            end
        end

        # backward pass

        for k in (nz-1):-1:1
            for j in 1:ny
                for i in 1:nx
                    s_pc[i, j, k] = s_pc[i, j, k] + q_pc[i, j, k] * s_pc(i, j, k + 1)
                end
            end
        end
    end

    # final result

    sOut .= s_pc

    return
end

function linOpr(sIn, Ls, opt, hortot, model)
    (; nx, ny, nz) = model.domain

    (; ac_b, acv_b, ach_b, al_b, ar_b, ab_b, af_b, ad_b, au_b, aru_b, ard_b, alu_b) = model.operator.op_fields
    (; ald_b, afu_b, afd_b, abu_b, abd_b, auu_b, add_b, aruu_b, ardd_b, aluu_b) = model.operator.op_fields
    (; aldd_b, afuu_b, afdd_b, abuu_b, abdd_b) = model.operator.op_fields


    # --------------------------------------
    #   Linear Operator in Poisson problem
    #   Functions as A*x
    # --------------------------------------

    # in/out variables
    # real, dimension[1:nx, 1:ny, 1:nz], intent(out) :: Ls
    # real, dimension[1:nx, 1:ny, 1:nz], intent(in) :: sIn

    # opt = expl =>
    # pressure solver for explicit problem and corresponding correction
    # of the winds
    # opt = impl =>
    # pressure solver for implicit problem and corresponding correction
    # of the winds and density fluctuations
    # character(len = *), intent(in) :: opt

    # hortot = tot =>
    # linear operator for total problem
    # hortot = hor =>
    # linear operator for horizontal problem
    # hortot = hnd =>
    # linear operator for horizontal problem without diagonal term
    # character(len = *), intent(in) :: hortot

    # local field (extended by ghost cells)
    # real, dimension(0:nx + 1, 0:ny + 1, 0:nz + 1) :: s

    # auxiliary fields for "dp"
    # real, dimension(0:ny + 1, 0:nz + 1) :: xSliceLeft_send, xSliceRight_send
    # real, dimension(0:ny + 1, 0:nz + 1) :: xSliceLeft_recv, xSliceRight_recv

    # real, dimension(0:nx + 1, 0:nz + 1) :: ySliceBack_send, ySliceForw_send
    # real, dimension(0:nx + 1, 0:nz + 1) :: ySliceBack_recv, ySliceForw_recv

    # local variables
    # integer :: i, j, k
    # real :: AL, AR, AB, AF, AD, AU, AC, ACH, ACV
    # real :: sL, sR, sB, sF, sD, sU, sC

    # real :: ARU, ARD, ALU, ALD, AFU, AFD, ABU, ABD
    # real :: AUU, ADD, ARUU, ARDD, ALUU, ALDD, AFUU, AFDD, ABUU, ABDD
    # real :: sRU, sRD, sLU, sLD, sFU, sFD, sBU, sBD
    # real :: sUU, sDD, sRUU, sRDD, sLUU, sLDD, sFUU, sFDD, sBUU, sBDD

    # MPI variables
    # integer :: dest, source, tag
    # integer :: sendcount, recvcount

    # work with auxiliary field s
    # s[1:nx, 1:ny, 1:nz] = sIn
    # sIn = zeros(nx, ny, nz)
    s = model.operator.cache.s_aux_field_lin_opr
    s[1:nx, 1:ny, 1:nz] .= sIn

    # Find neighbour procs
    # if (idim > 1) call mpi_cart_shift(comm, 0, 1, left, right, ierror)
    # if (jdim > 1) call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

    if model.parameters.model.model == "pseudo_incompressible"

        #----------------------------
        #   set Halo cells: xSlice
        #----------------------------

        #   if (idim > 1)
        #     # slice size
        #     sendcount = (ny + 2) * (nz + 2)
        #     recvcount = sendcount

        #     # read slice into contiguous array
        #     xSliceLeft_send[:, :] = s(1, :, :)
        #     xSliceRight_send[:, :] = s(nx, :, :)

        #     # left -> right
        #     source = left
        #     dest = right
        #     tag = 100

        #     call mpi_sendrecv(xSliceRight_send(0, 0), sendcount,mpi_double_precision, dest, tag, xSliceLeft_recv(0, 0),recvcount, mpi_double_precision, source, mpi_any_tag, comm,sts_left, ierror)

        #     # right -> left
        #     source = right
        #     dest = left
        #     tag = 100

        #     call mpi_sendrecv(xSliceLeft_send(0, 0), sendcount,mpi_double_precision, dest, tag, xSliceRight_recv(0, 0),recvcount, mpi_double_precision, source, mpi_any_tag, comm,sts_right, ierror)

        #     # right halos
        #     s[nx + 1, :, :] = xSliceRight_recv[:, :]

        #     # left halos
        #     s[0, :, :] = xSliceLeft_recv[:, :]
        #   else
        @. s[0, :, :] = s[nx, :, :]
        @. s[nx+1, :, :] = s[1, :, :]
        #   end

        #------------------------------
        #   set Halo cells: ySlice
        #------------------------------

        #   if (jdim > 1)
        #     # slice size
        #     sendcount = (nx + 2) * (nz + 2)
        #     recvcount = sendcount

        #     # read slice into contiguous array
        #     ySliceBack_send[:, :] = s(:, 1, :)
        #     ySliceForw_send[:, :] = s(:, ny, :)

        #     # back -> forw
        #     source = back
        #     dest = forw
        #     tag = 100

        #     call mpi_sendrecv(ySliceForw_send(0, 0), sendcount,mpi_double_precision, dest, tag, ySliceBack_recv(0, 0),recvcount, mpi_double_precision, source, mpi_any_tag, comm,sts_back, ierror)

        #     # forw -> back
        #     source = forw
        #     dest = back
        #     tag = 100

        #     call mpi_sendrecv(ySliceBack_send(0, 0), sendcount,mpi_double_precision, dest, tag, ySliceForw_recv(0, 0),recvcount, mpi_double_precision, source, mpi_any_tag, comm,sts_right, ierror)

        #     # forward halos
        #     s[:, ny + 1, :] = ySliceForw_recv[:, :]

        #     # backward halos
        #     s[:, 0, :] = ySliceBack_recv[:, :]
        #   else
        @. s[:, 0, :] = s[:, ny, :]
        @. s[:, ny+1, :] = s[:, 1, :]
        #   end

        #---------------------------------
        #         Loop over field
        #---------------------------------

        for k in 1:nz
            for j in 1:ny
                for i in 1:nx

                    # ------------------ A(i+1,j,k) ------------------

                    AR = ar_b[i, j, k]
                    sR = s[i+1, j, k]

                    # ------------------- A(i-1,j,k) --------------------

                    AL = al_b[i, j, k]
                    sL = s(i - 1, j, k)

                    # -------------------- A(i,j+1,k) ----------------------

                    AF = af_b[i, j, k]
                    sF = s(i, j + 1, k)

                    # --------------------- A(i,j-1,k) -----------------------

                    AB = ab_b[i, j, k]
                    sB = s(i, j - 1, k)

                    # --------------------- A(i,j,k+1) ------------------------

                    if (k < nz)
                        AU = au_b[i, j, k]
                        sU = s(i, j, k + 1)
                    else # k = nz -> upwad boundary (solid wall)
                        # A(i,j,nz+1) = 0
                        AU = 0.0
                        sU = 0.0
                    end

                    # --------------------- A(i,j,k-1) ------------------------

                    if (k > 1)
                        AD = ad_b[i, j, k]
                        sD = s(i, j, k - 1)
                    else # k = 1 -> downward boundary (solid wall)
                        # A(i,j,0) = 0
                        AD = 0.0
                        sD = 0.0
                    end

                    # -------------------- A(i,j,k) --------------------------

                    ACH = ach_b[i, j, k]
                    ACV = acv_b[i, j, k]

                    AC = ac_b[i, j, k]
                    sC = s[i, j, k]

                    # -------------------- apply Operator ---------------------

                    if hortot == "tot"
                        Ls[i, j, k] = AL * sL +
                                      AR * sR +
                                      AF * sF +
                                      AB * sB +
                                      AU * sU +
                                      AD * sD +
                                      AC * sC
                    elseif hortot == "hor"
                        Ls[i, j, k] = AL * sL + AR * sR + AF * sF + AB * sB + ACH * sC
                    elseif hortot == "hnd"
                        Ls[i, j, k] = AL * sL + AR * sR + AF * sF + AB * sB
                    else
                        @assert false "wrong hortot in linOpr"
                    end

                    # ----------------- A(i+1,j,k+1) -----------------

                    if (k < nz)
                        ARU = aru_b[i, j, k]
                        sRU = s(i + 1, j, k + 1)
                    else
                        ARU = 0.0
                        sRU = 0.0
                    end

                    # ----------------- A(i+1,j,k-1) -----------------

                    if (k > 1)
                        ARD = ard_b[i, j, k]
                        sRD = s(i + 1, j, k - 1)
                    else
                        ARD = 0.0
                        sRD = 0.0
                    end

                    # ----------------- A(i-1,j,k+1) -----------------

                    if (k < nz)
                        ALU = alu_b[i, j, k]
                        sLU = s(i - 1, j, k + 1)
                    else
                        ALU = 0.0
                        sLU = 0.0
                    end

                    # ----------------- A(i-1,j,k-1) -----------------

                    if (k > 1)
                        ALD = ald_b[i, j, k]
                        sLD = s(i - 1, j, k - 1)
                    else
                        ALD = 0.0
                        sLD = 0.0
                    end

                    # ----------------- A(i,j+1,k+1) -----------------

                    if (k < nz)
                        AFU = afu_b[i, j, k]
                        sFU = s(i, j + 1, k + 1)
                    else
                        AFU = 0.0
                        sFU = 0.0
                    end

                    # ----------------- A(i,j+1,k-1) -----------------

                    if (k > 1)
                        AFD = afd_b[i, j, k]
                        sFD = s(i, j + 1, k - 1)
                    else
                        AFD = 0.0
                        sFD = 0.0
                    end

                    # ----------------- A(i,j-1,k+1) -----------------

                    if (k < nz)
                        ABU = abu_b[i, j, k]
                        sBU = s(i, j - 1, k + 1)
                    else
                        ABU = 0.0
                        sBU = 0.0
                    end

                    # ----------------- A(i,j-1,k-1) -----------------

                    if (k > 1)
                        ABD = abd_b[i, j, k]
                        sBD = s(i, j - 1, k - 1)
                    else
                        ABD = 0.0
                        sBD = 0.0
                    end

                    # ------------------ A(i,j,k+2) -----------------

                    if (k < nz - 1)
                        AUU = auu_b[i, j, k]
                        sUU = s(i, j, k + 2)
                    else
                        AUU = 0.0
                        sUU = 0.0
                    end

                    # ------------------ A(i,j,k-2) -----------------

                    if (k > 2)
                        ADD = add_b[i, j, k]
                        sDD = s(i, j, k - 2)
                    else
                        ADD = 0.0
                        sDD = 0.0
                    end

                    # ----------------- A(i+1,j,k+2) -----------------

                    if (k < nz - 1)
                        ARUU = aruu_b[i, j, k]
                        sRUU = s(i + 1, j, k + 2)
                    else
                        ARUU = 0.0
                        sRUU = 0.0
                    end

                    # ----------------- A(i+1,j,k-2) -----------------

                    if (k > 2)
                        ARDD = ardd_b[i, j, k]
                        sRDD = s(i + 1, j, k - 2)
                    else
                        ARDD = 0.0
                        sRDD = 0.0
                    end

                    # ----------------- A(i-1,j,k+2) -----------------

                    if (k < nz - 1)
                        ALUU = aluu_b[i, j, k]
                        sLUU = s(i - 1, j, k + 2)
                    else
                        ALUU = 0.0
                        sLUU = 0.0
                    end

                    # ----------------- A(i-1,j,k-2) -----------------

                    if (k > 2)
                        ALDD = aldd_b[i, j, k]
                        sLDD = s(i - 1, j, k - 2)
                    else
                        ALDD = 0.0
                        sLDD = 0.0
                    end

                    # ----------------- A(i,j+1,k+2) -----------------

                    if (k < nz - 1)
                        AFUU = afuu_b[i, j, k]
                        sFUU = s(i, j + 1, k + 2)
                    else
                        AFUU = 0.0
                        sFUU = 0.0
                    end

                    # ----------------- A(i,j+1,k-2) -----------------

                    if (k > 2)
                        AFDD = afdd_b[i, j, k]
                        sFDD = s(i, j + 1, k - 2)
                    else
                        AFDD = 0.0
                        sFDD = 0.0
                    end

                    # ----------------- A(i,j-1,k+2) -----------------

                    if (k < nz - 1)
                        ABUU = abuu_b[i, j, k]
                        sBUU = s(i, j - 1, k + 2)
                    else
                        ABUU = 0.0
                        sBUU = 0.0
                    end

                    # ----------------- A(i,j-1,k-2) -----------------

                    if (k > 2)
                        ABDD = abdd_b[i, j, k]
                        sBDD = s(i, j - 1, k - 2)
                    else
                        ABDD = 0.0
                        sBDD = 0.0
                    end

                    # Update operator.
                    Ls[i, j, k] = Ls[i, j, k] +
                                  ARU * sRU +
                                  ARD * sRD +
                                  ALU * sLU +
                                  ALD * sLD +
                                  AFU * sFU +
                                  AFD * sFD +
                                  ABU * sBU +
                                  ABD * sBD +
                                  AUU * sUU +
                                  ADD * sDD +
                                  ARUU * sRUU +
                                  ARDD * sRDD +
                                  ALUU * sLUU +
                                  ALDD * sLDD +
                                  AFUU * sFUU +
                                  AFDD * sFDD +
                                  ABUU * sBUU +
                                  ABDD * sBDD
                end
            end
        end
    else
        @assert false "linOpr: unknown case model"
    end
end

function calc_RHS(b, model, dt)

    # TODO - These are not used
    # sum_local, sum_global = cache.sum_local_bicg, cache.sum_global_bicg
    #----------------------------------------
    #   calculates the RHS of the
    #   Poisson problem
    #----------------------------------------

    # in/out variables
    # type(var_type), intent(in) :: var
    # type(flux_type), intent(in) :: flux
    # real, intent(in) :: d

    # local vars
    # real :: uR, uL, vF, vB, wU, wD
    # real, dimension(1:nz) :: sum_local, sum_global

    # real :: pEdgeR, pEdgeL, pEdgeF, pEdgeB, pEdgeU, pEdgeD

    # integer :: i, j, k
    # real :: divSum

    # real :: bu, bv, bw, bl2loc, divL2_norm, divL2_norm_local

    # check L2-norm of divergence
    # real :: divL2, divMax

    # MPI stuff
    # real :: divL2_local, divSum_local
    # integer :: root

    # real :: fcscal

    #--------------------------------------------------
    #    Setup b = Ma^2 * P * u^*  (right hand side)
    #--------------------------------------------------

    divSum = 0.0
    divL2 = 0.0
    divMax = 0.0

    divSum_local = 0.0
    divL2_local = 0.0

    divL2_norm = 0.0
    divL2_norm_local = 0.0

    var = model.variables.prognostic_fields
    (; grid) = model
    (; nx, ny, nz) = model.domain
    (; jac, dx, dy, dz) = model.grid
    (; pstrattfc, rhostrattfc) = model.atmosphere

    c = model.constants
    if model.parameters.model.model == "pseudo_incompressible"
        # Calculate RHS for TFC.
        for k in 1:nz
            for j in 1:ny
                for i in 1:nx
                    # Calculate scaling factor.
                    fcscal = sqrt(pstrattfc[i, j, k]^2.0 / rhostrattfc[i, j, k])
                    # Store velocities at cell edges.
                    uR = var.u[i, j, k]
                    uL = var.u(i - 1, j, k)
                    vF = var.v[i, j, k]
                    vB = var.v(i, j - 1, k)
                    wU = var.w[i, j, k]
                    wD = var.w(i, j, k - 1)
                    # Calculate P at cell edges.
                    pEdgeR = 0.5 * (jac[i, j, k] * pstrattfc[i, j, k] +
                                    jac[i+1, j, k] * pstrattfc[i+1, j, k])
                    pEdgeL = 0.5 * (jac[i, j, k] * pstrattfc[i, j, k] +
                                    jac(i - 1, j, k) * pstrattfc(i - 1, j, k))
                    pEdgeF = 0.5 * (jac[i, j, k] * pstrattfc[i, j, k] +
                                    jac(i, j + 1, k) * pstrattfc(i, j + 1, k))
                    pEdgeB = 0.5 * (jac[i, j, k] * pstrattfc[i, j, k] +
                                    jac(i, j - 1, k) * pstrattfc(i, j - 1, k))
                    pEdgeU = jac[i, j, k] *
                             jac(i, j, k + 1) *
                             (pstrattfc[i, j, k] + pstrattfc(i, j, k + 1)) /
                             (jac[i, j, k] + jac(i, j, k + 1))
                    pEdgeD = jac[i, j, k] *
                             jac(i, j, k - 1) *
                             (pstrattfc[i, j, k] + pstrattfc(i, j, k - 1)) /
                             (jac[i, j, k] + jac(i, j, k - 1))
                    # Compute RHS.
                    bu = (pEdgeR * uR - pEdgeL * uL) / dx / jac[i, j, k] * c.ma^2.0 * c.kappa
                    bv = (pEdgeF * vF - pEdgeB * vB) / dy / jac[i, j, k] * c.ma^2.0 * c.kappa
                    bw = (pEdgeU * wU - pEdgeD * wD) / dz / jac[i, j, k] * c.ma^2.0 * c.kappa
                    # @show uR, vF, wU
                    # @show pEdgeR, pEdgeL, pEdgeF, dx, dy, dz, jac[i, j, k]
                    # @show bu, bv, bw
                    divSum_local = divSum_local + bu + bv + bw
                    bu = bu / fcscal
                    bv = bv / fcscal
                    bw = bw / fcscal
                    b[i, j, k] = bu + bv + bw
                    # Compute check sum for solvability criterion.
                    divL2_local = divL2_local + b[i, j, k]^2.0
                    bl2loc = bu^2.0 + bv^2.0 + bw^2.0
                    divL2_norm_local = divL2_norm_local + bl2loc
                    if (abs(b[i, j, k]) > divMax)
                        divMax = abs(b[i, j, k])
                    end
                end
            end
        end

        #   !MPI: sum divSum_local over all procs
        #   root = 0
        #   call mpi_reduce(divSum_local, divSum, 1, mpi_double_precision, mpi_sum, root, comm, ierror)

        #   call mpi_bcast(divSum, 1, mpi_double_precision, root, comm, ierror)

        #   !MPI: sum divL2_local over all procs
        #   root = 0
        #   call mpi_reduce(divL2_local, divL2, 1, mpi_double_precision, mpi_sum, root, comm, ierror)

        #   call mpi_bcast(divL2, 1, mpi_double_precision, root, comm, ierror)

        #   !MPI: sum divL2_norm_local over all procs
        #   root = 0
        #   call mpi_reduce(divL2_norm_local, divL2_norm, 1, mpi_double_precision, mpi_sum, root, comm, ierror)

        #   call mpi_bcast(divL2_norm, 1, mpi_double_precision, root, comm, ierror)

        # scale div
        divL2_local = sqrt(divL2_local / nx / ny / nz)
        sizeX = nx
        sizeY = ny
        sizeZ = nz
        divL2 = sqrt(divL2 / sizeX / sizeY / sizeZ)

        divL2_norm_local = sqrt(divL2_norm_local / nx / ny / nz)
        divL2_norm = sqrt(divL2_norm / sizeX / sizeY / sizeZ)

        b_norm = divL2

        if (divL2_norm != 0.0)
            tolref = divL2 / divL2_norm
        else
            if (divL2 == 0.0)
                tolref = 1.0
            else
                @assert false "ERROR: divL2_norm = 0 while divL2 != 0"
            end
        end

    else
        @assert false "poissonSolver: unknown case model."
    end
end

function poissonSolver(b, model, dt, errFlagBicg, nIter, opt, facray, facprs)
    (; nx, ny, nz) = model.domain
    # TODO: save on PoissonSolver
    sol = zeros(nx, ny, nz)
    # -------------------------------------------------
    # solves the Poisson problem with
    # application of linear operator L
    # -------------------------------------------------

    # in/out variables
    # type(var_type), intent(in) :: var
    # real, intent(in) :: dt, facray, facprs

    # logical, intent(out) :: errFlagBicg
    # integer, intent(out) :: nIter

    # real, dimension[1:nx, 1:ny, 1:nz], intent(in) :: b

    # facray multiplies the Rayleigh-damping terms so that they are only
    # handled in the implicit time stepping (sponge and immersed boundary)

    # facprs multiplies the time step so that the routine can be used
    # properly also in the implicit mode (where in sub-step 5 of the
    # semi-implicit scheme the pressure correction is over a full
    # time step, instead of half a time step)

    # opt = expl =>
    # pressure solver for explicit problem and corresponding correction
    # of the winds
    # opt = impl =>
    # pressure solver for implicit problem and corresponding correction
    # of the winds and density fluctuations
    # character(len = *), intent(in) :: opt

    # local vars
    # real, dimension[1:nx, 1:ny, 1:nz] :: sol # solution of Poisson problem

    # real :: res

    # real :: dtInv

    # integer :: i, j, k
    # real :: fcscal

    # Init
    if (dt == 0.0)
        @assert false "poissonSolver: dt = 0.0. Stopping."
    end
    dtInv = 1.0 / dt

    #--------------------------------
    #     Linear equation solver
    #     solve for dt * dp ...
    #--------------------------------

    if model.parameters.model.model == "pseudo_incompressible"
        val_PsIn(model, dt, opt, facray)
    else
        @assert false "linOpr: unknown case model"
    end

    bicgstab(b, dt, model, sol, nIter, errFlagBicg, opt)

    if (errFlagBicg)
        return
    end

    if model.parameters.model.model == "pseudo_incompressible"
        (; pstrattfc, rhostrattfc) = model.atmosphere
        for k in 1:nz
            for j in 1:ny
                for i in 1:nx
                    fcscal = sqrt(pstrattfc[i, j, k]^2 / rhostrattfc[i, j, k])
                    sol[i, j, k] = sol[i, j, k] / fcscal
                end
            end
        end
    end

    # now get dp from dt * dp ...
    # pass solution to pressure corrector
    @. model.atmosphere.dp[1:nx, 1:ny, 1:nz] .= dtInv / facprs * sol
end

function bicgstab(b_in, dt, model, sol, nIter, errFlag, opt)
    (; grid, parameters) = model
    (; matVec, v_pc, r_vm) = model.operator.cache
    (; tolcrit, tolpoisson, tolref, preconditioner) = parameters.poisson
    (; nx, ny, nz) = model.domain

    # --------------------------------------
    #    BiCGStab using linear operator
    #    preconditioner applied via A M^-1 M x = b
    #---------------------------------------

    # in/out variables
    # b_in = zeros(nx, ny, nz)
    # sol = zeros(nx, ny, nz)

    # FROM FORTRAN
    # real, intent(out) :: res # residual
    # integer, intent(out) :: nIter
    # logical, intent(out) :: errFlag

    # opt = expl =>
    # pressure solver for explicit problem and corresponding correction
    # of the winds
    # opt = impl =>
    # pressure solver for implicit problem and corresponding correction
    # of the winds and density fluctuations

    # FROM FORTRAN
    # character(len = *), intent(in) :: opt

    # Local parameters
    # FROM FORTRAN
    # integer :: maxIt

    # local variables
    # FROM FORTRAN
    # integer :: i, j, k, allocstat
    # integer :: j_b
    # real, dimension(:, :, :), allocatable :: p, r0, rOld, r, s, t, v, matVec, v_pc
    # real :: alpha, beta, omega

    # verbose
    giveInfo = true

    # MPI stuff
    # FROM FORTRAN
    # integer :: root
    # real :: res_local

    # real :: b_vm_norm, res_vm # FROM FORTRAN

    if (giveInfo)
        println("")
        println("(a)")
        println("BICGSTAB: solving linear system... ")
        println("(a)")
    end

    sol .= 0.0 # It was = 0.0 in fortran

    # Set parameters

    # modified convergence criterion so that iterations stop when either
    # (a) tolcrit = abs  =>  |Ax - b| < eps b_*
    #     with b_* a suitable norm deciding whether b (= the divergence
    #     criterion for the winds from the predictor) is small or not
    # (b) tolcrit = rel  =>  |Ax - b| < eps |b|
    # here eps = tolPoisson is the user-set convergence criterion
    # hypre has the criterion |Ax - b| < tol * |b|, hence, with
    # tolref = divL2/divL2_norm = |b|/b_*

    cache = model.operator.cache
    p = cache.p
    r0 = cache.r0
    rOld = cache.rOld
    r = cache.r
    s = cache.s
    b = cache.b
    t = cache.t
    v = cache.v
    matVec = cache.matVec
    v_pc = cache.v_pc
    if (tolcrit == "abs")
        tol = tolpoisson / tolref
    elseif (tolcrit == "rel")
        tol = tolpoisson
    end

    # error flag
    errFlag = false

    b .= b_in

    linOpr(sol, matVec, opt, "tot", model)
    # @assert false matVec[150, 1, 1],sol[150,1,1],b[150,1,1]
    r0 .= b - matVec
    p .= r0
    r .= r0

    res_local = 0.0
    for k in 1:nz
        for j in 1:ny
            for i in 1:nx
                res_local = res_local + r[i, j, k]^2
            end
        end
    end

    # MPI find global residual
    # root = 0
    # call mpi_reduce(res_local, res, 1, mpi_double_precision, mpi_sum, root, comm, ierror)

    # call mpi_bcast(res, 1, mpi_double_precision, root, comm, ierror)

    res = res_local

    sizeX = nx
    sizeY = ny
    sizeZ = nz
    res = sqrt(res / sizeX / sizeY / sizeZ)

    b_norm = res

    r_vm .= 0.0
    for k in 1:nz
        @. r_vm[:, :] = r_vm[:, :] + r[:, :, k]
    end
    r_vm .= r_vm / sizeZ

    res_local = 0.0
    for j in 1:ny
        for i in 1:nx
            res_local = res_local + r_vm[i, j]^2
        end
    end

    # root = 0
    # call mpi_reduce(res_local, res_vm, 1, mpi_double_precision, mpi_sum, root, comm, ierror)

    # call mpi_bcast(res_vm, 1, mpi_double_precision, root, comm, ierror)

    res_vm = res_local
    res_vm = sqrt(res_vm / sizeX / sizeY)

    b_vm_norm = res_vm

    if (res == 0.0 || res / b_norm <= tol)
        if (giveInfo)
            println(" ==> no iteration needed.")
        end
        nIter = 0
        return
    end

    # Loop

    println("starting loop over $(model.parameters.poisson.maxiter)")
    for j_b in 1:model.parameters.poisson.maxiter
        # @assert false p[150, 1, 1]
        # v = A*p
        if (preconditioner == "yes")
            preCond(p, v_pc, opt, model)
        else
            v_pc .= p
        end
        linOpr(v_pc, matVec, opt, "tot", model)

        v .= matVec

        # @assert false r[150, 1, 1], r0[150, 1, 1], v[150, 1, 1], p[150, 1, 1], dot(r,r0), dot(v,r0)
        # @show size(r)
        alpha = dot(r, r0) / dot(v, r0)
        @. s = r - alpha * v

        # t = A*s
        if (preconditioner == "yes")
            preCond(s, v_pc, opt, model)
        else
            v_pc .= s
        end
        linOpr(v_pc, matVec, opt, "tot", model)
        t .= matVec

        omega = dot(t, s) / dot(t, t)
        sol .= sol + alpha * p + omega * s

        rOld .= r
        r .= s - omega * t

        #-----------------------
        #   Abort criterion
        #-----------------------

        res_local = 0.0
        for k in 1:nz
            for j in 1:ny
                for i in 1:nx
                    res_local = res_local + r[i, j, k]^2
                end
            end
        end

        # MPI find global residual
        root = 0
        #   call mpi_reduce(res_local, res, 1, mpi_double_precision, mpi_sum, root, comm, ierror)

        #   call mpi_bcast(res, 1, mpi_double_precision, root, comm, ierror)

        res = res_local
        res = sqrt(res / sizeX / sizeY / sizeZ)

        r_vm .= 0.0
        for k in 1:nz
            r_vm[:, :] .= r_vm[:, :] + r(:, :, k)
        end
        r_vm .= r_vm ./ sizeZ

        res_local = 0.0
        for j in 1:ny
            for i in 1:nx
                res_local = res_local + r_vm[i, j]^2
            end
        end

        root = 0
        #   call mpi_reduce(res_local, res_vm, 1, mpi_double_precision, mpi_sum, root, comm, ierror)

        #   call mpi_bcast(res_vm, 1, mpi_double_precision, root, comm, ierror)

        res_vm = res_local
        res_vm = sqrt(res_vm / sizeX / sizeY)
        if (max(res / b_norm, res_vm / b_vm_norm) <= tol)
            if (giveInfo)
                println(" Nb.of iterations: j = ", j_b)
                println(" Final residual: res = ", res / b_norm)
                println(" Final residual v.m. = ", res_vm / b_vm_norm)
                println(" ")
            end

            nIter = j_b

            if (preconditioner == "yes")
                s .= sol
                preCond(s, sol, opt, semi)
            end

            return
        end

        beta = alpha / omega * dot(r, r0) / dot(rOld, r0)

        # @show p[150, 1, 1], r[150, 1, 1], v[150, 1, 1], beta, alpha, omega
        @. p = r + beta * (p - omega * v)

        # @assert false p[150, 1, 1], r[150, 1, 1], v[150, 1, 1], beta, alpha, omega

    end

    # max iteration

    errFlag = true
    println("BICGSTAB: max iteration reached.")
    # nIter = maxIt

end

function pressureBoundaryCondition(semi)
    (; grid, cache) = semi
    (; nx, ny, nz) = grid
    (; dp) = cache
    #--------------------------------------------------
    # set pressure correction dp in ghost cells for BC
    #--------------------------------------------------

    # auxiliary fields for "dp"
    # real, dimension(0:ny+1, 0:nz+1)::xSliceLeft_send, xSliceRight_send
    # real, dimension(0:ny+1, 0:nz+1)::xSliceLeft_recv, xSliceRight_recv

    # real, dimension(0:nx+1, 0:nz+1)::ySliceBack_send, ySliceForw_send
    # real, dimension(0:nx+1, 0:nz+1)::ySliceBack_recv, ySliceForw_recv

    # MPI STUFF (BIG TODO!)

    # xSliceLeft_send_ = zeros(ny + 2, nz + 2)
    # xSliceLeft_send = OffsetArray(xSliceLeft_send_, 0:ny+1, 0:nz+1)

    # xSliceRight_send,
    # xSliceLeft_recv,
    # xSliceRight_recv,
    # ySliceBack_send,
    # ySliceForw_send,
    # ySliceBack_recv,
    # ySliceForw_recv = (copy(xSliceLeft_send) for _ = 1:7)

    # MPI variables

    # MPI STUFF IGNORED
    # integer :: dest, source, tag
    # integer :: sendcount, recvcount

    # Find neighbour procs

    # MPI STUFF IGNORED
    # if (idim > 1) call mpi_cart_shift(comm, 0, 1, left, right, ierror)
    # if (jdim > 1) call mpi_cart_shift(comm, 1, 1, back, forw, ierror)

    #----------------------------
    #   set Halo cells: xSlice
    #----------------------------

    # if (idim > 1)
    #   # slice size
    #   sendcount = (ny + 2) * (nz + 2)
    #   recvcount = sendcount

    #   # read slice into contiguous array
    #   xSliceLeft_send[:, :] = dp(1, 0:ny + 1, :)
    #   xSliceRight_send[:, :] = dp(nx, 0:ny + 1, :)

    #   # left -> right
    #   source = left
    #   dest = right
    #   tag = 100

    #   call mpi_sendrecv(xSliceRight_send(0, 0), sendcount, mpi_double_precision, dest, tag, xSliceLeft_recv(0, 0), recvcount, mpi_double_precision, source, mpi_any_tag, comm, sts_left, ierror)

    #   # right -> left
    #   source = right
    #   dest = left
    #   tag = 100

    #   call mpi_sendrecv(xSliceLeft_send(0, 0), sendcount, mpi_double_precision, dest, tag, xSliceRight_recv(0, 0), recvcount, mpi_double_precision, source, mpi_any_tag, comm, sts_right, ierror)

    #   # right halos
    #   dp(nx + 1, 0:ny + 1, :) = xSliceRight_recv(0:ny + 1, :)

    #   # left halos
    #   dp(0, 0:ny + 1, :) = xSliceLeft_recv(0:ny + 1, :)

    # else

    @. dp[0, :, :] = dp[nx, :, :]
    @. dp[nx+1, :, :] = dp[1, :, :]

    # end

    #------------------------------
    #   set Halo cells: ySlice
    #------------------------------

    # if (jdim > 1)
    #   # slice size
    #   sendcount = (nx + 2) * (nz + 2)

    #   recvcount = sendcount

    #   # read slice into contiguous array
    #   ySliceBack_send[:, :] = dp(0:nx + 1, 1, :)

    #   ySliceForw_send[:, :] = dp(0:nx + 1, ny, :)

    #   # back -> forw
    #   source = back
    #   dest = forw
    #   tag = 100

    #   call mpi_sendrecv(ySliceForw_send(0, 0), sendcount, mpi_double_precision, dest, tag, ySliceBack_recv(0, 0), recvcount, mpi_double_precision, source, mpi_any_tag, comm, sts_back, ierror)

    #   # forw -> back
    #   source = forw
    #   dest = back
    #   tag = 100

    #   call mpi_sendrecv(ySliceBack_send(0, 0), sendcount, mpi_double_precision, dest, tag, ySliceForw_recv(0, 0), recvcount, mpi_double_precision, source, mpi_any_tag, comm, sts_right, ierror)

    #   # forward halos
    #   dp(0:nx + 1, ny + 1, :) = ySliceForw_recv(0:nx + 1, :)

    #   # backward halos
    #   dp(0:nx + 1, 0, :) = ySliceBack_recv(0:nx + 1, :)

    # else

    @. dp[:, 0, :] = dp[:, ny, :]
    @. dp[:, ny+1, :] = dp[:, 1, :]

    # end

    #----------------
    #   z-Boundary
    #----------------

    if zBoundary == "solid_wall"
        @. dp[:, :, 0] = dp[:, :, 1]
        @. dp[:, :, nz+1] = dp[:, :, nz]
    else
        @assert false "pressureBoundaryCondition: unknown case zBoundary."
    end
end

function correctorStep(semi, dt, opt, facray, facprs)
    (; cache, grid, met, equations, spongeLayer, sponge_uv) = semi
    (; var, dp, rhostrattfc, pstrattfc, jac, bvsStrat, kr_sp_tfc, kr_sp_w_tfc, corX, corY) = cache
    (; nx, ny, nz, dx, dy, dz) = grid
    (; kappaInv, MaInv2, g_ndim) = equations

    #------------------------------------------------
    #         correct pressure & velocity
    #------------------------------------------------

    # in/out variables

    # type(var_type), intent(inout) :: var (FROM FORTRAN)
    # real, intent(in) :: dt, facray, facprs (FROM FORTRAN)

    # facray multiplies the Rayleigh-damping terms so that they are only
    # handled in the implicit time stepping (sponge and immersed boundary)

    # facprs multiplies the time step so that the routine can be used
    # properly also in the implicit mode (where in sub-step 5 of the
    # semi-implicit scheme the pressure correction is over a full
    # time step, instead of half a time step)

    # opt = expl =>
    # pressure solver for explicit problem and corresponding correction
    # of the winds
    # opt = impl =>
    # pressure solver for implicit problem and corresponding correction
    # of the winds and density fluctuations
    # character(len = *), intent(in) :: opt (FROM FORTRAN)

    # local variables

    ## From Fortran
    # integer :: i, j, k
    # real :: rhoEdge, rhou, rhov, rho
    # real :: pGradX, pGradY, pGradZ
    # real :: du, dv, dw, db
    # real :: facu, facv, facw
    # real :: bvsstw

    # real :: rhow0, rhowm

    # real :: rhostrattfcEdgeU
    # real :: pEdgeR, pEdgeF, pEdgeU, pEdgeD
    # real :: dpEdgeR, dpUEdgeR, dpUUEdgeR, dpDEdgeR, dpDDEdgeR, dpEdgeF, dpUEdgeF, dpUUEdgeF, dpDEdgeF, dpDDEdgeF, dpREdgeU, dpLEdgeU, dpFEdgeU, dpBEdgeU, dpREdgeD, dpLEdgeD, dpFEdgeD, dpBEdgeD
    # real :: met13EdgeR, met23EdgeF, met13EdgeU, met23EdgeU, met33EdgeU, met13EdgeD, met23EdgeD, met33EdgeD
    # real :: pGradZEdgeU, pGradZEdgeD
    # real, dimension((- nbx):(nx + nbx), (- nby):(ny + nby), (- nbz):(nz + nbz)) :: corX, corY

    # integer :: k0, k1

    # --------------------------------------
    #             calc p + dp
    # --------------------------------------

    @. var.exner[0:(nx+1), 0:(ny+1), 0:(nz+1)] = var.exner[0:(nx+1), 0:(ny+1),
        0:(nz+1)] +
                                                 dp[0:(nx+1), 0:(ny+1),
        0:(nz+1)]

    if (opt == "impl")
        @. kr_sp_tfc = kr_sp_tfc * facray
        @. kr_sp_w_tfc = kr_sp_w_tfc * facray
    end

    # --------------------------------------
    #           calc du and u + du
    # --------------------------------------

    if (opt == "impl")
        for k in 1:nz
            for j in 1:ny
                for i in 0:nx
                    facu = 1.0

                    if (spongeLayer && sponge_uv)
                        facu = facu +
                               dt * 0.5 * (kr_sp_tfc[i, j, k] + kr_sp_tfc[i+1, j, k])
                    end

                    facv = facu

                    # Compute values at cell edges.
                    rhou = 0.5 * (var.rho[i, j, k] +
                                  var.rho[i+1, j, k] +
                                  rhostrattfc[i, j, k] +
                                  rhostrattfc[i+1, j, k])
                    pEdgeR = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i+1, j, k])
                    met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))
                    # Compute pressure difference gradient component.
                    if (k == 1 && zBoundary == "solid_wall")
                        dpUUEdgeR = 0.5 * (dp(i, j, k + 2) + dp(i + 1, j, k + 2))
                        dpUEdgeR = 0.5 * (dp(i, j, k + 1) + dp(i + 1, j, k + 1))
                        dpEdgeR = 0.5 * (dp[i, j, k] + dp[i+1, j, k])
                        pGradX = kappaInv * MaInv2 / rhou *
                                 pEdgeR *
                                 ((dp[i+1, j, k] - dp[i, j, k]) / dx +
                                  met13EdgeR *
                                  (-dpUUEdgeR + 4.0 * dpUEdgeR - 3.0 * dpEdgeR) *
                                  0.5 / dz)
                    elseif (k == nz && zBoundary == "solid_wall")
                        dpDDEdgeR = 0.5 * (dp(i, j, k - 2) + dp(i + 1, j, k - 2))
                        dpDEdgeR = 0.5 * (dp(i, j, k - 1) + dp(i + 1, j, k - 1))
                        dpEdgeR = 0.5 * (dp[i, j, k] + dp[i+1, j, k])
                        pGradX = kappaInv * MaInv2 / rhou *
                                 pEdgeR *
                                 ((dp[i+1, j, k] - dp[i, j, k]) / dx +
                                  met13EdgeR *
                                  (dpDDEdgeR - 4.0 * dpDEdgeR + 3.0 * dpEdgeR) *
                                  0.5 / dz)
                    else
                        dpUEdgeR = 0.5 * (dp(i, j, k + 1) + dp(i + 1, j, k + 1))
                        dpDEdgeR = 0.5 * (dp(i, j, k - 1) + dp(i + 1, j, k - 1))
                        pGradX = kappaInv * MaInv2 / rhou *
                                 pEdgeR *
                                 ((dp[i+1, j, k] - dp[i, j, k]) / dx +
                                  met13EdgeR * (dpUEdgeR - dpDEdgeR) * 0.5 / dz)
                    end
                    # Compute velocity correction.
                    corX[i, j, k] = facprs * dt / facu * pGradX
                    du = -corX[i, j, k]

                    var.u[i, j, k] = var.u[i, j, k] + du
                end
            end
        end
    elseif (opt == "expl")
        if (facprs != 1.0)
            "ERROR: wrong facprs in explicit sub-step"
        end
        for k in 1:nz
            for j in 1:ny
                for i in 0:nx
                    # Compute values at cell edges.
                    rhou = 0.5 * (var.rho[i, j, k] +
                                  var.rho[i+1, j, k] +
                                  rhostrattfc[i, j, k] +
                                  rhostrattfc[i+1, j, k])
                    pEdgeR = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i+1, j, k])
                    met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))
                    # Compute pressure difference gradient component.
                    # zBoundary = "solid_wall" # TODO - SERIOUS ISSUE. TREAT URGENTLY!
                    if (k == 1 && zBoundary == "solid_wall")
                        dpUUEdgeR = 0.5 * (dp(i, j, k + 2) + dp(i + 1, j, k + 2))
                        dpUEdgeR = 0.5 * (dp(i, j, k + 1) + dp(i + 1, j, k + 1))
                        dpEdgeR = 0.5 * (dp[i, j, k] + dp[i+1, j, k])
                        pGradX = kappaInv * MaInv2 / rhou *
                                 pEdgeR *
                                 ((dp[i+1, j, k] - dp[i, j, k]) / dx +
                                  met13EdgeR *
                                  (-dpUUEdgeR + 4.0 * dpUEdgeR - 3.0 * dpEdgeR) *
                                  0.5 / dz)
                    elseif (k == nz && zBoundary == "solid_wall")
                        dpDDEdgeR = 0.5 * (dp(i, j, k - 2) + dp(i + 1, j, k - 2))
                        dpDEdgeR = 0.5 * (dp(i, j, k - 1) + dp(i + 1, j, k - 1))
                        dpEdgeR = 0.5 * (dp[i, j, k] + dp[i+1, j, k])
                        pGradX = kappaInv * MaInv2 / rhou *
                                 pEdgeR *
                                 ((dp[i+1, j, k] - dp[i, j, k]) / dx +
                                  met13EdgeR *
                                  (dpDDEdgeR - 4.0 * dpDEdgeR + 3.0 * dpEdgeR) *
                                  0.5 / dz)
                    else
                        dpUEdgeR = 0.5 * (dp(i, j, k + 1) + dp(i + 1, j, k + 1))
                        dpDEdgeR = 0.5 * (dp(i, j, k - 1) + dp(i + 1, j, k - 1))
                        pGradX = kappaInv * MaInv2 / rhou *
                                 pEdgeR *
                                 ((dp[i+1, j, k] - dp[i, j, k]) / dx +
                                  met13EdgeR * (dpUEdgeR - dpDEdgeR) * 0.5 / dz)
                    end

                    du = dt * pGradX # TODO - BIGGEST TODO OF THIS CODE

                    var.u[i, j, k] = var.u[i, j, k] + du
                end
            end
        end
    else
        @assert "ERROR: wrong opt in correctorStep"
    end

    #-------------------------------------
    #         calc dv and v + dv
    #-------------------------------------

    if (opt == "impl")
        for k in 1:nz
            for j in 0:ny
                for i in 1:nx
                    facv = 1.0

                    if (spongeLayer && sponge_uv)
                        facv = facv +
                               dt * 0.5 * (kr_sp_tfc[i, j, k] + kr_sp_tfc(i, j + 1, k))
                    end

                    facu = facv

                    # Compute values at cell edges.
                    rhov = 0.5 * (var.rho[i, j, k] +
                                  var.rho(i, j + 1, k) +
                                  rhostrattfc[i, j, k] +
                                  rhostrattfc(i, j + 1, k))
                    pEdgeF = 0.5 * (pstrattfc[i, j, k] + pstrattfc(i, j + 1, k))
                    met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3))
                    # Compute pressure difference gradient component.
                    if (k == 1 && zBoundary == "solid_wall")
                        dpUUEdgeF = 0.5 * (dp(i, j, k + 2) + dp(i, j + 1, k + 2))
                        dpUEdgeF = 0.5 * (dp(i, j, k + 1) + dp(i, j + 1, k + 1))
                        dpEdgeF = 0.5 * (dp[i, j, k] + dp(i, j + 1, k))
                        pGradY = kappaInv * MaInv2 / rhov *
                                 pEdgeF *
                                 ((dp(i, j + 1, k) - dp[i, j, k]) / dy +
                                  met23EdgeF *
                                  (-dpUUEdgeF + 4.0 * dpUEdgeF - 3.0 * dpEdgeF) *
                                  0.5 / dz)
                    elseif (k == nz && zBoundary == "solid_wall")
                        dpDDEdgeF = 0.5 * (dp(i, j, k - 2) + dp(i, j + 1, k - 2))
                        dpDEdgeF = 0.5 * (dp(i, j, k - 1) + dp(i, j + 1, k - 1))
                        dpEdgeF = 0.5 * (dp[i, j, k] + dp(i, j + 1, k))
                        pGradY = kappaInv * MaInv2 / rhov *
                                 pEdgeF *
                                 ((dp(i, j + 1, k) - dp[i, j, k]) / dy +
                                  met23EdgeF *
                                  (dpDDEdgeF - 4.0 * dpDEdgeF + 3.0 * dpEdgeF) *
                                  0.5 / dz)
                    else
                        dpUEdgeF = 0.5 * (dp(i, j, k + 1) + dp(i, j + 1, k + 1))
                        dpDEdgeF = 0.5 * (dp(i, j, k - 1) + dp(i, j + 1, k - 1))
                        pGradY = kappaInv * MaInv2 / rhov *
                                 pEdgeF *
                                 ((dp(i, j + 1, k) - dp[i, j, k]) / dy +
                                  met23EdgeF * (dpUEdgeF - dpDEdgeF) * 0.5 / dz)
                    end
                    # Compute velocity correction.
                    corY[i, j, k] = facprs * dt / facv * pGradY
                    dv = -corY[i, j, k]

                    var.v[i, j, k] = var.v[i, j, k] + dv
                end
            end
        end
    elseif (opt == "expl")
        if (facprs != 1.0)
            @assert false "ERROR: wrong facprs in explicit sub-step"
        end
        for k in 1:nz
            for j in 0:ny
                for i in 1:nx
                    # Compute values at cell edges.
                    rhov = 0.5 * (var.rho[i, j, k] +
                                  var.rho(i, j + 1, k) +
                                  rhostrattfc[i, j, k] +
                                  rhostrattfc(i, j + 1, k))
                    pEdgeF = 0.5 * (pstrattfc[i, j, k] + pstrattfc(i, j + 1, k))
                    met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3))
                    # Compute pressure difference gradient component.
                    # zBoundary = "solid_wall" # TODO - SERIOUS ISSUE. TREAT URGENTLY!
                    if (k == 1 && zBoundary == "solid_wall")
                        dpUUEdgeF = 0.5 * (dp(i, j, k + 2) + dp(i, j + 1, k + 2))
                        dpUEdgeF = 0.5 * (dp(i, j, k + 1) + dp(i, j + 1, k + 1))
                        dpEdgeF = 0.5 * (dp[i, j, k] + dp(i, j + 1, k))
                        pGradY = kappaInv * MaInv2 / rhov *
                                 pEdgeF *
                                 ((dp(i, j + 1, k) - dp[i, j, k]) / dy +
                                  met23EdgeF *
                                  (-dpUUEdgeF + 4.0 * dpUEdgeF - 3.0 * dpEdgeF) *
                                  0.5 / dz)
                    elseif (k == nz && zBoundary == "solid_wall")
                        dpDDEdgeF = 0.5 * (dp(i, j, k - 2) + dp(i, j + 1, k - 2))
                        dpDEdgeF = 0.5 * (dp(i, j, k - 1) + dp(i, j + 1, k - 1))
                        dpEdgeF = 0.5 * (dp[i, j, k] + dp(i, j + 1, k))
                        pGradY = kappaInv * MaInv2 / rhov *
                                 pEdgeF *
                                 ((dp(i, j + 1, k) - dp[i, j, k]) / dy +
                                  met23EdgeF *
                                  (dpDDEdgeF - 4.0 * dpDEdgeF + 3.0 * dpEdgeF) *
                                  0.5 / dz)
                    else
                        dpUEdgeF = 0.5 * (dp(i, j, k + 1) + dp(i, j + 1, k + 1))
                        dpDEdgeF = 0.5 * (dp(i, j, k - 1) + dp(i, j + 1, k - 1))
                        pGradY = kappaInv * MaInv2 / rhov *
                                 pEdgeF *
                                 ((dp(i, j + 1, k) - dp[i, j, k]) / dy +
                                  met23EdgeF * (dpUEdgeF - dpDEdgeF) * 0.5 / dz)
                    end

                    dv = -dt * pGradY

                    var.v[i, j, k] = var.v[i, j, k] + dv
                end
            end
        end
    else
        @assert false "ERROR: wrong opt in correctorStep"
    end

    #-------------------------------------
    #         calc w and  w + dw
    #-------------------------------------

    # zBoundary = "solid_wall" # TODO - SERIOUS ISSUE. TREAT URGENTLY!
    if zBoundary == "solid_wall"
        k0 = 1
        k1 = nz - 1
    else
        @assert false "correctorStep: unknown case zBoundary."
    end

    if (opt == "impl")
        for k in k0:k1
            for j in 1:ny
                for i in 1:nx
                    facw = 1.0

                    if (spongeLayer)
                        facw = facw +
                               dt * (jac(i, j, k + 1) * kr_sp_w_tfc[i, j, k] +
                                     jac[i, j, k] * kr_sp_w_tfc(i, j, k + 1)) /
                               (jac[i, j, k] + jac(i, j, k + 1))
                    end

                    # Compute values at cell edges.
                    rhostrattfcEdgeU = (jac(i, j, k + 1) * rhostrattfc[i, j, k] +
                                        jac[i, j, k] * rhostrattfc(i, j, k + 1)) /
                                       (jac[i, j, k] + jac(i, j, k + 1))
                    rhoEdge = (jac(i, j, k + 1) * var.rho[i, j, k] +
                               jac[i, j, k] * var.rho(i, j, k + 1)) /
                              (jac[i, j, k] + jac(i, j, k + 1)) + rhostrattfcEdgeU
                    pEdgeU = (jac(i, j, k + 1) * pstrattfc[i, j, k] +
                              jac[i, j, k] * pstrattfc(i, j, k + 1)) /
                             (jac[i, j, k] + jac(i, j, k + 1))
                    bvsstw = (jac(i, j, k + 1) * bvsStrat[i, j, k] +
                              jac[i, j, k] * bvsStrat(i, j, k + 1)) /
                             (jac[i, j, k] + jac(i, j, k + 1))
                    met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 1, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 2, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 3, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    dpREdgeU = (jac(i + 1, j, k + 1) * dp[i+1, j, k] +
                                jac[i+1, j, k] * dp(i + 1, j, k + 1)) /
                               (jac[i+1, j, k] + jac(i + 1, j, k + 1))
                    dpLEdgeU = (jac(i - 1, j, k + 1) * dp(i - 1, j, k) +
                                jac(i - 1, j, k) * dp(i - 1, j, k + 1)) /
                               (jac(i - 1, j, k) + jac(i - 1, j, k + 1))
                    dpFEdgeU = (jac(i, j + 1, k + 1) * dp(i, j + 1, k) +
                                jac(i, j + 1, k) * dp(i, j + 1, k + 1)) /
                               (jac(i, j + 1, k) + jac(i, j + 1, k + 1))
                    dpBEdgeU = (jac(i, j - 1, k + 1) * dp(i, j - 1, k) +
                                jac(i, j - 1, k) * dp(i, j - 1, k + 1)) /
                               (jac(i, j - 1, k) + jac(i, j - 1, k + 1))
                    # Compute pressure difference gradient component.
                    pGradZ = kappaInv * MaInv2 / rhoEdge *
                             pEdgeU *
                             (met13EdgeU * (dpREdgeU - dpLEdgeU) * 0.5 / dx +
                              met23EdgeU * (dpFEdgeU - dpBEdgeU) * 0.5 / dy +
                              met33EdgeU * (dp(i, j, k + 1) - dp[i, j, k]) / dz)
                    # Compute velocity correction.
                    dw = -facprs * dt / (facw + rhostrattfcEdgeU / rhoEdge * bvsstw * dt^2.0) *
                         pGradZ -
                         1.0 / (facw + rhostrattfcEdgeU / rhoEdge * bvsstw * dt^2.0) *
                         rhostrattfcEdgeU / rhoEdge *
                         bvsstw *
                         dt^2.0 *
                         0.5 *
                         (jac(i, j, k + 1) *
                          (met(i, j, k, 1, 3) * (corX[i, j, k] + corX(i - 1, j, k)) +
                           met(i, j, k, 2, 3) * (corY[i, j, k] + corY(i, j - 1, k))) +
                          jac[i, j, k] * (met(i, j, k + 1, 1, 3) *
                                          (corX(i, j, k + 1) + corX(i - 1, j, k + 1)) +
                                          met(i, j, k + 1, 2, 3) *
                                          (corY(i, j, k + 1) + corY(i, j - 1, k + 1)))) /
                         (jac[i, j, k] + jac(i, j, k + 1))

                    var.w[i, j, k] = var.w[i, j, k] + dw
                end
            end
        end

    elseif (opt == "expl")
        if (facprs != 1.0)
            @assert false "ERROR: wrong facprs in explicit sub-step"
        end
        for k in k0:k1
            for j in 1:ny
                for i in 1:nx
                    # Compute values at cell edges.
                    rhostrattfcEdgeU = (jac(i, j, k + 1) * rhostrattfc[i, j, k] +
                                        jac[i, j, k] * rhostrattfc(i, j, k + 1)) /
                                       (jac[i, j, k] + jac(i, j, k + 1))
                    rhoEdge = (jac(i, j, k + 1) * var.rho[i, j, k] +
                               jac[i, j, k] * var.rho(i, j, k + 1)) /
                              (jac[i, j, k] + jac(i, j, k + 1)) + rhostrattfcEdgeU
                    pEdgeU = (jac(i, j, k + 1) * pstrattfc[i, j, k] +
                              jac[i, j, k] * pstrattfc(i, j, k + 1)) /
                             (jac[i, j, k] + jac(i, j, k + 1))
                    bvsstw = (jac(i, j, k + 1) * bvsStrat[i, j, k] +
                              jac[i, j, k] * bvsStrat(i, j, k + 1)) /
                             (jac[i, j, k] + jac(i, j, k + 1))
                    met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 1, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 2, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 3, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    dpREdgeU = (jac(i + 1, j, k + 1) * dp[i+1, j, k] +
                                jac[i+1, j, k] * dp(i + 1, j, k + 1)) /
                               (jac[i+1, j, k] + jac(i + 1, j, k + 1))
                    dpLEdgeU = (jac(i - 1, j, k + 1) * dp(i - 1, j, k) +
                                jac(i - 1, j, k) * dp(i - 1, j, k + 1)) /
                               (jac(i - 1, j, k) + jac(i - 1, j, k + 1))
                    dpFEdgeU = (jac(i, j + 1, k + 1) * dp(i, j + 1, k) +
                                jac(i, j + 1, k) * dp(i, j + 1, k + 1)) /
                               (jac(i, j + 1, k) + jac(i, j + 1, k + 1))
                    dpBEdgeU = (jac(i, j - 1, k + 1) * dp(i, j - 1, k) +
                                jac(i, j - 1, k) * dp(i, j - 1, k + 1)) /
                               (jac(i, j - 1, k) + jac(i, j - 1, k + 1))
                    # Compute pressure difference gradient component.
                    pGradZ = kappaInv * MaInv2 / rhoEdge *
                             pEdgeU *
                             (met13EdgeU * (dpREdgeU - dpLEdgeU) * 0.5 / dx +
                              met23EdgeU * (dpFEdgeU - dpBEdgeU) * 0.5 / dy +
                              met33EdgeU * (dp(i, j, k + 1) - dp[i, j, k]) / dz)
                    # Correct vertical velocity.
                    dw = -dt * pGradZ

                    var.w[i, j, k] = var.w[i, j, k] + dw
                end
            end
        end

    else
        @assert false "ERROR: wrong opt in correctorStep"
    end

    #-----------------------------------------------------------------
    #         calc rhop and rhop + drhop (only for implicit time step)
    #-----------------------------------------------------------------

    if (opt == "impl")
        for k in 1:nz
            for j in 1:ny
                for i in 1:nx
                    facw = 1.0

                    if (spongeLayer)
                        facw = facw + dt * kr_sp_w_tfc[i, j, k]
                    end

                    # Compute P coefficients.
                    pEdgeU = (jac(i, j, k + 1) * pstrattfc[i, j, k] +
                              jac[i, j, k] * pstrattfc(i, j, k + 1)) /
                             (jac[i, j, k] + jac(i, j, k + 1))
                    pEdgeD = (jac(i, j, k - 1) * pstrattfc[i, j, k] +
                              jac[i, j, k] * pstrattfc(i, j, k - 1)) /
                             (jac[i, j, k] + jac(i, j, k - 1))
                    # Compute density coefficients.
                    rhow0 = (jac(i, j, k + 1) * (var.rho[i, j, k] + rhostrattfc[i, j, k]) +
                             jac[i, j, k] * (var.rho(i, j, k + 1) + rhostrattfc(i, j, k + 1))) /
                            (jac[i, j, k] + jac(i, j, k + 1))
                    rhowm = (jac(i, j, k - 1) * (var.rho[i, j, k] + rhostrattfc[i, j, k]) +
                             jac[i, j, k] * (var.rho(i, j, k - 1) + rhostrattfc(i, j, k - 1))) /
                            (jac[i, j, k] + jac(i, j, k - 1))
                    rho = var.rho[i, j, k] + rhostrattfc[i, j, k]
                    # Interpolate metric tensor elements.
                    met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 1, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    met13EdgeD = (jac(i, j, k - 1) * met(i, j, k, 1, 3) +
                                  jac[i, j, k] * met(i, j, k - 1, 1, 3)) /
                                 (jac[i, j, k] + jac(i, j, k - 1))
                    met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 2, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    met23EdgeD = (jac(i, j, k - 1) * met(i, j, k, 2, 3) +
                                  jac[i, j, k] * met(i, j, k - 1, 2, 3)) /
                                 (jac[i, j, k] + jac(i, j, k - 1))
                    met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 3, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    met33EdgeD = (jac(i, j, k - 1) * met(i, j, k, 3, 3) +
                                  jac[i, j, k] * met(i, j, k - 1, 3, 3)) /
                                 (jac[i, j, k] + jac(i, j, k - 1))
                    # Interpolate pressure differences.
                    dpREdgeU = (jac(i + 1, j, k + 1) * dp[i+1, j, k] +
                                jac[i+1, j, k] * dp(i + 1, j, k + 1)) /
                               (jac[i+1, j, k] + jac(i + 1, j, k + 1))
                    dpLEdgeU = (jac(i - 1, j, k + 1) * dp(i - 1, j, k) +
                                jac(i - 1, j, k) * dp(i - 1, j, k + 1)) /
                               (jac(i - 1, j, k) + jac(i - 1, j, k + 1))
                    dpREdgeD = (jac(i + 1, j, k - 1) * dp[i+1, j, k] +
                                jac[i+1, j, k] * dp(i + 1, j, k - 1)) /
                               (jac[i+1, j, k] + jac(i + 1, j, k - 1))
                    dpLEdgeD = (jac(i - 1, j, k - 1) * dp(i - 1, j, k) +
                                jac(i - 1, j, k) * dp(i - 1, j, k - 1)) /
                               (jac(i - 1, j, k) + jac(i - 1, j, k - 1))
                    dpFEdgeU = (jac(i, j + 1, k + 1) * dp(i, j + 1, k) +
                                jac(i, j + 1, k) * dp(i, j + 1, k + 1)) /
                               (jac(i, j + 1, k) + jac(i, j + 1, k + 1))
                    dpBEdgeU = (jac(i, j - 1, k + 1) * dp(i, j - 1, k) +
                                jac(i, j - 1, k) * dp(i, j - 1, k + 1)) /
                               (jac(i, j - 1, k) + jac(i, j - 1, k + 1))
                    dpFEdgeD = (jac(i, j + 1, k - 1) * dp(i, j + 1, k) +
                                jac(i, j + 1, k) * dp(i, j + 1, k - 1)) /
                               (jac(i, j + 1, k) + jac(i, j + 1, k - 1))
                    dpBEdgeD = (jac(i, j - 1, k - 1) * dp(i, j - 1, k) +
                                jac(i, j - 1, k) * dp(i, j - 1, k - 1)) /
                               (jac(i, j - 1, k) + jac(i, j - 1, k - 1))
                    # Compute pressure difference gradients.
                    pGradZEdgeU = kappaInv * MaInv2 * pEdgeU / rhow0 *
                                  (0.5 * met13EdgeU * (dpREdgeU - dpLEdgeU) / dx +
                                   0.5 * met23EdgeU * (dpFEdgeU - dpBEdgeU) / dy +
                                   met33EdgeU * (dp(i, j, k + 1) - dp[i, j, k]) / dz)
                    pGradZEdgeD = kappaInv * MaInv2 * pEdgeD / rhowm *
                                  (0.5 * met13EdgeD * (dpREdgeD - dpLEdgeD) / dx +
                                   0.5 * met23EdgeD * (dpFEdgeD - dpBEdgeD) / dy +
                                   met33EdgeD * (dp[i, j, k] - dp(i, j, k - 1)) / dz)
                    # Adjust at boundaries.
                    if (k == 1 && zBoundary == "solid_wall")
                        pGradZEdgeD = 0.0
                    elseif (k == nz && zBoundary == "solid_wall")
                        pGradZEdgeU = 0.0
                    end
                    # Interpolate.
                    pGradZ = 0.5 * (pGradZEdgeU + pGradZEdgeD)
                    # Compute buoyancy correction.
                    db = -1.0 /
                         (facw + rhostrattfc[i, j, k] / rho * bvsStrat[i, j, k] * dt^2.0) *
                         (-rhostrattfc[i, j, k] / rho *
                          bvsStrat[i, j, k] *
                          facprs *
                          dt^2.0 *
                          jac[i, j, k] *
                          pGradZ +
                          rhostrattfc[i, j, k] / rho *
                          bvsStrat[i, j, k] *
                          dt *
                          jac[i, j, k] *
                          facw *
                          0.5 *
                          (met(i, j, k, 1, 3) * (corX[i, j, k] + corX(i - 1, j, k)) +
                           met(i, j, k, 2, 3) * (corY[i, j, k] + corY(i, j - 1, k))))

                    var.rhop[i, j, k] = var.rhop[i, j, k] - rho / g_ndim * db
                end
            end
        end
    end

    if (opt == "impl")
        @. kr_sp_tfc = kr_sp_tfc / facray
        @. kr_sp_w_tfc = kr_sp_w_tfc / facray
    end
end

function val_PsIn(model, dt, opt, facray)
    (; nx, ny, nz) = model.domain
    (; dx, dy, dz, jac, met) = model.grid
    (; pstrattfc, rhostrattfc) = model.atmosphere
    var = model.variables.prognostic_fields


    # Poisson solver cache
    (; ac_b, acv_b, ach_b, al_b, ar_b, ab_b, af_b, ad_b, au_b, aru_b, ard_b, alu_b) = model.operator.op_fields
    (; ald_b, afu_b, afd_b, abu_b, abd_b, auu_b, add_b, aruu_b, ardd_b, aluu_b) = model.operator.op_fields
    (; aldd_b, afuu_b, afdd_b, abuu_b, abdd_b) = model.operator.op_fields
    # Calculates the matrix values for the pressure solver
    # The solver solves for dt * dp, hence no dt in the matrix elements

    # facray multiplies the Rayleigh-damping terms so that they are only
    # handled in the implicit time stepping (sponge and immersed boundary)

    # opt = "expl" =>
    # Pressure solver for explicit problem and corresponding correction
    # of the winds
    # opt = "impl" =>
    # Pressure solver for implicit problem and corresponding correction
    # of the winds and density fluctuations

    if (opt == "expl")
        # Compute tensor elements for TFC.
        for k in 1:nz
            for j in 1:ny
                for i in 1:nx
                    # Compute scaling factors.
                    fcscal = sqrt(pstrattfc[i, j, k]^2.0 / rhostrattfc[i, j, k])
                    fcscal_r = sqrt(pstrattfc[i+1, j, k]^2.0 / rhostrattfc[i+1, j, k])
                    fcscal_l = sqrt(pstrattfc(i - 1, j, k)^2.0 / rhostrattfc(i - 1, j, k))
                    fcscal_f = sqrt(pstrattfc(i, j + 1, k)^2.0 / rhostrattfc(i, j + 1, k))
                    fcscal_b = sqrt(pstrattfc(i, j - 1, k)^2.0 / rhostrattfc(i, j - 1, k))
                    fcscal_u = sqrt(pstrattfc(i, j, k + 1)^2.0 / rhostrattfc(i, j, k + 1))
                    fcscal_d = sqrt(pstrattfc(i, j, k - 1)^2.0 / rhostrattfc(i, j, k - 1))
                    fcscal_ru = sqrt(pstrattfc(i + 1, j, k + 1)^2.0 /
                                     rhostrattfc(i + 1, j, k + 1))
                    fcscal_rd = sqrt(pstrattfc(i + 1, j, k - 1)^2.0 /
                                     rhostrattfc(i + 1, j, k - 1))
                    fcscal_lu = sqrt(pstrattfc(i - 1, j, k + 1)^2.0 /
                                     rhostrattfc(i - 1, j, k + 1))
                    fcscal_ld = sqrt(pstrattfc(i - 1, j, k - 1)^2.0 /
                                     rhostrattfc(i - 1, j, k - 1))
                    fcscal_fu = sqrt(pstrattfc(i, j + 1, k + 1)^2.0 /
                                     rhostrattfc(i, j + 1, k + 1))
                    fcscal_fd = sqrt(pstrattfc(i, j + 1, k - 1)^2.0 /
                                     rhostrattfc(i, j + 1, k - 1))
                    fcscal_bu = sqrt(pstrattfc(i, j - 1, k + 1)^2.0 /
                                     rhostrattfc(i, j - 1, k + 1))
                    fcscal_bd = sqrt(pstrattfc(i, j - 1, k - 1)^2.0 /
                                     rhostrattfc(i, j - 1, k - 1))
                    fcscal_uu = sqrt(pstrattfc(i, j, k + 2)^2.0 / rhostrattfc(i, j, k + 2))
                    fcscal_dd = sqrt(pstrattfc(i, j, k - 2)^2.0 / rhostrattfc(i, j, k - 2))
                    fcscal_ruu = sqrt(pstrattfc(i + 1, j, k + 2)^2.0 /
                                      rhostrattfc(i + 1, j, k + 2))
                    fcscal_rdd = sqrt(pstrattfc(i + 1, j, k - 2)^2.0 /
                                      rhostrattfc(i + 1, j, k - 2))
                    fcscal_luu = sqrt(pstrattfc(i - 1, j, k + 2)^2.0 /
                                      rhostrattfc(i - 1, j, k + 2))
                    fcscal_ldd = sqrt(pstrattfc(i - 1, j, k - 2)^2.0 /
                                      rhostrattfc(i - 1, j, k - 2))
                    fcscal_fuu = sqrt(pstrattfc(i, j + 1, k + 2)^2.0 /
                                      rhostrattfc(i, j + 1, k + 2))
                    fcscal_fdd = sqrt(pstrattfc(i, j + 1, k - 2)^2.0 /
                                      rhostrattfc(i, j + 1, k - 2))
                    fcscal_buu = sqrt(pstrattfc(i, j - 1, k + 2)^2.0 /
                                      rhostrattfc(i, j - 1, k + 2))
                    fcscal_bdd = sqrt(pstrattfc(i, j - 1, k - 2)^2.0 /
                                      rhostrattfc(i, j - 1, k - 2))

                    # Compute inverse Jacobian.
                    jacInv = 1.0 / jac[i, j, k]

                    # Compute P coefficients (divergence).
                    pEdgeRDiv = 0.5 * (jac[i, j, k] * pstrattfc[i, j, k] +
                                       jac[i+1, j, k] * pstrattfc[i+1, j, k])
                    pEdgeLDiv = 0.5 * (jac[i, j, k] * pstrattfc[i, j, k] +
                                       jac(i - 1, j, k) * pstrattfc(i - 1, j, k))
                    pEdgeFDiv = 0.5 * (jac[i, j, k] * pstrattfc[i, j, k] +
                                       jac(i, j + 1, k) * pstrattfc(i, j + 1, k))
                    pEdgeBDiv = 0.5 * (jac[i, j, k] * pstrattfc[i, j, k] +
                                       jac(i, j - 1, k) * pstrattfc(i, j - 1, k))
                    pEdgeUDiv = jac[i, j, k] *
                                jac(i, j, k + 1) *
                                (pstrattfc[i, j, k] + pstrattfc(i, j, k + 1)) /
                                (jac[i, j, k] + jac(i, j, k + 1))
                    pEdgeDDiv = jac[i, j, k] *
                                jac(i, j, k - 1) *
                                (pstrattfc[i, j, k] + pstrattfc(i, j, k - 1)) /
                                (jac[i, j, k] + jac(i, j, k - 1))

                    # Compute P coefficients (pressure gradient).
                    pEdgeRGra = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i+1, j, k])
                    pEdgeLGra = 0.5 * (pstrattfc[i, j, k] + pstrattfc(i - 1, j, k))
                    pEdgeFGra = 0.5 * (pstrattfc[i, j, k] + pstrattfc(i, j + 1, k))
                    pEdgeBGra = 0.5 * (pstrattfc[i, j, k] + pstrattfc(i, j - 1, k))
                    pEdgeUGra = (jac(i, j, k + 1) * pstrattfc[i, j, k] +
                                 jac[i, j, k] * pstrattfc(i, j, k + 1)) /
                                (jac[i, j, k] + jac(i, j, k + 1))
                    pEdgeDGra = (jac(i, j, k - 1) * pstrattfc[i, j, k] +
                                 jac[i, j, k] * pstrattfc(i, j, k - 1)) /
                                (jac[i, j, k] + jac(i, j, k - 1))

                    # Compute density coefficients.
                    rhoEdgeR = 0.5 * (var.rho[i, j, k] +
                                      var.rho[i+1, j, k] +
                                      rhostrattfc[i, j, k] +
                                      rhostrattfc[i+1, j, k])
                    rhoEdgeL = 0.5 * (var.rho[i, j, k] +
                                      var.rho(i - 1, j, k) +
                                      rhostrattfc[i, j, k] +
                                      rhostrattfc(i - 1, j, k))
                    rhoEdgeF = 0.5 * (var.rho[i, j, k] +
                                      var.rho(i, j + 1, k) +
                                      rhostrattfc[i, j, k] +
                                      rhostrattfc(i, j + 1, k))
                    rhoEdgeB = 0.5 * (var.rho[i, j, k] +
                                      var.rho(i, j - 1, k) +
                                      rhostrattfc[i, j, k] +
                                      rhostrattfc(i, j - 1, k))
                    rhoEdgeU = (jac(i, j, k + 1) * (var.rho[i, j, k] + rhostrattfc[i, j, k]) +
                                jac[i, j, k] *
                                (var.rho(i, j, k + 1) + rhostrattfc(i, j, k + 1))) /
                               (jac[i, j, k] + jac(i, j, k + 1))
                    rhoEdgeD = (jac(i, j, k - 1) * (var.rho[i, j, k] + rhostrattfc[i, j, k]) +
                                jac[i, j, k] *
                                (var.rho(i, j, k - 1) + rhostrattfc(i, j, k - 1))) /
                               (jac[i, j, k] + jac(i, j, k - 1))

                    # Interpolate metric-tensor elements.
                    met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))
                    met13EdgeL = 0.5 * (met(i, j, k, 1, 3) + met(i - 1, j, k, 1, 3))
                    met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3))
                    met23EdgeB = 0.5 * (met(i, j, k, 2, 3) + met(i, j - 1, k, 2, 3))
                    met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 1, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 2, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 3, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    met13EdgeD = (jac(i, j, k - 1) * met(i, j, k, 1, 3) +
                                  jac[i, j, k] * met(i, j, k - 1, 1, 3)) /
                                 (jac[i, j, k] + jac(i, j, k - 1))
                    met23EdgeD = (jac(i, j, k - 1) * met(i, j, k, 2, 3) +
                                  jac[i, j, k] * met(i, j, k - 1, 2, 3)) /
                                 (jac[i, j, k] + jac(i, j, k - 1))
                    met33EdgeD = (jac(i, j, k - 1) * met(i, j, k, 3, 3) +
                                  jac[i, j, k] * met(i, j, k - 1, 3, 3)) /
                                 (jac[i, j, k] + jac(i, j, k - 1))

                    # --------------------- A(i,j,k) ---------------------
                    # zBoundary = "solid_wall" # TODO - FIX THIS TEMPORARY HACK
                    if (k == 1 && zBoundary == "solid_wall")
                        AC = -jacInv / dx * (pEdgeRDiv / rhoEdgeR *
                                             pEdgeRGra *
                                             (1.0 / dx + 0.75 * met13EdgeR / dz) +
                                             pEdgeLDiv / rhoEdgeL *
                                             pEdgeLGra *
                                             (1.0 / dx - 0.75 * met13EdgeL / dz)) -
                             jacInv / dy * (pEdgeFDiv / rhoEdgeF *
                                            pEdgeFGra *
                                            (1.0 / dy + 0.75 * met23EdgeF / dz) +
                                            pEdgeBDiv / rhoEdgeB *
                                            pEdgeBGra *
                                            (1.0 / dy - 0.75 * met23EdgeB / dz)) -
                             jacInv / dz * pEdgeUDiv / rhoEdgeU * pEdgeUGra * met33EdgeU /
                             dz
                    elseif (k == nz && zBoundary == "solid_wall")
                        AC = -jacInv / dx * (pEdgeRDiv / rhoEdgeR *
                                             pEdgeRGra *
                                             (1.0 / dx - 0.75 * met13EdgeR / dz) +
                                             pEdgeLDiv / rhoEdgeL *
                                             pEdgeLGra *
                                             (1.0 / dx + 0.75 * met13EdgeL / dz)) -
                             jacInv / dy * (pEdgeFDiv / rhoEdgeF *
                                            pEdgeFGra *
                                            (1.0 / dy - 0.75 * met23EdgeF / dz) +
                                            pEdgeBDiv / rhoEdgeB *
                                            pEdgeBGra *
                                            (1.0 / dy + 0.75 * met23EdgeB / dz)) -
                             jacInv / dz * pEdgeDDiv / rhoEdgeD * pEdgeDGra * met33EdgeD /
                             dz
                    else
                        AC = -jacInv / dx * (pEdgeRDiv / rhoEdgeR * pEdgeRGra / dx +
                                             pEdgeLDiv / rhoEdgeL * pEdgeLGra / dx) -
                             jacInv / dy * (pEdgeFDiv / rhoEdgeF * pEdgeFGra / dy +
                                            pEdgeBDiv / rhoEdgeB * pEdgeBGra / dy) -
                             jacInv / dz *
                             (pEdgeUDiv / rhoEdgeU * pEdgeUGra * met33EdgeU / dz +
                              pEdgeDDiv / rhoEdgeD * pEdgeDGra * met33EdgeD / dz)
                    end

                    # -------------------- A(i+1,j,k) --------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AR = jacInv / dx * pEdgeRDiv / rhoEdgeR *
                             pEdgeRGra *
                             (1.0 / dx - 0.75 * met13EdgeR / dz) +
                             jacInv / dz * pEdgeUDiv / rhoEdgeU *
                             pEdgeUGra *
                             0.5 *
                             met13EdgeU / dx * jac(i + 1, j, k + 1) /
                             (jac[i+1, j, k] + jac(i + 1, j, k + 1))
                    elseif (k == nz && zBoundary == "solid_wall")
                        AR = jacInv / dx * pEdgeRDiv / rhoEdgeR *
                             pEdgeRGra *
                             (1.0 / dx + 0.75 * met13EdgeR / dz) -
                             jacInv / dz * pEdgeDDiv / rhoEdgeD *
                             pEdgeDGra *
                             0.5 *
                             met13EdgeD / dx * jac(i + 1, j, k - 1) /
                             (jac[i+1, j, k] + jac(i + 1, j, k - 1))
                    else
                        AR = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra / dx +
                             jacInv / dz *
                             (pEdgeUDiv / rhoEdgeU * pEdgeUGra * met13EdgeU * 0.5 / dx *
                              jac(i + 1, j, k + 1) /
                              (jac[i+1, j, k] + jac(i + 1, j, k + 1)) -
                              pEdgeDDiv / rhoEdgeD * pEdgeDGra * met13EdgeD * 0.5 / dx *
                              jac(i + 1, j, k - 1) /
                              (jac[i+1, j, k] + jac(i + 1, j, k - 1)))
                    end

                    # -------------------- A(i-1,j,k) --------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AL = jacInv / dx * pEdgeLDiv / rhoEdgeL *
                             pEdgeLGra *
                             (1.0 / dx + 0.75 * met13EdgeL / dz) -
                             jacInv / dz * pEdgeUDiv / rhoEdgeU *
                             pEdgeUGra *
                             0.5 *
                             met13EdgeU / dx * jac(i - 1, j, k + 1) /
                             (jac(i - 1, j, k) + jac(i - 1, j, k + 1))
                    elseif (k == nz && zBoundary == "solid_wall")
                        AL = jacInv / dx * pEdgeLDiv / rhoEdgeL *
                             pEdgeLGra *
                             (1.0 / dx - 0.75 * met13EdgeL / dz) +
                             jacInv / dz * pEdgeDDiv / rhoEdgeD *
                             pEdgeDGra *
                             0.5 *
                             met13EdgeD / dx * jac(i - 1, j, k - 1) /
                             (jac(i - 1, j, k) + jac(i - 1, j, k - 1))
                    else
                        AL = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra / dx -
                             jacInv / dz *
                             (pEdgeUDiv / rhoEdgeU * pEdgeUGra * met13EdgeU * 0.5 / dx *
                              jac(i - 1, j, k + 1) /
                              (jac(i - 1, j, k) + jac(i - 1, j, k + 1)) -
                              pEdgeDDiv / rhoEdgeD * pEdgeDGra * met13EdgeD * 0.5 / dx *
                              jac(i - 1, j, k - 1) /
                              (jac(i - 1, j, k) + jac(i - 1, j, k - 1)))
                    end

                    # -------------------- A(i,j+1,k) --------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AF = jacInv / dy * pEdgeFDiv / rhoEdgeF *
                             pEdgeFGra *
                             (1.0 / dy - 0.75 * met23EdgeF / dz) +
                             jacInv / dz * pEdgeUDiv / rhoEdgeU *
                             pEdgeUGra *
                             0.5 *
                             met23EdgeU / dy * jac(i, j + 1, k + 1) /
                             (jac(i, j + 1, k) + jac(i, j + 1, k + 1))
                    elseif (k == nz && zBoundary == "solid_wall")
                        AF = jacInv / dy * pEdgeFDiv / rhoEdgeF *
                             pEdgeFGra *
                             (1.0 / dy + 0.75 * met23EdgeF / dz) -
                             jacInv / dz * pEdgeDDiv / rhoEdgeD *
                             pEdgeDGra *
                             0.5 *
                             met23EdgeD / dy * jac(i, j + 1, k - 1) /
                             (jac(i, j + 1, k) + jac(i, j + 1, k - 1))
                    else
                        AF = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra / dy +
                             jacInv / dz *
                             (pEdgeUDiv / rhoEdgeU * pEdgeUGra * met23EdgeU * 0.5 / dy *
                              jac(i, j + 1, k + 1) /
                              (jac(i, j + 1, k) + jac(i, j + 1, k + 1)) -
                              pEdgeDDiv / rhoEdgeD * pEdgeDGra * met23EdgeD * 0.5 / dy *
                              jac(i, j + 1, k - 1) /
                              (jac(i, j + 1, k) + jac(i, j + 1, k - 1)))
                    end

                    # -------------------- A(i,j-1,k) --------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AB = jacInv / dy * pEdgeBDiv / rhoEdgeB *
                             pEdgeBGra *
                             (1.0 / dy + 0.75 * met23EdgeB / dz) -
                             jacInv / dz * pEdgeUDiv / rhoEdgeU *
                             pEdgeUGra *
                             0.5 *
                             met23EdgeU / dy * jac(i, j - 1, k + 1) /
                             (jac(i, j - 1, k) + jac(i, j - 1, k + 1))
                    elseif (k == nz && zBoundary == "solid_wall")
                        AB = jacInv / dy * pEdgeBDiv / rhoEdgeB *
                             pEdgeBGra *
                             (1.0 / dy - 0.75 * met23EdgeB / dz) +
                             jacInv / dz * pEdgeDDiv / rhoEdgeD *
                             pEdgeDGra *
                             0.5 *
                             met23EdgeD / dy * jac(i, j - 1, k - 1) /
                             (jac(i, j - 1, k) + jac(i, j - 1, k - 1))
                    else
                        AB = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra / dy -
                             jacInv / dz *
                             (pEdgeUDiv / rhoEdgeU * pEdgeUGra * met23EdgeU * 0.5 / dy *
                              jac(i, j - 1, k + 1) /
                              (jac(i, j - 1, k) + jac(i, j - 1, k + 1)) -
                              pEdgeDDiv / rhoEdgeD * pEdgeDGra * met23EdgeD * 0.5 / dy *
                              jac(i, j - 1, k - 1) /
                              (jac(i, j - 1, k) + jac(i, j - 1, k - 1)))
                    end

                    # -------------------- A(i,j,k+1) --------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AU = jacInv / dx *
                             (pEdgeRDiv / rhoEdgeR * pEdgeRGra * met13EdgeR / dz -
                              pEdgeLDiv / rhoEdgeL * pEdgeLGra * met13EdgeL / dz) +
                             jacInv / dy *
                             (pEdgeFDiv / rhoEdgeF * pEdgeFGra * met23EdgeF / dz -
                              pEdgeBDiv / rhoEdgeB * pEdgeBGra * met23EdgeB / dz) +
                             jacInv / dz * pEdgeUDiv / rhoEdgeU * pEdgeUGra * met33EdgeU /
                             dz
                    elseif (k == nz && zBoundary == "solid_wall")
                        AU = 0.0
                    else
                        AU = jacInv / dx *
                             (pEdgeRDiv / rhoEdgeR * pEdgeRGra * met13EdgeR * 0.25 / dz -
                              pEdgeLDiv / rhoEdgeL * pEdgeLGra * met13EdgeL * 0.25 / dz) +
                             jacInv / dy *
                             (pEdgeFDiv / rhoEdgeF * pEdgeFGra * met23EdgeF * 0.25 / dz -
                              pEdgeBDiv / rhoEdgeB * pEdgeBGra * met23EdgeB * 0.25 / dz) +
                             jacInv / dz * pEdgeUDiv / rhoEdgeU * pEdgeUGra * met33EdgeU /
                             dz
                    end

                    # -------------------- A(i,j,k-1) --------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AD = 0.0
                    elseif (k == nz && zBoundary == "solid_wall")
                        AD = -jacInv / dx *
                             (pEdgeRDiv / rhoEdgeR * pEdgeRGra * met13EdgeR / dz -
                              pEdgeLDiv / rhoEdgeL * pEdgeLGra * met13EdgeL / dz) -
                             jacInv / dy *
                             (pEdgeFDiv / rhoEdgeF * pEdgeFGra * met23EdgeF / dz -
                              pEdgeBDiv / rhoEdgeB * pEdgeBGra * met23EdgeB / dz) +
                             jacInv / dz * pEdgeDDiv / rhoEdgeD * pEdgeDGra * met33EdgeD /
                             dz
                    else
                        AD = -jacInv / dx *
                             (pEdgeRDiv / rhoEdgeR * pEdgeRGra * met13EdgeR * 0.25 / dz -
                              pEdgeLDiv / rhoEdgeL * pEdgeLGra * met13EdgeL * 0.25 / dz) -
                             jacInv / dy *
                             (pEdgeFDiv / rhoEdgeF * pEdgeFGra * met23EdgeF * 0.25 / dz -
                              pEdgeBDiv / rhoEdgeB * pEdgeBGra * met23EdgeB * 0.25 / dz) +
                             jacInv / dz * pEdgeDDiv / rhoEdgeD * pEdgeDGra * met33EdgeD /
                             dz
                    end

                    # ------------------- A(i+1,j,k+1) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ARU = jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * met13EdgeR /
                              dz +
                              jacInv / dz * pEdgeUDiv / rhoEdgeU *
                              pEdgeUGra *
                              met13EdgeU *
                              0.5 / dx * jac[i+1, j, k] /
                              (jac[i+1, j, k] + jac(i + 1, j, k + 1))
                    elseif (k == nz && zBoundary == "solid_wall")
                        ARU = 0.0
                    else
                        ARU = jacInv / dx * pEdgeRDiv / rhoEdgeR *
                              pEdgeRGra *
                              met13EdgeR *
                              0.25 / dz +
                              jacInv / dz * pEdgeUDiv / rhoEdgeU *
                              pEdgeUGra *
                              met13EdgeU *
                              0.5 / dx * jac[i+1, j, k] /
                              (jac[i+1, j, k] + jac(i + 1, j, k + 1))
                    end

                    # ------------------- A(i+1,j,k-1) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ARD = 0.0
                    elseif (k == nz && zBoundary == "solid_wall")
                        ARD = -jacInv / dx * pEdgeRDiv / rhoEdgeR * pEdgeRGra * met13EdgeR /
                              dz -
                              jacInv / dz * pEdgeDDiv / rhoEdgeD *
                              pEdgeDGra *
                              met13EdgeD *
                              0.5 / dx * jac[i+1, j, k] /
                              (jac[i+1, j, k] + jac(i + 1, j, k - 1))
                    else
                        ARD = -jacInv / dx * pEdgeRDiv / rhoEdgeR *
                              pEdgeRGra *
                              met13EdgeR *
                              0.25 / dz -
                              jacInv / dz * pEdgeDDiv / rhoEdgeD *
                              pEdgeDGra *
                              met13EdgeD *
                              0.5 / dx * jac[i+1, j, k] /
                              (jac[i+1, j, k] + jac(i + 1, j, k - 1))
                    end

                    # ------------------- A(i-1,j,k+1) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ALU = -jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * met13EdgeL /
                              dz -
                              jacInv / dz * pEdgeUDiv / rhoEdgeU *
                              pEdgeUGra *
                              met13EdgeU *
                              0.5 / dx * jac(i - 1, j, k) /
                              (jac(i - 1, j, k) + jac(i - 1, j, k + 1))
                    elseif (k == nz && zBoundary == "solid_wall")
                        ALU = 0.0
                    else
                        ALU = -jacInv / dx * pEdgeLDiv / rhoEdgeL *
                              pEdgeLGra *
                              met13EdgeL *
                              0.25 / dz -
                              jacInv / dz * pEdgeUDiv / rhoEdgeU *
                              pEdgeUGra *
                              met13EdgeU *
                              0.5 / dx * jac(i - 1, j, k) /
                              (jac(i - 1, j, k) + jac(i - 1, j, k + 1))
                    end

                    # ------------------- A(i-1,j,k-1) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ALD = 0.0
                    elseif (k == nz && zBoundary == "solid_wall")
                        ALD = jacInv / dx * pEdgeLDiv / rhoEdgeL * pEdgeLGra * met13EdgeL /
                              dz +
                              jacInv / dz * pEdgeDDiv / rhoEdgeD *
                              pEdgeDGra *
                              met13EdgeD *
                              0.5 / dx * jac(i - 1, j, k) /
                              (jac(i - 1, j, k) + jac(i - 1, j, k - 1))
                    else
                        ALD = jacInv / dx * pEdgeLDiv / rhoEdgeL *
                              pEdgeLGra *
                              met13EdgeL *
                              0.25 / dz +
                              jacInv / dz * pEdgeDDiv / rhoEdgeD *
                              pEdgeDGra *
                              met13EdgeD *
                              0.5 / dx * jac(i - 1, j, k) /
                              (jac(i - 1, j, k) + jac(i - 1, j, k - 1))
                    end

                    # ------------------- A(i,j+1,k+1) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AFU = jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * met23EdgeF /
                              dz +
                              jacInv / dz * pEdgeUDiv / rhoEdgeU *
                              pEdgeUGra *
                              met23EdgeU *
                              0.5 / dy * jac(i, j + 1, k) /
                              (jac(i, j + 1, k) + jac(i, j + 1, k + 1))
                    elseif (k == nz && zBoundary == "solid_wall")
                        AFU = 0.0
                    else
                        AFU = jacInv / dy * pEdgeFDiv / rhoEdgeF *
                              pEdgeFGra *
                              met23EdgeF *
                              0.25 / dz +
                              jacInv / dz * pEdgeUDiv / rhoEdgeU *
                              pEdgeUGra *
                              met23EdgeU *
                              0.5 / dy * jac(i, j + 1, k) /
                              (jac(i, j + 1, k) + jac(i, j + 1, k + 1))
                    end

                    # ------------------- A(i,j+1,k-1) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AFD = 0.0
                    elseif (k == nz && zBoundary == "solid_wall")
                        AFD = -jacInv / dy * pEdgeFDiv / rhoEdgeF * pEdgeFGra * met23EdgeF /
                              dz -
                              jacInv / dz * pEdgeDDiv / rhoEdgeD *
                              pEdgeDGra *
                              met23EdgeD *
                              0.5 / dy * jac(i, j + 1, k) /
                              (jac(i, j + 1, k) + jac(i, j + 1, k - 1))
                    else
                        AFD = -jacInv / dy * pEdgeFDiv / rhoEdgeF *
                              pEdgeFGra *
                              met23EdgeF *
                              0.25 / dz -
                              jacInv / dz * pEdgeDDiv / rhoEdgeD *
                              pEdgeDGra *
                              met23EdgeD *
                              0.5 / dy * jac(i, j + 1, k) /
                              (jac(i, j + 1, k) + jac(i, j + 1, k - 1))
                    end

                    # ------------------- A(i,j-1,k+1) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ABU = -jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * met23EdgeB /
                              dz -
                              jacInv / dz * pEdgeUDiv / rhoEdgeU *
                              pEdgeUGra *
                              met23EdgeU *
                              0.5 / dy * jac(i, j - 1, k) /
                              (jac(i, j - 1, k) + jac(i, j - 1, k + 1))
                    elseif (k == nz && zBoundary == "solid_wall")
                        ABU = 0.0
                    else
                        ABU = -jacInv / dy * pEdgeBDiv / rhoEdgeB *
                              pEdgeBGra *
                              met23EdgeB *
                              0.25 / dz -
                              jacInv / dz * pEdgeUDiv / rhoEdgeU *
                              pEdgeUGra *
                              met23EdgeU *
                              0.5 / dy * jac(i, j - 1, k) /
                              (jac(i, j - 1, k) + jac(i, j - 1, k + 1))
                    end

                    # ------------------- A(i,j-1,k-1) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ABD = 0.0
                    elseif (k == nz && zBoundary == "solid_wall")
                        ABD = jacInv / dy * pEdgeBDiv / rhoEdgeB * pEdgeBGra * met23EdgeB /
                              dz +
                              jacInv / dz * pEdgeDDiv / rhoEdgeD *
                              pEdgeDGra *
                              met23EdgeD *
                              0.5 / dy * jac(i, j - 1, k) /
                              (jac(i, j - 1, k) + jac(i, j - 1, k - 1))
                    else
                        ABD = jacInv / dy * pEdgeBDiv / rhoEdgeB *
                              pEdgeBGra *
                              met23EdgeB *
                              0.25 / dz +
                              jacInv / dz * pEdgeDDiv / rhoEdgeD *
                              pEdgeDGra *
                              met23EdgeD *
                              0.5 / dy * jac(i, j - 1, k) /
                              (jac(i, j - 1, k) + jac(i, j - 1, k - 1))
                    end

                    # ------------------- A(i,j,k+2) ---------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AUU = -jacInv / dx *
                              (pEdgeRDiv / rhoEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz -
                               pEdgeLDiv / rhoEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz) -
                              jacInv / dy *
                              (pEdgeFDiv / rhoEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz -
                               pEdgeBDiv / rhoEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz)
                    else
                        AUU = 0.0
                    end

                    # ------------------- A(i,j,k-2) ---------------------

                    if (k == nz && zBoundary == "solid_wall")
                        ADD = jacInv / dx *
                              (pEdgeRDiv / rhoEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz -
                               pEdgeLDiv / rhoEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz) +
                              jacInv / dy *
                              (pEdgeFDiv / rhoEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz -
                               pEdgeBDiv / rhoEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz)
                    else
                        ADD = 0.0
                    end

                    # ------------------ A(i+1,j,k+2) --------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ARUU = -jacInv / dx * pEdgeRDiv / rhoEdgeR *
                               pEdgeRGra *
                               0.25 *
                               met13EdgeR / dz
                    else
                        ARUU = 0.0
                    end

                    # ------------------ A(i+1,j,k-2) --------------------

                    if (k == nz && zBoundary == "solid_wall")
                        ARDD = jacInv / dx * pEdgeRDiv / rhoEdgeR *
                               pEdgeRGra *
                               0.25 *
                               met13EdgeR / dz
                    else
                        ARDD = 0.0
                    end

                    # ------------------ A(i-1,j,k+2) --------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ALUU = jacInv / dx * pEdgeLDiv / rhoEdgeL *
                               pEdgeLGra *
                               0.25 *
                               met13EdgeL / dz
                    else
                        ALUU = 0.0
                    end

                    # ------------------ A(i-1,j,k-2) --------------------

                    if (k == nz && zBoundary == "solid_wall")
                        ALDD = -jacInv / dx * pEdgeLDiv / rhoEdgeL *
                               pEdgeLGra *
                               0.25 *
                               met13EdgeL / dz
                    else
                        ALDD = 0.0
                    end

                    # ------------------ A(i,j+1,k+2) --------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AFUU = -jacInv / dy * pEdgeFDiv / rhoEdgeF *
                               pEdgeFGra *
                               0.25 *
                               met23EdgeF / dz
                    else
                        AFUU = 0.0
                    end

                    # ------------------ A(i,j+1,k-2) --------------------

                    if (k == nz && zBoundary == "solid_wall")
                        AFDD = jacInv / dy * pEdgeFDiv / rhoEdgeF *
                               pEdgeFGra *
                               0.25 *
                               met23EdgeF / dz
                    else
                        AFDD = 0.0
                    end

                    # ------------------ A(i,j-1,k+2) --------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ABUU = jacInv / dy * pEdgeBDiv / rhoEdgeB *
                               pEdgeBGra *
                               0.25 *
                               met23EdgeB / dz
                    else
                        ABUU = 0.0
                    end

                    # ------------------ A(i,j-1,k-2) --------------------

                    if (k == nz && zBoundary == "solid_wall")
                        ABDD = -jacInv / dy * pEdgeBDiv / rhoEdgeB *
                               pEdgeBGra *
                               0.25 *
                               met23EdgeB / dz
                    else
                        ABDD = 0.0
                    end

                    # Scale the tensor elements.
                    AC = AC / (fcscal^2.0)
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

                    # Set tensor elements for bicgstab.
                    ac_b[i, j, k] = AC
                    ar_b[i, j, k] = AR
                    al_b[i, j, k] = AL
                    af_b[i, j, k] = AF
                    ab_b[i, j, k] = AB
                    au_b[i, j, k] = AU
                    ad_b[i, j, k] = AD
                    aru_b[i, j, k] = ARU
                    ard_b[i, j, k] = ARD
                    alu_b[i, j, k] = ALU
                    ald_b[i, j, k] = ALD
                    afu_b[i, j, k] = AFU
                    afd_b[i, j, k] = AFD
                    abu_b[i, j, k] = ABU
                    abd_b[i, j, k] = ABD
                    auu_b[i, j, k] = AUU
                    add_b[i, j, k] = ADD
                    aruu_b[i, j, k] = ARUU
                    ardd_b[i, j, k] = ARDD
                    aluu_b[i, j, k] = ALUU
                    aldd_b[i, j, k] = ALDD
                    afuu_b[i, j, k] = AFUU
                    afdd_b[i, j, k] = AFDD
                    abuu_b[i, j, k] = ABUU
                    abdd_b[i, j, k] = ABDD

                    # Store horizontal and vertical components of AC (for
                    # preconditioner).
                    # # TODO: types
                    if (model.parameters.poisson.preconditioner == "yes")
                        ach_b[i, j, k] = -AR - AL - AF - AB
                        acv_b[i, j, k] = -AU - AD
                    end
                end
            end
        end
    elseif (opt == "impl")

        # Compute tensor elements for TFC.
        @. kr_sp_tfc = kr_sp_tfc * facray
        @. kr_sp_w_tfc = kr_sp_w_tfc * facray
        for k in 1:nz
            for j in 1:ny
                for i in 1:nx
                    # Compute scaling factors.
                    fcscal = sqrt(pstrattfc[i, j, k]^2.0 / rhostrattfc[i, j, k])
                    fcscal_r = sqrt(pstrattfc[i+1, j, k]^2.0 / rhostrattfc[i+1, j, k])
                    fcscal_l = sqrt(pstrattfc(i - 1, j, k)^2.0 / rhostrattfc(i - 1, j, k))
                    fcscal_f = sqrt(pstrattfc(i, j + 1, k)^2.0 / rhostrattfc(i, j + 1, k))
                    fcscal_b = sqrt(pstrattfc(i, j - 1, k)^2.0 / rhostrattfc(i, j - 1, k))
                    fcscal_u = sqrt(pstrattfc(i, j, k + 1)^2.0 / rhostrattfc(i, j, k + 1))
                    fcscal_d = sqrt(pstrattfc(i, j, k - 1)^2.0 / rhostrattfc(i, j, k - 1))
                    fcscal_ru = sqrt(pstrattfc(i + 1, j, k + 1)^2.0 /
                                     rhostrattfc(i + 1, j, k + 1))
                    fcscal_rd = sqrt(pstrattfc(i + 1, j, k - 1)^2.0 /
                                     rhostrattfc(i + 1, j, k - 1))
                    fcscal_lu = sqrt(pstrattfc(i - 1, j, k + 1)^2.0 /
                                     rhostrattfc(i - 1, j, k + 1))
                    fcscal_ld = sqrt(pstrattfc(i - 1, j, k - 1)^2.0 /
                                     rhostrattfc(i - 1, j, k - 1))
                    fcscal_fu = sqrt(pstrattfc(i, j + 1, k + 1)^2.0 /
                                     rhostrattfc(i, j + 1, k + 1))
                    fcscal_fd = sqrt(pstrattfc(i, j + 1, k - 1)^2.0 /
                                     rhostrattfc(i, j + 1, k - 1))
                    fcscal_bu = sqrt(pstrattfc(i, j - 1, k + 1)^2.0 /
                                     rhostrattfc(i, j - 1, k + 1))
                    fcscal_bd = sqrt(pstrattfc(i, j - 1, k - 1)^2.0 /
                                     rhostrattfc(i, j - 1, k - 1))
                    fcscal_uu = sqrt(pstrattfc(i, j, k + 2)^2.0 / rhostrattfc(i, j, k + 2))
                    fcscal_dd = sqrt(pstrattfc(i, j, k - 2)^2.0 / rhostrattfc(i, j, k - 2))
                    fcscal_ruu = sqrt(pstrattfc(i + 1, j, k + 2)^2.0 /
                                      rhostrattfc(i + 1, j, k + 2))
                    fcscal_rdd = sqrt(pstrattfc(i + 1, j, k - 2)^2.0 /
                                      rhostrattfc(i + 1, j, k - 2))
                    fcscal_luu = sqrt(pstrattfc(i - 1, j, k + 2)^2.0 /
                                      rhostrattfc(i - 1, j, k + 2))
                    fcscal_ldd = sqrt(pstrattfc(i - 1, j, k - 2)^2.0 /
                                      rhostrattfc(i - 1, j, k - 2))
                    fcscal_fuu = sqrt(pstrattfc(i, j + 1, k + 2)^2.0 /
                                      rhostrattfc(i, j + 1, k + 2))
                    fcscal_fdd = sqrt(pstrattfc(i, j + 1, k - 2)^2.0 /
                                      rhostrattfc(i, j + 1, k - 2))
                    fcscal_buu = sqrt(pstrattfc(i, j - 1, k + 2)^2.0 /
                                      rhostrattfc(i, j - 1, k + 2))
                    fcscal_bdd = sqrt(pstrattfc(i, j - 1, k - 2)^2.0 /
                                      rhostrattfc(i, j - 1, k - 2))

                    # Compute inverse Jacobian.
                    jacInv = 1.0 / jac[i, j, k]

                    # Compute P coefficients (divergence).
                    pEdgeRDiv = 0.5 * (jac[i, j, k] * pstrattfc[i, j, k] +
                                       jac[i+1, j, k] * pstrattfc[i+1, j, k])
                    pEdgeLDiv = 0.5 * (jac[i, j, k] * pstrattfc[i, j, k] +
                                       jac(i - 1, j, k) * pstrattfc(i - 1, j, k))
                    pEdgeFDiv = 0.5 * (jac[i, j, k] * pstrattfc[i, j, k] +
                                       jac(i, j + 1, k) * pstrattfc(i, j + 1, k))
                    pEdgeBDiv = 0.5 * (jac[i, j, k] * pstrattfc[i, j, k] +
                                       jac(i, j - 1, k) * pstrattfc(i, j - 1, k))
                    pEdgeUDiv = jac[i, j, k] *
                                jac(i, j, k + 1) *
                                (pstrattfc[i, j, k] + pstrattfc(i, j, k + 1)) /
                                (jac[i, j, k] + jac(i, j, k + 1))
                    pEdgeDDiv = jac[i, j, k] *
                                jac(i, j, k - 1) *
                                (pstrattfc[i, j, k] + pstrattfc(i, j, k - 1)) /
                                (jac[i, j, k] + jac(i, j, k - 1))

                    # Compute P coefficients (pressure gradient).
                    pEdgeRGra = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i+1, j, k])
                    pEdgeLGra = 0.5 * (pstrattfc[i, j, k] + pstrattfc(i - 1, j, k))
                    pEdgeFGra = 0.5 * (pstrattfc[i, j, k] + pstrattfc(i, j + 1, k))
                    pEdgeBGra = 0.5 * (pstrattfc[i, j, k] + pstrattfc(i, j - 1, k))
                    pEdgeUGra = (jac(i, j, k + 1) * pstrattfc[i, j, k] +
                                 jac[i, j, k] * pstrattfc(i, j, k + 1)) /
                                (jac[i, j, k] + jac(i, j, k + 1))
                    pEdgeDGra = (jac(i, j, k - 1) * pstrattfc[i, j, k] +
                                 jac[i, j, k] * pstrattfc(i, j, k - 1)) /
                                (jac[i, j, k] + jac(i, j, k - 1))
                    pUEdgeRGra = 0.5 * (pstrattfc(i, j, k + 1) + pstrattfc(i + 1, j, k + 1))
                    pUEdgeLGra = 0.5 * (pstrattfc(i, j, k + 1) + pstrattfc(i - 1, j, k + 1))
                    pUEdgeFGra = 0.5 * (pstrattfc(i, j, k + 1) + pstrattfc(i, j + 1, k + 1))
                    pUEdgeBGra = 0.5 * (pstrattfc(i, j, k + 1) + pstrattfc(i, j - 1, k + 1))
                    pDEdgeRGra = 0.5 * (pstrattfc(i, j, k - 1) + pstrattfc(i + 1, j, k - 1))
                    pDEdgeLGra = 0.5 * (pstrattfc(i, j, k - 1) + pstrattfc(i - 1, j, k - 1))
                    pDEdgeFGra = 0.5 * (pstrattfc(i, j, k - 1) + pstrattfc(i, j + 1, k - 1))
                    pDEdgeBGra = 0.5 * (pstrattfc(i, j, k - 1) + pstrattfc(i, j - 1, k - 1))

                    # Compute density coefficients.
                    rhostrattfcEdgeR = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc[i+1, j, k])
                    rhostrattfcEdgeL = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc(i - 1, j, k))
                    rhostrattfcEdgeF = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc(i, j + 1, k))
                    rhostrattfcEdgeB = 0.5 * (rhostrattfc[i, j, k] + rhostrattfc(i, j - 1, k))
                    rhostrattfcEdgeU = (jac(i, j, k + 1) * rhostrattfc[i, j, k] +
                                        jac[i, j, k] * rhostrattfc(i, j, k + 1)) /
                                       (jac[i, j, k] + jac(i, j, k + 1))
                    rhostrattfcEdgeD = (jac(i, j, k - 1) * rhostrattfc[i, j, k] +
                                        jac[i, j, k] * rhostrattfc(i, j, k - 1)) /
                                       (jac[i, j, k] + jac(i, j, k - 1))
                    rhoEdgeR = 0.5 * (var.rho[i, j, k] + var.rho[i+1, j, k]) +
                               rhostrattfcEdgeR
                    rhoEdgeL = 0.5 * (var.rho[i, j, k] + var.rho(i - 1, j, k)) +
                               rhostrattfcEdgeL
                    rhoEdgeF = 0.5 * (var.rho[i, j, k] + var.rho(i, j + 1, k)) +
                               rhostrattfcEdgeF
                    rhoEdgeB = 0.5 * (var.rho[i, j, k] + var.rho(i, j - 1, k)) +
                               rhostrattfcEdgeB
                    rhoEdgeU = (jac(i, j, k + 1) * var.rho[i, j, k] +
                                jac[i, j, k] * var.rho(i, j, k + 1)) /
                               (jac[i, j, k] + jac(i, j, k + 1)) + rhostrattfcEdgeU
                    rhoEdgeD = (jac(i, j, k - 1) * var.rho[i, j, k] +
                                jac[i, j, k] * var.rho(i, j, k - 1)) /
                               (jac[i, j, k] + jac(i, j, k - 1)) + rhostrattfcEdgeD

                    rhoUEdgeR = 0.5 * (var.rho(i, j, k + 1) +
                                       var.rho(i + 1, j, k + 1) +
                                       rhostrattfc(i, j, k + 1) +
                                       rhostrattfc(i + 1, j, k + 1))
                    rhoUEdgeL = 0.5 * (var.rho(i, j, k + 1) +
                                       var.rho(i - 1, j, k + 1) +
                                       rhostrattfc(i, j, k + 1) +
                                       rhostrattfc(i - 1, j, k + 1))
                    rhoUEdgeF = 0.5 * (var.rho(i, j, k + 1) +
                                       var.rho(i, j + 1, k + 1) +
                                       rhostrattfc(i, j, k + 1) +
                                       rhostrattfc(i, j + 1, k + 1))
                    rhoUEdgeB = 0.5 * (var.rho(i, j, k + 1) +
                                       var.rho(i, j - 1, k + 1) +
                                       rhostrattfc(i, j, k + 1) +
                                       rhostrattfc(i, j - 1, k + 1))
                    rhoDEdgeR = 0.5 * (var.rho(i, j, k - 1) +
                                       var.rho(i + 1, j, k - 1) +
                                       rhostrattfc(i, j, k - 1) +
                                       rhostrattfc(i + 1, j, k - 1))
                    rhoDEdgeL = 0.5 * (var.rho(i, j, k - 1) +
                                       var.rho(i - 1, j, k - 1) +
                                       rhostrattfc(i, j, k - 1) +
                                       rhostrattfc(i - 1, j, k - 1))
                    rhoDEdgeF = 0.5 * (var.rho(i, j, k - 1) +
                                       var.rho(i, j + 1, k - 1) +
                                       rhostrattfc(i, j, k - 1) +
                                       rhostrattfc(i, j + 1, k - 1))
                    rhoDEdgeB = 0.5 * (var.rho(i, j, k - 1) +
                                       var.rho(i, j - 1, k - 1) +
                                       rhostrattfc(i, j, k - 1) +
                                       rhostrattfc(i, j - 1, k - 1))

                    # Compute squared buoyancy frequency at edges.
                    bvsStratEdgeU = (jac(i, j, k + 1) * bvsStrat[i, j, k] +
                                     jac[i, j, k] * bvsStrat(i, j, k + 1)) /
                                    (jac[i, j, k] + jac(i, j, k + 1))
                    bvsStratEdgeD = (jac(i, j, k - 1) * bvsStrat[i, j, k] +
                                     jac[i, j, k] * bvsStrat(i, j, k - 1)) /
                                    (jac[i, j, k] + jac(i, j, k - 1))

                    # Interpolate metric-tensor elements.
                    met13EdgeR = 0.5 * (met(i, j, k, 1, 3) + met(i + 1, j, k, 1, 3))
                    met13EdgeL = 0.5 * (met(i, j, k, 1, 3) + met(i - 1, j, k, 1, 3))
                    met23EdgeF = 0.5 * (met(i, j, k, 2, 3) + met(i, j + 1, k, 2, 3))
                    met23EdgeB = 0.5 * (met(i, j, k, 2, 3) + met(i, j - 1, k, 2, 3))
                    met13EdgeU = (jac(i, j, k + 1) * met(i, j, k, 1, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 1, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    met23EdgeU = (jac(i, j, k + 1) * met(i, j, k, 2, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 2, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    met33EdgeU = (jac(i, j, k + 1) * met(i, j, k, 3, 3) +
                                  jac[i, j, k] * met(i, j, k + 1, 3, 3)) /
                                 (jac[i, j, k] + jac(i, j, k + 1))
                    met13EdgeD = (jac(i, j, k - 1) * met(i, j, k, 1, 3) +
                                  jac[i, j, k] * met(i, j, k - 1, 1, 3)) /
                                 (jac[i, j, k] + jac(i, j, k - 1))
                    met23EdgeD = (jac(i, j, k - 1) * met(i, j, k, 2, 3) +
                                  jac[i, j, k] * met(i, j, k - 1, 2, 3)) /
                                 (jac[i, j, k] + jac(i, j, k - 1))
                    met33EdgeD = (jac(i, j, k - 1) * met(i, j, k, 3, 3) +
                                  jac[i, j, k] * met(i, j, k - 1, 3, 3)) /
                                 (jac[i, j, k] + jac(i, j, k - 1))
                    met13UEdgeR = 0.5 *
                                  (met(i, j, k + 1, 1, 3) + met(i + 1, j, k + 1, 1, 3))
                    met13UEdgeL = 0.5 *
                                  (met(i, j, k + 1, 1, 3) + met(i - 1, j, k + 1, 1, 3))
                    met23UEdgeF = 0.5 *
                                  (met(i, j, k + 1, 2, 3) + met(i, j + 1, k + 1, 2, 3))
                    met23UEdgeB = 0.5 *
                                  (met(i, j, k + 1, 2, 3) + met(i, j - 1, k + 1, 2, 3))
                    met13DEdgeR = 0.5 *
                                  (met(i, j, k - 1, 1, 3) + met(i + 1, j, k - 1, 1, 3))
                    met13DEdgeL = 0.5 *
                                  (met(i, j, k - 1, 1, 3) + met(i - 1, j, k - 1, 1, 3))
                    met23DEdgeF = 0.5 *
                                  (met(i, j, k - 1, 2, 3) + met(i, j + 1, k - 1, 2, 3))
                    met23DEdgeB = 0.5 *
                                  (met(i, j, k - 1, 2, 3) + met(i, j - 1, k - 1, 2, 3))

                    # Compute Rayleigh damping terms.
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
                    if (spongeLayer)
                        if (sponge_uv)
                            facEdgeR = facEdgeR +
                                       dt * 0.5 *
                                       (kr_sp_tfc[i, j, k] + kr_sp_tfc[i+1, j, k])
                            facEdgeL = facEdgeL +
                                       dt * 0.5 *
                                       (kr_sp_tfc[i, j, k] + kr_sp_tfc(i - 1, j, k))
                            facEdgeF = facEdgeF +
                                       dt * 0.5 *
                                       (kr_sp_tfc[i, j, k] + kr_sp_tfc(i, j + 1, k))
                            facEdgeB = facEdgeB +
                                       dt * 0.5 *
                                       (kr_sp_tfc[i, j, k] + kr_sp_tfc(i, j - 1, k))
                            facUEdgeR = facUEdgeR +
                                        dt *
                                        0.5 *
                                        (kr_sp_tfc(i, j, k + 1) +
                                         kr_sp_tfc(i + 1, j, k + 1))
                            facUEdgeL = facUEdgeL +
                                        dt *
                                        0.5 *
                                        (kr_sp_tfc(i, j, k + 1) +
                                         kr_sp_tfc(i - 1, j, k + 1))
                            facUEdgeF = facUEdgeF +
                                        dt *
                                        0.5 *
                                        (kr_sp_tfc(i, j, k + 1) +
                                         kr_sp_tfc(i, j + 1, k + 1))
                            facUEdgeB = facUEdgeB +
                                        dt *
                                        0.5 *
                                        (kr_sp_tfc(i, j, k + 1) +
                                         kr_sp_tfc(i, j - 1, k + 1))
                            facDEdgeR = facDEdgeR +
                                        dt *
                                        0.5 *
                                        (kr_sp_tfc(i, j, k - 1) +
                                         kr_sp_tfc(i + 1, j, k - 1))
                            facDEdgeL = facDEdgeL +
                                        dt *
                                        0.5 *
                                        (kr_sp_tfc(i, j, k - 1) +
                                         kr_sp_tfc(i - 1, j, k - 1))
                            facDEdgeF = facDEdgeF +
                                        dt *
                                        0.5 *
                                        (kr_sp_tfc(i, j, k - 1) +
                                         kr_sp_tfc(i, j + 1, k - 1))
                            facDEdgeB = facDEdgeB +
                                        dt *
                                        0.5 *
                                        (kr_sp_tfc(i, j, k - 1) +
                                         kr_sp_tfc(i, j - 1, k - 1))
                        end
                        facEdgeU = facEdgeU +
                                   dt * (jac(i, j, k + 1) * kr_sp_w_tfc[i, j, k] +
                                         jac[i, j, k] * kr_sp_w_tfc(i, j, k + 1)) /
                                   (jac[i, j, k] + jac(i, j, k + 1))
                        facEdgeD = facEdgeD +
                                   dt * (jac(i, j, k - 1) * kr_sp_w_tfc[i, j, k] +
                                         jac[i, j, k] * kr_sp_w_tfc(i, j, k - 1)) /
                                   (jac[i, j, k] + jac(i, j, k - 1))
                    end

                    # Compute implicit coefficients.
                    impHorEdgeR = 1.0 / (facEdgeR^2.0)
                    impHorEdgeL = 1.0 / (facEdgeL^2.0)
                    impHorEdgeF = 1.0 / (facEdgeF^2.0)
                    impHorEdgeB = 1.0 / (facEdgeB^2.0)
                    impHorUEdgeR = 1.0 / (facUEdgeR^2.0)
                    impHorUEdgeL = 1.0 / (facUEdgeL^2.0)
                    impHorUEdgeF = 1.0 / (facUEdgeF^2.0)
                    impHorUEdgeB = 1.0 / (facUEdgeB^2.0)
                    impHorDEdgeR = 1.0 / (facDEdgeR^2.0)
                    impHorDEdgeL = 1.0 / (facDEdgeL^2.0)
                    impHorDEdgeF = 1.0 / (facDEdgeF^2.0)
                    impHorDEdgeB = 1.0 / (facDEdgeB^2.0)
                    impVerEdgeU = 1.0 / (facEdgeU +
                                         rhostrattfcEdgeU / rhoEdgeU * bvsStratEdgeU * dt^2.0)
                    impVerEdgeD = 1.0 / (facEdgeD +
                                         rhostrattfcEdgeD / rhoEdgeD * bvsStratEdgeD * dt^2.0)

                    # Compute gradient coefficients

                    # G(i + 1 / 2)
                    if (k == 1 && zBoundary == "solid_wall")
                        gEdgeR = jacInv / dx * pEdgeRDiv * impHorEdgeR * facEdgeR /
                                 rhoEdgeR +
                                 jacInv / dz * pEdgeUDiv * impVerEdgeU * rhostrattfcEdgeU /
                                 rhoEdgeU *
                                 bvsStratEdgeU *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 1, 3) *
                                 jac(i, j, k + 1) / (jac[i, j, k] + jac(i, j, k + 1)) *
                                 impHorEdgeR *
                                 facEdgeR / rhoEdgeR
                    elseif (k == nz && zBoundary == "solid_wall")
                        gEdgeR = jacInv / dx * pEdgeRDiv * impHorEdgeR * facEdgeR /
                                 rhoEdgeR -
                                 jacInv / dz * pEdgeDDiv * impVerEdgeD * rhostrattfcEdgeD /
                                 rhoEdgeD *
                                 bvsStratEdgeD *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 1, 3) *
                                 jac(i, j, k - 1) / (jac[i, j, k] + jac(i, j, k - 1)) *
                                 impHorEdgeR *
                                 facEdgeR / rhoEdgeR
                    else
                        gEdgeR = jacInv / dx * pEdgeRDiv * impHorEdgeR * facEdgeR /
                                 rhoEdgeR +
                                 jacInv / dz * pEdgeUDiv * impVerEdgeU * rhostrattfcEdgeU /
                                 rhoEdgeU *
                                 bvsStratEdgeU *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 1, 3) *
                                 jac(i, j, k + 1) / (jac[i, j, k] + jac(i, j, k + 1)) *
                                 impHorEdgeR *
                                 facEdgeR / rhoEdgeR -
                                 jacInv / dz * pEdgeDDiv * impVerEdgeD * rhostrattfcEdgeD /
                                 rhoEdgeD *
                                 bvsStratEdgeD *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 1, 3) *
                                 jac(i, j, k - 1) / (jac[i, j, k] + jac(i, j, k - 1)) *
                                 impHorEdgeR *
                                 facEdgeR / rhoEdgeR
                    end

                    # G(i - 1 / 2)
                    if (k == 1 && zBoundary == "solid_wall")
                        gEdgeL = -jacInv / dx * pEdgeLDiv * impHorEdgeL * facEdgeL /
                                 rhoEdgeL +
                                 jacInv / dz * pEdgeUDiv * impVerEdgeU * rhostrattfcEdgeU /
                                 rhoEdgeU *
                                 bvsStratEdgeU *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 1, 3) *
                                 jac(i, j, k + 1) / (jac[i, j, k] + jac(i, j, k + 1)) *
                                 impHorEdgeL *
                                 facEdgeL / rhoEdgeL
                    elseif (k == nz && zBoundary == "solid_wall")
                        gEdgeL = -jacInv / dx * pEdgeLDiv * impHorEdgeL * facEdgeL /
                                 rhoEdgeL -
                                 jacInv / dz * pEdgeDDiv * impVerEdgeD * rhostrattfcEdgeD /
                                 rhoEdgeD *
                                 bvsStratEdgeD *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 1, 3) *
                                 jac(i, j, k - 1) / (jac[i, j, k] + jac(i, j, k - 1)) *
                                 impHorEdgeL *
                                 facEdgeL / rhoEdgeL
                    else
                        gEdgeL = -jacInv / dx * pEdgeLDiv * impHorEdgeL * facEdgeL /
                                 rhoEdgeL +
                                 jacInv / dz * pEdgeUDiv * impVerEdgeU * rhostrattfcEdgeU /
                                 rhoEdgeU *
                                 bvsStratEdgeU *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 1, 3) *
                                 jac(i, j, k + 1) / (jac[i, j, k] + jac(i, j, k + 1)) *
                                 impHorEdgeL *
                                 facEdgeL / rhoEdgeL -
                                 jacInv / dz * pEdgeDDiv * impVerEdgeD * rhostrattfcEdgeD /
                                 rhoEdgeD *
                                 bvsStratEdgeD *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 1, 3) *
                                 jac(i, j, k - 1) / (jac[i, j, k] + jac(i, j, k - 1)) *
                                 impHorEdgeL *
                                 facEdgeL / rhoEdgeL
                    end

                    # G(j + 1 / 2)
                    if (k == 1 && zBoundary == "solid_wall")
                        gEdgeF = jacInv / dy * pEdgeFDiv * impHorEdgeF * facEdgeF /
                                 rhoEdgeF +
                                 jacInv / dz * pEdgeUDiv * impVerEdgeU * rhostrattfcEdgeU /
                                 rhoEdgeU *
                                 bvsStratEdgeU *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 2, 3) *
                                 jac(i, j, k + 1) / (jac[i, j, k] + jac(i, j, k + 1)) *
                                 impHorEdgeF *
                                 facEdgeF / rhoEdgeF
                    elseif (k == nz && zBoundary == "solid_wall")
                        gEdgeF = jacInv / dy * pEdgeFDiv * impHorEdgeF * facEdgeF /
                                 rhoEdgeF -
                                 jacInv / dz * pEdgeDDiv * impVerEdgeD * rhostrattfcEdgeD /
                                 rhoEdgeD *
                                 bvsStratEdgeD *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 2, 3) *
                                 jac(i, j, k - 1) / (jac[i, j, k] + jac(i, j, k - 1)) *
                                 impHorEdgeF *
                                 facEdgeF / rhoEdgeF
                    else
                        gEdgeF = jacInv / dy * pEdgeFDiv * impHorEdgeF * facEdgeF /
                                 rhoEdgeF +
                                 jacInv / dz * pEdgeUDiv * impVerEdgeU * rhostrattfcEdgeU /
                                 rhoEdgeU *
                                 bvsStratEdgeU *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 2, 3) *
                                 jac(i, j, k + 1) / (jac[i, j, k] + jac(i, j, k + 1)) *
                                 impHorEdgeF *
                                 facEdgeF / rhoEdgeF -
                                 jacInv / dz * pEdgeDDiv * impVerEdgeD * rhostrattfcEdgeD /
                                 rhoEdgeD *
                                 bvsStratEdgeD *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 2, 3) *
                                 jac(i, j, k - 1) / (jac[i, j, k] + jac(i, j, k - 1)) *
                                 impHorEdgeF *
                                 facEdgeF / rhoEdgeF
                    end

                    # G(j - 1 / 2)
                    if (k == 1 && zBoundary == "solid_wall")
                        gEdgeB = -jacInv / dy * pEdgeBDiv * impHorEdgeB * facEdgeB /
                                 rhoEdgeB +
                                 jacInv / dz * pEdgeUDiv * impVerEdgeU * rhostrattfcEdgeU /
                                 rhoEdgeU *
                                 bvsStratEdgeU *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 2, 3) *
                                 jac(i, j, k + 1) / (jac[i, j, k] + jac(i, j, k + 1)) *
                                 impHorEdgeB *
                                 facEdgeB / rhoEdgeB
                    elseif (k == nz && zBoundary == "solid_wall")
                        gEdgeB = -jacInv / dy * pEdgeBDiv * impHorEdgeB * facEdgeB /
                                 rhoEdgeB -
                                 jacInv / dz * pEdgeDDiv * impVerEdgeD * rhostrattfcEdgeD /
                                 rhoEdgeD *
                                 bvsStratEdgeD *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 2, 3) *
                                 jac(i, j, k - 1) / (jac[i, j, k] + jac(i, j, k - 1)) *
                                 impHorEdgeB *
                                 facEdgeB / rhoEdgeB
                    else
                        gEdgeB = -jacInv / dy * pEdgeBDiv * impHorEdgeB * facEdgeB /
                                 rhoEdgeB +
                                 jacInv / dz * pEdgeUDiv * impVerEdgeU * rhostrattfcEdgeU /
                                 rhoEdgeU *
                                 bvsStratEdgeU *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 2, 3) *
                                 jac(i, j, k + 1) / (jac[i, j, k] + jac(i, j, k + 1)) *
                                 impHorEdgeB *
                                 facEdgeB / rhoEdgeB -
                                 jacInv / dz * pEdgeDDiv * impVerEdgeD * rhostrattfcEdgeD /
                                 rhoEdgeD *
                                 bvsStratEdgeD *
                                 dt^2.0 *
                                 0.5 *
                                 met(i, j, k, 2, 3) *
                                 jac(i, j, k - 1) / (jac[i, j, k] + jac(i, j, k - 1)) *
                                 impHorEdgeB *
                                 facEdgeB / rhoEdgeB
                    end

                    # G(k + 1 / 2)
                    if (k == nz && zBoundary == "solid_wall")
                        gEdgeU = 0.0
                    else
                        gEdgeU = jacInv / dz * pEdgeUDiv * impVerEdgeU / rhoEdgeU
                    end

                    # G(k - 1 / 2)
                    if (k == 1 && zBoundary == "solid_wall")
                        gEdgeD = 0.0
                    else
                        gEdgeD = -jacInv / dz * pEdgeDDiv * impVerEdgeD / rhoEdgeD
                    end

                    # G(i + 1 / 2, k + 1)
                    if (k == nz && zBoundary == "solid_wall")
                        gUEdgeR = 0.0
                    else
                        gUEdgeR = jacInv / dz * pEdgeUDiv * impVerEdgeU * rhostrattfcEdgeU /
                                  rhoEdgeU *
                                  bvsStratEdgeU *
                                  dt^2.0 *
                                  0.5 *
                                  met(i, j, k + 1, 1, 3) *
                                  jac[i, j, k] / (jac[i, j, k] + jac(i, j, k + 1)) *
                                  impHorUEdgeR *
                                  facUEdgeR / rhoUEdgeR
                    end

                    # G(i - 1 / 2, k + 1)
                    if (k == nz && zBoundary == "solid_wall")
                        gUEdgeL = 0.0
                    else
                        gUEdgeL = jacInv / dz * pEdgeUDiv * impVerEdgeU * rhostrattfcEdgeU /
                                  rhoEdgeU *
                                  bvsStratEdgeU *
                                  dt^2.0 *
                                  0.5 *
                                  met(i, j, k + 1, 1, 3) *
                                  jac[i, j, k] / (jac[i, j, k] + jac(i, j, k + 1)) *
                                  impHorUEdgeL *
                                  facUEdgeL / rhoUEdgeL
                    end

                    # G(j + 1 / 2, k + 1)
                    if (k == nz && zBoundary == "solid_wall")
                        gUEdgeF = 0.0
                    else
                        gUEdgeF = jacInv / dz * pEdgeUDiv * impVerEdgeU * rhostrattfcEdgeU /
                                  rhoEdgeU *
                                  bvsStratEdgeU *
                                  dt^2.0 *
                                  0.5 *
                                  met(i, j, k + 1, 2, 3) *
                                  jac[i, j, k] / (jac[i, j, k] + jac(i, j, k + 1)) *
                                  impHorUEdgeF *
                                  facUEdgeF / rhoUEdgeF
                    end

                    # G(j - 1 / 2, k + 1)
                    if (k == nz && zBoundary == "solid_wall")
                        gUEdgeB = 0.0
                    else
                        gUEdgeB = jacInv / dz * pEdgeUDiv * impVerEdgeU * rhostrattfcEdgeU /
                                  rhoEdgeU *
                                  bvsStratEdgeU *
                                  dt^2.0 *
                                  0.5 *
                                  met(i, j, k + 1, 2, 3) *
                                  jac[i, j, k] / (jac[i, j, k] + jac(i, j, k + 1)) *
                                  impHorUEdgeB *
                                  facUEdgeB / rhoUEdgeB
                    end

                    # G(i + 1 / 2, k - 1)
                    if (k == 1 && zBoundary == "solid_wall")
                        gDEdgeR = 0.0
                    else
                        gDEdgeR = -jacInv / dz * pEdgeDDiv * impVerEdgeD * rhostrattfcEdgeD /
                                  rhoEdgeD *
                                  bvsStratEdgeD *
                                  dt^2.0 *
                                  0.5 *
                                  met(i, j, k - 1, 1, 3) *
                                  jac[i, j, k] / (jac[i, j, k] + jac(i, j, k - 1)) *
                                  impHorDEdgeR *
                                  facDEdgeR / rhoDEdgeR
                    end

                    # G(i - 1 / 2, k - 1)
                    if (k == 1 && zBoundary == "solid_wall")
                        gDEdgeL = 0.0
                    else
                        gDEdgeL = -jacInv / dz * pEdgeDDiv * impVerEdgeD * rhostrattfcEdgeD /
                                  rhoEdgeD *
                                  bvsStratEdgeD *
                                  dt^2.0 *
                                  0.5 *
                                  met(i, j, k - 1, 1, 3) *
                                  jac[i, j, k] / (jac[i, j, k] + jac(i, j, k - 1)) *
                                  impHorDEdgeL *
                                  facDEdgeL / rhoDEdgeL
                    end

                    # G(j + 1 / 2, k - 1)
                    if (k == 1 && zBoundary == "solid_wall")
                        gDEdgeF = 0.0
                    else
                        gDEdgeF = -jacInv / dz * pEdgeDDiv * impVerEdgeD * rhostrattfcEdgeD /
                                  rhoEdgeD *
                                  bvsStratEdgeD *
                                  dt^2.0 *
                                  0.5 *
                                  met(i, j, k - 1, 2, 3) *
                                  jac[i, j, k] / (jac[i, j, k] + jac(i, j, k - 1)) *
                                  impHorDEdgeF *
                                  facDEdgeF / rhoDEdgeF
                    end

                    # G(j - 1 / 2, k - 1)
                    if (k == 1 && zBoundary == "solid_wall")
                        gDEdgeB = 0.0
                    else
                        gDEdgeB = -jacInv / dz * pEdgeDDiv * impVerEdgeD * rhostrattfcEdgeD /
                                  rhoEdgeD *
                                  bvsStratEdgeD *
                                  dt^2.0 *
                                  0.5 *
                                  met(i, j, k - 1, 2, 3) *
                                  jac[i, j, k] / (jac[i, j, k] + jac(i, j, k - 1)) *
                                  impHorDEdgeB *
                                  facDEdgeB / rhoDEdgeB
                    end

                    # Compute tensor elements

                    # ------------------- A(i,j,k) --------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AC = -gEdgeR * pEdgeRGra * (1.0 / dx + 0.75 * met13EdgeR / dz) +
                             gEdgeL * pEdgeLGra * (1.0 / dx - 0.75 * met13EdgeL / dz) -
                             gEdgeF * pEdgeFGra * (1.0 / dy + 0.75 * met23EdgeF / dz) +
                             gEdgeB * pEdgeBGra * (1.0 / dy - 0.75 * met23EdgeB / dz) -
                             gEdgeU * pEdgeUGra * met33EdgeU / dz -
                             gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz -
                             gUEdgeL * pUEdgeLGra * 0.25 * met13UEdgeL / dz -
                             gUEdgeF * pUEdgeFGra * 0.25 * met23UEdgeF / dz -
                             gUEdgeB * pUEdgeBGra * 0.25 * met23UEdgeB / dz
                    elseif (k == 2 && zBoundary == "solid_wall")
                        AC = -gEdgeR * pEdgeRGra / dx + gEdgeL * pEdgeLGra / dx -
                             gEdgeF * pEdgeFGra / dy + gEdgeB * pEdgeBGra / dy -
                             gEdgeU * pEdgeUGra * met33EdgeU / dz +
                             gEdgeD * pEdgeDGra * met33EdgeD / dz -
                             gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz -
                             gUEdgeL * pUEdgeLGra * 0.25 * met13UEdgeL / dz -
                             gUEdgeF * pUEdgeFGra * 0.25 * met23UEdgeF / dz -
                             gUEdgeB * pUEdgeBGra * 0.25 * met23UEdgeB / dz +
                             gDEdgeR * pDEdgeRGra * met13DEdgeR / dz +
                             gDEdgeL * pDEdgeLGra * met13DEdgeL / dz +
                             gDEdgeF * pDEdgeFGra * met23DEdgeF / dz +
                             gDEdgeB * pDEdgeBGra * met23DEdgeB / dz
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        AC = -gEdgeR * pEdgeRGra / dx + gEdgeL * pEdgeLGra / dx -
                             gEdgeF * pEdgeFGra / dy + gEdgeB * pEdgeBGra / dy -
                             gEdgeU * pEdgeUGra * met33EdgeU / dz +
                             gEdgeD * pEdgeDGra * met33EdgeD / dz -
                             gUEdgeR * pUEdgeRGra * met13UEdgeR / dz -
                             gUEdgeL * pUEdgeLGra * met13UEdgeL / dz -
                             gUEdgeF * pUEdgeFGra * met23UEdgeF / dz -
                             gUEdgeB * pUEdgeBGra * met23UEdgeB / dz +
                             gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz +
                             gDEdgeL * pDEdgeLGra * 0.25 * met13DEdgeL / dz +
                             gDEdgeF * pDEdgeFGra * 0.25 * met23DEdgeF / dz +
                             gDEdgeB * pDEdgeBGra * 0.25 * met23DEdgeB / dz
                    elseif (k == nz && zBoundary == "solid_wall")
                        AC = -gEdgeR * pEdgeRGra * (1.0 / dx - 0.75 * met13EdgeR / dz) +
                             gEdgeL * pEdgeLGra * (1.0 / dx + 0.75 * met13EdgeL / dz) -
                             gEdgeF * pEdgeFGra * (1.0 / dy - 0.75 * met23EdgeF / dz) +
                             gEdgeB * pEdgeBGra * (1.0 / dy + 0.75 * met23EdgeB / dz) +
                             gEdgeD * pEdgeDGra * met33EdgeD / dz +
                             gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz +
                             gDEdgeL * pDEdgeLGra * 0.25 * met13DEdgeL / dz +
                             gDEdgeF * pDEdgeFGra * 0.25 * met23DEdgeF / dz +
                             gDEdgeB * pDEdgeBGra * 0.25 * met23DEdgeB / dz
                    else
                        AC = -gEdgeR * pEdgeRGra / dx + gEdgeL * pEdgeLGra / dx -
                             gEdgeF * pEdgeFGra / dy + gEdgeB * pEdgeBGra / dy -
                             gEdgeU * pEdgeUGra * met33EdgeU / dz +
                             gEdgeD * pEdgeDGra * met33EdgeD / dz -
                             gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz -
                             gUEdgeL * pUEdgeLGra * 0.25 * met13UEdgeL / dz -
                             gUEdgeF * pUEdgeFGra * 0.25 * met23UEdgeF / dz -
                             gUEdgeB * pUEdgeBGra * 0.25 * met23UEdgeB / dz +
                             gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz +
                             gDEdgeL * pDEdgeLGra * 0.25 * met13DEdgeL / dz +
                             gDEdgeF * pDEdgeFGra * 0.25 * met23DEdgeF / dz +
                             gDEdgeB * pDEdgeBGra * 0.25 * met23DEdgeB / dz
                    end

                    # ------------------ A(i+1,j,k) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AR = gEdgeR * pEdgeRGra * (1.0 / dx - 0.75 * met13EdgeR / dz) +
                             gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                             jac(i + 1, j, k + 1) /
                             (jac[i+1, j, k] + jac(i + 1, j, k + 1)) -
                             gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz
                    elseif (k == 2 && zBoundary == "solid_wall")
                        AR = gEdgeR * pEdgeRGra / dx +
                             gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                             jac(i + 1, j, k + 1) /
                             (jac[i+1, j, k] + jac(i + 1, j, k + 1)) +
                             gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                             jac(i + 1, j, k - 1) /
                             (jac[i+1, j, k] + jac(i + 1, j, k - 1)) -
                             gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz +
                             gDEdgeR * pDEdgeRGra * met13DEdgeR / dz
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        AR = gEdgeR * pEdgeRGra / dx +
                             gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                             jac(i + 1, j, k + 1) /
                             (jac[i+1, j, k] + jac(i + 1, j, k + 1)) +
                             gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                             jac(i + 1, j, k - 1) /
                             (jac[i+1, j, k] + jac(i + 1, j, k - 1)) -
                             gUEdgeR * pUEdgeRGra * met13UEdgeR / dz +
                             gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz
                    elseif (k == nz && zBoundary == "solid_wall")
                        AR = gEdgeR * pEdgeRGra * (1.0 / dx + 0.75 * met13EdgeR / dz) +
                             gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                             jac(i + 1, j, k - 1) /
                             (jac[i+1, j, k] + jac(i + 1, j, k - 1)) +
                             gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz
                    else
                        AR = gEdgeR * pEdgeRGra / dx +
                             gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                             jac(i + 1, j, k + 1) /
                             (jac[i+1, j, k] + jac(i + 1, j, k + 1)) +
                             gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                             jac(i + 1, j, k - 1) /
                             (jac[i+1, j, k] + jac(i + 1, j, k - 1)) -
                             gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz +
                             gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz
                    end

                    # ------------------ A(i-1,j,k) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AL = -gEdgeL * pEdgeLGra * (1.0 / dx + 0.75 * met13EdgeL / dz) -
                             gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                             jac(i - 1, j, k + 1) /
                             (jac(i - 1, j, k) + jac(i - 1, j, k + 1)) -
                             gUEdgeL * pUEdgeLGra * 0.25 * met13UEdgeL / dz
                    elseif (k == 2 && zBoundary == "solid_wall")
                        AL = -gEdgeL * pEdgeLGra / dx -
                             gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                             jac(i - 1, j, k + 1) /
                             (jac(i - 1, j, k) + jac(i - 1, j, k + 1)) -
                             gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                             jac(i - 1, j, k - 1) /
                             (jac(i - 1, j, k) + jac(i - 1, j, k - 1)) -
                             gUEdgeL * pUEdgeLGra * 0.25 * met13UEdgeL / dz +
                             gDEdgeL * pDEdgeLGra * met13DEdgeL / dz
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        AL = -gEdgeL * pEdgeLGra / dx -
                             gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                             jac(i - 1, j, k + 1) /
                             (jac(i - 1, j, k) + jac(i - 1, j, k + 1)) -
                             gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                             jac(i - 1, j, k - 1) /
                             (jac(i - 1, j, k) + jac(i - 1, j, k - 1)) -
                             gUEdgeL * pUEdgeLGra * met13UEdgeL / dz +
                             gDEdgeL * pDEdgeLGra * 0.25 * met13DEdgeL / dz
                    elseif (k == nz && zBoundary == "solid_wall")
                        AL = -gEdgeL * pEdgeLGra * (1.0 / dx - 0.75 * met13EdgeL / dz) -
                             gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                             jac(i - 1, j, k - 1) /
                             (jac(i - 1, j, k) + jac(i - 1, j, k - 1)) +
                             gDEdgeL * pDEdgeLGra * 0.25 * met13DEdgeL / dz
                    else
                        AL = -gEdgeL * pEdgeLGra / dx -
                             gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                             jac(i - 1, j, k + 1) /
                             (jac(i - 1, j, k) + jac(i - 1, j, k + 1)) -
                             gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                             jac(i - 1, j, k - 1) /
                             (jac(i - 1, j, k) + jac(i - 1, j, k - 1)) -
                             gUEdgeL * pUEdgeLGra * 0.25 * met13UEdgeL / dz +
                             gDEdgeL * pDEdgeLGra * 0.25 * met13DEdgeL / dz
                    end

                    # ------------------ A(i,j+1,k) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AF = gEdgeF * pEdgeFGra * (1.0 / dy - 0.75 * met23EdgeF / dz) +
                             gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                             jac(i, j + 1, k + 1) /
                             (jac(i, j + 1, k) + jac(i, j + 1, k + 1)) -
                             gUEdgeF * pUEdgeFGra * 0.25 * met23UEdgeF / dz
                    elseif (k == 2 && zBoundary == "solid_wall")
                        AF = gEdgeF * pEdgeFGra / dy +
                             gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                             jac(i, j + 1, k + 1) /
                             (jac(i, j + 1, k) + jac(i, j + 1, k + 1)) +
                             gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                             jac(i, j + 1, k - 1) /
                             (jac(i, j + 1, k) + jac(i, j + 1, k - 1)) -
                             gUEdgeF * pUEdgeFGra * 0.25 * met23UEdgeF / dz +
                             gDEdgeF * pDEdgeFGra * met23DEdgeF / dz
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        AF = gEdgeF * pEdgeFGra / dy +
                             gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                             jac(i, j + 1, k + 1) /
                             (jac(i, j + 1, k) + jac(i, j + 1, k + 1)) +
                             gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                             jac(i, j + 1, k - 1) /
                             (jac(i, j + 1, k) + jac(i, j + 1, k - 1)) -
                             gUEdgeF * pUEdgeFGra * met23UEdgeF / dz +
                             gDEdgeF * pDEdgeFGra * 0.25 * met23DEdgeF / dz
                    elseif (k == nz && zBoundary == "solid_wall")
                        AF = gEdgeF * pEdgeFGra * (1.0 / dy + 0.75 * met23EdgeF / dz) +
                             gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                             jac(i, j + 1, k - 1) /
                             (jac(i, j + 1, k) + jac(i, j + 1, k - 1)) +
                             gDEdgeF * pDEdgeFGra * 0.25 * met23DEdgeF / dz
                    else
                        AF = gEdgeF * pEdgeFGra / dy +
                             gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                             jac(i, j + 1, k + 1) /
                             (jac(i, j + 1, k) + jac(i, j + 1, k + 1)) +
                             gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                             jac(i, j + 1, k - 1) /
                             (jac(i, j + 1, k) + jac(i, j + 1, k - 1)) -
                             gUEdgeF * pUEdgeFGra * 0.25 * met23UEdgeF / dz +
                             gDEdgeF * pDEdgeFGra * 0.25 * met23DEdgeF / dz
                    end

                    # ------------------ A(i,j-1,k) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AB = -gEdgeB * pEdgeBGra * (1.0 / dy + 0.75 * met23EdgeB / dz) -
                             gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                             jac(i, j - 1, k + 1) /
                             (jac(i, j - 1, k) + jac(i, j - 1, k + 1)) -
                             gUEdgeB * pUEdgeBGra * 0.25 * met23UEdgeB / dz
                    elseif (k == 2 && zBoundary == "solid_wall")
                        AB = -gEdgeB * pEdgeBGra / dy -
                             gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                             jac(i, j - 1, k + 1) /
                             (jac(i, j - 1, k) + jac(i, j - 1, k + 1)) -
                             gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                             jac(i, j - 1, k - 1) /
                             (jac(i, j - 1, k) + jac(i, j - 1, k - 1)) -
                             gUEdgeB * pUEdgeBGra * 0.25 * met23UEdgeB / dz +
                             gDEdgeB * pDEdgeBGra * met23DEdgeB / dz
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        AB = -gEdgeB * pEdgeBGra / dy -
                             gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                             jac(i, j - 1, k + 1) /
                             (jac(i, j - 1, k) + jac(i, j - 1, k + 1)) -
                             gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                             jac(i, j - 1, k - 1) /
                             (jac(i, j - 1, k) + jac(i, j - 1, k - 1)) -
                             gUEdgeB * pUEdgeBGra * met23UEdgeB / dz +
                             gDEdgeB * pDEdgeBGra * 0.25 * met23DEdgeB / dz
                    elseif (k == nz && zBoundary == "solid_wall")
                        AB = -gEdgeB * pEdgeBGra * (1.0 / dy - 0.75 * met23EdgeB / dz) -
                             gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                             jac(i, j - 1, k - 1) /
                             (jac(i, j - 1, k) + jac(i, j - 1, k - 1)) +
                             gDEdgeB * pDEdgeBGra * 0.25 * met23DEdgeB / dz
                    else
                        AB = -gEdgeB * pEdgeBGra / dy -
                             gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                             jac(i, j - 1, k + 1) /
                             (jac(i, j - 1, k) + jac(i, j - 1, k + 1)) -
                             gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                             jac(i, j - 1, k - 1) /
                             (jac(i, j - 1, k) + jac(i, j - 1, k - 1)) -
                             gUEdgeB * pUEdgeBGra * 0.25 * met23UEdgeB / dz +
                             gDEdgeB * pDEdgeBGra * 0.25 * met23DEdgeB / dz
                    end

                    # ------------------ A(i,j,k+1) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AU = gEdgeR * pEdgeRGra * met13EdgeR / dz +
                             gEdgeL * pEdgeLGra * met13EdgeL / dz +
                             gEdgeF * pEdgeFGra * met23EdgeF / dz +
                             gEdgeB * pEdgeBGra * met23EdgeB / dz +
                             gEdgeU * pEdgeUGra * met33EdgeU / dz -
                             gUEdgeR * pUEdgeRGra / dx + gUEdgeL * pUEdgeLGra / dx -
                             gUEdgeF * pUEdgeFGra / dy + gUEdgeB * pUEdgeBGra / dy
                    elseif (k == 2 && zBoundary == "solid_wall")
                        AU = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz +
                             gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz +
                             gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz +
                             gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz +
                             gEdgeU * pEdgeUGra * met33EdgeU / dz -
                             gUEdgeR * pUEdgeRGra / dx + gUEdgeL * pUEdgeLGra / dx -
                             gUEdgeF * pUEdgeFGra / dy + gUEdgeB * pUEdgeBGra / dy -
                             gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz -
                             gDEdgeL * pDEdgeLGra * 0.25 * met13DEdgeL / dz -
                             gDEdgeF * pDEdgeFGra * 0.25 * met23DEdgeF / dz -
                             gDEdgeB * pDEdgeBGra * 0.25 * met23DEdgeB / dz
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        AU = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz +
                             gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz +
                             gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz +
                             gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz +
                             gEdgeU * pEdgeUGra * met33EdgeU / dz -
                             gUEdgeR * pUEdgeRGra * (1.0 / dx - 0.75 * met13UEdgeR / dz) +
                             gUEdgeL * pUEdgeLGra * (1.0 / dx + 0.75 * met13UEdgeL / dz) -
                             gUEdgeF * pUEdgeFGra * (1.0 / dy - 0.75 * met23UEdgeF / dz) +
                             gUEdgeB * pUEdgeBGra * (1.0 / dy + 0.75 * met23UEdgeB / dz)
                    elseif (k == nz && zBoundary == "solid_wall")
                        AU = 0.0
                    else
                        AU = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz +
                             gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz +
                             gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz +
                             gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz +
                             gEdgeU * pEdgeUGra * met33EdgeU / dz -
                             gUEdgeR * pUEdgeRGra / dx + gUEdgeL * pUEdgeLGra / dx -
                             gUEdgeF * pUEdgeFGra / dy + gUEdgeB * pUEdgeBGra / dy
                    end

                    # ------------------ A(i,j,k-1) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AD = 0.0
                    elseif (k == 2 && zBoundary == "solid_wall")
                        AD = -gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz -
                             gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz -
                             gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz -
                             gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz -
                             gEdgeD * pEdgeDGra * met33EdgeD / dz -
                             gDEdgeR * pDEdgeRGra * (1.0 / dx + 0.75 * met13DEdgeR / dz) +
                             gDEdgeL * pDEdgeLGra * (1.0 / dx - 0.75 * met13DEdgeL / dz) -
                             gDEdgeF * pDEdgeFGra * (1.0 / dy + 0.75 * met23DEdgeF / dz) +
                             gDEdgeB * pDEdgeBGra * (1.0 / dy - 0.75 * met23DEdgeB / dz)
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        AD = -gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz -
                             gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz -
                             gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz -
                             gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz -
                             gEdgeD * pEdgeDGra * met33EdgeD / dz -
                             gDEdgeR * pDEdgeRGra / dx + gDEdgeL * pDEdgeLGra / dx -
                             gDEdgeF * pDEdgeFGra / dy +
                             gDEdgeB * pDEdgeBGra / dy +
                             gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz +
                             gUEdgeL * pUEdgeLGra * 0.25 * met13UEdgeL / dz +
                             gUEdgeF * pUEdgeFGra * 0.25 * met23UEdgeF / dz +
                             gUEdgeB * pUEdgeBGra * 0.25 * met23UEdgeB / dz
                    elseif (k == nz && zBoundary == "solid_wall")
                        AD = -gEdgeR * pEdgeRGra * met13EdgeR / dz -
                             gEdgeL * pEdgeLGra * met13EdgeL / dz -
                             gEdgeF * pEdgeFGra * met23EdgeF / dz -
                             gEdgeB * pEdgeBGra * met23EdgeB / dz -
                             gEdgeD * pEdgeDGra * met33EdgeD / dz -
                             gDEdgeR * pDEdgeRGra / dx + gDEdgeL * pDEdgeLGra / dx -
                             gDEdgeF * pDEdgeFGra / dy + gDEdgeB * pDEdgeBGra / dy
                    else
                        AD = -gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz -
                             gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz -
                             gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz -
                             gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz -
                             gEdgeD * pEdgeDGra * met33EdgeD / dz -
                             gDEdgeR * pDEdgeRGra / dx + gDEdgeL * pDEdgeLGra / dx -
                             gDEdgeF * pDEdgeFGra / dy + gDEdgeB * pDEdgeBGra / dy
                    end

                    # ----------------- A(i+1,j,k+1) ------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ARU = gEdgeR * pEdgeRGra * met13EdgeR / dz +
                              gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                              jac[i+1, j, k] /
                              (jac[i+1, j, k] + jac(i + 1, j, k + 1)) +
                              gUEdgeR * pUEdgeRGra / dx
                    elseif (k == 2 && zBoundary == "solid_wall")
                        ARU = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz +
                              gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                              jac[i+1, j, k] /
                              (jac[i+1, j, k] + jac(i + 1, j, k + 1)) +
                              gUEdgeR * pUEdgeRGra / dx -
                              gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        ARU = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz +
                              gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                              jac[i+1, j, k] /
                              (jac[i+1, j, k] + jac(i + 1, j, k + 1)) +
                              gUEdgeR * pUEdgeRGra * (1.0 / dx + 0.75 * met13UEdgeR / dz)
                    elseif (k == nz && zBoundary == "solid_wall")
                        ARU = 0.0
                    else
                        ARU = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz +
                              gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                              jac[i+1, j, k] /
                              (jac[i+1, j, k] + jac(i + 1, j, k + 1)) +
                              gUEdgeR * pUEdgeRGra / dx
                    end

                    # ----------------- A(i+1,j,k-1) ------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ARD = 0.0
                    elseif (k == 2 && zBoundary == "solid_wall")
                        ARD = -gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz +
                              gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                              jac[i+1, j, k] /
                              (jac[i+1, j, k] + jac(i + 1, j, k - 1)) +
                              gDEdgeR * pDEdgeRGra * (1.0 / dx - 0.75 * met13DEdgeR / dz)
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        ARD = -gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz +
                              gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                              jac[i+1, j, k] /
                              (jac[i+1, j, k] + jac(i + 1, j, k - 1)) +
                              gDEdgeR * pDEdgeRGra / dx +
                              gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz
                    elseif (k == nz && zBoundary == "solid_wall")
                        ARD = -gEdgeR * pEdgeRGra * met13EdgeR / dz +
                              gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                              jac[i+1, j, k] /
                              (jac[i+1, j, k] + jac(i + 1, j, k - 1)) +
                              gDEdgeR * pDEdgeRGra / dx
                    else
                        ARD = -gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz +
                              gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                              jac[i+1, j, k] /
                              (jac[i+1, j, k] + jac(i + 1, j, k - 1)) +
                              gDEdgeR * pDEdgeRGra / dx
                    end

                    # ----------------- A(i-1,j,k+1) ------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ALU = gEdgeL * pEdgeLGra * met13EdgeL / dz -
                              gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                              jac(i - 1, j, k) /
                              (jac(i - 1, j, k) + jac(i - 1, j, k + 1)) -
                              gUEdgeL * pUEdgeLGra / dx
                    elseif (k == 2 && zBoundary == "solid_wall")
                        ALU = gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz -
                              gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                              jac(i - 1, j, k) /
                              (jac(i - 1, j, k) + jac(i - 1, j, k + 1)) -
                              gUEdgeL * pUEdgeLGra / dx -
                              gDEdgeL * pDEdgeLGra * 0.25 * met13DEdgeL / dz
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        ALU = gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz -
                              gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                              jac(i - 1, j, k) /
                              (jac(i - 1, j, k) + jac(i - 1, j, k + 1)) -
                              gUEdgeL * pUEdgeLGra * (1.0 / dx - 0.75 * met13UEdgeL / dz)
                    elseif (k == nz && zBoundary == "solid_wall")
                        ALU = 0.0
                    else
                        ALU = gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz -
                              gEdgeU * pEdgeUGra * 0.5 * met13EdgeU / dx *
                              jac(i - 1, j, k) /
                              (jac(i - 1, j, k) + jac(i - 1, j, k + 1)) -
                              gUEdgeL * pUEdgeLGra / dx
                    end

                    # ----------------- A(i-1,j,k-1) ------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ALD = 0.0
                    elseif (k == 2 && zBoundary == "solid_wall")
                        ALD = -gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz -
                              gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                              jac(i - 1, j, k) /
                              (jac(i - 1, j, k) + jac(i - 1, j, k - 1)) -
                              gDEdgeL * pDEdgeLGra * (1.0 / dx + 0.75 * met13DEdgeL / dz)
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        ALD = -gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz -
                              gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                              jac(i - 1, j, k) /
                              (jac(i - 1, j, k) + jac(i - 1, j, k - 1)) -
                              gDEdgeL * pDEdgeLGra / dx +
                              gUEdgeL * pUEdgeLGra * 0.25 * met13UEdgeL / dz
                    elseif (k == nz && zBoundary == "solid_wall")
                        ALD = -gEdgeL * pEdgeLGra * met13EdgeL / dz -
                              gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                              jac(i - 1, j, k) /
                              (jac(i - 1, j, k) + jac(i - 1, j, k - 1)) -
                              gDEdgeL * pDEdgeLGra / dx
                    else
                        ALD = -gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz -
                              gEdgeD * pEdgeDGra * 0.5 * met13EdgeD / dx *
                              jac(i - 1, j, k) /
                              (jac(i - 1, j, k) + jac(i - 1, j, k - 1)) -
                              gDEdgeL * pDEdgeLGra / dx
                    end

                    # ----------------- A(i,j+1,k+1) ------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AFU = gEdgeF * pEdgeFGra * met23EdgeF / dz +
                              gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                              jac(i, j + 1, k) /
                              (jac(i, j + 1, k) + jac(i, j + 1, k + 1)) +
                              gUEdgeF * pUEdgeFGra / dy
                    elseif (k == 2 && zBoundary == "solid_wall")
                        AFU = gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz +
                              gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                              jac(i, j + 1, k) /
                              (jac(i, j + 1, k) + jac(i, j + 1, k + 1)) +
                              gUEdgeF * pUEdgeFGra / dy -
                              gDEdgeF * pDEdgeFGra * 0.25 * met23DEdgeF / dz
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        AFU = gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz +
                              gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                              jac(i, j + 1, k) /
                              (jac(i, j + 1, k) + jac(i, j + 1, k + 1)) +
                              gUEdgeF * pUEdgeFGra * (1.0 / dy + 0.75 * met23UEdgeF / dz)
                    elseif (k == nz && zBoundary == "solid_wall")
                        AFU = 0.0
                    else
                        AFU = gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz +
                              gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                              jac(i, j + 1, k) /
                              (jac(i, j + 1, k) + jac(i, j + 1, k + 1)) +
                              gUEdgeF * pUEdgeFGra / dy
                    end

                    # ----------------- A(i,j+1,k-1) ------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AFD = 0.0
                    elseif (k == 2 && zBoundary == "solid_wall")
                        AFD = -gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz +
                              gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                              jac(i, j + 1, k) /
                              (jac(i, j + 1, k) + jac(i, j + 1, k - 1)) +
                              gDEdgeF * pDEdgeFGra * (1.0 / dy - 0.75 * met23DEdgeF / dz)
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        AFD = -gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz +
                              gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                              jac(i, j + 1, k) /
                              (jac(i, j + 1, k) + jac(i, j + 1, k - 1)) +
                              gDEdgeF * pDEdgeFGra / dy +
                              gUEdgeF * pUEdgeFGra * 0.25 * met23UEdgeF / dz
                    elseif (k == nz && zBoundary == "solid_wall")
                        AFD = -gEdgeF * pEdgeFGra * met23EdgeF / dz +
                              gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                              jac(i, j + 1, k) /
                              (jac(i, j + 1, k) + jac(i, j + 1, k - 1)) +
                              gDEdgeF * pDEdgeFGra / dy
                    else
                        AFD = -gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz +
                              gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                              jac(i, j + 1, k) /
                              (jac(i, j + 1, k) + jac(i, j + 1, k - 1)) +
                              gDEdgeF * pDEdgeFGra / dy
                    end

                    # ----------------- A(i,j-1,k+1) ------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ABU = gEdgeB * pEdgeBGra * met23EdgeB / dz -
                              gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                              jac(i, j - 1, k) /
                              (jac(i, j - 1, k) + jac(i, j - 1, k + 1)) -
                              gUEdgeB * pUEdgeBGra / dy
                    elseif (k == 2 && zBoundary == "solid_wall")
                        ABU = gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz -
                              gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                              jac(i, j - 1, k) /
                              (jac(i, j - 1, k) + jac(i, j - 1, k + 1)) -
                              gUEdgeB * pUEdgeBGra / dy -
                              gDEdgeB * pDEdgeBGra * 0.25 * met23DEdgeB / dz
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        ABU = gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz -
                              gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                              jac(i, j - 1, k) /
                              (jac(i, j - 1, k) + jac(i, j - 1, k + 1)) -
                              gUEdgeB * pUEdgeBGra * (1.0 / dy - 0.75 * met23UEdgeB / dz)
                    elseif (k == nz && zBoundary == "solid_wall")
                        ABU = 0.0
                    else
                        ABU = gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz -
                              gEdgeU * pEdgeUGra * 0.5 * met23EdgeU / dy *
                              jac(i, j - 1, k) /
                              (jac(i, j - 1, k) + jac(i, j - 1, k + 1)) -
                              gUEdgeB * pUEdgeBGra / dy
                    end

                    # ----------------- A(i,j-1,k-1) ------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ABD = 0.0
                    elseif (k == 2 && zBoundary == "solid_wall")
                        ABD = -gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz -
                              gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                              jac(i, j - 1, k) /
                              (jac(i, j - 1, k) + jac(i, j - 1, k - 1)) -
                              gDEdgeB * pDEdgeBGra * (1.0 / dy + 0.75 * met23DEdgeB / dz)
                    elseif (k == nz - 1 && zBoundary == "solid_wall")
                        ABD = -gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz -
                              gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                              jac(i, j - 1, k) /
                              (jac(i, j - 1, k) + jac(i, j - 1, k - 1)) -
                              gDEdgeB * pDEdgeBGra / dy +
                              gUEdgeB * pUEdgeBGra * 0.25 * met23UEdgeB / dz
                    elseif (k == nz && zBoundary == "solid_wall")
                        ABD = -gEdgeB * pEdgeBGra * met23EdgeB / dz -
                              gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                              jac(i, j - 1, k) /
                              (jac(i, j - 1, k) + jac(i, j - 1, k - 1)) -
                              gDEdgeB * pDEdgeBGra / dy
                    else
                        ABD = -gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz -
                              gEdgeD * pEdgeDGra * 0.5 * met23EdgeD / dy *
                              jac(i, j - 1, k) /
                              (jac(i, j - 1, k) + jac(i, j - 1, k - 1)) -
                              gDEdgeB * pDEdgeBGra / dy
                    end

                    # ------------------ A(i,j,k+2) -------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AUU = -gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz -
                              gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz -
                              gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz -
                              gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz +
                              gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz +
                              gUEdgeL * pUEdgeLGra * 0.25 * met13UEdgeL / dz +
                              gUEdgeF * pUEdgeFGra * 0.25 * met23UEdgeF / dz +
                              gUEdgeB * pUEdgeBGra * 0.25 * met23UEdgeB / dz
                    elseif ((k == nz - 1 || k == nz) && zBoundary == "solid_wall")
                        AUU = 0.0
                    else
                        AUU = gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz +
                              gUEdgeL * pUEdgeLGra * 0.25 * met13UEdgeL / dz +
                              gUEdgeF * pUEdgeFGra * 0.25 * met23UEdgeF / dz +
                              gUEdgeB * pUEdgeBGra * 0.25 * met23UEdgeB / dz
                    end

                    # ------------------ A(i,j,k-2) -------------------

                    if ((k == 1 || k == 2) && zBoundary == "solid_wall")
                        ADD = 0.0
                    elseif (k == nz && zBoundary == "solid_wall")
                        ADD = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz +
                              gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz +
                              gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz +
                              gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz -
                              gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz -
                              gDEdgeL * pDEdgeLGra * 0.25 * met13DEdgeL / dz -
                              gDEdgeF * pDEdgeFGra * 0.25 * met23DEdgeF / dz -
                              gDEdgeB * pDEdgeBGra * 0.25 * met23DEdgeB / dz
                    else
                        ADD = -gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz -
                              gDEdgeL * pDEdgeLGra * 0.25 * met13DEdgeL / dz -
                              gDEdgeF * pDEdgeFGra * 0.25 * met23DEdgeF / dz -
                              gDEdgeB * pDEdgeBGra * 0.25 * met23DEdgeB / dz
                    end

                    # ----------------- A(i+1,j,k+2) ------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ARUU = -gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz +
                               gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz
                    elseif ((k == nz - 1 || k == nz) && zBoundary == "solid_wall")
                        ARUU = 0.0
                    else
                        ARUU = gUEdgeR * pUEdgeRGra * 0.25 * met13UEdgeR / dz
                    end

                    # ----------------- A(i+1,j,k-2) ------------------

                    if ((k == 1 || k == 2) && zBoundary == "solid_wall")
                        ARDD = 0.0
                    elseif (k == nz && zBoundary == "solid_wall")
                        ARDD = gEdgeR * pEdgeRGra * 0.25 * met13EdgeR / dz -
                               gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz
                    else
                        ARDD = -gDEdgeR * pDEdgeRGra * 0.25 * met13DEdgeR / dz
                    end

                    # ----------------- A(i-1,j,k+2) ------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ALUU = -gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz +
                               gUEdgeL * pUEdgeLGra * 0.25 * met13UEdgeL / dz
                    elseif ((k == nz - 1 || k == nz) && zBoundary == "solid_wall")
                        ALUU = 0.0
                    else
                        ALUU = gUEdgeL * pUEdgeLGra * 0.25 * met13UEdgeL / dz
                    end

                    # ----------------- A(i-1,j,k-2) ------------------

                    if ((k == 1 || k == 2) && zBoundary == "solid_wall")
                        ALDD = 0.0
                    elseif (k == nz && zBoundary == "solid_wall")
                        ALDD = gEdgeL * pEdgeLGra * 0.25 * met13EdgeL / dz -
                               gDEdgeL * pDEdgeLGra * 0.25 * met13DEdgeL / dz
                    else
                        ALDD = -gDEdgeL * pDEdgeLGra * 0.25 * met13DEdgeL / dz
                    end

                    # ----------------- A(i,j+1,k+2) ------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        AFUU = -gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz +
                               gUEdgeF * pUEdgeFGra * 0.25 * met23UEdgeF / dz
                    elseif ((k == nz - 1 || k == nz) && zBoundary == "solid_wall")
                        AFUU = 0.0
                    else
                        AFUU = gUEdgeF * pUEdgeFGra * 0.25 * met23UEdgeF / dz
                    end

                    # ----------------- A(i,j+1,k-2) ------------------

                    if ((k == 1 || k == 2) && zBoundary == "solid_wall")
                        AFDD = 0.0
                    elseif (k == nz && zBoundary == "solid_wall")
                        AFDD = gEdgeF * pEdgeFGra * 0.25 * met23EdgeF / dz -
                               gDEdgeF * pDEdgeFGra * 0.25 * met23DEdgeF / dz
                    else
                        AFDD = -gDEdgeF * pDEdgeFGra * 0.25 * met23DEdgeF / dz
                    end

                    # ----------------- A(i,j-1,k+2) ------------------

                    if (k == 1 && zBoundary == "solid_wall")
                        ABUU = -gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz +
                               gUEdgeB * pUEdgeBGra * 0.25 * met23UEdgeB / dz
                    elseif ((k == nz - 1 || k == nz) && zBoundary == "solid_wall")
                        ABUU = 0.0
                    else
                        ABUU = gUEdgeB * pUEdgeBGra * 0.25 * met23UEdgeB / dz
                    end

                    # ----------------- A(i,j-1,k-2) ------------------

                    if ((k == 1 || k == 2) && zBoundary == "solid_wall")
                        ABDD = 0.0
                    elseif (k == nz && zBoundary == "solid_wall")
                        ABDD = gEdgeB * pEdgeBGra * 0.25 * met23EdgeB / dz -
                               gDEdgeB * pDEdgeBGra * 0.25 * met23DEdgeB / dz
                    else
                        ABDD = -gDEdgeB * pDEdgeBGra * 0.25 * met23DEdgeB / dz
                    end

                    # Scale the tensor elements.
                    AC = AC / (fcscal^2.0)
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

                    # Set matrix elements for bicgstab.
                    ac_b[i, j, k] = AC
                    ar_b[i, j, k] = AR
                    al_b[i, j, k] = AL
                    af_b[i, j, k] = AF
                    ab_b[i, j, k] = AB
                    au_b[i, j, k] = AU
                    ad_b[i, j, k] = AD
                    aru_b[i, j, k] = ARU
                    ard_b[i, j, k] = ARD
                    alu_b[i, j, k] = ALU
                    ald_b[i, j, k] = ALD
                    afu_b[i, j, k] = AFU
                    afd_b[i, j, k] = AFD
                    abu_b[i, j, k] = ABU
                    abd_b[i, j, k] = ABD
                    auu_b[i, j, k] = AUU
                    add_b[i, j, k] = ADD
                    aruu_b[i, j, k] = ARUU
                    ardd_b[i, j, k] = ARDD
                    aluu_b[i, j, k] = ALUU
                    aldd_b[i, j, k] = ALDD
                    afuu_b[i, j, k] = AFUU
                    afdd_b[i, j, k] = AFDD
                    abuu_b[i, j, k] = ABUU
                    abdd_b[i, j, k] = ABDD

                    # Store horizontal and vertical components of AC (for
                    # preconditioner).
                    if (preconditioner == "yes")
                        ach_b[i, j, k] = -AR - AL - AF - AB
                        acv_b[i, j, k] = -AU - AD
                    end
                end
            end
        end
        kr_sp_tfc .= kr_sp_tfc / facray
        kr_sp_w_tfc .= kr_sp_w_tfc / facray
    else
        @assert false "ERROR: wrong opt"
    end

    return nothing
end
