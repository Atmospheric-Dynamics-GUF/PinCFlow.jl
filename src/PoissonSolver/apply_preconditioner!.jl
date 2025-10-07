"""
```julia
apply_preconditioner!(
    sin::AbstractArray{<:AbstractFloat, 3},
    sout::AbstractArray{<:AbstractFloat, 3},
    state::State,
)
```

Apply a preconditioner to the Poisson problem.

This preconditioner integrates the auxiliary equation

```math
\\frac{\\mathrm{d} s}{\\mathrm{d} \\eta} = \\mathcal{L}_\\mathrm{h} \\left(s\\right) + \\mathcal{L}_\\mathrm{v} \\left(s\\right) - b,
```

where ``s`` is the iterative solution, ``\\eta`` is a pseudo-time variable, ``\\mathcal{L}_\\mathrm{v}`` contains the lower, center and upper diagonals of the linear operator, ``\\mathcal{L}_\\mathrm{h}`` contains all remaining elements, and ``b`` is the left-hand side. The integration is performed in a semi-implicit manner, following

```math
\\left(1 - \\Delta \\eta \\mathcal{L}_\\mathrm{v}\\right) \\left(s^{\\left(m + 1\\right)}\\right) = \\left(1 + \\Delta \\eta \\mathcal{L}_\\mathrm{h}\\right) \\left(s^{\\left(m\\right)}\\right) - \\Delta \\eta b,
```

where ``\\Delta \\eta = \\Delta \\tau / 2 \\left[\\left(\\Delta \\widehat{x}\\right)^{- 2} + \\left(\\Delta \\widehat{y}\\right)^{- 2}\\right]^{- 1}``, with ``\\Delta \\tau`` being a namelist parameter (`state.namelist.poisson.dtau`). Therein, the implicit problem is solved with the Thomas algorithm for tridiagonal matrices. The number of iterations is given by `state.namelist.poisson.preconditioner_iterations`. Since the Thomas algorithm consists of an upward elimination sweep and a downward pass, this method performs sequential one-way MPI communication if the domain is parallelized in the vertical.

# Arguments

  - `sin`: Residual array.

  - `sout`: Solution of the preconditioner.

  - `state`: Model state.

# See also

  - [`PinCFlow.PoissonSolver.apply_operator!`](@ref)
"""
function apply_preconditioner! end

function apply_preconditioner!(
    sin::AbstractArray{<:AbstractFloat, 3},
    sout::AbstractArray{<:AbstractFloat, 3},
    state::State,
)
    (; dtau, preconditioner_iterations) = state.namelists.poisson
    (; comm, zz_size, nz, nzz, ko, down, up) = state.domain
    (; dx, dy) = state.grid
    (; au_b, ac_b, ad_b) = state.poisson.tensor
    (; s_pc, q_pc, p_pc, s_pc_bc, q_pc_bc) = state.poisson.preconditioner

    # Initialize auxiliary fields.
    s_pc .= 0.0
    q_pc .= 0.0
    p_pc .= 0.0

    # Set pseudo-time step.
    deta = dtau / (2 * (1 / dx^2 + 1 / dy^2))

    # Iterate.
    @ivy for niter in 1:preconditioner_iterations
        apply_operator!(s_pc, q_pc, Horizontal(), state)
        s_pc .+= deta .* (q_pc .- sin)

        # Set the lower boundary.
        if ko == 0
            q_pc[:, :, 1] .=
                deta .* au_b[:, :, 1] ./ (1 .- deta .* ac_b[:, :, 1])
            s_pc[:, :, 1] ./= 1 .- deta .* ac_b[:, :, 1]
        else
            MPI.Recv!(q_pc_bc, comm; source = down, tag = 1)
            MPI.Recv!(s_pc_bc, comm; source = down, tag = 2)

            p_pc .=
                1 ./
                (1 .- deta .* ac_b[:, :, 1] .- deta .* ad_b[:, :, 1] .* q_pc_bc)
            q_pc[:, :, 1] .= deta .* au_b[:, :, 1] .* p_pc
            s_pc[:, :, 1] .=
                (s_pc[:, :, 1] .+ deta .* ad_b[:, :, 1] .* s_pc_bc) .* p_pc
        end

        # Perform upward sweep.
        for k in 2:nz
            p_pc .=
                1 ./ (
                    1 .- deta .* ac_b[:, :, k] .-
                    deta .* ad_b[:, :, k] .* q_pc[:, :, k - 1]
                )
            q_pc[:, :, k] .= deta .* au_b[:, :, k] .* p_pc
            s_pc[:, :, k] .=
                (s_pc[:, :, k] .+ deta .* ad_b[:, :, k] .* s_pc[:, :, k - 1]) .*
                p_pc
        end

        # Communicate the upper boundary and set it for the downward sweep.
        if ko + nzz != zz_size
            q_pc_bc .= q_pc[:, :, nz]
            s_pc_bc .= s_pc[:, :, nz]

            MPI.Send(q_pc_bc, comm; dest = up, tag = 1)
            MPI.Send(s_pc_bc, comm; dest = up, tag = 2)

            MPI.Recv!(s_pc_bc, comm; source = up)

            s_pc[:, :, nz] .+= q_pc[:, :, nz] .* s_pc_bc
        end

        # Perform downward sweep.
        for k in (nz - 1):-1:1
            s_pc[:, :, k] .+= q_pc[:, :, k] .* s_pc[:, :, k + 1]
        end

        # Communicate the lower boundary.
        if ko != 0
            s_pc_bc .= s_pc[:, :, 1]

            MPI.Send(s_pc_bc, comm; dest = down)
        end
    end

    sout .= s_pc

    return
end
