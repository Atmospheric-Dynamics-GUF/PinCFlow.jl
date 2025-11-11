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
    (; z_size) = state.namelists.domain
    (; dtau, preconditioner_iterations) = state.namelists.poisson
    (; dx, dy) = state.grid
    (; au_b, ac_b, ad_b) = state.poisson.tensor
    (; ath, bth, cth, fth, qth, pth, qth_bc, fth_bc) =
        state.variables.auxiliaries

    reset_thomas!(state)

    # Set pseudo-time step.
    deta = dtau / (2 * (1 / dx^2 + 1 / dy^2))

    ath .= .-deta .* ad_b
    bth .= 1.0 .- deta .* ac_b
    cth .= .-deta .* au_b

    # Iterate.
    @ivy for niter in 1:preconditioner_iterations
        apply_operator!(fth, qth, Horizontal(), state)
        fth .+= deta .* (qth .- sin)

        thomas_algorithm!(state, ath, bth, cth, fth, qth, pth, fth_bc, qth_bc)
    end

    sout .= fth

    return
end
