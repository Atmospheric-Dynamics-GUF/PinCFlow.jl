"""
```julia
solve_poisson!(
    state::State,
    dt::AbstractFloat,
    rayleigh_factor::AbstractFloat,
    tolref::AbstractFloat,
)::Tuple{Bool, <:Integer}
```

Solve the Poisson equation and return a tuple containing an error flag and the number of iterations.

Given a left-hand side and reference tolerance, this method computes the elements of the linear operator and solves the Poisson equation, using a preconditioned BicGStab algorithm. Both the Exner-pressure differences and the entire equation are scaled with ``\\sqrt{\\overline{\\rho}} / P`` in advance (the left-hand side has already been scaled at this point), so that the equation

```math
\\frac{\\sqrt{\\overline{\\rho}}}{P} \\mathrm{LHS} = \\frac{\\sqrt{\\overline{\\rho}}}{P} \\mathrm{RHS} \\left(\\frac{\\sqrt{\\overline{\\rho}}}{P} s\\right)
```

is solved for ``s``. The Exner-pressure differences are then given by ``\\Delta \\pi = \\left(\\sqrt{\\overline{\\rho}} / P\\right) \\left(s / \\Delta t\\right)``.

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `rayleigh_factor`: Factor by which the Rayleigh-damping coefficient is multiplied.

  - `tolref`: Reference tolerance for convergence criterion.

# See also

  - [`PinCFlow.PoissonSolver.compute_operator!`](@ref)

  - [`PinCFlow.PoissonSolver.apply_bicgstab!`](@ref)
"""
function solve_poisson! end

function solve_poisson!(
    state::State,
    dt::AbstractFloat,
    rayleigh_factor::AbstractFloat,
    tolref::AbstractFloat,
)::Tuple{Bool, <:Integer}
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; dpip) = state.variables.increments
    (; solution) = state.poisson

    solution .= 0.0

    if dt == 0.0
        error("Error in solve_poisson!: dt = 0.0!")
    end
    dtinv = 1.0 / dt

    compute_operator!(state, dt, rayleigh_factor)

    (errflagbicg, niterbicg) = apply_bicgstab!(state, tolref)

    if errflagbicg
        return (errflagbicg, niterbicg)
    end

    ii = i0:i1
    jj = j0:j1
    kk = k0:k1

    @ivy solution ./=
        sqrt.(pstrattfc[ii, jj, kk] .^ 2 ./ rhostrattfc[ii, jj, kk])

    # Pass solution to pressure correction.
    @ivy dpip[ii, jj, kk] .= dtinv .* solution

    return (errflagbicg, niterbicg)
end
