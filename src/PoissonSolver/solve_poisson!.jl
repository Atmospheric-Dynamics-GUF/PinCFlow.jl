"""
```julia
solve_poisson!(
    state::State,
    b::AbstractArray{<:AbstractFloat, 3},
    tolref::AbstractFloat,
    dt::AbstractFloat,
    rayleigh_factor::AbstractFloat,
)::Tuple{Bool, <:Integer}
```

Solve the Poisson equation and return a tuple containing an error flag and the number of iterations.

Given a right-hand side and reference tolerance, this method computes the elements of the linear operator and solves the Poisson equation, using a preconditioned BicGStab algorithm. Both the Exner-pressure differences and the entire equation are scaled with ``\\sqrt{\\overline{\\rho}} / P`` in advance (the right-hand side has already been scaled at this point), so that the equation

```math
\\frac{\\sqrt{\\overline{\\rho}}}{P} \\mathrm{LHS} \\left(\\frac{\\sqrt{\\overline{\\rho}}}{P} s\\right) = \\frac{\\sqrt{\\overline{\\rho}}}{P} \\mathrm{RHS}
```

is solved for ``s``. The Exner-pressure differnces are then given by ``\\Delta \\pi = \\left(\\sqrt{\\overline{\\rho}} / P\\right) \\left(s / \\Delta t\\right)``.

# Arguments

  - `state`: Model state.

  - `b`: Right-hand side.

  - `tolref`: Reference tolerance for convergence criterion.

  - `dt`: Time step.

  - `rayleigh_factor`: Factor by which the Rayleigh-damping coefficient is multiplied.

# See also

  - [`PinCFlow.PoissonSolver.compute_operator!`](@ref)

  - [`PinCFlow.PoissonSolver.apply_bicgstab!`](@ref)
"""
function solve_poisson! end

function solve_poisson!(
    state::State,
    b::AbstractArray{<:AbstractFloat, 3},
    tolref::AbstractFloat,
    dt::AbstractFloat,
    rayleigh_factor::AbstractFloat,
)::Tuple{Bool, <:Integer}
    (; namelists, domain, grid, poisson) = state
    (; model) = namelists.setting
    (; i0, i1, j0, j1, k0, k1) = domain
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; dpip) = state.variables.tendencies

    sol = state.poisson.solution
    sol .= 0.0

    if dt == 0.0
        error("Error in solve_poisson!: dt = 0.0!")
    end
    dtinv = 1.0 / dt

    compute_operator!(state, dt, rayleigh_factor)

    (errflagbicg, niterbicg) =
        apply_bicgstab!(b, tolref, sol, namelists, domain, grid, poisson)

    if errflagbicg
        return (errflagbicg, niterbicg)
    end

    for k in k0:k1, j in j0:j1, i in i0:i1
        fcscal = sqrt(pstrattfc[i, j, k]^2 / rhostrattfc[i, j, k])
        sol[i - i0 + 1, j - j0 + 1, k - k0 + 1] /= fcscal
    end

    # Pass solution to pressure correction.
    dpip[i0:i1, j0:j1, k0:k1] .= dtinv .* sol

    return (errflagbicg, niterbicg)
end
