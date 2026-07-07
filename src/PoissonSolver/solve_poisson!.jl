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

Given a left-hand side and reference tolerance, this method computes the elements of the linear operator and solves the Poisson equation, using a preconditioned BiCGSTAB algorithm. Both the Exner-pressure differences and the entire equation are scaled with ``\\sqrt{\\bar{\\rho}} / P`` in advance (the left-hand side has already been scaled at this point), so that the equation

```math
\\frac{\\sqrt{\\bar{\\rho}}}{P} \\mathrm{LHS} = \\frac{\\sqrt{\\bar{\\rho}}}{P} \\mathrm{RHS} \\left(\\frac{\\sqrt{\\bar{\\rho}}}{P} s\\right)
```

is solved for ``s``. The Exner-pressure differences are then given by ``\\Delta \\pi = \\left(\\sqrt{\\bar{\\rho}} / P\\right) \\left(s / \\Delta t\\right)``.

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

@ivy function solve_poisson!(
    state::State,
    dt::AbstractFloat,
    rayleigh_factor::AbstractFloat,
    tolref::AbstractFloat,
)::Tuple{Bool, <:Integer}
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; rhobar, pbar) = state.atmosphere
    (; dpip) = state.variables.increments
    (; solution) = state.poisson

    @share solution = 0.0

    if dt == 0.0
        error("Vanishing time step: dt = 0.0!")
    end
    dtinv = 1.0 / dt

    compute_operator!(state, dt, rayleigh_factor)

    (errflagbicg, niterbicg) = apply_bicgstab!(state, tolref)

    if errflagbicg
        return (errflagbicg, niterbicg)
    end

    @share for k in k0:k1, j in j0:j1, i in i0:i1
        is = i - i0 + 1
        js = j - j0 + 1
        ks = k - k0 + 1

        solution[is, js, ks] /= sqrt(pbar[i, j, k]^2 / rhobar[i, j, k])
        dpip[i, j, k] = dtinv * solution[is, js, ks]
    end

    return (errflagbicg, niterbicg)
end
