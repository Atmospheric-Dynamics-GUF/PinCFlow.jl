"""
```julia
solve_poisson!(
    state::State,
    b::AbstractArray{<:AbstractFloat, 3},
    tolref::AbstractFloat,
    dt::AbstractFloat,
    facray::AbstractFloat,
    facprs::AbstractFloat,
)
```

Solve the Poisson equation for pressure correction in the atmospheric flow solver.

This function performs the core pressure correction step by solving the discretized
Poisson equation using BiCGStab iterative method. The solution is scaled appropriately
for different equation models (Boussinesq, PseudoIncompressible, Compressible).

# Arguments

  - `state`: Complete simulation state containing all field variables and parameters
  - `b`: Right-hand side of the Poisson equation
  - `tolref`: Reference tolerance for convergence checking
  - `dt`: Time step size
  - `facray`: Rayleigh damping factor for sponge layers
  - `facprs`: Pressure correction factor

# Returns

  - `(errflagbicg, niterbicg)`: Error flag and number of iterations from BiCGStab solver

# Notes

  - The solution is stored in `state.poisson.solution` and then transferred to pressure correction field
  - For non-Boussinesq models, the solution is scaled by the compressibility factor
  - The function will error if `dt = 0.0` to prevent division by zero

# See also

  - [`PinCFlow.PoissonSolver.compute_operator!`](@ref)
  - [`PinCFlow.PoissonSolver.apply_bicgstab!`](@ref)
"""
function solve_poisson!(
    state::State,
    b::AbstractArray{<:AbstractFloat, 3},
    tolref::AbstractFloat,
    dt::AbstractFloat,
    facray::AbstractFloat,
    facprs::AbstractFloat,
)
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

    compute_operator!(state, dt, facray)

    (errflagbicg, niterbicg) =
        apply_bicgstab!(b, tolref, sol, namelists, domain, grid, poisson)

    if errflagbicg
        return (errflagbicg, niterbicg)
    end

    if model != Boussinesq()
        for k in k0:k1, j in j0:j1, i in i0:i1
            fcscal = sqrt(pstrattfc[i, j, k]^2 / rhostrattfc[i, j, k])
            sol[i - i0 + 1, j - j0 + 1, k - k0 + 1] /= fcscal
        end
    end

    # Pass solution to pressure correction.
    dpip[i0:i1, j0:j1, k0:k1] .= dtinv ./ facprs .* sol

    return (errflagbicg, niterbicg)
end
