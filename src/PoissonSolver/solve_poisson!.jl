function solve_poisson!(
    state::State,
    b::AbstractArray{<:AbstractFloat, 3},
    tolref::AbstractFloat,
    dt::AbstractFloat,
    opt::AbstractIntegration,
    model::PseudoIncompressible,
    facray::AbstractFloat,
    facprs::AbstractFloat,
)
    (; namelists, domain, grid, poisson) = state
    (; i0, i1, j0, j1, k0, k1) = domain
    (; rhostrattfc, pstrattfc) = state.atmosphere
    (; dpip) = state.variables.tendencies

    # Initialize solution and residual.
    sol = state.poisson.solution
    sol .= 0.0

    # Initialize.
    if dt == 0.0
        error("Error in solve_poisson!: dt = 0.0!")
    end
    dtinv = 1.0 / dt

    #--------------------------------
    #     Linear equation solver
    #     solve for dt * dp ...
    #--------------------------------

    compute_operator!(state, dt, opt, facray)

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
    dpip[i0:i1, j0:j1, k0:k1] .= dtinv ./ facprs .* sol

    return (errflagbicg, niterbicg)
end
