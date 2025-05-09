function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
    (; model) = state.namelists.setting
    return compute_compressible_buoyancy_factor(state, indices, variable, model)
end

function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::Compressible,
)
    return 1.0
end

function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::AbstractModel,
)
    (; rhostrattfc) = state.atmosphere
    (; rho) = state.variables.predictands
    (ix, jy, kz) = indices
    return rhostrattfc[ix, jy, kz] / (rho[ix, jy, kz] + rhostrattfc[ix, jy, kz])
end

function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::Compressible,
)
    return 1.0
end

function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::AbstractModel,
)
    (; jac) = state.grid
    (; rhostrattfc) = state.atmosphere
    (; rho) = state.variables.predictands
    (ix, jy, kz) = indices
    return (
        jac[ix, jy, kz + 1] * rhostrattfc[ix, jy, kz] +
        jac[ix, jy, kz] * rhostrattfc[ix, jy, kz + 1]
    ) / (
        jac[ix, jy, kz + 1] * (rho[ix, jy, kz] + rhostrattfc[ix, jy, kz]) +
        jac[ix, jy, kz] * (rho[ix, jy, kz + 1] + rhostrattfc[ix, jy, kz + 1])
    )
end
