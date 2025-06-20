"""
    compute_compressible_buoyancy_factor(state, indices, variable)

Compute compressible buoyancy factor for the given variable at specified grid indices.
Dispatches based on the model type from the state configuration.

# Arguments

  - `state::State`: Current simulation state
  - `indices::NTuple{3, <:Integer}`: Grid indices (ix, jy, kz)
  - `variable::AbstractVariable`: Variable type for which to compute the factor

# Returns

  - `AbstractFloat`: Buoyancy scaling factor
"""
function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
    (; model) = state.namelists.setting
    return compute_compressible_buoyancy_factor(state, indices, variable, model)
end

"""
    compute_compressible_buoyancy_factor(state, indices, variable::RhoP, model::Compressible)

Return unity buoyancy factor for density perturbations in fully compressible model.

# Returns

  - `Float64`: Always returns 1.0
"""
function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::Compressible,
)
    return 1.0
end

"""
    compute_compressible_buoyancy_factor(state, indices, variable::RhoP, model::AbstractModel)

Compute buoyancy factor for density perturbations in non-compressible models.
Factor represents the ratio of reference density to total density.

# Returns

  - `AbstractFloat`: ρ₀/(ρ + ρ₀) where ρ₀ is reference density and ρ is density perturbation
"""
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

"""
    compute_compressible_buoyancy_factor(state, indices, variable::W, model::Compressible)

Return unity buoyancy factor for vertical velocity in fully compressible model.

# Returns

  - `Float64`: Always returns 1.0
"""
function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::Compressible,
)
    return 1.0
end

"""
    compute_compressible_buoyancy_factor(state, indices, variable::W, model::AbstractModel)

Compute buoyancy factor for vertical velocity in non-compressible models.
Uses Jacobian-weighted interpolation between vertical levels.

# Returns

  - `AbstractFloat`: Jacobian-weighted factor accounting for grid stretching and density variation
"""
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
