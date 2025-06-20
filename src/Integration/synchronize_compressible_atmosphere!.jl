"""
    synchronize_compressible_atmosphere!(state::State, predictands::Predictands)

Synchronize compressible atmosphere fields based on model type.

This function dispatches to the appropriate model-specific implementation
for updating atmospheric stratification fields in compressible models.

# Arguments

  - `state::State`: Complete simulation state
  - `predictands::Predictands`: Current values of predicted variables
"""
function synchronize_compressible_atmosphere!(
    state::State,
    predictands::Predictands,
)
    (; model) = state.namelists.setting
    synchronize_compressible_atmosphere!(state, predictands, model)
    return
end

"""
    synchronize_compressible_atmosphere!(state::State, predictands::Predictands, model::AbstractModel)

No-op for non-compressible models.

Most model types don't require atmospheric synchronization as they don't
maintain stratified background atmosphere fields.

# Arguments

  - `state::State`: Simulation state (unused)
  - `predictands::Predictands`: Predicted variables (unused)
  - `model::AbstractModel`: Non-compressible model type
"""
function synchronize_compressible_atmosphere!(
    state::State,
    predictands::Predictands,
    model::AbstractModel,
)
    return
end

"""
    synchronize_compressible_atmosphere!(state::State, predictands::Predictands, model::Compressible)

Synchronize stratified atmosphere fields for compressible model.

Updates the background stratification fields based on current predictand values,
including pressure stratification and Brunt-Väisälä frequency computation.

# Arguments

  - `state::State`: Simulation state containing grid, domain, and atmosphere data
  - `predictands::Predictands`: Current predicted field values
  - `model::Compressible`: Compressible model type

# Details

  - Updates pressure stratification: `pstrattfc = p`
  - Computes Brunt-Väisälä frequency using current density and pressure fields
  - Applies vertical boundary conditions for the computed stratification
  - Handles special cases for boundary processes in MPI decomposition
"""
function synchronize_compressible_atmosphere!(
    state::State,
    predictands::Predictands,
    model::Compressible,
)
    (; namelists, domain) = state
    (; nbz) = namelists.domain
    (; zboundaries) = namelists.setting
    (; g_ndim) = state.constants
    (; sizezz, nzz, ko, k0, k1) = domain
    (; dz, jac) = state.grid
    (; pstrattfc, bvsstrattfc, thetastrattfc, rhostrattfc) = state.atmosphere
    (; rho, p) = predictands

    pstrattfc .= p

    for kz in k0:k1
        bvsstrattfc[:, :, kz] .=
            g_ndim .* p[:, :, kz] ./ (rho[:, :, kz] .+ rhostrattfc[:, :, kz]) ./
            (thetastrattfc[:, :, kz] .^ 2) ./ jac[:, :, kz] .*
            (thetastrattfc[:, :, kz + 1] .- thetastrattfc[:, :, kz - 1]) ./ 2 ./
            dz
    end

    set_vertical_boundaries_of_field!(
        bvsstrattfc,
        namelists,
        domain,
        zboundaries,
        +,
    )
    if ko == 0
        for k in 1:nbz
            bvsstrattfc[:, :, k] .=
                g_ndim .* p[:, :, k0 - 1] ./
                (rho[:, :, k0 - 1] .+ rhostrattfc[:, :, k0 - 1]) ./
                (thetastrattfc[:, :, k0 - 1] .^ 2) ./ jac[:, :, k0 - 1] .*
                (thetastrattfc[:, :, k0] .- thetastrattfc[:, :, k0 - 1]) ./ dz
        end
    end
    if ko + nzz == sizezz
        for k in 1:nbz
            bvsstrattfc[:, :, k1 + k] .=
                g_ndim .* p[:, :, k1 + 1] ./
                (rho[:, :, k1 + 1] .+ rhostrattfc[:, :, k1 + 1]) ./
                (thetastrattfc[:, :, k1 + 1] .^ 2) ./ jac[:, :, k1 + 1] .*
                (thetastrattfc[:, :, k1 + 1] .- thetastrattfc[:, :, k1]) ./ dz
        end
    end

    return
end
