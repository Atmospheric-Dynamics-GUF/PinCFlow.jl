"""
    synchronize_density_fluctuations!(state::State)

Synchronize density fluctuations based on the model type.

This function dispatches to the appropriate model-specific implementation
for synchronizing density fluctuations in the simulation state.

# Arguments

  - `state::State`: Complete simulation state containing model configuration
"""
function synchronize_density_fluctuations!(state::State)
    (; model) = state.namelists.setting
    synchronize_density_fluctuations!(state, model)
    return
end

"""
    synchronize_density_fluctuations!(state::State, model::Boussinesq)

No-op for Boussinesq model.

For Boussinesq approximation, density fluctuations don't require synchronization
as density is assumed constant except in buoyancy terms.

# Arguments

  - `state::State`: Simulation state (unused)
  - `model::Boussinesq`: Boussinesq model type
"""
function synchronize_density_fluctuations!(state::State, model::Boussinesq)
    return
end

"""
    synchronize_density_fluctuations!(state::State, model::PseudoIncompressible)

Synchronize density fluctuations for pseudo-incompressible model.

For pseudo-incompressible flow, the density fluctuation field `rhop` is
synchronized with the total density field `rho`.

# Arguments

  - `state::State`: Simulation state containing density fields
  - `model::PseudoIncompressible`: Pseudo-incompressible model type
"""
function synchronize_density_fluctuations!(
    state::State,
    model::PseudoIncompressible,
)
    (; rho, rhop) = state.variables.predictands

    rhop .= rho

    return
end

"""
    synchronize_density_fluctuations!(state::State, model::Compressible)

Synchronize density fluctuations for compressible model.

For compressible flow, the density fluctuation field `rhop` is computed
as the deviation from the hydrostatic balance, accounting for stratified
background atmosphere conditions.

# Arguments

  - `state::State`: Simulation state containing density and atmosphere fields
  - `model::Compressible`: Compressible model type

# Details

Computes: `rhop = rho + rhostrattfc - pstrattfc / thetastrattfc`
where the correction term represents deviations from hydrostatic equilibrium.
"""
function synchronize_density_fluctuations!(state::State, model::Compressible)
    (; rhostrattfc, thetastrattfc, pstrattfc) = state.atmosphere
    (; rho, rhop) = state.variables.predictands

    rhop .= rho .+ rhostrattfc .- pstrattfc ./ thetastrattfc

    return
end
