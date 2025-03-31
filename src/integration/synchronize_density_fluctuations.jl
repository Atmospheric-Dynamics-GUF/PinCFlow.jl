function synchronize_density_fluctuations!(state::State)
  (; model) = state.namelists.setting
  synchronize_density_fluctuations!(state, model)
  return
end

function synchronize_density_fluctuations!(
  state::State,
  model::PseudoIncompressible,
)

  # Get all necessary fields.
  (; rho, rhop) = state.variables.predictands

  # Synchronize density fluctuations.
  rhop .= rho

  # Return.
  return
end
