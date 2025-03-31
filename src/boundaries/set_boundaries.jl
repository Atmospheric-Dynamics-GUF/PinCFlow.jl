struct BoundaryPredictands end
struct BoundaryReconstructions end
struct BoundaryFluxes end

function set_boundaries!(state::State, variables::BoundaryPredictands)
  (; zboundaries) = state.namelists.setting
  set_zonal_boundaries!(state, variables)
  set_meridional_boundaries!(state, variables)
  set_vertical_boundaries!(state, variables, zboundaries)
  return
end

function set_boundaries!(state::State, variables::BoundaryReconstructions)
  set_zonal_boundaries!(state, variables)
  set_meridional_boundaries!(state, variables)
  return
end

function set_boundaries!(state::State, variables::BoundaryFluxes)
  (; zboundaries) = state.namelists.setting
  set_vertical_boundaries!(state, variables, zboundaries)
  return
end
