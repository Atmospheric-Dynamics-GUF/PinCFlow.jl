struct BoundaryPredictands end
struct BoundaryReconstructions end
struct BoundaryFluxes end

function set_boundaries!(state::State, variables::BoundaryPredictands)
  (; xboundaries, yboundaries, zboundaries) = state.namelists.boundaries
  set_zonal_boundaries!(state, variables, xboundaries)
  set_meridional_boundaries!(state, variables, yboundaries)
  set_vertical_boundaries!(state, variables, zboundaries)
  return
end

function set_boundaries!(state::State, variables::BoundaryReconstructions)
  (; xboundaries, yboundaries) = state.namelists.boundaries
  set_zonal_boundaries!(state, variables, xboundaries)
  set_meridional_boundaries!(state, variables, yboundaries)
  return
end

function set_boundaries!(state::State, variables::BoundaryFluxes)
  (; zboundaries) = state.namelists.boundaries
  set_vertical_boundaries!(state, variables, zboundaries)
  return
end
