abstract type AbstractBoundaryVariables end
struct BoundaryPredictands <: AbstractBoundaryVariables end
struct BoundaryReconstructions <: AbstractBoundaryVariables end
struct BoundaryFluxes <: AbstractBoundaryVariables end

function set_boundaries!(state::State, variables::BoundaryPredictands)
  (; xboundaries, yboundaries, zboundaries) = state.namelists.boundaries
  set_zonal_boundaries!(state, variables, xboundaries)
  set_meridional_boundaries!(state, variables, yboundaries)
  set_vertical_boundaries!(state, variables, zboundaries)
  return
end

function set_boundaries!(state::State, variables::BoundaryReconstructions)
  (; xboundaries, yboundaries, zboundaries) = state.namelists.boundaries
  set_zonal_boundaries!(state, variables, xboundaries)
  set_meridional_boundaries!(state, variables, yboundaries)
  return
end

function set_boundaries!(state::State, variables::BoundaryFluxes)
  (; xboundaries, yboundaries, zboundaries) = state.namelists.boundaries
  set_vertical_boundaries!(state, variables, zboundaries)
  return
end
