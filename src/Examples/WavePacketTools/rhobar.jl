# src/Examples/WavePacketTools/rhobar.jl

function rhobar(state::State, x::Real, y::Real, z::Real)::Real
    (; atmosphere) = state
    (; rhoref) = state.constants

    @ivy return atmosphere.rhobar[ijk(state, x, y, z)] .* rhoref
end
