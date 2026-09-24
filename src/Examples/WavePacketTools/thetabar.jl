# src/Examples/WavePacketTools/thetabar.jl

@ivy function thetabar(state::State, x::Real, y::Real, z::Real)::Real
    (; atmosphere) = state
    (; thetaref) = state.constants

    atmosphere.thetabar[ijk(state, x, y, z)] .* thetaref
end
