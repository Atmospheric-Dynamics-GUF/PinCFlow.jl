# src/Examples/WavePacketTools/thetabar.jl

function thetabar(state::State, x::Real, y::Real, z::Real)::Real
    (; atmosphere) = state
    (; thetaref) = state.constants

    @ivy return atmosphere.thetabar[ijk(state, x, y, z)] .* thetaref
end
