# src/Examples/WavePacketTools/n2.jl

@ivy function n2(state::State, x::Real, y::Real, z::Real)::Real
    (; atmosphere) = state
    (; tref) = state.constants

    atmosphere.n2[ijk(state, x, y, z)] ./ tref .^ 2
end
