# src/Examples/WavePacketTools/thetabar.jl

"""
```julia
thetabar(state::State, x::Real, y::Real, z::Real)::Real
```

Return the reference potential temperature of `state` at ``\\left(x, y, z\\right)``.

# Arguments

  - `state`: Auxiliary model state.

  - `x`: Zonal position.

  - `y`: Meridional position.

  - `z`: Vertical position.

# See also

  - [`PinCFlow.Examples.WavePacketTools.ijk`](@ref)
"""
function thetabar end

@ivy function thetabar(state::State, x::Real, y::Real, z::Real)::Real
    (; atmosphere) = state
    (; thetaref) = state.constants

    return atmosphere.thetabar[ijk(state, x, y, z)] .* thetaref
end
