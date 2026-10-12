# src/Examples/WavePacketTools/n2.jl

"""
```julia
n2(state::State, x::Real, y::Real, z::Real)::Real
```

Return the squared buoyancy frequency of `state` at ``\\left(x, y, z\\right)``.

# Arguments

  - `state`: Auxiliary model state.

  - `x`: Zonal position.

  - `y`: Meridional position.

  - `z`: Vertical position.

# See also

  - [`PinCFlow.Examples.WavePacketTools.ijk`](@ref)
"""
function n2 end

@ivy function n2(state::State, x::Real, y::Real, z::Real)::Real
    (; atmosphere) = state
    (; tref) = state.constants

    return atmosphere.n2[ijk(state, x, y, z)] ./ tref .^ 2
end
