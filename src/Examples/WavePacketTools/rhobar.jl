# src/Examples/WavePacketTools/rhobar.jl

"""
```julia
rhobar(state::State, x::Real, y::Real, z::Real)::Real
```

Return the reference density of `state` at ``\\left(x, y, z\\right)``.

# Arguments

  - `state`: Auxiliary model state.

  - `x`: Zonal position.

  - `y`: Meridional position.

  - `z`: Vertical position.

# See also

  - [`PinCFlow.Examples.WavePacketTools.ijk`](@ref)
"""
function rhobar end

@ivy function rhobar(state::State, x::Real, y::Real, z::Real)::Real
    (; atmosphere) = state
    (; rhoref) = state.constants

    return atmosphere.rhobar[ijk(state, x, y, z)] .* rhoref
end
