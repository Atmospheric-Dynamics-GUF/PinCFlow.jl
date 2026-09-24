# src/Examples/WavePacketTools/ijk.jl

"""
```julia
ijk(state::State, x::Real, y::Real, z::Real)::CartesianIndex
```

Return the indices of the grid point in `state` that is nearest to ``\\left(x, y, z\\right)``.

# Arguments

  - `state`: Auxiliary model state.

  - `x`: Zonal position.

  - `y`: Meridional position.

  - `z`: Vertical position.
"""
function ijk end

@ivy function ijk(state::State, x::Real, y::Real, z::Real)::CartesianIndex
    (; lref) = state.constants
    (; grid) = state

    i = argmin(abs.(x .- grid.x .* lref))
    j = argmin(abs.(y .- grid.y .* lref))
    k = argmin(abs.(z .- grid.zc[i, j, :] .* lref))

    return CartesianIndex(i, j, k)
end
