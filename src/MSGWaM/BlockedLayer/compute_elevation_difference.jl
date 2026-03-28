"""
```julia
compute_elevation_difference(
    state::State,
    i::Integer,
    j::Integer,
)::AbstractFloat
```

Return ``\\Delta h = \\sum_\\alpha \\left|h_{\\mathrm{w}, \\alpha}\\right|`` as an approximation of the elevation difference between the local background orography and local summits of the true orography.

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid index.

  - `j`: Meridional grid index.

!!! danger "Experimental"
    The blocked-layer scheme is an experimental feature that hasn't been validated yet.
"""
function compute_elevation_difference end

function compute_elevation_difference(
    state::State,
    i::Integer,
    j::Integer,
)::AbstractFloat
    (; hw) = state.grid
    @ivy return sum(abs, hw[:, i, j])
end
