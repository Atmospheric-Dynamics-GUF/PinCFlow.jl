"""
```julia
compute_slope(state::State, i::Integer, j::Integer)::NTuple{2, <:AbstractFloat}
```

Return ``\\boldsymbol{k}_h = \\left(\\sum_\\alpha \\left|h_{\\mathrm{w}, \\alpha}\\right| \\boldsymbol{k}_{h, \\alpha}\\right) / \\left(\\sum_\\alpha \\left|h_{\\mathrm{w}, \\alpha}\\right|\\right)`` as an approximation of the local slope of the true orography.

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid index.

  - `j`: Meridional grid index.

!!! danger "Experimental"
    The blocked-layer scheme is an experimental feature that hasn't been validated yet.
"""
function compute_slope end

@ivy function compute_slope(
    state::State,
    i::Integer,
    j::Integer,
)::NTuple{2, <:AbstractFloat}
    (; hw, kh, lh) = state.grid

    return (
        mapreduce((a, b) -> abs(a) * b, +, hw[:, i, j], kh[:, i, j]) /
        sum(abs, hw[:, i, j]),
        mapreduce((a, b) -> abs(a) * b, +, hw[:, i, j], lh[:, i, j]) /
        sum(abs, hw[:, i, j]),
    )
end
