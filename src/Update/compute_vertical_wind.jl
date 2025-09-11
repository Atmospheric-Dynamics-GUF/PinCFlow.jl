"""
```julia
compute_vertical_wind(
    i::Integer,
    j::Integer,
    k::Integer,
    state::State,
)::AbstractFloat
```

Compute and return the Cartesian vertical wind at the grid point `(i, j, k + 1 / 2)`.

# Arguments

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `state`: Model state.

# See also

  - [`PinCFlow.Update.transform`](@ref)
"""
function compute_vertical_wind end

function compute_vertical_wind(
    i::Integer,
    j::Integer,
    k::Integer,
    state::State,
)::AbstractFloat
    (; u, v, w) = state.variables.predictands

    @ivy uedger = u[i, j, k]
    @ivy uuedger = u[i, j, k + 1]
    @ivy uedgel = u[i - 1, j, k]
    @ivy uuedgel = u[i - 1, j, k + 1]
    @ivy vedgef = v[i, j, k]
    @ivy vuedgef = v[i, j, k + 1]
    @ivy vedgeb = v[i, j - 1, k]
    @ivy vuedgeb = v[i, j - 1, k + 1]
    @ivy wedgeu = w[i, j, k]

    return transform(
        i,
        j,
        k,
        uedger,
        uuedger,
        uedgel,
        uuedgel,
        vedgef,
        vuedgef,
        vedgeb,
        vuedgeb,
        wedgeu,
        Cartesian(),
        state,
    )
end
