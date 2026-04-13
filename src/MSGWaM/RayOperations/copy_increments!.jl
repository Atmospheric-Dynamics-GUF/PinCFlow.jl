"""
```julia
copy_increments!(
    increments::WKBIncrements,
    r::Pair{<:Integer, <:Integer},
    i::Pair{<:Integer, <:Integer},
    j::Pair{<:Integer, <:Integer},
    k::Pair{<:Integer, <:Integer},
)
```

Copy all increments of the ray volume specified by the first components of the index pairs to that specified by the second components.

# Arguments

  - `increments`: Collection of increments for ray-volume properties.

  - `r`: Ray-volume indices.

  - `i`: Zonal grid-cell indices.

  - `j`: Meridional grid-cell indices.

  - `k`: Vertical grid-cell indices.
"""
function copy_increments! end

function copy_increments!(
    increments::WKBIncrements,
    r::Pair{<:Integer, <:Integer},
    i::Pair{<:Integer, <:Integer},
    j::Pair{<:Integer, <:Integer},
    k::Pair{<:Integer, <:Integer},
)
    @ivy for field in fieldnames(WKBIncrements)
        getfield(increments, field)[r[2], i[2], j[2], k[2]] =
            getfield(increments, field)[r[1], i[1], j[1], k[1]]
    end

    return
end
