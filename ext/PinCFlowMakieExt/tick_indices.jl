"""
```julia
tick_indices(levels::AbstractVector{<:Real})::AbstractVector{<:Integer}
```

Return three symmetric indices for `levels`.

# Arguments

  - `levels`: Abstract input vector for which the symmetric indices are determined.
"""
function tick_indices end

function tick_indices(levels::AbstractVector{<:Real})::AbstractVector{<:Integer}
    n = length(levels)
    return n % 2 == 0 ? [1, fld(n, 3) + 1, n - fld(n, 3), n] : [1, cld(n, 2), n]
end
