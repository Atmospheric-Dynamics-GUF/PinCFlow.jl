"""
```julia
get_next_level(
    i::Integer,
    j::Integer,
    z::AbstractFloat,
    domain::Domain,
    grid::Grid,
)
```

# Arguments

  - `i::Integer`: Grid index in x-direction
  - `j::Integer`: Grid index in y-direction
  - `z::AbstractFloat`: Target vertical coordinate
  - `domain::Domain`: Domain configuration
  - `grid::Grid`: Grid structure
"""
function get_next_level(
    i::Integer,
    j::Integer,
    z::AbstractFloat,
    domain::Domain,
    grid::Grid,
)
    (; sizezz, nzz, ko, k0, k1) = domain
    (; ztfc) = grid

    @views k = argmin(abs.(ztfc[i, j, :] .- z))
    if ztfc[i, j, k] < z
        k += 1
    end

    if ko == 0
        k = max(k, k0)
    else
        if k < 3
            error("Error in get_next_level: k = ", k, " < 3")
        end
    end

    if ko + nzz == sizezz
        k = min(k, k1 + 1)
    else
        if k > nzz - 1
            error(
                "Error in get_next_level: k = ",
                k,
                " > ",
                nzz - 1,
                " = nzz - 1",
            )
        end
    end

    return k
end
