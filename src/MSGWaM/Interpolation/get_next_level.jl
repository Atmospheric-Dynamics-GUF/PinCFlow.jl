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

Determine the index of the next level above `z` at the horizontal position `(i, j)`.

This method is heavily used for interpolation to ray-volume positions. To ensure that the vertical boundary conditions are met and no out-of-bounds errors occur, the following constraints are set.

  - In MPI processes at the lower boundary of the domain, the returned index cannot be smaller than `domain.k0`, in other processes, it cannot be smaller than 3.
  - In MPI processes at the upper boundary of the domain, the returned index cannot be larger than `domain.k1 + 1`, in other processes, it cannot be larger than `domain.nzz - 1`.

# Arguments

  - `i`: Zonal index.
  - `j`: Meridional index.
  - `z`: Vertical position.
  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
  - `grid`: Collection of parameters and fields that describe the grid.

# Returns

  - `::Integer`: Index of the next level.
"""
function get_next_level end

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
