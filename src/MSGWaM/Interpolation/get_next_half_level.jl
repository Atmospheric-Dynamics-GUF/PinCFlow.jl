"""
```julia
get_next_half_level(
    i::Integer,
    j::Integer,
    z::AbstractFloat,
    state::State,
)::Integer
```

Determine and return the index of the next half-level above `z` at the horizontal position `(i, j)`.

This method is heavily used for interpolation to ray-volume positions. To ensure that the vertical boundary conditions are met and no out-of-bounds errors occur, the following constraints are set.

  - In MPI processes at the lower boundary of the domain, the returned index cannot be smaller than `state.domain.k0`, in other processes, it cannot be smaller than 3.

  - In MPI processes at the upper boundary of the domain, the returned index cannot be larger than `state.domain.k1`, in other processes, it cannot be larger than `state.domain.nzz - 1`.

# Arguments

  - `i`: Zonal index.

  - `j`: Meridional index.

  - `z`: Vertical position.

  - `state`: Model state.
"""
function get_next_half_level end

function get_next_half_level(
    i::Integer,
    j::Integer,
    z::AbstractFloat,
    state::State,
)::Integer
    (; sizezz, nzz, ko, k0, k1) = state.domain
    (; ztildetfc) = state.grid

    @ivy k = argmin(@. abs(ztildetfc[i, j, :] - z))
    @ivy if ztildetfc[i, j, k] < z
        k += 1
    end

    if ko == 0
        k = max(k, k0)
    else
        if k < 3
            error("Error in get_next_half_level: k = ", k, " < 3")
        end
    end

    if ko + nzz == sizezz
        k = min(k, k1)
    else
        if k > nzz - 1
            error(
                "Error in get_next_half_level: k = ",
                k,
                " > ",
                nzz - 1,
                " = nzz - 1",
            )
        end
    end

    return k
end
