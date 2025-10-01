"""
```julia
get_next_level(
    i::Integer,
    j::Integer,
    z::AbstractFloat,
    state::State;
    dkd::Integer = 0,
    dku::Integer = 0,
)::Integer
```

Determine and return the index of the next level above `z` at the horizontal position ``\\left(i, j\\right)``.

This method is heavily used for interpolation to ray-volume positions. To ensure that the vertical boundary conditions are met and no out-of-bounds errors occur, the following constraints are set.

  - In MPI processes at the lower boundary of the domain, the computed index has the lower bound `state.domain.k0`. In other processes, an error is thrown if it is below `1 + dkd`.

  - In MPI processes at the upper boundary of the domain, the computed index has the upper bound `state.domain.k1 + 1`. In other processes, an error is thrown if it is above `state.domain.nzz - dku`.

In case an error is thrown, the parameter `cfl_wave` of the discretization namelist should be set to a smaller value.

# Arguments

  - `i`: Zonal index.

  - `j`: Meridional index.

  - `z`: Vertical position.

  - `state`: Model state.

# Keywords

  - `dkd`: Number of levels needed below.

  - `dku`: Number of levels needed above.
"""
function get_next_level end

function get_next_level(
    i::Integer,
    j::Integer,
    z::AbstractFloat,
    state::State;
    dkd::Integer = 0,
    dku::Integer = 0,
)::Integer
    (; sizezz, nzz, ko, k0, k1) = state.domain
    (; ztfc) = state.grid

    @ivy k = argmin(abs.(ztfc[i, j, :] .- z))
    @ivy if ztfc[i, j, k] < z
        k += 1
    end

    if ko == 0
        k = max(k, k0)
    else
        if k < 1 + dkd
            error(
                "Error in get_next_level: k = ",
                k,
                " < ",
                1 + dkd,
                " = 1 + dkd",
                "\nPlease choose a smaller WKB-CFL number.",
            )
        end
    end

    if ko + nzz == sizezz
        k = min(k, k1 + 1)
    else
        if k > nzz - dku
            error(
                "Error in get_next_level: k = ",
                k,
                " > ",
                nzz - dku,
                " = nzz - dku",
                "\nPlease choose a smaller WKB-CFL number.",
            )
        end
    end

    return k
end
