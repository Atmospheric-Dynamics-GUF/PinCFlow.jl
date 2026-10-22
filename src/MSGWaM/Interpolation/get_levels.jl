"""
```julia
get_levels(
    i::Integer,
    j::Integer,
    z::AbstractFloat,
    state::State;
    di::Integer = false,
    dj::Integer = false,
    dkd::Integer = 0,
    dku::Integer = 0,
)::NTuple{2, <:Integer}
```

Determine and return the indices of the levels that enclose `z` at the horizontal position ``\\left(i, j\\right)``, ``\\left(i + 1 / 2, j\\right)``, ``\\left(i, j + 1 / 2\\right)``, or ``\\left(i + 1 / 2, j + 1 / 2\\right)``, depending on `di` and `dj`.

This method is heavily used for interpolation to ray-volume positions. To ensure that the vertical boundary conditions are met and no out-of-bounds errors occur, the following constraints are set.

  - In MPI processes at the lower boundary of the domain, the computed index has the lower bound `state.domain.k0`. In other processes, an error is thrown if it is below `1 + dkd`.

  - In MPI processes at the upper boundary of the domain, the computed index has the upper bound `state.domain.k1 + 1`. In other processes, an error is thrown if it is above `state.domain.nzz - dku`.

In case an error is thrown, the parameter `wkb_cfl_number` of the discretization namelist should be set to a smaller value.

# Arguments

  - `i`: Zonal index.

  - `j`: Meridional index.

  - `z`: Vertical position.

  - `state`: Model state.

# Keywords

  - `di`: Switch for shifting the grid to ``i + 1 / 2``.

  - `dj`: Switch for shifting the grid to ``j + 1 / 2``.

  - `dkd`: Number of levels needed below.

  - `dku`: Number of levels needed above.

  - `old_grid`: Switch for using the old grid (`state.grid.zcold`).
"""
function get_levels end

@ivy function get_levels(
    i::Integer,
    j::Integer,
    z::AbstractFloat,
    state::State;
    di::Integer = false,
    dj::Integer = false,
    dkd::Integer = 0,
    dku::Integer = 0,
    old_grid::Bool = false,
)::NTuple{2, <:Integer}
    (; grid) = state
    (; z_size) = state.namelists.domain
    (; nz, nzz, ko, k0, k1) = state.domain

    zc = old_grid ? grid.zcold : grid.zc

    ku = 1
    while (
        zc[i, j, ku] +
        zc[i + di, j, ku] +
        zc[i, j + dj, ku] +
        zc[i + di, j + dj, ku]
    ) / 4 < z
        ku += 1
    end
    kd = ku - 1

    if ko == 0
        kd = max(kd, k0)
        ku = max(ku, k0)
    else
        if kd < 1 + dkd
            error(
                "Vertical index is too small: kd = ",
                kd,
                " < ",
                1 + dkd,
                " = 1 + dkd",
            )
        end
    end

    if ko + nz == z_size
        kd = min(kd, k1)
        ku = min(ku, k1)
    else
        if ku > nzz - dku
            error(
                "Vertical index is too large: ku = ",
                ku,
                " > ",
                nzz - dku,
                " = nzz - dku",
            )
        end
    end

    return (kd, ku)
end
