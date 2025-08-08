"""
```julia
set_vertical_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries,
    mode::Function;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
    staggered = false,
)
```

Enforce vertical boundary conditions for a 3D array (assuming solid-wall boundaries).

Halo exchange is used for multi-process domains (`npz > 1`). Use `mode = +` (`mode = -`) for line-reflected (point-reflected) ghost-cell values.

```julia
set_vertical_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Exchange halo values of a 5D array if multiple processes are used in the vertical (`npz > 1`).

This method is applied to reconstruction arrays. Vertical boundary conditions are not enforced for these but for the fluxes determined from them.

# Arguments

- `field`: Input array.
- `namelists`: Namelists with all model parameters.
- `domain`: Collection of domain-decomposition and MPI-communication parameters.
- `zboundaries`: Vertical boundary conditions.
- `mode`: Method used for setting the boundary-cell values.

# Keywords

- `layers`: The number of boundary layers in each dimension. Use `-1` for the default values from `namelists`.
- `staggered`: A switch for whether or not the field is on the staggered vertical grid.

# See also

- [`PinCFlow.MPIOperations.set_vertical_halos_of_field!`](@ref)
"""
function set_vertical_boundaries_of_field! end

function set_vertical_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries,
    mode::Function;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
    staggered = false,
)
    (; npz) = namelists.domain
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    if npz > 1
        set_vertical_halos_of_field!(
            field,
            namelists,
            domain,
            zboundaries;
            layers,
        )
    end

    i = (i0 - nbx):(i1 + nbx)
    j = (j0 - nby):(j1 + nby)

    if ko == 0
        if staggered
            field[i, j, k0 - 1] .= 0.0
            for k in 1:nbz
                @views field[i, j, k0 - k] .= mode(field[i, j, k0 + k - 2])
            end
        else
            for k in 1:nbz
                @views field[i, j, k0 - k] .= mode(field[i, j, k0 + k - 1])
            end
        end
    end

    if ko + nzz == sizezz
        if staggered
            field[i, j, k1] .= 0.0
            for k in 1:nbz
                @views field[i, j, k1 + k] .= mode(field[i, j, k1 - k])
            end
        else
            for k in 1:nbz
                @views field[i, j, k1 + k] .= mode(field[i, j, k1 - k + 1])
            end
        end
    end

    return
end

function set_vertical_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; npz) = namelists.domain

    if npz > 1
        set_vertical_halos_of_field!(
            field,
            namelists,
            domain,
            zboundaries;
            layers,
        )
    end

    return
end
