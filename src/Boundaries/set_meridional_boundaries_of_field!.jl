"""
```julia
set_meridional_boundaries_of_field!(
    field::AbstractMatrix{<:AbstractFloat},
    namelists::Namelists,
    domain::Domain,
)
```

Enforce meridional boundary conditions for a matrix.

Halo exchange is used for multi-process domains (`npy > 1`), otherwise periodic boundaries are set by copying values from opposite domain edges.

```julia
set_meridional_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Enforce meridional boundary conditions for a 3D array.

Halo exchange is used in the same manner as in the method for matrices.

```julia
set_meridional_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Enforce meridional boundary conditions for a 5D array.

Halo exchange is used in the same manner as in the methods for matrices and 3D arrays. The first three dimensions of the array are assumed to represent the dimensions of physical space.

# Arguments

  - `field`: Input array.

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

# Keywords

  - `layers`: The number of boundary layers in each dimension. Use `-1` for the default values from `namelists`.

# See also

  - [`PinCFlow.MPIOperations.set_meridional_halos_of_field!`](@ref)
"""
function set_meridional_boundaries_of_field! end

@ivy function set_meridional_boundaries_of_field!(
    field::AbstractMatrix{<:AbstractFloat},
    namelists::Namelists,
    domain::Domain,
)
    (; y_size, nby) = namelists.domain
    (; j0, j1, nxx) = domain

    if y_size > 1
        set_meridional_halos_of_field!(field, namelists, domain)
    else
        for j in 1:nby
            @share for i in 1:nxx
                field[i, j0 - j] = field[i, j1 - j + 1]
                field[i, j1 + j] = field[i, j0 + j - 1]
            end
        end
    end

    return
end

@ivy function set_meridional_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; y_size) = namelists.domain
    (; i0, i1, j0, j1, k0, k1) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    if y_size > 1
        set_meridional_halos_of_field!(field, namelists, domain; layers)
    else
        @share vector = false for k in (k0 - nbz):(k1 + nbz)
            for j in 1:nby
                @share thread = false for i in (i0 - nbx):(i1 + nbx)
                    field[i, j0 - j, k] = field[i, j1 - j + 1, k]
                    field[i, j1 + j, k] = field[i, j0 + j - 1, k]
                end
            end
        end
    end

    return
end

@ivy function set_meridional_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; y_size) = namelists.domain
    (; i0, i1, j0, j1, k0, k1) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    if y_size > 1
        set_meridional_halos_of_field!(field, namelists, domain; layers)
    else
        @share vector = false for m in 1:2, l in 1:3, k in (k0 - nbz):(k1 + nbz)
            for j in 1:nby
                @share thread = false for i in (i0 - nbx):(i1 + nbx)
                    field[i, j0 - j, k, l, m] = field[i, j1 - j + 1, k, l, m]
                    field[i, j1 + j, k, l, m] = field[i, j0 + j - 1, k, l, m]
                end
            end
        end
    end

    return
end
