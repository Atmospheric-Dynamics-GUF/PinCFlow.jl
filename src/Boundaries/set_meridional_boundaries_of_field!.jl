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

function set_meridional_boundaries_of_field!(
    field::AbstractMatrix{<:AbstractFloat},
    namelists::Namelists,
    domain::Domain,
)
    (; npy, nby) = namelists.domain
    (; j0, j1) = domain

    @ivy if npy > 1
        set_meridional_halos_of_field!(field, namelists, domain)
    else
        for j in 1:nby
            @. field[:, j0 - j] = field[:, j1 - j + 1]
            @. field[:, j1 + j] = field[:, j0 + j - 1]
        end
    end

    return
end

function set_meridional_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; npy) = namelists.domain
    (; i0, i1, j0, j1, k0, k1) = domain

    @ivy nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    @ivy nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    @ivy nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    @ivy if npy > 1
        set_meridional_halos_of_field!(field, namelists, domain; layers)
    else
        ii = (i0 - nbx):(i1 + nbx)
        kk = (k0 - nbz):(k1 + nbz)

        for j in 1:nby
            @. field[ii, j0 - j, kk] = field[ii, j1 - j + 1, kk]
            @. field[ii, j1 + j, kk] = field[ii, j0 + j - 1, kk]
        end
    end

    return
end

function set_meridional_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; npy) = namelists.domain
    (; i0, i1, j0, j1, k0, k1) = domain

    @ivy nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    @ivy nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    @ivy nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    @ivy if npy > 1
        set_meridional_halos_of_field!(field, namelists, domain; layers)
    else
        ii = (i0 - nbx):(i1 + nbx)
        kk = (k0 - nbz):(k1 + nbz)

        for j in 1:nby
            @. field[ii, j0 - j, kk, :, :] = field[ii, j1 - j + 1, kk, :, :]
            @. field[ii, j1 + j, kk, :, :] = field[ii, j0 + j - 1, kk, :, :]
        end
    end

    return
end
