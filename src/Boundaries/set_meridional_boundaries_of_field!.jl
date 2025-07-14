"""
```julia
set_meridional_boundaries_of_field!(
    field::AbstractMatrix{<:AbstractFloat},
    namelists::Namelists,
    domain::Domain,
)
```

Enforce meridional boundary conditions for 2D fields.

Halo exchange is used for multi-process domains (`npy > 1`), otherwise periodic boundaries are set by copying values from opposite domain edges.

# Arguments

  - `field`: 2D input field.
  - `namelists`: Namelists with all model parameters.
  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

# See also

  - [`PinCFlow.MPIOperations.set_meridional_halos_of_field!`](@ref)
"""
function set_meridional_boundaries_of_field!(
    field::AbstractMatrix{<:AbstractFloat},
    namelists::Namelists,
    domain::Domain,
)
    (; npy, nby) = namelists.domain
    (; j0, j1) = domain

    if npy > 1
        set_meridional_halos_of_field!(field, namelists, domain)
    else
        for j in 1:nby
            @views field[:, j0 - j] .= field[:, j1 - j + 1]
            @views field[:, j1 + j] .= field[:, j0 + j - 1]
        end
    end

    return
end

"""
```julia
set_meridional_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Enforce meridional boundary conditions for 3D fields.

Halo exchange is used for multi-process domains (`npy > 1`), otherwise periodic boundaries are set by copying values from opposite domain edges.

# Arguments

  - `field`: 3D input field.
  - `namelists`: Namelists with all model parameters.
  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
  - `layers`: The number of boundary layers in each dimension. Use `-1` for the default values from `namelists`.

# See also

  - [`PinCFlow.MPIOperations.set_meridional_halos_of_field!`](@ref)
"""
function set_meridional_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; npy) = namelists.domain
    (; i0, i1, j0, j1, k0, k1) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    if npy > 1
        set_meridional_halos_of_field!(field, namelists, domain; layers)
    else
        i = (i0 - nbx):(i1 + nbx)
        k = (k0 - nbz):(k1 + nbz)

        for j in 1:nby
            @views field[i, j0 - j, k] .= field[i, j1 - j + 1, k]
            @views field[i, j1 + j, k] .= field[i, j0 + j - 1, k]
        end
    end

    return
end

"""
```julia
set_meridional_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Enforce meridional boundary conditions for 5D fields.

Halo exchange is used for multi-process domains (`npy > 1`), otherwise periodic boundaries are set by copying values from opposite domain edges. Boundary conditions are enforced across all elements in dimensions 4 and 5.

# Arguments

  - `field`: 5D input field.
  - `namelists`: Namelists with all model parameters.
  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
  - `layers`: The number of boundary layers in each dimension. Use `-1` for the default values from `namelists`.

# See also

  - [`PinCFlow.MPIOperations.set_meridional_halos_of_field!`](@ref)
"""
function set_meridional_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; npy) = namelists.domain
    (; i0, i1, j0, j1, k0, k1) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    if npy > 1
        set_meridional_halos_of_field!(field, namelists, domain; layers)
    else
        i = (i0 - nbx):(i1 + nbx)
        k = (k0 - nbz):(k1 + nbz)

        for j in 1:nby
            @views field[i, j0 - j, k, :, :] .= field[i, j1 - j + 1, k, :, :]
            @views field[i, j1 + j, k, :, :] .= field[i, j0 + j - 1, k, :, :]
        end
    end

    return
end
