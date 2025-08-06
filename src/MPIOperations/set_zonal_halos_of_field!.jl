"""
```julia
set_zonal_halos_of_field!(
    field::AbstractMatrix{<:AbstractFloat},
    namelists::Namelists,
    domain::Domain,
)
```

Exchange all zonal halo values of a matrix by performing bidirectional MPI communication between left and right neighbor processes.

```julia
set_zonal_halos_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Exchange a specified number of zonal halo values of a 3D array with an algorithm similar to that implemented in the method for matrices.

```julia
set_zonal_halos_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Exchange a specified number of zonal halo values of a 5D array with an algorithm similar to that implemented in the method for 3D arrays.

The first three dimensions of the array are assumed to represent the dimensions of physical space.

# Arguments

  - `field`: Input array.
  - `namelists`: Namelists with all model parameters.
  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

# Keywords

  - `layers`: The number of halo layers in each dimension. Use `-1` for the default values from `namelists`.
"""
function set_zonal_halos_of_field! end

function set_zonal_halos_of_field!(
    field::AbstractMatrix{<:AbstractFloat},
    namelists::Namelists,
    domain::Domain,
)
    (; nbx) = namelists.domain
    (; comm, i0, i1, left, right) = domain

    @views MPI.Sendrecv!(
        field[(i1 - nbx + 1):i1, :],
        field[(i0 - nbx):(i0 - 1), :],
        comm;
        dest = right,
        source = left,
    )

    @views MPI.Sendrecv!(
        field[i0:(i0 + nbx - 1), :],
        field[(i1 + 1):(i1 + nbx), :],
        comm;
        dest = left,
        source = right,
    )

    return
end

function set_zonal_halos_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; comm, i0, i1, j0, j1, k0, k1, left, right) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    j = (j0 - nby):(j1 + nby)
    k = (k0 - nbz):(k1 + nbz)

    @views MPI.Sendrecv!(
        field[(i1 - nbx + 1):i1, j, k],
        field[(i0 - nbx):(i0 - 1), j, k],
        comm;
        dest = right,
        source = left,
    )

    @views MPI.Sendrecv!(
        field[i0:(i0 + nbx - 1), j, k],
        field[(i1 + 1):(i1 + nbx), j, k],
        comm;
        dest = left,
        source = right,
    )

    return
end

function set_zonal_halos_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; comm, i0, i1, j0, j1, k0, k1, left, right) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    j = (j0 - nby):(j1 + nby)
    k = (k0 - nbz):(k1 + nbz)

    @views MPI.Sendrecv!(
        field[(i1 - nbx + 1):i1, j, k, :, :],
        field[(i0 - nbx):(i0 - 1), j, k, :, :],
        comm;
        dest = right,
        source = left,
    )

    @views MPI.Sendrecv!(
        field[i0:(i0 + nbx - 1), j, k, :, :],
        field[(i1 + 1):(i1 + nbx), j, k, :, :],
        comm;
        dest = left,
        source = right,
    )

    return
end
