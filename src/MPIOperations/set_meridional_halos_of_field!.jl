"""
```julia
set_meridional_halos_of_field!(
    field::AbstractMatrix{<:AbstractFloat},
    namelists::Namelists,
    domain::Domain,
)
```

Exchange meridional halo values of a matrix by performing bidirectional MPI communication between backward and forward neighbor processes.

```julia
set_meridional_halos_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Exchange meridional halo values of a 3D array by performing bidirectional MPI communication between backward and forward neighbor processes.

```julia
set_meridional_halos_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Exchange meridional halo values of a 5D array by performing bidirectional MPI communication between backward and forward neighbor processes.

The first three dimensions of the array are assumed to represent the dimensions of physical space.

# Arguments

  - `field`: Input array.
  - `namelists`: Namelists with all model parameters.
  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

# Keywords

  - `layers`: The number of halo layers in each dimension. Use `-1` for the default values from `namelists`.
"""
function set_meridional_halos_of_field! end

function set_meridional_halos_of_field!(
    field::AbstractMatrix{<:AbstractFloat},
    namelists::Namelists,
    domain::Domain,
)
    (; nby) = namelists.domain
    (; comm, j0, j1, backward, forward) = domain

    @views MPI.Sendrecv!(
        field[:, (j1 - nby + 1):j1],
        field[:, (j0 - nby):(j0 - 1)],
        comm;
        dest = forward,
        source = backward,
    )

    @views MPI.Sendrecv!(
        field[:, j0:(j0 + nby - 1)],
        field[:, (j1 + 1):(j1 + nby)],
        comm;
        dest = backward,
        source = forward,
    )

    return
end

function set_meridional_halos_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; comm, i0, i1, j0, j1, k0, k1, backward, forward) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    i = (i0 - nbx):(i1 + nbx)
    k = (k0 - nbz):(k1 + nbz)

    @views MPI.Sendrecv!(
        field[i, (j1 - nby + 1):j1, k],
        field[i, (j0 - nby):(j0 - 1), k],
        comm;
        dest = forward,
        source = backward,
    )

    @views MPI.Sendrecv!(
        field[i, j0:(j0 + nby - 1), k],
        field[i, (j1 + 1):(j1 + nby), k],
        comm;
        dest = backward,
        source = forward,
    )

    return
end

function set_meridional_halos_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; comm, i0, i1, j0, j1, k0, k1, backward, forward) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    i = (i0 - nbx):(i1 + nbx)
    k = (k0 - nbz):(k1 + nbz)

    @views MPI.Sendrecv!(
        field[i, (j1 - nby + 1):j1, k, :, :],
        field[i, (j0 - nby):(j0 - 1), k, :, :],
        comm;
        dest = forward,
        source = backward,
    )

    @views MPI.Sendrecv!(
        field[i, j0:(j0 + nby - 1), k, :, :],
        field[i, (j1 + 1):(j1 + nby), k, :, :],
        comm;
        dest = backward,
        source = forward,
    )

    return
end
