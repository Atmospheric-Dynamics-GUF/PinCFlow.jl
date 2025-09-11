"""
```julia
set_vertical_halos_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Exchange a specified number of vertical halo values of a 3D array by performing MPI communication between downward and upward neighbor processes.

Solid walls are assumed at the vertical boundaries of the domain. The corresponding ghost-cell values are not changed.

```julia
set_vertical_halos_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Exchange a specified number of vertical halo values of a 5D array with an algorithm similar to that implemented in the above method.

The vertical domain boundaries are treated as described above. The first three dimensions of the array are assumed to represent the dimensions of physical space.

# Arguments

  - `field`: Input array.

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

# Keywords

  - `layers`: The number of halo layers in each dimension. Use `-1` for the default values from `namelists`.
"""
function set_vertical_halos_of_field! end

function set_vertical_halos_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; nbz) = namelists.domain
    (; comm, sizezz, nzz, ko, i0, i1, j0, j1, k0, k1, down, up) = domain

    @ivy nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    @ivy nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    @ivy nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    ii = (i0 - nbx):(i1 + nbx)
    jj = (j0 - nby):(j1 + nby)

    @ivy if ko == 0
        MPI.Sendrecv!(
            field[ii, jj, (k1 - nbz + 1):k1],
            field[ii, jj, (k1 + 1):(k1 + nbz)],
            comm;
            dest = up,
            source = up,
        )
    elseif ko + nzz == sizezz
        MPI.Sendrecv!(
            field[ii, jj, k0:(k0 + nbz - 1)],
            field[ii, jj, (k0 - nbz):(k0 - 1)],
            comm;
            dest = down,
            source = down,
        )
    else
        MPI.Sendrecv!(
            field[ii, jj, (k1 - nbz + 1):k1],
            field[ii, jj, (k0 - nbz):(k0 - 1)],
            comm;
            dest = up,
            source = down,
        )

        MPI.Sendrecv!(
            field[ii, jj, k0:(k0 + nbz - 1)],
            field[ii, jj, (k1 + 1):(k1 + nbz)],
            comm;
            dest = down,
            source = up,
        )
    end

    return
end

function set_vertical_halos_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; nbz) = namelists.domain
    (; comm, sizezz, nzz, ko, i0, i1, j0, j1, k0, k1, down, up) = domain

    @ivy nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    @ivy nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    @ivy nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    ii = (i0 - nbx):(i1 + nbx)
    jj = (j0 - nby):(j1 + nby)

    @ivy if ko == 0
        MPI.Sendrecv!(
            field[ii, jj, (k1 - nbz + 1):k1, :, :],
            field[ii, jj, (k1 + 1):(k1 + nbz), :, :],
            comm;
            dest = up,
            source = up,
        )
    elseif ko + nzz == sizezz
        MPI.Sendrecv!(
            field[ii, jj, k0:(k0 + nbz - 1), :, :],
            field[ii, jj, (k0 - nbz):(k0 - 1), :, :],
            comm;
            dest = down,
            source = down,
        )
    else
        MPI.Sendrecv!(
            field[ii, jj, (k1 - nbz + 1):k1, :, :],
            field[ii, jj, (k0 - nbz):(k0 - 1), :, :],
            comm;
            dest = up,
            source = down,
        )

        MPI.Sendrecv!(
            field[ii, jj, k0:(k0 + nbz - 1), :, :],
            field[ii, jj, (k1 + 1):(k1 + nbz), :, :],
            comm;
            dest = down,
            source = up,
        )
    end

    return
end
