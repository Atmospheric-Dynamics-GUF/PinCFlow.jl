function set_vertical_halos_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; nbz) = namelists.domain
    (; comm, sizezz, nzz, ko, i0, i1, j0, j1, k0, k1, down, up) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    i = (i0 - nbx):(i1 + nbx)
    j = (j0 - nby):(j1 + nby)

    if ko == 0
        @views MPI.Sendrecv!(
            field[i, j, (k1 - nbz + 1):k1],
            field[i, j, (k1 + 1):(k1 + nbz)],
            comm;
            dest = up,
            source = up,
        )
    elseif ko + nzz == sizezz
        @views MPI.Sendrecv!(
            field[i, j, k0:(k0 + nbz - 1)],
            field[i, j, (k0 - nbz):(k0 - 1)],
            comm;
            dest = down,
            source = down,
        )
    else
        @views MPI.Sendrecv!(
            field[i, j, (k1 - nbz + 1):k1],
            field[i, j, (k0 - nbz):(k0 - 1)],
            comm;
            dest = up,
            source = down,
        )

        @views MPI.Sendrecv!(
            field[i, j, k0:(k0 + nbz - 1)],
            field[i, j, (k1 + 1):(k1 + nbz)],
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
    domain::Domain,
    zboundaries::SolidWallBoundaries;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; nbz) = namelists.domain
    (; comm, sizezz, nzz, ko, i0, i1, j0, j1, k0, k1, down, up) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    i = (i0 - nbx):(i1 + nbx)
    j = (j0 - nby):(j1 + nby)

    if ko == 0
        @views MPI.Sendrecv!(
            field[i, j, (k1 - nbz + 1):k1, :, :],
            field[i, j, (k1 + 1):(k1 + nbz), :, :],
            comm;
            dest = up,
            source = up,
        )
    elseif ko + nzz == sizezz
        @views MPI.Sendrecv!(
            field[i, j, k0:(k0 + nbz - 1), :, :],
            field[i, j, (k0 - nbz):(k0 - 1), :, :],
            comm;
            dest = down,
            source = down,
        )
    else
        @views MPI.Sendrecv!(
            field[i, j, (k1 - nbz + 1):k1, :, :],
            field[i, j, (k0 - nbz):(k0 - 1), :, :],
            comm;
            dest = up,
            source = down,
        )

        @views MPI.Sendrecv!(
            field[i, j, k0:(k0 + nbz - 1), :, :],
            field[i, j, (k1 + 1):(k1 + nbz), :, :],
            comm;
            dest = down,
            source = up,
        )
    end

    return
end
