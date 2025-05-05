function set_zonal_halos_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; comm, ny, nz, i0, i1, j0, j1, k0, k1, left, right) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    (send_right, send_left, receive_right, receive_left) =
        (zeros(nbx, ny + 2 * nby, nz + 2 * nbz) for i in 1:4)

    j = (j0 - nby):(j1 + nby)
    k = (k0 - nbz):(k1 + nbz)

    for i in 1:nbx
        @views send_right[i, :, :] .= field[i1 - i + 1, j, k]
        @views send_left[i, :, :] .= field[i0 + i - 1, j, k]
    end

    MPI.Sendrecv!(
        send_right,
        receive_left,
        comm;
        dest = right,
        source = left,
    )

    MPI.Sendrecv!(
        send_left,
        receive_right,
        comm;
        dest = left,
        source = right,
    )

    for i in 1:nbx
        @views field[i0 - i, j, k] .= receive_left[i, :, :]
        @views field[i1 + i, j, k] .= receive_right[i, :, :]
    end

    return
end

function set_zonal_halos_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; comm, ny, nz, i0, i1, j0, j1, k0, k1, left, right) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    (send_right, send_left, receive_right, receive_left) =
        (zeros(nbx, ny + 2 * nby, nz + 2 * nbz, 3, 2) for i in 1:4)

    j = (j0 - nby):(j1 + nby)
    k = (k0 - nbz):(k1 + nbz)

    for i in 1:nbx
        @views send_right[i, :, :, :, :] .= field[i1 - i + 1, j, k, :, :]
        @views send_left[i, :, :, :, :] .= field[i0 + i - 1, j, k, :, :]
    end

    MPI.Sendrecv!(
        send_right,
        receive_left,
        comm;
        dest = right,
        source = left,
    )

    MPI.Sendrecv!(
        send_left,
        receive_right,
        comm;
        dest = left,
        source = right,
    )

    for i in 1:nbx
        @views field[i0 - i, j, k, :, :] .= receive_left[i, :, :, :, :]
        @views field[i1 + i, j, k, :, :] .= receive_right[i, :, :, :, :]
    end

    return
end
