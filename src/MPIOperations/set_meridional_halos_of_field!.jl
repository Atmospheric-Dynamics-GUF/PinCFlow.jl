function set_meridional_halos_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; comm, nx, nz, i0, i1, j0, j1, k0, k1, back, forw) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    (send_forw, send_back, receive_forw, receive_back) =
        (zeros(nx + 2 * nbx, nby, nz + 2 * nbz) for i in 1:4)

    i = (i0 - nbx):(i1 + nbx)
    k = (k0 - nbz):(k1 + nbz)

    for j in 1:nby
        @views send_forw[:, j, :] .= field[i, j1 - j + 1, k]
        @views send_back[:, j, :] .= field[i, j0 + j - 1, k]
    end

    MPI.Sendrecv!(
        send_forw,
        receive_back,
        comm;
        dest = forw,
        source = back,
    )

    MPI.Sendrecv!(
        send_back,
        receive_forw,
        comm;
        dest = back,
        source = forw,
    )

    for j in 1:nby
        @views field[i, j0 - j, k] .= receive_back[:, j, :]
        @views field[i, j1 + j, k] .= receive_forw[:, j, :]
    end

    return
end

function set_meridional_halos_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; comm, nx, nz, i0, i1, j0, j1, k0, k1, back, forw) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    (send_forw, send_back, receive_forw, receive_back) =
        (zeros(nx + 2 * nbx, nby, nz + 2 * nbz, 3, 2) for i in 1:4)

    i = (i0 - nbx):(i1 + nbx)
    k = (k0 - nbz):(k1 + nbz)

    for j in 1:nby
        @views send_forw[:, j, :, :, :] .= field[i, j1 - j + 1, k, :, :]
        @views send_back[:, j, :, :, :] .= field[i, j0 + j - 1, k, :, :]
    end

    MPI.Sendrecv!(
        send_forw,
        receive_back,
        comm;
        dest = forw,
        source = back,
    )

    MPI.Sendrecv!(
        send_back,
        receive_forw,
        comm;
        dest = back,
        source = forw,
    )

    for j in 1:nby
        @views field[i, j0 - j, k, :, :] .= receive_back[:, j, :, :, :]
        @views field[i, j1 + j, k, :, :] .= receive_forw[:, j, :, :, :]
    end

    return
end
