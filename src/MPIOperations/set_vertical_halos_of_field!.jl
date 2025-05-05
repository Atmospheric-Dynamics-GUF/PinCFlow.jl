function set_vertical_halos_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; nbz) = namelists.domain
    (; comm, sizezz, nzz, nx, ny, ko, i0, i1, j0, j1, k0, k1, down, up) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    (send_up, send_down, receive_up, receive_down) =
        (zeros(nx + 2 * nbx, ny + 2 * nby, nbz) for i in 1:4)

    i = (i0 - nbx):(i1 + nbx)
    j = (j0 - nby):(j1 + nby)

    if ko == 0
        for k in 1:nbz
            @views send_up[:, :, k] .= field[i, j, k1 - k + 1]
        end

        MPI.Sendrecv!(send_up, receive_up, comm; dest = up, source = up)

        for k in 1:nbz
            @views field[i, j, k1 + k] .= receive_up[:, :, k]
        end
    elseif ko + nzz == sizezz
        for k in 1:nbz
            @views send_down[:, :, k] .= field[i, j, k0 + k - 1]
        end

        MPI.Sendrecv!(send_down, receive_down, comm; dest = down, source = down)

        for k in 1:nbz
            @views field[i, j, k0 - k] .= receive_down[:, :, k]
        end
    else
        for k in 1:nbz
            @views send_up[:, :, k] .= field[i, j, k1 - k + 1]
            @views send_down[:, :, k] .= field[i, j, k0 + k - 1]
        end

        MPI.Sendrecv!(send_up, receive_down, comm; dest = up, source = down)

        MPI.Sendrecv!(send_down, receive_up, comm; dest = down, source = up)

        for k in 1:nbz
            @views field[i, j, k0 - k] .= receive_down[:, :, k]
            @views field[i, j, k1 + k] .= receive_up[:, :, k]
        end
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
    (; comm, sizezz, nzz, nx, ny, ko, i0, i1, j0, j1, k0, k1, down, up) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    (send_up, send_down, receive_up, receive_down) =
        (zeros(nx + 2 * nbx, ny + 2 * nby, nbz, 3, 2) for i in 1:4)

    i = (i0 - nbx):(i1 + nbx)
    j = (j0 - nby):(j1 + nby)

    if ko == 0
        for k in 1:nbz
            @views send_up[:, :, k, :, :] .= field[i, j, k1 - k + 1, :, :]
        end

        MPI.Sendrecv!(send_up, receive_up, comm; dest = up, source = up)

        for k in 1:nbz
            @views field[i, j, k1 + k, :, :] .= receive_up[:, :, k, :, :]
        end
    elseif ko + nzz == sizezz
        for k in 1:nbz
            @views send_down[:, :, k, :, :] .= field[i, j, k0 + k - 1, :, :]
        end

        MPI.Sendrecv!(send_down, receive_down, comm; dest = down, source = down)

        for k in 1:nbz
            @views field[i, j, k0 - k, :, :] .= receive_down[:, :, k, :, :]
        end
    else
        for k in 1:nbz
            @views send_up[:, :, k, :, :] .= field[i, j, k1 - k + 1, :, :]
            @views send_down[:, :, k, :, :] .= field[i, j, k0 + k - 1, :, :]
        end

        MPI.Sendrecv!(send_up, receive_down, comm; dest = up, source = down)

        MPI.Sendrecv!(send_down, receive_up, comm; dest = down, source = up)

        for k in 1:nbz
            @views field[i, j, k0 - k, :, :] .= receive_down[:, :, k, :, :]
            @views field[i, j, k1 + k, :, :] .= receive_up[:, :, k, :, :]
        end
    end

    return
end
