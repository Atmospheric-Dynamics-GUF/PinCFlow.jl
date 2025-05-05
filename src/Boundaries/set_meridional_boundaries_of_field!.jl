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
