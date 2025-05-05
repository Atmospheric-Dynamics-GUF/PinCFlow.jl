function set_zonal_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1)
)
    (; npx) = namelists.domain
    (; i0, i1, j0, j1, k0, k1) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    if npx > 1
        set_zonal_halos_of_field!(field, namelists, domain; layers)
    else
        j = (j0 - nby):(j1 + nby)
        k = (k0 - nbz):(k1 + nbz)

        for i in 1:nbx
            @views field[i0 - i, j, k] .= field[i1 - i + 1, j, k]
            @views field[i1 + i, j, k] .= field[i0 + i - 1, j, k]
        end
    end

    return
end

function set_zonal_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1)
)
    (; npx) = namelists.domain
    (; i0, i1, j0, j1, k0, k1) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    if npx > 1
        set_zonal_halos_of_field!(field, namelists, domain; layers)
    else
        j = (j0 - nby):(j1 + nby)
        k = (k0 - nbz):(k1 + nbz)

        for i in 1:nbx
            @views field[i0 - i, j, k, :, :] .= field[i1 - i + 1, j, k, :, :]
            @views field[i1 + i, j, k, :, :] .= field[i0 + i - 1, j, k, :, :]
        end
    end

    return
end
