"""
    set_zonal_boundaries_of_field!(field::AbstractMatrix{<:AbstractFloat}, namelists::Namelists, domain::Domain)

Enforce zonal boundary conditions for 2D fields.

Halo exchange is used for multi-process domains (`npx > 1`), otherwise periodic boundaries are set by copying values from opposite domain edges.
"""
function set_zonal_boundaries_of_field!(
    field::AbstractMatrix{<:AbstractFloat},
    namelists::Namelists,
    domain::Domain,
)
    (; npx, nbx) = namelists.domain
    (; i0, i1) = domain

    if npx > 1
        set_zonal_halos_of_field!(field, namelists, domain)
    else
        for i in 1:nbx
            @views field[i0 - i, :] .= field[i1 - i + 1, :]
            @views field[i1 + i, :] .= field[i0 + i - 1, :]
        end
    end

    return
end

"""
    set_zonal_boundaries_of_field!(field::AbstractArray{<:Real, 3}, namelists::Namelists, domain::Domain; <keyword arguments>)

Enforce zonal boundary conditions for 3D fields.

Halo exchange is used for multi-process domains (`npx > 1`), otherwise periodic boundaries are set by copying values from opposite domain edges.

# Arguments

- `layers::NTuple{3, <:Integer} = (-1, -1, -1)`: The number of boundary layers in each dimension. Use `-1` for the default values from `namelists`.
"""
function set_zonal_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
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

"""
    set_zonal_boundaries_of_field!(field::AbstractArray{<:AbstractFloat, 5}, namelists::Namelists, domain::Domain; <keyword arguments>)

Enforce zonal boundary conditions for 5D fields.

Halo exchange is used for multi-process domains (`npx > 1`), otherwise periodic boundaries are set by copying values from opposite domain edges. Boundary conditions are enforced across all elements in dimensions 4 and 5.

# Arguments

- `layers::NTuple{3, <:Integer} = (-1, -1, -1)`: The number of boundary layers in each dimension. Use `-1` for the default values from `namelists`.
"""
function set_zonal_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
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
