"""
    set_zonal_boundaries_of_field!(field::AbstractMatrix, namelists, domain)

Set zonal boundary conditions for 2D fields. Uses halo exchange for multi-process domains
(`npx > 1`), otherwise applies periodic boundaries.
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
    set_zonal_boundaries_of_field!(field::AbstractArray{<:Real, 3}, namelists, domain; layers)

Set zonal boundary conditions for 3D fields.

# Arguments

  - `layers::NTuple{3, <:Integer}`: Boundary layer sizes (nbx, nby, nbz). Use -1 for defaults.
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
    set_zonal_boundaries_of_field!(field::AbstractArray{<:AbstractFloat, 5}, namelists, domain; layers)

Set zonal boundary conditions for 5D fields. Applies boundaries to all elements in
dimensions 4 and 5.
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
