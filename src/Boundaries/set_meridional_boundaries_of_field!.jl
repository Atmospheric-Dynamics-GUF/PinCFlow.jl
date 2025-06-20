"""
    set_meridional_boundaries_of_field!(field::AbstractMatrix, namelists, domain)

Set meridional boundary conditions for 2D fields. Uses halo exchange for multi-process domains (`npy > 1`),
otherwise applies periodic boundaries by copying values from opposite domain edges.

# Arguments

  - `field::AbstractMatrix{<:AbstractFloat}`: 2D field array to apply boundaries to
  - `namelists::Namelists`: Configuration containing domain parameters (`npy`, `nby`)
  - `domain::Domain`: Domain indices containing `j0`, `j1`
"""
function set_meridional_boundaries_of_field!(
    field::AbstractMatrix{<:AbstractFloat},
    namelists::Namelists,
    domain::Domain,
)
    (; npy, nby) = namelists.domain
    (; j0, j1) = domain

    if npy > 1
        set_meridional_halos_of_field!(field, namelists, domain)
    else
        for j in 1:nby
            @views field[:, j0 - j] .= field[:, j1 - j + 1]
            @views field[:, j1 + j] .= field[:, j0 + j - 1]
        end
    end

    return
end

"""
    set_meridional_boundaries_of_field!(field::AbstractArray{<:Real, 3}, namelists, domain; layers)

Set meridional boundary conditions for 3D fields.

# Arguments

  - `field::AbstractArray{<:Real, 3}`: 3D field array to apply boundaries to
  - `namelists::Namelists`: Configuration containing domain parameters (`npy`, `nbx`, `nby`, `nbz`)
  - `domain::Domain`: Domain indices containing `i0`, `i1`, `j0`, `j1`, `k0`, `k1`
  - `layers::NTuple{3, <:Integer}`: Tuple `(nbx, nby, nbz)` specifying boundary layer thickness. Use `-1` for default values from namelists.
"""
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

"""
    set_meridional_boundaries_of_field!(field::AbstractArray{<:AbstractFloat, 5}, namelists, domain; layers)

Set meridional boundary conditions for 5D fields. Applies boundary conditions across all elements
in dimensions 4 and 5.

# Arguments

  - `field::AbstractArray{<:AbstractFloat, 5}`: 5D field array to apply boundaries to
  - `namelists::Namelists`: Configuration containing domain parameters (`npy`, `nbx`, `nby`, `nbz`)
  - `domain::Domain`: Domain indices containing `i0`, `i1`, `j0`, `j1`, `k0`, `k1`
  - `layers::NTuple{3, <:Integer}`: Tuple `(nbx, nby, nbz)` specifying boundary layer thickness. Use `-1` for default values from namelists.
"""
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
