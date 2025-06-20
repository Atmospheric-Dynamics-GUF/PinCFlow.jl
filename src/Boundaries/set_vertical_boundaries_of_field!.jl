"""
    set_vertical_boundaries_of_field!(field, namelists, domain, zboundaries::SolidWallBoundaries, mode; layers, staggered)

Set vertical boundary conditions for 3D fields at solid walls.

# Arguments

  - `mode::Function`: Boundary condition mode (+ for symmetric, - for antisymmetric)
  - `staggered::Bool`: Whether field is on staggered vertical grid (sets boundary values to zero)
  - `layers::NTuple{3, <:Integer}`: Boundary layer sizes. Use -1 for defaults.
"""
function set_vertical_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries,
    mode::Function;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
    staggered = false,
)
    (; npz) = namelists.domain
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    if npz > 1
        set_vertical_halos_of_field!(
            field,
            namelists,
            domain,
            zboundaries;
            layers,
        )
    end

    i = (i0 - nbx):(i1 + nbx)
    j = (j0 - nby):(j1 + nby)

    if ko == 0
        if staggered
            field[i, j, k0 - 1] .= 0.0
            for k in 1:nbz
                @views field[i, j, k0 - k] .= mode(field[i, j, k0 + k - 2])
            end
        else
            for k in 1:nbz
                @views field[i, j, k0 - k] .= mode(field[i, j, k0 + k - 1])
            end
        end
    end

    if ko + nzz == sizezz
        if staggered
            field[i, j, k1] .= 0.0
            for k in 1:nbz
                @views field[i, j, k1 + k] .= mode(field[i, j, k1 - k])
            end
        else
            for k in 1:nbz
                @views field[i, j, k1 + k] .= mode(field[i, j, k1 - k + 1])
            end
        end
    end

    return
end

"""
    set_vertical_boundaries_of_field!(field::AbstractArray{<:AbstractFloat, 5}, namelists, domain, zboundaries; layers)

Set vertical boundary conditions for 5D fields. Uses halo exchange for multi-process domains.
"""
function set_vertical_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; npz) = namelists.domain

    if npz > 1
        set_vertical_halos_of_field!(
            field,
            namelists,
            domain,
            zboundaries;
            layers,
        )
    end

    return
end
