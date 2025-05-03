function set_vertical_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 3},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries,
    mode::Function;
    staggered = false,
)
    (; nbz, npz) = namelists.domain
    (; sizezz, nzz, ko, k0, k1) = domain

    if npz > 1
        set_vertical_halos_of_field!(field, namelists, domain, zboundaries)
    end

    if ko == 0
        if staggered
            field[:, :, k0 - 1] .= 0.0
            for k in 1:nbz
                @views field[:, :, k0 - k] .= mode(field[:, :, k0 + k - 2])
            end
        else
            for k in 1:nbz
                @views field[:, :, k0 - k] .= mode(field[:, :, k0 + k - 1])
            end
        end
    end

    if ko + nzz == sizezz
        if staggered
            field[:, :, k1] .= 0.0
            for k in 1:nbz
                @views field[:, :, k1 + k] .= mode(field[:, :, k1 - k])
            end
        else
            for k in 1:nbz
                @views field[:, :, k1 + k] .= mode(field[:, :, k1 - k + 1])
            end
        end
    end

    return
end

function set_vertical_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries,
)
    (; npz) = namelists.domain

    if npz > 1
        set_vertical_halos_of_field!(field, namelists, domain, zboundaries)
    end

    return
end
