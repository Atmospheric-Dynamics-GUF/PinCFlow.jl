function set_vertical_boundaries_of_reduced_field!(
    field::AbstractArray{<:AbstractFloat, 3},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries,
    mode::Function;
)
    (; npz) = namelists.domain
    (; sizezz, nzz, ko, i0, i1, j0, j1, k0, k1) = domain

    if npz > 1
        set_vertical_halos_of_reduced_field!(
            field,
            namelists,
            domain,
            zboundaries,
        )
    end

    ix = (i0 - 1):(i1 + 1)
    jy = (j0 - 1):(j1 + 1)

    if ko == 0
        @views field[ix, jy, k0 - 1] .= mode(field[ix, jy, k0])
    end

    if ko + nzz == sizezz
        @views field[ix, jy, k1 + 1] .= mode(field[ix, jy, k1])
    end

    return
end
