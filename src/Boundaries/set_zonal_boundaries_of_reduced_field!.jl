function set_zonal_boundaries_of_reduced_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain,
)
    (; npx) = namelists.domain
    (; i0, i1, j0, j1, k0, k1) = domain

    if npx > 1
        set_zonal_halos_of_reduced_field!(field, domain)
    else
        jy = (j0 - 1):(j1 + 1)
        kz = (k0 - 1):(k1 + 1)

        @views field[i0 - 1, jy, kz] .= field[i1, jy, kz]
        @views field[i1 + 1, jy, kz] .= field[i0, jy, kz]
    end

    return
end
