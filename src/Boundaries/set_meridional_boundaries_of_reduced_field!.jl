function set_meridional_boundaries_of_reduced_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain,
)
    (; npy) = namelists.domain
    (; i0, i1, j0, j1, k0, k1) = domain

    if npy > 1
        set_meridional_halos_of_reduced_field!(field, domain)
    else
        ix = (i0 - 1):(i1 + 1)
        kz = (k0 - 1):(k1 + 1)

        @views field[ix, j0 - 1, kz] .= field[ix, j1, kz]
        @views field[ix, j1 + 1, kz] .= field[ix, j0, kz]
    end

    return
end
