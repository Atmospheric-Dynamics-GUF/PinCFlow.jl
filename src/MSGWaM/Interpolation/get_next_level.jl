function get_next_level(
    i::Integer,
    j::Integer,
    z::AbstractFloat,
    domain::Domain,
    grid::Grid,
)
    (; sizezz, nzz, nz, ko, k0, k1) = domain
    (; ztfc) = grid

    @views k = argmin(abs.(ztfc[i, j, :] .- z))
    if ztfc[i, j, k] < z && k < nzz
        k += 1
    end
    if ko + k < k0
        k = k0
    end
    if ko + k > sizezz - (nzz - nz) / 2 + 1
        k = k1 + 1
    end

    return k
end
