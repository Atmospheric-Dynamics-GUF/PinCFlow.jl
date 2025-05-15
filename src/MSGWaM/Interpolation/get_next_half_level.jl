function get_next_half_level(
    i::Integer,
    j::Integer,
    z::AbstractFloat,
    domain::Domain,
    grid::Grid,
)
    (; sizezz, nzz, nz, ko, k0, k1) = domain
    (; ztildetfc) = grid

    @views k = argmin(abs.(ztildetfc[i, j, :] .- z))
    if ztildetfc[i, j, k] < z && k < nzz
        k += 1
    end
    if ko + k < k0
        k = k0
    end
    if ko + k > sizezz - (nzz - nz) / 2
        k = k1
    end

    return k
end
