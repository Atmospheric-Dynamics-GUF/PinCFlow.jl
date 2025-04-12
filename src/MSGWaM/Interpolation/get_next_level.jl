function get_next_level(
    i::Integer,
    j::Integer,
    z::AbstractFloat,
    domain::Domain,
    grid::Grid,
)

    # Get all necessary fields.
    (; nzz) = domain
    (; ztfc) = grid

    # Compute the vertical index.
    @views k = argmin(abs.(ztfc[i, j, :] .- z))
    if ztfc[i, j, k] < z
        k += 1
    end
    if k < 2
        k = 2
    end
    if k > nzz
        k = nzz
    end
    return k
end
