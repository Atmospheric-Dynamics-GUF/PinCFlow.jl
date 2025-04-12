function compute_spectral_bounds(wavenumbers::AbstractVector{<:AbstractFloat})

    # Compute minima and maxima.
    @views min_p = minimum(wavenumbers[wavenumbers .> 0.0]; init = 0.0)
    @views max_p = maximum(wavenumbers[wavenumbers .> 0.0]; init = 0.0)
    @views min_n = maximum(wavenumbers[wavenumbers .< 0.0]; init = 0.0)
    @views max_n = minimum(wavenumbers[wavenumbers .< 0.0]; init = 0.0)

    # Adjust bounds that are zero.
    if min_p == max_p == min_n == max_n == 0.0
        min_p = 1.0
        max_p = 2.0
        min_n = 1.0
        max_n = 2.0
    elseif min_p == max_p == 0.0
        min_p = min_n
        max_p = max_n
    elseif min_n == max_n == 0.0
        min_n = min_p
        max_n = max_p
    end

    # Prevent zero-width intervals.
    if min_n == max_n
        min_n = 0.5 * min_n
        max_n = 2.0 * max_n
    end
    if min_p == max_p
        min_p = 0.5 * min_p
        max_p = 2.0 * max_p
    end

    # Return.
    return (min_p, max_p, min_n, max_n)
end
