"""
    compute_spectral_bounds(wavenumbers::AbstractVector{<:AbstractFloat})

# Arguments

  - `wavenumbers::AbstractVector{<:AbstractFloat}`: Vector of wavenumber values to analyze
"""
function compute_spectral_bounds(wavenumbers::AbstractVector{<:AbstractFloat})

    # Initialize minima and maxima.
    krminp = 0.0
    krmaxp = 0.0
    krminn = 0.0
    krmaxn = 0.0

    # Compute minima and maxima.
    for kr in wavenumbers
        if kr > 0
            if krminp == 0
                krminp = kr
            else
                krminp = min(krminp, kr)
            end
            krmaxp = max(krmaxp, kr)
        elseif kr < 0
            if krminn == 0
                krminn = -kr
            else
                krminn = min(krminn, -kr)
            end
            krmaxn = max(krmaxn, -kr)
        end
    end

    # Adjust bounds that are zero.
    if krminp == krmaxp == krminn == krmaxn == 0
        krminp = 1.0
        krmaxp = 2.0
        krminn = 1.0
        krmaxn = 2.0
    elseif krminp == krmaxp == 0
        krminp = krminn
        krmaxp = krmaxn
    elseif krminn == krmaxn == 0
        krminn = krminp
        krmaxn = krmaxp
    end

    # Prevent zero-width intervals.
    if krminn == krmaxn
        krminn /= 2
        krmaxn *= 2
    end
    if krminp == krmaxp
        krminp /= 2
        krmaxp *= 2
    end

    # Return.
    return (krminp, krmaxp, krminn, krmaxn)
end
