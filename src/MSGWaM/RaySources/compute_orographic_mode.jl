function compute_orographic_mode(
    displm::AbstractFloat,
    wnk::AbstractFloat,
    wnl::AbstractFloat,
    uavg::AbstractFloat,
    vavg::AbstractFloat,
    rhoavg::AbstractFloat,
    bvsavg::AbstractFloat,
    f_cor_nd::AbstractFloat,
    branch::Integer,
)

    # Compute horizontal wavenumber.
    wnh = sqrt(wnk^2 + wnl^2)

    # Compute intrinsic frequency from orographic wavenumbers.
    omi = -uavg * wnk - vavg * wnl

    # Adjust the signs to be consistent with the chosen frequency
    # branch.
    if omi * branch < 0
        omi = -omi
        wnk = -wnk
        wnl = -wnl
    end

    # Compute vertical wavenumber and wave-action density.
    if omi^2 > f_cor_nd^2 && omi^2 < bvsavg

        # Compute vertical wavenumber.
        wnm = -branch * sqrt(wnh^2 * (bvsavg - omi^2) / (omi^2 - f_cor_nd^2))

        # Compute wave-action density.
        wad = 0.5 * rhoavg * displm^2 * omi * (wnh^2 + wnm^2) / wnh^2

        # Set to zero if something went wrong.
        if wad != wad || wnm != wnm
            wad = 0.0
            wnm = 0.0
        end

        # Account for critical and reflecting levels.
    else
        wad = 0.0
        wnm = 0.0
    end

    # Return.
    return (omi, wnk, wnl, wnm, wad)
end
