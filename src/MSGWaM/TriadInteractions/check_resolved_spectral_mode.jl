function check_resolved_spectral_mode end

function check_resolved_spectral_mode(
    spec_tend::TriadTendencies,
    kpvalue::AbstractFloat,
    mvalue::AbstractFloat,
    ::Triad2D,
)::Bool
    (; kpc, mc, ml) = spec_tend.spec_grid

    iz2 = ml ÷ 2 + 1

    return kpc[1] <= kpvalue <= kpc[end] &&
           ((mc[1] <= mvalue <= mc[iz2]) || (mc[iz2 + 1] <= mvalue <= mc[end]))

end