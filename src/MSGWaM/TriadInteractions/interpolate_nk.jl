function interpolate_nk end

function interpolate_nk(
    spec_tend::TriadTendencies,
    kpvalue::AbstractFloat,
    mvalue::AbstractFloat,
    triad_mode::Triad2D,
)::AbstractFloat

    (; kp, m, kpc, mc, kpl, ml, loglkp, loglm) = spec_tend.spec_grid
    (; c_o, alphakp, alpham, beta) = spec_tend.interp_coef

    if m[1] > 0
        return interpolate_nk(spec_tend, kpvalue, mvalue, Triad3DIso())
    end

    mhalf = ml ÷ 2
    iz2 = mhalf + 1

    # Positive m branch
    if kpc[1] <= kpvalue <= kpc[end] && mc[iz2 + 1] <= mvalue <= mc[end]

        kpi = ceil(Int, 1 + log(kpvalue / kp[1]) / loglkp)
        kpi = clamp(kpi, 1, kpl)

        mi = ceil(Int, 1 + log(mvalue / m[iz2]) / loglm)
        mi = clamp(mi, 1, mhalf)
        mi += mhalf

        nkvalue = c_o[kpi, mi] -
                  alphakp[kpi, mi] * kpvalue -
                  alpham[kpi, mi] * mvalue -
                  beta[kpi, mi] * kpvalue * mvalue

    # Negative m branch
    elseif kpc[1] <= kpvalue <= kpc[end] && mc[1] <= mvalue <= mc[iz2]

        kpi = ceil(Int, 1 + log(kpvalue / kp[1]) / loglkp)
        kpi = clamp(kpi, 1, kpl)

        mabs = abs(mvalue)

        mi = ceil(Int, 1 + log(mabs / m[iz2]) / loglm)
        mi = clamp(mi, 1, mhalf)
        mi = mhalf - mi + 1

        nkvalue = c_o[kpi, mi] -
                  alphakp[kpi, mi] * kpvalue -
                  alpham[kpi, mi] * mabs -
                  beta[kpi, mi] * kpvalue * mabs

    else
        nkvalue = 0.0
    end

    return max(nkvalue, 0.0)
end

function interpolate_nk(
    spec_tend::TriadTendencies,
    kpi::Integer,
    mvalue::AbstractFloat,
    ::Triad2D,
    )::AbstractFloat

    (; m, mc, ml, loglm) = spec_tend.spec_grid
    (; c_o, alpham) = spec_tend.interp_coef

    if m[1] > 0.0
        error("Discrete-k interpolation is currently implemented for signed-m Triad2D.")
    end

    mhalf = ml ÷ 2
    iz2 = mhalf + 1

    #----------------------------------------------------------
    # Check whether mvalue lies inside the represented
    # signed-m spectral domain.
    #----------------------------------------------------------

    if mvalue > 0.0
        if !(mc[iz2 + 1] <= mvalue <= mc[end])
            return 0.0
        end
    elseif mvalue < 0.0
        if !(mc[1] <= mvalue <= mc[iz2])
            return 0.0
        end
    else
        return 0.0
    end

    #----------------------------------------------------------
    # Locate the interpolation interval on the logarithmic
    # |m| grid.
    #----------------------------------------------------------

    mabs = abs(mvalue)

    mj = ceil(Int, 1 + log(mabs / m[iz2]) / loglm)
    mj = clamp(mj, 1, mhalf)

    if mvalue > 0.0
        mi = mhalf + mj
    else
        mi = mhalf - mj + 1
    end

    #----------------------------------------------------------
    # Linear interpolation/extrapolation in physical |m|.
    #
    # The coefficients were computed beforehand from the
    # current wave spectrum.
    #----------------------------------------------------------

    nkvalue = c_o[kpi, mi] - alpham[kpi, mi] * mabs

    return max(nkvalue, 0.0)
end

function interpolate_nk(spec_tend::TriadTendencies,
    kpvalue::AbstractFloat,
    mvalue::AbstractFloat,
    triad_mode::Triad3DIso)::AbstractFloat

    (; kp, m, kpc, mc, kpl, ml, loglkp, loglm) = spec_tend.spec_grid
    (; c_o, alphakp, alpham, beta) = spec_tend.interp_coef

    @ivy if kpc[1] <= kpvalue <= kpc[end] && mc[1] <= mvalue <= mc[end]
        kpi = ceil(Int, 1 + log(kpvalue / kp[1]) / loglkp)
        kpi = clamp(kpi, 1, kpl)

        mi = ceil(Int, 1 + log(mvalue / m[1]) / loglm)
        mi = clamp(mi, 1, ml)

        nkvalue = c_o[kpi, mi] - alphakp[kpi, mi] * kpvalue -
                  alpham[kpi, mi] * mvalue - beta[kpi, mi] * kpvalue * mvalue
    else
        nkvalue = 0.0
    end

    return max(nkvalue, 0.0)
end

