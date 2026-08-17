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

