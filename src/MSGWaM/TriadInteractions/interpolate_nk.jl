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
    nk::AbstractMatrix{<:AbstractFloat},
    kpi::Integer,
    mvalue::AbstractFloat,
    triad_mode::Triad2D,
)::AbstractFloat

    (; m, mc, ml, loglm) = spec_tend.spec_grid

    if m[1] > 0
        error("Discrete-k interpolation is currently implemented for signed-m Triad2D.")
    end

    mhalf = ml ÷ 2
    iz2 = mhalf + 1

    @assert mhalf >= 2

    #----------------------------------------------------------
    # Check whether mvalue lies inside the represented
    # signed-m spectral domain.
    #----------------------------------------------------------

    if mvalue > 0
        if !(mc[iz2 + 1] <= mvalue <= mc[end])
            return 0.0
        end
    elseif mvalue < 0
        if !(mc[1] <= mvalue <= mc[iz2])
            return 0.0
        end
    else
        return 0.0
    end

    #----------------------------------------------------------
    # Locate mvalue on the logarithmic |m| grid.
    #----------------------------------------------------------

    mabs = abs(mvalue)

    mj = ceil(Int, 1 + log(mabs / m[iz2]) / loglm)
    mj = clamp(mj, 1, mhalf)

    # For values below the first centre, use the first two
    # centres for linear extrapolation.
    #
    # Otherwise interpolate between the neighbouring centres
    # mj - 1 and mj.
    if mj == 1
        mj1 = 1
        mj2 = 2
    else
        mj1 = mj - 1
        mj2 = mj
    end

    #----------------------------------------------------------
    # Convert the positive-|m| indices to indices of the
    # signed-m spectral grid.
    #----------------------------------------------------------

    if mvalue > 0
        mi1 = mhalf + mj1
        mi2 = mhalf + mj2
    else
        mi1 = mhalf - mj1 + 1
        mi2 = mhalf - mj2 + 1
    end

    m1 = abs(m[mi1])
    m2 = abs(m[mi2])

    n1 = nk[kpi, mi1]
    n2 = nk[kpi, mi2]

    #----------------------------------------------------------
    # Linear interpolation/extrapolation in physical |m|.
    #----------------------------------------------------------

    nkvalue = n1 + (n2 - n1) * (mabs - m1) / (m2 - m1)

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

