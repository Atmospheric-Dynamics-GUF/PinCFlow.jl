function update_interpolation_coef_discrete! end

function update_interpolation_coef_discrete!(
    spec_tend::TriadTendencies,
    nk::AbstractMatrix{<:AbstractFloat},
    ::Triad2D,
)
    (; m, ml) = spec_tend.spec_grid
    (; c_o, alpham) = spec_tend.interp_coef

    if m[1] > 0.0
        error("Discrete-k interpolation is currently implemented for signed-m Triad2D.")
    end

    mhalf = ml ÷ 2
    iz1 = mhalf
    iz2 = mhalf + 1

    @assert mhalf >= 2

    @ivy for kpi in axes(nk, 1)

        #------------------------------------------------------
        # Positive-m branch
        #------------------------------------------------------

        # Inner extrapolation interval.
        mi1 = iz2
        mi2 = iz2 + 1

        m1 = abs(m[mi1])
        m2 = abs(m[mi2])

        slope = (nk[kpi, mi2] - nk[kpi, mi1]) / (m2 - m1)

        alpham[kpi, iz2] = -slope
        c_o[kpi, iz2] = nk[kpi, mi1] - slope * m1

        # Interpolation intervals between positive-m centres.
        for mi in (iz2 + 1):ml
            mi1 = mi - 1
            mi2 = mi

            m1 = abs(m[mi1])
            m2 = abs(m[mi2])

            slope = (nk[kpi, mi2] - nk[kpi, mi1]) / (m2 - m1)

            alpham[kpi, mi] = -slope
            c_o[kpi, mi] = nk[kpi, mi1] - slope * m1
        end

        #------------------------------------------------------
        # Negative-m branch
        #------------------------------------------------------

        # Inner extrapolation interval.
        mi1 = iz1
        mi2 = iz1 - 1

        m1 = abs(m[mi1])
        m2 = abs(m[mi2])

        slope = (nk[kpi, mi2] - nk[kpi, mi1]) / (m2 - m1)

        alpham[kpi, iz1] = -slope
        c_o[kpi, iz1] = nk[kpi, mi1] - slope * m1

        # Interpolation intervals between negative-m centres.
        for mi in 1:(iz1 - 1)
            mi1 = mi + 1
            mi2 = mi

            m1 = abs(m[mi1])
            m2 = abs(m[mi2])

            slope = (nk[kpi, mi2] - nk[kpi, mi1]) / (m2 - m1)

            alpham[kpi, mi] = -slope
            c_o[kpi, mi] = nk[kpi, mi1] - slope * m1
        end
    end

    return nothing
end