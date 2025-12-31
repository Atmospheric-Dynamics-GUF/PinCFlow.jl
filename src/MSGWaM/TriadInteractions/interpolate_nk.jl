function interpolate_nk end

function interpolate_nk(spec_tend::TriadTendencies,
     kpvalue::AbstractFloat, 
     mvalue::AbstractFloat)::AbstractFloat

    (; kp, m, kpl, ml, loglkp, loglm) = spec_tend.spec_grid
    (; c_o, alphakp, alpham, beta) = spec_tend.interp_coef

    if kpvalue > 0 && kpvalue <= kp[end] && mvalue > 0 && mvalue <= m[end]
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