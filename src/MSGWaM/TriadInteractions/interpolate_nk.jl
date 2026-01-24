function interpolate_nk end

function interpolate_nk(spec_tend::TriadTendencies,
     kpvalue::AbstractFloat, 
     mvalue::AbstractFloat,
     triad_mode::Triad3DIso)::AbstractFloat

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

function interpolate_nk(spec_tend::TriadTendencies,
     kpvalue::AbstractFloat, 
     mvalue::AbstractFloat,
     triad_mode::Triad2D)::AbstractFloat

    (; kp, m, kpl, ml, loglkp, loglm) = spec_tend.spec_grid
    (; c_o, alphakp, alpham, beta) = spec_tend.interp_coef

    if m[1] > 0
        return interpolate_nk(spec_tend, kpvalue, mvalue, Triad3DISo())
    end

    if kpvalue > 0 && kpvalue <= kp[end] && mvalue > 0 && mvalue <= m[end] 
        kpi = ceil(Int, 1 + log(kpvalue / kp[1]) / loglkp)
        kpi = clamp(kpi, 1, kpl)
        mi = ceil(Int, 1 + log(mvalue / m[Int(ml / 2 + 1)]) / loglm)
        mi = clamp(mi, 1, Int(ml / 2))
        mi =  Int(ml / 2 + mi) #adjusting the correct index for complete non positive m grid
        nkvalue = c_o[kpi, mi] - alphakp[kpi, mi] * kpvalue - 
                alpham[kpi, mi] * mvalue - beta[kpi, mi] * kpvalue * mvalue
    
    elseif kpvalue > 0 && kpvalue <= kp[end] && mvalue < 0 && mvalue >= m[1]
        kpi = ceil(Int, 1 + log(kpvalue / kp[1]) / loglkp)
        kpi = clamp(kpi, 1, kpl)
        mi = ceil(Int, 1 + log(abs(mvalue) / m[Int(ml / 2 + 1)]) / loglm)
        mi = clamp(mi, 1, ml)
        mi = Int(ml / 2 - mi + 1)  #reflection for negative mvalue
        nkvalue = c_o[kpi, mi] - alphakp[kpi, mi] * kpvalue - 
                alpham[kpi, mi] * abs(mvalue) - beta[kpi, mi] * kpvalue * abs(mvalue)
    else
        nkvalue = 0.0
    end
    return max(nkvalue, 0.0)
end
