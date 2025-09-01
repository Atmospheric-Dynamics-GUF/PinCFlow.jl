
function sat_ratio(Qv::AbstractFloat, p::AbstractFloat, p_si::AbstractFloat, rhoMean::AbstractFloat, cons::IceConstants)
    (; epsil0hat) = cons
    return Qv / rhoMean * p / p_si / epsil0hat
end     

function dot_qv(S::AbstractFloat, N::AbstractFloat, T::AbstractFloat, p::AbstractFloat, p_si::AbstractFloat, cons::IceConstants)
    (; Dep) = cons
    return - Dep * (S - 1.0) * T * N * p_si / p
end    

function dot_n(S::AbstractFloat, rho::AbstractFloat, cons::IceConstants)
    (; J_nuc, B_nuc, S_c) = cons
    return rho * J_nuc * exp(B_nuc * (S - S_c))
end

function psat_ice(T::AbstractFloat, iceconstants::IceConstants)
    # compute saturation pressure over ice
    # NB: temperature nondimensionalized with thetaRef_tropopause
    # bc. Psat_ice_ref = Psat_ice_ref (thetaRef_tropopause)

    (; L_hat, thetaRefRatio) = iceconstants

    return exp(L_hat * (1. - 1. / (T * thetaRefRatio)))

  end

 