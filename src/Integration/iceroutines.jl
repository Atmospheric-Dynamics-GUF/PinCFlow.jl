
function sat_ratio(Qv::AbstractFloat, p::AbstractFloat, p_si::AbstractFloat, rhoMean::AbstractFloat, cons::IceConstants)
    (; epsil0) = cons # episl0hat = epsil0 / pRef
    return Qv / rhoMean * p / p_si / epsil0
end     

function dot_qv(S::AbstractFloat, N::AbstractFloat, T::AbstractFloat, p::AbstractFloat, p_si::AbstractFloat, cons::IceConstants)
    (; DepS) = cons # Dep = DepS / pRef
    return - DepS * (S - 1.0) * T * N * p_si / p
end    

function dot_n(S::AbstractFloat, rho::AbstractFloat, cons::IceConstants)
    (; J_nuc, B_nuc, S_c) = cons
    return rho * J_nuc * exp(B_nuc * (S - S_c))
end


function psat_ice(T::AbstractFloat, iceconstants::IceConstants)
    # compute saturation pressure over ice
    # NB: temperature nondimensionalized with thetaRef_tropopause
    # bc. Psat_ice_ref = Psat_ice_ref (thetaRef_tropopause)

    (; L_hat, thetaRefRatio, pRef) = iceconstants

    return 1/pRef * exp(L_hat * (1. - 1. / (T * thetaRefRatio)))
    # changes: added 1/pRef to correctly non dimensionalize the pressure
  end

 
#=
function psat_ice(T::AbstractFloat, iceconstants::IceConstants)
    # compute saturation pressure over ice

    (; Li_hat, es0_hat, T0_hat) = iceconstants

    return es0_hat*exp(Li_hat * (1. / T0_hat - 1. / T))
  end
=#

  

 function nIce_param_nuc(N::AbstractFloat, S::AbstractFloat, dotS::AbstractFloat, T::AbstractFloat, cons::IceConstants)
    # compute nucleated number of ice crystals
    # using the asymptotic solution
    (; S_c, DepS, epsil0) = cons

    #CHANGES DepS --> Dep
    #NIce_param_nuc = 2. * dotS / (Dep*1.0e7 * (S - 1.) * T) - N
    NIce_param_nuc = 2. * dotS / (DepS/epsil0 * (S - 1.) * T) - N
    #
    if NIce_param_nuc < N
    # CHANGES  
    #  NIce_param_nuc = NIce_param_nuc - N
      NIce_param_nuc = N
    end
    return NIce_param_nuc

  end 

  # for heterogeneous nucleation
  function ni(T::AbstractFloat, rho:: AbstractFloat, cons::IceConstants)
    # computes the activated IN at a given temperature (per kg)
    # temperature nondimensionalized with tRef (0.2 has units of K^-1)
    # rho is also nondimensional
    (; Ni0_hat, T0_hat, thetaRef) = cons
    return Ni0_hat/rho * exp(0.2*thetaRef*(T0_hat - T))
  end 

  function psat_water(T::AbstractFloat, cons::IceConstants)
    # computes equilibrium vapour pressure over water (COSMO documentation eq. 5.33)
    # pressure nondimensionalized with p0v_hat
    # temperature dimensionalized with tRef
    (; p0v_hat, aw, bw_hat, Tr_hat) = cons
    return p0v_hat * exp(aw * (T - Tr_hat) / (T - bw_hat))
  end

  function psat_ice2(T::AbstractFloat, cons::IceConstants)
    # computes equilibrium vapour pressure over ice (COSMO documentation eq. 5.35)
    # pressure nondimensionalized with p0v_hat
    # temperature dimensionalized with tRef
    (; p0v_hat, ai, bi_hat, Tr_hat) = cons
    return p0v_hat * exp(ai * (T - Tr_hat) / (T - bi_hat))
  end

  function qv_sw(p::AbstractFloat, p_sw::AbstractFloat, cons::IceConstants)
    # computes saturation mass fraction over water (Lohmann - An Introduction to Clouds 2.80 // COSMO documentation eq. 5.32)
    (; epsil0) = cons 
    return epsil0 * p_sw / (p - (1 - epsil0) * p_sw)
  end

  function qv_si(p::AbstractFloat, p_si::AbstractFloat, cons::IceConstants)
    # computes saturation mass fraction over ice (Lohmann - An Introduction to Clouds 2.80 // COSMO documentation eq. 5.34)
    # Note: p_si is not consistent with the used formula in the COSMO model (COSMO documentation eq. 5.35)
    # COSMO uses: pv_si = pv_0 * exp(ai * (T-Tr)/(T-bi)); ai=21.875, Tr = 273.16 K, bi = 35.86 K
    (; epsil0) = cons
    return epsil0 * p_si / (p - (1 - epsil0) * p_si)
  end