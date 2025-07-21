struct IceConstants{A <: AbstractFloat}
    
    PSatIceRef::A

    thetaRef_trp:: A
    thetaRefRatio:: A
    
    epsil0:: A
    epsil0hat:: A
    
    Dep0:: A
    J_nuc:: A
    B_nuc:: A
    S_c:: A
    
    L_ice:: A
    R_v:: A
    mRef:: A
    meanMassIce:: A
    
    Dep:: A
    DepS:: A
    L_hat:: A
    Li_hat:: A

    n :: A
    q :: A
    qv :: A

end    

function IceConstants(constants::Constants)

    (; thetaref, tref, pref, rhoref, lref) =  constants

    thetaRef = thetaref
    tRef = tref
    pRef = pref
    rhoRef = rhoref
    lRef = lref
    PsatIceRef = 1.               # reference saturation pressure [Pa]
    thetaRef_trp = 210.          # reference temperature in the tropopause
    thetaRefRatio = thetaRef / thetaRef_trp

    epsil0 = 0.62                # approx Mole_mass_water/Mole_mass_dryAir
    epsil0hat = epsil0 * PsatIceRef / pRef
    Dep0 = 4.3E-8                # C_0 coeff paper 
    J_nuc = 4.9E4                # nucleation rate
    B_nuc = 337.                  # nucleation exponent
    S_c = 1.5

    L_ice = 2.8E6                # constant latent heat ice [J/kg]
    R_v  = 461.                  # specific gas constant for water vapor [J/kg/K]
    mRef = rhoRef * lRef ^ 3     # reference mass
    meanMassIce = 1.E-12         # mean mass ice crystals [kg]
    Dep = Dep0 * thetaRef * tRef * meanMassIce ^ (1. / 3.) * PsatIceRef / pRef / mRef
    DepS = Dep0 * thetaRef * tRef * meanMassIce ^ (1. / 3.) / mRef
    L_hat = L_ice / R_v / thetaRef_trp
    Li_hat = L_ice / R_v / thetaRef
    
    J_nuc = J_nuc * tRef * mRef
    
    # units for dimensionalization predictands
    n = 1. / mRef 
    q = 1.
    qv = 1.

    return IceConstants(
    PsatIceRef,
    thetaRef_trp,
    thetaRefRatio,
    epsil0,
    epsil0hat,
    Dep0,
    J_nuc,
    B_nuc,
    S_c,
    L_ice,
    R_v,
    mRef,
    meanMassIce,
    Dep,
    DepS,
    L_hat,
    Li_hat,
    n,
    q,
    qv
    )
end