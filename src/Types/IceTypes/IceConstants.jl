struct IceConstants{A <: AbstractFloat}
    
    PSatIceRef::A

    theta_ref_ratio:: A
    theta_ref_trp:: A
    
    epsil0::AbstractArray
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

end    

function IceConstants(stateconstant::Constants)

    (; thetaRef, tRef, pRef, rhoRef, lRef) =  stateconstants

    return IceConstants( 
        PsatIceRef = 1, #reference saturation pressure [Pa]
        thetaRef_trp = 210., # reference temperature in the tropopause
        thetaRefRatio = thetaRef / thetaRef_trp,

        epsil0 = 0.62, #\approx Mole_mass_water/Mole_mass_dryAir
        epsil0hat = epsil0 * PsatIceRef / pRef,
        Dep0 = 4.3E-8, # C_0 coeff paper 
        J_nuc = 4.9E4, #nucleation rate
        B_nuc = 337, #nucleation exponent
        S_c = 1.5,

        L_ice = 2.8E6, # constant latent heat  ice [J/kg]
        R_v  = 461., # specific gas constant for water vapor [J/kg/K]
        mRef = rhoRef * lRef ^ 3, #reverence mass
        meanMassIce = 1.E-12, # mean mass ice crystals [kg]
        Dep = Dep0 * thetaRef * tRef * meanMassIce ^ (1. / 3.) * PsatIceRef / pRef / mRef,
        DepS =  Dep0 * thetaRef * tRef * meanMassIce ^ (1. / 3.) / mRef,
        L_hat = L_ice / R_v / theta_ref_trp, 
        Li_hat = L_ice / R_v / thetaRef )

end