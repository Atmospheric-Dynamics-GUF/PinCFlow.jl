struct IceConstants{A <: AbstractFloat}
    
    #PSatIceRef::A # ,maybe change to PsatIceRef::A
    PSatIceRef::A

    thetaRef_trp:: A
    thetaRefRatio:: A
    tRef :: A
    thetaRef :: A
    
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

    n_hom :: A
    q_hom :: A
    qv :: A
    n_in :: A
    n_het :: A
    q_het :: A

    # constants for heterogeneous nucleation

    mi0 :: A
    Td_hat :: A
    Tnuc_hat :: A 
    T0_hat :: A
    Ni0 :: A 
    p0v :: A 
    aw :: A 
    bw_hat :: A
    ai :: A 
    bi_hat :: A
    Tr_hat :: A
    Nref :: A
    Ni0_hat :: A
    p0v_hat :: A
    es0_hat :: A
    mi0_hat :: A
    rhoRef :: A
    nRef :: A
    lRef :: A
    pRef :: A
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

    # constants for heterogeneous nucleation
    
    mi0 = 1.e-12                # initial mass of the heterogeneously nucleated cloud ice crystals [kg]
    #Td = 248.15                # general threshold temperature [K]#
    Td_hat = 248.15 / thetaRef
    Tnuc = 267.15               # nucleation threshold temperature [K]
    Tnuc_hat = Tnuc / thetaRef
    T0 = 273.15                   # freezing point temperature [K]
    T0_hat = T0 / thetaRef
    Ni0 = 1.0e2                 # number density of cloud ice crystals [m^-3]
    p0v = 610.78                # pressure constants for vapour pressure over a plane surface of ice/water [Pa]
    es0 = 611.2                 # pressure constants for vapour pressure over a plane surface of ice [Pa] (Lohmann 2.60)
    aw = 17.27                  # dimensionless
    bw = 35.86                  # [K]
    bw_hat = bw / thetaRef
    ai = 21.875                 # dimensionsless
    bi = 7.66                   # [K]
    bi_hat = bi / thetaRef
    Tr = 273.16                 # [K]
    Tr_hat = Tr / thetaRef
    Nref = 1.0 / (lRef^3)
    nRef = 1.0 / mRef

    # Dimensionless version of constants
    Ni0_hat = Ni0 / Nref
    p0v_hat = p0v / pRef
    mi0_hat = mi0 / mRef
    es0_hat = es0 / pRef
    # units for dimensionalization predictands
    n_hom = 1. / mRef 
    q_hom = 1.
    qv = 1.
    n_in = 1. / mRef
    n_het = 1. / mRef
    q_het = 1.

    return IceConstants(
        PsatIceRef,

        thetaRef_trp,
        thetaRefRatio,
        tRef,
        thetaRef,

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

        n_hom,
        q_hom,
        qv,
        n_in,
        n_het,
        q_het,
        
        mi0,
        Td_hat,
        Tnuc_hat,
        T0_hat,
        Ni0,
        p0v,
        aw,
        bw_hat,
        ai,
        bi_hat,
        Tr_hat,
        Nref,
        Ni0_hat,
        p0v_hat,
        es0_hat,
        mi0_hat,
        rhoRef,
        nRef,
        lRef,
        pRef
    )
end