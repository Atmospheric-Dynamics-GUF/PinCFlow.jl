struct Constants{F<:AbstractFloat}
    # Natural constants.
    gamma::F
    gammainv::F
    kappa::F
    kappainv::F
    rsp::F
    g::F

    # Reference quantities.
    rhoref::F
    pref::F
    aref::F
    uref::F
    lref::F
    tref::F
    thetaref::F
    fref::F

    # Non-dimensionalized gravitational acceleration.
    g_ndim::F

    # Flow parameters.
    re::F
    ma::F
    mainv2::F
    ma2::F
    fr::F
    frinv2::F
    fr2::F
    sig::F

    # Small float.
    small::F
end

function Constants(p::Parameters)
    # Constants
    gamma = 1.4
    gamma_1 = gamma - 1.0
    kappa = (gamma - 1.0) / gamma
    # TODO: do we really want to keep all these constants? can be calculated on the fly
    kappaInv = 1.0 / kappa

    gammaInv = 1.0 / gamma

    # Reference quantities
    Rsp = 287.0
    g = 9.81
    rhoRef = 1.184
    pRef = 101325.0
    aRef = sqrt(pRef / rhoRef)
    uRef = aRef
    lRef = pRef / rhoRef / g
    tRef = lRef / aRef
    thetaRef = aRef^2.0 / Rsp
    Ma = uRef / aRef # == 1?
    MaInv2 = 1.0 / Ma^2
    Fr = uRef / sqrt(g * lRef)
    sig = Ma^2 / Fr^2

    FRef = rhoRef * uRef^2 / lRef

    dt = p.discretization.dtmax_dim / tRef

    # isothermal atmosphere


    # Reynolds number
    if p.atmosphere.specifyreynolds
        ReInv = p.atmosphere.mu_viscous_dim / (uRef * lRef)
    end

    Re = if ReInv < 1.0e-20
        1.0e20
    else
        1.0 / ReInv
    end

    g_ndim = g / (uRef^2 / lRef)

    return Constants(gamma, gammaInv, kappa, kappaInv, Rsp, g,
        rhoRef, pRef, aRef, uRef, lRef, tRef, thetaRef, FRef,
        g * tRef^2 / lRef,
        g_ndim,
        Re, # Re (placeholder)
        Ma, MaInv2, Ma^2, Fr, 1 / Fr^2, Fr^2, sig,
        1e-14) # small

end



# struct Model
#     parameters::Parameters
#     constants::Constants
#     time::Time
#     domain::Domain
#     grid::Grid
#     atmosphere::Atmosphere
#     operator::Operator
#     variables::Variables
#     fluxes::Fluxes
# end

# function Model(p::Parameters)

# end
