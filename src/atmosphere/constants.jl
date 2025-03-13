
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
    kappainv = 1.0 / kappa
    gammainv = 1.0 / gamma

    # Reference quantities
    rsp = 287.0
    g = 9.81
    rhoref = 1.184
    pref = 101325.0
    aref = sqrt(pref / rhoref)
    uref = aref
    lref = pref / rhoref / g
    tref = lref / aref
    thetaref = aref^2.0 / rsp
    ma = uref / aref # == 1?
    mainv2 = 1.0 / ma^2
    fr = uref / sqrt(g * lref)
    sig = ma^2 / fr^2

    fref = rhoref * uref^2 / lref
    dt = p.discretization.dtmax_dim / tref

    # isothermal atmosphere


    # Reynolds number
    if p.atmosphere.specifyreynolds
        reinv = p.atmosphere.mu_viscous_dim / (uref * lref)
    else
        reinv = 0
    end

    re = if reinv < 1.0e-20
        1.0e20
    else
        1.0 / reinv
    end
    g_ndim = g / (uref^2 / lref)

    # TODO:
    ma2 = 0.0
    frinv2 = 0.0
    fr2 = 0.0
    small = nextfloat(0.0)
    return Constants(
        gamma,
        gammainv,
        kappa,
        kappainv,
        rsp,
        g,
        rhoref,
        pref,
        aref,
        uref,
        lref,
        tref,
        thetaref,
        fref,
        g_ndim,
        re,
        ma,
        mainv2,
        ma2,
        fr,
        frinv2,
        fr2,
        sig,
        small
    )

end
