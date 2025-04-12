struct Atmosphere{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractVector{<:AbstractFloat},
}
    pstrattfc::A
    thetastrattfc::A
    rhostrattfc::A
    bvsstrattfc::A
    f_cor_nd::B
end

function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
)
    (; model) = namelists.setting
    (; background, corset) = namelists.atmosphere
    return Atmosphere(
        namelists,
        constants,
        domain,
        grid,
        model,
        background,
        corset,
    )
end

function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::PseudoIncompressible,
    background::Isothermal,
    corset::ConstantCoriolis,
)
    # Get parameters.
    (; temp0_dim, press0_dim, f_coriolis_dim) = namelists.atmosphere
    (; tref, thetaref, pref, kappa, sig, gamma, g_ndim) = constants
    (; nxx, nyy, nzz, k0, k1) = domain
    (; ztfc, jac, dz) = grid

    # Initialize the background fields.
    (pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc) =
        (zeros(nxx, nyy, nzz) for i in 1:4)

    t0 = temp0_dim / thetaref
    p0 = press0_dim / pref

    # Define 3D background fields.
    for k in (k0 - 2):(k1 + 2)
        # Define pStratTFC.
        pstrattfc[:, :, k] .= p0 .* exp.(-sig .* ztfc[:, :, k] ./ gamma ./ t0)
        # Define thetaStratTFC.
        thetastrattfc[:, :, k] .=
            t0 .* exp.(kappa .* sig ./ t0 .* ztfc[:, :, k])
        # Define rhoStratTFC.
        rhostrattfc[:, :, k] .= pstrattfc[:, :, k] ./ thetastrattfc[:, :, k]
    end

    # Define bvsStratTFC.
    bvsstrattfc .= 0.0
    # Lower boundary.
    bvsstrattfc[:, :, k0 - 2] .=
        g_ndim ./ thetastrattfc[:, :, k0 - 1] ./ jac[:, :, k0 - 1] .*
        (thetastrattfc[:, :, k0] .- thetastrattfc[:, :, k0 - 1]) ./ dz
    bvsstrattfc[:, :, k0 - 1] .= bvsstrattfc[:, :, k0 - 2]
    # Between boundaries.
    for k in k0:k1
        bvsstrattfc[:, :, k] .=
            g_ndim ./ thetastrattfc[:, :, k] ./ jac[:, :, k] .* 0.5 .*
            (thetastrattfc[:, :, k + 1] .- thetastrattfc[:, :, k - 1]) ./ dz
    end
    # Upper boundary.
    bvsstrattfc[:, :, k1 + 1] .=
        g_ndim ./ thetastrattfc[:, :, k1 + 1] ./ jac[:, :, k1 + 1] .*
        (thetastrattfc[:, :, k1 + 1] .- thetastrattfc[:, :, k1]) ./ dz
    bvsstrattfc[:, :, k1 + 2] .= bvsstrattfc[:, :, k1 + 1]

    # Set Coriolis parameter.
    f_cor_nd = zeros(nyy)
    f_cor_nd .= f_coriolis_dim .* tref

    # Return an Atmosphere instance.
    return Atmosphere(
        pstrattfc,
        thetastrattfc,
        rhostrattfc,
        bvsstrattfc,
        f_cor_nd,
    )
end
