struct Atmosphere{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractVector{<:AbstractFloat},
}
    pstrattfc::A
    thetastrattfc::A
    rhostrattfc::A
    bvsstrattfc::A
    fc::B
end

function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
)
    (; model) = namelists.setting
    (; background) = namelists.atmosphere
    return Atmosphere(namelists, constants, domain, grid, model, background)
end

function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Boussinesq,
    background::UniformBoussinesq,
)
    (; theta0_dim, coriolis_mode) = namelists.atmosphere
    (; thetaref) = constants
    (; nxx, nyy, nzz) = domain

    # Set the background fields.
    rhostrattfc = ones(nxx, nyy, nzz)
    thetastrattfc = theta0_dim ./ thetaref .* ones(nxx, nyy, nzz)
    pstrattfc = rhostrattfc .* thetastrattfc
    bvsstrattfc = zeros(nxx, nyy, nzz)

    # Set the Coriolis parameter.
    fc = compute_coriolis_parameter(
        namelists,
        constants,
        domain,
        grid,
        coriolis_mode,
    )

    # Return an Atmosphere instance.
    return Atmosphere(pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc, fc)
end

function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Boussinesq,
    background::StratifiedBoussinesq,
)
    (; buoyancy_frequency, theta0_dim, coriolis_mode) = namelists.atmosphere
    (; tref, thetaref) = constants
    (; nxx, nyy, nzz) = domain

    # Set the background fields.
    rhostrattfc = ones(nxx, nyy, nzz)
    thetastrattfc = theta0_dim ./ thetaref .* ones(nxx, nyy, nzz)
    pstrattfc = rhostrattfc .* thetastrattfc
    bvsstrattfc = (buoyancy_frequency .* tref) .^ 2 .* ones(nxx, nyy, nzz)

    # Set the Coriolis parameter.
    fc = compute_coriolis_parameter(
        namelists,
        constants,
        domain,
        grid,
        coriolis_mode,
    )

    # Return an Atmosphere instance.
    return Atmosphere(pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc, fc)
end

function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::AbstractModel,
    background::Isothermal,
)
    (; nbz) = namelists.domain
    (; zboundaries) = namelists.setting
    (; temp0_dim, press0_dim, coriolis_mode) = namelists.atmosphere
    (; thetaref, pref, kappa, sig, gamma, g_ndim) = constants
    (; sizezz, nxx, nyy, nzz, ko, k0, k1) = domain
    (; ztfc, jac, dz) = grid

    # Initialize the background fields.
    (pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc) =
        (zeros(nxx, nyy, nzz) for i in 1:4)

    t0 = temp0_dim / thetaref
    p0 = press0_dim / pref

    # Compute the background fields.
    pstrattfc .= p0 .* exp.(-sig .* ztfc ./ gamma ./ t0)
    thetastrattfc .= t0 .* exp.(kappa .* sig ./ t0 .* ztfc)
    rhostrattfc .= pstrattfc ./ thetastrattfc

    # Compute the squared buoyancy frequency.
    bvsstrattfc .= 0.0
    for k in k0:k1
        bvsstrattfc[:, :, k] .=
            g_ndim ./ thetastrattfc[:, :, k] ./ jac[:, :, k] .* 0.5 .*
            (thetastrattfc[:, :, k + 1] .- thetastrattfc[:, :, k - 1]) ./ dz
    end

    # Compute the squared buoyancy frequency at the boundaries.
    set_vertical_boundaries_of_field!(
        bvsstrattfc,
        namelists,
        domain,
        zboundaries,
        +,
    )
    if ko == 0
        for k in 1:nbz
            bvsstrattfc[:, :, k] .=
                g_ndim ./ thetastrattfc[:, :, k0 - 1] ./ jac[:, :, k0 - 1] .*
                (thetastrattfc[:, :, k0] .- thetastrattfc[:, :, k0 - 1]) ./ dz
        end
    end
    if ko + nzz == sizezz
        for k in 1:nbz
            bvsstrattfc[:, :, k1 + k] .=
                g_ndim ./ thetastrattfc[:, :, k1 + 1] ./ jac[:, :, k1 + 1] .*
                (thetastrattfc[:, :, k1 + 1] .- thetastrattfc[:, :, k1]) ./ dz
        end
    end

    # Set the Coriolis parameter.
    fc = compute_coriolis_parameter(
        namelists,
        constants,
        domain,
        grid,
        coriolis_mode,
    )

    # Return an Atmosphere instance.
    return Atmosphere(pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc, fc)
end
