"""
```julia
Atmosphere{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractVector{<:AbstractFloat},
}
```

A type representing the atmospheric background state in terrain-following coordinates.

# Fields

  - `pstrattfc::A`: Background pressure ``p/p_\\mathrm{ref}``
  - `thetastrattfc::A`: Background potential temperature ``\\theta/\\theta_\\mathrm{ref}``
  - `rhostrattfc::A`: Background density ``\\rho/\\rho_\\mathrm{ref}``
  - `bvsstrattfc::A`: Brunt-Väisälä frequency squared ``N^2/N_\\mathrm{ref}^2``
  - `fc::B`: Coriolis parameter ``f/f_\\mathrm{ref}``

# Model Types

  - `Boussinesq`: Boussinesq approximation
  - Subtypes of `AbstractModel`: Other atmospheric models

# Background Types

  - `Isothermal`: Isothermal background state
  - `UniformBoussinesq`: Uniform background for Boussinesq model
  - `StratifiedBoussinesq`: Stratified background for Boussinesq model

# Notes

  - All fields are non-dimensionalized using reference values from `Constants`
  - Coordinates are in terrain-following system
  - Grid decomposition for parallel computation is handled through `Domain`
"""
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

"""
```julia
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
)
```

Creates an `Atmosphere` instance with model and background types specified in namelists.

# Arguments

  - `namelists`: Configuration settings including:

      + `atmosphere`: Background state parameters
      + `domain`: Domain configuration
      + `setting`: Model settings

  - `constants`: Physical constants and reference values
  - `domain`: Grid domain information including parallel decomposition
  - `grid`: Grid metrics and coordinate information

# Returns

  - `::Atmosphere`: `Atmosphere` instance.
"""
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

"""
```julia
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Boussinesq,
    background::UniformBoussinesq,
)
```

Creates a uniform Boussinesq atmosphere with:

  - Constant density (``\\rho = \\rho_\\mathrm{ref}``)
  - Uniform potential temperature (``\\theta = \\theta_0``)
  - Zero buoyancy frequency (``N^2 = 0``)
  - Coriolis parameter based on specified mode

# Returns

  - `::Atmosphere`: `Atmosphere` instance.

# See also

  - [`PinCFlow.Types.FoundationalTypes.compute_coriolis_parameter`](@ref)
"""
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

"""
```julia
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Boussinesq,
    background::StratifiedBoussinesq,
)
```

Creates a stratified Boussinesq atmosphere with:

  - Constant density: ``\\rho = \\rho_\\mathrm{ref}``
  - Uniform potential temperature: ``\\theta = \\theta_0``
  - Constant pressure from equation of state: ``p = \\rho \\theta``
  - Constant buoyancy frequency (``N^2 = N_0^2``)
  - Coriolis parameter based on specified mode

# Returns

  - `::Atmosphere`: `Atmosphere` instance.

# See also

  - [`PinCFlow.Types.FoundationalTypes.compute_coriolis_parameter`](@ref)
"""
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

"""
```julia
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::AbstractModel,
    background::Isothermal,
)
```

Creates an isothermal atmosphere with:

  - Exponential pressure profile ``p(z) = p_0 \\exp(-\\sigma z/\\gamma T_0)``
  - Temperature profile ``\\theta(z) = T_0 \\exp(\\kappa \\sigma z/T_0)``
  - Density from equation of state ``\\rho = p/\\theta``
  - ``N^2`` computed from vertical ``\\theta`` gradient
  - Handles boundary conditions for ``N^2`` calculation

# Returns

  - `::Atmosphere`: `Atmosphere` instance.

# See also

  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)
  - [`PinCFlow.Types.FoundationalTypes.compute_coriolis_parameter`](@ref)
"""
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
