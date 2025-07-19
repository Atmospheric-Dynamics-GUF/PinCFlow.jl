"""
```julia
Atmosphere{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractVector{<:AbstractFloat},
}
```

Composite type for atmospheric background fields and Coriolis frequency.

# Fields

  - `pstrattfc::A`: Mass-weighted potential temperature ``P \\left(z\\right)`` (``P \\left(x, y, z, t\\right)`` in compressible mode).
  - `thetastrattfc::A`: Background potential temperature ``\\overline{\\theta} \\left(z\\right)``.
  - `rhostrattfc::A`: Background density ``\\overline{\\rho} \\left(z\\right)``.
  - `bvsstrattfc::A`: Squared buoyancy frequency ``N^2 \\left(z\\right)``.
  - `fc::B`: Coriolis frequency ``f \\left(y\\right)``.
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

Creates an `Atmosphere` instance, depending on the background and dynamic equations specified in `namelists`.

Dispatches to specific methods based on the background and dynamic equations specified in `namelists`.

# Arguments

  - `namelists`: Namelists with all model parameters.
  - `constants`: Physical constants and reference values.
  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
  - `grid`: Collection of parameters and fields that describe the grid.

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

Returns an Atmosphere instance with background fields describing a uniform (i.e. neutral) Boussinesq atmosphere.

The background fields are given by

```math
\\begin{align*}
    \\overline{\\rho} & = \\rho_0,\\\\
    \\overline{\\theta} & = \\theta_0,\\\\
    P & = \\overline{\\rho} \\overline{\\theta},\\\\
    N^2 & = 0,
\\end{align*}
```

where ``\\rho_0`` and ``\\theta_0`` are given by `constants.rhoref` and `namelists.atmosphere.theta0_dim`, respectively. The Coriolis frequency is computed with `compute_coriolis_frequency`, depending on `namelists.atmosphere.coriolis_mode`.

# Arguments

  - `namelists`: Namelists with all model parameters.
  - `constants`: Physical constants and reference values.
  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
  - `grid`: Collection of parameters and fields that describe the grid.
  - `model`: Dynamic equations.
  - `background`: Atmospheric background.

# Returns

  - `::Atmosphere`: `Atmosphere` instance.

# See also

  - [`PinCFlow.Types.FoundationalTypes.compute_coriolis_frequency`](@ref)
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
    fc = compute_coriolis_frequency(
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

Returns an Atmosphere instance with background fields describing a stratified Boussinesq atmosphere.

The background fields are given by

```math
\\begin{align*}
    \\overline{\\rho} & = \\rho_0,\\\\
    \\overline{\\theta} & = \\theta_0,\\\\
    P & = \\overline{\\rho} \\overline{\\theta},\\\\
    N^2 & = N_0^2,
\\end{align*}
```

where ``\\rho_0``, ``\\theta_0`` and ``N_0`` are given by `constants.rhoref`, `namelists.atmosphere.theta0_dim` and `namelists.atmosphere.buoyancy_frequency`, respectively. The Coriolis frequency is computed with `compute_coriolis_frequency`, depending on `namelists.atmosphere.coriolis_mode`.

# Arguments

  - `namelists`: Namelists with all model parameters.
  - `constants`: Physical constants and reference values.
  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
  - `grid`: Collection of parameters and fields that describe the grid.
  - `model`: Dynamic equations.
  - `background`: Atmospheric background.

# Returns

  - `::Atmosphere`: `Atmosphere` instance.

# See also

  - [`PinCFlow.Types.FoundationalTypes.compute_coriolis_frequency`](@ref)
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
    fc = compute_coriolis_frequency(
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

Returns an Atmosphere instance with background fields describing an isothermal atmosphere.

The background fields are given by

```math
\\begin{align*}
    P \\left(z\\right) & = p_0 \\exp \\left(- \\frac{\\sigma z}{\\gamma T_0}\\right),\\\\
    \\overline{\\theta} \\left(z\\right) & = T_0 \\exp \\left(\\frac{\\kappa \\sigma z}{T_0}\\right),\\\\
    \\overline{\\rho} \\left(z\\right) & = \\frac{P \\left(z\\right)}{\\overline{\\theta} \\left(z\\right)},\\\\
    N^2 & = \\frac{g}{\\overline{\\theta}} \\frac{\\overline{\\theta}_{k + 1} - \\overline{\\theta}_{k - 1}}{2 J \\Delta \\widehat{z}},
\\end{align*}
```

where ``p_0``, ``T_0``, ``\\sigma``, ``\\gamma`` and ``\\kappa`` are given by `namelists.atmosphere.press0_dim`, `namelists.atmosphere.temp0_dim`, `constants.sig`, `constants.gamma` and `constants.kappa`, respectively. The Coriolis frequency is computed with `compute_coriolis_frequency`, depending on `namelists.atmosphere.coriolis_mode`.

# Arguments

  - `namelists`: Namelists with all model parameters.
  - `constants`: Physical constants and reference values.
  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
  - `grid`: Collection of parameters and fields that describe the grid.
  - `model`: Dynamic equations.
  - `background`: Atmospheric background.

# Returns

  - `::Atmosphere`: `Atmosphere` instance.

# See also

  - [`PinCFlow.Boundaries.set_vertical_boundaries_of_field!`](@ref)
  - [`PinCFlow.Types.FoundationalTypes.compute_coriolis_frequency`](@ref)
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
    fc = compute_coriolis_frequency(
        namelists,
        constants,
        domain,
        grid,
        coriolis_mode,
    )

    # Return an Atmosphere instance.
    return Atmosphere(pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc, fc)
end
