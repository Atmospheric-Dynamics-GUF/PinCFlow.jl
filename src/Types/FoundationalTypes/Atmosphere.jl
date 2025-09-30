"""
```julia
Atmosphere{A <: AbstractArray{<:AbstractFloat, 3}}
```

Composite type for atmospheric background fields.

```julia
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
)::Atmosphere
```

Create an `Atmosphere` instance by dispatching to a method specific for the background and dynamic equations set in `namelists`.

```julia
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Boussinesq,
    background::UniformBoussinesq,
)::Atmosphere
```

Create an `Atmosphere` instance with background fields describing a uniform (i.e. neutral) Boussinesq atmosphere.

The background fields are given by

```math
\\begin{align*}
    \\overline{\\rho} & = \\rho_0, & \\overline{\\theta} & = \\theta_0, & P & = \\overline{\\rho} \\overline{\\theta}, & N^2 & = 0,
\\end{align*}
```

where ``\\rho_0`` and ``\\theta_0`` are given by `constants.rhoref` and `namelists.atmosphere.potential_temperature`, respectively.

```julia
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Boussinesq,
    background::StratifiedBoussinesq,
)::Atmosphere
```

Create an `Atmosphere` instance with background fields describing a stratified Boussinesq atmosphere.

The background fields are given by

```math
\\begin{align*}
    \\overline{\\rho} & = \\rho_0, & \\overline{\\theta} & = \\theta_0, & P & = \\overline{\\rho} \\overline{\\theta}, & N^2 & = N_0^2,
\\end{align*}
```

where ``\\rho_0``, ``\\theta_0`` and ``N_0`` are given by `constants.rhoref`, `namelists.atmosphere.potential_temperature` and `namelists.atmosphere.buoyancy_frequency`, respectively.

```julia
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::AbstractModel,
    background::Isothermal,
)::Atmosphere
```

Create an `Atmosphere` instance with background fields describing an isothermal atmosphere.

The background fields are given by

```math
\\begin{align*}
    P \\left(z\\right) & = p_0 \\exp \\left(- \\frac{\\sigma z}{\\gamma T_0}\\right),\\\\
    \\overline{\\theta} \\left(z\\right) & = T_0 \\exp \\left(\\frac{\\kappa \\sigma z}{T_0}\\right),\\\\
    \\overline{\\rho} \\left(z\\right) & = \\frac{P \\left(z\\right)}{\\overline{\\theta} \\left(z\\right)},\\\\
    N^2 & = \\frac{g}{\\overline{\\theta}} \\frac{\\overline{\\theta}_{k + 1} - \\overline{\\theta}_{k - 1}}{2 J \\Delta \\widehat{z}},
\\end{align*}
```

where ``p_0``, ``T_0``, ``\\sigma``, ``\\gamma`` and ``\\kappa`` are given by `namelists.atmosphere.ground_pressure`, `namelists.atmosphere.temperature`, `constants.sig`, `constants.gamma` and `constants.kappa`, respectively.

# Fields

  - `pstrattfc::A`: Mass-weighted potential temperature ``P \\left(z\\right)`` (``P \\left(x, y, z, t\\right)`` in compressible mode).

  - `thetastrattfc::A`: Background potential temperature ``\\overline{\\theta} \\left(z\\right)``.

  - `rhostrattfc::A`: Background density ``\\overline{\\rho} \\left(z\\right)``.

  - `bvsstrattfc::A`: Squared buoyancy frequency ``N^2 \\left(z\\right)``.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `grid`: Collection of parameters and fields that describe the grid.

  - `model`: Dynamic equations.

  - `background`: Atmospheric background.

# See also

  - [`PinCFlow.Types.FoundationalTypes.set_vertical_boundaries_of_field!`](@ref)
"""
struct Atmosphere{A <: AbstractArray{<:AbstractFloat, 3}}
    pstrattfc::A
    thetastrattfc::A
    rhostrattfc::A
    bvsstrattfc::A
end

function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
)::Atmosphere
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
)::Atmosphere
    (; potential_temperature) = namelists.atmosphere
    (; thetaref) = constants
    (; nxx, nyy, nzz) = domain

    # Set the background fields.
    rhostrattfc = ones(nxx, nyy, nzz)
    thetastrattfc = potential_temperature ./ thetaref .* ones(nxx, nyy, nzz)
    pstrattfc = rhostrattfc .* thetastrattfc
    bvsstrattfc = zeros(nxx, nyy, nzz)

    return Atmosphere(pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc)
end

function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Boussinesq,
    background::StratifiedBoussinesq,
)::Atmosphere
    (; buoyancy_frequency, potential_temperature) = namelists.atmosphere
    (; tref, thetaref) = constants
    (; nxx, nyy, nzz) = domain

    # Set the background fields.
    rhostrattfc = ones(nxx, nyy, nzz)
    thetastrattfc = potential_temperature ./ thetaref .* ones(nxx, nyy, nzz)
    pstrattfc = rhostrattfc .* thetastrattfc
    bvsstrattfc = (buoyancy_frequency .* tref) .^ 2 .* ones(nxx, nyy, nzz)

    return Atmosphere(pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc)
end

function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::AbstractModel,
    background::Isothermal,
)::Atmosphere
    (; nbz) = namelists.domain
    (; temperature, ground_pressure) = namelists.atmosphere
    (; thetaref, pref, kappa, sig, gamma, g_ndim) = constants
    (; ndzz, nxx, nyy, nzz, ko, k0, k1) = domain
    (; ztfc, jac, dz) = grid

    # Initialize the background fields.
    (pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc) =
        (zeros(nxx, nyy, nzz) for i in 1:4)

    t0 = temperature / thetaref
    p0 = ground_pressure / pref

    # Compute the background fields.
    pstrattfc .= p0 .* exp.(.-sig .* ztfc ./ gamma ./ t0)
    thetastrattfc .= t0 .* exp.(kappa .* sig ./ t0 .* ztfc)
    rhostrattfc .= pstrattfc ./ thetastrattfc

    # Compute the squared buoyancy frequency.
    bvsstrattfc .= 0.0
    @ivy for k in k0:k1
        bvsstrattfc[:, :, k] .=
            g_ndim ./ thetastrattfc[:, :, k] ./ jac[:, :, k] .* 0.5 .*
            (thetastrattfc[:, :, k + 1] .- thetastrattfc[:, :, k - 1]) ./ dz
    end

    # Compute the squared buoyancy frequency at the boundaries.
    set_vertical_boundaries_of_field!(bvsstrattfc, namelists, domain, +)
    @ivy if ko == 0
        for k in 1:nbz
            bvsstrattfc[:, :, k] .=
                g_ndim ./ thetastrattfc[:, :, k0 - 1] ./ jac[:, :, k0 - 1] .*
                (thetastrattfc[:, :, k0] .- thetastrattfc[:, :, k0 - 1]) ./ dz
        end
    end
    @ivy if ko + nzz == ndzz
        for k in 1:nbz
            bvsstrattfc[:, :, k1 + k] .=
                g_ndim ./ thetastrattfc[:, :, k1 + 1] ./ jac[:, :, k1 + 1] .*
                (thetastrattfc[:, :, k1 + 1] .- thetastrattfc[:, :, k1]) ./ dz
        end
    end

    return Atmosphere(pstrattfc, thetastrattfc, rhostrattfc, bvsstrattfc)
end
