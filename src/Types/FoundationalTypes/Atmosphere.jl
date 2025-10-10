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
    model::Union{PseudoIncompressible, Compressible},
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

```julia 
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Union{PseudoIncompressible, Compressible},
    background::Isentropic,
)::Atmosphere
```

Create an `Atmosphere` instance with background fields describing an isentropic atmosphere.

The background fields are given by 

```math
\\begin{align*}
    P \\left(z\\right) & = p_0 \\left( 1 - \\frac{\\kappa\\sigma z}{\\theta_0}\\right)^{\\frac{1}{\\gamma - 1}} \\;, \\\\
    \\overline{\\theta} & = \\theta_0 \\;, \\\\
    \\overline{\\rho}\\left(z\\right) & = \\frac{P \\left(z\\right)}{\\overline{\\theta} \\left(z\\right)}\\;,\\\\
    N^2 & = 0 \\;, 
\\end{align*}
```

where ``p_0``, ``\\theta_0``, ``\\sigma``, ``\\gamma`` and ``\\kappa`` are given by `namelists.atmosphere.ground_pressure`, `namelists.atmosphere.potential_temperature`, `constants.sig`, `constants.gamma` and `constants.kappa`, respectively.

```julia 
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Union{PseudoIncompressible, Compressible},
    background::Realistic,
)::Atmosphere
```

Create an `Atmosphere` instance with background fields describing a realistic atmosphere with an isentropic troposphere, an isothermal stratosphere, and a tropopause located at the altitude ``z_{\\mathrm{TP}}``.

The background fields are given by 

```math 
\\begin{align*}
    P \\left(z \\right) & = 
    \\begin{cases} 
        p_0 \\left( 1 - \\frac{\\kappa\\sigma z}{\\theta_0}\\right)^{\\frac{1}{\\gamma - 1}} & z \\leq z_{\\mathrm{TP}}\\;, \\\\
        p_0^{\\kappa} p_{\\mathrm{TP}}^{1/\\gamma}\\exp\\left[-\\frac{\\sigma(z-z_{\\mathrm{TP}})}{\\gamma T_{\\mathrm{TP}}}\\right] & z > z_{\\mathrm{TP}} \\;,
    \\end{cases} \\\\
    \\overline{\\theta}\\left(z\\right) & = 
    \\begin{cases}
        \\theta_0 & z \\leq z_{\\mathrm{TP}} \\;, \\\\
        \\theta_0 \\exp\\left[\\frac{\\kappa\\sigma(z-z_{\\mathrm{TP}})}{T_{\\mathrm{TP}}}\\right] & z > z_{\\mathrm{TP}} \\;,
    \\end{cases} \\\\
    \\overline{\\rho}\\left(z\\right) & = \\frac{P \\left(z\\right)}{\\overline{\\theta} \\left(z\\right)}\\;,\\\\
    N^2 & = \\frac{g}{\\overline{\\theta}} \\frac{\\overline{\\theta}_{k + 1} - \\overline{\\theta}_{k - 1}}{2 J \\Delta \\widehat{z}}\\;, 
\\end{align*}
```
and
```math
\\begin{align*}
    p_{\\mathrm{TP}} & = p_0 \\left(1 - \\frac{\\kappa\\sigma z_{\\mathrm{TP}}}{\\theta_0}\\right)^{\\frac{1}{\\gamma - 1}} \\;, \\\\
    T_{\\mathrm{TP}} & = \\theta_0 \\left(\\frac{p_{\\mathrm{TP}}}{p_0}\\right)^{\\kappa}\\;,
\\end{align*}
```

where ``p_0``, ``\\theta_0``, ``z_{\\mathrm{TP}}``, ``\\sigma``, ``\\gamma`` and ``\\kappa`` are given by `namelists.atmosphere.ground_pressure`, `namelists.atmosphere.potential_temperature`, `namelists.atmosphere.tropopause_height`, `constants.sig`, `constants.gamma` and `constants.kappa`, respectively.

```julia 
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Union{PseudoIncompressible, Compressible},
    background::LapseRates,
)::Atmosphere
```

Create an `Atmosphere` instance with background fields describing a troposphere and a stratosphere with lapse rates ``\\Gamma_{\\mathrm{TS}}`` and ``\\Gamma_{\\mathrm{TS}}``, respectively, and a tropopause located at the altitude ``z_{\\mathrm{TP}}``.

The background fields are given by 
```math 
\\begin{align*}
    T\\left(z\\right) & = 
    \\begin{cases}
        T_0 - \\Gamma_{\\mathrm{TS}} z & z \\leq z_{\\mathrm{TP}} \\;, \\\\
        T_0 - \\Gamma_{\\mathrm{TS}} z_{\\mathrm{TP}} - \\Gamma_{\\mathrm{SS}} \\left(z - z_{\\mathrm{TP}}\\right) & z > z_{\\mathrm{TP}} \\;,
    \\end{cases} \\\\
    P\\left(z\\right) & = 
    \\begin{cases}
        p_0 \\left(1 - \\frac{\\Gamma_{\\mathrm{TS}} z}{T_0}\\right)^{\\frac{g}{R \\Gamma_{\\mathrm{TS}}}} & z \\leq z_{\\mathrm{TP}} \\; \\& \\; \\Gamma_{\\mathrm{TS}} \\neq 0 \\;, \\\\
        p_0 \\exp\\left(- \\frac{z \\sigma}{\\gamma T_0} \\right) & z \\leq z_{\\mathrm{TP}} \\; \\& \\; \\Gamma_{\\mathrm{TS}} = 0 \\;, \\\\
        p_{\\mathrm{TP}}\\left[1 - \\frac{\\Gamma_{\\mathrm{SS}} \\left(z - z_{\\mathrm{TP}} \\right)}{T_{\\mathrm{TP}}} \\right]^{\\frac{g}{R\\Gamma_{\\mathrm{SS}}}} & z > z_{\\mathrm{TP}} \\; \\& \\; \\Gamma_{\\mathrm{SS}} \\neq 0 \\;, \\\\
        p_{\\mathrm{TP}}\\exp\\left[- \\frac{\\left(z - z_{\\mathrm{TP}} \\right)\\sigma}{\\gamma T_{\\mathrm{TP}}} \\right] & z > z_{\\mathrm{TP}} \\; \\& \\; \\Gamma_{\\mathrm{SS}} = 0 \\;,
    \\end{cases} \\\\
    \\overline{\\theta}\\left(z\\right) & =
    \\begin{cases}
        T\\left(z\\right) \\left[\\frac{p_0}{P\\left(z\\right)}\\right]^{\\frac{R\\Gamma_{\\mathrm{TS}}}{g}} & z \\leq z_{\\mathrm{TP}} \\; \\& \\; \\Gamma_{\\mathrm{TS}} \\neq 0 \\;, \\\\
        T_0 \\exp\\left(\\frac{\\kappa\\sigma z}{T_0}\\right) & z \\leq z_{\\mathrm{TP}} \\;\\&\\; \\Gamma_{\\mathrm{TS}} = 0 \\;, \\\\
        T\\left(z\\right)\\left[\\frac{p_{\\mathrm{TP}}}{P\\left(z\\right)}\\right]^{\\frac{R\\Gamma_{\\mathrm{SS}}}{g}} & z > z_{\\mathrm{TP}} \\; \\& \\; \\Gamma_{\\mathrm{SS}} \\neq 0 \\;, \\\\
        \\theta_{\\mathrm{TP}}\\exp\\left[\\frac{\\kappa\\sigma \\left(z-z_{\\mathrm{TP}}\\right)}{T\\left(z_{\\mathrm{TP}}\\right)}\\right] & z > z_{\\mathrm{TP}} \\;\\&\\; \\Gamma_{\\mathrm{SS}} = 0 \\;, 
    \\end{cases} \\\\
    \\overline{\\rho}\\left(z\\right) & = \\frac{P \\left(z\\right)}{\\overline{\\theta} \\left(z\\right)}\\;,\\\\
    N^2 & = \\frac{g}{\\overline{\\theta}} \\frac{\\overline{\\theta}_{k + 1} - \\overline{\\theta}_{k - 1}}{2 J \\Delta \\widehat{z}}\\;, 
\\end{align*}
```
and
```math
\\begin{align*}
    p_{\\mathrm{TP}} & = 
    \\begin{cases}
        p_0 \\left(1 - \\frac{\\Gamma_{\\mathrm{TS}} z_{\\mathrm{TP}}}{T_0}\\right)^{\\frac{g}{R \\Gamma_{\\mathrm{TS}}}} & \\Gamma_{\\mathrm{TS}} \\neq 0 \\;, \\\\
        p_0 \\exp\\left(- \\frac{z_{\\mathrm{TP}}\\sigma}{\\gamma T_0} \\right) & \\Gamma_{\\mathrm{TS}} = 0 \\;, 
    \\end{cases} \\\\
    \\theta_{\\mathrm{TP}} &= T_0 \\exp\\left(\\frac{\\kappa\\sigma z_{\\mathrm{TP}}}{T_0} \\right)\\;,
\\end{align*}
```

where ``p_0``, ``T_0``, ``z_{\\mathrm{TP}}``, ``\\Gamma_{\\mathrm{TS}}``, ``\\Gamma_{\\mathrm{SS}}``, ``\\sigma``, ``\\gamma`` and ``\\kappa`` are given by `namelists.atmosphere.ground_pressure`, `namelists.atmosphere.temperature`, `namelists.atmosphere.tropopause_height`, `namelists.atmosphere.troposphere_lapse_rate`, `namelists.atmosphere.stratosphere_lapse_rate`, `constants.sig`, `constants.gamma` and `constants.kappa`, respectively.

# Fields

  - `pbar::A`: Mass-weighted potential temperature ``P \\left(z\\right)`` (``P \\left(x, y, z, t\\right)`` in compressible mode).

  - `thetabar::A`: Background potential temperature ``\\overline{\\theta} \\left(z\\right)``.

  - `rhobar::A`: Background density ``\\overline{\\rho} \\left(z\\right)``.

  - `n2::A`: Squared buoyancy frequency ``N^2 \\left(z\\right)``.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `grid`: Collection of parameters and fields that describe the grid.

  - `model`: Dynamic equations.

  - `background`: Atmospheric background.

# See also

  - [`PinCFlow.Types.FoundationalTypes.compute_n2!`](@ref)
"""
struct Atmosphere{A <: AbstractArray{<:AbstractFloat, 3}}
    pbar::A
    thetabar::A
    rhobar::A
    n2::A
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
    rhobar = ones(nxx, nyy, nzz)
    thetabar = potential_temperature ./ thetaref .* ones(nxx, nyy, nzz)
    pbar = rhobar .* thetabar
    n2 = zeros(nxx, nyy, nzz)

    return Atmosphere(pbar, thetabar, rhobar, n2)
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
    rhobar = ones(nxx, nyy, nzz)
    thetabar = potential_temperature ./ thetaref .* ones(nxx, nyy, nzz)
    pbar = rhobar .* thetabar
    n2 = (buoyancy_frequency .* tref) .^ 2 .* ones(nxx, nyy, nzz)

    return Atmosphere(pbar, thetabar, rhobar, n2)
end

function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Union{PseudoIncompressible, Compressible},
    background::Isothermal,
)::Atmosphere
    (; nbz) = namelists.domain
    (; temperature, ground_pressure) = namelists.atmosphere
    (; thetaref, pref, kappa, sig, gamma, g_ndim) = constants
    (; zz_size, nxx, nyy, nzz, ko, k0, k1) = domain
    (; zc, jac, dz) = grid

    # Initialize the background fields.
    (pbar, thetabar, rhobar, n2) = (zeros(nxx, nyy, nzz) for i in 1:4)

    t0 = temperature / thetaref
    p0 = ground_pressure / pref

    # Compute the background fields.
    pbar .= p0 .* exp.(.-sig .* zc ./ gamma ./ t0)
    thetabar .= t0 .* exp.(kappa .* sig ./ t0 .* zc)
    rhobar .= pbar ./ thetabar

    compute_n2!(namelists, constants, domain, grid, thetabar, n2)

    return Atmosphere(pbar, thetabar, rhobar, n2)
end

function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Union{PseudoIncompressible, Compressible},
    background::Isentropic,
)::Atmosphere
    (; lz) = namelists.domain
    (; potential_temperature, ground_pressure) = namelists.atmosphere
    (; thetaref, pref, kappa, sig, rsp, gamma, g) = constants
    (; nxx, nyy, nzz) = domain
    (; zc) = grid

    # Initialize the background fields.
    (pbar, thetabar, rhobar, n2) = (zeros(nxx, nyy, nzz) for i in 1:4)

    min_potential_temperature = kappa * g / rsp * lz
    if potential_temperature < min_potential_temperature
        error(
            "Potential temperature too low for given configuration: potential_temperature = ",
            potential_temperature,
            "K",
            " < minimum potential_temperature = ",
            min_potential_temperature,
            "K",
        )
    end

    theta0 = potential_temperature / thetaref
    p0 = ground_pressure / pref

    # Compute the background fields.
    n2 .= 0.0
    thetabar .= theta0
    pbar .= p0 .* (1.0 .- kappa .* sig ./ theta0 .* zc) .^ (1.0 ./ (gamma .- 1.0))
    rhobar .= pbar ./ thetabar

    return Atmosphere(pbar, thetabar, rhobar, n2)
end

function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Union{PseudoIncompressible, Compressible},
    background::Realistic,
)::Atmosphere
    (; nbx, nby, nbz) = namelists.domain
    (; potential_temperature, ground_pressure, tropopause_height) =
        namelists.atmosphere
    (; thetaref, lref, pref, kappa, sig, rsp, gamma, g, g_ndim) = constants
    (; zz_size, nxx, nyy, nzz, ko, k0, k1, j0, j1, i0, i1) = domain
    (; zc, jac, dz) = grid

    # Initialize the background fields.
    (pbar, thetabar, rhobar, n2) = (zeros(nxx, nyy, nzz) for i in 1:4)

    p0 = ground_pressure / pref
    ztrop = tropopause_height / lref
    theta0 = potential_temperature / thetaref
    ptrop = p0 * (1.0 - kappa * sig / theta0 * ztrop)^(1.0 / (gamma - 1.0))
    ttrop = theta0 * (ptrop / p0)^kappa

    min_potential_temperature = kappa * g / rsp * tropopause_height
    if potential_temperature < min_potential_temperature
        error(
            "Potential temperature in isentropic troposphere too low for given configuration: potential_temperature = ",
            potenential_temperature,
            "K",
            " < minimum potenential_temperature = ",
            min_potential_temperature,
            "K",
        )
    end

    @ivy for k in 1:nzz, j in 1:nyy, i in 1:nxx
        if zc[i, j, k] <= ztrop
            thetabar[i, j, k] = theta0
            pbar[i, j, k] =
                p0 *
                (1.0 - kappa * sig / theta0 * zc[i, j, k])^(1.0 / (gamma - 1.0))
        else
            thetabar[i, j, k] =
                theta0 * exp(kappa * sig / ttrop * (zc[i, j, k] - ztrop))
            pbar[i, j, k] =
                p0^kappa *
                ptrop^gammainv *
                exp(-sig * gammainv / ttrop * (zc[i, j, k] - ztrop))
        end
    end
    rhobar .= pbar ./ thetabar

    compute_n2!(namelists, constants, domain, grid, thetabar, n2)

    return Atmosphere(pbar, thetabar, rhobar, n2)
end

function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
    model::Union{PseudoIncompressible, Compressible},
    background::LapseRates,
)::Atmosphere
    (; nbx, nby, nbz) = namelists.domain
    (;
        ground_pressure,
        troposphere_lapse_rate,
        stratosphere_lapse_rate,
        tropopause_height,
        temperature,
    ) = namelists.atmosphere
    (; thetaref, lref, pref, kappa, sig, rsp, g, g_ndim, gamma) = constants
    (; zz_size, nxx, nyy, nzz, ko, k0, k1, j0, j1, i0, i1) = domain
    (; zc, jac, dz) = grid

    # Initialize the background fields.
    (pbar, thetabar, rhobar, n2) = (zeros(nxx, nyy, nzz) for i in 1:4)

    gamma_t = troposphere_lapse_rate / thetaref * lref
    gamma_s = stratosphere_lapse_rate / thetaref * lref

    p0 = ground_pressure / pref
    t0 = temperature / thetaref
    ztrop = tropopause_height / lref

    if gamma_t != 0.0
        power_t = g / (rsp * troposphere_lapse_rate)
        ptrop = p0 * (1.0 - gamma_t * ztrop / t0)^power_t
    else
        ptrop = p0 * exp(-ztrop * sig / gamma / t0)
        pttrop = t0 * exp(kappa * sig / t0 * ztrop)
    end
    if gamma_s != 0.0
        power_s = g / (rsp * stratosphere_lapse_rate)
    end

    ttrop = t0 - gamma_t * ztrop

    @ivy for k in 1:nzz, j in 1:nyy, i in 1:nxx
        if zc[i, j, k] <= ztrop
            tbar = t0 - gamma_t * zc[i, j, k]

            if gamma_t != 0.0
                pbar[i, j, k] = p0 * (1.0 - gamma_t * zc[i, j, k] / t0)^power_t
                thetabar[i, j, k] = tbar * (p0 / pbar[i, j, k])^(1 / power_t)
            else
                pbar[i, j, k] = p0 * exp(-zc[i, j, k] * sig / gamma / t0)
                thetabar[i, j, k] = t0 * exp(kappa * sig / t0 * zc[i, j, k])
            end

        else
            tbar = ttrop - gamma_s * (zc[i, j, k] - ztrop)

            if gamma_s != 0.0
                pbar[i, j, k] =
                    ptrop *
                    (1.0 - gamma_s * (zc[i, j, k] - ztrop) / ttrop)^power_s
                thetabar[i, j, k] = tbar * (ptrop / pbar[i, j, k])^(1 / power_s)
            else
                pbar[i, j, k] =
                    ptrop * exp(-(zc[i, j, k] - ztrop) * sig / gamma / ttrop)
                thetabar[i, j, k] =
                    pttrop * exp(kappa * sig / ttrop * (zc[i, j, k] - ztrop))
            end
        end
    end

    rhobar .= pbar ./ thetabar

    compute_n2!(namelists, constants, domain, grid, thetabar, n2)

    return Atmosphere(pbar, thetabar, rhobar, n2)
end
