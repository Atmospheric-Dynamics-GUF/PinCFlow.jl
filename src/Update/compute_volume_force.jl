"""
```julia
compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Union{U, V, W, P, Chi},
)::AbstractFloat
```

Return the volume force in the equation specified by `variable`, by dispatching to an equation-and-WKB-mode specific method.

```julia
compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Union{U, V, W, Chi},
    wkb_mode::Val{:NoWKB},
)::AbstractFloat
```

Return ``0`` as the volume force in non-WKB configurations (for all variables except the mass-weighted potential temperature).

```julia
compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::AbstractFloat
```

Return the gravity-wave drag on the zonal momentum, interpolated to ``\\left(i + 1 / 2, j, k\\right)``.

```julia
compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::AbstractFloat
```

Return the gravity-wave drag on the meridional momentum, interpolated to ``\\left(i, j + 1 / 2, k\\right)``.

```julia
compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::AbstractFloat
```

Return the gravity-wave drag on the transformed vertical momentum, interpolated to ``\\left(i, j, k + 1 / 2\\right)``, as given by

```math
\\left(\\frac{\\partial \\hat{w}}{\\partial t}\\right)_\\mathrm{w} = \\left[G^{1 3} \\left(\\frac{\\partial u}{\\partial t}\\right)_\\mathrm{w}\\right]_{k + 1 / 2} + \\left[G^{2 3} \\left(\\frac{\\partial v}{\\partial t}\\right)_\\mathrm{w}\\right]_{k + 1 / 2}.
```

```julia
compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::P,
    wkb_mode::Val{:NoWKB},
)::AbstractFloat
```

Return the conductive heating.

```julia
compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::P,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::AbstractFloat
```

Return the sum of gravity-wave impact on the mass-weighted potential temperature and conductive heating.

```julia
compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::Chi,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::AbstractFloat
```

Return the tracer flux convergence due to gravity waves.

```julia
compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::TKE,
)::AbstractFloat
```

Return the mass-weighted impact of shear ``\\mathcal{S}`` and buoyancy ``\\mathcal{B}`` on the TKE, given by

```math
\\left(\\frac{\\partial \\rho e_\\mathrm{k}}{\\partial t}\\right) = \\rho\\mathcal{S} + \\rho\\mathcal{B}
```

where

```math
\\begin{align*}
\\mathcal{S} &= K_\\mathrm{M}\\left[\\left(\\frac{\\partial u}{\\partial \\hat{z}}\\right)^2 + \\left(\\frac{\\partial v}{\\partial \\hat{z}}\\right)^2\\right] \\;, \\\\
\\mathcal{B} &= -K_\\mathrm{H}\\left(N^2 + \\frac{\\partial b}{\\partial \\hat{z}}\\right) \\;,
\\end{align*}
```

and ``K_\\mathrm{M}`` and ``K_\\mathrm{H}`` represent the eddy diffusion coefficients for momentum and heat, respectively. The buoyancy term is calculated by dispatching to the model-specific method.

```julia
compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::TKE,
    model::Union{Val{:PseudoIncompressible}, Val{:Compressible}},
)::AbstractFloat
```

Return the buoyancy forcing on the TKE for configurations in pseudo-incompressible and compressible mode.

```julia
compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::TKE,
    model::Val{:Boussinesq},
)::AbstractFloat
```

Return the buoyancy forcing on the TKE for configurations in Boussinesq mode.

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `variable`: Variable (equation) of choice.

  - `wkb_mode`: Approximations used by MS-GWaM.

  - `model`: Dynamic equations.

# See also

  - [`PinCFlow.Update.conductive_heating`](@ref)

  - [`PinCFlow.Update.compute_momentum_diffusion_terms`](@ref)
"""
function compute_volume_force end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Union{U, V, W, P, Chi},
)::AbstractFloat
    (; wkb_mode) = state.namelists.wkb

    @dispatch_wkb_mode return compute_volume_force(
        state,
        i,
        j,
        k,
        variable,
        Val(wkb_mode),
    )
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Union{U, V, W, Chi},
    wkb_mode::Val{:NoWKB},
)::AbstractFloat
    return 0.0
end

@ivy function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::AbstractFloat
    (; dudt) = state.wkb.tendencies

    return (dudt[i, j, k] + dudt[i + 1, j, k]) / 2
end

@ivy function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::AbstractFloat
    (; dvdt) = state.wkb.tendencies

    return (dvdt[i, j, k] + dvdt[i, j + 1, k]) / 2
end

@ivy function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::AbstractFloat
    (; jac, met) = state.grid
    (; dudt, dvdt) = state.wkb.tendencies

    return (
        jac[i, j, k + 1] * (
            met[i, j, k, 1, 3] * dudt[i, j, k] +
            met[i, j, k, 2, 3] * dvdt[i, j, k]
        ) +
        jac[i, j, k] * (
            met[i, j, k + 1, 1, 3] * dudt[i, j, k + 1] +
            met[i, j, k + 1, 2, 3] * dvdt[i, j, k + 1]
        )
    ) / (jac[i, j, k] + jac[i, j, k + 1])
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::P,
    wkb_mode::Val{:NoWKB},
)::AbstractFloat
    return conductive_heating(state, i, j, k)
end

@ivy function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::P,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::AbstractFloat
    (; dthetadt) = state.wkb.tendencies

    return dthetadt[i, j, k] + conductive_heating(state, i, j, k)
end

@ivy function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::Chi,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)::AbstractFloat
    (; leading_order_impact) = state.namelists.tracer
    (; dchidt0) = state.tracer.tracerwkbtendencies
    (; model) = state.namelists.atmosphere

    impact = 0.0

    if leading_order_impact && model == :Compressible
        impact += dchidt0[i, j, k]
    end
    return impact
end

@ivy function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::TKE,
)::AbstractFloat
    (; shear_production, buoyancy_production) =
        state.turbulence.turbulenceauxiliaries
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere
    (; model) = state.namelists.atmosphere

    shear =
        turbulence_diffusion_coefficient(state, i, j, k, KM()) * (
            compute_momentum_diffusion_terms(state, i, j, k, U(), Z())^2.0 +
            compute_momentum_diffusion_terms(state, i, j, k, V(), Z())^2.0
        )

    shear_production[i, j, k] = shear

    @dispatch_model buoyancy =
        compute_volume_force(state, i, j, k, variable, Val(model))

    buoyancy_production[i, j, k] = buoyancy

    return (rho[i, j, k] + rhobar[i, j, k]) * (shear + buoyancy)
end

@ivy function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::TKE,
    model::Union{Val{:PseudoIncompressible}, Val{:Compressible}},
)::AbstractFloat
    (; rho) = state.variables.predictands
    (; rhobar, n2) = state.atmosphere
    (; jac, dz) = state.grid
    (; g_ndim) = state.constants

    bu = g_ndim * (1 / (rho[i, j, k + 1] / rhobar[i, j, k + 1] + 1) - 1)
    bd = g_ndim * (1 / (rho[i, j, k - 1] / rhobar[i, j, k - 1] + 1) - 1)

    buoyancy =
        -turbulence_diffusion_coefficient(state, i, j, k, KH()) *
        (n2[i, j, k] + (bu - bd) / (jac[i, j, k] * 2.0 * dz))

    return buoyancy
end

@ivy function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::TKE,
    model::Val{:Boussinesq},
)::AbstractFloat
    (; rhop) = state.variables.predictands
    (; rhobar, n2) = state.atmosphere
    (; jac, dz) = state.grid
    (; g_ndim) = state.constants

    bu = g_ndim * (1 / (rhop[i, j, k + 1] / rhobar[i, j, k + 1] + 1) - 1)
    bd = g_ndim * (1 / (rhop[i, j, k - 1] / rhobar[i, j, k - 1] + 1) - 1)

    buoyancy =
        -turbulence_diffusion_coefficient(state, i, j, k, KH()) *
        (n2[i, j, k] + (bu - bd) / (jac[i, j, k] * 2.0 * dz))

    return buoyancy
end
