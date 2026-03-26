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
    wkb_mode::NoWKB,
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
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
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
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
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
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
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
    wkb_mode::NoWKB,
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
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
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
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)::AbstractFloat
```

Return the tracer flux convergence due to gravity waves.

```julia
compute_volume_force(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::TKE,
)::AbstractFloat
```

Return the mass-weighted impact of shear and buoyancy on the TKE by dispatching to the appropriate method.

```julia
compute_volume_force(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::TKE,
    model::Union{PseudoIncompressible, Compressible},
)::AbstractFloat
```

Return the mass-weighted impact of shear and buoyancy on the TKE in pseudo-incompressible and compressible mode.

```julia
compute_volume_force(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::TKE,
    model::Boussinesq,
)::AbstractFloat
```

Return the mass-weighted impact of shear and buoyancy on the TKE in Boussinesq mode.

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `variable`: Variable (equation) of choice.

  - `wkb_mode`: Approximations used by MS-GWaM.

  - `p0`: Predictands

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

    return compute_volume_force(state, i, j, k, variable, wkb_mode)
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Union{U, V, W, Chi},
    wkb_mode::NoWKB,
)::AbstractFloat
    return 0.0
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)::AbstractFloat
    (; dudt) = state.wkb.tendencies

    @ivy return (dudt[i, j, k] + dudt[i + 1, j, k]) / 2
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)::AbstractFloat
    (; dvdt) = state.wkb.tendencies

    @ivy return (dvdt[i, j, k] + dvdt[i, j + 1, k]) / 2
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)::AbstractFloat
    (; jac, met) = state.grid
    (; dudt, dvdt) = state.wkb.tendencies

    @ivy return (
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
    wkb_mode::NoWKB,
)::AbstractFloat
    return conductive_heating(state, i, j, k)
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::P,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)::AbstractFloat
    (; dthetadt) = state.wkb.tendencies

    @ivy return dthetadt[i, j, k] + conductive_heating(state, i, j, k)
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::Chi,
    wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
)::AbstractFloat
    (; leading_order_impact, next_order_impact, turbulence_impact) = state.namelists.tracer
    (; dchidt0, dchidt1, dchidtq) = state.tracer.tracerwkbtendencies
    (; model) = state.namelists.atmosphere

    impact = 0.0

    @ivy if leading_order_impact && model == Compressible()
        impact += dchidt0[i, j, k]
    end
    @ivy if next_order_impact
        impact += dchidt1[i, j, k]
    end
    @ivy if turbulence_impact
        impact += dchidtq[i, j, k]
    end
    return impact
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::TKE,
)::AbstractFloat
    (; model) = state.namelists.atmosphere

    return compute_volume_force(state, i, j, k, variables, model)
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::TKE,
    model::Union{PseudoIncompressible, Compressible},
)::AbstractFloat
    (; shearproduction, buoyancyproduction) =
        state.turbulence.turbulenceauxiliaries
    (; km, kh) = state.turbulence.turbulencediffusioncoefficients
    (; rho) = state.variables.predictands
    (; rhobar, n2) = state.atmosphere
    (; jac, dz) = state.grid
    (; g_ndim) = state.constants

    shear =
        km[i, j, k] * (
            compute_momentum_diffusion_terms(state, i, j, k, U(), Z())^2.0 +
            compute_momentum_diffusion_terms(state, i, j, k, V(), Z())^2.0
        )

    shearproduction[i, j, k] = shear

    bu = g_ndim * (1 / (rho[i, j, k + 1] / rhobar[i, j, k + 1] + 1) - 1)
    bd = g_ndim * (1 / (rho[i, j, k - 1] / rhobar[i, j, k - 1] + 1) - 1)

    buoyancy =
        -kh[i, j, k] * (n2[i, j, k] + (bu - bd) / (jac[i, j, k] * 2.0 * dz))

    buoyancyproduction[i, j, k] = buoyancy

    return (rho[i, j, k] + rhobar[i, j, k]) * (shear + buoyancy)
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::TKE,
    model::Boussinesq,
)::AbstractFloat
    (; shearproduction, buoyancyproduction) =
        state.turbulence.turbulenceauxiliaries
    (; km, kh) = state.turbulence.turbulencediffusioncoefficients
    (; rhop) = state.variables.predictands
    (; rhobar, n2) = state.atmosphere
    (; jac, dz) = state.grid
    (; g_ndim) = state.constants

    shear =
        km[i, j, k] * (
            compute_momentum_diffusion_terms(state, i, j, k, U(), Z())^2.0 +
            compute_momentum_diffusion_terms(state, i, j, k, V(), Z())^2.0
        )

    shearproduction[i, j, k] = shear

    bu = g_ndim * (1 / (rhop[i, j, k + 1] / rhobar[i, j, k + 1] + 1) - 1)
    bd = g_ndim * (1 / (rhop[i, j, k - 1] / rhobar[i, j, k - 1] + 1) - 1)

    buoyancy =
        -kh[i, j, k] * (n2[i, j, k] + (bu - bd) / (jac[i, j, k] * 2 * dz))

    buoyancyproduction[i, j, k] = buoyancy

    return (rhop[i, j, k] + rhobar[i, j, k]) * (shear + buoyancy)
end
