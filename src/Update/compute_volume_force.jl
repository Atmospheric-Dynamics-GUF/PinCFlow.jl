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
\\left(\\frac{\\partial \\widehat{w}}{\\partial t}\\right)_\\mathrm{w} = \\left[G^{1 3} \\left(\\frac{\\partial u}{\\partial t}\\right)_\\mathrm{w}\\right]_{k + 1 / 2} + \\left[G^{2 3} \\left(\\frac{\\partial v}{\\partial t}\\right)_\\mathrm{w}\\right]_{k + 1 / 2}.
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

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `variable`: Variable (equation) of choice.

  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.Update.conductive_heating`](@ref)
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
    (; leading_order_impact) = state.namelists.tracer
    (; chiq0) = state.tracer.tracerforcings
    (; model) = state.namelists.atmosphere

    @ivy if leading_order_impact && model == Compressible()
        return chiq0.dchidt[i, j, k]
    else
        return 0.0
    end
end
