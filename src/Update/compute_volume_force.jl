"""
```julia
compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::AbstractVariable,
)::AbstractFloat
```

Return the volume force in the equation specified by `variable`.

```julia
compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::AbstractVariable,
    test_case::AbstractTestCase,
)::AbstractFloat
```

Return ``0`` as the volume force in non-WKB test cases (for all variables except the mass-weighted potential temperature).

```julia
compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    test_case::AbstractWKBTestCase,
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
    test_case::AbstractWKBTestCase,
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
    test_case::AbstractWKBTestCase,
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
    test_case::AbstractTestCase,
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
    test_case::AbstractWKBTestCase,
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
    test_case::AbstractWKBTestCase,
)::AbstractFloat
```

Return the tracer flux convergence due to gravity waves.

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `variable`: Variable (equation) of choice.

  - `test_case`: Test case on which the current simulation is based.

# See also

  - [`PinCFlow.Update.conductive_heating`](@ref)
"""
function compute_volume_force end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::AbstractVariable,
)::AbstractFloat
    (; test_case) = state.namelists.setting

    return compute_volume_force(state, i, j, k, variable, test_case)
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::AbstractVariable,
    test_case::AbstractTestCase,
)::AbstractFloat
    return 0.0
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    test_case::AbstractWKBTestCase,
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
    test_case::AbstractWKBTestCase,
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
    test_case::AbstractWKBTestCase,
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
    test_case::AbstractTestCase,
)::AbstractFloat
    return conductive_heating(state, i, j, k)
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::P,
    test_case::AbstractWKBTestCase,
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
    test_case::AbstractWKBTestCase,
)::AbstractFloat
    (; leading_order_impact) = state.namelists.tracer
    (; chiq0) = state.tracer.tracerforcings
    (; model) = state.namelists.setting

    @ivy if leading_order_impact && model == Compressible()
        return chiq0.dchidt[i, j, k]
    else
        return 0.0
    end
end

function compute_volume_force(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
    variables::TKE,
)::AbstractFloat
    (; km, kh, shearproduction, buoyancyproduction) =
        state.turbulence.turbulenceauxiliaries
    (; rho) = p0
    (; rhobar, n2) = state.atmosphere
    (; jac, dz) = state.grid
    (; g_ndim, lref, tref) = state.constants
    (; ko, k0)

    shear =
        km[i, j, k] * (
            compute_momentum_diffusion_terms(state, p0, i, j, k, U(), Z())^2.0 +
            compute_momentum_diffusion_terms(state, p0, i, j, k, V(), Z())^2.0
        )

    shearproduction[i, j, k] = shear

    bu = g_ndim * rho[i, j, k + 1] / (rho[i, j, k + 1] + rhobar[i, j, k + 1])
    bd = g_ndim * rho[i, j, k - 1] / (rho[i, j, k - 1] + rhobar[i, j, k - 1])

    buoyancy =
        -kh[i, j, k] * (n2[i, j, k]  + (bu - bd) / (jac[i, j, k] * 2.0 * dz))

    buoyancyproduction[i, j, k] = buoyancy

    return (rho[i, j, k] + rhobar[i, j, k]) * (shear + buoyancy)
end