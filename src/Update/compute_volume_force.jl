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
    testcase::AbstractTestCase,
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
    testcase::AbstractWKBTestCase,
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
    testcase::AbstractWKBTestCase,
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
    testcase::AbstractWKBTestCase,
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
    testcase::AbstractTestCase,
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
    testcase::AbstractWKBTestCase,
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
    testcase::AbstractWKBTestCase
)::AbstractFloat
```

Return the tracer flux convergence due to gravity waves.

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `variable`: Variable (equation) of choice.

  - `testcase`: Test case on which the current simulation is based.

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
    (; testcase) = state.namelists.setting

    return compute_volume_force(state, i, j, k, variable, testcase)
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::AbstractVariable,
    testcase::AbstractTestCase,
)::AbstractFloat
    return 0.0
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    testcase::AbstractWKBTestCase,
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
    testcase::AbstractWKBTestCase,
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
    testcase::AbstractWKBTestCase,
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
    testcase::AbstractTestCase,
)::AbstractFloat
    return conductive_heating(state, i, j, k)
end

function compute_volume_force(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::P,
    testcase::AbstractWKBTestCase,
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
    testcase::AbstractWKBTestCase,
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
