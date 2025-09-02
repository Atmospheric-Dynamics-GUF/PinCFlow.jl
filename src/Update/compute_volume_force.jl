"""
```julia
compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)::AbstractFloat
```

Return the volume force in the equation specified by `variable`.

```julia
compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
    testcase::AbstractTestCase,
)::AbstractFloat
```

Return ``0`` as the volume force in non-WKB test cases.

```julia
compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::U,
    testcase::AbstractWKBTestCase,
)::AbstractFloat
```

Return the gravity-wave drag on the zonal momentum, interpolated to ``\\left(i + 1 / 2, j, k\\right)``.

```julia
compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::V,
    testcase::AbstractWKBTestCase,
)::AbstractFloat
```

Return the gravity-wave drag on the meridional momentum, interpolated to ``\\left(i, j + 1 / 2, k\\right)``.

```julia
compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
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
    indices::NTuple{3, <:Integer},
    variable::P,
    testcase::AbstractWKBTestCase,
)::AbstractFloat
```

Return the gravity-wave impact on the mass-weighted potential temperature (diabatic heating).

# Arguments

  - `state`: Model state.

  - `indices`: Grid-cell indices.

  - `variable`: Variable (equation) of choice.

  - `testcase`: Test case on which the current simulation is based.
"""
function compute_volume_force end

function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)::AbstractFloat
    (; testcase) = state.namelists.setting

    force = compute_volume_force(state, indices, variable, testcase) +
        conductive_heating(state, indices)

    return force
end

function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
    testcase::AbstractTestCase,
)::AbstractFloat
    return 0.0
end

function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::U,
    testcase::AbstractWKBTestCase,
)::AbstractFloat
    (; dudt) = state.wkb.tendencies
    (ix, jy, kz) = indices

    return (dudt[ix, jy, kz] + dudt[ix + 1, jy, kz]) / 2
end

function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::V,
    testcase::AbstractWKBTestCase,
)::AbstractFloat
    (; dvdt) = state.wkb.tendencies
    (ix, jy, kz) = indices

    return (dvdt[ix, jy, kz] + dvdt[ix, jy + 1, kz]) / 2
end

function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    testcase::AbstractWKBTestCase,
)::AbstractFloat
    (; jac, met) = state.grid
    (; dudt, dvdt) = state.wkb.tendencies
    (ix, jy, kz) = indices

    return (
        jac[ix, jy, kz + 1] * (
            met[ix, jy, kz, 1, 3] * dudt[ix, jy, kz] +
            met[ix, jy, kz, 2, 3] * dvdt[ix, jy, kz]
        ) +
        jac[ix, jy, kz] * (
            met[ix, jy, kz + 1, 1, 3] * dudt[ix, jy, kz + 1] +
            met[ix, jy, kz + 1, 2, 3] * dvdt[ix, jy, kz + 1]
        )
    ) / (jac[ix, jy, kz] + jac[ix, jy, kz + 1])
end

function compute_volume_force(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::P,
    testcase::AbstractWKBTestCase,
)::AbstractFloat
    (; dthetadt) = state.wkb.tendencies
    (ix, jy, kz) = indices

    return dthetadt[ix, jy, kz]
end
