"""
```julia
compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
```

Compute the factor by which the wind should be multiplied at ``\\left(i + 1 / 2, j, k\\right)``, ``\\left(i, j + 1 / 2, k\\right)`` or ``\\left(i, j, k + 1 / 2\\right)`` by dispatching to a method specific for the dynamic equations and the component indicated by `variable`.

In compressible mode, the Euler steps that are used to integrate the right-hand side of the momentum equation update ``\\left(J P\\right)_{i + 1 / 2} u_{i + 1 / 2}``, ``\\left(J P\\right)_{j + 1 / 2} v_{j + 1 / 2}`` and ``\\left(J P\\right)_{k + 1 / 2} \\widehat{w}_{k + 1 / 2}`` instead of ``u_{i + 1 / 2}``, ``v_{j + 1 / 2}`` and ``\\widehat{w}_{k + 1 / 2}``.

```julia
compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
    model::AbstractModel,
)
```

Return ``1`` as the factor by which the wind should be multiplied in non-compressible mode.

```julia
compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::U,
    model::Compressible,
)
```

Return ``\\left(J P\\right)_{i + 1 / 2}`` as the factor by which the zonal wind should be multiplied in compressible mode.

```julia
compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::V,
    model::Compressible,
)
```

Return ``\\left(J P\\right)_{j + 1 / 2}`` as the factor by which the meridional wind should be multiplied in compressible mode.

```julia
compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::Compressible,
)
```

Return ``\\left(J P\\right)_{k + 1 / 2}`` as the factor by which the transformed vertical wind should be multiplied in compressible mode.

The interpolation is given by

```math
\\left(J P\\right)_{k + 1 / 2} = \\frac{J J_{k + 1} \\left(P + P_{k + 1}\\right)}{J + J_{k + 1}}.
```

# Arguments

  - `state`: Model state.

  - `indices`: Grid-cell indices.

  - `variable`: Variable for which the factor is needed.

  - `model`: Dynamic equations.

# Returns

  - `::AbstractFloat`: Wind factor.
"""
function compute_compressible_wind_factor end

function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
    (; model) = state.namelists.setting
    return compute_compressible_wind_factor(state, indices, variable, model)
end

function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
    model::AbstractModel,
)
    return 1.0
end

function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::U,
    model::Compressible,
)
    (; jac) = state.grid
    (; p) = state.variables.predictands
    (ix, jy, kz) = indices
    return (
        jac[ix, jy, kz] * p[ix, jy, kz] +
        jac[ix + 1, jy, kz] * p[ix + 1, jy, kz]
    ) / 2
end

function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::V,
    model::Compressible,
)
    (; jac) = state.grid
    (; p) = state.variables.predictands
    (ix, jy, kz) = indices
    return (
        jac[ix, jy, kz] * p[ix, jy, kz] +
        jac[ix, jy + 1, kz] * p[ix, jy + 1, kz]
    ) / 2
end

function compute_compressible_wind_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::Compressible,
)
    (; jac) = state.grid
    (; p) = state.variables.predictands
    (ix, jy, kz) = indices
    return jac[ix, jy, kz] *
           jac[ix, jy, kz + 1] *
           (p[ix, jy, kz] + p[ix, jy, kz + 1]) /
           (jac[ix, jy, kz] + jac[ix, jy, kz + 1])
end
