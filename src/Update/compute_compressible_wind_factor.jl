"""
```julia
compute_compressible_wind_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Union{U, V, W},
)::AbstractFloat
```

Compute the factor by which the wind should be multiplied at ``\\left(i + 1 / 2, j, k\\right)``, ``\\left(i, j + 1 / 2, k\\right)`` or ``\\left(i, j, k + 1 / 2\\right)``, by dispatching to a method specific for the dynamic equations and the component indicated by `variable`, and return the result.

In compressible mode, the Euler steps that are used to integrate the right-hand side of the momentum equation update ``\\left(J P\\right)_{i + 1 / 2} u_{i + 1 / 2}``, ``\\left(J P\\right)_{j + 1 / 2} v_{j + 1 / 2}`` and ``\\left(J P\\right)_{k + 1 / 2} \\widehat{w}_{k + 1 / 2}`` instead of ``u_{i + 1 / 2}``, ``v_{j + 1 / 2}`` and ``\\widehat{w}_{k + 1 / 2}``.

```julia
compute_compressible_wind_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Union{U, V, W},
    model::Union{Boussinesq, PseudoIncompressible},
)::AbstractFloat
```

Return ``1`` as the factor by which the wind should be multiplied in non-compressible mode.

```julia
compute_compressible_wind_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    model::Compressible,
)::AbstractFloat
```

Return ``\\left(J P\\right)_{i + 1 / 2}`` as the factor by which the zonal wind should be multiplied in compressible mode.

```julia
compute_compressible_wind_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    model::Compressible,
)::AbstractFloat
```

Return ``\\left(J P\\right)_{j + 1 / 2}`` as the factor by which the meridional wind should be multiplied in compressible mode.

```julia
compute_compressible_wind_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    model::Compressible,
)::AbstractFloat
```

Return ``\\left(J P\\right)_{k + 1 / 2}`` as the factor by which the transformed vertical wind should be multiplied in compressible mode.

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `variable`: Variable for which the factor is needed.

  - `model`: Dynamic equations.
"""
function compute_compressible_wind_factor end

function compute_compressible_wind_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Union{U, V, W},
)::AbstractFloat
    (; model) = state.namelists.setting
    return compute_compressible_wind_factor(state, i, j, k, variable, model)
end

function compute_compressible_wind_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Union{U, V, W},
    model::Union{Boussinesq, PseudoIncompressible},
)::AbstractFloat
    return 1.0
end

function compute_compressible_wind_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    model::Compressible,
)::AbstractFloat
    (; jac) = state.grid
    (; p) = state.variables.predictands
    @ivy return (
        jac[i, j, k] * p[i, j, k] + jac[i + 1, j, k] * p[i + 1, j, k]
    ) / 2
end

function compute_compressible_wind_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    model::Compressible,
)::AbstractFloat
    (; jac) = state.grid
    (; p) = state.variables.predictands
    @ivy return (
        jac[i, j, k] * p[i, j, k] + jac[i, j + 1, k] * p[i, j + 1, k]
    ) / 2
end

function compute_compressible_wind_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    model::Compressible,
)::AbstractFloat
    (; jac) = state.grid
    (; p) = state.variables.predictands
    @ivy return jac[i, j, k] *
                jac[i, j, k + 1] *
                (p[i, j, k] + p[i, j, k + 1]) /
                (jac[i, j, k] + jac[i, j, k + 1])
end
