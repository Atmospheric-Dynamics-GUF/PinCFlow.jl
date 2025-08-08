"""
```julia
compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
```

Compute the factor by which the buoyancy term should be multiplied at ``\\left(i, j, k\\right)`` or ``\\left(i, j, k + 1 / 2\\right)`` by dispatching to a method specific for the dynamic equations and `variable`.

In pseudo-incompressible mode, the squared buoyancy frequency used by PinCFlow is

```math
N^2 = \\frac{g}{\\overline{\\theta}} \\frac{\\mathrm{d} \\overline{\\theta}}{\\mathrm{d} z},
```

whereas in compressible mode, it is

```math
N^2 = \\frac{g P}{\\rho \\overline{\\theta}^2} \\frac{\\mathrm{d} \\overline{\\theta}}{\\mathrm{d} z}.
```

In both modes, the buoyancy term is expressed in terms of ``N^2``. Thus, one has

```math
\\left(\\frac{\\partial b'}{\\partial t}\\right)_{N^2} = f_{b'} N^2 w,
```

with ``f_{b'} = 1`` in compressible mode and ``f_{b'} = \\overline{\\rho} / \\rho`` in pseudo-incompressible mode. This method returns either ``f_{b'}`` at ``\\left(i, j, k\\right)`` or ``f_w``, which is the interpolation of ``f_{b'}`` to ``\\left(i, j, k + 1 / 2\\right)``, based on the type of `variable`. Note that both values of ``f_{b'}`` are equivalent in Boussinesq mode, where ``\\rho = \\overline{\\rho} = \\rho_0`` (since the density fluctuations are treated separately).

```julia
compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::Compressible,
)
```

Return ``f_{b'} = 1`` as the factor by which the buoyancy term should be multiplied in compressible mode.

```julia
compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::AbstractModel,
)
```

Return ``f_{b'} = \\overline{\\rho} / \\rho`` as the factor by which the buoyancy term should be multiplied at ``\\left(i, j, k\\right)`` in pseudo-incompressible mode (this method is also used in Boussinesq mode).

```julia
compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::Compressible,
)
```

Return ``f_w = 1`` as the factor by which the buoyancy term should be multiplied in compressible mode.

```julia
compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::AbstractModel,
)
```

Return ``f_w = \\overline{\\rho}_{k + 1 / 2} / \\rho_{k + 1 / 2}`` as the factor by which the buoyancy term should be multiplied at ``\\left(i, j, k + 1 / 2\\right)`` in pseudo-incompressible mode (this method is also used in Boussinesq mode).

The interpolation to ``\\left(i, j, k + 1 / 2\\right)`` follows

```math
f_w = \\frac{\\overline{\\rho}_{k + 1 / 2}}{\\rho_{k + 1 / 2}} = \\frac{J_{k + 1} \\overline{\\rho} + J \\overline{\\rho}_{k + 1}}{J_{k + 1} \\rho + J \\rho_{k + 1}}.
```

# Arguments

- `state`: Model state.
- `indices`: Grid-cell indices.
- `variable`: Variable for which the factor is needed.
- `model`: Dynamic equations.

# Returns

- `::AbstractFloat`: Buoyancy factor.
"""
function compute_compressible_buoyancy_factor end

function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)
    (; model) = state.namelists.setting
    return compute_compressible_buoyancy_factor(state, indices, variable, model)
end

function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::Compressible,
)
    return 1.0
end

function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::AbstractModel,
)
    (; rhostrattfc) = state.atmosphere
    (; rho) = state.variables.predictands
    (ix, jy, kz) = indices
    return rhostrattfc[ix, jy, kz] / (rho[ix, jy, kz] + rhostrattfc[ix, jy, kz])
end

function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::Compressible,
)
    return 1.0
end

function compute_compressible_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::AbstractModel,
)
    (; jac) = state.grid
    (; rhostrattfc) = state.atmosphere
    (; rho) = state.variables.predictands
    (ix, jy, kz) = indices
    return (
        jac[ix, jy, kz + 1] * rhostrattfc[ix, jy, kz] +
        jac[ix, jy, kz] * rhostrattfc[ix, jy, kz + 1]
    ) / (
        jac[ix, jy, kz + 1] * (rho[ix, jy, kz] + rhostrattfc[ix, jy, kz]) +
        jac[ix, jy, kz] * (rho[ix, jy, kz + 1] + rhostrattfc[ix, jy, kz + 1])
    )
end
