"""
```julia
compute_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)::AbstractFloat
```

Compute the factor by which the buoyancy term should be multiplied at ``\\left(i, j, k\\right)`` or ``\\left(i, j, k + 1 / 2\\right)``, by dispatching to a method specific for the dynamic equations and `variable`, and return the result.

In pseudo-incompressible mode, ``\\rho'`` are deviations of the total density from ``\\overline{\\rho}``, which describes the reference atmosphere. However, in compressible mode, ``\\rho' = \\rho - P / \\overline{\\theta}`` does not reduce to this, i.e. the density background has a spatiotemporal dependence. As a consequence, the right-hand side of the prognostic equation for ``\\rho'`` is given by

```math
\\left(\\frac{\\partial \\rho'}{\\partial t}\\right)_{N^2} = f_{\\rho'} \\frac{N^2 \\rho w}{g},
```

with ``f_{\\rho'} = \\overline{\\rho} / \\rho`` in pseudo-incompressible mode and ``f_{\\rho'} = P / \\left(\\rho \\overline{\\theta}\\right)`` in compressible mode. This method returns either ``f_{\\rho'}`` at ``\\left(i, j, k\\right)`` or ``f_w``, which is the interpolation of ``f_{\\rho'}`` to ``\\left(i, j, k + 1 / 2\\right)``, based on the type of `variable`.

```julia
compute_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::Compressible,
)::AbstractFloat
```

Return ``f_{\\rho'} = P / \\left(\\rho \\overline{\\theta}\\right)`` as the factor by which the buoyancy term should be multiplied at ``\\left(i, j, k\\right)`` in compressible mode.

```julia
compute_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::AbstractModel,
)::AbstractFloat
```

Return ``f_{\\rho'} = \\overline{\\rho} / \\rho`` as the factor by which the buoyancy term should be multiplied at ``\\left(i, j, k\\right)`` in pseudo-incompressible mode (this method is also used in Boussinesq mode).

```julia
compute_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::Compressible,
)::AbstractFloat
```

Return ``f_w = \\left(P / \\overline{\\theta}\\right)_{k + 1 / 2} / \\rho_{k + 1 / 2}`` as the factor by which the buoyancy term should be multiplied in compressible mode.

```julia
compute_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::AbstractModel,
)::AbstractFloat
```

Return ``f_w = \\overline{\\rho}_{k + 1 / 2} / \\rho_{k + 1 / 2}`` as the factor by which the buoyancy term should be multiplied at ``\\left(i, j, k + 1 / 2\\right)`` in pseudo-incompressible mode (this method is also used in Boussinesq mode).

# Arguments

  - `state`: Model state.

  - `indices`: Grid-cell indices ``\\left(i, j, k\\right)``.

  - `variable`: Variable for which the factor is needed.

  - `model`: Dynamic equations.
"""
function compute_buoyancy_factor end

function compute_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::AbstractVariable,
)::AbstractFloat
    (; model) = state.namelists.setting
    return compute_buoyancy_factor(state, indices, variable, model)
end

function compute_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::Compressible,
)::AbstractFloat
    (; rhostrattfc, thetastrattfc, pstrattfc) = state.atmosphere
    (; rho) = state.variables.predictands
    (ix, jy, kz) = indices
    return pstrattfc[ix, jy, kz] / thetastrattfc[ix, jy, kz] /
           (rho[ix, jy, kz] + rhostrattfc[ix, jy, kz])
end

function compute_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::RhoP,
    model::AbstractModel,
)::AbstractFloat
    (; rhostrattfc) = state.atmosphere
    (; rho) = state.variables.predictands
    (ix, jy, kz) = indices
    return rhostrattfc[ix, jy, kz] / (rho[ix, jy, kz] + rhostrattfc[ix, jy, kz])
end

function compute_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::Compressible,
)::AbstractFloat
    (; jac) = state.grid
    (; rhostrattfc, thetastrattfc, pstrattfc) = state.atmosphere
    (; rho) = state.variables.predictands
    (ix, jy, kz) = indices
    return (
        jac[ix, jy, kz + 1] * pstrattfc[ix, jy, kz] /
        thetastrattfc[ix, jy, kz] +
        jac[ix, jy, kz] * pstrattfc[ix, jy, kz + 1] /
        thetastrattfc[ix, jy, kz + 1]
    ) / (
        jac[ix, jy, kz + 1] * (rho[ix, jy, kz] + rhostrattfc[ix, jy, kz]) +
        jac[ix, jy, kz] * (rho[ix, jy, kz + 1] + rhostrattfc[ix, jy, kz + 1])
    )
end

function compute_buoyancy_factor(
    state::State,
    indices::NTuple{3, <:Integer},
    variable::W,
    model::AbstractModel,
)::AbstractFloat
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
