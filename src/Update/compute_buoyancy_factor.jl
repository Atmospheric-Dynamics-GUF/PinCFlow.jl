"""
```julia
compute_buoyancy_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Union{RhoP, W},
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
    i::Integer,
    j::Integer,
    k::Integer,
    variable::RhoP,
    model::Compressible,
)::AbstractFloat
```

Return ``f_{\\rho'} = P / \\left(\\rho \\overline{\\theta}\\right)`` as the factor by which the buoyancy term should be multiplied at ``\\left(i, j, k\\right)`` in compressible mode.

```julia
compute_buoyancy_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::RhoP,
    model::Union{Boussinesq, PseudoIncompressible},
)::AbstractFloat
```

Return ``f_{\\rho'} = \\overline{\\rho} / \\rho`` as the factor by which the buoyancy term should be multiplied at ``\\left(i, j, k\\right)`` in pseudo-incompressible mode (this method is also used in Boussinesq mode).

```julia
compute_buoyancy_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    model::Compressible,
)::AbstractFloat
```

Return ``f_w = \\left(P / \\overline{\\theta}\\right)_{k + 1 / 2} / \\rho_{k + 1 / 2}`` as the factor by which the buoyancy term should be multiplied in compressible mode.

```julia
compute_buoyancy_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    model::Union{Boussinesq, PseudoIncompressible},
)::AbstractFloat
```

Return ``f_w = \\overline{\\rho}_{k + 1 / 2} / \\rho_{k + 1 / 2}`` as the factor by which the buoyancy term should be multiplied at ``\\left(i, j, k + 1 / 2\\right)`` in pseudo-incompressible mode (this method is also used in Boussinesq mode).

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `variable`: Variable for which the factor is needed.

  - `model`: Dynamic equations.
"""
function compute_buoyancy_factor end

function compute_buoyancy_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::Union{RhoP, W},
)::AbstractFloat
    (; model) = state.namelists.setting
    return compute_buoyancy_factor(state, i, j, k, variable, model)
end

function compute_buoyancy_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::RhoP,
    model::Compressible,
)::AbstractFloat
    (; rhobar, thetabar, pbar) = state.atmosphere
    (; rho) = state.variables.predictands
    @ivy return pbar[i, j, k] / thetabar[i, j, k] /
                (rho[i, j, k] + rhobar[i, j, k])
end

function compute_buoyancy_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::RhoP,
    model::Union{Boussinesq, PseudoIncompressible},
)::AbstractFloat
    (; rhobar) = state.atmosphere
    (; rho) = state.variables.predictands
    @ivy return rhobar[i, j, k] / (rho[i, j, k] + rhobar[i, j, k])
end

function compute_buoyancy_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    model::Compressible,
)::AbstractFloat
    (; jac) = state.grid
    (; rhobar, thetabar, pbar) = state.atmosphere
    (; rho) = state.variables.predictands
    @ivy return (
        jac[i, j, k + 1] * pbar[i, j, k] / thetabar[i, j, k] +
        jac[i, j, k] * pbar[i, j, k + 1] / thetabar[i, j, k + 1]
    ) / (
        jac[i, j, k + 1] * (rho[i, j, k] + rhobar[i, j, k]) +
        jac[i, j, k] * (rho[i, j, k + 1] + rhobar[i, j, k + 1])
    )
end

function compute_buoyancy_factor(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    model::Union{Boussinesq, PseudoIncompressible},
)::AbstractFloat
    (; jac) = state.grid
    (; rhobar) = state.atmosphere
    (; rho) = state.variables.predictands
    @ivy return (
        jac[i, j, k + 1] * rhobar[i, j, k] + jac[i, j, k] * rhobar[i, j, k + 1]
    ) / (
        jac[i, j, k + 1] * (rho[i, j, k] + rhobar[i, j, k]) +
        jac[i, j, k] * (rho[i, j, k + 1] + rhobar[i, j, k + 1])
    )
end
