"""
```julia
conductive_heating(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
)::AbstractFloat
```

Compute and return the conductive heating by dispatching to specialized methods dependent on the model.

```julia
conductive_heating(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    model::Boussinesq,
)::AbstractFloat
```

Return ``0`` as conductive heating in Boussinesq mode.

```julia
conductive_heating(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    model::PseudoIncompressible,
)::AbstractFloat
```

Return ``0`` as conductive heating in PseudoIncompressible mode.

```julia
conductive_heating(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    model::Compressible,
)::AbstractFloat
```

Compute and return the conductive heating as the convergence of potential temperature fluxes (weighted by the density), i.e.

```math
\\left(\\frac{\\partial P}{\\partial t}\\right)_\\lambda = - \\frac{\\rho}{J} \\left(\\frac{\\mathcal{F}_{i + 1 / 2}^{\\theta, \\widehat{x}} - \\mathcal{F}_{i - 1 / 2}^{\\theta, \\widehat{x}}}{\\Delta \\widehat{x}} + \\frac{\\mathcal{F}_{j + 1 / 2}^{\\theta, \\widehat{y}} - \\mathcal{F}_{j - 1 / 2}^{\\theta, \\widehat{y}}}{\\Delta \\widehat{y}} + \\frac{\\mathcal{F}_{k + 1 / 2}^{\\theta, \\widehat{z}} - \\mathcal{F}_{k - 1 / 2}^{\\theta, \\widehat{z}}}{\\Delta \\widehat{z}}\\right).
```

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `model`: Dynamic equations.
"""
function conductive_heating end

function conductive_heating(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
)::AbstractFloat
    (; model) = state.namelists.atmosphere

    return conductive_heating(state, i, j, k, model)
end

function conductive_heating(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    model::Boussinesq,
)::AbstractFloat
    return 0.0
end

function conductive_heating(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    model::PseudoIncompressible,
)::AbstractFloat
    return 0.0
end

function conductive_heating(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    model::Compressible,
)::AbstractFloat
    (; phitheta) = state.variables.fluxes
    (; rho) = state.variables.predictands
    (; jac, dx, dy, dz) = state.grid
    (; rhobar) = state.atmosphere

    @ivy rhotot = (rho[i, j, k] + rhobar[i, j, k]) / jac[i, j, k]

    @ivy return -rhotot * (
        (phitheta[i, j, k, 1] - phitheta[i - 1, j, k, 1]) / dx +
        (phitheta[i, j, k, 2] - phitheta[i, j - 1, k, 2]) / dy +
        (phitheta[i, j, k, 3] - phitheta[i, j, k - 1, 3]) / dz
    )
end
