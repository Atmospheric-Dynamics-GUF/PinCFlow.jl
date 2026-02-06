"""
```julia 
buoyancy_term(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
)::AbstractFloat
```

Compute the buoyancy term in the prognostic equation for the TKE by dispatching to the appropriate model-specific method.

The buoyancy term is given by 

```math 
\\mathcal{B} = -K_H \\left(N^2 + \\frac{\\partial b}{\\partial z}\\right)
```

```julia 
buoyancy_term(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
    model::Boussinesq,
)::AbstractFloat
```

Compute the buoyancy term in Boussinesq mode.
```

```julia 
buoyancy_term(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
    model::Union{PseudoIncompressible, Compressible},
)::AbstractFloat
```

Compute the buoyancy term in pseudo-incompressible and compressible mode.

# Arguments

  - `state`: Model state.

  - `p0`: The predictands used to compute the buoyancy term.

  - `i`: Zonal grid-cell index. 

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `model`: Dynamic equations.
"""
function buoyancy_term end

function buoyancy_term(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
)::AbstractFloat
    (; model) = state.namelists.atmosphere

    return buoyancy_term(state, p0, i, j, k, model)
end

function buoyancy_term(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
    model::Boussinesq,
)::AbstractFloat
    (; rhop) = p0
    (; pbar, rhobar) = state.atmosphere
    (; kh) = state.turbulence.turbulencediffusioncoefficients
    (; jac, dz) = state.grid
    (; g_ndim) = state.constants

    bu = g_ndim * (1 / (rhop[i, j, k + 1] / rhobar[i, j, k + 1] + 1) - 1)
    bd = g_ndim * (1 / (rhop[i, j, k - 1] / rhobar[i, j, k - 1] + 1) - 1)

    buoyancy =
        -kh[i, j, k] * (n2[i, j, k] + (bu - bd) / (jac[i, j, k] * 2.0 * dz))

    return buoyancy
end

function buoyancy_term(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
    model::Union{PseudoIncompressible, Compressible},
)::AbstractFloat
    (; rho) = p0
    (; pbar, rhobar, n2) = state.atmosphere
    (; kh) = state.turbulence.turbulencediffusioncoefficients
    (; jac, dz) = state.grid
    (; g_ndim) = state.constants

    bu = g_ndim * (1 / (rho[i, j, k + 1] / rhobar[i, j, k + 1] + 1) - 1)
    bd = g_ndim * (1 / (rho[i, j, k - 1] / rhobar[i, j, k - 1] + 1) - 1)

    buoyancy =
        -kh[i, j, k] * (n2[i, j, k] + (bu - bd) / (jac[i, j, k] * 2 * dz))

    return buoyancy
end
