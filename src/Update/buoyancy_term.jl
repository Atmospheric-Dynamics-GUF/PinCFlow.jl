"""
```julia 
buoyancy_term(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
)
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
)
```

Compute the buoyancy term in Boussinesq mode.

"""
function buoyancy_term end

function buoyancy_term(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
)
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
)
    (; rhop) = p0
    (; pbar, rhobar) = state.atmosphere
    (; kh) = state.turbulence.turbulencediffusioncoefficients
    (; jac, dz) = state.grid
    (; g_ndim) = state.constants

    bu = g_ndim * (1 / (rhop[i, j, k + 1] / rhobar[i, j, k + 1] + 1) - 1)
    bd = g_ndim * (1 / (rhop[i, j, k - 1] / rhobar[i, j, k - 1] + 1) - 1)

    buoyancy =
        -kh[i, j, k] * (n2[i, j, k] + (bu - bd) / (jac[i, j, k] * 2.0 * dz))

    # thetau = pbar[i, j, k + 1] / (rhop[i, j, k + 1] + rhobar[i, j, k + 1])
    # theta = pbar[i, j, k] / (rhop[i, j, k] + rhobar[i, j, k])
    # thetad = pbar[i, j, k - 1] / (rhop[i, j, k - 1] + rhobar[i, j, k - 1])

    # wthetap = (thetau - thetad) / (jac[i, j, k] * 2 * dz)
    # buoyancy = -g_ndim * wthetap / theta

    return buoyancy
end

function buoyancy_term(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
    model::Union{PseudoIncompressible, Compressible},
)
    (; rho) = p0
    (; pbar, rhobar, n2) = state.atmosphere
    (; kh) = state.turbulence.turbulencediffusioncoefficients
    (; jac, dz) = state.grid
    (; g_ndim) = state.constants

    bu = g_ndim * (1 / (rho[i, j, k + 1] / rhobar[i, j, k + 1] + 1) - 1)
    bd = g_ndim * (1 / (rho[i, j, k - 1] / rhobar[i, j, k - 1] + 1) - 1)

    buoyancy =
        -kh[i, j, k] * (n2[i, j, k] + 0 * (bu - bd) / (jac[i, j, k] * 2 * dz))

    # thetau = pbar[i, j, k + 1] / (rho[i, j, k + 1] + rhobar[i, j, k + 1])
    # theta = pbar[i, j, k] / (rho[i, j, k] + rhobar[i, j, k])
    # thetad = pbar[i, j, k - 1] / (rho[i, j, k - 1] + rhobar[i, j, k - 1])

    # wthetap = (thetau - thetad) / (jac[i, j, k] * 2 * dz)
    # buoyancy = -g_ndim * wthetap / theta

    return buoyancy
end