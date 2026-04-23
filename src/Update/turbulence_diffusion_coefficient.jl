"""
```julia 
turbulence_diffusion_coefficient(state::State, i::Integer, j::Integer, k::Integer, variable::KM)
```

Compute the eddy diffusion coefficient for momentum at ``(i, j, k)``. 

The eddy diffusion coefficient for momentum is given by 

```math 
    K_M = l_v \\sqrt{2 e_\\mathrm{k}} \\;, 
```

with turbulence mixing length `l_v` stored in `state.turbulence.turbulenceconstants.lv`.

```julia 
turbulence_diffusion_coefficient(state::State, i::Integer, j::Integer, k::Integer, variable::KH)
```

Compute the eddy diffusion coefficient for heat at ``(i, j, k)``. 

The eddy diffusion coefficient for heat is given by 

```math 
    K_H  = l_b \\sqrt{2 e_\\mathrm{k}} \\;, 
```

with turbulence mixing length `l_b` stored in `state.turbulence.turbulenceconstants.lb`.

```julia 
turbulence_diffusion_coefficient(state::State, i::Integer, j::Integer, k::Integer, variable::KEK)
```

Compute the eddy diffusion coefficient for turbulent kinetic energy at ``(i, j, k)``. 

The eddy diffusion coefficient for turbulent kinetic energy is given by 

```math 
    K_{e_\\mathrm{k}}  = l_t \\sqrt{2 e_\\mathrm{k}} \\;, 
```

with turbulence mixing length `l_t` stored in `state.turbulence.turbulenceconstants.lt`.

# Arguments:

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `variable`: Eddy diffusion coefficient to be computed.
"""
function turbulence_diffusion_coefficient end

function turbulence_diffusion_coefficient(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::KM,
)::AbstractFloat
    (; tke) = state.turbulence.turbulencepredictands
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere
    (; lv) = state.turbulence.turbulenceconstants

    return lv * sqrt(2 * tke[i, j, k] / (rho[i, j, k] + rhobar[i, j, k]))
end

function turbulence_diffusion_coefficient(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::KH,
)::AbstractFloat
    (; tke) = state.turbulence.turbulencepredictands
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere
    (; lb) = state.turbulence.turbulenceconstants

    return lb * sqrt(2 * tke[i, j, k] / (rho[i, j, k] + rhobar[i, j, k]))
end

function turbulence_diffusion_coefficient(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::KEK,
)::AbstractFloat
    (; tke) = state.turbulence.turbulencepredictands
    (; rho) = state.variables.predictands
    (; rhobar) = state.atmosphere
    (; lt) = state.turbulence.turbulenceconstants

    return lt * sqrt(2 * tke[i, j, k] / (rho[i, j, k] + rhobar[i, j, k]))
end