"""
```julia
turbulence_computation!(
    state::State,
    p2::Predictands,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    process::Dissipation,
)
```

Performs the dissipation step of the turbulence computation.

```julia
turbulence_computation!(
    state::State,
    p2::Predictands,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    process::Diffusion,
)
```

Performs the diffusion step of the turbulence computation.

# Arguments

  - `state`: Model state.

  - `p2`: The predictands that are used to compute the transporting velocities in the computation of the fluxes.

  - `dtstage`: Fractional time step.

  - `time`: Simulation time.

  - `process`: Denotes the process of the turbulence computation
"""
function turbulence_computation! end

function turbulence_computation!(
    state::State,
    p2::Predictands,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    process::Dissipation,
)
    (; cepsilon) = state.turbulence.turbulenceconstants
    (; tke) = state.turbulence.turbulencepredictands
    
    tke = cepsilon

    return
end

function turbulence_computation!(
    state::State,
    p2::Predictands,
    dtstage::AbstractFloat,
    time::AbstractFloat,
    process::Diffusion,
)
    
    # Function Definition

    return
end
