"""
```julia
compute_mean_flow_effect!(state::State, dt::AbstractFloat)
```

Calculate the mean-flow impact of unresolved gravity waves by dispatching to a WKB-mode-specific method.

```julia
compute_mean_flow_effect!(
    state::State,
    dt::AbstractFloat,
    wkb_mode::Val{:NoWKB},
)
```

Return for non-WKB configurations.

```julia
compute_mean_flow_effect!(
    state::State,
    dt::AbstractFloat,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)
```

Calculate the mean-flow impact of unresolved gravity waves.

This method consists of the following steps:

  1. Computation of several spectral integrals (using `compute_gw_integrals!`), most of which represent gravity-wave fluxes. 
  
  1. Enforcement of boundary conditions of the spectral integrals (using `set_boundaries!`). 
  
  1. Removal of small-scale features in tracer-specific integrals that may occur due to a coarse ray-volume distribution (using `smoothing`).
  
  1. Further computation of tracer-specific spectral integrals based on the previously computed integrals.

  1. Enforcement of boundary conditions of the spectral integrals. 

  1. Computation of the corresponding tendencies (using `compute_gw_tendencies!`).

  1. Enforce boundary conditions for the tendencies.

  1. Smoothing of the tencencies.

  1. If MS-GWaM parameterizes mountain waves, the tendencies are adjusted to account for the formation of blocked layers (using `include_blocked_flow_drag!`)
  
  1. Enforce boundary conditions again.

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `wkb_mode`: Approximations used by MS-GWaM.

# See also

  - [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_integrals!`](@ref)

  - [`PinCFlow.Boundaries.set_boundaries!`](@ref)

  - [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_tendencies!`](@ref)

  - [`PinCFlow.MSGWaM.MeanFlowEffect.smoothing!`](@ref)

  - [`PinCFlow.MSGWaM.BlockedLayer.include_blocked_flow_drag!`](@ref)
"""
function compute_mean_flow_effect! end

function compute_mean_flow_effect!(state::State, dt::AbstractFloat)
    (; wkb_mode) = state.namelists.wkb
    @dispatch_wkb_mode compute_mean_flow_effect!(state, dt, Val(wkb_mode))
    return
end

function compute_mean_flow_effect!(
    state::State,
    dt::AbstractFloat,
    wkb_mode::Val{:NoWKB},
)
    return
end

function compute_mean_flow_effect!(
    state::State,
    dt::AbstractFloat,
    wkb_mode::Union{Val{:SteadyState}, Val{:SingleColumn}, Val{:MultiColumn}},
)
    compute_gw_integrals!(state)

    set_boundaries!(state, BoundaryWKBIntegrals())

    smoothing!(state, Integrals())

    compute_gw_integrals!(state, dt)

    set_boundaries!(state, BoundaryWKBIntegrals())

    compute_gw_tendencies!(state)

    set_boundaries!(state, BoundaryWKBTendencies())

    smoothing!(state, Tendencies())

    include_blocked_flow_drag!(state)

    set_boundaries!(state, BoundaryWKBTendencies())

    return
end
