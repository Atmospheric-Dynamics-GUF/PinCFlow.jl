"""
```julia 
wkb_initialization!(state::State)
```

Initialize MS-GWaM by dispatching to a configuration-specific method.

```julia 
wkb_initialization!(state::State, wkb_mode::Val{:NoWKB})
```

Return for configurations without MS-GWaM.

```julia 
wkb_initialization!(
    state::State,
    wkb_mode::Union{Val{:MultiColumn}, Val{:SingleColumn}, Val{:SteadyState}},
)
```

Initialize MS-GWaM by calling `initialize_rays!`. 

For configurations with tracer transport and next-order impact calculation, also compute the leading-order gravity-wave integrals, enforce boundary conditions, and smooth the leading-order wave amplitudes.

# Arguments 

  - `state`: Model state.

  - `wkb_mode`: Approximations used by MS-GWaM.

# See also

  - [`PinCFlow.MSGWaM.RayUpdate.initialize_rays!`](@ref)

  - [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_integrals!`](@ref)

  - [`PinCFlow.Boundaries.set_boundaries!`](@ref)

  - [`PinCFlow.MSGWaM.MeanFlowEffect.smoothing!`](@ref)
"""
function wkb_initialization! end

function wkb_initialization!(state::State)
    (; wkb_mode) = state.namelists.wkb

    @dispatch_wkb_mode wkb_initialization!(state, Val(wkb_mode))
    return
end

function wkb_initialization!(state::State, wkb_mode::Val{:NoWKB})
    return
end

function wkb_initialization!(
    state::State,
    wkb_mode::Union{Val{:MultiColumn}, Val{:SingleColumn}, Val{:SteadyState}},
)
    (; tracer_setup, next_order_impact) = state.namelists.tracer

    initialize_rays!(state)

    if tracer_setup === :TracerOn && next_order_impact
        compute_gw_integrals!(state)
        set_boundaries!(state, BoundaryWKBIntegrals())
        smoothing!(state, Integrals())
    end
    return
end
