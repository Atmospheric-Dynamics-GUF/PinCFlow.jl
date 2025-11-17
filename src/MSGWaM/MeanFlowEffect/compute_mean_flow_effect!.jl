"""
```julia
compute_mean_flow_effect!(state::State)
```

Calculate the mean-flow impact of unresolved gravity waves by dispatching to a test-case-specific method.

```julia
compute_mean_flow_effect!(state::State, test_case::AbstractTestCase)
```

Return for non-WKB test cases.

```julia
compute_mean_flow_effect!(state::State, test_case::AbstractWKBTestCase)
```

Calculate the mean-flow impact of unresolved gravity waves.

This method first computes several spectral integrals (using `compute_gw_integrals!`), most of which represent gravity-wave fluxes. After the boundary conditions for these have been enforced (using `set_boundaries!`), the corresponding tendencies are calculated (using `compute_gw_tendencies!`). These also have boundary conditions that need to be enforced (once again using `set_boundaries!`) before they are smoothed to remove small-scale features that may occur due to a coarse ray-volume distribution (using `smooth_gw_tendencies!`). Afterwards, if MSGWaM parameterizes mountain waves, the tendencies are adjusted to account for the formation of blocked layers (using `apply_blocked_layer_scheme!`), before the boundary conditions are enforced again.

# Arguments

  - `state`: Model state.

  - `test_case`: Test case on which the current simulation is based.

# See also

  - [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_integrals!`](@ref)

  - [`PinCFlow.Boundaries.set_boundaries!`](@ref)

  - [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_tendencies!`](@ref)

  - [`PinCFlow.MSGWaM.MeanFlowEffect.smooth_gw_tendencies!`](@ref)

  - [`PinCFlow.MSGWaM.MeanFlowEffect.apply_blocked_layer_scheme!`](@ref)
"""
function compute_mean_flow_effect! end

function compute_mean_flow_effect!(state::State)
    (; test_case) = state.namelists.setting
    compute_mean_flow_effect!(state, test_case)
    return
end

function compute_mean_flow_effect!(state::State, test_case::AbstractTestCase)
    return
end

function compute_mean_flow_effect!(state::State, test_case::AbstractWKBTestCase)
    compute_gw_integrals!(state)

    set_boundaries!(state, BoundaryWKBIntegrals())

    compute_gw_tendencies!(state)

    set_boundaries!(state, BoundaryWKBTendencies())

    smooth_gw_tendencies!(state)

    apply_blocked_layer_scheme!(state)

    set_boundaries!(state, BoundaryWKBTendencies())

    return
end
