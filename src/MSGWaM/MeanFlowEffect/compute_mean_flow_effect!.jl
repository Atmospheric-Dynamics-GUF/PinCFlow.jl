"""
```julia
compute_mean_flow_effect!(state::State)
```

Calculate the mean-flow impact of unresolved gravity waves, based on the test case.

Dispatches to the appropriate method, depending on test case.

# Arguments

  - `state`: Model state.
"""
function compute_mean_flow_effect!(state::State)
    (; testcase) = state.namelists.setting
    compute_mean_flow_effect!(state, testcase)
    return
end

"""
```julia
compute_mean_flow_effect!(state::State, testcase::AbstractTestCase)
```

Return for non-WKB test cases.

# Arguments

  - `state`: Model state.
  - `testcase`: Test case on which the current simulation is based.
"""
function compute_mean_flow_effect!(state::State, testcase::AbstractTestCase)
    return
end

"""
```julia
compute_mean_flow_effect!(state::State, testcase::AbstractWKBTestCase)
```

Calculate the mean-flow impact of unresolved gravity waves.

This method first computes several spectral integrals (using `compute_gw_integrals!`), most of which represent gravity-wave fluxes. After the boundary conditions for these have been enforced (using `set_boundaries!`), the corresponding tendencies are calculated (using `compute_gw_tendencies!`). These also have boundary conditions that need to be enforced (once again using `set_boundaries!`) before they are smoothed to remove small-scale features that may occur due to a coarse ray-volume distribution (using `smooth_gw_tendencies!`). Afterwards, if MSGWaM parameterizes mountain waves, the tendencies are adjusted to account for the formation of blocked layers (using `apply_blocked_layer_scheme!`), before the boundary conditions are enforced again.

# Arguments

  - `state`: Model state.
  - `testcase`: Test case on which the current simulation is based.

# See also

  - [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_integrals!`](@ref)
  - [`PinCFlow.Boundaries.set_boundaries!`](@ref)
  - [`PinCFlow.MSGWaM.MeanFlowEffect.compute_gw_tendencies!`](@ref)
  - [`PinCFlow.MSGWaM.MeanFlowEffect.smooth_gw_tendencies!`](@ref)
  - [`PinCFlow.MSGWaM.MeanFlowEffect.apply_blocked_layer_scheme!`](@ref)
"""
function compute_mean_flow_effect!(state::State, testcase::AbstractWKBTestCase)
    (; wkb_mode) = state.namelists.wkb

    compute_gw_integrals!(state, wkb_mode)

    set_boundaries!(state, BoundaryGWIntegrals())

    compute_gw_tendencies!(state)

    set_boundaries!(state, BoundaryGWTendencies())

    smooth_gw_tendencies!(state)

    apply_blocked_layer_scheme!(state)

    set_boundaries!(state, BoundaryGWTendencies())

    return
end
