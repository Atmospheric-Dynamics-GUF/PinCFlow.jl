"""
```julia
compute_mean_flow_effect!(state::State)
```

Entry point for computing gravity wave effects on the mean flow.

Dispatches to the appropriate implementation based on test case type.

# Arguments

  - `state::State`: Complete simulation state
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

No-op for non-WKB test cases.

Standard test cases do not use gravity wave-mean flow interactions.

# Arguments

  - `state::State`: Simulation state (unused)
  - `testcase::AbstractTestCase`: Non-WKB test case
"""
function compute_mean_flow_effect!(state::State, testcase::AbstractTestCase)
    return
end

"""
```julia
compute_mean_flow_effect!(state::State, testcase::AbstractWKBTestCase)
```

Compute gravity wave effects on mean flow for WKB test cases.

Implements the complete calculation of gravity wave momentum and heat transport
effects on the resolved mean flow, including wave stress divergence and
Eliassen-Palm flux contributions.

# Arguments

  - `state::State`: Complete simulation state
  - `testcase::AbstractWKBTestCase`: WKB test case specification

# Algorithm Overview

 1. **Compute Integrals**: Calculate momentum flux and EP flux components
 2. **Set Boundaries**: Apply boundary conditions to integral fields
 3. **Compute Tendencies**: Calculate drag and heating from flux divergence
 4. **Set Boundaries**: Apply boundary conditions to tendency fields
 5. **Smooth Tendencies**: Apply spatial filtering for numerical stability
 6. **Apply Blocking**: Add blocked layer drag for mountain waves
 7. **Final Boundaries**: Ensure boundary consistency

# Physical Process Chain

## Wave-Mean Flow Interaction

  - **Wave propagation**: Rays transport momentum and energy
  - **Local deposition**: Waves deposit momentum where they dissipate
  - **Mean flow acceleration**: Momentum convergence accelerates mean flow
  - **Feedback**: Modified mean flow affects subsequent wave propagation

## Momentum Transport

  - **Reynolds stress**: `⟨u'w'⟩`, `⟨v'w'⟩` from vertical momentum flux
  - **Form stress**: `⟨u'u'⟩`, `⟨v'v'⟩`, `⟨u'v'⟩` from horizontal correlations
  - **Pressure forces**: Implicit in wave momentum budget

## Heat Transport

  - **Eliassen-Palm fluxes**: `⟨u'θ'⟩`, `⟨v'θ'⟩` modify temperature gradients
  - **Rotational effects**: Heat fluxes coupled to momentum through geostrophy

# Boundary Conditions

Applied at each stage to ensure:

  - **Physical consistency**: Proper flux conservation at boundaries
  - **Numerical stability**: Smooth transitions in halo regions
  - **Domain coupling**: Correct MPI boundary exchange

# Spatial Filtering

Optional smoothing reduces:

  - **Numerical noise**: From discrete ray sampling
  - **Grid-scale oscillations**: Unphysical small-scale features
  - **Stability issues**: In regions with few rays

# Blocked Layer Integration

For mountain waves, combines:

  - **Wave drag**: From resolved gravity wave propagation
  - **Blocked flow drag**: From parameterized sub-grid blocking

# Conservation Properties

Maintains conservation of:

  - **Total momentum**: Momentum lost by waves gained by mean flow
  - **Total energy**: Energy dissipated by waves heats the atmosphere
  - **Angular momentum**: Consistent with rotation effects

# Output

Updates gravity wave tendency fields:

  - `tendencies.dudt`: Zonal wind tendency [m/s²]
  - `tendencies.dvdt`: Meridional wind tendency [m/s²]
  - `tendencies.dthetadt`: Potential temperature tendency [K/s]
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
