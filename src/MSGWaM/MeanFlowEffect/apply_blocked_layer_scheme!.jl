"""
```julia
apply_blocked_layer_scheme!(state::State)
```

Entry point for applying blocked layer drag scheme based on test case.

Dispatches to the appropriate blocked layer implementation depending on
the test case type.

# Arguments

  - `state::State`: Complete simulation state
"""
function apply_blocked_layer_scheme!(state::State)
    (; testcase) = state.namelists.setting
    apply_blocked_layer_scheme!(state, testcase)
    return
end

"""
```julia
apply_blocked_layer_scheme!(state::State, testcase::AbstractWKBTestCase)
```

Default implementation for general WKB test cases. No blocked layer scheme applied.

Most WKB test cases do not require flow blocking parameterization.

# Arguments

  - `state::State`: Complete simulation state (unused)
  - `testcase::AbstractWKBTestCase`: General WKB test case
"""
function apply_blocked_layer_scheme!(
    state::State,
    testcase::AbstractWKBTestCase,
)
    return
end

"""
```julia
apply_blocked_layer_scheme!(state::State, testcase::WKBMountainWave)
```

Apply blocked layer drag scheme for mountain wave simulations.

Implements a parameterization of blocked flow drag that occurs when stable
stratification and strong winds create a blocked layer near the surface,
modifying the gravity wave drag profile.

# Arguments

  - `state::State`: Complete simulation state
  - `testcase::WKBMountainWave`: Mountain wave test case

# Physical Background

## Flow Blocking

When air flows over topography with strong stratification:

  - **Blocked layer**: Region where flow is deflected around rather than over terrain
  - **Critical height**: `h_c = N·h/U` where blocking transitions to wave generation
  - **Drag distribution**: Modified to account for both wave drag and blocked flow drag

## Blocking Criteria

Blocking occurs when the linearity parameter exceeds threshold:
`L = N·h/U > L_threshold`

where:

  - `N`: Brunt-Väisälä frequency
  - `h`: Terrain height
  - `U`: Background wind speed
  - `L_threshold`: Critical linearity parameter

# Algorithm

 1. **Blocking Check**: Only apply if blocking is enabled in configuration
 2. **Grid Cell Loop**: Process each grid cell in computational domain
 3. **Blocked Fraction**: Compute fraction of cell within blocked layer
 4. **Average Wavenumbers**: Weight by topographic spectrum amplitudes
 5. **Perpendicular Wind**: Component normal to terrain wave vector
 6. **Drag Calculation**: Compute blocked flow drag using drag coefficient
 7. **Tendency Blending**: Blend wave drag and blocked drag based on height

# Drag Formulation

Blocked flow drag per unit volume:
`F_drag = -C_d · ρ · k_h / (2π) · |U_⊥| · U_⊥`

where:

  - `C_d`: Drag coefficient
  - `ρ`: Total density (perturbation + background)
  - `k_h`: Horizontal wavenumber magnitude
  - `U_⊥`: Wind component perpendicular to terrain waves

# Vertical Distribution

  - **Below blocked layer**: Pure blocked flow drag
  - **Above blocked layer**: Pure gravity wave drag
  - **Transition zone**: Weighted combination based on height fraction

# Momentum Budget

Conserves total momentum while redistributing drag between:

  - Gravity wave mechanisms (resolved by ray tracing)
  - Blocked flow mechanisms (parameterized drag)

# Configuration

Controlled by:

  - `blocking`: Enable/disable blocking scheme
  - `drag_coefficient`: Strength of blocked flow drag
  - `long_threshold`: Critical linearity parameter for blocking onset
"""
function apply_blocked_layer_scheme!(state::State, testcase::WKBMountainWave)
    (; blocking, drag_coefficient) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dz, jac, ztildetfc, topography_spectrum, k_spectrum, l_spectrum) =
        state.grid
    (; rhostrattfc) = state.atmosphere
    (; rho, u, v) = state.variables.predictands
    (; zb) = state.wkb
    (; dudt, dvdt, dthetadt) = state.wkb.tendencies

    if !blocking
        return
    end

    # Initialize arrays for blocked-flow drag computation.
    (kavg, uperp, drag) = (zeros(2) for i in 1:3)

    # Adjust the drag to account for blocking.
    for kz in k0:k1, jy in j0:j1, ix in i0:i1
        fraction =
            (
                min(zb[ix, jy], ztildetfc[ix, jy, kz]) -
                ztildetfc[ix, jy, kz - 1]
            ) / jac[ix, jy, kz] / dz
        @views if fraction <= 0
            continue
        else
            kavg[1] =
                sum(
                    abs.(topography_spectrum[:, ix, jy]) .*
                    k_spectrum[:, ix, jy],
                ) / sum(abs.(topography_spectrum[:, ix, jy]))
            kavg[2] =
                sum(
                    abs.(topography_spectrum[:, ix, jy]) .*
                    l_spectrum[:, ix, jy],
                ) / sum(abs.(topography_spectrum[:, ix, jy]))

            uperp .=
                (
                    (u[ix, jy, kz] .+ u[ix - 1, jy, kz]) .* kavg[1] .+
                    (v[ix, jy, kz] .+ v[ix, jy - 1, kz]) .* kavg[2]
                ) ./ 2 .* kavg ./ dot(kavg, kavg)
            drag .=
                -drag_coefficient .*
                (rho[ix, jy, kz] .+ rhostrattfc[ix, jy, kz]) .*
                sqrt(dot(kavg, kavg)) ./ (2 .* pi) .* sqrt(dot(uperp, uperp)) .*
                uperp
            dudt[ix, jy, kz] =
                fraction * drag[1] + (1 - fraction) * dudt[ix, jy, kz]
            dvdt[ix, jy, kz] =
                fraction * drag[2] + (1 - fraction) * dvdt[ix, jy, kz]
            dthetadt[ix, jy, kz] = (1 - fraction) * dthetadt[ix, jy, kz]
        end
    end
end
