"""
```julia
apply_saturation_scheme!(state::State, dt::AbstractFloat)
```

Entry point for wave saturation scheme based on test case type.

Dispatches to the appropriate saturation implementation depending on
simulation configuration.

# Arguments

  - `state::State`: Complete simulation state
  - `dt::AbstractFloat`: Time step for saturation calculation
"""
function apply_saturation_scheme!(state::State, dt::AbstractFloat)
    (; testcase) = state.namelists.setting
    apply_saturation_scheme!(state, dt, testcase)
    return
end

"""
```julia
apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    testcase::AbstractTestCase,
)
```

No-op for non-WKB test cases.

Standard test cases don't use ray tracing or wave saturation.

# Arguments

  - `state::State`: Simulation state (unused)
  - `dt::AbstractFloat`: Time step (unused)
  - `testcase::AbstractTestCase`: Non-WKB test case
"""
function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    testcase::AbstractTestCase,
)
    return
end

"""
```julia
apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    testcase::AbstractWKBTestCase,
)
```

Apply saturation scheme for WKB test cases based on WKB mode.

Dispatches to the specific WKB mode implementation.

# Arguments

  - `state::State`: Simulation state containing WKB configuration
  - `dt::AbstractFloat`: Time step for saturation calculation
  - `testcase::AbstractWKBTestCase`: WKB test case specification
"""
function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    testcase::AbstractWKBTestCase,
)
    (; wkb_mode) = state.namelists.wkb
    apply_saturation_scheme!(state, dt, wkb_mode)
    return
end

"""
```julia
apply_saturation_scheme!(state::State, dt::AbstractFloat, wkb_mode::SteadyState)
```

No-op for steady-state WKB mode.

Steady-state mode handles saturation differently, typically within
the propagation step rather than as a separate scheme.

# Arguments

  - `state::State`: Simulation state (unused)
  - `dt::AbstractFloat`: Time step (unused)
  - `wkb_mode::SteadyState`: Steady-state WKB mode
"""
function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    wkb_mode::SteadyState,
)
    return
end

"""
```julia
apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    wkb_mode::AbstractWKBMode,
)
```

Apply wave breaking saturation scheme to limit wave amplitudes.

Implements a parameterization of wave breaking that prevents wave amplitudes
from exceeding critical values, representing the onset of convective instability
and turbulent mixing in the atmosphere.

# Arguments

  - `state::State`: Complete simulation state
  - `dt::AbstractFloat`: Time step for diffusion calculation
  - `wkb_mode::AbstractWKBMode`: WKB mode (MultiColumn, SingleColumn, etc.)

# Saturation Theory

Wave breaking occurs when wave-induced potential temperature perturbations
exceed a critical fraction of the background stratification:

`|θ'| > α_sat · θ₀ · N/g`

where:

  - `α_sat`: Saturation parameter (typically 0.1-1.0)
  - `θ₀`: Background potential temperature
  - `N`: Brunt-Väisälä frequency
  - `g`: Gravitational acceleration

# Algorithm

 1. **Saturation Integral**: Compute `mb2 = ∫ (wave momentum flux)`
 2. **Critical Value**: `mb2_crit = α_sat² · N²`
 3. **Diffusion Coefficient**: `κ = (mb2 - mb2_crit) / (2·dt·mb2k2)` if mb2 > mb2_crit
 4. **Wave Action Reduction**: `A_new = A_old · max(0, 1 - 2κ·dt·|k|²)`
 5. **Ray Removal**: Remove rays with zero wave action

# Saturation Integrals

  - `mb2`: Total squared momentum flux `∫ ρ·A·ω·(k²+l²+m²)`
  - `mb2k2`: Weighted momentum flux for diffusion `∫ ρ·A·ω·|k|²·(vertical transport)`

# Physical Interpretation

  - **Pre-breaking**: Waves propagate without significant damping
  - **At saturation**: Turbulent diffusion removes excess wave energy
  - **Post-breaking**: Wave field relaxes toward saturated state

# Grid Cell Processing

Applied independently to each grid cell, allowing for spatially varying
saturation based on local wave field and background conditions.

# Diagnostics

  - Reports saturation violations if residual exceeds tolerance
  - Removes rays with negligible wave action density

# Conservation

Total wave action is not conserved during saturation - the "lost" action
represents conversion to turbulent kinetic energy and mixing.
"""
function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    wkb_mode::AbstractWKBMode,
)
    (; domain, grid) = state
    (; nray, rays, diffusion) = state.wkb
    (; sizex, sizey) = state.namelists.domain
    (; lsaturation, alpha_sat) = state.namelists.wkb
    (; io, jo, i0, i1, j0, j1, k0, k1) = state.domain
    (; lx, ly, dx, dy, ztfc) = state.grid

    if !lsaturation
        return
    end

    for kz in k0:k1, jy in j0:j1, ix in i0:i1

        # Compute saturation integrals for wave-action reduction.
        (mb2, mb2k2) = compute_saturation_integrals(state, (ix, jy, kz))

        # Calculate the turbulent eddy diffusivity.
        n2r = interpolate_stratification(ztfc[ix, jy, kz], state, N2())
        if mb2k2 == 0 || mb2 < alpha_sat^2 * n2r^2
            diffusion[ix, jy, kz] = 0
        else
            diffusion[ix, jy, kz] =
                (mb2 - alpha_sat^2 * n2r^2) / (2 * dt * mb2k2)
        end

        # Reduce the wave-action density.
        for iray in 1:nray[ix, jy, kz]
            xr = rays.x[iray, ix, jy, kz]
            yr = rays.y[iray, ix, jy, kz]
            zr = rays.z[iray, ix, jy, kz]

            if sizex > 1
                ix = floor(Int, (xr - lx[1]) / dx) + i0 - io
            else
                ix = i0
            end

            if sizey > 1
                jy = floor(Int, (yr - ly[1]) / dy) + j0 - jo
            else
                jy = j0
            end

            kz = get_next_half_level(ix, jy, zr, domain, grid)

            wnrk = rays.k[iray, ix, jy, kz]
            wnrl = rays.l[iray, ix, jy, kz]
            wnrm = rays.m[iray, ix, jy, kz]

            kappa = diffusion[ix, jy, kz]

            rays.dens[iray, ix, jy, kz] *=
                max(0, 1 - dt * 2 * kappa * (wnrk^2 + wnrl^2 + wnrm^2))
        end

        # Compute the saturation integrals again for diagnostics.
        (mb2, mb2k2) = compute_saturation_integrals(state, (ix, jy, kz))

        # Check if saturation is violated.
        n2r = interpolate_stratification(ztfc[ix, jy, kz], state, N2())
        if mb2 - alpha_sat^2 * n2r^2 > 1.0E-3 * alpha_sat^2 * n2r^2
            println("Saturation violated at (ix, jy, kz) = ", (ix, jy, kz))
            println("mb2[ix, jy, kz] = ", mb2)
            println("alpha_sat^2 * n2r^2 = ", alpha_sat^2 * n2r^2)
            println("")
        end
    end

    # Rmove rays with zero wave-action density.
    remove_rays!(state)

    return
end
