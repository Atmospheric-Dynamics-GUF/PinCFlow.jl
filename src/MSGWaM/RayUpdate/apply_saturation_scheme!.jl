"""
```julia
apply_saturation_scheme!(state::State, dt::AbstractFloat)
```

Apply the saturation scheme by dispatching to a test-case-specific method.

```julia
apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    testcase::AbstractTestCase,
)
```

Return for non-WKB test cases.

```julia
apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    testcase::AbstractWKBTestCase,
)
```

Apply the saturation scheme by dispatching to a WKB-mode-specific method.

```julia
apply_saturation_scheme!(state::State, dt::AbstractFloat, wkb_mode::SteadyState)
```

Return for steady-state configurations.

In steady-state mode, saturation is handled by [`PinCFlow.MSGWaM.RayUpdate.propagate_rays!`](@ref).

```julia
apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    wkb_mode::AbstractWKBMode,
)
```

Apply the saturation scheme.

Saturation is assumed to occur when the static-instability criterion

```math
\\sum\\limits_\\alpha \\left(m_\\alpha \\left|b_{\\mathrm{w}, \\alpha}\\right|\\right)^2 f_\\alpha \\geq \\alpha_\\mathrm{S}^2 N^4
```

is locally fulfilled (i.e. within a grid cell). Therein, ``\\boldsymbol{k}_\\alpha = \\left(k_\\alpha, l_\\alpha, m_\\alpha\\right)^\\mathrm{T}`` is the wavevector, ``\\left|b_{\\mathrm{w}, \\alpha}\\right|^2`` is the squared gravity-wave amplitude of the buoyancy, ``f_\\alpha`` is the maximum grid-cell fraction each ray volume can cover, ``N`` is the buoyancy frequency and ``\\alpha_\\mathrm{S}`` is a saturation coefficient that represents the uncertainties of the criterion. The phase-space wave-action density is then reduced in accordance with

```math
\\frac{\\Delta \\left|b_{\\mathrm{w}, \\alpha}\\right|^2}{\\Delta t} = - 2 K \\left|\\boldsymbol{k}_\\alpha\\right|^2 \\left|b_{\\mathrm{w}, \\alpha}\\right|^2,
```

which is based on the assumption that wave breaking leads to turbulent fluxes that may be parameterized with a flux-gradient ansatz. The turbulent viscosity and diffusivity

```math
K = \\left[2 \\Delta t\\sum\\limits_\\alpha \\left(m_\\alpha \\left|b_{\\mathrm{w}, \\alpha}\\right| \\left|\\boldsymbol{k}_\\alpha\\right|\\right)^2 f_\\alpha\\right]^{- 1} \\max \\left[0, \\sum_\\alpha \\left(m_\\alpha \\left|b_{\\mathrm{w}, \\alpha}\\right|\\right)^2 f_\\alpha - \\alpha_\\mathrm{S}^2 N^4\\right]
```

is such that wave action is reduced exactly to the saturation threshold. The two sums involved in this scheme (discretizations of spectral integrals) are computed with `compute_spectral_integrals`.

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `testcase`: Test case on which the current simulation is based.

  - `wkb_mode`: Approximations used by MSGWaM.

# See also

  - [`PinCFlow.MSGWaM.RayOperations.compute_saturation_integrals`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_stratification`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.get_next_half_level`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.remove_rays!`](@ref)
"""
function apply_saturation_scheme! end

function apply_saturation_scheme!(state::State, dt::AbstractFloat)
    (; testcase) = state.namelists.setting
    apply_saturation_scheme!(state, dt, testcase)
    return
end

function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    testcase::AbstractTestCase,
)
    return
end

function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    testcase::AbstractWKBTestCase,
)
    (; wkb_mode) = state.namelists.wkb
    apply_saturation_scheme!(state, dt, wkb_mode)
    return
end

function apply_saturation_scheme!(
    state::State,
    dt::AbstractFloat,
    wkb_mode::SteadyState,
)
    return
end

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

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1

        # Compute saturation integrals for wave-action reduction.
        (mb2, mb2k2) = compute_saturation_integrals(state, i, j, k)

        # Calculate the turbulent eddy diffusivity.
        n2r = interpolate_stratification(ztfc[i, j, k], state, N2())
        if mb2k2 == 0 || mb2 < alpha_sat^2 * n2r^2
            diffusion[i, j, k] = 0
        else
            diffusion[i, j, k] = (mb2 - alpha_sat^2 * n2r^2) / (2 * dt * mb2k2)
        end

        # Reduce the wave-action density.
        for r in 1:nray[i, j, k]
            xr = rays.x[r, i, j, k]
            yr = rays.y[r, i, j, k]
            zr = rays.z[r, i, j, k]

            if sizex > 1
                ix = floor(Int, (xr + lx / 2) / dx) + i0 - io
            else
                ix = i0
            end

            if sizey > 1
                jy = floor(Int, (yr + ly / 2) / dy) + j0 - jo
            else
                jy = j0
            end

            kz = get_next_half_level(ix, jy, zr, state)

            wnrk = rays.k[r, i, j, k]
            wnrl = rays.l[r, i, j, k]
            wnrm = rays.m[r, i, j, k]

            kappa = diffusion[ix, jy, kz]

            rays.dens[r, i, j, k] *=
                max(0, 1 - dt * 2 * kappa * (wnrk^2 + wnrl^2 + wnrm^2))
        end

        # Compute the saturation integrals again for diagnostics.
        (mb2, mb2k2) = compute_saturation_integrals(state, i, j, k)

        # Check if saturation is violated.
        n2r = interpolate_stratification(ztfc[i, j, k], state, N2())
        if mb2 - alpha_sat^2 * n2r^2 > 1.0E-3 * alpha_sat^2 * n2r^2
            println("Saturation violated at (i, j, k) = ", (i, j, k))
            println("mb2[i, j, k] = ", mb2)
            println("alpha_sat^2 * n2r^2 = ", alpha_sat^2 * n2r^2)
            println("")
        end
    end

    # Rmove rays with zero wave-action density.
    remove_rays!(state)

    return
end
