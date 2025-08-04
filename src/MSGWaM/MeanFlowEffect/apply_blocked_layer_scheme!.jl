"""
```julia
apply_blocked_layer_scheme!(state::State)
```

Compute the blocked-flow drag and adjust the mean-flow impact accordingly, based on the test case.

Dispatches to the appropriate method depending on the test case.

# Arguments

  - `state`: Model state.
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

Return for non-mountain-wave test cases.

# Arguments

  - `state`: Model state.
  - `testcase`: Test case on which the current simulation is based.
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

Compute the blocked-flow drag and adjust the mean-flow impact accordingly.

The blocked-flow drag is given by

```math
\\left(\\frac{\\partial \\boldsymbol{u}_\\mathrm{b}}{\\partial t}\\right)_\\mathrm{B} = - D \\frac{\\left|\\boldsymbol{k}_h\\right| \\left|\\boldsymbol{u}_\\mathrm{p}\\right|}{2 \\pi} \\boldsymbol{u}_\\mathrm{p},
```

where ``\\boldsymbol{u}_\\mathrm{b} = \\left(u_\\mathrm{b}, v_\\mathrm{b}, 0\\right)^\\mathrm{T}`` is the resolved horizontal wind, ``D`` is a dimensionless drag coefficient (represented by `state.namelists.wkb.drag_coefficient`),

```math
\\boldsymbol{k}_h = \\frac{\\sum_\\alpha \\left|h_{\\mathrm{w}, \\alpha}\\right| \\boldsymbol{k}_{h, \\alpha}}{\\sum_\\alpha \\left|h_{\\mathrm{w}, \\alpha}\\right|}
```

is a weighted average of the orography's wavevectors ``\\boldsymbol{k}_{h, \\alpha}`` (with ``h_{\\mathrm{w}, \\alpha}`` being the corresponding spectral modes) and

```math
\\boldsymbol{u}_\\mathrm{p} = \\frac{\\left(\\boldsymbol{u}_\\mathrm{b} \\cdot \\boldsymbol{k}_h\\right) \\boldsymbol{k}_h}{\\left|\\boldsymbol{k}_h\\right|^2}
```

is the projection of ``\\boldsymbol{u}_\\mathrm{b}`` onto ``\\boldsymbol{k}_h``. This drag replaces the drag due to gravity waves below ``z_\\mathrm{B}``, the upper edge of the blocked layer that has been determined by [`PinCFlow.MSGWaM.RaySources.activate_orographic_source!`](@ref). In grid cells that contain this upper edge, blocking and gravity waves both contribute to the total drag, weighted by the corresponding grid-cell fractions. The gravity-wave heating is treated similarly, with the "blocking contribution" being zero.

# Arguments

  - `state`: Complete simulation state
  - `testcase`: Test case on which the current simulation is based.
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
        if fraction <= 0
            continue
        else
            @views kavg[1] =
                sum(
                    abs.(topography_spectrum[:, ix, jy]) .*
                    k_spectrum[:, ix, jy],
                ) / sum(abs.(topography_spectrum[:, ix, jy]))
            @views kavg[2] =
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
