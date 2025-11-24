"""
```julia
apply_blocked_layer_scheme!(state::State)
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

  - `state`: Model state.
"""
function apply_blocked_layer_scheme! end

function apply_blocked_layer_scheme!(state::State)
    (; blocking, drag_coefficient) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dz, jac, zctilde, hw, kh, lh) = state.grid
    (; rhobar) = state.atmosphere
    (; rho, u, v) = state.variables.predictands
    (; zb) = state.wkb
    (; dudt, dvdt, dthetadt) = state.wkb.tendencies

    if !blocking
        return
    end

    # Initialize arrays for blocked-flow drag computation.
    (kavg, uperp, drag) = (zeros(2) for i in 1:3)

    # Adjust the drag to account for blocking.
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        fraction =
            (min(zb[i, j], zctilde[i, j, k]) - zctilde[i, j, k - 1]) /
            jac[i, j, k] / dz
        if fraction <= 0
            continue
        else
            hsum = sum(abs, hw[:, i, j])
            ksum = mapreduce((a, b) -> abs(a) * b, +, hw[:, i, j], kh[:, i, j])
            lsum = mapreduce((a, b) -> abs(a) * b, +, hw[:, i, j], lh[:, i, j])
            kavg[1] = ksum / hsum
            kavg[2] = lsum / hsum

            kavg2 = sum(a -> a^2, kavg)
            uperp .=
                (
                    (u[i, j, k] .+ u[i - 1, j, k]) .* kavg[1] .+
                    (v[i, j, k] .+ v[i, j - 1, k]) .* kavg[2]
                ) ./ 2 .* kavg ./ kavg2

            uperp2 = sum(a -> a^2, uperp)
            drag .=
                .-drag_coefficient .* (rho[i, j, k] .+ rhobar[i, j, k]) .*
                sqrt.(kavg2) ./ (2 .* pi) .* sqrt.(uperp2) .* uperp

            dudt[i, j, k] = fraction * drag[1] + (1 - fraction) * dudt[i, j, k]
            dvdt[i, j, k] = fraction * drag[2] + (1 - fraction) * dvdt[i, j, k]
            dthetadt[i, j, k] = (1 - fraction) * dthetadt[i, j, k]
        end
    end
end
