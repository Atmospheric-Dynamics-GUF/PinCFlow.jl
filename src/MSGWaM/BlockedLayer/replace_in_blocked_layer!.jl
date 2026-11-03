"""
```julia
replace_in_blocked_layer!(
    replacement::Function,
    state::State,
    array::AbstractArray{<:AbstractFloat, 3},
)::Nothing
```

In each grid cell `(i, j, k)`, replace the corresponding element in `array` with a weighted average between its value and `replacement(state, i, j, k)`, where the weights are given by the grid-cell fractions that are below and above the upper edge of the blocked layer.

The operation can be written as

```math
a \\rightarrow \\frac{\\delta z_\\mathrm{B}}{J \\Delta \\hat{z}} r + \\left(1 - \\frac{\\delta z_\\mathrm{B}}{J \\Delta \\hat{z}}\\right) a,
```

where ``a`` and ``r`` represent `array[i, j, k]` and `replacement(state, i, j, k)`, respectively, and ``\\delta z_\\mathrm{B}`` is defined by

```math
\\delta z_\\mathrm{B} = \\max \\left[0, \\min \\left(h_\\mathrm{b} + \\frac{\\Delta z_\\mathrm{B}}{2}, z_{k + 1 / 2}\\right) - z_{k - 1 / 2}\\right].
```

# Arguments

  - `state`: Model state.

# See also

  - [`PinCFlow.MSGWaM.BlockedLayer.compute_slope`](@ref)
"""
function replace_in_blocked_layer! end

@ivy function replace_in_blocked_layer!(
    replacement::Function,
    state::State,
    array::AbstractArray{<:AbstractFloat, 3},
)::Nothing
    (; blocking) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dz, jac, zctilde, hb) = state.grid
    (; deltazb) = state.wkb

    if !blocking
        return nothing
    end

    for k in k0:k1, j in j0:j1, i in i0:i1
        fraction =
            max(
                0,
                min(hb[i, j] + deltazb[i, j] / 2, zctilde[i, j, k]) -
                zctilde[i, j, k - 1],
            ) / jac[i, j, k] / dz
        if fraction > 0
            array[i, j, k] =
                fraction * replacement(state, i, j, k) +
                (1 - fraction) * array[i, j, k]
        end
    end

    nothing
end
