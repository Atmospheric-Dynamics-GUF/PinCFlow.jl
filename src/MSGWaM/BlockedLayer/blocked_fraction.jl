"""
```julia
blocked_fraction(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
)::AbstractFloat
```

Return the blocked fraction of the grid cell at `(i, j, k)`.

The blocked fraction is given by ``\\delta z_\\mathrm{B} / \\left(J \\Delta \\hat{z}\\right)``, where ``\\delta z_\\mathrm{B}`` is defined by

```math
\\delta z_\\mathrm{B} = \\max \\left[0, \\min \\left(h_\\mathrm{b} + \\frac{\\Delta z_\\mathrm{B}}{2}, z_{k + 1 / 2}\\right) - z_{k - 1 / 2}\\right].
```

If `state.namelists.wkb.blocking` is set to `false`, ``0`` is returned instead.

# Arguments

  - `state`: Model state.

  - `i`: Zonal index.

  - `j`: Meridional index.

  - `k`: Vertical index.
"""
function blocked_fraction end

@ivy function blocked_fraction(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
)::AbstractFloat
    (; blocking) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dz, jac, zctilde, hb) = state.grid
    (; deltazb) = state.wkb

    blocking ?
    max(
        0,
        min(hb[i, j] + deltazb[i, j] / 2, zctilde[i, j, k]) -
        zctilde[i, j, k - 1],
    ) / jac[i, j, k] / dz : 0.0
end
