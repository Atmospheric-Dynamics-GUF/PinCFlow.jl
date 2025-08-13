"""
```julia
compute_stress_tensor(
    i::Integer,
    j::Integer,
    k::Integer,
    mu::Integer,
    nu::Integer,
    predictands::Predictands,
    grid::Grid,
)
```

Compute the element `(mu, nu)` of the Cartesian viscous stress tensor at the grid point `(i, j, k)`, divided by the dynamic viscosity.

The discretized elements of the Cartesian viscous stress tensor are given by

```math
\\begin{align*}
    \\Pi^{1 1} & = \\eta \\left[\\frac{2}{\\Delta \\widehat{x}} \\left(u_{i + 1 / 2} - u_{i - 1 / 2}\\right) + \\frac{G^{1 3}}{\\Delta \\widehat{z}} \\left(u_{k + 1} - u_{k - 1}\\right) - \\frac{2}{3} \\delta\\right],\\\\
    \\Pi^{1 2} & = \\eta \\left[\\frac{1}{2 \\Delta \\widehat{y}} \\left(u_{j + 1} - u_{j - 1}\\right) + \\frac{G^{2 3}}{2 \\Delta \\widehat{z}} \\left(u_{k + 1} - u_{k - 1}\\right) + \\frac{1}{2 \\Delta \\widehat{x}} \\left(v_{i + 1} - v_{i - 1}\\right) + \\frac{G^{1 3}}{2 \\Delta \\widehat{z}} \\left(v_{k + 1} - v_{k - 1}\\right)\\right],\\\\
    \\Pi^{1 3} & = \\eta \\left[\\frac{1}{2 J \\Delta \\widehat{z}} \\left(u_{k + 1} - u_{k - 1}\\right) + \\frac{1}{2 \\Delta \\widehat{x}} \\left(w_{i + 1} - w_{i - 1}\\right) + \\frac{G^{1 3}}{\\Delta \\widehat{z}} \\left(w_{k + 1 / 2} - w_{k - 1 / 2}\\right)\\right],\\\\
    \\Pi^{2 2} & = \\eta \\left[\\frac{2}{\\Delta \\widehat{y}} \\left(v_{j + 1 / 2} - v_{j - 1 / 2}\\right) + \\frac{G^{2 3}}{\\Delta \\widehat{z}} \\left(v_{k + 1} - v_{k - 1}\\right) - \\frac{2}{3} \\delta\\right],\\\\
    \\Pi^{2 3} & = \\eta \\left[\\frac{1}{2 J \\Delta \\widehat{z}} \\left(v_{k + 1} - v_{k - 1}\\right) + \\frac{1}{2 \\Delta \\widehat{y}} \\left(w_{j + 1} - w_{j - 1}\\right) + \\frac{G^{2 3}}{\\Delta \\widehat{z}} \\left(w_{k + 1 / 2} - w_{k - 1 / 2}\\right)\\right],\\\\
    \\Pi^{3 3} & = \\eta \\left[\\frac{2}{J \\Delta \\widehat{z}} \\left(w_{k + 1 / 2} - w_{k - 1 / 2}\\right) - \\frac{2}{3} \\delta\\right],
\\end{align*}
```

where

```math
\\delta = \\frac{1}{J} \\left[\\frac{1}{\\Delta \\widehat{x}} \\left(J_{i + 1 / 2} u_{i + 1 / 2} - J_{i - 1 / 2} u_{i - 1 / 2}\\right) + \\frac{1}{\\Delta \\widehat{y}} \\left(J_{j + 1 / 2} v_{j + 1 / 2} - J_{j - 1 / 2} v_{j - 1 / 2}\\right) + \\frac{1}{\\Delta \\widehat{z}} \\left(J_{k + 1 / 2} \\widehat{w}_{k + 1 / 2} - J_{k - 1 / 2} \\widehat{w}_{k - 1 / 2}\\right)\\right].
```

# Arguments

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

  - `mu`: First contravariant tensor index.

  - `nu`: Second contravariant tensor index.

  - `predictands`: Prognostic variables.

  - `grid`: Collection of parameters and fields that describe the grid.

# Returns

  - `::AbstractFloat`: Element `(mu, nu)` of the Cartesian stress tensor at `(i, j, k)`.

# See also

  - [`PinCFlow.Update.compute_vertical_wind`](@ref)
"""
function compute_stress_tensor end

function compute_stress_tensor(
    i::Integer,
    j::Integer,
    k::Integer,
    mu::Integer,
    nu::Integer,
    predictands::Predictands,
    grid::Grid,
)
    (; u, v, w) = predictands
    (; dx, dy, dz, jac, met) = grid

    jacedger = 0.5 * (jac[i, j, k] + jac[i + 1, j, k])
    jacedgel = 0.5 * (jac[i, j, k] + jac[i - 1, j, k])
    jacedgef = 0.5 * (jac[i, j, k] + jac[i, j + 1, k])
    jacedgeb = 0.5 * (jac[i, j, k] + jac[i, j - 1, k])
    jacedgeu =
        2.0 * jac[i, j, k] * jac[i, j, k + 1] /
        (jac[i, j, k] + jac[i, j, k + 1])
    jacedged =
        2.0 * jac[i, j, k] * jac[i, j, k - 1] /
        (jac[i, j, k] + jac[i, j, k - 1])

    uf = 0.5 * (u[i, j + 1, k] + u[i - 1, j + 1, k])
    ub = 0.5 * (u[i, j - 1, k] + u[i - 1, j - 1, k])
    uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1])
    ud = 0.5 * (u[i, j, k - 1] + u[i - 1, j, k - 1])
    vr = 0.5 * (v[i + 1, j, k] + v[i + 1, j - 1, k])
    vl = 0.5 * (v[i - 1, j, k] + v[i - 1, j - 1, k])
    vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1])
    vd = 0.5 * (v[i, j, k - 1] + v[i, j - 1, k - 1])
    wr =
        0.5 * (
            compute_vertical_wind(i + 1, j, k, predictands, grid) +
            compute_vertical_wind(i + 1, j, k - 1, predictands, grid)
        )
    wl =
        0.5 * (
            compute_vertical_wind(i - 1, j, k, predictands, grid) +
            compute_vertical_wind(i - 1, j, k - 1, predictands, grid)
        )
    wf =
        0.5 * (
            compute_vertical_wind(i, j + 1, k, predictands, grid) +
            compute_vertical_wind(i, j + 1, k - 1, predictands, grid)
        )
    wb =
        0.5 * (
            compute_vertical_wind(i, j - 1, k, predictands, grid) +
            compute_vertical_wind(i, j - 1, k - 1, predictands, grid)
        )

    if mu == 1 && nu == 1
        stress_tensor =
            2.0 * (u[i, j, k] - u[i - 1, j, k]) / dx +
            met[i, j, k, 1, 3] * (uu - ud) / dz -
            2.0 / 3.0 * (
                (jacedger * u[i, j, k] - jacedgel * u[i - 1, j, k]) / dx +
                (jacedgef * v[i, j, k] - jacedgeb * v[i, j - 1, k]) / dy +
                (jacedgeu * w[i, j, k] - jacedged * w[i, j, k - 1]) / dz
            ) / jac[i, j, k]
    elseif (mu == 1 && nu == 2) || (mu == 2 && nu == 1)
        stress_tensor =
            0.5 * (uf - ub) / dy +
            0.5 * met[i, j, k, 2, 3] * (uu - ud) / dz +
            0.5 * (vr - vl) / dx +
            0.5 * met[i, j, k, 1, 3] * (vu - vd) / dz
    elseif (mu == 1 && nu == 3) || (mu == 3 && nu == 1)
        stress_tensor =
            0.5 * (uu - ud) / dz / jac[i, j, k] +
            0.5 * (wr - wl) / dx +
            met[i, j, k, 1, 3] * (
                compute_vertical_wind(i, j, k, predictands, grid) -
                compute_vertical_wind(i, j, k - 1, predictands, grid)
            ) / dz
    elseif mu == 2 && nu == 2
        stress_tensor =
            2.0 * (v[i, j, k] - v[i, j - 1, k]) / dy +
            met[i, j, k, 2, 3] * (vu - vd) / dz -
            2.0 / 3.0 * (
                (jacedger * u[i, j, k] - jacedgel * u[i - 1, j, k]) / dx +
                (jacedgef * v[i, j, k] - jacedgeb * v[i, j - 1, k]) / dy +
                (jacedgeu * w[i, j, k] - jacedged * w[i, j, k - 1]) / dz
            ) / jac[i, j, k]
    elseif (mu == 2 && nu == 3) || (mu == 3 && nu == 2)
        stress_tensor =
            0.5 * (vu - vd) / dz / jac[i, j, k] +
            0.5 * (wf - wb) / dy +
            met[i, j, k, 2, 3] * (
                compute_vertical_wind(i, j, k, predictands, grid) -
                compute_vertical_wind(i, j, k - 1, predictands, grid)
            ) / dz
    elseif mu == 3 && nu == 3
        stress_tensor =
            2.0 * (
                compute_vertical_wind(i, j, k, predictands, grid) -
                compute_vertical_wind(i, j, k - 1, predictands, grid)
            ) / dz / jac[i, j, k] -
            2.0 / 3.0 * (
                (jacedger * u[i, j, k] - jacedgel * u[i - 1, j, k]) / dx +
                (jacedgef * v[i, j, k] - jacedgeb * v[i, j - 1, k]) / dy +
                (jacedgeu * w[i, j, k] - jacedged * w[i, j, k - 1]) / dz
            ) / jac[i, j, k]
    end

    return stress_tensor
end
