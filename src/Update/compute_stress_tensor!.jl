"""
```julia
compute_stress_tensor!(state::State)
```

Compute the Cartesian stress tensor elements.

The discretized elements of the Cartesian stress tensor are given by

```math
\\begin{align*}
    \\Pi^{1 1} & = \\frac{2}{\\Delta \\hat{x}} \\left(u_{i + 1 / 2} - u_{i - 1 / 2}\\right) + \\frac{G^{1 3}}{\\Delta \\hat{z}} \\left(u_{k + 1} - u_{k - 1}\\right) - \\frac{2}{3} \\delta,\\\\
    \\Pi^{1 2} & = \\frac{1}{2 \\Delta \\hat{y}} \\left(u_{j + 1} - u_{j - 1}\\right) + \\frac{G^{2 3}}{2 \\Delta \\hat{z}} \\left(u_{k + 1} - u_{k - 1}\\right) + \\frac{1}{2 \\Delta \\hat{x}} \\left(v_{i + 1} - v_{i - 1}\\right) + \\frac{G^{1 3}}{2 \\Delta \\hat{z}} \\left(v_{k + 1} - v_{k - 1}\\right),\\\\
    \\Pi^{1 3} & = \\frac{1}{2 J \\Delta \\hat{z}} \\left(u_{k + 1} - u_{k - 1}\\right) + \\frac{1}{2 \\Delta \\hat{x}} \\left(w_{i + 1} - w_{i - 1}\\right) + \\frac{G^{1 3}}{\\Delta \\hat{z}} \\left(w_{k + 1 / 2} - w_{k - 1 / 2}\\right),\\\\
    \\Pi^{2 2} & = \\frac{2}{\\Delta \\hat{y}} \\left(v_{j + 1 / 2} - v_{j - 1 / 2}\\right) + \\frac{G^{2 3}}{\\Delta \\hat{z}} \\left(v_{k + 1} - v_{k - 1}\\right) - \\frac{2}{3} \\delta,\\\\
    \\Pi^{2 3} & = \\frac{1}{2 J \\Delta \\hat{z}} \\left(v_{k + 1} - v_{k - 1}\\right) + \\frac{1}{2 \\Delta \\hat{y}} \\left(w_{j + 1} - w_{j - 1}\\right) + \\frac{G^{2 3}}{\\Delta \\hat{z}} \\left(w_{k + 1 / 2} - w_{k - 1 / 2}\\right),\\\\
    \\Pi^{3 3} & = \\frac{2}{J \\Delta \\hat{z}} \\left(w_{k + 1 / 2} - w_{k - 1 / 2}\\right) - \\frac{2}{3} \\delta,
\\end{align*}
```

where

```math
\\begin{align*}
    \\delta & = \\frac{1}{J} \\left[\\frac{1}{\\Delta \\hat{x}} \\left(J_{i + 1 / 2} u_{i + 1 / 2} - J_{i - 1 / 2} u_{i - 1 / 2}\\right) + \\frac{1}{\\Delta \\hat{y}} \\left(J_{j + 1 / 2} v_{j + 1 / 2} - J_{j - 1 / 2} v_{j - 1 / 2}\\right)\\right.\\\\
    & \\qquad \\quad + \\left.\\frac{1}{\\Delta \\hat{z}} \\left(J_{k + 1 / 2} \\hat{w}_{k + 1 / 2} - J_{k - 1 / 2} \\hat{w}_{k - 1 / 2}\\right)\\right].
\\end{align*}
```

``\\Pi^{1 1}``, ``\\Pi^{1 2}``, ``\\Pi^{1 3}``, ``\\Pi^{2 2}``, ``\\Pi^{2 3}``, and ``\\Pi^{3 3}`` are stored in `state.variables.auxiliaries.stress_tensor_11`, `state.variables.auxiliaries.stress_tensor_12`, `state.variables.auxiliaries.stress_tensor_13`, `state.variables.auxiliaries.stress_tensor_22`, `state.variables.auxiliaries.stress_tensor_23`, and `state.variables.auxiliaries.stress_tensor_33`, respectively.

# Arguments

  - `state`: Model state.

# See also

  - [`PinCFlow.Update.compute_vertical_wind`](@ref)
"""
function compute_stress_tensor! end

@ivy function compute_stress_tensor!(state::State)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (;
        stress_tensor_11,
        stress_tensor_12,
        stress_tensor_13,
        stress_tensor_22,
        stress_tensor_23,
        stress_tensor_33,
    ) = state.variables.auxiliaries
    (; u, v, w) = state.variables.predictands
    (; dx, dy, dz, jac, met) = state.grid
    (; re) = state.constants

    if 1 / re <= eps()
        return
    end

    for k in (k0 - 2):(k1 + 1), j in (j0 - 2):(j1 + 1), i in (i0 - 2):(i1 + 1)
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
                compute_vertical_wind(i + 1, j, k, state) +
                compute_vertical_wind(i + 1, j, k - 1, state)
            )
        wl =
            0.5 * (
                compute_vertical_wind(i - 1, j, k, state) +
                compute_vertical_wind(i - 1, j, k - 1, state)
            )
        wf =
            0.5 * (
                compute_vertical_wind(i, j + 1, k, state) +
                compute_vertical_wind(i, j + 1, k - 1, state)
            )
        wb =
            0.5 * (
                compute_vertical_wind(i, j - 1, k, state) +
                compute_vertical_wind(i, j - 1, k - 1, state)
            )

        stress_tensor_11[i, j, k] =
            2.0 * (u[i, j, k] - u[i - 1, j, k]) / dx +
            met[i, j, k, 1, 3] * (uu - ud) / dz -
            2.0 / 3.0 * (
                (jacedger * u[i, j, k] - jacedgel * u[i - 1, j, k]) / dx +
                (jacedgef * v[i, j, k] - jacedgeb * v[i, j - 1, k]) / dy +
                (jacedgeu * w[i, j, k] - jacedged * w[i, j, k - 1]) / dz
            ) / jac[i, j, k]

        stress_tensor_12[i, j, k] =
            0.5 * (uf - ub) / dy +
            0.5 * met[i, j, k, 2, 3] * (uu - ud) / dz +
            0.5 * (vr - vl) / dx +
            0.5 * met[i, j, k, 1, 3] * (vu - vd) / dz

        stress_tensor_13[i, j, k] =
            0.5 * (uu - ud) / dz / jac[i, j, k] +
            0.5 * (wr - wl) / dx +
            met[i, j, k, 1, 3] * (
                compute_vertical_wind(i, j, k, state) -
                compute_vertical_wind(i, j, k - 1, state)
            ) / dz

        stress_tensor_22[i, j, k] =
            2.0 * (v[i, j, k] - v[i, j - 1, k]) / dy +
            met[i, j, k, 2, 3] * (vu - vd) / dz -
            2.0 / 3.0 * (
                (jacedger * u[i, j, k] - jacedgel * u[i - 1, j, k]) / dx +
                (jacedgef * v[i, j, k] - jacedgeb * v[i, j - 1, k]) / dy +
                (jacedgeu * w[i, j, k] - jacedged * w[i, j, k - 1]) / dz
            ) / jac[i, j, k]

        stress_tensor_23[i, j, k] =
            0.5 * (vu - vd) / dz / jac[i, j, k] +
            0.5 * (wf - wb) / dy +
            met[i, j, k, 2, 3] * (
                compute_vertical_wind(i, j, k, state) -
                compute_vertical_wind(i, j, k - 1, state)
            ) / dz

        stress_tensor_33[i, j, k] =
            2.0 * (
                compute_vertical_wind(i, j, k, state) -
                compute_vertical_wind(i, j, k - 1, state)
            ) / dz / jac[i, j, k] -
            2.0 / 3.0 * (
                (jacedger * u[i, j, k] - jacedgel * u[i - 1, j, k]) / dx +
                (jacedgef * v[i, j, k] - jacedgeb * v[i, j - 1, k]) / dy +
                (jacedgeu * w[i, j, k] - jacedged * w[i, j, k - 1]) / dz
            ) / jac[i, j, k]
    end

    return
end