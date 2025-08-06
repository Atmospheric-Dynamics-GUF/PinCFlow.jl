"""
```julia
compute_derivatives(state::State, indices::NTuple{4, <:Integer}, phitype::DUDX)
```

Compute the zonal derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial x``) at two specified positions on the grid.

The derivative is given by

```math
\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial x}\\right) = \\frac{u_{\\mathrm{b}, i + 1 / 2} - u_{\\mathrm{b}, i - 1 / 2}}{\\Delta \\widehat{x}} + G^{13} \\frac{u_{\\mathrm{b}, i + 1 / 2, k + 1} + u_{\\mathrm{b}, i - 1 / 2, k + 1} - u_{\\mathrm{b}, i + 1 / 2, k - 1} - u_{\\mathrm{b}, i - 1 / 2, k - 1}}{4 \\Delta \\widehat{z}}.
```

```julia
compute_derivatives(state::State, indices::NTuple{4, <:Integer}, phitype::DUDY)
```

Compute the meridional derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial y``) at two specified positions on the grid.

The derivative is given by

```math
\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial y}\\right)_{i + 1 / 2, j + 1 / 2} = \\frac{u_{\\mathrm{b}, i + 1 / 2, j + 1} - u_{\\mathrm{b}, i + 1 / 2}}{\\Delta \\widehat{y}} + \\frac{G^{23} + G^{23}_{i + 1} + G^{23}_{j + 1} + G^{23}_{i + 1, j + 1}}{4} \\frac{u_{\\mathrm{b}, i + 1 / 2, k + 1} + u_{\\mathrm{b}, i + 1 / 2, j + 1, k + 1} - u_{\\mathrm{b}, i + 1 / 2, k - 1} - u_{\\mathrm{b}, i + 1 / 2, j + 1, k - 1}}{4 \\Delta \\widehat{z}}.
```

```julia
compute_derivatives(state::State, indices::NTuple{4, <:Integer}, phitype::DUDZ)
```

Compute the vertical derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial z``) at two specified positions on the grid.

The derivative is given by

```math
\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial z}\\right)_{i + 1 / 2, k + 1 / 2} = \\frac{u_{\\mathrm{b}, i + 1 / 2, k + 1} - u_{\\mathrm{b}, i + 1 / 2}}{\\Delta \\widehat{z}} \\left(\\frac{J J_{k + 1}}{J + J_{k + 1}} + \\frac{J_{i + 1} J_{i + 1, k + 1}}{J_{i + 1} + J_{i + 1, k + 1}}\\right)^{- 1}.
```

At grid points beyond the vertical boundaries, it is set to zero.

```julia
compute_derivatives(state::State, indices::NTuple{4, <:Integer}, phitype::DVDX)
```

Compute the zonal derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial x``) at two specified positions on the grid.

The derivative is given by

```math
\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial x}\\right)_{i + 1 / 2, j + 1 / 2} = \\frac{v_{\\mathrm{b}, i + 1, j + 1 / 2} - v_{\\mathrm{b}, j + 1 / 2}}{\\Delta \\widehat{x}} + \\frac{G^{13} + G^{13}_{i + 1} + G^{13}_{j + 1} + G^{13}_{i + 1, j + 1}}{4} \\frac{v_{\\mathrm{b}, j + 1 / 2, k + 1} + v_{\\mathrm{b}, i + 1, j + 1 / 2, k + 1} - v_{\\mathrm{b}, j + 1 / 2, k - 1} - v_{\\mathrm{b}, i + 1, j + 1 / 2, k - 1}}{4 \\Delta \\widehat{z}}.
```

```julia
compute_derivatives(state::State, indices::NTuple{4, <:Integer}, phitype::DVDY)
```

Compute the meridional derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial y``) at two specified positions on the grid.

The derivative is given by

```math
\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial y}\\right) = \\frac{v_{\\mathrm{b}, j + 1 / 2} - v_{\\mathrm{b}, j - 1 / 2}}{\\Delta \\widehat{y}} + G^{23} \\frac{v_{\\mathrm{b}, j + 1 / 2, k + 1} + v_{\\mathrm{b}, j - 1 / 2, k + 1} - v_{\\mathrm{b}, j + 1 / 2, k - 1} - v_{\\mathrm{b}, j - 1 / 2, k - 1}}{4 \\Delta \\widehat{z}}.
```

```julia
compute_derivatives(state::State, indices::NTuple{4, <:Integer}, phitype::DVDZ)
```

Compute the vertical derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial z``) at two specified positions on the grid.

The derivative is given by

```math
\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial z}\\right)_{j + 1 / 2, k + 1 / 2} = \\frac{v_{\\mathrm{b}, j + 1 / 2, k + 1} - v_{\\mathrm{b}, j + 1 / 2}}{\\Delta \\widehat{z}} \\left(\\frac{J J_{k + 1}}{J + J_{k + 1}} + \\frac{J_{j + 1} J_{j + 1, k + 1}}{J_{j + 1} + J_{j + 1, k + 1}}\\right)^{- 1}.
```

At grid points beyond the vertical boundaries, it is set to zero.

# Arguments

  - `state`: Model state.
  - `indices`: Grid indices `(ix, jy, kzd, kzu)` of the two positions at which to compute the derivative.
  - `phitype`: Type of derivative to compute.

# Returns

  - `::AbstractFloat`: Specified derivative near the grid point `(ix, jy, kzd)`.
  - `::AbstractFloat`: Specified derivative near the grid point `(ix, jy, kzu)`.
"""
function compute_derivatives end

function compute_derivatives(
    state::State,
    indices::NTuple{4, <:Integer},
    phitype::DUDX,
)
    (; dx, dz, met) = state.grid
    (; u) = state.variables.predictands

    (ix, jy, kzd, kzu) = indices

    phid =
        (u[ix, jy, kzd] - u[ix - 1, jy, kzd]) / dx +
        met[ix, jy, kzd, 1, 3] *
        0.25 *
        (
            u[ix, jy, kzd + 1] + u[ix - 1, jy, kzd + 1] - u[ix, jy, kzd - 1] -
            u[ix - 1, jy, kzd - 1]
        ) / dz
    phiu =
        (u[ix, jy, kzu] - u[ix - 1, jy, kzu]) / dx +
        met[ix, jy, kzu, 1, 3] *
        0.25 *
        (
            u[ix, jy, kzu + 1] + u[ix - 1, jy, kzu + 1] - u[ix, jy, kzu - 1] -
            u[ix - 1, jy, kzu - 1]
        ) / dz

    return (phid, phiu)
end

function compute_derivatives(
    state::State,
    indices::NTuple{4, <:Integer},
    phitype::DUDY,
)
    (; dy, dz, met) = state.grid
    (; u) = state.variables.predictands

    (ix, jy, kzd, kzu) = indices

    phid =
        (u[ix, jy + 1, kzd] - u[ix, jy, kzd]) / dy +
        0.25 *
        (
            met[ix, jy, kzd, 2, 3] +
            met[ix + 1, jy, kzd, 2, 3] +
            met[ix, jy + 1, kzd, 2, 3] +
            met[ix + 1, jy + 1, kzd, 2, 3]
        ) *
        0.25 *
        (
            u[ix, jy, kzd + 1] + u[ix, jy + 1, kzd + 1] - u[ix, jy, kzd - 1] -
            u[ix, jy + 1, kzd - 1]
        ) / dz
    phiu =
        (u[ix, jy + 1, kzu] - u[ix, jy, kzu]) / dy +
        0.25 *
        (
            met[ix, jy, kzu, 2, 3] +
            met[ix + 1, jy, kzu, 2, 3] +
            met[ix, jy + 1, kzu, 2, 3] +
            met[ix + 1, jy + 1, kzu, 2, 3]
        ) *
        0.25 *
        (
            u[ix, jy, kzu + 1] + u[ix, jy + 1, kzu + 1] - u[ix, jy, kzu - 1] -
            u[ix, jy + 1, kzu - 1]
        ) / dz

    return (phid, phiu)
end

function compute_derivatives(
    state::State,
    indices::NTuple{4, <:Integer},
    phitype::DUDZ,
)
    (; lz, dz, ztildetfc, jac, topography_surface) = state.grid
    (; u) = state.variables.predictands

    (ix, jy, kzd, kzu) = indices

    if ztildetfc[ix, jy, kzu] < topography_surface[ix, jy]
        phid = 0.0
        phiu = 0.0
    elseif ztildetfc[ix, jy, kzd] < topography_surface[ix, jy]
        phid = 0.0
        phiu =
            (u[ix, jy, kzu + 1] - u[ix, jy, kzu]) / dz / (
                jac[ix, jy, kzu] * jac[ix, jy, kzu + 1] /
                (jac[ix, jy, kzu] + jac[ix, jy, kzu + 1]) +
                jac[ix + 1, jy, kzu] * jac[ix + 1, jy, kzu + 1] /
                (jac[ix + 1, jy, kzu] + jac[ix + 1, jy, kzu + 1])
            )
    else
        if ztildetfc[ix, jy, kzu] < lz[2]
            phid =
                (u[ix, jy, kzd + 1] - u[ix, jy, kzd]) / dz / (
                    jac[ix, jy, kzd] * jac[ix, jy, kzd + 1] /
                    (jac[ix, jy, kzd] + jac[ix, jy, kzd + 1]) +
                    jac[ix + 1, jy, kzd] * jac[ix + 1, jy, kzd + 1] /
                    (jac[ix + 1, jy, kzd] + jac[ix + 1, jy, kzd + 1])
                )
            phiu =
                (u[ix, jy, kzu + 1] - u[ix, jy, kzu]) / dz / (
                    jac[ix, jy, kzu] * jac[ix, jy, kzu + 1] /
                    (jac[ix, jy, kzu] + jac[ix, jy, kzu + 1]) +
                    jac[ix + 1, jy, kzu] * jac[ix + 1, jy, kzu + 1] /
                    (jac[ix + 1, jy, kzu] + jac[ix + 1, jy, kzu + 1])
                )
        elseif ztildetfc[ix, jy, kzd] < lz[2]
            phid =
                (u[ix, jy, kzd + 1] - u[ix, jy, kzd]) / dz / (
                    jac[ix, jy, kzd] * jac[ix, jy, kzd + 1] /
                    (jac[ix, jy, kzd] + jac[ix, jy, kzd + 1]) +
                    jac[ix + 1, jy, kzd] * jac[ix + 1, jy, kzd + 1] /
                    (jac[ix + 1, jy, kzd] + jac[ix + 1, jy, kzd + 1])
                )
            phiu = 0.0
        else
            phid = 0.0
            phiu = 0.0
        end
    end

    return (phid, phiu)
end

function compute_derivatives(
    state::State,
    indices::NTuple{4, <:Integer},
    phitype::DVDX,
)
    (; dx, dz, met) = state.grid
    (; v) = state.variables.predictands

    (ix, jy, kzd, kzu) = indices

    phid =
        (v[ix + 1, jy, kzd] - v[ix, jy, kzd]) / dx +
        0.25 *
        (
            met[ix, jy, kzd, 1, 3] +
            met[ix + 1, jy, kzd, 1, 3] +
            met[ix, jy + 1, kzd, 1, 3] +
            met[ix + 1, jy + 1, kzd, 1, 3]
        ) *
        0.25 *
        (
            v[ix, jy, kzd + 1] + v[ix + 1, jy, kzd + 1] - v[ix, jy, kzd - 1] -
            v[ix + 1, jy, kzd - 1]
        ) / dz
    phiu =
        (v[ix + 1, jy, kzu] - v[ix, jy, kzu]) / dx +
        0.25 *
        (
            met[ix, jy, kzu, 1, 3] +
            met[ix + 1, jy, kzu, 1, 3] +
            met[ix, jy + 1, kzu, 1, 3] +
            met[ix + 1, jy + 1, kzu, 1, 3]
        ) *
        0.25 *
        (
            v[ix, jy, kzu + 1] + v[ix + 1, jy, kzu + 1] - v[ix, jy, kzu - 1] -
            v[ix + 1, jy, kzu - 1]
        ) / dz

    return (phid, phiu)
end

function compute_derivatives(
    state::State,
    indices::NTuple{4, <:Integer},
    phitype::DVDY,
)
    (; dy, dz, met) = state.grid
    (; v) = state.variables.predictands

    (ix, jy, kzd, kzu) = indices

    phid =
        (v[ix, jy, kzd] - v[ix, jy - 1, kzd]) / dy +
        met[ix, jy, kzd, 2, 3] *
        0.25 *
        (
            v[ix, jy, kzd + 1] + v[ix, jy - 1, kzd + 1] - v[ix, jy, kzd - 1] -
            v[ix, jy - 1, kzd - 1]
        ) / dz
    phiu =
        (v[ix, jy, kzu] - v[ix, jy - 1, kzu]) / dy +
        met[ix, jy, kzu, 2, 3] *
        0.25 *
        (
            v[ix, jy, kzu + 1] + v[ix, jy - 1, kzu + 1] - v[ix, jy, kzu - 1] -
            v[ix, jy - 1, kzu - 1]
        ) / dz

    return (phid, phiu)
end

function compute_derivatives(
    state::State,
    indices::NTuple{4, <:Integer},
    phitype::DVDZ,
)
    (; lz, dz, ztildetfc, jac, topography_surface) = state.grid
    (; v) = state.variables.predictands

    (ix, jy, kzd, kzu) = indices

    if ztildetfc[ix, jy, kzu] < topography_surface[ix, jy]
        phid = 0.0
        phiu = 0.0
    elseif ztildetfc[ix, jy, kzd] < topography_surface[ix, jy]
        phid = 0.0
        phiu =
            (v[ix, jy, kzu + 1] - v[ix, jy, kzu]) / dz / (
                jac[ix, jy, kzu] * jac[ix, jy, kzu + 1] /
                (jac[ix, jy, kzu] + jac[ix, jy, kzu + 1]) +
                jac[ix, jy + 1, kzu] * jac[ix, jy + 1, kzu + 1] /
                (jac[ix, jy + 1, kzu] + jac[ix, jy + 1, kzu + 1])
            )
    else
        if ztildetfc[ix, jy, kzu] < lz[2]
            phid =
                (v[ix, jy, kzd + 1] - v[ix, jy, kzd]) / dz / (
                    jac[ix, jy, kzd] * jac[ix, jy, kzd + 1] /
                    (jac[ix, jy, kzd] + jac[ix, jy, kzd + 1]) +
                    jac[ix, jy + 1, kzd] * jac[ix, jy + 1, kzd + 1] /
                    (jac[ix, jy + 1, kzd] + jac[ix, jy + 1, kzd + 1])
                )
            phiu =
                (v[ix, jy, kzu + 1] - v[ix, jy, kzu]) / dz / (
                    jac[ix, jy, kzu] * jac[ix, jy, kzu + 1] /
                    (jac[ix, jy, kzu] + jac[ix, jy, kzu + 1]) +
                    jac[ix, jy + 1, kzu] * jac[ix, jy + 1, kzu + 1] /
                    (jac[ix, jy + 1, kzu] + jac[ix, jy + 1, kzu + 1])
                )
        elseif ztildetfc[ix, jy, kzd] < lz[2]
            phid =
                (v[ix, jy, kzd + 1] - v[ix, jy, kzd]) / dz / (
                    jac[ix, jy, kzd] * jac[ix, jy, kzd + 1] /
                    (jac[ix, jy, kzd] + jac[ix, jy, kzd + 1]) +
                    jac[ix, jy + 1, kzd] * jac[ix, jy + 1, kzd + 1] /
                    (jac[ix, jy + 1, kzd] + jac[ix, jy + 1, kzd + 1])
                )
            phiu = 0.0
        else
            phid = 0.0
            phiu = 0.0
        end
    end

    return (phid, phiu)
end
