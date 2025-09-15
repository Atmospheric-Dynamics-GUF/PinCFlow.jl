"""
```julia
compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DUDX,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the zonal derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial x``) at ``\\left(i, j, k_\\mathrm{D}\\right)`` and ``\\left(i, j, k_\\mathrm{U}\\right)``.

The derivative is given by

```math
\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial x}\\right) = \\frac{u_{\\mathrm{b}, i + 1 / 2} - u_{\\mathrm{b}, i - 1 / 2}}{\\Delta \\widehat{x}} + G^{13} \\frac{u_{\\mathrm{b}, i + 1 / 2, k + 1} + u_{\\mathrm{b}, i - 1 / 2, k + 1} - u_{\\mathrm{b}, i + 1 / 2, k - 1} - u_{\\mathrm{b}, i - 1 / 2, k - 1}}{4 \\Delta \\widehat{z}}.
```

```julia
compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DUDY,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the meridional derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial y``) at ``\\left(i + 1 / 2, j + 1 / 2, k_\\mathrm{D}\\right)`` and ``\\left(i + 1 / 2, j + 1 / 2, k_\\mathrm{U}\\right)``.

The derivative is given by

```math
\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial y}\\right)_{i + 1 / 2, j + 1 / 2} = \\frac{u_{\\mathrm{b}, i + 1 / 2, j + 1} - u_{\\mathrm{b}, i + 1 / 2}}{\\Delta \\widehat{y}} + \\frac{G^{23} + G^{23}_{i + 1} + G^{23}_{j + 1} + G^{23}_{i + 1, j + 1}}{4} \\frac{u_{\\mathrm{b}, i + 1 / 2, k + 1} + u_{\\mathrm{b}, i + 1 / 2, j + 1, k + 1} - u_{\\mathrm{b}, i + 1 / 2, k - 1} - u_{\\mathrm{b}, i + 1 / 2, j + 1, k - 1}}{4 \\Delta \\widehat{z}}.
```

```julia
compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DUDZ,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the vertical derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial z``) at ``\\left(i + 1 / 2, j, k_\\mathrm{D} + 1 / 2\\right)`` and ``\\left(i + 1 / 2, j, k_\\mathrm{U} + 1 / 2\\right)``.

The derivative is given by

```math
\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial z}\\right)_{i + 1 / 2, k + 1 / 2} = \\frac{u_{\\mathrm{b}, i + 1 / 2, k + 1} - u_{\\mathrm{b}, i + 1 / 2}}{\\Delta \\widehat{z}} \\left(\\frac{J J_{k + 1}}{J + J_{k + 1}} + \\frac{J_{i + 1} J_{i + 1, k + 1}}{J_{i + 1} + J_{i + 1, k + 1}}\\right)^{- 1}.
```

At grid points beyond the vertical boundaries, it is set to zero.

```julia
compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DVDX,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the zonal derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial x``) at ``\\left(i + 1 / 2, j + 1 / 2, k_\\mathrm{D}\\right)`` and ``\\left(i + 1 / 2, j + 1 / 2, k_\\mathrm{U}\\right)``.

The derivative is given by

```math
\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial x}\\right)_{i + 1 / 2, j + 1 / 2} = \\frac{v_{\\mathrm{b}, i + 1, j + 1 / 2} - v_{\\mathrm{b}, j + 1 / 2}}{\\Delta \\widehat{x}} + \\frac{G^{13} + G^{13}_{i + 1} + G^{13}_{j + 1} + G^{13}_{i + 1, j + 1}}{4} \\frac{v_{\\mathrm{b}, j + 1 / 2, k + 1} + v_{\\mathrm{b}, i + 1, j + 1 / 2, k + 1} - v_{\\mathrm{b}, j + 1 / 2, k - 1} - v_{\\mathrm{b}, i + 1, j + 1 / 2, k - 1}}{4 \\Delta \\widehat{z}}.
```

```julia
compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DVDY,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the meridional derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial y``) at ``\\left(i, j, k_\\mathrm{D}\\right)`` and ``\\left(i, j, k_\\mathrm{U}\\right)``.

The derivative is given by

```math
\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial y}\\right) = \\frac{v_{\\mathrm{b}, j + 1 / 2} - v_{\\mathrm{b}, j - 1 / 2}}{\\Delta \\widehat{y}} + G^{23} \\frac{v_{\\mathrm{b}, j + 1 / 2, k + 1} + v_{\\mathrm{b}, j - 1 / 2, k + 1} - v_{\\mathrm{b}, j + 1 / 2, k - 1} - v_{\\mathrm{b}, j - 1 / 2, k - 1}}{4 \\Delta \\widehat{z}}.
```

```julia
compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DVDZ,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the vertical derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial z``) at ``\\left(i, j + 1 / 2, k_\\mathrm{D} + 1 / 2\\right)`` and ``\\left(i, j + 1 / 2, k_\\mathrm{U} + 1 / 2\\right)``.

The derivative is given by

```math
\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial z}\\right)_{j + 1 / 2, k + 1 / 2} = \\frac{v_{\\mathrm{b}, j + 1 / 2, k + 1} - v_{\\mathrm{b}, j + 1 / 2}}{\\Delta \\widehat{z}} \\left(\\frac{J J_{k + 1}}{J + J_{k + 1}} + \\frac{J_{j + 1} J_{j + 1, k + 1}}{J_{j + 1} + J_{j + 1, k + 1}}\\right)^{- 1}.
```

At grid points beyond the vertical boundaries, it is set to zero.

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `kd`: Lower vertical grid-cell index.

  - `ku`: Upper vertical grid-cell index.

  - `phitype`: Type of derivative to compute.
"""
function compute_derivatives end

function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DUDX,
)::NTuple{2, <:AbstractFloat}
    (; dx, dz, met) = state.grid
    (; u) = state.variables.predictands

    @ivy phid =
        (u[i, j, kd] - u[i - 1, j, kd]) / dx +
        met[i, j, kd, 1, 3] *
        0.25 *
        (
            u[i, j, kd + 1] + u[i - 1, j, kd + 1] - u[i, j, kd - 1] -
            u[i - 1, j, kd - 1]
        ) / dz
    @ivy phiu =
        (u[i, j, ku] - u[i - 1, j, ku]) / dx +
        met[i, j, ku, 1, 3] *
        0.25 *
        (
            u[i, j, ku + 1] + u[i - 1, j, ku + 1] - u[i, j, ku - 1] -
            u[i - 1, j, ku - 1]
        ) / dz

    return (phid, phiu)
end

function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DUDY,
)::NTuple{2, <:AbstractFloat}
    (; dy, dz, met) = state.grid
    (; u) = state.variables.predictands

    @ivy phid =
        (u[i, j + 1, kd] - u[i, j, kd]) / dy +
        0.25 *
        (
            met[i, j, kd, 2, 3] +
            met[i + 1, j, kd, 2, 3] +
            met[i, j + 1, kd, 2, 3] +
            met[i + 1, j + 1, kd, 2, 3]
        ) *
        0.25 *
        (
            u[i, j, kd + 1] + u[i, j + 1, kd + 1] - u[i, j, kd - 1] -
            u[i, j + 1, kd - 1]
        ) / dz
    @ivy phiu =
        (u[i, j + 1, ku] - u[i, j, ku]) / dy +
        0.25 *
        (
            met[i, j, ku, 2, 3] +
            met[i + 1, j, ku, 2, 3] +
            met[i, j + 1, ku, 2, 3] +
            met[i + 1, j + 1, ku, 2, 3]
        ) *
        0.25 *
        (
            u[i, j, ku + 1] + u[i, j + 1, ku + 1] - u[i, j, ku - 1] -
            u[i, j + 1, ku - 1]
        ) / dz

    return (phid, phiu)
end

function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DUDZ,
)::NTuple{2, <:AbstractFloat}
    (; lz, dz, ztildetfc, jac, topography_surface) = state.grid
    (; u) = state.variables.predictands

    @ivy if ztildetfc[i, j, ku] < topography_surface[i, j]
        phid = 0.0
        phiu = 0.0
    elseif ztildetfc[i, j, kd] < topography_surface[i, j]
        phid = 0.0
        phiu =
            (u[i, j, ku + 1] - u[i, j, ku]) / dz / (
                jac[i, j, ku] * jac[i, j, ku + 1] /
                (jac[i, j, ku] + jac[i, j, ku + 1]) +
                jac[i + 1, j, ku] * jac[i + 1, j, ku + 1] /
                (jac[i + 1, j, ku] + jac[i + 1, j, ku + 1])
            )
    else
        if ztildetfc[i, j, ku] < lz
            phid =
                (u[i, j, kd + 1] - u[i, j, kd]) / dz / (
                    jac[i, j, kd] * jac[i, j, kd + 1] /
                    (jac[i, j, kd] + jac[i, j, kd + 1]) +
                    jac[i + 1, j, kd] * jac[i + 1, j, kd + 1] /
                    (jac[i + 1, j, kd] + jac[i + 1, j, kd + 1])
                )
            phiu =
                (u[i, j, ku + 1] - u[i, j, ku]) / dz / (
                    jac[i, j, ku] * jac[i, j, ku + 1] /
                    (jac[i, j, ku] + jac[i, j, ku + 1]) +
                    jac[i + 1, j, ku] * jac[i + 1, j, ku + 1] /
                    (jac[i + 1, j, ku] + jac[i + 1, j, ku + 1])
                )
        elseif ztildetfc[i, j, kd] < lz
            phid =
                (u[i, j, kd + 1] - u[i, j, kd]) / dz / (
                    jac[i, j, kd] * jac[i, j, kd + 1] /
                    (jac[i, j, kd] + jac[i, j, kd + 1]) +
                    jac[i + 1, j, kd] * jac[i + 1, j, kd + 1] /
                    (jac[i + 1, j, kd] + jac[i + 1, j, kd + 1])
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
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DVDX,
)::NTuple{2, <:AbstractFloat}
    (; dx, dz, met) = state.grid
    (; v) = state.variables.predictands

    @ivy phid =
        (v[i + 1, j, kd] - v[i, j, kd]) / dx +
        0.25 *
        (
            met[i, j, kd, 1, 3] +
            met[i + 1, j, kd, 1, 3] +
            met[i, j + 1, kd, 1, 3] +
            met[i + 1, j + 1, kd, 1, 3]
        ) *
        0.25 *
        (
            v[i, j, kd + 1] + v[i + 1, j, kd + 1] - v[i, j, kd - 1] -
            v[i + 1, j, kd - 1]
        ) / dz
    @ivy phiu =
        (v[i + 1, j, ku] - v[i, j, ku]) / dx +
        0.25 *
        (
            met[i, j, ku, 1, 3] +
            met[i + 1, j, ku, 1, 3] +
            met[i, j + 1, ku, 1, 3] +
            met[i + 1, j + 1, ku, 1, 3]
        ) *
        0.25 *
        (
            v[i, j, ku + 1] + v[i + 1, j, ku + 1] - v[i, j, ku - 1] -
            v[i + 1, j, ku - 1]
        ) / dz

    return (phid, phiu)
end

function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DVDY,
)::NTuple{2, <:AbstractFloat}
    (; dy, dz, met) = state.grid
    (; v) = state.variables.predictands

    @ivy phid =
        (v[i, j, kd] - v[i, j - 1, kd]) / dy +
        met[i, j, kd, 2, 3] *
        0.25 *
        (
            v[i, j, kd + 1] + v[i, j - 1, kd + 1] - v[i, j, kd - 1] -
            v[i, j - 1, kd - 1]
        ) / dz
    @ivy phiu =
        (v[i, j, ku] - v[i, j - 1, ku]) / dy +
        met[i, j, ku, 2, 3] *
        0.25 *
        (
            v[i, j, ku + 1] + v[i, j - 1, ku + 1] - v[i, j, ku - 1] -
            v[i, j - 1, ku - 1]
        ) / dz

    return (phid, phiu)
end

function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DVDZ,
)::NTuple{2, <:AbstractFloat}
    (; lz, dz, ztildetfc, jac, topography_surface) = state.grid
    (; v) = state.variables.predictands

    @ivy if ztildetfc[i, j, ku] < topography_surface[i, j]
        phid = 0.0
        phiu = 0.0
    elseif ztildetfc[i, j, kd] < topography_surface[i, j]
        phid = 0.0
        phiu =
            (v[i, j, ku + 1] - v[i, j, ku]) / dz / (
                jac[i, j, ku] * jac[i, j, ku + 1] /
                (jac[i, j, ku] + jac[i, j, ku + 1]) +
                jac[i, j + 1, ku] * jac[i, j + 1, ku + 1] /
                (jac[i, j + 1, ku] + jac[i, j + 1, ku + 1])
            )
    else
        if ztildetfc[i, j, ku] < lz
            phid =
                (v[i, j, kd + 1] - v[i, j, kd]) / dz / (
                    jac[i, j, kd] * jac[i, j, kd + 1] /
                    (jac[i, j, kd] + jac[i, j, kd + 1]) +
                    jac[i, j + 1, kd] * jac[i, j + 1, kd + 1] /
                    (jac[i, j + 1, kd] + jac[i, j + 1, kd + 1])
                )
            phiu =
                (v[i, j, ku + 1] - v[i, j, ku]) / dz / (
                    jac[i, j, ku] * jac[i, j, ku + 1] /
                    (jac[i, j, ku] + jac[i, j, ku + 1]) +
                    jac[i, j + 1, ku] * jac[i, j + 1, ku + 1] /
                    (jac[i, j + 1, ku] + jac[i, j + 1, ku + 1])
                )
        elseif ztildetfc[i, j, kd] < lz
            phid =
                (v[i, j, kd + 1] - v[i, j, kd]) / dz / (
                    jac[i, j, kd] * jac[i, j, kd + 1] /
                    (jac[i, j, kd] + jac[i, j, kd + 1]) +
                    jac[i, j + 1, kd] * jac[i, j + 1, kd + 1] /
                    (jac[i, j + 1, kd] + jac[i, j + 1, kd + 1])
                )
            phiu = 0.0
        else
            phid = 0.0
            phiu = 0.0
        end
    end

    return (phid, phiu)
end
