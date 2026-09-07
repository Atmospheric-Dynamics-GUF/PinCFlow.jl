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
\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial x}\\right) = \\frac{u_{\\mathrm{b}, i + 1 / 2} - u_{\\mathrm{b}, i - 1 / 2}}{\\Delta \\hat{x}} + G^{13} \\frac{u_{\\mathrm{b}, k + 1} - u_{\\mathrm{b}, k - 1}}{2 \\Delta \\hat{z}}.
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
\\begin{align*}
    \\left(\\frac{\\partial u_\\mathrm{b}}{\\partial y}\\right)_{i + 1 / 2, j + 1 / 2} & = \\frac{u_{\\mathrm{b}, i + 1 / 2, j + 1} - u_{\\mathrm{b}, i + 1 / 2}}{\\Delta \\hat{y}} + G_{i + 1 / 2, j + 1 / 2}^{23} \\frac{u_{\\mathrm{b}, i + 1 / 2, j + 1 / 2, k + 1} - u_{\\mathrm{b}, i + 1 / 2, j + 1 / 2, k - 1}}{2 \\Delta \\hat{z}}.
\\end{align*}
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
\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial z}\\right)_{i + 1 / 2, k + 1 / 2} = \\frac{u_{\\mathrm{b}, i + 1 / 2, k + 1} - u_{\\mathrm{b}, i + 1 / 2}}{J_{i + 1 / 2, k + 1 / 2} \\Delta \\hat{z}}.
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
\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial x}\\right)_{i + 1 / 2, j + 1 / 2} = \\frac{v_{\\mathrm{b}, i + 1, j + 1 / 2} - v_{\\mathrm{b}, j + 1 / 2}}{\\Delta \\hat{x}} + G_{i + 1 / 2, j + 1 / 2}^{13} \\frac{v_{\\mathrm{b}, i + 1 / 2, j + 1 / 2, k + 1} - v_{\\mathrm{b}, i + 1 / 2, j + 1 / 2, k - 1}}{2 \\Delta \\hat{z}}.
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
\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial y}\\right) = \\frac{v_{\\mathrm{b}, j + 1 / 2} - v_{\\mathrm{b}, j - 1 / 2}}{\\Delta \\hat{y}} + G^{23} \\frac{v_{\\mathrm{b}, k + 1} - v_{\\mathrm{b}, k - 1} }{2 \\Delta \\hat{z}}.
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
\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial z}\\right)_{j + 1 / 2, k + 1 / 2} = \\frac{v_{\\mathrm{b}, j + 1 / 2, k + 1} - v_{\\mathrm{b}, j + 1 / 2}}{J_{j + 1 / 2, k + 1 / 2} \\Delta \\hat{z}}.
```

At grid points beyond the vertical boundaries, it is set to zero.

```julia
compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    field::AbstractArray{<:AbstractFloat, 3},
    phitype::DX,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the zonal derivative of a scalar field (``\\partial \\psi_\\mathrm{b} / \\partial x``) at ``\\left(i + 1 / 2, j, k_\\mathrm{D}\\right)`` and ``\\left(i + 1 / 2, j, k_\\mathrm{U}\\right)``.

The derivative is given by

```math
\\left(\\frac{\\partial \\psi_\\mathrm{b}}{\\partial x}\\right)_{i + 1 / 2} = \\frac{\\psi_{\\mathrm{b}, i + 1} - \\psi_\\mathrm{b}}{\\Delta \\hat{x}} + G_{i + 1 / 2}^{13} \\frac{\\psi_{\\mathrm{b}, i + 1 / 2, k + 1} - \\psi_{\\mathrm{b}, i + 1 / 2, k - 1}}{2 \\Delta \\hat{z}}.
```

```julia
compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    field::AbstractArray{<:AbstractFloat, 3},
    phitype::DDY,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the meridional derivative of a scalar field (``\\partial \\psi_\\mathrm{b} / \\partial y``) at ``\\left(i, j + 1 / 2, k_\\mathrm{D}\\right)`` and ``\\left(i, j + 1 / 2, k_\\mathrm{U}\\right)``.

The derivative is given by

```math
\\left(\\frac{\\partial \\psi_\\mathrm{b}}{\\partial y}\\right)_{j + 1 / 2} = \\frac{\\psi_{\\mathrm{b}, j + 1} - \\psi_\\mathrm{b}}{\\Delta \\hat{y}} + G_{j + 1 / 2}^{23} \\frac{\\psi_{\\mathrm{b}, j + 1 / 2, k + 1} - \\psi_{\\mathrm{b}, j + 1 / 2, k - 1}}{2 \\Delta \\hat{z}}.
```

```julia
compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    field::AbstractArray{<:AbstractFloat, 3}phitype::DChiDZ,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the vertical derivative of a scalar field (``\\partial \\psi_\\mathrm{b} / \\partial z``) at ``\\left(i, j, k_\\mathrm{D} + 1 / 2\\right)`` and ``\\left(i, j, k_\\mathrm{U} + 1 / 2\\right)``.

The derivative is given by

```math
\\left(\\frac{\\partial \\psi_\\mathrm{b}}{\\partial z}\\right)_{k + 1 / 2} = \\frac{\\psi_{\\mathrm{b}, k + 1} - \\psi_\\mathrm{b}}{J_{k + 1 / 2} \\Delta \\hat{z}}.
```

At grid points beyond the vertical boundaries, it is set to zero.

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `kd`: Lower vertical grid-cell index.

  - `ku`: Upper vertical grid-cell index.

  - `field`: Input array.

  - `phitype`: Type of derivative to compute.
"""
function compute_derivatives end

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DUDX,
)::NTuple{2, <:AbstractFloat}
    (; dx, dz, met) = state.grid
    (; u) = state.variables.predictands

    phid =
        (u[i, j, kd] - u[i - 1, j, kd]) / dx +
        met[i, j, kd, 1, 3] *
        0.25 *
        (
            u[i, j, kd + 1] + u[i - 1, j, kd + 1] - u[i, j, kd - 1] -
            u[i - 1, j, kd - 1]
        ) / dz
    phiu =
        (u[i, j, ku] - u[i - 1, j, ku]) / dx +
        met[i, j, ku, 1, 3] *
        0.25 *
        (
            u[i, j, ku + 1] + u[i - 1, j, ku + 1] - u[i, j, ku - 1] -
            u[i - 1, j, ku - 1]
        ) / dz

    return (phid, phiu)
end

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DUDY,
)::NTuple{2, <:AbstractFloat}
    (; dy, dz, met) = state.grid
    (; u) = state.variables.predictands

    phid =
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
    phiu =
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

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DUDZ,
)::NTuple{2, <:AbstractFloat}
    (; lz, dz, zctilde, jac, hb) = state.grid
    (; u) = state.variables.predictands

    if (zctilde[i, j, ku] + zctilde[i + 1, j, ku]) / 2 <
       (hb[i, j] + hb[i + 1, j]) / 2
        phid = 0.0
        phiu = 0.0
    elseif (zctilde[i, j, kd] + zctilde[i + 1, j, kd]) / 2 <
           (hb[i, j] + hb[i + 1, j]) / 2
        phid = 0.0
        phiu =
            (u[i, j, ku + 1] - u[i, j, ku]) / dz / (
                jac[i, j, ku] * jac[i, j, ku + 1] /
                (jac[i, j, ku] + jac[i, j, ku + 1]) +
                jac[i + 1, j, ku] * jac[i + 1, j, ku + 1] /
                (jac[i + 1, j, ku] + jac[i + 1, j, ku + 1])
            )
    else
        if (zctilde[i, j, ku] + zctilde[i + 1, j, ku]) / 2 < lz
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
        elseif (zctilde[i, j, kd] + zctilde[i + 1, j, kd]) / 2 < lz
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

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DVDX,
)::NTuple{2, <:AbstractFloat}
    (; dx, dz, met) = state.grid
    (; v) = state.variables.predictands

    phid =
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
    phiu =
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

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DVDY,
)::NTuple{2, <:AbstractFloat}
    (; dy, dz, met) = state.grid
    (; v) = state.variables.predictands

    phid =
        (v[i, j, kd] - v[i, j - 1, kd]) / dy +
        met[i, j, kd, 2, 3] *
        0.25 *
        (
            v[i, j, kd + 1] + v[i, j - 1, kd + 1] - v[i, j, kd - 1] -
            v[i, j - 1, kd - 1]
        ) / dz
    phiu =
        (v[i, j, ku] - v[i, j - 1, ku]) / dy +
        met[i, j, ku, 2, 3] *
        0.25 *
        (
            v[i, j, ku + 1] + v[i, j - 1, ku + 1] - v[i, j, ku - 1] -
            v[i, j - 1, ku - 1]
        ) / dz

    return (phid, phiu)
end

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    phitype::DVDZ,
)::NTuple{2, <:AbstractFloat}
    (; lz, dz, zctilde, jac, hb) = state.grid
    (; v) = state.variables.predictands

    if (zctilde[i, j, ku] + zctilde[i, j + 1, ku]) / 2 <
       (hb[i, j] + hb[i, j + 1]) / 2
        phid = 0.0
        phiu = 0.0
    elseif (zctilde[i, j, kd] + zctilde[i, j + 1, kd]) / 2 <
           (hb[i, j] + hb[i, j + 1]) / 2
        phid = 0.0
        phiu =
            (v[i, j, ku + 1] - v[i, j, ku]) / dz / (
                jac[i, j, ku] * jac[i, j, ku + 1] /
                (jac[i, j, ku] + jac[i, j, ku + 1]) +
                jac[i, j + 1, ku] * jac[i, j + 1, ku + 1] /
                (jac[i, j + 1, ku] + jac[i, j + 1, ku + 1])
            )
    else
        if (zctilde[i, j, ku] + zctilde[i, j + 1, ku]) / 2 < lz
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
        elseif (zctilde[i, j, kd] + zctilde[i, j + 1, kd]) / 2 < lz
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

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    field::AbstractArray{<:AbstractFloat, 3},
    phitype::DX,
)::NTuple{2, <:AbstractFloat}
    (; dx, dz, met) = state.grid

    cc = field[i, j, kd]
    cr = field[i + 1, j, kd]
    cu = field[i, j, kd + 1]
    cd = field[i, j, kd - 1]
    cru = field[i + 1, j, kd + 1]
    crd = field[i + 1, j, kd - 1]

    phid =
        (cr - cc) / dx +
        0.5 *
        (met[i, j, kd, 1, 3] + met[i + 1, j, kd, 1, 3]) *
        0.25 *
        (cu + cru - cd - crd) / dz

    cc = field[i, j, ku]
    cr = field[i + 1, j, ku]
    cu = field[i, j, ku + 1]
    cd = field[i, j, ku - 1]
    cru = field[i + 1, j, ku + 1]
    crd = field[i + 1, j, ku - 1]

    phiu =
        (cr - cc) / dx +
        0.5 *
        (met[i, j, ku, 1, 3] + met[i + 1, j, ku, 1, 3]) *
        0.25 *
        (cu + cru - cd - crd) / dz

    return (phid, phiu)
end

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    field::AbstractArray{<:AbstractFloat, 3},
    phitype::DY,
)::NTuple{2, <:AbstractFloat}
    (; dy, dz, met) = state.grid

    cc = field[i, j, kd]
    cf = field[i, j + 1, kd]
    cu = field[i, j, kd + 1]
    cd = field[i, j, kd - 1]
    cfu = field[i, j + 1, kd + 1]
    cfd = field[i, j + 1, kd - 1]

    phid =
        (cf - cc) / dy +
        0.5 *
        (met[i, j, kd, 2, 3] + met[i, j + 1, kd, 2, 3]) *
        0.25 *
        (cu + cfu - cd - cfd) / dz

    cc = field[i, j, ku]
    cf = field[i, j + 1, ku]
    cu = field[i, j, ku + 1]
    cd = field[i, j, ku - 1]
    cfu = field[i, j + 1, ku + 1]
    cfd = field[i, j + 1, ku - 1]

    phiu =
        (cf - cc) / dy +
        0.5 *
        (met[i, j, ku, 2, 3] + met[i, j + 1, ku, 2, 3]) *
        0.25 *
        (cu + cfu - cd - cfd) / dz

    return (phid, phiu)
end

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    field::AbstractArray{<:AbstractFloat, 3},
    phitype::DZ,
)::NTuple{2, <:AbstractFloat}
    (; lz, dz, zctilde, jac, hb) = state.grid

    cuc = field[i, j, ku]
    cdc = field[i, j, kd]
    cuu = field[i, j, ku + 1]
    cdu = field[i, j, kd + 1]

    if zctilde[i, j, ku] < hb[i, j]
        phid = 0.0
        phiu = 0.0
    elseif zctilde[i, j, kd] < hb[i, j]
        phid = 0.0
        phiu =
            (cuu - cuc) / dz / (
                2.0 * jac[i, j, ku] * jac[i, j, ku + 1] /
                (jac[i, j, ku] + jac[i, j, ku + 1])
            )
    else
        if zctilde[i, j, ku] < lz
            phid =
                (cdu - cdc) / dz / (
                    2.0 * jac[i, j, kd] * jac[i, j, kd + 1] /
                    (jac[i, j, kd] + jac[i, j, kd + 1])
                )
            phiu =
                (cuu - cuc) / dz / (
                    2.0 * jac[i, j, ku] * jac[i, j, ku + 1] /
                    (jac[i, j, ku] + jac[i, j, ku + 1])
                )
        elseif zctilde[i, j, kd] < lz
            phid =
                (cdu - cdc) / dz / (
                    2.0 * jac[i, j, kd] * jac[i, j, kd + 1] /
                    (jac[i, j, kd] + jac[i, j, kd + 1])
                )
            phiu = 0.0
        else
            phid = 0.0
            phiu = 0.0
        end
    end

    return (phid, phiu)
end
