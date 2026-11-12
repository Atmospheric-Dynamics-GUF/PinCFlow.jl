"""
```julia
compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    derivative::DUDY,
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
    derivative::DUDZ,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the vertical derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial z``) at ``\\left(i + 1 / 2, j, k_\\mathrm{D} + 1 / 2\\right)`` and ``\\left(i + 1 / 2, j, k_\\mathrm{U} + 1 / 2\\right)``.

The derivative is given by

```math
\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial z}\\right)_{i + 1 / 2, k + 1 / 2} = \\frac{u_{\\mathrm{b}, i + 1 / 2, k + 1} - u_{\\mathrm{b}, i + 1 / 2}}{J_{i + 1 / 2, k + 1 / 2} \\Delta \\hat{z}}.
```

In the case of solid-wall vertical boundary conditions, the derivative is set to zero at grid points beyond the vertical boundaries.

```julia
compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    derivative::DVDX,
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
    derivative::DVDZ,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the vertical derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial z``) at ``\\left(i, j + 1 / 2, k_\\mathrm{D} + 1 / 2\\right)`` and ``\\left(i, j + 1 / 2, k_\\mathrm{U} + 1 / 2\\right)``.

The derivative is given by

```math
\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial z}\\right)_{j + 1 / 2, k + 1 / 2} = \\frac{v_{\\mathrm{b}, j + 1 / 2, k + 1} - v_{\\mathrm{b}, j + 1 / 2}}{J_{j + 1 / 2, k + 1 / 2} \\Delta \\hat{z}}.
```

In the case of solid-wall vertical boundary conditions, the derivative is set to zero at grid points beyond the vertical boundaries.

```julia
compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    field::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    derivative::DX,
)::Union{
    NTuple{2, T},
    NTuple{2, Complex{T}},
    Tuple{T, Complex{T}},
    Tuple{Complex{T}, T},
} where {T <: Real}
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
    field::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    derivative::DY,
)::Union{
    NTuple{2, T},
    NTuple{2, Complex{T}},
    Tuple{T, Complex{T}},
    Tuple{Complex{T}, T},
} where {T <: Real}
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
    field::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    derivative::DZ,
)::Union{
    NTuple{2, T},
    NTuple{2, Complex{T}},
    Tuple{T, Complex{T}},
    Tuple{Complex{T}, T},
} where {T <: Real}
```

Compute and return the vertical derivative of a scalar field (``\\partial \\psi_\\mathrm{b} / \\partial z``) at ``\\left(i, j, k_\\mathrm{D} + 1 / 2\\right)`` and ``\\left(i, j, k_\\mathrm{U} + 1 / 2\\right)``.

The derivative is given by

```math
\\left(\\frac{\\partial \\psi_\\mathrm{b}}{\\partial z}\\right)_{k + 1 / 2} = \\frac{\\psi_{\\mathrm{b}, k + 1} - \\psi_\\mathrm{b}}{J_{k + 1 / 2} \\Delta \\hat{z}}.
```

In the case of solid-wall vertical boundary conditions, the derivative is set to zero at grid points beyond the vertical boundaries.

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `kd`: Lower vertical grid-cell index.

  - `ku`: Upper vertical grid-cell index.

  - `field`: Input array.

  - `derivative`: Type of derivative to compute.
"""
function compute_derivatives end

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    derivative::DUDY,
)::NTuple{2, <:AbstractFloat}
    (; dy, dz, met23) = state.grid
    (; u) = state.variables.predictands

    phid =
        (u[i, j + 1, kd] - u[i, j, kd]) / dy +
        0.25 *
        (
            met23[i, j, kd] +
            met23[i + 1, j, kd] +
            met23[i, j + 1, kd] +
            met23[i + 1, j + 1, kd]
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
            met23[i, j, ku] +
            met23[i + 1, j, ku] +
            met23[i, j + 1, ku] +
            met23[i + 1, j + 1, ku]
        ) *
        0.25 *
        (
            u[i, j, ku + 1] + u[i, j + 1, ku + 1] - u[i, j, ku - 1] -
            u[i, j + 1, ku - 1]
        ) / dz

    (phid, phiu)
end

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    derivative::DUDZ,
)::NTuple{2, <:AbstractFloat}
    (; lz, dz, zctilde, jac, hb) = state.grid
    (; u) = state.variables.predictands
    (; vertical_boundary_condition) = state.namelists.domain

    if vertical_boundary_condition === :Periodic
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
    else
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
    end

    (phid, phiu)
end

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    derivative::DVDX,
)::NTuple{2, <:AbstractFloat}
    (; dx, dz, met13) = state.grid
    (; v) = state.variables.predictands

    phid =
        (v[i + 1, j, kd] - v[i, j, kd]) / dx +
        0.25 *
        (
            met13[i, j, kd] +
            met13[i + 1, j, kd] +
            met13[i, j + 1, kd] +
            met13[i + 1, j + 1, kd]
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
            met13[i, j, ku] +
            met13[i + 1, j, ku] +
            met13[i, j + 1, ku] +
            met13[i + 1, j + 1, ku]
        ) *
        0.25 *
        (
            v[i, j, ku + 1] + v[i + 1, j, ku + 1] - v[i, j, ku - 1] -
            v[i + 1, j, ku - 1]
        ) / dz

    (phid, phiu)
end

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    derivative::DVDZ,
)::NTuple{2, <:AbstractFloat}
    (; lz, dz, zctilde, jac, hb) = state.grid
    (; v) = state.variables.predictands
    (; vertical_boundary_condition) = state.namelists.domain

    if vertical_boundary_condition === :Periodic
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
    else
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
    end

    (phid, phiu)
end

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    field::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    derivative::DX,
)::Union{
    NTuple{2, T},
    NTuple{2, Complex{T}},
    Tuple{T, Complex{T}},
    Tuple{Complex{T}, T},
} where {T <: Real}
    (; dx, dz, met13) = state.grid

    cc = field[i, j, kd]
    cr = field[i + 1, j, kd]
    cu = field[i, j, kd + 1]
    cd = field[i, j, kd - 1]
    cru = field[i + 1, j, kd + 1]
    crd = field[i + 1, j, kd - 1]

    phid =
        (cr - cc) / dx +
        0.5 *
        (met13[i, j, kd] + met13[i + 1, j, kd]) *
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
        (met13[i, j, ku] + met13[i + 1, j, ku]) *
        0.25 *
        (cu + cru - cd - crd) / dz

    (phid, phiu)
end

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    field::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    derivative::DY,
)::Union{
    NTuple{2, T},
    NTuple{2, Complex{T}},
    Tuple{T, Complex{T}},
    Tuple{Complex{T}, T},
} where {T <: Real}
    (; dy, dz, met23) = state.grid

    cc = field[i, j, kd]
    cf = field[i, j + 1, kd]
    cu = field[i, j, kd + 1]
    cd = field[i, j, kd - 1]
    cfu = field[i, j + 1, kd + 1]
    cfd = field[i, j + 1, kd - 1]

    phid =
        (cf - cc) / dy +
        0.5 *
        (met23[i, j, kd] + met23[i, j + 1, kd]) *
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
        (met23[i, j, ku] + met23[i, j + 1, ku]) *
        0.25 *
        (cu + cfu - cd - cfd) / dz

    (phid, phiu)
end

@ivy function compute_derivatives(
    state::State,
    i::Integer,
    j::Integer,
    kd::Integer,
    ku::Integer,
    field::Union{AbstractArray{T, 3}, AbstractArray{Complex{T}, 3}},
    derivative::DZ,
)::Union{
    NTuple{2, T},
    NTuple{2, Complex{T}},
    Tuple{T, Complex{T}},
    Tuple{Complex{T}, T},
} where {T <: Real}
    (; lz, dz, zctilde, jac, hb) = state.grid
    (; vertical_boundary_condition) = state.namelists.domain

    cuc = field[i, j, ku]
    cdc = field[i, j, kd]
    cuu = field[i, j, ku + 1]
    cdu = field[i, j, kd + 1]
    if vertical_boundary_condition === :Periodic
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
    else
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
    end

    (phid, phiu)
end
