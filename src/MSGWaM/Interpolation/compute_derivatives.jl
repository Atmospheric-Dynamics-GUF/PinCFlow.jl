"""
```julia
compute_derivatives(
    state::State,
    indices::NTuple{4, <:Integer},
    phitype::DUDX,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the zonal derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial x``) near the two grid points specified by `indices`.

The derivative is given by

```math
\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial x}\\right) = \\frac{u_{\\mathrm{b}, i + 1 / 2} - u_{\\mathrm{b}, i - 1 / 2}}{\\Delta \\widehat{x}} + G^{13} \\frac{u_{\\mathrm{b}, i + 1 / 2, k + 1} + u_{\\mathrm{b}, i - 1 / 2, k + 1} - u_{\\mathrm{b}, i + 1 / 2, k - 1} - u_{\\mathrm{b}, i - 1 / 2, k - 1}}{4 \\Delta \\widehat{z}}.
```

```julia
compute_derivatives(
    state::State,
    indices::NTuple{4, <:Integer},
    phitype::DUDY,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the meridional derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial y``) near the two grid points specified by `indices`.

The derivative is given by

```math
\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial y}\\right)_{i + 1 / 2, j + 1 / 2} = \\frac{u_{\\mathrm{b}, i + 1 / 2, j + 1} - u_{\\mathrm{b}, i + 1 / 2}}{\\Delta \\widehat{y}} + \\frac{G^{23} + G^{23}_{i + 1} + G^{23}_{j + 1} + G^{23}_{i + 1, j + 1}}{4} \\frac{u_{\\mathrm{b}, i + 1 / 2, k + 1} + u_{\\mathrm{b}, i + 1 / 2, j + 1, k + 1} - u_{\\mathrm{b}, i + 1 / 2, k - 1} - u_{\\mathrm{b}, i + 1 / 2, j + 1, k - 1}}{4 \\Delta \\widehat{z}}.
```

```julia
compute_derivatives(
    state::State,
    indices::NTuple{4, <:Integer},
    phitype::DUDZ,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the vertical derivative of the zonal wind (``\\partial u_\\mathrm{b} / \\partial z``) near the two grid points specified by `indices`.

The derivative is given by

```math
\\left(\\frac{\\partial u_\\mathrm{b}}{\\partial z}\\right)_{i + 1 / 2, k + 1 / 2} = \\frac{u_{\\mathrm{b}, i + 1 / 2, k + 1} - u_{\\mathrm{b}, i + 1 / 2}}{\\Delta \\widehat{z}} \\left(\\frac{J J_{k + 1}}{J + J_{k + 1}} + \\frac{J_{i + 1} J_{i + 1, k + 1}}{J_{i + 1} + J_{i + 1, k + 1}}\\right)^{- 1}.
```

At grid points beyond the vertical boundaries, it is set to zero.

```julia
compute_derivatives(
    state::State,
    indices::NTuple{4, <:Integer},
    phitype::DVDX,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the zonal derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial x``) near the two grid points specified by `indices`.

The derivative is given by

```math
\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial x}\\right)_{i + 1 / 2, j + 1 / 2} = \\frac{v_{\\mathrm{b}, i + 1, j + 1 / 2} - v_{\\mathrm{b}, j + 1 / 2}}{\\Delta \\widehat{x}} + \\frac{G^{13} + G^{13}_{i + 1} + G^{13}_{j + 1} + G^{13}_{i + 1, j + 1}}{4} \\frac{v_{\\mathrm{b}, j + 1 / 2, k + 1} + v_{\\mathrm{b}, i + 1, j + 1 / 2, k + 1} - v_{\\mathrm{b}, j + 1 / 2, k - 1} - v_{\\mathrm{b}, i + 1, j + 1 / 2, k - 1}}{4 \\Delta \\widehat{z}}.
```

```julia
compute_derivatives(
    state::State,
    indices::NTuple{4, <:Integer},
    phitype::DVDY,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the meridional derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial y``) near the two grid points specified by `indices`.

The derivative is given by

```math
\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial y}\\right) = \\frac{v_{\\mathrm{b}, j + 1 / 2} - v_{\\mathrm{b}, j - 1 / 2}}{\\Delta \\widehat{y}} + G^{23} \\frac{v_{\\mathrm{b}, j + 1 / 2, k + 1} + v_{\\mathrm{b}, j - 1 / 2, k + 1} - v_{\\mathrm{b}, j + 1 / 2, k - 1} - v_{\\mathrm{b}, j - 1 / 2, k - 1}}{4 \\Delta \\widehat{z}}.
```

```julia
compute_derivatives(
    state::State,
    indices::NTuple{4, <:Integer},
    phitype::DVDZ,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the vertical derivative of the meridional wind (``\\partial v_\\mathrm{b} / \\partial z``) near the two grid points specified by `indices`.

The derivative is given by

```math
\\left(\\frac{\\partial v_\\mathrm{b}}{\\partial z}\\right)_{j + 1 / 2, k + 1 / 2} = \\frac{v_{\\mathrm{b}, j + 1 / 2, k + 1} - v_{\\mathrm{b}, j + 1 / 2}}{\\Delta \\widehat{z}} \\left(\\frac{J J_{k + 1}}{J + J_{k + 1}} + \\frac{J_{j + 1} J_{j + 1, k + 1}}{J_{j + 1} + J_{j + 1, k + 1}}\\right)^{- 1}.
```

At grid points beyond the vertical boundaries, it is set to zero.

```julia 
compute_derivatives(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDX,
)
```

Compute and return the zonal derivative of the tracer field (``\\partial \\chi_\\mathrm{b} / \\partial x``) near the location specified by `xlc`, `ylc`, and `zlc`.

```julia 
compute_derivatives(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDY,
)
```

Compute and return the meridional derivative of the tracer field (``\\partial \\chi_\\mathrm{b} / \\partial y``) near the location specified by `xlc`, `ylc`, and `zlc`.

```julia 
compute_derivatives(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDZ,
)
```

Compute and return the vertical derivative of the tracer field (``\\partial \\chi_\\mathrm{b} / \\partial z``) near the location specified by `xlc`, `ylc`, and `zlc`.

# Arguments

  - `state`: Model state.

  - `indices`: Indices `(ix, jy, kzd, kzu)` of the two grid points at which to compute the derivative.

  - `phitype`: Type of derivative to compute.

  - `xlc`: Location in `\\widehat{x}`-direction.

  - `ylc`: Location in `\\widehat{y}`-direction.

  - `zlc`: Location in `\\widehat{z}`-direction.
"""
function compute_derivatives end

function compute_derivatives(
    state::State,
    indices::NTuple{4, <:Integer},
    phitype::DUDX,
)::NTuple{2, <:AbstractFloat}
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
)::NTuple{2, <:AbstractFloat}
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
)::NTuple{2, <:AbstractFloat}
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
)::NTuple{2, <:AbstractFloat}
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
)::NTuple{2, <:AbstractFloat}
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
)::NTuple{2, <:AbstractFloat}
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

function compute_derivatives(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDX,
)
    (; sizex, sizey) = state.namelists.domain
    (; lx, ly, lz, dx, dy, dz, x, y, ztfc) = state.grid
    (; nxx, nyy, nzz, io, jo, ko, i0, j0, k0) = state.domain
    (; chi) = state.tracer.tracerpredictands
    (; rho) = state.variables.predictands
    (; rhostrattfc) = state.atmosphere
    (; domain, grid) = state

    if sizex == 1
        return 0.0
    else
        # find closest points in zonal direction 
        ixl = floor(Int, (xlc - lx[1]) / dx) + i0 - 1 - io
        if (ixl < 1)
            error("Error in compute_derivatives (DChiDX): ixl = ", ixl, " < 1")
        end
        ixr = ixl + 1
        if ixr > nxx
            error(
                "Error in compute_derivatives (DChiDX): ixr = ",
                ixr,
                "> nxx = ",
                nxx,
            )
        end

        xr = x[ixr + io]
        xl = x[ixl + io]

        if abs(xr - xlc) > abs(xlc - xl)
            ix = ixl
        else
            ix = ixr
        end

        # find closest points in meridional direction
        if sizey == 1
            jyb = j0
            jyf = j0
        else
            jyb = floor(Int, (ylc - ly[1] - dy / 2) / dy) + j0 - jo
            if (jyb < 1)
                error(
                    "Error in compute_derivatives (DChiDX): jyl = ",
                    jyl,
                    " < 1",
                )
            end
            jyf = jyb + 1
            if jyf > nyy
                error(
                    "Error in compute_derivatives (DChiDX): jyr = ",
                    jyr,
                    " > nyy = ",
                    nyy,
                )
            end
        end
        yf = y[jyf + jo]
        yb = y[jyb + jo]

        if abs(yf - ylc) > abs(ylc - yb)
            jy = jyb
        else
            jy = jyf
        end

        # find closest index on z-axis 
        kzu = get_next_level(ix, jy, zlc, domain, grid)
        kzd = kzu - 1

        zu = ztfc[ix, jy, kzu]
        zd = ztfc[ix, jy, kzd]

        if abs(zu - zlc) > abs(zlc - zd)
            kz = kzd
        else
            kz = kzu
        end

        dchil =
            (
                chi[ixl + 1, jy, kz] /
                (rho[ixl + 1, jy, kz] + rhostrattfc[ixl + 1, jy, kz]) -
                chi[ixl, jy, kz] /
                (rho[ixl, jy, kz] + rhostrattfc[ixl, jy, kz])
            ) / dx
        dchir =
            (
                chi[ixr + 1, jy, kz] /
                (rho[ixr + 1, jy, kz] + rhostrattfc[ixr + 1, jy, kz]) -
                chi[ixr, jy, kz] /
                (rho[ixr, jy, kz] + rhostrattfc[ixr, jy, kz])
            ) / dx

        if xr < xl
            then
            error(
                "Error in compute_derivatives (DChiDX): xr = ",
                xr,
                " < xl = ",
                xl,
            )
        elseif xr == xl
            factor = 0.0
        elseif xlc > xr
            factor = 0.0
        elseif xlc > xl
            factor = (xr - xlc) / dx
        else
            factor = 1.0
        end

        dchidx = factor * dchil + (1.0 - factor) * dchir

        return dchidx
    end
end

function compute_derivatives(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDY,
)
    (; sizex, sizey) = state.namelists.domain
    (; lx, ly, lz, dx, dy, dz, x, y, ztfc) = state.grid
    (; nxx, nyy, nzz, io, jo, ko, i0, j0, k0) = state.domain
    (; chi) = state.tracer.tracerpredictands
    (; rho) = state.variables.predictands
    (; rhostrattfc) = state.atmosphere
    (; domain, grid) = state

    if sizey == 1
        return 0.0
    else
        # find closest points in zonal direction 
        if sizex == 1
            ixl = i0
            ixr = i0
        else
            ixl = floor(Int, (xlc - lx[1]) / dx) + i0 - 1 - io
            if (ixl < 1)
                error(
                    "Error in compute_derivatives (DChiDY): ixl = ",
                    ixl,
                    " < 1",
                )
            end
            ixr = ixl + 1
            if ixr > nxx
                error(
                    "Error in compute_derivatives (DChiDY): ixr = ",
                    ixr,
                    "> nxx = ",
                    nxx,
                )
            end
        end
        xr = x[ixr + io]
        xl = x[ixl + io]

        if abs(xr - xlc) > abs(xlc - xl)
            ix = ixl
        else
            ix = ixr
        end

        # find closest points in meridional direction
        jyb = floor(Int, (ylc - ly[1] - dy / 2) / dy) + j0 - jo
        if (jyb < 1)
            error("Error in compute_derivatives (DChiDY): jyl = ", jyl, " < 1")
        end
        jyf = jyb + 1
        if jyf > nyy
            error(
                "Error in compute_derivatives (DChiDY): jyr = ",
                jyr,
                " > nyy = ",
                nyy,
            )
        end

        yf = y[jyf + jo]
        yb = y[jyb + jo]

        if abs(yf - ylc) > abs(ylc - yb)
            jy = jyb
        else
            jy = jyf
        end

        # find closest index on z-axis 
        kzu = get_next_level(ix, jy, zlc, domain, grid)
        kzd = kzu - 1

        zu = ztfc[ix, jy, kzu]
        zd = ztfc[ix, jy, kzd]

        if abs(zu - zlc) > abs(zlc - zd)
            kz = kzd
        else
            kz = kzu
        end

        dchif =
            (
                chi[ix, jyf + 1, kz] /
                (rho[ix, jyf + 1, kz] + rhostrattfc[ix, jyf + 1, kz]) -
                chi[ix, jyf, kz] /
                (rho[ix, jyf, kz] + rhostrattfc[ix, jyf, kz])
            ) / dy
        dchib =
            (
                chi[ix, jyb + 1, kz] /
                (rho[ix, jyb + 1, kz] + rhostrattfc[ix, jyb + 1, kz]) -
                chi[ix, jyb, kz] /
                (rho[ix, jyb, kz] + rhostrattfc[ix, jyb, kz])
            ) / dy

        if yf < yb
            then
            error(
                "Error in compute_derivatives (DChiDY): yf = ",
                yf,
                " < yb = ",
                yb,
            )
        elseif yf == yb
            factor = 0.0
        elseif ylc > yf
            factor = 0.0
        elseif ylc > yb
            factor = (yf - ylc) / dy
        else
            factor = 1.0
        end

        dchidy = factor * dchib + (1.0 - factor) * dchif

        return dchidy
    end
end

function compute_derivatives(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
    phitype::DChiDZ,
)
    (; sizex, sizey) = state.namelists.domain
    (; lx, ly, lz, dx, dy, dz, x, y, ztfc) = state.grid
    (; nxx, nyy, nzz, io, jo, ko, i0, j0, k0) = state.domain
    (; chi) = state.tracer.tracerpredictands
    (; rho) = state.variables.predictands
    (; rhostrattfc) = state.atmosphere
    (; domain, grid) = state

    # find closest points in zonal direction 
    if sizex == 1
        ixl = i0
        ixr = i0
    else
        ixl = floor(Int, (xlc - lx[1]) / dx) + i0 - 1 - io
        if (ixl < 1)
            error("Error in compute_derivatives (DChiDZ): ixl = ", ixl, " < 1")
        end
        ixr = ixl + 1
        if ixr > nxx
            error(
                "Error in compute_derivatives (DChiDZ): ixr = ",
                ixr,
                "> nxx = ",
                nxx,
            )
        end
    end
    xr = x[ixr + io]
    xl = x[ixl + io]

    if abs(xr - xlc) > abs(xlc - xl)
        ix = ixl
    else
        ix = ixr
    end

    # find closest points in meridional direction
    if sizey == 1
        jyb = j0
        jyf = j0
    else
        jyb = floor(Int, (ylc - ly[1] - dy / 2) / dy) + j0 - jo
        if (jyb < 1)
            error("Error in compute_derivatives (DChiDZ): jyl = ", jyl, " < 1")
        end
        jyf = jyb + 1
        if jyf > nyy
            error(
                "Error in compute_derivatives (DChiDZ): jyr = ",
                jyr,
                " > nyy = ",
                nyy,
            )
        end
    end
    yf = y[jyf + jo]
    yb = y[jyb + jo]

    if abs(yf - ylc) > abs(ylc - yb)
        jy = jyb
    else
        jy = jyf
    end

    # find closest index on z-axis 
    kzu = get_next_level(ix, jy, zlc, domain, grid)
    kzd = kzu - 1

    zu = ztfc[ix, jy, kzu]
    zd = ztfc[ix, jy, kzd]

    if abs(zu - zlc) > abs(zlc - zd)
        kz = kzd
    else
        kz = kzu
    end

    dchiu =
        (
            chi[ix, jy, kzu + 1] /
            (rho[ix, jy, kzu + 1] + rhostrattfc[ix, jy, kzu + 1]) -
            chi[ix, jy, kzu] / (rho[ix, jy, kzu] + rhostrattfc[ix, jy, kzu])
        ) / dz
    dchid =
        (
            chi[ix, jy, kzd + 1] /
            (rho[ix, jy, kzd + 1] + rhostrattfc[ix, jy, kzd + 1]) -
            chi[ix, jy, kzd] / (rho[ix, jy, kzd] + rhostrattfc[ix, jy, kzd])
        ) / dz

    if zu < zd
        then
        error(
            "Error in compute_derivatives (DChiDZ): zu = ",
            zu,
            " < zd = ",
            zd,
        )
    elseif zu == zd
        factor = 0.0
    elseif zlc > zu
        factor = 0.0
    elseif zlc > zd
        factor = (zu - zlc) / dz
    else
        factor = 1.0
    end

    dchidz = factor * dchid + (1.0 - factor) * dchiu

    return dchidz
end
