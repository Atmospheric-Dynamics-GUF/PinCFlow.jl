"""
```julia
compute_saturation_integrals(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
)::NTuple{2, <:AbstractFloat}
```

Compute and return the two spectral integrals ``S_1`` and ``S_2``, as needed by the saturation scheme (in the grid cell ``\\left(i, j, k\\right)``).

Computes the sums

```math
\\begin{align*}
    S_1 & \\approx \\sum\\limits_r \\left(m_r \\left|b_{\\mathrm{w}, r}\\right|\\right)^2 f_r,\\\\
    S_2 & \\approx \\sum\\limits_r \\left(m_r \\left|b_{\\mathrm{w}, r}\\right| \\left|\\boldsymbol{k}_r\\right|\\right)^2 f_r,
\\end{align*}
```

where

```math
f_r = \\max \\left(1, \\frac{\\Delta x_r}{\\Delta \\widehat{x}}\\right) \\max \\left(1, \\frac{\\Delta y_r}{\\Delta \\widehat{y}}\\right) \\max \\left(1, \\frac{\\Delta z_r}{J \\Delta \\widehat{z}}\\right)
```

is the maximum grid-cell fraction that can be covered by each ray volume (with ``\\left(\\Delta x_r, \\Delta y_r, \\Delta z_r\\right)`` being the ray-volume extents in physical space) and

```math
\\left|b_{\\mathrm{w}, r}\\right|^2 = \\frac{2}{\\overline{\\rho}} \\frac{N_r^4 \\left(k_r^2 + l_r^2\\right)}{\\widehat{\\omega}_r \\left|\\boldsymbol{k}_r\\right|^2} \\mathcal{N}_r \\Delta k_r \\Delta l_r \\Delta m_r
```

is the squared gravity-wave amplitude of the buoyancy. Therein, ``N_r^2`` is the squared buoyancy frequency interpolated to the ray-volume position (using `interpolate_stratification`) and ``\\left(\\Delta k_r, \\Delta l_r, \\Delta m_r\\right)`` are the ray-volume extents in spectral space.

# Arguments

  - `state`: Model state.

  - `i`: Zonal grid-cell index.

  - `j`: Meridional grid-cell index.

  - `k`: Vertical grid-cell index.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_half_level`](@ref)

  - [`PinCFlow.MSGWaM.RayOperations.compute_intrinsic_frequency`](@ref)

  - [`PinCFlow.MSGWaM.Interpolation.interpolate_stratification`](@ref)
"""
function compute_saturation_integrals end

function compute_saturation_integrals(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
)::NTuple{2, <:AbstractFloat}
    (; domain, grid) = state
    (; x_size, y_size) = state.namelists.domain
    (; io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, dz, jac) = grid
    (; rhobar) = state.atmosphere
    (; nray, rays) = state.wkb

    # Initialize Integrals.
    mb2 = 0.0
    mb2k2 = 0.0

    # Loop over ray volumes.
    @ivy for r in 1:nray[i, j, k]

        # Skip ray volumes with zero wave-action density.
        if rays.dens[r, i, j, k] == 0
            continue
        end

        xr = rays.x[r, i, j, k]
        yr = rays.y[r, i, j, k]
        zr = rays.z[r, i, j, k]

        dxr = rays.dxray[r, i, j, k]
        dyr = rays.dyray[r, i, j, k]
        dzr = rays.dzray[r, i, j, k]

        if x_size > 1
            iray = floor(Int, (xr + lx / 2) / dx) + i0 - io
        else
            iray = i0
        end

        if y_size > 1
            jray = floor(Int, (yr + ly / 2) / dy) + j0 - jo
        else
            jray = j0
        end

        kray = get_next_half_level(iray, jray, zr, state)

        n2r = interpolate_stratification(zr, state, N2())

        wnrk = rays.k[r, i, j, k]
        wnrl = rays.l[r, i, j, k]
        wnrm = rays.m[r, i, j, k]

        wnrhs = wnrk^2 + wnrl^2

        dwnrk = rays.dkray[r, i, j, k]
        dwnrl = rays.dlray[r, i, j, k]
        dwnrm = rays.dmray[r, i, j, k]

        omir = compute_intrinsic_frequency(state, r, i, j, k)

        densr = rays.dens[r, i, j, k]

        dzi = min(dzr, jac[iray, jray, kray] * dz)
        facpsp = dzi / jac[iray, jray, kray] / dz * dwnrm

        if x_size > 1
            dxi = min(dxr, dx)
            facpsp = facpsp * dxi / dx * dwnrk
        end

        if y_size > 1
            dyi = min(dyr, dy)
            facpsp = facpsp * dyi / dy * dwnrl
        end

        integral1 = wnrhs * wnrm^2 / ((wnrhs + wnrm^2) * omir) * facpsp

        mb2 += 2 * n2r^2 / rhobar[iray, jray, kray] * densr * integral1

        integral2 = wnrhs * wnrm^2 / omir * facpsp

        mb2k2 += 2 * n2r^2 / rhobar[iray, jray, kray] * densr * integral2
    end

    return (mb2, mb2k2)
end
