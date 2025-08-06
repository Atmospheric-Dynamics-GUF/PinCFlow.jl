"""
```julia
compute_saturation_integrals(state::State, indices::NTuple{3, <:Integer})
```

Compute two spectral integrals needed by the saturation scheme (in the grid cell specified by `indices`).

Computes the sums

```math
\\begin{align*}
    S_1 & \\approx \\sum\\limits_\\alpha \\left(m_\\alpha \\left|b_{\\mathrm{w}, \\alpha}\\right|\\right)^2 f_\\alpha,\\\\
    S_2 & \\approx \\sum\\limits_\\alpha \\left(m_\\alpha \\left|b_{\\mathrm{w}, \\alpha}\\right| \\left|\\boldsymbol{k}_\\alpha\\right|\\right)^2 f_\\alpha,
\\end{align*}
```

where ``\\boldsymbol{k}_\\alpha = \\left(k_\\alpha, l_\\alpha, m_\\alpha\\right)^\\mathrm{T}`` is the wavevector,

```math
f_\\alpha = \\max \\left(1, \\frac{\\Delta x_\\alpha}{\\Delta \\widehat{x}}\\right) \\max \\left(1, \\frac{\\Delta y_\\alpha}{\\Delta \\widehat{y}}\\right) \\max \\left(1, \\frac{\\Delta z_\\alpha}{J \\Delta \\widehat{z}}\\right)
```

is the maximum grid-cell fraction that can be covered by each ray volume (with ``\\left(\\Delta x_\\alpha, \\Delta y_\\alpha, \\Delta z_\\alpha\\right)`` being the ray-volume extents in physical space) and

```math
\\left|b_{\\mathrm{w}, \\alpha}\\right|^2 = \\frac{2}{\\overline{\\rho}} \\frac{N_\\alpha^4 \\left(k_\\alpha^2 + l_\\alpha^2\\right)}{\\widehat{\\omega}_\\alpha \\left|\\boldsymbol{k}_\\alpha\\right|^2} \\mathcal{N}_\\alpha \\Delta k_\\alpha \\Delta l_\\alpha \\Delta m_\\alpha
```

is the squared gravity-wave amplitude of the buoyancy. Therein, ``\\overline{\\rho}`` is the background density, ``N_\\alpha^2`` is the squared buoyancy frequency interpolated to the ray-volume position (using `interpolate_stratification`), ``\\widehat{\\omega}_\\alpha`` is the intrinsic frequency (calculated with `compute_intrinsic_frequency`), ``\\mathcal{N}_\\alpha`` is the phase-space wave-action density and ``\\left(\\Delta k_\\alpha, \\Delta l_\\alpha, \\Delta m_\\alpha\\right)`` are the ray-volume extents in spectral space.

# Arguments

  - `state`: Model state.
  - `indices`: Grid-cell indices.

# Returns

  - `::AbstractFloat`: Discretized saturation integral ``S_1``.
  - `::AbstractFloat`: Discretized saturation integral ``S_2``.

# See also

  - [`PinCFlow.MSGWaM.Interpolation.get_next_half_level`](@ref)
  - [`PinCFlow.MSGWaM.RayOperations.compute_intrinsic_frequency`](@ref)
  - [`PinCFlow.MSGWaM.Interpolation.interpolate_stratification`](@ref)
"""
function compute_saturation_integrals end

function compute_saturation_integrals(
    state::State,
    indices::NTuple{3, <:Integer},
)
    (; domain, grid) = state
    (; sizex, sizey) = state.namelists.domain
    (; io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, dz, jac) = grid
    (; rhostrattfc) = state.atmosphere
    (; nray, rays) = state.wkb

    # Get indices.
    (ixrv, jyrv, kzrv) = indices

    # Initialize Integrals.
    mb2 = 0.0
    mb2k2 = 0.0

    # Loop over ray volumes.
    for iray in 1:nray[ixrv, jyrv, kzrv]

        # Skip ray volumes with zero wave-action density.
        if rays.dens[iray, ixrv, jyrv, kzrv] == 0
            continue
        end

        xr = rays.x[iray, ixrv, jyrv, kzrv]
        yr = rays.y[iray, ixrv, jyrv, kzrv]
        zr = rays.z[iray, ixrv, jyrv, kzrv]

        dxr = rays.dxray[iray, ixrv, jyrv, kzrv]
        dyr = rays.dyray[iray, ixrv, jyrv, kzrv]
        dzr = rays.dzray[iray, ixrv, jyrv, kzrv]

        if sizex > 1
            ix = floor(Int, (xr - lx[1]) / dx) + i0 - io
        else
            ix = i0
        end

        if sizey > 1
            jy = floor(Int, (yr - ly[1]) / dy) + j0 - jo
        else
            jy = j0
        end

        kz = get_next_half_level(ix, jy, zr, domain, grid)

        n2r = interpolate_stratification(zr, state, N2())

        wnrk = rays.k[iray, ixrv, jyrv, kzrv]
        wnrl = rays.l[iray, ixrv, jyrv, kzrv]
        wnrm = rays.m[iray, ixrv, jyrv, kzrv]

        wnrhs = wnrk^2 + wnrl^2

        dwnrk = rays.dkray[iray, ixrv, jyrv, kzrv]
        dwnrl = rays.dlray[iray, ixrv, jyrv, kzrv]
        dwnrm = rays.dmray[iray, ixrv, jyrv, kzrv]

        omir = compute_intrinsic_frequency(state, (iray, ixrv, jyrv, kzrv))

        densr = rays.dens[iray, ixrv, jyrv, kzrv]

        dzi = min(dzr, jac[ix, jy, kz] * dz)
        facpsp = dzi / jac[ix, jy, kz] / dz * dwnrm

        if sizex > 1
            dxi = min(dxr, dx)
            facpsp = facpsp * dxi / dx * dwnrk
        end

        if sizey > 1
            dyi = min(dyr, dy)
            facpsp = facpsp * dyi / dy * dwnrl
        end

        integral1 = wnrhs * wnrm^2 / ((wnrhs + wnrm^2) * omir) * facpsp

        mb2 += 2 * n2r^2 / rhostrattfc[ix, jy, kz] * densr * integral1

        integral2 = wnrhs * wnrm^2 / omir * facpsp

        mb2k2 += 2 * n2r^2 / rhostrattfc[ix, jy, kz] * densr * integral2
    end

    # Return the results.
    return (mb2, mb2k2)
end
