"""
```julia
compute_saturation_integrals(state::State, indices::NTuple{3, <:Integer}) ->
    Tuple{AbstractFloat, AbstractFloat}
```

Compute wave saturation integrals for a grid cell.

Calculates the momentum flux and diffusion integrals needed for the wave
saturation scheme, which parameterizes wave breaking when wave amplitudes
exceed critical values.

# Arguments

  - `state::State`: Complete simulation state
  - `indices::NTuple{3, <:Integer}`: Grid cell indices (ix, jy, kz)

# Returns

  - `Tuple{AbstractFloat, AbstractFloat}`: Saturation integrals (mb2, mb2k2)

      + `mb2`: Total momentum flux integral `∫ ρ·A·ω·(k²+l²+m²)`
      + `mb2k2`: Diffusion-weighted integral for saturation calculation

# Physical Theory

## Wave Breaking Criterion

Convective instability occurs when potential temperature perturbations exceed:
`|θ'| > α_sat · θ₀ · N/g`

## Momentum Flux Integral

`mb2 = ∫ (2N²/ρ₀) · A · (k²m²/(k²+l²+m²)) · (1/ω) · dV_phasespace`

where:

  - `N²`: Brunt-Väisälä frequency squared
  - `ρ₀`: Background density
  - `A`: Wave action density
  - `k²,l²,m²`: Wavenumber components
  - `ω`: Intrinsic frequency
  - `dV_phasespace`: Phase space volume element

## Diffusion Integral

`mb2k2`: Similar to mb2 but weighted for diffusion calculations

# Algorithm

 1. **Ray Loop**: Iterate over all ray volumes in the grid cell
 2. **Position Mapping**: Determine which grid cell each ray occupies
 3. **Phase Space Factors**: Compute spatial and spectral overlap fractions
 4. **Local Properties**: Interpolate stratification and compute frequency
 5. **Integral Contributions**: Add weighted contributions from each ray

# Saturation Scheme

Used to compute turbulent diffusion coefficient:

```julia
if mb2 > α_sat² · N²
    κ = (mb2 - α_sat² · N²) / (2 · dt · mb2k2)
else
    κ = 0
end
```

# Phase Space Weighting

Accounts for:

  - Ray volume overlap with grid cell
  - Spectral resolution factors
  - Coordinate system metrics (terrain-following)

# Applications

  - Wave breaking parameterization
  - Turbulent mixing calculations
  - Energy dissipation estimates
  - Momentum flux saturation

# Grid Cell Association

Rays are associated with grid cells based on their center positions,
with appropriate interpolation for rays spanning multiple cells.
"""
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
