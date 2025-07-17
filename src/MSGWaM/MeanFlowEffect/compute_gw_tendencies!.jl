"""
    compute_gw_tendencies!(state)

Compute gravity wave momentum and heating tendencies from wave stress divergence.

Calculates the acceleration and heating of the mean flow due to gravity wave
momentum transport and Eliassen-Palm flux divergence, based on previously
computed wave momentum flux integrals.

# Arguments

  - `state`: Complete simulation state containing GW integrals and atmospheric fields

# Physical Theory

## Momentum Equation

The gravity wave momentum tendencies represent the deceleration of the mean flow
due to momentum flux divergence:

`∂u/∂t = -∇ · (momentum flux) + Eliassen-Palm force`

For each momentum component:

  - **Zonal**: `∂u/∂t = -∂(uw)/∂z - ∂(uu)/∂x - ∂(uv)/∂y + Fx`
  - **Meridional**: `∂v/∂t = -∂(vw)/∂z - ∂(uv)/∂x - ∂(vv)/∂y + Fy`

## Eliassen-Palm Forces

When rotation is important (`fc ≠ 0`):

  - `Fx = ρ · fc² · N² · k · m / (ρ₀ · g · ω · (k² + l² + m²))`
  - `Fy = ρ · fc² · N² · l · m / (ρ₀ · g · ω · (k² + l² + m²))`

## Heating Terms

Potential temperature tendencies from EP flux divergence:

  - `∂θ/∂t = θ₀/fc · [∂(vθ)/∂y + ∂(uθ)/∂x + met₁₃·∂(uθ)/∂z + met₂₃·∂(vθ)/∂z]`

# Coordinate System

Accounts for terrain-following coordinates using metric tensor:

  - `met[i,j,k,1,3]`: X-Z metric component
  - `met[i,j,k,2,3]`: Y-Z metric component
  - `met[i,j,k,3,3]`: Z-Z metric component

# Grid Staggering

  - Momentum fluxes: Computed at appropriate staggered locations
  - Derivatives: Centered differences with appropriate metric factors
  - Density weighting: Uses total density `ρ = ρ' + ρ₀`

# Vertical Derivatives

Uses terrain-following coordinate system:

  - `∂/∂z` → `(1/jac) · ∂/∂ζ` where `ζ` is the terrain-following coordinate
  - Grid spacing corrections for stretched coordinates

# Domain Restrictions

Only computes tendencies above minimum WKB altitude:

  - `z > lz[1] + zmin_wkb_dim / lref`
  - Below this level, tendencies remain zero

# Output Fields

Updates the following tendency fields:

  - `tendencies.dudt[i,j,k]`: Zonal wind tendency

  - `tendencies.dvdt[i,j,k]`: Meridional wind tendency
  - `tendencies.dthetadt[i,j,k]`: Potential temperature tendency

    # Set the Coriolis parameter.

# Conservation Properties

  - Momentum: Conserved through proper flux divergence calculation
  - Energy: Heating terms maintain thermodynamic consistency
  - Angular momentum: Preserved when rotation effects included
"""
function compute_gw_tendencies!(state)
    (; sizex, sizey) = state.namelists.domain
    (; coriolis_frequency) = state.namelists.atmosphere
    (; zmin_wkb_dim) = state.namelists.wkb
    (; tref, lref) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; lz, dx, dy, dz, ztfc, jac, met) = state.grid
    (; rhostrattfc) = state.atmosphere
    (; rho) = state.variables.predictands
    (; integrals, tendencies) = state.wkb

    # Set the Coriolis parameter.
    fc = coriolis_frequency * tref

    for field in fieldnames(GWTendencies)
        getfield(tendencies, field) .= 0.0
    end

    for kz in k0:k1, jy in j0:j1, ix in i0:i1
        if ztfc[ix, jy, kz] < lz[1] + zmin_wkb_dim / lref
            continue
        end

        rhotot = rho[ix, jy, kz] + rhostrattfc[ix, jy, kz]

        # Compute the drag on the zonal wind.

        tendencies.dudt[ix, jy, kz] =
            -rhotot / rhostrattfc[ix, jy, kz] / jac[ix, jy, kz] *
            (integrals.uw[ix, jy, kz + 1] - integrals.uw[ix, jy, kz - 1]) /
            (2.0 * dz)

        if sizex > 1
            tendencies.dudt[ix, jy, kz] -=
                rhotot / rhostrattfc[ix, jy, kz] * (
                    (
                        integrals.uu[ix + 1, jy, kz] -
                        integrals.uu[ix - 1, jy, kz]
                    ) / (2.0 * dx) +
                    met[ix, jy, kz, 1, 3] * (
                        integrals.uu[ix, jy, kz + 1] -
                        integrals.uu[ix, jy, kz - 1]
                    ) / (2.0 * dz)
                )
        end

        if sizey > 1
            tendencies.dudt[ix, jy, kz] -=
                rhotot / rhostrattfc[ix, jy, kz] * (
                    (
                        integrals.uv[ix, jy + 1, kz] -
                        integrals.uv[ix, jy - 1, kz]
                    ) / (2.0 * dy) +
                    met[ix, jy, kz, 2, 3] * (
                        integrals.uv[ix, jy, kz + 1] -
                        integrals.uv[ix, jy, kz - 1]
                    ) / (2.0 * dz)
                )
        end

        tendencies.dudt[ix, jy, kz] += rhotot * integrals.etx[ix, jy, kz]

        # Compute the drag on the meridional wind.

        tendencies.dvdt[ix, jy, kz] =
            -rhotot / rhostrattfc[ix, jy, kz] / jac[ix, jy, kz] *
            (integrals.vw[ix, jy, kz + 1] - integrals.vw[ix, jy, kz - 1]) /
            (2.0 * dz)

        if sizex > 1
            tendencies.dvdt[ix, jy, kz] -=
                rhotot / rhostrattfc[ix, jy, kz] * (
                    (
                        integrals.uv[ix + 1, jy, kz] -
                        integrals.uv[ix - 1, jy, kz]
                    ) / (2.0 * dx) +
                    met[ix, jy, kz, 1, 3] * (
                        integrals.uv[ix, jy, kz + 1] -
                        integrals.uv[ix, jy, kz - 1]
                    ) / (2.0 * dz)
                )
        end

        if sizey > 1
            tendencies.dvdt[ix, jy, kz] -=
                rhotot / rhostrattfc[ix, jy, kz] * (
                    (
                        integrals.vv[ix, jy + 1, kz] -
                        integrals.vv[ix, jy - 1, kz]
                    ) / (2.0 * dy) +
                    met[ix, jy, kz, 2, 3] * (
                        integrals.vv[ix, jy, kz + 1] -
                        integrals.vv[ix, jy, kz - 1]
                    ) / (2.0 * dz)
                )
        end

        tendencies.dvdt[ix, jy, kz] += rhotot * integrals.ety[ix, jy, kz]

        # Compute the heating.

        if fc != 0.0 && (sizex > 1 || sizey > 1)
            if sizex > 1
                tendencies.dthetadt[ix, jy, kz] +=
                    rhotot * (
                        (
                            integrals.utheta[ix + 1, jy, kz] -
                            integrals.utheta[ix - 1, jy, kz]
                        ) / (2.0 * dx) +
                        met[ix, jy, kz, 1, 3] * (
                            integrals.utheta[ix, jy, kz + 1] -
                            integrals.utheta[ix, jy, kz - 1]
                        ) / (2.0 * dz)
                    )
            end

            if sizey > 1
                tendencies.dthetadt[ix, jy, kz] +=
                    rhotot * (
                        (
                            integrals.vtheta[ix, jy + 1, kz] -
                            integrals.vtheta[ix, jy - 1, kz]
                        ) / (2.0 * dy) +
                        met[ix, jy, kz, 2, 3] * (
                            integrals.vtheta[ix, jy, kz + 1] -
                            integrals.vtheta[ix, jy, kz - 1]
                        ) / (2.0 * dz)
                    )
            end
        end

        compute_leading_order_tracer_forcing!(
            state,
            ix,
            jy,
            kz,
            state.namelists.tracer.tracersetup,
        )
    end
end