"""
    compute_derivatives(state::State, indices::NTuple{4, <:Integer}, phitype::DUDX) -> Tuple{AbstractFloat, AbstractFloat}

Compute zonal derivative of zonal wind (∂u/∂x) at specified grid levels.

Calculates the horizontal gradient of zonal wind using finite differences
with terrain-following coordinate corrections. This derivative is essential
for gravity wave refraction calculations and represents horizontal wind shear.

# Arguments

  - `state::State`: Complete simulation state containing wind fields and grid
  - `indices::NTuple{4, <:Integer}`: Grid indices (ix, jy, kzd, kzu) for upper/lower levels
  - `phitype::DUDX`: Type dispatch for zonal wind, zonal derivative

# Returns

  - `Tuple{AbstractFloat, AbstractFloat}`: (∂u/∂x_lower, ∂u/∂x_upper) at the two levels

# Physical Significance

∂u/∂x represents:

  - **Horizontal wind shear**: Rate of zonal wind change in east-west direction
  - **Flow deformation**: Stretching/compression of fluid elements
  - **Wave refraction**: Causes bending of gravity wave ray paths
  - **Vorticity**: Contributes to horizontal vorticity component

# Finite Difference Formula

For terrain-following coordinates:

```
∂u/∂x = (u[i] - u[i-1])/dx + met[1,3] * (∂u/∂ζ)
```

where:

  - First term: Standard horizontal difference
  - Second term: Terrain-following correction using metric tensor
  - `ζ`: Terrain-following vertical coordinate

# Metric Tensor Correction

The `met[ix,jy,kz,1,3]` term accounts for:

  - **Coordinate coupling**: x-z coordinate system interactions
  - **Terrain influence**: Topographic effects on horizontal derivatives
  - **Grid stretching**: Non-orthogonal coordinate transformations

# Grid Staggering

  - **U-grid points**: Zonal wind located at cell faces in x-direction
  - **Horizontal spacing**: Uses `dx` between adjacent u-points
  - **Vertical averaging**: Terrain correction uses average over vertical levels

# Applications in Wave Dynamics

  - **Ray refraction**: ∂u/∂x appears in wavenumber evolution equation dk/dt
  - **Group velocity**: Modifies horizontal wave propagation characteristics
  - **Critical layers**: Large shears can create wave absorption regions
  - **Instability**: Strong shears trigger wave breaking mechanisms
"""
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

"""
    compute_derivatives(state::State, indices::NTuple{4, <:Integer}, phitype::DUDY) -> Tuple{AbstractFloat, AbstractFloat}

Compute meridional derivative of zonal wind (∂u/∂y).

Calculates the north-south gradient of zonal wind using finite differences
with appropriate metric corrections for terrain-following coordinates.

# Arguments

  - `state::State`: Simulation state with wind and grid data
  - `indices::NTuple{4, <:Integer}`: Grid level indices (ix, jy, kzd, kzu)
  - `phitype::DUDY`: Type dispatch for zonal wind, meridional derivative

# Returns

  - `Tuple{AbstractFloat, AbstractFloat}`: (∂u/∂y_lower, ∂u/∂y_upper)

# Physical Interpretation

∂u/∂y represents:

  - **Cross-flow shear**: Zonal wind variation in meridional direction
  - **Baroclinic shear**: Often associated with temperature gradients
  - **Coriolis coupling**: Interacts with Coriolis force in momentum equations
  - **Wave propagation**: Affects meridional wave refraction

# Finite Difference Implementation

```
∂u/∂y = (u[j+1] - u[j])/dy + <met[2,3]> * (∂u/∂ζ)
```

where `<met[2,3]>` is spatially averaged metric tensor component.

# Metric Averaging

Since u-points don't coincide with cell centers, requires averaging:

  - **Four-point average**: Over surrounding cell centers
  - **Spatial weighting**: Equal weights for regular grids
  - **Coordinate consistency**: Maintains second-order accuracy
"""
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

"""
    compute_derivatives(state::State, indices::NTuple{4, <:Integer}, phitype::DUDZ) -> Tuple{AbstractFloat, AbstractFloat}

Compute vertical derivative of zonal wind (∂u/∂z).

Calculates vertical wind shear using finite differences in terrain-following
coordinates with proper handling of topographic boundaries and coordinate
system metrics.

# Arguments

  - `state::State`: Simulation state
  - `indices::NTuple{4, <:Integer}`: Grid indices (ix, jy, kzd, kzu)
  - `phitype::DUDZ`: Type dispatch for vertical zonal wind derivative

# Returns

  - `Tuple{AbstractFloat, AbstractFloat}`: (∂u/∂z_lower, ∂u/∂z_upper)

# Physical Significance

∂u/∂z represents:

  - **Vertical wind shear**: Primary driver of atmospheric turbulence
  - **Thermal wind**: Related to horizontal temperature gradients
  - **Wave breaking**: Strong shears trigger gravity wave dissipation
  - **Momentum transport**: Key for understanding drag forces

# Terrain-Following Coordinates

Vertical derivative in terrain coordinates:

```
∂u/∂z = (1/jac) * ∂u/∂ζ
```

where:

  - `jac`: Jacobian of coordinate transformation
  - `ζ`: Terrain-following coordinate
  - Accounts for grid stretching near topography

# Boundary Conditions

  - **Below ground**: Set derivative to zero
  - **At surface**: Handle partial grid cells appropriately
  - **Domain top**: Extrapolate using available levels
  - **Topography**: Ensure physical consistency near terrain

# Jacobian Averaging

Uses harmonic mean of adjacent Jacobians:

```
jac_effective = 2 * jac1 * jac2 / (jac1 + jac2)
```

Provides stable derivative calculation in stretched coordinates.
"""
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

"""
    compute_derivatives(state::State, indices::NTuple{4, <:Integer}, phitype::DVDX) -> Tuple{AbstractFloat, AbstractFloat}

Compute zonal derivative of meridional wind (∂v/∂x).

Similar to ∂u/∂y but for meridional wind component, representing how
north-south winds vary in the east-west direction.

# Physical Significance

  - **Cross-flow gradient**: Meridional wind change in zonal direction
  - **Geostrophic balance**: Related to pressure gradients via Coriolis force
  - **Wave dynamics**: Affects zonal propagation of meridional wave modes
  - **Vorticity**: Contributes to vertical vorticity component
"""
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

"""
    compute_derivatives(state::State, indices::NTuple{4, <:Integer}, phitype::DVDY) -> Tuple{AbstractFloat, AbstractFloat}

Compute meridional derivative of meridional wind (∂v/∂y).

Represents the primary meridional wind shear component, analogous to
∂u/∂x but for the meridional flow direction.

# Physical Significance

  - **Meridional stretching**: Flow convergence/divergence in y-direction
  - **Temperature advection**: Related to heat transport processes
  - **Storm dynamics**: Important for cyclone and anticyclone development
  - **Wave propagation**: Modifies meridional wave characteristics
"""
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

"""
    compute_derivatives(state::State, indices::NTuple{4, <:Integer}, phitype::DVDZ) -> Tuple{AbstractFloat, AbstractFloat}

Compute vertical derivative of meridional wind (∂v/∂z).

Calculates vertical shear of meridional wind component, representing how
north-south winds change with altitude.

# Physical Significance

  - **Vertical wind shear**: Meridional component of 3D shear vector
  - **Baroclinic instability**: Driver of extratropical storm development
  - **Jet stream structure**: Characterizes upper-level wind patterns
  - **Wave-mean flow interaction**: Couples waves to background meridional flow

# Implementation Notes

Similar structure to ∂u/∂z but uses v-wind grid points and appropriate
metric tensor components for meridional wind staggering.
"""
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
