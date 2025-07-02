"""
    interpolate_sponge(xlc::AbstractFloat, ylc::AbstractFloat, zlc::AbstractFloat, state::State) -> AbstractFloat

Interpolate sponge layer damping coefficient to arbitrary 3D position.

Computes the wave damping coefficient for sponge layers at any point in the
computational domain through trilinear interpolation. Sponge layers prevent
spurious wave reflections from domain boundaries and provide controlled
wave absorption in specified regions.

# Arguments

  - `xlc, ylc, zlc::AbstractFloat`: Target coordinates for interpolation
  - `state::State`: Complete simulation state containing sponge configuration

# Returns

  - `AbstractFloat`: Sponge damping coefficient [dimensionless] at target position

# Physical Background

Sponge layers implement exponential damping:

  - **Wave amplitude**: `A(t) = A₀ * exp(-α * t)` where α is damping coefficient
  - **Damping timescale**: `τ = 1/α` determines absorption rate
  - **Spatial variation**: Coefficient varies smoothly across domain

# Sponge Layer Design

## Purpose

  - **Boundary absorption**: Prevent wave reflection from domain edges
  - **Numerical stability**: Remove spurious high-frequency modes
  - **Physical realism**: Mimic atmospheric dissipation processes
  - **Computational efficiency**: Avoid expensive boundary conditions

## Implementation

  - **Gradual onset**: Smooth spatial transition to avoid numerical artifacts
  - **Direction-dependent**: Different strengths in x, y, z directions
  - **Unified formulation**: Single coefficient for all wave components

# Algorithm

 1. **Grid Point Location**: Find 8 surrounding grid points in 3D space
 2. **Coefficient Extraction**: Get sponge values from `alphaunifiedsponge` array
 3. **Coordinate Determination**: Compute vertical levels for each horizontal point
 4. **Trilinear Interpolation**: Use general `interpolate` function for final value

# Grid Point Identification

For each dimension:

  - **X-direction**: Floor operation to find left/right grid indices
  - **Y-direction**: Floor operation to find backward/forward indices
  - **Z-direction**: Use `get_next_level` for proper vertical level assignment

# Coordinate System Handling

  - **Domain dimensions**: Automatically handles 1D, 2D, and 3D configurations
  - **Terrain-following**: Uses `ztfc` coordinates for vertical interpolation
  - **Staggered grids**: Accounts for different grid staggering patterns

# Applications

## Wave Damping

Used in ray propagation for exponential amplitude reduction:

```julia
A_new = A_old * exp(-2 * α_sponge * dt)
```

## Boundary Zones

  - **Upper atmosphere**: Prevent wave reflection from model top
  - **Lateral boundaries**: Absorb outgoing waves in limited-area models
  - **Near-surface**: Optional damping of unrealistic wave modes

# Performance Optimization

  - **Efficient interpolation**: Reuses general interpolation infrastructure
  - **Minimal overhead**: Fast coefficient lookup during ray propagation
  - **Memory efficient**: Single storage array for all sponge coefficients

# Numerical Considerations

    # Dermine closest points in horizontal direction.

  - **Smooth transitions**: Avoids discontinuous jumps in damping rate
  - **Conservation**: Maintains overall wave energy budget accounting
  - **Stability**: Prevents numerical instabilities from excessive damping
"""
function interpolate_sponge(
    xlc::AbstractFloat,
    ylc::AbstractFloat,
    zlc::AbstractFloat,
    state::State,
)
    (; namelists, domain, grid) = state
    (; sizex, sizey) = namelists.domain
    (; io, jo, i0, j0) = domain
    (; lx, ly, dx, dy, x, y, ztfc) = grid
    (; alphaunifiedsponge) = state.sponge

    # Dermine closest points in horizontal direction.
    if sizex > 1
        ixl = floor(Int, (xlc - lx[1] - dx / 2) / dx) + i0 - io
        ixr = ixl + 1
    else
        ixl = i0
        ixr = i0
    end
    xl = x[io + ixl]
    xr = x[io + ixr]

    # Determine closest points in meridional direction.
    if sizey > 1
        jyb = floor(Int, (ylc - ly[1] - dy / 2) / dy) + j0 - jo
        jyf = jyb + 1
    else
        jyb = j0
        jyf = j0
    end
    yb = y[jo + jyb]
    yf = y[jo + jyf]

    # Determine closest points in vertical direction and set interpolation
    # values.

    kzlbu = get_next_level(ixl, jyb, zlc, domain, grid)
    kzlbd = kzlbu - 1
    zlbd = ztfc[ixl, jyb, kzlbd]
    zlbu = ztfc[ixl, jyb, kzlbu]

    kzlfu = get_next_level(ixl, jyf, zlc, domain, grid)
    kzlfd = kzlfu - 1
    zlfd = ztfc[ixl, jyf, kzlfd]
    zlfu = ztfc[ixl, jyf, kzlfu]

    kzrbu = get_next_level(ixr, jyb, zlc, domain, grid)
    kzrbd = kzrbu - 1
    zrbd = ztfc[ixr, jyb, kzrbd]
    zrbu = ztfc[ixr, jyb, kzrbu]

    kzrfu = get_next_level(ixr, jyf, zlc, domain, grid)
    kzrfd = kzrfu - 1
    zrfd = ztfc[ixr, jyf, kzrfd]
    zrfu = ztfc[ixr, jyf, kzrfu]

    philbd = alphaunifiedsponge[ixl, jyb, kzlbd]
    philbu = alphaunifiedsponge[ixl, jyb, kzlbu]

    philfd = alphaunifiedsponge[ixl, jyf, kzlfd]
    philfu = alphaunifiedsponge[ixl, jyf, kzlfu]

    phirbd = alphaunifiedsponge[ixr, jyb, kzrbd]
    phirbu = alphaunifiedsponge[ixr, jyb, kzrbu]

    phirfd = alphaunifiedsponge[ixr, jyf, kzrfd]
    phirfu = alphaunifiedsponge[ixr, jyf, kzrfu]

    # Interpolate.
    phi = interpolate(
        namelists;
        philbd,
        philbu,
        philfd,
        philfu,
        phirbd,
        phirbu,
        phirfd,
        phirfu,
        zlbd,
        zlbu,
        zlfd,
        zlfu,
        zrbd,
        zrbu,
        zrfd,
        zrfu,
        zlc,
        yb,
        yf,
        ylc,
        xl,
        xr,
        xlc,
    )

    return phi
end
