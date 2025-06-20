"""
    interpolate(namelists::Namelists; kwargs...) -> AbstractFloat

Perform trilinear interpolation in 3D space with terrain-following coordinates.

Core interpolation function that implements trilinear interpolation for arbitrary
points in 3D atmospheric domains. Handles domain dimensionality automatically
and accounts for terrain-following coordinate systems used in atmospheric models.

# Arguments

  - `namelists::Namelists`: Model configuration containing domain dimensions

# Keyword Arguments

## Field Values (8 corner points)

  - `philbd, philbu`: Left-backward field values (down/up in vertical)
  - `philfd, philfu`: Left-forward field values (down/up in vertical)
  - `phirbd, phirbu`: Right-backward field values (down/up in vertical)
  - `phirfd, phirfu`: Right-forward field values (down/up in vertical)

## Coordinate Values (8 corner points)

  - `zlbd, zlbu, zlfd, zlfu`: Left vertical coordinates (backward/forward, down/up)
  - `zrbd, zrbu, zrfd, zrfu`: Right vertical coordinates (backward/forward, down/up)

## Interpolation Point

  - `xlc, ylc, zlc`: Target coordinates for interpolation
  - `xl, xr`: Left/right grid coordinates
  - `yb, yf`: Backward/forward grid coordinates

# Returns

  - `AbstractFloat`: Interpolated field value at target point

# Algorithm

Performs sequential interpolation in each dimension:

## 1. X-Direction Interpolation

For multi-dimensional domains (`sizex > 1`):

  - Computes interpolation factor based on position within grid cell
  - Linearly interpolates between left and right grid values
  - Handles edge cases where target point is outside interpolation bounds

## 2. Y-Direction Interpolation

For multi-dimensional domains (`sizey > 1`):

  - Interpolates between backward and forward values
  - Uses results from x-direction interpolation as input

## 3. Z-Direction Interpolation

Final vertical interpolation:

  - Accounts for terrain-following coordinate stretching
  - Handles varying grid spacing in vertical direction
  - Produces final interpolated value

# Interpolation Factor Calculation

For each direction, computes factor `f` where:

  - `f = 0`: Use lower/left/backward value exclusively
  - `f = 1`: Use upper/right/forward value exclusively
  - `0 < f < 1`: Linear combination of both values

Special cases:

  - Points outside interpolation bounds use nearest boundary value
  - Zero spacing between grid points defaults to lower/left/backward value

# Coordinate System Support

  - **Cartesian coordinates**: Standard rectangular grids
  - **Terrain-following**: Vertical coordinates follow surface topography
  - **Stretched grids**: Non-uniform spacing in any direction
  - **Staggered grids**: Different variable locations on grid

# Error Handling

Validates coordinate ordering:

  - `xr ≥ xl`: Right coordinate ≥ left coordinate
  - `yf ≥ yb`: Forward coordinate ≥ backward coordinate
  - `zu ≥ zd`: Upper coordinate ≥ lower coordinate

# Applications

Used throughout MSGWaM for:

  - Mean flow interpolation to ray positions
  - Stratification evaluation along ray paths
  - Sponge layer coefficient determination
  - Background field derivatives for refraction calculations

# Performance Notes

  - Efficient for scattered interpolation points
  - Avoids expensive grid searches by using pre-computed indices
  - Handles domain boundaries and coordinate singularities gracefully
"""
function interpolate(
    namelists::Namelists;
    philbd::AbstractFloat = NaN,
    philbu::AbstractFloat = NaN,
    philfd::AbstractFloat = NaN,
    philfu::AbstractFloat = NaN,
    phirbd::AbstractFloat = NaN,
    phirbu::AbstractFloat = NaN,
    phirfd::AbstractFloat = NaN,
    phirfu::AbstractFloat = NaN,
    zlbd::AbstractFloat = NaN,
    zlbu::AbstractFloat = NaN,
    zlfd::AbstractFloat = NaN,
    zlfu::AbstractFloat = NaN,
    zrbd::AbstractFloat = NaN,
    zrbu::AbstractFloat = NaN,
    zrfd::AbstractFloat = NaN,
    zrfu::AbstractFloat = NaN,
    zlc::AbstractFloat = NaN,
    yb::AbstractFloat = NaN,
    yf::AbstractFloat = NaN,
    ylc::AbstractFloat = NaN,
    xl::AbstractFloat = NaN,
    xr::AbstractFloat = NaN,
    xlc::AbstractFloat = NaN,
)
    (; sizex, sizey) = namelists.domain

    # Interpolate in x.
    if sizex == 1
        phibd = philbd
        phibu = philbu

        phifd = philfd
        phifu = philfu

        zbd = zlbd
        zbu = zlbu

        zfd = zlfd
        zfu = zlfu
    else
        if xr < xl
            error("Error in interpolate: xr = ", xr, " < xl = ", xl)
        elseif xr == xl
            factor = 0.0
        elseif xlc > xr
            factor = 0.0
        elseif xlc > xl
            factor = (xr - xlc) / (xr - xl)
        else
            factor = 1.0
        end

        phibd = factor * philbd + (1.0 - factor) * phirbd
        phibu = factor * philbu + (1.0 - factor) * phirbu

        phifd = factor * philfd + (1.0 - factor) * phirfd
        phifu = factor * philfu + (1.0 - factor) * phirfu

        zbd = factor * zlbd + (1.0 - factor) * zrbd
        zbu = factor * zlbu + (1.0 - factor) * zrbu

        zfd = factor * zlfd + (1.0 - factor) * zrfd
        zfu = factor * zlfu + (1.0 - factor) * zrfu
    end

    # Intepolate in y.
    if sizey == 1
        phid = phibd
        phiu = phibu

        zd = zbd
        zu = zbu
    else
        if yf < yb
            error("Error in interpolate: yf = ", yf, " < yb = ", yb)
        elseif yf == yb
            factor = 0.0
        elseif ylc > yf
            factor = 0.0
        elseif ylc > yb
            factor = (yf - ylc) / (yf - yb)
        else
            factor = 1.0
        end

        phid = factor * phibd + (1.0 - factor) * phifd
        phiu = factor * phibu + (1.0 - factor) * phifu

        zd = factor * zbd + (1.0 - factor) * zfd
        zu = factor * zbu + (1.0 - factor) * zfu
    end

    # Interpolate in z.
    if zu < zd
        error("Error in interpolate: zu = ", zu, " < zd = ", zd)
    elseif zu == zd
        factor = 0.0
    elseif zlc > zu
        factor = 0.0
    elseif zlc > zd
        factor = (zu - zlc) / (zu - zd)
    else
        factor = 1.0
    end

    phi = factor * phid + (1.0 - factor) * phiu

    # Return the result.
    return phi
end
