"""
    interpolate_stratification(zlc::AbstractFloat, state::State, strtype::N2) -> AbstractFloat

Interpolate Brunt-Väisälä frequency squared (N²) to a specified vertical level.

Computes the atmospheric stratification parameter N² at arbitrary vertical
positions through linear interpolation between model grid levels. N² determines
gravity wave propagation characteristics and is fundamental to wave dynamics.

# Arguments

  - `zlc::AbstractFloat`: Target vertical coordinate for interpolation
  - `state::State`: Complete simulation state containing atmospheric data
  - `strtype::N2`: Type dispatch for N² interpolation

# Returns

  - `AbstractFloat`: Brunt-Väisälä frequency squared [s⁻²] at target level

# Physical Background

The Brunt-Väisälä frequency N² characterizes atmospheric stability:

  - `N² = (g/θ) * (dθ/dz)` where g is gravity, θ potential temperature
  - **N² > 0**: Stable stratification (gravity waves can propagate)
  - **N² = 0**: Neutral stratification (marginal stability)
  - **N² < 0**: Unstable stratification (convective instability)

# Algorithm

 1. **Level Location**: Find grid levels bracketing target height using `get_next_level`
 2. **Value Extraction**: Get N² values from atmospheric state at bracket levels
 3. **Linear Interpolation**: Compute weighted average based on vertical position
 4. **Boundary Handling**: Extrapolate using nearest values for points outside domain

# Interpolation Formula

```
N²(z) = f * N²_lower + (1-f) * N²_upper
```

where `f` is the interpolation factor based on relative position between levels.

# Grid Considerations

  - Uses full model levels (`ztfc`) for coordinate reference
  - Accounts for terrain-following coordinate system
  - Handles variable vertical grid spacing

# Applications in Wave Dynamics

  - **Dispersion relation**: N² appears directly in gravity wave frequency
  - **Vertical wavelength**: Determines local vertical scale of waves
  - **Critical levels**: N² = 0 levels where waves are absorbed
  - **Wave breaking**: Large N² gradients can trigger instabilities

# Performance Notes

  - Efficient single-level interpolation for ray tracing
  - Uses pre-computed atmospheric profiles
  - Minimal computational overhead per ray evaluation
"""
function interpolate_stratification(
    zlc::AbstractFloat,
    state::State,
    strtype::N2,
)
    (; domain, grid) = state
    (; bvsstrattfc) = state.atmosphere
    (; i0, j0) = domain
    (; ztfc) = grid

    kzu = get_next_level(i0, j0, zlc, domain, grid)
    kzd = kzu - 1

    zd = ztfc[i0, j0, kzd]
    zu = ztfc[i0, j0, kzu]
    strd = bvsstrattfc[i0, j0, kzd]
    stru = bvsstrattfc[i0, j0, kzu]

    if zu < zd
        error(
            "Error in interpolate_stratification (N2): zu = ",
            zu,
            " < zd = ",
            zd,
        )
    elseif zu == zd
        factor = 0.0
    elseif zlc > zu
        factor = 0.0
    elseif zlc > zd
        factor = (zu - zlc) / (zu - zd)
    else
        factor = 1.0
    end

    str = factor * strd + (1.0 - factor) * stru

    return str
end

"""
    interpolate_stratification(zlc::AbstractFloat, state::State, strtype::DN2DZ) -> AbstractFloat

Interpolate vertical derivative of Brunt-Väisälä frequency squared (dN²/dz).

Computes the vertical gradient of atmospheric stratification at arbitrary
vertical positions. This quantity is crucial for gravity wave refraction
and affects wave propagation characteristics.

# Arguments

  - `zlc::AbstractFloat`: Target vertical coordinate for interpolation
  - `state::State`: Complete simulation state containing atmospheric data
  - `strtype::DN2DZ`: Type dispatch for dN²/dz interpolation

# Returns

  - `AbstractFloat`: Vertical derivative of N² [s⁻² m⁻¹] at target level

# Physical Significance

The vertical gradient dN²/dz affects:

  - **Wave refraction**: Changes in N² bend wave ray paths
  - **Critical level formation**: Sharp N² gradients create wave absorption regions
  - **Vertical group velocity**: Modified by stratification variations
  - **Wave breaking**: Strong gradients can trigger convective instability

# Algorithm

 1. **Level Location**: Find half-levels bracketing target using `get_next_half_level`
 2. **Derivative Computation**: Calculate dN²/dz using centered differences at each level
 3. **Terrain-Following Correction**: Account for coordinate system metrics (`jac`)
 4. **Linear Interpolation**: Weighted average of derivatives at bracket levels

# Derivative Calculation

At each half-level:

```
dN²/dz = (N²_above - N²_below) / (effective_spacing)
```

where effective spacing accounts for terrain-following coordinates:

```
effective_spacing = 2 * jac_lower * jac_upper / (jac_lower + jac_upper) * dz
```

# Terrain-Following Coordinates

  - Uses half-levels (`ztildetfc`) for derivative evaluation
  - Jacobian factors (`jac`) account for coordinate stretching
  - Maintains accuracy in regions with steep topography

# Applications

  - **Ray propagation**: Refraction calculations in `propagate_rays!`
  - **Wavenumber evolution**: Contributes to dk/dt, dl/dt, dm/dt equations
  - **Wave-mean flow interaction**: Background flow modification effects
  - **Instability analysis**: Convective breakdown prediction

# Numerical Considerations

  - Second-order accurate centered differences
  - Proper metric factor weighting for terrain-following grids
  - Stable interpolation even with steep coordinate stretching
"""
function interpolate_stratification(
    zlc::AbstractFloat,
    state::State,
    strtype::DN2DZ,
)
    (; domain, grid) = state
    (; bvsstrattfc) = state.atmosphere
    (; i0, j0) = domain
    (; dz, ztildetfc, jac) = grid

    kzu = get_next_half_level(i0, j0, zlc, domain, grid)
    kzd = kzu - 1

    zd = ztildetfc[i0, j0, kzd]
    zu = ztildetfc[i0, j0, kzu]

    strd =
        (bvsstrattfc[i0, j0, kzd + 1] - bvsstrattfc[i0, j0, kzd]) / (
            2.0 * jac[i0, j0, kzd] * jac[i0, j0, kzd + 1] /
            (jac[i0, j0, kzd] + jac[i0, j0, kzd + 1])
        ) / dz
    stru =
        (bvsstrattfc[i0, j0, kzu + 1] - bvsstrattfc[i0, j0, kzu]) / (
            2.0 * jac[i0, j0, kzu] * jac[i0, j0, kzu + 1] /
            (jac[i0, j0, kzu] + jac[i0, j0, kzu + 1])
        ) / dz

    if zu < zd
        error(
            "Error in interpolate_stratification (DN2DZ): zu = ",
            zu,
            " < zd = ",
            zd,
        )
    elseif zu == zd
        factor = 0.0
    elseif zlc > zu
        factor = 0.0
    elseif zlc > zd
        factor = (zu - zlc) / (zu - zd)
    else
        factor = 1.0
    end

    str = factor * strd + (1.0 - factor) * stru

    return str
end
