"""
    compute_intrinsic_frequency(state::State, indices::NTuple{4, <:Integer}) -> AbstractFloat

Compute the intrinsic frequency of a gravity wave ray volume.

Calculates the intrinsic (Doppler-shifted) frequency of a gravity wave based on
its wavenumber vector and the local atmospheric conditions, using the full
gravity wave dispersion relation including rotational effects.

# Arguments

  - `state::State`: Complete simulation state containing atmospheric data
  - `indices::NTuple{4, <:Integer}`: Ray volume indices (iray, ix, jy, kz)

# Returns

  - `AbstractFloat`: Intrinsic frequency ω [s⁻¹]

# Dispersion Relation

Uses the full gravity wave dispersion relation with rotation:

`ω² = N² · k_h² / (k_h² + m²) + f² · m² / (k_h² + m²)`

where:

  - `ω`: Intrinsic frequency
  - `N²`: Brunt-Väisälä frequency squared (local stratification)
  - `k_h`: Horizontal wavenumber magnitude `√(k² + l²)`
  - `m`: Vertical wavenumber
  - `f`: Coriolis parameter

# Physical Interpretation

## Frequency Components

  - **Buoyancy term**: `N² · k_h² / (k_h² + m²)` - restoring force from stratification
  - **Inertial term**: `f² · m² / (k_h² + m²)` - restoring force from rotation

## Limiting Cases

  - **High frequency limit** (`ω ≫ f`): `ω ≈ N · k_h / √(k_h² + m²)` (internal gravity waves)
  - **Low frequency limit** (`ω ≈ f`): Inertia-gravity waves
  - **f → 0**: `ω = N · k_h / √(k_h² + m²)` (non-rotating case)

# Branch Selection

The `branchr` parameter determines the frequency branch:

  - `branchr = +1`: Positive frequency branch (upward phase propagation)
  - `branchr = -1`: Negative frequency branch (downward phase propagation)

# Local Conditions

  - **Stratification**: Interpolated to ray position from background profile
  - **Rotation**: Uses domain-wide Coriolis parameter
  - **Coordinates**: All calculations in terrain-following coordinates

# Applications

Used for:

  - Group velocity calculations
  - Wave action conservation
  - Dispersion relation enforcement
  - Critical level detection
  - Energy/momentum flux computations

# Critical Levels

Wave propagation is only possible when:

  - `|ω| > f` (above inertial frequency)
  - `|ω| < N` (below buoyancy frequency)

Outside this range, waves become evanescent or are absorbed.

# Units

Returns frequency in units consistent with the time normalization,
typically [s⁻¹] or dimensionless if time is scaled.
"""
function compute_intrinsic_frequency(
    state::State,
    indices::NTuple{4, <:Integer},
)
    (; coriolis_frequency) = state.namelists.atmosphere
    (; branchr) = state.namelists.wkb
    (; tref) = state.constants
    (; rays) = state.wkb

    zr = rays.z[indices...]
    kr = rays.k[indices...]
    lr = rays.l[indices...]
    mr = rays.m[indices...]
    khr = sqrt(kr^2 + lr^2)

    n2r = interpolate_stratification(zr, state, N2())
    fc = coriolis_frequency * tref

    return branchr * sqrt(n2r * khr^2 + fc^2 * mr^2) / sqrt(khr^2 + mr^2)
end
