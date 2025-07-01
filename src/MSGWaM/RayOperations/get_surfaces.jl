"""
```julia
get_surfaces(rays::Rays, indices::NTuple{4, <:Integer}) ->
    Tuple{AbstractFloat, AbstractFloat, AbstractFloat}
```

Compute phase space surface areas for a ray volume.

Calculates the areas of the three phase space surfaces that bound a ray volume
in the 6D phase space (x, y, z, k, l, m). These surfaces represent the
action-angle canonical conjugate pairs.

# Arguments

  - `rays::Rays`: Ray volume data structure
  - `indices::NTuple{4, <:Integer}`: Ray indices (iray, ix, jy, kz)

# Returns

  - `Tuple{AbstractFloat, AbstractFloat, AbstractFloat}`: Surface areas (axk, ayl, azm)

      + `axk`: X-K surface area = dx × dk
      + `ayl`: Y-L surface area = dy × dl
      + `azm`: Z-M surface area = dz × dm

# Phase Space Interpretation

In Hamiltonian mechanics, these surfaces represent:

  - **X-K surface**: Position-wavenumber conjugate pair in x-direction
  - **Y-L surface**: Position-wavenumber conjugate pair in y-direction
  - **Z-M surface**: Position-wavenumber conjugate pair in z-direction

# Conservation Properties

These surface areas are conserved quantities (adiabatic invariants) during
ray propagation in slowly varying media, making them important for:

  - Wave action conservation
  - Energy transport calculations
  - Phase space volume preservation
  - Ray merging/splitting operations

# Applications

Used in:

  - Ray volume initialization
  - Wave action density calculations
  - Momentum flux computations
  - Ray merging algorithms
  - Phase space volume conservation checks

# Example

For a ray with spatial extent (1km, 2km, 0.5km) and spectral extent
(0.001, 0.0005, 0.002) rad/m, the surface areas would be:

  - axk = 1000 × 0.001 = 1.0
  - ayl = 2000 × 0.0005 = 1.0
  - azm = 500 × 0.002 = 1.0
"""
function get_surfaces(rays::Rays, indices::NTuple{4, <:Integer})
    (dxr, dyr, dzr) = get_physical_extent(rays, indices)
    (dkr, dlr, dmr) = get_spectral_extent(rays, indices)
    return (dxr * dkr, dyr * dlr, dzr * dmr)
end
