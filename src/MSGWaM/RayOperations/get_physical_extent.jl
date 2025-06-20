"""
    get_physical_extent(rays::Rays, indices::NTuple{4, <:Integer}) -> Tuple{AbstractFloat, AbstractFloat, AbstractFloat}

Get the physical spatial extents of a ray volume.

Retrieves the spatial dimensions of a ray volume in physical coordinates,
representing the size of the wave packet in each spatial direction.

# Arguments

  - `rays::Rays`: Ray volume data structure containing all ray properties
  - `indices::NTuple{4, <:Integer}`: Ray indices (iray, ix, jy, kz)

# Returns

  - `Tuple{AbstractFloat, AbstractFloat, AbstractFloat}`: Physical extents (dxr, dyr, dzr)

      + `dxr`: Extent in x-direction (zonal) [length units]
      + `dyr`: Extent in y-direction (meridional) [length units]
      + `dzr`: Extent in z-direction (vertical) [length units]

# Physical Interpretation

These extents represent the spatial "size" of the wave packet:

  - **Horizontal extents**: Lateral spreading of wave energy
  - **Vertical extent**: Vertical wavelength scale
  - **Evolution**: Extents can change due to dispersion and shear

# Ray Volume Geometry

Ray volumes are represented as rectangular boxes in physical space:

  - Center position: `(xr, yr, zr)`
  - Extent: `(dxr, dyr, dzr)`
  - Boundaries: `[xr±dxr/2, yr±dyr/2, zr±dzr/2]`

# Applications

Used for:

  - Spatial overlap calculations with grid cells
  - Ray-grid interpolation weights
  - Splitting criteria (when extent exceeds grid spacing)
  - Boundary condition applications
  - Phase space volume computations

# Units

Typically in the same units as the computational grid:

  - Horizontal: meters or kilometers
  - Vertical: meters (may be stretched in terrain-following coordinates)

# Evolution

Extents evolve according to:

  - `d(dx)/dt = ∂cgx/∂x · dx` (differential stretching)
  - Similar equations for dy and dz
  - Can grow or shrink depending on background flow gradients
"""
function get_physical_extent(rays::Rays, indices::NTuple{4, <:Integer})
    return (
        rays.dxray[indices...],
        rays.dyray[indices...],
        rays.dzray[indices...],
    )
end
