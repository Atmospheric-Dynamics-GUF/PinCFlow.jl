"""
```julia
Grid{
    A <: AbstractVector{<:AbstractFloat},
    B <: AbstractFloat,
    C <: AbstractMatrix{<:AbstractFloat},
    D <: AbstractArray{<:AbstractFloat, 3},
    E <: AbstractArray{<:AbstractFloat, 3},
    F <: AbstractArray{<:AbstractFloat, 5},
}
```

Grid information.

# Fields

## Domain Boundaries

  - `lx`: Non-dimensional domain boundaries in x-direction [x_min, x_max]
  - `ly`: Non-dimensional domain boundaries in y-direction [y_min, y_max]
  - `lz`: Non-dimensional domain boundaries in z-direction [z_min, z_max]

## Grid Spacing

  - `dx`: Uniform grid spacing in x-direction
  - `dy`: Uniform grid spacing in y-direction
  - `dz`: Uniform grid spacing in z-direction (before vertical stretching)

## Coordinate Arrays

  - `x`: Cell-centered x-coordinates for the local domain
  - `y`: Cell-centered y-coordinates for the local domain
  - `z`: Cell-centered z-coordinates for the local domain

## Topography

  - `topography_surface`: 2D surface topography field
  - `topography_spectrum`: 3D spectral representation of topography
  - `k_spectrum`: Wavenumber array in x-direction for spectral methods
  - `l_spectrum`: Wavenumber array in y-direction for spectral methods

## Coordinate Transformation

  - `jac`: 3D Jacobian of the terrain-following coordinate transformation
  - `met`: 5D metric tensor (3×3 at each grid point) for coordinate transformation

## Physical Coordinates

  - `ztfc`: Physical z-coordinates at cell centers (terrain-following)
  - `ztildetfc`: Physical z-coordinates at cell faces (terrain-following)

# Terrain-Following Coordinates

The grid uses terrain-following coordinates where:

  - Horizontal coordinates (x,y) remain Cartesian
  - Vertical coordinate ζ follows the terrain: ζ = (z - h(x,y))/(H - h(x,y))
  - Physical height: z_phys = h(x,y) + ζ * (H - h(x,y))

where h(x,y) is surface topography and H is domain height.

# See also

  - [`PinCFlow.Types.FoundationalTypes.compute_topography`](@ref)
  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)
  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)
"""
struct Grid{
    A <: AbstractVector{<:AbstractFloat},
    B <: AbstractFloat,
    C <: AbstractMatrix{<:AbstractFloat},
    D <: AbstractArray{<:AbstractFloat, 3},
    E <: AbstractArray{<:AbstractFloat, 3},
    F <: AbstractArray{<:AbstractFloat, 5},
}

    # Scaled domain.
    lx::A
    ly::A
    lz::A

    # Grid spacings.
    dx::B
    dy::B
    dz::B

    # Coordinates.
    x::A
    y::A
    z::A

    # Topography.
    topography_surface::C
    topography_spectrum::D
    k_spectrum::D
    l_spectrum::D

    # Jacobian and metric tensor.
    jac::E
    met::F

    # Vertical layers.
    ztfc::E
    ztildetfc::E
end

"""
```julia
Grid(namelists::Namelists, constants::Constants, domain::Domain) -> Grid
```

Construct a `Grid` instance for mesh generation with terrain-following coordinates.

This constructor creates a 3D computational grid with terrain-following coordinates (TFC),
including topography handling, coordinate transformations, and metric tensor computations.

# Arguments

  - `namelists`: Configuration object containing grid and domain parameters
  - `constants`: Physical constants including reference length scale
  - `domain`: MPI domain decomposition information

# Returns

  - `::Grid`: `Grid` instance.

# Additional Information

## Grid Generation Process

 1. Non-dimensionalizes domain boundaries using reference length scale
 2. Computes uniform grid spacings in each direction
 3. Generates coordinate arrays for the local MPI subdomain
 4. Applies vertical stretching transformation based on `stretch_exponent`
 5. Computes topography using specified test case
 6. Calculates Jacobian and metric tensor for terrain-following coordinates
 7. Transforms computational coordinates to physical terrain-following coordinates
"""
function Grid(namelists::Namelists, constants::Constants, domain::Domain)
    (; sizex, sizey, sizez, lx_dim, ly_dim, lz_dim, nbz) = namelists.domain
    (; testcase) = namelists.setting
    (; stretch_exponent,) = namelists.grid
    (; nxx, nyy, nzz, sizexx, sizeyy, sizezz, ko, i0, i1, j0, j1, k0) = domain
    (; lref) = constants

    # Non-dimensionalize domain boundaries.
    lx = [lx_dim...] ./ lref
    ly = [ly_dim...] ./ lref
    lz = [lz_dim...] ./ lref

    # Compute grid spacings.
    dx = (lx[2] - lx[1]) / sizex
    dy = (ly[2] - ly[1]) / sizey
    dz = (lz[2] - lz[1]) / sizez

    # Compute x-coordinate.
    x = zeros(sizexx)
    for ix in 1:sizexx
        x[ix] = lx[1] + (ix - i0) * dx + dx / 2
    end

    # Compute y-coordinate.
    y = zeros(sizeyy)
    for jy in 1:sizeyy
        y[jy] = ly[1] + (jy - j0) * dy + dy / 2
    end

    # Compute z-coordinate.
    z = zeros(sizezz)
    for kz in 1:sizezz
        z[kz] = lz[1] + (kz - k0) * dz + dz / 2
    end

    # Initialize the stretched vertical grid.
    (ztildes, zs) = (zeros(sizezz) for i in 1:2)

    # Compute the stretched vertical grid.
    for kz in 1:sizezz
        level = z[kz] + 0.5 * dz
        if level < 0
            ztildes[kz] = -lz[2] * (-level / lz[2])^stretch_exponent
        elseif level > lz[2]
            ztildes[kz] =
                2 * lz[2] -
                lz[2] * ((2 * lz[2] - level) / lz[2])^stretch_exponent
        else
            ztildes[kz] = lz[2] * (level / lz[2])^stretch_exponent
        end
    end
    for kz in 2:sizezz
        zs[kz] = 0.5 * (ztildes[kz] + ztildes[kz - 1])
    end
    zs[1] = ztildes[1] - 0.5 * (ztildes[2 * nbz] - ztildes[2 * nbz - 1])

    if lz[1] != 0.0
        error("Error in Grid: lz[1] must be zero for TFC!")
    end

    # Compute the topography.
    (topography_surface, topography_spectrum, k_spectrum, l_spectrum) =
        compute_topography(namelists, constants, domain, x, y, testcase)

    # Initialize Jacobian and metric tensor.
    jac = zeros(nxx, nyy, nzz)
    (met13, met23, met33) = (zeros(nxx, nyy, nzz) for i in 1:3)
    met = zeros(nxx, nyy, nzz, 3, 3)

    # Set the start index for the computation of the Jacobian and metric tensor.
    kz0 = ko == 0 ? 2 : 1

    # Compute the Jacobian.
    for kz in kz0:nzz
        jac[:, :, kz] .=
            (lz[2] .- topography_surface) ./ lz[2] .*
            (ztildes[ko + kz] .- ztildes[ko + kz - 1]) ./ dz
    end
    ko == 0 && @views jac[:, :, 1] .= jac[:, :, 2 * nbz]

    # Compute the metric tensor.

    met[:, :, :, 1, 2] .= 0.0
    met[:, :, :, 2, 1] .= 0.0
    met[:, :, :, 1, 1] .= 1.0
    met[:, :, :, 2, 2] .= 1.0

    for kz in kz0:nzz, jy in 1:nyy, ix in i0:i1
        met13[ix, jy, kz] =
            (topography_surface[ix + 1, jy] - topography_surface[ix - 1, jy]) /
            (2.0 * dx) * (zs[ko + kz] - lz[2]) /
            (lz[2] - topography_surface[ix, jy]) * dz /
            (ztildes[ko + kz] - ztildes[ko + kz - 1])
    end
    set_zonal_boundaries_of_field!(met13, namelists, domain)
    ko == 0 && @views met13[:, :, 1] .=
        met13[:, :, 2 * nbz] .* (zs[1] .- lz[2]) ./ (zs[2 * nbz] .- lz[2])
    met[:, :, :, 1, 3] .= met13
    met[:, :, :, 3, 1] .= met13

    for kz in 2:nzz, jy in j0:j1, ix in 1:nxx
        met23[ix, jy, kz] =
            (topography_surface[ix, jy + 1] - topography_surface[ix, jy - 1]) /
            (2.0 * dy) * (zs[ko + kz] - lz[2]) /
            (lz[2] - topography_surface[ix, jy]) * dz /
            (ztildes[ko + kz] - ztildes[ko + kz - 1])
    end
    set_meridional_boundaries_of_field!(met23, namelists, domain)
    ko == 0 && @views met23[:, :, 1] .=
        met23[:, :, 2 * nbz] .* (zs[1] .- lz[2]) ./ (zs[2 * nbz] .- lz[2])
    met[:, :, :, 2, 3] .= met23
    met[:, :, :, 3, 2] .= met23

    for kz in kz0:nzz, jy in j0:j1, ix in i0:i1
        met33[ix, jy, kz] =
            (
                (lz[2] / (lz[2] - topography_surface[ix, jy]))^2.0 +
                (
                    (zs[ko + kz] - lz[2]) /
                    (lz[2] - topography_surface[ix, jy])
                )^2.0 * (
                    (
                        (
                            topography_surface[ix + 1, jy] -
                            topography_surface[ix - 1, jy]
                        ) / (2.0 * dx)
                    )^2.0 +
                    (
                        (
                            topography_surface[ix, jy + 1] -
                            topography_surface[ix, jy - 1]
                        ) / (2.0 * dy)
                    )^2.0
                )
            ) * (dz / (ztildes[ko + kz] - ztildes[ko + kz - 1]))^2.0
    end
    ko == 0 && for jy in j0:j1, ix in i0:i1
        met33[ix, jy, 1] =
            (
                (lz[2] / (lz[2] - topography_surface[ix, jy]))^2.0 +
                ((zs[1] - lz[2]) / (lz[2] - topography_surface[ix, jy]))^2.0 * (
                    (
                        (
                            topography_surface[ix + 1, jy] -
                            topography_surface[ix - 1, jy]
                        ) / (2.0 * dx)
                    )^2.0 +
                    (
                        (
                            topography_surface[ix, jy + 1] -
                            topography_surface[ix, jy - 1]
                        ) / (2.0 * dy)
                    )^2.0
                )
            ) * (dz / (ztildes[2 * nbz] - ztildes[2 * nbz - 1]))^2.0
    end
    set_zonal_boundaries_of_field!(met33, namelists, domain)
    set_meridional_boundaries_of_field!(met33, namelists, domain)
    met[:, :, :, 3, 3] .= met33

    # Initialize the physical layers.
    (ztildetfc, ztfc) = (zeros(nxx, nyy, nzz) for i in 1:2)

    # Compute the physical layers.
    for kz in 1:nzz
        ztildetfc[:, :, kz] .=
            (lz[2] .- topography_surface) ./ lz[2] .* ztildes[ko + kz] .+
            topography_surface
        ztfc[:, :, kz] .=
            (lz[2] .- topography_surface) ./ lz[2] .* zs[ko + kz] .+
            topography_surface
    end

    # Return a Grid instance.
    return Grid(
        lx,
        ly,
        lz,
        dx,
        dy,
        dz,
        x,
        y,
        z,
        topography_surface,
        topography_spectrum,
        k_spectrum,
        l_spectrum,
        jac,
        met,
        ztfc,
        ztildetfc,
    )
end
