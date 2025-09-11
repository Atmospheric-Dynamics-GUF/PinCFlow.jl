"""
```julia
Grid{
    A <: AbstractFloat,
    B <: AbstractVector{<:AbstractFloat},
    C <: AbstractMatrix{<:AbstractFloat},
    D <: AbstractArray{<:AbstractFloat, 3},
    E <: AbstractArray{<:AbstractFloat, 3},
    F <: AbstractArray{<:AbstractFloat, 5},
}
```

Collection of parameters and fields that describe the grid.

```julia
Grid(namelists::Namelists, constants::Constants, domain::Domain)::Grid
```

Construct a `Grid` instance, using the specifications in `namelists.grid` and the MPI decomposition described by `domain`.

This constructor creates a 3D parallelized grid for a terrain-following, vertically stretched coordinate system. The global computational grid is defined by

```math
\\begin{align*}
    \\widehat{x}_i & = - L_x / 2 + \\left(i - i_0 + \\frac{1}{2}\\right) \\Delta \\widehat{x}, & i & \\in \\left\\{1, ..., N_{\\mathrm{t}, \\widehat{x}}\\right\\},\\\\
    \\widehat{y}_j & = - L_y / 2 + \\left(j - j_0 + \\frac{1}{2}\\right) \\Delta \\widehat{y}, & j & \\in \\left\\{1, ..., N_{\\mathrm{t}, \\widehat{y}}\\right\\},\\\\
    \\widehat{z}_k & = \\left(k - k_0 + \\frac{1}{2}\\right) \\Delta \\widehat{z}, & k & \\in \\left\\{1, ..., N_{\\mathrm{t}, \\widehat{z}}\\right\\},\\\\
\\end{align*}
```

where ``\\left(L_x, L_y\\right)``, ``\\left(i_0, j_0, k_0\\right)`` and ``\\left(\\Delta \\widehat{x}, \\Delta \\widehat{y}, \\Delta \\widehat{z}\\right)`` are the horizontal extents of the domain, the lower index bounds of the MPI subdomains and the grid spacings (determined from the total extents and grid-point counts of the domain), respectively. The total number of grid points in each dimension (``N_{\\mathrm{t}, \\widehat{x}}``, ``N_{\\mathrm{t}, \\widehat{y}}`` or ``N_{\\mathrm{t}, \\widehat{z}}``) is defined as the respective number of physical grid points (``N_{\\widehat{x}}``, ``N_{\\widehat{y}}`` or ``N_{\\widehat{z}}``) plus two times the number of ghost cells beyond each boundary (``N_{\\mathrm{b}, \\widehat{x}}``, ``N_{\\mathrm{b}, \\widehat{y}}`` or ``N_{\\mathrm{b}, \\widehat{z}}``). Throughout the documentation, the position of any variable on this grid is indicated with the indices ``\\left(i, j, k\\right)`` in its subscript. Therein, unshifted indices are omitted for the sake of brevity. The grid is staggered, i.e. the wind components are defined at the midpoints of those cell surfaces that are orthogonal to their respective directions. Interpolations are therefore necessary in many places. These are indicated as in

```math
\\begin{align*}
    \\rho_{i + 1 / 2} & = \\frac{\\rho + \\rho_{i + 1}}{2}, & \\rho_{j + 1 / 2} & = \\frac{\\rho + \\rho_{j + 1}}{2}, & \\rho_{k + 1 / 2} & = \\frac{J_{k + 1} \\rho + J \\rho_{k + 1}}{J_{k + 1} + J},\\\\
    u & = \\frac{u_{i - 1 / 2} + u_{i + 1 / 2}}{2}, & u_{j + 1 / 2} & = \\frac{u + u_{j + 1}}{2}, & u_{k + 1 / 2} & = \\frac{J_{k + 1} u + J u_{k + 1}}{J_{k + 1} + J},\\\\
    v_{i + 1 / 2} & = \\frac{v + v_{i + 1}}{2}, & v & = \\frac{v_{j - 1 / 2} + v_{j + 1 / 2}}{2}, & v_{k + 1 / 2} & = \\frac{J_{k + 1} v + J v_{k + 1}}{J_{k + 1} + J},\\\\
    \\widehat{w}_{i + 1 / 2} & = \\frac{\\widehat{w} + \\widehat{w}_{i + 1}}{2}, & \\widehat{w}_{j + 1 / 2} & = \\frac{\\widehat{w} + \\widehat{w}_{j + 1}}{2}, & \\widehat{w} & = \\frac{\\widehat{w}_{k - 1 / 2} + \\widehat{w}_{k + 1 / 2}}{2}.
\\end{align*}
```

The vertical layer centers and edges of the stretched and physical grids are given by

```math
\\begin{align*}
    \\widetilde{z}_{k + 1 / 2} & = L_z \\left(\\frac{\\widehat{z}_{k + 1 / 2}}{L_z}\\right)^s, & z_{k + 1 / 2} & = \\frac{L_z - h}{L_z} \\widetilde{z}_{k + 1 / 2} + h,\\\\
    \\widetilde{z} & = \\frac{\\widetilde{z}_{k + 1 / 2} + \\widetilde{z}_{k - 1 / 2}}{2}, & z & = \\frac{L_z - h}{L_z} \\widetilde{z} + h,
\\end{align*}
```

where ``L_z``, ``s`` and ``h`` are the vertical extent of the domain (`namelists.domain.lz_dim`), the vertical-stretching parameter (`namelists.grid.stretch_exponent`) and the surface topography (as returned by `compute_topography`), respectively. Finally, the Jacobian is

```math
J = \\frac{L_z - h}{L_z} \\frac{\\widetilde{z}_{k + 1 / 2} - \\widetilde{z}_{k - 1 / 2}}{\\Delta \\widehat{z}}
```

and the non-Cartesian elements of the metric tensor are

```math
\\begin{align*}
    G^{1 3} & = \\frac{h_{\\mathrm{b}, i + 1} - h_{\\mathrm{b}, i - 1}}{2 \\Delta \\widehat{x}} \\frac{\\widetilde{z} - L_z}{L_z - h} \\frac{\\Delta \\widehat{z}}{\\widetilde{z}_{k + 1 / 2} - \\widetilde{z}_{k - 1 / 2}},\\\\
    G^{2 3} & = \\frac{h_{\\mathrm{b}, j + 1} - h_{\\mathrm{b}, j - 1}}{2 \\Delta \\widehat{y}} \\frac{\\widetilde{z} - L_z}{L_z - h} \\frac{\\Delta \\widehat{z}}{\\widetilde{z}_{k + 1 / 2} - \\widetilde{z}_{k - 1 / 2}},\\\\
    G^{3 3} & = \\left\\{\\left(\\frac{L_z}{L_z - h}\\right)^2 + \\left(\\frac{\\widetilde{z} - L_z}{L_z - h}\\right)^2 \\left[\\left(\\frac{h_{\\mathrm{b}, i + 1} - h_{\\mathrm{b}, i - 1}}{2 \\Delta \\widehat{x}}\\right)^2 + \\left(\\frac{h_{\\mathrm{b}, j + 1} - h_{\\mathrm{b}, j - 1}}{2 \\Delta \\widehat{y}}\\right)^2\\right]\\right\\}\\\\
    & \\quad \\times \\left(\\frac{\\Delta \\widehat{z}}{\\widetilde{z}_{k + 1 / 2} - \\widetilde{z}_{k - 1 / 2}}\\right)^2.
\\end{align*}
```

# Fields

Domain extent:

  - `lx::A`: Non-dimensional domain extent in ``\\widehat{x}``-direction.

  - `ly::A`: Non-dimensional domain extent in ``\\widehat{y}``-direction.

  - `lz::A`: Non-dimensional domain extent in ``\\widehat{z}``-direction.

Grid spacing:

  - `dx::A`: Grid spacing ``\\Delta \\widehat{x}``.

  - `dy::A`: Grid spacing ``\\Delta \\widehat{y}``.

  - `dz::A`: Grid spacing ``\\Delta \\widehat{z}``.

Coordinate arrays:

  - `x::B`: Cell-centered ``\\widehat{x}``-coordinate of the entire domain.

  - `y::B`: Cell-centered ``\\widehat{y}``-coordinate of the entire domain.

  - `z::B`: Cell-centered ``\\widehat{z}``-coordinate of the entire domain.

Topography:

  - `topography_surface::C`: Resolved surface topography.

  - `topography_spectrum::D`: Spectrum of the unresolved surface topography.

  - `k_spectrum::D`: Zonal wavenumbers of the spectrum.

  - `l_spectrum::D`: Meridional wavenumbers of the spectrum.

Coordinate transformation.

  - `jac::E`: Jacobian.

  - `met::F`: Metric tensor.

Physical coordinates:

  - `ztfc::E`: Physical height at cell centers.

  - `ztildetfc::E`: Physical height at vertical cell edges.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

# See also

  - [`PinCFlow.Types.FoundationalTypes.compute_topography`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.Types.FoundationalTypes.set_meridional_boundaries_of_field!`](@ref)
"""
struct Grid{
    A <: AbstractFloat,
    B <: AbstractVector{<:AbstractFloat},
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
    dx::A
    dy::A
    dz::A

    # Coordinates.
    x::B
    y::B
    z::B

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

function Grid(namelists::Namelists, constants::Constants, domain::Domain)::Grid
    (; sizex, sizey, sizez, lx_dim, ly_dim, lz_dim, nbz) = namelists.domain
    (; testcase) = namelists.setting
    (; stretch_exponent,) = namelists.grid
    (; nxx, nyy, nzz, sizexx, sizeyy, sizezz, ko, i0, i1, j0, j1, k0) = domain
    (; lref) = constants

    # Non-dimensionalize domain boundaries.
    lx = lx_dim / lref
    ly = ly_dim / lref
    lz = lz_dim / lref

    # Compute grid spacings.
    dx = lx / sizex
    dy = ly / sizey
    dz = lz / sizez

    # Compute x-coordinate.
    x = zeros(sizexx)
    @ivy for ix in 1:sizexx
        x[ix] = -lx / 2 + (ix - i0) * dx + dx / 2
    end

    # Compute y-coordinate.
    y = zeros(sizeyy)
    @ivy for jy in 1:sizeyy
        y[jy] = -ly / 2 + (jy - j0) * dy + dy / 2
    end

    # Compute z-coordinate.
    z = zeros(sizezz)
    @ivy for kz in 1:sizezz
        z[kz] = (kz - k0) * dz + dz / 2
    end

    # Initialize the stretched vertical grid.
    (ztildes, zs) = (zeros(sizezz) for i in 1:2)

    # Compute the stretched vertical grid.
    @ivy for kz in 1:sizezz
        level = z[kz] + 0.5 * dz
        if level < 0
            ztildes[kz] = -lz * (-level / lz)^stretch_exponent
        elseif level > lz
            ztildes[kz] = 2 * lz - lz * ((2 * lz - level) / lz)^stretch_exponent
        else
            ztildes[kz] = lz * (level / lz)^stretch_exponent
        end
    end
    @ivy for kz in 2:sizezz
        zs[kz] = 0.5 * (ztildes[kz] + ztildes[kz - 1])
    end
    @ivy zs[1] = ztildes[1] - 0.5 * (ztildes[2 * nbz] - ztildes[2 * nbz - 1])

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
    @ivy for kz in kz0:nzz
        jac[:, :, kz] .=
            (lz .- topography_surface) ./ lz .*
            (ztildes[ko + kz] .- ztildes[ko + kz - 1]) ./ dz
    end
    @ivy ko == 0 && jac[:, :, 1] .= jac[:, :, 2 * nbz]

    # Compute the metric tensor.

    @ivy met[:, :, :, 1, 2] .= 0.0
    @ivy met[:, :, :, 2, 1] .= 0.0
    @ivy met[:, :, :, 1, 1] .= 1.0
    @ivy met[:, :, :, 2, 2] .= 1.0

    @ivy for kz in kz0:nzz, jy in 1:nyy, ix in i0:i1
        met13[ix, jy, kz] =
            (topography_surface[ix + 1, jy] - topography_surface[ix - 1, jy]) /
            (2.0 * dx) * (zs[ko + kz] - lz) /
            (lz - topography_surface[ix, jy]) * dz /
            (ztildes[ko + kz] - ztildes[ko + kz - 1])
    end
    set_zonal_boundaries_of_field!(met13, namelists, domain)
    @ivy ko == 0 && met13[:, :, 1] .=
        met13[:, :, 2 * nbz] .* (zs[1] .- lz) ./ (zs[2 * nbz] .- lz)
    @ivy met[:, :, :, 1, 3] .= met13
    @ivy met[:, :, :, 3, 1] .= met13

    @ivy for kz in 2:nzz, jy in j0:j1, ix in 1:nxx
        met23[ix, jy, kz] =
            (topography_surface[ix, jy + 1] - topography_surface[ix, jy - 1]) /
            (2.0 * dy) * (zs[ko + kz] - lz) /
            (lz - topography_surface[ix, jy]) * dz /
            (ztildes[ko + kz] - ztildes[ko + kz - 1])
    end
    set_meridional_boundaries_of_field!(met23, namelists, domain)
    @ivy ko == 0 && met23[:, :, 1] .=
        met23[:, :, 2 * nbz] .* (zs[1] .- lz) ./ (zs[2 * nbz] .- lz)
    @ivy met[:, :, :, 2, 3] .= met23
    @ivy met[:, :, :, 3, 2] .= met23

    @ivy for kz in kz0:nzz, jy in j0:j1, ix in i0:i1
        met33[ix, jy, kz] =
            (
                (lz / (lz - topography_surface[ix, jy]))^2.0 +
                ((zs[ko + kz] - lz) / (lz - topography_surface[ix, jy]))^2.0 *
                (
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
    @ivy ko == 0 && for jy in j0:j1, ix in i0:i1
        met33[ix, jy, 1] =
            (
                (lz / (lz - topography_surface[ix, jy]))^2.0 +
                ((zs[1] - lz) / (lz - topography_surface[ix, jy]))^2.0 * (
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
    @ivy met[:, :, :, 3, 3] .= met33

    # Initialize the physical layers.
    (ztildetfc, ztfc) = (zeros(nxx, nyy, nzz) for i in 1:2)

    # Compute the physical layers.
    @ivy for kz in 1:nzz
        ztildetfc[:, :, kz] .=
            (lz .- topography_surface) ./ lz .* ztildes[ko + kz] .+
            topography_surface
        ztfc[:, :, kz] .=
            (lz .- topography_surface) ./ lz .* zs[ko + kz] .+
            topography_surface
    end

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
