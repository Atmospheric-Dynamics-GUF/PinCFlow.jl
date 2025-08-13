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

Collection of parameters and fields that describe the grid.

```julia
Grid(namelists::Namelists, constants::Constants, domain::Domain)
```

Construct a `Grid` instance, using the specifications in `namelists.grid` and the MPI decomposition described by `domain`.

This constructor creates a 3D parallelized grid for a terrain-following, vertically stretched coordinate system. The global computational grid is defined by

```math
\\begin{align*}
    \\widehat{x} & = L_x^{\\left(0\\right)} + \\left(i - i_0 + \\frac{1}{2}\\right) \\Delta \\widehat{x},\\\\
    \\widehat{y} & = L_y^{\\left(0\\right)} + \\left(j - j_0 + \\frac{1}{2}\\right) \\Delta \\widehat{y},\\\\
    \\widehat{z} & = L_z^{\\left(0\\right)} + \\left(k - k_0 + \\frac{1}{2}\\right) \\Delta \\widehat{z},\\\\
\\end{align*}
```

where ``\\left(L_x^{\\left(0\\right)}, L_y^{\\left(0\\right)}, L_z^{\\left(0\\right)}\\right)``, ``\\left(i_0, j_0, k_0\\right)`` and ``\\left(\\Delta \\widehat{x}, \\Delta \\widehat{y}, \\Delta \\widehat{z}\\right)`` are the lower bounds of the domain, the lower index bounds of the MPI subdomains and the grid spacings (determined from the total extents and grid-point counts of the domain), respectively. The vertical layer centers and edges of the stretched and physical grids are given by

```math
\\begin{align*}
    \\widetilde{z}_{k + 1 / 2} & = L_z \\left(\\frac{\\widehat{z}_{k + 1 / 2}}{L_z}\\right)^s, & z_{k + 1 / 2} & = \\frac{L_z - h_\\mathrm{b}}{L_z} \\widetilde{z}_{k + 1 / 2} + h_\\mathrm{b},\\\\
    \\widetilde{z} & = \\frac{\\widetilde{z}_{k + 1 / 2} + \\widetilde{z}_{k - 1 / 2}}{2}, & z & = \\frac{L_z - h_\\mathrm{b}}{L_z} \\widetilde{z} + h_\\mathrm{b},
\\end{align*}
```

where ``L_z``, ``s`` and ``h_\\mathrm{b}`` are the vertical extent of the domain (`diff(namelists.domain.lz)`), the vertical-stretching parameter (`namelists.grid.stretch_exponent`) and the resolved surface topography (as returned by `compute_topography`), respectively. Finally, the Jacobian is

```math
J = \\frac{L_z - h_\\mathrm{b}}{L_z} \\frac{\\widetilde{z}_{k + 1 / 2} - \\widetilde{z}_{k - 1 / 2}}{\\Delta \\widehat{z}}
```

and the non-Cartesian elements of the metric tensor are

```math
\\begin{align*}
    G^{1 3} & = \\frac{h_{\\mathrm{b}, i + 1} - h_{\\mathrm{b}, i - 1}}{2 \\Delta \\widehat{x}} \\frac{\\widetilde{z} - L_z}{L_z - h_\\mathrm{b}} \\frac{\\Delta \\widehat{z}}{\\widetilde{z}_{k + 1 / 2} - \\widetilde{z}_{k - 1 / 2}},\\\\
    G^{2 3} & = \\frac{h_{\\mathrm{b}, j + 1} - h_{\\mathrm{b}, j - 1}}{2 \\Delta \\widehat{y}} \\frac{\\widetilde{z} - L_z}{L_z - h_\\mathrm{b}} \\frac{\\Delta \\widehat{z}}{\\widetilde{z}_{k + 1 / 2} - \\widetilde{z}_{k - 1 / 2}},\\\\
    G^{3 3} & = \\left\\{\\left(\\frac{L_z}{L_z - h_\\mathrm{b}}\\right)^2 + \\left(\\frac{\\widetilde{z} - L_z}{L_z - h_\\mathrm{b}}\\right)^2 \\left[\\left(\\frac{h_{\\mathrm{b}, i + 1} - h_{\\mathrm{b}, i - 1}}{2 \\Delta \\widehat{x}}\\right)^2 + \\left(\\frac{h_{\\mathrm{b}, j + 1} - h_{\\mathrm{b}, j - 1}}{2 \\Delta \\widehat{y}}\\right)^2\\right]\\right\\} \\left(\\frac{\\Delta \\widehat{z}}{\\widetilde{z}_{k + 1 / 2} - \\widetilde{z}_{k - 1 / 2}}\\right)^2.
\\end{align*}
```

# Fields

Domain boundaries:

  - `lx::A`: Non-dimensional domain boundaries in ``\\widehat{x}``-direction.

  - `ly::A`: Non-dimensional domain boundaries in ``\\widehat{y}``-direction.

  - `lz::A`: Non-dimensional domain boundaries in ``\\widehat{z}``-direction.

Grid spacing:

  - `dx::B`: Grid spacing ``\\Delta \\widehat{x}``.

  - `dy::B`: Grid spacing ``\\Delta \\widehat{y}``.

  - `dz::B`: Grid spacing ``\\Delta \\widehat{z}``.

Coordinate arrays:

  - `x::A`: Cell-centered ``\\widehat{x}``-coordinate of the entire domain.

  - `y::A`: Cell-centered ``\\widehat{y}``-coordinate of the entire domain.

  - `z::A`: Cell-centered ``\\widehat{z}``-coordinate of the entire domain.

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
