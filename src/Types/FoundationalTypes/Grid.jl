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
    \\widehat{x}_i & = - L_x / 2 + \\left(i - i_0 + \\frac{1}{2}\\right) \\Delta \\widehat{x},\\\\
    \\widehat{y}_j & = - L_y / 2 + \\left(j - j_0 + \\frac{1}{2}\\right) \\Delta \\widehat{y},\\\\
    \\widehat{z}_k & = \\left(k - k_0 + \\frac{1}{2}\\right) \\Delta \\widehat{z},\\\\
\\end{align*}
```

where ``\\left(L_x, L_y\\right)``, ``\\left(i_0, j_0, k_0\\right)`` and ``\\left(\\Delta \\widehat{x}, \\Delta \\widehat{y}, \\Delta \\widehat{z}\\right)`` are the horizontal extents of the domain, the lower index bounds of the MPI subdomains and the grid spacings (determined from the total extents and grid-point counts of the domain), respectively. Throughout the documentation, the position of any variable on this grid is indicated with the indices ``\\left(i, j, k\\right)`` in its subscript. Therein, unshifted indices are omitted for the sake of brevity. The grid is staggered, i.e. the wind components are defined at the midpoints of those cell surfaces that are orthogonal to their respective directions. Interpolations are therefore necessary in many places. These are indicated as in

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

where ``L_z``, ``s`` and ``h`` are the vertical extent of the domain (`namelists.domain.lz`), the vertical-stretching parameter (`namelists.grid.stretch_exponent`) and the surface topography (as returned by `compute_topography`), respectively. Finally, the Jacobian is

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

Horizontal coordinates:

  - `x::B`: Cell-centered ``\\widehat{x}``-coordinate of the entire domain.

  - `y::B`: Cell-centered ``\\widehat{y}``-coordinate of the entire domain.

Topography:

  - `hb::C`: Resolved surface topography.

  - `hw::D`: Spectrum of the unresolved surface topography.

  - `kh::D`: Zonal wavenumbers of the spectrum.

  - `lh::D`: Meridional wavenumbers of the spectrum.

Coordinate transformation.

  - `jac::E`: Jacobian.

  - `met::F`: Metric tensor.

Vertical coordinates:

  - `zc::E`: Physical height at cell centers.

  - `zctilde::E`: Physical height at vertical cell edges.

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

    # Horizontal coordinates.
    x::B
    y::B

    # Topography.
    hb::C
    hw::D
    kh::D
    lh::D

    # Jacobian and metric tensor.
    jac::E
    met::F

    # Vertical coordinates.
    zc::E
    zctilde::E
end

function Grid(namelists::Namelists, constants::Constants, domain::Domain)::Grid
    (; x_size, y_size, z_size, nbz) = namelists.domain
    (; test_case) = namelists.setting
    (; stretch_exponent,) = namelists.grid
    (; nxx, nyy, nzz, xx_size, yy_size, zz_size, ko, i0, i1, j0, j1, k0) =
        domain
    (; lref) = constants

    # Non-dimensionalize domain boundaries.
    lx = namelists.domain.lx / lref
    ly = namelists.domain.ly / lref
    lz = namelists.domain.lz / lref

    # Compute grid spacings.
    dx = lx / x_size
    dy = ly / y_size
    dz = lz / z_size

    # Compute x-coordinate.
    x = zeros(xx_size)
    @ivy for i in 1:xx_size
        x[i] = -lx / 2 + (i - i0) * dx + dx / 2
    end

    # Compute y-coordinate.
    y = zeros(yy_size)
    @ivy for j in 1:yy_size
        y[j] = -ly / 2 + (j - j0) * dy + dy / 2
    end

    # Compute z-coordinate.
    z = zeros(zz_size)
    @ivy for k in 1:zz_size
        z[k] = (k - k0) * dz + dz / 2
    end

    # Allocate the stretched vertical grid.
    (ztildes, zs) = (zeros(zz_size) for i in 1:2)

    # Compute the stretched vertical grid.
    @ivy for k in 1:zz_size
        level = z[k] + 0.5 * dz
        if level < 0
            ztildes[k] = -lz * (-level / lz)^stretch_exponent
        elseif level > lz
            ztildes[k] = 2 * lz - lz * ((2 * lz - level) / lz)^stretch_exponent
        else
            ztildes[k] = lz * (level / lz)^stretch_exponent
        end
    end
    @ivy for k in 2:zz_size
        zs[k] = 0.5 * (ztildes[k] + ztildes[k - 1])
    end
    @ivy zs[1] = ztildes[1] - 0.5 * (ztildes[2 * nbz] - ztildes[2 * nbz - 1])

    # Compute the topography.
    (hb, hw, kh, lh) =
        compute_topography(namelists, constants, domain, x, y, test_case)

    # Allocate Jacobian and metric tensor.
    jac = zeros(nxx, nyy, nzz)
    (met13, met23, met33) = (zeros(nxx, nyy, nzz) for i in 1:3)
    met = zeros(nxx, nyy, nzz, 3, 3)

    # Set the start index for the computation of the Jacobian and metric tensor.
    kmin = ko == 0 ? 2 : 1
    kmax = nzz

    # Compute the Jacobian.
    @ivy for k in kmin:kmax
        jac[:, :, k] .=
            (lz .- hb) ./ lz .* (ztildes[ko + k] .- ztildes[ko + k - 1]) ./ dz
    end
    @ivy ko == 0 && (jac[:, :, 1] .= jac[:, :, 2 * nbz])

    # Compute the metric tensor.

    @ivy met[:, :, :, 1, 2] .= 0.0
    @ivy met[:, :, :, 2, 1] .= 0.0
    @ivy met[:, :, :, 1, 1] .= 1.0
    @ivy met[:, :, :, 2, 2] .= 1.0

    @ivy for k in kmin:kmax, j in 1:nyy, i in i0:i1
        met13[i, j, k] =
            (hb[i + 1, j] - hb[i - 1, j]) / (2.0 * dx) * (zs[ko + k] - lz) /
            (lz - hb[i, j]) * dz / (ztildes[ko + k] - ztildes[ko + k - 1])
    end
    set_zonal_boundaries_of_field!(met13, namelists, domain)
    @ivy ko == 0 && (
        met13[:, :, 1] .=
            met13[:, :, 2 * nbz] .* (zs[1] .- lz) ./ (zs[2 * nbz] .- lz)
    )
    @ivy met[:, :, :, 1, 3] .= met13
    @ivy met[:, :, :, 3, 1] .= met13

    @ivy for k in 2:nzz, j in j0:j1, i in 1:nxx
        met23[i, j, k] =
            (hb[i, j + 1] - hb[i, j - 1]) / (2.0 * dy) * (zs[ko + k] - lz) /
            (lz - hb[i, j]) * dz / (ztildes[ko + k] - ztildes[ko + k - 1])
    end
    set_meridional_boundaries_of_field!(met23, namelists, domain)
    @ivy ko == 0 && (
        met23[:, :, 1] .=
            met23[:, :, 2 * nbz] .* (zs[1] .- lz) ./ (zs[2 * nbz] .- lz)
    )
    @ivy met[:, :, :, 2, 3] .= met23
    @ivy met[:, :, :, 3, 2] .= met23

    @ivy for k in kmin:kmax, j in j0:j1, i in i0:i1
        met33[i, j, k] =
            (
                (lz / (lz - hb[i, j]))^2.0 +
                ((zs[ko + k] - lz) / (lz - hb[i, j]))^2.0 * (
                    ((hb[i + 1, j] - hb[i - 1, j]) / (2.0 * dx))^2.0 +
                    ((hb[i, j + 1] - hb[i, j - 1]) / (2.0 * dy))^2.0
                )
            ) * (dz / (ztildes[ko + k] - ztildes[ko + k - 1]))^2.0
    end
    @ivy ko == 0 && for j in j0:j1, i in i0:i1
        met33[i, j, 1] =
            (
                (lz / (lz - hb[i, j]))^2.0 +
                ((zs[1] - lz) / (lz - hb[i, j]))^2.0 * (
                    ((hb[i + 1, j] - hb[i - 1, j]) / (2.0 * dx))^2.0 +
                    ((hb[i, j + 1] - hb[i, j - 1]) / (2.0 * dy))^2.0
                )
            ) * (dz / (ztildes[2 * nbz] - ztildes[2 * nbz - 1]))^2.0
    end
    set_zonal_boundaries_of_field!(met33, namelists, domain)
    set_meridional_boundaries_of_field!(met33, namelists, domain)
    @ivy met[:, :, :, 3, 3] .= met33

    # Allocate the physical layers.
    (zctilde, zc) = (zeros(nxx, nyy, nzz) for i in 1:2)

    # Compute the physical layers.
    @ivy for k in 1:nzz
        zctilde[:, :, k] .= (lz .- hb) ./ lz .* ztildes[ko + k] .+ hb
        zc[:, :, k] .= (lz .- hb) ./ lz .* zs[ko + k] .+ hb
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
        hb,
        hw,
        kh,
        lh,
        jac,
        met,
        zc,
        zctilde,
    )
end
