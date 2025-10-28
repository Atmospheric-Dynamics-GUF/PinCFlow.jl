"""
```julia
Sponge{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractFloat,
    D <: AbstractFloat,
    E <: AbstractFloat,
    F <: AbstractFloat,
    G <: AbstractFloat,
    H <: AbstractFloat,
    I <: AbstractFloat,
    J <: AbstractFloat,
    K <: AbstractVector{<:AbstractFloat},
}
```

Composite type for Rayleigh-damping coefficients, sponge bounds and extents, as well as an auxiliary array for the computation of horizontal means.

```julia
Sponge(namelists::Namelists, domain::Domain, grid::Grid)::Sponge
```

Construct a `Sponge` instance, using the model parameters in `namelists`.

The vertical extent of the sponge is set to the fraction `namelists.sponge.sponge_extent` of the vertical extent of the domain. The horizontal extents of the LHS sponge are computed similarly, using the same parameter multiplied by `0.5` (since the sponge is centered at the horizontal boundaries).

# Fields

Rayleigh-damping coefficients:

  - `alphar::A`: Coefficient of the LHS sponge (used in all prognostic equations).

  - `betar::B`: Coefficient of the RHS sponge (used in the momentum equation).

Vertical sponge extent:

  - `zsponge::C`: Lower edge of both sponges.

  - `dzsponge::D`: Vertical extent of both sponges.

Horizontal sponge extent:

  - `xsponge0::E`: Right edge of the LHS sponge.

  - `xsponge1::F`: Left edge of the LHS sponge.

  - `ysponge0::G`: Forward edge of the LHS sponge.

  - `ysponge1::H`: Backward edge of the LHS sponge.

  - `dxsponge::I`: Halved zonal extent of the LHS sponge.

  - `dysponge::J`: Halved meridional extent of the LHS sponge.

Auxiliary array:

  - `horizontal_mean::K`: Auxiliary array for the computation of horizontal means.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `grid`: Collection of parameters and fields that describe the grid.
"""
struct Sponge{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
    C <: AbstractFloat,
    D <: AbstractFloat,
    E <: AbstractFloat,
    F <: AbstractFloat,
    G <: AbstractFloat,
    H <: AbstractFloat,
    I <: AbstractFloat,
    J <: AbstractFloat,
    K <: AbstractVector{<:AbstractFloat},
}

    # Damping coefficients.
    alphar::A
    betar::B

    # Vertical sponge extent.
    zsponge::C
    dzsponge::D

    # Lateral sponge extent.
    xsponge0::E
    xsponge1::F
    ysponge0::G
    ysponge1::H
    dxsponge::I
    dysponge::J

    # Auxiliary array for horizontal means.
    horizontal_mean::K
end

function Sponge(namelists::Namelists, domain::Domain, grid::Grid)::Sponge
    (; sponge_extent) = namelists.sponge
    (; nxx, nyy, nzz, nz) = domain
    (; lx, ly, lz) = grid

    (betar, alphar) = (zeros(nxx, nyy, nzz) for i in 1:3)

    dzsponge = sponge_extent * lz
    zsponge = lz - dzsponge
    dxsponge = sponge_extent * lx / 2
    dysponge = sponge_extent * ly / 2
    xsponge0 = -lx / 2 + dxsponge
    ysponge0 = -ly / 2 + dysponge
    xsponge1 = lx / 2 - dxsponge
    ysponge1 = ly / 2 - dysponge

    horizontal_mean = zeros(nz)

    return Sponge(
        alphar,
        betar,
        zsponge,
        dzsponge,
        xsponge0,
        xsponge1,
        ysponge0,
        ysponge1,
        dxsponge,
        dysponge,
        horizontal_mean,
    )
end
