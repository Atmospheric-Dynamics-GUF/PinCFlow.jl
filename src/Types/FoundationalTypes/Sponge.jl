"""
```julia
Sponge{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractFloat,
    C <: AbstractVector{<:AbstractFloat},
}
```

Composite type for Rayleigh-damping coefficients, sponge bounds and extents, as well as an auxiliary array for the computation of horizontal means.

```julia
Sponge(namelists::Namelists, domain::Domain, grid::Grid)::Sponge
```

Construct a `Sponge` instance, using the model parameters in `namelists`.

The vertical extent of the sponge is set to the fraction `namelists.sponge.spongeheight` of the vertical extent of the domain. The horizontal extents of the LHS sponge are computed similarly, using the same parameter multiplied by `0.5` (since the sponge is centered at the horizontal boundaries).

# Fields

Rayleigh-damping coefficients:

  - `alphar::A`: Coefficient of the LHS sponge (used in all prognostic equations).

  - `betar::A`: Coefficient of the RHS sponge (used in the momentum equation).

Vertical sponge extent:

  - `zsponge::B`: Lower edge of both sponges.

  - `dzsponge::B`: Vertical extent of both sponges.

Horizontal sponge extent:

  - `xsponge0::B`: Right edge of the LHS sponge.

  - `xsponge1::B`: Left edge of the LHS sponge.

  - `ysponge0::B`: Forward edge of the LHS sponge.

  - `ysponge1::B`: Backward edge of the LHS sponge.

  - `dxsponge::B`: Halved zonal extent of the LHS sponge.

  - `dysponge::B`: Halved meridional extent of the LHS sponge.

Auxiliary array:

  - `horizontal_mean::C`: Auxiliary array for the computation of horizontal means.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `grid`: Collection of parameters and fields that describe the grid.
"""
struct Sponge{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractFloat,
    C <: AbstractVector{<:AbstractFloat},
}

    # Damping coefficients.
    alphar::A
    betar::A

    # Vertical sponge extent.
    zsponge::B
    dzsponge::B

    # Lateral sponge extent.
    xsponge0::B
    xsponge1::B
    ysponge0::B
    ysponge1::B
    dxsponge::B
    dysponge::B

    # Auxiliary array for horizontal means.
    horizontal_mean::C
end

function Sponge(namelists::Namelists, domain::Domain, grid::Grid)::Sponge
    (; spongeheight) = namelists.sponge
    (; nxx, nyy, nzz, nz) = domain
    (; lx, ly, lz) = grid

    # Initialize the sponge coefficients.
    (betar, alphar) = (zeros(nxx, nyy, nzz) for i in 1:3)

    # Set up the sponges.
    dzsponge = spongeheight * lz
    zsponge = lz - dzsponge
    dxsponge = spongeheight * lx / 2
    dysponge = spongeheight * ly / 2
    xsponge0 = -lx / 2 + dxsponge
    ysponge0 = -ly / 2 + dysponge
    xsponge1 = lx / 2 - dxsponge
    ysponge1 = ly / 2 - dysponge

    # Initialize the auxiliary array for horizontal means.
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
