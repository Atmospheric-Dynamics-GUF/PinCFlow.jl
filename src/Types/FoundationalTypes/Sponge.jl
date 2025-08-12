"""
```julia
Sponge{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractFloat,
    C <: AbstractVector{<:AbstractFloat},
}
```

Composite type for Rayleigh-damping coefficients, sponge-layer bounds and extents, as well as an auxiliary array for the computation of horizontal means.

```julia
Sponge(namelists::Namelists, domain::Domain, grid::Grid)
```

Construct a `Sponge` instance, using the model parameters in `namelists`.

The vertical extent of the sponge is set to the fraction `namelists.sponge.spongeheight` of the vertical extent of the domain. The horizontal extents of the unified sponge are computed similarly, using the same parameter multiplied by `0.5` (since the sponge is centered at the horizontal boundaries).

# Fields

Rayleigh-damping coefficients:

  - `alphaunifiedsponge::A`: Coefficient of the unified sponge (used in all prognostic equations).
  - `kr_sp_tfc::A`: Coefficient of the non-unified sponge (used in the horizontal-momentum equation).
  - `kr_sp_w_tfc::A`: Coefficient of the non-unified sponge (used in the transformed-vertical-momentum equation).

Vertical sponge extent:

  - `zsponge::B`: Lower edge of the sponge.
  - `dzsponge::B`: Vertical extent of the sponge.

Horizontal sponge extent:

  - `xsponge0::B`: Right edge of the unified sponge.
  - `xsponge1::B`: Left edge of the unified sponge.
  - `ysponge0::B`: Forward edge of the unified sponge.
  - `ysponge1::B`: Backward edge of the unified sponge.
  - `dxsponge::B`: Halved zonal extent of the unified sponge.
  - `dysponge::B`: Halved meridional extent of the unified sponge.

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
    alphaunifiedsponge::A
    kr_sp_tfc::A
    kr_sp_w_tfc::A

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

function Sponge(namelists::Namelists, domain::Domain, grid::Grid)
    (; spongeheight) = namelists.sponge
    (; nxx, nyy, nzz, nz) = domain
    (; lx, ly, lz) = grid

    # Initialize the sponge layer coefficients.
    (kr_sp_tfc, kr_sp_w_tfc, alphaunifiedsponge) =
        (zeros(nxx, nyy, nzz) for i in 1:3)

    # Set up the sponge layer.
    dzsponge = spongeheight * (lz[2] - lz[1])
    zsponge = lz[2] - dzsponge
    dxsponge = 0.5 * spongeheight * (lx[2] - lx[1])
    dysponge = 0.5 * spongeheight * (ly[2] - ly[1])
    xsponge0 = lx[1] + dxsponge
    ysponge0 = ly[1] + dysponge
    xsponge1 = lx[2] - dxsponge
    ysponge1 = ly[2] - dysponge

    # Initialize the auxiliary array for horizontal means.
    horizontal_mean = zeros(nz)

    # Return a Sponge instance.
    return Sponge(
        alphaunifiedsponge,
        kr_sp_tfc,
        kr_sp_w_tfc,
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
