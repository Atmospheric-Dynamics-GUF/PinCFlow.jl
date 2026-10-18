"""
```julia
Atmosphere{A <: AbstractArray{<:AbstractFloat, 3}}
```

Composite type for atmospheric background fields.

```julia
Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
)::Atmosphere
```

Create an `Atmosphere` instance, with the background fields computed via `compute_background_fields!`.

# Fields

  - `pbar::A`: Mass-weighted potential temperature ``P \\left(z\\right)`` (``P \\left(x, y, z, t\\right)`` in compressible mode).

  - `thetabar::A`: Background potential temperature ``\\bar{\\theta} \\left(z\\right)``.

  - `rhobar::A`: Background density ``\\bar{\\rho} \\left(z\\right)``.

  - `n2::A`: Squared buoyancy frequency ``N^2 \\left(z\\right)``.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `grid`: Collection of parameters and fields that describe the grid.

# See also

  - [`PinCFlow.Types.FoundationalTypes.compute_background_fields!`](@ref)
"""
struct Atmosphere{A <: AbstractArray{<:AbstractFloat, 3}}
    pbar::A
    thetabar::A
    rhobar::A
    n2::A
end

function Atmosphere(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    grid::Grid,
)::Atmosphere
    (; nxx, nyy, nzz) = domain

    (pbar, thetabar, rhobar, n2) = (zeros(nxx, nyy, nzz) for i in 1:4)

    atmosphere = Atmosphere(pbar, thetabar, rhobar, n2)

    compute_background_fields!(namelists, constants, domain, grid, atmosphere)

    return atmosphere
end
