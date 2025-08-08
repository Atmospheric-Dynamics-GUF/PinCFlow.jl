"""
```julia
IcePredictands{A <: AbstractArray{<:AbstractFloat, 3}}
```

Arrays for prognostic ice variables.

```julia
IcePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)
```

Construct an `IcePredictands` instance with dimensions and initial values depending on the general configuration of ice physics, by dispatching to the appropriate method.

```julia
IcePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    icesetup::NoIce,
    variables::Variables,
)
```

Construct an `IcePredictands` instance with zero-size arrays for configurations without ice physics.

```julia
IcePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    icesetup::AbstractIce,
    variables::Variables,
)
```

Construct an `IcePredictands` instance with all arrays initialized as ``z \\rho`` (non-dimensionalized).

# Fields

- `n::A`: Ice-crystal number concentration.
- `q::A`: Ice mixing ratio.
- `qv::A`: Water-vapor mixing ratio.

# Arguments

- `namelists`: Namelists with all model parameters.
- `constants`: Physical constants and reference values.
- `domain`: Collection of domain-decomposition and MPI-communication parameters.
- `atmosphere`: Atmospheric-background fields.
- `grid`: Collection of parameters and fields describing the grid.
- `icesetup`: General ice-physics configuration.
- `variables`: Container for arrays needed for the prediction of the prognostic variables.
"""
struct IcePredictands{A <: AbstractArray{<:AbstractFloat, 3}}
    n::A
    q::A
    qv::A
end

function IcePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)
    (; icesetup) = namelists.ice

    return IcePredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        icesetup,
        variables,
    )
end

function IcePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    icesetup::NoIce,
    variables::Variables,
)
    n = zeros(0, 0, 0)
    q = zeros(0, 0, 0)
    qv = zeros(0, 0, 0)

    return IcePredictands(n, q, qv)
end

function IcePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    icesetup::AbstractIce,
    variables::Variables,
)
    (; nxx, nyy, nzz) = domain
    (; ztfc) = grid
    (; rhostrattfc) = atmosphere
    (; rho) = variables.predictands

    n = zeros(nxx, nyy, nzz)
    q = zeros(nxx, nyy, nzz)
    qv = zeros(nxx, nyy, nzz)

    n .= ztfc .* (rho .+ rhostrattfc)
    q .= ztfc .* (rho .+ rhostrattfc)
    qv .= ztfc .* (rho .+ rhostrattfc)

    return IcePredictands(n, q, qv)
end
