"""
```julia
TurbulencePredictands{A <: AbstractArray{<:AbstractFloat, 3}}
```

Arrays for prognostic turbulence variables.

```julia
TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)
```

Construct a `TurbulencePredictands` instance with dimensions and initial values depending on the general configuration of turbulence physics, by dispatching to the appropriate method.

```julia
TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulencesetup::NoTurbulence,
    variables::Variables,
)
```

Construct a `TurbulencePredictands` instance with zero-size arrays for configurations without turbulence physics.

```julia
TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulencesetup::AbstractTurbulence,
    variables::Variables,
)
```

Construct a `TurbulencePredictands` instance with both arrays initialized as ``t_\\mathrm{ref}^2 / \\left(10 L_\\mathrm{ref}^2\\right) \\rho`` (non-dimensionalized), where ``t_\\mathrm{ref}`` and ``L_\\mathrm{ref}`` are given by the properties `tref` and `lref` of `constants`, respectively.

# Fields

  - `tke::A`: Turbulent kinetic energy.
  - `tte::A`: Total turbulent energy.

# Arguments

  - `namelists`: Namelists with all model parameters.
  - `constants`: Physical constants and reference values.
  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
  - `atmosphere`: Atmospheric-background fields.
  - `grid`: Collection of parameters and fields describing the grid.
  - `turbulencesetup`: General turbulence-physics configuration.
  - `variables`: Container for arrays needed for the prediction of the prognostic variables.
"""
struct TurbulencePredictands{A <: AbstractArray{<:AbstractFloat, 3}}
    tke::A
    tte::A
end

function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)
    (; turbulencesetup) = namelists.turbulence

    return TurbulencePredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        turbulencesetup,
        variables,
    )
end

function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulencesetup::NoTurbulence,
    variables::Variables,
)
    tke = zeros(0, 0, 0)
    tte = zeros(0, 0, 0)

    return TurbulencePredictands(tke, tte)
end

function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulencesetup::AbstractTurbulence,
    variables::Variables,
)
    (; nxx, nyy, nzz) = domain
    (; ztfc) = grid
    (; rhostrattfc) = atmosphere
    (; rho) = variables.predictands
    (; lref, tref) = constants

    tke =
        ones(nxx, nyy, nzz) .* 0.1 .* (tref ^ 2.0) ./ (lref ^ 2.0) .*
        (rho .+ rhostrattfc)
    tte =
        ones(nxx, nyy, nzz) .* 0.1 .* (tref ^ 2.0) ./ (lref ^ 2.0) .*
        (rho .+ rhostrattfc)

    return TurbulencePredictands(tke, tte)
end
