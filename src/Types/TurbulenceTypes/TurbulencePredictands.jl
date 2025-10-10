"""
```julia
TurbulencePredictands{A <: AbstractArray{<:AbstractFloat, 3}}
```

Arrays for turbulence energies.

```julia
TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)::TurbulencePredictands
```

Construct a `TurbulencePredictands` instance with dimensions and initial values depending on the general configuration of turbulence transport, by dispatching to the appropriate method.

```julia
TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulence_setup::AbstractTurbulence,
    variables::Variables,
)::TurbulencePredictands
```

Construct a `TurbulencePredictands` instance with zero-size arrays for configurations without turbulence transport.

```julia
TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulence_setup::LinearTurbulence,
    variables::Variables,
)::TurbulencePredictands
```

Construct a `TurbulencePredictands` instance with an initialized non-dimensional turbulence linearly increasing with altitude. The turbulence field is multiplied by the density.

# Fields

  - `chi::A`: Non-dimensional turbulence.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `atmosphere`: Atmospheric-background fields.

  - `grid`: Collection of parameters and fields describing the grid.

  - `turbulence_setup`: General turbulence-transport configuration.

  - `variables`: Container for arrays needed for the prediction of the prognostic variables.
"""
struct TurbulencePredictands{A <: AbstractArray{<:AbstractFloat, 3}}
    chi::A
end

function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)::TurbulencePredictands
    (; turbulence_setup) = namelists.turbulence

    return TurbulencePredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        turbulence_setup,
        variables,
    )
end

function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulence_setup::AbstractTurbulence,
    variables::Variables,
)::TurbulencePredictands
    return TurbulencePredictands(
        [zeros(0, 0, 0) for field in fieldnames(TurbulencePredictands)]...,
    )
end

function TurbulencePredictands(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    turbulence_setup::LinearTurbulence,
    variables::Variables,
)::TurbulencePredictands
    (; nxx, nyy, nzz) = domain
    (; zc) = grid
    (; lref) = constants
    (; rhobar) = atmosphere
    (; rho) = variables.predictands
    (; lref) = constants
    (; alphaturbulence) = namelists.turbulence

    chi = zeros(nxx, nyy, nzz)
    chi .= alphaturbulence .* lref .* zc

    chi .= chi .* (rho .+ rhobar)

    return TurbulencePredictands(chi)
end
