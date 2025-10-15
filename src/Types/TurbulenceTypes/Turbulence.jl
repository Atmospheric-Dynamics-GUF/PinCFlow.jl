"""
```julia
Turbulence{
    A <: TurbulencePredictands,
    B <: TurbulenceIncrements,
    C <: TurbulenceAuxiliaries,
    D <: TurbulenceReconstructions,
    E <: TurbulenceFluxes,
    F <: TurbulenceConstants,
}
```

Container for arrays and constants needed by the turbulence scheme.

```julia
Turbulence(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)::Turbulence
```

Construct a `Turbulence` instance, with array dimensions, initial values and constants set according to the model configuration.

# Fields

  - `turbulencepredictands::A`: Prognostic variables of the turbulence scheme.

  - `turbulenceincrements::B`: Runge-Kutta updates of the turbulence variables.

  - `turbulenceauxiliaries::C`: Background values.

  - `turbulencereconstructions::D`: Reconstructions of the turbulence variables.

  - `turbulencefluxes::E`: Fluxes of the turbulence variables.

  - `turbulenceconstants::F`: Constants used in the turbulence calculation

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `atmosphere`: Atmospheric-background fields.

  - `grid`: Collection of parameters and fields describing the grid.

  - `variables`: Container for arrays needed for the prediction of the prognostic variables.

# See also

  - [`PinCFlow.Types.TurbulenceTypes.TurbulencePredictands`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceIncrements`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceAuxiliaries`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceReconstructions`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceFluxes`](@ref)

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceConstants`](@ref)
"""
struct Turbulence{
    A <: TurbulencePredictands,
    B <: TurbulenceIncrements,
    C <: TurbulenceAuxiliaries,
    D <: TurbulenceReconstructions,
    E <: TurbulenceFluxes,
    F <: TurbulenceConstants,
}
    turbulencepredictands::A
    turbulenceincrements::B
    turbulenceauxiliaries::C
    turbulencereconstructions::D
    turbulencefluxes::E
    turbulenceconstants::F
end

function Turbulence(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)::Turbulence
    turbulencepredictands = TurbulencePredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        variables,
    )
    turbulenceincrements = TurbulenceIncrements(namelists, domain)
    turbulenceauxiliaries =
        TurbulenceAuxiliaries(turbulencepredictands, constants, domain)
    turbulencereconstructions = TurbulenceReconstructions(namelists, domain)
    turbulencefluxes = TurbulenceFluxes(namelists, domain)
    turbulenceconstants = TurbulenceConstants(namelists, constants)

    return Turbulence(
        turbulencepredictands,
        turbulenceincrements,
        turbulenceauxiliaries,
        turbulencereconstructions,
        turbulencefluxes,
        turbulenceconstants,
    )
end
