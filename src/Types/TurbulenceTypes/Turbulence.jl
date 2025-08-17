"""
```julia
Turbulence{
    A <: TurbulencePredictands,
    B <: TurbulenceIncrements,
    C <: TurbulenceAuxiliaries,
    D <: TurbulenceReconstructions,
    E <: TurbulenceFluxes,
}
```

Container for arrays needed by the turbulence scheme.

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

Construct a `Turbulence` instance, with array dimensions and initial values set according to the model configuration.

# Fields

  - `turbulencepredictands::A`: Prognostic variables of the turbulence scheme.

  - `turbulenceincrements::B`: Runge-Kutta updates of the turbulence variables.

  - `turbulenceauxiliaries::C`: Background values.

  - `turbulencereconstructions::D`: Reconstructions of the turbulence variables.

  - `turbulencefluxes::E`: Fluxes of the turbulence variables.

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
"""
struct Turbulence{
    A <: TurbulencePredictands,
    B <: TurbulenceIncrements,
    C <: TurbulenceAuxiliaries,
    D <: TurbulenceReconstructions,
    E <: TurbulenceFluxes,
}
    turbulencepredictands::A
    turbulenceincrements::B
    turbulenceauxiliaries::C
    turbulencereconstructions::D
    turbulencefluxes::E
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
    turbulenceauxiliaries = TurbulenceAuxiliaries(constants)
    turbulencereconstructions = TurbulenceReconstructions(namelists, domain)
    turbulencefluxes = TurbulenceFluxes(namelists, domain)

    return Turbulence(
        turbulencepredictands,
        turbulenceincrements,
        turbulenceauxiliaries,
        turbulencereconstructions,
        turbulencefluxes,
    )
end
