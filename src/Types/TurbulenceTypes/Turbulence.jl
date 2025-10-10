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

Container for arrays needed for turbulence transport.

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

  - `turbulencepredictands::A`: Turbulence energies.

  - `turbulenceincrements::B`: Runge-Kutta updates of the turbulence energies.

  - `turbulenceauxiliaries::C`: Initial states of the turbulence energies.

  - `turbulencereconstructions::D`: Reconstructions of the turbulence energies.

  - `turbulencefluxes::E`: Fluxes of the turbulence energies.

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

  - [`PinCFlow.Types.TurbulenceTypes.TurbulenceForcings`](@ref)
"""
struct Turbulence{
    A <: TurbulencePredictands,
    B <: TurbulenceIncrements,
    C <: TurbulenceAuxiliaries,
    D <: TurbulenceReconstructions,
    E <: TurbulenceFluxes,
    F <: TurbulenceForcings,
}
    turbulencepredictands::A
    turbulenceincrements::B
    turbulenceauxiliaries::C
    turbulencereconstructions::D
    turbulencefluxes::E
    turbulenceforcings::F
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
    turbulenceauxiliaries = TurbulenceAuxiliaries(turbulencepredictands)
    turbulencereconstructions = TurbulenceReconstructions(namelists, domain)
    turbulencefluxes = TurbulenceFluxes(namelists, domain)
    turbulenceforcings = TurbulenceForcings(namelists, domain)

    return Turbulence(
        turbulencepredictands,
        turbulenceincrements,
        turbulenceauxiliaries,
        turbulencereconstructions,
        turbulencefluxes,
        turbulenceforcings,
    )
end
