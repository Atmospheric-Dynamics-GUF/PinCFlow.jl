"""
```julia
Ice{
    A <: IcePredictands,
    B <: IceIncrements,
    C <: IceAuxiliaries,
    D <: IceReconstructions,
    E <: IceFluxes,
}
```

Container for arrays needed by the ice-physics scheme.

```julia
Ice(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)::Ice
```

Construct an `Ice` instance, with array dimensions and initial values set according to the model configuration.

# Fields

  - `icepredictands::A`: Prognostic variables of the ice-physics scheme.

  - `iceincrements::B`: Runge-Kutta updates of the ice variables.

  - `iceauxiliaries::C`: Initial states of the ice variables.

  - `icereconstructions::D`: Reconstructions of the ice variables.

  - `icefluxes::E`: Fluxes of the ice variables.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `atmosphere`: Atmospheric-background fields.

  - `grid`: Collection of parameters and fields describing the grid.

  - `variables`: Container for arrays needed for the prediction of the prognostic variables.

# See also

  - [`PinCFlow.Types.IceTypes.IcePredictands`](@ref)

  - [`PinCFlow.Types.IceTypes.IceIncrements`](@ref)

  - [`PinCFlow.Types.IceTypes.IceAuxiliaries`](@ref)

  - [`PinCFlow.Types.IceTypes.IceReconstructions`](@ref)

  - [`PinCFlow.Types.IceTypes.IceFluxes`](@ref)
"""
struct Ice{
    A <: IcePredictands,
    B <: IceIncrements,
    C <: IceAuxiliaries,
    D <: IceReconstructions,
    E <: IceFluxes,
}
    icepredictands::A
    iceincrements::B
    iceauxiliaries::C
    icereconstructions::D
    icefluxes::E
end

function Ice(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)::Ice
    icepredictands = IcePredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        variables,
    )
    iceincrements = IceIncrements(namelists, domain)
    iceauxiliaries = IceAuxiliaries(icepredictands)
    icereconstructions = IceReconstructions(namelists, domain)
    icefluxes = IceFluxes(namelists, domain)

    return Ice(
        icepredictands,
        iceincrements,
        iceauxiliaries,
        icereconstructions,
        icefluxes,
    )
end
