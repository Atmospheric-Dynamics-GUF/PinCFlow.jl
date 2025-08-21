"""
```julia
Tracer{
    A <: TracerPredictands,
    B <: TracerIncrements,
    C <: TracerAuxiliaries,
    D <: TracerReconstructions,
    E <: TracerFluxes,
}
```

Container for arrays needed for tracer transport.

```julia
Tracer(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)::Tracer
```

Construct a `Tracer` instance, with array dimensions and initial values set according to the model configuration.

# Fields

  - `tracerpredictands::A`: Tracers.

  - `tracerincrements::B`: Runge-Kutta updates of the tracers.

  - `tracerauxiliaries::C`: Initial states of the tracers.

  - `tracerreconstructions::D`: Reconstructions of the tracers.

  - `tracerfluxes::E`: Fluxes of the tracers.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `constants`: Physical constants and reference values.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `atmosphere`: Atmospheric-background fields.

  - `grid`: Collection of parameters and fields describing the grid.

  - `variables`: Container for arrays needed for the prediction of the prognostic variables.

# See also

  - [`PinCFlow.Types.TracerTypes.TracerPredictands`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerIncrements`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerAuxiliaries`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerReconstructions`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerFluxes`](@ref)
"""
struct Tracer{
    A <: TracerPredictands,
    B <: TracerIncrements,
    C <: TracerAuxiliaries,
    D <: TracerReconstructions,
    E <: TracerFluxes,
}
    tracerpredictands::A
    tracerincrements::B
    tracerauxiliaries::C
    tracerreconstructions::D
    tracerfluxes::E
end

function Tracer(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)::Tracer
    tracerpredictands = TracerPredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        variables,
    )
    tracerincrements = TracerIncrements(namelists, domain)
    tracerauxiliaries = TracerAuxiliaries(tracerpredictands)
    tracerreconstructions = TracerReconstructions(namelists, domain)
    tracerfluxes = TracerFluxes(namelists, domain)

    return Tracer(
        tracerpredictands,
        tracerincrements,
        tracerauxiliaries,
        tracerreconstructions,
        tracerfluxes,
    )
end
