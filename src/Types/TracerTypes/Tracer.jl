"""
```julia
Tracer{
    A <: TracerPredictands,
    B <: TracerIncrements,
    C <: TracerReconstructions,
    D <: TracerFluxes,
    E <: TracerWKBIntegrals,
    F <: TracerWKBTendencies,
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

  - `tracerreconstructions::C`: Reconstructions of the tracers.

  - `tracerfluxes::D`: Fluxes of the tracers.

  - `tracerwkbintegrals::E`: Integrals of gravity-wave induced tracer fluxes.

  - `tracerwkbtendencies::F`: Tracer impact of unresolved gravity waves.

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

  - [`PinCFlow.Types.TracerTypes.TracerReconstructions`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerFluxes`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerWKBIntegrals`](@ref)

  - [`PinCFlow.Types.TracerTypes.TracerWKBTendencies`](@ref)
"""
struct Tracer{
    A <: TracerPredictands,
    B <: TracerIncrements,
    C <: TracerReconstructions,
    D <: TracerFluxes,
    E <: TracerWKBIntegrals,
    F <: TracerWKBTendencies,
}
    tracerpredictands::A
    tracerincrements::B
    tracerreconstructions::C
    tracerfluxes::D
    tracerwkbintegrals::E
    tracerwkbtendencies::F
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
    tracerreconstructions = TracerReconstructions(namelists, domain)
    tracerfluxes = TracerFluxes(namelists, domain)
    tracerwkbintegrals = TracerWKBIntegrals(namelists, domain)
    tracerwkbtendencies = TracerWKBTendencies(namelists, domain)

    return Tracer(
        tracerpredictands,
        tracerincrements,
        tracerreconstructions,
        tracerfluxes,
        tracerwkbintegrals,
        tracerwkbtendencies,
    )
end
