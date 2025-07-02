struct Tracer{
    A <: TracerPredictands,
    B <: TracerTendencies,
    C <: TracerAuxiliaries,
    D <: TracerReconstructions,
    E <: TracerFluxes,
}
    tracerpredictands::A
    tracertendencies::B
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
)
    tracerpredictands = TracerPredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        variables,
    )
    tracertendencies = TracerTendencies(namelists, domain)
    tracerauxiliaries = TracerAuxiliaries(tracerpredictands)
    tracerreconstructions = TracerReconstructions(namelists, domain)
    tracerfluxes = TracerFluxes(namelists, domain)

    return Tracer(
        tracerpredictands,
        tracertendencies,
        tracerauxiliaries,
        tracerreconstructions,
        tracerfluxes,
    )
end
