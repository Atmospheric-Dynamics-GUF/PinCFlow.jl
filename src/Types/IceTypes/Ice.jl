struct Ice{
    A <: IcePredictands,
    B <: IceTendencies,
    C <: IceAuxiliaries,
    D <: IceReconstructions,
    E <: IceFluxes,
    F <: IceSource,
    }
    
    icepredictands::A
    icetendencies::B
    iceauxiliaries::C
    icereconstructions::D
    icefluxes::E
    icesource::F
end

function Ice(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)
    icepredictands = IcePredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        variables,
    )
    icetendencies = IceTendencies(namelists, domain)
    iceauxiliaries = IceAuxiliaries(icepredictands)
    icereconstructions = IceReconstructions(namelists, domain)
    icefluxes = IceFluxes(namelists, domain)
    icesource = IceSource(domain)

    return Ice(
        icepredictands,
        icetendencies,
        iceauxiliaries,
        icereconstructions,
        icefluxes,
        icesource
    )
end
