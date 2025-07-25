struct Ice{
    A <: IcePredictands,
    B <: IceTendencies,
    C <: IceAuxiliaries,
    D <: IceReconstructions,
    E <: IceFluxes,
    F <: IceSource,
    G <: IceConstants,
    H <: GW,
    I <: SgsGW,
    J <: SgsPredictands
    }
    
    icepredictands::A
    icetendencies::B
    iceauxiliaries::C
    icereconstructions::D
    icefluxes::E
    icesource::F
    iceconstants::G
    gw::H
    sgs::I
    sgspredictands::J
end

function Ice(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)
    iceconstants = IceConstants(constants)

    icepredictands = IcePredictands(
        namelists,
        constants,
        domain,
        atmosphere,
        grid,
        variables,
        iceconstants 
    )
    icetendencies = IceTendencies(namelists, domain)
    iceauxiliaries = IceAuxiliaries(icepredictands)
    icereconstructions = IceReconstructions(namelists, domain)
    icefluxes = IceFluxes(namelists, domain)
    icesource = IceSource(domain)
    gw = GW(namelists, domain)
    sgs = SgsGW(namelists, domain)
    # sgs_predictands = SgsPredictands(
    #     namelists,
    #     constants,
    #     domain,
    #     atmosphere,
    #     grid)
    sgspredictands = SgsPredictands(icepredictands)

    return Ice(
        icepredictands,
        icetendencies,
        iceauxiliaries,
        icereconstructions,
        icefluxes,
        icesource,
        iceconstants,
        gw,
        sgs,
        sgspredictands
    )
end
