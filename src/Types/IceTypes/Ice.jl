"""
```julia
Ice{
    A <: IcePredictands,
    B <: IceTendencies,
    C <: IceAuxiliaries,
    D <: IceReconstructions,
    E <: IceFluxes,
}
```

```julia
Ice(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
    variables::Variables,
)
```
"""
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
    J <: SgsPredictands,
    K <: SubGrid
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
    subgrid::K
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
    icesource = IceSource(namelists, domain)
    gw = GW(namelists, domain)
    subgrid = SubGrid(namelists, domain, grid)
    sgs = SgsGW(namelists, domain, subgrid)
    sgspredictands = SgsPredictands(namelists, domain, icepredictands, subgrid)

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
        sgspredictands,
        subgrid
    )
end
