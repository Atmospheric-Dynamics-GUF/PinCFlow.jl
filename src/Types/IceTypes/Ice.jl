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
    F <: IceSource,
    G <: IceConstants,
    H <: GW,
    I <: SgsGW,
    J <: SgsPredictands,
    K <: SubGrid
    }
    
    icepredictands::A
    iceincrements::B
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
    iceincrements = IceIncrements(namelists, domain)
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
        iceincrements,
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
