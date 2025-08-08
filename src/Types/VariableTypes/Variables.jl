"""
```julia
Variables{
    A <: Predictands,
    B <: Tendencies,
    C <: Backups,
    D <: Auxiliaries,
    E <: Reconstructions,
    F <: Fluxes,
}
```

Container for arrays needed for the prediction of the prognostic variables.

```julia
Variables(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
)
```

Construct a `Variables` instance, with array dimensions and initial values set according to the model configuration.

# Fields

- `predictands::A`: Prognostic variables.
- `tendencies::B`: Runge-Kutta updates and pressure correction.
- `backups::C`: Backups of the prognostic variables needed in the semi-implicit time scheme.
- `auxiliaries::D`: Auxiliary array needed in the reconstruction.
- `reconstructions::E`: Reconstructions of the prognostic variables.
- `fluxes::F`: Fluxes of the prognostic variables.

# Arguments

- `namelists`: Namelists with all model parameters.
- `constants`: Physical constants and reference values.
- `domain`: Collection of domain-decomposition and MPI-communication parameters.
- `atmosphere`: Atmospheric-background fields.

# See also

- [`PinCFlow.Types.VariableTypes.Predictands`](@ref)
- [`PinCFlow.Types.VariableTypes.Tendencies`](@ref)
- [`PinCFlow.Types.VariableTypes.Backups`](@ref)
- [`PinCFlow.Types.VariableTypes.Auxiliaries`](@ref)
- [`PinCFlow.Types.VariableTypes.Reconstructions`](@ref)
- [`PinCFlow.Types.VariableTypes.Fluxes`](@ref)
"""
struct Variables{
    A <: Predictands,
    B <: Tendencies,
    C <: Backups,
    D <: Auxiliaries,
    E <: Reconstructions,
    F <: Fluxes,
}
    predictands::A
    tendencies::B
    backups::C
    auxiliaries::D
    reconstructions::E
    fluxes::F
end

function Variables(
    namelists::Namelists,
    constants::Constants,
    domain::Domain,
    atmosphere::Atmosphere,
    grid::Grid,
)

    # Initialize all fields.
    predictands = Predictands(namelists, constants, domain, atmosphere, grid)
    tendencies = Tendencies(namelists, domain)
    backups = Backups(domain)
    auxiliaries = Auxiliaries(domain)
    reconstructions = Reconstructions(domain)
    fluxes = Fluxes(namelists, domain)

    # Return a Variables instance.
    return Variables(
        predictands,
        tendencies,
        backups,
        auxiliaries,
        reconstructions,
        fluxes,
    )
end
