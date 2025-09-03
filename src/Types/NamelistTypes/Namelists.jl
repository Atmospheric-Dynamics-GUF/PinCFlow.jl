"""
```julia
Namelists{
    A <: DomainNamelist,
    B <: OutputNamelist,
    C <: SettingNamelist,
    D <: DiscretizationNamelist,
    E <: PoissonNamelist,
    F <: AtmosphereNamelist,
    G <: GridNamelist,
    H <: SpongeNamelist,
    I <: WKBNamelist,
    J <: TracerNamelist,
}
```

Collection of all configurable model parameters.

```julia
Namelists(;
    domain::DomainNamelist = DomainNamelist(),
    output::OutputNamelist = OutputNamelist(),
    setting::SettingNamelist = SettingNamelist(),
    discretization::DiscretizationNamelist = DiscretizationNamelist(),
    poisson::PoissonNamelist = PoissonNamelist(),
    atmosphere::AtmosphereNamelist = AtmosphereNamelist(),
    grid::GridNamelist = GridNamelist(),
    sponge::SpongeNamelist = SpongeNamelist(),
    wkb::WKBNamelist = WKBNamelist(),
    tracer::TracerNamelist = TracerNamelist(),
)::Namelists
```

Construct a `Namelists` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `domain::A`: Namelist for parameters describing the model domain.

  - `output::B`: Namelist for I/O parameters.

  - `setting::C`: Namelist for parameters describing the general model setting.

  - `discretization::D`: Namelist for parameters describing discretization.

  - `poisson::E`: Namelist for parameters used by the Poisson solver.

  - `atmosphere::F`: Namelist for parameters describing the atmospheric background.

  - `grid::G`: Namelist for parameters describing the grid.

  - `sponge::H`: Namelist for parameters describing the sponge.

  - `wkb::I`: Namelist for parameters used by MSGWaM.

  - `tracer::J`: Namelist for parameters configuring tracer transport.

# See also

  - [`PinCFlow.Types.NamelistTypes.DomainNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.OutputNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.SettingNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.DiscretizationNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.PoissonNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.AtmosphereNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.GridNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.SpongeNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.WKBNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.TracerNamelist`](@ref)
"""
struct Namelists{
    A <: DomainNamelist,
    B <: OutputNamelist,
    C <: SettingNamelist,
    D <: DiscretizationNamelist,
    E <: PoissonNamelist,
    F <: AtmosphereNamelist,
    G <: GridNamelist,
    H <: SpongeNamelist,
    I <: WKBNamelist,
    J <: TracerNamelist,
}
    domain::A
    output::B
    setting::C
    discretization::D
    poisson::E
    atmosphere::F
    grid::G
    sponge::H
    wkb::I
    tracer::J
end

function Namelists(;
    domain::DomainNamelist = DomainNamelist(),
    output::OutputNamelist = OutputNamelist(),
    setting::SettingNamelist = SettingNamelist(),
    discretization::DiscretizationNamelist = DiscretizationNamelist(),
    poisson::PoissonNamelist = PoissonNamelist(),
    atmosphere::AtmosphereNamelist = AtmosphereNamelist(),
    grid::GridNamelist = GridNamelist(),
    sponge::SpongeNamelist = SpongeNamelist(),
    wkb::WKBNamelist = WKBNamelist(),
    tracer::TracerNamelist = TracerNamelist(),
)::Namelists
    return Namelists(
        domain,
        output,
        setting,
        discretization,
        poisson,
        atmosphere,
        grid,
        sponge,
        wkb,
        tracer,
    )
end
