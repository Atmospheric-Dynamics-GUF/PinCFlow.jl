"""
```julia
Namelists{
    A <: DomainNamelist,
    B <: OutputNamelist,
    C <: DiscretizationNamelist,
    D <: PoissonNamelist,
    E <: AtmosphereNamelist,
    F <: GridNamelist,
    G <: SpongeNamelist,
    H <: WKBNamelist,
    I <: TracerNamelist,
}
```

Collection of all configurable model parameters.

```julia
Namelists(;
    domain::DomainNamelist = DomainNamelist(),
    output::OutputNamelist = OutputNamelist(),
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

  - `discretization::C`: Namelist for parameters describing discretization.

  - `poisson::D`: Namelist for parameters used by the Poisson solver.

  - `atmosphere::E`: Namelist for parameters describing the atmospheric background.

  - `grid::F`: Namelist for parameters describing the grid.

  - `sponge::G`: Namelist for parameters describing the sponge.

  - `wkb::H`: Namelist for parameters used by MSGWaM.

  - `tracer::I`: Namelist for parameters configuring tracer transport.

# See also

  - [`PinCFlow.Types.NamelistTypes.DomainNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.OutputNamelist`](@ref)

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
    C <: DiscretizationNamelist,
    D <: PoissonNamelist,
    E <: AtmosphereNamelist,
    F <: GridNamelist,
    G <: SpongeNamelist,
    H <: WKBNamelist,
    I <: TracerNamelist,
}
    domain::A
    output::B
    discretization::C
    poisson::D
    atmosphere::E
    grid::F
    sponge::G
    wkb::H
    tracer::I
end

function Namelists(;
    domain::DomainNamelist = DomainNamelist(),
    output::OutputNamelist = OutputNamelist(),
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
        discretization,
        poisson,
        atmosphere,
        grid,
        sponge,
        wkb,
        tracer,
    )
end
