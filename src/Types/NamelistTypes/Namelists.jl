"""
```julia
Namelists{
    A <: DomainNamelist,
    B <: OutputNamelist,
<<<<<<< HEAD
    C <: SettingNamelist,
    D <: DiscretizationNamelist,
    E <: PoissonNamelist,
    F <: AtmosphereNamelist,
    G <: GridNamelist,
    H <: SpongeNamelist,
    I <: WKBNamelist,
    J <: TracerNamelist,
    K <: IceNamelist,
    M <: WavePacketNamelist,
=======
    C <: DiscretizationNamelist,
    D <: PoissonNamelist,
    E <: AtmosphereNamelist,
    F <: GridNamelist,
    G <: SpongeNamelist,
    H <: WKBNamelist,
    I <: TracerNamelist,
>>>>>>> cf395edbf2
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
    ice::IceNamelist = IceNamelist(),
    wavepacket::WavePacketNamelist = WavePacketNamelist(),
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

  - `ice::K`: Namelist for parameters configuring ice physics.

  - `wavepacket::M`: Namelist for parameters used for the `WavePacket` test case.

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

  - [`PinCFlow.Types.NamelistTypes.IceNamelist`](@ref)

  - [`PinCFlow.Types.NamelistTypes.WavePacketNamelist`](@ref)
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
    K <: IceNamelist,
    M <: WavePacketNamelist,
    N <: MultiWavePacketNamelist,
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
    ice::K
    wavepacket::M
    multiwavepackets::N
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
    ice::IceNamelist = IceNamelist(),
    wavepacket::WavePacketNamelist = WavePacketNamelist(),
    multiwavepackets::MultiWavePacketNamelist = MultiWavePacketNamelist(),
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
        ice,
        wavepacket,
        multiwavepackets,
    )
end
