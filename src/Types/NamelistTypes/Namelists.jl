"""
```julia
Namelists
```

Represents all configurable parameters for the simulation.

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
  - [`PinCFlow.Types.NamelistTypes.IceNamelist`](@ref)
  - [`PinCFlow.Types.NamelistTypes.TurbulenceNamelist`](@ref)
  - [`PinCFlow.Types.NamelistTypes.WavePacketNamelist`](@ref)
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
    K <: IceNamelist,
    L <: TurbulenceNamelist,
    M <: WavePacketNamelist,
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
    ice::K
    turbulence::L
    wavepacket::M
end

"""
```julia
Namelists(;
    domain = DomainNamelist(),
    output = OutputNamelist(),
    setting = SettingNamelist(),
    discretization = DiscretizationNamelist(),
    poisson = PoissonNamelist(),
    atmosphere = AtmosphereNamelist(),
    grid = GridNamelist(),
    sponge = SpongeNamelist(),
    wkb = WKBNamelist(),
)
```

Create a new `Namelist` instance with the provided namelists. Omitted namelists will be initialized with default values.

# Returns

  - `::Namelists`: `Namelists` instance.
"""
function Namelists(;
    domain = DomainNamelist(),
    output = OutputNamelist(),
    setting = SettingNamelist(),
    discretization = DiscretizationNamelist(),
    poisson = PoissonNamelist(),
    atmosphere = AtmosphereNamelist(),
    grid = GridNamelist(),
    sponge = SpongeNamelist(),
    wkb = WKBNamelist(),
    tracer = TracerNamelist(),
    ice = IceNamelist(),
    turbulence = TurbulenceNamelist(),
    wavepacket = WavePacketNamelist(),
)
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
        ice,
        turbulence,
        wavepacket,
    )
end
