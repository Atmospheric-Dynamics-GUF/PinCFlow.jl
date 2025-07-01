"""
```julia
Namelists
```

Represents all configurable parameters for the simulation.
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
    )
end
