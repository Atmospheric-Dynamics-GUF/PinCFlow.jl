struct State{
    A <: Namelists,
    B <: Time,
    C <: Constants,
    D <: Domain,
    E <: Grid,
    F <: Atmosphere,
    G <: Sponge,
    H <: Poisson,
    I <: Variables,
    J <: WKB,
}
    namelists::A
    time::B
    constants::C
    domain::D
    grid::E
    atmosphere::F
    sponge::G
    poisson::H
    variables::I
    wkb::J
end

function State(namelists::Namelists)

    # Initialize everything.
    constants = Constants(namelists)
    time = Time()
    domain = Domain(namelists)
    grid = Grid(namelists, constants, domain)
    atmosphere = Atmosphere(namelists, constants, domain, grid)
    sponge = Sponge(namelists, domain, grid)
    poisson = Poisson(namelists, domain)
    variables = Variables(namelists, constants, domain, atmosphere)
    wkb = WKB(namelists, constants, domain, grid)

    # Return a State instance.
    return State(
        namelists,
        time,
        constants,
        domain,
        grid,
        atmosphere,
        sponge,
        poisson,
        variables,
        wkb,
    )
end
