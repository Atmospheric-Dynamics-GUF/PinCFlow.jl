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
end

function State(namelists::Namelists)

  # Get parameters.
  (; model) = namelists.setting
  (; background) = namelists.atmosphere

  # Initialize everything.
  constants = Constants(namelists)
  time = Time()
  domain = Domain(namelists)
  grid = Grid(namelists, constants, domain)
  atmosphere = Atmosphere(namelists, constants, domain, grid, model, background)
  sponge = Sponge(namelists, domain, grid)
  poisson = Poisson(namelists, domain)
  variables = Variables(namelists, constants, domain)

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
  )
end
