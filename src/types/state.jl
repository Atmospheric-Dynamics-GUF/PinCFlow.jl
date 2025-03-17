struct State{
  A <: Namelists,
  B <: Time,
  C <: Constants,
  D <: Domain,
  E <: Grid,
  F <: Atmosphere,
  G <: Operator,
  H <: Variables,
}
  namelists::A
  time::B
  constants::C
  domain::D
  grid::E
  atmosphere::F
  operator::G
  variables::H
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
  operator = Operator(domain)
  variables = Variables(namelists, constants, domain)

  # Return a State instance.
  return State(
    namelists,
    time,
    constants,
    domain,
    grid,
    atmosphere,
    operator,
    variables,
  )
end
