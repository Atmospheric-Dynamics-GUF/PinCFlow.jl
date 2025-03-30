struct Namelists{
  A <: DomainNamelist,
  B <: OutputNamelist,
  C <: SettingNamelist,
  D <: DiscretizationNamelist,
  E <: PoissonNamelist,
  F <: AtmosphereNamelist,
  G <: GridNamelist,
  H <: SpongeNamelist,
}
  domain::A
  output::B
  setting::C
  discretization::D
  poisson::E
  atmosphere::F
  grid::G
  sponge::H
end

function Namelists(;
  domain = DomainNamelist(),
  output = OutputNamelist(),
  setting = SettingNamelist(),
  discretization = DiscretizationNamelist(),
  poisson = PoissonNamelist(),
  atmosphere = AtmosphereNamelist(),
  grid = GridNamelist(),
  sponge = SpongeNamelist(),
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
  )
end
