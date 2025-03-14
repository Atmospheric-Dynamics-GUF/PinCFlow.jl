struct Namelists{
    A<:DomainNamelist,
    B<:OutputNamelist,
    C<:SettingNamelist,
    D<:DiscretizationNamelist,
    E<:PoissonNamelist,
    F<:AtmosphereNamelist,
    G<:GridNamelist,
    H<:SpongeNamelist,
    I<:BoundariesNamelist,
}
    domain::A
    output::B
    setting::C
    discretization::D
    poisson::E
    atmosphere::F
    grid::G
    sponge::H
    boundaries::I
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
  boundaries = BoundariesNamelist(),
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
    boundaries,
  )
end
