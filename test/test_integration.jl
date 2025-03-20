
using PinCFlow
domain = PinCFlow.DomainNamelist(;
  sizex = 30,
  sizey = 1,
  sizez = 10,
  lx_dim = [0.0, 60000.0],
  ly_dim = [0.0, 40000.0],
  lz_dim = [0.0, 20000.0],
)

atmosphere =
  PinCFlow.AtmosphereNamelist(; backgroundflow_dim = [10.0, 0.0, 0.0])
grid = PinCFlow.GridNamelist(;
  mountain_case = 3,
  mountainheight_dim = 400.0,
  mountainwidth_dim = 1000.0,
)
namelists = Namelists(; domain = domain, grid = grid, atmosphere = atmosphere)

# this breaks after two iterations!
state = PinCFlow.integrate(namelists)

@test get_norms(state.variables.predictands.u) ==
      (171.9769087732426, 2.4244442121203518, 0.03759161012798296)
@test get_norms(state.variables.predictands.v) == (0.0, 0.0, 0.0)
@test get_norms(state.variables.predictands.w) ==
      (0.46842763331095266, 0.020501897078321726, 0.0020527458908644945)
