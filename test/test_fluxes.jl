if isdefined(@__MODULE__, :LanguageServer)
  include("../src/PinCFlow.jl")
end
using Setfield

@testset_broken "Test fluxes.jl" begin
  namelists = PinCFlow.Namelists()
  domain = PinCFlow.DomainNamelist(; sizex = 30, sizey = 1, sizez = 10)
  @set namelists.domain = domain
  state = PinCFlow.State(namelists)

  normalizer = 0.01 # To get smaller values where arithmetic is more accurate
  set_lin_array_normalizer!(arr) = set_lin_array!(arr; normalizer = normalizer)
  set_lin_array_normalizer!.((
    state.atmosphere.pstrattfc,
    state.atmosphere.rhostrattfc,
    state.grid.topography_surface,
  ))

  var = state.variables.predictands
  set_lin_array_normalizer!(var.rho)
  set_lin_array_normalizer!(var.rhop)
  set_lin_array_normalizer!(var.u)
  set_lin_array_normalizer!(var.v)
  set_lin_array_normalizer!(var.w)

  PinCFlow.reconstruct!(state)

  test_arr(
    state.variables.reconstructions.rhotilde,
    15120.0,
    122.96340919151518,
    1.0;
    tol = 1e-14,
  )
  test_arr(
    state.variables.reconstructions.rhoptilde,
    15120.0,
    122.96340919151518,
    1.0;
    tol = 1e-14,
  )
  test_arr(
    state.variables.utilde,
    9374.400000000001,
    80.74848357709314,
    1.1300000000000001;
    tol = 1e-14,
  )
  test_arr(
    state.variables.vtilde,
    9374.400000000001,
    80.74805260809762,
    1.1300000000000001;
    tol = 1e-14,
  )
  test_arr(
    state.variables.wtilde,
    8703.234242111952,
    73.30783802039744,
    0.8492267931815118;
    tol = 1e-14,
  )

  set_lin_array_normalizer!(state.variables.backups.rho)
  set_lin_array_normalizer!(state.variables.backups.rhop)
  set_lin_array_normalizer!(state.variables.backups.u)
  set_lin_array_normalizer!(state.variables.backups.v)
  set_lin_array_normalizer!(state.variables.backups.w)

  PinCFlow.compute_fluxes!(semi)

  flux = state.variables.fluxes
  test_arr(
    flux.phirho,
    234.48977709288997,
    7.465169904723465,
    0.4313585705042774;
    tol = 1e-14,
  )
  test_arr(
    flux.phirhop,
    117.24488854644498,
    3.7325849523617327,
    0.2156792852521387;
    tol = 1e-14,
  )
  test_arr(
    flux.phiu,
    88.4864783397629,
    3.0875380464248465,
    0.22591343710336997;
    tol = 1e-14,
  )
  test_arr(
    flux.phiv,
    149.20376149456763,
    3.964286453926213,
    0.22591343710336997;
    tol = 1e-14,
  )
  test_arr(
    flux.phiw,
    81.05895692784512,
    2.705278309488443,
    0.18549395063610127;
    tol = 1e-14,
  )
end
