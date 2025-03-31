
@testset "Atmosphere at rest" begin
  namelists = Namelists()
  state = PinCFlow.integrate(namelists)
  @test all(state.variables.predictands.u .== 0)
  @test all(state.variables.predictands.v .== 0)
  @test all(state.variables.predictands.w .== 0)
end
