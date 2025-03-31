@testset "Namelists" begin
    @testset "Domainnamelist" begin
  namelist = PinCFlow.Namelists()
  @test namelist.domain.sizex == 4 && namelist.domain.sizey == 4 && namelist.domain.sizez == 4
  @test namelist.domain.nprocx == 1 && namelist.domain.nprocy == 1
  end
end
