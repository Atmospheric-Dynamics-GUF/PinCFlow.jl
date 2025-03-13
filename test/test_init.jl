using Setfield
@testset "Initialization" begin
    p = default_parameters()
    domain = DomainParameters(sizex=10)
    p = @set p.domain = domain
    @test p.domain.sizex == 10
    @test p.domain.sizey == 1
    m = Model(p)
end
