module TestExamples

include("test_pincflow.jl")

@testset "periodic_hill" begin
    @test_pinc_include(
        joinpath(EXAMPLES_DIR, "periodic_hill.jl"),
        l2 = [
            0.006555150635456577,
            4.159469879942826,
            0.0,
            0.004301707226655797,
            0.00015665630037031673,
            0.006531937247859811,
        ],
        linf = [
            0.000164891523374835,
            0.03460092539380474,
            0.0,
            0.00011599587150920556,
            3.21250293586638e-6,
            0.0001643172503966307,
        ]
    )
end

@testset "mountain wave" begin
    @test_pinc_include(
        joinpath(EXAMPLES_DIR, "mountain_wave.jl"),
        l2 = [
            0.001031058600826577,
            1.7906430623509666,
            0.0011272119805721063,
            0.001497573424424862,
            2.1426668808783215e-5,
            0.0010290483234598807,
        ],
        linf = [
            0.00024712277854250565,
            0.03446943529063215,
            0.00010273999644613423,
            0.000259944482327785,
            2.3151923637276867e-6,
            0.00024665016295258605,
        ],
        domain = DomainNamelist(;
            sizex = 8,
            sizey = 8,
            sizez = 8,
            lx_dim = 2.0E+4,
            ly_dim = 2.0E+4,
            lz_dim = 2.0E+4,
        )
    )
end

end # end module
