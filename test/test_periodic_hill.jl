module TestExamples

include("test_pincflow.jl")

@testset "periodic_hill" begin
    @test_example(
        joinpath(submit_dir, "periodic_hill.jl"),
        l2 = Float32[
            0.0025576192,
            1509.1896,
            4.2984476,
            3600.0,
            3439.5984,
            0.06355105,
        ],
        linf = Float32[
            0.00031970255,
            314.72037,
            1.0072244,
            3600.0,
            552.3491,
            0.012468175,
        ],
        sizex = 8,
        sizey = 1,
        sizez = 8,
    )
end

end # end module
