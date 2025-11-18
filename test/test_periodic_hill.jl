l2 = (
    n2 = 0.0031957547f0,
    p = 1888.8247f0,
    pip = 1.5724572f-5,
    rhobar = 5.38609f0,
    rhop = 0.00034377302f0,
    t = 3600.0f0,
    thetabar = 4300.3516f0,
    us = 141.42049f0,
    vs = 0.0f0,
    wts = 0.1385239f0,
    x = 18165.902f0,
    y = 0.0f0,
    z = 115340.16f0,
    ztilde = 124110.04f0,
)
linf = (
    n2 = 0.00031957557f0,
    p = 321.19247f0,
    pip = 4.693946f-6,
    rhobar = 1.0363418f0,
    rhop = 7.6056924f-5,
    t = 3600.0f0,
    thetabar = 556.86066f0,
    us = 10.052516f0,
    vs = 0.0f0,
    wts = 0.024587423f0,
    x = 9000.0f0,
    y = 0.0f0,
    z = 19000.488f0,
    ztilde = 20000.0f0,
)
reference = (l2, linf)

@testset "periodic_hill" begin
<<<<<<< HEAD
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
        x_size = 8,
        y_size = 1,
        z_size = 8,
=======
    test_example(
        joinpath(submit_directory, "periodic_hill.jl"),
        reference,
        :x_size => 10,
        :z_size => 10,
        :output => OutputNamelist(; prepare_restart = true),
>>>>>>> afc93468
    )
end
