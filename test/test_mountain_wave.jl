l2 = (
    n2 = 0.010105865f0,
    p = 5973.465f0,
    pip = 3.281434f-5,
    rhobar = 17.034365f0,
    rhop = 0.0011904717f0,
    t = 3600.0f0,
    thetabar = 13598.672f0,
    us = 447.23755f0,
    vs = 0.3299895f0,
    wts = 0.47778645f0,
    x = 18165.902f0,
    y = 18165.902f0,
    z = 364725.28f0,
    ztilde = 392459.12f0,
)
linf = (
    n2 = 0.00031957557f0,
    p = 321.18332f0,
    pip = 3.477205f-6,
    rhobar = 1.0363004f0,
    rhop = 0.00040524942f0,
    t = 3600.0f0,
    thetabar = 556.882f0,
    us = 10.135668f0,
    vs = 0.050961703f0,
    wts = 0.13209467f0,
    x = 9000.0f0,
    y = 9000.0f0,
    z = 19001.666f0,
    ztilde = 20000.0f0,
)
reference = (l2, linf)

@testset "mountain_wave" begin
    test_example(
        joinpath(submit_directory, "mountain_wave.jl"),
        reference,
        :x_size => 10,
        :y_size => 10,
        :z_size => 10,
        :output => OutputNamelist(; prepare_restart = true),
    )
end
