l2 = (
    n2 = 0.010105865f0,
    p = 5973.465f0,
    pip = 4.3203294f-5,
    rhobar = 17.034365f0,
    rhop = 0.0011906031f0,
    t = 3600.0f0,
    thetabar = 13598.672f0,
    us = 447.23752f0,
    vs = 0.33409667f0,
    wts = 0.5572854f0,
    x = 18165.902f0,
    y = 18165.902f0,
    z = 364725.28f0,
    ztilde = 392459.12f0,
)
linf = (
    n2 = 0.00031957557f0,
    p = 321.18332f0,
    pip = 1.0270111f-5,
    rhobar = 1.0363004f0,
    rhop = 0.0004052551f0,
    t = 3600.0f0,
    thetabar = 556.882f0,
    us = 10.1355295f0,
    vs = 0.05098579f0,
    wts = 0.13203365f0,
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
