l2 = (
    n2 = 0.0f0,
    p = 1756.7067f0,
    pip = 0.00038125255f0,
    rhobar = 5.855689f0,
    rhop = 0.009143501f0,
    t = 3600.0f0,
    thetabar = 3000.0f0,
    us = 8.429902f0,
    vs = 2.9628015f0,
    wts = 9.879402f0,
    x = 18165.902f0,
    y = 0.0f0,
    z = 115325.625f0,
    ztilde = 124096.734f0,
)
linf = (
    n2 = 0.0f0,
    p = 320.76392f0,
    pip = 9.277299f-5,
    rhobar = 1.069213f0,
    rhop = 0.0058253715f0,
    t = 3600.0f0,
    thetabar = 300.0f0,
    us = 2.1983979f0,
    vs = 0.6476124f0,
    wts = 3.5191288f0,
    x = 9000.0f0,
    y = 0.0f0,
    z = 19000.0f0,
    ztilde = 20000.0f0,
)
reference = (l2, linf)

@testset "Cold bubble" begin
    test_example(
        joinpath(scripts_directory, "cold_bubble.jl"),
        reference,
        :x_size => 10,
        :z_size => 10,
        :output => OutputNamelist(; prepare_restart = true);
        update_references,
    )
end
