l2 = (
    n2 = 0.0f0,
    p = 1756.7067f0,
    pip = 0.0013329095f0,
    rhobar = 5.855689f0,
    rhop = 0.010418067f0,
    t = 3600.0f0,
    thetabar = 3000.0f0,
    us = 37.253185f0,
    vs = 5.496259f0,
    wts = 5.7381635f0,
    x = 363318.03f0,
    y = 0.0f0,
    z = 115325.625f0,
    ztilde = 124096.734f0,
)
linf = (
    n2 = 0.0f0,
    p = 320.76392f0,
    pip = 0.00034201576f0,
    rhobar = 1.069213f0,
    rhop = 0.0058253715f0,
    t = 3600.0f0,
    thetabar = 300.0f0,
    us = 9.882739f0,
    vs = 1.0850737f0,
    wts = 1.4321839f0,
    x = 180000.0f0,
    y = 0.0f0,
    z = 19000.0f0,
    ztilde = 20000.0f0,
)
reference = (l2, linf)

@testset "Cold bubble" begin
    test_example(
        joinpath(submit_directory, "cold_bubble.jl"),
        reference,
        :x_size => 10,
        :z_size => 10,
        :output => OutputNamelist(; prepare_restart = true);
        update_references,
    )
end
