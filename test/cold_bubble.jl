l2 = (
    n2 = 0.0f0,
    p = 1756.7067f0,
    pip = 0.00038126486f0,
    rhobar = 5.855689f0,
    rhop = 0.0091435f0,
    t = 3600.0f0,
    thetabar = 3000.0f0,
    us = 8.428999f0,
    vs = 2.9697578f0,
    wts = 9.878598f0,
    x = 18165.902f0,
    y = 0.0f0,
    z = 115325.625f0,
    ztilde = 124096.734f0,
)
linf = (
    n2 = 0.0f0,
    p = 320.76392f0,
    pip = 9.276894f-5,
    rhobar = 1.069213f0,
    rhop = 0.0058253715f0,
    t = 3600.0f0,
    thetabar = 300.0f0,
    us = 2.1978178f0,
    vs = 0.6512958f0,
    wts = 3.5186996f0,
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
