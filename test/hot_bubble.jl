l2 = (
    n2 = 0.0f0,
    p = 2484.3525f0,
    pip = 0.002410898f0,
    rhobar = 5.855689f0,
    rhop = 0.009429613f0,
    t = 3600.0f0,
    thetabar = 3000.0f0,
    us = 61.343918f0,
    vs = 5.8135195f0,
    wts = 5.4669003f0,
    x = 363318.03f0,
    y = 0.0f0,
    z = 115325.625f0,
    ztilde = 124096.734f0,
)
linf = (
    n2 = 0.0f0,
    p = 320.76422f0,
    pip = 0.0010347428f0,
    rhobar = 1.069213f0,
    rhop = 0.0058253715f0,
    t = 3600.0f0,
    thetabar = 300.0f0,
    us = 22.151293f0,
    vs = 1.9208302f0,
    wts = 1.5948994f0,
    x = 180000.0f0,
    y = 0.0f0,
    z = 19000.0f0,
    ztilde = 20000.0f0,
)
reference = (l2, linf)

@testset "Hot bubble" begin
    test_example(
        joinpath(submit_directory, "hot_bubble.jl"),
        reference,
        :x_size => 10,
        :z_size => 10,
        :output => OutputNamelist(; prepare_restart = true);
        update_references,
    )
end
