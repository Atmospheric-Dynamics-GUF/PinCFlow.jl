l2 = (
    n2 = 0.0f0,
    p = 2484.3005f0,
    pip = 0.0017557857f0,
    rhobar = 5.855689f0,
    rhop = 0.008923059f0,
    t = 3600.0f0,
    thetabar = 3000.0f0,
    tke = 0.019977719f0,
    us = 22.346428f0,
    vs = 0.8279369f0,
    wts = 25.47495f0,
    x = 18165.902f0,
    y = 0.0f0,
    z = 115325.625f0,
    ztilde = 124096.734f0,
)
linf = (
    n2 = 0.0f0,
    p = 320.76392f0,
    pip = 0.0005637321f0,
    rhobar = 1.069213f0,
    rhop = 0.0058253715f0,
    t = 3600.0f0,
    thetabar = 300.0f0,
    tke = 0.007874703f0,
    us = 6.7488136f0,
    vs = 0.29772216f0,
    wts = 6.1135535f0,
    x = 9000.0f0,
    y = 0.0f0,
    z = 19000.0f0,
    ztilde = 20000.0f0,
)
reference = (l2, linf)

@testset "Hot bubble" begin
    test_example(
        joinpath(scripts_directory, "hot_bubble.jl"),
        reference,
        :x_size => 10,
        :z_size => 10,
        :output => OutputNamelist(; prepare_restart = true);
        update_references,
    )
end
