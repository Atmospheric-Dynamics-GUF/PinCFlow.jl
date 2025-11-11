l2 = (
    n2 = 0.0f0,
    p = 2484.3528f0,
    pip = 0.0018367001f0,
    rhobar = 5.855689f0,
    rhop = 0.00893393f0,
    t = 3600.0f0,
    thetabar = 3000.0f0,
    us = 20.658684f0,
    vs = 0.5071614f0,
    wts = 23.970451f0,
    x = 18165.902f0,
    y = 0.0f0,
    z = 115325.625f0,
    ztilde = 124096.734f0,
)
linf = (
    n2 = 0.0f0,
    p = 320.76392f0,
    pip = 0.0005893867f0,
    rhobar = 1.069213f0,
    rhop = 0.0058253715f0,
    t = 3600.0f0,
    thetabar = 300.0f0,
    us = 5.9817967f0,
    vs = 0.17109516f0,
    wts = 5.1726527f0,
    x = 9000.0f0,
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
