l2 = (
    n2 = 0.0031957547f0,
    p = 1888.8247f0,
    pip = 1.3406276f-5,
    rhobar = 5.38609f0,
    rhop = 0.000343484f0,
    t = 3600.0f0,
    thetabar = 4300.3516f0,
    us = 141.42049f0,
    vs = 0.0f0,
    wts = 0.13482842f0,
    x = 18165.902f0,
    y = 0.0f0,
    z = 115340.16f0,
    ztilde = 124110.04f0,
)
linf = (
    n2 = 0.00031957557f0,
    p = 321.19247f0,
    pip = 1.8913985f-6,
    rhobar = 1.0363418f0,
    rhop = 7.596874f-5,
    t = 3600.0f0,
    thetabar = 556.86066f0,
    us = 10.053132f0,
    vs = 0.0f0,
    wts = 0.02500533f0,
    x = 9000.0f0,
    y = 0.0f0,
    z = 19000.488f0,
    ztilde = 20000.0f0,
)
reference = (l2, linf)

@testset "periodic_hill" begin
    test_example(
        joinpath(submit_directory, "periodic_hill.jl"),
        reference,
        "x_size = 10",
        "z_size = 10",
        "output_variables = ()",
        "prepare_restart = true",
    )
end
