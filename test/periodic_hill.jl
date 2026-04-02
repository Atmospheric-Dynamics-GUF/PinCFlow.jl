l2 = (
    buoyancyproduction = 0.00021221871f0,
    kek = 2.1212764f0,
    kh = 2.1212764f0,
    km = 2.1212764f0,
    pip = 0.00062133715f0,
    rhop = 0.008697031f0,
    shearproduction = 2.202516f-7,
    t = 3600.0f0,
    tke = 0.0007070922f0,
    us = 140.49748f0,
    vs = 0.0f0,
    wts = 4.2872596f0,
    x = 18165.902f0,
    y = 0.0f0,
    z = 116062.984f0,
    ztilde = 124774.234f0,
)
linf = (
    buoyancyproduction = 2.3234028f-5,
    kek = 0.2123454f0,
    kh = 0.2123454f0,
    km = 0.2123454f0,
    pip = 0.0002016583f0,
    rhop = 0.0023945163f0,
    shearproduction = 1.009385f-7,
    t = 3600.0f0,
    tke = 5.0101324f-5,
    us = 10.880933f0,
    vs = 0.0f0,
    wts = 1.1136314f0,
    x = 9000.0f0,
    y = 0.0f0,
    z = 19024.389f0,
    ztilde = 20000.0f0,
)
reference = (l2, linf)

@testset "Periodic hill" begin
    test_example(
        joinpath(scripts_directory, "periodic_hill.jl"),
        reference,
        :x_size => 10,
        :z_size => 10,
        :output => OutputNamelist(; prepare_restart = true);
        update_references,
    )
end
