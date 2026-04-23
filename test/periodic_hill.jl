l2 = (
    pip = 0.00062332454f0,
    rhop = 0.0086944355f0,
    t = 3600.0f0,
    tke = 0.0007071068f0,
    us = 140.49619f0,
    vs = 0.0f0,
    wts = 4.285317f0,
    x = 18165.902f0,
    y = 0.0f0,
    z = 116062.984f0,
    ztilde = 124774.234f0,
)
linf = (
    pip = 0.0002016583f0,
    rhop = 0.0023938704f0,
    t = 3600.0f0,
    tke = 5.0f-5,
    us = 10.87999f0,
    vs = 0.0f0,
    wts = 1.1131233f0,
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
