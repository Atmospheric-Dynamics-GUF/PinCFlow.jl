l2 = (
    pip = 0.00055369554f0,
    rhop = 0.005398585f0,
    t = 3600.0f0,
    us = 140.5774f0,
    vs = 0.0f0,
    wts = 4.074997f0,
    x = 18165.902f0,
    y = 0.0f0,
    z = 116062.984f0,
    ztilde = 124774.234f0,
)
linf = (
    pip = 0.0002016583f0,
    rhop = 0.0016297717f0,
    t = 3600.0f0,
    us = 10.618719f0,
    vs = 0.0f0,
    wts = 0.98847735f0,
    x = 9000.0f0,
    y = 0.0f0,
    z = 19024.389f0,
    ztilde = 20000.0f0,
)
reference = (l2, linf)

@testset "Periodic hill" begin
    test_example(
        joinpath(submit_directory, "periodic_hill.jl"),
        reference,
        :x_size => 10,
        :z_size => 10,
        :output => OutputNamelist(; prepare_restart = true);
        update_references,
    )
end
