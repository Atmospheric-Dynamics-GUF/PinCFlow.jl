l2 = (
    chi = 3.642599f0,
    pip = 0.0076706754f0,
    rhop = 0.0f0,
    t = 3600.0f0,
    us = 30.65904f0,
    vs = 30.659004f0,
    wts = 0.0f0,
    x = 363318.03f0,
    y = 363318.03f0,
    z = 5000.0f0,
    ztilde = 10000.0f0,
)
linf = (
    chi = 0.9045085f0,
    pip = 0.0018946502f0,
    rhop = 0.0f0,
    t = 3600.0f0,
    us = 9.283495f0,
    vs = 9.283495f0,
    wts = 0.0f0,
    x = 180000.0f0,
    y = 180000.0f0,
    z = 500.0f0,
    ztilde = 1000.0f0,
)
reference = (l2, linf)

@testset "Periodic hill" begin
    test_example(
        joinpath(submit_directory, "vortex.jl"),
        reference,
        :x_size => 10,
        :y_size => 10,
        :output => OutputNamelist(; prepare_restart = true);
        update_references,
    )
end
