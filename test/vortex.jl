l2 = (
    buoyancyproduction = 0.0f0,
    chi = 3.6041312f0,
    kek = 2.1213202f0,
    kh = 2.1213202f0,
    km = 2.1213202f0,
    pip = 3.849339f-5,
    rhop = 0.0f0,
    shearproduction = 0.0f0,
    t = 3600.0f0,
    tke = 0.0007071068f0,
    us = 2.9554124f0,
    vs = 2.9554124f0,
    wts = 0.0f0,
    x = 18165.902f0,
    y = 18165.902f0,
    z = 5000.0f0,
    ztilde = 10000.0f0,
)
linf = (
    buoyancyproduction = 0.0f0,
    chi = 0.9045085f0,
    kek = 0.21213204f0,
    kh = 0.21213204f0,
    km = 0.21213204f0,
    pip = 9.473251f-6,
    rhop = 0.0f0,
    shearproduction = 0.0f0,
    t = 3600.0f0,
    tke = 5.0f-5,
    us = 0.9283495f0,
    vs = 0.9283495f0,
    wts = 0.0f0,
    x = 9000.0f0,
    y = 9000.0f0,
    z = 500.0f0,
    ztilde = 1000.0f0,
)
reference = (l2, linf)

@testset "Vortex" begin
    test_example(
        joinpath(scripts_directory, "vortex.jl"),
        reference,
        :x_size => 10,
        :y_size => 10,
        :output => OutputNamelist(; prepare_restart = true);
        update_references,
    )
end
