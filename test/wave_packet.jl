l2 = (
    n2 = 0.011484171f0,
    p = 4077.7273f0,
    pip = 1.9591602f-5,
    rhobar = 13.219808f0,
    rhop = 0.00014683898f0,
    t = 3600.0f0,
    thetabar = 18475.896f0,
    us = 0.15009244f0,
    vs = 0.14708927f0,
    wts = 0.25232828f0,
    x = 18165.902f0,
    y = 18165.902f0,
    z = 729383.3f0,
    ztilde = 784856.7f0,
)
linf = (
    n2 = 0.0004251976f0,
    p = 294.45767f0,
    pip = 7.2291814f-6,
    rhobar = 0.98152554f0,
    rhop = 7.640168f-5,
    t = 3600.0f0,
    thetabar = 1003.64667f0,
    us = 0.043309662f0,
    vs = 0.043617133f0,
    wts = 0.07246128f0,
    x = 9000.0f0,
    y = 9000.0f0,
    z = 38000.0f0,
    ztilde = 40000.0f0,
)
reference = (l2, linf)

@testset "Wave packet" begin
    cp(
        joinpath(scripts_directory, "wave_packet_tools.jl"),
        "wave_packet_tools.jl",
    )
    test_example(
        joinpath(scripts_directory, "wave_packet.jl"),
        reference,
        :x_size => 10,
        :y_size => 10,
        :z_size => 10,
        :output => OutputNamelist(; prepare_restart = true);
        update_references,
    )
    rm("wave_packet_tools.jl")
end
