l2 = (
    dkr = 0.008907088f0,
    dlr = 0.008907009f0,
    dmr = 0.006298247f0,
    dxr = 34610.06f0,
    dyr = 33954.152f0,
    dzr = 53814.496f0,
    kr = 0.06297686f0,
    lr = 0.062975086f0,
    mr = 0.063003115f0,
    n2 = 0.011484171f0,
    nr = 4.4345897f10,
    p = 5766.7773f0,
    pip = 9.787384f-7,
    rhobar = 13.219808f0,
    rhop = 8.9572374f-8,
    t = 3600.0f0,
    thetabar = 18475.896f0,
    us = 0.016690716f0,
    vs = 0.015832756f0,
    wts = 5.3174186f-5,
    x = 18165.902f0,
    xr = 89528.06f0,
    y = 18165.902f0,
    yr = 88949.2f0,
    z = 729383.3f0,
    zr = 624391.0f0,
    ztilde = 784856.7f0,
)
linf = (
    dkr = 0.00035570745f0,
    dlr = 0.0003557073f0,
    dmr = 0.0002513274f0,
    dxr = 2000.0f0,
    dyr = 2000.0f0,
    dzr = 4000.0f0,
    kr = 0.0025186455f0,
    lr = 0.0025186369f0,
    mr = 0.002517031f0,
    n2 = 0.0004251976f0,
    nr = 7.197886f9,
    p = 294.45767f0,
    pip = 4.9705655f-8,
    rhobar = 0.98152554f0,
    rhop = 7.2660042f-9,
    t = 3600.0f0,
    thetabar = 1003.64667f0,
    us = 0.0039137662f0,
    vs = 0.0037157412f0,
    wts = 1.2830066f-5,
    x = 9000.0f0,
    xr = 7522.5044f0,
    y = 9000.0f0,
    yr = 7522.472f0,
    z = 38000.0f0,
    zr = 31040.926f0,
    ztilde = 40000.0f0,
)
reference = (l2, linf)

@testset "WKB Wave packet" begin
    cp(
        joinpath(scripts_directory, "wave_packet_tools.jl"),
        "wave_packet_tools.jl";
        force = true,
    )
    test_example(
        joinpath(scripts_directory, "wkb_wave_packet.jl"),
        reference,
        :x_size => 10,
        :y_size => 10,
        :z_size => 10,
        :output => OutputNamelist(; prepare_restart = true);
        update_references,
    )
    rm("wave_packet_tools.jl")
end
