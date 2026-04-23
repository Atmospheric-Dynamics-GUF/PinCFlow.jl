l2 = (
    dkr = 0.0093903635f0,
    dlr = 0.009390376f0,
    dmr = 0.0066399924f0,
    dxr = 34234.277f0,
    dyr = 34539.58f0,
    dzr = 56356.01f0,
    kr = 0.066400856f0,
    lr = 0.06640136f0,
    mr = 0.06639635f0,
    n2 = 0.011484171f0,
    nr = 4.60334f10,
    p = 5766.73f0,
    pip = 5.167412f-5,
    rhobar = 13.219808f0,
    rhop = 1.735658f-6,
    t = 3600.0f0,
    thetabar = 18475.896f0,
    tke = 0.0022360676f0,
    us = 0.006486799f0,
    vs = 0.006070223f0,
    wts = 0.00051546056f0,
    x = 18165.902f0,
    xr = 99434.016f0,
    y = 18165.902f0,
    yr = 96132.88f0,
    z = 729383.3f0,
    zr = 678006.44f0,
    ztilde = 784856.7f0,
)
linf = (
    dkr = 0.0003554692f0,
    dlr = 0.00035546927f0,
    dmr = 0.0002513274f0,
    dxr = 2000.0f0,
    dyr = 2000.0f0,
    dzr = 4000.0f0,
    kr = 0.002513892f0,
    lr = 0.0025139165f0,
    mr = 0.0025135484f0,
    n2 = 0.0004251976f0,
    nr = 7.197886f9,
    p = 294.45767f0,
    pip = 2.7410977f-6,
    rhobar = 0.98152554f0,
    rhop = 1.5865861f-7,
    t = 3600.0f0,
    thetabar = 1003.64667f0,
    tke = 5.0000013f-5,
    us = 0.0018015321f0,
    vs = 0.0016603299f0,
    wts = 9.308403f-5,
    x = 9000.0f0,
    xr = 7519.3354f0,
    y = 9000.0f0,
    yr = 7519.3037f0,
    z = 38000.0f0,
    zr = 31038.988f0,
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
