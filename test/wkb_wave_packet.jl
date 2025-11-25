l2 = (
    dkr = 0.008720753f0,
    dlr = 0.008720691f0,
    dmr = 0.0061664907f0,
    dxr = 35774.8f0,
    dyr = 34668.01f0,
    dzr = 52839.38f0,
    kr = 0.0616586f0,
    lr = 0.061655745f0,
    mr = 0.061686628f0,
    n2 = 0.011484171f0,
    nr = 4.3560223f10,
    p = 5766.7773f0,
    pip = 9.799567f-7,
    rhobar = 13.219808f0,
    rhop = 9.1163585f-8,
    t = 3600.0f0,
    thetabar = 18475.896f0,
    us = 0.018543681f0,
    vs = 0.017528562f0,
    wts = 6.602516f-5,
    x = 18165.902f0,
    xr = 90936.35f0,
    y = 18165.902f0,
    yr = 88669.05f0,
    z = 729383.3f0,
    zr = 623523.4f0,
    ztilde = 784856.7f0,
)
linf = (
    dkr = 0.0003556991f0,
    dlr = 0.00035570015f0,
    dmr = 0.0002513274f0,
    dxr = 2000.0f0,
    dyr = 2000.0f0,
    dzr = 4000.0f0,
    kr = 0.0025185477f0,
    lr = 0.0025185319f0,
    mr = 0.0025169451f0,
    n2 = 0.0004251976f0,
    nr = 7.197886f9,
    p = 294.45767f0,
    pip = 4.7752987f-8,
    rhobar = 0.98152554f0,
    rhop = 7.609013f-9,
    t = 3600.0f0,
    thetabar = 1003.64667f0,
    us = 0.0045344345f0,
    vs = 0.004240775f0,
    wts = 1.5284726f-5,
    x = 9000.0f0,
    xr = 7522.246f0,
    y = 9000.0f0,
    yr = 7522.056f0,
    z = 38000.0f0,
    zr = 31040.598f0,
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
