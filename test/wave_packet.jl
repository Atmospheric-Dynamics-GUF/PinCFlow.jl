l2 = (
    buoyancyproduction = 0.0024361617f0,
    kek = 6.7082057f0,
    kh = 6.7082057f0,
    km = 6.7082057f0,
    n2 = 0.011484171f0,
    p = 4077.7273f0,
    pip = 1.9592067f-5,
    rhobar = 13.219808f0,
    rhop = 0.00014683888f0,
    shearproduction = 1.2759932f-11,
    t = 3600.0f0,
    thetabar = 18475.896f0,
    tke = 0.0022360685f0,
    us = 0.15008229f0,
    vs = 0.1470851f0,
    wts = 0.25232783f0,
    x = 18165.902f0,
    y = 18165.902f0,
    z = 729383.3f0,
    ztilde = 784856.7f0,
)
linf = (
    buoyancyproduction = 9.0260655f-5,
    kek = 0.21215414f0,
    kh = 0.21215414f0,
    km = 0.21215414f0,
    n2 = 0.0004251976f0,
    p = 294.45767f0,
    pip = 7.229181f-6,
    rhobar = 0.98152554f0,
    rhop = 7.640168f-5,
    shearproduction = 3.2465233f-12,
    t = 3600.0f0,
    thetabar = 1003.64667f0,
    tke = 5.0010898f-5,
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
        "wave_packet_tools.jl";
        force = true,
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
