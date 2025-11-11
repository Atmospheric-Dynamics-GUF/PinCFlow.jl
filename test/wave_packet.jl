l2 = (
    n2 = 0.009136622f0,
    p = 5837.2305f0,
    pip = 4.0536876f0,
    rhobar = 19.042942f0,
    rhop = 0.3046972f0,
    t = 3600.0f0,
    thetabar = 10786.326f0,
    us = 49.6171f0,
    vs = 34.40523f0,
    wts = 64.71526f0,
    x = 18165.902f0,
    y = 18165.902f0,
    z = 364691.66f0,
    ztilde = 392428.34f0,
)
linf = (
    n2 = 0.00042362066f0,
    p = 320.76392f0,
    pip = 0.84784794f0,
    rhobar = 1.069213f0,
    rhop = 0.06671554f0,
    t = 3600.0f0,
    thetabar = 442.28027f0,
    us = 5.9756207f0,
    vs = 5.5485463f0,
    wts = 11.815212f0,
    x = 9000.0f0,
    y = 9000.0f0,
    z = 19000.0f0,
    ztilde = 20000.0f0,
)
reference = (l2, linf)

@testset "Wave packet" begin
    test_example(
        joinpath(submit_directory, "wave_packet.jl"),
        reference,
        :x_size => 10,
        :y_size => 10,
        :z_size => 10,
        :output => OutputNamelist(; prepare_restart = true);
        update_references,
    )
end
