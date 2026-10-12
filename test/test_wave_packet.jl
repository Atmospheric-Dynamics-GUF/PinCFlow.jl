function test_wave_packet()
    l2 = (
        n2 = 0.03744983f0,
        p = 5563.24f0,
        pip = 5.3333635f-5,
        rhobar = 15.7452135f0,
        rhop = 0.0009089816f0,
        t = 600.0f0,
        thetabar = 43872.85f0,
        tke = 0.0063245553f0,
        u = 1.3253986f0,
        us = 1.6192237f0,
        v = 1.325949f0,
        vs = 1.620395f0,
        w = 2.736039f0,
        wts = 4.841161f0,
        x = 25787.594f0,
        y = 25787.594f0,
        z = 1.8617196f6,
        ztilde = 1.9152024f6,
    )
    linf = (
        n2 = 0.00042322697f0,
        p = 123.33977f0,
        pip = 9.576407f-6,
        rhobar = 0.4023616f0,
        rhop = 0.00015550917f0,
        t = 600.0f0,
        thetabar = 695.6166f0,
        tke = 5.003752f-5,
        u = 0.19352245f0,
        us = 0.24351631f0,
        v = 0.19409175f0,
        vs = 0.24438018f0,
        w = 0.38159657f0,
        wts = 0.7121469f0,
        x = 9500.0f0,
        y = 9500.0f0,
        z = 29500.0f0,
        ztilde = 30000.0f0,
    )
    reference = (l2, linf)

    @testset "Wave packet" begin
        test_example(wave_packet, keywords, reference; update)
    end

    return
end
