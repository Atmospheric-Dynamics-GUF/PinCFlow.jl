function test_wave_packet()
    l2 = (
        n2 = 0.011484171f0,
        p = 3961.7478f0,
        pip = 1.9827688f-5,
        rhobar = 12.988534f0,
        rhop = 0.00011797931f0,
        t = 900.0f0,
        thetabar = 18475.896f0,
        tke = 0.0022360685f0,
        u = 0.1385172f0,
        us = 0.15143847f0,
        v = 0.13760449f0,
        vs = 0.15059878f0,
        w = 0.1974145f0,
        wts = 0.2463817f0,
        x = 18165.902f0,
        y = 18165.902f0,
        z = 729383.3f0,
        ztilde = 784856.7f0,
    )
    linf = (
        n2 = 0.0004251976f0,
        p = 294.45767f0,
        pip = 7.2657526f-6,
        rhobar = 0.98152554f0,
        rhop = 5.7667487f-5,
        t = 900.0f0,
        thetabar = 1003.64667f0,
        tke = 5.001521f-5,
        u = 0.037085615f0,
        us = 0.043309957f0,
        v = 0.036703926f0,
        vs = 0.04361742f0,
        w = 0.051753767f0,
        wts = 0.072462775f0,
        x = 9000.0f0,
        y = 9000.0f0,
        z = 38000.0f0,
        ztilde = 40000.0f0,
    )
    reference = (l2, linf)

    keywords = (
        x_size = 10,
        y_size = 10,
        z_size = 10,
        prepare_restart = true,
        visualize = false,
    )

    @testset "Wave packet" begin
        test_example(wave_packet, keywords, reference; update_references)
    end

    return
end
