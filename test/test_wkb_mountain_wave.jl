function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0018781909f0,
        dlr = 0.0018785506f0,
        dmr = 0.0025569731f0,
        dxr = 788494.2f0,
        dyr = 788364.25f0,
        dzr = 18042.904f0,
        kr = 0.018782577f0,
        lr = 3.1936274f-6,
        mr = 0.041188624f0,
        n2 = 0.012520827f0,
        nr = 2.9261713f14,
        p = 5733.0596f0,
        pip = 0.0005967469f0,
        rhobar = 17.83844f0,
        rhop = 0.00033051078f0,
        t = 3600.0f0,
        thetabar = 12638.998f0,
        us = 447.22772f0,
        vs = 0.20059822f0,
        wts = 0.03101894f0,
        x = 363318.03f0,
        xr = 606065.4f0,
        y = 363318.03f0,
        yr = 651951.75f0,
        z = 364706.28f0,
        zr = 220968.7f0,
        ztilde = 392441.75f0,
    )
    linf = (
        dkr = 6.300326f-5,
        dlr = 6.292205f-5,
        dmr = 0.00039536806f0,
        dxr = 40000.0f0,
        dyr = 40000.0f0,
        dzr = 1996.0245f0,
        kr = 0.0006288605f0,
        lr = 2.8112373f-7,
        mr = 0.0023443382f0,
        n2 = 0.00057970756f0,
        nr = 1.1042815f13,
        p = 320.91135f0,
        pip = 6.432183f-5,
        rhobar = 1.0579953f0,
        rhop = 0.00012087278f0,
        t = 3600.0f0,
        thetabar = 569.07f0,
        us = 10.119688f0,
        vs = 0.05485425f0,
        wts = 0.007681421f0,
        x = 180000.0f0,
        xr = 29974.87f0,
        y = 180000.0f0,
        yr = 30014.275f0,
        z = 19001.988f0,
        zr = 13316.297f0,
        ztilde = 20000.0f0,
    )
    reference = (l2, linf)

    keywords = (
        x_size = 10,
        y_size = 10,
        z_size = 10,
        prepare_restart = true,
        visualize = false,
    )

    @testset "WKB mountain wave" begin
        test_example(wkb_mountain_wave, keywords, reference; update_references)
    end

    return
end
