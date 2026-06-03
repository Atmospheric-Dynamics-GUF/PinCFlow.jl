function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0017588523f0,
        dlr = 0.0017592126f0,
        dmr = 0.0024346884f0,
        dxr = 761689.94f0,
        dyr = 761510.94f0,
        dzr = 16030.756f0,
        kr = 0.017589377f0,
        lr = 3.115409f-6,
        mr = 0.03486739f0,
        n2 = 0.012520827f0,
        nr = 2.7448467f14,
        p = 5733.0596f0,
        pip = 0.00059674797f0,
        rhobar = 17.83844f0,
        rhop = 0.00033060482f0,
        t = 3600.0f0,
        thetabar = 12638.998f0,
        tke = 0.002236068f0,
        us = 447.22742f0,
        vs = 0.2006127f0,
        wts = 0.03101719f0,
        x = 363318.03f0,
        xr = 583176.44f0,
        y = 363318.03f0,
        yr = 608133.3f0,
        z = 364706.28f0,
        zr = 188553.98f0,
        ztilde = 392441.75f0,
    )
    linf = (
        dkr = 6.30285f-5,
        dlr = 6.297082f-5,
        dmr = 0.0006420181f0,
        dxr = 41264.863f0,
        dyr = 40015.688f0,
        dzr = 2883.5125f0,
        kr = 0.00062886043f0,
        lr = 2.8110185f-7,
        mr = 0.0023430092f0,
        n2 = 0.00057970756f0,
        nr = 1.1042542f13,
        p = 320.91135f0,
        pip = 6.432183f-5,
        rhobar = 1.0579953f0,
        rhop = 0.00012087787f0,
        t = 3600.0f0,
        thetabar = 569.07f0,
        tke = 5.0005714f-5,
        us = 10.119257f0,
        vs = 0.05485646f0,
        wts = 0.0076800627f0,
        x = 180000.0f0,
        xr = 29975.04f0,
        y = 180000.0f0,
        yr = 30014.275f0,
        z = 19001.988f0,
        zr = 13394.341f0,
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
