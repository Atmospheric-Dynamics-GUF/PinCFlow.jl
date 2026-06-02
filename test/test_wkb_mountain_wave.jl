function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0018781829f0,
        dlr = 0.0018785471f0,
        dmr = 0.0025567617f0,
        dxr = 788496.2f0,
        dyr = 788365.1f0,
        dzr = 18050.56f0,
        kr = 0.0187825f0,
        lr = 3.1978634f-6,
        mr = 0.04119062f0,
        n2 = 0.012520827f0,
        nr = 2.9261314f14,
        p = 5733.0596f0,
        pip = 0.00059674797f0,
        rhobar = 17.83844f0,
        rhop = 0.0003306045f0,
        t = 3600.0f0,
        thetabar = 12638.998f0,
        tke = 0.002236068f0,
        us = 447.22742f0,
        vs = 0.20061252f0,
        wts = 0.031017218f0,
        x = 363318.03f0,
        xr = 608370.9f0,
        y = 363318.03f0,
        yr = 651953.06f0,
        z = 364706.28f0,
        zr = 221021.48f0,
        ztilde = 392441.75f0,
    )
    linf = (
        dkr = 6.300329f-5,
        dlr = 6.292205f-5,
        dmr = 0.00039532495f0,
        dxr = 40000.0f0,
        dyr = 40000.0f0,
        dzr = 1996.0245f0,
        kr = 0.00062886043f0,
        lr = 2.8110222f-7,
        mr = 0.002344334f0,
        n2 = 0.00057970756f0,
        nr = 1.1042542f13,
        p = 320.91135f0,
        pip = 6.432183f-5,
        rhobar = 1.0579953f0,
        rhop = 0.000120877856f0,
        t = 3600.0f0,
        thetabar = 569.07f0,
        tke = 5.0005714f-5,
        us = 10.119257f0,
        vs = 0.05485631f0,
        wts = 0.0076800436f0,
        x = 180000.0f0,
        xr = 29975.04f0,
        y = 180000.0f0,
        yr = 30014.275f0,
        z = 19001.988f0,
        zr = 13316.349f0,
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
