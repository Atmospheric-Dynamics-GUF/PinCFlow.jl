function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0013637786f0,
        dlr = 0.0013628941f0,
        dmr = 0.0014700559f0,
        dxr = 678107.44f0,
        dyr = 675763.8f0,
        dzr = 23403.05f0,
        kr = 0.013622303f0,
        lr = 1.9720014f-6,
        mr = 0.033773884f0,
        n2 = 0.012520827f0,
        nr = 2.3978588f14,
        p = 5733.0596f0,
        pip = 0.00059674704f0,
        rhobar = 17.83844f0,
        rhop = 0.0003305382f0,
        t = 3600.0f0,
        thetabar = 12638.998f0,
        tke = 0.002236068f0,
        us = 447.22742f0,
        vs = 0.20060988f0,
        wts = 0.03100347f0,
        x = 363318.03f0,
        xr = 470110.94f0,
        y = 363318.03f0,
        yr = 461071.66f0,
        z = 364706.28f0,
        zr = 181062.89f0,
        ztilde = 392441.75f0,
    )
    linf = (
        dkr = 6.306611f-5,
        dlr = 6.2921834f-5,
        dmr = 9.414942f-5,
        dxr = 40000.0f0,
        dyr = 40000.0f0,
        dzr = 1996.0245f0,
        kr = 0.0006288618f0,
        lr = 2.0336576f-7,
        mr = 0.0024997713f0,
        n2 = 0.00057970756f0,
        nr = 1.5666467f13,
        p = 320.91135f0,
        pip = 6.432183f-5,
        rhobar = 1.0579953f0,
        rhop = 0.00012086283f0,
        t = 3600.0f0,
        thetabar = 569.07f0,
        tke = 5.0005714f-5,
        us = 10.119353f0,
        vs = 0.054848224f0,
        wts = 0.0076817446f0,
        x = 180000.0f0,
        xr = 29920.406f0,
        y = 180000.0f0,
        yr = 30015.31f0,
        z = 19001.988f0,
        zr = 14204.129f0,
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
