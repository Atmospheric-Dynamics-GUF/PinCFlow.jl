function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0018782009f0,
        dlr = 0.0018785988f0,
        dmr = 0.0025586833f0,
        dxr = 788487.6f0,
        dyr = 788364.75f0,
        dzr = 18072.197f0,
        kr = 0.018782482f0,
        lr = 3.9171714f-6,
        mr = 0.041188233f0,
        n2 = 0.012520827f0,
        nr = 2.9218515f14,
        p = 5733.0596f0,
        pip = 0.00059673685f0,
        rhobar = 17.83844f0,
        rhop = 0.0003344956f0,
        t = 3600.0f0,
        thetabar = 12638.998f0,
        tke = 0.002236068f0,
        us = 447.2314f0,
        vs = 0.19802675f0,
        wts = 0.030618837f0,
        x = 363318.03f0,
        xr = 608365.0f0,
        y = 363318.03f0,
        yr = 651938.6f0,
        z = 364706.28f0,
        zr = 221077.8f0,
        ztilde = 392441.75f0,
    )
    linf = (
        dkr = 6.304347f-5,
        dlr = 6.2991516f-5,
        dmr = 0.00039608765f0,
        dxr = 40000.0f0,
        dyr = 40000.0f0,
        dzr = 1996.0245f0,
        kr = 0.0006288686f0,
        lr = 2.9519356f-7,
        mr = 0.0023464868f0,
        n2 = 0.00057970756f0,
        nr = 1.1042474f13,
        p = 320.91135f0,
        pip = 6.432183f-5,
        rhobar = 1.0579953f0,
        rhop = 0.00012305923f0,
        t = 3600.0f0,
        thetabar = 569.07f0,
        tke = 5.0005816f-5,
        us = 10.117606f0,
        vs = 0.054088634f0,
        wts = 0.0076158717f0,
        x = 180000.0f0,
        xr = 29975.475f0,
        y = 180000.0f0,
        yr = 30014.271f0,
        z = 19001.988f0,
        zr = 13313.045f0,
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
