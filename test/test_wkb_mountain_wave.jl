function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0011751519f0,
        dlr = 0.0011753674f0,
        dmr = 0.0012153187f0,
        dxr = 485775.0f0,
        dyr = 485779.44f0,
        dzr = 18207.799f0,
        kr = 0.011752316f0,
        lr = 1.6774112f-6,
        mr = 0.028595012f0,
        n2 = 0.012520827f0,
        nr = 1.7993104f14,
        p = 5733.0596f0,
        pip = 0.00059674f0,
        rhobar = 17.83844f0,
        rhop = 0.00033116358f0,
        t = 3600.0f0,
        thetabar = 12638.998f0,
        tke = 0.002236068f0,
        us = 447.22928f0,
        vs = 0.20151146f0,
        wts = 0.030926727f0,
        x = 363318.03f0,
        xr = 372741.78f0,
        y = 363318.03f0,
        yr = 408672.75f0,
        z = 364706.28f0,
        zr = 153482.83f0,
        ztilde = 392441.75f0,
    )
    linf = (
        dkr = 6.288036f-5,
        dlr = 6.285317f-5,
        dmr = 8.345443f-5,
        dxr = 40000.0f0,
        dyr = 40000.0f0,
        dzr = 1996.0245f0,
        kr = 0.00062880444f0,
        lr = 2.1376468f-7,
        mr = 0.002339892f0,
        n2 = 0.00057970756f0,
        nr = 1.1041821f13,
        p = 320.91135f0,
        pip = 6.432183f-5,
        rhobar = 1.0579953f0,
        rhop = 0.000121180805f0,
        t = 3600.0f0,
        thetabar = 569.07f0,
        tke = 5.000573f-5,
        us = 10.117885f0,
        vs = 0.055415075f0,
        wts = 0.007690605f0,
        x = 180000.0f0,
        xr = 29990.355f0,
        y = 180000.0f0,
        yr = 30011.508f0,
        z = 19001.988f0,
        zr = 13223.245f0,
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
        test_example(wkb_mountain_wave, keywords, reference; update)
    end

    return
end
