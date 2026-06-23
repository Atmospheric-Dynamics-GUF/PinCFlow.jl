function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0011751519f0,
        dlr = 0.0011753674f0,
        dmr = 0.0012153698f0,
        dxr = 485775.03f0,
        dyr = 485779.44f0,
        dzr = 18207.76f0,
        kr = 0.011752316f0,
        lr = 1.677306f-6,
        mr = 0.028594991f0,
        n2 = 0.012520827f0,
        nr = 1.7992885f14,
        p = 5733.0596f0,
        pip = 0.0005967409f0,
        rhobar = 17.83844f0,
        rhop = 0.00033124947f0,
        t = 3600.0f0,
        thetabar = 12638.998f0,
        tke = 0.002236068f0,
        us = 447.22897f0,
        vs = 0.20153242f0,
        wts = 0.030925881f0,
        x = 363318.03f0,
        xr = 372742.75f0,
        y = 363318.03f0,
        yr = 408672.75f0,
        z = 364706.28f0,
        zr = 153482.69f0,
        ztilde = 392441.75f0,
    )
    linf = (
        dkr = 6.288036f-5,
        dlr = 6.285317f-5,
        dmr = 8.345443f-5,
        dxr = 40000.0f0,
        dyr = 40000.0f0,
        dzr = 1996.0245f0,
        kr = 0.0006288044f0,
        lr = 2.1374763f-7,
        mr = 0.0023398877f0,
        n2 = 0.00057970756f0,
        nr = 1.1041572f13,
        p = 320.91135f0,
        pip = 6.432183f-5,
        rhobar = 1.0579953f0,
        rhop = 0.00012118336f0,
        t = 3600.0f0,
        thetabar = 569.07f0,
        tke = 5.000573f-5,
        us = 10.117493f0,
        vs = 0.055419385f0,
        wts = 0.007689721f0,
        x = 180000.0f0,
        xr = 29990.668f0,
        y = 180000.0f0,
        yr = 30011.508f0,
        z = 19001.988f0,
        zr = 13223.299f0,
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
