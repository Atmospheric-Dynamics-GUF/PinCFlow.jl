function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0017247543f0,
        dlr = 0.0017251098f0,
        dmr = 0.0017544362f0,
        dxr = 704247.75f0,
        dyr = 704247.25f0,
        dzr = 25139.852f0,
        kr = 0.017248727f0,
        lr = 2.711103f-6,
        mr = 0.043083366f0,
        n2 = 0.012520827f0,
        nr = 2.4701057f14,
        p = 5733.0596f0,
        pip = 0.0005967452f0,
        rhobar = 17.83844f0,
        rhop = 0.0003305276f0,
        t = 3600.0f0,
        thetabar = 12638.998f0,
        us = 447.22766f0,
        vs = 0.19881797f0,
        wts = 0.03097478f0,
        x = 363318.03f0,
        xr = 539040.94f0,
        y = 363318.03f0,
        yr = 600703.7f0,
        z = 364706.28f0,
        zr = 231974.25f0,
        ztilde = 392441.75f0,
    )
    linf = (
        dkr = 6.2885956f-5,
        dlr = 6.285545f-5,
        dmr = 8.345443f-5,
        dxr = 40000.0f0,
        dyr = 40000.0f0,
        dzr = 1996.0245f0,
        kr = 0.0006288603f0,
        lr = 2.7030126f-7,
        mr = 0.0023649868f0,
        n2 = 0.00057970756f0,
        nr = 1.104297f13,
        p = 320.91135f0,
        pip = 6.432183f-5,
        rhobar = 1.0579953f0,
        rhop = 0.000120980316f0,
        t = 3600.0f0,
        thetabar = 569.07f0,
        us = 10.120789f0,
        vs = 0.054217454f0,
        wts = 0.007657155f0,
        x = 180000.0f0,
        xr = 29828.967f0,
        y = 180000.0f0,
        yr = 30014.717f0,
        z = 19001.988f0,
        zr = 13710.953f0,
        ztilde = 20000.0f0,
    )
    reference = (l2, linf)

    keywords = (
        x_size = 10,
        y_size = 10,
        z_size = 10,
        npx = 1,
        npy = 1,
        npz = 1,
        output = OutputNamelist(; prepare_restart = true),
        visualize = false,
    )

    @testset "WKB mountain wave" begin
        test_example(wkb_mountain_wave, keywords, reference; update_references)
    end

    return
end
