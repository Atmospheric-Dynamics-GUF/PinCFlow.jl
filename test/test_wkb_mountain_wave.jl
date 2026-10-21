function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.00140488f0,
        dlr = 0.001398942f0,
        dmr = 0.0037843657f0,
        dxr = 221713.22f0,
        dyr = 222395.48f0,
        dzr = 5569.3833f0,
        kr = 0.013965767f0,
        lr = 3.1620802f-6,
        mr = 0.036955386f0,
        n2 = 0.040395368f0,
        nr = 2.7704423f15,
        p = 36435.88f0,
        pip = 0.012910796f0,
        rhobar = 113.39746f0,
        rhop = 0.006651708f0,
        t = 300.0f0,
        thetabar = 41262.887f0,
        us = 1265.1171f0,
        uw = 3.2433107f0,
        vs = 3.5869687f0,
        wts = 0.11078721f0,
        x = 257875.94f0,
        xr = 564839.9f0,
        y = 257875.94f0,
        yr = 552401.5f0,
        z = 365470.78f0,
        zr = 11380.941f0,
        ztilde = 379233.8f0,
    )
    linf = (
        dkr = 6.450755f-5,
        dlr = 6.340685f-5,
        dmr = 0.00019519957f0,
        dxr = 10354.15f0,
        dyr = 10026.717f0,
        dzr = 431.71576f0,
        kr = 0.00062975334f0,
        lr = 4.505031f-7,
        mr = 0.0016969205f0,
        n2 = 0.00031935345f0,
        nr = 3.9182366f14,
        p = 344.90543f0,
        pip = 0.00042407282f0,
        rhobar = 1.145016f0,
        rhop = 0.0006181189f0,
        t = 300.0f0,
        thetabar = 351.64978f0,
        us = 10.2784605f0,
        uw = 0.549437f0,
        vs = 0.1291953f0,
        wts = 0.0062398165f0,
        x = 95000.0f0,
        xr = 47774.977f0,
        y = 95000.0f0,
        yr = 45007.86f0,
        z = 4879.7573f0,
        zr = 1107.1731f0,
        ztilde = 5000.0f0,
    )
    reference = (l2, linf)

    @testset "WKB mountain wave" begin
        test_example(wkb_mountain_wave, keywords, reference; update)
    end

    return
end
