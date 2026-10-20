function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0014157434f0,
        dlr = 0.0014080313f0,
        dmr = 0.0038572003f0,
        dxr = 222763.83f0,
        dyr = 223758.31f0,
        dzr = 5601.6045f0,
        kr = 0.014050493f0,
        lr = 3.2636462f-6,
        mr = 0.03718075f0,
        n2 = 0.02856384f0,
        nr = 2.6667583f15,
        p = 25764.057f0,
        pip = 0.012910779f0,
        rhobar = 80.18411f0,
        rhop = 0.006651344f0,
        t = 300.0f0,
        thetabar = 29177.266f0,
        us = 1265.1172f0,
        uw = 3.235108f0,
        vs = 3.5866296f0,
        wts = 0.11084713f0,
        x = 257875.94f0,
        xr = 566026.3f0,
        y = 257875.94f0,
        yr = 552538.0f0,
        z = 258426.88f0,
        zr = 11310.329f0,
        ztilde = 268158.8f0,
    )
    linf = (
        dkr = 6.470947f-5,
        dlr = 6.3525426f-5,
        dmr = 0.00019786698f0,
        dxr = 10368.592f0,
        dyr = 10027.303f0,
        dzr = 300.43854f0,
        kr = 0.0006294487f0,
        lr = 4.3280124f-7,
        mr = 0.0016969579f0,
        n2 = 0.00031935345f0,
        nr = 3.918286f14,
        p = 344.90543f0,
        pip = 0.00042407282f0,
        rhobar = 1.145016f0,
        rhop = 0.00061810936f0,
        t = 300.0f0,
        thetabar = 351.64978f0,
        us = 10.27837f0,
        uw = 0.5430088f0,
        vs = 0.12918442f0,
        wts = 0.0062386347f0,
        x = 95000.0f0,
        xr = 47774.54f0,
        y = 95000.0f0,
        yr = 45007.812f0,
        z = 4879.7573f0,
        zr = 1070.2614f0,
        ztilde = 5000.0f0,
    )
    reference = (l2, linf)

    @testset "WKB mountain wave" begin
        test_example(wkb_mountain_wave, keywords, reference; update)
    end

    return
end
