function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0014157423f0,
        dlr = 0.0014080118f0,
        dmr = 0.0038459997f0,
        dxr = 222765.27f0,
        dyr = 223758.55f0,
        dzr = 5601.7344f0,
        kr = 0.014050498f0,
        lr = 3.2488942f-6,
        mr = 0.03718997f0,
        n2 = 0.02856384f0,
        nr = 2.6640047f15,
        p = 25764.057f0,
        pip = 0.01291078f0,
        rhobar = 80.18411f0,
        rhop = 0.006651466f0,
        t = 300.0f0,
        thetabar = 29177.266f0,
        us = 1265.1171f0,
        uw = 3.239377f0,
        vs = 3.5867424f0,
        wts = 0.11082995f0,
        x = 257875.94f0,
        xr = 566030.2f0,
        y = 257875.94f0,
        yr = 552537.9f0,
        z = 258426.88f0,
        zr = 11302.54f0,
        ztilde = 268158.8f0,
    )
    linf = (
        dkr = 6.470928f-5,
        dlr = 6.35242f-5,
        dmr = 0.00019871292f0,
        dxr = 10367.958f0,
        dyr = 10027.305f0,
        dzr = 301.23404f0,
        kr = 0.0006294488f0,
        lr = 4.32204f-7,
        mr = 0.0016969448f0,
        n2 = 0.00031935345f0,
        nr = 3.9182795f14,
        p = 344.90543f0,
        pip = 0.00042407282f0,
        rhobar = 1.145016f0,
        rhop = 0.00061810896f0,
        t = 300.0f0,
        thetabar = 351.64978f0,
        us = 10.278389f0,
        uw = 0.5435354f0,
        vs = 0.12919033f0,
        wts = 0.0062374487f0,
        x = 95000.0f0,
        xr = 47775.195f0,
        y = 95000.0f0,
        yr = 45007.816f0,
        z = 4879.7573f0,
        zr = 1070.9567f0,
        ztilde = 5000.0f0,
    )
    reference = (l2, linf)

    @testset "WKB mountain wave" begin
        test_example(wkb_mountain_wave, keywords, reference; update)
    end

    return
end
