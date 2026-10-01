function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0014157484f0,
        dlr = 0.0014080136f0,
        dmr = 0.0038586515f0,
        dxr = 222766.36f0,
        dyr = 223758.64f0,
        dzr = 5602.1216f0,
        kr = 0.014050499f0,
        lr = 3.247451f-6,
        mr = 0.037179865f0,
        n2 = 0.02856384f0,
        nr = 2.6691635f15,
        p = 25764.057f0,
        pip = 0.01291078f0,
        rhobar = 80.18411f0,
        rhop = 0.006651403f0,
        t = 300.0f0,
        thetabar = 29177.266f0,
        us = 1265.1172f0,
        uw = 3.2371957f0,
        vs = 3.5866885f0,
        wts = 0.11083743f0,
        x = 257875.94f0,
        xr = 566029.5f0,
        y = 257875.94f0,
        yr = 552537.8f0,
        z = 258426.88f0,
        zr = 11301.149f0,
        ztilde = 268158.8f0,
    )
    linf = (
        dkr = 6.4709246f-5,
        dlr = 6.3524145f-5,
        dmr = 0.00019787297f0,
        dxr = 10368.59f0,
        dyr = 10027.303f0,
        dzr = 300.68198f0,
        kr = 0.0006294487f0,
        lr = 4.3327887f-7,
        mr = 0.001696954f0,
        n2 = 0.00031935345f0,
        nr = 3.918283f14,
        p = 344.90543f0,
        pip = 0.00042407282f0,
        rhobar = 1.145016f0,
        rhop = 0.00061810936f0,
        t = 300.0f0,
        thetabar = 351.64978f0,
        us = 10.27838f0,
        uw = 0.5433125f0,
        vs = 0.12918764f0,
        wts = 0.0062380005f0,
        x = 95000.0f0,
        xr = 47774.457f0,
        y = 95000.0f0,
        yr = 45007.812f0,
        z = 4879.7573f0,
        zr = 1070.3611f0,
        ztilde = 5000.0f0,
    )
    reference = (l2, linf)

    @testset "WKB mountain wave" begin
        test_example(wkb_mountain_wave, keywords, reference; update)
    end

    return
end
