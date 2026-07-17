function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0014157874f0,
        dlr = 0.0014080801f0,
        dmr = 0.0038584f0,
        dxr = 222761.88f0,
        dyr = 223758.61f0,
        dzr = 5603.218f0,
        kr = 0.014050533f0,
        lr = 3.3125214f-6,
        mr = 0.037176125f0,
        n2 = 0.02856384f0,
        nr = 2.6683832f15,
        p = 25764.057f0,
        pip = 0.012910779f0,
        rhobar = 80.18411f0,
        rhop = 0.0066514253f0,
        t = 300.0f0,
        thetabar = 29177.266f0,
        us = 1265.1171f0,
        uw = 3.237847f0,
        vs = 3.5867116f0,
        wts = 0.11083427f0,
        x = 257875.94f0,
        xr = 566029.7f0,
        y = 257875.94f0,
        yr = 552537.7f0,
        z = 258426.88f0,
        zr = 11305.939f0,
        ztilde = 268158.8f0,
    )
    linf = (
        dkr = 6.47789f-5,
        dlr = 6.3524225f-5,
        dmr = 0.00019791294f0,
        dxr = 10367.9375f0,
        dyr = 10027.304f0,
        dzr = 301.23303f0,
        kr = 0.0006295199f0,
        lr = 4.3339182f-7,
        mr = 0.0016969462f0,
        n2 = 0.00031935345f0,
        nr = 3.9182795f14,
        p = 344.90543f0,
        pip = 0.00042407282f0,
        rhobar = 1.145016f0,
        rhop = 0.0006181091f0,
        t = 300.0f0,
        thetabar = 351.64978f0,
        us = 10.27838f0,
        uw = 0.54324675f0,
        vs = 0.12918915f0,
        wts = 0.0062375395f0,
        x = 95000.0f0,
        xr = 47774.48f0,
        y = 95000.0f0,
        yr = 45007.816f0,
        z = 4879.7573f0,
        zr = 1070.9481f0,
        ztilde = 5000.0f0,
    )
    reference = (l2, linf)

    @testset "WKB mountain wave" begin
        test_example(wkb_mountain_wave, keywords, reference; update)
    end

    return
end
