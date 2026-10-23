function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.001332497f0,
        dlr = 0.0013306037f0,
        dmr = 0.0036350584f0,
        dxr = 198547.92f0,
        dyr = 211686.33f0,
        dzr = 4897.639f0,
        kr = 0.01329964f0,
        lr = 1.3212675f-6,
        mr = 0.03532602f0,
        n2 = 0.04039537f0,
        nr = 3.2408218f15,
        p = 36448.062f0,
        pip = 0.00041995512f0,
        rhobar = 113.4528f0,
        rhop = 0.0027996362f0,
        t = 300.0f0,
        thetabar = 41258.24f0,
        us = 1264.9857f0,
        uw = 3.071055f0,
        vs = 1.8981813f0,
        wts = 0.025909295f0,
        x = 257875.94f0,
        xr = 536958.0f0,
        y = 257875.94f0,
        yr = 516137.44f0,
        z = 365197.1f0,
        zr = 12142.003f0,
        ztilde = 378970.5f0,
    )
    linf = (
        dkr = 6.333f-5,
        dlr = 6.297292f-5,
        dmr = 0.00018435801f0,
        dxr = 10286.999f0,
        dyr = 10007.421f0,
        dzr = 392.53662f0,
        kr = 0.0006288009f0,
        lr = 2.1093372f-7,
        mr = 0.00168578f0,
        n2 = 0.00031935345f0,
        nr = 8.8445624f14,
        p = 344.90543f0,
        pip = 2.0716772f-5,
        rhobar = 1.145016f0,
        rhop = 0.00025560477f0,
        t = 300.0f0,
        thetabar = 351.63614f0,
        us = 10.190774f0,
        uw = 0.8809092f0,
        vs = 0.07595289f0,
        wts = 0.0019564042f0,
        x = 95000.0f0,
        xr = 47762.24f0,
        y = 95000.0f0,
        yr = 45003.887f0,
        z = 4878.568f0,
        zr = 956.95337f0,
        ztilde = 5000.0f0,
    )
    reference = (l2, linf)

    @testset "WKB mountain wave" begin
        test_example(wkb_mountain_wave, keywords, reference; update)
    end

    return
end
