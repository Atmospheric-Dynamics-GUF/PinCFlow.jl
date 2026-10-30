function test_wkb_wave_vertical_periodicity()
    l2 = (
        dkr = 0.0007446156f0,
        dlr = 0.0f0,
        dmr = 0.022338469f0,
        dxr = 533291.7f0,
        dyr = 355527.78f0,
        dzr = 11757.976f0,
        kr = 0.0074461563f0,
        lr = 0.0f0,
        mr = 0.2233847f0,
        nr = 3.348891f10,
        pip = 1.4950321f-7,
        rhop = 2.3718956f-6,
        t = 3600.0f0,
        tke = 0.0014142136f0,
        u = 0.0041400203f0,
        us = 0.004167613f0,
        vs = 0.0008043565f0,
        wts = 9.516912f-5,
        x = 386813.9f0,
        xr = 2.6784802f6,
        y = 0.0f0,
        yr = 0.0f0,
        z = 115433.96f0,
        zr = 198639.28f0,
        ztilde = 119791.484f0,
    )
    linf = (
        dkr = 2.0943951f-5,
        dlr = 0.0f0,
        dmr = 0.0006283185f0,
        dxr = 15000.0f0,
        dyr = 10000.0f0,
        dzr = 500.0f0,
        kr = 0.00020943952f0,
        lr = 0.0f0,
        mr = 0.0062831854f0,
        nr = 3.4838643f9,
        pip = 1.9886539f-8,
        rhop = 2.4372275f-7,
        t = 3600.0f0,
        tke = 5.0f-5,
        u = 0.0005357763f0,
        us = 0.0005449238f0,
        vs = 0.000106871805f0,
        wts = 9.071602f-6,
        x = 142500.0f0,
        xr = 147978.28f0,
        y = 0.0f0,
        yr = 0.0f0,
        z = 9750.0f0,
        zr = 9807.608f0,
        ztilde = 10000.0f0,
    )
    reference = (l2, linf)

    @testset "WKB Wave packet vertical periodicity" begin
        test_example(wkb_wave_vertical_periodicity, keywords, reference; update)
    end

    return
end
