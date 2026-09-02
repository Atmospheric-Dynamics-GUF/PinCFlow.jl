function test_wkb_wave_vertical_periodicity()
    l2 = (
        dkr = 0.0006943187f0,
        dlr = 0.0f0,
        dmr = 0.020829568f0,
        dxr = 458972.38f0,
        dyr = 331511.7f0,
        dzr = 13511.57f0,
        kr = 0.006943169f0,
        lr = 0.0f0,
        mr = 0.20829494f0,
        nr = 3.1329065f10,
        pip = 1.4885413f-7,
        rhop = 2.3603181f-6,
        t = 3600.0f0,
        tke = 0.0014142136f0,
        u = 0.004121158f0,
        us = 0.00414913f0,
        vs = 0.00080253463f0,
        wts = 9.438114f-5,
        x = 386813.9f0,
        xr = 2.5600755f6,
        y = 0.0f0,
        yr = 0.0f0,
        z = 115433.96f0,
        zr = 180287.47f0,
        ztilde = 119791.484f0,
    )
    linf = (
        dkr = 2.0945137f-5,
        dlr = 0.0f0,
        dmr = 0.0006283811f0,
        dxr = 22500.176f0,
        dyr = 10000.0f0,
        dzr = 750.0002f0,
        kr = 0.00020944173f0,
        lr = 0.0f0,
        mr = 0.006283339f0,
        nr = 3.4838643f9,
        pip = 1.9769615f-8,
        rhop = 2.4253916f-7,
        t = 3600.0f0,
        tke = 5.0f-5,
        u = 0.00053904904f0,
        us = 0.0005472084f0,
        vs = 0.00010716307f0,
        wts = 8.967339f-6,
        x = 142500.0f0,
        xr = 147978.28f0,
        y = 0.0f0,
        yr = 0.0f0,
        z = 9750.0f0,
        zr = 9807.609f0,
        ztilde = 10000.0f0,
    )
    reference = (l2, linf)

    @testset "WKB Wave packet vertical periodicity" begin
        test_example(wkb_wave_vertical_periodicity, keywords, reference; update)
    end

    return
end
