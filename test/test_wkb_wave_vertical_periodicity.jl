function test_wkb_wave_vertical_periodicity()
    l2 = (
        dkr = 0.0013395981f0,
        dlr = 0.0f0,
        dmr = 0.040187888f0,
        dxr = 463731.12f0,
        dyr = 639609.25f0,
        dzr = 12423.517f0,
        kr = 0.013395963f0,
        lr = 0.0f0,
        mr = 0.40187797f0,
        nr = 5.9232735f10,
        pip = 3.1354924f-7,
        rhop = 5.1562433f-6,
        t = 3600.0f0,
        tke = 0.0028284271f0,
        u = 0.008616386f0,
        us = 0.00863608f0,
        vs = 0.0016766292f0,
        wts = 0.0001973366f0,
        x = 547551.4f0,
        xr = 4.955065f6,
        y = 0.0f0,
        yr = 0.0f0,
        z = 230922.06f0,
        zr = 361828.12f0,
        ztilde = 235265.81f0,
    )
    linf = (
        dkr = 2.094473f-5,
        dlr = 0.0f0,
        dmr = 0.00062834675f0,
        dxr = 11250.044f0,
        dyr = 10000.0f0,
        dzr = 250.00024f0,
        kr = 0.00020944198f0,
        lr = 0.0f0,
        mr = 0.0062833494f0,
        nr = 3.5491046f9,
        pip = 2.1395728f-8,
        rhop = 2.6830966f-7,
        t = 3600.0f0,
        tke = 5.0f-5,
        u = 0.0005892883f0,
        us = 0.00061238976f0,
        vs = 0.00011413504f0,
        wts = 1.14418945f-5,
        x = 146250.0f0,
        xr = 149184.78f0,
        y = 0.0f0,
        yr = 0.0f0,
        z = 9875.0f0,
        zr = 9995.109f0,
        ztilde = 10000.0f0,
    )
    reference = (l2, linf)

    @testset "WKB Wave packet vertical periodicity" begin
        test_example(wkb_wave_vertical_periodicity, keywords, reference; update)
    end

    return
end
