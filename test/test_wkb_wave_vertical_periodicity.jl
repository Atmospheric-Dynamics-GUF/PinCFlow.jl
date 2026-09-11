function test_wkb_wave_vertical_periodicity()
    l2 = (
        chi = 60573.37f0,
        dchidt1 = 1.0416784f-5,
        dkr = 0.00031136567f0,
        dlr = 0.0f0,
        dmr = 0.009340969f0,
        dxr = 416221.12f0,
        dyr = 148660.69f0,
        dzr = 12569.804f0,
        kr = 0.0031135902f0,
        lr = 0.0f0,
        mr = 0.093405746f0,
        nr = 3.8520124f11,
        pip = 1.5695243f-6,
        rhop = 2.1579326f-5,
        t = 3600.0f0,
        tke = 0.0007071068f0,
        us = 0.04556574f0,
        vs = 0.00856454f0,
        wts = 0.0010427983f0,
        x = 272488.53f0,
        xr = 1.1348975f6,
        y = 0.0f0,
        yr = 0.0f0,
        z = 57662.812f0,
        zr = 83570.94f0,
        ztilde = 62048.367f0,
    )
    linf = (
        chi = 4975.0625f0,
        dchidt1 = 3.343046f-6,
        dkr = 2.0947757f-5,
        dlr = 0.0f0,
        dmr = 0.00062856125f0,
        dxr = 30007.625f0,
        dyr = 10000.0f0,
        dzr = 1000.00244f0,
        kr = 0.00020947757f0,
        lr = 0.0f0,
        mr = 0.0062856483f0,
        nr = 8.081998f10,
        pip = 3.699721f-7,
        rhop = 3.97466f-6,
        t = 3600.0f0,
        tke = 5.0f-5,
        us = 0.011027801f0,
        vs = 0.0021181218f0,
        wts = 0.00018608589f0,
        x = 135000.0f0,
        xr = 140479.75f0,
        y = 0.0f0,
        yr = 0.0f0,
        z = 9500.0f0,
        zr = 9932.626f0,
        ztilde = 10000.0f0,
    )
    reference = (l2, linf)

    @testset "WKB Wave packet vertical periodicity" begin
        test_example(wkb_wave_vertical_periodicity, keywords, reference; update)
    end

    return
end
