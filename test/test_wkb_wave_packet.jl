function test_wkb_wave_packet()
    l2 = (
        dkr = 0.006197144f0,
        dlr = 0.0061971447f0,
        dmr = 0.004382043f0,
        dxr = 28425.285f0,
        dyr = 28213.412f0,
        dzr = 40000.0f0,
        kr = 0.043820158f0,
        lr = 0.043820117f0,
        mr = 0.043821115f0,
        n2 = 0.011484171f0,
        nr = 3.7145567f10,
        p = 5766.7773f0,
        pip = 2.4653392f-7,
        rhobar = 13.219808f0,
        rhop = 3.4271082f-8,
        t = 900.0f0,
        thetabar = 18475.896f0,
        us = 0.0030097133f0,
        vs = 0.0029441144f0,
        wts = 5.8236925f-5,
        x = 18165.902f0,
        xr = 40935.902f0,
        y = 18165.902f0,
        yr = 38780.145f0,
        z = 729383.3f0,
        zr = 381800.12f0,
        ztilde = 784856.7f0,
    )
    linf = (
        dkr = 0.0003554358f0,
        dlr = 0.00035543618f0,
        dmr = 0.0002513274f0,
        dxr = 2000.0f0,
        dyr = 2000.0f0,
        dzr = 4000.0f0,
        kr = 0.0025134152f0,
        lr = 0.0025134177f0,
        mr = 0.0025133824f0,
        n2 = 0.0004251976f0,
        nr = 7.197886f9,
        p = 294.45767f0,
        pip = 1.118148f-8,
        rhobar = 0.98152554f0,
        rhop = 5.432737f-9,
        t = 900.0f0,
        thetabar = 1003.64667f0,
        us = 0.0007413181f0,
        vs = 0.00072259636f0,
        wts = 1.017745f-5,
        x = 9000.0f0,
        xr = 4504.834f0,
        y = 9000.0f0,
        yr = 4504.8833f0,
        z = 38000.0f0,
        zr = 25009.656f0,
        ztilde = 40000.0f0,
    )
    reference = (l2, linf)

    keywords = (
        x_size = 10,
        y_size = 10,
        z_size = 10,
        prepare_restart = true,
        visualize = false,
    )

    @testset "WKB Wave packet" begin
        test_example(wkb_wave_packet, keywords, reference; update_references)
    end

    return
end
