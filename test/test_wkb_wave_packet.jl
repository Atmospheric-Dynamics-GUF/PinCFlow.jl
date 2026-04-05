function test_wkb_wave_packet()
    l2 = (
        dkr = 0.0061971447f0,
        dlr = 0.0061971447f0,
        dmr = 0.004382043f0,
        dxr = 27622.453f0,
        dyr = 29206.162f0,
        dzr = 40000.0f0,
        kr = 0.043820426f0,
        lr = 0.043820426f0,
        mr = 0.043820444f0,
        n2 = 0.011484171f0,
        nr = 3.6936524f10,
        p = 5766.7773f0,
        pip = 2.430678f-7,
        rhobar = 13.219808f0,
        rhop = 4.3287912f-8,
        t = 900.0f0,
        thetabar = 18475.896f0,
        us = 0.003316415f0,
        vs = 0.0031345782f0,
        wts = 2.8396018f-5,
        x = 18165.902f0,
        xr = 41103.484f0,
        y = 18165.902f0,
        yr = 38171.094f0,
        z = 729383.3f0,
        zr = 381946.75f0,
        ztilde = 784856.7f0,
    )
    linf = (
        dkr = 0.00035543073f0,
        dlr = 0.00035543076f0,
        dmr = 0.0002513274f0,
        dxr = 2000.0f0,
        dyr = 2000.0f0,
        dzr = 4000.0f0,
        kr = 0.0025132773f0,
        lr = 0.0025132773f0,
        mr = 0.0025132766f0,
        n2 = 0.0004251976f0,
        nr = 7.197886f9,
        p = 294.45767f0,
        pip = 1.1186584f-8,
        rhobar = 0.98152554f0,
        rhop = 8.370002f-9,
        t = 900.0f0,
        thetabar = 1003.64667f0,
        us = 0.0007332119f0,
        vs = 0.00067978573f0,
        wts = 5.2112564f-6,
        x = 9000.0f0,
        xr = 4504.8203f0,
        y = 9000.0f0,
        yr = 4504.822f0,
        z = 38000.0f0,
        zr = 25009.64f0,
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
