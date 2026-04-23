function test_wkb_wave_packet()
    l2 = (
        dkr = 0.006207327f0,
        dlr = 0.0062073325f0,
        dmr = 0.0043892446f0,
        dxr = 27477.207f0,
        dyr = 29274.488f0,
        dzr = 40049.97f0,
        kr = 0.043892622f0,
        lr = 0.043892775f0,
        mr = 0.04389183f0,
        n2 = 0.011484171f0,
        nr = 3.7145567f10,
        p = 5766.7656f0,
        pip = 1.0205832f-5,
        rhobar = 13.219808f0,
        rhop = 3.4633177f-7,
        t = 900.0f0,
        thetabar = 18475.896f0,
        tke = 0.0022360678f0,
        us = 0.00072021f0,
        vs = 0.000708969f0,
        wts = 0.00018230997f0,
        x = 18165.902f0,
        xr = 41577.207f0,
        y = 18165.902f0,
        yr = 39917.16f0,
        z = 729383.3f0,
        zr = 382879.62f0,
        ztilde = 784856.7f0,
    )
    linf = (
        dkr = 0.00035543632f0,
        dlr = 0.00035543588f0,
        dmr = 0.0002513274f0,
        dxr = 2000.0f0,
        dyr = 2000.0f0,
        dzr = 4000.0f0,
        kr = 0.002513444f0,
        lr = 0.0025134473f0,
        mr = 0.0025132827f0,
        n2 = 0.0004251976f0,
        nr = 7.197886f9,
        p = 294.45767f0,
        pip = 5.375786f-7,
        rhobar = 0.98152554f0,
        rhop = 3.1430723f-8,
        t = 900.0f0,
        thetabar = 1003.64667f0,
        tke = 5.0000002f-5,
        us = 0.00016283107f0,
        vs = 0.000161258f0,
        wts = 3.4135137f-5,
        x = 9000.0f0,
        xr = 4504.803f0,
        y = 9000.0f0,
        yr = 4504.7983f0,
        z = 38000.0f0,
        zr = 25009.652f0,
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
