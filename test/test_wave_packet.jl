function test_wave_packet()
    l2 = (
        buoyancyproduction = 0.0024361643f0,
        kek = 6.708211f0,
        kh = 6.708211f0,
        km = 6.708211f0,
        n2 = 0.011484171f0,
        p = 4077.7273f0,
        pip = 1.9590077f-5,
        rhobar = 13.219808f0,
        rhop = 0.00015629096f0,
        shearproduction = 1.3830034f-11,
        t = 900.0f0,
        thetabar = 18475.896f0,
        tke = 0.0022360685f0,
        u = 0.13855109f0,
        us = 0.15149634f0,
        v = 0.13763934f0,
        vs = 0.15065783f0,
        w = 0.19741125f0,
        wts = 0.24634102f0,
        x = 18165.902f0,
        y = 18165.902f0,
        z = 729383.3f0,
        ztilde = 784856.7f0,
    )
    linf = (
        buoyancyproduction = 9.031374f-5,
        kek = 0.21217921f0,
        kh = 0.21217921f0,
        km = 0.21217921f0,
        n2 = 0.0004251976f0,
        p = 294.45767f0,
        pip = 7.229181f-6,
        rhobar = 0.98152554f0,
        rhop = 7.640168f-5,
        shearproduction = 3.9158342f-12,
        t = 900.0f0,
        thetabar = 1003.64667f0,
        tke = 5.001521f-5,
        u = 0.037073277f0,
        us = 0.043309662f0,
        v = 0.03669157f0,
        vs = 0.043617133f0,
        w = 0.051705882f0,
        wts = 0.07246128f0,
        x = 9000.0f0,
        y = 9000.0f0,
        z = 38000.0f0,
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

    @testset "Wave packet" begin
        test_example(wave_packet, keywords, reference; update_references)
    end

    return
end
