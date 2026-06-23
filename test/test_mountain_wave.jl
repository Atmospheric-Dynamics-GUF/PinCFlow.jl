function test_mountain_wave()
    l2 = (
        n2 = 0.010105865f0,
        p = 5973.465f0,
        pip = 4.321387f-5,
        rhobar = 17.034365f0,
        rhop = 0.0011904156f0,
        t = 3600.0f0,
        thetabar = 13598.672f0,
        tke = 0.002236067f0,
        us = 447.2375f0,
        vs = 0.33395347f0,
        w = 0.22903982f0,
        wts = 0.55782604f0,
        x = 18165.902f0,
        y = 18165.902f0,
        z = 364725.28f0,
        ztilde = 392459.12f0,
    )
    linf = (
        n2 = 0.00031957557f0,
        p = 321.18332f0,
        pip = 1.0270111f-5,
        rhobar = 1.0363004f0,
        rhop = 0.00040521479f0,
        t = 3600.0f0,
        thetabar = 556.882f0,
        tke = 5.001706f-5,
        us = 10.13566f0,
        vs = 0.050962564f0,
        w = 0.059119597f0,
        wts = 0.13209596f0,
        x = 9000.0f0,
        y = 9000.0f0,
        z = 19001.666f0,
        ztilde = 20000.0f0,
    )
    reference = (l2, linf)

    keywords = (
        x_size = 10,
        y_size = 10,
        z_size = 10,
        prepare_restart = true,
        visualize = false,
    )

    @testset "Mountain wave" begin
        test_example(mountain_wave, keywords, reference; update)
    end

    return
end
