function test_cold_bubble()
    l2 = (
        buoyancyproduction = 1.0340406f-5,
        kek = 4.2448754f0,
        kh = 4.2448754f0,
        km = 4.2448754f0,
        n2 = 0.0f0,
        p = 1756.7067f0,
        pip = 0.00038125078f0,
        rhobar = 5.855689f0,
        rhop = 0.009143527f0,
        shearproduction = 1.13948f-6,
        t = 3600.0f0,
        thetabar = 3000.0f0,
        thetap = 4.445302f0,
        tke = 0.004332921f0,
        us = 8.4282f0,
        vs = 2.9626446f0,
        wts = 9.877679f0,
        x = 18165.902f0,
        y = 0.0f0,
        z = 115325.625f0,
        ztilde = 124096.734f0,
    )
    linf = (
        buoyancyproduction = 5.283655f-6,
        kek = 1.2896109f0,
        kh = 1.2896109f0,
        km = 1.2896109f0,
        n2 = 0.0f0,
        p = 320.76392f0,
        pip = 9.277032f-5,
        rhobar = 1.069213f0,
        rhop = 0.0058253715f0,
        shearproduction = 6.185649f-7,
        t = 3600.0f0,
        thetabar = 300.0f0,
        thetap = 2.8450513f0,
        tke = 0.0018419594f0,
        us = 2.1976736f0,
        vs = 0.647581f0,
        wts = 3.5185103f0,
        x = 9000.0f0,
        y = 0.0f0,
        z = 19000.0f0,
        ztilde = 20000.0f0,
    )
    reference = (l2, linf)

    keywords =
        (x_size = 10, z_size = 10, prepare_restart = true, visualize = false)

    @testset "Cold bubble" begin
        test_example(cold_bubble, keywords, reference; update_references)
    end

    return
end
