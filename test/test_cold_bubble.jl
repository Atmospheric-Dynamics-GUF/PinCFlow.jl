function test_cold_bubble()
    l2 = (
        n2 = 0.0f0,
        p = 1756.7067f0,
        pip = 0.00038125357f0,
        rhobar = 5.855689f0,
        rhop = 0.0091435f0,
        t = 3600.0f0,
        thetabar = 3000.0f0,
        thetap = 4.445302f0,
        us = 8.429899f0,
        vs = 2.962789f0,
        wts = 9.879302f0,
        x = 18165.902f0,
        y = 0.0f0,
        z = 115325.625f0,
        ztilde = 124096.734f0,
    )
    linf = (
        n2 = 0.0f0,
        p = 320.76392f0,
        pip = 9.277293f-5,
        rhobar = 1.069213f0,
        rhop = 0.0058253715f0,
        t = 3600.0f0,
        thetabar = 300.0f0,
        thetap = 2.8450513f0,
        us = 2.198411f0,
        vs = 0.64761156f0,
        wts = 3.5191293f0,
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
