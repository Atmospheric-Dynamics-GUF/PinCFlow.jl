function test_hot_bubble()
    l2 = (
        n2 = 0.0f0,
        p = 2484.3513f0,
        pip = 0.0018407463f0,
        rhobar = 5.855689f0,
        rhop = 0.008935174f0,
        t = 3600.0f0,
        thetabar = 3000.0f0,
        thetap = 9.650698f0,
        tke = 0.020130103f0,
        us = 20.533022f0,
        vs = 0.8346823f0,
        wts = 23.84519f0,
        x = 18165.902f0,
        y = 0.0f0,
        z = 115325.625f0,
        ztilde = 124096.734f0,
    )
    linf = (
        n2 = 0.0f0,
        p = 320.76392f0,
        pip = 0.0005915232f0,
        rhobar = 1.069213f0,
        rhop = 0.0058253715f0,
        t = 3600.0f0,
        thetabar = 300.0f0,
        thetap = 6.075528f0,
        tke = 0.00893059f0,
        us = 5.930701f0,
        vs = 0.30192947f0,
        wts = 5.097805f0,
        x = 9000.0f0,
        y = 0.0f0,
        z = 19000.0f0,
        ztilde = 20000.0f0,
    )
    reference = (l2, linf)

    keywords =
        (x_size = 10, z_size = 10, prepare_restart = true, visualize = false)

    @testset "Hot bubble" begin
        test_example(hot_bubble, keywords, reference; update_references)
    end

    return
end
