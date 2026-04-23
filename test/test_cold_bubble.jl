function test_cold_bubble()
    l2 = (
        n2 = 0.0f0,
        p = 1756.7067f0,
        pip = 0.00038122848f0,
        rhobar = 5.855689f0,
        rhop = 0.009143565f0,
        t = 3600.0f0,
        thetabar = 3000.0f0,
        thetap = 4.4452977f0,
        tke = 0.0043104338f0,
        us = 8.426216f0,
        vs = 2.962375f0,
        wts = 9.875353f0,
        x = 18165.902f0,
        y = 0.0f0,
        z = 115325.625f0,
        ztilde = 124096.734f0,
    )
    linf = (
        n2 = 0.0f0,
        p = 320.76392f0,
        pip = 9.276486f-5,
        rhobar = 1.069213f0,
        rhop = 0.0058253715f0,
        t = 3600.0f0,
        thetabar = 300.0f0,
        thetap = 2.8450513f0,
        tke = 0.0018422024f0,
        us = 2.19673f0,
        vs = 0.64749074f0,
        wts = 3.5177119f0,
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
