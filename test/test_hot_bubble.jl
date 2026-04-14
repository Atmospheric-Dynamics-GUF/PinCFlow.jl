function test_hot_bubble()
    l2 = (
        buoyancyproduction = 0.00011643528f0,
        kek = 7.818534f0,
        kh = 7.818534f0,
        km = 7.818534f0,
        n2 = 0.0f0,
        p = 2484.3425f0,
        pip = 0.0018536093f0,
        rhobar = 5.855689f0,
        rhop = 0.008940074f0,
        shearproduction = 6.189122f-7,
        t = 3600.0f0,
        thetabar = 3000.0f0,
        thetap = 9.747072f0,
        tke = 0.01908623f0,
        us = 20.135872f0,
        vs = 0.8207722f0,
        wts = 23.656435f0,
        x = 18165.902f0,
        y = 0.0f0,
        z = 115325.625f0,
        ztilde = 124096.734f0,
    )
    linf = (
        buoyancyproduction = 6.497864f-5,
        kek = 2.9067802f0,
        kh = 2.9067802f0,
        km = 2.9067802f0,
        n2 = 0.0f0,
        p = 320.76392f0,
        pip = 0.0006050003f0,
        rhobar = 1.069213f0,
        rhop = 0.0058253715f0,
        shearproduction = 2.0827716f-7,
        t = 3600.0f0,
        thetabar = 300.0f0,
        thetap = 6.075528f0,
        tke = 0.009418048f0,
        us = 5.664044f0,
        vs = 0.30034977f0,
        wts = 4.707085f0,
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
