function test_hot_bubble()
    l2 = (
        n2 = 0.0f0,
        p = 2484.3513f0,
        pip = 0.0018407352f0,
        rhobar = 5.855689f0,
        rhop = 0.008935173f0,
        t = 3600.0f0,
        thetabar = 3000.0f0,
        thetap = 9.650683f0,
        tke = 0.020087982f0,
        us = 20.53315f0,
        vs = 0.8346811f0,
        wts = 23.84529f0,
        x = 18165.902f0,
        y = 0.0f0,
        z = 115325.625f0,
        ztilde = 124096.734f0,
    )
    linf = (
        n2 = 0.0f0,
        p = 320.76392f0,
        pip = 0.00059151894f0,
        rhobar = 1.069213f0,
        rhop = 0.0058253715f0,
        t = 3600.0f0,
        thetabar = 300.0f0,
        thetap = 6.075528f0,
        tke = 0.0089118695f0,
        us = 5.930776f0,
        vs = 0.30192846f0,
        wts = 5.0978894f0,
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
