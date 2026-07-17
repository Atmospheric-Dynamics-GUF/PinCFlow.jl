function test_cold_bubble()
    l2 = (
        n2 = 0.0f0,
        p = 21220.572f0,
        pip = 0.0018712163f0,
        rhobar = 70.73524f0,
        rhop = 0.09722272f0,
        t = 300.0f0,
        thetabar = 26832.816f0,
        thetap = 29.200014f0,
        tke = 0.1826098f0,
        us = 60.441563f0,
        vs = 60.441563f0,
        wts = 120.458176f0,
        x = 12893.797f0,
        y = 12893.797f0,
        z = 516236.38f0,
        ztilde = 535723.8f0,
    )
    linf = (
        n2 = 0.0f0,
        p = 341.38608f0,
        pip = 0.00042946145f0,
        rhobar = 1.1379536f0,
        rhop = 0.009421057f0,
        t = 300.0f0,
        thetabar = 300.0f0,
        thetap = 2.9994044f0,
        tke = 0.030105574f0,
        us = 8.299103f0,
        vs = 8.299103f0,
        wts = 14.577757f0,
        x = 4750.0f0,
        y = 4750.0f0,
        z = 9750.0f0,
        ztilde = 10000.0f0,
    )
    reference = (l2, linf)

    @testset "Cold bubble" begin
        test_example(cold_bubble, keywords, reference; update)
    end

    return
end
