function test_cold_bubble()
    l2 = (
        n2 = 0.0f0,
        p = 30010.42f0,
        pip = 0.0018711103f0,
        rhobar = 100.03474f0,
        rhop = 0.09722274f0,
        t = 300.0f0,
        thetabar = 37947.332f0,
        thetap = 29.200098f0,
        tke = 0.293085f0,
        us = 60.43407f0,
        vs = 60.43407f0,
        wts = 120.44959f0,
        x = 12893.797f0,
        y = 12893.797f0,
        z = 730068.5f0,
        ztilde = 757627.9f0,
    )
    linf = (
        n2 = 0.0f0,
        p = 341.38608f0,
        pip = 0.0004293566f0,
        rhobar = 1.1379536f0,
        rhop = 0.009421078f0,
        t = 300.0f0,
        thetabar = 300.0f0,
        thetap = 2.9994044f0,
        tke = 0.04835543f0,
        us = 8.29613f0,
        vs = 8.29613f0,
        wts = 14.577412f0,
        x = 4750.0f0,
        y = 4750.0f0,
        z = 9750.0f0,
        ztilde = 10000.0f0,
    )
    reference = (l2, linf)

    @testset "Cold bubble" begin
        test_example(cold_bubble, keywords, reference; update)
    end

    nothing
end
