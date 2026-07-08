function test_hot_bubble()
    l2 = (
        n2 = 0.0f0,
        p = 30010.416f0,
        pip = 0.0044959094f0,
        rhobar = 70.73524f0,
        rhop = 0.08934707f0,
        t = 300.0f0,
        thetabar = 26832.816f0,
        thetap = 52.045544f0,
        tke = 0.7188447f0,
        us = 119.88153f0,
        vs = 119.88153f0,
        wts = 171.93753f0,
        x = 12893.797f0,
        y = 12893.797f0,
        z = 516236.38f0,
        ztilde = 535723.8f0,
    )
    linf = (
        n2 = 0.0f0,
        p = 341.38608f0,
        pip = 0.0007016741f0,
        rhobar = 1.1379536f0,
        rhop = 0.009277864f0,
        t = 300.0f0,
        thetabar = 300.0f0,
        thetap = 5.0405197f0,
        tke = 0.08789462f0,
        us = 17.522087f0,
        vs = 17.522087f0,
        wts = 15.267322f0,
        x = 4750.0f0,
        y = 4750.0f0,
        z = 9750.0f0,
        ztilde = 10000.0f0,
    )
    reference = (l2, linf)

    @testset "Hot bubble" begin
        test_example(hot_bubble, keywords, reference; update)
    end

    return
end
