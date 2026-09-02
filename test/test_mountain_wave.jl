function test_mountain_wave()
    l2 = (
        n2 = 0.02856384f0,
        p = 25773.348f0,
        pip = 0.00031950124f0,
        rhobar = 80.22628f0,
        rhop = 0.008588406f0,
        t = 600.0f0,
        thetabar = 29173.713f0,
        tke = 0.006324555f0,
        us = 1265.0928f0,
        vs = 2.291456f0,
        w = 1.9963831f0,
        wts = 2.7518039f0,
        x = 25787.594f0,
        y = 25787.594f0,
        z = 258216.62f0,
        ztilde = 267956.44f0,
    )
    linf = (
        n2 = 0.00031935345f0,
        p = 344.89035f0,
        pip = 3.042257f-5,
        rhobar = 1.144946f0,
        rhop = 0.0017714053f0,
        t = 600.0f0,
        thetabar = 351.61438f0,
        tke = 5.004641f-5,
        us = 10.529754f0,
        vs = 0.32042885f0,
        w = 0.29883713f0,
        wts = 0.31943735f0,
        x = 9500.0f0,
        y = 9500.0f0,
        z = 4876.6665f0,
        ztilde = 5000.0f0,
    )
    reference = (l2, linf)

    @testset "Mountain wave" begin
        test_example(mountain_wave, keywords, reference; update)
    end

    return
end
