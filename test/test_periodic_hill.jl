function test_periodic_hill()
    l2 = (
        pip = 0.0012510886f0,
        rhop = 0.021196544f0,
        t = 3600.0f0,
        tke = 0.0014142136f0,
        us = 281.0013f0,
        vs = 0.0f0,
        w = 6.4545374f0,
        wts = 9.853274f0,
        x = 25787.594f0,
        y = 0.0f0,
        z = 232335.95f0,
        ztilde = 240994.72f0,
    )
    linf = (
        pip = 0.00024125425f0,
        rhop = 0.0024787914f0,
        t = 3600.0f0,
        tke = 5.0f-5,
        us = 11.487826f0,
        vs = 0.0f0,
        w = 0.7097565f0,
        wts = 1.3313593f0,
        x = 9500.0f0,
        y = 0.0f0,
        z = 19512.424f0,
        ztilde = 20000.0f0,
    )
    reference = (l2, linf)

    @testset "Periodic hill" begin
        test_example(periodic_hill, keywords, reference; update)
    end

    return
end
