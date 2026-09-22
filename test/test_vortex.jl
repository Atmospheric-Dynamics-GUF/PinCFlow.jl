function test_vortex()
    l2 = (
        chi = 7.2657404f0,
        pip = 2.5491812f-5,
        rhop = 0.0f0,
        t = 3600.0f0,
        tke = 0.0014142136f0,
        us = 6.676539f0,
        vs = 6.676539f0,
        wts = 0.0f0,
        x = 25787.594f0,
        y = 25787.594f0,
        z = 10000.0f0,
        ztilde = 20000.0f0,
    )
    linf = (
        chi = 0.9514584f0,
        pip = 5.5970127f-6,
        rhop = 0.0f0,
        t = 3600.0f0,
        tke = 5.0f-5,
        us = 1.2179157f0,
        vs = 1.2179157f0,
        wts = 0.0f0,
        x = 9500.0f0,
        y = 9500.0f0,
        z = 500.0f0,
        ztilde = 1000.0f0,
    )
    reference = (l2, linf)

    @testset "Vortex" begin
        test_example(vortex, keywords, reference; update)
    end

    return
end
