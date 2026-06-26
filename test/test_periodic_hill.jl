function test_periodic_hill()
    l2 = (
        pip = 0.00055369455f0,
        rhop = 0.0053984374f0,
        t = 3600.0f0,
        tke = 0.0007071068f0,
        us = 140.5774f0,
        vs = 0.0f0,
        w = 2.124289f0,
        wts = 4.0749683f0,
        x = 18165.902f0,
        y = 0.0f0,
        z = 116062.984f0,
        ztilde = 124774.234f0,
    )
    linf = (
        pip = 0.0002016583f0,
        rhop = 0.0016297551f0,
        t = 3600.0f0,
        tke = 5.0f-5,
        us = 10.618661f0,
        vs = 0.0f0,
        w = 0.5624462f0,
        wts = 0.988464f0,
        x = 9000.0f0,
        y = 0.0f0,
        z = 19024.389f0,
        ztilde = 20000.0f0,
    )
    reference = (l2, linf)

    keywords =
        (x_size = 10, z_size = 10, prepare_restart = true, visualize = false)

    @testset "Periodic hill" begin
        test_example(periodic_hill, keywords, reference; update)
    end

    return
end
