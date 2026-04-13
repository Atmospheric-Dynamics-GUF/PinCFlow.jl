function test_periodic_hill()
    l2 = (
        buoyancyproduction = 0.00021213203f0,
        kek = 2.1213202f0,
        kh = 2.1213202f0,
        km = 2.1213202f0,
        pip = 0.000621337f0,
        rhop = 0.008697031f0,
        shearproduction = 2.2032053f-7,
        t = 3600.0f0,
        tke = 0.0007071068f0,
        us = 140.49748f0,
        vs = 0.0f0,
        w = 2.7645657f0,
        wts = 4.28726f0,
        x = 18165.902f0,
        y = 0.0f0,
        z = 116062.984f0,
        ztilde = 124774.234f0,
    )
    linf = (
        buoyancyproduction = 2.1213204f-5,
        kek = 0.21213204f0,
        kh = 0.21213204f0,
        km = 0.21213204f0,
        pip = 0.0002016583f0,
        rhop = 0.0023945163f0,
        shearproduction = 1.01032846f-7,
        t = 3600.0f0,
        tke = 5.0f-5,
        us = 10.880933f0,
        vs = 0.0f0,
        w = 0.70323414f0,
        wts = 1.1136314f0,
        x = 9000.0f0,
        y = 0.0f0,
        z = 19024.389f0,
        ztilde = 20000.0f0,
    )
    reference = (l2, linf)

    keywords =
        (x_size = 10, z_size = 10, prepare_restart = true, visualize = false)

    @testset "Periodic hill" begin
        test_example(periodic_hill, keywords, reference; update_references)
    end

    return
end
