function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0014379211f0,
        dlr = 0.0014362237f0,
        dmr = 0.003909416f0,
        dxr = 218585.31f0,
        dyr = 228499.08f0,
        dzr = 5355.4595f0,
        kr = 0.014355713f0,
        lr = 1.4668492f-6,
        mr = 0.038147204f0,
        n2 = 0.04039537f0,
        nr = 3.2829122f15,
        p = 36448.062f0,
        pip = 0.000411167f0,
        rhobar = 113.4528f0,
        rhop = 0.0027354418f0,
        t = 300.0f0,
        thetabar = 41258.24f0,
        us = 1264.9856f0,
        uw = 3.0415967f0,
        vs = 1.8948448f0,
        wts = 0.025369892f0,
        x = 257875.94f0,
        xr = 583925.7f0,
        y = 257875.94f0,
        yr = 573979.3f0,
        z = 365197.1f0,
        zr = 12168.835f0,
        ztilde = 378970.5f0,
    )
    linf = (
        dkr = 6.341313f-5,
        dlr = 6.298655f-5,
        dmr = 0.00018481903f0,
        dxr = 10343.222f0,
        dyr = 10007.269f0,
        dzr = 360.3735f0,
        kr = 0.0006288331f0,
        lr = 2.2346704f-7,
        mr = 0.0016861538f0,
        n2 = 0.00031935345f0,
        nr = 8.845047f14,
        p = 344.90543f0,
        pip = 2.0523112f-5,
        rhobar = 1.145016f0,
        rhop = 0.0002496933f0,
        t = 300.0f0,
        thetabar = 351.63614f0,
        us = 10.190177f0,
        uw = 0.8631419f0,
        vs = 0.07521223f0,
        wts = 0.0019595667f0,
        x = 95000.0f0,
        xr = 47787.3f0,
        y = 95000.0f0,
        yr = 45004.29f0,
        z = 4878.568f0,
        zr = 939.103f0,
        ztilde = 5000.0f0,
    )
    reference = (l2, linf)

    @testset "WKB mountain wave" begin
        test_example(wkb_mountain_wave, keywords, reference; update)
    end

    nothing
end
