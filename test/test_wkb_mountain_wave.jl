function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0018781841f0,
        dlr = 0.0018785482f0,
        dmr = 0.002556423f0,
        dxr = 788492.9f0,
        dyr = 788363.44f0,
        dzr = 18077.58f0,
        kr = 0.018782519f0,
        lr = 3.2470582f-6,
        mr = 0.041186426f0,
        n2 = 0.012520827f0,
        nr = 2.925965f14,
        p = 5733.0596f0,
        pip = 0.00059674983f0,
        rhobar = 17.83844f0,
        rhop = 0.00033393208f0,
        t = 3600.0f0,
        thetabar = 12638.998f0,
        tke = 0.002236068f0,
        us = 447.22708f0,
        vs = 0.19900095f0,
        wts = 0.030781437f0,
        x = 363318.03f0,
        xr = 608380.4f0,
        y = 363318.03f0,
        yr = 651952.6f0,
        z = 364706.28f0,
        zr = 221004.02f0,
        ztilde = 392441.75f0,
    )
    linf = (
        dkr = 6.300841f-5,
        dlr = 6.29236f-5,
        dmr = 0.0003950022f0,
        dxr = 40000.0f0,
        dyr = 40000.0f0,
        dzr = 1996.0245f0,
        kr = 0.0006288675f0,
        lr = 2.8359972f-7,
        mr = 0.0023444057f0,
        n2 = 0.00057970756f0,
        nr = 1.1042035f13,
        p = 320.91135f0,
        pip = 6.432183f-5,
        rhobar = 1.0579953f0,
        rhop = 0.00012263547f0,
        t = 3600.0f0,
        thetabar = 569.07f0,
        tke = 5.0005798f-5,
        us = 10.116935f0,
        vs = 0.054423608f0,
        wts = 0.0075844564f0,
        x = 180000.0f0,
        xr = 29975.828f0,
        y = 180000.0f0,
        yr = 30014.432f0,
        z = 19001.988f0,
        zr = 13316.544f0,
        ztilde = 20000.0f0,
    )
    reference = (l2, linf)

    keywords = (
        x_size = 10,
        y_size = 10,
        z_size = 10,
        prepare_restart = true,
        visualize = false,
    )

    @testset "WKB mountain wave" begin
        test_example(wkb_mountain_wave, keywords, reference; update_references)
    end

    return
end
