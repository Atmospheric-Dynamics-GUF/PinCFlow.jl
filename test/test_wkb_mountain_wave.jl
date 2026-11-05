function test_wkb_mountain_wave()
    l2 = (
        dkr = 0.0015521273f0,
        dlr = 0.0015486742f0,
        dmr = 0.0041514807f0,
        dxr = 247336.0f0,
        dyr = 246244.7f0,
        dzr = 5450.598f0,
        kr = 0.015469856f0,
        lr = 3.2606072f-6,
        mr = 0.04092164f0,
        n2 = 0.02856384f0,
        nr = 3.0979288f15,
        p = 25764.057f0,
        pip = 0.012910715f0,
        rhobar = 80.18411f0,
        rhop = 0.006657116f0,
        t = 300.0f0,
        thetabar = 29177.266f0,
        us = 1265.1134f0,
        uw = 3.2381527f0,
        vs = 3.5894208f0,
        wts = 0.1101422f0,
        x = 257875.94f0,
        xr = 614493.5f0,
        y = 257875.94f0,
        yr = 580276.9f0,
        z = 258426.88f0,
        zr = 13577.741f0,
        ztilde = 268158.8f0,
    )
    linf = (
        dkr = 6.437767f-5,
        dlr = 6.3473184f-5,
        dmr = 0.00019561134f0,
        dxr = 10348.517f0,
        dyr = 10025.849f0,
        dzr = 354.03897f0,
        kr = 0.0006298297f0,
        lr = 4.160275f-7,
        mr = 0.0016961319f0,
        n2 = 0.00031935345f0,
        nr = 3.917712f14,
        p = 344.90543f0,
        pip = 0.00042407282f0,
        rhobar = 1.145016f0,
        rhop = 0.0006181314f0,
        t = 300.0f0,
        thetabar = 351.64978f0,
        us = 10.279384f0,
        uw = 0.5600827f0,
        vs = 0.12932862f0,
        wts = 0.0061857565f0,
        x = 95000.0f0,
        xr = 47795.52f0,
        y = 95000.0f0,
        yr = 45007.598f0,
        z = 4879.7573f0,
        zr = 1161.098f0,
        ztilde = 5000.0f0,
    )
    reference = (l2, linf)

    @testset "WKB mountain wave" begin
        test_example(wkb_mountain_wave, keywords, reference; update)
    end

    nothing
end
