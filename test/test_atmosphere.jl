@testset "Atmosphere" begin
    semi = initialize_values(3, 1, 1, 3, 3, 3, 0, 60000, 0, 40000, 0, 20000)
    (; equations, cache, grid) = semi

    @test grid.nx == 3
    @test grid.nz == 1
    @test grid.ny == 1
    @test grid.lx[0] ≈ 0.0
    @test grid.lx[1] ≈ 6.877891931902294
    @test grid.ly[0] ≈ 0.0
    @test grid.ly[1] ≈ 4.585261287934863
    @test grid.lz[0] ≈ 0.0
    @test grid.lz[1] ≈ 2.2926306439674313
    @test grid.dx ≈ 2.2926306439674313
    @test grid.dy ≈ 4.585261287934863
    @test grid.dz ≈ 2.2926306439674313

    (; gamma,
    gamma_1,
    kappaInv,
    gammaInv,
    Rsp,
    g,
    rhoRef,
    pRef,
    aRef,
    uRef,
    lRef,
    tRef,
    thetaRef,
    Ma,
    Fr,
    kappa,
    sig,
    press0_dim,
    Temp0_dim,
    T0,
    N2,
    NN,
    mu_viscous_dim,
    ReInv,
    Re,
    g_ndim,
    MaInv2,
    maxIter,
    maxTime,
    tStepChoice,
    dtMax_dim,
    cfl,
    dt,
    model,
    backgroundFlow_dim,
    f_Coriolis_dim,
    corset,
    background) = equations

    tol = 1e-14

    @test isapprox(equations.gamma, 1.4; atol = tol, rtol = tol)
    @test isapprox(equations.gamma_1, 0.4; atol = tol, rtol = tol)
    @test isapprox(equations.kappaInv, 1.4 / 0.4; atol = tol, rtol = tol)
    @test isapprox(equations.gammaInv, 1.0 / 1.4; atol = tol, rtol = tol)
    @test isapprox(equations.Rsp, 287.0; atol = tol, rtol = tol)
    @test isapprox(equations.g, 9.81; atol = tol, rtol = tol)
    @test isapprox(equations.pRef, 101325.0; atol = tol, rtol = tol)
    @test isapprox(gamma, 1.4, atol = 1e-14, rtol = 1e-14)
    @test isapprox(g, 9.81, atol = 1e-14, rtol = 1e-14)
    @test isapprox(rhoRef, 1.184, atol = 1e-14, rtol = 1e-14)
    @test isapprox(pRef, 101325.0, atol = 1e-14, rtol = 1e-14)
    @test isapprox(aRef, 292.5381125550948, atol = 1e-14, rtol = 1e-14)
    @test isapprox(uRef, 292.5381125550948, atol = 1e-14, rtol = 1e-14)
    @test isapprox(lRef, 8723.60319034631, atol = 1e-14, rtol = 1e-14)
    @test isapprox(tRef, 29.820398833342992, atol = 1e-14, rtol = 1e-14)
    @test isapprox(thetaRef, 298.1830916282137, atol = 1e-14, rtol = 1e-14)
    @test isapprox(Ma, 1.0, atol = 1e-14, rtol = 1e-14)
    @test isapprox(Fr, 1.0, atol = 1e-14, rtol = 1e-14)
    @test isapprox(kappa, 0.28571428571428564, atol = 1e-14, rtol = 1e-14)
    @test isapprox(sig, 1.0, atol = 1e-14, rtol = 1e-14)
    @test isapprox(press0_dim, 100000.0, atol = 1e-14, rtol = 1e-14)
    @test isapprox(Temp0_dim, 300.0, atol = 1e-14, rtol = 1e-14)
    @test isapprox(T0, 1.0060932642487044, atol = 1e-14, rtol = 1e-14)
    @test isapprox(N2, 0.28398389678877484, atol = 1e-14, rtol = 1e-14)
    @test isapprox(NN, 0.5329013949960864, atol = 1e-14, rtol = 1e-14)
    @test isapprox(mu_viscous_dim, 0.0, atol = 1e-14, rtol = 1e-14)
    @test isapprox(ReInv, 0.0, atol = 1e-14, rtol = 1e-14)
    @test isapprox(Re, 1.0e20, atol = 1e-14, rtol = 1e-14)
    @test isapprox(g_ndim, 0.9999999999999998, atol = 1e-14, rtol = 1e-14)
    @test isapprox(MaInv2, 1.0, atol = 1e-14, rtol = 1e-14)
    @test isapprox(maxIter, 5000, atol = 1e-14, rtol = 1e-14)
    @test isapprox(maxTime, 3600.0, atol = 1e-14, rtol = 1e-14)
    @test tStepChoice == "cfl"
    @test isapprox(dtMax_dim, 1.0, atol = 1e-14, rtol = 1e-14)
    @test isapprox(cfl, 0.05, atol = 1e-14, rtol = 1e-14)
    @test isapprox(dt, 0.033534092068610195, atol = 1e-14, rtol = 1e-14)
    @test model == "pseudo_incompressible"
    @test sum(backgroundFlow_dim .≈ (10.0, 0.0, 0.0)) == 3
    @test isapprox(f_Coriolis_dim, 0.0, atol = 1e-14, rtol = 1e-14)
    @test corset == "constant"
    @test background == "isothermal"

    initialize_atmosphere!(semi)

    initialize_variables!(semi)

    (; topography_surface, zTFC, zS) = grid

    # TODO - Improve these tests when jac, met are arrays
    @test cache.jac[2, 2, 2] ≈ 0.98
    @test semi.met[1, 1, 1, 1, 3] ≈ -0.004987779939149084
    @test semi.met[3, 1, 1, 2, 3] ≈ 0.0
    @test semi.met[1, 1, 1, 3, 3] ≈ 1.0001246360352993

    test_arr(cache.pStrat, 1107.533828869411, 101.40469577297178, 11.340157693052404,
             tol = 1e-14)
    test_arr(cache.rhoStrat, 2595.1704734926416, 261.8137328357221, 29.930580142611795,
             tol = 1e-14)
    test_arr(cache.bvsStrat, 89.74505372656458, 5.21192418228228, 0.4002353847198065,
             tol = 1e-14)
    test_arr(cache.thetaStrat, 821.7649286701387, 53.45708015679106, 5.12305138831322,
             tol = 1e-14)

    test_arr(topography_surface, 0.9677073885235098, 0.21012437605443837,
             0.04585261287934863, tol = 1e-14)
    test_arr(zTFC, 2919.0839060619323, 132.19881726855255, 8.02420725388601, tol = 1e-14)
    test_arr(zS, 36.6820903034789, 14.857944720776754, 8.02420725388601, tol = 1e-14)

    sum_vert_wind = sum(PinCFlow_dev.vertWind(1, 1, k, semi) for k in -2:3)
    @test sum_vert_wind ≈ 0.00102295002764069
end
