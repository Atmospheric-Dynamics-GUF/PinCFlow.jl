@testset "Initialization" begin
  domain = PinCFlow.DomainNamelist(;
    sizex = 3,
    sizey = 1,
    sizez = 1,
    lx_dim = [0.0, 60000.0],
    ly_dim = [0.0, 40000.0],
    lz_dim = [0.0, 20000.0],
  )
  grid = PinCFlow.GridNamelist(;
    mountain_case = 3,
    mountainheight_dim = 400.0,
    mountainwidth_dim = 1000.0,
  )
  namelists = Namelists(; domain = domain, grid = grid)
  state = PinCFlow.State(namelists)
  @testset "Domain" begin
    @test state.domain.nx == 3
    @test state.domain.nz == 1
    @test state.domain.ny == 1
    @test state.domain.nxx == 10
    @test state.domain.nyy == 8
    @test state.domain.nzz == 8
  end
  @testset "Grid" begin
    @test state.grid.lx[0] ≈ 0.0
    @test state.grid.lx[1] ≈ 6.877891931902294
    @test state.grid.ly[0] ≈ 0.0
    @test state.grid.ly[1] ≈ 4.585261287934863
    @test state.grid.lz[0] ≈ 0.0
    @test state.grid.lz[1] ≈ 2.2926306439674313
    @test state.grid.dx ≈ 2.2926306439674313
    @test state.grid.dy ≈ 4.585261287934863
    @test state.grid.dz ≈ 2.2926306439674313
  end

  tol = 1e-14
  @testset "Constants" begin
    @test isapprox(state.constants.gamma, 1.4; atol = tol, rtol = tol)
    @test isapprox(state.constants.kappainv, 1.4 / 0.4; atol = tol, rtol = tol)
    @test isapprox(state.constants.gammainv, 1.0 / 1.4; atol = tol, rtol = tol)
    @test isapprox(state.constants.rsp, 287.0; atol = tol, rtol = tol)
    @test isapprox(state.constants.g, 9.81; atol = tol, rtol = tol)
    @test isapprox(state.constants.pref, 101325.0; atol = tol, rtol = tol)
    @test isapprox(state.constants.rhoref, 1.184, atol = 1e-14, rtol = tol)
    @test isapprox(
      state.constants.aref,
      292.5381125550948,
      atol = tol,
      rtol = tol,
    )
    @test isapprox(
      state.constants.uref,
      292.5381125550948,
      atol = tol,
      rtol = tol,
    )
    @test isapprox(
      state.constants.lref,
      8723.60319034631,
      atol = tol,
      rtol = tol,
    )
    @test isapprox(
      state.constants.tref,
      29.820398833342992,
      atol = tol,
      rtol = tol,
    )
    @test isapprox(
      state.constants.thetaref,
      298.1830916282137,
      atol = tol,
      rtol = tol,
    )
    @test isapprox(state.constants.ma, 1.0, atol = tol, rtol = tol)
    @test isapprox(state.constants.fr, 1.0, atol = tol, rtol = tol)
    @test isapprox(
      state.constants.kappa,
      0.28571428571428564,
      atol = tol,
      rtol = tol,
    )
    @test isapprox(state.constants.sig, 1.0, atol = tol, rtol = tol)
    @test isapprox(
      state.atmosphere.t0,
      1.0060932642487044,
      atol = tol,
      rtol = tol,
    )
    @test isapprox(
      state.atmosphere.n2,
      0.28398389678877484,
      atol = tol,
      rtol = tol,
    )
    @test isapprox(
      state.atmosphere.nn,
      0.5329013949960864,
      atol = tol,
      rtol = tol,
    )
    @test isapprox(state.constants.re, 1.0e20, atol = tol, rtol = tol)
    @test isapprox(
      state.constants.g_ndim,
      0.9999999999999998,
      atol = tol,
      rtol = tol,
    )
  end

  @testset "Jacobibian" begin
    @test state.grid.jac[2, 2, 2] ≈ 0.98
  end
  @testset "Metric tensor" begin
    @test state.grid.met[1, 1, 1, 1, 3] ≈ -0.004987779939149084
    @test state.grid.met[3, 1, 1, 2, 3] ≈ 0.0
    @test state.grid.met[1, 1, 1, 3, 3] ≈ 1.0001246360352993
  end
end