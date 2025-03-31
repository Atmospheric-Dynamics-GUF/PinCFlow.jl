@testset "Boundary" begin
  domain = PinCFlow.DomainNamelist(;
    sizex = 30,
    sizey = 30,
    sizez = 10,
    lx_dim = [0.0, 60000.0],
    ly_dim = [0.0, 40000.0],
    lz_dim = [0.0, 20000],
  )
  atmosphere =
    PinCFlow.AtmosphereNamelist(; backgroundflow_dim = [10.0, 0.0, 0.0])
  namelists = PinCFlow.Namelists(; domain = domain, atmosphere = atmosphere)
  state = PinCFlow.State(namelists)

  var = state.variables.predictands

  @testset "Zonal" begin
    set_lin_array!(var.u)

    @test var.u[31, :, :] != var.u[1, :, :]
    @test var.u[32, :, :] != var.u[2, :, :]
    @test var.u[33, :, :] != var.u[3, :, :]

    PinCFlow.set_zonal_boundaries_of_field!(var.u, namelists, state.domain)

    @test var.u[31, :, :] == var.u[1, :, :]
    @test var.u[32, :, :] == var.u[2, :, :]
    @test var.u[33, :, :] == var.u[3, :, :]

    @test var.u[0, :, :] == var.u[30, :, :]
    @test var.u[-1, :, :] == var.u[29, :, :]
    @test var.u[-2, :, :] == var.u[28, :, :]
  end

  @testset "Meridional" begin
    set_lin_array!(var.v)

    PinCFlow.set_meridional_boundaries_of_field!(var.v, namelists, state.domain)

    @test var.v[:, 31, :] == var.v[:, 1, :]
    @test var.v[:, 32, :] == var.v[:, 2, :]
    @test var.v[:, 33, :] == var.v[:, 3, :]

    @test var.v[:, 0, :] == var.v[:, 30, :]
    @test var.v[:, -1, :] == var.v[:, 29, :]
    @test var.v[:, -2, :] == var.v[:, 28, :]
  end

  @testset "Vertical" begin
    set_lin_array!(var.w)
    PinCFlow.set_vertical_boundaries!(
      state,
      PinCFlow.BoundaryPredictands(),
      PinCFlow.SolidWallBoundaries(),
    )

    @test all(var.w[:, :, 0] .== 0)
    @test all(var.w[:, :, 10] .== 0)

    @test_broken all(var.v[:, :, 1] .== 0)
    # TODO: is this excpected? 
    @test var.w[:, :, -1] == -var.w[:, :, 1]
    @test var.w[:, :, -2] == -var.w[:, :, 2]
    @test var.w[:, :, -3] == -var.w[:, :, 3]
    @test var.w[:, :, 11] == -var.w[:, :, 9]
    @test var.w[:, :, 12] == -var.w[:, :, 8]
    @test var.w[:, :, 13] == -var.w[:, :, 7]
  end
end
