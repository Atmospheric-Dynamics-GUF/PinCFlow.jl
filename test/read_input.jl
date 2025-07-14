using Random
using Test
using Dates

domain = DomainNamelist(;
    sizex = 2,
    sizey = 3,
    sizez = 4,
    nbx = 3,
    nby = 3,
    nbz = 3,
    lx_dim = (0.0, 4.0E+5),
    ly_dim = (0.0, 4.0E+5),
    lz_dim = (0.0, 2.0E+4),
    npx = 1,
    npy = 1,
)
output = OutputNamelist(;
    output_variables = (),
    prepare_restart = true,
    restart = true,
    iin = -1,
    output_steps = true,
    noutput = 1,
    maxiter = 1,
    outputtimediff = 3.6E+3,
    maxtime = 3.6E+3,
    input_file = "./test/pincflow_input.h5",
    output_file = "./test/pincflow_output.h5",
)
setting = SettingNamelist(;
    model = PseudoIncompressible(),
    testcase = WKBMountainWave(),
    zboundaries = SolidWallBoundaries(),
)
discretization = DiscretizationNamelist()
poisson = PoissonNamelist()
atmosphere = AtmosphereNamelist()
grid = GridNamelist()
sponge = SpongeNamelist()
wkb = WKBNamelist(;
    xrmin_dim = 0.0E+0,
    xrmax_dim = 4.0E+5,
    yrmin_dim = 0.0,
    yrmax_dim = 4.0E+5,
    zrmin_dim = 0.0,
    zrmax_dim = 2.0E+4,
    nrxl = 1,
    nryl = 1,
    nrzl = 1,
    nrk_init = 1,
    nrl_init = 1,
    nrm_init = 1,
    nray_fac = 4,
    fac_dk_init = 1.0E-1,
    fac_dl_init = 1.0E-1,
    fac_dm_init = 1.0E-1,
    branchr = -1,
    merge_mode = ConstantWaveAction(),
    nsmth_wkb = 2,
    lsmth_wkb = true,
    sm_filter = Shapiro(),
    zmin_wkb_dim = 0.0E+0,
    lsaturation = true,
    alpha_sat = 1.0E+0,
    wkb_mode = MultiColumn(),
    blocking = false,
    nwm = 1,
)

namelists = Namelists(;
    domain = domain,
    output = output,
    setting = setting,
    discretization = discretization,
    poisson = poisson,
    atmosphere = atmosphere,
    grid = grid,
    sponge = sponge,
    wkb = wkb,
)

state = State(namelists)
(; k0, k1, i0, i1, j0, j1) = state.domain
(; sizex, sizey, sizez) = state.namelists.domain

cpu_start_time = now()
PinCFlow.Output.create_output(state)

state.variables.predictands.u .=
    rand(Float32, size(state.variables.predictands.u))
state.variables.predictands.v .=
    rand(Float32, size(state.variables.predictands.v))
state.variables.predictands.w .=
    rand(Float32, size(state.variables.predictands.w))

# set rays
#
nray = 5
state.wkb.rays.dens .= 0.0
state.wkb.rays.x .= zeros(size(state.wkb.rays.x))
# make sure that rays.dens has some "holes", i.e. zero values
for kz in k0:k1, jy in j0:j1, ix in i0:i1
    indices = randperm(state.wkb.nray_max)[1:nray]
    state.wkb.rays.dens[indices, ix, jy, kz] .= 1.0
    # mark some x values so that we can look for them later
    state.wkb.rays.x[indices, ix, jy, kz] .= 1.0
end

dens_sum = sum(state.wkb.rays.dens)
@assert isapprox(dens_sum, sizex * sizey * sizez * nray)

PinCFlow.Output.write_output(state, 0.0, 1, cpu_start_time)

restart_output = OutputNamelist(;
    output_variables = (),
    restart = true,
    iin = 2,
    input_file = "./test/pincflow_output.h5",
    output_file = "./test/pincflow_output.h5",
)
in_state = PinCFlow.State(
    Namelists(;
        domain = domain,
        output = restart_output,
        setting = setting,
        discretization = discretization,
        poisson = poisson,
        atmosphere = atmosphere,
        grid = grid,
        sponge = sponge,
        wkb = wkb,
    ),
)

t = PinCFlow.Output.read_input!(in_state)

# check that fields match **inside** grid
(; i0, i1, j0, j1, k0, k1) = in_state.domain
@test all(
    isapprox.(
        in_state.variables.predictands.u[i0:i1, j0:j1, k0:k1],
        state.variables.predictands.u[i0:i1, j0:j1, k0:k1],
        rtol = 1e-6,
    ),
)
@test all(
    isapprox.(
        in_state.variables.predictands.v[i0:i1, j0:j1, k0:k1],
        state.variables.predictands.v[i0:i1, j0:j1, k0:k1],
        rtol = 1e-6,
    ),
)
@test all(
    isapprox.(
        in_state.variables.predictands.w[i0:i1, j0:j1, k0:k1],
        state.variables.predictands.w[i0:i1, j0:j1, k0:k1],
        rtol = 1e-6,
    ),
)

@test isapprox(dens_sum, sum(in_state.wkb.rays.dens), rtol = 1e-6)
@test isapprox(sum(in_state.wkb.nray), dens_sum)

# find the marked x values and test that the shuffeling is correct
marked_x_values = findall(isapprox.(in_state.wkb.rays.x, 1.0, rtol = 1e-6))
foreach(marked_x_values) do i
    @test isapprox(in_state.wkb.rays.dens[i], 1.0, rtol = 1e-6)
end

rm("./test/pincflow_output.h5")
