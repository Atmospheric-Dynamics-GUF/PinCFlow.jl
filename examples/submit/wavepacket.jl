# examples/submit/wavepacket.jl

using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

include("./wavepacket_ini.jl")

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npz = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

atmosphere = AtmosphereNamelist(;
    temperature = t0,
    ground_pressure = p0,
    coriolis_frequency = f,
    initial_u = uprime,
    initial_v = vprime,
    initial_w = wprime,
    initial_pip = piprime,
    initial_rhop = rhoprime,
    initial_thetap = thetaprime,
)
domain = DomainNamelist(;
    x_size = 32,
    y_size = 1,
    z_size = 854,
    lx = 30.0E+3,
    ly = 1.0E+3,
    lz = 80.E+3,
    npx,
    npz,
)
output = OutputNamelist(;
    output_variables = (:u, :w),
    output_file = "STIH_tke_with-coupling.h5",
    output_interval = 3600.0,
    tmax = 3600.0,
)
sponge = SpongeNamelist(;
    alpharmax = 1.0,
    damp_horizontal_wind_on_rhs = true,
    relax_to_mean = false,
)
turbulence = TurbulenceNamelist(;
    turbulence_scheme = TKEScheme(),
    momentum_coupling = true,
)

namelists = Namelists(; atmosphere, domain, output, sponge, turbulence)

integrate(namelists)
