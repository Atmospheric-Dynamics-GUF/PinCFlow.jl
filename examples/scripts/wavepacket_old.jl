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
    model = PseudoIncompressible(),
) 
domain = DomainNamelist(;
    x_size = 16,
    y_size = 1,
    z_size = 427,
    lx = 30.0E+3,
    ly = 1.0E+3,
    lz = 80.E+3,
    npx,
    npz,
)
output = OutputNamelist(;
    output_variables = (:u, :w),
    output_file = "wavepacket.h5",
    output_interval = 360.0,
    tmax = 3600.0,
)

turbulence = TurbulenceNamelist(;
    turbulence_scheme = TKEScheme(),
)
tracer = TracerNamelist(;
    tracer_setup = TracerOn(),
    initial_tracer = chitotal)

namelists = Namelists(; atmosphere, domain, output, turbulence, tracer)

integrate(namelists)
