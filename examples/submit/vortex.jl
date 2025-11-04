# examples/submit/vortex.jl

using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

lx = 400000.0
ly = 400000.0

rx = lx / 4
ry = ly / 4

atmosphere = AtmosphereNamelist(;
    model = Boussinesq(),
    background = NeutralStratification(),
    initial_u = (x, y, z) -> begin
        r = sqrt((x / rx)^2 + (y / ry)^2)
        if r <= 1
            return -50 * y / ry * (1 + cos(pi * r)) / 2
        else
            return 0.0
        end
    end,
    initial_v = (x, y, z) -> begin
        r = sqrt((x / rx)^2 + (y / ry)^2)
        if r <= 1
            return 50 * x / rx * (1 + cos(pi * r)) / 2
        else
            return 0.0
        end
    end,
)

domain =
    DomainNamelist(; x_size = 40, y_size = 40, z_size = 1, lx, ly, npx, npy)

output = OutputNamelist(;
    output_variables = (:chi, :u, :v),
    output_file = "vortex.h5",
)

tracer = TracerNamelist(;
    tracer_setup = TracerOn(),
    initial_tracer = (x, y, z) -> begin
        r = sqrt(((abs(x) - rx) / rx)^2 + (y / ry)^2)
        if r <= 1
            return sign(x) * (1 + cos(pi * r)) / 2
        else
            return 0.0
        end
    end
)

integrate(Namelists(; atmosphere, domain, output, tracer))
