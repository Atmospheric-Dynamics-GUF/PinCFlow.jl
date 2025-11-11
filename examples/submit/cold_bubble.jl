# examples/submit/cold_bubble.jl

using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npz = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

lx = 20000.0
lz = 20000.0

rx = lx / 8
rz = lz / 8

atmosphere = AtmosphereNamelist(;
    background = Isentropic(),
    initial_rhop = (x, y, z) -> begin
        r = sqrt((x / rx)^2 + ((z - 3 * rz) / rz)^2)
        if r <= 1
            return 0.005 * (1 + cos(pi * r))
        else
            return 0.0
        end
    end,
)
discretization = DiscretizationNamelist(; dtmax = 60.0)
domain = DomainNamelist(; x_size = 40, z_size = 40, lx, lz, npx, npz)
output = OutputNamelist(;
    output_variables = (:thetap,),
    output_file = "cold_bubble.h5",
)

integrate(Namelists(; atmosphere, discretization, domain, output))
