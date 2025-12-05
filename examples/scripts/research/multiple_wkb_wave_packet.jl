# examples/scripts/wkb_wave_packet.jl

using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
npz = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

x_size = 16 #x_size = 16
y_size = 16 #y_size = 16
z_size = 32 #z_size = 32

lx = 20000.0
ly = 20000.0
lz = 40000.0

wave_modes = 2


x_c = [0.0, 10000.0]
y_c = [0.0, 10000.0]
z_c = [30000.0, 20000.0 ]

sigma_xc = [lx / 6, lx / 6]
sigma_yc = [ly / 6, ly / 6]
sigma_zc = [lz / 12, lz / 12]

a0 = [0.05, 0.1]

k = [20 * pi / lx, 15 * pi / lx]
l = [20 * pi / ly, 15 * pi / ly]
m = [32 * pi / lz, -25 * pi / lz]

k_perp = sqrt.(k.^2 + l.^2)

lkp = 4 * max(abs.(k_perp)...)
lm = 4 * max(abs.(m)...)
kp_size = 16
m_size = 16



model = Compressible()
background = Realistic()
coriolis_frequency = 0.0001

atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)
domain = DomainNamelist(;
    x_size,
    y_size,
    z_size,
    lx,
    ly,
    lz,
    base_comm = MPI.COMM_SELF,
)

triad = TriadNamelist(;
    kp_size,
    m_size,
    lkp,
    lm,
    triad_int = true,
)
auxiliary_state = State(Namelists(; atmosphere, domain, triad))
(; g, kappa, rsp, lref, tref, rhoref, thetaref) = auxiliary_state.constants

include("multiple_wave_packet_tools.jl")

atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)
output = OutputNamelist(;
    output_variables = (:u, :w, :wavespectrum),
    save_ray_volumes = true,
    output_file = "multiple_wkb_wave_packet.h5",
    tmax = 600.0,
    output_interval = 30.0,
)
wkb = WKBNamelist(;
    wkb_mode = MultiColumn(),
    wave_modes,
    initial_wave_field = (alpha, x, y, z) ->
        (k[alpha], l[alpha], m[alpha], omega(alpha, x, y, z), wave_action_density(alpha, x, y, z)),
)

integrate(Namelists(; atmosphere, domain, output, wkb, triad))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("wkb_wave_packet.h5") do data
        plot_output(
            "examples/results/multiple_wkb_wave_packet.svg",
            data,
            ("nr", 8, 8, 16, 2);
            time_unit = "min",
        )
        return
    end
end
