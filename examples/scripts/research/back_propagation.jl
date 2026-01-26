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

#model parameters

x_size = 2 #x_size = 16
y_size = 1 #y_size = 16
z_size = 229 #z_size = 32

lx = 50000.0
ly = 50000.0
lz = 40000.0

wave_modes = 2

x_c = [0.0, 0.0]
y_c = [0.0, 0.0]
z_c = [20000.0, 20000.0]

sigma_xc = [lx, lx]
sigma_yc = [1.0, 1.0]
sigma_zc = [1000.0, 1000.0]

a0 = [0.5, 0.5]

k = [4 * pi / 50000.0, 4 * pi / 50000.0]
l = [0.0, 0.0]
#m = [5 * 2 * pi / 5000.0, -2.07 * 2 * pi / 5000.0]
#postive wave number
#m = [0.03575356081779317]

#negative wave number
m = [0.03575356081779317, -0.014809609793861218]

"""
x_c = [0.0]
y_c = [0.0]
z_c = [20000.0]

sigma_xc = [lx]
sigma_yc = [ly]
sigma_zc = [lz / 20]

a0 = [0.05]

k = [2 * pi / lx]
l = [0.0]
m = [-2.07 * 2 * pi / 5000.0]

"""

k_perp = sqrt.(k.^2 + l.^2)

lkp = 4 * max(abs.(k_perp)...)
lm = 4 * max(abs.(m)...)
k_size = 16
m_size = 32

lk = 4 * max(abs.(k)...)

# defining the shear
shear_s = 1.0
backg_shear = (x, y, z) -> shear_s * sin(2 * pi * z / lz)


model = Boussinesq()
background = StableStratification()
coriolis_frequency = 0.0000

atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)
domain = DomainNamelist(;
    x_size,
    y_size,
    z_size,
    lx,
    ly,
    lz,
    npx = 1,
    npy = 1,
    npz = 1,
    base_comm = MPI.COMM_SELF,
)

triad = TriadNamelist(;
    k_size,
    m_size,
    lk,
    lm,

    triad_mode = NoTriad(),
)

wkb = WKBNamelist(;
    wkb_mode = SingleColumn(),
    branch = +1,
    impact_altitude = 0.0,
)
auxiliary_state = State(Namelists(; atmosphere, domain, wkb, triad))
(; g, kappa, rsp, lref, tref, rhoref, thetaref) = auxiliary_state.constants

include("multiple_wave_packet_tools.jl")

atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency,
initial_u = (x, y, z) -> backg_shear(x, y, z))
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)
output = OutputNamelist(;
    output_variables = (:u, :w, :wavespectrum),
    save_ray_volumes = true,
    output_steps = true,
    nout = 10,
    iterations = 700,
    output_file = "back_prop_pos_neg.h5",
    tmax = 1000.0,
    output_interval = 100.0,
)
wkb = WKBNamelist(;
    wkb_mode = SingleColumn(),
    branch = +1,
    impact_altitude = 0.0,
    wave_modes,
    initial_wave_field = (alpha, x, y, z) ->
        (k[alpha], l[alpha], m[alpha], omega(alpha, x, y, z), wave_action_density(alpha, x, y, z)),
)

integrate(Namelists(; atmosphere, domain, output, wkb, triad))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("back_prop_pos_neg.h5") do data
        plot_output(
            "examples/results/back_prop_pos_neg.svg",
            data,
            ("nr", 2, 1, 67, 71);
            time_unit = "min",
        )
        return
    end
end
