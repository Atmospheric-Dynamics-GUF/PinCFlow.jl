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
z_size = 120 #z_size = 176 #z_size = 32

lx = 50000.0
ly = 50000.0
lz = 40000.0

wave_modes = 2

#First set of initial conditi
x_c = [0.0, 0.0]
y_c = [0.0, 0.0]
z_c = [20769.375, 18079.259765625]

sigma_xc = [lx / 2, lx / 2]
sigma_yc = [1.0, 1.0]
sigma_zc = [lz / 40, lz / 40]

a0 = [0.8, 0.8]

k = [2 * pi / 25000.0, 2 * pi / 25000.0]
l = [0.0, 0.0]
#m = [5 * 2 * pi / 5000.0, -2.07 * 2 * pi / 5000.0]
m = [0.013695018365979195, -0.011669567786157131]

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

lk = 3.5 * max(abs.(k)...)
lm = 2.5 * max(abs.(m)...)
k_size = 32
m_size = 42



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
    base_comm = MPI.COMM_SELF,
)

triad = TriadNamelist(;
    k_size,
    m_size,
    lk,
    lm,
    triad_mode = Triad2D(),
    time_scheme = EulerMethod(),
)

wkb = WKBNamelist(;
    wkb_mode = SingleColumn(),
    branch = +1,
    impact_altitude = 0.0,
)
auxiliary_state = State(Namelists(; atmosphere, domain, wkb, triad))
(; g, kappa, rsp, lref, tref, rhoref, thetaref) = auxiliary_state.constants

include("multiple_hyd_wave_packet_tools.jl")

atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency,
initial_u = (x, y, z) -> backg_shear(x, y, z))
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)
output = OutputNamelist(;
    output_variables = (:u, :w, :wavespectrum),
    save_ray_volumes = true,
    output_steps = true,
    nout = 3,  #5,
    iterations = 250,   #130,
    output_file = "triad_wave_propagation_250_hyds.h5",
    tmax = 80,
    output_interval = 20.0,
)
wkb = WKBNamelist(;
    wkb_mode = SingleColumn(),
    branch = +1,
    impact_altitude = 0.0,
    wave_modes,
    multiplication_factor = 6,
    
    initial_wave_field = (alpha, x, y, z) ->
        (k[alpha], l[alpha], m[alpha], omega(alpha, x, y, z), wave_action_density(alpha, x, y, z)),
)

integrate(Namelists(; atmosphere, domain, output, wkb, triad))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("triad_wave_propagation_250_hyds.h5") do data
        plot_output(
            "examples/results/triad_wave_propagation_250_hyds.svg",
            data,
            ("nr", 2, 1, 67, 1);
            time_unit = "min",
        )
        return
    end
end
