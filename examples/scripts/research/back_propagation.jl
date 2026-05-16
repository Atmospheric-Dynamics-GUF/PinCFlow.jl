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

nthreads_triad = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : Threads.nthreads()


#model parameters

x_size = 2 #x_size = 16
y_size = 1 #y_size = 16
z_size = 115 #z_size = 176 #z_size = 32

lx = 50000.0
ly = 50000.0
lz = 40000.0

wave_modes = 1

#Third set of initial conditions
x_c = [0.0, 0.0]
y_c = [0.0, 0.0]
z_c = [16000.0, 16000.0]

sigma_xc = [lx / 2, lx / 2]
sigma_yc = [1.0, 1.0]
sigma_zc = [3000.0, 3000.0]

#test case: 1
a0 = [0.1, 0.05]
#test case:2
#a0 = [0.05, 0.05]


k = [2 * pi / 25000.0, 2 * pi / 25000.0]
l = [0.0, 0.0]
#m = [5 * 2 * pi / 5000.0, -2.07 * 2 * pi / 5000.0]
m = [0.005363034122668976]
#m =  [-0.0077385419819957415]

k_max = 8.0 * max(abs.(k)...)
m_max = 3.5 * max(abs.(m)...)
k_min = 0.25 * (2π / lx)
k_size = 22
m_size = 22


# defining the shear
shear_s = 0.2
#backg_shear = (x, y, z) -> shear_s * sin(2 * pi * z / lz)
backg_shear = (x, y, z) -> shear_s


model = Boussinesq()
background = RadiatedBoussinesq()
#background = StableStratification()
coriolis_frequency = 0.0000

atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)
domain = DomainNamelist(;
    x_size,
    y_size,
    z_size,
    lx,
    ly,
    lz,
    npx,
    npy,
    npz,
    base_comm = MPI.COMM_SELF,
)

triad = TriadNamelist(;
    k_size,
    m_size,
    k_max,
    m_max,
    k_min,
    triad_mode = NoTriad(),
    time_scheme = EulerMethod(),
    rm_index = (7, 8),
    nthreads_triad,
)

wkb = WKBNamelist(;
    wkb_mode = SingleColumn(),
    branch = +1,
    impact_altitude = 0.0,
    #use_saturation = false,
)
auxiliary_state = State(Namelists(; atmosphere, domain, wkb, triad))
(; g, kappa, rsp, lref, tref, rhoref, thetaref) = auxiliary_state.constants

include("multiple_wave_packet_tools.jl")

atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency,
initial_u = (x, y, z) -> backg_shear(x, y, z))
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)
output = OutputNamelist(;
    output_variables = (:u, :w, :rhop, :thetap, :wavespectrum),
    save_ray_volumes = true,
    output_steps = true,
    nout = 10,  #5,
    iterations = 300,   #130,
    output_file = "back_propagation_pos.h5",
    tmax = 80,
    output_interval = 20.0,
)
wkb = WKBNamelist(;
    wkb_mode = SingleColumn(),
    branch = +1,
    impact_altitude = 0.0,
    wave_modes,
    multiplication_factor = 8,
    initial_wave_field = (alpha, x, y, z) ->
        (k[alpha], l[alpha], m[alpha], omega(alpha, x, y, z), wave_action_density(alpha, x, y, z)),
)

integrate(Namelists(; atmosphere, domain, output, wkb, triad))