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
nbx = 2
nby = 1
nbz = 3


#model parameters

x_size = 80 #x_size = 16
y_size = 1 #y_size = 16
z_size = 2400 #z_size = 32

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


# defining the shear
shear_s = 1.0
backg_shear = (x, y, z) -> shear_s * sin(2 * pi * z / lz)

model = Boussinesq()
background = StableStratification()
coriolis_frequency = 0.0

####End of model parameter############





atmosphere = AtmosphereNamelist(; background, model, coriolis_frequency)

#auxiliary atmosphere with shear
atmosphere_1 = AtmosphereNamelist(;
    coriolis_frequency,
    initial_u = (x, y, z) -> backg_shear(x, y, z),
    model,
    background,
)

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
    triad_mode = NoTriad(),
)

wkb = WKBNamelist(;
    wkb_mode = NoWKB(),
    branch = +1,
    impact_altitude = 0.0,
)
auxiliary_state = State(Namelists(; atmosphere, domain, wkb, triad))
(; g, kappa, rsp, lref, tref, rhoref, thetaref) = auxiliary_state.constants

include("multiple_wave_packet_tools.jl")

alphas = [1, 2]

atmosphere = AtmosphereNamelist(;
    background,
    model,
    coriolis_frequency,
    initial_rhop = (x, y, z) ->
        sum(rhobar(x, y, z) *
        (1 / (1 + real(bhat(alpha, x, y, z) * exp(1im * phi(alpha, x, y, z))) / g) - 1) for alpha in alphas),
    initial_u = (x, y, z) -> sum(real(uhat(alpha, x, y, z) * exp(1im * phi(alpha, x, y, z))) for alpha in alphas),
    initial_v = (x, y, z) -> sum(real(vhat(alpha, x, y, z) * exp(1im * phi(alpha, x, y, z))) for alpha in alphas),
    initial_w = (x, y, z) -> sum(real(what(alpha, x, y, z) * exp(1im * phi(alpha, x, y, z))) for alpha in alphas),
    initial_pip = (x, y, z) ->
        sum(real(pihat(alpha, x, y, z) * exp(1im * phi(alpha, x, y, z))) for alpha in alphas),
)
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)
output = OutputNamelist(;
    output_variables = (:u, :w, :rhop, :thetap),
    #save_ray_volumes = true,
    output_steps = true,
    nout = 3,  #5,
    iterations = 400, 
    output_file = "resolved_flow_test_1.h5",
    tmax = 100.0,
    output_interval = 50.0,
)
wkb = WKBNamelist(;
    wkb_mode = NoWKB(),
    branch = +1,
    impact_altitude = 0.0,
    wave_modes,
    initial_wave_field = (alpha, x, y, z) ->
        (k[alpha], l[alpha], m[alpha], omega(alpha, x, y, z), wave_action_density(alpha, x, y, z)),
)

#Method-1 By creating the auxiliary state for the the shear and using it to only calculate the fluxes
#creating the auxiliary state with shear
(state1) = State(Namelists(; atmosphere = atmosphere_1, domain, output, wkb, triad))
shear1 = deepcopy(state1.variables.predictands.u)


#calling the integrate function suitable for the propagation in shear
integrate(Namelists(; atmosphere, domain, output, wkb, triad), shear1)

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("resolved_flow_test_1.h5") do data
        plot_output(
            "examples/results/resolved_flow_test_1.svg",
            data,
            ("w", 26, 1, 160, 1);
            time_unit = "min",
        )
        return
    end
end

