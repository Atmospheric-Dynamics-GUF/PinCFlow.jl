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

x_size = 50 #x_size = 16
y_size = 1 #y_size = 16
z_size = 320 #z_size = 32

lx = 50000.0
ly = 50000.0
lz = 40000.0

wave_modes = 2


x_c = [0.0, 0.0]
y_c = [0.0, 0.0]
z_c = [16000.0, 22000.0 ]

sigma_xc = [lx / 6, lx / 6]
sigma_yc = [ly / 6, ly / 6]
sigma_zc = [lz / 20, lz / 20]

a0 = [0.5, 0.5]

k = [2 * pi / lx, 2 * pi / lx]
l = [2 * pi / ly, 2 * pi / ly]
m = [32 * pi / lz, -25 * pi / lz]

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
    initial_u = (x, y, z) -> backg_shear(x, y, z) + sum(real(uhat(alpha, x, y, z) * exp(1im * phi(alpha, x, y, z))) for alpha in alphas),
    initial_v = (x, y, z) -> sum(real(vhat(alpha, x, y, z) * exp(1im * phi(alpha, x, y, z))) for alpha in alphas),
    initial_w = (x, y, z) -> sum(real(what(alpha, x, y, z) * exp(1im * phi(alpha, x, y, z))) for alpha in alphas),
    initial_pip = (x, y, z) ->
        sum(real(pihat(alpha, x, y, z) * exp(1im * phi(alpha, x, y, z))) for alpha in alphas),
)
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)
output = OutputNamelist(;
    output_variables = (:u, :w),
    #save_ray_volumes = true,
    output_file = "resolved_wave_prop_shear_method_2.h5",
    tmax = 2000.0,
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

#Method-2 By adding the shear with the initial condition and hopping that shear and waves do not interact with each other


integrate(Namelists(; atmosphere, domain, output, wkb, triad))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("resolved_wave_prop_shear_method_2.h5") do data
        plot_output(
            "examples/results/resolved_wave_prop_shear_method_2.svg",
            data,
            ("w", 26, 1, 160, 1);
            time_unit = "min",
        )
        return
    end
end

