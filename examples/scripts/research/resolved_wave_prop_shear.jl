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
z_c = [21937.755, 16851.213]

sigma_xc = [lx, lx]
sigma_yc = [1.0, 1.0]
sigma_zc = [lz / 20, lz / 20]

a0 = [0.5, 0.5]

k = [2 * pi / lx, 2 * pi / lx]
l = [0.0, 0.0]
#m = [5 * 2 * pi / 5000.0, -2.07 * 2 * pi / 5000.0]
m = [0.0043791985000000005, -0.0033200091999999997]

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
    initial_u = (x, y, z) -> sum(real(uhat(alpha, x, y, z) * exp(1im * phi(alpha, x, y, z))) for alpha in alphas),
    initial_v = (x, y, z) -> sum(real(vhat(alpha, x, y, z) * exp(1im * phi(alpha, x, y, z))) for alpha in alphas),
    initial_w = (x, y, z) -> sum(real(what(alpha, x, y, z) * exp(1im * phi(alpha, x, y, z))) for alpha in alphas),
    initial_pip = (x, y, z) ->
        sum(real(pihat(alpha, x, y, z) * exp(1im * phi(alpha, x, y, z))) for alpha in alphas),
)
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)
output = OutputNamelist(;
    output_variables = (:u, :w),
    #save_ray_volumes = true,
    output_file = "resolved_wave_prop_shear.h5",
    tmax = 5000.0,
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
    h5open("resolved_wave_prop_shear.h5") do data
        plot_output(
            "examples/results/resolved_wave_prop_shear.svg",
            data,
            ("w", 26, 1, 160, 1);
            time_unit = "min",
        )
        return
    end
end

