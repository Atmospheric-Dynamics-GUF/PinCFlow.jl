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

x_size = 40
y_size = 40
z_size = 80

lx = 40000.0
ly = 40000.0
lz = 40000.0

rx = 0.0
ry = 0.5
rz = 0.1

sigma_x = 8000
sigma_y = 8000.0
sigma_z = 4000

x0 = 0.0
y0 = 0.0
z0 = lz / 2

k = (2 * pi * 10  / lx) 
l = (2 * pi * 10  / ly) 
m = (2 * pi * 20 / lz) 

#k_2 = (2 * pi/lx) 
#l_2 = (2 * pi/lx) 
#m_2 =  0.05
a0 = 0.05


#lz_1 = lz*1/2
#lz_2 = lz/5
#sigma_sq = 4.0E+3
#vertical_wind_k_1 = (x, y, z) -> amp * exp(-(z - lz_1)^2/(4*sigma_sq)) * cos(k_1*x + l_1*y + m_1*z)

#vertical_wind_k_2 = (x, y, z) -> amp * exp(-(z - lz_2)^2/(4*sigma_sq)) * cos(k_2*x + l_2*y + m_2*z)

#out_of_phase_w_k_1 = (x, y, z) -> amp * exp(-(z - lz_1)^2/(4*sigma_sq)) * sin(k_1*x + l_1*y + m_1*z)

#out_of_phase_w_k_2 = (x, y, z) -> amp * exp(-(z - lz_2)^2/(4*sigma_sq)) * sin(k_2*x + l_2*y + m_2*z)



#M = 2 * pi/1E+4
#z1 = 1.0E+4 
#N_Z = (x, y, z) -> z1 <= z && z <= z1 +  (2 * pi / M) ?
#               sqrt(buoyancy_f^2 + (buoyancy_f^2 * 0.8) * (sin(M * (z - z1)))) : buoyancy_f

#####



#atmosphere = AtmosphereNamelist(;
#    coriolis_frequency = 0.0E+0,
#    initial_rhop = (x, y, z) -> (-rho_o / gravity) * ((N_Z(x,y,z) * sqrt(k_1^2+l_1^2+ m_1^2) / sqrt(k_1^2+l_1^2)) * out_of_phase_w_k_1(x,y,z)
#                                + (N_Z(x,y,z) * sqrt(k_2^2+l_2^2+ m_2^2) / sqrt(k_2^2+l_2^2)) * out_of_phase_w_k_2(x,y,z)),
#    initial_u = (x, y, z) -> - (k_1 * m_1 / (k_1^2+l_1^2)) * vertical_wind_k_1(x,y,z) - (k_2 * m_2 / (k_2^2+l_2^2)) * vertical_wind_k_2(x,y,z),
#    initial_v = (x, y, z) -> - (l_1 * m_1 / (k_1^2+l_1^2)) * vertical_wind_k_1(x,y,z) - (l_2 * m_2 / (k_2^2+l_2^2)) * vertical_wind_k_2(x,y,z),
#    initial_w = (x, y, z) -> vertical_wind_k_1(x,y,z) + vertical_wind_k_2(x,y,z),
#    initial_pip = (x, y, z) -> - (1/c_p/theta_o) * (m_1 * N_Z(x,y,z) * vertical_wind_k_1(x,y,z) / (sqrt(k_1^2+l_1^2+ m_1^2) * sqrt(k_1^2+l_1^2))
#                                +m_2 * N_Z(x,y,z) * vertical_wind_k_2(x,y,z) / (sqrt(k_2^2+l_2^2+ m_2^2) * sqrt(k_2^2+l_2^2))) ,
    
    #initial_u = (x, y, z) -> amp * exp(-(z - lz_1)^2/(4*sigma_sq)) * cos(k_1*x + l_1*y + m_1*z)
    #            + amp * exp(-(z - lz_2)^2/(4*sigma_sq)) * cos(k_2*x + l_2*y + m_2*z),
    #initial_v = (x, y, z) -> -(k_1 * m_1 / (k_1^2+l_1^2)) * amp * exp(-(z - lz_1)^2/(4*sigma_sq)) * cos(k_1*x + l_1*y + m_1*z)
    #            -(k_2 * m_2 / (k_2^2+l_2^2)) * amp * exp(-(z - lz_2)^2/(4*sigma_sq)) * cos(k_2*x + l_2*y + m_2*z),
    #initial_w = (x, y, z) -> -(l_1 * m_1 / (k_1^2+l_1^2)) * amp * exp(-(z - lz_1)^2/(4*sigma_sq)) * cos(k_1*x + l_1*y + m_1*z)
    #            -(l_2 * m_2 / (k_2^2+l_2^2)) * amp * exp(-(z - lz_2)^2/(4*sigma_sq)) * cos(k_2*x + l_2*y + m_2*z),
#    model = Boussinesq(),
#    background = RadiatedBoussinesq(),
#)


background = Realistic()
coriolis_frequency = 0.0

atmosphere = AtmosphereNamelist(; background, coriolis_frequency)
domain = DomainNamelist(;
    x_size,
    y_size,
    z_size,
    lx,
    ly,
    lz,
    base_comm = MPI.COMM_SELF,
)
auxiliary_state = State(Namelists(; atmosphere, domain))
(; g, kappa, rsp) = auxiliary_state.constants

include("wave_packet_tools.jl")

atmosphere = AtmosphereNamelist(;
    background,
    coriolis_frequency,
    initial_rhop = (x, y, z) ->
        rhobar(x, y, z) *
        (1 / (1 + real(bhat(x, y, z) * exp(1im * phi(x, y, z))) / g) - 1),
    initial_u = (x, y, z) -> real(uhat(x, y, z) * exp(1im * phi(x, y, z))),
    #initial_v = (x, y, z) -> real(vhat(x, y, z) * exp(1im * phi(x, y, z))),
    initial_w = (x, y, z) -> real(what(x, y, z) * exp(1im * phi(x, y, z))),
    initial_pip = (x, y, z) ->
        real(pihat(x, y, z) * exp(1im * phi(x, y, z))),
)
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)
output = OutputNamelist(;
    output_variables = (:u, :v, :w),
    output_file = "wave_packet1.h5",
    tmax = 50.0,
    output_interval = 10.0,
)

integrate(Namelists(; atmosphere, domain, output))

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    h5open("wave_packet1.h5") do data
        plot_output(
            "examples/results/wave_packet1.svg",
            data,
            ("u", 20, 20, 40, 2),
            ("v", 20, 20, 40, 2),
            ("w", 20, 20, 40, 2);
            time_unit = "min",
        )
        return
    end
end
