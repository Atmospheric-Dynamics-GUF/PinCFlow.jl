using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npz = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

#local parameteres for wavepackets
##natural parametere
rho_o = 1.0
gravity = 9.8
buoyancy_f = 0.01
c_p = 1005.0 #doubtfull
theta_o = 3.0E+2 #/thetaref
#####
lx = 5.0E+4
ly = 5.0E+4
lz = 4.0E+4

k_1 = (2 * pi*20/lx)   #in m^{-1}
l_1 = (2 * pi*20/lx) 
m_1 = 0.05 
amp = 10
lz_1 = lz/2
sigma_sq = 10.0E+3

vertical_wind_k_1 = (x, y, z) -> amp * exp(-(z - lz_1)^2/(4*sigma_sq)) 

vertical_wind_k_2 = (x, y, z) -> amp * exp(-(z - lz_2)^2/(4*sigma_sq)) 

out_of_phase_w_k_1 = (x, y, z) -> amp * exp(-(z - lz_1)^2/(4*sigma_sq)) 

out_of_phase_w_k_2 = (x, y, z) -> amp * exp(-(z - lz_2)^2/(4*sigma_sq)) 

M = 2 * pi/1E+4
z1 = 1.0E+4 
N_Z = (x, y, z) -> z1 <= z && z <= z1 +  (2 * pi / M) ?
               #sqrt(buoyancy_f^2 + (buoyancy_f^2 * 0.8) * (sin(M * (z - z1)))) 
               buoyancy_f : buoyancy_f

wave_density = (alpha, x, y, z) -> lz_1 - 4*sigma_sq <= z && z <= lz_1 + 4*sigma_sq ?
              (1/2) * (k_1^2 + l_1^2 + m_1^2)^(3/2) * vertical_wind_k_1(x,y,z)^2 /(N_Z(x, y, z)*(k_1^2 + l_1^2)^(3/2))  : 0


atmosphere = AtmosphereNamelist(;
    coriolis_frequency = 0.0E+0,
    #initial_rhop = (x, y, z) -> (-rho_o / gravity) * ((N_Z(x,y,z) * sqrt(k_1^2+l_1^2+ m_1^2) / sqrt(k_1^2+l_1^2)) * out_of_phase_w_k_1(x,y,z)),
    #                            + (N_Z(x,y,z) * sqrt(k_2^2+l_2^2+ m_2^2) / sqrt(k_2^2+l_2^2)) * out_of_phase_w_k_2(x,y,z)),
    initial_u = (x, y, z) -> - (k_1 * m_1 / (k_1^2+l_1^2)) * vertical_wind_k_1(x,y,z), 
                            #- (k_2 * m_2 / (k_2^2+l_2^2)) * vertical_wind_k_2(x,y,z),
    initial_v = (x, y, z) -> - (l_1 * m_1 / (k_1^2+l_1^2)) * vertical_wind_k_1(x,y,z), 
                                #- (l_2 * m_2 / (k_2^2+l_2^2)) * vertical_wind_k_2(x,y,z),
    #initial_w = (x, y, z) -> vertical_wind_k_1(x,y,z), 
                        #+ vertical_wind_k_2(x,y,z),
    #initial_pip = (x, y, z) -> - (1/c_p/theta_o) * (m_1 * N_Z(x,y,z) * vertical_wind_k_1(x,y,z) / (sqrt(k_1^2+l_1^2+ m_1^2) * sqrt(k_1^2+l_1^2))),
    #                            +m_2 * N_Z(x,y,z) * vertical_wind_k_2(x,y,z) / (sqrt(k_2^2+l_2^2+ m_2^2) * sqrt(k_2^2+l_2^2))) ,
    model = Boussinesq(),
    background = RadiatedBoussinesq(),
)
domain = DomainNamelist(;
    x_size = 40,
    y_size = 40,
    z_size = 40,
    lx,
    ly,
    lz,
    npx,
    npz,
)
grid = GridNamelist(;)


sponge = SpongeNamelist(;
    #sponge_extent = 1.0E-1,
    #alpharmax = 1.79E-2,
    #lateral_sponge = true,
    #sponge_type = ExponentialSponge(),
    #relax_to_mean = false,
    #relaxation_wind = (1.0E+1, 0.0E+0, 0.0E+0),
)
output =
    OutputNamelist(; output_variables = (:u, :dudt,), #:dudt 
    output_file = "wkb_prop_perturb.h5", 
    output_interval = 1.0E+2,
    save_ray_volumes = true,
    tmax = 1.0E+3,
)
sponge = SpongeNamelist(; 
#betarmax = 1.0E+0
)

wkb = WKBNamelist(; wkb_mode = MultiColumn(),
        wave_modes = 1,
        initial_wave_field = (alpha, x, y, z) -> lz_1 - 4*sigma_sq <= z && z <= lz_1 + 4*sigma_sq ?
         (k_1, l_1, m_1, 
        N_Z(x, y, z) * sqrt(k_1^2+l_1^2)/sqrt(k_1^2+l_1^2+m_1^2), wave_density(alpha, x, y, z)
        ) : (0.0, 0.0, 0.0, 0.0, 0.0),
)
integrate(Namelists(; atmosphere, domain, grid, output, sponge, wkb))
