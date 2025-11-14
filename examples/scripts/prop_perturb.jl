using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npz = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

#local parameteres for wavepackets
##natural parametere
# Set natural constants.
    gamma = 1.4
    gammainv = 1.0 / gamma
    kappa = (gamma - 1.0) / gamma
    kappainv = 1.0 / kappa
    rsp = 287.0
    g = 9.81

    # Set reference quantities.
    rhoref = 1.184 # in kg/m^3
    pref = 101325.0 # in Pa = kg/m/s^2
    aref = sqrt(pref / rhoref) # in m/s
    uref = aref # in m/s
    lref = pref / rhoref / g # in m
    tref = lref / aref # in s
    thetaref = aref^2 / rsp # in K
    fref = rhoref * uref^2 / lref # in N/m^3
rho_o = 1.0
gravity = 9.8
buoyancy_f = 0.01
c_p = 1005.0 #doubtfull
theta_o = 3.0E+2/thetaref
 
#####
lx = 5.0E+4
ly = 5.0E+4
lz = 4.0E+4

k_1 = (2 * pi*10/lx) 
l_1 = (2 * pi*10/lx) 
m_1 = -0.0297
k_2 = (2 * pi/lx) 
l_2 = (2 * pi/lx) 
m_2 =  0.05
amp = 0.5
lz_1 = lz*1/2
lz_2 = lz/5
sigma_sq = 4.0E+3
vertical_wind_k_1 = (x, y, z) -> amp * exp(-(z - lz_1)^2/(4*sigma_sq)) * cos(k_1*x + l_1*y + m_1*z)

vertical_wind_k_2 = (x, y, z) -> amp * exp(-(z - lz_2)^2/(4*sigma_sq)) * cos(k_2*x + l_2*y + m_2*z)

out_of_phase_w_k_1 = (x, y, z) -> amp * exp(-(z - lz_1)^2/(4*sigma_sq)) * sin(k_1*x + l_1*y + m_1*z)

out_of_phase_w_k_2 = (x, y, z) -> amp * exp(-(z - lz_2)^2/(4*sigma_sq)) * sin(k_2*x + l_2*y + m_2*z)



M = 2 * pi/1E+4
z1 = 1.0E+4 
N_Z = (x, y, z) -> z1 <= z && z <= z1 +  (2 * pi / M) ?
               sqrt(buoyancy_f^2 + (buoyancy_f^2 * 0.8) * (sin(M * (z - z1)))) : buoyancy_f

#####



atmosphere = AtmosphereNamelist(;
    coriolis_frequency = 0.0E+0,
    initial_rhop = (x, y, z) -> (-rho_o / gravity) * ((N_Z(x,y,z) * sqrt(k_1^2+l_1^2+ m_1^2) / sqrt(k_1^2+l_1^2)) * out_of_phase_w_k_1(x,y,z)
                                + (N_Z(x,y,z) * sqrt(k_2^2+l_2^2+ m_2^2) / sqrt(k_2^2+l_2^2)) * out_of_phase_w_k_2(x,y,z)),
    initial_u = (x, y, z) -> - (k_1 * m_1 / (k_1^2+l_1^2)) * vertical_wind_k_1(x,y,z) - (k_2 * m_2 / (k_2^2+l_2^2)) * vertical_wind_k_2(x,y,z),
    initial_v = (x, y, z) -> - (l_1 * m_1 / (k_1^2+l_1^2)) * vertical_wind_k_1(x,y,z) - (l_2 * m_2 / (k_2^2+l_2^2)) * vertical_wind_k_2(x,y,z),
    initial_w = (x, y, z) -> vertical_wind_k_1(x,y,z) + vertical_wind_k_2(x,y,z),
    initial_pip = (x, y, z) -> - (1/c_p/theta_o) * (m_1 * N_Z(x,y,z) * vertical_wind_k_1(x,y,z) / (sqrt(k_1^2+l_1^2+ m_1^2) * sqrt(k_1^2+l_1^2))
                                +m_2 * N_Z(x,y,z) * vertical_wind_k_2(x,y,z) / (sqrt(k_2^2+l_2^2+ m_2^2) * sqrt(k_2^2+l_2^2))) ,
    
    #initial_u = (x, y, z) -> amp * exp(-(z - lz_1)^2/(4*sigma_sq)) * cos(k_1*x + l_1*y + m_1*z)
    #            + amp * exp(-(z - lz_2)^2/(4*sigma_sq)) * cos(k_2*x + l_2*y + m_2*z),
    #initial_v = (x, y, z) -> -(k_1 * m_1 / (k_1^2+l_1^2)) * amp * exp(-(z - lz_1)^2/(4*sigma_sq)) * cos(k_1*x + l_1*y + m_1*z)
    #            -(k_2 * m_2 / (k_2^2+l_2^2)) * amp * exp(-(z - lz_2)^2/(4*sigma_sq)) * cos(k_2*x + l_2*y + m_2*z),
    #initial_w = (x, y, z) -> -(l_1 * m_1 / (k_1^2+l_1^2)) * amp * exp(-(z - lz_1)^2/(4*sigma_sq)) * cos(k_1*x + l_1*y + m_1*z)
    #            -(l_2 * m_2 / (k_2^2+l_2^2)) * amp * exp(-(z - lz_2)^2/(4*sigma_sq)) * cos(k_2*x + l_2*y + m_2*z),
    model = Boussinesq(),
    background = RadiatedBoussinesq(),
)
domain = DomainNamelist(;
    x_size = 40,
    y_size = 40,
    z_size = 60,
    lx,
    ly,
    lz,
    npx,
    npz,
)
grid = GridNamelist(;
    #resolved_topography = (x, y) -> 5 * (1 + cos(pi / 10000 * x)),
)
output =
    OutputNamelist(; output_variables = (:w, :dudt,), 
    output_file = "prop_perturb.h5", 
    output_interval = 0.5E+2,
    tmax = 0.5E+3,
)
sponge = SpongeNamelist(;
    #rhs_sponge = (x, y, z, t, dt) ->
    #    z >= zr ? sin(pi / 2 * (z - zr) / (lz - zr))^2 / dt : 0.0,
)

wkb = WKBNamelist(; wkb_mode = NoWKB())

integrate(Namelists(; atmosphere, domain, grid, output, sponge, wkb))
