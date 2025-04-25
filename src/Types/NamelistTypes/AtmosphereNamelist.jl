struct AtmosphereNamelist{
    A <: Bool,
    B <: AbstractFloat,
    C <: AbstractBackground,
    D <: NTuple{3, <:AbstractFloat},
    E <: AbstractCoriolisMode,
}
    specifyreynolds::A
    reinv::B
    mu_viscous_dim::B
    background::C
    buoyancy_frequency::B
    theta0_dim::B
    temp0_dim::B
    press0_dim::B
    backgroundflow_dim::D
    coriolis_frequency::B
    coriolis_mode::E
end

function AtmosphereNamelist(;
    specifyreynolds = false,
    reinv = 0.0E+0,
    mu_viscous_dim = 0.0E+0,
    background = Isothermal(),
    buoyancy_frequency = 1.0E-2,
    theta0_dim = 3.0E+2,
    temp0_dim = 3.0E+2,
    press0_dim = 1.0E+5,
    backgroundflow_dim = (0.0E+0, 0.0E+0, 0.0E+0),
    coriolis_frequency = 0.0E+0,
    coriolis_mode = FPlane(),
)
    return AtmosphereNamelist(
        specifyreynolds,
        reinv,
        mu_viscous_dim,
        background,
        buoyancy_frequency,
        theta0_dim,
        temp0_dim,
        press0_dim,
        backgroundflow_dim,
        coriolis_frequency,
        coriolis_mode,
    )
end
