"""
```julia
AtmosphereNamelist{
    A <: Bool,
    B <: AbstractFloat,
    C <: AbstractBackground,
    D <: NTuple{3, <:AbstractFloat},
    E <: AbstractCoriolisMode,
}
```

Namelist for the atmosphere (see constructor for parameter descriptions).
"""
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

"""
```julia
AtmosphereNamelist(;
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
```

Constructor for the `AtmosphereNamelist` type, which defines parameters for the background atmosphere.

# Arguments

  - `specifyreynolds`: Flag to specify Reynolds number instead of viscosity
  - `reinv`: Inverse Reynolds number
  - `mu_viscous_dim`: Dimensional viscosity
  - `background`: Background atmosphere model. Can be `Isothermal()`, `StratifiedBoussinesq()`, or `Stratified()`
  - `buoyancy_frequency`: Buoyancy frequency. Used only if `background` is `StratifiedBoussinesq()`
  - `theta0_dim`: Reference potential temperature
  - `temp0_dim`: Reference temperature
  - `press0_dim`: Reference pressure
  - `backgroundflow_dim`: Background flow in 3D
  - `coriolis_frequency`: Coriolis frequency. Used only if `background` is `Stratified()`
  - `coriolis_mode`: Coriolis mode. Used only if `background` is `Stratified()`

# Returns

  - `AtmosphereNamelist`: Initialized namelist with specified parameters
"""
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
