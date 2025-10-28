"""
```julia
AtmosphereNamelist{
    A <: AbstractModel,
    B <: Bool,
    C <: AbstractFloat,
    D <: AbstractFloat,
    E <: AbstractFloat,
    F <: AbstractFloat,
    G <: AbstractBackground,
    H <: AbstractFloat,
    I <: AbstractFloat,
    J <: AbstractFloat,
    K <: AbstractFloat,
    L <: AbstractFloat,
    M <: AbstractFloat,
    N <: AbstractFloat,
    O <: AbstractFloat,
    P <: Function,
    Q <: Function,
    R <: Function,
    S <: Function,
    T <: Function,
    U <: Function,
}
```

Namelist for parameters used in the definition of the atmospheric background and the initialization of prognostic variables.

```julia
AtmosphereNamelist(;
    model::AbstractModel = PseudoIncompressible(),
    specify_reynolds_number::Bool = false,
    inverse_reynolds_number::AbstractFloat = 0.0E+0,
    kinematic_viscosity::AbstractFloat = 1.5E-5,
    thermal_conductivity::AbstractFloat = 3.0E-5,
    kinematic_diffusivity::AbstractFloat = 0.0E+0,
    background::AbstractBackground = Isothermal(),
    buoyancy_frequency::AbstractFloat = 1.0E-2,
    potential_temperature::AbstractFloat = 3.0E+2,
    temperature::AbstractFloat = 3.0E+2,
    ground_pressure::AbstractFloat = 1.0E+5,
    coriolis_frequency::AbstractFloat = 1.0E-4,
    tropopause_height::AbstractFloat = 1.0E+4,
    troposphere_lapse_rate::AbstractFloat = 6.5E-3,
    stratosphere_lapse_rate::AbstractFloat = -5.0E-3,
    initial_rhop::Function = (x, y, z) -> 0.0,
    initial_thetap::Function = (x, y, z) -> 0.0,
    initial_u::Function = (x, y, z) -> 0.0,
    initial_v::Function = (x, y, z) -> 0.0,
    initial_w::Function = (x, y, z) -> 0.0,
    initial_pip::Function = (x, y, z) -> 0.0,
)::AtmosphereNamelist
```

Construct an `AtmosphereNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `model::A`: Dynamic equations.

  - `specify_reynolds_number::B`: Flag to specify inverse Reynolds number instead of viscosity.

  - `inverse_reynolds_number::C`: Inverse Reynolds number.

  - `kinematic_viscosity::D`: Kinematic viscosity at the surface.

  - `thermal_conductivity::E`: Thermal conductivity at the surface.

  - `kinematic_diffusivity::F`: Kinematic diffusivity at the surface.

  - `background::G`: Atmospheric background.

  - `buoyancy_frequency::H`: Buoyancy frequency if `background == StratifiedBoussinesq()`.

  - `potential_temperature::I`: Reference potential temperature.

  - `temperature::J`: Reference temperature.

  - `ground_pressure::K`: Reference pressure.

  - `coriolis_frequency::L`: Coriolis frequency of the ``f``-plane.

  - `tropopause_height::M`: Height of the tropopause for `background == Realistic()` or `background == LapseRates()`.

  - `troposphere_lapse_rate::N`: Lapse rate in the troposphere for `background == LapseRates()`.

  - `stratosphere_lapse_rate::O`: Lapse rate in the stratosphere for `background == LapseRates()`.

  - `initial_rhop::P`: Function used to initialize the density fluctuations.

  - `initial_thetap::Q`: Function used to initialize the potential-temperature fluctuations (only relevant in compressible mode).

  - `initial_u::R`: Function used to initialize the zonal wind.

  - `initial_v::S`: Function used to initialize the meridional wind.

  - `initial_w::T`: Function used to initialize the vertical wind.

  - `initial_pip::U`: Function used to initialize the Exner-pressure fluctuations.
"""
struct AtmosphereNamelist{
    A <: AbstractModel,
    B <: Bool,
    C <: AbstractFloat,
    D <: AbstractFloat,
    E <: AbstractFloat,
    F <: AbstractFloat,
    G <: AbstractBackground,
    H <: AbstractFloat,
    I <: AbstractFloat,
    J <: AbstractFloat,
    K <: AbstractFloat,
    L <: AbstractFloat,
    M <: AbstractFloat,
    N <: AbstractFloat,
    O <: AbstractFloat,
    P <: Function,
    Q <: Function,
    R <: Function,
    S <: Function,
    T <: Function,
    U <: Function,
}
    model::A
    specify_reynolds_number::B
    inverse_reynolds_number::C
    kinematic_viscosity::D
    thermal_conductivity::E
    kinematic_diffusivity::F
    background::G
    buoyancy_frequency::H
    potential_temperature::I
    temperature::J
    ground_pressure::K
    coriolis_frequency::L
    tropopause_height::M
    troposphere_lapse_rate::N
    stratosphere_lapse_rate::O
    initial_rhop::P
    initial_thetap::Q
    initial_u::R
    initial_v::S
    initial_w::T
    initial_pip::U
end

function AtmosphereNamelist(;
    model::AbstractModel = PseudoIncompressible(),
    specify_reynolds_number::Bool = false,
    inverse_reynolds_number::AbstractFloat = 0.0E+0,
    kinematic_viscosity::AbstractFloat = 1.5E-5,
    thermal_conductivity::AbstractFloat = 3.0E-5,
    kinematic_diffusivity::AbstractFloat = 0.0E+0,
    background::AbstractBackground = Isothermal(),
    buoyancy_frequency::AbstractFloat = 1.0E-2,
    potential_temperature::AbstractFloat = 3.0E+2,
    temperature::AbstractFloat = 3.0E+2,
    ground_pressure::AbstractFloat = 1.0E+5,
    coriolis_frequency::AbstractFloat = 1.0E-4,
    tropopause_height::AbstractFloat = 1.0E+4,
    troposphere_lapse_rate::AbstractFloat = 6.5E-3,
    stratosphere_lapse_rate::AbstractFloat = -5.0E-3,
    initial_rhop::Function = (x, y, z) -> 0.0,
    initial_thetap::Function = (x, y, z) -> 0.0,
    initial_u::Function = (x, y, z) -> 0.0,
    initial_v::Function = (x, y, z) -> 0.0,
    initial_w::Function = (x, y, z) -> 0.0,
    initial_pip::Function = (x, y, z) -> 0.0,
)::AtmosphereNamelist
    return AtmosphereNamelist(
        model,
        specify_reynolds_number,
        inverse_reynolds_number,
        kinematic_viscosity,
        thermal_conductivity,
        kinematic_diffusivity,
        background,
        buoyancy_frequency,
        potential_temperature,
        temperature,
        ground_pressure,
        coriolis_frequency,
        tropopause_height,
        troposphere_lapse_rate,
        stratosphere_lapse_rate,
        initial_rhop,
        initial_thetap,
        initial_u,
        initial_v,
        initial_w,
        initial_pip,
    )
end
