"""
```julia
AtmosphereNamelist{
    A <: Bool,
    B <: AbstractFloat,
    C <: AbstractBackground,
    D <: Function,
    E <: Function,
    F <: Function,
    G <: Function,
    H <: Function,
    I <: Function,
}
```

Namelist for parameters used in the definition of the atmospheric background and the initialization of prognostic variables.

```julia
AtmosphereNamelist(;
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

  - `specify_reynolds_number::A`: Flag to specify inverse Reynolds number instead of viscosity.

  - `inverse_reynolds_number::B`: Inverse Reynolds number.

  - `kinematic_viscosity::B`: Kinematic viscosity at the surface.

  - `thermal_conductivity::B`: Thermal conductivity at the surface.

  - `kinematic_diffusivity::B`: Kinematic diffusivity at the surface.

  - `background::C`: Atmospheric background.

  - `buoyancy_frequency::B`: Buoyancy frequency if `background == StratifiedBoussinesq()`.

  - `potential_temperature::B`: Reference potential temperature.

  - `temperature::B`: Reference temperature.

  - `ground_pressure::B`: Reference pressure.

  - `coriolis_frequency::B`: Coriolis frequency of the ``f``-plane.

  - `tropopause_height::B`: Height of the tropopause for `background == Realistic()` or `background == LapseRates()`.

  - `troposphere_lapse_rate::B`: Lapse rate in the troposphere for `background == LapseRates()`.

  - `stratosphere_lapse_rate::B`: Lapse rate in the stratosphere for `background == LapseRates()`.

  - `initial_rhop::D`: Function used to initialize the density fluctuations.

  - `initial_thetap::E`: Function used to initialize the potential-temperature fluctuations (only relevant in compressible mode).

  - `initial_u::F`: Function used to initialize the zonal wind.

  - `initial_v::G`: Function used to initialize the meridional wind.

  - `initial_w::H`: Function used to initialize the vertical wind.

  - `initial_pip::I`: Function used to initialize the Exner-pressure fluctuations.
"""
struct AtmosphereNamelist{
    A <: Bool,
    B <: AbstractFloat,
    C <: AbstractBackground,
    D <: Function,
    E <: Function,
    F <: Function,
    G <: Function,
    H <: Function,
    I <: Function,
}
    specify_reynolds_number::A
    inverse_reynolds_number::B
    kinematic_viscosity::B
    thermal_conductivity::B
    kinematic_diffusivity::B
    background::C
    buoyancy_frequency::B
    potential_temperature::B
    temperature::B
    ground_pressure::B
    coriolis_frequency::B
    tropopause_height::B
    troposphere_lapse_rate::B
    stratosphere_lapse_rate::B
    initial_rhop::D
    initial_thetap::E
    initial_u::F
    initial_v::G
    initial_w::H
    initial_pip::I
end

function AtmosphereNamelist(;
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
