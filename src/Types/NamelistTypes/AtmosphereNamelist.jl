"""
```julia
AtmosphereNamelist{
    A <: AbstractModel,
    B <: Bool,
    C <: Float64,
    D <: AbstractBackground,
    E <: Function,
    F <: Function,
    G <: Function,
    H <: Function,
    I <: Function,
    J <: Function,
}
```

Namelist for parameters used in the definition of the atmospheric background and the initialization of prognostic variables.

```julia
AtmosphereNamelist(;
    model::AbstractModel = PseudoIncompressible(),
    specify_reynolds_number::Bool = false,
    inverse_reynolds_number::Real = 0.0E+0,
    kinematic_viscosity::Real = 1.5E-5,
    thermal_conductivity::Real = 3.0E-5,
    kinematic_diffusivity::Real = 0.0E+0,
    background::AbstractBackground = Isothermal(),
    buoyancy_frequency::Real = 1.0E-2,
    potential_temperature::Real = 3.0E+2,
    temperature::Real = 3.0E+2,
    ground_pressure::Real = 1.0E+5,
    coriolis_frequency::Real = 1.0E-4,
    tropopause_height::Real = 1.0E+4,
    troposphere_lapse_rate::Real = 6.5E-3,
    stratosphere_lapse_rate::Real = -5.0E-3,
    initial_rhop::Function = (x, y, z) -> 0.0,
    initial_u::Function = (x, y, z) -> 0.0,
    initial_v::Function = (x, y, z) -> 0.0,
    initial_w::Function = (x, y, z) -> 0.0,
    initial_pip::Function = (x, y, z) -> 0.0,
)::AtmosphereNamelist
```

Construct an `AtmosphereNamelist` instance with the given keyword arguments as properties, converting them to meet the type constraints.

# Fields/Keywords

  - `model::A`: Dynamic equations.

  - `specify_reynolds_number::B`: Flag to specify inverse Reynolds number instead of viscosity.

  - `inverse_reynolds_number::C`: Inverse Reynolds number.

  - `kinematic_viscosity::C`: Kinematic viscosity at the surface.

  - `thermal_conductivity::C`: Thermal conductivity at the surface.

  - `kinematic_diffusivity::C`: Kinematic diffusivity at the surface.

  - `background::D`: Atmospheric background.

  - `buoyancy_frequency::C`: Buoyancy frequency if `background == StableStratification()`.

  - `potential_temperature::C`: Reference potential temperature.

  - `temperature::C`: Reference temperature.

  - `ground_pressure::C`: Reference pressure.

  - `coriolis_frequency::C`: Coriolis frequency of the ``f``-plane.

  - `tropopause_height::C`: Height of the tropopause for `background == Realistic()` or `background == LapseRates()`.

  - `troposphere_lapse_rate::C`: Lapse rate in the troposphere for `background == LapseRates()`.

  - `stratosphere_lapse_rate::C`: Lapse rate in the stratosphere for `background == LapseRates()`.

  - `initial_rhop::E`: Function used to initialize the density fluctuations.

  - `initial_u::F`: Function used to initialize the zonal wind.

  - `initial_v::G`: Function used to initialize the meridional wind.

  - `initial_w::H`: Function used to initialize the vertical wind.

  - `initial_pip::I`: Function used to initialize the Exner-pressure fluctuations.
"""
struct AtmosphereNamelist{
    A <: AbstractModel,
    B <: Bool,
    C <: Float64,
    D <: AbstractBackground,
    E <: Function,
    F <: Function,
    G <: Function,
    H <: Function,
    I <: Function,
}
    model::A
    specify_reynolds_number::B
    inverse_reynolds_number::C
    kinematic_viscosity::C
    thermal_conductivity::C
    kinematic_diffusivity::C
    background::D
    buoyancy_frequency::C
    potential_temperature::C
    temperature::C
    ground_pressure::C
    coriolis_frequency::C
    tropopause_height::C
    troposphere_lapse_rate::C
    stratosphere_lapse_rate::C
    initial_rhop::E
    initial_u::F
    initial_v::G
    initial_w::H
    initial_pip::I
end

function AtmosphereNamelist(;
    model::AbstractModel = PseudoIncompressible(),
    specify_reynolds_number::Bool = false,
    inverse_reynolds_number::Real = 0.0E+0,
    kinematic_viscosity::Real = 1.5E-5,
    thermal_conductivity::Real = 3.0E-5,
    kinematic_diffusivity::Real = 0.0E+0,
    background::AbstractBackground = Isothermal(),
    buoyancy_frequency::Real = 1.0E-2,
    potential_temperature::Real = 3.0E+2,
    temperature::Real = 3.0E+2,
    ground_pressure::Real = 1.0E+5,
    coriolis_frequency::Real = 1.0E-4,
    tropopause_height::Real = 1.0E+4,
    troposphere_lapse_rate::Real = 6.5E-3,
    stratosphere_lapse_rate::Real = -5.0E-3,
    initial_rhop::Function = (x, y, z) -> 0.0,
    initial_u::Function = (x, y, z) -> 0.0,
    initial_v::Function = (x, y, z) -> 0.0,
    initial_w::Function = (x, y, z) -> 0.0,
    initial_pip::Function = (x, y, z) -> 0.0,
)::AtmosphereNamelist
    return AtmosphereNamelist(
        model,
        specify_reynolds_number,
        Float64(inverse_reynolds_number),
        Float64(kinematic_viscosity),
        Float64(thermal_conductivity),
        Float64(kinematic_diffusivity),
        background,
        Float64(buoyancy_frequency),
        Float64(potential_temperature),
        Float64(temperature),
        Float64(ground_pressure),
        Float64(coriolis_frequency),
        Float64(tropopause_height),
        Float64(troposphere_lapse_rate),
        Float64(stratosphere_lapse_rate),
        initial_rhop,
        initial_u,
        initial_v,
        initial_w,
        initial_pip,
    )
end
