"""
```julia
AtmosphereNamelist
```

Namelist for parameters used in the definition of the atmospheric background and the initialization of prognostic variables.

```julia
AtmosphereNamelist(;
    model::Symbol = :PseudoIncompressible,
    specify_reynolds_number::Bool = false,
    inverse_reynolds_number::Real = 0.0E+0,
    kinematic_viscosity::Real = 1.5E-5,
    thermal_conductivity::Real = 3.0E-5,
    kinematic_diffusivity::Real = 0.0E+0,
    background::Symbol = :Isothermal,
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
    initial_thetap::Function = (x, y, z) -> 0.0,
    buoyancy_initialization::Symbol = :initial_rhop,
)::AtmosphereNamelist
```

Construct an `AtmosphereNamelist` instance with the given keyword arguments as properties, converting them to meet the type constraints.

# Fields/Keywords

  - `model::Symbol`: Dynamic equations.

  - `specify_reynolds_number::Bool`: Flag to specify inverse Reynolds number instead of viscosity.

  - `inverse_reynolds_number::Float64`: Inverse Reynolds number.

  - `kinematic_viscosity::Float64`: Kinematic viscosity at the surface.

  - `thermal_conductivity::Float64`: Thermal conductivity at the surface.

  - `kinematic_diffusivity::Float64`: Kinematic diffusivity at the surface.

  - `background::Symbol`: Atmospheric background.

  - `buoyancy_frequency::Float64`: Buoyancy frequency if `background == StableStratification()`.

  - `potential_temperature::Float64`: Reference potential temperature.

  - `temperature::Float64`: Reference temperature.

  - `ground_pressure::Float64`: Reference pressure.

  - `coriolis_frequency::Float64`: Coriolis frequency of the ``f``-plane.

  - `tropopause_height::Float64`: Height of the tropopause for `background == Realistic()` or `background == LapseRates()`.

  - `troposphere_lapse_rate::Float64`: Lapse rate in the troposphere for `background == LapseRates()`.

  - `stratosphere_lapse_rate::Float64`: Lapse rate in the stratosphere for `background == LapseRates()`.

  - `initial_rhop::FunctionWrapper{Float64, NTuple{3, Float64}}`: Function used to initialize the density fluctuations if `buoyancy_initialization == :initial_rhop`.

  - `initial_u::FunctionWrapper{Float64, NTuple{3, Float64}}`: Function used to initialize the zonal wind.

  - `initial_v::FunctionWrapper{Float64, NTuple{3, Float64}}`: Function used to initialize the meridional wind.

  - `initial_w::FunctionWrapper{Float64, NTuple{3, Float64}}`: Function used to initialize the vertical wind.

  - `initial_pip::FunctionWrapper{Float64, NTuple{3, Float64}}`: Function used to initialize the Exner-pressure fluctuations.

  - `initial_thetap::FunctionWrapper{Float64, NTuple{3, Float64}}`: Function used to initialize the density fluctuations if `buoyancy_initialization == :initial_thetap`.

  - `buoyancy_initialization::Symbol`: Switch for using `initial_rhop` (if set to `:initial_rhop`) or `initial_thetap` (if set to `:initial_thetap`) to initialize the buoyancy.
"""
struct AtmosphereNamelist
    model::Symbol
    specify_reynolds_number::Bool
    inverse_reynolds_number::Float64
    kinematic_viscosity::Float64
    thermal_conductivity::Float64
    kinematic_diffusivity::Float64
    background::Symbol
    buoyancy_frequency::Float64
    potential_temperature::Float64
    temperature::Float64
    ground_pressure::Float64
    coriolis_frequency::Float64
    tropopause_height::Float64
    troposphere_lapse_rate::Float64
    stratosphere_lapse_rate::Float64
    initial_rhop::FunctionWrapper{Float64, NTuple{3, Float64}}
    initial_u::FunctionWrapper{Float64, NTuple{3, Float64}}
    initial_v::FunctionWrapper{Float64, NTuple{3, Float64}}
    initial_w::FunctionWrapper{Float64, NTuple{3, Float64}}
    initial_pip::FunctionWrapper{Float64, NTuple{3, Float64}}
    initial_thetap::FunctionWrapper{Float64, NTuple{3, Float64}}
    buoyancy_initialization::Symbol
end

function AtmosphereNamelist(;
    model::Symbol = :PseudoIncompressible,
    specify_reynolds_number::Bool = false,
    inverse_reynolds_number::Real = 0.0E+0,
    kinematic_viscosity::Real = 1.5E-5,
    thermal_conductivity::Real = 3.0E-5,
    kinematic_diffusivity::Real = 0.0E+0,
    background::Symbol = :Isothermal,
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
    initial_thetap::Function = (x, y, z) -> 0.0,
    buoyancy_initialization::Symbol = :initial_rhop,
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
        initial_thetap,
        buoyancy_initialization,
    )
end
