"""
```julia
AtmosphereNamelist{
    A <: Bool,
    B <: AbstractFloat,
    C <: AbstractBackground,
    D <: NTuple{3, <:AbstractFloat},
}
```

Namelist for parameters describing the atmospheric background.

```julia
AtmosphereNamelist(;
    specify_reynolds_number::Bool = false,
    inverse_reynolds_number::AbstractFloat = 0.0E+0,
    kinematic_viscosity::AbstractFloat = 0.0E+0,
    thermal_conductivity::AbstractFloat = 0.0E+0,
    kinematic_diffusivity::AbstractFloat = 0.0E+0,
    background::AbstractBackground = Isothermal(),
    buoyancy_frequency::AbstractFloat = 1.0E-2,
    potential_temperature::AbstractFloat = 3.0E+2,
    temperature::AbstractFloat = 3.0E+2,
    ground_pressure::AbstractFloat = 1.0E+5,
    initial_wind::NTuple{3, <:AbstractFloat} = (0.0E+0, 0.0E+0, 0.0E+0),
    coriolis_frequency::AbstractFloat = 0.0E+0,
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

  - `initial_wind::D`: Initial wind.

  - `coriolis_frequency::B`: Coriolis frequency of the ``f``-plane.
"""
struct AtmosphereNamelist{
    A <: Bool,
    B <: AbstractFloat,
    C <: AbstractBackground,
    D <: NTuple{3, <:AbstractFloat},
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
    initial_wind::D
    coriolis_frequency::B
    tropopause_height::B
    lapse_rate_troposphere::B
    lapse_rate_stratosphere::B
end

function AtmosphereNamelist(;
    specify_reynolds_number::Bool = false,
    inverse_reynolds_number::AbstractFloat = 0.0E+0,
    kinematic_viscosity::AbstractFloat = 0.0E+0,
    thermal_conductivity::AbstractFloat = 0.0E+0,
    kinematic_diffusivity::AbstractFloat = 0.0E+0,
    background::AbstractBackground = Isothermal(),
    buoyancy_frequency::AbstractFloat = 1.0E-2,
    potential_temperature::AbstractFloat = 3.0E+2,
    temperature::AbstractFloat = 3.0E+2,
    ground_pressure::AbstractFloat = 1.0E+5,
    initial_wind::NTuple{3, <:AbstractFloat} = (0.0E+0, 0.0E+0, 0.0E+0),
    coriolis_frequency::AbstractFloat = 0.0E+0,
    tropopause_height = 1.0E+4,
    lapse_rate_troposphere = 6.5E-3,
    lapse_rate_stratosphere = -5.0E-3,
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
        initial_wind,
        coriolis_frequency,
        tropopause_height,
        lapse_rate_troposphere,
        lapse_rate_stratosphere,
    )
end
