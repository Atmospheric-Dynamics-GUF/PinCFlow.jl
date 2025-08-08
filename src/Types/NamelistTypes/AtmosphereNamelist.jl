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

Namelist for parameters describing the atmospheric background.

```julia
AtmosphereNamelist(;
    specifyreynolds::Bool = false,
    reinv::AbstractFloat = 0.0E+0,
    mu_viscous_dim::AbstractFloat = 0.0E+0,
    background::AbstractBackground = Isothermal(),
    buoyancy_frequency::AbstractFloat = 1.0E-2,
    theta0_dim::AbstractFloat = 3.0E+2,
    temp0_dim::AbstractFloat = 3.0E+2,
    press0_dim::AbstractFloat = 1.0E+5,
    backgroundflow_dim::NTuple{3, <:AbstractFloat} = (0.0E+0, 0.0E+0, 0.0E+0),
    coriolis_frequency::AbstractFloat = 0.0E+0,
    coriolis_mode::AbstractCoriolisMode = FPlane(),
)
```

Construct an `AtmosphereNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

- `specifyreynolds::A`: Flag to specify inverse Reynolds number instead of viscosity.
- `reinv::B`: Inverse Reynolds number.
- `mu_viscous_dim::B`: Kinematic viscosity at the surface.
- `background::C`: Atmospheric background.
- `buoyancy_frequency::B`: Buoyancy frequency if `background == StratifiedBoussinesq()`.
- `theta0_dim::B`: Reference potential temperature.
- `temp0_dim::B`: Reference temperature.
- `press0_dim::B`: Reference pressure.
- `backgroundflow_dim::D`: Initial wind.
- `coriolis_frequency::B`: Coriolis frequency if `coriolis_mode == FPlane()`
- `coriolis_mode::E`: Approximation used for the Coriolis frequency.
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

function AtmosphereNamelist(;
    specifyreynolds::Bool = false,
    reinv::AbstractFloat = 0.0E+0,
    mu_viscous_dim::AbstractFloat = 0.0E+0,
    background::AbstractBackground = Isothermal(),
    buoyancy_frequency::AbstractFloat = 1.0E-2,
    theta0_dim::AbstractFloat = 3.0E+2,
    temp0_dim::AbstractFloat = 3.0E+2,
    press0_dim::AbstractFloat = 1.0E+5,
    backgroundflow_dim::NTuple{3, <:AbstractFloat} = (0.0E+0, 0.0E+0, 0.0E+0),
    coriolis_frequency::AbstractFloat = 0.0E+0,
    coriolis_mode::AbstractCoriolisMode = FPlane(),
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
