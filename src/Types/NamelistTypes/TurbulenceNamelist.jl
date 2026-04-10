"""
```julia
TurbulenceNamelist
```

Namelist for the inclusion of a turbulence parameterization and the turbulent diffusion of momentum, mass-weighted potential temperature, and tracers.

```julia
TurbulenceNamelist(;
    turbulence_scheme::Symbol = :TKEScheme,
    momentum_coupling::Bool = true,
    entropy_coupling::Bool = true,
    tracer_coupling::Bool = true,
    initial_tke::Function = (x, y, z) -> 5e-5,
)::TurbulenceNamelist
```

Construct a `TurbulenceNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `turbulence_scheme::A`: General turbulence parameterization configuration.

  - `momentum_coupling::B`: Switch for turbulent diffusion of momentum.

  - `entropy_coupling::B`: Switch for turbulent diffusion of the the mass-specific potential temperature.

  - `tracer_coupling::B`: Switch for turbulent diffusion of tracers.

  - `initial_tke::C`: Function used to initialize the mass-specific turbulent kinetic energy.
"""
struct TurbulenceNamelist
    turbulence_scheme::Symbol
    momentum_coupling::Bool
    entropy_coupling::Bool
    tracer_coupling::Bool
    initial_tke::FunctionWrapper{Float64, NTuple{3, Float64}}
end

function TurbulenceNamelist(;
    turbulence_scheme::Symbol = :TKEScheme,
    momentum_coupling::Bool = true,
    entropy_coupling::Bool = true,
    tracer_coupling::Bool = true,
    initial_tke::Function = (x, y, z) -> 5e-5,
)::TurbulenceNamelist
    return TurbulenceNamelist(
        turbulence_scheme,
        momentum_coupling,
        entropy_coupling,
        tracer_coupling,
        initial_tke,
    )
end
