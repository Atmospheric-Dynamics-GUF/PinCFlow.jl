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
    wave_impact::Bool = true,
)::TurbulenceNamelist
```

Construct a `TurbulenceNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `turbulence_scheme::Symbol`: General turbulence parameterization configuration.

  - `momentum_coupling::Bool`: Switch for turbulent diffusion of momentum.

  - `entropy_coupling::Bool`: Switch for turbulent diffusion of the mass-specific potential temperature.

  - `tracer_coupling::Bool`: Switch for turbulent diffusion of tracers.

  - `initial_tke::FunctionWrapper{Float64, NTuple{3, Float64}}`: Function used to initialize the mass-specific turbulent kinetic energy.

  - `wave_impact::Bool`: Switch for turbulence production due to unresolved gravity wave shear.
"""
struct TurbulenceNamelist
    turbulence_scheme::Symbol
    momentum_coupling::Bool
    entropy_coupling::Bool
    tracer_coupling::Bool
    initial_tke::FunctionWrapper{Float64, NTuple{3, Float64}}
    wave_impact::Bool
    turbulence_filter_order::Integer
    turbulence_filter_type::Symbol
    smooth_turbulence::Bool
end

function TurbulenceNamelist(;
    turbulence_scheme::Symbol = :TKEScheme,
    momentum_coupling::Bool = true,
    entropy_coupling::Bool = true,
    tracer_coupling::Bool = true,
    initial_tke::Function = (x, y, z) -> 5e-5,
    wave_impact::Bool = true,
    turbulence_filter_order::Integer = 3,
    turbulence_filter_type::Symbol = :BoxFilter,
    smooth_turbulence::Bool = true,
)::TurbulenceNamelist
    return TurbulenceNamelist(
        turbulence_scheme,
        momentum_coupling,
        entropy_coupling,
        tracer_coupling,
        initial_tke,
        wave_impact,
        turbulence_filter_order,
        turbulence_filter_type,
        smooth_turbulence,
    )
end
