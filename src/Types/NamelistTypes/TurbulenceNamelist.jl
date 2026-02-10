"""
```julia
TurbulenceNamelist{A <: AbstractTurbulence}
```

Namelist for the inclusion of a turbulence parameterization.

```julia
TurbulenceNamelist(;
    turbulence_scheme::AbstractTurbulence = NoTurbulence(),
)::TurbulenceNamelist
```

Construct a `TurbulenceNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `turbulence_scheme::A`: Turbulence parameterization scheme.
"""
struct TurbulenceNamelist{A <: AbstractTurbulence, B <: Bool, C <: Function}
    turbulence_scheme::A
    momentum_coupling::B
    entropy_coupling::B
    tracer_coupling::B
    initial_tke::C
end

function TurbulenceNamelist(;
    turbulence_scheme::AbstractTurbulence = NoTurbulence(),
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
