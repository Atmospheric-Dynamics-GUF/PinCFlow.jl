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
struct TurbulenceNamelist{A <: AbstractTurbulence, B <: Bool}
    turbulence_scheme::A
    momentum_coupling::B
    entropy_coupling::B
end

function TurbulenceNamelist(;
    turbulence_scheme::AbstractTurbulence = NoTurbulence(),
    momentum_coupling::Bool = true,
    entropy_coupling::Bool = true,
)::TurbulenceNamelist
    return TurbulenceNamelist(
        turbulence_scheme,
        momentum_coupling,
        entropy_coupling,
    )
end
