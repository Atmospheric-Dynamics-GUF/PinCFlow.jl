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
struct TurbulenceNamelist{A <: AbstractTurbulence}
    turbulence_scheme::A
end

function TurbulenceNamelist(;
    turbulence_scheme::AbstractTurbulence = NoTurbulence(),
)::TurbulenceNamelist
    return TurbulenceNamelist(turbulence_scheme)
end
