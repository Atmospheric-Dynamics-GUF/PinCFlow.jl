"""
```julia
TurbulenceNamelist{A <: AbstractTurbulence}
```

Namelist for the inclusion of turbulence physics.

```julia
TurbulenceNamelist(;
    turbulencesetup::AbstractTurbulence = NoTurbulence(),
)::TurbulenceNamelist
```

Construct a `TurbulenceNamelist` instance with the given keyword arguments as properties.

# Fields/Keywords

  - `turbulencesetup::A`: General turbulence-physics configuration.
"""
struct TurbulenceNamelist{A <: AbstractTurbulence}
    turbulencesetup::A
end

function TurbulenceNamelist(;
    turbulencesetup::AbstractTurbulence = NoTurbulence(),
)::TurbulenceNamelist
    return TurbulenceNamelist(turbulencesetup)
end
