"""
```julia
TurbulenceNamelist{A <: AbstractTurbulence}
```

```julia
TurbulenceNamelist(; turbulencesetup::AbstractTurbulence = NoTurbulence())
```
"""
struct TurbulenceNamelist{A <: AbstractTurbulence}
    turbulencesetup::A
end

function TurbulenceNamelist(;
    turbulencesetup::AbstractTurbulence = NoTurbulence(),
)
    return TurbulenceNamelist(turbulencesetup)
end
