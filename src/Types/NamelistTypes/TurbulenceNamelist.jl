"""
```julia
TurbulenceNamelist{A <: AbstractTurbulence}
```
"""
struct TurbulenceNamelist{A <: AbstractTurbulence}
    turbulencesetup::A
end

"""
```julia
TurbulenceNamelist(; turbulencesetup::AbstractTurbulence = NoTurbulence())
```
"""
function TurbulenceNamelist(;
    turbulencesetup::AbstractTurbulence = NoTurbulence(),
)
    return TurbulenceNamelist(turbulencesetup)
end
