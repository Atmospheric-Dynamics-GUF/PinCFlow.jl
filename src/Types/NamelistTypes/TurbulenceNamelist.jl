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
TurbulenceNamelist(; turbulencesetup = NoTurbulence())
```
"""
function TurbulenceNamelist(; turbulencesetup = NoTurbulence())
    return TurbulenceNamelist(turbulencesetup)
end
