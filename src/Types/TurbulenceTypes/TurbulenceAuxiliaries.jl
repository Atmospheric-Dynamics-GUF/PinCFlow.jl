"""
```julia
TurbulenceAuxiliaries{A <: AbstractFloat}
```

Background values for turbulence variables.

```julia
TurbulenceAuxiliaries(constants::Constants)::TurbulenceAuxiliaries
```

Construct a `TurbulenceAuxiliaries` instance with both fields set to ``t_\\mathrm{ref}^2 / \\left(10 L_\\mathrm{ref}^2\\right)``, where ``t_\\mathrm{ref}`` and ``L_\\mathrm{ref}`` are given by the properties `tref` and `lref` of `constants`, respectively.

# Fields

# Arguments

  - `constants`: Physical constants and reference values.
"""
struct TurbulenceAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    shearproduction::A
    buoyancyproduction::A
end

function TurbulenceAuxiliaries(
    namelists::Namelists,
    domain::Domain,
)::TurbulenceAuxiliaries
    (; turbulence_scheme) = namelists.turbulence

    return TurbulenceAuxiliaries(domain, turbulence_scheme)
end

function TurbulenceAuxiliaries(
    domain::Domain,
    turbulence_scheme::NoTurbulence,
)::TurbulenceAuxiliaries
    return TurbulenceAuxiliaries([zeros(0, 0, 0) for i in 1:2]...)
end

function TurbulenceAuxiliaries(
    domain::Domain,
    turbulence_scheme::TKEScheme,
)::TurbulenceAuxiliaries
    (; nxx, nyy, nzz, nx, ny, nz) = domain
    return TurbulenceAuxiliaries([zeros(nxx, nyy, nzz) for i in 1:2]...)
end
