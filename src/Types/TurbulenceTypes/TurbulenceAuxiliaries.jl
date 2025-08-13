"""
```julia
TurbulenceAuxiliaries{A <: AbstractFloat}
```

Background values for turbulence variables.

```julia
TurbulenceAuxiliaries(constants::Constants)
```

Construct a `TurbulenceAuxiliaries` instance with both fields set to ``t_\\mathrm{ref}^2 / \\left(10 L_\\mathrm{ref}^2\\right)``, where ``t_\\mathrm{ref}`` and ``L_\\mathrm{ref}`` are given by the properties `tref` and `lref` of `constants`, respectively.

# Fields

  - `tkebg::A`: Background value of the turbulent kinetic energy.

  - `ttebg::A`: Background value of the total turbulent energy.

# Arguments

  - `constants`: Physical constants and reference values.
"""
struct TurbulenceAuxiliaries{A <: AbstractFloat}
    tkebg::A
    ttebg::A
end

function TurbulenceAuxiliaries(constants::Constants)
    (; lref, tref) = constants

    tkebg = 0.1 * tref ^ 2.0 / lref ^ 2.0
    ttebg = 0.1 * tref ^ 2.0 / lref ^ 2.0

    return TurbulenceAuxiliaries(tkebg, ttebg)
end
