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

  - `tkebg::A`: Background value of the turbulent kinetic energy.

  - `ttebg::A`: Background value of the total turbulent energy.

# Arguments

  - `constants`: Physical constants and reference values.
"""
struct TurbulenceAuxiliaries{
    A <: AbstractFloat,
    B <: AbstractArray{<:AbstractFloat, 3},
}
    tkebg::A
    ttebg::A
    tkebackup::B
end

function TurbulenceAuxiliaries(
    turbulencepredictands::TurbulencePredictands,
    constants::Constants,
)::TurbulenceAuxiliaries
    (; lref, tref) = constants
    (; tke) = turbulencepredictands

    tkebg = 0.1 * tref^2.0 / lref^2.0
    ttebg = 0.1 * tref^2.0 / lref^2.0
    tkebackup = deepcopy(tke)

    return TurbulenceAuxiliaries(tkebg, ttebg, tkebackup)
end
