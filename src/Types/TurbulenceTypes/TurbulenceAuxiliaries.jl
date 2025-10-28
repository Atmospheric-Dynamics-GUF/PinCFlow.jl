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
struct TurbulenceAuxiliaries{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 1},
}
    athomas::B
    bthomas::B
    cthomas::B
    fthomas::B
    qthomas::B
    shearproduction::A
    buoyancyproduction::A
end

function TurbulenceAuxiliaries(
    turbulencepredictands::TurbulencePredictands,
    constants::Constants,
    domain::Domain,
)::TurbulenceAuxiliaries
    (; lref, tref) = constants
    (; tke) = turbulencepredictands

    athomas = zeros(size(tke)[3])
    bthomas = zeros(size(tke)[3])
    cthomas = zeros(size(tke)[3])
    fthomas = zeros(size(tke)[3])
    qthomas = zeros(size(tke)[3])
    shearproduction = zeros(size(tke))
    buoyancyproduction = zeros(size(tke))

    return TurbulenceAuxiliaries(
        athomas,
        bthomas,
        cthomas,
        fthomas,
        qthomas,
        shearproduction,
        buoyancyproduction,
    )
end
