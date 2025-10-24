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

# Arguments

  - `constants`: Physical constants and reference values.
"""
struct TurbulenceAuxiliaries{
    A <: AbstractFloat,
    B <: AbstractArray{<:AbstractFloat, 3},
}
    tkebg::A
    km::B
    kh::B
    kek::B
    athomas::B
    bthomas::B
    cthomas::B
    fthomas::B
    qthomas::B 
    shearproduction::B 
    buoyancyproduction::B
end

function TurbulenceAuxiliaries(
    turbulencepredictands::TurbulencePredictands,
    constants::Constants,
    domain::Domain,
)::TurbulenceAuxiliaries
    (; lref, tref) = constants
    (; tke) = turbulencepredictands

    tkebg = 1.0E-8
    km = zeros(size(tke))
    kh = zeros(size(tke))
    kek = zeros(size(tke))
    athomas = zeros(size(tke))
    bthomas = zeros(size(tke))
    cthomas = zeros(size(tke))
    fthomas = zeros(size(tke))
    qthomas = zeros(size(tke))
    shearproduction = zeros(size(tke))
    buoyancyproduction = zeros(size(tke))

    return TurbulenceAuxiliaries(
        tkebg,
        km,
        kh,
        kek,
        athomas,
        bthomas,
        cthomas,
        fthomas,
        qthomas,
        shearproduction,
        buoyancyproduction,
    )
end
