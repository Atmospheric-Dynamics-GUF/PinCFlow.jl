"""
```julia
TurbulenceConstants{A <: AbstractFloat}
```

Composite type for the dimensional and non-dimensional turbulence length scales and the minimum turbulent kinetic energy.

```julia
TurbulenceConstants(
    namelists::Namelists,
    constants::Constants,
)::TurbulenceConstants
```

Construct a `TurbulenceConstants` instance.

# Fields

  - `lturb::A`: Characteristic turbulence length scale ``L = 30 \\ \\mathrm{m}``.

  - `ld::A`: Non-dimensional turbulent mixing length for dissipation ``l_\\mathrm{d} = \\sqrt{2} L/L_\\mathrm{ref}``.

  - `lv::A`: Non-dimensional turbulent mixing length for momentum diffusion ``l_\\mathrm{d} = L/L_\\mathrm{ref}/\\sqrt{2}``.

  - `lb::A`: Non-dimensional turbulent mixing length for entropy diffusion ``l_\\mathrm{b} =  L/L_\\mathrm{ref}/\\sqrt{2}``.

  - `lt::A`: Non-dimensional turbulent mixing length for turbulence diffusion ``l_\\mathrm{t} = L/L_\\mathrm{ref}/\\sqrt{2}``.

  - `tkemin::A`: Minimum turbulent kinetic energy value ``e_\\mathrm{k} = 5\\times 10^{-5} \\ \\mathrm{m^2 \\ s^{-2}}``.
  
# Arguments

  - `constants`: Physical constsants and reference values.
"""
struct TurbulenceConstants{A <: AbstractFloat}

    # Dissipation constant
    lturb::A
    ld::A
    lv::A
    lb::A
    lt::A
    tkemin::A
end

function TurbulenceConstants(
    constants::Constants,
)::TurbulenceConstants
    (; lref, tref) = constants

    lturb = 30.0
    ld = sqrt(2) * lturb / lref
    lv = lb = lt = lturb / lref / sqrt(2)
    tkemin = 5.E-5 * tref^2.0 / lref^2.0

    return TurbulenceConstants(lturb, ld, lv, lb, lt, tkemin)
end
