"""
```julia
TurbulenceConstants{
    A <: AbstractFloat,
}
```

Composite type for the turbulence constants like the dissipation and diffusion co-efficients, closure and calibration constants, and the stability functions.

```julia
TurbulenceConstants(namelists::Namelists)::TurbulenceConstants
```

Construct a `TurbulenceConstants` instance, using the model parameters in `namelists`.

This constructor creates an instance with all the turbulence constants.

# Fields

Free parameters:

  - `cepsilon::A`: Closure constant controlling the intensity of turbulent dissipation.
  

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `grid`: Collection of parameters and fields that describe the grid.
"""
struct TurbulenceConstants{A <: AbstractFloat}

    # Dissipation constant
    ck::A
    lturb::A
    tkemin::A
end

function TurbulenceConstants(
    namelists::Namelists,
    constants::Constants,
)::TurbulenceConstants
    (; lref, tref, rhoref) = constants

    ck = 0.5 # according to Therry and Lacarrere (1983)
    lturb = 1.0E+2
    tkemin = 5.E-5 * tref^2.0 / lref^2.0

    return TurbulenceConstants(ck, lturb, tkemin)
end
