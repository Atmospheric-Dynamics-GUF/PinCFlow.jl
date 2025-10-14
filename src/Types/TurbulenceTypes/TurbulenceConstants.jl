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
    cepsilon::A
    ck::A
    c3::A
    lturb::A
end

function TurbulenceConstants(namelists::Namelists)::TurbulenceConstants
    cepsilon = 8.71E-1
    ck = 0.084
    c3 = 1.185
    lturb = 1.0E+2

    return TurbulenceConstants(cepsilon, ck, c3, lturb)
end
