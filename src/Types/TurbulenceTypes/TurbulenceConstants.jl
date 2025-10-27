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
    lturb::A
    lturb_ndim::A
    tkemin::A
    prandtl::A
    prandtlinv::A
end

function TurbulenceConstants(
    namelists::Namelists,
    constants::Constants,
)::TurbulenceConstants
    (; lref, tref) = constants

    lturb = 30.0
    lturb_ndim = lturb / lref
    tkemin = 5.E-5 * tref^2.0 / lref^2.0
    prandtl = 0.85
    prandtlinv = 1.0 / prandtl

    return TurbulenceConstants(lturb, lturb_ndim, tkemin, prandtl, prandtlinv)
end
