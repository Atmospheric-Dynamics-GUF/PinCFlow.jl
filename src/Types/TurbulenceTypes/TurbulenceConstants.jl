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

  - `c3::A`: The inverse Prandtl number at neutrality

  - `cK::A`: Free parameter 

Exchange co-efficients:

  - `KM::A`: Exchange co-efficient for momentum ``K_M = C_K L \\chi_3 \\sqrt(e_K)``.

  - `KH::A`: Exchange co-efficient for heat ``K_H = C_3 C_K L \\phi_3 \\sqrt(e_K)``.

Stability functions:

  - `chi3::A`: Stability function used in the computation of the exchange co-efficient for momentum.

  - `phi3::A`: Stability function used in the computation of the exchange co-efficient for heat.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `grid`: Collection of parameters and fields that describe the grid.
"""
struct TurbulenceConstants{
    A <: AbstractFloat,
}

    # Dissipation constant
    cepsilon::A
end

function TurbulenceConstants(namelists::Namelists)::TurbulenceConstants
    (; spongeheight) = namelists.sponge

    # Initialize the dissipation coefficients.
    cepsilon = 0.85

    # Return a TurbulenceConstants instance.
    return TurbulenceConstants(
        cepsilon,
    )
end
