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

  - `turbulencelengthscale::A`: The turbulence length scale L

  - `c3::A`: The inverse Prandtl number at neutrality

  - `cK::A`: Free parameter
  
Constants:

  - `alphaturb::A`: The dissipation constant ``\\alpha = \\frac{C_\\epsilon}{L}``.

  - `k_diff::A`: Pre-factor used to calculate the diffusion constant ``K = \\frac{C_{ek}}{L C_\\epsilon} <w'w'> ``, where
  the diffusion co-efficient ``K_ek = \\frac{C_{ek}}{C_\\epsilon} <w'w'> \\frac{\\sqrt{e_K}}{L}``

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
    turbulencelengthscale::A
    alphaturb::A
    k_diff::A
end

function TurbulenceConstants(namelists::Namelists)::TurbulenceConstants
    (; spongeheight) = namelists.sponge

    # Initialize the dissipation coefficients.
    cepsilon = 8.5E-1
    turbulencelengthscale = 1.0E-1
    alphaturb = cepsilon/turbulencelengthscale

    # Calcualte the pre-factor for the diffusion co-efficient
    c_ek = 1.0E-1
    rif = 1
    F = 6.4E-1
    lambdaturb = 4.0E-1
    lambda_2 = (5 * lambdaturb - 6 * sqrt(F) * lambdaturb)/25
    lambda_3 = (5 * lambdaturb + 14 * sqrt(F) * lambdaturb)/75
    lambda_4 = lambdaturb/4
    wpwp = (2 * (1 - 3 * lambda_3 + lambda_2 - 4 * lambda_4 * rif))/ (3 * (1- rif))
    k_diff = (c_ek * wpwp)/(cepsilon * turbulencelengthscale)

    # Return a TurbulenceConstants instance.
    return TurbulenceConstants(
        cepsilon,
        turbulencelengthscale,
        alphaturb,
        k_diff
    )
end
