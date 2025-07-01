"""
```julia
compute_flux(usurf::AbstractFloat, phiup::AbstractFloat, phidown::AbstractFloat)
```

Compute upwind flux based on surface velocity direction.

Selects the appropriate scalar value (`phiup` or `phidown`) based on the sign of the surface velocity `usurf` to ensure numerical stability through upwind differencing.

# Arguments

  - `usurf::AbstractFloat`: Surface velocity (positive = rightward/upward flow)
  - `phiup::AbstractFloat`: Scalar value upstream when flow is negative
  - `phidown::AbstractFloat`: Scalar value upstream when flow is positive

# Returns

  - `AbstractFloat`: Flux value `usurf * phi` where `phi` is chosen by flow direction
"""
function compute_flux(
    usurf::AbstractFloat,
    phiup::AbstractFloat,
    phidown::AbstractFloat,
)
    if usurf > 0.0
        return usurf * phiup
    else
        return usurf * phidown
    end
end
