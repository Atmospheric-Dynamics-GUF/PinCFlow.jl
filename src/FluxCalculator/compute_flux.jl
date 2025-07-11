"""
```julia
compute_flux(usurf::AbstractFloat, phiup::AbstractFloat, phidown::AbstractFloat)
```

Compute the upstream flux from reconstructed values, based on the sign of the transporting velocity.

# Arguments

  - `usurf`: Transporting velocity.
  - `phiup`: Upstream reconstruction for `usurf > 0`.
  - `phidown`: Downstream reconstruction for `usurf > 0`.

# Returns

  - `::Float64`: Product of `usurf` and the upstream reconstruction.
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
