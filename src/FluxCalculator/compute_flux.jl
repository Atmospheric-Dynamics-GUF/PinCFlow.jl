"""
```julia
compute_flux(usurf::AbstractFloat, phiup::AbstractFloat, phidown::AbstractFloat)::AbstractFloat
```

Compute and return the upstream flux from reconstructed values, based on the sign of the transporting velocity.

# Arguments

  - `usurf`: Transporting velocity.

  - `phiup`: Upstream reconstruction for `usurf > 0`.

  - `phidown`: Downstream reconstruction for `usurf > 0`.
"""
function compute_flux end

function compute_flux(
    usurf::AbstractFloat,
    phiup::AbstractFloat,
    phidown::AbstractFloat,
)::AbstractFloat
    if usurf > 0.0
        return usurf * phiup
    else
        return usurf * phidown
    end
end
