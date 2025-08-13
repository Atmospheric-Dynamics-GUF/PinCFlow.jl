"""
```julia
Correction{A <: AbstractArray{<:AbstractFloat, 3}}
```

Correction terms used to update the horizontal wind in the corrector step.

```julia
Correction(domain::Domain)
```

Create a `Correction` instance with zero-initialized arrays sized according to the dimensions of the MPI subdomain.

# Fields

  - `corx::A`: Correction term for the zonal wind.

  - `cory::A`: Correction term for the meridional wind.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct Correction{A <: AbstractArray{<:AbstractFloat, 3}}
    corx::A
    cory::A
end

function Correction(domain::Domain)

    # Get all necessary fields.
    (; nxx, nyy, nzz) = domain

    # Return a Correction instance.
    return Correction([zeros(nxx, nyy, nzz) for i in 1:2]...)
end
