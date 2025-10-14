"""
```julia
Auxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```

Auxiliary array used in the reconstruction of prognostic variables.

```julia
Auxiliaries(domain::Domain)::Auxiliaries
```

Construct an `Auxiliaries` instance with a zero-initialized auxiliary array sized according to the MPI subdomain dimensions.

# Fields

  - `phi::A`: Auxiliary array used as input for [`PinCFlow.FluxCalculator.apply_3d_muscl!`](@ref).

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct Auxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    phi::A
end

function Auxiliaries(domain::Domain)::Auxiliaries
    (; nxx, nyy, nzz) = domain

    return Auxiliaries(zeros(nxx, nyy, nzz))
end
