"""
```julia
Auxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```

Auxiliary array used in the reconstruction of prognostic variables.

# Fields

  - `phi::A`: Auxiliary array used as input for [`PinCFlow.FluxCalculator.apply_3d_muscl!`](@ref).
"""
struct Auxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    phi::A
end

"""
```julia
Auxiliaries(domain::Domain)
```

Construct an `Auxiliaries` instance with a zero-initialized auxiliary array sized according to the MPI subdomain dimensions.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

# Returns

  - `::Auxiliaries`: `Auxiliaries` instance with a zero-initialized array.
"""
function Auxiliaries(domain::Domain)
    (; nxx, nyy, nzz) = domain

    # Initialize the auxiliaries.
    phi = zeros(nxx, nyy, nzz)

    # Return an Auxiliaries instance.
    return Auxiliaries(phi)
end
