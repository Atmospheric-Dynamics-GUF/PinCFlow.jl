"""
```julia
Auxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
```

Auxiliary array used in the reconstruction of prognostic variables.

```julia
Auxiliaries(namelists::Namelists, domain::Domain)::Auxiliaries
```

Construct an `Auxiliaries` instance with a zero-initialized auxiliary array sized according to the MPI subdomain dimensions.

# Fields

  - `phi::A`: Auxiliary array used as input for [`PinCFlow.FluxCalculator.apply_3d_muscl!`](@ref).

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct Auxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    phi::A
end

function Auxiliaries(namelists::Namelists, domain::Domain)::Auxiliaries
    (; float_type) = namelists.discretization
    (; nxx, nyy, nzz) = domain

    return Auxiliaries(zeros(float_type, nxx, nyy, nzz))
end
