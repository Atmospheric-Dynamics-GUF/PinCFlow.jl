"""
```julia
Operator{A <: AbstractArray{<:AbstractFloat, 3}}
```

Workspace array for applying the linear operator of the Poisson solver.

```julia
Operator(namelists::Namelists, domain::Domain)::Operator
```

Create an `Operator` instance with a zero-initialized array sized according to the dimensions of the MPI subdomain.

# Fields

  - `s::A`: Auxiliary array for enforcing boundary conditions and performing MPI communication prior to the application of the linear operator.

# Arguments

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
struct Operator{A <: AbstractArray{<:AbstractFloat, 3}}
    s::A
end

function Operator(namelists::Namelists, domain::Domain)::Operator
    (; float_type) = namelists.discretization
    (; nxx, nyy, nzz) = domain

    return Operator(zeros(float_type, nxx, nyy, nzz))
end
