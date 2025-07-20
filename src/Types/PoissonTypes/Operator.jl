"""
```julia
Operator{A <: AbstractArray{<:AbstractFloat, 3}}
```

Workspace array for applying the linear operator of the Poisson solver.

# Fields

  - `s::A`: Auxiliary array for enforcing boundary conditions and performing MPI communication prior to the application of the linear operator.
"""
struct Operator{A <: AbstractArray{<:AbstractFloat, 3}}
    s::A
end

"""
```julia
Operator(domain::Domain)
```

Initialize operator-workspace array sized according to the dimensions of the MPI subdomain.

# Arguments

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

# Returns

  - `::Operator`: `Operator` instance with a zero-initialized array.
"""
function Operator(domain::Domain)

    # Get all necessary fields.
    (; nxx, nyy, nzz) = domain

    # Return an Operator instance.
    return Operator(zeros(nxx, nyy, nzz))
end
