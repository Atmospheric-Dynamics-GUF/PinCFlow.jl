"""
```julia
compute_global_dot_product(
    a::AbstractArray{<:AbstractFloat, 3},
    b::AbstractArray{<:AbstractFloat, 3},
    domain::Domain,
)::AbstractFloat
```

Compute and return the dot product ``\\boldsymbol{a} \\cdot \\boldsymbol{b} = \\sum_i a_i \\cdot b_i`` of two 3D arrays distributed across MPI processes.

# Arguments

  - `a`: First input array.

  - `b`: Second input array (must have the same shape as `a`).

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.
"""
function compute_global_dot_product end

function compute_global_dot_product(
    a::AbstractArray{<:AbstractFloat, 3},
    b::AbstractArray{<:AbstractFloat, 3},
    domain::Domain,
)::AbstractFloat

    # Get parameters.
    (; comm) = domain

    # Compute local dot product.
    local_dot_product = mapreduce((a, b) -> a * b, +, a, b)

    # Sum over all processes.
    global_dot_product = MPI.Allreduce(local_dot_product, +, comm)

    # Return the result.
    return global_dot_product
end
