"""
```julia
compute_global_dot_product(
    a::AbstractArray{<:AbstractFloat, 3},
    b::AbstractArray{<:AbstractFloat, 3},
    domain::Domain,
)
```

Compute dot product of two 3D arrays distributed across MPI processes.

# Arguments

  - `a, b::AbstractArray{<:AbstractFloat, 3}`: Input arrays (must have identical shapes)
  - `domain::Domain`: MPI domain decomposition information

# Returns

  - `global_dot_product::AbstractFloat`: Global dot product ``∑ᵢ aᵢ·bᵢ`` across all MPI processes.

# Implementation

 1. Validates array shapes match
 2. Computes local dot product using `LinearAlgebra.dot`
 3. Reduces via `MPI.Allreduce` with sum operation

# Use Cases

  - Iterative solver convergence testing
  - Distributed array norm computation
  - Global scalar reductions
"""
function compute_global_dot_product(
    a::AbstractArray{<:AbstractFloat, 3},
    b::AbstractArray{<:AbstractFloat, 3},
    domain::Domain,
)

    # Get parameters.
    (; comm) = domain

    # Get shapes.
    asize = size(a)
    bsize = size(b)

    # Check if shapes agree.
    for i in 1:3
        if asize[i] != bsize[i]
            error("Error in compute_global_dot_product: Shapes disagree!")
        end
    end

    # Compute local dot product.
    local_dot_product = 0.0
    for k in 1:asize[3], j in 1:asize[2]
        @views local_dot_product += dot(a[:, j, k], b[:, j, k])
    end

    # Sum over all processes.
    global_dot_product = MPI.Allreduce(local_dot_product, +, comm)

    # Return the result.
    return global_dot_product
end
