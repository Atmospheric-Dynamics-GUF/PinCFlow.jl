"""
```julia
compute_global_dot_product(
    a::AbstractArray{<:AbstractFloat, 3},
    b::AbstractArray{<:AbstractFloat, 3},
    state::State,
)::AbstractFloat
```

Compute and return the dot product ``\\boldsymbol{a} \\cdot \\boldsymbol{b} = \\sum_i a_i \\cdot b_i`` of two 3D arrays distributed across MPI processes.

# Arguments

  - `a`: First input array.

  - `b`: Second input array (must have the same shape as `a`).

  - `state`: Model state.
"""
function compute_global_dot_product end

function compute_global_dot_product(
    a::AbstractArray{<:AbstractFloat, 3},
    b::AbstractArray{<:AbstractFloat, 3},
    state::State,
)::AbstractFloat
    (; comm, nx, ny, nz) = state.domain

    # Compute local dot product.
    local_dot_product = 0.0
    @share (+, local_dot_product) for k in 1:nz, j in 1:ny, i in 1:nx
        local_dot_product += a[i, j, k] * b[i, j, k]
    end

    # Sum over all processes.
    global_dot_product = MPI.Allreduce(local_dot_product, +, comm)

    return global_dot_product
end
