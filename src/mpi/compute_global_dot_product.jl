function compute_global_dot_product(
  a::AbstractArray{AbstractFloat, 3},
  b::AbstractArray{AbstractFloat, 3},
  domain::Domain,
)

  # Get parameters.
  (; comm, root) = domain

  # Get shapes.
  asize = size(a)
  bsize = size(b)

  # Check if shapes agree.
  for i in 1:3
    if asize[i] != bsize[i]
      println("Error in global_dot_product: Shapes disagree!")
      exit()
    end
  end

  # Compute local dot product.
  local_dot_product = 0.0
  for k in 1:asize[3]
    for j in 1:asize[2]
      local_dot_product = local_dot_product + dot(a[:, j, k], b[:, j, k])
    end
  end

  # Sum over all processes.
  global_dot_product = MPI.Allreduce(local_dot_product, +, comm)

  # Return the result.
  return global_dot_product
end
