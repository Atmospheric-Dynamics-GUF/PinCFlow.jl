struct Operator{A <: AbstractArray{<:AbstractFloat, 3}}
  s::A
end

function Operator(domain::Domain)

  # Get all necessary fields.
  (; nx, ny, nz) = domain

  # Initialize s.
  s = OffsetArray(
    zeros((nx + 2, ny + 2, nz + 4)),
    0:(nx + 1),
    0:(ny + 1),
    -1:(nz + 2),
  )

  # Return an Operator instance.
  return Operator(s)
end
