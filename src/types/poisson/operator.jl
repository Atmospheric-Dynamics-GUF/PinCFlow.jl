struct Operator{A <: AbstractArray{<:AbstractFloat, 3}}
  s::A
end

function Operator(domain::Domain)

  # Get all necessary fields.
  (; nx, ny, nz) = domain

  # Initialize s.
  s = OffsetArray(
    zeros((nx + 2, ny + 2, nz + 2)),
    0:(nx + 1),
    0:(ny + 1),
    0:(nz + 1),
  )

  # Return an Operator instance.
  return Operator(s)
end
