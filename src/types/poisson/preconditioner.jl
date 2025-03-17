struct Preconditioner{A <: Array{<:AbstractFloat, 3}}
  s::A
  q::A
  p::A
end

function Preconditioner(domain)

  # Get all necessary fields.
  (; nx, ny, nz)

  # Initialize the preconditioner fields.
  (s, q, p) = (zeros((nx, ny, nz)) for i in 1:3)

  # Return a Preconditioner instance.
  return Preconditioner(s, q, p)
end
