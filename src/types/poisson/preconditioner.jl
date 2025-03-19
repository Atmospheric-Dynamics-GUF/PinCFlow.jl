struct Preconditioner{A <: AbstractArray{<:AbstractFloat, 3}}
  s_pc::A
  q_pc::A
  p_pc::A
end

function Preconditioner(domain)

  # Get all necessary fields.
  (; nx, ny, nz) = domain

  # Initialize the preconditioner fields.
  (s_pc, q_pc, p_pc) = (zeros((nx, ny, nz)) for i in 1:3)

  # Return a Preconditioner instance.
  return Preconditioner(s_pc, q_pc, p_pc)
end
