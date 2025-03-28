struct Preconditioner{
  A <: AbstractArray{<:AbstractFloat, 3},
  B <: AbstractMatrix{<:AbstractFloat},
}
  s_pc::A
  q_pc::A
  p_pc::B
end

function Preconditioner(domain::Domain)

  # Get all necessary fields.
  (; nx, ny, nz) = domain

  # Initialize the preconditioner fields.
  (s_pc, q_pc) = (zeros((nx, ny, nz)) for i in 1:2)
  p_pc = zeros(nx, ny)

  # Return a Preconditioner instance.
  return Preconditioner(s_pc, q_pc, p_pc)
end
