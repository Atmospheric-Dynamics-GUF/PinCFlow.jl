struct BicGStab{
  A <: AbstractMatrix{<:AbstractFloat},
  B <: AbstractArray{<:AbstractFloat, 3},
}
  r_vm::A
  p::B
  r0::B
  rold::B
  r::B
  s::B
  t::B
  v::B
  matvec::B
  v_pc::B
end

function BicGStab(domain::Domain)

  # Get all necessary fields.
  (; nx, ny, nz) = domain

  # Initialize BicGStab fields.
  r_vm = zeros(nx, ny)
  (p, r0, rold, r, s, t, v, matvec, v_pc) = (zeros(nx, ny, nz) for i in 1:9)

  # Return a BicGStab instance.
  return BicGStab(r_vm, p, r0, rold, r, s, t, v, matvec, v_pc)
end
