struct Correction{A <: AbstractArray{<:AbstractFloat, 3}}
  corx::A
  cory::A
end

function Correction(domain::Domain)

  # Get all necessary fields.
  (; nxx, nyy, nzz) = domain

  # Return a Correction instance.
  return Correction([zeros(nxx, nyy, nzz) for i in 1:2]...)
end
