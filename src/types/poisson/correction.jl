struct Correction{A <: AbstractArray{<:AbstractFloat, 3}}
  corx::A
  cory::A
end

function Correction(domain::Domain)

  # Get all necessary fields.
  (; nxx, nyy, nzz) = domain

  # Initialize correction fields.
  (corx, cory) = (zeros(nxx, nyy, nzz) for i in 1:2)

  # Return a Correction instance.
  return Correction(corx, cory)
end
