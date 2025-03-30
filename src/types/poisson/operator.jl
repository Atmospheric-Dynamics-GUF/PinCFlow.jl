struct Operator{A <: AbstractArray{<:AbstractFloat, 3}}
  s::A
end

function Operator(domain::Domain)

  # Get all necessary fields.
  (; nxx, nyy, nzz) = domain

  # Initialize s.
  s = zeros(nxx, nyy, nzz)

  # Return an Operator instance.
  return Operator(s)
end
