struct Operator{A <: AbstractArray{<:AbstractFloat, 3}}
  s::A
end

function Operator(domain::Domain)

  # Get all necessary fields.
  (; nxx, nyy, nzz) = domain

  # Return an Operator instance.
  return Operator(zeros(nxx, nyy, nzz))
end
