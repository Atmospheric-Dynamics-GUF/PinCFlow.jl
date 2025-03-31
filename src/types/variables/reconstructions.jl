struct Reconstructions{A <: AbstractArray{<:AbstractFloat, 5}}
  rhotilde::A
  rhoptilde::A
  utilde::A
  vtilde::A
  wtilde::A
end

function Reconstructions(domain::Domain)

  # Get parameters.
  (; nxx, nyy, nzz) = domain

  # Initialize the reconstructed variables.
  (rhotilde, rhoptilde, utilde, vtilde, wtilde) =
    (zeros(nxx, nyy, nzz, 3, 2) for i in 1:5)

  # Return a Reconstructions instance.
  return Reconstructions(rhotilde, rhoptilde, utilde, vtilde, wtilde)
end
