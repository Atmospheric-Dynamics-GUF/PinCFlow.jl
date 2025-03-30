struct Predictands{A <: AbstractArray{<:AbstractFloat, 3}}
  rho::A
  rhop::A
  u::A
  v::A
  w::A
  pip::A
end

function Predictands(namelists::Namelists, constants::Constants, domain::Domain)
  (; model, testcase) = namelists.setting
  return Predictands(namelists, constants, domain, model, testcase)
end

function Predictands(
  namelists::Namelists,
  constants::Constants,
  domain::Domain,
  model::PseudoIncompressible,
  testcase::MountainWave,
)

  # Get parameters.
  (; backgroundflow_dim) = namelists.atmosphere
  (; uref) = constants
  (; nxx, nyy, nzz) = domain

  # Initialize the predictands.
  (rho, rhop, u, v, w, pip) = (zeros(nxx, nyy, nzz) for i in 1:6)

  # Initial winds.
  u .= backgroundflow_dim[1] / uref
  v .= backgroundflow_dim[2] / uref
  w .= backgroundflow_dim[3] / uref

  # Return a Predictands instance.
  return Predictands(rho, rhop, u, v, w, pip)
end
