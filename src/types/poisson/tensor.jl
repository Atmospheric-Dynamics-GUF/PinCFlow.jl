struct Tensor{A <: AbstractArray{<:AbstractFloat, 3}}
  ac_b::A
  acv_b::A
  ach_b::A
  al_b::A
  ar_b::A
  ab_b::A
  af_b::A
  ad_b::A
  au_b::A
  ald_b::A
  alu_b::A
  ard_b::A
  aru_b::A
  abd_b::A
  abu_b::A
  afd_b::A
  afu_b::A
  add_b::A
  auu_b::A
  aldd_b::A
  aluu_b::A
  ardd_b::A
  aruu_b::A
  abdd_b::A
  abuu_b::A
  afdd_b::A
  afuu_b::A
end

function Tensor(domain::Domain)

  # Get parameters.
  (; nx, ny, nz) = domain

  # Initialize the elements of the operator.
  (
    ac_b,
    acv_b,
    ach_b,
    al_b,
    ar_b,
    ab_b,
    af_b,
    ad_b,
    au_b,
    ald_b,
    alu_b,
    ard_b,
    aru_b,
    abd_b,
    abu_b,
    afd_b,
    afu_b,
    add_b,
    auu_b,
    aldd_b,
    aluu_b,
    ardd_b,
    aruu_b,
    abdd_b,
    abuu_b,
    afdd_b,
    afuu_b,
  ) = (zeros((nx, ny, nz)) for i in 1:27)

  # Return an Tensor instance.
  return Tensor(
    ac_b,
    acv_b,
    ach_b,
    al_b,
    ar_b,
    ab_b,
    af_b,
    ad_b,
    au_b,
    ald_b,
    alu_b,
    ard_b,
    aru_b,
    abd_b,
    abu_b,
    afd_b,
    afu_b,
    add_b,
    auu_b,
    aldd_b,
    aluu_b,
    ardd_b,
    aruu_b,
    abdd_b,
    abuu_b,
    afdd_b,
    afuu_b,
  )
end
