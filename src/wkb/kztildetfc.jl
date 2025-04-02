function kztildetfc(
  i::Integer,
  j::Integer,
  z::AbstractFloat,
  domain::Domain,
  grid::Grid,
)

  # Get all necessary fields.
  (; nzz) = domain
  (; ztildetfc) = grid

  # Compute the vertical index.
  @views k = argmin(abs.(ztildetfc[i, j, :] .- z))
  if ztildetfc[i, j, k] < z
    k += 1
  end
  if k < 2
    k = 2
  end
  if k > nzz
    k = nzz
  end
  return k
end