function apply_3d_muscl!(
  phi::AbstractArray{<:AbstractFloat, 3},
  phitilde::AbstractArray{<:AbstractFloat, 5},
  nxx::Integer,
  nyy::Integer,
  nzz::Integer,
  limitertype::MCVariant,
)

  # Reconstruct in x.
  for kz in 2:(nzz - 1), jy in 2:(nyy - 1)
    @views apply_1d_muscl!(phi[:, jy, kz], phitilde[:, jy, kz, 1, :], nxx)
  end

  # Reconstruct in y.
  for kz in 2:(nzz - 1), ix in 2:(nxx - 1)
    @views apply_1d_muscl!(phi[ix, :, kz], phitilde[ix, :, kz, 2, :], nyy)
  end

  # Reconstruct in z.
  for jy in 2:(nyy - 1), ix in 2:(nxx - 1)
    @views apply_1d_muscl!(phi[ix, jy, :], phitilde[ix, jy, :, 3, :], nzz)
  end

  # Return.
  return
end
