function apply_3d_muscl!(
  phi::AbstractArray{<:AbstractFloat, 3},
  phitilde::AbstractArray{<:AbstractFloat, 5},
  domain::Domain,
  auxiliaries::Auxiliaries,
  limitertype::MCVariant,
)
  (; nxx, nyy, nzz) = domain
  (; phix, phiy, phiz, phitildex, phitildey, phitildez) = auxiliaries

  for kz in 2:(nzz - 1), jy in 2:(nyy - 1)
    @views phix .= phi.parent[:, jy, kz]
    apply_1d_muscl!(phix, phitildex, nxx)
    phitilde.parent[:, jy, kz, 1, :] .= phitildex
  end

  for kz in 2:(nzz - 1), ix in 2:(nxx - 1)
    @views phiy .= phi.parent[ix, :, kz]
    apply_1d_muscl!(phiy, phitildey, nyy)
    phitilde.parent[ix, :, kz, 2, :] .= phitildey
  end

  for jy in 2:(nyy - 1), ix in 2:(nxx - 1)
    @views phiz .= phi.parent[ix, jy, :]
    apply_1d_muscl!(phiz, phitildez, nzz)
    phitilde.parent[ix, jy, :, 3, :] .= phitildez
  end

  return
end
