function apply_3d_muscl!(
  phi::AbstractArray{<:AbstractFloat, 3},
  phitilde::AbstractArray{<:AbstractFloat, 5},
  sizex::Integer,
  sizey::Integer,
  sizez::Integer,
  limitertype::MCVariant,
)
  phix = zeros(sizex)
  phiy = zeros(sizey)
  phiz = zeros(sizez)

  phitildex = zeros(sizex, 2)
  phitildey = zeros(sizey, 2)
  phitildez = zeros(sizez, 2)

  for kz in 2:(sizez - 1), jy in 2:(sizey - 1)
    phix .= phi.parent[:, jy, kz]
    apply_1d_muscl!(phix, phitildex, sizex)
    phitilde.parent[:, jy, kz, 1, :] .= phitildex
  end

  for kz in 2:(sizez - 1), ix in 2:(sizex - 1)
    phiy .= phi.parent[ix, :, kz]
    apply_1d_muscl!(phiy, phitildey, sizey)
    phitilde.parent[ix, :, kz, 2, :] .= phitildey
  end

  for jy in 2:(sizey - 1), ix in 2:(sizex - 1)
    phiz .= phi.parent[ix, jy, :]
    apply_1d_muscl!(phiz, phitildez, sizez)
    phitilde.parent[ix, jy, :, 3, :] .= phitildez
  end

  return
end
