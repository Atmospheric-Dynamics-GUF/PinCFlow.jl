function apply_3d_muscl!(
  phi::OffsetArray{<:AbstractFloat, 3},
  phitilde::OffsetArray{<:AbstractFloat, 5},
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

  for kz in 2:(sizez - 1)
    for jy in 2:(sizey - 1)
      phix[:] .= phi.parent[:, jy, kz]
      muscle_reconstruct1D_mcvariant!(phix, phitildex, sizex)
      phitilde.parent[:, jy, kz, 1, :] .= phitildex.parent
    end
  end

  for kz in 2:(sizez - 1)
    for ix in 2:(sizex - 1)
      phiy .= phi.parent[ix, :, kz]
      muscle_reconstruct1D_mcvariant!(phiY, phitildey, sizey)
      phitilde.parent[ix, :, kz, 2, :] .= phitildey.parent
    end
  end

  for jy in 2:(sizey - 1)
    for ix in 2:(sizex - 1)
      phiz .= phi.parent[ix, jy, :]
      muscle_reconstruct1D_mcvariant!(phiz, phitildez, sizez)
      phitilde.parent[ix, jy, :, 3, :] .= phitildez.parent
    end
  end

  return
end
