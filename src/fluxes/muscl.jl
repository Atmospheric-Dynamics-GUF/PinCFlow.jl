function apply_1d_muscl!(
  phi::Vector{<:AbstractFloat},
  phitilde::Matrix{<:AbstractFloat},
  phisize::Integer,
)

  # Initialize phitilde.
  phitilde .= 1000.0

  for i in 2:(phisize - 1)
    deltal = phi[i] - phi[i - 1]
    deltar = phi[i + 1] - phi[i]

    if deltar == 0.0
      phitilde[i, 2] = phi[i]
      phitilde[i, 1] = phi[i]
    else
      if deltal == 0.0
        theta = deltal / deltar
        s = (2.0 + theta) / 3.0
        sigmal = max(0.0, min(2 * theta, s, 2.0))

        phitilde[i, 2] = phi[i] + 0.5 * sigmal * deltar
        phitilde[i, 1] = phi[i]
      else
        theta = deltal / deltar

        s = (2.0 + theta) / 3.0
        sigmal = max(0.0, min(2 * theta, s, 2.0))

        s = (2.0 + 1.0 / theta) / 3.0
        sigmar = max(0.0, min(2 / theta, s, 2.0))

        phitilde[i, 2] = phi[i] + 0.5 * sigmal * deltar
        phitilde[i, 1] = phi[i] - 0.5 * sigmar * deltal
      end
    end
  end

  return
end

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