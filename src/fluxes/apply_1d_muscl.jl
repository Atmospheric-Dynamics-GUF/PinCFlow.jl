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
