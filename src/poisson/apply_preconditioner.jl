function apply_preconditioner!(
  sin::AbstractArray{<:AbstractFloat, 3},
  sout::AbstractArray{<:AbstractFloat, 3},
  namelists::Namelists,
  domain::Domain,
  grid::Grid,
  poisson::Poisson,
)
  (; dtau, maxiteradi) = namelists.poisson
  (; nx, ny, nz) = domain
  (; dx, dy) = grid
  (; au_b, ac_b, ad_b) = poisson.tensor
  (; s_pc, q_pc, p_pc) = poisson.preconditioner

  # Initialize auxiliary fields.
  s_pc .= 0.0
  q_pc .= 0.0
  p_pc .= 0.0

  # Set pseudo-time step.
  deta = dtau / (2.0 * (1.0 / dx^2 + 1.0 / dy^2))

  # Iterate.
  for niter in 1:maxiteradi
    if niter == 0
      s_pc .= sin
    else
      # Treat all diagonal elements implicitly.
      apply_operator!(s_pc, q_pc, Horizontal(), namelists, domain, poisson)
      s_pc .= s_pc .+ deta .* (q_pc .- sin)
    end

    # Perform upward sweep.
    for j in 1:ny
      for i in 1:nx
        au_b[i, j, nz] = 0.0
      end
    end

    for j in 1:ny, i in 1:nx
      if niter == 0
        q_pc[i, j, 1] = -au_b[i, j, 1] / ac_b[i, j, 1]
        s_pc[i, j, 1] = s_pc[i, j, 1] / ac_b[i, j, 1]
      else
        # Treat all diagonal elements implicity.
        q_pc[i, j, 1] = deta * au_b[i, j, 1] / (1.0 - deta * ac_b[i, j, 1])
        s_pc[i, j, 1] = s_pc[i, j, 1] / (1.0 - deta * ac_b[i, j, 1])
      end
    end

    for k in 2:nz, j in 1:ny, i in 1:nx
      if niter == 0
        p_pc[i, j] = 1.0 / (ac_b[i, j, k] + ad_b[i, j, k] * q_pc[i, j, k - 1])

        q_pc[i, j, k] = -au_b[i, j, k] * p_pc[i, j]

        s_pc[i, j, k] =
          (s_pc[i, j, k] - ad_b[i, j, k] * s_pc[i, j, k - 1]) * p_pc[i, j]
      else
        # Treat all diagonal elements implicitly.
        p_pc[i, j] =
          1.0 / (
            1.0 - deta * ac_b[i, j, k] -
            deta * ad_b[i, j, k] * q_pc[i, j, k - 1]
          )

        q_pc[i, j, k] = deta * au_b[i, j, k] * p_pc[i, j]

        s_pc[i, j, k] =
          (s_pc[i, j, k] + deta * ad_b[i, j, k] * s_pc[i, j, k - 1]) *
          p_pc[i, j]
      end
    end

    # Perform backward pass.
    for k in (nz - 1):-1:1, j in 1:ny, i in 1:nx
      s_pc[i, j, k] = s_pc[i, j, k] + q_pc[i, j, k] * s_pc[i, j, k + 1]
    end
  end

  # Set final result.
  sout .= s_pc

  # Return.
  return
end
