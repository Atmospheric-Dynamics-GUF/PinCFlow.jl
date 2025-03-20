function solve_poisson!(
  state::State,
  b::AbstractArray{<:AbstractFloat, 3},
  tolref::AbstractFloat,
  dt::AbstractFloat,
  opt::AbstractIntegration,
  model::PseudoIncompressible,
  facray::AbstractFloat,
  facprs::AbstractFloat,
)
  (; namelists, domain, grid, poisson) = state
  (; nx, ny, nz) = domain
  (; rhostrattfc, pstrattfc) = state.atmosphere
  (; dpip) = state.variables.tendencies

  # Initialize solution and residual.
  sol = state.poisson.solution
  sol .= 0.0
  res = 0.0

  # Initialize.
  if dt == 0.0
    error("Error in solve_poisson!: dt = 0.0!")
  end
  dtinv = 1.0 / dt

  #--------------------------------
  #     Linear equation solver
  #     solve for dt * dp ...
  #--------------------------------

  compute_operator!(state, dt, opt, facray)

  (errflagbicg, niterbicg) = apply_bicgstab!(
    b,
    tolref,
    dt,
    sol,
    res,
    namelists,
    domain,
    grid,
    poisson,
  )

  if errflagbicg
    return
  end

  for k in 1:nz
    for j in 1:ny
      for i in 1:nx
        fcscal = sqrt(pstrattfc[i, j, k]^2 / rhostrattfc[i, j, k])
        sol[i, j, k] /= fcscal
      end
    end
  end

  # Pass solution to pressure correction.
  dpip[1:nx, 1:ny, 1:nz] = dtinv ./ facprs .* sol

  return (errflagbicg, niterbicg)
end
