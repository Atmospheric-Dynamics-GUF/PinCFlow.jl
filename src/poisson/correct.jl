function correct!(
  state::State,
  dt::AbstractFloat,
  opt::EXPL,
  facray::AbstractFloat,
  facprs::AbstractFloat,
)
  correct!(state, dt, U())
  correct!(state, dt, V())
  correct!(state, dt, W())
  correct!(state, PiP())
  return
end

function correct!(
  state::State,
  dt::AbstractFloat,
  opt::IMPL,
  facray::AbstractFloat,
  facprs::AbstractFloat,
)
  correct!(state, dt, U(), facray, facprs)
  correct!(state, dt, V(), facray, facprs)
  correct!(state, dt, W(), facray, facprs)
  correct!(state, dt, RhoP(), facray, facprs)
  correct!(state, PiP())
  return
end

function correct!(state::State, dt::AbstractFloat, variable::U)
  (; zboundaries) = state.namelists.boundaries
  (; kappainv, mainv2) = state.constants
  (; nx, ny, nz) = state.domain
  (; dx, dz, met) = state.grid
  (; rhostrattfc, pstrattfc) = state.atmosphere
  (; dpip) = state.variables.tendencies
  (; rho, u) = state.variables.predictands

  for k in 1:nz, j in 1:ny, i in 0:nx
    # Compute values at cell edges.
    rhou =
      0.5 * (
        rho[i, j, k] +
        rho[i + 1, j, k] +
        rhostrattfc[i, j, k] +
        rhostrattfc[i + 1, j, k]
      )
    pedger = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i + 1, j, k])
    met13edger = 0.5 * (met[i, j, k, 1, 3] + met[i + 1, j, k, 1, 3])

    # Compute pressure difference gradient component.
    if k == 1 && zboundaries == SolidWallBoundaries()
      dpuuedger = 0.5 * (dpip[i, j, k + 2] + dpip[i + 1, j, k + 2])
      dpuedger = 0.5 * (dpip[i, j, k + 1] + dpip[i + 1, j, k + 1])
      dpedger = 0.5 * (dpip[i, j, k] + dpip[i + 1, j, k])
      pgradx =
        kappainv * mainv2 / rhou *
        pedger *
        (
          (dpip[i + 1, j, k] - dpip[i, j, k]) / dx +
          met13edger * (-dpuuedger + 4.0 * dpuedger - 3.0 * dpedger) * 0.5 / dz
        )
    elseif k == nz && zboundaries == SolidWallBoundaries()
      dpddedger = 0.5 * (dpip[i, j, k - 2] + dpip[i + 1, j, k - 2])
      dpdedger = 0.5 * (dpip[i, j, k - 1] + dpip[i + 1, j, k - 1])
      dpedger = 0.5 * (dpip[i, j, k] + dpip[i + 1, j, k])
      pgradx =
        kappainv * mainv2 / rhou *
        pedger *
        (
          (dpip[i + 1, j, k] - dpip[i, j, k]) / dx +
          met13edger * (dpddedger - 4.0 * dpdedger + 3.0 * dpedger) * 0.5 / dz
        )
    else
      dpuedger = 0.5 * (dpip[i, j, k + 1] + dpip[i + 1, j, k + 1])
      dpdedger = 0.5 * (dpip[i, j, k - 1] + dpip[i + 1, j, k - 1])
      pgradx =
        kappainv * mainv2 / rhou *
        pedger *
        (
          (dpip[i + 1, j, k] - dpip[i, j, k]) / dx +
          met13edger * (dpuedger - dpdedger) * 0.5 / dz
        )
    end

    du = -dt * pgradx

    u[i, j, k] += du
  end

  return
end

function correct!(
  state::State,
  dt::AbstractFloat,
  variable::U,
  facray::AbstractFloat,
  facprs::AbstractFloat,
)
  (; spongelayer, sponge_uv) = state.namelists.sponge
  (; zboundaries) = state.namelists.boundaries
  (; kappainv, mainv2) = state.constants
  (; nx, ny, nz) = state.domain
  (; dx, dz, met) = state.grid
  (; rhostrattfc, pstrattfc) = state.atmosphere
  (; kr_sp_tfc) = state.sponge
  (; corx) = state.poisson.correction
  (; dpip) = state.variables.tendencies
  (; rho, u) = state.variables.predictands

  for k in 1:nz, j in 1:ny, i in 0:nx
    facu = 1.0

    if spongelayer && sponge_uv
      facu += dt * 0.5 * (kr_sp_tfc[i, j, k] + kr_sp_tfc[i + 1, j, k]) * facray
    end

    # Compute values at cell edges.
    rhou =
      0.5 * (
        rho[i, j, k] +
        rho[i + 1, j, k] +
        rhostrattfc[i, j, k] +
        rhostrattfc[i + 1, j, k]
      )
    pedger = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i + 1, j, k])
    met13edger = 0.5 * (met[i, j, k, 1, 3] + met[i + 1, j, k, 1, 3])

    # Compute pressure difference gradient component.
    if k == 1 && zboundaries == SolidWallBoundaries()
      dpuuedger = 0.5 * (dpip[i, j, k + 2] + dpip[i + 1, j, k + 2])
      dpuedger = 0.5 * (dpip[i, j, k + 1] + dpip[i + 1, j, k + 1])
      dpedger = 0.5 * (dpip[i, j, k] + dpip[i + 1, j, k])
      pgradx =
        kappainv * mainv2 / rhou *
        pedger *
        (
          (dpip[i + 1, j, k] - dpip[i, j, k]) / dx +
          met13edger * (-dpuuedger + 4.0 * dpuedger - 3.0 * dpedger) * 0.5 / dz
        )
    elseif k == nz && zboundaries == SolidWallBoundaries()
      dpddedger = 0.5 * (dpip[i, j, k - 2] + dpip[i + 1, j, k - 2])
      dpdedger = 0.5 * (dpip[i, j, k - 1] + dpip[i + 1, j, k - 1])
      dpedger = 0.5 * (dpip[i, j, k] + dpip[i + 1, j, k])
      pgradx =
        kappainv * mainv2 / rhou *
        pedger *
        (
          (dpip[i + 1, j, k] - dpip[i, j, k]) / dx +
          met13edger * (dpddedger - 4.0 * dpdedger + 3.0 * dpedger) * 0.5 / dz
        )
    else
      dpuedger = 0.5 * (dpip[i, j, k + 1] + dpip[i + 1, j, k + 1])
      dpdedger = 0.5 * (dpip[i, j, k - 1] + dpip[i + 1, j, k - 1])
      pgradx =
        kappainv * mainv2 / rhou *
        pedger *
        (
          (dpip[i + 1, j, k] - dpip[i, j, k]) / dx +
          met13edger * (dpuedger - dpdedger) * 0.5 / dz
        )
    end

    # Compute velocity correction.
    corx[i, j, k] = facprs * dt / facu * pgradx
    du = -corx[i, j, k]

    u[i, j, k] += du
  end

  return
end

function correct!(state::State, dt::AbstractFloat, variable::V)
  (; zboundaries) = state.namelists.boundaries
  (; kappainv, mainv2) = state.constants
  (; nx, ny, nz) = state.domain
  (; dy, dz, met) = state.grid
  (; rhostrattfc, pstrattfc) = state.atmosphere
  (; dpip) = state.variables.tendencies
  (; rho, v) = state.variables.predictands

  for k in 1:nz, j in 0:ny, i in 1:nx
    # Compute values at cell edges.
    rhov =
      0.5 * (
        rho[i, j, k] +
        rho[i, j + 1, k] +
        rhostrattfc[i, j, k] +
        rhostrattfc[i, j + 1, k]
      )
    pedgef = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j + 1, k])
    met23edgef = 0.5 * (met[i, j, k, 2, 3] + met[i, j + 1, k, 2, 3])

    # Compute pressure difference gradient component.
    if k == 1 && zboundaries == SolidWallBoundaries()
      dpuuedgef = 0.5 * (dpip[i, j, k + 2] + dpip[i, j + 1, k + 2])
      dpuedgef = 0.5 * (dpip[i, j, k + 1] + dpip[i, j + 1, k + 1])
      dpedgef = 0.5 * (dpip[i, j, k] + dpip[i, j + 1, k])
      pgrady =
        kappainv * mainv2 / rhov *
        pedgef *
        (
          (dpip[i, j + 1, k] - dpip[i, j, k]) / dy +
          met23edgef * (-dpuuedgef + 4.0 * dpuedgef - 3.0 * dpedgef) * 0.5 / dz
        )
    elseif k == nz && zboundaries == SolidWallBoundaries()
      dpddedgef = 0.5 * (dpip[i, j, k - 2] + dpip[i, j + 1, k - 2])
      dpdedgef = 0.5 * (dpip[i, j, k - 1] + dpip[i, j + 1, k - 1])
      dpedgef = 0.5 * (dpip[i, j, k] + dpip[i, j + 1, k])
      pgrady =
        kappainv * mainv2 / rhov *
        pedgef *
        (
          (dpip[i, j + 1, k] - dpip[i, j, k]) / dy +
          met23edgef * (dpddedgef - 4.0 * dpdedgef + 3.0 * dpedgef) * 0.5 / dz
        )
    else
      dpuedgef = 0.5 * (dpip[i, j, k + 1] + dpip[i, j + 1, k + 1])
      dpdedgef = 0.5 * (dpip[i, j, k - 1] + dpip[i, j + 1, k - 1])
      pgrady =
        kappainv * mainv2 / rhov *
        pedgef *
        (
          (dpip[i, j + 1, k] - dpip[i, j, k]) / dy +
          met23edgef * (dpuedgef - dpdedgef) * 0.5 / dz
        )
    end

    dv = -dt * pgrady

    v[i, j, k] += dv
  end

  return
end

function correct!(
  state::State,
  dt::AbstractFloat,
  variable::V,
  facray::AbstractFloat,
  facprs::AbstractFloat,
)
  (; spongelayer, sponge_uv) = state.namelists.sponge
  (; zboundaries) = state.namelists.boundaries
  (; kappainv, mainv2) = state.constants
  (; nx, ny, nz) = state.domain
  (; dy, dz, met) = state.grid
  (; rhostrattfc, pstrattfc) = state.atmosphere
  (; kr_sp_tfc) = state.sponge
  (; cory) = state.poisson.correction
  (; dpip) = state.variables.tendencies
  (; rho, v) = state.variables.predictands

  for k in 1:nz, j in 0:ny, i in 1:nx
    facv = 1.0

    if spongelayer && sponge_uv
      facv += dt * 0.5 * (kr_sp_tfc[i, j, k] + kr_sp_tfc[i, j + 1, k]) * facray
    end

    # Compute values at cell edges.
    rhov =
      0.5 * (
        rho[i, j, k] +
        rho[i, j + 1, k] +
        rhostrattfc[i, j, k] +
        rhostrattfc[i, j + 1, k]
      )
    pedgef = 0.5 * (pstrattfc[i, j, k] + pstrattfc[i, j + 1, k])
    met23edgef = 0.5 * (met[i, j, k, 2, 3] + met[i, j + 1, k, 2, 3])

    # Compute pressure difference gradient component.
    if k == 1 && zboundaries == SolidWallBoundaries()
      dpuuedgef = 0.5 * (dpip[i, j, k + 2] + dpip[i, j + 1, k + 2])
      dpuedgef = 0.5 * (dpip[i, j, k + 1] + dpip[i, j + 1, k + 1])
      dpedgef = 0.5 * (dpip[i, j, k] + dpip[i, j + 1, k])
      pgrady =
        kappainv * mainv2 / rhov *
        pedgef *
        (
          (dpip[i, j + 1, k] - dpip[i, j, k]) / dy +
          met23edgef * (-dpuuedgef + 4.0 * dpuedgef - 3.0 * dpedgef) * 0.5 / dz
        )
    elseif k == nz && zboundaries == SolidWallBoundaries()
      dpddedgef = 0.5 * (dpip[i, j, k - 2] + dpip[i, j + 1, k - 2])
      dpdedgef = 0.5 * (dpip[i, j, k - 1] + dpip[i, j + 1, k - 1])
      dpedgef = 0.5 * (dpip[i, j, k] + dpip[i, j + 1, k])
      pgrady =
        kappainv * mainv2 / rhov *
        pedgef *
        (
          (dpip[i, j + 1, k] - dpip[i, j, k]) / dy +
          met23edgef * (dpddedgef - 4.0 * dpdedgef + 3.0 * dpedgef) * 0.5 / dz
        )
    else
      dpuedgef = 0.5 * (dpip[i, j, k + 1] + dpip[i, j + 1, k + 1])
      dpdedgef = 0.5 * (dpip[i, j, k - 1] + dpip[i, j + 1, k - 1])
      pgrady =
        kappainv * mainv2 / rhov *
        pedgef *
        (
          (dpip[i, j + 1, k] - dpip[i, j, k]) / dy +
          met23edgef * (dpuedgef - dpdedgef) * 0.5 / dz
        )
    end

    # Compute velocity correction.
    cory[i, j, k] = facprs * dt / facv * pgrady
    dv = -cory[i, j, k]

    v[i, j, k] += dv
  end

  return
end

function correct!(state::State, dt::AbstractFloat, variable::W)
  (; zboundaries) = state.namelists.boundaries
  (; kappainv, mainv2) = state.constants
  (; nx, ny, nz) = state.domain
  (; dx, dy, dz, jac, met) = state.grid
  (; rhostrattfc, pstrattfc, bvsstrattfc) = state.atmosphere
  (; dpip) = state.variables.tendencies
  (; rho, w) = state.variables.predictands

  if zboundaries == SolidWallBoundaries()
    k0 = 1
    k1 = nz - 1
  else
    error("Error in correct!: Unknown zboundaries!")
  end

  for k in k0:k1, j in 1:ny, i in 1:nx
    # Compute values at cell edges.
    rhostratedgeu =
      (
        jac[i, j, k + 1] * rhostrattfc[i, j, k] +
        jac[i, j, k] * rhostrattfc[i, j, k + 1]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    rhoedge =
      (jac[i, j, k + 1] * rho[i, j, k] + jac[i, j, k] * rho[i, j, k + 1]) /
      (jac[i, j, k] + jac[i, j, k + 1]) + rhostratedgeu
    pedgeu =
      (
        jac[i, j, k + 1] * pstrattfc[i, j, k] +
        jac[i, j, k] * pstrattfc[i, j, k + 1]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    bvsstw =
      (
        jac[i, j, k + 1] * bvsstrattfc[i, j, k] +
        jac[i, j, k] * bvsstrattfc[i, j, k + 1]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met13edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 1, 3] +
        jac[i, j, k] * met[i, j, k + 1, 1, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met23edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 2, 3] +
        jac[i, j, k] * met[i, j, k + 1, 2, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met33edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 3, 3] +
        jac[i, j, k] * met[i, j, k + 1, 3, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    dpredgeu =
      (
        jac[i + 1, j, k + 1] * dpip[i + 1, j, k] +
        jac[i + 1, j, k] * dpip[i + 1, j, k + 1]
      ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
    dpledgeu =
      (
        jac[i - 1, j, k + 1] * dpip[i - 1, j, k] +
        jac[i - 1, j, k] * dpip[i - 1, j, k + 1]
      ) / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
    dpfedgeu =
      (
        jac[i, j + 1, k + 1] * dpip[i, j + 1, k] +
        jac[i, j + 1, k] * dpip[i, j + 1, k + 1]
      ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
    dpbedgeu =
      (
        jac[i, j - 1, k + 1] * dpip[i, j - 1, k] +
        jac[i, j - 1, k] * dpip[i, j - 1, k + 1]
      ) / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])

    # Compute pressure difference gradient component.
    pgradz =
      kappainv * mainv2 / rhoedge *
      pedgeu *
      (
        met13edgeu * (dpredgeu - dpledgeu) * 0.5 / dx +
        met23edgeu * (dpfedgeu - dpbedgeu) * 0.5 / dy +
        met33edgeu * (dpip[i, j, k + 1] - dpip[i, j, k]) / dz
      )

    # Correct vertical velocity.
    dw = -dt * pgradz

    w[i, j, k] += dw
  end

  return
end

function correct!(
  state::State,
  dt::AbstractFloat,
  variable::W,
  facray::AbstractFloat,
  facprs::AbstractFloat,
)
  (; spongelayer, sponge_uv) = state.namelists.sponge
  (; zboundaries) = state.namelists.boundaries
  (; kappainv, mainv2) = state.constants
  (; nx, ny, nz) = state.domain
  (; dx, dy, dz, jac, met) = state.grid
  (; rhostrattfc, pstrattfc, bvsstrattfc) = state.atmosphere
  (; kr_sp_w_tfc) = state.sponge
  (; corx, cory) = state.poisson.correction
  (; dpip) = state.variables.tendencies
  (; rho, w) = state.variables.predictands

  if zboundaries == SolidWallBoundaries()
    k0 = 1
    k1 = nz - 1
  else
    error("Error in correct!: Unknown zboundaries!")
  end

  for k in k0:k1, j in 1:ny, i in 1:nx
    facw = 1.0

    if spongelayer
      facw +=
        dt * (
          jac[i, j, k + 1] * kr_sp_w_tfc[i, j, k] +
          jac[i, j, k] * kr_sp_w_tfc[i, j, k + 1]
        ) / (jac[i, j, k] + jac[i, j, k + 1]) * facray
    end

    # Compute values at cell edges.
    rhostratedgeu =
      (
        jac[i, j, k + 1] * rhostrattfc[i, j, k] +
        jac[i, j, k] * rhostrattfc[i, j, k + 1]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    rhoedge =
      (jac[i, j, k + 1] * rho[i, j, k] + jac[i, j, k] * rho[i, j, k + 1]) /
      (jac[i, j, k] + jac[i, j, k + 1]) + rhostratedgeu
    pedgeu =
      (
        jac[i, j, k + 1] * pstrattfc[i, j, k] +
        jac[i, j, k] * pstrattfc[i, j, k + 1]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    bvsstw =
      (
        jac[i, j, k + 1] * bvsstrattfc[i, j, k] +
        jac[i, j, k] * bvsstrattfc[i, j, k + 1]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met13edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 1, 3] +
        jac[i, j, k] * met[i, j, k + 1, 1, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met23edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 2, 3] +
        jac[i, j, k] * met[i, j, k + 1, 2, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met33edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 3, 3] +
        jac[i, j, k] * met[i, j, k + 1, 3, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    dpredgeu =
      (
        jac[i + 1, j, k + 1] * dpip[i + 1, j, k] +
        jac[i + 1, j, k] * dpip[i + 1, j, k + 1]
      ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
    dpledgeu =
      (
        jac[i - 1, j, k + 1] * dpip[i - 1, j, k] +
        jac[i - 1, j, k] * dpip[i - 1, j, k + 1]
      ) / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
    dpfedgeu =
      (
        jac[i, j + 1, k + 1] * dpip[i, j + 1, k] +
        jac[i, j + 1, k] * dpip[i, j + 1, k + 1]
      ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
    dpbedgeu =
      (
        jac[i, j - 1, k + 1] * dpip[i, j - 1, k] +
        jac[i, j - 1, k] * dpip[i, j - 1, k + 1]
      ) / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])

    # Compute pressure difference gradient component.
    pgradz =
      kappainv * mainv2 / rhoedge *
      pedgeu *
      (
        met13edgeu * (dpredgeu - dpledgeu) * 0.5 / dx +
        met23edgeu * (dpfedgeu - dpbedgeu) * 0.5 / dy +
        met33edgeu * (dpip[i, j, k + 1] - dpip[i, j, k]) / dz
      )

    # Compute velocity correction.
    dw =
      -facprs * dt / (facw + rhostratedgeu / rhoedge * bvsstw * dt^2.0) *
      pgradz -
      1.0 / (facw + rhostratedgeu / rhoedge * bvsstw * dt^2.0) * rhostratedgeu /
      rhoedge *
      bvsstw *
      dt^2.0 *
      0.5 *
      (
        jac[i, j, k + 1] * (
          met[i, j, k, 1, 3] * (corx[i, j, k] + corx[i - 1, j, k]) +
          met[i, j, k, 2, 3] * (cory[i, j, k] + cory[i, j - 1, k])
        ) +
        jac[i, j, k] * (
          met[i, j, k + 1, 1, 3] * (corx[i, j, k + 1] + corx[i - 1, j, k + 1]) +
          met[i, j, k + 1, 2, 3] * (cory[i, j, k + 1] + cory[i, j - 1, k + 1])
        )
      ) / (jac[i, j, k] + jac[i, j, k + 1])

    w[i, j, k] += dw
  end

  return
end

function correct!(
  state::State,
  dt::AbstractFloat,
  variable::RhoP,
  facray::AbstractFloat,
  facprs::AbstractFloat,
)
  (; spongelayer, sponge_uv) = state.namelists.sponge
  (; zboundaries) = state.namelists.boundaries
  (; kappainv, mainv2, g_ndim) = state.constants
  (; nx, ny, nz) = state.domain
  (; dx, dy, dz, jac, met) = state.grid
  (; rhostrattfc, pstrattfc, bvsstrattfc) = state.atmosphere
  (; kr_sp_w_tfc) = state.sponge
  (; corx, cory) = state.poisson.correction
  (; dpip) = state.variables.tendencies
  (; rho, rhop, w) = state.variables.predictands

  for k in 1:nz, j in 1:ny, i in 1:nx
    facw = 1.0

    if spongelayer
      facw += dt * kr_sp_w_tfc[i, j, k] * facray
    end

    # Compute P coefficients.
    pedgeu =
      (
        jac[i, j, k + 1] * pstrattfc[i, j, k] +
        jac[i, j, k] * pstrattfc[i, j, k + 1]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    pedged =
      (
        jac[i, j, k - 1] * pstrattfc[i, j, k] +
        jac[i, j, k] * pstrattfc[i, j, k - 1]
      ) / (jac[i, j, k] + jac[i, j, k - 1])

    # Compute density coefficients.
    rhow0 =
      (
        jac[i, j, k + 1] * (rho[i, j, k] + rhostrattfc[i, j, k]) +
        jac[i, j, k] * (rho[i, j, k + 1] + rhostrattfc[i, j, k + 1])
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    rhowm =
      (
        jac[i, j, k - 1] * (rho[i, j, k] + rhostrattfc[i, j, k]) +
        jac[i, j, k] * (rho[i, j, k - 1] + rhostrattfc[i, j, k - 1])
      ) / (jac[i, j, k] + jac[i, j, k - 1])

    # Interpolate metric tensor elements.
    met13edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 1, 3] +
        jac[i, j, k] * met[i, j, k + 1, 1, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met13edged =
      (
        jac[i, j, k - 1] * met[i, j, k, 1, 3] +
        jac[i, j, k] * met[i, j, k - 1, 1, 3]
      ) / (jac[i, j, k] + jac[i, j, k - 1])
    met23edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 2, 3] +
        jac[i, j, k] * met[i, j, k + 1, 2, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met23edged =
      (
        jac[i, j, k - 1] * met[i, j, k, 2, 3] +
        jac[i, j, k] * met[i, j, k - 1, 2, 3]
      ) / (jac[i, j, k] + jac[i, j, k - 1])
    met33edgeu =
      (
        jac[i, j, k + 1] * met[i, j, k, 3, 3] +
        jac[i, j, k] * met[i, j, k + 1, 3, 3]
      ) / (jac[i, j, k] + jac[i, j, k + 1])
    met33edged =
      (
        jac[i, j, k - 1] * met[i, j, k, 3, 3] +
        jac[i, j, k] * met[i, j, k - 1, 3, 3]
      ) / (jac[i, j, k] + jac[i, j, k - 1])

    # Interpolate pressure differences.
    dpredgeu =
      (
        jac[i + 1, j, k + 1] * dpip[i + 1, j, k] +
        jac[i + 1, j, k] * dpip[i + 1, j, k + 1]
      ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
    dpledgeu =
      (
        jac[i - 1, j, k + 1] * dpip[i - 1, j, k] +
        jac[i - 1, j, k] * dpip[i - 1, j, k + 1]
      ) / (jac[i - 1, j, k] + jac[i - 1, j, k + 1])
    dpredged =
      (
        jac[i + 1, j, k - 1] * dpip[i + 1, j, k] +
        jac[i + 1, j, k] * dpip[i + 1, j, k - 1]
      ) / (jac[i + 1, j, k] + jac[i + 1, j, k - 1])
    dpledged =
      (
        jac[i - 1, j, k - 1] * dpip[i - 1, j, k] +
        jac[i - 1, j, k] * dpip[i - 1, j, k - 1]
      ) / (jac[i - 1, j, k] + jac[i - 1, j, k - 1])
    dpfedgeu =
      (
        jac[i, j + 1, k + 1] * dpip[i, j + 1, k] +
        jac[i, j + 1, k] * dpip[i, j + 1, k + 1]
      ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
    dpbedgeu =
      (
        jac[i, j - 1, k + 1] * dpip[i, j - 1, k] +
        jac[i, j - 1, k] * dpip[i, j - 1, k + 1]
      ) / (jac[i, j - 1, k] + jac[i, j - 1, k + 1])
    dpfedged =
      (
        jac[i, j + 1, k - 1] * dpip[i, j + 1, k] +
        jac[i, j + 1, k] * dpip[i, j + 1, k - 1]
      ) / (jac[i, j + 1, k] + jac[i, j + 1, k - 1])
    dpbedged =
      (
        jac[i, j - 1, k - 1] * dpip[i, j - 1, k] +
        jac[i, j - 1, k] * dpip[i, j - 1, k - 1]
      ) / (jac[i, j - 1, k] + jac[i, j - 1, k - 1])

    # Compute pressure difference gradients.
    pgradzedgeu =
      kappainv * mainv2 * pedgeu / rhow0 * (
        0.5 * met13edgeu * (dpredgeu - dpledgeu) / dx +
        0.5 * met23edgeu * (dpfedgeu - dpbedgeu) / dy +
        met33edgeu * (dpip[i, j, k + 1] - dpip[i, j, k]) / dz
      )
    pgradzedged =
      kappainv * mainv2 * pedged / rhowm * (
        0.5 * met13edged * (dpredged - dpledged) / dx +
        0.5 * met23edged * (dpfedged - dpbedged) / dy +
        met33edged * (dpip[i, j, k] - dpip[i, j, k - 1]) / dz
      )

    # Adjust at boundaries.
    if k == 1 && zboundaries == SolidWallBoundaries()
      pgradzedged = 0.0
    elseif k == nz && zboundaries == SolidWallBoundaries()
      pgradzedgeu = 0.0
    end

    # Interpolate.
    pgradz = 0.5 * (pgradzedgeu + pgradzedged)

    # Compute buoyancy correction.
    db =
      -1.0 / (
        facw +
        rhostrattfc[i, j, k] / (rho[i, j, k] + rhostrattfc[i, j, k]) *
        bvsstrattfc[i, j, k] *
        dt^2.0
      ) * (
        -rhostrattfc[i, j, k] / (rho[i, j, k] + rhostrattfc[i, j, k]) *
        bvsstrattfc[i, j, k] *
        facprs *
        dt^2.0 *
        jac[i, j, k] *
        pgradz +
        rhostrattfc[i, j, k] / (rho[i, j, k] + rhostrattfc[i, j, k]) *
        bvsstrattfc[i, j, k] *
        dt *
        jac[i, j, k] *
        facw *
        0.5 *
        (
          met[i, j, k, 1, 3] * (corx[i, j, k] + corx[i - 1, j, k]) +
          met[i, j, k, 2, 3] * (cory[i, j, k] + cory[i, j - 1, k])
        )
      )

    rhop[i, j, k] -= (rho[i, j, k] + rhostrattfc[i, j, k]) / g_ndim * db
  end

  return
end

function correct!(state::State, variable::PiP)
  (; nx, ny, nz) = state.domain
  (; pip) = state.variables.predictands
  (; dpip) = state.variables.tendencies

  @views pip[0:(nx + 1), 0:(ny + 1), 0:(nz + 1)] .=
    pip[0:(nx + 1), 0:(ny + 1), 0:(nz + 1)] .+
    dpip[0:(nx + 1), 0:(ny + 1), 0:(nz + 1)]

  return
end
