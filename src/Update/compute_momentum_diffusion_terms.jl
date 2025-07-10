function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    direction::X,
)
    (; u, v, w) = state.variables.predictands
    (; dx, dy, dz, jac, met) = state.grid

    uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1]) # u_{i, j, k+1}
    ud = 0.5 * (u[i, j, k - 1] + u[i - 1, j, k - 1]) # u_{i, j, k-1}

    diffux =
        (u[i, j, k] - u[i - i, j, k]) / dx +
        met[i, j, k, 1, 3] * (uu - ud) / (2.0 * dz)

    return diffux
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    direction::Y,
)
    (; u, v, w) = state.variables.predictands
    (; dx, dy, dz, jac, met) = state.grid

    uf = 0.5 * (u[i, j + 1, k] + u[i - 1, j + 1, k]) # u_{i, j+1, k}
    ub = 0.5 * (u[i, j - 1, k] + u[i - 1, j - 1, k]) # u_{i, j-1, k}

    uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1]) # u_{i, j, k+1}
    ud = 0.5 * (u[i, j, k - 1] + u[i - 1, j, k - 1]) # u_{i, j, k-1}

    diffuy =
        (uf - ub) / (2.0 * dy) + met[i, j, k, 2, 3] * (uu - ud) / (2.0 * dz)

    return diffuy
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::U,
    direction::Z,
)
    (; u, v, w) = state.variables.predictands
    (; dx, dy, dz, jac, met) = state.grid

    uf = 0.5 * (u[i, j + 1, k] + u[i - 1, j + 1, k]) # u_{i, j+1, k}
    ub = 0.5 * (u[i, j - 1, k] + u[i - 1, j - 1, k]) # u_{i, j-1, k}

    uu = 0.5 * (u[i, j, k + 1] + u[i - 1, j, k + 1]) # u_{i, j, k+1}
    ud = 0.5 * (u[i, j, k - 1] + u[i - 1, j, k - 1]) # u_{i, j, k-1}

    diffuz =
        met[i, j, k, 1, 3] * (u[i, j, k] - u[i - 1, j, k]) / dx +
        met[i, j, k, 2, 3] * (uf - ub) / (2.0 * dy) +
        met[i, j, k, 3, 3] * (uu - ud) / (2.0 * dz)

    return diffuz
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    direction::X,
)
    (; u, v, w) = state.variables.predictands
    (; dx, dy, dz, jac, met) = state.grid

    vr = 0.5 * (v[i + 1, j, k] + v[i + 1, j - 1, k]) # v_{i+1, j, k}
    vl = 0.5 * (v[i - 1, j, k] + v[i - 1, j - 1, k]) # v_{i-1, j, k}
    vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1]) # v_{i, j, k+1}
    vd = 0.5 * (v[i, j, k - 1] + v[i, j - 1, k - 1]) # v_{i, j, k-1}

    diffvx =
        (vr - vl) / (2.0 * dx) + met[i, j, k, 1, 3] * (vu - vd) / (2.0 * dz)

    return diffvx
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    direction::Y,
)
    (; u, v, w) = state.variables.predictands
    (; dx, dy, dz, jac, met) = state.grid

    vr = 0.5 * (v[i + 1, j, k] + v[i + 1, j - 1, k]) # v_{i+1, j, k}
    vl = 0.5 * (v[i - 1, j, k] + v[i - 1, j - 1, k]) # v_{i-1, j, k}
    vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1]) # v_{i, j, k+1}
    vd = 0.5 * (v[i, j, k - 1] + v[i, j - 1, k - 1]) # v_{i, j, k-1}

    diffvy =
        (v[i, j, k] - v[i, j - 1, k]) / dy +
        met[i, j, k, 2, 3] * (vu - vd) / (2.0 * dz)

    return diffvy
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::V,
    direction::Z,
)
    (; u, v, w) = state.variables.predictands
    (; dx, dy, dz, jac, met) = state.grid

    vr = 0.5 * (v[i + 1, j, k] + v[i + 1, j - 1, k]) # v_{i+1, j, k}
    vl = 0.5 * (v[i - 1, j, k] + v[i - 1, j - 1, k]) # v_{i-1, j, k}
    vu = 0.5 * (v[i, j, k + 1] + v[i, j - 1, k + 1]) # v_{i, j, k+1}
    vd = 0.5 * (v[i, j, k - 1] + v[i, j - 1, k - 1]) # v_{i, j, k-1}

    diffvz =
        met[i, j, k, 1, 3] * (vr - vl) / (2 * dx) +
        met[i, j, k, 2, 3] * (v[i, j, k] - v[i, j - 1, k]) / dy +
        met[i, j, k, 3, 3] * (vu - vd) / (2.0 * dz)

    return diffvz
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    direction::X,
)
    (; u, v, w) = state.variables.predictands
    (; dx, dy, dz, jac, met) = state.grid

    wr = 0.5 * (w[i + 1, j, k] + w[i + 1, j, k - 1])
    wl = 0.5 * (w[i - 1, j, k] + w[i - 1, j, k - 1])
    wf = 0.5 * (w[i, j + 1, k] + w[i, j + 1, k - 1])
    wb = 0.5 * (w[i, j - 1, k] + w[i, j - 1, k - 1])

    diffwx =
        (wr - wl) / (2.0 * dx) +
        met[i, j, k, 1, 3] * (w[i, j, k] - w[i, j, k - 1]) / dz

    return diffwx
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    direction::Y,
)
    (; u, v, w) = state.variables.predictands
    (; dx, dy, dz, jac, met) = state.grid

    wr = 0.5 * (w[i + 1, j, k] + w[i + 1, j, k - 1])
    wl = 0.5 * (w[i - 1, j, k] + w[i - 1, j, k - 1])
    wf = 0.5 * (w[i, j + 1, k] + w[i, j + 1, k - 1])
    wb = 0.5 * (w[i, j - 1, k] + w[i, j - 1, k - 1])

    diffwy =
        (wf - wb) / (2.0 * dy) +
        met[i, j, k, 2, 3] * (w[i, j, k] - w[i, j, k - 1]) / dz

    return diffwy
end

function compute_momentum_diffusion_terms(
    state::State,
    i::Integer,
    j::Integer,
    k::Integer,
    variable::W,
    direction::Z,
)
    (; u, v, w) = state.variables.predictands
    (; dx, dy, dz, jac, met) = state.grid

    wr = 0.5 * (w[i + 1, j, k] + w[i + 1, j, k - 1])
    wl = 0.5 * (w[i - 1, j, k] + w[i - 1, j, k - 1])
    wf = 0.5 * (w[i, j + 1, k] + w[i, j + 1, k - 1])
    wb = 0.5 * (w[i, j - 1, k] + w[i, j - 1, k - 1])

    diffwz =
        met[i, j, k, 1, 3] * (wr - wl) / (2.0 * dx) +
        met[i, j, k, 2, 3] * (wf - wb) / (2.0 * dy) +
        met[i, j, k, 3, 3] * (w[i, j, k] - w[i, j, k - 1]) / dz

    return diffwz
end
