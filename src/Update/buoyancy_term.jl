function buoyancy_term end

function buoyancy_term(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
)
    (; model) = state.namelists.atmosphere

    return buoyancy_term(state, p0, i, j, k, model)
end

function buoyancy_term(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
    model::Boussinesq,
)
    (; kh) = state.turbulence.turbulencediffusioncoefficients
    (; rhop) = p0
    (; pbar, rhobar) = state.atmosphere
    (; jac, dz) = state.grid
    (; g_ndim) = state.constants

    thetau = pbar[i, j, k + 1] / (rhop[i, j, k + 1] + rhobar[i, j, k + 1])
    theta = pbar[i, j, k] / (rhop[i, j, k] + rhobar[i, j, k])
    thetad = pbar[i, j, k - 1] / (rhop[i, j, k - 1] + rhobar[i, j, k - 1])

    wthetap = -kh[i, j, k] * (thetau - thetad) / (jac[i, j, k] * 2.0 * dz)

    buoyancy = wthetap * g_ndim / theta

    return buoyancy
end

function buoyancy_term(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
    model::PseudoIncompressible,
)
    (; kh) = state.turbulence.turbulencediffusioncoefficients
    (; rho) = p0
    (; pbar, rhobar) = state.atmosphere
    (; jac, dz) = state.grid
    (; g_ndim) = state.constants

    thetau = pbar[i, j, k + 1] / (rho[i, j, k + 1] + rhobar[i, j, k + 1])
    theta = pbar[i, j, k] / (rho[i, j, k] + rhobar[i, j, k])
    thetad = pbar[i, j, k - 1] / (rho[i, j, k - 1] + rhobar[i, j, k - 1])

    wthetap = -kh[i, j, k] * (thetau - thetad) / (jac[i, j, k] * 2.0 * dz)

    buoyancy = wthetap * g_ndim / theta

    return buoyancy
end

function buoyancy_term(
    state::State,
    p0::Predictands,
    i::Integer,
    j::Integer,
    k::Integer,
    model::Compressible,
)
    (; kh) = state.turbulence.turbulencediffusioncoefficients
    (; rho, p) = p0
    (; rhobar) = state.atmosphere
    (; jac, dz) = state.grid
    (; g_ndim) = state.constants

    thetau = p[i, j, k + 1] / (rho[i, j, k + 1] + rhobar[i, j, k + 1])
    theta = p[i, j, k] / (rho[i, j, k] + rhobar[i, j, k])
    thetad = p[i, j, k - 1] / (rho[i, j, k - 1] + rhobar[i, j, k - 1])

    wthetap = -kh[i, j, k] * (thetau - thetad) / (jac[i, j, k] * 2.0 * dz)

    buoyancy = wthetap * g_ndim / theta

    return buoyancy
end