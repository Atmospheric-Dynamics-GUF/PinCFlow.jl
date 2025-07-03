function conductive_heating(
    state::State,
    indices::NTuple{3, <:Integer},
)
    
    (; model) = state.namelists.setting

    return conductive_heating(state, indices, model)
end

function conductive_heating(
    state::State, 
    indices::NTuple{3, <:Integer},
    model::Boussinesq,
)

    return 0.0
end

function conductive_heating(
    state::State, 
    indices::NTuple{3, <:Integer},
    model::PseudoIncompressible,
)

    return 0.0
end

function conductive_heating(
    state::State, 
    indices::NTuple{3, <:Integer},
    model::Compressible,
)

    (; phitheta) = state.variables.fluxes
    (; rho) = state.variables.predictands
    (; jac, dx, dy, dz) = state.grid
    (; rhostrattfc) = state.atmosphere
    (i, j, k) = indices

    rhotot = (rho[i, j, k] + rhostrattfc[i, j, k]) / jac[i, j, k]

    heating = rhotot * ((phitheta[i, j, k, 1] - 
        phitheta[i - 1, j, k, 1]) / dx + (phitheta[i, j, k, 2] - 
        phitheta[i, j - 1, k, 2]) / dy + (phitheta[i, j, k, 3] - 
        phitheta[i, j, k - 1, 3]) / dz)


    return heating
end

