function set_compressible_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    model::AbstractModel,
)
    return
end

function set_compressible_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Compressible,
)
    (; nbz) = state.namelists.domain
    (; k0, k1) = state.domain
    (; p) = state.variables.predictands
    for k in 1:nbz
        @views p[:, :, k0 - k] .= p[:, :, k0 + k - 1]
        @views p[:, :, k1 + k] .= p[:, :, k1 - k + 1]
    end
    return
end

function set_compressible_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    model::Compressible,
)
    (; k0, k1) = state.domain
    (; phip) = state.variables.fluxes
    phip[:, :, k0 - 1, 3] .= 0.0
    phip[:, :, k1, 3] .= 0.0
    return
end