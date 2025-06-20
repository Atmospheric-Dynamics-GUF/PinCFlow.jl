"""
    set_compressible_vertical_boundaries!(state, variables, model::AbstractModel)

No-op for non-compressible models.
"""
function set_compressible_vertical_boundaries!(
    state::State,
    variables::AbstractBoundaryVariables,
    model::AbstractModel,
)
    return
end

"""
    set_compressible_vertical_boundaries!(state, variables::BoundaryPredictands, model::Compressible)

Set vertical boundaries for pressure field in compressible model with symmetric conditions.
"""
function set_compressible_vertical_boundaries!(
    state::State,
    variables::BoundaryPredictands,
    model::Compressible,
)
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; p) = state.variables.predictands

    set_vertical_boundaries_of_field!(p, namelists, domain, zboundaries, +)

    return
end

"""
    set_compressible_vertical_boundaries!(state, variables::BoundaryFluxes, model::Compressible)

Set vertical pressure flux boundaries to zero at solid walls.
"""
function set_compressible_vertical_boundaries!(
    state::State,
    variables::BoundaryFluxes,
    model::Compressible,
)
    (; sizezz, nzz, ko, k0, k1)
    (; phip) = state.variables.fluxes

    if ko == 0
        phip[:, :, k0 - 1, 3] .= 0.0
    end

    if ko + nzz == sizezz
        phip[:, :, k1, 3] .= 0.0
    end

    return
end
