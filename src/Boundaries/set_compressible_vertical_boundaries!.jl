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
    (; namelists, domain) = state
    (; zboundaries) = namelists.setting
    (; p) = state.variables.predictands

    set_vertical_boundaries_of_field!(p, namelists, domain, zboundaries, +)

    return
end

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
